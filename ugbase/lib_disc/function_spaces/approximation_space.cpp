/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include "approximation_space.h"
#include "lib_disc/domain.h"
#include "lib_disc/common/groups_util.h"
#include "common/util/string_util.h"
#include "common/profiler/profiler.h"

#ifdef UG_PARALLEL
	#include "pcl/pcl.h"
#endif

#include "lib_disc/dof_manager/dof_distribution.h"
#include "grid_function.h"

#include <algorithm> // std::sort
#include <sstream> // std::stringstream
using namespace std;

//	for debugging only:
//#include "lib_grid/file_io/file_io.h"
//#define APPROX_SPACE_PERFORM_CHANGED_GRID_DEBUG_SAVES
//#define APPROX_SPACE_PERFORM_DISTRIBUTED_GRID_DEBUG_SAVES

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// IApproximationSpace
////////////////////////////////////////////////////////////////////////////////

IApproximationSpace::
IApproximationSpace(SmartPtr<subset_handler_type> spMGSH,
                    SmartPtr<grid_type> spMG,
                    const AlgebraType& algebraType)
{
	init(spMGSH, spMG, algebraType);
}

IApproximationSpace::
IApproximationSpace(SmartPtr<subset_handler_type> spMGSH,
                    SmartPtr<grid_type> spMG)
{
	init(spMGSH, spMG, DefaultAlgebra::get());
}

void IApproximationSpace::
init(SmartPtr<subset_handler_type> spMGSH,
     SmartPtr<grid_type> spMG,
     const AlgebraType& algebraType)
{
	m_spMG = spMG;
	m_spMGSH = spMGSH;
	m_spDoFDistributionInfo = SmartPtr<DoFDistributionInfo>(new DoFDistributionInfo(spMGSH));
	m_algebraType = algebraType;
	m_bAdaptionIsActive = false;
	m_RevCnt = RevisionCounter(this);

	this->set_dof_distribution_info(m_spDoFDistributionInfo);

//	get blocksize of algebra
	const int blockSize = m_algebraType.blocksize();

//	a)	If blocksize fixed and > 1, we need grouping in dof manager. Thus,
//		the dofmanager hopefully fits (i.e. same number of grouped
//		dofs everywhere.)
	if(blockSize > 1) m_bGrouped = true;
//	b) 	If blocksize flexible, we group
	else if (blockSize == AlgebraType::VariableBlockSize) m_bGrouped = true;
//	c)	If blocksize == 1, we do not group. This will allow us to handle
//		this case for any problem.
	else if (blockSize == 1) m_bGrouped = false;
	else
		UG_THROW("Cannot determine blocksize of Algebra.");

//	this class listens to the grid-adaption-messages
	register_at_adaption_msg_hub();
}


IApproximationSpace::
~IApproximationSpace()
{
	if(m_spSurfaceView.valid())
		m_spSurfaceView = SmartPtr<SurfaceView>(nullptr);
}

template <typename TElem>
bool MightContainGhosts(const GridLayoutMap& layoutMap, int lvl)
{
	if(!layoutMap.has_layout<TElem>(INT_V_MASTER)) return false;
	if(lvl >= (int)layoutMap.get_layout<TElem>(INT_V_MASTER).num_levels()) return false;
	if(layoutMap.get_layout<TElem>(INT_V_MASTER).empty(lvl)) return false;

	typedef typename GridLayoutMap::Types<TElem>::Layout::LevelLayout TLayout;
	typedef typename TLayout::const_iterator InterfaceIterator;

	const TLayout& elemLayout = layoutMap.get_layout<TElem>(INT_V_MASTER).layout_on_level(lvl);
	for(InterfaceIterator iIter = elemLayout.begin(); iIter != elemLayout.end(); ++iIter){
		if(!elemLayout.interface(iIter).empty())
			return true;
	}

	return false;
}

bool IApproximationSpace::might_contain_ghosts(int lvl) const
{
	if(lvl < 0 || lvl > (int)num_levels()-1)
		UG_THROW("ApproximationSpace: Level not contained.");

	const DistributedGridManager* pDistGridMgr = m_spMG->distributed_grid_manager();
	if(!pDistGridMgr) return false;

	bool bGhosts = false;
	const GridLayoutMap& layoutMap = pDistGridMgr->grid_layout_map();
	if(max_dofs(VERTEX)) bGhosts |=  MightContainGhosts<Vertex>(layoutMap, lvl);
	if(max_dofs(EDGE)) bGhosts |=  MightContainGhosts<Edge>(layoutMap, lvl);
	if(max_dofs(FACE)) bGhosts |=  MightContainGhosts<Face>(layoutMap, lvl);
	if(max_dofs(VOLUME)) bGhosts |=  MightContainGhosts<Volume>(layoutMap, lvl);

	return bGhosts;
}

bool IApproximationSpace::might_contain_ghosts() const
{
	bool bGhosts = false;
	for(int lvl = 0; lvl < (int)num_levels(); ++lvl)
		bGhosts |= might_contain_ghosts(lvl);

	return bGhosts;
}

////////////////////////////////////////////////////////////////////////////////
// add
////////////////////////////////////////////////////////////////////////////////

void IApproximationSpace::
add(const std::vector<std::string>& vName, const char* fetype, int order)
{
	const int dim = DimensionOfSubsets(*m_spMGSH);
	if(dim == DIM_SUBSET_EMPTY_GRID)
		UG_THROW("ApproximationSpace: Cannot find dimension of grid. Maybe your grid is empty?");

	add(vName, ConvertStringToLFEID(fetype, dim, order));
}

void IApproximationSpace::
add(const std::vector<std::string>& vName, const char* fetype)
{
	const int dim = DimensionOfSubsets(*m_spMGSH);
	if(dim == DIM_SUBSET_EMPTY_GRID)
		UG_THROW("ApproximationSpace: Cannot find dimension of grid. Maybe your grid is empty?");

	add(vName, ConvertStringToLFEID(fetype, dim));
}

void IApproximationSpace::
add(const char* name, const char* fetype, int order)
{
	add(TokenizeTrimString(name), fetype, order);
}

void IApproximationSpace::
add(const char* name, const char* fetype)
{
	add(TokenizeTrimString(name), fetype);
}


void IApproximationSpace::
add(const std::vector<std::string>& vName, const char* fetype, int order,
    const std::vector<std::string>& vSubsets)
{
	SubsetGroup ssGrp(m_spMGSH, vSubsets);
	const int dim = ssGrp.get_highest_subset_dimension();

//	check
	if(dim == DIM_SUBSET_EMPTY_GRID)
		UG_THROW("ApproximationSpace: Cannot find dimension for new function on"
				"the subsets. Maybe your grid is empty?");

	add(vName, ConvertStringToLFEID(fetype, dim, order), vSubsets);
}

void IApproximationSpace::
add(const std::vector<std::string>& vName, const char* fetype,
    const std::vector<std::string>& vSubsets)
{
	SubsetGroup ssGrp(m_spMGSH, vSubsets);
	const int dim = ssGrp.get_highest_subset_dimension();

//	check
	if(dim == DIM_SUBSET_EMPTY_GRID)
		UG_THROW("ApproximationSpace: Cannot find dimension for new function on"
				"the subsets. Maybe your grid is empty?");

	add(vName, ConvertStringToLFEID(fetype, dim), vSubsets);
}

void IApproximationSpace::
add(const char* name, const char* fetype, int order, const char* subsets)
{
	add(TokenizeTrimString(name), fetype, order, TokenizeTrimString(subsets));
}

void IApproximationSpace::
add(const char* name, const char* fetype, const char* subsets)
{
	add(TokenizeTrimString(name), fetype, TokenizeTrimString(subsets));
}

////////////////////////////////////////////////////////////////////////////////
// DoFDistributions
////////////////////////////////////////////////////////////////////////////////

SmartPtr<DoFDistribution>
IApproximationSpace::dof_distribution(const GridLevel& gl, bool bCreate)
{
	for(size_t i = 0; i < m_vDD.size(); ++i)
		if(m_vDD[i]->grid_level() == gl)
			return m_vDD[i];

	if(!bCreate)
		UG_THROW("ApproxSpace: Could not create the DoFDistribution to GridLevel "<<gl);

	create_dof_distribution(gl);

	return dof_distribution(gl, false);
}

SmartPtr<DoFDistribution>
IApproximationSpace::dd(const GridLevel& gl, bool bCreate)
{
	return dof_distribution(gl, bCreate);
}

ConstSmartPtr<DoFDistribution>
IApproximationSpace::dof_distribution(const GridLevel& gl, bool bCreate) const
{
	return const_cast<IApproximationSpace*>(this)->dof_distribution(gl, bCreate);
}

ConstSmartPtr<DoFDistribution>
IApproximationSpace::dd(const GridLevel& gl, bool bCreate) const
{
	return dof_distribution(gl, bCreate);
}

std::vector<SmartPtr<DoFDistribution> >
IApproximationSpace::dof_distributions() const
{
	return m_vDD;
}


void IApproximationSpace::init_levels()
{
	PROFILE_FUNC();
	for(size_t lvl = 0; lvl < num_levels(); ++lvl){
		dof_distribution(GridLevel(lvl, GridLevel::LEVEL, false));
		dof_distribution(GridLevel(lvl, GridLevel::LEVEL, true));
	}
}

void IApproximationSpace::init_surfaces()
{
	PROFILE_FUNC();
	for(size_t lvl = 0; lvl < num_levels(); ++lvl)
		dof_distribution(GridLevel(lvl, GridLevel::SURFACE, false));

	init_top_surface();
}

void IApproximationSpace::init_top_surface()
{
	PROFILE_FUNC();
	dof_distribution(GridLevel(GridLevel::TOP, GridLevel::SURFACE, false));
}

////////////////////////////////////////////////////////////////////////////////
// DoFDistribution Creation
////////////////////////////////////////////////////////////////////////////////

bool SortDD(SmartPtr<DoFDistribution> spDD1, SmartPtr<DoFDistribution> spDD2){
	return spDD1->grid_level() < spDD2->grid_level();
}

void IApproximationSpace::create_dof_distribution(const GridLevel& gl)
{

	dof_distribution_info_required();
	surface_view_required();

//	get DoFIndexStorage if it is reusable
	SmartPtr<DoFIndexStorage> spIndexStrg;
	if(gl.is_level()){
		if(gl.ghosts()){
			if(m_spDoFIndexStrgForLevelWithGhost.invalid())
				m_spDoFIndexStrgForLevelWithGhost = SmartPtr<DoFIndexStorage>(
						new DoFIndexStorage(m_spMG, m_spDoFDistributionInfo));
			spIndexStrg = m_spDoFIndexStrgForLevelWithGhost;
		}
		else{
			if(m_spDoFIndexStrgForLevelNoGhost.invalid())
				m_spDoFIndexStrgForLevelNoGhost = SmartPtr<DoFIndexStorage>(
						new DoFIndexStorage(m_spMG, m_spDoFDistributionInfo));
			spIndexStrg = m_spDoFIndexStrgForLevelNoGhost;
		}
	}

//	create DoFDistribution
	SmartPtr<DoFDistribution> spDD = SmartPtr<DoFDistribution>(new
		DoFDistribution(m_spMG, m_spMGSH, m_spDoFDistributionInfo,
						m_spSurfaceView, gl, m_bGrouped, spIndexStrg));

//	add to list and sort
	m_vDD.push_back(spDD);
	std::sort(m_vDD.begin(), m_vDD.end(), SortDD);
}

void IApproximationSpace::surface_view_required()
{
//	allocate surface view if needed
	if(!m_spSurfaceView.valid())
		m_spSurfaceView = SmartPtr<SurfaceView>(new SurfaceView(m_spMGSH));
}

void IApproximationSpace::dof_distribution_info_required()
{
//	init dd-info (and fix the function pattern by that)
	m_spDoFDistributionInfo->init();

//	check that used algebra-type matches requirements
//	get blocksize of algebra
	const int blockSize = m_algebraType.blocksize();

//	if blockSize is 1, we're fine if dd is non-grouped
	if(blockSize == 1){
		if(m_bGrouped == true)
			UG_THROW("ApproximationSpace: Using grouped DD, but Algebra is 1x1.")
	}

//	if variable block algebra
	else if(blockSize == AlgebraType::VariableBlockSize){
		UG_THROW("ApproximationSpace: Variable algebra currently not supported.")
	}

//	if block algebra, check that number of sub-elements is zero or == blockSize
	else if(blockSize > 1){
		for(int r = 0; r < NUM_REFERENCE_OBJECTS; ++r){
			const ReferenceObjectID roid = (ReferenceObjectID)r;

			for(int si = 0; si < m_spDDI->num_subsets(); ++si){
				const int  numDoFs = m_spDDI->num_dofs(roid, si);

				if(numDoFs != 0 && numDoFs != blockSize)
					UG_THROW("ApproximationSpace: Using Block-Algebra with "
							"Blocksize "<<blockSize<<". Therefore, the number of"
							" dofs on each ReferenceObject must equal the blocksize"
							" or be zero. But number of dofs on "<<roid<<" in "
							"subset "<<si<<" is "<<numDoFs<<".");
			}
		}
	}

//	catch other (invalid) settings
	else
		UG_THROW("Cannot determine blocksize of Algebra.");
}

////////////////////////////////////////////////////////////////////////////////
// Grid-Change Handling
////////////////////////////////////////////////////////////////////////////////

void IApproximationSpace::reinit()
{
	PROFILE_FUNC();
//	update surface view
	if(m_spSurfaceView.valid())
		m_spSurfaceView->refresh_surface_states();

//	reinit all existing dof distributions
	for(size_t i = 0; i < m_vDD.size(); ++i){
		m_vDD[i]->reinit();
	}

//	increase revision counter
	++m_RevCnt;
}

void IApproximationSpace::register_at_adaption_msg_hub()
{
//	register function for grid adaption
	SPMessageHub msgHub = m_spMGSH->multi_grid()->message_hub();
	m_spGridAdaptionCallbackID =
		msgHub->register_class_callback(this,
		&ug::IApproximationSpace::grid_changed_callback);

	m_spGridDistributionCallbackID =
		msgHub->register_class_callback(this,
		&ug::IApproximationSpace::grid_distribution_callback);
}

void IApproximationSpace::
grid_changed_callback(const GridMessage_Adaption& msg)
{
	if(msg.adaption_begins())
		m_bAdaptionIsActive = true;

	else if(m_bAdaptionIsActive){
			if(msg.adaption_ends())
			{
				reinit();
				m_bAdaptionIsActive = false;

				#ifdef APPROX_SPACE_PERFORM_CHANGED_GRID_DEBUG_SAVES
					{
						static int counter = 0;
						std::stringstream ss;
						ss << "grid-changed-surface-view" << counter << "-p" << pcl::ProcRank() << ".ugx";
						UG_LOG("PERFORMING SURFACE VIEW DEBUG SAVE IN IApproximationSpace::grid_changed_callback: " << ss.str() << "\n");
						SaveSurfaceViewTransformed(*m_spMG, *m_spSurfaceView, ss.str().c_str(), 0.1);
						++counter;
					}
					{
						#ifdef UG_PARALLEL
							static int counter = 0;
							std::stringstream ss;
							ss << "grid-changed-parallel-layout-" << counter << "-p" << pcl::ProcRank() << ".ugx";
							UG_LOG("PERFORMING GRID LAYOUT DEBUG SAVE IN IApproximationSpace::grid_changed_callback: " << ss.str() << "\n");
							SaveParallelGridLayout(*m_spMG, ss.str().c_str(), 0.1);
							++counter;
						#endif
					}
				#endif
			}
	}

	else{
		UG_THROW("Before any grid-adaption may be performed, the approximation"
				" space has to be informed that grid-adaption shall begin. "
				"You may use IRefiner::grid_adaption_begins() or schedule "
				"an appropriate message to the associated grids message-hub.");
	}
}

void IApproximationSpace::
grid_distribution_callback(const GridMessage_Distribution& msg)
{
	PROFILE_FUNC();
	switch(msg.msg()){
		case GMDT_DISTRIBUTION_STARTS:
			break;

		case GMDT_DISTRIBUTION_STOPS:
			reinit();
			#ifdef APPROX_SPACE_PERFORM_DISTRIBUTED_GRID_DEBUG_SAVES
				{
					static int counter = 0;
					std::stringstream ss;
					ss << "grid-distributed-surface-view" << counter << "-p" << pcl::ProcRank() << ".ugx";
					UG_LOG("PERFORMING SURFACE VIEW DEBUG SAVE IN IApproximationSpace::grid_distribution_callback: " << ss.str() << "\n");
					SaveSurfaceViewTransformed(*m_spMG, *m_spSurfaceView, ss.str().c_str(), 0.1);
					++counter;
				}
				{
					#ifdef UG_PARALLEL
						static int counter = 0;
						std::stringstream ss;
						ss << "grid-distributed-parallel-layout-" << counter << "-p" << pcl::ProcRank() << ".ugx";
						UG_LOG("PERFORMING GRID LAYOUT DEBUG SAVE IN IApproximationSpace::grid_distribution_callback: " << ss.str() << "\n");
						SaveParallelGridLayout(*m_spMG, ss.str().c_str(), 0.1);
						++counter;
					#endif
				}
			#endif
			break;

		default:
			break;
	}
}

////////////////////////////////////////////////////////////////////////////////
// Statistic
////////////////////////////////////////////////////////////////////////////////

void PrintDoFCount(const vector<DoFCount>& vDC,
                   const string& sInfo,
                   const string& sAlgebra,
                   const string& sflags)
{
	const bool bPrintCmps = (sflags.find("component") != string::npos);
	const bool bPrintInterface = (sflags.find("interface") != string::npos);
	const bool bPrintSurface = (sflags.find("surface") != string::npos);
	const bool bPrintSubset = (sflags.find("subset") != string::npos);

//	check for output
	if(vDC.size() == 0)
		UG_THROW("Expected something to print.")

//	constants for size of output
	static const int LEVEL = 14;
	static const int COMPONENT = 7;
	static const int INTERFACE = 8;
	static const int SURFACE = 12;
	static const int NUMBER = 12;
	static const char* sSep = " | ";
	static const char* sLeft = "| ";
	static const char* sRight = " |";

//	constants for selection
	static const int ALL_FCT = DoFCount::ALL_FCT;
	static const int ALL_SUBSET = DoFCount::ALL_SUBSET;
	static const byte ALL_ES = DoFCount::ALL_ES;
	static const byte ALL_SS = DoFCount::ALL_SS;
	static const byte UNIQUE_ES = DoFCount::UNIQUE_ES;
	static const byte UNIQUE_SS = DoFCount::UNIQUE_SS;

//	Table Header
	stringstream ssHead;
	ssHead << setw(LEVEL) << "GridLevel  " << sSep;

//	Components
	vector<pair<string,int> > vCmp;
	vCmp.push_back(pair<string,int>("all", ALL_FCT));
	if(bPrintCmps){
		ssHead << setw(COMPONENT) << "Comps" << sSep;
		for(int fct = 0; fct < (int)vDC[0].num_fct(); ++fct){
			stringstream name; name << fct << ": "<< vDC[0].name(fct);
			vCmp.push_back(pair<string,int>(SnipString(name.str(), COMPONENT, 2), fct));
		}
	}

//	Interface
	if(bPrintInterface) {
		ssHead << setw(INTERFACE) << "Parallel" << sSep;
	}

//	Surface
	if(bPrintSurface) {
		ssHead << setw(SURFACE) << "Surface" << sSep;
	}

//	Subsets
	ssHead << setw(NUMBER) << "Domain";
	vector<int> vSubset; vSubset.push_back(ALL_SUBSET);
	if(bPrintSubset) {
		for(int si = 0; si < vDC[0].num_subsets(); ++si){
			stringstream name; name << si << ": "<< vDC[0].subset_name(si);
			ssHead << sSep << setw(NUMBER) << SnipString(name.str(), NUMBER, 2);
			vSubset.push_back(si);
		}
	}

//	size of a line
	int LINE = ssHead.str().size();
	if(LINE < 76) LINE = 76;

	UG_LOG(sLeft << repeat('-', LINE) << sRight << endl);
	UG_LOG(sLeft << left << setw(LINE) << sInfo << right << sRight << endl);
	UG_LOG(sLeft << left << setw(LINE) << sAlgebra << right << sRight << endl);
	UG_LOG(sLeft << setw(LINE) << "" << sRight << endl);
	UG_LOG(sLeft << setw(LINE) << left << ssHead.str() << right << sRight << endl);
	UG_LOG(sLeft << repeat('-', LINE) << sRight << endl);


//	Loop Level
	for(size_t i = 0; i < vDC.size(); ++i){

		const DoFCount& dc = vDC[i];
		const GridLevel gl = dc.grid_level();
		stringstream ssGL; ssGL << gl;

	//	always print unique (w.r.t interface) number
		vector<pair<string, byte> > vInIS;
		vInIS.push_back(pair<string,byte>("unique",UNIQUE_ES));
		vector<pair<string, byte> > vContainsIS;

	//	if PrintInterface: add more output
		if(bPrintInterface){
			vContainsIS.push_back(pair<string,byte>("m (&)",ES_H_MASTER));
			vContainsIS.push_back(pair<string,byte>("s (&)",ES_H_SLAVE));
			vInIS.push_back(pair<string,byte>("all",ALL_ES));

		//	if grid level with ghost: add more output
			if(gl.is_level() && gl.ghosts()){
				vContainsIS.push_back(pair<string,byte>("vm (&)",ES_V_MASTER));
				vContainsIS.push_back(pair<string,byte>("vs (&)",ES_V_SLAVE));
				vInIS.push_back(pair<string,byte>("no (x)",ES_NONE));
				vInIS.push_back(pair<string,byte>("m (x)",ES_H_MASTER));
				vInIS.push_back(pair<string,byte>("s (x)",ES_H_SLAVE));
				vInIS.push_back(pair<string,byte>("vm (x)",ES_V_MASTER));
				vInIS.push_back(pair<string,byte>("vs (x)", ES_V_SLAVE));
				vInIS.push_back(pair<string,byte>("m+vm (x)", ES_H_MASTER | ES_V_MASTER));
				vInIS.push_back(pair<string,byte>("m+vs (x)", ES_H_MASTER | ES_V_SLAVE));
				vInIS.push_back(pair<string,byte>("s+vm (x)", ES_H_SLAVE | ES_V_MASTER));
				vInIS.push_back(pair<string,byte>("s+vs (x)", ES_H_SLAVE | ES_V_SLAVE));
			}
		}

	//	always print unique (w.r.t. surface) number
		vector<pair<string, byte> > vInSS;
		if(gl.is_surface())
			vInSS.push_back(pair<string,byte>("unique",UNIQUE_SS));
		else
			vInSS.push_back(pair<string,byte>("---",UNIQUE_SS));

	//	if PrintSurface and a surface level: add more output
		if(bPrintSurface && gl.is_surface()){
			vInSS.push_back(pair<string,byte>("all",ALL_SS));
			vInSS.push_back(pair<string,byte>("pure",SurfaceView::MG_SURFACE_PURE));
			vInSS.push_back(pair<string,byte>("shadowing",SurfaceView::MG_SURFACE_RIM));
			vInSS.push_back(pair<string,byte>("shadow-cpy",SurfaceView::MG_SHADOW_RIM_COPY));
			vInSS.push_back(pair<string,byte>("shadow-nocpy",SurfaceView::MG_SHADOW_RIM_NONCOPY));
		}

		UG_LOG(sLeft<<setw(LEVEL) << left << ssGL.str() << right);
		stringstream ss; ss << sLeft<<setw(LEVEL)<<"";
		string LvlBegin(ss.str());

	// 	Loop Comps
		for(size_t cmp = 0; cmp < vCmp.size(); ++cmp){

			string LineBegin(LvlBegin);

		//	write component at first appearance
			const int fct = vCmp[cmp].second;
			if(bPrintCmps) {
				if(cmp > 0) UG_LOG(LvlBegin);
				UG_LOG(sSep<<setw(COMPONENT) << left << vCmp[cmp].first << right);
				stringstream ss; ss <<sSep << setw(COMPONENT)<<"";
				LineBegin.append(ss.str());
			}

			bool bPrintBegin = false;

		//	print interface numbers
			for(size_t is = 0; is < vInIS.size(); ++is){
				for(size_t ss = 0; ss < vInSS.size(); ++ss){
					if(bPrintBegin) {UG_LOG(LineBegin);} else bPrintBegin = true;
					if(bPrintInterface) UG_LOG(sSep << setw(INTERFACE) << vInIS[is].first);
					if(bPrintSurface) UG_LOG(sSep << setw(SURFACE) << vInSS[ss].first);
					for(size_t si = 0; si < vSubset.size(); ++si)
						UG_LOG(sSep << setw(NUMBER) << ConvertNumber(dc.num(fct,vSubset[si],vInSS[ss].second,vInIS[is].second),NUMBER,4));
					UG_LOG(sRight << endl);
				}
			}

			for(size_t is = 0; is < vContainsIS.size(); ++is){
				for(size_t ss = 0; ss < vInSS.size(); ++ss){
					if(bPrintBegin) {UG_LOG(LineBegin);} else bPrintBegin = true;
					if(bPrintInterface) UG_LOG(sSep << setw(INTERFACE) << vContainsIS[is].first);
					if(bPrintSurface) UG_LOG(sSep << setw(SURFACE) << vInSS[ss].first);
					for(size_t si = 0; si < vSubset.size(); ++si)
						UG_LOG(sSep << setw(NUMBER) << ConvertNumber(dc.num_contains(fct,vSubset[si],vInSS[ss].second,vContainsIS[is].second),NUMBER,4));
					UG_LOG(sRight << endl);
				}
			}
		}
	}

	UG_LOG(sLeft << repeat('-', LINE) << sRight << endl);

	UG_LOG(left);
	if(sflags.find("legend") != string::npos){
		UG_LOG(sLeft<<setw(LINE)<<" GridLevel: underlying grid part"<<sRight<<endl);
		UG_LOG(sLeft<<setw(LINE)<<"            lev  = level view (all elems in a grid level)"<<sRight<<endl);
		UG_LOG(sLeft<<setw(LINE)<<"            surf = surface view of a level"<<sRight<<endl);
		UG_LOG(sLeft<<setw(LINE)<<"                 = all elems in level + elems without child in lower levels"<<sRight<<endl);
		UG_LOG(sLeft<<setw(LINE)<<"            top  = top surface (leaf elems)"<<sRight<<endl);
		UG_LOG(sLeft<<setw(LINE)<<"            g    = with ghost elems"<<sRight<<endl);

		if(bPrintCmps){
			UG_LOG(sLeft<<setw(LINE)<<" Comps: DoFs in single components"<<sRight<<endl);
		}

		if(bPrintInterface){
			UG_LOG(sLeft<<setw(LINE)<<" Parallel:  DoFs in parallel interfaces"<<sRight<<endl);
			UG_LOG(sLeft<<setw(LINE)<<"            (x) = DoFs exactly matching parallel state"<<sRight<<endl);
			UG_LOG(sLeft<<setw(LINE)<<"            (&) = DoFs containing parallel state"<<sRight<<endl);
			UG_LOG(sLeft<<setw(LINE)<<"            m = (horiz.) master, s = (horiz.) slave"<<sRight<<endl);
			UG_LOG(sLeft<<setw(LINE)<<"            vm = vert. master, vs = vert. slave"<<sRight<<endl);
			UG_LOG(sLeft<<setw(LINE)<<"            no = not contained in interface"<<sRight<<endl);
			UG_LOG(sLeft<<setw(LINE)<<"            all = all DoFs"<<sRight<<endl);
			UG_LOG(sLeft<<setw(LINE)<<"            unique = neglecting DoF copies (as if serial run)"<<sRight<<endl);
			UG_LOG(sLeft<<setw(LINE)<<"                   = no + m (&) + vs (x)"<<sRight<<endl);
		}

		if(bPrintSurface){
			UG_LOG(sLeft<<setw(LINE)<<" Surface:  DoFs in surface states (matching state exactly)"<<sRight<<endl);
			UG_LOG(sLeft<<setw(LINE)<<"           pure = inner surface DoF"<<sRight<<endl);
			UG_LOG(sLeft<<setw(LINE)<<"           shadowing = Shadowing"<<sRight<<endl);
			UG_LOG(sLeft<<setw(LINE)<<"           shadow-cpy = Shadow Copy (i.e. has same type shadowing)"<<sRight<<endl);
			UG_LOG(sLeft<<setw(LINE)<<"           shadow-nocpy = Shadow Non-Copy (i.e. has not same type shadowing)"<<sRight<<endl);
			UG_LOG(sLeft<<setw(LINE)<<"           all = all DoFs"<<sRight<<endl);
			UG_LOG(sLeft<<setw(LINE)<<"           unique = number of uniquely numbered DoFs"<<sRight<<endl);
			UG_LOG(sLeft<<setw(LINE)<<"                  = pure + shadowing + shadow-nocpy"<<sRight<<endl);
		}

		UG_LOG(sLeft<<setw(LINE)<<" Domain: DoFs on whole domain"<<sRight<<endl);
		if(bPrintSubset){
			UG_LOG(sLeft<<setw(LINE)<<" Subset: DoFs on subset only"<<sRight<<endl);
		}

		UG_LOG(sLeft<<setw(LINE)<<""<<sRight<<endl);
		UG_LOG(sLeft<<setw(LINE)<<" Call options: print_statistic(\"opt1, opt2, ...\")"<<sRight<<endl);
		UG_LOG(sLeft<<setw(LINE)<<"   proc:      show DoFs for single proc"<<sRight<<endl);
		UG_LOG(sLeft<<setw(LINE)<<"   subset:    show DoFs per subset"<<sRight<<endl);
		UG_LOG(sLeft<<setw(LINE)<<"   interface: show DoFs per parallel interface"<<sRight<<endl);
		UG_LOG(sLeft<<setw(LINE)<<"   surface:   show DoFs per surface state"<<sRight<<endl);
		UG_LOG(sLeft<<setw(LINE)<<"   component: show DoFs per component"<<sRight<<endl);
		UG_LOG(sLeft<<setw(LINE)<<"   legend:    show this legend"<<sRight<<endl);
		UG_LOG(sLeft<<setw(LINE)<<"   all:       enable all options"<<sRight<<endl);
	} else {
		UG_LOG(sLeft << setw(LINE) << "For Legend and Options: print_statistic(\"legend\")."<< sRight << endl);
	}
	UG_LOG(right);

	UG_LOG(sLeft << repeat('-', LINE) << sRight << endl);
}

void IApproximationSpace::print_statistic() const
{
	print_statistic("subset");
}

void IApproximationSpace::print_statistic(std::string flags) const
{
	PROFILE_FUNC();

//	if nothing printed
	if(m_vDD.empty()){
		static const char* sLeft = " | ";
		static const char* sRight = " | ";
		const int LINE = 60;
		UG_LOG(" --" << repeat('-', LINE) << "-- " << endl);
		UG_LOG(left);
		UG_LOG(sLeft << setw(LINE) << "No DoFDistributions created."<<sRight<<endl);
		UG_LOG(sLeft << setw(LINE) << "NOTE: DoFDistributions are created only on request."<<sRight<<endl);
		UG_LOG(sLeft << setw(LINE) << "      However, you may force creation using:"<<sRight<<endl);
		UG_LOG(sLeft << setw(LINE) << "       - ApproximationSpace::init_levels()"<<sRight<<endl);
		UG_LOG(sLeft << setw(LINE) << "       - ApproximationSpace::init_surfaces()"<<sRight<<endl);
		UG_LOG(sLeft << setw(LINE) << "       - ApproximationSpace::init_top_surface()"<<sRight<<endl);
		UG_LOG(right);
		UG_LOG(" --" << repeat('-', LINE) << "-- " << endl);
		return;
	}

//	Get DoF Counts
	PROFILE_BEGIN(CountLocalDoFStatistic);
	vector<DoFCount> vDC(m_vDD.size());
	for(size_t i = 0; i < vDC.size(); ++i)
		vDC[i] = m_vDD[i]->dof_count();
	PROFILE_END();

	string sflags = ToLower(flags);
	if(sflags.find("all") != string::npos)
		sflags = string("proc, subset, interface, surface, component, legend");

	stringstream ssDDOneProc; ssDDOneProc<<" Number of DoFs";
	int numProcs = 1;
#ifdef UG_PARALLEL
	numProcs = pcl::NumProcs();
	ssDDOneProc<<" (Proc: "<<pcl::ProcRank()<<" of "<< pcl::NumProcs()<<")";
#endif

	bool bPrintOneProc = false;
	if(sflags.find("proc") != string::npos) bPrintOneProc = true;
	if(numProcs == 1) bPrintOneProc = false;

//	Algebra Info
	const int blockSize = DefaultAlgebra::get().blocksize();
	stringstream ssAlgebra; ssAlgebra << " Algebra: ";
	if(blockSize != AlgebraType::VariableBlockSize)
		ssAlgebra<<"Block "<<blockSize<<" (divide by "<<blockSize<<" for #Index)";
	else ssAlgebra <<"Flex";

//	Print infos
	PROFILE_BEGIN(PrintLocalDoFStatistic);
	if(bPrintOneProc)
		PrintDoFCount(vDC, ssDDOneProc.str(), ssAlgebra.str(), sflags);
	PROFILE_END();

	PROFILE_BEGIN(CountGlobalDoFStatistic);
	for(size_t i = 0; i < vDC.size(); ++i)
		vDC[i].sum_values_over_procs(ug::GetLogAssistant().get_output_process());
	PROFILE_END();

	PROFILE_BEGIN(PrintGlobalDoFStatistic);
	PrintDoFCount(vDC, " Number of DoFs (All Procs)", ssAlgebra.str(), sflags);
	PROFILE_END();
}


#ifdef UG_PARALLEL
static size_t NumIndices(const IndexLayout& Layout)
{
	size_t sum = 0;
	for(IndexLayout::const_iterator iter = Layout.begin();
			iter != Layout.end(); ++iter)
		sum += Layout.interface(iter).size();
	return sum;
}
#endif

void IApproximationSpace::print_layout_statistic() const
{
#ifdef UG_PARALLEL
	static const int LEVEL = 14;
	static const int NUMBER = 12;
	static const int SEP = 3;
	static const int LINE = LEVEL + 4*NUMBER + 4*SEP;
	static const char* sSep = " | ";

//	Write header line
	UG_LOG(" --" << repeat('-', LINE) << "-- " << endl);
	stringstream ss; ss << " Index Layouts on Proc " <<
	GetLogAssistant().get_output_process() << " of "<< pcl::NumProcs()
	<< " Procs: " << repeat(' ', 15);
	UG_LOG(sSep << setw(LINE)<<ss.str() << sSep << endl);

	UG_LOG(sSep << setw(LEVEL) << "GridLevel  " << sSep);
	UG_LOG(setw(NUMBER) << "Master  " << sSep);
	UG_LOG(setw(NUMBER) << "Slave  " << sSep);
	UG_LOG(setw(NUMBER) << "vert. Master" << sSep);
	UG_LOG(setw(NUMBER) << "vert. Slave" << sSep << endl);
	UG_LOG(" |-" << repeat('-', LINE) << "-| " << endl);

//	Write Infos for Levels
	for(size_t i = 0; i < m_vDD.size(); ++i){
		stringstream ss; ss << m_vDD[i]->grid_level();
		UG_LOG(sSep << setw(LEVEL) << left << ss.str() << right << sSep);
		UG_LOG(setw(NUMBER) << NumIndices(m_vDD[i]->layouts()->master()) << sSep);
		UG_LOG(setw(NUMBER) << NumIndices(m_vDD[i]->layouts()->slave()) << sSep);
		UG_LOG(setw(NUMBER) << NumIndices(m_vDD[i]->layouts()->vertical_master()) << sSep);
		UG_LOG(setw(NUMBER) << NumIndices(m_vDD[i]->layouts()->vertical_slave()) << sSep << endl);
	}
	UG_LOG(" --" << repeat('-', LINE) << "-- " << endl);

#else
	UG_LOG(" No Layouts in sequential code.\n");
#endif
}

////////////////////////////////////////////////////////////////////////////////
// ApproximationSpace
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain>
ApproximationSpace<TDomain>::
ApproximationSpace(SmartPtr<domain_type> domain)
	: IApproximationSpace(domain->subset_handler(), domain->grid()),
	  m_spDomain(domain)
{
	if(!m_spDomain.valid())
		UG_THROW("Domain, passed to ApproximationSpace, is invalid.");
	if(!m_spMGSH.valid())
		UG_THROW("SubsetHandler, passed to ApproximationSpace, is invalid.");
};

template <typename TDomain>
ApproximationSpace<TDomain>::
ApproximationSpace(SmartPtr<domain_type> domain, const AlgebraType& algebraType)
	: IApproximationSpace(domain->subset_handler(), domain->grid(), algebraType),
	  m_spDomain(domain)
{
	if(!m_spDomain.valid())
		UG_THROW("Domain, passed to ApproximationSpace, is invalid.");
	if(!m_spMGSH.valid())
		UG_THROW("SubsetHandler, passed to ApproximationSpace, is invalid.");
};

} // end namespace ug

#ifdef UG_DIM_1
template class ug::ApproximationSpace<ug::Domain1d>;
#endif
#ifdef UG_DIM_2
template class ug::ApproximationSpace<ug::Domain2d>;
#endif
#ifdef UG_DIM_3
template class ug::ApproximationSpace<ug::Domain3d>;
#endif
