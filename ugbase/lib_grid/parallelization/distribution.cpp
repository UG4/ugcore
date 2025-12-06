/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#include <sstream>
#include "common/static_assert.h"
#include "common/util/table.h"
#include "distribution.h"
#include "distributed_grid.h"
#include "lib_grid/tools/selector_multi_grid.h"
#include "lib_grid/algorithms/selection_util.h"
//ø #include "lib_grid/algorithms/subset_util.h"
#include "lib_grid/algorithms/attachment_util.h"
#include "lib_grid/algorithms/debug_util.h"
#include "parallelization_util.h"
#include "lib_grid/file_io/file_io.h"
#include "lib_grid/global_attachments.h"

//#define LG_DISTRIBUTION_DEBUG
//#define LG_DISTRIBUTION_Z_OUTPUT_TRANSFORM 40
// #define UG_DLOG(id, idx, msg)	UG_LOG("DLOG: " << msg)

using namespace std;

namespace ug{

static DebugID LG_DIST("LG_DIST");


struct TargetProcInfo
{
	TargetProcInfo() = default;
	TargetProcInfo(int pID, byte_t intfcState) :
		procID(pID), interfaceState(intfcState) {}

	int procID;
	byte_t interfaceState; // or-combinations of constants from InterfaceStates
};

using ADistInfo = Attachment<vector<TargetProcInfo> >;

///	Automatically attaches ADistInfo to all elements of a grid.
/**	On destruction, the attached dist info attachments are removed again.
 * You may access the dist-info in an element through the get method.
 * The get method returns a reference to the attached DistInfo object.
 *
 * Make sure that the given grid is valid while the DistInfoSupplier exists.
 */
class DistInfoSupplier{
	public:
		DistInfoSupplier(Grid& grid) : m_grid(grid), m_aDistInfo("distribution-info")
		{
			m_grid.attach_to_all(m_aDistInfo);
			m_aaDistInfoVRT.access(grid, m_aDistInfo);
			m_aaDistInfoEDGE.access(grid, m_aDistInfo);
			m_aaDistInfoFACE.access(grid, m_aDistInfo);
			m_aaDistInfoVOL.access(grid, m_aDistInfo);
		}

		~DistInfoSupplier()
		{
			m_grid.detach_from_all(m_aDistInfo);
		}

		vector<TargetProcInfo>& get(Vertex* vrt)	{return m_aaDistInfoVRT[vrt];}
		vector<TargetProcInfo>& get(Edge* edge)		{return m_aaDistInfoEDGE[edge];}
		vector<TargetProcInfo>& get(Face* face)			{return m_aaDistInfoFACE[face];}
		vector<TargetProcInfo>& get(Volume* vol)		{return m_aaDistInfoVOL[vol];}
		vector<TargetProcInfo>& get(GridObject* obj)
		{
			int objType = obj->base_object_id();
			switch(objType){
				case GridBaseObjectId::VERTEX:	return get(static_cast<Vertex*>(obj));
				case GridBaseObjectId::EDGE:		return get(static_cast<Edge*>(obj));
				case GridBaseObjectId::FACE:		return get(static_cast<Face*>(obj));
				case GridBaseObjectId::VOLUME:	return get(static_cast<Volume*>(obj));
				default:	UG_THROW("Unknown geometric object base type."); 
			}
		}

		ADistInfo dist_info_attachment()	{return m_aDistInfo;}

		template <typename TElem>
		StringStreamTable get_debug_info(TElem* e)
		{
			const vector<TargetProcInfo>& di = get(e);
			StringStreamTable t;
			for(size_t i = 0; i < di.size(); ++i){
				t(0, i+1) << "p" << di[i].procID;
			}

			size_t ri = 1;
			t(ri, 0) << "normal";
			for(size_t i = 0; i < di.size(); ++i)
				t(ri, i+1) << ((di[i].interfaceState & InterfaceStates::IS_NORMAL) != 0);

			ri = 2;
			t(ri, 0) << "vmaster";
			for(size_t i = 0; i < di.size(); ++i)
				t(ri, i+1) << ((di[i].interfaceState & InterfaceStates::IS_VMASTER) != 0);

			ri = 3;
			t(ri, 0) << "vslave";
			for(size_t i = 0; i < di.size(); ++i)
				t(ri, i+1) << ((di[i].interfaceState & InterfaceStates::IS_VSLAVE) != 0);

			ri = 4;
			t(ri, 0) << "dummy";
			for(size_t i = 0; i < di.size(); ++i)
				t(ri, i+1) << ((di[i].interfaceState & InterfaceStates::IS_DUMMY) != 0);

			ri = 5;
			t(ri, 0) << "has parent";
			for(size_t i = 0; i < di.size(); ++i)
				t(ri, i+1) << ((di[i].interfaceState & InterfaceStates::HAS_PARENT) != 0);

			return t;
		}

	private:
	//	copy construction unsupported.
		DistInfoSupplier(const DistInfoSupplier& di) : m_grid(di.m_grid) {}

		Grid& 		m_grid;
		ADistInfo	m_aDistInfo;
		Grid::AttachmentAccessor<Vertex, ADistInfo> m_aaDistInfoVRT;
		Grid::AttachmentAccessor<Edge, ADistInfo> m_aaDistInfoEDGE;
		Grid::AttachmentAccessor<Face, ADistInfo> m_aaDistInfoFACE;
		Grid::AttachmentAccessor<Volume, ADistInfo> m_aaDistInfoVOL;
};


////////////////////////////////////////////////////////////////////////////////
///	Communicates the distribution infos through existing interfaces
/**	Distribution infos are packed into the send buffer for each node and are
 * either merged with existing entries or existing entries are simply overwritten.
 * The merge/overwrite behavior can be chosen through the member method enable_merge.
 */
template <typename TLayout>
class ComPol_SynchronizeDistInfos : public pcl::ICommunicationPolicy<TLayout>
{
	public:
		using Layout = TLayout;
		using GeomObj = typename Layout::Type;
		using Element = typename Layout::Element;
		using Interface = typename Layout::Interface;
		using InterfaceIter = typename Interface::const_iterator;

		ComPol_SynchronizeDistInfos(DistInfoSupplier& distInfos, bool merge) :
			m_distInfos(distInfos), m_mergeEnabled(merge) {}

		~ComPol_SynchronizeDistInfos() override = default;

		void enable_merge(bool enable) {m_mergeEnabled = enable;}
		bool merge_enabled() {return m_mergeEnabled;}

		int get_required_buffer_size(const Interface& interface) override {return -1;}

	///	write target processes and move-flag
		bool collect(BinaryBuffer& buff, const Interface& intfc) override
		{
			for(auto iter = intfc.begin(); iter != intfc.end(); ++iter){
				Element elem = intfc.get_element(iter);
				Serialize(buff, m_distInfos.get(elem));
			}
			return true;
		}

	///	read target processes and move-flag
		bool extract(BinaryBuffer& buff, const Interface& intfc) override
		{
			if(m_mergeEnabled){
				vector<TargetProcInfo> tpInfo;
				for(auto iter = intfc.begin(); iter != intfc.end(); ++iter){
					tpInfo.clear();
					Deserialize(buff, tpInfo);

					Element elem = intfc.get_element(iter);
					vector<TargetProcInfo>& tpInfoDest = m_distInfos.get(elem);
					size_t initialInfoSize = tpInfoDest.size();

					for(size_t i_src = 0; i_src < tpInfo.size(); ++i_src){
						int procID = tpInfo[i_src].procID;
						bool gotOne = false;

					//	we only have to check entries up to initialInfoSize, since
					//	all following entries have been added during this operation.
					//	Since there are no double entries in tpInfo, there's no
					//	need to check against those new entries in tpInfoDest.
						for(size_t i = 0; i < initialInfoSize; ++i){
							if(procID == tpInfoDest[i].procID){
								tpInfoDest[i].interfaceState
												|= tpInfo[i_src].interfaceState;

								gotOne = true;
								break;
							}
						}

						if(!gotOne)
							tpInfoDest.push_back(tpInfo[i_src]);
					}
				}
			}
			else{
				for(auto iter = intfc.begin(); iter != intfc.end(); ++iter){
					Element elem = intfc.get_element(iter);
					vector<TargetProcInfo>& tpInfo = m_distInfos.get(elem);
					tpInfo.clear();
					Deserialize(buff, tpInfo);
				}
			}

			return true;
		}

	protected:
		DistInfoSupplier& 	m_distInfos;
		bool				m_mergeEnabled;
};


////////////////////////////////////////////////////////////////////////////////
template <typename TElem>
static void SynchronizeDistInfos(MultiGrid& mg, DistInfoSupplier& distInfos)
{
	using ElemLayout = typename GridLayoutMap::Types<TElem>::Layout;
	GridLayoutMap& glm = mg.distributed_grid_manager()->grid_layout_map();
	pcl::InterfaceCommunicator<ElemLayout>	com;
	ComPol_SynchronizeDistInfos<ElemLayout>	compolSync(distInfos, true);

	compolSync.enable_merge(true);
	com.exchange_data(glm, InterfaceNodeTypes::INT_H_SLAVE, InterfaceNodeTypes::INT_H_MASTER, compolSync);
	com.communicate();

	compolSync.enable_merge(false);
	com.exchange_data(glm, InterfaceNodeTypes::INT_H_MASTER, InterfaceNodeTypes::INT_H_SLAVE, compolSync);
	com.communicate();

//	if an element is a multi-v-master, all v-master copies are connected to
//	all v-slave copies. However, multi-v-masters are not connected with one another
//	through h-interfaces. That's a pity and requires this triple communication
	compolSync.enable_merge(true);
	com.exchange_data(glm, InterfaceNodeTypes::INT_V_MASTER, InterfaceNodeTypes::INT_V_SLAVE, compolSync);
	com.communicate();

	compolSync.enable_merge(true);
	com.exchange_data(glm, InterfaceNodeTypes::INT_V_SLAVE, InterfaceNodeTypes::INT_V_MASTER, compolSync);
	com.communicate();

	compolSync.enable_merge(false);
	com.exchange_data(glm, InterfaceNodeTypes::INT_V_MASTER, InterfaceNodeTypes::INT_V_SLAVE, compolSync);
	com.communicate();
}

#ifdef LG_DISTRIBUTION_DEBUG
////////////////////////////////////////////////////////////////////////////////
static void SaveDistSelectorToFile(MGSelector& msel, const char* filename)
{
//	create a subset handler which holds different subsets for the different selection states
	MultiGrid& mg = *msel.multi_grid();
	SubsetHandler sh(mg);

	for(size_t lvl = 0; lvl < msel.num_levels(); ++lvl){
		for(MGSelector::traits<Volume>::level_iterator iter = msel.begin<Volume>(lvl);
			iter != msel.end<Volume>(lvl); ++iter)
		{
			sh.assign_subset(*iter, msel.get_selection_status(*iter));
		}

		for(MGSelector::traits<Face>::level_iterator iter = msel.begin<Face>(lvl);
			iter != msel.end<Face>(lvl); ++iter)
		{
			sh.assign_subset(*iter, msel.get_selection_status(*iter));
		}

		for(MGSelector::traits<Edge>::level_iterator iter = msel.begin<Edge>(lvl);
			iter != msel.end<Edge>(lvl); ++iter)
		{
			sh.assign_subset(*iter, msel.get_selection_status(*iter));
		}

		for(MGSelector::traits<Vertex>::level_iterator iter = msel.begin<Vertex>(lvl);
			iter != msel.end<Vertex>(lvl); ++iter)
		{
			sh.assign_subset(*iter, msel.get_selection_status(*iter));
		}
	}

	const char* subsetNames[] = {"unassigned", "normal", "vmaster", "normal+vmaster",
								 "vslave", "normal+vslave", "vmaster+vslave",
								 "normal+vmaster+vslave", "dummy", "normal+dummy",
								 "vmaster+dummy", "normal+vmaster+dummy", "vslave+dummy",
								 "normal+vslave+dummy", "vmaster+vslave+dummy",
								 "normal+vmaster+vslave+dummy"};

	for(int i = 0; i < 16; ++i)
		sh.subset_info(i).name = subsetNames[i];

	AssignSubsetColors(sh);
	EraseEmptySubsets(sh);
	//SaveGridHierarchyTransformed(mg, sh, filename, LG_DISTRIBUTION_Z_OUTPUT_TRANSFORM);
    SaveGridToFile(mg, sh, filename);
}

////////////////////////////////////////////////////////////////////////////////
static void SaveDistInfosToFile(MultiGrid& mg, DistInfoSupplier& infoSupplier,
								const char* filename)
{
//	create a subset handler which holds different subsets for the different selection states
	SubsetHandler sh(mg);

//	write a file for each adressed process
	for(int pi = 0; pi < pcl::NumProcs(); ++pi){
		sh.clear();

		for(MultiGrid::traits<Volume>::iterator iter = mg.begin<Volume>();
			iter != mg.end<Volume>(); ++iter)
		{
			vector<TargetProcInfo>& infos = infoSupplier.get(*iter);
			for(size_t i = 0; i < infos.size(); ++i){
				if(infos[i].procID == pi)
					sh.assign_subset(*iter, infos[i].interfaceState);
			}
		}

		for(MultiGrid::traits<Face>::iterator iter = mg.begin<Face>();
			iter != mg.end<Face>(); ++iter)
		{
			vector<TargetProcInfo>& infos = infoSupplier.get(*iter);
			for(size_t i = 0; i < infos.size(); ++i){
				if(infos[i].procID == pi)
					sh.assign_subset(*iter, infos[i].interfaceState);
			}
		}

		for(MultiGrid::traits<Edge>::iterator iter = mg.begin<Edge>();
			iter != mg.end<Edge>(); ++iter)
		{
			vector<TargetProcInfo>& infos = infoSupplier.get(*iter);
			for(size_t i = 0; i < infos.size(); ++i){
				if(infos[i].procID == pi)
					sh.assign_subset(*iter, infos[i].interfaceState);
			}
		}

		for(MultiGrid::traits<Vertex>::iterator iter = mg.begin<Vertex>();
			iter != mg.end<Vertex>(); ++iter)
		{
			vector<TargetProcInfo>& infos = infoSupplier.get(*iter);
			for(size_t i = 0; i < infos.size(); ++i){
				if(infos[i].procID == pi)
					sh.assign_subset(*iter, infos[i].interfaceState);
			}
		}

		const char* subsetNames[] = {"unassigned", "normal", "vmaster", "normal+vmaster",
									 "vslave", "normal+vslave", "vmaster+vslave",
									 "normal+vmaster+vslave", "dummy", "normal+dummy",
									 "vmaster+dummy", "normal+vmaster+dummy", "vslave+dummy",
									 "normal+vslave+dummy", "vmaster+vslave+dummy",
									 "normal+vmaster+vslave+dummy"};

		for(int i = 0; i < 16; ++i)
			sh.subset_info(i).name = subsetNames[i];

		AssignSubsetColors(sh);
		EraseEmptySubsets(sh);
		if(sh.num_subsets() > 0){
			stringstream ss;
			ss << filename << "_p" << pcl::ProcRank() << "_for_p" << pi << ".ugx";
			//SaveGridHierarchyTransformed(mg, sh, ss.str().c_str(), LG_DISTRIBUTION_Z_OUTPUT_TRANSFORM);
		        SaveGridToFile(mg, sh,ss.str().c_str());
		}
	}
}

template <typename TElem>
static void WriteDistInfosToTextFile(MultiGrid& mg, DistInfoSupplier& infoSupplier,
									 const char* filename)
{
	using TElemIter = typename MultiGrid::traits<TElem>::iterator;

	Table<std::stringstream> table(mg.num<TElem>() + 1, 3);
	table(0, 0) << "lvl";	table(0, 1) << "center";	table(0, 2) << "interface states";

	int row = 1;
	for(size_t lvl = 0; lvl < mg.num_levels(); ++lvl){
		for(TElemIter iter = mg.begin<TElem>(lvl); iter != mg.end<TElem>(lvl);
			++iter, ++row)
		{
			TElem* e = *iter;

			table(row, 0) << lvl;
			table(row, 1) << GetGridObjectCenter(mg, e);

			vector<TargetProcInfo>& infos = infoSupplier.get(e);

			for(size_t i = 0; i < infos.size(); ++i){
				table(row, 2) << "p" << infos[i].procID << ": ";
				byte is = infos[i].interfaceState;
				if(is & IS_NORMAL)	table(row, 2) << "normal ";
				if(is & IS_VMASTER)	table(row, 2) << "vmaster ";
				if(is & IS_VSLAVE)	table(row, 2) << "vslave ";
				if(is & IS_DUMMY)	table(row, 2) << "dummy ";
			}

			table(row, 2) << "| ";
		}
	}

	ofstream out(filename);
	if(!out){
		UG_THROW("Couldn't open file " << filename << " for output.");
	}
	out << table;
	out.close();
}

template <typename TElem>
static string LocateElement(MultiGrid& mg, TElem* e)
{
	stringstream ssLocator;
	ssLocator << "at " << GetGridObjectCenter(mg, e)
			  << " on level " << mg.get_level(e);
	return ssLocator.str();
}

template <typename TElem>
static bool PerformValidityCheck(DistributedGridManager& dgm)
{
	using TElemIter = typename Grid::traits<TElem>::iterator;

	bool isValid = true;

	UG_LOG("DEBUG: Performing validity check on distributed grid ");
	switch(TElem::BASE_OBJECT_ID){
	case VERTEX:	UG_LOG("for vertices:\n"); break;
	case EDGE:		UG_LOG("for edges:\n"); break;
	case FACE:		UG_LOG("for faces:\n"); break;
	case VOLUME:	UG_LOG("for volumes:\n"); break;
	}

	MultiGrid& mg = *dgm.get_assigned_grid();
	for(TElemIter iter = mg.begin<TElem>(); iter != mg.end<TElem>(); ++iter)
	{
		TElem* e = *iter;
	//	Make sure that pure vertical masters (ghosts) do not have children
		if(dgm.is_ghost(e)){
			if(mg.has_children(e)){
				UG_LOG("  Ghost has child " << LocateElement(mg, e) << endl);
				isValid = false;
			}
		}
	//	Make sure that pure vertical slaves do not have parents
		else if(dgm.contains_status(e, ES_V_SLAVE)
				&& (!dgm.is_in_horizontal_interface(e)))
		{
			if(mg.get_parent(e)){
				UG_LOG("  Pure vertical slave has parent " << LocateElement(mg, e) << endl);
				isValid = false;
			}
		}
	}

	UG_LOG("DEBUG: Validity check done with result ");
	if(isValid){
		UG_LOG("SUCCESS\n");
	}
	else{
		UG_LOG("FAIL\n");
	}

	return isValid;
}

static bool PerformValidityCheck(DistributedGridManager& dgm)
{
	bool isValid = true;
	isValid &= PerformValidityCheck<Vertex>(dgm);
	isValid &= PerformValidityCheck<Edge>(dgm);
	isValid &= PerformValidityCheck<Face>(dgm);
	isValid &= PerformValidityCheck<Volume>(dgm);
	return isValid;
}

#endif

////////////////////////////////////////////////////////////////////////////////
// ATTENTION - THIS DOESN'T REALLY WORK!
//	mpirun -n 4 ugshell -ex adaptive_mg/moving_front.lua -redistributionSteps 2 -redistributionProcs 2
//template <typename TElem>
//void AdjustGhostSelection(MGSelector& msel, ISelector::status_t status)
//{
//	DistributedGridManager& dgm = *msel.grid()->distributed_grid_manager();
//	GridLayoutMap& glm = dgm.grid_layout_map();
//
//	if(!glm.has_layout<TElem>(INT_V_MASTER))
//		return;
//
// using Layout = typename GridLayoutMap::Types<TElem>::Layout;
// using LIter = typename Layout::iterator;
// using Interface = typename Layout::Interface;
// using IIter = typename Interface::iterator;
//
//	Layout& layout = glm.get_layout<TElem>(INT_V_MASTER);
//	for(size_t lvl = 0; lvl < layout.num_levels(); ++lvl){
//		for(LIter liter = layout.begin(lvl); liter != layout.end(lvl); ++liter){
//			Interface& intfc = layout.interface(liter);
//			for(IIter iiter = intfc.begin(); iiter != intfc.end(); ++iiter){
//				TElem* e = intfc.get_element(iiter);
//				if(dgm.is_ghost(e)){
//					msel.select(e, status);
//				}
//			}
//		}
//	}
//}
//
//void AdjustGhostSelection(MGSelector& msel, ISelector::status_t status)
//{
//	UG_LOG("DEBUG: AdjustGhostSelection, #sel-vrts before: " << msel.num<Vertex>() << endl);
//	Grid& g = *msel.grid();
//	if(g.num_vertices())
//		AdjustGhostSelection<Vertex>(msel, status);
//	if(g.num_edges())
//		AdjustGhostSelection<Edge>(msel, status);
//	if(g.num_faces())
//		AdjustGhostSelection<Face>(msel, status);
//	if(g.num_volumes())
//		AdjustGhostSelection<Volume>(msel, status);
//	UG_LOG("DEBUG: AdjustGhostSelection, #sel-vrts after: " << msel.num<Vertex>() << endl);
//}


////////////////////////////////////////////////////////////////////////////////
///	Recursively selects unselected sides.
template <typename TElem>
static void SelectAssociatedSides(MGSelector& msel, TElem* e,
								  ISelector::status_t status = ISelector::SELECTED)
{
	//UG_DLOG(LG_DIST, 3, "dist-start: SelectAssociatedSides\n");
	GDIST_PROFILE_FUNC();

	UG_ASSERT(msel.multi_grid(), "");
	MultiGrid& mg = *msel.multi_grid();

	using TSide = typename TElem::side;
	typename MultiGrid::traits<TSide>::secure_container sides;

	mg.associated_elements(sides, e);
	for(size_t i = 0; i < sides.size(); ++i){
		TSide* s = sides[i];
		//if(!msel.is_selected(sides[i])){
			ISelector::status_t nstate = status | msel.get_selection_status(s);
			msel.select(s, nstate);
			if(TElem::HAS_SIDES)
				SelectAssociatedSides(msel, s, nstate);
		//}
	}

	//UG_DLOG(LG_DIST, 3, "dist-stop: SelectAssociatedSides\n");
}


////////////////////////////////////////////////////////////////////////////////
/**	selects unselected constrained elements of all selected constraining elements
 * and associated unselected low-dim elems. An exception is made for constraining
 * elements which are pure vertical masters. Associated constrained elements won't
 * be selected in this case, since pure vertical masters mustn't have children.*/
static void SelectAssociatedConstrainedElements(MGSelector& msel,
								ISelector::status_t status = ISelector::SELECTED)
{
	UG_DLOG(LG_DIST, 3, "dist-start: SelectAssociatedConstrainedElements\n");
	GDIST_PROFILE_FUNC();

	const bool selectAll = true;

//	constraining triangles
	{
		using TElem = ConstrainingTriangle;
		using TIter = MGSelector::traits<TElem>::level_iterator;
		for(size_t lvl = 0; lvl < msel.num_levels(); ++lvl){
			for(TIter iter = msel.begin<TElem>(lvl);
				iter != msel.end<TElem>(lvl); ++iter)
			{
				ConstrainingFace* e = *iter;
			//	we won't select constrained elements of pure v-masters!
				if(msel.get_selection_status(e) == InterfaceStates::IS_VMASTER)
					continue;

				for(size_t i = 0; i < e->num_constrained_vertices(); ++i){
					Vertex* cd = e->constrained_vertex(i);
					ISelector::status_t nstate = status | msel.get_selection_status(cd);
					if(selectAll || !msel.is_selected(cd)){
						msel.select(cd, nstate);
					}
				}
				for(size_t i = 0; i < e->num_constrained_edges(); ++i){
					Edge* cd = e->constrained_edge(i);
					ISelector::status_t nstate = status | msel.get_selection_status(cd);
					if(selectAll || !msel.is_selected(cd)){
						msel.select(cd, nstate);
						SelectAssociatedSides(msel, cd, nstate);
					}
				}
				for(size_t i = 0; i < e->num_constrained_faces(); ++i){
					Face* cd = e->constrained_face(i);
					ISelector::status_t nstate = status | msel.get_selection_status(cd);
					if(selectAll || !msel.is_selected(cd)){
						msel.select(cd, nstate);
						SelectAssociatedSides(msel, cd, nstate);
					}
				}
			}
		}
	}

//	constraining quadrilaterals
	{
		using TElem = ConstrainingQuadrilateral;
		using TIter = MGSelector::traits<TElem>::level_iterator;
		for(size_t lvl = 0; lvl < msel.num_levels(); ++lvl){
			for(TIter iter = msel.begin<TElem>(lvl);
				iter != msel.end<TElem>(lvl); ++iter)
			{
				ConstrainingFace* e = *iter;
			//	we won't select constrained elements of pure v-masters!
				if(msel.get_selection_status(e) == InterfaceStates::IS_VMASTER)
					continue;

				for(size_t i = 0; i < e->num_constrained_vertices(); ++i){
					Vertex* cd = e->constrained_vertex(i);
					ISelector::status_t nstate = status | msel.get_selection_status(cd);
					if(selectAll || !msel.is_selected(cd)){
						msel.select(cd, nstate);
					}
				}
				for(size_t i = 0; i < e->num_constrained_edges(); ++i){
					Edge* cd = e->constrained_edge(i);
					ISelector::status_t nstate = status | msel.get_selection_status(cd);
					if(selectAll || !msel.is_selected(cd)){
						msel.select(cd, nstate);
						SelectAssociatedSides(msel, cd, nstate);
					}
				}
				for(size_t i = 0; i < e->num_constrained_faces(); ++i){
					Face* cd = e->constrained_face(i);
					ISelector::status_t nstate = status | msel.get_selection_status(cd);
					if(selectAll || !msel.is_selected(cd)){
						msel.select(cd, nstate);
						SelectAssociatedSides(msel, cd, nstate);
					}
				}
			}
		}
	}

//	constraining edges
	{
		using TElem = ConstrainingEdge;
		using TIter = MGSelector::traits<TElem>::level_iterator;
		for(size_t lvl = 0; lvl < msel.num_levels(); ++lvl){
			for(TIter iter = msel.begin<TElem>(lvl);
				iter != msel.end<TElem>(lvl); ++iter)
			{
				ConstrainingEdge* e = *iter;
			//	we won't select constrained elements of pure v-masters!
				if(msel.get_selection_status(e) == InterfaceStates::IS_VMASTER)
					continue;

				for(size_t i = 0; i < e->num_constrained_vertices(); ++i){
					Vertex* cd = e->constrained_vertex(i);
					ISelector::status_t nstate = status | msel.get_selection_status(cd);
					if(selectAll || !msel.is_selected(cd)){
						msel.select(cd, nstate);
					}
				}
				for(size_t i = 0; i < e->num_constrained_edges(); ++i){
					Edge* cd = e->constrained_edge(i);
					ISelector::status_t nstate = status | msel.get_selection_status(cd);
					if(selectAll || !msel.is_selected(cd)){
						msel.select(cd, nstate);
						SelectAssociatedSides(msel, cd, nstate);
					}
				}
			}
		}
	}
	UG_DLOG(LG_DIST, 3, "dist-stop: SelectAssociatedConstrainedElements\n");
}

////////////////////////////////////////////////////////////////////////////////
/**	selects unselected constraining elements of all selected constrained elements
 * and associated unselected low-dim elems.*/
//static void SelectAssociatedConstrainingElements(MGSelector& msel,
//								ISelector::status_t status = ISelector::SELECTED)
//{
//	UG_DLOG(LG_DIST, 3, "dist-start: SelectAssociatedConstrainingElements\n");
//	GDIST_PROFILE_FUNC();
//
//	const bool selectAll = true;
//
////	constrained triangles
//	{
// using TElem = ConstrainedTriangle;
// using TIter = MGSelector::traits<TElem>::level_iterator;
//		for(size_t lvl = 0; lvl < msel.num_levels(); ++lvl){
//			for(TIter iter = msel.begin<TElem>(lvl);
//				iter != msel.end<TElem>(lvl); ++iter)
//			{
//				ConstrainedFace* e = *iter;
//				if(GridObject* cg = e->get_constraining_object()){
//				//	we won't select pure v-masters!
////					if(msel.get_selection_status(cg) == IS_VMASTER)
////						continue;
//
//					ISelector::status_t nstate = status | msel.get_selection_status(cg);
//					if(selectAll || !msel.is_selected(cg)){
//						msel.select(cg, nstate);
//						UG_ASSERT(dynamic_cast<ConstrainingFace*>(cg),
//								  "constraining object of a face has to be a "
//								  "ConstrainingFace!");
//						SelectAssociatedSides(msel, static_cast<Face*>(cg), nstate);
//					}
//				}
//			}
//		}
//	}
//
////	constrained quadrilaterals
//	{
// using TElem = ConstrainedQuadrilateral;
// using TIter = MGSelector::traits<TElem>::level_iterator;
//		for(size_t lvl = 0; lvl < msel.num_levels(); ++lvl){
//			for(TIter iter = msel.begin<TElem>(lvl);
//				iter != msel.end<TElem>(lvl); ++iter)
//			{
//				ConstrainedFace* e = *iter;
//				if(GridObject* cg = e->get_constraining_object()){
//				//	we won't select pure v-masters!
////					if(msel.get_selection_status(cg) == IS_VMASTER)
////						continue;
//
//					ISelector::status_t nstate = status | msel.get_selection_status(cg);
//					if(selectAll || !msel.is_selected(cg)){
//						msel.select(cg, nstate);
//						UG_ASSERT(dynamic_cast<Face*>(cg),
//								  "constraining object of a face has to be a "
//								  "Face!");
//						SelectAssociatedSides(msel, static_cast<Face*>(cg), nstate);
//					}
//				}
//			}
//		}
//	}
//
////	constrained edges
//	{
//		using TElem = ConstrainedEdge;
//		using TIter = MGSelector::traits<TElem>::level_iterator;
//		for(size_t lvl = 0; lvl < msel.num_levels(); ++lvl){
//			for(TIter iter = msel.begin<TElem>(lvl);
//				iter != msel.end<TElem>(lvl); ++iter)
//			{
//				ConstrainedEdge* e = *iter;
//				if(GridObject* cg = e->get_constraining_object()){
//				//	we won't select pure v-masters!
////					if(msel.get_selection_status(cg) == IS_VMASTER)
////						continue;
//
//					ISelector::status_t nstate = status | msel.get_selection_status(cg);
//					if(selectAll || !msel.is_selected(cg)){
//						msel.select(cg, nstate);
//						switch(cg->base_object_id()){
//						case EDGE:
//							SelectAssociatedSides(msel, static_cast<Edge*>(cg), nstate);
//							break;
//						case FACE:
//							SelectAssociatedSides(msel, static_cast<Face*>(cg), nstate);
//							break;
//						}
//					}
//				}
//			}
//		}
//	}
//
////	constrained vertices
//	{
// using TElem = ConstrainedVertex;
// using TIter = MGSelector::traits<TElem>::level_iterator;
//		for(size_t lvl = 0; lvl < msel.num_levels(); ++lvl){
//			for(TIter iter = msel.begin<TElem>(lvl);
//				iter != msel.end<TElem>(lvl); ++iter)
//			{
//				ConstrainedVertex* e = *iter;
//				if(GridObject* cg = e->get_constraining_object()){
//				//	we won't select pure v-masters!
////					if(msel.get_selection_status(cg) == IS_VMASTER)
////						continue;
//
//					ISelector::status_t nstate = status | msel.get_selection_status(cg);
//					if(selectAll || !msel.is_selected(cg)){
//						msel.select(cg, nstate);
//						switch(cg->base_object_id()){
//						case EDGE:
//							SelectAssociatedSides(msel, static_cast<Edge*>(cg), nstate);
//							break;
//						case FACE:
//							SelectAssociatedSides(msel, static_cast<Face*>(cg), nstate);
//							break;
//						}
//					}
//				}
//			}
//		}
//	}
//	UG_DLOG(LG_DIST, 3, "dist-stop: SelectAssociatedConstrainingElements\n");
//}




static void SelectChildrenOfSelectedShadowRimFaces
(
	MGSelector& msel,
	ISelector::status_t status = ISelector::SELECTED
)
{
	UG_DLOG(LG_DIST, 3, "dist-start: SelectChildrenOfSelectedShadowVertices\n");
	GDIST_PROFILE_FUNC();

	UG_ASSERT(msel.multi_grid(), "The selector has to operate on a grid!");
	MultiGrid& mg = *msel.multi_grid();
	DistributedGridManager& distGridMgr = *mg.distributed_grid_manager();

	Grid::volume_traits::secure_container vols;

	for (size_t lvl = 0; lvl < msel.num_levels(); ++lvl)
	{
		MGSelector::traits<Face>::level_iterator iter = msel.begin<Face>(lvl);
		for (; iter != msel.end<Face>(lvl); ++iter)
		{
			Face* face = *iter;
			if (msel.get_selection_status(face) == InterfaceStates::IS_VMASTER)
				continue;

			size_t numCh = mg.num_child_faces(face);
			if (!numCh) continue;

			// check whether face has an associated volume which does not have children
			// TODO: This is not a correct criterion for surface rim.
			mg.associated_elements(vols, face);
			for (size_t i = 0; i < vols.size(); ++i)
			{
				if (!(distGridMgr.is_ghost(vols[i]) || mg.has_children(vols[i])))
				{
					for (size_t ch = 0; ch < numCh; ++ch)
					{
						Face* child = mg.get_child_face(face, ch);
						msel.select(child, msel.get_selection_status(child) | status);
					}
					break;
				}
			}
		}
	}
	UG_DLOG(LG_DIST, 3, "dist-stop: SelectChildrenOfSelectedShadowVertices\n");
}


static void SelectChildrenOfSelectedShadowRimEdges
(
	MGSelector& msel,
	ISelector::status_t status = ISelector::SELECTED
)
{
	UG_DLOG(LG_DIST, 3, "dist-start: SelectChildrenOfSelectedShadowVertices\n");
	GDIST_PROFILE_FUNC();

	UG_ASSERT(msel.multi_grid(), "The selector has to operate on a grid!");
	MultiGrid& mg = *msel.multi_grid();
	DistributedGridManager& distGridMgr = *mg.distributed_grid_manager();

	Grid::face_traits::secure_container faces;

	for (size_t lvl = 0; lvl < msel.num_levels(); ++lvl)
	{
		MGSelector::traits<Edge>::level_iterator iter = msel.begin<Edge>(lvl);
		for (; iter != msel.end<Edge>(lvl); ++iter)
		{
			Edge* edge = *iter;
			if (msel.get_selection_status(edge) == InterfaceStates::IS_VMASTER)
				continue;

			size_t numCh = mg.num_child_edges(edge);
			if (!numCh) continue;

			// check whether edge has an associated face which does not have children
			// TODO: This is not a correct criterion for surface rim.
			mg.associated_elements(faces, edge);
			for (size_t i = 0; i < faces.size(); ++i)
			{
				if (!(distGridMgr.is_ghost(faces[i]) || mg.has_children(faces[i])))
				{
					for (size_t ch = 0; ch < numCh; ++ch)
					{
						Edge* child = mg.get_child_edge(edge, ch);
						msel.select(child, msel.get_selection_status(child) | status);
					}
					break;
				}
			}
		}
	}
	UG_DLOG(LG_DIST, 3, "dist-stop: SelectChildrenOfSelectedShadowVertices\n");
}


////////////////////////////////////////////////////////////////////////////////
///	Recursively selects all children of selected vertices
/**	This method is required, since if a distributed vertex has a child which
 * is not connected to other distributed elements, then the child wouldn't be selected.
 * Note that for edges and faces the methods SelectAssociatedConstrainedElements
 * takes care of this. (This is only true if there is no anisotropic refinement!)
 * Children of pure vertical masters won't be selected, since those mustn't have
 * children.
 */
static void SelectChildrenOfSelectedShadowVertices(MGSelector& msel,
								ISelector::status_t status = ISelector::SELECTED)
{
	UG_DLOG(LG_DIST, 3, "dist-start: SelectChildrenOfSelectedShadowVertices\n");
	GDIST_PROFILE_FUNC();

	UG_ASSERT(msel.multi_grid(), "The selector has to operate on a grid!");
	MultiGrid& mg = *msel.multi_grid();
	DistributedGridManager& distGridMgr = *mg.distributed_grid_manager();

	Grid::edge_traits::secure_container edges;

	for(size_t lvl = 0; lvl < msel.num_levels(); ++lvl){
		for(MGSelector::traits<Vertex>::level_iterator iter = msel.begin<Vertex>(lvl);
			iter != msel.end<Vertex>(lvl); ++iter)
		{
			Vertex* vrt = *iter;
			if(msel.get_selection_status(vrt) == InterfaceStates::IS_VMASTER)
				continue;

			Vertex* child = mg.get_child_vertex(vrt);
			if(!child)
				continue;

		//	check whether vrt has an associated edge which does not have children
			mg.associated_elements(edges, vrt);
			for(size_t i = 0; i < edges.size(); ++i){
				if(!(distGridMgr.is_ghost(edges[i]) || mg.has_children(edges[i]))){
					msel.select(child, msel.get_selection_status(child) | status);
					break;
				}
			}
		}
	}
	UG_DLOG(LG_DIST, 3, "dist-stop: SelectChildrenOfSelectedShadowVertices\n");
}


////////////////////////////////////////////////////////////////////////////////
/**	The method operates on selected entries only. Make sure that all elements
 * of type TElem which are being sent to a process are selected.
 *
 * If a selected element has no children and if it is a vertical master, it will
 * be marked as vertical master again.
 *
 * If a selected element has unselected children, then those children will be
 * selected as vertical master.
 *
 * This method only works correctly if called for the elements of highest dimension.
 */
template <typename TElem>
static void AssignVerticalMasterAndSlaveStates(MGSelector& msel, bool partitionForLocalProc)
{
	UG_DLOG(LG_DIST, 3, "dist-start: AssignVerticalMasterAndSlaveStates\n");
	GDIST_PROFILE_FUNC();

	UG_ASSERT(msel.multi_grid(), "Selector has to operate on a MultiGrid");
	MultiGrid& mg = *msel.multi_grid();
	DistributedGridManager& distGridMgr = *mg.distributed_grid_manager();

//	we start on the highest level and go downwards to avoid side
//	effects from repeated selection adjustment.
	using TIter = typename MGSelector::traits<TElem>::level_iterator;
	for(int lvl = (int)msel.num_levels() - 1; lvl >= 0 ; --lvl){
		for(TIter iter = msel.begin<TElem>(lvl);
			iter != msel.end<TElem>(lvl);)
		{
			TElem* e = *iter;
			++iter;

		//TODO: Check this! if e is VSLAVE and is sent to another proc,
		//		it's children have to be sent there too, since they will
		//		be VMASTER on this new proc!
			if((msel.get_selection_status(e) & InterfaceStates::IS_VMASTER)
				|| (msel.get_selection_status(e) & InterfaceStates::IS_VSLAVE))
			{
			//	nothing to do here...
				continue;
			}

		//	assign vertical master states first
			size_t numChildren = mg.num_children<TElem>(e);
			GridObject* parent = mg.get_parent(e);
			bool parentIsSelected = false;
			if(parent)
				parentIsSelected = msel.is_selected(parent);

			if(numChildren){
				for(size_t i = 0; i < numChildren; ++i){
					auto* c = mg.get_child<TElem>(e, i);
					if(!msel.is_selected(c))
						msel.select(c, InterfaceStates::IS_VMASTER);
				}
			}
			else if(distGridMgr.contains_status(e, ElementStatusTypes::ES_V_MASTER)){
				if(parentIsSelected || (partitionForLocalProc && (lvl == 0))){
					msel.select(e, InterfaceStates::IS_VMASTER);
					continue;
				}
				else{
					msel.deselect(e);
					continue;
				}
			}

		//	and now slave states
			if(parent){
				if(!msel.is_selected(parent))
					msel.select(e, InterfaceStates::IS_VSLAVE);
			}
			else{
				if(distGridMgr.contains_status(e, ElementStatusTypes::ES_V_SLAVE))
					msel.select(e, InterfaceStates::IS_VSLAVE);
			}
		}
	}

	UG_DLOG(LG_DIST, 3, "dist-stop: AssignVerticalMasterAndSlaveStates\n");
}

////////////////////////////////////////////////////////////////////////////////
/**	VSlaves will be ignored.*/
template <typename TElem>
static void SelectUnselectedRootElementsAsVMasters(MGSelector& msel)
{
	UG_DLOG(LG_DIST, 3, "dist-start: SelectUnselectedRootElementsAsVMasters\n");
	GDIST_PROFILE_FUNC();

	using TIter = typename Grid::traits<TElem>::iterator;

	UG_ASSERT(msel.multi_grid(), "Selector has to operate on a MultiGrid");
	MultiGrid& mg = *msel.multi_grid();
	DistributedGridManager& distGridMgr = *mg.distributed_grid_manager();

	for(TIter iter = mg.begin<TElem>(0); iter != mg.end<TElem>(0); ++iter){
		if(!msel.is_selected(*iter)){
			if(!distGridMgr.contains_status(*iter, ElementStatusTypes::ES_V_SLAVE))
				msel.select(*iter, InterfaceStates::IS_VMASTER);
		}
	}
	UG_DLOG(LG_DIST, 3, "dist-stop: SelectUnselectedRootElementsAsVMasters\n");
}

////////////////////////////////////////////////////////////////////////////////
/**	VMasters will be ignored.*/
template <typename TElem>
static void SelectSelectedRootElementsAsVSlaves(MGSelector& msel)
{
	UG_DLOG(LG_DIST, 3, "dist-start: SelectSelectedRootElementsAsVSlaves\n");
	GDIST_PROFILE_FUNC();

	using TIter = typename Grid::traits<TElem>::iterator;

	UG_ASSERT(msel.multi_grid(), "Selector has to operate on a MultiGrid");
	MultiGrid& mg = *msel.multi_grid();
	DistributedGridManager& distGridMgr = *mg.distributed_grid_manager();

	for(TIter iter = msel.begin<TElem>(0); iter != msel.end<TElem>(0); ++iter){
		if(!distGridMgr.contains_status(*iter, ElementStatusTypes::ES_V_MASTER))
			msel.select(*iter, InterfaceStates::IS_VSLAVE);
	}
	UG_DLOG(LG_DIST, 3, "dist-stop: SelectSelectedRootElementsAsVSlaves\n");
}

////////////////////////////////////////////////////////////////////////////////
static void SelectElementsForTargetPartition(MGSelector& msel,
								SubsetHandler& shPartition, int partitionIndex,
								bool partitionForLocalProc,
								bool createVerticalInterfaces)
{
	UG_DLOG(LG_DIST, 3, "dist-start: SelectElementsForTargetPartition\n");
	GDIST_PROFILE_FUNC();

//	elements which do not have parents (so called root-elements), and which are
//	not v-slaves have to have a copy on the local proc.
//	If they are not contained in the partition for the local proc, we'll add a
//	copy on the local proc and make it the v-master copy.
//	Note that this should only affect elements in the base-level.
//	The assignment is performed in SelectUnselectedRootElementsAsVMasters.
//todo	assert that only elements in the base-level do not have parents
//		(regarding the global grid)

	UG_ASSERT(msel.multi_grid(), "Selector has to operate on a MultiGrid");
	MultiGrid& mg = *msel.multi_grid();

	if(mg.num<Volume>() > 0){
		if(partitionIndex >= 0)
			SelectSubsetElements<Volume>(msel, shPartition, partitionIndex, InterfaceStates::IS_NORMAL);
		if(createVerticalInterfaces){
			if(partitionForLocalProc)
				SelectUnselectedRootElementsAsVMasters<Volume>(msel);
			else
				SelectSelectedRootElementsAsVSlaves<Volume>(msel);
		}
	}
	else if(mg.num<Face>() > 0){
		if(partitionIndex >= 0)
			SelectSubsetElements<Face>(msel, shPartition, partitionIndex, InterfaceStates::IS_NORMAL);
		if(createVerticalInterfaces){
			if(partitionForLocalProc)
				SelectUnselectedRootElementsAsVMasters<Face>(msel);
			else
				SelectSelectedRootElementsAsVSlaves<Face>(msel);
		}
	}
	else if(mg.num<Edge>() > 0){
		if(partitionIndex >= 0)
			SelectSubsetElements<Edge>(msel, shPartition, partitionIndex, InterfaceStates::IS_NORMAL);
		if(createVerticalInterfaces){
			if(partitionForLocalProc)
				SelectUnselectedRootElementsAsVMasters<Edge>(msel);
			else
				SelectSelectedRootElementsAsVSlaves<Edge>(msel);
		}
	}
	else if(mg.num<Vertex>() > 0){
		if(partitionIndex >= 0)
			SelectSubsetElements<Vertex>(msel, shPartition, partitionIndex, InterfaceStates::IS_NORMAL);
		if(createVerticalInterfaces){
			if(partitionForLocalProc)
				SelectUnselectedRootElementsAsVMasters<Vertex>(msel);
			else
				SelectSelectedRootElementsAsVSlaves<Vertex>(msel);
		}
	}

	if(mg.num<Volume>() > 0){
		if(createVerticalInterfaces)
			AssignVerticalMasterAndSlaveStates<Volume>(msel, partitionForLocalProc);
		AssignSelectionStateToSides<Volume>(msel, true);
	}
	else if(mg.num<Face>() > 0){
		if(createVerticalInterfaces)
			AssignVerticalMasterAndSlaveStates<Face>(msel, partitionForLocalProc);
		AssignSelectionStateToSides<Face>(msel, true);
	}
	else if(mg.num<Edge>() > 0){
		if(createVerticalInterfaces)
			AssignVerticalMasterAndSlaveStates<Edge>(msel, partitionForLocalProc);
		AssignSelectionStateToSides<Edge>(msel, true);
	}
	else if(mg.num<Vertex>() > 0){
		if(createVerticalInterfaces)
			AssignVerticalMasterAndSlaveStates<Vertex>(msel, partitionForLocalProc);
	//	no sides to assign...
	}

	// adjust distribution
	DistributedGridManager& dgm = *mg.distributed_grid_manager();
	SmartPtr<DistroAdjuster> spDA = dgm.distro_adjuster();
	if (spDA.valid()) spDA->adjust(msel, partitionForLocalProc, createVerticalInterfaces);

//	select associated constraining elements first, since they may reference
//	additional unselected constrained elements.
//	UG_LOG("DEBUG: SELECTING CONSTRAINING ELEMENTS...\n");
//	SelectAssociatedConstrainingElements(msel, IS_DUMMY);
	SelectAssociatedConstrainedElements(msel, InterfaceStates::IS_DUMMY | InterfaceStates::HAS_PARENT);
	SelectChildrenOfSelectedShadowRimFaces(msel, InterfaceStates::IS_DUMMY | InterfaceStates::HAS_PARENT);
	SelectChildrenOfSelectedShadowRimEdges(msel, InterfaceStates::IS_DUMMY | InterfaceStates::HAS_PARENT);
	SelectChildrenOfSelectedShadowVertices(msel, InterfaceStates::IS_DUMMY | InterfaceStates::HAS_PARENT);
	UG_DLOG(LG_DIST, 3, "dist-stop: SelectElementsForTargetPartition\n");
}

////////////////////////////////////////////////////////////////////////////////
template <typename TElem>
static void AddTargetProcToDistInfos(MGSelector& msel,
									DistInfoSupplier& distInfos, int targetProc)
{
	UG_DLOG(LG_DIST, 3, "dist-start: AddTargetProcToDistInfos\n");
	GDIST_PROFILE_FUNC();

	using TElemIter = typename Grid::traits<TElem>::iterator;


	for(size_t lvl = 0; lvl < msel.num_levels(); ++lvl){
		for(TElemIter iter = msel.begin<TElem>(lvl);
			iter != msel.end<TElem>(lvl); ++iter)
		{
			TElem* e = *iter;
			byte_t selState = msel.get_selection_status(e);

			distInfos.get(e).push_back(
					TargetProcInfo(targetProc, selState));
		}
	}

	UG_DLOG(LG_DIST, 3, "dist-stop: AddTargetProcToDistInfos\n");
}


////////////////////////////////////////////////////////////////////////////////
///	DistInfos are post-processed and some values are adjusted (primarily missing vslaves-marks are added)
/**	In some situations a copy of an element may be marked as vmaster but some
 * associated copies are neither marked as vmaster or vslave. This would be invalid
 * and we have to mark those copies as vslaves in those situations.
 *
 * This occurs in situations where a low-dim element with copies
 * on p1 and p2 (no v-interface) is distributed from p1 to a third process.
 *
 * This method should be called after the distribution infos have been synchronized.
 * Since it performs the exactly same actions on all processes for synchronized
 * dist-infos, no further communication is required afterwards.
 */
template <typename TElem>
static void PostProcessDistInfos(MultiGrid& mg, DistInfoSupplier& distInfos)
{
//	iterate over all elements and check for each whether a copy is marked as vmaster.
//	If this is the case, all other elements have to be in v-interfaces, too.
//	If a copy isn't in a v-interface, it will be marked as vslave.
	for(typename MultiGrid::traits<TElem>::iterator iter = mg.begin<TElem>();
		iter != mg.end<TElem>(); ++iter)
	{
		TElem* e = *iter;
		vector<TargetProcInfo>& di = distInfos.get(e);
		if(di.size() < 2)
			continue;

		bool gotVMaster = false;
		bool gotNeither = false;
		for(size_t i = 0; i < di.size(); ++i){
			TargetProcInfo& tpi = di[i];
			if(tpi.interfaceState & InterfaceStates::IS_VMASTER){
				gotVMaster = true;
			}
			else if(!(tpi.interfaceState & InterfaceStates::IS_VSLAVE)){
				gotNeither = true;
			}
		}

		if(gotVMaster && gotNeither){
		//	those which have neither a vmaster nor vslave mark have to be marked as vslaves.
			for(size_t i = 0; i < di.size(); ++i){
				TargetProcInfo& tpi = di[i];
				if(!(tpi.interfaceState & (InterfaceStates::IS_VMASTER | InterfaceStates::IS_VSLAVE))){
					tpi.interfaceState |= InterfaceStates::IS_VSLAVE;
				}
			}
		}
	}

	#ifdef LG_DISTRIBUTION_DEBUG
		UG_LOG("DEBUG: DUMMY CHECK\n");
		for(typename MultiGrid::traits<TElem>::iterator iter = mg.begin<TElem>();
			iter != mg.end<TElem>(); ++iter)
		{
			TElem* e = *iter;
			vector<TargetProcInfo>& di = distInfos.get(e);

			if(di.size() < 2)
				continue;

			bool allDummies = true;
			for(size_t i = 0; i < di.size(); ++i){
				TargetProcInfo& tpi = di[i];
				bool isNormal = ((tpi.interfaceState & IS_NORMAL) != 0);
				bool isVMaster = ((tpi.interfaceState & IS_VMASTER) != 0);
				bool isVSlave = ((tpi.interfaceState & IS_VSLAVE) != 0);
				bool isDummy = ((tpi.interfaceState & IS_DUMMY) != 0);
				if(isNormal || isVMaster || isVSlave){
					allDummies = false;
					break;
				}
				else{
					if(!isDummy){
						UG_THROW("Element doesn't have a valid interface state: "
								 << ElementDebugInfo(mg, e));
					}
				}
			}
			if(allDummies){
				UG_THROW("The element (" << ElementDebugInfo(mg, e) << ") has only dummy marks:\n"
						<< distInfos.get_debug_info(e));
			}
		}
	#endif
}


////////////////////////////////////////////////////////////////////////////////
/**	\param partitionIsEmpty	If no elements are selected for a target-partition, the
 * 							corresponding entry is set to false. This can happen even
 * 							if shPartition contains elements for that partition:
 * 							vmaster elements whose parents are not sent to the same
 * 							partition don't have to be sent either and are thus ignored...
 */
static void FillDistInfos(MultiGrid& mg, SubsetHandler& shPartition, MGSelector& msel,
						DistInfoSupplier& distInfos, const std::vector<int>* processMap,
						const pcl::ProcessCommunicator& procComm,
						bool createVerticalInterfaces,
						vector<bool>& partitionIsEmpty)
{
	UG_DLOG(LG_DIST, 3, "dist-start: FillDistInfos\n");
	GDIST_PROFILE_FUNC();

	partitionIsEmpty.resize(shPartition.num_subsets());

	for(int i_part = 0; i_part < shPartition.num_subsets(); ++i_part){

		int targetProc = i_part;
		if(processMap)
			targetProc = (*processMap)[i_part];

		bool localPartition = (targetProc == pcl::ProcRank());

		msel.clear();
		SelectElementsForTargetPartition(msel, shPartition, i_part,
									 localPartition, createVerticalInterfaces);

		partitionIsEmpty[i_part] = msel.empty();

		if(!partitionIsEmpty[i_part]){
		//DEBUG:	temporarily save selection to a file
			#ifdef LG_DISTRIBUTION_DEBUG
			{
				stringstream ss;
				ss << "dist-selection-p" << pcl::ProcRank() << "for-p"<< i_part << ".ugx";
				SaveDistSelectorToFile(msel, ss.str().c_str());
			}
			#endif

			AddTargetProcToDistInfos<Volume>(msel, distInfos, targetProc);
			AddTargetProcToDistInfos<Face>(msel, distInfos, targetProc);
			AddTargetProcToDistInfos<Edge>(msel, distInfos, targetProc);
			AddTargetProcToDistInfos<Vertex>(msel, distInfos, targetProc);
		}
	}

#ifdef LG_DISTRIBUTION_DEBUG
	{
		//stringstream ss;
		//ss << "dist_infos_vrt_before_sync_p" << pcl::ProcRank() << ".ugx";
		//WriteDistInfosToTextFile<Vertex>(mg, distInfos, ss.str().c_str());
		SaveDistInfosToFile(mg, distInfos, "dist_infos_before_sync");
	}
#endif

	SynchronizeDistInfos<Vertex>(mg, distInfos);
	SynchronizeDistInfos<Edge>(mg, distInfos);
	SynchronizeDistInfos<Face>(mg, distInfos);
	SynchronizeDistInfos<Volume>(mg, distInfos);

	PostProcessDistInfos<Vertex>(mg, distInfos);
	PostProcessDistInfos<Edge>(mg, distInfos);
	PostProcessDistInfos<Face>(mg, distInfos);
	PostProcessDistInfos<Volume>(mg, distInfos);

#ifdef LG_DISTRIBUTION_DEBUG
	{
		//stringstream ss;
		//ss << "dist_infos_vrt_after_sync_p" << pcl::ProcRank() << ".ugx";
		//WriteDistInfosToTextFile<Vertex>(mg, distInfos, ss.str().c_str());
		SaveDistInfosToFile(mg, distInfos, "dist_infos_after_sync");
	}
#endif

	UG_DLOG(LG_DIST, 3, "dist-stop: FillDistInfos\n");
}

////////////////////////////////////////////////////////////////////////////////
/**	Based on the list of target processes given by DistInfoSupplier, layouts and
 * interfaces are generated in the given GridLayoutMap.
 *
 * \todo	Think about caching interfaces to speed up this method.
 */
template <typename TElem>
static void CreateLayoutsFromDistInfos(MultiGrid& mg, GridLayoutMap& glm,
										DistInfoSupplier& distInfos,
										AGeomObjID& aGID)
{
	UG_DLOG(LG_DIST, 3, "dist-start: CreateLayoutsFromDistInfos\n");
	GDIST_PROFILE_FUNC();

	using TIter = typename MultiGrid::traits<TElem>::iterator;


	int localProcID = pcl::ProcRank();

	for(size_t lvl = 0; lvl < mg.num_levels(); ++lvl){
		for(TIter iter = mg.begin<TElem>(lvl); iter != mg.end<TElem>(lvl); ++iter)
		{
			TElem* e = *iter;
			vector<TargetProcInfo>& di = distInfos.get(e);

			if(di.size() < 2)
				continue;

		//	get the process with the lowest rank, on which a normal copy of this
		//	element lies (ignore pure vertical masters)
		//	this lowest rank is required to decide, which process a horizontal
		//	master should reside on
			byte_t localInterfaceState = 0;
			int minProc = pcl::NumProcs();
			int minVMasterProc = pcl::NumProcs();
			int minVMasterNoVSlave = pcl::NumProcs();
			int minNormalProc = pcl::NumProcs();

		//	the lowest proc which holds a v-slave or a normal entry.
		//	dummies are ignored here, since we don't want them to be h-masters.
			int minRegularHMasterProc = pcl::NumProcs();
			bool isVMaster = false;
			bool isVSlave = false;
			bool isNormal = false;
			bool isDummy = false;

			bool vMasterExists = false;
			bool dummyExists = false;
			//bool isDummy = false;
			//bool isNormal = false;
			bool createNormalHInterface = false;
//			bool forceHInterface = false;

			int numVSlaveProcs = 0;
			for(size_t i = 0; i < di.size(); ++i){
				TargetProcInfo& tpi = di[i];
				if(tpi.procID == localProcID)
					localInterfaceState = tpi.interfaceState;

				if(tpi.interfaceState & InterfaceStates::IS_VMASTER){
				//	if there is more than one vmaster, then we have to build
				//	h interfaces
//					if(vMasterExists){
//						createNormalHInterface = true;
//						forceHInterface = true;
//					}

					vMasterExists = true;
					if(tpi.procID == localProcID){
						isVMaster = true;
					}
					if(tpi.procID < minVMasterProc)
						minVMasterProc = tpi.procID;
					if(!(tpi.interfaceState & InterfaceStates::IS_VSLAVE)){
						if(tpi.procID < minVMasterNoVSlave)
							minVMasterNoVSlave = tpi.procID;
					}
				}

				if(tpi.interfaceState & InterfaceStates::IS_VSLAVE){
					if(tpi.procID == localProcID)
						isVSlave = true;
					if(tpi.procID < minRegularHMasterProc)
						minRegularHMasterProc = tpi.procID;
					++numVSlaveProcs;
				}

				if(tpi.interfaceState & (InterfaceStates::IS_NORMAL)){
					createNormalHInterface = true;
					if(tpi.procID < minRegularHMasterProc)
						minRegularHMasterProc = tpi.procID;
					if(tpi.procID < minNormalProc)
						minNormalProc = tpi.procID;
					if(tpi.procID == localProcID)
						isNormal = true;
				}

				if(tpi.interfaceState & (InterfaceStates::IS_DUMMY)){
					createNormalHInterface = true;
					dummyExists = true;
					if(tpi.procID == localProcID)
						isDummy = true;

				//	if you don't want to have dummies, which are h-masters, then
				//	remove the following lines
//					if(tpi.procID < minRegularHMasterProc)
//						minRegularHMasterProc = tpi.procID;

				}
				if(tpi.procID < minProc)
					minProc = tpi.procID;
			}

			UG_ASSERT((!createNormalHInterface)
					  || (minRegularHMasterProc < pcl::NumProcs()),
					  "invalid h-master process. The local node (" << ElementDebugInfo(mg, e)
					  << ") has the following flags:\n"
					  << distInfos.get_debug_info(e) << "\n");

		//	if one process is marked as vmaster but not as a vslave, then we have
		//	to be careful if we adjust states on processes which are marked as
		//	both vmaster and vslave.
			if(minVMasterNoVSlave < pcl::NumProcs())
				minVMasterProc = minVMasterNoVSlave;

		//	in some situations, lower dimensional elements can be marked as a
		//	vmaster and as vslave at the same time. We currently allow for elements
		//	to be part of both interfaces at the same time.
			if(isVMaster && isVSlave){
//				if(localProcID == minVMasterProc)
//					isVSlave = false;
//				else
//					isVMaster = false;

//				isVMaster = isVSlave = false;

			//	adjacent normal full-dimensional elements should thus exist and a
			//	horizontal interface has to be built.
				createNormalHInterface = true;
				UG_ASSERT(minRegularHMasterProc < pcl::NumProcs(), "invalid h-master process");
			}
//			else if((!(isVMaster || isVSlave)) && vMasterExists){
//			//	check whether a vmaster copy exists. If this is the case,
//			//	the element itself has to be a vslave.
//			//	This occurs in situations where a low-dim element with copies
//			//	on p1 and p2 (no v-interface) is distributed from p1 to a third process.
//			//	a v-interface on p1 will then be created and
//				isVSlave = true;
//			}

		//	dummies are only required where no normal or slave state is set
			if(isDummy && (isNormal || isVSlave))
				isDummy = false;

		//	there only may be one v-master copy
//			if(isVMaster && (localProcID != minVMasterProc)){
//				isVMaster = false;
//				isVSlave = true;
//			}

		//	if this condition is fulfilled, some kind of h-interface will be built
			//bool createHInterface = createNormalHInterface || (numVSlaveProcs > 1);

			for(size_t i = 0; i < di.size(); ++i){
				TargetProcInfo& tpi = di[i];
				if(tpi.procID == localProcID)
					continue;

				bool tpIsVMaster = (tpi.interfaceState & InterfaceStates::IS_VMASTER);
				bool tpIsVSlave = (tpi.interfaceState & InterfaceStates::IS_VSLAVE);
				bool tpIsNormal = (tpi.interfaceState & InterfaceStates::IS_NORMAL);
				bool tpIsDummy = (tpi.interfaceState & InterfaceStates::IS_DUMMY);


				if(tpIsVMaster && tpIsVSlave){
//					if(tpi.procID == minVMasterProc)
//						tpIsVSlave = false;
//					else
//						tpIsVMaster = false;

//					tpIsVMaster = tpIsVSlave = false;

					createNormalHInterface = true;
					UG_ASSERT(minRegularHMasterProc < pcl::NumProcs(), "invalid h-master process");
				}
//				else if((!(tpIsVMaster || tpIsVSlave)) && vMasterExists){
//					tpIsVSlave = true;
//				}

				if(tpIsDummy && (tpIsNormal || tpIsVSlave))
					tpIsDummy = false;
//				if(tpIsVMaster && (tpi.procID != minVMasterProc)){
//					tpIsVMaster = false;
//					tpIsVSlave = true;
//				}

				bool interfaceCreated = false;

			//	add entry to vertical interface if necessary
				//if(isVSlave && (tpi.procID == minVMasterProc)){
				if(isVSlave && tpIsVMaster){
					glm.get_layout<TElem>(InterfaceNodeTypes::INT_V_SLAVE).
						interface(tpi.procID, lvl).push_back(e);
				}
				if(isVMaster && tpIsVSlave){
					glm.get_layout<TElem>(InterfaceNodeTypes::INT_V_MASTER).
						interface(tpi.procID, lvl).push_back(e);

				//	we have to destroy parent-child relations to possibly existing
				//	children. Those children are now indirectly connected through
				//	copies on other processes.
					//if(!createHInterface && mg.has_children(e))
					if(!createNormalHInterface && mg.has_children(e))
						mg.clear_child_connections(e);
				}
				if(isVSlave && tpIsVSlave){
					UG_ASSERT(minRegularHMasterProc < pcl::NumProcs(), "invalid h-master process");
				//	we still have to build a horizontal interface, this time
				//	however only between vertical slaves
//					if(tpIsVSlave && (!tpWasVMaster)){
//						if(!(isVMaster || wasVMaster)){
//					if(isVSlave && tpIsVSlave){
//					if(tpIsVSlave){
//						if(!(isVMaster)){

					if(localProcID == minRegularHMasterProc){
					//	horizontal master
						interfaceCreated = true;
						glm.get_layout<TElem>(InterfaceNodeTypes::INT_H_MASTER).
							interface(tpi.procID, lvl).push_back(e);
					}
					else if(tpi.procID == minRegularHMasterProc){
					//	horizontal slave
						interfaceCreated = true;
						glm.get_layout<TElem>(InterfaceNodeTypes::INT_H_SLAVE).
							interface(tpi.procID, lvl).push_back(e);
					}

//						}
//					}
				}

			//	add entry to horizontal interface if necessary
				if(!interfaceCreated && createNormalHInterface){
					UG_ASSERT(minRegularHMasterProc < pcl::NumProcs(), "invalid h-master process");

				//	check whether the target process would also create a normal h interface
					if(localProcID == minRegularHMasterProc){
					//	horizontal master
					//	only build the interface if the process is not a pure
					//	v-master
						if(tpi.interfaceState != InterfaceStates::IS_VMASTER){
						//if((tpi.procID != minVMasterProc) || (tpi.interfaceState != IS_VMASTER)){
						//if(forceHInterface || (tpi.interfaceState != IS_VMASTER)){
							glm.get_layout<TElem>(InterfaceNodeTypes::INT_H_MASTER).
								interface(tpi.procID, lvl).push_back(e);
						}
					}
					else if(tpi.procID == minRegularHMasterProc){
					//	horizontal slave
					//	only build the interface if the process is not a pure
					//	v-master
						if(localInterfaceState != InterfaceStates::IS_VMASTER){
						//if((localProcID != minVMasterProc) || (localInterfaceState != IS_VMASTER)){
						//if(forceHInterface || (localInterfaceState != IS_VMASTER)){
							glm.get_layout<TElem>(InterfaceNodeTypes::INT_H_SLAVE).
								interface(tpi.procID, lvl).push_back(e);
						}
					}
				}

			//	finally we have to make sure, that dummies which do not have a parent
			//	are v-slaves.
				if(dummyExists && (!vMasterExists)){
				//	the lowest normal process will be transformed to a v-master
					UG_ASSERT(minNormalProc < pcl::NumProcs(), "invalid minNormalProc!");

					if(isDummy && (!(localInterfaceState & InterfaceStates::HAS_PARENT))
						&& (tpi.procID == minNormalProc))
					{
						glm.get_layout<TElem>(InterfaceNodeTypes::INT_V_SLAVE).
								interface(tpi.procID, lvl).push_back(e);
					}
					else if(tpIsDummy && (!(tpi.interfaceState & InterfaceStates::HAS_PARENT))
							&& (localProcID == minNormalProc))
					{
						glm.get_layout<TElem>(InterfaceNodeTypes::INT_V_MASTER).
								interface(tpi.procID, lvl).push_back(e);
					}
				}
			}
		}
	}

//	Now sort the interface entries in the different layouts
	CompareByAttachment<TElem, AGeomObjID> gidCmp(mg, aGID);
	if(glm.has_layout<TElem>(InterfaceNodeTypes::INT_H_MASTER))
		glm.get_layout<TElem>(InterfaceNodeTypes::INT_H_MASTER).sort_interface_entries(gidCmp);
	if(glm.has_layout<TElem>(InterfaceNodeTypes::INT_H_SLAVE))
		glm.get_layout<TElem>(InterfaceNodeTypes::INT_H_SLAVE).sort_interface_entries(gidCmp);
	if(glm.has_layout<TElem>(InterfaceNodeTypes::INT_V_MASTER))
		glm.get_layout<TElem>(InterfaceNodeTypes::INT_V_MASTER).sort_interface_entries(gidCmp);
	if(glm.has_layout<TElem>(InterfaceNodeTypes::INT_V_SLAVE))
		glm.get_layout<TElem>(InterfaceNodeTypes::INT_V_SLAVE).sort_interface_entries(gidCmp);

	UG_DLOG(LG_DIST, 3, "dist-stop: CreateLayoutsFromDistInfos\n");
}

///	Adds serializers for all registered global attachments.
/**	Make sure that the same global attachments are attached to the given grid
 * on all processes before calling this method. Use 'SynchronizeAttachedGlobalAttachments'
 * to achieve this.*/
template <typename TElem>
static
void AddGlobalAttachmentsToSerializer (
		GridDataSerializationHandler& handler,
		Grid& grid)
{
	const vector<string>& attachmentNames = GlobalAttachments::declared_attachment_names();
	for(size_t i = 0; i < attachmentNames.size(); ++i){
		const string& name = attachmentNames[i];
		if(GlobalAttachments::is_attached<TElem>(grid, name)){
			GlobalAttachments::add_data_serializer<TElem>(handler, grid, name);
		}
	}
}

/**	Attaches global attachments to 'g' that are attached at 'g' on some other process,
 * thus synchronizing the set of attached global attachments of 'g'.*/
static void SynchronizeAttachedGlobalAttachments (
				Grid& g,
				const pcl::ProcessCommunicator& procComm)
{
	UG_DLOG(LG_DIST, 1, "SynchronizeAttachedGlobalAttachments start\n");

//	now make sure that the same global attachments are attached everywhere
	const vector<string>&  names = GlobalAttachments::declared_attachment_names();

//	make sure that the same number of attachments is declared on all processes
	const size_t maxNumAttachmentNames = procComm.allreduce(names.size(), PCL_RO_MAX);
	UG_COND_THROW(pcl::OneProcTrue(maxNumAttachmentNames != names.size(), procComm),
				  "Different number of global attachments declared on different processes!");

	if(maxNumAttachmentNames == 0)
		return;

	vector<string>	attachedNamesAndTypes;

	vector<byte_t> locAttached(names.size(), 0);
	for(size_t i = 0; i < names.size(); ++i){
		byte_t& b = locAttached[i];
		if(GlobalAttachments::is_attached<Vertex>(g, names[i]))
			b |= 1;
		if(GlobalAttachments::is_attached<Edge>(g, names[i]))
			b |= 1<<1;
		if(GlobalAttachments::is_attached<Face>(g, names[i]))
			b |= 1<<2;
		if(GlobalAttachments::is_attached<Volume>(g, names[i]))
			b |= 1<<3;
	}

	vector<byte_t> globAttached(names.size());
	procComm.allreduce(locAttached, globAttached, PCL_RO_BOR);

	for(size_t i = 0; i < names.size(); ++i){
		byte_t& b = globAttached[i];
		if((b & 1) && !GlobalAttachments::is_attached<Vertex>(g, names[i]))
			GlobalAttachments::attach<Vertex>(g, names[i]);

		if((b & 1<<1) && !GlobalAttachments::is_attached<Edge>(g, names[i]))
			GlobalAttachments::attach<Edge>(g, names[i]);

		if((b & 1<<2) && !GlobalAttachments::is_attached<Face>(g, names[i]))
			GlobalAttachments::attach<Face>(g, names[i]);

		if((b & 1<<3) && !GlobalAttachments::is_attached<Volume>(g, names[i]))
			GlobalAttachments::attach<Volume>(g, names[i]);
	}
	UG_DLOG(LG_DIST, 1, "SynchronizeAttachedGlobalAttachments end\n");
}


////////////////////////////////////////////////////////////////////////////////
bool DistributeGrid(MultiGrid& mg,
					SubsetHandler& shPartition,
					GridDataSerializationHandler& serializer,
					bool createVerticalInterfaces,
					const std::vector<int>* processMap,
					const pcl::ProcessCommunicator& procComm)
{
	GDIST_PROFILE_FUNC();
	PCL_DEBUG_BARRIER(procComm);
	GDIST_PROFILE(performRedistribution);
	UG_STATIC_ASSERT(IS_DUMMY < 256, RedistributeGrid_IS_DUMMY_too_big);

	UG_DLOG(LG_DIST, 3, "dist-start: DistributeGrid\n");
	const char* errprefix = "ERROR in DistributeGrid: ";

	if(!mg.is_parallel()){
		UG_THROW(errprefix << "Can't distribute a serial grid! Compile ug with -DPARALLEL=ON");
	}


	UG_DLOG(LG_DIST, 2, "dist: Informing msg-hub that distribution starts\n");
	GDIST_PROFILE(gdist_distStartsCallback);
	GridDataSerializationHandler	userDataSerializer;
//	add global attachments to the user-data-serializer
	SynchronizeAttachedGlobalAttachments(mg, procComm);
	AddGlobalAttachmentsToSerializer<Vertex>(userDataSerializer, mg);
	AddGlobalAttachmentsToSerializer<Edge>(userDataSerializer, mg);
	AddGlobalAttachmentsToSerializer<Face>(userDataSerializer, mg);
	AddGlobalAttachmentsToSerializer<Volume>(userDataSerializer, mg);

	SPMessageHub msgHub = mg.message_hub();
	msgHub->post_message(GridMessage_Distribution(GridMessageDistributionType::GMDT_DISTRIBUTION_STARTS, userDataSerializer));

	PCL_DEBUG_BARRIER(procComm);
	GDIST_PROFILE_END();


	DistributedGridManager& distGridMgr = *mg.distributed_grid_manager();
	GridLayoutMap& glm = distGridMgr.grid_layout_map();

//	The selector will be of frequent use to speed up some algorithms
	MGSelector msel(mg);

	#ifdef LG_DISTRIBUTION_DEBUG
		PerformValidityCheck(distGridMgr);
	#endif

//	Since we will change huge parts of the underlying grid and the grid-layout-map,
//	we'll disable auto-insertion of elements in the distributed-grid-manager.
//	This means we carefully have to take care of all interface changes.
	distGridMgr.enable_interface_management(false);

////////////////////////////////
//	GLOBAL IDS
//todo:	only create global ids if they aren't already present
	GDIST_PROFILE(gdist_CreateGlobalIDs);
	UG_DLOG(LG_DIST, 2, "dist-DistributeGrid: Create global vertex ids\n");
	CreateAndDistributeGlobalIDs<Vertex>(mg, glm);
	CreateAndDistributeGlobalIDs<Edge>(mg, glm);
	CreateAndDistributeGlobalIDs<Face>(mg, glm);
	CreateAndDistributeGlobalIDs<Volume>(mg, glm);
	MultiElementAttachmentAccessor<AGeomObjID> aaID(mg, aGeomObjID);
	PCL_DEBUG_BARRIER(procComm);
	GDIST_PROFILE_END();


	#ifdef LG_DISTRIBUTION_DEBUG
	{
		UG_LOG("DEBUG: WRITING GLOBAL VERTEX IDS TO FILE\n");
		stringstream ss;
		ss << "global_ids_vrt_p" << pcl::ProcRank() << ".txt";
		WriteDebugValuesToFile<Vertex>(ss.str().c_str(), mg, aGeomObjID, false);
	}
	{
		UG_LOG("DEBUG: WRITING GLOBAL EDGE IDS TO FILE\n");
		stringstream ss;
		ss << "global_ids_edge_p" << pcl::ProcRank() << ".txt";
		WriteDebugValuesToFile<Edge>(ss.str().c_str(), mg, aGeomObjID, false);
	}
	{
		UG_LOG("DEBUG: WRITING GLOBAL FACE IDS TO FILE\n");
		stringstream ss;
		ss << "global_ids_face_p" << pcl::ProcRank() << ".txt";
		WriteDebugValuesToFile<Face>(ss.str().c_str(), mg, aGeomObjID, false);
	}
	#endif

////////////////////////////////
//	FILL THE DISTRIBUTION INFOS (INVOLVES COMMUNICATION...)
	GDIST_PROFILE(gdist_FillDistInfos);
	UG_DLOG(LG_DIST, 2, "dist-DistributeGrid: Fill distribution infos\n");
	vector<bool> partitionIsEmpty;
	DistInfoSupplier distInfos(mg);
	FillDistInfos(mg, shPartition, msel, distInfos, processMap, procComm,
				  createVerticalInterfaces, partitionIsEmpty);
	PCL_DEBUG_BARRIER(procComm);
	GDIST_PROFILE_END();

//	DEBUG: output distInfos...
	#ifdef LG_DISTRIBUTION_DEBUG
	{
		SaveDistInfosToFile(mg, distInfos, "dist_infos_before_distribution");
	}
	#endif

////////////////////////////////
//	COMMUNICATE INVOLVED PROCESSES
	GDIST_PROFILE(gdist_CommunicateInvolvedProcs);
	UG_DLOG(LG_DIST, 2, "dist-DistributeGrid: CommunicateInvolvedProcesses\n");

//	each process has to know with which other processes it
//	has to communicate.
	vector<int> sendToRanks, recvFromRanks, sendPartitionInds;

	if(processMap && (shPartition.num_subsets() > (int)processMap->size())){
		UG_THROW("process-map is too small for the given number of partitions!");
	}

//	for each subset which is not emtpy we'll have to send data to
//	the associated process.
	for(int si = 0; si < shPartition.num_subsets(); ++si){
	//	instead of simply querying shPartition.empty(si), we'll check partitionIsEmpty[si],
	//	since this array tells whether data is actually sent to a partition.
	//	E.g. vmasters which are contained in shPartition are not necessarily sent
	//	to the associated target process...
		if(!partitionIsEmpty[si]){
			int toProc = si;
		//	if a process map exists, we'll use the associated process
			if(processMap)
				toProc = processMap->at(si);

			sendToRanks.push_back(toProc);
			sendPartitionInds.push_back(si);
		}
	}

	pcl::CommunicateInvolvedProcesses(recvFromRanks, sendToRanks, procComm);

	PCL_DEBUG_BARRIER(procComm);
	GDIST_PROFILE_END();


////////////////////////////////
//	SERIALIZE THE GRID, THE GLOBAL IDS AND THE DISTRIBUTION INFOS.
	GDIST_PROFILE(gdist_Serialization);
	UG_DLOG(LG_DIST, 2, "dist-DistributeGrid: Serialization\n");
	AInt aLocalInd("distribution-tmp-local-index");
	mg.attach_to_all(aLocalInd);
	MultiElementAttachmentAccessor<AInt> aaInt(mg, aLocalInd);

//	outBufs will be used to serialize and distribute the grid.
//	don't resize outBufs later on! Would be expensive!
	std::vector<BinaryBuffer> outBufs(sendToRanks.size());

//	the magic number is used for debugging to make sure that the stream is read correctly
	int magicNumber1 = 75234587;
	int magicNumber2 = 560245;

	ADistInfo aDistInfo = distInfos.dist_info_attachment();

	GridDataSerializationHandler distInfoSerializer;
	distInfoSerializer.add(GeomObjAttachmentSerializer<Vertex, ADistInfo>::create(mg, aDistInfo));
	distInfoSerializer.add(GeomObjAttachmentSerializer<Edge, ADistInfo>::create(mg, aDistInfo));
	distInfoSerializer.add(GeomObjAttachmentSerializer<Face, ADistInfo>::create(mg, aDistInfo));
	distInfoSerializer.add(GeomObjAttachmentSerializer<Volume, ADistInfo>::create(mg, aDistInfo));

//	now perform the serialization
	int localPartitionInd = -1;
	for(size_t i_to = 0; i_to < sendPartitionInds.size(); ++i_to){
		int partInd = sendPartitionInds[i_to];
		bool localPartition = (sendToRanks[i_to] == pcl::ProcRank());
		if(localPartition)
			localPartitionInd = partInd;

	//	don't serialize the local partition since we'll keep it here on the local
	//	process anyway.
		if(!localPartition){
			BinaryBuffer& out = outBufs[i_to];

		//	write a magic number for debugging purposes
			out.write((char*)&magicNumber1, sizeof(int));

		//	select the elements of the current partition
			msel.clear();
			SelectElementsForTargetPartition(msel, shPartition, partInd,
										 localPartition, createVerticalInterfaces);
			//AdjustGhostSelection(msel, ISelector::DESELECTED);

			SerializeMultiGridElements(mg, msel.get_grid_objects(), aaInt, out, &aaID);


		//	serialize associated data
			distInfoSerializer.write_infos(out);
			distInfoSerializer.serialize(out, msel.get_grid_objects());
			serializer.write_infos(out);
			serializer.serialize(out, msel.get_grid_objects());
			userDataSerializer.write_infos(out);
			userDataSerializer.serialize(out, msel.get_grid_objects());

		//	write a magic number for debugging purposes
			out.write((char*)&magicNumber2, sizeof(int));
		}
	}
	PCL_DEBUG_BARRIER(procComm);
	GDIST_PROFILE_END();



////////////////////////////////
//	COMMUNICATE SERIALIZED DATA
	GDIST_PROFILE(gdist_CommunicateSerializedData);
	UG_DLOG(LG_DIST, 2, "dist-DistributeGrid: Distribute data\n");
//	now distribute the packs between involved processes
	std::vector<BinaryBuffer> inBufs(recvFromRanks.size());

	procComm.distribute_data(GetDataPtr(inBufs), GetDataPtr(recvFromRanks),
							(int)recvFromRanks.size(),
							GetDataPtr(outBufs), GetDataPtr(sendToRanks),
							(int)sendToRanks.size());

//	clear out-buffers, since they are no longer needed
	for(size_t i = 0; i < outBufs.size(); ++i)
		outBufs[i] = BinaryBuffer();

	PCL_DEBUG_BARRIER(procComm);
	GDIST_PROFILE_END();


////////////////////////////////
//	INTERMEDIATE CLEANUP
	GDIST_PROFILE(gdist_IntermediateCleanup);
	UG_DLOG(LG_DIST, 2, "dist-DistributeGrid: Intermediate cleanup\n");

	msgHub->post_message(GridMessage_Creation(GridMessageCreationType::GMCT_CREATION_STARTS));

//	we have to remove all elements which won't stay on the local process.
//	To do so, we'll first select all elements that stay, invert that selection
//	and erase all elements which are selected thereafter.
	if(createVerticalInterfaces || (localPartitionInd != -1)){
		msel.clear();
		SelectElementsForTargetPartition(msel, shPartition, localPartitionInd,
									 	 true, createVerticalInterfaces);
		InvertSelection(msel);

	//	make sure that constrained/constraining connections won't be harmed
	//	this is a little cumbersome in the moment. Ideally constrained/constraining
	//	elements should unregister from each other automatically on destruction.
		// for(size_t lvl = 0; lvl < msel.num_levels(); ++lvl){
		// 	for(ConstrainedVertexIterator iter = msel.begin<ConstrainedVertex>(lvl);
		// 		iter != msel.end<ConstrainedVertex>(lvl); ++iter)
		// 	{
		// 		GridObject* co = (*iter)->get_constraining_object();
		// 		if(co && !msel.is_selected(co)){
		// 			switch(co->base_object_id()){
		// 				case EDGE:{
		// 					if(ConstrainingEdge* ce = dynamic_cast<ConstrainingEdge*>(co))
		// 						ce->unconstrain_object(*iter);
		// 				}break;
		// 				case FACE:{
		// 					if(co->reference_object_id() == ROID_TRIANGLE){
		// 						if(ConstrainingTriangle* ce = dynamic_cast<ConstrainingTriangle*>(co))
		// 							ce->unconstrain_object(*iter);
		// 					}
		// 					else{
		// 						if(ConstrainingQuadrilateral* ce = dynamic_cast<ConstrainingQuadrilateral*>(co))
		// 							ce->unconstrain_object(*iter);
		// 					}
		// 				}break;
		// 				default: break;
		// 			}
		// 		}
		// 	}

		// 	for(ConstrainedEdgeIterator iter = msel.begin<ConstrainedEdge>(lvl);
		// 		iter != msel.end<ConstrainedEdge>(lvl); ++iter)
		// 	{
		// 		GridObject* co = (*iter)->get_constraining_object();
		// 		if(co && !msel.is_selected(co)){
		// 			switch(co->base_object_id()){
		// 				case EDGE:{
		// 					if(ConstrainingEdge* ce = dynamic_cast<ConstrainingEdge*>(co))
		// 						ce->unconstrain_object(*iter);
		// 				}break;
		// 				case FACE:{
		// 					if(co->reference_object_id() == ROID_TRIANGLE){
		// 						if(ConstrainingTriangle* ce = dynamic_cast<ConstrainingTriangle*>(co))
		// 							ce->unconstrain_object(*iter);
		// 					}
		// 					else{
		// 						if(ConstrainingQuadrilateral* ce = dynamic_cast<ConstrainingQuadrilateral*>(co))
		// 							ce->unconstrain_object(*iter);
		// 					}
		// 				}break;
		// 				default: break;
		// 			}
		// 		}
		// 	}

		// 	for(ConstrainingEdgeIterator iter = msel.begin<ConstrainingEdge>(lvl);
		// 		iter != msel.end<ConstrainingEdge>(lvl); ++iter)
		// 	{
		// 		ConstrainingEdge* e = *iter;
		// 		for(size_t i = 0; i < e->num_constrained_vertices(); ++i){
		// 			ConstrainedVertex* cv = dynamic_cast<ConstrainedVertex*>(e->constrained_vertex(i));
		// 			UG_ASSERT(cv, "Constrained vertices have to be of the type ConstrainedVertex");
		// 			cv->set_constraining_object(nullptr);
		// 		}

		// 		for(size_t i = 0; i < e->num_constrained_edges(); ++i){
		// 			ConstrainedEdge* cde = dynamic_cast<ConstrainedEdge*>(e->constrained_edge(i));
		// 			UG_ASSERT(cde, "Constrained edges have to be of the type ConstrainedEdge");
		// 			cde->set_constraining_object(nullptr);
		// 		}
		// 	}


		// 	for(ConstrainedTriangleIterator iter = msel.begin<ConstrainedTriangle>(lvl);
		// 		iter != msel.end<ConstrainedTriangle>(lvl); ++iter)
		// 	{
		// 		GridObject* co = (*iter)->get_constraining_object();
		// 		if(co && !msel.is_selected(co)){
		// 			if(ConstrainingTriangle* ce = dynamic_cast<ConstrainingTriangle*>(co))
		// 				ce->unconstrain_object(*iter);
		// 		}
		// 	}

		// 	for(ConstrainedQuadrilateralIterator iter = msel.begin<ConstrainedQuadrilateral>(lvl);
		// 		iter != msel.end<ConstrainedQuadrilateral>(lvl); ++iter)
		// 	{
		// 		GridObject* co = (*iter)->get_constraining_object();
		// 		if(co && !msel.is_selected(co)){
		// 			if(ConstrainingQuadrilateral* ce = dynamic_cast<ConstrainingQuadrilateral*>(co))
		// 				ce->unconstrain_object(*iter);
		// 		}
		// 	}

		// 	for(ConstrainingTriangleIterator iter = msel.begin<ConstrainingTriangle>(lvl);
		// 		iter != msel.end<ConstrainingTriangle>(lvl); ++iter)
		// 	{
		// 		ConstrainingFace* e = *iter;
		// 		for(size_t i = 0; i < e->num_constrained_vertices(); ++i){
		// 			ConstrainedVertex* cv = dynamic_cast<ConstrainedVertex*>(e->constrained_vertex(i));
		// 			UG_ASSERT(cv, "Constrained vertices have to be of the type ConstrainedVertex");
		// 			cv->set_constraining_object(nullptr);
		// 		}

		// 		for(size_t i = 0; i < e->num_constrained_edges(); ++i){
		// 			ConstrainedEdge* cde = dynamic_cast<ConstrainedEdge*>(e->constrained_edge(i));
		// 			UG_ASSERT(cde, "Constrained edges have to be of the type ConstrainedEdge");
		// 			cde->set_constraining_object(nullptr);
		// 		}

		// 		for(size_t i = 0; i < e->num_constrained_faces(); ++i){
		// 			ConstrainedFace* cdf = dynamic_cast<ConstrainedFace*>(e->constrained_face(i));
		// 			UG_ASSERT(cdf, "Constrained faces have to be of the type ConstrainedFace");
		// 			cdf->set_constraining_object(nullptr);
		// 		}
		// 	}

		// 	for(ConstrainingQuadrilateralIterator iter = msel.begin<ConstrainingQuadrilateral>(lvl);
		// 		iter != msel.end<ConstrainingQuadrilateral>(lvl); ++iter)
		// 	{
		// 		ConstrainingFace* e = *iter;
		// 		for(size_t i = 0; i < e->num_constrained_vertices(); ++i){
		// 			ConstrainedVertex* cv = dynamic_cast<ConstrainedVertex*>(e->constrained_vertex(i));
		// 			UG_ASSERT(cv, "Constrained vertices have to be of the type ConstrainedVertex");
		// 			cv->set_constraining_object(nullptr);
		// 		}

		// 		for(size_t i = 0; i < e->num_constrained_edges(); ++i){
		// 			ConstrainedEdge* cde = dynamic_cast<ConstrainedEdge*>(e->constrained_edge(i));
		// 			UG_ASSERT(cde, "Constrained edges have to be of the type ConstrainedEdge");
		// 			cde->set_constraining_object(nullptr);
		// 		}

		// 		for(size_t i = 0; i < e->num_constrained_faces(); ++i){
		// 			ConstrainedFace* cdf = dynamic_cast<ConstrainedFace*>(e->constrained_face(i));
		// 			UG_ASSERT(cdf, "Constrained faces have to be of the type ConstrainedFace");
		// 			cdf->set_constraining_object(nullptr);
		// 		}
		// 	}
		// }

		GDIST_PROFILE(gdist_ErasingObjects);
		EraseSelectedObjects(msel);
		GDIST_PROFILE_END();
	}
	else{
	//	nothing remains on the local process...
		GDIST_PROFILE(gdist_ClearGeometry);
		mg.clear_geometry();
		GDIST_PROFILE_END();
	}

	{
		GDIST_PROFILE(gdist_ClearLayoutMap);
	//	the grid layout map will be rebuilt from scratch
		glm.clear();
		GDIST_PROFILE_END();
	}
	PCL_DEBUG_BARRIER(procComm);
	GDIST_PROFILE_END();

//	DEBUGGING...
	// {
	// 	static int counter = 0;
	// 	stringstream ss;
	// 	ss << "parallel-grid-layout-before-redist-(cleared)" << counter << "-p" << pcl::ProcRank() << ".ugx";
	// 	UG_LOG("DEBUG SAVE OF PARALLEL GRID LAYOUT IN DistributeGrid\n");
	// 	SaveParallelGridLayout(mg, ss.str().c_str(), 2);
	// 	++counter;
	// }

////////////////////////////////
//	DESERIALIZE INCOMING GRIDS
	GDIST_PROFILE(gdist_Deserialize);
	distInfoSerializer.deserialization_starts();
	serializer.deserialization_starts();
	userDataSerializer.deserialization_starts();

	vector<Vertex*>	vrts;
	vector<Edge*> edges;
	vector<Face*> faces;
	vector<Volume*> vols;

	for(size_t i = 0; i < recvFromRanks.size(); ++i){
	//	there is nothing to serialize from the local rank
		if(recvFromRanks[i] == pcl::ProcRank())
			continue;

		BinaryBuffer& in = inBufs[i];

		UG_DLOG(LG_DIST, 2, "Deserializing from rank " << recvFromRanks[i] << "\n");

	//	read the magic number and make sure that it matches our magicNumber
		int tmp = 0;
		in.read((char*)&tmp, sizeof(int));
		if(tmp != magicNumber1){
			UG_THROW("ERROR in RedistributeGrid: "
					 "Magic number mismatch before deserialization.\n");
		}

		DeserializeMultiGridElements(mg, in, &vrts, &edges, &faces, &vols, &aaID);

	//	deserialize the associated data (global ids have already been deserialized)
		distInfoSerializer.read_infos(in);
		distInfoSerializer.deserialize(in, vrts.begin(), vrts.end());
		distInfoSerializer.deserialize(in, edges.begin(), edges.end());
		distInfoSerializer.deserialize(in, faces.begin(), faces.end());
		distInfoSerializer.deserialize(in, vols.begin(), vols.end());

		serializer.read_infos(in);
		serializer.deserialize(in, vrts.begin(), vrts.end());
		serializer.deserialize(in, edges.begin(), edges.end());
		serializer.deserialize(in, faces.begin(), faces.end());
		serializer.deserialize(in, vols.begin(), vols.end());

		userDataSerializer.read_infos(in);
		userDataSerializer.deserialize(in, vrts.begin(), vrts.end());
		userDataSerializer.deserialize(in, edges.begin(), edges.end());
		userDataSerializer.deserialize(in, faces.begin(), faces.end());
		userDataSerializer.deserialize(in, vols.begin(), vols.end());

	//	read the magic number and make sure that it matches our magicNumber
		tmp = 0;
		in.read((char*)&tmp, sizeof(int));
		if(tmp != magicNumber2){
			UG_THROW("ERROR in RedistributeGrid: "
					 "Magic number mismatch after deserialization.\n");
		}

		UG_DLOG(LG_DIST, 2, "Deserialization from rank " << recvFromRanks[i] << " done\n");
	}

	//	clear in-buffers, since they are no longer needed
	for(size_t i = 0; i < inBufs.size(); ++i)
		inBufs[i] = BinaryBuffer();

	PCL_DEBUG_BARRIER(procComm);
	GDIST_PROFILE_END();

//	DEBUG: output distInfos...
	#ifdef LG_DISTRIBUTION_DEBUG
	{
		SaveDistInfosToFile(mg, distInfos, "dist_infos_after_distribution");
	}
	#endif

////////////////////////////////
//	CREATE LAYOUTS
	GDIST_PROFILE(gdist_CreateLayouts);
	CreateLayoutsFromDistInfos<Vertex>(mg, glm, distInfos, aGeomObjID);
	CreateLayoutsFromDistInfos<Edge>(mg, glm, distInfos, aGeomObjID);
	CreateLayoutsFromDistInfos<Face>(mg, glm, distInfos, aGeomObjID);
	CreateLayoutsFromDistInfos<Volume>(mg, glm, distInfos, aGeomObjID);
	PCL_DEBUG_BARRIER(procComm);
	GDIST_PROFILE_END();

////////////////////////////////
//	UPDATE THE DISTRIBUTED GRID MANAGER
	GDIST_PROFILE(gdist_UpdateDistGridManager);
	UG_DLOG(LG_DIST, 2, "dist-DistributeGrid: Update DistributedGridManager\n");
	glm.remove_empty_interfaces();
	distGridMgr.enable_interface_management(true);
	distGridMgr.grid_layouts_changed(false);
	PCL_DEBUG_BARRIER(procComm);
	GDIST_PROFILE_END();

	mg.detach_from_all(aLocalInd);

	#ifdef LG_DISTRIBUTION_DEBUG
		PerformValidityCheck(distGridMgr);
	#endif

//	DEBUGGING...
	// {
	// 	static int counter = 0;
	// 	stringstream ss;
	// 	ss << "parallel-grid-layout-after-redist-" << counter << "-p" << pcl::ProcRank() << ".ugx";
	// 	UG_LOG("DEBUG SAVE OF PARALLEL GRID LAYOUT IN DistributeGrid\n");
	// 	SaveParallelGridLayout(mg, ss.str().c_str(), 2);
	// 	++counter;

	// 	if(!TestGridLayoutMap(mg, glm)){
	// 		UG_THROW("TestGridLayoutMap failed after redistribution!");
	// 	}
	// }


	GDIST_PROFILE_END_(performRedistribution);

//	execute callbacks for external postprocessing
	GDIST_PROFILE(gdist_ExternalPostProcessing);
	UG_DLOG(LG_DIST, 2, "dist: Informing msg-hub that distribution stops\n");

	msgHub->post_message(GridMessage_Creation(GridMessageCreationType::GMCT_CREATION_STOPS));
	//msgHub->post_message(GridMessage_Distribution(GMDT_GRID_SERIALIZATION_DONE));

//	we'll inform deserializers now, that deserialization is complete.
	distInfoSerializer.deserialization_done();
	serializer.deserialization_done();
	userDataSerializer.deserialization_done();
	msgHub->post_message(GridMessage_Distribution(GridMessageDistributionType::GMDT_DISTRIBUTION_STOPS, userDataSerializer));
	//msgHub->post_message(GridMessage_Distribution(GMDT_DATA_SERIALIZATION_DONE));
	PCL_DEBUG_BARRIER(procComm);
	GDIST_PROFILE_END();

	UG_DLOG(LG_DIST, 3, "dist-stop: DistributeGrid\n");
	return true;
}

}// end of namespace
