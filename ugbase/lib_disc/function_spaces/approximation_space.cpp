/*
 * approximation_space.cpp
 *
 *  Created on: 13.12.2011
 *      Author: andreasvogel
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
		m_spSurfaceView = SmartPtr<SurfaceView>(NULL);
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

ConstSmartPtr<DoFDistribution>
IApproximationSpace::dof_distribution(const GridLevel& gl, bool bCreate) const
{
	return const_cast<IApproximationSpace*>(this)->dof_distribution(gl, bCreate);
}

void IApproximationSpace::init_levels()
{
	PROFILE_FUNC();
	for(size_t lvl = 0; lvl < num_levels(); ++lvl)
		create_dof_distribution(GridLevel(lvl, GridLevel::LEVEL, true));
}

void IApproximationSpace::init_surfaces()
{
	PROFILE_FUNC();
	for(size_t lvl = 0; lvl < num_levels(); ++lvl)
		create_dof_distribution(GridLevel(lvl, GridLevel::SURFACE, false));

	init_top_surface();
}

void IApproximationSpace::init_top_surface()
{
	PROFILE_FUNC();
	create_dof_distribution(GridLevel(GridLevel::TOPLEVEL, GridLevel::SURFACE, false));
}

////////////////////////////////////////////////////////////////////////////////
// DoFDistribution Creation
////////////////////////////////////////////////////////////////////////////////

void IApproximationSpace::create_dof_distribution(const GridLevel& gl)
{

	dof_distribution_info_required();
	surface_view_required();

	// todo: reuse DoFIndexStorage in case of Level
	SmartPtr<DoFDistribution> spDD = SmartPtr<DoFDistribution>(new
		DoFDistribution(m_spMG, m_spMGSH, m_spDoFDistributionInfo,
						m_spSurfaceView, gl, m_bGrouped));

	m_vDD.push_back(spDD);
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
						ss << "grid-changed-surface-view" << counter << "-p" << pcl::GetProcRank() << ".ugx";
						UG_LOG("PERFORMING SURFACE VIEW DEBUG SAVE IN IApproximationSpace::grid_changed_callback: " << ss.str() << "\n");
						SaveSurfaceViewTransformed(*m_spMG, *m_spSurfaceView, ss.str().c_str(), 0.1);
						++counter;
					}
					{
						#ifdef UG_PARALLEL
							static int counter = 0;
							std::stringstream ss;
							ss << "grid-changed-parallel-layout-" << counter << "-p" << pcl::GetProcRank() << ".ugx";
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
					ss << "grid-distributed-surface-view" << counter << "-p" << pcl::GetProcRank() << ".ugx";
					UG_LOG("PERFORMING SURFACE VIEW DEBUG SAVE IN IApproximationSpace::grid_distribution_callback: " << ss.str() << "\n");
					SaveSurfaceViewTransformed(*m_spMG, *m_spSurfaceView, ss.str().c_str(), 0.1);
					++counter;
				}
				{
					#ifdef UG_PARALLEL
						static int counter = 0;
						std::stringstream ss;
						ss << "grid-distributed-parallel-layout-" << counter << "-p" << pcl::GetProcRank() << ".ugx";
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

void IApproximationSpace::
print_statistic(ConstSmartPtr<DoFDistribution> dd, int verboseLev) const
{
//	Total number of DoFs
	UG_LOG(std::setw(10) << ConvertNumber(dd->num_indices(),10,6) << " | ");

	const int blockSize = DefaultAlgebra::get().blocksize();

//	Overall block size
	if(blockSize != AlgebraType::VariableBlockSize) {UG_LOG(std::setw(8)  << blockSize);}
	else {UG_LOG("variable");}
	UG_LOG("  | " );

//	Subset informations
	if(verboseLev >= 1)
	{
		for(int si = 0; si < dd->num_subsets(); ++si)
		{
			UG_LOG( " (" << dd->subset_name(si) << ",");
			UG_LOG(std::setw(8) << ConvertNumber(dd->num_indices(si),8,4) << ") ");

		}
	}
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

void IApproximationSpace::
print_parallel_statistic(ConstSmartPtr<DoFDistribution> dd, int verboseLev) const
{
#ifdef UG_PARALLEL
//	Get Process communicator;
	const pcl::ProcessCommunicator& pCom = dd->layouts()->proc_comm();

//	hack since pcl does not support much constness
	DoFDistribution* nonconstDD = const_cast<DoFDistribution*>(dd.get());

//	compute local dof numbers of all masters; this the number of all dofs
//	minus the number of all slave dofs (double counting of slaves can not
//	occur, since each slave is only slave of one master), minus the number of
//	all vertical master dofs (double counting can occure, since a vertical master
//	can be in a horizontal slave interface.
//	Therefore: 	if there are no vert. master dofs, we can simply calculate the
//				number of DoFs by counting. If there are vert. master dofs
//				we need to remove doubles when counting.
	int numMasterDoF;
	if(NumIndices(dd->layouts()->vertical_master()) > 0)
	{
	//	create vector of vertical masters and horiz. slaves, since those are
	//	not reguarded as masters on the dd. All others are masters. (Especially,
	//	also vertical slaves are masters w.r.t. computing like smoothing !!)
		std::vector<IndexLayout::Element> vIndex;
		CollectElements(vIndex, nonconstDD->layouts()->vertical_master());
		CollectElements(vIndex, nonconstDD->layouts()->slave(), false);

	//	sort vector and remove doubles
		std::sort(vIndex.begin(), vIndex.end());
		vIndex.erase(std::unique(vIndex.begin(), vIndex.end()),
		                         vIndex.end());

	//	now remove all horizontal masters, since they are still master though
	//	in the vert master interface
		std::vector<IndexLayout::Element> vIndex2;
		CollectElements(vIndex2, nonconstDD->layouts()->master());
		for(size_t i = 0; i < vIndex2.size(); ++i)
			vIndex.erase(std::remove(vIndex.begin(), vIndex.end(), vIndex2[i]),
			             vIndex.end());

	//	the remaining DoFs are the "slaves"
		numMasterDoF = dd->num_indices() - vIndex.size();
	}
	else
	{
	//	easy case: only subtract masters from slaves
		numMasterDoF = dd->num_indices() - NumIndices(dd->layouts()->slave());
	}

//	global and local values
//	we use uint64 here instead of size_t since the communication via
//	mpi does not support the usage of size_t.
	std::vector<uint64> tNumGlobal, tNumLocal;

//	write number of Masters on this process
	tNumLocal.push_back(numMasterDoF);

//	write local number of dof in a subset for all subsets. For simplicity, we
//	only communicate the total number of dofs (i.e. master+slave)
//	\todo: count slaves in subset and subtract them to get only masters
	for(int si = 0; si < dd->num_subsets(); ++si)
		tNumLocal.push_back(dd->num_indices(si));

//	resize receive array
	tNumGlobal.resize(tNumLocal.size());

//	sum up over processes
	if(!pCom.empty())
	{
		pCom.allreduce(&tNumLocal[0], &tNumGlobal[0], tNumGlobal.size(),
		               PCL_DT_UNSIGNED_LONG_LONG, PCL_RO_SUM);
	}
	else if (pcl::GetNumProcesses() == 1)
	{
		for(size_t i = 0; i < tNumGlobal.size(); ++i)
		{
			tNumGlobal[i] = tNumLocal[i];
		}
	}
	else
	{
		UG_LOG(" Unable to compute informations.");
		return;
	}

//	Total number of DoFs (last arg of 'ConvertNumber()' ("precision") is total
//  width - 4 (two for binary prefix, two for space and decimal point)
	UG_LOG(std::setw(10) << ConvertNumber(tNumGlobal[0],10,6) <<" | ");

//	Overall block size
	const int blockSize = DefaultAlgebra::get().blocksize();
	if(blockSize != AlgebraType::VariableBlockSize) {UG_LOG(std::setw(8)  << blockSize);}
	else {UG_LOG("variable");}
	UG_LOG("  | " );

//	Subset informations
	if(verboseLev>=1)
	{
		for(int si = 0; si < dd->num_subsets(); ++si)
		{
			UG_LOG( " (" << dd->subset_name(si) << ",");
			UG_LOG(std::setw(8) << ConvertNumber(tNumGlobal[si+1],8,4) << ") ");
		}
	}
#endif
}


void IApproximationSpace::print_statistic(int verboseLev) const
{
//	Write info
	UG_LOG("Statistic on DoF-Distribution");
#ifdef UG_PARALLEL
	UG_LOG(" on process " << pcl::GetProcRank() <<
	       " of " << pcl::GetNumProcesses() << " processes");
#endif
	UG_LOG(":\n");

//	check, what to print
	bool bPrintSurface = false, bPrintLevel = false;
	for(size_t i = 0; i < m_vDD.size(); ++i){
		if(m_vDD[i]->grid_level().is_level()) bPrintLevel = true;
		if(m_vDD[i]->grid_level().is_surface()) bPrintSurface = true;
	}

//	Write header line
	if(bPrintLevel || bPrintSurface)
	{
		UG_LOG("         |   Total   | BlockSize | ");
		if(verboseLev >= 1) UG_LOG("(Subset, DoFs) ");
		UG_LOG("\n");
	}

//	Write Infos for Levels
	if(bPrintLevel)
	{
		UG_LOG("  Level  |\n");
		for(size_t i = 0; i < m_vDD.size(); ++i){
			if(!m_vDD[i]->grid_level().is_level()) continue;

			UG_LOG("  " << std::setw(5) << m_vDD[i]->grid_level().level()  << "  |");
			print_statistic(m_vDD[i], verboseLev);
			UG_LOG(std::endl);
		}
	}

//	Write Infos for Surface Grid
	if(bPrintSurface)
	{
		UG_LOG(" Surface |\n");
		for(size_t i = 0; i < m_vDD.size(); ++i){
			if(!m_vDD[i]->grid_level().is_surface()) continue;
			if(m_vDD[i]->grid_level().level() == GridLevel::TOPLEVEL) continue;

			UG_LOG("  " << std::setw(5) << m_vDD[i]->grid_level().level() << "  |");
			print_statistic(m_vDD[i], verboseLev);
			UG_LOG(std::endl);
		}

		for(size_t i = 0; i < m_vDD.size(); ++i){
			if(!m_vDD[i]->grid_level().is_surface()) continue;
			if(m_vDD[i]->grid_level().level() != GridLevel::TOPLEVEL) continue;

			UG_LOG("     top |");
			print_statistic(m_vDD[i], verboseLev);
			UG_LOG(std::endl);
		}
	}

//	if nothing printed
	if(!bPrintLevel && ! bPrintSurface)
	{
		UG_LOG(	"   No DoFs distributed yet (done automatically). \n"
				"   In order to force DoF creation use \n"
				"   ApproximationSpace::init_levels() or "
				"ApproximationSpace::init_surfaces().\n\n");
	}

#ifdef UG_PARALLEL
//	only in case of more than 1 proc
	if(pcl::GetNumProcesses() < 2) return;

//	header
	UG_LOG("Statistic on DoFDistribution on all Processes (m= master, s=slave):\n");

//	Write header line
	if(bPrintLevel || bPrintSurface)
	{
		UG_LOG("         |   Total   | BlockSize | ");
		if(verboseLev >= 1) UG_LOG("(Subset, DoFs) ");
		UG_LOG("\n");
	}

//	Write Infos for Levels
	if(bPrintLevel)
	{
		UG_LOG("  Level  |\n");
		for(size_t i = 0; i < m_vDD.size(); ++i){
			if(!m_vDD[i]->grid_level().is_level()) continue;

			UG_LOG("  " << std::setw(5) << m_vDD[i]->grid_level().level()  << "  |");
			print_parallel_statistic(m_vDD[i], verboseLev);
			UG_LOG(std::endl);
		}
	}

//	Write Infos for Surface Grid
	if(bPrintSurface)
	{
		UG_LOG(" Surface |\n");
		for(size_t i = 0; i < m_vDD.size(); ++i){
			if(!m_vDD[i]->grid_level().is_surface()) continue;
			if(m_vDD[i]->grid_level().level() == GridLevel::TOPLEVEL) continue;

			UG_LOG("  " << std::setw(5) << m_vDD[i]->grid_level().level() << "  |");
			print_parallel_statistic(m_vDD[i], verboseLev);
			UG_LOG(std::endl);
		}

		for(size_t i = 0; i < m_vDD.size(); ++i){
			if(!m_vDD[i]->grid_level().is_surface()) continue;
			if(m_vDD[i]->grid_level().level() != GridLevel::TOPLEVEL) continue;

			UG_LOG("     top |");
			print_parallel_statistic(m_vDD[i], verboseLev);
			UG_LOG(std::endl);
		}
	}
#endif
}

#ifdef UG_PARALLEL
static void PrintLayoutStatistic(ConstSmartPtr<DoFDistribution> dd)
{
	UG_LOG(std::setw(8) << NumIndices(dd->layouts()->master()) <<" | ");
	UG_LOG(std::setw(8) << NumIndices(dd->layouts()->slave()) <<" | ");
	UG_LOG(std::setw(12) << NumIndices(dd->layouts()->vertical_master()) <<" | ");
	UG_LOG(std::setw(12) << NumIndices(dd->layouts()->vertical_slave()));
}
#endif

void IApproximationSpace::print_layout_statistic(int verboseLev) const
{
#ifdef UG_PARALLEL
//	Write info
	UG_LOG("Layouts on Process " <<  GetLogAssistant().get_output_process() << ":\n");

//	Write header line
	UG_LOG(" Level |  Master  |  Slave   | vert. Master | vert. Slave\n");
	UG_LOG("----------------------------------------------------------\n");

//	Write Infos for Levels
	for(size_t i = 0; i < m_vDD.size(); ++i){
		if(!m_vDD[i]->grid_level().is_level()) continue;

		UG_LOG(" " << std::setw(5)<< m_vDD[i]->grid_level().level() << " | ");
		PrintLayoutStatistic(m_vDD[i]);
		UG_LOG(std::endl);
	}

//	Write Infos for Surface Grid
	UG_LOG(" Surf  |\n");
	for(size_t i = 0; i < m_vDD.size(); ++i){
		if(!m_vDD[i]->grid_level().is_surface()) continue;

		UG_LOG(" " << std::setw(5)<< m_vDD[i]->grid_level().level() << " | ");
		PrintLayoutStatistic(m_vDD[i]);
		UG_LOG(std::endl);
	}

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
