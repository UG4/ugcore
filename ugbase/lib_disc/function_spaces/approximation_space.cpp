/*
 * approximation_space.cpp
 *
 *  Created on: 13.12.2011
 *      Author: andreasvogel
 */

#include "approximation_space.h"
#include "lib_disc/domain.h"

namespace ug{

IApproximationSpace::
IApproximationSpace(SmartPtr<subset_handler_type> spMGSH, bool bGroup)
	: m_spMGSH(spMGSH), m_spFunctionPattern(new FunctionPattern(spMGSH)),
	  m_bGrouped(bGroup)
#ifdef UG_PARALLEL
	, m_pDistGridMgr(NULL)
#endif
{}

IApproximationSpace::
IApproximationSpace(SmartPtr<subset_handler_type> spMGSH)
	: m_spMGSH(spMGSH),  m_spFunctionPattern(new FunctionPattern(spMGSH))
#ifdef UG_PARALLEL
	, m_pDistGridMgr(NULL)
#endif
{
//	get blocksize of algebra
	const int blockSize = DefaultAlgebra::get().blocksize();

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
		UG_THROW_FATAL("Cannot determine blocksize of Algebra.");
}

std::vector<ConstSmartPtr<SurfaceDoFDistribution> >
IApproximationSpace::surface_dof_distributions() const
{
	std::vector<ConstSmartPtr<SurfaceDoFDistribution> > vDD;
	if(num_levels() == 0) return vDD;

	vDD.resize(m_vSurfDD.size());
	const_cast<IApproximationSpace*>(this)->surf_dd_required(0,num_levels()-1);
	for(size_t i = 0; i < m_vSurfDD.size(); ++i)
		vDD[i] = m_vSurfDD[i];
	return vDD;
}

SmartPtr<SurfaceDoFDistribution>
IApproximationSpace::surface_dof_distribution(int level)
{
	if(level != GridLevel::TOPLEVEL){
		surf_dd_required(level,level); return m_vSurfDD[level];
	}

	else{
		top_surf_dd_required();
		return m_spTopSurfDD;
	}
}

ConstSmartPtr<SurfaceDoFDistribution>
IApproximationSpace::surface_dof_distribution(int level) const
{
	return const_cast<IApproximationSpace*>(this)->surface_dof_distribution(level);
}

std::vector<ConstSmartPtr<LevelDoFDistribution> >
IApproximationSpace::level_dof_distributions() const
{
	std::vector<ConstSmartPtr<LevelDoFDistribution> > vDD;
	if(num_levels() == 0) return vDD;

	vDD.resize(m_vLevDD.size());
	const_cast<IApproximationSpace*>(this)->level_dd_required(0,num_levels()-1);
	for(size_t i = 0; i < m_vLevDD.size(); ++i)
		vDD[i] = m_vLevDD[i];
	return vDD;
}

SmartPtr<LevelDoFDistribution>
IApproximationSpace::level_dof_distribution(int level)
{
	level_dd_required(level,level); return m_vLevDD[level];
}

ConstSmartPtr<LevelDoFDistribution>
IApproximationSpace::level_dof_distribution(int level) const
{
	const_cast<IApproximationSpace*>(this)->level_dd_required(level,level);
	return m_vLevDD[level];
}

void IApproximationSpace::init_level()
{
	if(num_levels() > 0)
		level_dd_required(0, num_levels()-1);
}

void IApproximationSpace::init_surface()
{
	int numLevGlobal = num_levels();
#ifdef UG_PARALLEL
	int numLevLocal = num_levels();
	pcl::ProcessCommunicator commWorld;
	commWorld.allreduce(&numLevLocal, &numLevGlobal, 1, PCL_DT_INT, PCL_RO_MAX);
#endif

	if(numLevGlobal > 0){
		surf_dd_required(0, numLevGlobal-1);
		top_surf_dd_required();
	}
}

void IApproximationSpace::defragment()
{
//	update surface view
	if(m_spSurfaceView.is_valid())
		m_spSurfaceView->mark_shadows();

//	defragment level dd
	if(num_levels() > m_vLevDD.size())
		level_dd_required(m_vLevDD.size(), num_levels()-1);

	for(size_t lev = 0; lev < m_vLevDD.size(); ++lev)
		if(m_vLevDD[lev].is_valid()) m_vLevDD[lev]->defragment();

//	defragment top dd
	if(m_spTopSurfDD.is_valid()) m_spTopSurfDD->defragment();

//	defragment surface dd
	for(size_t lev = 0; lev < m_vSurfDD.size(); ++lev)
		if(m_vSurfDD[lev].is_valid()) m_vSurfDD[lev]->defragment();
}

template <typename TDD>
void
IApproximationSpace::
print_statistic(ConstSmartPtr<TDD> dd, int verboseLev) const
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

void IApproximationSpace::print_statistic(int verboseLev) const
{
//	Write info
	UG_LOG("DoFDistribution");
#ifdef UG_PARALLEL
	UG_LOG(" on Process " << pcl::GetProcRank());
#endif
	UG_LOG(":\n");

//	Write header line
	UG_LOG("         |   Total   | BlockSize | ");
	if(verboseLev >= 1) UG_LOG("(Subset, DoFs) ");
	UG_LOG("\n");

//	Write Infos for Levels
	UG_LOG("  Level  |\n");
	for(size_t l = 0; l < m_vLevDD.size(); ++l)
	{
		UG_LOG("  " << std::setw(5) << l << "  |");
		print_statistic(level_dof_distribution(l), verboseLev);
		UG_LOG(std::endl);
	}

//	done ?!
	if(m_vSurfDD.empty() && !m_spTopSurfDD.is_valid()) return;

//	Write Infos for Surface Grid
	UG_LOG(" Surface |\n");

	if(!m_vSurfDD.empty())
	{
		for(size_t l = 0; l < m_vSurfDD.size(); ++l)
		{
			UG_LOG("  " << std::setw(5) << l << "  |");
			print_statistic(surface_dof_distribution(l), verboseLev);
			UG_LOG(std::endl);
		}
	}

	if(m_spTopSurfDD.is_valid())
	{
		UG_LOG("     top |");
		print_statistic(surface_dof_distribution(GridLevel::TOPLEVEL), verboseLev);
		UG_LOG(std::endl);
	}

}

void IApproximationSpace::
print_local_dof_statistic(ConstSmartPtr<MGDoFDistribution> dd, int verboseLev) const
{
//	Subset informations
	UG_LOG(dd->num_subsets() << " Subset(s) used (Subset Name, dim): ");
	for(int si = 0; si < dd->num_subsets(); ++si)
	{
		if(si > 0) UG_LOG(", ");
		UG_LOG("(" << dd->subset_name(si) << ", " << dd->dim_subset(si) << ")");
	}
	UG_LOG("\n");

//	Function informations
	UG_LOG(dd->num_fct() << " Function(s) defined (Symbolic Name, dim): ");
	for(size_t fct = 0; fct < dd->num_fct(); ++fct)
	{
		if(fct > 0) UG_LOG(", ");
		UG_LOG("(" << dd->name(fct) << ", " << dd->dim(fct) << ")");
	}
	UG_LOG("\n");

//	print subsets of functions
	if(verboseLev >= 2)
	{
		UG_LOG("Function definition on subsets: \n");
		for(size_t fct = 0; fct < dd->num_fct(); ++fct)
		{
			UG_LOG("   "<<dd->name(fct) << ": ");
			if(dd->is_def_everywhere(fct)) UG_LOG("[everywhere] ");
			bool bFirst = true;
			for(int si = 0; si < dd->num_subsets(); ++si)
			{
				if(bFirst) bFirst = false; else UG_LOG(", ");
				if(!dd->is_def_in_subset(fct, si)) continue;
				UG_LOG(dd->subset_name(si));
			}
		}
		UG_LOG("\n");
	}

//	print info of dofdistribution
	dd->print_local_dof_statistic(verboseLev);
}

void IApproximationSpace::print_local_dof_statistic(int verboseLev) const
{
	UG_LOG("\nLocal DoF Statistic for DoFDistribution:");
	UG_LOG("\n-----------------------------------------------\n");

	if(m_spLevMGDD.is_valid())
	{
		print_local_dof_statistic(m_spLevMGDD, verboseLev);
		return;
	}
	if(!m_vSurfDD.empty())
	{
		print_local_dof_statistic(m_vSurfDD[0], verboseLev);
		return;
	}
	if(m_spTopSurfDD.is_valid())
	{
		print_local_dof_statistic(m_spTopSurfDD, verboseLev);
		return;
	}

	UG_LOG(" not avaible .\n");
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

template <typename TDD>
static void PrintLayoutStatistic(ConstSmartPtr<TDD> dd)
{
	UG_LOG(std::setw(8) << NumIndices(dd->master_layout()) <<" | ");
	UG_LOG(std::setw(8) << NumIndices(dd->slave_layout()) <<" | ");
	UG_LOG(std::setw(12) << NumIndices(dd->vertical_master_layout()) <<" | ");
	UG_LOG(std::setw(12) << NumIndices(dd->vertical_slave_layout()));
}
#endif

void IApproximationSpace::print_layout_statistic(int verboseLev) const
{
#ifdef UG_PARALLEL
//	Write info
	UG_LOG("Layouts on Process " <<  pcl::GetOutputProcRank() << ":\n");

//	Write header line
	UG_LOG(" Level |  Master  |  Slave   | vert. Master | vert. Slave\n");
	UG_LOG("----------------------------------------------------------\n");

//	Write Infos for Levels
	for(size_t l = 0; l < m_vLevDD.size(); ++l)
	{
		UG_LOG(" " << std::setw(5)<< l << " | ");
		PrintLayoutStatistic(level_dof_distribution(l));
		UG_LOG(std::endl);
	}

//	Write Infos for Surface Grid
	if(!m_vSurfDD.empty())
	{
		UG_LOG(" Surf  |\n");
		for(size_t l = 0; l < m_vLevDD.size(); ++l)
		{
			UG_LOG(" " << std::setw(5)<< l << " | ");
			PrintLayoutStatistic(surface_dof_distribution(l));
			UG_LOG(std::endl);
		}
	}

#else
	UG_LOG(" No Layouts in sequential code.\n");
#endif
}

void IApproximationSpace::level_dd_required(size_t fromLevel, size_t toLevel)
{
//	check correct arguments
	if(fromLevel > toLevel)
		UG_THROW_FATAL("fromLevel must be small than toLevel");

//	check level
	if(toLevel >= this->num_levels())
		UG_THROW_FATAL("Required Level "<<toLevel<<", but only "<<
					   this->num_levels()<<" in the MultiGrid.");

//	if not yet MGLevelDD allocated
	if(!m_spLevMGDD.is_valid())
		m_spLevMGDD = SmartPtr<LevelMGDoFDistribution>
					(new LevelMGDoFDistribution(this->m_spMGSH, *m_spFunctionPattern
#ifdef UG_PARALLEL
					                              ,m_pDistGridMgr
#endif
					));

//	resize level
	if(m_vLevDD.size() < toLevel+1) m_vLevDD.resize(toLevel+1, NULL);

//	allocate Level DD if needed
	for(size_t lvl = fromLevel; lvl <= toLevel; ++lvl)
	{
		if(!m_vLevDD[lvl].is_valid()){
			m_vLevDD[lvl] = SmartPtr<LevelDoFDistribution>
							(new LevelDoFDistribution(m_spLevMGDD, m_spMGSH, lvl));
		}
	}
}

void IApproximationSpace::surf_dd_required(size_t fromLevel, size_t toLevel)
{
//	check correct arguments
	if(fromLevel > toLevel)
		UG_THROW_FATAL("fromLevel must be small than toLevel");

//	resize level
	if(m_vSurfDD.size() < toLevel+1) m_vSurfDD.resize(toLevel+1, NULL);

	surface_level_view_required(fromLevel, toLevel);

//	allocate Level DD if needed
	for(size_t lvl = fromLevel; lvl <= toLevel; ++lvl)
	{
		if(!m_vSurfDD[lvl].is_valid())
			m_vSurfDD[lvl] = SmartPtr<SurfaceDoFDistribution>
							(new SurfaceDoFDistribution(
									this->m_spMGSH, *m_spFunctionPattern, m_vSurfLevView[lvl], lvl
#ifdef UG_PARALLEL
					                              ,m_pDistGridMgr
#endif
									));
	}
}

void IApproximationSpace::top_surf_dd_required()
{
	top_surface_level_view_required();

//	allocate Level DD if needed
	if(!m_spTopSurfDD.is_valid())
		m_spTopSurfDD = SmartPtr<SurfaceDoFDistribution>
							(new SurfaceDoFDistribution(
									this->m_spMGSH, *m_spFunctionPattern, m_spTopSurfLevView, GridLevel::TOPLEVEL
#ifdef UG_PARALLEL
					                              ,m_pDistGridMgr
#endif
									));
}

void IApproximationSpace::top_surface_level_view_required()
{
//	allocate surface view if needed
	if(!m_spSurfaceView.is_valid())
		m_spSurfaceView = SmartPtr<SurfaceView>
						 (new SurfaceView(m_spMGSH
#ifdef UG_PARALLEL
					                              ,m_pDistGridMgr
#endif
						                       ));

//	allocate Level DD if needed
	if(!m_spTopSurfLevView.is_valid())
		m_spTopSurfLevView = SmartPtr<SurfaceLevelView>
							(new SurfaceLevelView(m_spSurfaceView, GridLevel::TOPLEVEL));
}

void IApproximationSpace::surface_level_view_required(size_t fromLevel, size_t toLevel)
{
//	check correct arguments
	if(fromLevel > toLevel)
		UG_THROW_FATAL("fromLevel must be small than toLevel");

//	allocate surface view if needed
	if(!m_spSurfaceView.is_valid())
		m_spSurfaceView = SmartPtr<SurfaceView>
						 (new SurfaceView(m_spMGSH
#ifdef UG_PARALLEL
					                              ,m_pDistGridMgr
#endif
						                       ));

//	resize level
	if(m_vSurfLevView.size() < toLevel+1) m_vSurfLevView.resize(toLevel+1, NULL);

//	allocate Level DD if needed
	for(size_t lvl = fromLevel; lvl <= toLevel; ++lvl)
	{
		if(!m_vSurfLevView[lvl].is_valid())
			m_vSurfLevView[lvl] = SmartPtr<SurfaceLevelView>
							(new SurfaceLevelView(m_spSurfaceView, lvl));
	}
}


template <typename TDomain>
ApproximationSpace<TDomain>::
ApproximationSpace(SmartPtr<domain_type> domain)
	: IApproximationSpace(domain->subset_handler()),
	  m_spDomain(domain)
{
	if(!m_spDomain.is_valid())
		UG_THROW_FATAL("Domain, passed to ApproximationSpace, is invalid.");
	if(!m_spMGSH.is_valid())
		UG_THROW_FATAL("SubsetHandler, passed to ApproximationSpace, is invalid.");

#ifdef UG_PARALLEL
	this->set_dist_grid_mgr(domain->distributed_grid_manager());
#endif
};

} // end namespace ug

template class ug::ApproximationSpace<ug::Domain1d>;
template class ug::ApproximationSpace<ug::Domain2d>;
template class ug::ApproximationSpace<ug::Domain3d>;
