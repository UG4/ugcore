/*
 * approximation_space_impl.h
 *
 *  Created on: 19.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__FUNCTION_SPACE__APPROXIMATION_SPACE_IMPL__
#define __H__UG__LIB_DISC__FUNCTION_SPACE__APPROXIMATION_SPACE_IMPL__

#include "approximation_space.h"
#include "../../common/common.h"

namespace ug{

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
void
ApproximationSpace<TDomain, TDoFDistribution, TAlgebra>::init(bool bInitDoFs)
{
//	check if already init
	if(m_bInit) return;

//	lock function pattern
	this->lock();

//	set subsethandler to DofManager
	if(!m_MGDoFManager.assign_multi_grid_subset_handler(*(this->m_pMGSH)))
		UG_THROW_FATAL("In 'ApproximationSpace::init':"
						" Cannot assign multi grid subset handler.");

//	set the function pattern for dofmanager
	if(!m_MGDoFManager.assign_function_pattern(*this))
		UG_THROW_FATAL("In 'ApproximationSpace::init':"
						" Cannot assign Function Pattern.");

#ifdef UG_PARALLEL
//	set distributed grid manager
	m_MGDoFManager.set_distributed_grid_manager(
			*this->m_pDomain->get_distributed_grid_manager());
#endif

	if(bInitDoFs)
		if(!m_MGDoFManager.enable_indices())
			UG_THROW_FATAL("In 'ApproximationSpace::init':"
							" Cannot distribute dofs.");

//	remember init flag
	m_bInit = true;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
typename ApproximationSpace<TDomain, TDoFDistribution, TAlgebra>::function_type*
ApproximationSpace<TDomain, TDoFDistribution, TAlgebra>::create_level_function(size_t level)
{
//	init space
	init();

//	enable level dofs
	if(!m_bLevelDoFInit){
		if(!m_MGDoFManager.enable_level_indices()){
			UG_THROW_FATAL( "ApproximationSpace: Cannot distribute level dofs.");
		}
		else{
			m_bLevelDoFInit = true;
		}
	}

//	get level dof distribution
	dof_distribution_type* dofDistr = m_MGDoFManager.get_level_dof_distribution(level);

//	check distribution
	if(dofDistr == NULL)
		UG_THROW_FATAL( "ApproximationSpace: No level DoFDistribution created.");

//	create new function
	return new function_type(*this, *dofDistr);
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
typename ApproximationSpace<TDomain, TDoFDistribution, TAlgebra>::function_type*
ApproximationSpace<TDomain, TDoFDistribution, TAlgebra>::create_surface_function()
{
//	init space
	init();

//	enable surface dofs
	if(!m_bSurfDoFInit){
		if(!m_MGDoFManager.enable_surface_indices()){
			UG_THROW_FATAL( "ApproximationSpace: Cannot distribute surface dofs.");
		}
		else{
			m_bSurfDoFInit = true;
		}
	}

//	get surface dof distribution
	dof_distribution_type* dofDistr = m_MGDoFManager.get_surface_dof_distribution();

//	check distribution
	if(dofDistr == NULL)
		UG_THROW_FATAL( "ApproximationSpace: No surface DoFDistribution created.");

//	create new function
	return new function_type(*this, *dofDistr);
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
typename ApproximationSpace<TDomain, TDoFDistribution, TAlgebra>::dof_distribution_type&
ApproximationSpace<TDomain, TDoFDistribution, TAlgebra>::get_surface_dof_distribution()
{
//	init space
	init();

//	enable surface dofs
	if(!m_bSurfDoFInit){
		if(!m_MGDoFManager.enable_surface_indices()){
			UG_THROW_FATAL( "ApproximationSpace: Cannot distribute surface dofs.");
		}
		else{
			m_bSurfDoFInit = true;
		}
	}

	dof_distribution_type* dofDistr = m_MGDoFManager.get_surface_dof_distribution();

	if(dofDistr == NULL)
		UG_THROW_FATAL( "ApproximationSpace: No surface DoFDistribution created.");

	return *dofDistr;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
const typename ApproximationSpace<TDomain, TDoFDistribution, TAlgebra>::dof_distribution_type&
ApproximationSpace<TDomain, TDoFDistribution, TAlgebra>::get_surface_dof_distribution() const
{
	ApproximationSpace<TDomain, TDoFDistribution, TAlgebra>* This =
			const_cast<ApproximationSpace<TDomain, TDoFDistribution, TAlgebra>*>(this);

	return This->get_surface_dof_distribution();
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
typename ApproximationSpace<TDomain, TDoFDistribution, TAlgebra>::dof_distribution_type&
ApproximationSpace<TDomain, TDoFDistribution, TAlgebra>::get_level_dof_distribution(size_t level)
{
//	init space
	init();

//	enable surface dofs
	if(!m_bLevelDoFInit){
		if(!m_MGDoFManager.enable_level_indices()){
			UG_THROW_FATAL( "ApproximationSpace: Cannot distribute level dofs.");
		}
		else{
			m_bLevelDoFInit = true;
		}
	}

	return *(m_MGDoFManager.get_level_dof_distribution(level));
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
const typename ApproximationSpace<TDomain, TDoFDistribution, TAlgebra>::dof_distribution_type&
ApproximationSpace<TDomain, TDoFDistribution, TAlgebra>::get_level_dof_distribution(size_t level) const
{
	ApproximationSpace<TDomain, TDoFDistribution, TAlgebra>* This =
			const_cast<ApproximationSpace<TDomain, TDoFDistribution, TAlgebra>*>(this);

	return This->get_level_dof_distribution(level);
}

}


#endif /* __H__UG__LIB_DISC__FUNCTION_SPACE__APPROXIMATION_SPACE_IMPL__ */
