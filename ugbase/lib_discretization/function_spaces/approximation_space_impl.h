/*
 * approximation_space_impl.h
 *
 *  Created on: 19.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__FUNCTION_SPACE__APPROXIMATION_SPACE_IMPL__
#define __H__LIBDISCRETIZATION__FUNCTION_SPACE__APPROXIMATION_SPACE_IMPL__

#include "approximation_space.h"

namespace ug{

struct UG_ERROR_DoFDistributionMissing{};

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool
ApproximationSpace<TDomain, TDoFDistribution, TAlgebra>::init()
{
//	Check, that domain is given
	if(this->m_pMGSH == NULL)
	{
		UG_LOG("ERROR in 'ApproximationSpace::init':"
				" No Domain assigned to Approximation Space.\n");
		return false;
	}

//	check, if already initialized
	if(m_bInit)
	{
		UG_LOG("WARNING in 'ApproximationSpace::init':"
				" Approximation Space already initialized. You cannot alter"
				" the pattern. This call was useless.\n");
		return true;
	}

//	lock function pattern
	this->lock();

//	set subsethandler to DofManager
	if(!m_MGDoFManager.assign_multi_grid_subset_handler(*(this->m_pMGSH)))
	{
		UG_LOG("In 'ApproximationSpace::init':"
				" Cannot assign multi grid subset handler.\n");
		return false;
	}

//	set the function pattern for dofmanager
	if(!m_MGDoFManager.assign_function_pattern(*this))
	{
		UG_LOG("In 'ApproximationSpace::init':"
				" Cannot assign Function Pattern.\n");
		return false;
	}

#ifdef UG_PARALLEL
//	set distributed grid manager
	m_MGDoFManager.set_distributed_grid_manager(
			*this->m_pDomain->get_distributed_grid_manager());
#endif

//	enable all dofs
	if(!m_MGDoFManager.enable_dofs())
	{
		UG_LOG("In 'ApproximationSpace::init':"
				" Cannot distribute dofs.\n");
		return false;
	}

//	remember init flag
	m_bInit = true;

//	we're done
	return true;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
typename ApproximationSpace<TDomain, TDoFDistribution, TAlgebra>::function_type*
ApproximationSpace<TDomain, TDoFDistribution, TAlgebra>::create_level_function(size_t level)
{
	if(!m_bInit)
	{
		UG_LOG("Approximation Space not initialized.\n");
		return NULL;
	}

	dof_distribution_type* dofDistr = m_MGDoFManager.get_level_dof_distribution(level);
	if(dofDistr == NULL)
	{
		throw(UG_ERROR_DoFDistributionMissing());
	}

	function_type* gridFct = new function_type(*this, *dofDistr);
	return gridFct;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
typename ApproximationSpace<TDomain, TDoFDistribution, TAlgebra>::function_type*
ApproximationSpace<TDomain, TDoFDistribution, TAlgebra>::create_surface_function()
{
	if(!m_bInit)
	{
		UG_LOG("Approximation Space not initialized.\n");
		return NULL;
	}

	dof_distribution_type* dofDistr = m_MGDoFManager.get_surface_dof_distribution();
	if(dofDistr == NULL)
	{
		throw(UG_ERROR_DoFDistributionMissing());
	}

	function_type* gridFct = new function_type(*this, *dofDistr);

	return gridFct;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
typename ApproximationSpace<TDomain, TDoFDistribution, TAlgebra>::dof_distribution_type&
ApproximationSpace<TDomain, TDoFDistribution, TAlgebra>::get_surface_dof_distribution()
{
	if(!m_bInit)
		throw(UG_ERROR_DoFDistributionMissing());

	dof_distribution_type* dofDistr = m_MGDoFManager.get_surface_dof_distribution();

	if(dofDistr == NULL)
		throw(UG_ERROR_DoFDistributionMissing());

	return *dofDistr;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
const typename ApproximationSpace<TDomain, TDoFDistribution, TAlgebra>::dof_distribution_type&
ApproximationSpace<TDomain, TDoFDistribution, TAlgebra>::get_surface_dof_distribution() const
{
	if(!m_bInit)
		throw(UG_ERROR_DoFDistributionMissing());

	const dof_distribution_type* dofDistr = m_MGDoFManager.get_surface_dof_distribution();

	if(dofDistr == NULL)
		throw(UG_ERROR_DoFDistributionMissing());

	return *dofDistr;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
typename ApproximationSpace<TDomain, TDoFDistribution, TAlgebra>::dof_distribution_type&
ApproximationSpace<TDomain, TDoFDistribution, TAlgebra>::get_level_dof_distribution(size_t level)
{
	if(!m_bInit)
		throw(UG_ERROR_DoFDistributionMissing());

	return *(m_MGDoFManager.get_level_dof_distribution(level));
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
const typename ApproximationSpace<TDomain, TDoFDistribution, TAlgebra>::dof_distribution_type&
ApproximationSpace<TDomain, TDoFDistribution, TAlgebra>::get_level_dof_distribution(size_t level) const
{
	if(!m_bInit)
		throw(UG_ERROR_DoFDistributionMissing());

	return *(m_MGDoFManager.get_level_dof_distribution(level));
}

}


#endif /* __H__LIBDISCRETIZATION__FUNCTION_SPACE__APPROXIMATION_SPACE_IMPL__ */
