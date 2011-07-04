/*
 * dirichet_post_process_interface_impl.h
 *
 *  Created on: 04.08.2010
 *      Author: andreasvogel
 */


#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__POST_PROCESS__DIRICHLET_BOUNDARY__DIRICHLET_POST_PROCESS_INTERFACE_IMPL__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__POST_PROCESS__DIRICHLET_BOUNDARY__DIRICHLET_POST_PROCESS_INTERFACE_IMPL__

#include "./post_process_interface.h"

namespace ug{


////////////////////////////////////////////////
// Set geometric object
////////////////////////////////////////////////

template<typename TAlgebra>
bool
IDirichletPostProcess<TAlgebra>::
set_geometric_object_type(int id, IDirichletPostProcessNeed need)
{
	if(function_registered(id, need))
	{
		m_id = id;
		return true;
	}
	else
	{
		m_id = -1;
		UG_LOG("No or not all functions registered "
				"for object with reference object id " << id << ".\n");
	}
	return false;
};

template<typename TAlgebra>
bool
IDirichletPostProcess<TAlgebra>::
function_registered(int id, IDirichletPostProcessNeed need)
{
	// if nothing required, return true
	if(need == IDPPN_NONE) return true;

	// loop functions must exist in any case
	if(!prepare_element_loop_function_registered(id))
	{
		UG_LOG("prepare_element_loop() function not registered for id " << id <<".\n");
		return false;
	}
	if(!prepare_element_function_registered(id))
	{
		UG_LOG("prepare_element(...) function not registered for id " << id <<".\n");
		return false;
	}
	if(!finish_element_loop_function_registered(id))
	{
		UG_LOG("finish_element_loop() function not registered for id " << id <<".\n");
		return false;
	}

	// check if functions for defect exist
	if(need == IDPPN_DEFECT)
	{
		if(!post_process_d_function_registered(id))
		{
			UG_LOG("post_process_d(...) function not registered for id " << id <<".\n");
			return false;
		}

		if(!post_process_f_function_registered(id))
		{
			UG_LOG("post_process_f(...) function not registered for id " << id <<".\n");
			return false;
		}
	}

	// check if functions for jacobian exist
	if(need == IDPPN_JACOBIAN)
	{
		if(!post_process_J_function_registered(id))
		{
			UG_LOG("post_process_J(...) function not registered for id " << id <<".\n");
			return false;
		}
	}

	// check if functions for jacobian exist
	if(need == IDPPN_LINEAR)
	{
		if(!post_process_J_function_registered(id))
		{
			UG_LOG("post_process_J(...) function not registered for id " << id <<".\n");
			return false;
		}
		if(!post_process_f_function_registered(id))
		{
			UG_LOG("post_process_f(...) function not registered for id " << id <<".\n");
			return false;
		}
	}

	// check if functions for jacobian exist
	if(need == IDPPN_SOLUTION)
	{
		if(!set_solution_function_registered(id))
		{
			UG_LOG("set_solution(...) function not registered for id " << id <<".\n");
			return false;
		}
	}

	return true;
};

////////////////////////////////////////////////
// Prepare Element Loop
////////////////////////////////////////////////

template<typename TAlgebra>
template<typename TAssFunc>
void
IDirichletPostProcess<TAlgebra>::
reg_prepare_vol_loop_fct(int id, TAssFunc func)
{
//	make sure that there is enough space
	if((size_t)id >= m_vPrepareElementLoopFunc.size())
		m_vPrepareElementLoopFunc.resize(id+1, 0);

	m_vPrepareElementLoopFunc[id] = (PrepareElementLoopFunc)func;
};

template<typename TAlgebra>
bool
IDirichletPostProcess<TAlgebra>::
prepare_element_loop_function_registered(int id)
{
	if(id >= 0 && (size_t)id < m_vPrepareElementLoopFunc.size())
	{
		if(m_vPrepareElementLoopFunc[id] != 0)
			return true;
	}
	return false;
}


////////////////////////////////////////////////
// Prepare Element
////////////////////////////////////////////////

template<typename TAlgebra>
template<typename TAssFunc>
void
IDirichletPostProcess<TAlgebra>::
reg_prepare_vol_fct(int id, TAssFunc func)
{
//	make sure that there is enough space
	if((size_t)id >= m_vPrepareElementFunc.size())
		m_vPrepareElementFunc.resize(id+1, 0);

	m_vPrepareElementFunc[id] = (PrepareElementFunc)func;
};

template<typename TAlgebra>
bool
IDirichletPostProcess<TAlgebra>::
prepare_element_function_registered(int id)
{
	if(id >= 0 && (size_t)id < m_vPrepareElementFunc.size())
	{
		if(m_vPrepareElementFunc[id] != 0)
			return true;
	}
	return false;
}


////////////////////////////////////////////////
// Finish Element Loop
////////////////////////////////////////////////

template<typename TAlgebra>
template<typename TAssFunc>
void
IDirichletPostProcess<TAlgebra>::
reg_finish_vol_loop_fct(int id, TAssFunc func)
{
//	make sure that there is enough space
	if((size_t)id >= m_vFinishElementLoopFunc.size())
		m_vFinishElementLoopFunc.resize(id+1, 0);

	m_vFinishElementLoopFunc[id] = (FinishElementLoopFunc)func;
};

template<typename TAlgebra>
bool
IDirichletPostProcess<TAlgebra>::
finish_element_loop_function_registered(int id)
{
	if(id >= 0 && (size_t)id < m_vFinishElementLoopFunc.size())
	{
		if(m_vFinishElementLoopFunc[id] != 0)
			return true;
	}
	return false;
}


////////////////////////////////////////////////
// Assemble Jacobian
////////////////////////////////////////////////

template<typename TAlgebra>
template<typename TAssFunc>
void
IDirichletPostProcess<TAlgebra>::
register_post_process_J_function(int id, TAssFunc func)
{
//	make sure that there is enough space
	if((size_t)id >= m_vPostProcessJFunc.size())
		m_vPostProcessJFunc.resize(id+1, 0);

	m_vPostProcessJFunc[id] = (PostProcessJFunc)func;
};

template<typename TAlgebra>
bool
IDirichletPostProcess<TAlgebra>::
post_process_J_function_registered(int id)
{
	if(id >= 0 && (size_t)id < m_vPostProcessJFunc.size())
	{
		if(m_vPostProcessJFunc[id] != 0)
			return true;
	}
	return false;
}


////////////////////////////////////////////////
// Assemble Defect
////////////////////////////////////////////////

template<typename TAlgebra>
template<typename TAssFunc>
void
IDirichletPostProcess<TAlgebra>::
register_post_process_d_function(int id, TAssFunc func)
{
//	make sure that there is enough space
	if((size_t)id >= m_vPostProcessDFunc.size())
		m_vPostProcessDFunc.resize(id+1, 0);

	m_vPostProcessDFunc[id] = (PostProcessDFunc)func;
};

template<typename TAlgebra>
bool
IDirichletPostProcess<TAlgebra>::
post_process_d_function_registered(int id)
{
	if(id >= 0 && (size_t)id < m_vPostProcessDFunc.size())
	{
		if(m_vPostProcessDFunc[id] != 0)
			return true;
	}
	return false;
}


////////////////////////////////////////////////
// Assemble rhs
////////////////////////////////////////////////

template<typename TAlgebra>
template<typename TAssFunc>
void
IDirichletPostProcess<TAlgebra>::
register_post_process_f_function(int id, TAssFunc func)
{
//	make sure that there is enough space
	if((size_t)id >= m_vPostProcessFFunc.size())
		m_vPostProcessFFunc.resize(id+1, 0);

	m_vPostProcessFFunc[id] = (PostProcessFFunc)func;
};

template<typename TAlgebra>
bool
IDirichletPostProcess<TAlgebra>::
post_process_f_function_registered(int id)
{
	if(id >= 0 && (size_t)id < m_vPostProcessFFunc.size())
	{
		if(m_vPostProcessFFunc[id] != 0)
			return true;
	}
	return false;
}


////////////////////////////////////////////////
// set solution
////////////////////////////////////////////////

template<typename TAlgebra>
template<typename TAssFunc>
void
IDirichletPostProcess<TAlgebra>::
register_set_solution_function(int id, TAssFunc func)
{
//	make sure that there is enough space
	if((size_t)id >= m_vSetSolutionFunc.size())
		m_vSetSolutionFunc.resize(id+1, 0);

	m_vSetSolutionFunc[id] = (SetSolutionFunc)func;
};

template<typename TAlgebra>
bool
IDirichletPostProcess<TAlgebra>::
set_solution_function_registered(int id)
{
	if(id >= 0 && (size_t)id < m_vSetSolutionFunc.size())
	{
		if(m_vSetSolutionFunc[id] != 0)
			return true;
	}
	return false;
}


}

#endif /* __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__POST_PROCESS__DIRICHLET_BOUNDARY__DIRICHLET_POST_PROCESS_INTERFACE_IMPL__ */
