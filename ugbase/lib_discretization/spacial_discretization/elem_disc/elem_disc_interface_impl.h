/*
 * elem_disc_interface_impl.h
 *
 *  Created on: 07.07.2010
 *      Author: andreasvogel
 */


#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__ELEM_DISC_INTERFACE_IMPL__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__ELEM_DISC_INTERFACE_IMPL__

#include "elem_disc_interface.h"

namespace ug{


////////////////////////////////////////////////
// Set geometric object
////////////////////////////////////////////////

template<typename TAlgebra>
bool
IElemDisc<TAlgebra>::
set_geometric_object_type(int id, IElemDiscNeed need)
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
IElemDisc<TAlgebra>::
function_registered(int id, IElemDiscNeed need)
{
	// if nothing required, return true
	if(need == IEDN_NONE) return true;

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
	if(need == IEDN_DEFECT)
	{
		if(!assemble_A_function_registered(id))
		{
			UG_LOG("assemble_A(...) function not registered for id " << id <<".\n");
			return false;
		}

		if(!assemble_M_function_registered(id))
		{
			UG_LOG("assemble_M(...) function not registered for id " << id <<".\n");
			return false;
		}
		if(!assemble_f_function_registered(id))
		{
			UG_LOG("assemble_f(...) function not registered for id " << id <<".\n");
			return false;
		}
	}

	// check if functions for jacobian exist
	if(need == IEDN_JACOBIAN)
	{
		if(!assemble_JA_function_registered(id))
		{
			UG_LOG("assemble_JA(...) function not registered for id " << id <<".\n");
			return false;
		}

		if(!assemble_JM_function_registered(id))
		{
			UG_LOG("assemble_JM(...) function not registered for id " << id <<".\n");
			return false;
		}
	}

	// check if functions for linear exist
	if(need == IEDN_LINEAR)
	{
		if(!assemble_JA_function_registered(id))
		{
			UG_LOG("assemble_JA(...) function not registered for id " << id <<".\n");
			return false;
		}

		if(!assemble_JM_function_registered(id))
		{
			UG_LOG("assemble_JM(...) function not registered for id " << id <<".\n");
			return false;
		}
		if(!assemble_f_function_registered(id))
		{
			UG_LOG("assemble_f(...) function not registered for id " << id <<".\n");
			return false;
		}
	}

	// check if functions for stiffness matrix exist
	if(need == IEDN_STIFFNESS)
	{
		if(!assemble_JA_function_registered(id))
		{
			UG_LOG("assemble_JA(...) function not registered for id " << id <<".\n");
			return false;
		}
	}

	// check if functions for stiffness matrix exist
	if(need == IEDN_MASS)
	{
		if(!assemble_JM_function_registered(id))
		{
			UG_LOG("assemble_JM(...) function not registered for id " << id <<".\n");
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
IElemDisc<TAlgebra>::
register_prepare_element_loop_function(int id, TAssFunc func)
{
//	make sure that there is enough space
	if((size_t)id >= m_vPrepareElementLoopFunc.size())
		m_vPrepareElementLoopFunc.resize(id+1, 0);

	m_vPrepareElementLoopFunc[id] = (PrepareElementLoopFunc)func;
};

template<typename TAlgebra>
bool
IElemDisc<TAlgebra>::
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
IElemDisc<TAlgebra>::
register_prepare_element_function(int id, TAssFunc func)
{
//	make sure that there is enough space
	if((size_t)id >= m_vPrepareElementFunc.size())
		m_vPrepareElementFunc.resize(id+1, 0);

	m_vPrepareElementFunc[id] = (PrepareElementFunc)func;
};

template<typename TAlgebra>
bool
IElemDisc<TAlgebra>::
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
IElemDisc<TAlgebra>::
register_finish_element_loop_function(int id, TAssFunc func)
{
//	make sure that there is enough space
	if((size_t)id >= m_vFinishElementLoopFunc.size())
		m_vFinishElementLoopFunc.resize(id+1, 0);

	m_vFinishElementLoopFunc[id] = (FinishElementLoopFunc)func;
};

template<typename TAlgebra>
bool
IElemDisc<TAlgebra>::
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
IElemDisc<TAlgebra>::
register_assemble_JA_function(int id, TAssFunc func)
{
//	make sure that there is enough space
	if((size_t)id >= m_vAssembleJAFunc.size())
		m_vAssembleJAFunc.resize(id+1, 0);

	m_vAssembleJAFunc[id] = (AssembleJAFunc)func;
};

template<typename TAlgebra>
bool
IElemDisc<TAlgebra>::
assemble_JA_function_registered(int id)
{
	if(id >= 0 && (size_t)id < m_vAssembleJAFunc.size())
	{
		if(m_vAssembleJAFunc[id] != 0)
			return true;
	}
	return false;
}


template<typename TAlgebra>
template<typename TAssFunc>
void
IElemDisc<TAlgebra>::
register_assemble_JM_function(int id, TAssFunc func)
{
//	make sure that there is enough space
	if((size_t)id >= m_vAssembleJMFunc.size())
		m_vAssembleJMFunc.resize(id+1, 0);

	m_vAssembleJMFunc[id] = (AssembleJMFunc)func;
};

template<typename TAlgebra>
bool
IElemDisc<TAlgebra>::
assemble_JM_function_registered(int id)
{
	if(id >= 0 && (size_t)id < m_vAssembleJMFunc.size())
	{
		if(m_vAssembleJMFunc[id] != 0)
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
IElemDisc<TAlgebra>::
register_assemble_A_function(int id, TAssFunc func)
{
//	make sure that there is enough space
	if((size_t)id >= m_vAssembleAFunc.size())
		m_vAssembleAFunc.resize(id+1, 0);

	m_vAssembleAFunc[id] = (AssembleAFunc)func;
};

template<typename TAlgebra>
bool
IElemDisc<TAlgebra>::
assemble_A_function_registered(int id)
{
	if(id >= 0 && (size_t)id < m_vAssembleAFunc.size())
	{
		if(m_vAssembleAFunc[id] != 0)
			return true;
	}
	return false;
}


template<typename TAlgebra>
template<typename TAssFunc>
void
IElemDisc<TAlgebra>::
register_assemble_M_function(int id, TAssFunc func)
{
//	make sure that there is enough space
	if((size_t)id >= m_vAssembleMFunc.size())
		m_vAssembleMFunc.resize(id+1, 0);

	m_vAssembleMFunc[id] = (AssembleMFunc)func;
};

template<typename TAlgebra>
bool
IElemDisc<TAlgebra>::
assemble_M_function_registered(int id)
{
	if(id >= 0 && (size_t)id < m_vAssembleMFunc.size())
	{
		if(m_vAssembleMFunc[id] != 0)
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
IElemDisc<TAlgebra>::
register_assemble_f_function(int id, TAssFunc func)
{
//	make sure that there is enough space
	if((size_t)id >= m_vAssembleFFunc.size())
		m_vAssembleFFunc.resize(id+1, 0);

	m_vAssembleFFunc[id] = (AssembleFFunc)func;
};

template<typename TAlgebra>
bool
IElemDisc<TAlgebra>::
assemble_f_function_registered(int id)
{
	if(id >= 0 && (size_t)id < m_vAssembleFFunc.size())
	{
		if(m_vAssembleFFunc[id] != 0)
			return true;
	}
	return false;
}


}

#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__ELEM_DISC_INTERFACE_IMPL__ */
