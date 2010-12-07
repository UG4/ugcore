/*
 * elem_disc_interface_impl.h
 *
 *  Created on: 07.07.2010
 *      Author: andreasvogel
 */


#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__ELEM_DISC_INTERFACE_IMPL__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__ELEM_DISC_INTERFACE_IMPL__

#include "elem_disc_interface.h"
#include "lib_discretization/common/groups_util.h"

namespace ug{

/////////////////////////////////////////////////////
// Functions and Subsets (default implementation)
/////////////////////////////////////////////////////

template<typename TAlgebra>
bool
IElemDisc<TAlgebra>::
set_functions(const char* functions)
{
//	check that function pattern exists
	if(m_pPattern == NULL)
	{
		UG_LOG("IElemDisc::set_functions: Function Pattern not set.\n");
		return false;
	}

//	create Function Group
	FunctionGroup functionGroup;

//	convert string
	if(!ConvertStringToFunctionGroup(functionGroup, *m_pPattern, functions))
	{
		UG_LOG("ERROR while parsing Functions.\n");
		return false;
	}

//	forward request
	return set_functions(functionGroup);
}

template<typename TAlgebra>
bool
IElemDisc<TAlgebra>::
set_functions(const FunctionGroup& funcGroup)
{
//	check that pattern is set
	if(m_pPattern == NULL)
	{
		UG_LOG("IElemDisc::set_functions: No underlying Function Pattern set.\n");
		return false;
	}

//	check, that pattern is the same
	if(m_pPattern != funcGroup.get_function_pattern())
	{
		UG_LOG("IElemDisc::set_functions: Given Function group does "
				"not have correct underlying Function Pattern.\n");
		return false;
	}

//	check number of functions
	if(funcGroup.num_fct() != this->num_fct())
	{
		UG_LOG("IElemDisc::set_functions: Exactly " << this->num_fct() <<
				" functions needed, but given "<< funcGroup.num_fct() <<" functions.\n");
		return false;
	}

//	check trial spaces
	for(size_t i = 0; i < this->num_fct(); ++i)
	{
		if(funcGroup.local_shape_function_set_id(i) !=
				this->local_shape_function_set_id(i))
		{
			UG_LOG("IElemDisc::set_functions: Function " << i << " has incorrect"
					" local shape function type.\n");
			return false;
		}
	}

//	get Dimension of subset
	int dim = funcGroup.get_function_dimension();

	if(dim == -1)
	{
		UG_LOG("IElemDisc::set_functions: Given function group does not have a unique dimension.\n");
		return false;
	}

	if(!m_SubsetGroup.empty())
	{
	//	check Dimension of function group
		int dim_subsets = m_SubsetGroup.get_subset_dimension();
		if(dim != dim_subsets)
		{
			UG_LOG("IElemDisc::set_functions: Dimension of already set subset "
					"group has dimension " << dim_subsets << ", but passed "
					"function lives in " << dim << "d. Cannot choose function.\n");
			return false;
		}

	// 	check that function is defined for segment
		for(size_t i = 0; i <funcGroup.num_fct(); ++i)
		{
			const size_t fct = funcGroup[i];
			for(size_t si = 0; si < m_SubsetGroup.num_subsets(); ++si)
			{
				const int subsetIndex = m_SubsetGroup[si];
				if(!m_pPattern->is_def_in_subset(fct, subsetIndex))
				{
					UG_LOG("IElemDisc::set_subsets: Function "
							<< fct << " not defined in subset " << subsetIndex << ".\n");
					return false;
				}
			}
		}
	}

//	remember group (copy)
	m_FunctionGroup = funcGroup;
	return true;
}

template<typename TAlgebra>
bool
IElemDisc<TAlgebra>::
set_subsets(const char* subsets)
{
//	check that function pattern exists
	if(m_pPattern == NULL)
	{
		UG_LOG("IElemDisc::set_subsets: Function Pattern not set.\n");
		return false;
	}

//	create Function Group
	SubsetGroup subsetGroup;

//	convert string
	if(!ConvertStringToSubsetGroup(subsetGroup, *m_pPattern, subsets))
	{
		UG_LOG("ERROR while parsing Subsets.\n");
		return false;
	}

//	forward request
	return set_subsets(subsetGroup);
}

template<typename TAlgebra>
bool
IElemDisc<TAlgebra>::
set_subsets(const SubsetGroup& subsetGroup)
{
//	check that pattern is set
	if(m_pPattern == NULL)
	{
		UG_LOG("IElemDisc::set_subsets: No underlying Function Pattern set.\n");
		return false;
	}

//	check, that subset handler is the same
	if(m_pPattern->get_subset_handler() != subsetGroup.get_subset_handler())
	{
		UG_LOG("IElemDisc::set_subsets: Given subset group does "
				"not have correct underlying subset handler.\n");
		return false;
	}

//	get Dimension of subset
	int dim = subsetGroup.get_subset_dimension();

	if(dim == -1)
	{
		UG_LOG("IElemDisc::set_subsets: Given subset group does not have a unique dimension.\n");
		return false;
	}

	if(!m_FunctionGroup.empty())
	{
		//	check Dimension of function group
		int dim_functions = m_FunctionGroup.get_function_dimension();
		if(dim != dim_functions)
		{
			UG_LOG("IElemDisc::set_subsets: Dimension of already set function group"
					" does not have same dimension as passed subsets.\n");
			return false;
		}

	// 	check that function is defined for segment
		for(size_t i = 0; i < m_FunctionGroup.num_fct(); ++i)
		{
			const size_t fct = m_FunctionGroup[i];
			for(size_t si = 0; si < subsetGroup.num_subsets(); ++si)
			{
				const int subsetIndex = subsetGroup[si];
				if(!m_pPattern->is_def_in_subset(fct, subsetIndex))
				{
					UG_LOG("IElemDisc::set_subsets: Function "
							<< fct << " not defined in subset " << subsetIndex << ".\n");
					return false;
				}
			}
		}
	}

//	remember group (copy)
	m_SubsetGroup = subsetGroup;
	return true;
}


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
