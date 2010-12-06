/*
 * domain_discretization_impl.h
 *
 *  Created on: 29.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DOMAIN_DISCRETIZATION__IMPL__
#define __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DOMAIN_DISCRETIZATION__IMPL__

#include "domain_discretization_impl.h"
#include "lib_discretization/common/groups_util.h"

namespace ug{

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//  Time Independent
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

//////////////////////////
// assemble Jacobian
//////////////////////////
template <typename TDoFDistribution, typename TAlgebra>
IAssembleReturn
DomainDiscretization<TDoFDistribution, TAlgebra>::
assemble_jacobian(matrix_type& J, const vector_type& u,
                  const dof_distribution_type& dofDistr)
{
//	assemble normal element discretizations
	for(size_t i = 0; i < m_vElemDisc.size(); ++i)
	{
	//	get subset group
		const SubsetGroup& subsetGroup = m_vElemDisc[i]->get_subset_group();

	//	loop subsets
		for(size_t j = 0; j < subsetGroup.num_subsets(); ++j)
		{
		//	get subset index and dimesion
			const int subsetIndex = subsetGroup[j];
			const int dim = subsetGroup.get_subset_dimension(j);

		//	assemble
			if(!AssembleJacobian<IElemDisc<TAlgebra>, TDoFDistribution, TAlgebra>
			(*m_vElemDisc[i], J, u, dofDistr, subsetIndex, dim))
				return IAssemble_ERROR;
		}
	}

//	assemble coupled element discretizations
	for(size_t i = 0; i < m_vCoupledDisc.size(); ++i)
	{
	//	loop subsets
		for(size_t j = 0; j < m_vCoupledDisc[i].subsetGroup.num_subsets(); ++j)
		{
		//	get subset index and dimesion
			const int subsetIndex = m_vCoupledDisc[i].subsetGroup[j];
			const int dim = m_vCoupledDisc[i].subsetGroup.get_subset_dimension(j);

			//	assemble
			if(!AssembleJacobian<CoupledSystem<TAlgebra>, TDoFDistribution, TAlgebra>
				(*m_vCoupledDisc[i].disc, J, u, dofDistr, subsetIndex, dim))
					return IAssemble_ERROR;
		}
	}

//	post process
	for(size_t type = 0; type < PPT_NUM_POST_PROCESS_TYPES; ++type)
	{
		for(size_t i = 0; i < m_vvPostProcess[type].size(); ++i)
		{
			if(m_vvPostProcess[type][i]->post_process_jacobian(J, u, dofDistr)
					!= IAssemble_OK)
				return IAssemble_ERROR;
		}
	}

//	done
	return IAssemble_OK;
}


//////////////////////////
// assemble Defect
//////////////////////////
template <typename TDoFDistribution, typename TAlgebra>
IAssembleReturn
DomainDiscretization<TDoFDistribution, TAlgebra>::
assemble_defect(vector_type& d, const vector_type& u,
                const dof_distribution_type& dofDistr)
{
//	assemble normal element discretizations
	for(size_t i = 0; i < m_vElemDisc.size(); ++i)
	{
	//	get subset group
		const SubsetGroup& subsetGroup = m_vElemDisc[i]->get_subset_group();

	//	loop subsets
		for(size_t j = 0; j < subsetGroup.num_subsets(); ++j)
		{
		//	get subset index and dimesion
			const int subsetIndex = subsetGroup[j];
			const int dim = subsetGroup.get_subset_dimension(j);

			if(!AssembleDefect<IElemDisc<TAlgebra>, TDoFDistribution, TAlgebra>
				(*m_vElemDisc[i], d, u, dofDistr, subsetIndex, dim))
					return IAssemble_ERROR;
		}
	}

//	assemble coupled element discretizations
	for(size_t i = 0; i < m_vCoupledDisc.size(); ++i)
	{
	//	loop subsets
		for(size_t j = 0; j < m_vCoupledDisc[i].subsetGroup.num_subsets(); ++j)
		{
		//	get subset index and dimesion
			const int subsetIndex = m_vCoupledDisc[i].subsetGroup[j];
			const int dim = m_vCoupledDisc[i].subsetGroup.get_subset_dimension(j);

			//	assemble
			if(!AssembleDefect<CoupledSystem<TAlgebra>, TDoFDistribution, TAlgebra>
				(*m_vCoupledDisc[i].disc, d, u, dofDistr, subsetIndex, dim))
					return IAssemble_ERROR;
		}
	}

//	post process
	for(size_t type = 0; type < PPT_NUM_POST_PROCESS_TYPES; ++type)
	{
		for(size_t i = 0; i < m_vvPostProcess[type].size(); ++i)
		{
			if(m_vvPostProcess[type][i]->post_process_defect(d, u, dofDistr)
					!= IAssemble_OK)
				return IAssemble_ERROR;
		}
	}

//	done
	return IAssemble_OK;
}

//////////////////////////
// assemble Matrix and RHS
//////////////////////////
template <typename TDoFDistribution, typename TAlgebra>
IAssembleReturn
DomainDiscretization<TDoFDistribution, TAlgebra>::
assemble_linear(matrix_type& mat, vector_type& rhs, const vector_type& u,
                const dof_distribution_type& dofDistr)
{
//	assemble normal element discretizations
	for(size_t i = 0; i < m_vElemDisc.size(); ++i)
	{
	//	get subset group
		const SubsetGroup& subsetGroup = m_vElemDisc[i]->get_subset_group();

	//	loop subsets
		for(size_t j = 0; j < subsetGroup.num_subsets(); ++j)
		{
		//	get subset index and dimesion
			const int subsetIndex = subsetGroup[j];
			const int dim = subsetGroup.get_subset_dimension(j);

			if(!AssembleLinear<IElemDisc<TAlgebra>, TDoFDistribution, TAlgebra>
			(*m_vElemDisc[i], mat, rhs, u, dofDistr, subsetIndex, dim))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
						" Cannot assemble linear for Elem Disc.\n");
				return IAssemble_ERROR;
			}
		}
	}

//	assemble coupled element discretizations
	for(size_t i = 0; i < m_vCoupledDisc.size(); ++i)
	{
	//	loop subsets
		for(size_t j = 0; j < m_vCoupledDisc[i].subsetGroup.num_subsets(); ++j)
		{
		//	get subset index and dimesion
			const int subsetIndex = m_vCoupledDisc[i].subsetGroup[j];
			const int dim = m_vCoupledDisc[i].subsetGroup.get_subset_dimension(j);

			//	assemble
			if(!AssembleLinear<CoupledSystem<TAlgebra>, TDoFDistribution, TAlgebra>
				(*m_vCoupledDisc[i].disc, mat, rhs, u, dofDistr, subsetIndex, dim))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
						" Cannot assemble linear for Cpl Elem Disc.\n");
				return IAssemble_ERROR;
			}
		}
	}

//	post process
	for(size_t type = 0; type < PPT_NUM_POST_PROCESS_TYPES; ++type)
	{
		for(size_t i = 0; i < m_vvPostProcess[type].size(); ++i)
		{
			if(m_vvPostProcess[type][i]->post_process_linear(mat, rhs, u, dofDistr)
					!= IAssemble_OK)
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
						" Cannot post process.\n");
				return IAssemble_ERROR;
			}
		}
	}

//	done
	return IAssemble_OK;
}

//////////////////////////////////
// set solution in dirichlet dofs
//////////////////////////////////
template <typename TDoFDistribution, typename TAlgebra>
IAssembleReturn
DomainDiscretization<TDoFDistribution, TAlgebra>::
assemble_solution(vector_type& u, const dof_distribution_type& dofDistr)
{
//	post process dirichlet
	for(size_t i = 0; i < m_vvPostProcess[PPT_DIRICHLET].size(); ++i)
	{
		if(m_vvPostProcess[PPT_DIRICHLET][i]->post_process_solution(u, dofDistr)
				!= IAssemble_OK)
			return IAssemble_ERROR;
	}

//	done
	return IAssemble_OK;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//  Time Independent
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////
// Assemble Jacobian
//////////////////////////////////
template <typename TDoFDistribution, typename TAlgebra>
IAssembleReturn
DomainDiscretization<TDoFDistribution, TAlgebra>::
assemble_jacobian(matrix_type& J, const vector_type& u,
                  const dof_distribution_type& dofDistr,
                  number time, number s_m, number s_a)
{
//	assemble normal element discretizations
	for(size_t i = 0; i < m_vElemDisc.size(); ++i)
	{
	//	get subset group
		const SubsetGroup& subsetGroup = m_vElemDisc[i]->get_subset_group();

	//	loop subsets
		for(size_t j = 0; j < subsetGroup.num_subsets(); ++j)
		{
		//	get subset index and dimesion
			const int subsetIndex = subsetGroup[j];
			const int dim = subsetGroup.get_subset_dimension(j);

			if(!AssembleJacobian<IElemDisc<TAlgebra>, TDoFDistribution, TAlgebra>
				(*m_vElemDisc[i], J, u, dofDistr, subsetIndex, dim, time, s_m, s_a))
					return IAssemble_ERROR;
		}
	}

//	assemble coupled element discretizations
	for(size_t i = 0; i < m_vCoupledDisc.size(); ++i)
	{
	//	loop subsets
		for(size_t j = 0; j < m_vCoupledDisc[i].subsetGroup.num_subsets(); ++j)
		{
		//	get subset index and dimesion
			const int subsetIndex = m_vCoupledDisc[i].subsetGroup[j];
			const int dim = m_vCoupledDisc[i].subsetGroup.get_subset_dimension(j);

			//	assemble
			if(!AssembleJacobian<CoupledSystem<TAlgebra>, TDoFDistribution, TAlgebra>
				(*m_vCoupledDisc[i].disc, J, u, dofDistr, subsetIndex, dim, time, s_m, s_a))
					return IAssemble_ERROR;
		}
	}

//	post process
	for(size_t type = 0; type < PPT_NUM_POST_PROCESS_TYPES; ++type)
	{
		for(size_t i = 0; i < m_vvPostProcess[type].size(); ++i)
		{
			if(m_vvPostProcess[type][i]->post_process_jacobian(J, u, dofDistr, time)
					!= IAssemble_OK)
				return IAssemble_ERROR;
		}
	}

//	done
	return IAssemble_OK;
}

//////////////////////////////////
// Assemble Defect
//////////////////////////////////
template <typename TDoFDistribution, typename TAlgebra>
IAssembleReturn
DomainDiscretization<TDoFDistribution, TAlgebra>::
assemble_defect(vector_type& d, const vector_type& u,
                const dof_distribution_type& dofDistr,
                number time, number s_m, number s_a)
{
//	assemble normal element discretizations
	for(size_t i = 0; i < m_vElemDisc.size(); ++i)
	{
	//	get subset group
		const SubsetGroup& subsetGroup = m_vElemDisc[i]->get_subset_group();

	//	loop subsets
		for(size_t j = 0; j < subsetGroup.num_subsets(); ++j)
		{
		//	get subset index and dimesion
			const int subsetIndex = subsetGroup[j];
			const int dim = subsetGroup.get_subset_dimension(j);

			if(!AssembleDefect<IElemDisc<TAlgebra>, TDoFDistribution, TAlgebra>
				(*m_vElemDisc[i], d, u, dofDistr, subsetIndex, dim, time, s_m, s_a))
					return IAssemble_ERROR;
		}
	}

//	assemble coupled element discretizations
	for(size_t i = 0; i < m_vCoupledDisc.size(); ++i)
	{
	//	loop subsets
		for(size_t j = 0; j < m_vCoupledDisc[i].subsetGroup.num_subsets(); ++j)
		{
		//	get subset index and dimesion
			const int subsetIndex = m_vCoupledDisc[i].subsetGroup[j];
			const int dim = m_vCoupledDisc[i].subsetGroup.get_subset_dimension(j);

			//	assemble
			if(!AssembleDefect<CoupledSystem<TAlgebra>, TDoFDistribution, TAlgebra>
				(*m_vCoupledDisc[i].disc, d, u, dofDistr, subsetIndex, dim, time, s_m, s_a))
					return IAssemble_ERROR;
		}
	}

//	post process
	for(size_t type = 0; type < PPT_NUM_POST_PROCESS_TYPES; ++type)
	{
		for(size_t i = 0; i < m_vvPostProcess[type].size(); ++i)
		{
			if(m_vvPostProcess[type][i]->post_process_defect(d, u, dofDistr, time)
					!= IAssemble_OK)
				return IAssemble_ERROR;
		}
	}

//	done
	return IAssemble_OK;
}

//////////////////////////////////
// Assemble Matrix and RHS
//////////////////////////////////
template <typename TDoFDistribution, typename TAlgebra>
IAssembleReturn
DomainDiscretization<TDoFDistribution, TAlgebra>::
assemble_linear(matrix_type& mat, vector_type& rhs, const vector_type& u,
                const dof_distribution_type& dofDistr,
                number time, number s_m, number s_a)
{
//	assemble normal element discretizations
	for(size_t i = 0; i < m_vElemDisc.size(); ++i)
	{
	//	get subset group
		const SubsetGroup& subsetGroup = m_vElemDisc[i]->get_subset_group();

	//	loop subsets
		for(size_t j = 0; j < subsetGroup.num_subsets(); ++j)
		{
		//	get subset index and dimesion
			const int subsetIndex = subsetGroup[j];
			const int dim = subsetGroup.get_subset_dimension(j);

			if(!AssembleLinear<IElemDisc<TAlgebra>, TDoFDistribution, TAlgebra>
				(*m_vElemDisc[i], mat, rhs, u, dofDistr, subsetIndex, dim, time, s_m, s_a))
					return IAssemble_ERROR;
		}
	}

//	assemble coupled element discretizations
	for(size_t i = 0; i < m_vCoupledDisc.size(); ++i)
	{
	//	loop subsets
		for(size_t j = 0; j < m_vCoupledDisc[i].subsetGroup.num_subsets(); ++j)
		{
		//	get subset index and dimesion
			const int subsetIndex = m_vCoupledDisc[i].subsetGroup[j];
			const int dim = m_vCoupledDisc[i].subsetGroup.get_subset_dimension(j);

			//	assemble
			if(!AssembleLinear<CoupledSystem<TAlgebra>, TDoFDistribution, TAlgebra>
				(*m_vCoupledDisc[i].disc, mat, rhs, u, dofDistr, subsetIndex, dim, time, s_m, s_a))
					return IAssemble_ERROR;
		}
	}

//	post process
	for(size_t type = 0; type < PPT_NUM_POST_PROCESS_TYPES; ++type)
	{
		for(size_t i = 0; i < m_vvPostProcess[type].size(); ++i)
		{
			if(m_vvPostProcess[type][i]->post_process_linear(mat, rhs, u, dofDistr, time)
					!= IAssemble_OK)
				return IAssemble_ERROR;
		}
	}

//	done
	return IAssemble_OK;
}

//////////////////////////////////
// set dirichlet values
//////////////////////////////////
template <typename TDoFDistribution, typename TAlgebra>
IAssembleReturn
DomainDiscretization<TDoFDistribution, TAlgebra>::
assemble_solution(vector_type& u, const dof_distribution_type& dofDistr, number time)
{
//	post process
	for(size_t i = 0; i < m_vvPostProcess[PPT_DIRICHLET].size(); ++i)
	{
		if(m_vvPostProcess[PPT_DIRICHLET][i]->post_process_solution(u, dofDistr, time)
				!= IAssemble_OK)
			return IAssemble_ERROR;
	}

//	done
	return IAssemble_OK;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//  Infrastructure
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////
// addcoupled system
//////////////////////////////////
template <typename TDoFDistribution, typename TAlgebra>
bool
DomainDiscretization<TDoFDistribution, TAlgebra>::
add(CoupledSystem<TAlgebra>& coupledSystem,
    const FunctionGroup& functionGroup,
    const SubsetGroup& subsetGroup)
{
// 	check if number of functions match
	if(coupledSystem.num_fct() != functionGroup.num_fct())
	{
		UG_LOG("Wrong number of local functions for Elemet Discretization.\n");
		UG_LOG("Needed: " << coupledSystem.num_fct() << ", given: "
		       	  << functionGroup.num_fct() << ".\n");
		UG_LOG("Cannot add element discretization.\n");
		return false;
	}

//	check that type of function match
	for(size_t i = 0; i < coupledSystem.num_fct(); ++i)
	{
		if(coupledSystem.local_shape_function_set_id(functionGroup[i]) !=
				functionGroup.local_shape_function_set_id(i))
		{
			UG_LOG("Function " << functionGroup.get_function_name(i)
			       << " has wrong Shape Function.\n");
			UG_LOG("Solution does not match requirements of discretization.\n");
			return false;
		}
	}

//	get Dimension of subset
	int dim = subsetGroup.get_subset_dimension();

	if(dim == -1)
	{
		UG_LOG("Given subset group does not have a unique dimension.\n");
		return false;
	}

//	get Dimension of function group
	int dim_functions = functionGroup.get_function_dimension();

	if(dim_functions == -1)
	{
		UG_LOG("Given function group does not have a unique dimension.\n");
		return false;
	}

//	check that dimension of subsets and dimension of function group is equal
	if(dim != dim_functions)
	{
		UG_LOG("Dimension of Function Group and dimension "
				"of subset group does not match.\n");
		return false;
	}

//	remember (copy)
	m_vCoupledDisc.push_back(CoupledDisc(coupledSystem, functionGroup, subsetGroup));
	return true;
}

//////////////////////////////////
// add coupled system (by string)
//////////////////////////////////
template <typename TDoFDistribution, typename TAlgebra>
bool
DomainDiscretization<TDoFDistribution, TAlgebra>::
add(CoupledSystem<TAlgebra>& coupledSystem, const FunctionPattern& pattern,
    const char* functions, const char* subsets)
{
//	create Function Group and Subset Group
	FunctionGroup functionGroup;
	SubsetGroup subsetGroup;

//	convert strings
	if(!ConvertStringToSubsetGroup(subsetGroup, pattern, subsets))
	{
		UG_LOG("ERROR while parsing Subsets.\n");
		return false;
	}
	if(!ConvertStringToFunctionGroup(functionGroup, pattern, functions))
	{
		UG_LOG("ERROR while parsing Functions.\n");
		return false;
	}

//	forward request
	return add(coupledSystem, functionGroup, subsetGroup);
}

}

#endif /*__H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DOMAIN_DISCRETIZATION__IMPL__*/
