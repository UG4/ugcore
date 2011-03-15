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

//////////////////////////
// assemble Mass Matrix
//////////////////////////
template <typename TDoFDistribution, typename TAlgebra>
IAssembleReturn
DomainDiscretization<TDoFDistribution, TAlgebra>::
assemble_mass_matrix(matrix_type& M, const vector_type& u,
                     const dof_distribution_type& dofDistr)
{
//	Union of Subsets
	SubsetGroup unionSubsets;

//	create list of all subsets
	if(!CreateUnionOfSubsets(unionSubsets, m_vElemDisc))
	{
		UG_LOG("ERROR: Can not create union of subsets.\n");
		return IAssemble_ERROR;
	}

//	loop subsets
	for(size_t i = 0; i < unionSubsets.num_subsets(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = unionSubsets.dim(i);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(m_bForceRegGrid) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc<TAlgebra>*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, si);

	//	assemble on suitable elements
		switch(dim)
		{
		case 1:
			if(!AssembleMassMatrix<Edge>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									     M, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_mass_matrix':"
						"Cannot assemble Edges.\n");
				return IAssemble_ERROR;
			}
			break;
		case 2:
			if(!AssembleMassMatrix<Triangle>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									     M, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_mass_matrix':"
						"Cannot assemble Triangles.\n");
				return IAssemble_ERROR;
			}
			if(!AssembleMassMatrix<Quadrilateral>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									     M, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_mass_matrix':"
						"Cannot assemble Quadrilaterals.\n");
				return IAssemble_ERROR;
			}
			break;
		case 3:
			if(!AssembleMassMatrix<Tetrahedron>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									     M, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_mass_matrix':"
						"Cannot assemble Tetrahedrons.\n");
				return IAssemble_ERROR;
			}
			if(!AssembleMassMatrix<Pyramid>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									     M, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_mass_matrix':"
						"Cannot assemble Pyramids.\n");
				return IAssemble_ERROR;
			}
			if(!AssembleMassMatrix<Prism>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									     M, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_mass_matrix':"
						"Cannot assemble Prisms.\n");
				return IAssemble_ERROR;
			}
			if(!AssembleMassMatrix<Hexahedron>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									     M, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_mass_matrix':"
						"Cannot assemble Hexahedrons.\n");
				return IAssemble_ERROR;
			}
			break;
		default:UG_LOG("ERROR in 'DomainDiscretization::assemble_mass_matrix':"
				"Dimension " << dim << " not supported.\n");
				return IAssemble_ERROR;
		}
	}

//	post process
	for(size_t type = 0; type < PPT_NUM_POST_PROCESS_TYPES; ++type)
	{
		for(size_t i = 0; i < m_vvPostProcess[type].size(); ++i)
		{
			if(m_vvPostProcess[type][i]->post_process_jacobian(M, u, dofDistr)
					!= IAssemble_OK)
				return IAssemble_ERROR;
		}
	}

//	done
	return IAssemble_OK;
}

template <typename TDoFDistribution, typename TAlgebra>
IAssembleReturn
DomainDiscretization<TDoFDistribution, TAlgebra>::
assemble_stiffness_matrix(matrix_type& A, const vector_type& u,
                          const dof_distribution_type& dofDistr)
{
//	Union of Subsets
	SubsetGroup unionSubsets;

//	create list of all subsets
	if(!CreateUnionOfSubsets(unionSubsets, m_vElemDisc))
	{
		UG_LOG("ERROR: Can not create union of subsets.\n");
		return IAssemble_ERROR;
	}

//	loop subsets
	for(size_t i = 0; i < unionSubsets.num_subsets(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = unionSubsets.dim(i);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(m_bForceRegGrid) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc<TAlgebra>*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, si);

	//	assemble on suitable elements
		switch(dim)
		{
		case 1:
			if(!AssembleStiffnessMatrix<Edge>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
										         A, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_stiffness_matrix':"
						"Cannot assemble Edges.\n");
				return IAssemble_ERROR;
			}
			break;
		case 2:
			if(!AssembleStiffnessMatrix<Triangle>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
										         A, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_stiffness_matrix':"
						"Cannot assemble Triangles.\n");
				return IAssemble_ERROR;
			}
			if(!AssembleStiffnessMatrix<Quadrilateral>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
										         A, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_stiffness_matrix':"
						"Cannot assemble Quadrilaterals.\n");
				return IAssemble_ERROR;
			}
			break;
		case 3:
			if(!AssembleStiffnessMatrix<Tetrahedron>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
										         A, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_stiffness_matrix':"
						"Cannot assemble Tetrahedrons.\n");
				return IAssemble_ERROR;
			}
			if(!AssembleStiffnessMatrix<Pyramid>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
										         A, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_stiffness_matrix':"
						"Cannot assemble Pyramids.\n");
				return IAssemble_ERROR;
			}
			if(!AssembleStiffnessMatrix<Prism>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
										         A, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_stiffness_matrix':"
						"Cannot assemble Prisms.\n");
				return IAssemble_ERROR;
			}
			if(!AssembleStiffnessMatrix<Hexahedron>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
										         A, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_stiffness_matrix':"
						"Cannot assemble Hexahedrons.\n");
				return IAssemble_ERROR;
			}
			break;
		default:UG_LOG("ERROR in 'DomainDiscretization::assemble_stiffness_matrix':"
				"Dimension " << dim << " not supported.\n");
				return IAssemble_ERROR;
		}
	}

//	post process
	for(size_t type = 0; type < PPT_NUM_POST_PROCESS_TYPES; ++type)
	{
		for(size_t i = 0; i < m_vvPostProcess[type].size(); ++i)
		{
			if(m_vvPostProcess[type][i]->post_process_jacobian(A, u, dofDistr)
					!= IAssemble_OK)
				return IAssemble_ERROR;
		}
	}

//	done
	return IAssemble_OK;
}

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
assemble_jacobian(matrix_type& J,
                  const vector_type& u,
                  const dof_distribution_type& dofDistr)
{
//	Union of Subsets
	SubsetGroup unionSubsets;

//	create list of all subsets
	if(!CreateUnionOfSubsets(unionSubsets, m_vElemDisc))
	{
		UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
				" Can not create union of subsets.\n");
		return IAssemble_ERROR;
	}

//	loop subsets
	for(size_t i = 0; i < unionSubsets.num_subsets(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = unionSubsets.dim(i);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(m_bForceRegGrid) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc<TAlgebra>*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, si);

	//	assemble on suitable elements
		switch(dim)
		{
		case 1:
			if(!AssembleJacobian<Edge>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
			                           J, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
						"Cannot assemble Edges.\n");
				return IAssemble_ERROR;
			}
			break;
		case 2:
			if(!AssembleJacobian<Triangle>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
			                           J, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
						"Cannot assemble Triangles.\n");
				return IAssemble_ERROR;
			}
			if(!AssembleJacobian<Quadrilateral>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
			                           J, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
						"Cannot assemble Quadrilaterals.\n");
				return IAssemble_ERROR;
			}
			break;
		case 3:
			if(!AssembleJacobian<Tetrahedron>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
			                           J, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
						"Cannot assemble Tetrahedrons.\n");
				return IAssemble_ERROR;
			}
			if(!AssembleJacobian<Pyramid>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
			                           J, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
						"Cannot assemble Pyramids.\n");
				return IAssemble_ERROR;
			}
			if(!AssembleJacobian<Prism>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
			                           J, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
						"Cannot assemble Prisms.\n");
				return IAssemble_ERROR;
			}
			if(!AssembleJacobian<Hexahedron>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
			                           J, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
						"Cannot assemble Hexahedrons.\n");
				return IAssemble_ERROR;
			}
			break;
		default:UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
				"Dimension " << dim << " not supported.\n");
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
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
						" Cannot execute post process " << i << ".\n");
				return IAssemble_ERROR;
			}
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
assemble_defect(vector_type& d,
                const vector_type& u,
                const dof_distribution_type& dofDistr)
{
//	Union of Subsets
	SubsetGroup unionSubsets;

//	create list of all subsets
	if(!CreateUnionOfSubsets(unionSubsets, m_vElemDisc))
	{
		UG_LOG("ERROR: Can not create union of subsets.\n");
		return IAssemble_ERROR;
	}

//	loop subsets
	for(size_t i = 0; i < unionSubsets.num_subsets(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = unionSubsets.dim(i);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(m_bForceRegGrid) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc<TAlgebra>*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, si);

	//	assemble on suitable elements
		switch(dim)
		{
		case 1:
			if(!AssembleDefect<Edge>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   d, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_defect':"
						"Cannot assemble Edges.\n");
				return IAssemble_ERROR;
			}
			break;
		case 2:
			if(!AssembleDefect<Triangle>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   d, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_defect':"
						"Cannot assemble Triangles.\n");
				return IAssemble_ERROR;
			}
			if(!AssembleDefect<Quadrilateral>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   d, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_defect':"
						"Cannot assemble Quadrilaterals.\n");
				return IAssemble_ERROR;
			}
			break;
		case 3:
			if(!AssembleDefect<Tetrahedron>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   d, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_defect':"
						"Cannot assemble Tetrahedrons.\n");
				return IAssemble_ERROR;
			}
			if(!AssembleDefect<Pyramid>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   d, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_defect':"
						"Cannot assemble Pyramids.\n");
				return IAssemble_ERROR;
			}
			if(!AssembleDefect<Prism>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   d, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_defect':"
						"Cannot assemble Prisms.\n");
				return IAssemble_ERROR;
			}
			if(!AssembleDefect<Hexahedron>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   d, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_defect':"
						"Cannot assemble Hexahedrons.\n");
				return IAssemble_ERROR;
			}
			break;
		default:UG_LOG("ERROR in 'DomainDiscretization::assemble_defect':"
				"Dimension " << dim << " not supported.\n");
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
assemble_linear(matrix_type& mat, vector_type& rhs,
                const vector_type& u,
                const dof_distribution_type& dofDistr)
{
//	Union of Subsets
	SubsetGroup unionSubsets;

//	create list of all subsets
	if(!CreateUnionOfSubsets(unionSubsets, m_vElemDisc))
	{
		UG_LOG("ERROR: Can not create union of subsets.\n");
		return IAssemble_ERROR;
	}

//	loop subsets
	for(size_t i = 0; i < unionSubsets.num_subsets(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = unionSubsets.dim(i);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(m_bForceRegGrid) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc<TAlgebra>*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, si);

	//	assemble on suitable elements
		switch(dim)
		{
		case 1:
			if(!AssembleLinear<Edge>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   mat, rhs, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
						"Cannot assemble Edges.\n");
				return IAssemble_ERROR;
			}
			break;
		case 2:
			if(!AssembleLinear<Triangle>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   mat, rhs, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
						"Cannot assemble Triangles.\n");
				return IAssemble_ERROR;
			}
			if(!AssembleLinear<Quadrilateral>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   mat, rhs, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
						"Cannot assemble Quadrilaterals.\n");
				return IAssemble_ERROR;
			}
			break;
		case 3:
			if(!AssembleLinear<Tetrahedron>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   mat, rhs, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
						"Cannot assemble Tetrahedrons.\n");
				return IAssemble_ERROR;
			}
			if(!AssembleLinear<Pyramid>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   mat, rhs, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
						"Cannot assemble Pyramids.\n");
				return IAssemble_ERROR;
			}
			if(!AssembleLinear<Prism>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   mat, rhs, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
						"Cannot assemble Prisms.\n");
				return IAssemble_ERROR;
			}
			if(!AssembleLinear<Hexahedron>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   mat, rhs, u))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
						"Cannot assemble Hexahedrons.\n");
				return IAssemble_ERROR;
			}
			break;
		default:UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
				"Dimension " << dim << " not supported.\n");
				return IAssemble_ERROR;
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
assemble_jacobian(matrix_type& J,
                  const vector_type& u, number time,
                  const PreviousSolutions<vector_type>& prevSol,
                  const dof_distribution_type& dofDistr,
                  number s_m0, number s_a0)
{
//	check that s_m is 1.0
	if(s_m0 != 1.0)
	{
		UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
				" It is assumed to have s_m == 1\n");
	}

//	Union of Subsets
	SubsetGroup unionSubsets;

//	create list of all subsets
	if(!CreateUnionOfSubsets(unionSubsets, m_vElemDisc))
	{
		UG_LOG("ERROR: Can not create union of subsets.\n");
		return IAssemble_ERROR;
	}

//	loop subsets
	for(size_t i = 0; i < unionSubsets.num_subsets(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = unionSubsets.dim(i);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(m_bForceRegGrid) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc<TAlgebra>*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, si);

	//	assemble on suitable elements
		switch(dim)
		{
		case 1:
			if(!AssembleJacobian<Edge>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   J, u, s_a0, time))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
						"Cannot assemble Edges.\n");
				return IAssemble_ERROR;
			}
			break;
		case 2:
			if(!AssembleJacobian<Triangle>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   J, u, s_a0, time))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
						"Cannot assemble Triangles.\n");
				return IAssemble_ERROR;
			}
			if(!AssembleJacobian<Quadrilateral>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   J, u, s_a0, time))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
						"Cannot assemble Quadrilaterals.\n");
				return IAssemble_ERROR;
			}
			break;
		case 3:
			if(!AssembleJacobian<Tetrahedron>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   J, u, s_a0, time))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
						"Cannot assemble Tetrahedrons.\n");
				return IAssemble_ERROR;
			}
			if(!AssembleJacobian<Pyramid>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   J, u, s_a0, time))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
						"Cannot assemble Pyramids.\n");
				return IAssemble_ERROR;
			}
			if(!AssembleJacobian<Prism>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   J, u, s_a0, time))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
						"Cannot assemble Prisms.\n");
				return IAssemble_ERROR;
			}
			if(!AssembleJacobian<Hexahedron>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   J, u, s_a0, time))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
						"Cannot assemble Hexahedrons.\n");
				return IAssemble_ERROR;
			}
			break;
		default:UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
				"Dimension " << dim << " not supported.\n");
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
assemble_defect(vector_type& d,
                const vector_type& u, number time,
                const PreviousSolutions<vector_type>& prevSol,
                const dof_distribution_type& dofDistr,
                number s_m, number s_a)
{
//	Union of Subsets
	SubsetGroup unionSubsets;

//	create list of all subsets
	if(!CreateUnionOfSubsets(unionSubsets, m_vElemDisc))
	{
		UG_LOG("ERROR: Can not create union of subsets.\n");
		return IAssemble_ERROR;
	}

//	loop subsets
	for(size_t i = 0; i < unionSubsets.num_subsets(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = unionSubsets.dim(i);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(m_bForceRegGrid) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc<TAlgebra>*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, si);

	//	assemble on suitable elements
		switch(dim)
		{
		case 1:
			if(!AssembleDefect<Edge>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   d, u, s_m, s_a, time))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_defect':"
						"Cannot assemble Edges.\n");
				return IAssemble_ERROR;
			}
			break;
		case 2:
			if(!AssembleDefect<Triangle>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   d, u, s_m, s_a, time))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_defect':"
						"Cannot assemble Triangles.\n");
				return IAssemble_ERROR;
			}
			if(!AssembleDefect<Quadrilateral>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   d, u, s_m, s_a, time))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_defect':"
						"Cannot assemble Quadrilaterals.\n");
				return IAssemble_ERROR;
			}
			break;
		case 3:
			if(!AssembleDefect<Tetrahedron>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   d, u, s_m, s_a, time))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_defect':"
						"Cannot assemble Tetrahedrons.\n");
				return IAssemble_ERROR;
			}
			if(!AssembleDefect<Pyramid>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   d, u, s_m, s_a, time))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_defect':"
						"Cannot assemble Pyramids.\n");
				return IAssemble_ERROR;
			}
			if(!AssembleDefect<Prism>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   d, u, s_m, s_a, time))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_defect':"
						"Cannot assemble Prisms.\n");
				return IAssemble_ERROR;
			}
			if(!AssembleDefect<Hexahedron>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   d, u, s_m, s_a, time))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_defect':"
						"Cannot assemble Hexahedrons.\n");
				return IAssemble_ERROR;
			}
			break;
		default:UG_LOG("ERROR in 'DomainDiscretization::assemble_defect':"
				"Dimension " << dim << " not supported.\n");
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
assemble_linear(matrix_type& mat, vector_type& rhs,
                const vector_type& u, number time,
                const PreviousSolutions<vector_type>& prevSol,
                const dof_distribution_type& dofDistr,
                number s_m, number s_a)
{
//	Union of Subsets
	SubsetGroup unionSubsets;

//	create list of all subsets
	if(!CreateUnionOfSubsets(unionSubsets, m_vElemDisc))
	{
		UG_LOG("ERROR: Can not create union of subsets.\n");
		return IAssemble_ERROR;
	}

//	loop subsets
	for(size_t i = 0; i < unionSubsets.num_subsets(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = unionSubsets.dim(i);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(m_bForceRegGrid) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc<TAlgebra>*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, si);

	//	assemble on suitable elements
		switch(dim)
		{
		case 1:
			if(!AssembleLinear<Edge>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
			                         mat, rhs, u, s_m, s_a, time))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
						"Cannot assemble Edges.\n");
				return IAssemble_ERROR;
			}
			break;
		case 2:
			if(!AssembleLinear<Triangle>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   mat, rhs, u, s_m, s_a, time))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
						"Cannot assemble Triangles.\n");
				return IAssemble_ERROR;
			}
			if(!AssembleLinear<Quadrilateral>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   mat, rhs, u, s_m, s_a, time))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
						"Cannot assemble Quadrilaterals.\n");
				return IAssemble_ERROR;
			}
			break;
		case 3:
			if(!AssembleLinear<Tetrahedron>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   mat, rhs, u, s_m, s_a, time))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
						"Cannot assemble Tetrahedrons.\n");
				return IAssemble_ERROR;
			}
			if(!AssembleLinear<Pyramid>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   mat, rhs, u, s_m, s_a, time))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
						"Cannot assemble Pyramids.\n");
				return IAssemble_ERROR;
			}
			if(!AssembleLinear<Prism>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   mat, rhs, u, s_m, s_a, time))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
						"Cannot assemble Prisms.\n");
				return IAssemble_ERROR;
			}
			if(!AssembleLinear<Hexahedron>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   mat, rhs, u, s_m, s_a, time))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
						"Cannot assemble Hexahedrons.\n");
				return IAssemble_ERROR;
			}
			break;
		default:UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
				"Dimension " << dim << " not supported.\n");
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
assemble_solution(vector_type& u, number time, const dof_distribution_type& dofDistr)
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

}

#endif /*__H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DOMAIN_DISCRETIZATION__IMPL__*/
