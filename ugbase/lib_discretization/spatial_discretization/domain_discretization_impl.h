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
bool
DomainDiscretization<TDoFDistribution, TAlgebra>::
assemble_mass_matrix(matrix_type& M, const vector_type& u,
                     const dof_distribution_type& dofDistr)
{
//	reset matrix to zero and resize
	const size_t numDoFs = dofDistr.num_dofs();
	M.resize(0,0);
	M.resize(numDoFs, numDoFs);
	M.set(0.0);

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	const ISubsetHandler& sh = *dofDistr.get_function_pattern().get_subset_handler();
	if(!CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, sh))
	{
		UG_LOG("ERROR in 'DomainDiscretization':"
				" Can not Subset Groups and union.\n");
		return false;
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
		std::vector<IElemDisc*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	assemble on suitable elements
		switch(dim)
		{
		case 1:
			if(!AssembleMassMatrix<Edge,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									     M, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_mass_matrix':"
						"Cannot assemble Edges.\n");
				return false;
			}
			break;
		case 2:
			if(!AssembleMassMatrix<Triangle,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									     M, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_mass_matrix':"
						"Cannot assemble Triangles.\n");
				return false;
			}
			if(!AssembleMassMatrix<Quadrilateral,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									     M, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_mass_matrix':"
						"Cannot assemble Quadrilaterals.\n");
				return false;
			}
			break;
		case 3:
			if(!AssembleMassMatrix<Tetrahedron,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									     M, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_mass_matrix':"
						"Cannot assemble Tetrahedrons.\n");
				return false;
			}
			if(!AssembleMassMatrix<Pyramid,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									     M, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_mass_matrix':"
						"Cannot assemble Pyramids.\n");
				return false;
			}
			if(!AssembleMassMatrix<Prism,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									     M, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_mass_matrix':"
						"Cannot assemble Prisms.\n");
				return false;
			}
			if(!AssembleMassMatrix<Hexahedron,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									     M, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_mass_matrix':"
						"Cannot assemble Hexahedrons.\n");
				return false;
			}
			break;
		default:UG_LOG("ERROR in 'DomainDiscretization::assemble_mass_matrix':"
				"Dimension " << dim << " not supported.\n");
				return false;
		}
	}

//	post process
	for(size_t type = 0; type < PPT_NUM_POST_PROCESS_TYPES; ++type)
	{
		for(size_t i = 0; i < m_vvPostProcess[type].size(); ++i)
		{
			if(m_vvPostProcess[type][i]->post_process_jacobian(M, u, dofDistr)
					!= true)
				return false;
		}
	}

//	Remember parallel storage type
#ifdef UG_PARALLEL
	M.set_storage_type(PST_ADDITIVE);
	dof_distribution_type* dist = const_cast<dof_distribution_type*>(&dofDistr);
	CopyLayoutsAndCommunicatorIntoMatrix(M, *dist);
#endif

//	done
	return true;
}

template <typename TDoFDistribution, typename TAlgebra>
bool
DomainDiscretization<TDoFDistribution, TAlgebra>::
assemble_stiffness_matrix(matrix_type& A, const vector_type& u,
                          const dof_distribution_type& dofDistr)
{
//	reset matrix to zero and resize
	const size_t numDoFs = dofDistr.num_dofs();
	A.resize(0,0);
	A.resize(numDoFs, numDoFs);
	A.set(0.0);

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	const ISubsetHandler& sh = *dofDistr.get_function_pattern().get_subset_handler();
	if(!CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, sh))
	{
		UG_LOG("ERROR in 'DomainDiscretization':"
				" Can not Subset Groups and union.\n");
		return false;
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
		std::vector<IElemDisc*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	assemble on suitable elements
		switch(dim)
		{
		case 1:
			if(!AssembleStiffnessMatrix<Edge,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
										         A, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_stiffness_matrix':"
						"Cannot assemble Edges.\n");
				return false;
			}
			break;
		case 2:
			if(!AssembleStiffnessMatrix<Triangle,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
										         A, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_stiffness_matrix':"
						"Cannot assemble Triangles.\n");
				return false;
			}
			if(!AssembleStiffnessMatrix<Quadrilateral,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
										         A, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_stiffness_matrix':"
						"Cannot assemble Quadrilaterals.\n");
				return false;
			}
			break;
		case 3:
			if(!AssembleStiffnessMatrix<Tetrahedron,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
										         A, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_stiffness_matrix':"
						"Cannot assemble Tetrahedrons.\n");
				return false;
			}
			if(!AssembleStiffnessMatrix<Pyramid,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
										         A, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_stiffness_matrix':"
						"Cannot assemble Pyramids.\n");
				return false;
			}
			if(!AssembleStiffnessMatrix<Prism,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
										         A, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_stiffness_matrix':"
						"Cannot assemble Prisms.\n");
				return false;
			}
			if(!AssembleStiffnessMatrix<Hexahedron,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
										         A, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_stiffness_matrix':"
						"Cannot assemble Hexahedrons.\n");
				return false;
			}
			break;
		default:UG_LOG("ERROR in 'DomainDiscretization::assemble_stiffness_matrix':"
				"Dimension " << dim << " not supported.\n");
				return false;
		}
	}

//	post process
	for(size_t type = 0; type < PPT_NUM_POST_PROCESS_TYPES; ++type)
	{
		for(size_t i = 0; i < m_vvPostProcess[type].size(); ++i)
		{
			if(m_vvPostProcess[type][i]->post_process_jacobian(A, u, dofDistr)
					!= true)
				return false;
		}
	}

//	Remember parallel storage type
#ifdef UG_PARALLEL
	A.set_storage_type(PST_ADDITIVE);
	dof_distribution_type* dist = const_cast<dof_distribution_type*>(&dofDistr);
	CopyLayoutsAndCommunicatorIntoMatrix(A, *dist);
#endif

//	done
	return true;
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
bool
DomainDiscretization<TDoFDistribution, TAlgebra>::
assemble_jacobian(matrix_type& J,
                  const vector_type& u,
                  const dof_distribution_type& dofDistr)
{
//	reset matrix to zero and resize
	const size_t numDoFs = dofDistr.num_dofs();
	J.resize(0,0);
	J.resize(numDoFs, numDoFs);
	J.set(0.0);

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	const ISubsetHandler& sh = *dofDistr.get_function_pattern().get_subset_handler();
	if(!CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, sh))
	{
		UG_LOG("ERROR in 'DomainDiscretization':"
				" Can not Subset Groups and union.\n");
		return false;
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
		std::vector<IElemDisc*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	assemble on suitable elements
		switch(dim)
		{
		case 1:
			if(!AssembleJacobian<Edge,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
			                           J, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
						"Cannot assemble Edges.\n");
				return false;
			}
			break;
		case 2:
			if(!AssembleJacobian<Triangle,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
			                           J, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
						"Cannot assemble Triangles.\n");
				return false;
			}
			if(!AssembleJacobian<Quadrilateral,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
			                           J, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
						"Cannot assemble Quadrilaterals.\n");
				return false;
			}
			break;
		case 3:
			if(!AssembleJacobian<Tetrahedron,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
			                           J, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
						"Cannot assemble Tetrahedrons.\n");
				return false;
			}
			if(!AssembleJacobian<Pyramid,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
			                           J, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
						"Cannot assemble Pyramids.\n");
				return false;
			}
			if(!AssembleJacobian<Prism,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
			                           J, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
						"Cannot assemble Prisms.\n");
				return false;
			}
			if(!AssembleJacobian<Hexahedron,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
			                           J, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
						"Cannot assemble Hexahedrons.\n");
				return false;
			}
			break;
		default:UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
				"Dimension " << dim << " not supported.\n");
				return false;
		}
	}

//	post process
	for(size_t type = 0; type < PPT_NUM_POST_PROCESS_TYPES; ++type)
	{
		for(size_t i = 0; i < m_vvPostProcess[type].size(); ++i)
		{
			if(m_vvPostProcess[type][i]->post_process_jacobian(J, u, dofDistr)
					!= true)
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
						" Cannot execute post process " << i << ".\n");
				return false;
			}
		}
	}

//	Remember parallel storage type
#ifdef UG_PARALLEL
	J.set_storage_type(PST_ADDITIVE);
	dof_distribution_type* dist = const_cast<dof_distribution_type*>(&dofDistr);
	CopyLayoutsAndCommunicatorIntoMatrix(J, *dist);
#endif

//	done
	return true;
}


//////////////////////////
// assemble Defect
//////////////////////////
template <typename TDoFDistribution, typename TAlgebra>
bool
DomainDiscretization<TDoFDistribution, TAlgebra>::
assemble_defect(vector_type& d,
                const vector_type& u,
                const dof_distribution_type& dofDistr)
{
//	reset matrix to zero and resize
	const size_t numDoFs = dofDistr.num_dofs();
	d.resize(numDoFs);
	d.set(0.0);

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	const ISubsetHandler& sh = *dofDistr.get_function_pattern().get_subset_handler();
	if(!CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, sh))
	{
		UG_LOG("ERROR in 'DomainDiscretization':"
				" Can not Subset Groups and union.\n");
		return false;
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
		std::vector<IElemDisc*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	assemble on suitable elements
		switch(dim)
		{
		case 1:
			if(!AssembleDefect<Edge,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   d, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_defect':"
						"Cannot assemble Edges.\n");
				return false;
			}
			break;
		case 2:
			if(!AssembleDefect<Triangle,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   d, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_defect':"
						"Cannot assemble Triangles.\n");
				return false;
			}
			if(!AssembleDefect<Quadrilateral,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   d, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_defect':"
						"Cannot assemble Quadrilaterals.\n");
				return false;
			}
			break;
		case 3:
			if(!AssembleDefect<Tetrahedron,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   d, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_defect':"
						"Cannot assemble Tetrahedrons.\n");
				return false;
			}
			if(!AssembleDefect<Pyramid,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   d, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_defect':"
						"Cannot assemble Pyramids.\n");
				return false;
			}
			if(!AssembleDefect<Prism,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   d, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_defect':"
						"Cannot assemble Prisms.\n");
				return false;
			}
			if(!AssembleDefect<Hexahedron,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   d, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_defect':"
						"Cannot assemble Hexahedrons.\n");
				return false;
			}
			break;
		default:UG_LOG("ERROR in 'DomainDiscretization::assemble_defect':"
				"Dimension " << dim << " not supported.\n");
				return false;
		}
	}

//	post process
	for(size_t type = 0; type < PPT_NUM_POST_PROCESS_TYPES; ++type)
	{
		for(size_t i = 0; i < m_vvPostProcess[type].size(); ++i)
		{
			if(m_vvPostProcess[type][i]->post_process_defect(d, u, dofDistr)
					!= true)
				return false;
		}
	}

//	Remember parallel storage type
#ifdef UG_PARALLEL
	d.set_storage_type(PST_ADDITIVE);
#endif

//	done
	return true;
}

//////////////////////////
// assemble Matrix and RHS
//////////////////////////
template <typename TDoFDistribution, typename TAlgebra>
bool
DomainDiscretization<TDoFDistribution, TAlgebra>::
assemble_linear(matrix_type& mat, vector_type& rhs,
                const vector_type& u,
                const dof_distribution_type& dofDistr)
{
//	reset matrix to zero and resize
	const size_t numDoFs = dofDistr.num_dofs();
	mat.resize(0,0);
	mat.resize(numDoFs, numDoFs);
	mat.set(0.0);

	rhs.resize(numDoFs);
	rhs.set(0.0);

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	const ISubsetHandler& sh = *dofDistr.get_function_pattern().get_subset_handler();
	if(!CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, sh))
	{
		UG_LOG("ERROR in 'DomainDiscretization':"
				" Can not Subset Groups and union.\n");
		return false;
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
		std::vector<IElemDisc*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	assemble on suitable elements
		switch(dim)
		{
		case 1:
			if(!AssembleLinear<Edge,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   mat, rhs, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
						"Cannot assemble Edges.\n");
				return false;
			}
			break;
		case 2:
			if(!AssembleLinear<Triangle,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   mat, rhs, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
						"Cannot assemble Triangles.\n");
				return false;
			}
			if(!AssembleLinear<Quadrilateral,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   mat, rhs, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
						"Cannot assemble Quadrilaterals.\n");
				return false;
			}
			break;
		case 3:
			if(!AssembleLinear<Tetrahedron,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   mat, rhs, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
						"Cannot assemble Tetrahedrons.\n");
				return false;
			}
			if(!AssembleLinear<Pyramid,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   mat, rhs, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
						"Cannot assemble Pyramids.\n");
				return false;
			}
			if(!AssembleLinear<Prism,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   mat, rhs, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
						"Cannot assemble Prisms.\n");
				return false;
			}
			if(!AssembleLinear<Hexahedron,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   mat, rhs, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
						"Cannot assemble Hexahedrons.\n");
				return false;
			}
			break;
		default:UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
				"Dimension " << dim << " not supported.\n");
				return false;
		}
	}

//	post process
	for(size_t type = 0; type < PPT_NUM_POST_PROCESS_TYPES; ++type)
	{
		for(size_t i = 0; i < m_vvPostProcess[type].size(); ++i)
		{
			if(m_vvPostProcess[type][i]->post_process_linear(mat, rhs, u, dofDistr)
					!= true)
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
						" Cannot post process.\n");
				return false;
			}
		}
	}

//	Remember parallel storage type
#ifdef UG_PARALLEL
	mat.set_storage_type(PST_ADDITIVE);
	dof_distribution_type* dist = const_cast<dof_distribution_type*>(&dofDistr);
	CopyLayoutsAndCommunicatorIntoMatrix(mat, *dist);

	rhs.set_storage_type(PST_ADDITIVE);
#endif

//	done
	return true;
}

/////////////////////////
// assemble Matrix and RHS
//////////////////////////
template <typename TDoFDistribution, typename TAlgebra>
bool
DomainDiscretization<TDoFDistribution, TAlgebra>::
assemble_rhs(vector_type& rhs,
			const vector_type& u,
			const dof_distribution_type& dofDistr)
{
//	reset matrix to zero and resize
	const size_t numDoFs = dofDistr.num_dofs();
	rhs.resize(numDoFs);
	rhs.set(0.0);

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	const ISubsetHandler& sh = *dofDistr.get_function_pattern().get_subset_handler();
	if(!CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, sh))
	{
		UG_LOG("ERROR in 'DomainDiscretization':"
				" Can not Subset Groups and union.\n");
		return false;
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
		std::vector<IElemDisc*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	assemble on suitable elements
		switch(dim)
		{
		case 1:
			if(!AssembleRhs<Edge,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   rhs, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_rhs':"
						"Cannot assemble Edges.\n");
				return false;
			}
			break;
		case 2:
			if(!AssembleRhs<Triangle,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   rhs, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_rhs':"
						"Cannot assemble Triangles.\n");
				return false;
			}
			if(!AssembleRhs<Quadrilateral,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   rhs, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_rhs':"
						"Cannot assemble Quadrilaterals.\n");
				return false;
			}
			break;
		case 3:
			if(!AssembleRhs<Tetrahedron,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   rhs, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_rhs':"
						"Cannot assemble Tetrahedrons.\n");
				return false;
			}
			if(!AssembleRhs<Pyramid,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   rhs, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_rhs':"
						"Cannot assemble Pyramids.\n");
				return false;
			}
			if(!AssembleRhs<Prism,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   rhs, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_rhs':"
						"Cannot assemble Prisms.\n");
				return false;
			}
			if(!AssembleRhs<Hexahedron,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   rhs, u, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_rhs':"
						"Cannot assemble Hexahedrons.\n");
				return false;
			}
			break;
		default:UG_LOG("ERROR in 'DomainDiscretization::assemble_rhs':"
				"Dimension " << dim << " not supported.\n");
				return false;
		}
	}

//	post process
	for(size_t type = 0; type < PPT_NUM_POST_PROCESS_TYPES; ++type)
	{
		for(size_t i = 0; i < m_vvPostProcess[type].size(); ++i)
		{
			if(m_vvPostProcess[type][i]->post_process_rhs(rhs, u, dofDistr)
					!= true)
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_rhs':"
						" Cannot post process.\n");
				return false;
			}
		}
	}

//	Remember parallel storage type
#ifdef UG_PARALLEL
	rhs.set_storage_type(PST_ADDITIVE);
#endif

//	done
	return true;
}

//////////////////////////////////
// set solution in dirichlet dofs
//////////////////////////////////
template <typename TDoFDistribution, typename TAlgebra>
bool
DomainDiscretization<TDoFDistribution, TAlgebra>::
assemble_solution(vector_type& u, const dof_distribution_type& dofDistr)
{
//	post process dirichlet
	for(size_t i = 0; i < m_vvPostProcess[PPT_DIRICHLET].size(); ++i)
	{
		if(m_vvPostProcess[PPT_DIRICHLET][i]->post_process_solution(u, dofDistr)
				!= true)
			return false;
	}

//	post process constraints
	for(size_t i = 0; i < m_vvPostProcess[PPT_CONSTRAINTS].size(); ++i)
	{
		if(m_vvPostProcess[PPT_CONSTRAINTS][i]->post_process_solution(u, dofDistr)
				!= true)
			return false;
	}

//	done
	return true;
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
bool
DomainDiscretization<TDoFDistribution, TAlgebra>::
assemble_jacobian(matrix_type& J,
                  const vector_type& u, number time,
                  const SolutionTimeSeries<vector_type>& solList,
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
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	const ISubsetHandler& sh = *dofDistr.get_function_pattern().get_subset_handler();
	if(!CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, sh))
	{
		UG_LOG("ERROR in 'DomainDiscretization':"
				" Can not Subset Groups and union.\n");
		return false;
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
		std::vector<IElemDisc*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	assemble on suitable elements
		switch(dim)
		{
		case 1:
			if(!AssembleJacobian<Edge,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   J, u, time, solList, s_a0, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
						"Cannot assemble Edges.\n");
				return false;
			}
			break;
		case 2:
			if(!AssembleJacobian<Triangle,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   J, u, time, solList, s_a0, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
						"Cannot assemble Triangles.\n");
				return false;
			}
			if(!AssembleJacobian<Quadrilateral,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   J, u, time, solList, s_a0, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
						"Cannot assemble Quadrilaterals.\n");
				return false;
			}
			break;
		case 3:
			if(!AssembleJacobian<Tetrahedron,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   J, u, time, solList, s_a0, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
						"Cannot assemble Tetrahedrons.\n");
				return false;
			}
			if(!AssembleJacobian<Pyramid,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   J, u, time, solList, s_a0, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
						"Cannot assemble Pyramids.\n");
				return false;
			}
			if(!AssembleJacobian<Prism,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   J, u, time, solList, s_a0, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
						"Cannot assemble Prisms.\n");
				return false;
			}
			if(!AssembleJacobian<Hexahedron,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   J, u, time, solList, s_a0, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
						"Cannot assemble Hexahedrons.\n");
				return false;
			}
			break;
		default:UG_LOG("ERROR in 'DomainDiscretization::assemble_jacobian':"
				"Dimension " << dim << " not supported.\n");
				return false;
		}
	}

//	post process
	for(size_t type = 0; type < PPT_NUM_POST_PROCESS_TYPES; ++type)
	{
		for(size_t i = 0; i < m_vvPostProcess[type].size(); ++i)
		{
			if(m_vvPostProcess[type][i]->post_process_jacobian(J, u, dofDistr, time)
					!= true)
				return false;
		}
	}

//	Remember parallel storage type
#ifdef UG_PARALLEL
	J.set_storage_type(PST_ADDITIVE);
	dof_distribution_type* dist = const_cast<dof_distribution_type*>(&dofDistr);
	CopyLayoutsAndCommunicatorIntoMatrix(J, *dist);
#endif

//	done
	return true;
}

//////////////////////////////////
// Assemble Defect
//////////////////////////////////
template <typename TDoFDistribution, typename TAlgebra>
bool
DomainDiscretization<TDoFDistribution, TAlgebra>::
assemble_defect(vector_type& d,
                const vector_type& u, number time,
                const SolutionTimeSeries<vector_type>& solList,
                const dof_distribution_type& dofDistr,
                number s_m, number s_a)
{
//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	const ISubsetHandler& sh = *dofDistr.get_function_pattern().get_subset_handler();
	if(!CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, sh))
	{
		UG_LOG("ERROR in 'DomainDiscretization':"
				" Can not Subset Groups and union.\n");
		return false;
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
		std::vector<IElemDisc*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	assemble on suitable elements
		switch(dim)
		{
		case 1:
			if(!AssembleDefect<Edge,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   d, u, time, solList, s_m, s_a, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_defect':"
						"Cannot assemble Edges.\n");
				return false;
			}
			break;
		case 2:
			if(!AssembleDefect<Triangle,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   d, u, time, solList, s_m, s_a, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_defect':"
						"Cannot assemble Triangles.\n");
				return false;
			}
			if(!AssembleDefect<Quadrilateral,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   d, u, time, solList, s_m, s_a, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_defect':"
						"Cannot assemble Quadrilaterals.\n");
				return false;
			}
			break;
		case 3:
			if(!AssembleDefect<Tetrahedron,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   d, u, time, solList, s_m, s_a, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_defect':"
						"Cannot assemble Tetrahedrons.\n");
				return false;
			}
			if(!AssembleDefect<Pyramid,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   d, u, time, solList, s_m, s_a, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_defect':"
						"Cannot assemble Pyramids.\n");
				return false;
			}
			if(!AssembleDefect<Prism,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   d, u, time, solList, s_m, s_a, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_defect':"
						"Cannot assemble Prisms.\n");
				return false;
			}
			if(!AssembleDefect<Hexahedron,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   d, u, time, solList, s_m, s_a, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_defect':"
						"Cannot assemble Hexahedrons.\n");
				return false;
			}
			break;
		default:UG_LOG("ERROR in 'DomainDiscretization::assemble_defect':"
				"Dimension " << dim << " not supported.\n");
				return false;
		}
	}

//	post process
	for(size_t type = 0; type < PPT_NUM_POST_PROCESS_TYPES; ++type)
	{
		for(size_t i = 0; i < m_vvPostProcess[type].size(); ++i)
		{
			if(m_vvPostProcess[type][i]->post_process_defect(d, u, dofDistr, time)
					!= true)
				return false;
		}
	}

//	Remember parallel storage type
#ifdef UG_PARALLEL
	d.set_storage_type(PST_ADDITIVE);
#endif

//	done
	return true;
}

//////////////////////////////////
// Assemble Matrix and RHS
//////////////////////////////////
template <typename TDoFDistribution, typename TAlgebra>
bool
DomainDiscretization<TDoFDistribution, TAlgebra>::
assemble_linear(matrix_type& mat, vector_type& rhs,
                const vector_type& u, number time,
                const SolutionTimeSeries<vector_type>& solList,
                const dof_distribution_type& dofDistr,
                number s_m, number s_a)
{
//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	const ISubsetHandler& sh = *dofDistr.get_function_pattern().get_subset_handler();
	if(!CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, sh))
	{
		UG_LOG("ERROR in 'DomainDiscretization':"
				" Can not Subset Groups and union.\n");
		return false;
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
		std::vector<IElemDisc*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	assemble on suitable elements
		switch(dim)
		{
		case 1:
			if(!AssembleLinear<Edge,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
			                         mat, rhs, u, time, solList, s_m, s_a, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
						"Cannot assemble Edges.\n");
				return false;
			}
			break;
		case 2:
			if(!AssembleLinear<Triangle,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   mat, rhs, u, time, solList, s_m, s_a, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
						"Cannot assemble Triangles.\n");
				return false;
			}
			if(!AssembleLinear<Quadrilateral,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   mat, rhs, u, time, solList, s_m, s_a, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
						"Cannot assemble Quadrilaterals.\n");
				return false;
			}
			break;
		case 3:
			if(!AssembleLinear<Tetrahedron,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   mat, rhs, u, time, solList, s_m, s_a, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
						"Cannot assemble Tetrahedrons.\n");
				return false;
			}
			if(!AssembleLinear<Pyramid,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   mat, rhs, u, time, solList, s_m, s_a, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
						"Cannot assemble Pyramids.\n");
				return false;
			}
			if(!AssembleLinear<Prism,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   mat, rhs, u, time, solList, s_m, s_a, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
						"Cannot assemble Prisms.\n");
				return false;
			}
			if(!AssembleLinear<Hexahedron,TDoFDistribution,TAlgebra>(vSubsetElemDisc, dofDistr, si, bNonRegularGrid,
									   mat, rhs, u, time, solList, s_m, s_a, m_pSelector))
			{
				UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
						"Cannot assemble Hexahedrons.\n");
				return false;
			}
			break;
		default:UG_LOG("ERROR in 'DomainDiscretization::assemble_linear':"
				"Dimension " << dim << " not supported.\n");
				return false;
		}
	}


//	post process
	for(size_t type = 0; type < PPT_NUM_POST_PROCESS_TYPES; ++type)
	{
		for(size_t i = 0; i < m_vvPostProcess[type].size(); ++i)
		{
			if(m_vvPostProcess[type][i]->post_process_linear(mat, rhs, u, dofDistr, time)
					!= true)
				return false;
		}
	}

//	Remember parallel storage type
#ifdef UG_PARALLEL
	mat.set_storage_type(PST_ADDITIVE);
	dof_distribution_type* dist = const_cast<dof_distribution_type*>(&dofDistr);
	CopyLayoutsAndCommunicatorIntoMatrix(mat, *dist);

	rhs.set_storage_type(PST_ADDITIVE);
#endif


//	done
	return true;
}

//////////////////////////////////
// set dirichlet values
//////////////////////////////////
template <typename TDoFDistribution, typename TAlgebra>
bool
DomainDiscretization<TDoFDistribution, TAlgebra>::
assemble_solution(vector_type& u, number time, const dof_distribution_type& dofDistr)
{
//	dirichlet
	for(size_t i = 0; i < m_vvPostProcess[PPT_DIRICHLET].size(); ++i)
	{
		if(m_vvPostProcess[PPT_DIRICHLET][i]->post_process_solution(u, dofDistr, time)
				!= true)
			return false;
	}

//	constraints
	for(size_t i = 0; i < m_vvPostProcess[PPT_CONSTRAINTS].size(); ++i)
	{
		if(m_vvPostProcess[PPT_CONSTRAINTS][i]->post_process_solution(u, dofDistr, time)
				!= true)
			return false;
	}

//	done
	return true;
}

}

#endif /*__H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DOMAIN_DISCRETIZATION__IMPL__*/
