/*
 * domain_disc_impl.h
 *
 *  Created on: 29.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DOMAIN_DISC_IMPL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DOMAIN_DISC_IMPL__

#include "domain_disc.h"
#include "lib_disc/common/groups_util.h"

namespace ug{

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool DomainDiscretization<TDomain, TDoFDistribution, TAlgebra>::
update_elem_discs()
{
//	check Approximation space
	if(m_pApproxSpace == NULL)
	{
		UG_LOG("ERROR in DomainDiscretization: Before using the "
				"DomainDiscretization an ApproximationSpace must be set to it. "
				"Please use DomainDiscretization:set_approximation_space to "
				"set an appropriate Space.\n");
		return false;
	}

//	set approximation space and extract IElemDiscs
	m_vElemDisc.clear();
	for(size_t i = 0; i < m_vDomainElemDisc.size(); ++i)
	{
		m_vDomainElemDisc[i]->set_approximation_space(*m_pApproxSpace);
		m_vElemDisc.push_back(m_vDomainElemDisc[i]);
	}

//	ok
	return true;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool DomainDiscretization<TDomain, TDoFDistribution, TAlgebra>::
update_constraints()
{
//	check Approximation space
	if(m_pApproxSpace == NULL)
	{
		UG_LOG("ERROR in DomainDiscretization: Before using the "
				"DomainDiscretization an ApproximationSpace must be set to it. "
				"Please use DomainDiscretization:set_approximation_space to "
				"set an appropriate Space.\n");
		return false;
	}


	for(size_t type = 0; type < NUM_CONSTRAINT_TYPES; ++type)
	{
		for(size_t i = 0; i < m_vvConstraints[type].size(); ++i)
		{
			m_vvConstraints[type][i]->set_approximation_space(*m_pApproxSpace);
		}
	}

//	ok
	return true;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool DomainDiscretization<TDomain, TDoFDistribution, TAlgebra>::
update_disc_items()
{
//	return flag
	bool bRet = true;

	bRet &= update_elem_discs();
	bRet &= update_constraints();

//	ok
	return bRet;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
typename DomainDiscretization<TDomain, TDoFDistribution, TAlgebra>::dof_distribution_type&
DomainDiscretization<TDomain, TDoFDistribution, TAlgebra>::
get_surface_dd()
{
//	check Approximation space
	if(m_pApproxSpace == NULL)
	{
		UG_LOG("ERROR in DomainDiscretization: Before using the "
				"DomainDiscretization an ApproximationSpace must be set to it. "
				"Please use DomainDiscretization:set_approximation_space to "
				"set an appropriate Space.\n");
		throw(UGFatalError("ApproximationSpace missing in DomainDiscretization."));
	}

	return m_pApproxSpace->surface_dof_distribution();
}


///////////////////////////////////////////////////////////////////////////////
// Mass Matrix
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
void DomainDiscretization<TDomain, TDoFDistribution, TAlgebra>::
assemble_mass_matrix(matrix_type& M, const vector_type& u,
                     const dof_distribution_type& dd)
{
//	update the elem discs
	if(!update_disc_items())
		UG_THROW_FATAL("DomainDisc: Cannot update disc items.")

//	reset matrix to zero and resize
	const size_t numIndex = dd.num_indices();
	M.resize(0,0);
	M.resize(numIndex, numIndex);
	M.set(0.0);

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	const ISubsetHandler& sh = dd.get_function_pattern().subset_handler();
	if(!CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, sh))
		UG_THROW_FATAL("ERROR in 'DomainDiscretization':"
					   " Can not Subset Groups and union.\n");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.num_subsets(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = subset_dim(unionSubsets, i);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(m_bForceRegGrid) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	success flag
		bool bSuc = true;

	//	assemble on suitable elements
		switch(dim)
		{
		case 1:
			bSuc &= AssembleMassMatrix<Edge,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u, m_pSelector);
			break;
		case 2:
			bSuc &= AssembleMassMatrix<Triangle,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u, m_pSelector);
			bSuc &= AssembleMassMatrix<Quadrilateral,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u, m_pSelector);
			break;
		case 3:
			bSuc &= AssembleMassMatrix<Tetrahedron,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u, m_pSelector);
			bSuc &= AssembleMassMatrix<Pyramid,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u, m_pSelector);
			bSuc &= AssembleMassMatrix<Prism,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u, m_pSelector);
			bSuc &= AssembleMassMatrix<Hexahedron,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u, m_pSelector);
			break;
		default:
			UG_THROW_FATAL("ERROR in 'DomainDiscretization::assemble_mass_matrix':"
							"Dimension " << dim << " (subset="<<si<<") not supported.\n");
		}

	//	check success
		if(!bSuc)
			UG_THROW_FATAL("ERROR in 'DomainDiscretization::assemble_mass_matrix':"
						" Assembling of elements of Dimension " << dim << " in "
						" subset "<<si<< " failed.\n");
	}

//	post process
	for(size_t type = 0; type < NUM_CONSTRAINT_TYPES; ++type)
	{
		for(size_t i = 0; i < m_vvConstraints[type].size(); ++i)
		{
			if(!m_vvConstraints[type][i]->adjust_jacobian(M, u, dd))
				UG_THROW_FATAL("DomainDiscretization: Cannot adjust jacobian.");
		}
	}

//	Remember parallel storage type
#ifdef UG_PARALLEL
	M.set_storage_type(PST_ADDITIVE);
	dof_distribution_type* dist = const_cast<dof_distribution_type*>(&dd);
	CopyLayoutsAndCommunicatorIntoMatrix(M, *dist);
#endif
}

///////////////////////////////////////////////////////////////////////////////
// Stiffness Matrix
///////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
void DomainDiscretization<TDomain, TDoFDistribution, TAlgebra>::
assemble_stiffness_matrix(matrix_type& A, const vector_type& u,
                          const dof_distribution_type& dd)
{
//	update the elem discs
	if(!update_disc_items())
		UG_THROW_FATAL("Cannot update disc items.");

//	reset matrix to zero and resize
	const size_t numIndex = dd.num_indices();
	A.resize(0,0);
	A.resize(numIndex, numIndex);
	A.set(0.0);

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	const ISubsetHandler& sh = dd.get_function_pattern().subset_handler();
	if(!CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, sh))
		UG_THROW_FATAL("ERROR in 'DomainDiscretization':"
						" Can not Subset Groups and union.\n");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.num_subsets(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
	//	get dimension of the subset
		const int dim = subset_dim(unionSubsets, i);


	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(m_bForceRegGrid) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	success flag
		bool bSuc = true;

	//	assemble on suitable elements
		switch(dim)
		{
		case 1:
			bSuc &= AssembleStiffnessMatrix<Edge,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u, m_pSelector);
			break;
		case 2:
			bSuc &= AssembleStiffnessMatrix<Triangle,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u, m_pSelector);
			bSuc &= AssembleStiffnessMatrix<Quadrilateral,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u, m_pSelector);
			break;
		case 3:
			bSuc &= AssembleStiffnessMatrix<Tetrahedron,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u, m_pSelector);
			bSuc &= AssembleStiffnessMatrix<Pyramid,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u, m_pSelector);
			bSuc &= AssembleStiffnessMatrix<Prism,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u, m_pSelector);
			bSuc &= AssembleStiffnessMatrix<Hexahedron,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u, m_pSelector);
			break;
		default:
			UG_THROW_FATAL("ERROR in 'DomainDiscretization::assemble_stiffness_matrix':"
							"Dimension " << dim << " (subset="<<si<<") not supported.\n");
		}

	//	check success
		if(!bSuc)
			UG_THROW_FATAL("ERROR in 'DomainDiscretization::assemble_stiffness_matrix':"
					" Assembling of elements of Dimension " << dim << " in "
					" subset "<<si<< " failed.\n");
	}

//	post process
	for(size_t type = 0; type < NUM_CONSTRAINT_TYPES; ++type)
	{
		for(size_t i = 0; i < m_vvConstraints[type].size(); ++i)
		{
			if(!m_vvConstraints[type][i]->adjust_jacobian(A, u, dd))
				UG_THROW_FATAL("Cannot adjust jacobian.");
		}
	}

//	Remember parallel storage type
#ifdef UG_PARALLEL
	A.set_storage_type(PST_ADDITIVE);
	dof_distribution_type* dist = const_cast<dof_distribution_type*>(&dd);
	CopyLayoutsAndCommunicatorIntoMatrix(A, *dist);
#endif
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//  Time Independent (stationary)
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// Jacobian (stationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
void DomainDiscretization<TDomain, TDoFDistribution, TAlgebra>::
assemble_jacobian(matrix_type& J,
                  const vector_type& u,
                  const dof_distribution_type& dd)
{
//	update the elem discs
	if(!update_disc_items())
		UG_THROW_FATAL("Cannot update disc items.");

//	reset matrix to zero and resize
	const size_t numIndex = dd.num_indices();
	J.resize(0,0);
	J.resize(numIndex, numIndex);
	J.set(0.0);

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	const ISubsetHandler& sh = dd.get_function_pattern().subset_handler();
	if(!CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, sh))
		UG_THROW_FATAL("ERROR in 'DomainDiscretization':"
						" Can not Subset Groups and union.\n");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.num_subsets(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = subset_dim(unionSubsets, i);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(m_bForceRegGrid) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	success flag
		bool bSuc = true;

	//	assemble on suitable elements
		switch(dim)
		{
		case 1:
			bSuc &= AssembleJacobian<Edge,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, u, m_pSelector);
			break;
		case 2:
			bSuc &= AssembleJacobian<Triangle,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, u, m_pSelector);
			bSuc &= AssembleJacobian<Quadrilateral,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, u, m_pSelector);
			break;
		case 3:
			bSuc &= AssembleJacobian<Tetrahedron,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, u, m_pSelector);
			bSuc &= AssembleJacobian<Pyramid,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, u, m_pSelector);
			bSuc &= AssembleJacobian<Prism,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, u, m_pSelector);
			bSuc &= AssembleJacobian<Hexahedron,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, u, m_pSelector);
			break;
		default:
			UG_THROW_FATAL("ERROR in 'DomainDiscretization::assemble_jacobian (stationary)':"
							"Dimension " << dim << " (subset="<<si<<") not supported.\n");
		}

	//	check success
		if(!bSuc)
			UG_THROW_FATAL("ERROR in 'DomainDiscretization::assemble_jacobian (stationary)':"
							" Assembling of elements of Dimension " << dim << " in "
							" subset "<<si<< " failed.\n");
	}

//	post process
	for(size_t type = 0; type < NUM_CONSTRAINT_TYPES; ++type)
	{
		for(size_t i = 0; i < m_vvConstraints[type].size(); ++i)
		{
			if(!m_vvConstraints[type][i]->adjust_jacobian(J, u, dd))
				UG_THROW_FATAL("ERROR in 'DomainDiscretization::assemble_jacobian':"
								" Cannot execute post process " << i << ".\n");
		}
	}
//	Remember parallel storage type
#ifdef UG_PARALLEL
	J.set_storage_type(PST_ADDITIVE);
	dof_distribution_type* dist = const_cast<dof_distribution_type*>(&dd);
	CopyLayoutsAndCommunicatorIntoMatrix(J, *dist);
#endif
}


///////////////////////////////////////////////////////////////////////////////
// Defect (stationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
void DomainDiscretization<TDomain, TDoFDistribution, TAlgebra>::
assemble_defect(vector_type& d,
                const vector_type& u,
                const dof_distribution_type& dd)
{
//	update the elem discs
	if(!update_disc_items())
		UG_THROW_FATAL("Cannot update disc items.");

//	reset matrix to zero and resize
	const size_t numIndex = dd.num_indices();
	d.resize(numIndex);
	d.set(0.0);

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	const ISubsetHandler& sh = dd.get_function_pattern().subset_handler();
	if(!CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, sh))
		UG_THROW_FATAL("ERROR in 'DomainDiscretization':"
						" Can not Subset Groups and union.\n");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.num_subsets(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = subset_dim(unionSubsets, i);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(m_bForceRegGrid) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	success flag
		bool bSuc = true;

	//	assemble on suitable elements
		switch(dim)
		{
		case 1:
			bSuc &= AssembleDefect<Edge,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, u, m_pSelector);
			break;
		case 2:
			bSuc &= AssembleDefect<Triangle,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, u, m_pSelector);
			bSuc &= AssembleDefect<Quadrilateral,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, u, m_pSelector);
			break;
		case 3:
			bSuc &= AssembleDefect<Tetrahedron,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, u, m_pSelector);
			bSuc &= AssembleDefect<Pyramid,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, u, m_pSelector);
			bSuc &= AssembleDefect<Prism,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, u, m_pSelector);
			bSuc &= AssembleDefect<Hexahedron,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, u, m_pSelector);
			break;
		default:
			UG_THROW_FATAL("ERROR in 'DomainDiscretization::assemble_defect (stationary)':"
						"Dimension " << dim << " (subset="<<si<<") not supported.\n");
		}

	//	check success
		if(!bSuc)
			UG_THROW_FATAL("ERROR in 'DomainDiscretization::assemble_defect (stationary)':"
							" Assembling of elements of Dimension " << dim << " in "
							" subset "<<si<< " failed.\n");
	}

//	post process
	for(size_t type = 0; type < NUM_CONSTRAINT_TYPES; ++type)
	{
		for(size_t i = 0; i < m_vvConstraints[type].size(); ++i)
		{
			if(!m_vvConstraints[type][i]->adjust_defect(d, u, dd))
				UG_THROW_FATAL("Cannot adjust defect.");
		}
	}

//	Remember parallel storage type
#ifdef UG_PARALLEL
	d.set_storage_type(PST_ADDITIVE);
#endif
}

///////////////////////////////////////////////////////////////////////////////
// Matrix and RHS (stationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
void DomainDiscretization<TDomain, TDoFDistribution, TAlgebra>::
assemble_linear(matrix_type& mat, vector_type& rhs,
                const dof_distribution_type& dd)
{
//	update the elem discs
	if(!update_disc_items())
		UG_THROW_FATAL("Cannot update disc items.");

//	reset matrix to zero and resize
	const size_t numIndex = dd.num_indices();
	mat.resize(0,0);
	mat.resize(numIndex, numIndex);
	mat.set(0.0);

	rhs.resize(numIndex);
	rhs.set(0.0);

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	const ISubsetHandler& sh = dd.get_function_pattern().subset_handler();
	if(!CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, sh))
		UG_THROW_FATAL("ERROR in 'DomainDiscretization':"
						" Can not Subset Groups and union.\n");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.num_subsets(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = subset_dim(unionSubsets, i);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(m_bForceRegGrid) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	success flag
		bool bSuc = true;

	//	assemble on suitable elements
		switch(dim)
		{
		case 1:
			bSuc &= AssembleLinear<Edge,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, m_pSelector);
			break;
		case 2:
			bSuc &= AssembleLinear<Triangle,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, m_pSelector);
			bSuc &= AssembleLinear<Quadrilateral,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, m_pSelector);
			break;
		case 3:
			bSuc &= AssembleLinear<Tetrahedron,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, m_pSelector);
			bSuc &= AssembleLinear<Pyramid,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, m_pSelector);
			bSuc &= AssembleLinear<Prism,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, m_pSelector);
			bSuc &= AssembleLinear<Hexahedron,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, m_pSelector);
			break;
		default:
			UG_THROW_FATAL("ERROR in 'DomainDiscretization::assemble_linear (stationary)':"
							"Dimension "<<dim<<" (subset="<<si<<") not supported.\n");
		}

	//	check success
		if(!bSuc)
			UG_THROW_FATAL("ERROR in 'DomainDiscretization::assemble_linear (stationary)':"
							" Assembling of elements of Dimension " << dim << " in "
							" subset "<<si<< " failed.\n");
	}

//	post process
	for(size_t type = 0; type < NUM_CONSTRAINT_TYPES; ++type)
	{
		for(size_t i = 0; i < m_vvConstraints[type].size(); ++i)
		{
			if(!m_vvConstraints[type][i]->adjust_linear(mat, rhs, dd))
				UG_THROW_FATAL("ERROR in 'DomainDiscretization::assemble_linear':"
						" Cannot post process.\n");
		}
	}

//	Remember parallel storage type
#ifdef UG_PARALLEL
	mat.set_storage_type(PST_ADDITIVE);
	dof_distribution_type* dist = const_cast<dof_distribution_type*>(&dd);
	CopyLayoutsAndCommunicatorIntoMatrix(mat, *dist);

	rhs.set_storage_type(PST_ADDITIVE);
#endif
}

///////////////////////////////////////////////////////////////////////////////
// RHS (stationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
void DomainDiscretization<TDomain, TDoFDistribution, TAlgebra>::
assemble_rhs(vector_type& rhs,
			const vector_type& u,
			const dof_distribution_type& dd)
{
//	update the elem discs
	if(!update_disc_items())
		UG_THROW_FATAL("Cannot update disc items.");

//	reset matrix to zero and resize
	const size_t numIndex = dd.num_indices();
	rhs.resize(numIndex);
	rhs.set(0.0);

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	const ISubsetHandler& sh = dd.get_function_pattern().subset_handler();
	if(!CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, sh))
		UG_THROW_FATAL("ERROR in 'DomainDiscretization':"
						" Can not Subset Groups and union.\n");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.num_subsets(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = subset_dim(unionSubsets, i);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(m_bForceRegGrid) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	success flag
		bool bSuc = true;

	//	assemble on suitable elements
		switch(dim)
		{
		case 1:
			bSuc &= AssembleRhs<Edge,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u, m_pSelector);
			break;
		case 2:
			bSuc &= AssembleRhs<Triangle,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u, m_pSelector);
			bSuc &= AssembleRhs<Quadrilateral,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u, m_pSelector);
			break;
		case 3:
			bSuc &= AssembleRhs<Tetrahedron,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u, m_pSelector);
			bSuc &= AssembleRhs<Pyramid,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u, m_pSelector);
			bSuc &= AssembleRhs<Prism,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u, m_pSelector);
			bSuc &= AssembleRhs<Hexahedron,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u, m_pSelector);
			break;
		default:
			UG_THROW_FATAL("ERROR in 'DomainDiscretization::assemble_rhs (stationary)':"
						"Dimension "<<dim<<" (subset="<<si<<") not supported.\n");
		}

	//	check success
		if(!bSuc)
			UG_THROW_FATAL("ERROR in 'DomainDiscretization::assemble_rhs (stationary)':"
							" Assembling of elements of Dimension " << dim << " in "
							" subset "<<si<< " failed.\n");
	}

//	post process
	for(size_t type = 0; type < NUM_CONSTRAINT_TYPES; ++type)
	{
		for(size_t i = 0; i < m_vvConstraints[type].size(); ++i)
		{
			if(!m_vvConstraints[type][i]->adjust_rhs(rhs, u, dd))
				UG_THROW_FATAL("ERROR in 'DomainDiscretization::assemble_rhs':"
							" Cannot post process.\n");
		}
	}

//	Remember parallel storage type
#ifdef UG_PARALLEL
	rhs.set_storage_type(PST_ADDITIVE);
#endif
}

///////////////////////////////////////////////////////////////////////////////
// set constraints (stationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
void DomainDiscretization<TDomain, TDoFDistribution, TAlgebra>::
adjust_solution(vector_type& u, const dof_distribution_type& dd)
{
	if(!update_constraints())
		UG_THROW_FATAL("Cannot update constraints.");

//	post process dirichlet
	for(size_t i = 0; i < m_vvConstraints[CT_DIRICHLET].size(); ++i)
	{
		if(!m_vvConstraints[CT_DIRICHLET][i]->adjust_solution(u, dd))
			UG_THROW_FATAL("Cannot adjust solution.");
	}

//	post process constraints
	for(size_t i = 0; i < m_vvConstraints[CT_CONSTRAINTS].size(); ++i)
	{
		if(!m_vvConstraints[CT_CONSTRAINTS][i]->adjust_solution(u, dd))
			UG_THROW_FATAL("Cannot adjust solution.");
	}
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//  Time Dependent (instationary)
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Prepare Timestep (instationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
void DomainDiscretization<TDomain, TDoFDistribution, TAlgebra>::
prepare_timestep(const VectorTimeSeries<vector_type>& vSol,
                const dof_distribution_type& dd)
{
//	update the elem discs
	if(!update_disc_items())
		UG_THROW_FATAL("Cannot update disc items.");

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	const ISubsetHandler& sh = dd.get_function_pattern().subset_handler();
	if(!CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, sh))
		UG_THROW_FATAL("ERROR in 'DomainDiscretization':"
				" Can not Subset Groups and union.\n");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.num_subsets(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = subset_dim(unionSubsets, i);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(m_bForceRegGrid) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	success flag
		bool bSuc = true;

	//	assemble on suitable elements
		switch(dim)
		{
		case 1:
			bSuc &= PrepareTimestep<Edge,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol, m_pSelector);
			break;
		case 2:
			bSuc &= PrepareTimestep<Triangle,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol, m_pSelector);
			bSuc &= PrepareTimestep<Quadrilateral,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol, m_pSelector);
			break;
		case 3:
			bSuc &= PrepareTimestep<Tetrahedron,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol, m_pSelector);
			bSuc &= PrepareTimestep<Pyramid,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol, m_pSelector);
			bSuc &= PrepareTimestep<Prism,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol, m_pSelector);
			bSuc &= PrepareTimestep<Hexahedron,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol, m_pSelector);
			break;
		default:
			UG_THROW_FATAL("ERROR in 'DomainDiscretization::prepare_timestep (instationary)':"
						"Dimension "<<dim<<" (subset="<<si<<") not supported.\n");
		}

	//	check success
		if(!bSuc)
			UG_THROW_FATAL("ERROR in 'DomainDiscretization::prepare_timestep (instationary)':"
							" Assembling of elements of Dimension " << dim << " in "
							" subset "<<si<< " failed.\n");
	}

	//hier Dirichlet-Constraints abfragen?

}

///////////////////////////////////////////////////////////////////////////////
// Jacobian (instationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
void DomainDiscretization<TDomain, TDoFDistribution, TAlgebra>::
assemble_jacobian(matrix_type& J,
                  const VectorTimeSeries<vector_type>& vSol,
                  const number s_a0,
                  const dof_distribution_type& dd)
{
//	update the elem discs
	if(!update_disc_items())
		UG_THROW_FATAL("Cannot update disc items.");

//	get current time
	const number time = vSol.time(0);

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	const ISubsetHandler& sh = dd.get_function_pattern().subset_handler();
	if(!CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, sh))
		UG_THROW_FATAL("ERROR in 'DomainDiscretization':"
					" Can not Subset Groups and union.\n");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.num_subsets(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = subset_dim(unionSubsets, i);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(m_bForceRegGrid) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	success flag
		bool bSuc = true;

	//	assemble on suitable elements
		switch(dim)
		{
		case 1:
			bSuc &= AssembleJacobian<Edge,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, vSol, s_a0, m_pSelector);
			break;
		case 2:
			bSuc &= AssembleJacobian<Triangle,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, vSol, s_a0, m_pSelector);
			bSuc &= AssembleJacobian<Quadrilateral,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, vSol, s_a0, m_pSelector);
			break;
		case 3:
			bSuc &= AssembleJacobian<Tetrahedron,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, vSol, s_a0, m_pSelector);
			bSuc &= AssembleJacobian<Pyramid,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, vSol, s_a0, m_pSelector);
			bSuc &= AssembleJacobian<Prism,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, vSol, s_a0, m_pSelector);
			bSuc &= AssembleJacobian<Hexahedron,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, vSol, s_a0, m_pSelector);
			break;
		default:
			UG_THROW_FATAL("ERROR in 'DomainDiscretization::assemble_jacobian (instationary)':"
				"Dimension " << dim << " (subset="<<si<<") not supported.\n");
		}

	//	check success
		if(!bSuc)
			UG_THROW_FATAL("ERROR in 'DomainDiscretization::assemble_jacobian (instationary)':"
					" Assembling of elements of Dimension " << dim << " in "
					" subset "<<si<< " failed.\n");
	}

//	post process
	for(size_t type = 0; type < NUM_CONSTRAINT_TYPES; ++type)
	{
		for(size_t i = 0; i < m_vvConstraints[type].size(); ++i)
		{
			if(!m_vvConstraints[type][i]->adjust_jacobian(J, vSol.solution(0), dd, time))
				UG_THROW_FATAL("Cannot adjust jacobian.");
		}
	}

//	Remember parallel storage type
#ifdef UG_PARALLEL
	J.set_storage_type(PST_ADDITIVE);
	dof_distribution_type* dist = const_cast<dof_distribution_type*>(&dd);
	CopyLayoutsAndCommunicatorIntoMatrix(J, *dist);
#endif
}

///////////////////////////////////////////////////////////////////////////////
// Defect (instationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
void DomainDiscretization<TDomain, TDoFDistribution, TAlgebra>::
assemble_defect(vector_type& d,
                const VectorTimeSeries<vector_type>& vSol,
                const std::vector<number>& vScaleMass,
                const std::vector<number>& vScaleStiff,
                const dof_distribution_type& dd)
{
//	update the elem discs
	if(!update_disc_items())
		UG_THROW_FATAL("Cannot update disc items.");

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	const ISubsetHandler& sh = dd.get_function_pattern().subset_handler();
	if(!CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, sh))
		UG_THROW_FATAL("ERROR in 'DomainDiscretization':"
				" Can not Subset Groups and union.\n");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.num_subsets(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = subset_dim(unionSubsets, i);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(m_bForceRegGrid) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	success flag
		bool bSuc = true;

	//	assemble on suitable elements
		switch(dim)
		{
		case 1:
			bSuc &= AssembleDefect<Edge,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, vSol, vScaleMass, vScaleStiff, m_pSelector);
			break;
		case 2:
			bSuc &= AssembleDefect<Triangle,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, vSol, vScaleMass, vScaleStiff, m_pSelector);
			bSuc &= AssembleDefect<Quadrilateral,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, vSol, vScaleMass, vScaleStiff, m_pSelector);
			break;
		case 3:
			bSuc &= AssembleDefect<Tetrahedron,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, vSol, vScaleMass, vScaleStiff, m_pSelector);
			bSuc &= AssembleDefect<Pyramid,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, vSol, vScaleMass, vScaleStiff, m_pSelector);
			bSuc &= AssembleDefect<Prism,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, vSol, vScaleMass, vScaleStiff, m_pSelector);
			bSuc &= AssembleDefect<Hexahedron,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, vSol, vScaleMass, vScaleStiff, m_pSelector);
			break;
		default:
			UG_THROW_FATAL("ERROR in 'DomainDiscretization::assemble_defect (instationary)':"
						"Dimension "<<dim<<" (subset="<<si<<") not supported.\n");
		}

	//	check success
		if(!bSuc)
			UG_THROW_FATAL("ERROR in 'DomainDiscretization::assemble_defect (instationary)':"
							" Assembling of elements of Dimension " << dim << " in "
							" subset "<<si<< " failed.\n");
	}

//	post process
	for(size_t type = 0; type < NUM_CONSTRAINT_TYPES; ++type)
	{
		for(size_t i = 0; i < m_vvConstraints[type].size(); ++i)
		{
			if(!m_vvConstraints[type][i]->adjust_defect(d, vSol.solution(0), dd, vSol.time(0)))
				UG_THROW_FATAL("Cannot adjust defect.");
		}
	}

//	Remember parallel storage type
#ifdef UG_PARALLEL
	d.set_storage_type(PST_ADDITIVE);
#endif
}

///////////////////////////////////////////////////////////////////////////////
// Matrix and RHS (instationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
void DomainDiscretization<TDomain, TDoFDistribution, TAlgebra>::
assemble_linear(matrix_type& mat, vector_type& rhs,
                const VectorTimeSeries<vector_type>& vSol,
                const std::vector<number>& vScaleMass,
                const std::vector<number>& vScaleStiff,
                const dof_distribution_type& dd)
{
//	update the elem discs
	if(!update_disc_items())
		UG_THROW_FATAL("Cannot update disc items.");

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	const ISubsetHandler& sh = dd.get_function_pattern().subset_handler();
	if(!CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, sh))
		UG_THROW_FATAL("ERROR in 'DomainDiscretization':"
						" Can not Subset Groups and union.\n");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.num_subsets(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = subset_dim(unionSubsets, i);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(m_bForceRegGrid) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	success flag
		bool bSuc = true;

	//	assemble on suitable elements
		switch(dim)
		{
		case 1:
			bSuc &= AssembleLinear<Edge,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff, m_pSelector);
			break;
		case 2:
			bSuc &= AssembleLinear<Triangle,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff, m_pSelector);
			bSuc &= AssembleLinear<Quadrilateral,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff, m_pSelector);
			break;
		case 3:
			bSuc &= AssembleLinear<Tetrahedron,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff, m_pSelector);
			bSuc &= AssembleLinear<Pyramid,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff, m_pSelector);
			bSuc &= AssembleLinear<Prism,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff, m_pSelector);
			bSuc &= AssembleLinear<Hexahedron,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff, m_pSelector);
			break;
		default:
			UG_THROW_FATAL("ERROR in 'DomainDiscretization::assemble_linear (instationary)':"
						"Dimension "<<dim<<" (subset="<<si<<") not supported.\n");
		}

	//	check success
		if(!bSuc)
			UG_THROW_FATAL("ERROR in 'DomainDiscretization::assemble_linear (instationary)':"
					" Assembling of elements of Dimension " << dim << " in "
					" subset "<<si<< " failed.\n");
	}


//	post process
	for(size_t type = 0; type < NUM_CONSTRAINT_TYPES; ++type)
	{
		for(size_t i = 0; i < m_vvConstraints[type].size(); ++i)
		{
			if(!m_vvConstraints[type][i]->adjust_linear(mat, rhs, dd, vSol.time(0)))
				UG_THROW_FATAL("Cannot adjust linear.");
		}
	}

//	Remember parallel storage type
#ifdef UG_PARALLEL
	mat.set_storage_type(PST_ADDITIVE);
	dof_distribution_type* dist = const_cast<dof_distribution_type*>(&dd);
	CopyLayoutsAndCommunicatorIntoMatrix(mat, *dist);

	rhs.set_storage_type(PST_ADDITIVE);
#endif
}

///////////////////////////////////////////////////////////////////////////////
// set constraint values (instationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
void DomainDiscretization<TDomain, TDoFDistribution, TAlgebra>::
adjust_solution(vector_type& u, number time, const dof_distribution_type& dd)
{
	if(!update_constraints())
		UG_THROW_FATAL("Cannot update constraints.");

//	dirichlet
	for(size_t i = 0; i < m_vvConstraints[CT_DIRICHLET].size(); ++i)
	{
		if(!m_vvConstraints[CT_DIRICHLET][i]->adjust_solution(u, dd, time))
			UG_THROW_FATAL("Cannot adjust solution.");
	}

//	constraints
	for(size_t i = 0; i < m_vvConstraints[CT_CONSTRAINTS].size(); ++i)
	{
		if(!m_vvConstraints[CT_CONSTRAINTS][i]->adjust_solution(u, dd, time))
			UG_THROW_FATAL("Cannot adjust solution.");
	}
}

///////////////////////////////////////////////////////////////////////////////
// Finish Timestep (instationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
void DomainDiscretization<TDomain, TDoFDistribution, TAlgebra>::
finish_timestep(const VectorTimeSeries<vector_type>& vSol,
                const dof_distribution_type& dd)
{
//	update the elem discs
	if(!update_disc_items())
		UG_THROW_FATAL("Cannot update disc items.");

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	const ISubsetHandler& sh = dd.get_function_pattern().subset_handler();
	if(!CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, sh))
		UG_THROW_FATAL("ERROR in 'DomainDiscretization':"
				" Can not Subset Groups and union.\n");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.num_subsets(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = subset_dim(unionSubsets, i);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(m_bForceRegGrid) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	success flag
		bool bSuc = true;

	//	assemble on suitable elements
		switch(dim)
		{
		case 1:
			bSuc &= FinishTimestep<Edge,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol, m_pSelector);
			break;
		case 2:
			bSuc &= FinishTimestep<Triangle,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol, m_pSelector);
			bSuc &= FinishTimestep<Quadrilateral,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol, m_pSelector);
			break;
		case 3:
			bSuc &= FinishTimestep<Tetrahedron,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol, m_pSelector);
			bSuc &= FinishTimestep<Pyramid,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol, m_pSelector);
			bSuc &= FinishTimestep<Prism,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol, m_pSelector);
			bSuc &= FinishTimestep<Hexahedron,TDoFDistribution,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol, m_pSelector);
			break;
		default:
			UG_THROW_FATAL("ERROR in 'DomainDiscretization::finish_timestep (instationary)':"
						"Dimension "<<dim<<" (subset="<<si<<") not supported.\n");
		}

	//	check success
		if(!bSuc)
			UG_THROW_FATAL("ERROR in 'DomainDiscretization::finish_timestep (instationary)':"
							" Assembling of elements of Dimension " << dim << " in "
							" subset "<<si<< " failed.\n");
	}

}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
int DomainDiscretization<TDomain, TDoFDistribution, TAlgebra>::
subset_dim(SubsetGroup& subsetGrp, int subsetInd)
{
	#ifdef UG_PARALLEL
		return subsetGrp.dim(subsetInd, &m_procCom);
	#else
		return subsetGrp.dim(subsetInd);
	#endif
}

}

#endif /*__H__UG__LIB_DISC__SPATIAL_DISC__DOMAIN_DISC_IMPL__*/
