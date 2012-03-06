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
#include "lib_disc/spatial_disc/elem_disc/elem_disc_assemble_util.h"
#ifdef UG_PARALLEL
#include "lib_disc/parallelization/parallelization_util.h"
#endif

namespace ug{

template <typename TDomain, typename TAlgebra>
bool DomainDiscretization<TDomain, TAlgebra>::update_elem_discs()
{
//	check Approximation space
	if(!m_spApproxSpace.is_valid())
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
		m_vDomainElemDisc[i]->set_approximation_space(m_spApproxSpace);
		m_vElemDisc.push_back(m_vDomainElemDisc[i]);
	}

//	ok
	return true;
}

template <typename TDomain, typename TAlgebra>
bool DomainDiscretization<TDomain, TAlgebra>::update_constraints()
{
//	check Approximation space
	if(!m_spApproxSpace.is_valid())
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
			m_vvConstraints[type][i]->set_approximation_space(m_spApproxSpace);
		}
	}

//	ok
	return true;
}

template <typename TDomain, typename TAlgebra>
bool DomainDiscretization<TDomain, TAlgebra>::update_disc_items()
{
//	return flag
	bool bRet = true;

	bRet &= update_elem_discs();
	bRet &= update_constraints();

//	ok
	return bRet;
}


///////////////////////////////////////////////////////////////////////////////
// Mass Matrix
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra>
template <typename TDD>
void DomainDiscretization<TDomain, TAlgebra>::
assemble_mass_matrix(matrix_type& M, const vector_type& u,
                     ConstSmartPtr<TDD> dd)
{
//	update the elem discs
	if(!update_disc_items())
		UG_THROW_FATAL("DomainDisc: Cannot update disc items.")

//	reset matrix to zero and resize
	const size_t numIndex = dd->num_indices();
	M.resize(0,0);
	M.resize(numIndex, numIndex);
	M.set(0.0);

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	if(!CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, dd->subset_handler()))
		UG_THROW_FATAL("ERROR in 'DomainDiscretization':"
					   " Can not Subset Groups and union.\n");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.num_subsets(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = DimensionOfSubset(*dd->subset_handler(), si);

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
			bSuc &= AssembleMassMatrix<Edge,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u, m_pBoolMarker);
			break;
		case 2:
			bSuc &= AssembleMassMatrix<Triangle,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u, m_pBoolMarker);
			bSuc &= AssembleMassMatrix<Quadrilateral,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u, m_pBoolMarker);
			break;
		case 3:
			bSuc &= AssembleMassMatrix<Tetrahedron,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u, m_pBoolMarker);
			bSuc &= AssembleMassMatrix<Pyramid,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u, m_pBoolMarker);
			bSuc &= AssembleMassMatrix<Prism,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u, m_pBoolMarker);
			bSuc &= AssembleMassMatrix<Hexahedron,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u, m_pBoolMarker);
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
		if(type == CT_CONSTRAINTS && !m_bConstraintsEnabled) continue;
		for(size_t i = 0; i < m_vvConstraints[type].size(); ++i)
		{
			if(!m_vvConstraints[type][i]->adjust_jacobian(M, u, dd))
				UG_THROW_FATAL("DomainDiscretization: Cannot adjust jacobian.");
		}
	}

//	Remember parallel storage type
#ifdef UG_PARALLEL
	M.set_storage_type(PST_ADDITIVE);
	TDD* dist = const_cast<TDD*>(&dd);
	CopyLayoutsAndCommunicatorIntoMatrix(M, *dist);
#endif
}

///////////////////////////////////////////////////////////////////////////////
// Stiffness Matrix
///////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
template <typename TDD>
void DomainDiscretization<TDomain, TAlgebra>::
assemble_stiffness_matrix(matrix_type& A, const vector_type& u,
                          ConstSmartPtr<TDD> dd)
{
//	update the elem discs
	if(!update_disc_items())
		UG_THROW_FATAL("Cannot update disc items.");

//	reset matrix to zero and resize
	const size_t numIndex = dd->num_indices();
	A.resize(0,0);
	A.resize(numIndex, numIndex);
	A.set(0.0);

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	if(!CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, dd->subset_handler()))
		UG_THROW_FATAL("ERROR in 'DomainDiscretization':"
						" Can not Subset Groups and union.\n");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.num_subsets(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = DimensionOfSubset(*dd->subset_handler(), si);

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
			bSuc &= AssembleStiffnessMatrix<Edge,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u, m_pBoolMarker);
			break;
		case 2:
			bSuc &= AssembleStiffnessMatrix<Triangle,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u, m_pBoolMarker);
			bSuc &= AssembleStiffnessMatrix<Quadrilateral,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u, m_pBoolMarker);
			break;
		case 3:
			bSuc &= AssembleStiffnessMatrix<Tetrahedron,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u, m_pBoolMarker);
			bSuc &= AssembleStiffnessMatrix<Pyramid,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u, m_pBoolMarker);
			bSuc &= AssembleStiffnessMatrix<Prism,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u, m_pBoolMarker);
			bSuc &= AssembleStiffnessMatrix<Hexahedron,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u, m_pBoolMarker);
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
		if(type == CT_CONSTRAINTS && !m_bConstraintsEnabled) continue;
		for(size_t i = 0; i < m_vvConstraints[type].size(); ++i)
		{
			if(!m_vvConstraints[type][i]->adjust_jacobian(A, u, dd))
				UG_THROW_FATAL("Cannot adjust jacobian.");
		}
	}

//	Remember parallel storage type
#ifdef UG_PARALLEL
	A.set_storage_type(PST_ADDITIVE);
	TDD* dist = const_cast<TDD*>(&dd);
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
template <typename TDomain, typename TAlgebra>
template <typename TDD>
void DomainDiscretization<TDomain, TAlgebra>::
assemble_jacobian(matrix_type& J,
                  const vector_type& u,
                  ConstSmartPtr<TDD> dd)
{
//	update the elem discs
	if(!update_disc_items())
		UG_THROW_FATAL("Cannot update disc items.");

//	reset matrix to zero and resize
	const size_t numIndex = dd->num_indices();
	J.resize(0,0);
	J.resize(numIndex, numIndex);
	J.set(0.0);

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	try{
		CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, dd->subset_handler());
	}UG_CATCH_THROW("'DomainDiscretization': Can not create Subset Groups and Union.");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.num_subsets(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = DimensionOfSubset(*dd->subset_handler(), si);

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
			bSuc &= AssembleJacobian<Edge,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, u, m_pBoolMarker);
			break;
		case 2:
			bSuc &= AssembleJacobian<Triangle,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, u, m_pBoolMarker);
			bSuc &= AssembleJacobian<Quadrilateral,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, u, m_pBoolMarker);
			break;
		case 3:
			bSuc &= AssembleJacobian<Tetrahedron,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, u, m_pBoolMarker);
			bSuc &= AssembleJacobian<Pyramid,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, u, m_pBoolMarker);
			bSuc &= AssembleJacobian<Prism,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, u, m_pBoolMarker);
			bSuc &= AssembleJacobian<Hexahedron,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, u, m_pBoolMarker);
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
	try{
	for(size_t type = 0; type < NUM_CONSTRAINT_TYPES; ++type){
		if(type == CT_CONSTRAINTS && !m_bConstraintsEnabled) continue;
		for(size_t i = 0; i < m_vvConstraints[type].size(); ++i)
			m_vvConstraints[type][i]->adjust_jacobian(J, u, dd->grid_level());
	}
	}UG_CATCH_THROW("DomainDiscretization::assemble_jacobian:"
					" Cannot execute post process.");

//	Remember parallel storage type
#ifdef UG_PARALLEL
	J.set_storage_type(PST_ADDITIVE);
	TDD* pDD = const_cast<TDD*>(dd.get_impl());
	CopyLayoutsAndCommunicatorIntoMatrix(J, *pDD);
#endif
}


///////////////////////////////////////////////////////////////////////////////
// Defect (stationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra>
template <typename TDD>
void DomainDiscretization<TDomain, TAlgebra>::
assemble_defect(vector_type& d,
                const vector_type& u,
                ConstSmartPtr<TDD> dd)
{
//	update the elem discs
	if(!update_disc_items())
		UG_THROW_FATAL("Cannot update disc items.");

//	reset matrix to zero and resize
	const size_t numIndex = dd->num_indices();
	d.resize(numIndex);
	d.set(0.0);

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	try{
		CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, dd->subset_handler());
	}UG_CATCH_THROW("'DomainDiscretization': Can not create Subset Groups and Union.");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.num_subsets(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = DimensionOfSubset(*dd->subset_handler(), si);

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
			bSuc &= AssembleDefect<Edge,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, u, m_pBoolMarker);
			break;
		case 2:
			bSuc &= AssembleDefect<Triangle,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, u, m_pBoolMarker);
			bSuc &= AssembleDefect<Quadrilateral,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, u, m_pBoolMarker);
			break;
		case 3:
			bSuc &= AssembleDefect<Tetrahedron,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, u, m_pBoolMarker);
			bSuc &= AssembleDefect<Pyramid,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, u, m_pBoolMarker);
			bSuc &= AssembleDefect<Prism,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, u, m_pBoolMarker);
			bSuc &= AssembleDefect<Hexahedron,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, u, m_pBoolMarker);
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
	try{
	for(size_t type = 0; type < NUM_CONSTRAINT_TYPES; ++type){
		if(type == CT_CONSTRAINTS && !m_bConstraintsEnabled) continue;
		for(size_t i = 0; i < m_vvConstraints[type].size(); ++i)
			m_vvConstraints[type][i]->adjust_defect(d, u, dd->grid_level());
	}
	} UG_CATCH_THROW("Cannot adjust defect.");

//	Remember parallel storage type
#ifdef UG_PARALLEL
	d.set_storage_type(PST_ADDITIVE);
#endif
}

///////////////////////////////////////////////////////////////////////////////
// Matrix and RHS (stationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra>
template <typename TDD>
void DomainDiscretization<TDomain, TAlgebra>::
assemble_linear(matrix_type& mat, vector_type& rhs,
                ConstSmartPtr<TDD> dd)
{
//	update the elem discs
	if(!update_disc_items())
		UG_THROW_FATAL("Cannot update disc items.");

//	reset matrix to zero and resize
	const size_t numIndex = dd->num_indices();
	mat.resize(0,0);
	mat.resize(numIndex, numIndex);
	mat.set(0.0);

	rhs.resize(numIndex);
	rhs.set(0.0);

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	try{
		CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, dd->subset_handler());
	}UG_CATCH_THROW("'DomainDiscretization': Can not create Subset Groups and Union.");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.num_subsets(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = DimensionOfSubset(*dd->subset_handler(), si);

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
			bSuc &= AssembleLinear<Edge,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, m_pBoolMarker);
			break;
		case 2:
			bSuc &= AssembleLinear<Triangle,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, m_pBoolMarker);
			bSuc &= AssembleLinear<Quadrilateral,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, m_pBoolMarker);
			break;
		case 3:
			bSuc &= AssembleLinear<Tetrahedron,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, m_pBoolMarker);
			bSuc &= AssembleLinear<Pyramid,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, m_pBoolMarker);
			bSuc &= AssembleLinear<Prism,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, m_pBoolMarker);
			bSuc &= AssembleLinear<Hexahedron,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, m_pBoolMarker);
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
	try{
	for(size_t type = 0; type < NUM_CONSTRAINT_TYPES; ++type){
		if(type == CT_CONSTRAINTS && !m_bConstraintsEnabled) continue;
		for(size_t i = 0; i < m_vvConstraints[type].size(); ++i)
			m_vvConstraints[type][i]->adjust_linear(mat, rhs, dd->grid_level());
	}
	}UG_CATCH_THROW("DomainDiscretization::assemble_linear: Cannot post process.");

//	Remember parallel storage type
#ifdef UG_PARALLEL
	mat.set_storage_type(PST_ADDITIVE);
	TDD* pDD = const_cast<TDD*>(dd.get_impl());
	CopyLayoutsAndCommunicatorIntoMatrix(mat, *pDD);
	rhs.set_storage_type(PST_ADDITIVE);
#endif
}

///////////////////////////////////////////////////////////////////////////////
// RHS (stationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra>
template <typename TDD>
void DomainDiscretization<TDomain, TAlgebra>::
assemble_rhs(vector_type& rhs,
			const vector_type& u,
			ConstSmartPtr<TDD> dd)
{
//	update the elem discs
	if(!update_disc_items())
		UG_THROW_FATAL("Cannot update disc items.");

//	reset matrix to zero and resize
	const size_t numIndex = dd->num_indices();
	rhs.resize(numIndex);
	rhs.set(0.0);

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	if(!CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, dd->subset_handler()))
		UG_THROW_FATAL("ERROR in 'DomainDiscretization':"
						" Can not Subset Groups and union.\n");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.num_subsets(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = DimensionOfSubset(*dd->subset_handler(), si);

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
			bSuc &= AssembleRhs<Edge,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u, m_pBoolMarker);
			break;
		case 2:
			bSuc &= AssembleRhs<Triangle,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u, m_pBoolMarker);
			bSuc &= AssembleRhs<Quadrilateral,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u, m_pBoolMarker);
			break;
		case 3:
			bSuc &= AssembleRhs<Tetrahedron,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u, m_pBoolMarker);
			bSuc &= AssembleRhs<Pyramid,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u, m_pBoolMarker);
			bSuc &= AssembleRhs<Prism,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u, m_pBoolMarker);
			bSuc &= AssembleRhs<Hexahedron,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u, m_pBoolMarker);
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
		if(type == CT_CONSTRAINTS && !m_bConstraintsEnabled) continue;
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
template <typename TDomain, typename TAlgebra>
template <typename TDD>
void DomainDiscretization<TDomain, TAlgebra>::
adjust_solution(vector_type& u, ConstSmartPtr<TDD> dd)
{
	if(!update_constraints())
		UG_THROW_FATAL("Cannot update constraints.");

	try{
//	post process dirichlet
	for(size_t i = 0; i < m_vvConstraints[CT_DIRICHLET].size(); ++i)
		m_vvConstraints[CT_DIRICHLET][i]->adjust_solution(u, dd->grid_level());

//	post process constraints
	for(size_t i = 0; i < m_vvConstraints[CT_CONSTRAINTS].size(); ++i)
		m_vvConstraints[CT_CONSTRAINTS][i]->adjust_solution(u, dd->grid_level());

	} UG_CATCH_THROW("Cannot adjust solution.");
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//  Time Dependent (instationary)
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Prepare Timestep (instationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra>
template <typename TDD>
void DomainDiscretization<TDomain, TAlgebra>::
prepare_timestep(const VectorTimeSeries<vector_type>& vSol,
                ConstSmartPtr<TDD> dd)
{
//	update the elem discs
	if(!update_disc_items())
		UG_THROW_FATAL("Cannot update disc items.");

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	try{
		CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, dd->subset_handler());
	}UG_CATCH_THROW("'DomainDiscretization': Can not create Subset Groups and Union.");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.num_subsets(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = DimensionOfSubset(*dd->subset_handler(), si);

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
			bSuc &= PrepareTimestep<Edge,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol, m_pBoolMarker);
			break;
		case 2:
			bSuc &= PrepareTimestep<Triangle,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol, m_pBoolMarker);
			bSuc &= PrepareTimestep<Quadrilateral,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol, m_pBoolMarker);
			break;
		case 3:
			bSuc &= PrepareTimestep<Tetrahedron,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol, m_pBoolMarker);
			bSuc &= PrepareTimestep<Pyramid,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol, m_pBoolMarker);
			bSuc &= PrepareTimestep<Prism,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol, m_pBoolMarker);
			bSuc &= PrepareTimestep<Hexahedron,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol, m_pBoolMarker);
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
template <typename TDomain, typename TAlgebra>
template <typename TDD>
void DomainDiscretization<TDomain, TAlgebra>::
assemble_jacobian(matrix_type& J,
                  const VectorTimeSeries<vector_type>& vSol,
                  const number s_a0,
                  ConstSmartPtr<TDD> dd)
{
//	update the elem discs
	if(!update_disc_items())
		UG_THROW_FATAL("Cannot update disc items.");

//	reset matrix to zero and resize
	const size_t numIndex = dd->num_indices();
	J.resize(0,0);
	J.resize(numIndex, numIndex);
	J.set(0.0);

//	get current time
	const number time = vSol.time(0);

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	try{
		CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, dd->subset_handler());
	}UG_CATCH_THROW("'DomainDiscretization': Can not create Subset Groups and Union.");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.num_subsets(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = DimensionOfSubset(*dd->subset_handler(), si);

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
			bSuc &= AssembleJacobian<Edge,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, vSol, s_a0, m_pBoolMarker);
			break;
		case 2:
			bSuc &= AssembleJacobian<Triangle,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, vSol, s_a0, m_pBoolMarker);
			bSuc &= AssembleJacobian<Quadrilateral,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, vSol, s_a0, m_pBoolMarker);
			break;
		case 3:
			bSuc &= AssembleJacobian<Tetrahedron,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, vSol, s_a0, m_pBoolMarker);
			bSuc &= AssembleJacobian<Pyramid,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, vSol, s_a0, m_pBoolMarker);
			bSuc &= AssembleJacobian<Prism,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, vSol, s_a0, m_pBoolMarker);
			bSuc &= AssembleJacobian<Hexahedron,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, vSol, s_a0, m_pBoolMarker);
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
	try{
	for(size_t type = 0; type < NUM_CONSTRAINT_TYPES; ++type){
		if(type == CT_CONSTRAINTS && !m_bConstraintsEnabled) continue;
		for(size_t i = 0; i < m_vvConstraints[type].size(); ++i)
			m_vvConstraints[type][i]->adjust_jacobian(J, vSol.solution(0), dd->grid_level(), time);
	}
	}UG_CATCH_THROW("Cannot adjust jacobian.");

//	Remember parallel storage type
#ifdef UG_PARALLEL
	J.set_storage_type(PST_ADDITIVE);
	TDD* pDD = const_cast<TDD*>(dd.get_impl());
	CopyLayoutsAndCommunicatorIntoMatrix(J, *pDD);
#endif
}

///////////////////////////////////////////////////////////////////////////////
// Defect (instationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra>
template <typename TDD>
void DomainDiscretization<TDomain, TAlgebra>::
assemble_defect(vector_type& d,
                const VectorTimeSeries<vector_type>& vSol,
                const std::vector<number>& vScaleMass,
                const std::vector<number>& vScaleStiff,
                ConstSmartPtr<TDD> dd)
{
//	update the elem discs
	if(!update_disc_items())
		UG_THROW_FATAL("Cannot update disc items.");

//	reset matrix to zero and resize
	const size_t numIndex = dd->num_indices();
	d.resize(numIndex);
	d.set(0.0);

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	try{
		CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, dd->subset_handler());
	}UG_CATCH_THROW("'DomainDiscretization': Can not create Subset Groups and Union.");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.num_subsets(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = DimensionOfSubset(*dd->subset_handler(), si);

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
			bSuc &= AssembleDefect<Edge,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, vSol, vScaleMass, vScaleStiff, m_pBoolMarker);
			break;
		case 2:
			bSuc &= AssembleDefect<Triangle,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, vSol, vScaleMass, vScaleStiff, m_pBoolMarker);
			bSuc &= AssembleDefect<Quadrilateral,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, vSol, vScaleMass, vScaleStiff, m_pBoolMarker);
			break;
		case 3:
			bSuc &= AssembleDefect<Tetrahedron,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, vSol, vScaleMass, vScaleStiff, m_pBoolMarker);
			bSuc &= AssembleDefect<Pyramid,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, vSol, vScaleMass, vScaleStiff, m_pBoolMarker);
			bSuc &= AssembleDefect<Prism,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, vSol, vScaleMass, vScaleStiff, m_pBoolMarker);
			bSuc &= AssembleDefect<Hexahedron,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, vSol, vScaleMass, vScaleStiff, m_pBoolMarker);
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
	try{
	for(size_t type = 0; type < NUM_CONSTRAINT_TYPES; ++type){
		if(type == CT_CONSTRAINTS && !m_bConstraintsEnabled) continue;
		for(size_t i = 0; i < m_vvConstraints[type].size(); ++i)
			m_vvConstraints[type][i]->adjust_defect(d, vSol.solution(0), dd->grid_level(), vSol.time(0));
	}
	} UG_CATCH_THROW("Cannot adjust defect.");

//	Remember parallel storage type
#ifdef UG_PARALLEL
	d.set_storage_type(PST_ADDITIVE);
#endif
}

///////////////////////////////////////////////////////////////////////////////
// Matrix and RHS (instationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra>
template <typename TDD>
void DomainDiscretization<TDomain, TAlgebra>::
assemble_linear(matrix_type& mat, vector_type& rhs,
                const VectorTimeSeries<vector_type>& vSol,
                const std::vector<number>& vScaleMass,
                const std::vector<number>& vScaleStiff,
                ConstSmartPtr<TDD> dd)
{
//	update the elem discs
	if(!update_disc_items())
		UG_THROW_FATAL("Cannot update disc items.");

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	try{
		CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, dd->subset_handler());
	}UG_CATCH_THROW("'DomainDiscretization': Can not create Subset Groups and Union.");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.num_subsets(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = DimensionOfSubset(*dd->subset_handler(), si);

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
			bSuc &= AssembleLinear<Edge,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff, m_pBoolMarker);
			break;
		case 2:
			bSuc &= AssembleLinear<Triangle,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff, m_pBoolMarker);
			bSuc &= AssembleLinear<Quadrilateral,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff, m_pBoolMarker);
			break;
		case 3:
			bSuc &= AssembleLinear<Tetrahedron,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff, m_pBoolMarker);
			bSuc &= AssembleLinear<Pyramid,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff, m_pBoolMarker);
			bSuc &= AssembleLinear<Prism,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff, m_pBoolMarker);
			bSuc &= AssembleLinear<Hexahedron,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff, m_pBoolMarker);
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
	try{
	for(size_t type = 0; type < NUM_CONSTRAINT_TYPES; ++type){
		if(type == CT_CONSTRAINTS && !m_bConstraintsEnabled) continue;
		for(size_t i = 0; i < m_vvConstraints[type].size(); ++i)
			m_vvConstraints[type][i]->adjust_linear(mat, rhs, dd->grid_level(), vSol.time(0));
	}
	} UG_CATCH_THROW("Cannot adjust linear.");

//	Remember parallel storage type
#ifdef UG_PARALLEL
	mat.set_storage_type(PST_ADDITIVE);
	TDD* pDD = const_cast<TDD*>(dd.get_impl());
	CopyLayoutsAndCommunicatorIntoMatrix(mat, *pDD);

	rhs.set_storage_type(PST_ADDITIVE);
#endif
}

///////////////////////////////////////////////////////////////////////////////
// set constraint values (instationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra>
template <typename TDD>
void DomainDiscretization<TDomain, TAlgebra>::
adjust_solution(vector_type& u, number time, ConstSmartPtr<TDD> dd)
{
	if(!update_constraints())
		UG_THROW_FATAL("Cannot update constraints.");

	try{

//	dirichlet
	for(size_t i = 0; i < m_vvConstraints[CT_DIRICHLET].size(); ++i)
		m_vvConstraints[CT_DIRICHLET][i]->adjust_solution(u, dd->grid_level(), time);

//	constraints
	for(size_t i = 0; i < m_vvConstraints[CT_CONSTRAINTS].size(); ++i)
		m_vvConstraints[CT_CONSTRAINTS][i]->adjust_solution(u, dd->grid_level(), time);

	} UG_CATCH_THROW(" Cannot adjust solution.");
}

///////////////////////////////////////////////////////////////////////////////
// Finish Timestep (instationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra>
template <typename TDD>
void DomainDiscretization<TDomain, TAlgebra>::
finish_timestep(const VectorTimeSeries<vector_type>& vSol,
                ConstSmartPtr<TDD> dd)
{
//	update the elem discs
	if(!update_disc_items())
		UG_THROW_FATAL("Cannot update disc items.");

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	try{
		CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, dd->subset_handler());
	}UG_CATCH_THROW("'DomainDiscretization': Can not create Subset Groups and Union.");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.num_subsets(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = DimensionOfSubset(*dd->subset_handler(), si);

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
			bSuc &= FinishTimestep<Edge, TDD, TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol, m_pBoolMarker);
			break;
		case 2:
			bSuc &= FinishTimestep<Triangle,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol, m_pBoolMarker);
			bSuc &= FinishTimestep<Quadrilateral,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol, m_pBoolMarker);
			break;
		case 3:
			bSuc &= FinishTimestep<Tetrahedron,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol, m_pBoolMarker);
			bSuc &= FinishTimestep<Pyramid,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol, m_pBoolMarker);
			bSuc &= FinishTimestep<Prism,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol, m_pBoolMarker);
			bSuc &= FinishTimestep<Hexahedron,TDD,TAlgebra>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol, m_pBoolMarker);
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

} // end namespace ug

#endif /*__H__UG__LIB_DISC__SPATIAL_DISC__DOMAIN_DISC_IMPL__*/
