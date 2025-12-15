/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DOMAIN_DISC_IMPL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DOMAIN_DISC_IMPL__

#include "domain_disc.h"

#include "common/profiler/profiler.h"
//#include "lib_disc/common/groups_util.h"
//#include "lib_disc/function_spaces/error_indicator_util.h"
#include "lib_disc/spatial_disc/subset_assemble_util.h"
#include "lib_grid/algorithms/debug_util.h"

#ifdef UG_PARALLEL
#include "lib_disc/parallelization/parallelization_util.h"
#endif


namespace ug {

template <typename TElemDisc>
static void prep_assemble_loop(std::vector<TElemDisc*> vElemDisc)
{
	for(size_t i = 0; i < vElemDisc.size(); ++i)
	{
		vElemDisc[i]->prep_assemble_loop();
	}
}

template <typename TElemDisc>
static void post_assemble_loop(std::vector<TElemDisc*> vElemDisc)
{
	for(size_t i = 0; i < vElemDisc.size(); ++i)
	{
		vElemDisc[i]->post_assemble_loop();
	}
}


template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::update_elem_discs()
{
//	check Approximation space
	if(!m_spApproxSpace.valid())
		UG_THROW("DomainDiscretization: Before using the "
				"DomainDiscretization an ApproximationSpace must be set to it. "
				"Please use DomainDiscretization:set_approximation_space to "
				"set an appropriate Space.");

//	set approximation space and extract IElemDiscs
	m_vElemDisc.clear();
	for(size_t i = 0; i < m_vDomainElemDisc.size(); ++i)
	{
		m_vDomainElemDisc[i]->set_approximation_space(m_spApproxSpace);

		if(!(m_spAssTuner->elem_disc_type_enabled(m_vDomainElemDisc[i]->type()))) continue;
		m_vElemDisc.push_back(m_vDomainElemDisc[i].get());
	}
}

template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::update_elem_errors()
{
//	check Approximation space
	if(!m_spApproxSpace.valid())
		UG_THROW("DomainDiscretization: Before using the "
				"DomainDiscretization an ApproximationSpace must be set to it. "
				"Please use DomainDiscretization:set_approximation_space to "
				"set an appropriate Space.");

//	set approximation space and extract IElemDiscs
	m_vElemError.clear();
	for(size_t i = 0; i < m_vDomainElemError.size(); ++i)
	{
		m_vDomainElemError[i]->set_approximation_space(m_spApproxSpace);

		if(!(m_spAssTuner->elem_disc_type_enabled(m_vDomainElemError[i]->type()))) continue;
		m_vElemError.push_back(m_vDomainElemError[i].get());
	}

}


template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::update_constraints()
{
//	check Approximation space
	if(!m_spApproxSpace.valid())
		UG_THROW("DomainDiscretization: Before using the "
				"DomainDiscretization an ApproximationSpace must be set to it. "
				"Please use DomainDiscretization:set_approximation_space to "
				"set an appropriate Space.");


	for(size_t i = 0; i < m_vConstraint.size(); ++i)
		m_vConstraint[i]->set_approximation_space(m_spApproxSpace);
}

template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::update_disc_items()
{
	update_elem_discs();
	update_constraints();
}

template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::update_error_items()
{
	update_elem_errors();
	update_constraints();
}

///////////////////////////////////////////////////////////////////////////////
// Mass Matrix
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
assemble_mass_matrix(matrix_type& M, const vector_type& u,
                     ConstSmartPtr<DoFDistribution> dd)
{
	PROFILE_FUNC_GROUP("discretization");
//	update the elem discs
	update_disc_items();
	prep_assemble_loop(m_vElemDisc);

//	reset matrix to zero and resize
	m_spAssTuner->resize(dd, M);

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	try{
		CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, dd->subset_handler());
	}UG_CATCH_THROW("'DomainDiscretization': Can not create Subset Groups and Union.");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.size(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = DimensionOfSubset(*dd->subset_handler(), si);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(m_spAssTuner->regular_grid_forced()) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc<TDomain>*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	assemble on suitable elements
		try
		{
		switch(dim)
		{
		case 0:
			this->AssembleMassMatrix<RegularVertex>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u);
			break;
		case 1:
			this->AssembleMassMatrix<RegularEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->AssembleMassMatrix<ConstrainingEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u);
			break;
		case 2:
			this->AssembleMassMatrix<Triangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u);
			this->AssembleMassMatrix<Quadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->AssembleMassMatrix<ConstrainingTriangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u);
			this->AssembleMassMatrix<ConstrainingQuadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u);
			break;
		case 3:
			this->AssembleMassMatrix<Tetrahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u);
			this->AssembleMassMatrix<Pyramid>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u);
			this->AssembleMassMatrix<Prism>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u);
			this->AssembleMassMatrix<Hexahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u);
			this->AssembleMassMatrix<Octahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u);
			break;
		default:
			UG_THROW("DomainDiscretization::assemble_mass_matrix:"
							"Dimension "<<dim<<" (subset="<<si<<") not supported.");
		}
		}
		UG_CATCH_THROW("DomainDiscretization::assemble_mass_matrix:"
						" Assembling of elements of Dimension " << dim << " in "
						" subset "<<si<< " failed.");
	}

//	post process
	try{
	for(int type = 1; type < CT_ALL; type = type << 1){
		if(!(m_spAssTuner->constraint_type_enabled(type))) continue;
		for(size_t i = 0; i < m_vConstraint.size(); ++i)
			if(m_vConstraint[i]->type() & type)
			{
				m_vConstraint[i]->set_ass_tuner(m_spAssTuner);
				m_vConstraint[i]->adjust_jacobian(M, u, dd, type);
			}
	}
	post_assemble_loop(m_vElemDisc);
	}UG_CATCH_THROW("DomainDiscretization::assemble_mass_matrix:"
					" Cannot execute post process.");

//	Remember parallel storage type
#ifdef UG_PARALLEL
	M.set_storage_type(PST_ADDITIVE);
	M.set_layouts(dd->layouts());
#endif
}

/**
 * This function adds the contributions of all passed element discretizations
 * on one given subset to the global Mass matrix.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd				DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	M				Mass matrix
 * \param[in]		u				solution
 */
template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
template <typename TElem>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
AssembleMassMatrix( const std::vector<IElemDisc<domain_type>*>& vElemDisc,
					ConstSmartPtr<DoFDistribution> dd,
					int si, bool bNonRegularGrid,
					matrix_type& M,
					const vector_type& u)
{
	//	check if only some elements are selected
	if(m_spAssTuner->selected_elements_used())
	{
		std::vector<TElem*> vElem;
		m_spAssTuner->collect_selected_elements(vElem, dd, si);

		//	assembling is carried out only over those elements
		//	which are selected and in subset si
		gass_type::template AssembleMassMatrix<TElem>
			(vElemDisc, m_spApproxSpace->domain(), dd, vElem.begin(), vElem.end(), si,
			 bNonRegularGrid, M, u, m_spAssTuner);
	}
	else
	{
		//	general case: assembling over all elements in subset si
		gass_type::template AssembleMassMatrix<TElem>
			(vElemDisc, m_spApproxSpace->domain(), dd,
				dd->begin<TElem>(si), dd->end<TElem>(si), si,
					bNonRegularGrid, M, u, m_spAssTuner);
	}
}

///////////////////////////////////////////////////////////////////////////////
// Stiffness Matrix
///////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
assemble_stiffness_matrix(matrix_type& A, const vector_type& u,
                          ConstSmartPtr<DoFDistribution> dd)
{
	PROFILE_FUNC_GROUP("discretization");
//	update the elem discs
	update_disc_items();
	prep_assemble_loop(m_vElemDisc);

//	reset matrix to zero and resize
	m_spAssTuner->resize(dd, A);

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	try{
		CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, dd->subset_handler());
	}UG_CATCH_THROW("'DomainDiscretization': Can not create Subset Groups and Union.");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.size(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = DimensionOfSubset(*dd->subset_handler(), si);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(m_spAssTuner->regular_grid_forced()) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc<TDomain>*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	assemble on suitable elements
		try
		{
		switch(dim)
		{
		case 0:
			this->AssembleStiffnessMatrix<RegularVertex>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u);
			break;
		case 1:
			this->AssembleStiffnessMatrix<RegularEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->AssembleStiffnessMatrix<ConstrainingEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u);
			break;
		case 2:
			this->AssembleStiffnessMatrix<Triangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u);
			this->AssembleStiffnessMatrix<Quadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->AssembleStiffnessMatrix<ConstrainingTriangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u);
			this->AssembleStiffnessMatrix<ConstrainingQuadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u);
			break;
		case 3:
			this->AssembleStiffnessMatrix<Tetrahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u);
			this->AssembleStiffnessMatrix<Pyramid>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u);
			this->AssembleStiffnessMatrix<Prism>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u);
			this->AssembleStiffnessMatrix<Hexahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u);
			this->AssembleStiffnessMatrix<Octahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u);
			break;
		default:
			UG_THROW("DomainDiscretization::assemble_stiffness_matrix:"
							"Dimension "<<dim<<" (subset="<<si<<") not supported.");
		}
		}
		UG_CATCH_THROW("DomainDiscretization::assemble_stiffness_matrix:"
					" Assembling of elements of Dimension " << dim << " in "
					" subset "<<si<< " failed.");
	}

//	post process
	try{
	for(int type = 1; type < CT_ALL; type = type << 1){
		if(!(m_spAssTuner->constraint_type_enabled(type))) continue;
		for(size_t i = 0; i < m_vConstraint.size(); ++i)
			if(m_vConstraint[i]->type() & type)
			{
				m_vConstraint[i]->set_ass_tuner(m_spAssTuner);
				m_vConstraint[i]->adjust_jacobian(A, u, dd, type);
			}
	}
	post_assemble_loop(m_vElemDisc);
	}UG_CATCH_THROW("DomainDiscretization::assemble_stiffness_matrix:"
					" Cannot execute post process.");

//	Remember parallel storage type
#ifdef UG_PARALLEL
	A.set_storage_type(PST_ADDITIVE);
	A.set_layouts(dd->layouts());
#endif
}

/**
 * This function adds the contributions of all passed element discretizations
 * on one given subset to the global Stiffness matrix.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd				DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	A				Stiffness matrix
 * \param[in]		u				solution
 */
template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
template <typename TElem>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
AssembleStiffnessMatrix(	const std::vector<IElemDisc<domain_type>*>& vElemDisc,
							ConstSmartPtr<DoFDistribution> dd,
							int si, bool bNonRegularGrid,
							matrix_type& A,
							const vector_type& u)
{
	//	check if only some elements are selected
	if(m_spAssTuner->selected_elements_used())
	{
		std::vector<TElem*> vElem;
		m_spAssTuner->collect_selected_elements(vElem, dd, si);

		//	assembling is carried out only over those elements
		//	which are selected and in subset si
		gass_type::template AssembleStiffnessMatrix<TElem>
			(vElemDisc, m_spApproxSpace->domain(), dd, vElem.begin(), vElem.end(), si,
			 bNonRegularGrid, A, u, m_spAssTuner);
	}
	else
	{
		//	general case: assembling over all elements in subset si
		gass_type::template AssembleStiffnessMatrix<TElem>
			(vElemDisc, m_spApproxSpace->domain(), dd,
				dd->begin<TElem>(si), dd->end<TElem>(si), si,
					bNonRegularGrid, A, u, m_spAssTuner);
	}
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//  Time Independent (stationary)
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// Jacobian (stationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
assemble_jacobian(matrix_type& J,
                  const vector_type& u,
                  ConstSmartPtr<DoFDistribution> dd)
{
	PROFILE_FUNC_GROUP("discretization");
//	update the elem discs
	update_disc_items();
	prep_assemble_loop(m_vElemDisc);

//	reset matrix to zero and resize
	m_spAssTuner->resize(dd, J);

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	pre process -  modifies the solution, used for computing the defect
	const vector_type* pModifyU = &u;
	SmartPtr<vector_type> pModifyMemory;
	if( m_spAssTuner->modify_solution_enabled() ){
		pModifyMemory = u.clone();
		pModifyU = pModifyMemory.get();
		try{
		for(int type = 1; type < CT_ALL; type = type << 1){
			if(!(m_spAssTuner->constraint_type_enabled(type))) continue;
			for(size_t i = 0; i < m_vConstraint.size(); ++i)
				if(m_vConstraint[i]->type() & type)
					m_vConstraint[i]->modify_solution(*pModifyMemory, u, dd, type);
		}
		} UG_CATCH_THROW("Cannot modify solution.");
	}



//	create list of all subsets
	try{
		CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, dd->subset_handler());
	}UG_CATCH_THROW("'DomainDiscretization': Can not create Subset Groups and Union.");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.size(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = DimensionOfSubset(*dd->subset_handler(), si);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(m_spAssTuner->regular_grid_forced()) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc<TDomain>*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	assemble on suitable elements
		try
		{
		switch(dim)
		{
		case 0:
			this->AssembleJacobian<RegularVertex>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, *pModifyU);
			break;
		case 1:
			this->AssembleJacobian<RegularEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, *pModifyU);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->AssembleJacobian<ConstrainingEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, *pModifyU);
			break;
		case 2:
			this->AssembleJacobian<Triangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, *pModifyU);
			this->AssembleJacobian<Quadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, *pModifyU);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->AssembleJacobian<ConstrainingTriangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, *pModifyU);
			this->AssembleJacobian<ConstrainingQuadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, *pModifyU);
			break;
		case 3:
			this->AssembleJacobian<Tetrahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, *pModifyU);
			this->AssembleJacobian<Pyramid>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, *pModifyU);
			this->AssembleJacobian<Prism>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, *pModifyU);
			this->AssembleJacobian<Hexahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, *pModifyU);
			this->AssembleJacobian<Octahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, *pModifyU);
			break;
		default:
			UG_THROW("DomainDiscretization::assemble_jacobian (stationary):"
							"Dimension "<<dim<<"(subset="<<si<<") not supported");
		}
		}
		UG_CATCH_THROW("DomainDiscretization::assemble_jacobian (stationary):"
						" Assembling of elements of Dimension " << dim << " in "
						" subset "<<si<< " failed.");
	}

//	post process
	try{
	for(int type = 1; type < CT_ALL; type = type << 1){
		if(!(m_spAssTuner->constraint_type_enabled(type))) continue;
		for(size_t i = 0; i < m_vConstraint.size(); ++i)
			if(m_vConstraint[i]->type() & type)
			{
				m_vConstraint[i]->set_ass_tuner(m_spAssTuner);
				m_vConstraint[i]->adjust_jacobian(J, *pModifyU, dd, type);
			}
	}
	post_assemble_loop(m_vElemDisc);
	}UG_CATCH_THROW("DomainDiscretization::assemble_jacobian:"
					" Cannot execute post process.");

//	Remember parallel storage type
#ifdef UG_PARALLEL
	J.set_storage_type(PST_ADDITIVE);
	J.set_layouts(dd->layouts());
#endif
}

/**
 * This function adds the contributions of all passed element discretizations
 * on one given subset to the global Jacobian in the stationary case.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	J				jacobian
 * \param[in]		u				solution
 */
template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
template <typename TElem>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
AssembleJacobian(	const std::vector<IElemDisc<domain_type>*>& vElemDisc,
					ConstSmartPtr<DoFDistribution> dd,
					int si, bool bNonRegularGrid,
					matrix_type& J,
					const vector_type& u)
{
	//	check if only some elements are selected
	if(m_spAssTuner->selected_elements_used())
	{
		std::vector<TElem*> vElem;
		m_spAssTuner->collect_selected_elements(vElem, dd, si);

		//	assembling is carried out only over those elements
		//	which are selected and in subset si
		gass_type::template AssembleJacobian<TElem>
			(vElemDisc, m_spApproxSpace->domain(), dd, vElem.begin(), vElem.end(), si,
			 bNonRegularGrid, J, u, m_spAssTuner);
	}
	else
	{
		//	general case: assembling over all elements in subset si
		gass_type::template AssembleJacobian<TElem>
			(vElemDisc, m_spApproxSpace->domain(), dd,
				dd->begin<TElem>(si), dd->end<TElem>(si), si,
					bNonRegularGrid, J, u, m_spAssTuner);
	}
}

///////////////////////////////////////////////////////////////////////////////
// Defect (stationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
assemble_defect(vector_type& d,
                const vector_type& u,
                ConstSmartPtr<DoFDistribution> dd)
{
	PROFILE_FUNC_GROUP("discretization");
//	update the elem discs
	update_disc_items();
	prep_assemble_loop(m_vElemDisc);

//	reset matrix to zero and resize
	m_spAssTuner->resize(dd, d);

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	pre process -  modifies the solution, used for computing the defect
	const vector_type* pModifyU = &u;
	SmartPtr<vector_type> pModifyMemory;
	if( m_spAssTuner->modify_solution_enabled() ){
		pModifyMemory = u.clone();
		pModifyU = pModifyMemory.get();
		try{
		for(int type = 1; type < CT_ALL; type = type << 1){
			if(!(m_spAssTuner->constraint_type_enabled(type))) continue;
			for(size_t i = 0; i < m_vConstraint.size(); ++i)
				if(m_vConstraint[i]->type() & type)
					m_vConstraint[i]->modify_solution(*pModifyMemory, u, dd, type);
		}
		} UG_CATCH_THROW("Cannot modify solution.");
	}

//	create list of all subsets
	try{
		CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, dd->subset_handler());
	}UG_CATCH_THROW("'DomainDiscretization': Can not create Subset Groups and Union.");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.size(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = DimensionOfSubset(*dd->subset_handler(), si);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(m_spAssTuner->regular_grid_forced()) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc<TDomain>*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	assemble on suitable elements
		try
		{
		switch(dim)
		{
		case 0:
			this->AssembleDefect<RegularVertex>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, *pModifyU);
			break;
		case 1:
			this->AssembleDefect<RegularEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, *pModifyU);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->AssembleDefect<ConstrainingEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, *pModifyU);
			break;
		case 2:
			this->AssembleDefect<Triangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, *pModifyU);
			this->AssembleDefect<Quadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, *pModifyU);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->AssembleDefect<ConstrainingTriangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, *pModifyU);
			this->AssembleDefect<ConstrainingQuadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, *pModifyU);
			break;
		case 3:
			this->AssembleDefect<Tetrahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, *pModifyU);
			this->AssembleDefect<Pyramid>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, *pModifyU);
			this->AssembleDefect<Prism>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, *pModifyU);
			this->AssembleDefect<Hexahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, *pModifyU);
			this->AssembleDefect<Octahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, *pModifyU);
			break;
		default:
			UG_THROW("DomainDiscretization::assemble_defect (stationary):"
							"Dimension "<<dim<<" (subset="<<si<<") not supported.");
		}
		}
		UG_CATCH_THROW("DomainDiscretization::assemble_defect (stationary):"
						" Assembling of elements of Dimension " << dim << " in "
						" subset "<<si<< " failed.");
	}

//	post process
	try{

	// Dirichlet first, since hanging nodes might be constrained by Dirichlet nodes
	if (m_spAssTuner->constraint_type_enabled(CT_DIRICHLET))
	{
		for (size_t i = 0; i < m_vConstraint.size(); ++i)
		{
			if (m_vConstraint[i]->type() & CT_DIRICHLET)
			{
				m_vConstraint[i]->set_ass_tuner(m_spAssTuner);
				m_vConstraint[i]->adjust_defect(d, *pModifyU, dd, CT_DIRICHLET);
			}
		}
	}

	for(int type = 1; type < CT_ALL; type = type << 1){
		if(!(m_spAssTuner->constraint_type_enabled(type))) continue;
		for(size_t i = 0; i < m_vConstraint.size(); ++i)
			if(m_vConstraint[i]->type() & type)
			{
				m_vConstraint[i]->set_ass_tuner(m_spAssTuner);
				m_vConstraint[i]->adjust_defect(d, *pModifyU, dd, type);
			}
	}
	post_assemble_loop(m_vElemDisc);
	} UG_CATCH_THROW("Cannot adjust defect.");


//	Remember parallel storage type
#ifdef UG_PARALLEL
	d.set_storage_type(PST_ADDITIVE);
#endif
}

/**
 * This function adds the contributions of all passed element discretizations
 * on one given subset to the global Defect in the stationary case.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd				DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	d				defect
 * \param[in]		u				solution
 */
template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
template <typename TElem>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
AssembleDefect( const std::vector<IElemDisc<domain_type>*>& vElemDisc,
				ConstSmartPtr<DoFDistribution> dd,
				int si, bool bNonRegularGrid,
				vector_type& d,
				const vector_type& u)
{
	//	check if only some elements are selected
	if(m_spAssTuner->selected_elements_used())
	{
		std::vector<TElem*> vElem;
		m_spAssTuner->collect_selected_elements(vElem, dd, si);

		//	assembling is carried out only over those elements
		//	which are selected and in subset si
		gass_type::template AssembleDefect<TElem>
			(vElemDisc, m_spApproxSpace->domain(), dd, vElem.begin(), vElem.end(), si,
			 bNonRegularGrid, d, u, m_spAssTuner);
	}
	else
	{
		//	general case: assembling over all elements in subset si
		gass_type::template AssembleDefect<TElem>
			(vElemDisc, m_spApproxSpace->domain(), dd,
				dd->begin<TElem>(si), dd->end<TElem>(si), si,
					bNonRegularGrid, d, u, m_spAssTuner);
	}
}

///////////////////////////////////////////////////////////////////////////////
// Matrix and RHS (stationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
assemble_linear(matrix_type& mat, vector_type& rhs,
                ConstSmartPtr<DoFDistribution> dd)
{
	PROFILE_FUNC_GROUP("discretization");
//	update the elem discs
	update_disc_items();
	prep_assemble_loop(m_vElemDisc);

//	reset matrix to zero and resize
	m_spAssTuner->resize(dd, mat);
	m_spAssTuner->resize(dd, rhs);

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	try{
		CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, dd->subset_handler());
	}UG_CATCH_THROW("DomainDiscretization: Can not create Subset Groups and Union.");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.size(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = DimensionOfSubset(*dd->subset_handler(), si);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(m_spAssTuner->regular_grid_forced()) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc<TDomain>*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	assemble on suitable elements
		try
		{
		switch(dim)
		{
		case 0:
			this->AssembleLinear<RegularVertex>
					(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs);
			break;
		case 1:
			this->AssembleLinear<RegularEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->AssembleLinear<ConstrainingEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs);
			break;
		case 2:
			this->AssembleLinear<Triangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs);
			this->AssembleLinear<Quadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->AssembleLinear<ConstrainingTriangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs);
			this->AssembleLinear<ConstrainingQuadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs);
			break;
		case 3:
			this->AssembleLinear<Tetrahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs);
			this->AssembleLinear<Pyramid>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs);
			this->AssembleLinear<Prism>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs);
			this->AssembleLinear<Hexahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs);
			this->AssembleLinear<Octahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs);
			break;
		default:
			UG_THROW("DomainDiscretization::assemble_linear (stationary):"
							"Dimension "<<dim<<" (subset="<<si<<") not supported.");
		}
		}
		UG_CATCH_THROW("DomainDiscretization::assemble_linear (stationary):"
						" Assembling of elements of Dimension " << dim << " in "
						" subset "<<si<< " failed.");
	}

//	post process
	try{
	for(int type = 1; type < CT_ALL; type = type << 1){
		if(!(m_spAssTuner->constraint_type_enabled(type))) continue;
		for(size_t i = 0; i < m_vConstraint.size(); ++i)
			if(m_vConstraint[i]->type() & type)
			{
				m_vConstraint[i]->set_ass_tuner(m_spAssTuner);
				m_vConstraint[i]->adjust_linear(mat, rhs, dd, type);
			}
	}
	post_assemble_loop(m_vElemDisc);
	}UG_CATCH_THROW("DomainDiscretization::assemble_linear: Cannot post process.");

//	Remember parallel storage type
#ifdef UG_PARALLEL
	mat.set_storage_type(PST_ADDITIVE);
	mat.set_layouts(dd->layouts());
	rhs.set_storage_type(PST_ADDITIVE);
#endif
}

/**
 * This function adds the contributions of all passed element discretizations
 * on one given subset to the global Matrix and the global Right-Hand Side
 * of the Linear problem in the stationary case.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd				DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	A				Matrix
 * \param[in,out]	rhs				Right-hand side
 */
template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
template <typename TElem>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
AssembleLinear( const std::vector<IElemDisc<domain_type>*>& vElemDisc,
				ConstSmartPtr<DoFDistribution> dd,
				int si, bool bNonRegularGrid,
				matrix_type& A,
				vector_type& rhs)
{
	//	check if only some elements are selected
	if(m_spAssTuner->selected_elements_used())
	{
		std::vector<TElem*> vElem;
		m_spAssTuner->collect_selected_elements(vElem, dd, si);

		//	assembling is carried out only over those elements
		//	which are selected and in subset si
		gass_type::template AssembleLinear<TElem>
			(vElemDisc, m_spApproxSpace->domain(), dd, vElem.begin(), vElem.end(), si,
			 bNonRegularGrid, A, rhs, m_spAssTuner);
	}
	else
	{
		//	general case: assembling over all elements in subset si
		gass_type::template AssembleLinear<TElem>
			(vElemDisc, m_spApproxSpace->domain(), dd,
				dd->begin<TElem>(si), dd->end<TElem>(si), si,
					bNonRegularGrid, A, rhs, m_spAssTuner);
	}
}

///////////////////////////////////////////////////////////////////////////////
// RHS (stationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
assemble_rhs(vector_type& rhs,
			const vector_type& u,
			ConstSmartPtr<DoFDistribution> dd)
{
	PROFILE_FUNC_GROUP("discretization");
//	update the elem discs
	update_disc_items();
	prep_assemble_loop(m_vElemDisc);

//	reset matrix to zero and resize
	m_spAssTuner->resize(dd, rhs);

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	try{
		CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, dd->subset_handler());
	}UG_CATCH_THROW("'DomainDiscretization': Can not create Subset Groups and Union.");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.size(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = DimensionOfSubset(*dd->subset_handler(), si);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(m_spAssTuner->regular_grid_forced()) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc<TDomain>*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	assemble on suitable elements
		try
		{
		switch(dim)
		{
		case 0:
			this->AssembleRhs<RegularVertex>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u);
			break;
		case 1:
			this->AssembleRhs<RegularEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->AssembleRhs<ConstrainingEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u);
			break;
		case 2:
			this->AssembleRhs<Triangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u);
			this->AssembleRhs<Quadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->AssembleRhs<ConstrainingTriangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u);
			this->AssembleRhs<ConstrainingQuadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u);
			break;
		case 3:
			this->AssembleRhs<Tetrahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u);
			this->AssembleRhs<Pyramid>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u);
			this->AssembleRhs<Prism>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u);
			this->AssembleRhs<Hexahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u);
			this->AssembleRhs<Octahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u);
			break;
		default:
			UG_THROW("DomainDiscretization::assemble_rhs (stationary):"
							"Dimension "<<dim<<" (subset="<<si<<") not supported.");
		}
		}
		UG_CATCH_THROW("DomainDiscretization::assemble_rhs (stationary):"
						" Assembling of elements of Dimension " << dim << " in "
						" subset "<<si<< " failed.");
	}

//	post process
	try{
	for(int type = 1; type < CT_ALL; type = type << 1){
		if(!(m_spAssTuner->constraint_type_enabled(type))) continue;
		for(size_t i = 0; i < m_vConstraint.size(); ++i)
			if(m_vConstraint[i]->type() & type)
			{
				m_vConstraint[i]->set_ass_tuner(m_spAssTuner);
				m_vConstraint[i]->adjust_rhs(rhs, u, dd, type);
			}
	}
	post_assemble_loop(m_vElemDisc);
	}UG_CATCH_THROW("DomainDiscretization::assemble_rhs:"
					" Cannot execute post process.");

//	Remember parallel storage type
#ifdef UG_PARALLEL
	rhs.set_storage_type(PST_ADDITIVE);
#endif
}

/**
 * This function adds the contributions of all passed element discretizations
 * on one given subset to the global Right-Hand Side.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd				DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	rhs				Right-hand side
 * \param[in]		u				solution
 */
template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
template <typename TElem>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
AssembleRhs(	const std::vector<IElemDisc<domain_type>*>& vElemDisc,
				ConstSmartPtr<DoFDistribution> dd,
				int si, bool bNonRegularGrid,
				vector_type& rhs,
				const vector_type& u)
{
	//	check if only some elements are selected
	if(m_spAssTuner->selected_elements_used())
	{
		std::vector<TElem*> vElem;
		m_spAssTuner->collect_selected_elements(vElem, dd, si);

		//	assembling is carried out only over those elements
		//	which are selected and in subset si
		gass_type::template AssembleRhs<TElem>
			(vElemDisc, m_spApproxSpace->domain(), dd, vElem.begin(), vElem.end(), si,
			 bNonRegularGrid, rhs, u, m_spAssTuner);
	}
	else
	{
		//	general case: assembling over all elements in subset si
		gass_type::template AssembleRhs<TElem>
			(vElemDisc, m_spApproxSpace->domain(), dd, dd->begin<TElem>(si), dd->end<TElem>(si), si,
					bNonRegularGrid, rhs, u, m_spAssTuner);
	}
}

///////////////////////////////////////////////////////////////////////////////
// set constraints (stationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
adjust_solution(vector_type& u, ConstSmartPtr<DoFDistribution> dd)
{
	PROFILE_FUNC_GROUP("discretization");
	update_constraints();

	// NOTE: it is crucial, that dirichlet pp are processed before constraints.
	// 	 	 otherwise we may start with an inconsistent solution in the solvers
	std::vector<int> vType(6);
	vType[0] = CT_DIRICHLET;	// hanging (or other constrained) nodes might depend on Dirichlet nodes
	vType[1] = CT_CONSTRAINTS;
	vType[2] = CT_HANGING;
	vType[3] = CT_MAY_DEPEND_ON_HANGING;
	vType[4] = CT_ASSEMBLED;
	vType[5] = CT_DIRICHLET;	// hanging DoFs might be Dirichlet (Dirichlet overrides)

	// if assembling is carried out at one DoF only, u needs to be resized
	if (m_spAssTuner->single_index_assembling_enabled()) u.resize(1);

	try{
	for(size_t i = 0; i < vType.size(); ++i){
		int type = vType[i];
		if(!(m_spAssTuner->constraint_type_enabled(type))) continue;
		for(size_t i = 0; i < m_vConstraint.size(); ++i)
			if(m_vConstraint[i]->type() & type)
			{
				m_vConstraint[i]->set_ass_tuner(m_spAssTuner);
				m_vConstraint[i]->adjust_solution(u, dd, type);
			}
	}

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
template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
prepare_timestep
(
	ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
	number future_time,
	ConstSmartPtr<DoFDistribution> dd
)
{
	PROFILE_FUNC_GROUP("discretization");
//	update the elem discs
	update_disc_items();
	prep_assemble_loop(m_vElemDisc);

//	find out whether grid is regular
	ConstSmartPtr<ISubsetHandler> sh = dd->subset_handler();
	size_t num_subsets = sh->num_subsets();
	bool bNonRegularGrid = false;
	for (size_t si = 0; si < num_subsets; ++si)
		bNonRegularGrid |= !SubsetIsRegularGrid(*sh, si);

//	overrule by regular grid if required
	if(m_spAssTuner->regular_grid_forced()) bNonRegularGrid = false;

//	call assembler's PrepareTimestep
	try
	{
		gass_type::PrepareTimestep(m_vElemDisc, dd, bNonRegularGrid, vSol, future_time, m_spAssTuner);
	}
	UG_CATCH_THROW("DomainDiscretization::prepare_timestep (instationary):" <<
				   " Preparing time step failed.");

	post_assemble_loop(m_vElemDisc);
}

///////////////////////////////////////////////////////////////////////////////
// Prepare Timestep Elem (instationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
prepare_timestep_elem(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
                ConstSmartPtr<DoFDistribution> dd)
{
	PROFILE_FUNC_GROUP("discretization");
//	update the elem discs
	update_disc_items();
	prep_assemble_loop(m_vElemDisc);

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	try{
		CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, dd->subset_handler());
	}UG_CATCH_THROW("'DomainDiscretization': Can not create Subset Groups and Union.");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.size(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = DimensionOfSubset(*dd->subset_handler(), si);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(m_spAssTuner->regular_grid_forced()) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc<TDomain>*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	assemble on suitable elements
		try
		{
		switch(dim)
		{
		case 0:
			this->PrepareTimestepElem<RegularVertex>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			break;
		case 1:
			this->PrepareTimestepElem<RegularEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->PrepareTimestepElem<ConstrainingEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			break;
		case 2:
			this->PrepareTimestepElem<Triangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			this->PrepareTimestepElem<Quadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->PrepareTimestepElem<ConstrainingTriangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			this->PrepareTimestepElem<ConstrainingQuadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			break;
		case 3:
			this->PrepareTimestepElem<Tetrahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			this->PrepareTimestepElem<Pyramid>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			this->PrepareTimestepElem<Prism>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			this->PrepareTimestepElem<Hexahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			this->PrepareTimestepElem<Octahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			break;
		default:
			UG_THROW("DomainDiscretization::prepare_timestep_elem (instationary):"
							"Dimension "<<dim<<" (subset="<<si<<") not supported.");
		}
		}
		UG_CATCH_THROW("DomainDiscretization::prepare_timestep_elem (instationary):"
						" Assembling of elements of Dimension " << dim << " in "
						" subset "<<si<< " failed.");
	}
	post_assemble_loop(m_vElemDisc);
}

/**
 * This function prepares the global discretization for a time-stepping scheme
 * by calling the "prepare_timestep_elem" methods of all passed element
 * discretizations on one given subset.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd				DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in]		vSol			current and previous solutions
 */
template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
template <typename TElem>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
PrepareTimestepElem(const std::vector<IElemDisc<domain_type>*>& vElemDisc,
				ConstSmartPtr<DoFDistribution> dd,
				int si, bool bNonRegularGrid,
				ConstSmartPtr<VectorTimeSeries<vector_type> > vSol)
{
	//	check if only some elements are selected
	if(m_spAssTuner->selected_elements_used())
	{
		std::vector<TElem*> vElem;
		m_spAssTuner->collect_selected_elements(vElem, dd, si);

		//	assembling is carried out only over those elements
		//	which are selected and in subset si
		gass_type::template PrepareTimestepElem<TElem>
			(vElemDisc, m_spApproxSpace->domain(), dd, vElem.begin(), vElem.end(), si,
			 bNonRegularGrid, vSol, m_spAssTuner);
	}
	else
	{
		//	general case: assembling over all elements in subset si
		gass_type::template PrepareTimestepElem<TElem>
			(vElemDisc, m_spApproxSpace->domain(), dd,
				dd->begin<TElem>(si), dd->end<TElem>(si), si,
					bNonRegularGrid, vSol, m_spAssTuner);
	}
}

///////////////////////////////////////////////////////////////////////////////
// Jacobian (instationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
assemble_jacobian(matrix_type& J,
                  ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
                  const number s_a0,
                  ConstSmartPtr<DoFDistribution> dd)
{
	// do not do anything if matrix is supposed to be const
	if (m_spAssTuner->matrix_is_const()) return;

	PROFILE_FUNC_GROUP("discretization");
//	update the elem discs
	update_disc_items();
	prep_assemble_loop(m_vElemDisc);

//	reset matrix to zero and resize
	m_spAssTuner->resize(dd, J);

//	get current time
	const number time = vSol->time(0);

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	try{
		CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, dd->subset_handler());
	}UG_CATCH_THROW("'DomainDiscretization': Can not create Subset Groups and Union.");

//	preprocess -  modifies the solution, used for computing the defect
	ConstSmartPtr<VectorTimeSeries<vector_type> > pModifyU = vSol;
	SmartPtr<VectorTimeSeries<vector_type> > pModifyMemory;
	if( m_spAssTuner->modify_solution_enabled() ){
		pModifyMemory = vSol->clone();
		pModifyU = pModifyMemory;
		try{
		for(int type = 1; type < CT_ALL; type = type << 1){
			if(!(m_spAssTuner->constraint_type_enabled(type))) continue;
			for(size_t i = 0; i < m_vConstraint.size(); ++i)
				if(m_vConstraint[i]->type() & type)
					m_vConstraint[i]->modify_solution(pModifyMemory, vSol, dd, type);
		}
		} UG_CATCH_THROW("'DomainDiscretization': Cannot modify solution.");
	}

//	loop subsets
	for(size_t i = 0; i < unionSubsets.size(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = DimensionOfSubset(*dd->subset_handler(), si);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(m_spAssTuner->regular_grid_forced()) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc<TDomain>*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	assemble on suitable elements
		try
		{
		switch(dim)
		{
		case 0:
			this->AssembleJacobian<RegularVertex>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, pModifyU, s_a0);
			break;
		case 1:
			this->AssembleJacobian<RegularEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, pModifyU, s_a0);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->AssembleJacobian<ConstrainingEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, pModifyU, s_a0);
			break;
		case 2:
			this->AssembleJacobian<Triangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, pModifyU, s_a0);
			this->AssembleJacobian<Quadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, pModifyU, s_a0);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->AssembleJacobian<ConstrainingTriangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, pModifyU, s_a0);
			this->AssembleJacobian<ConstrainingQuadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, pModifyU, s_a0);
			break;
		case 3:
			this->AssembleJacobian<Tetrahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, pModifyU, s_a0);
			this->AssembleJacobian<Pyramid>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, pModifyU, s_a0);
			this->AssembleJacobian<Prism>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, pModifyU, s_a0);
			this->AssembleJacobian<Hexahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, pModifyU, s_a0);
			this->AssembleJacobian<Octahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, pModifyU, s_a0);
			break;
		default:
			UG_THROW("DomainDiscretization::assemble_jacobian (instationary):"
							"Dimension "<<dim<<" (subset="<<si<<") not supported.");
		}
		}
		UG_CATCH_THROW("DomainDiscretization::assemble_jacobian (instationary):"
						" Assembling of elements of Dimension " << dim << " in "
						" subset "<<si<< " failed.");
	}

//	post process
	try{
	for(int type = 1; type < CT_ALL; type = type << 1){
		if(!(m_spAssTuner->constraint_type_enabled(type))) continue;
		for(size_t i = 0; i < m_vConstraint.size(); ++i)
			if(m_vConstraint[i]->type() & type)
			{
				m_vConstraint[i]->set_ass_tuner(m_spAssTuner);
				m_vConstraint[i]->adjust_jacobian(J, *pModifyU->solution(0), dd, type, time, pModifyU,s_a0);
			}
	}
	post_assemble_loop(m_vElemDisc);
	}UG_CATCH_THROW("Cannot adjust jacobian.");

//	Remember parallel storage type
#ifdef UG_PARALLEL
	J.set_storage_type(PST_ADDITIVE);
	J.set_layouts(dd->layouts());
#endif
}

/**
 * This function adds the contributions of all passed element discretizations
 * on one given subset to the global Jacobian in the time-dependent case.
 * Note, that s_m0 == 1
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd				DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	J				jacobian
 * \param[in]		vSol			current and previous solutions
 * \param[in]		s_a0			scaling factor for stiffness part
 */
template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
template <typename TElem>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
AssembleJacobian(	const std::vector<IElemDisc<domain_type>*>& vElemDisc,
					ConstSmartPtr<DoFDistribution> dd,
					int si, bool bNonRegularGrid,
					matrix_type& J,
					ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
					number s_a0)
{
	//	check if only some elements are selected
	if(m_spAssTuner->selected_elements_used())
	{
		std::vector<TElem*> vElem;
		m_spAssTuner->collect_selected_elements(vElem, dd, si);
	
		//	assembling is carried out only over those elements
		//	which are selected and in subset si
		gass_type::template AssembleJacobian<TElem>
			(vElemDisc, m_spApproxSpace->domain(), dd, vElem.begin(), vElem.end(), si,
			 bNonRegularGrid, J, vSol, s_a0, m_spAssTuner);
	}
	else
	{
		gass_type::template AssembleJacobian<TElem>
			(vElemDisc, m_spApproxSpace->domain(), dd,
				dd->begin<TElem>(si), dd->end<TElem>(si), si,
					bNonRegularGrid, J, vSol, s_a0, m_spAssTuner);
	}
}

///////////////////////////////////////////////////////////////////////////////
// Defect (instationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
assemble_defect(vector_type& d,
                ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
                const std::vector<number>& vScaleMass,
                const std::vector<number>& vScaleStiff,
                ConstSmartPtr<DoFDistribution> dd)
{
	PROFILE_FUNC_GROUP("discretization");
//	update the elem discs
	update_disc_items();
	prep_assemble_loop(m_vElemDisc);

//	reset vector to zero and resize
	m_spAssTuner->resize(dd, d);

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	try{
		CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, dd->subset_handler());
	}UG_CATCH_THROW("'DomainDiscretization': Can not create Subset Groups and Union.");

//	pre process -  modifies the solution, used for computing the defect
	ConstSmartPtr<VectorTimeSeries<vector_type> > pModifyU = vSol;
	SmartPtr<VectorTimeSeries<vector_type> > pModifyMemory;
	if( m_spAssTuner->modify_solution_enabled() ){
		pModifyMemory = vSol->clone();
		pModifyU = pModifyMemory;
		try{
		for(int type = 1; type < CT_ALL; type = type << 1){
			if(!(m_spAssTuner->constraint_type_enabled(type))) continue;
			for(size_t i = 0; i < m_vConstraint.size(); ++i)
				if(m_vConstraint[i]->type() & type)
					m_vConstraint[i]->modify_solution(pModifyMemory, vSol, dd, type);
		}
		} UG_CATCH_THROW("'DomainDiscretization: Cannot modify solution.");
	}

//	loop subsets
	for(size_t i = 0; i < unionSubsets.size(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = DimensionOfSubset(*dd->subset_handler(), si);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(m_spAssTuner->regular_grid_forced()) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc<TDomain>*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	assemble on suitable elements
		try
		{
		switch(dim)
		{
		case 0:
			this->AssembleDefect<RegularVertex>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff);
			break;
		case 1:
			this->AssembleDefect<RegularEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->AssembleDefect<ConstrainingEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff);
			break;
		case 2:
			this->AssembleDefect<Triangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff);
			this->AssembleDefect<Quadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->AssembleDefect<ConstrainingTriangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff);
			this->AssembleDefect<ConstrainingQuadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff);
			break;
		case 3:
			this->AssembleDefect<Tetrahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff);
			this->AssembleDefect<Pyramid>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff);
			this->AssembleDefect<Prism>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff);
			this->AssembleDefect<Hexahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff);
			this->AssembleDefect<Octahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff);
			break;
		default:
			UG_THROW("DomainDiscretization::assemble_defect (instationary):"
							"Dimension "<<dim<<" (subset="<<si<<") not supported.");
		}
		}
		UG_CATCH_THROW("DomainDiscretization::assemble_defect (instationary):"
						" Assembling of elements of Dimension " << dim << " in"
						" subset "<< si << " failed.");
	}

//	post process
	try{
	for(int type = 1; type < CT_ALL; type = type << 1){
		if(!(m_spAssTuner->constraint_type_enabled(type))) continue;
		for(size_t i = 0; i < m_vConstraint.size(); ++i)
			if(m_vConstraint[i]->type() & type)
			{
				m_vConstraint[i]->set_ass_tuner(m_spAssTuner);
				m_vConstraint[i]->adjust_defect(d, *pModifyU->solution(0), dd, type, pModifyU->time(0), pModifyU, &vScaleMass, &vScaleStiff);
			}
	}
	post_assemble_loop(m_vElemDisc);
	} UG_CATCH_THROW("Cannot adjust defect.");

//	Remember parallel storage type
#ifdef UG_PARALLEL
	d.set_storage_type(PST_ADDITIVE);
#endif
}

/*
 * This function adds the contributions of all passed element discretizations
 * on one given subset to the global Defect in the instationary case.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd				DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	d				defect
 * \param[in]		vSol			current and previous solutions
 * \param[in]		vScaleMass		scaling factors for mass part
 * \param[in]		vScaleStiff		scaling factors for stiffness part
 */
template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
template <typename TElem>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
AssembleDefect( const std::vector<IElemDisc<domain_type>*>& vElemDisc,
				ConstSmartPtr<DoFDistribution> dd,
				int si, bool bNonRegularGrid,
				vector_type& d,
				ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
				const std::vector<number>& vScaleMass,
				const std::vector<number>& vScaleStiff)
{
	//	check if only some elements are selected
	if(m_spAssTuner->selected_elements_used())
	{
		std::vector<TElem*> vElem;
		m_spAssTuner->collect_selected_elements(vElem, dd, si);
		
		//	assembling is carried out only over those elements
		//	which are selected and in subset si
		gass_type::template AssembleDefect<TElem>
			(vElemDisc, m_spApproxSpace->domain(), dd, vElem.begin(), vElem.end(), si,
			 bNonRegularGrid, d, vSol, vScaleMass, vScaleStiff, m_spAssTuner);
	}
	else
	{
		//	general case: assembling over all elements in subset si
		gass_type::template AssembleDefect<TElem>
			(vElemDisc, m_spApproxSpace->domain(), dd,
				dd->begin<TElem>(si), dd->end<TElem>(si), si,
					bNonRegularGrid, d, vSol, vScaleMass, vScaleStiff, m_spAssTuner);
	}
}

///////////////////////////////////////////////////////////////////////////////
// Matrix and RHS (instationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
assemble_linear(matrix_type& mat, vector_type& rhs,
                ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
                const std::vector<number>& vScaleMass,
                const std::vector<number>& vScaleStiff,
                ConstSmartPtr<DoFDistribution> dd)
{
	PROFILE_FUNC_GROUP("discretization");
//	update the elem discs
	update_disc_items();

//	reset matrix to zero and resize
	if (!m_spAssTuner->matrix_is_const())
		m_spAssTuner->resize(dd, mat);
	m_spAssTuner->resize(dd, rhs);

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	try{
		CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, dd->subset_handler());
	}UG_CATCH_THROW("'DomainDiscretization': Can not create Subset Groups and Union.");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.size(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = DimensionOfSubset(*dd->subset_handler(), si);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(m_spAssTuner->regular_grid_forced()) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc<TDomain>*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	assemble on suitable elements
		try
		{
		switch(dim)
		{
		case 0:
			this->AssembleLinear<RegularVertex>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff);
			break;
		case 1:
			this->AssembleLinear<RegularEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->AssembleLinear<ConstrainingEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff);
			break;
		case 2:
			this->AssembleLinear<Triangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff);
			this->AssembleLinear<Quadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->AssembleLinear<ConstrainingTriangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff);
			this->AssembleLinear<ConstrainingQuadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff);
			break;
		case 3:
			this->AssembleLinear<Tetrahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff);
			this->AssembleLinear<Pyramid>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff);
			this->AssembleLinear<Prism>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff);
			this->AssembleLinear<Hexahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff);
			this->AssembleLinear<Octahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff);
			break;
		default:
			UG_THROW("DomainDiscretization::assemble_linear (instationary):"
							"Dimension "<<dim<<" (subset="<<si<<") not supported.");
		}
		}
		UG_CATCH_THROW("DomainDiscretization::assemble_linear (instationary):"
						" Assembling of elements of Dimension " << dim << " in "
						" subset "<<si<< " failed.");
	}

//	post process
	try{
	for(int type = 1; type < CT_ALL; type = type << 1){
		if(!(m_spAssTuner->constraint_type_enabled(type))) continue;
		for(size_t i = 0; i < m_vConstraint.size(); ++i)
			if(m_vConstraint[i]->type() & type)
			{
				m_vConstraint[i]->set_ass_tuner(m_spAssTuner);
				m_vConstraint[i]->adjust_linear(mat, rhs, dd, type, vSol->time(0));
			}
	}
	} UG_CATCH_THROW("Cannot adjust linear.");

//	Remember parallel storage type
#ifdef UG_PARALLEL
	mat.set_storage_type(PST_ADDITIVE);
	mat.set_layouts(dd->layouts());

	rhs.set_storage_type(PST_ADDITIVE);
#endif
}

/**
 * This function adds the contributions of all passed element discretizations
 * on one given subset to the global Matrix and the global Right-Hand Side
 * of the Linear problem in the stationary case.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd				DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	A				Matrix
 * \param[in,out]	rhs				Right-hand side
 * \param[in]		vSol			current and previous solutions
 * \param[in]		vScaleMass		scaling factors for mass part
 * \param[in]		vScaleStiff		scaling factors for stiffness part
 */
template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
template <typename TElem>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
AssembleLinear( const std::vector<IElemDisc<domain_type>*>& vElemDisc,
				ConstSmartPtr<DoFDistribution> dd,
				int si, bool bNonRegularGrid,
				matrix_type& A,
				vector_type& rhs,
				ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
				const std::vector<number>& vScaleMass,
				const std::vector<number>& vScaleStiff)
{
	//	check if only some elements are selected
	if(m_spAssTuner->selected_elements_used())
	{
		std::vector<TElem*> vElem;
		m_spAssTuner->collect_selected_elements(vElem, dd, si);

		//	assembling is carried out only over those elements
		//	which are selected and in subset si
		gass_type::template AssembleLinear<TElem>
			(vElemDisc, m_spApproxSpace->domain(), dd, vElem.begin(), vElem.end(), si,
			 bNonRegularGrid, A, rhs, vSol, vScaleMass, vScaleStiff, m_spAssTuner);
	}
	else
	{
		//	general case: assembling over all elements in subset si
		gass_type::template AssembleLinear<TElem>
			(vElemDisc, m_spApproxSpace->domain(), dd,
				dd->begin<TElem>(si), dd->end<TElem>(si), si,
					bNonRegularGrid, A, rhs, vSol, vScaleMass, vScaleStiff, m_spAssTuner);
	}
}

///////////////////////////////////////////////////////////////////////////////
// RHS (instationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
assemble_rhs(vector_type& rhs,
             ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
             const std::vector<number>& vScaleMass,
             const std::vector<number>& vScaleStiff,
             ConstSmartPtr<DoFDistribution> dd)
{
	PROFILE_FUNC_GROUP("discretization");
//	update the elem discs
	update_disc_items();

//	reset vector to zero and resize
	m_spAssTuner->resize(dd, rhs);

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	try{
		CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, dd->subset_handler());
	}UG_CATCH_THROW("'DomainDiscretization': Can not create Subset Groups and Union.");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.size(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = DimensionOfSubset(*dd->subset_handler(), si);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(m_spAssTuner->regular_grid_forced()) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc<TDomain>*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	assemble on suitable elements
		try
		{
		switch(dim)
		{
		case 0:
			this->AssembleRhs<RegularVertex>
			 	 (vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff);
			break;
		case 1:
			this->AssembleRhs<RegularEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->AssembleRhs<ConstrainingEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff);
			break;
		case 2:
			this->AssembleRhs<Triangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff);
			this->AssembleRhs<Quadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->AssembleRhs<ConstrainingTriangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff);
			this->AssembleRhs<ConstrainingQuadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff);
			break;
		case 3:
			this->AssembleRhs<Tetrahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff);
			this->AssembleRhs<Pyramid>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff);
			this->AssembleRhs<Prism>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff);
			this->AssembleRhs<Hexahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff);
			this->AssembleRhs<Octahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff);
			break;
		default:
			UG_THROW("DomainDiscretization::assemble_rhs (instationary):"
							"Dimension "<<dim<<" (subset="<<si<<") not supported.");
		}
		}
		UG_CATCH_THROW("DomainDiscretization::assemble_rhs (instationary):"
						" Assembling of elements of Dimension " << dim << " in "
						" subset "<<si<< " failed.");
	}


//	post process
	try{
	for(int type = 1; type < CT_ALL; type = type << 1){
		if(!(m_spAssTuner->constraint_type_enabled(type))) continue;
		for(size_t i = 0; i < m_vConstraint.size(); ++i)
			if(m_vConstraint[i]->type() & type)
			{
				m_vConstraint[i]->set_ass_tuner(m_spAssTuner);
				m_vConstraint[i]->adjust_rhs(rhs, *(vSol->solution(0)), dd, type, vSol->time(0));
			}
	}
	} UG_CATCH_THROW("Cannot adjust linear.");

//	Remember parallel storage type
#ifdef UG_PARALLEL
	rhs.set_storage_type(PST_ADDITIVE);
#endif
}

/**
 * This function adds the contributions of all passed element discretizations
 * on one given subset to the global Right-Hand Side in the time-dependent case.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd				DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in,out]	rhs				Right-hand side
 * \param[in]		vSol			current and previous solutions
 * \param[in]		vScaleMass		scaling factors for mass part
 * \param[in]		vScaleStiff		scaling factors for stiffness part
 */
template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
template <typename TElem>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
AssembleRhs(	const std::vector<IElemDisc<domain_type>*>& vElemDisc,
				ConstSmartPtr<DoFDistribution> dd,
				int si, bool bNonRegularGrid,
				vector_type& rhs,
				ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
				const std::vector<number>& vScaleMass,
				const std::vector<number>& vScaleStiff)
{
	//	check if only some elements are selected
	if(m_spAssTuner->selected_elements_used())
	{
		std::vector<TElem*> vElem;
		m_spAssTuner->collect_selected_elements(vElem, dd, si);

		//	assembling is carried out only over those elements
		//	which are selected and in subset si
		gass_type::template AssembleRhs<TElem>
			(vElemDisc, m_spApproxSpace->domain(), dd, vElem.begin(), vElem.end(), si,
			 bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff, m_spAssTuner);
	}
	else
	{
		//	general case: assembling over all elements in subset si
		gass_type::template AssembleRhs<TElem>
			(vElemDisc, m_spApproxSpace->domain(), dd,
				dd->begin<TElem>(si), dd->end<TElem>(si), si,
					bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff, m_spAssTuner);
	}
}

///////////////////////////////////////////////////////////////////////////////
// set constraint values (instationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
adjust_solution(vector_type& u, number time, ConstSmartPtr<DoFDistribution> dd)
{
	PROFILE_FUNC_GROUP("discretization");
	update_constraints();

	// NOTE: it is crucial, that dirichlet pp are processed before constraints.
	// 	 	 otherwise we may start with an inconsistent solution in the solvers
	std::vector<int> vType(5);
	vType[0] = CT_DIRICHLET;	// hanging (or other constrained) nodes might depend on Dirichlet nodes
	vType[1] = CT_CONSTRAINTS;
	vType[2] = CT_HANGING;
	vType[3] = CT_MAY_DEPEND_ON_HANGING;
	vType[4] = CT_DIRICHLET;	// hanging DoFs might be Dirichlet (Dirichlet overrides)

	// if assembling is carried out at one DoF only, u needs to be resized
	if (m_spAssTuner->single_index_assembling_enabled()) u.resize(1);

	try{

//	constraints
	for(size_t i = 0; i < vType.size(); ++i){
		int type = vType[i];
		if(!(m_spAssTuner->constraint_type_enabled(type))) continue;
		for(size_t i = 0; i < m_vConstraint.size(); ++i)
			if(m_vConstraint[i]->type() & type)
			{
				m_vConstraint[i]->set_ass_tuner(m_spAssTuner);
				m_vConstraint[i]->adjust_solution(u, dd, type, time);
			}
	}
	} UG_CATCH_THROW(" Cannot adjust solution.");
}

///////////////////////////////////////////////////////////////////////////////
// Initialization of the exports (optional)
///////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
init_all_exports(ConstSmartPtr<DoFDistribution> dd, bool bAsTimeDependent)
{
	PROFILE_FUNC_GROUP("discretization");
//	update the elem discs
	update_disc_items();

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	try{
		CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, dd->subset_handler());
	}UG_CATCH_THROW("'DomainDiscretization': Can not create Subset Groups and Union.");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.size(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = DimensionOfSubset(*dd->subset_handler(), si);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(m_spAssTuner->regular_grid_forced()) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc<TDomain>*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	assemble on suitable elements
		try
		{
		switch(dim)
		{
		case 0:
			this->InitAllExports<RegularVertex>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, bAsTimeDependent);
			break;
		case 1:
			this->InitAllExports<RegularEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, bAsTimeDependent);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->InitAllExports<ConstrainingEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, bAsTimeDependent);
			break;
		case 2:
			this->InitAllExports<Triangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, bAsTimeDependent);
			this->InitAllExports<Quadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, bAsTimeDependent);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->InitAllExports<ConstrainingTriangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, bAsTimeDependent);
			this->InitAllExports<ConstrainingQuadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, bAsTimeDependent);
			break;
		case 3:
			this->InitAllExports<Tetrahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, bAsTimeDependent);
			this->InitAllExports<Pyramid>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, bAsTimeDependent);
			this->InitAllExports<Prism>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, bAsTimeDependent);
			this->InitAllExports<Hexahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, bAsTimeDependent);
			this->InitAllExports<Octahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, bAsTimeDependent);
			break;
		default:
			UG_THROW("DomainDiscretization::init_all_exports:"
							"Dimension "<<dim<<" (subset="<<si<<") not supported.");
		}
		}
		UG_CATCH_THROW("DomainDiscretization::assemble_stiffness_matrix:"
					" Assembling of elements of Dimension " << dim << " in "
					" subset "<<si<< " failed.");
	}

}

/**
 * This function initalizes the exports and does not do anything else.
 *
 * It transfers the dof distributions to the exports and creates the function
 * maps there.
 *
 * This initizalization is optional but after that, the exports can be used
 * by other objects (like in the plotting).
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd				DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in]		bAsTimeDependent flag to simulate the instationary case for the initialization
 */
template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
template <typename TElem>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
InitAllExports(	const std::vector<IElemDisc<domain_type>*>& vElemDisc,
							ConstSmartPtr<DoFDistribution> dd,
							int si, bool bNonRegularGrid, bool bAsTimeDependent)
{
	//	check if only some elements are selected
	if(m_spAssTuner->selected_elements_used())
	{
		std::vector<TElem*> vElem;
		m_spAssTuner->collect_selected_elements(vElem, dd, si);

		//	assembling is carried out only over those elements
		//	which are selected and in subset si
		gass_type::template InitAllExports<TElem>
			(vElemDisc, dd, vElem.begin(), vElem.end(), si, bNonRegularGrid, bAsTimeDependent);
	}
	else
	{
		//	general case: assembling over all elements in subset si
		gass_type::template InitAllExports<TElem>
			(vElemDisc, dd,
				dd->begin<TElem>(si), dd->end<TElem>(si), si,
					bNonRegularGrid, bAsTimeDependent);
	}
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Error estimator
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename T>
void CollectIErrEstData(std::vector<IErrEstData<TDomain>*> &vErrEstData, const T &vElemDisc)
{
	for (std::size_t i = 0; i < vElemDisc.size(); ++i)
	{
		SmartPtr<IErrEstData<TDomain> > sp_err_est_data = vElemDisc[i]->err_est_data();
		IErrEstData<TDomain>* err_est_data = sp_err_est_data.get();
			if (err_est_data == nullptr) continue; // no data specified
			if (std::find (vErrEstData.begin(), vErrEstData.end(), err_est_data) != vErrEstData.end())
				continue; // this one is already in the array
			if (err_est_data->consider_me()) vErrEstData.push_back(err_est_data);
	}

}

///////////////////////////////////////////////////////////////////////////////
// Error estimator (stationary)
///////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
calc_error
(
	const vector_type& u,
	ConstSmartPtr<DoFDistribution> dd,
	error_vector_type* u_vtk
)
{
	PROFILE_FUNC_GROUP("error_estimator");

//	get multigrid
	SmartPtr<MultiGrid> pMG = const_cast<DoFDistribution *>(dd.get())->multi_grid();

// check, whether separate error data exists
	const bool useErrorData = !m_vDomainElemError.empty();
	// UG_LOG("useErrorData=" << (useErrorData) << std::endl);

//	update the elem discs
	if (useErrorData) { update_error_items();}
	else { update_disc_items();}

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	try
	{
		if (useErrorData) { CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemError, dd->subset_handler());}
		else { CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, dd->subset_handler());}
	}
	UG_CATCH_THROW("'DomainDiscretization': Can not create Subset Groups and Union.");

//	get the error estimator data for all the discretizations
	std::vector<IErrEstData<TDomain>*> vErrEstData;
	if (useErrorData) CollectIErrEstData(vErrEstData, m_vElemError);
	else  CollectIErrEstData(vErrEstData, m_vElemDisc);
	/*for(std::size_t i = 0; i < m_vElemDisc.size(); ++i)
	{
		SmartPtr<IErrEstData<TDomain> > sp_err_est_data = m_vElemDisc[i]->err_est_data();
		IErrEstData<TDomain>* err_est_data = sp_err_est_data.get();
		if (err_est_data == nullptr) continue; // no data specified
		if (std::find (vErrEstData.begin(), vErrEstData.end(), err_est_data) != vErrEstData.end())
			continue; // this one is already in the array
		if (err_est_data->consider_me()) vErrEstData.push_back(err_est_data);
	}*/

//	preprocess the error estimator data in the discretizations
	try
	{
		for (size_t i = 0; i < vErrEstData.size(); ++i)
			vErrEstData[i]->alloc_err_est_data(dd->surface_view(), dd->grid_level());
	}
	UG_CATCH_THROW("DomainDiscretization::calc_error: Cannot prepare the error estimator");



//	loop subsets to assemble the estimators
	for (size_t i = 0; i < unionSubsets.size(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = DimensionOfSubset(*dd->subset_handler(), si);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	Elem Disc on the subset
		using error_vector_type = std::vector<IElemError<TDomain>*>;
		error_vector_type vSubsetElemError;

	//	get all element discretizations that work on the subset
		GetElemDiscItemOnSubset<IElemError<TDomain>, IElemDisc<TDomain> > (vSubsetElemError, m_vElemDisc, vSSGrp, si);
		if (useErrorData)
		{
			GetElemDiscItemOnSubset<IElemError<TDomain>, IElemError<TDomain> >
			(vSubsetElemError, m_vElemError, vSSGrp, si);
		}
		else
		{
			GetElemDiscItemOnSubset<IElemError<TDomain>, IElemDisc<TDomain> >
			(vSubsetElemError, m_vElemDisc, vSSGrp, si);
		}

		/*
		UG_LOG("m_vElemDisc.size="<<m_vElemDisc.size()<< std::endl);
		UG_LOG("m_vElemError.size="<<m_vElemError.size()<< std::endl);
		UG_LOG("vSubsetElemError.size="<<vSubsetElemError.size()<< std::endl);
		*/

	//	remove from elemDisc list those with !err_est_enabled()
		typename error_vector_type::iterator it = vSubsetElemError.begin();
		while (it != vSubsetElemError.end())
		{
			if (!(*it)->err_est_enabled())
				it = vSubsetElemError.erase(it);
			else ++it;
		}
	//	UG_LOG("vSubsetElemError.size="<<vSubsetElemError.size()<< std::endl);

	//	assemble on suitable elements
		try
		{
			switch (dim)
			{
			case 1:
				this->AssembleErrorEstimator<RegularEdge>
					(vSubsetElemError, dd, si, bNonRegularGrid, u);
				// When assembling over lower-dim manifolds that contain hanging nodes:
				this->AssembleErrorEstimator<ConstrainingEdge>
					(vSubsetElemError, dd, si, bNonRegularGrid, u);
				break;
			case 2:
				this->AssembleErrorEstimator<Triangle>
					(vSubsetElemError, dd, si, bNonRegularGrid, u);
				this->AssembleErrorEstimator<Quadrilateral>
					(vSubsetElemError, dd, si, bNonRegularGrid, u);
				// When assembling over lower-dim manifolds that contain hanging nodes:
				this->AssembleErrorEstimator<ConstrainingTriangle>
					(vSubsetElemError, dd, si, bNonRegularGrid, u);
				this->AssembleErrorEstimator<ConstrainingQuadrilateral>
					(vSubsetElemError, dd, si, bNonRegularGrid, u);
				break;
			case 3:
				this->AssembleErrorEstimator<Tetrahedron>
					(vSubsetElemError, dd, si, bNonRegularGrid, u);
				this->AssembleErrorEstimator<Pyramid>
					(vSubsetElemError, dd, si, bNonRegularGrid, u);
				this->AssembleErrorEstimator<Prism>
					(vSubsetElemError, dd, si, bNonRegularGrid, u);
				this->AssembleErrorEstimator<Hexahedron>
					(vSubsetElemError, dd, si, bNonRegularGrid, u);
				this->AssembleErrorEstimator<Octahedron>
					(vSubsetElemError, dd, si, bNonRegularGrid, u);
				break;
			default:
				UG_THROW("DomainDiscretization::calc_error:"
								" Dimension "<<dim<<" (subset="<<si<<") not supported.");
			}
		}
		UG_CATCH_THROW("DomainDiscretization::calc_error:"
						" Assembling of elements of Dimension " << dim << " in "
						" subset "<<si<< " failed.");
	}

//	apply constraints
	try
	{
		for (int type = 1; type < CT_ALL; type = type << 1)
		{
			if (!(m_spAssTuner->constraint_type_enabled(type))) continue;
			for (size_t i = 0; i < m_vConstraint.size(); ++i)
				if ((m_vConstraint[i]->type() & type) && m_vConstraint[i]->err_est_enabled())
				{
					m_vConstraint[i]->set_ass_tuner(m_spAssTuner);
					m_vConstraint[i]->adjust_error(u, dd, type);
				}
		}
	}
	UG_CATCH_THROW("Cannot adjust error assemblings.");

//	summarize the error estimator data in the discretizations
	try
	{
		for (std::size_t i = 0; i < vErrEstData.size(); ++i)
			vErrEstData[i]->summarize_err_est_data(m_spApproxSpace->domain());
	}
	UG_CATCH_THROW("DomainDiscretization::calc_error: Cannot summarize the error estimator");

// perform integrations for error estimators and mark elements
	using elem_type = typename domain_traits<dim>::element_type;
	using elem_iter_type = typename SurfaceView::traits<elem_type>::const_iterator;

	m_mgElemErrors.attach_indicators(pMG);

	// loop surface elements
	ConstSmartPtr<SurfaceView> sv = dd->surface_view();
	const GridLevel& gl = dd->grid_level();
	elem_iter_type elem_iter_end = sv->end<elem_type> (gl, SurfaceView::ALL);
	for (elem_iter_type elem = sv->begin<elem_type> (gl, SurfaceView::ALL); elem != elem_iter_end; ++elem)
	{
		// clear attachment (to be on the safe side)
		m_mgElemErrors.error(*elem) = 0.0;

		// get corner coordinates
		std::vector<MathVector<dim> > vCornerCoords = std::vector<MathVector<dim> >(0);
		CollectCornerCoordinates(vCornerCoords, *elem, m_spApproxSpace->domain()->position_accessor(), false);

		// integrate for all estimators, then add up
		for (std::size_t ee = 0; ee < vErrEstData.size(); ++ee)
			m_mgElemErrors.error(*elem) += vErrEstData[ee]->scaling_factor()
								* vErrEstData[ee]->get_elem_error_indicator(*elem, &vCornerCoords[0]);
	}

//	write error estimator values to vtk
	if (u_vtk)
	{
		// local indices and solution
		LocalIndices ind; LocalVector locU;

		// cast u_vtk to grid_function
		auto* uVTK = dynamic_cast<GridFunction<TDomain, CPUAlgebra>*>(u_vtk);
		if (!uVTK)
		{
			UG_THROW("Argument passed as output for error function is not a GridFunction of suitable type.");
		}

		// clear previous values
		uVTK->set(0.0);

		// map attachments to grid function
		ConstSmartPtr<SurfaceView> sv = uVTK->approx_space()->dof_distribution(gl)->surface_view();
		elem_iter_type elem_iter_end = sv->end<elem_type> (gl, SurfaceView::ALL);
		for (elem_iter_type elem = sv->begin<elem_type> (gl, SurfaceView::ALL); elem != elem_iter_end; ++elem)
		{
			// 	get global indices
			uVTK->approx_space()->dof_distribution(gl)->indices(*elem, ind, false);

			UG_COND_THROW(ind.num_fct() != 1,
				"Number of functions in grid function passed for error indicator values is not 1 on "
				<< ElementDebugInfo(*uVTK->domain()->grid(), *elem) << ".");

			UG_COND_THROW(ind.num_dof(0) != 1,
				"Number of DoFs in grid function passed for error indicator values is not 1 on "
				<< ElementDebugInfo(*uVTK->domain()->grid(), *elem) << ".");

			// 	adapt local algebra
			locU.resize(ind);

			// 	read local values of u
			GetLocalVector(locU, *uVTK);

			// assign error value
			locU(0,0) = m_mgElemErrors.error(*elem);

			// add to grid function
			AddLocalVector(*uVTK, locU);
		}
	}

	this->m_bErrorCalculated = true;

//	postprocess the error estimators in the discretizations
	try{
		for(std::size_t i = 0; i < vErrEstData.size(); ++i)
			vErrEstData[i]->release_err_est_data();
	}
	UG_CATCH_THROW("DomainDiscretization::calc_error: Cannot release the error estimator");
}

/**
 * This function assembles the error estimator associated with all the
 * element discretizations in the internal data structure for one given subset.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd				DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in]		u				solution
 */
template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
template <typename TElem>
inline void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
AssembleErrorEstimator(	const std::vector<IElemError<domain_type>*>& vElemDisc,
						ConstSmartPtr<DoFDistribution> dd,
						int si, bool bNonRegularGrid,
						const vector_type& u)
{
	//	general case: assembling over all elements in subset si
	gass_type::template AssembleErrorEstimator<TElem>
		(vElemDisc, m_spApproxSpace->domain(), dd,
			dd->begin<TElem>(si), dd->end<TElem>(si),
				si, bNonRegularGrid, u);
}

///////////////////////////////////////////////////////////////////////////////
// Error estimator (instationary)
///////////////////////////////////////////////////////////////////////////////



template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
calc_error(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		   ConstSmartPtr<DoFDistribution> dd,
		   const std::vector<number>& vScaleMass,
		   const std::vector<number>& vScaleStiff,
		   error_vector_type* u_vtk)
{
	PROFILE_FUNC_GROUP("error_estimator");

//	get multigrid
	SmartPtr<MultiGrid> pMG = const_cast<DoFDistribution *>(dd.get())->multi_grid();

// check, whether separate error data exists
	const bool useErrorData = !m_vDomainElemError.empty();

//	update elem items
	if (useErrorData) { update_error_items();}
	else { update_disc_items();}

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	try
	{
		if (useErrorData) {CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemError, dd->subset_handler());}
		else {CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, dd->subset_handler());}
	}
	UG_CATCH_THROW("'DomainDiscretization': Can not create Subset Groups and Union.");

//	collect the error estimator data (for all discretizations)
	std::vector<IErrEstData<TDomain>*> vErrEstData;

	if (useErrorData) CollectIErrEstData(vErrEstData, m_vElemError);
	else CollectIErrEstData(vErrEstData, m_vElemDisc);

//	preprocess the error estimator data in the discretizations
	try
	{
		for (std::size_t i = 0; i < vErrEstData.size(); ++i)
			vErrEstData[i]->alloc_err_est_data(dd->surface_view(), dd->grid_level());
	}
	UG_CATCH_THROW("DomainDiscretization::calc_error: Cannot prepare the error estimator");

//	loop subsets to assemble the estimators
	for (std::size_t i = 0; i < unionSubsets.size(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = DimensionOfSubset(*dd->subset_handler(), si);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	collect all element discretizations that work on the subset
	using error_vector_type = std::vector<IElemError<TDomain>*>;
		error_vector_type vSubsetElemError;

		if (useErrorData) {
			 GetElemDiscItemOnSubset<IElemError<TDomain>, IElemError<TDomain> >
							(vSubsetElemError, m_vElemError, vSSGrp, si);
		}
		else
		{
			 GetElemDiscItemOnSubset<IElemError<TDomain>, IElemDisc<TDomain> >
										(vSubsetElemError, m_vElemDisc, vSSGrp, si);
		}

	//	remove from elemDisc list those with !err_est_enabled()
		typename error_vector_type::iterator it = vSubsetElemError.begin();
		while (it != vSubsetElemError.end())
		{
			if (!(*it)->err_est_enabled())
				it = vSubsetElemError.erase(it);
			else ++it;
		}

	//	assemble on suitable elements
		try
		{
			switch (dim)
			{
			case 1:
				this->AssembleErrorEstimator<RegularEdge>
					(vSubsetElemError, dd, si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
				// When assembling over lower-dim manifolds that contain hanging nodes:
				this->AssembleErrorEstimator<ConstrainingEdge>
					(vSubsetElemError, dd, si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
				break;
			case 2:
				this->AssembleErrorEstimator<Triangle>
					(vSubsetElemError, dd, si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
				this->AssembleErrorEstimator<Quadrilateral>
					(vSubsetElemError, dd, si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
				// When assembling over lower-dim manifolds that contain hanging nodes:
				this->AssembleErrorEstimator<ConstrainingTriangle>
					(vSubsetElemError, dd, si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
				this->AssembleErrorEstimator<ConstrainingQuadrilateral>
					(vSubsetElemError, dd, si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
				break;
			case 3:
				this->AssembleErrorEstimator<Tetrahedron>
					(vSubsetElemError, dd, si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
				this->AssembleErrorEstimator<Pyramid>
					(vSubsetElemError, dd, si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
				this->AssembleErrorEstimator<Prism>
					(vSubsetElemError, dd, si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
				this->AssembleErrorEstimator<Hexahedron>
					(vSubsetElemError, dd, si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
				this->AssembleErrorEstimator<Octahedron>
					(vSubsetElemError, dd, si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
				break;
			default:
				UG_THROW("DomainDiscretization::calc_error:"
								" Dimension "<<dim<<" (subset="<<si<<") not supported.");
			}
		}
		UG_CATCH_THROW("DomainDiscretization::calc_error:"
						" Assembling of elements of Dimension " << dim << " in "
						" subset "<< si << "failed.");
	}

//	apply constraints
	try
	{
		for (int type = 1; type < CT_ALL; type = type << 1)
		{
			if (!(m_spAssTuner->constraint_type_enabled(type))) continue;
			for (size_t i = 0; i < m_vConstraint.size(); ++i)
				if ((m_vConstraint[i]->type() & type) && m_vConstraint[i]->err_est_enabled())
				{
					m_vConstraint[i]->set_ass_tuner(m_spAssTuner);
					m_vConstraint[i]->adjust_error(*vSol->solution(0), dd, type, vSol->time(0),
						vSol, &vScaleMass, &vScaleStiff);
				}
		}
	}
	UG_CATCH_THROW("Cannot adjust error assemblings.");

//	summarize the error estimator data in the discretizations
	try
	{
		for (std::size_t i = 0; i < vErrEstData.size(); ++i)
			vErrEstData[i]->summarize_err_est_data(m_spApproxSpace->domain());
	}
	UG_CATCH_THROW("DomainDiscretization::calc_error: Cannot summarize the error estimator.");

	// perform integrations for error estimators and mark elements
	using elem_type = typename domain_traits<dim>::element_type;
	using elem_iter_type = typename SurfaceView::traits<elem_type>::const_iterator;

	m_mgElemErrors.attach_indicators(pMG);

	// loop surface elements
	ConstSmartPtr<SurfaceView> sv = dd->surface_view();
	const GridLevel& gl = dd->grid_level();
	elem_iter_type elem_iter_end = sv->end<elem_type> (gl, SurfaceView::ALL);
	for (elem_iter_type elem = sv->begin<elem_type> (gl, SurfaceView::ALL); elem != elem_iter_end; ++elem)
	{
		// clear attachment
		m_mgElemErrors.error(*elem) = 0.0;

		// get corner coordinates
		std::vector<MathVector<dim> > vCornerCoords = std::vector<MathVector<dim> >(0);
		CollectCornerCoordinates(vCornerCoords, *elem, m_spApproxSpace->domain()->position_accessor(), false);

		// integrate for all estimators, then add up
		for (std::size_t ee = 0; ee < vErrEstData.size(); ++ee)
			m_mgElemErrors.error(*elem) += vErrEstData[ee]->scaling_factor()
								* vErrEstData[ee]->get_elem_error_indicator(*elem, &vCornerCoords[0]);
	}

//	write error estimator values to vtk
	if (u_vtk)
	{
		// local indices and solution
		LocalIndices ind; LocalVector locU;

		// cast u_vtk to grid_function
		auto* uVTK = dynamic_cast<GridFunction<TDomain,TAlgebra>*>(u_vtk);
		if (!uVTK)
		{
			UG_THROW("Argument passed as output for error function is not a GridFunction.");
		}

		// clear previous values
		uVTK->set(0.0);

		// map attachments to grid function
		ConstSmartPtr<SurfaceView> sv = uVTK->approx_space()->dof_distribution(gl)->surface_view();
		elem_iter_type elem_iter_end = sv->end<elem_type> (gl, SurfaceView::ALL);
		for (elem_iter_type elem = sv->begin<elem_type> (gl, SurfaceView::ALL); elem != elem_iter_end; ++elem)
		{
			// 	get global indices
			uVTK->approx_space()->dof_distribution(gl)->indices(*elem, ind, false);

			// 	adapt local algebra
			locU.resize(ind);

			// 	read local values of u
			GetLocalVector(locU, *uVTK);

			// assign error value
			locU(0,0) = m_mgElemErrors.error(*elem);

			// add to grid function
			AddLocalVector(*uVTK, locU);
		}
	}

	this->m_bErrorCalculated = true;

//	postprocess the error estimators in the discretizations
	try{
		for(std::size_t i = 0; i < vErrEstData.size(); ++i)
			vErrEstData[i]->release_err_est_data();
	}
	UG_CATCH_THROW("DomainDiscretization::calc_error: Cannot release the error estimator");
}

/**
 * This function assembles the error estimator associated with all the
 * element discretizations in the internal data structure for one given subset.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd				DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in]		vScaleMass		scaling factors for mass part
 * \param[in]		vScaleStiff		scaling factors for stiffness part
 * \param[in]		vSol				solution
 */
template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
template <typename TElem>
inline void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
AssembleErrorEstimator(	const std::vector<IElemError<domain_type>*>& vElemDisc,
						ConstSmartPtr<DoFDistribution> dd,
						int si, bool bNonRegularGrid,
						const std::vector<number>& vScaleMass,
						const std::vector<number>& vScaleStiff,
						ConstSmartPtr<VectorTimeSeries<vector_type> > vSol)
{
	//	general case: assembling over all elements in subset si
	gass_type::template AssembleErrorEstimator<TElem>
		(vElemDisc, m_spApproxSpace->domain(), dd,
			dd->begin<TElem>(si), dd->end<TElem>(si),
				si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
}

template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
mark_with_strategy
(	IRefiner& refiner,
	SmartPtr<IElementMarkingStrategy<TDomain> >spMarkingStrategy
)
{
	// check that error indicators have been calculated
	if (!this->m_bErrorCalculated)
	{
		UG_THROW("Error indicators have to be calculated first by a call to 'calc_error'.");
	}

	// mark elements for refinement
	if (spMarkingStrategy.valid())
	{
		spMarkingStrategy->mark( m_mgElemErrors, refiner, this->dd(GridLevel(GridLevel::TOP, GridLevel::ViewType::SURFACE)));
	}
}

template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
invalidate_error()
{
	using elem_type = typename domain_traits<dim>::element_type;

	// check that error indicators have been calculated
	if (m_bErrorCalculated)
	{
		m_bErrorCalculated = false;
		m_mgElemErrors.detach_indicators();
		// this->m_pMG->detach_from<elem_type>(this->m_aError);
	}
}

template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
bool DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
is_error_valid()
{
	// check that error indicators have been calculated
	return this->m_bErrorCalculated;
}

///////////////////////////////////////////////////////////////////////////////
// Finish Timestep (instationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
finish_timestep
(
	ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
	ConstSmartPtr<DoFDistribution> dd
)
{
	PROFILE_FUNC_GROUP("discretization");
//	update the elem discs
	update_disc_items();
	prep_assemble_loop(m_vElemDisc);

//	find out whether grid is regular
	ConstSmartPtr<ISubsetHandler> sh = dd->subset_handler();
	size_t num_subsets = sh->num_subsets();
	bool bNonRegularGrid = false;
	for (size_t si = 0; si < num_subsets; ++si)
		bNonRegularGrid |= !SubsetIsRegularGrid(*sh, si);

//	overrule by regular grid if required
	if(m_spAssTuner->regular_grid_forced()) bNonRegularGrid = false;

//	call assembler's FinishTimestep
	try
	{
		gass_type::FinishTimestep(m_vElemDisc, dd, bNonRegularGrid, vSol, m_spAssTuner);
	}
	UG_CATCH_THROW("DomainDiscretization::finish_timestep (instationary):" <<
				   " Finishing time step failed.");

	post_assemble_loop(m_vElemDisc);
}

///////////////////////////////////////////////////////////////////////////////
// Finish Timestep Elem (instationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
finish_timestep_elem(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
                ConstSmartPtr<DoFDistribution> dd)
{
	PROFILE_FUNC_GROUP("discretization");
//	update the elem discs
	update_disc_items();

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	try{
		CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, dd->subset_handler());
	}UG_CATCH_THROW("'DomainDiscretization': Can not create Subset Groups and Union.");

//	loop subsets
	for(size_t i = 0; i < unionSubsets.size(); ++i)
	{
	//	get subset
		const int si = unionSubsets[i];

	//	get dimension of the subset
		const int dim = DimensionOfSubset(*dd->subset_handler(), si);

	//	request if subset is regular grid
		bool bNonRegularGrid = !unionSubsets.regular_grid(i);

	//	overrule by regular grid if required
		if(m_spAssTuner->regular_grid_forced()) bNonRegularGrid = false;

	//	Elem Disc on the subset
		std::vector<IElemDisc<TDomain>*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	assemble on suitable elements
		try
		{
		switch(dim)
		{
		case 0:
			this->FinishTimestepElem<RegularVertex>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			break;
		case 1:
			this->FinishTimestepElem<RegularEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->FinishTimestepElem<ConstrainingEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			break;
		case 2:
			this->FinishTimestepElem<Triangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			this->FinishTimestepElem<Quadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->FinishTimestepElem<ConstrainingTriangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			this->FinishTimestepElem<ConstrainingQuadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			break;
		case 3:
			this->FinishTimestepElem<Tetrahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			this->FinishTimestepElem<Pyramid>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			this->FinishTimestepElem<Prism>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			this->FinishTimestepElem<Hexahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			this->FinishTimestepElem<Octahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			break;
		default:
			UG_THROW("DomainDiscretization::finish_timestep_elem (instationary):"
							"Dimension "<<dim<<" (subset="<<si<<") not supported.");
		}
		}
		UG_CATCH_THROW("DomainDiscretization::finish_timestep_elem (instationary):"
						" Assembling of elements of Dimension " << dim << " in "
						" subset "<<si<< " failed.");
	}

}

/**
 * This function finalizes the global discretization in a time-stepping scheme
 * by calling the "finish_timestep_elem" methods of all passed element
 * discretizations on one given subset.
 *
 * \param[in]		vElemDisc		element discretizations
 * \param[in]		dd				DoF Distribution
 * \param[in]		si				subset index
 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
 * \param[in]		vSol			current and previous solutions
 */
template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
template <typename TElem>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
FinishTimestepElem(const std::vector<IElemDisc<domain_type>*>& vElemDisc,
			   ConstSmartPtr<DoFDistribution> dd,
			   int si, bool bNonRegularGrid,
			   ConstSmartPtr<VectorTimeSeries<vector_type> > vSol)
{
	//	check if only some elements are selected
	if(m_spAssTuner->selected_elements_used())
	{
		std::vector<TElem*> vElem;
		m_spAssTuner->collect_selected_elements(vElem, dd, si);

		//	assembling is carried out only over those elements
		//	which are selected and in subset si
		gass_type::template FinishTimestepElem<TElem>
			(vElemDisc, m_spApproxSpace->domain(), dd, vElem.begin(), vElem.end(), si,
			 bNonRegularGrid, vSol, m_spAssTuner);
	}
	else
	{
		//	general case: assembling over all elements in subset si
		gass_type::template FinishTimestepElem<TElem>
			(vElemDisc, m_spApproxSpace->domain(), dd,
				dd->begin<TElem>(si), dd->end<TElem>(si), si,
					bNonRegularGrid, vSol, m_spAssTuner);
	}
}

} // end namespace ug

#endif