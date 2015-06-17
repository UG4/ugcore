/*
 * domain_disc_impl.h
 *
 *  Created on: 29.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DOMAIN_DISC_IMPL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DOMAIN_DISC_IMPL__

#include "common/profiler/profiler.h"
#include "domain_disc.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/function_spaces/error_indicator_util.h"
#ifdef UG_PARALLEL
#include "lib_disc/parallelization/parallelization_util.h"
#endif

namespace ug{

template <class TElemDisc>
static void prep_assemble_loop(std::vector<TElemDisc*> vElemDisc)
{
	for(size_t i = 0; i < vElemDisc.size(); ++i)
	{
		vElemDisc[i]->prep_assemble_loop();
	}
}

template <class TElemDisc>
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
		case 1:
			this->template AssembleMassMatrix<RegularEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->template AssembleMassMatrix<ConstrainingEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u);
			break;
		case 2:
			this->template AssembleMassMatrix<Triangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u);
			this->template AssembleMassMatrix<Quadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->template AssembleMassMatrix<ConstrainingTriangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u);
			this->template AssembleMassMatrix<ConstrainingQuadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u);
			break;
		case 3:
			this->template AssembleMassMatrix<Tetrahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u);
			this->template AssembleMassMatrix<Pyramid>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u);
			this->template AssembleMassMatrix<Prism>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u);
			this->template AssembleMassMatrix<Hexahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, M, u);
			this->template AssembleMassMatrix<Octahedron>
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
				m_vConstraint[i]->adjust_jacobian(M, u, dd);
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
				dd->template begin<TElem>(si), dd->template end<TElem>(si), si,
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
		case 1:
			this->template AssembleStiffnessMatrix<RegularEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->template AssembleStiffnessMatrix<ConstrainingEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u);
			break;
		case 2:
			this->template AssembleStiffnessMatrix<Triangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u);
			this->template AssembleStiffnessMatrix<Quadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->template AssembleStiffnessMatrix<ConstrainingTriangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u);
			this->template AssembleStiffnessMatrix<ConstrainingQuadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u);
			break;
		case 3:
			this->template AssembleStiffnessMatrix<Tetrahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u);
			this->template AssembleStiffnessMatrix<Pyramid>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u);
			this->template AssembleStiffnessMatrix<Prism>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u);
			this->template AssembleStiffnessMatrix<Hexahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, A, u);
			this->template AssembleStiffnessMatrix<Octahedron>
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
				m_vConstraint[i]->adjust_jacobian(A, u, dd);
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
				dd->template begin<TElem>(si), dd->template end<TElem>(si), si,
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
					m_vConstraint[i]->modify_solution(*pModifyMemory, u, dd);
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
		case 1:
			this->template AssembleJacobian<RegularEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, *pModifyU);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->template AssembleJacobian<ConstrainingEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, *pModifyU);
			break;
		case 2:
			this->template AssembleJacobian<Triangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, *pModifyU);
			this->template AssembleJacobian<Quadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, *pModifyU);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->template AssembleJacobian<ConstrainingTriangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, *pModifyU);
			this->template AssembleJacobian<ConstrainingQuadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, *pModifyU);
			break;
		case 3:
			this->template AssembleJacobian<Tetrahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, *pModifyU);
			this->template AssembleJacobian<Pyramid>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, *pModifyU);
			this->template AssembleJacobian<Prism>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, *pModifyU);
			this->template AssembleJacobian<Hexahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, *pModifyU);
			this->template AssembleJacobian<Octahedron>
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
				m_vConstraint[i]->adjust_jacobian(J, *pModifyU, dd);
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
				dd->template begin<TElem>(si), dd->template end<TElem>(si), si,
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
					m_vConstraint[i]->modify_solution(*pModifyMemory, u, dd);
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
		case 1:
			this->template AssembleDefect<RegularEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, *pModifyU);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->template AssembleDefect<ConstrainingEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, *pModifyU);
			break;
		case 2:
			this->template AssembleDefect<Triangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, *pModifyU);
			this->template AssembleDefect<Quadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, *pModifyU);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->template AssembleDefect<ConstrainingTriangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, *pModifyU);
			this->template AssembleDefect<ConstrainingQuadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, *pModifyU);
			break;
		case 3:
			this->template AssembleDefect<Tetrahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, *pModifyU);
			this->template AssembleDefect<Pyramid>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, *pModifyU);
			this->template AssembleDefect<Prism>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, *pModifyU);
			this->template AssembleDefect<Hexahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, *pModifyU);
			this->template AssembleDefect<Octahedron>
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
	for(int type = 1; type < CT_ALL; type = type << 1){
		if(!(m_spAssTuner->constraint_type_enabled(type))) continue;
		for(size_t i = 0; i < m_vConstraint.size(); ++i)
			if(m_vConstraint[i]->type() & type)
			{
				m_vConstraint[i]->set_ass_tuner(m_spAssTuner);
				m_vConstraint[i]->adjust_defect(d, *pModifyU, dd);
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
				dd->template begin<TElem>(si), dd->template end<TElem>(si), si,
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
		case 1:
			this->template AssembleLinear<RegularEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->template AssembleLinear<ConstrainingEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs);
			break;
		case 2:
			this->template AssembleLinear<Triangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs);
			this->template AssembleLinear<Quadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->template AssembleLinear<ConstrainingTriangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs);
			this->template AssembleLinear<ConstrainingQuadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs);
			break;
		case 3:
			this->template AssembleLinear<Tetrahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs);
			this->template AssembleLinear<Pyramid>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs);
			this->template AssembleLinear<Prism>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs);
			this->template AssembleLinear<Hexahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs);
			this->template AssembleLinear<Octahedron>
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
				m_vConstraint[i]->adjust_linear(mat, rhs, dd);
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
				dd->template begin<TElem>(si), dd->template end<TElem>(si), si,
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
		case 1:
			this->template AssembleRhs<RegularEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->template AssembleRhs<ConstrainingEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u);
			break;
		case 2:
			this->template AssembleRhs<Triangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u);
			this->template AssembleRhs<Quadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->template AssembleRhs<ConstrainingTriangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u);
			this->template AssembleRhs<ConstrainingQuadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u);
			break;
		case 3:
			this->template AssembleRhs<Tetrahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u);
			this->template AssembleRhs<Pyramid>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u);
			this->template AssembleRhs<Prism>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u);
			this->template AssembleRhs<Hexahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, u);
			this->template AssembleRhs<Octahedron>
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
				m_vConstraint[i]->adjust_rhs(rhs, u, dd);
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
			(vElemDisc, m_spApproxSpace->domain(), dd, dd->template begin<TElem>(si), dd->template end<TElem>(si), si,
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
	std::vector<int> vType(2);
	vType[0] = CT_DIRICHLET;
	vType[1] = CT_CONSTRAINTS;

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
				m_vConstraint[i]->adjust_solution(u, dd);
			}
	}

	} UG_CATCH_THROW("Cannot adjust solution.");
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
	vector_type* u_vtk
)
{
	PROFILE_FUNC_GROUP("error_estimator");

//	get multigrid
	SmartPtr<MultiGrid> pMG = ((DoFDistribution *) dd.get())->multi_grid();

//	update the elem discs
	update_disc_items();

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	try
	{
		CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, dd->subset_handler());
	}
	UG_CATCH_THROW("'DomainDiscretization': Can not create Subset Groups and Union.");

//	get the error estimator data for all the discretizations
	std::vector<IErrEstData<TDomain>*> vErrEstData;
	for(std::size_t i = 0; i < m_vElemDisc.size(); ++i)
	{
		SmartPtr<IErrEstData<TDomain> > sp_err_est_data = m_vElemDisc[i]->err_est_data();
		IErrEstData<TDomain>* err_est_data = sp_err_est_data.get();
		if (err_est_data == NULL) continue; // no data specified
		if (std::find (vErrEstData.begin(), vErrEstData.end(), err_est_data) != vErrEstData.end())
			continue; // this one is already in the array
		if (err_est_data->consider_me()) vErrEstData.push_back(err_est_data);
	}

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
		std::vector<IElemDisc<TDomain>*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	remove from elemDisc list those with !err_est_enabled()
		typename std::vector<IElemDisc<TDomain>*>::iterator it = vSubsetElemDisc.begin();
		while (it != vSubsetElemDisc.end())
		{
			if (!(*it)->err_est_enabled())
				it = vSubsetElemDisc.erase(it);
			else ++it;
		}

	//	assemble on suitable elements
		try
		{
			switch (dim)
			{
			case 1:
				this->template AssembleErrorEstimator<RegularEdge>
					(vSubsetElemDisc, dd, si, bNonRegularGrid, u);
				// When assembling over lower-dim manifolds that contain hanging nodes:
				this->template AssembleErrorEstimator<ConstrainingEdge>
					(vSubsetElemDisc, dd, si, bNonRegularGrid, u);
				break;
			case 2:
				this->template AssembleErrorEstimator<Triangle>
					(vSubsetElemDisc, dd, si, bNonRegularGrid, u);
				this->template AssembleErrorEstimator<Quadrilateral>
					(vSubsetElemDisc, dd, si, bNonRegularGrid, u);
				// When assembling over lower-dim manifolds that contain hanging nodes:
				this->template AssembleErrorEstimator<ConstrainingTriangle>
					(vSubsetElemDisc, dd, si, bNonRegularGrid, u);
				this->template AssembleErrorEstimator<ConstrainingQuadrilateral>
					(vSubsetElemDisc, dd, si, bNonRegularGrid, u);
				break;
			case 3:
				this->template AssembleErrorEstimator<Tetrahedron>
					(vSubsetElemDisc, dd, si, bNonRegularGrid, u);
				this->template AssembleErrorEstimator<Pyramid>
					(vSubsetElemDisc, dd, si, bNonRegularGrid, u);
				this->template AssembleErrorEstimator<Prism>
					(vSubsetElemDisc, dd, si, bNonRegularGrid, u);
				this->template AssembleErrorEstimator<Hexahedron>
					(vSubsetElemDisc, dd, si, bNonRegularGrid, u);
				this->template AssembleErrorEstimator<Octahedron>
					(vSubsetElemDisc, dd, si, bNonRegularGrid, u);
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
					m_vConstraint[i]->adjust_error(u, dd);
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
	typedef typename domain_traits<dim>::element_type elem_type;
	typedef typename SurfaceView::traits<elem_type>::const_iterator elem_iter_type;

	// default value negative in order to distinguish between newly added elements (e.g. after refinement)
	// and elements which an error indicator is known for
	pMG->template attach_to_dv<elem_type>(m_aError, -1.0);
	m_pMG = pMG;
	m_aaError = aa_type(*pMG, m_aError);

	// loop surface elements
	ConstSmartPtr<SurfaceView> sv = dd->surface_view();
	const GridLevel& gl = dd->grid_level();
	elem_iter_type elem_iter_end = sv->template end<elem_type> (gl, SurfaceView::ALL);
	for (elem_iter_type elem = sv->template begin<elem_type> (gl, SurfaceView::ALL); elem != elem_iter_end; ++elem)
	{
		// clear attachment (to be on the safe side)
		m_aaError[*elem] = 0.0;

		// get corner coordinates
		std::vector<MathVector<dim> > vCornerCoords = std::vector<MathVector<dim> >(0);
		CollectCornerCoordinates(vCornerCoords, *elem, m_spApproxSpace->domain()->position_accessor(), false);

		// integrate for all estimators, then add up
		for (std::size_t ee = 0; ee < vErrEstData.size(); ++ee)
			m_aaError[*elem] += vErrEstData[ee]->scaling_factor()
								* vErrEstData[ee]->get_elem_error_indicator(*elem, &vCornerCoords[0]);
	}

//	write error estimator values to vtk
	if (u_vtk)
	{
		// local indices and solution
		LocalIndices ind; LocalVector locU;

		// cast u_vtk to grid_function
		GridFunction<TDomain,TAlgebra>* uVTK = dynamic_cast<GridFunction<TDomain,TAlgebra>*>(u_vtk);
		if (!uVTK)
		{
			UG_THROW("Argument passed as output for error function is not a GridFunction.");
		}

		// clear previous values
		uVTK->set(0.0);

		// map attachments to grid function
		ConstSmartPtr<SurfaceView> sv = uVTK->approx_space()->dof_distribution(gl)->surface_view();
		elem_iter_type elem_iter_end = sv->template end<elem_type> (gl, SurfaceView::ALL);
		for (elem_iter_type elem = sv->template begin<elem_type> (gl, SurfaceView::ALL); elem != elem_iter_end; ++elem)
		{
			// 	get global indices
			uVTK->approx_space()->dof_distribution(gl)->indices(*elem, ind, false);

			// 	adapt local algebra
			locU.resize(ind);

			// 	read local values of u
			GetLocalVector(locU, *uVTK);

			// assign error value
			locU(0,0) = m_aaError[*elem];

			// add to grid function
			AddLocalVector(*uVTK, locU);
		}
	}

	m_bErrorCalculated = true;

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
AssembleErrorEstimator(	const std::vector<IElemDisc<domain_type>*>& vElemDisc,
						ConstSmartPtr<DoFDistribution> dd,
						int si, bool bNonRegularGrid,
						const vector_type& u)
{
	//	general case: assembling over all elements in subset si
	gass_type::template AssembleErrorEstimator<TElem>
		(vElemDisc, m_spApproxSpace->domain(), dd,
			dd->template begin<TElem>(si), dd->template end<TElem>(si),
				si, bNonRegularGrid, u);
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
		gass_type::PrepareTimestep(m_vElemDisc, dd, bNonRegularGrid, vSol, m_spAssTuner);
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
		case 1:
			this->template PrepareTimestepElem<RegularEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->template PrepareTimestepElem<ConstrainingEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			break;
		case 2:
			this->template PrepareTimestepElem<Triangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			this->template PrepareTimestepElem<Quadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->template PrepareTimestepElem<ConstrainingTriangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			this->template PrepareTimestepElem<ConstrainingQuadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			break;
		case 3:
			this->template PrepareTimestepElem<Tetrahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			this->template PrepareTimestepElem<Pyramid>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			this->template PrepareTimestepElem<Prism>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			this->template PrepareTimestepElem<Hexahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			this->template PrepareTimestepElem<Octahedron>
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
				dd->template begin<TElem>(si), dd->template end<TElem>(si), si,
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
					m_vConstraint[i]->modify_solution(pModifyMemory, vSol, dd);
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
		case 1:
			this->template AssembleJacobian<RegularEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, pModifyU, s_a0);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->template AssembleJacobian<ConstrainingEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, pModifyU, s_a0);
			break;
		case 2:
			this->template AssembleJacobian<Triangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, pModifyU, s_a0);
			this->template AssembleJacobian<Quadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, pModifyU, s_a0);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->template AssembleJacobian<ConstrainingTriangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, pModifyU, s_a0);
			this->template AssembleJacobian<ConstrainingQuadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, pModifyU, s_a0);
			break;
		case 3:
			this->template AssembleJacobian<Tetrahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, pModifyU, s_a0);
			this->template AssembleJacobian<Pyramid>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, pModifyU, s_a0);
			this->template AssembleJacobian<Prism>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, pModifyU, s_a0);
			this->template AssembleJacobian<Hexahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, J, pModifyU, s_a0);
			this->template AssembleJacobian<Octahedron>
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
				m_vConstraint[i]->adjust_jacobian(J, *pModifyU->solution(0), dd, time, pModifyU,s_a0);
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
				dd->template begin<TElem>(si), dd->template end<TElem>(si), si,
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
					m_vConstraint[i]->modify_solution(pModifyMemory, vSol, dd);
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
		case 1:
			this->template AssembleDefect<RegularEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->template AssembleDefect<ConstrainingEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff);
			break;
		case 2:
			this->template AssembleDefect<Triangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff);
			this->template AssembleDefect<Quadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->template AssembleDefect<ConstrainingTriangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff);
			this->template AssembleDefect<ConstrainingQuadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff);
			break;
		case 3:
			this->template AssembleDefect<Tetrahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff);
			this->template AssembleDefect<Pyramid>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff);
			this->template AssembleDefect<Prism>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff);
			this->template AssembleDefect<Hexahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff);
			this->template AssembleDefect<Octahedron>
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
				m_vConstraint[i]->adjust_defect(d, *pModifyU->solution(0), dd, pModifyU->time(0), pModifyU, &vScaleMass, &vScaleStiff);
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
				dd->template begin<TElem>(si), dd->template end<TElem>(si), si,
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
		case 1:
			this->template AssembleLinear<RegularEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->template AssembleLinear<ConstrainingEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff);
			break;
		case 2:
			this->template AssembleLinear<Triangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff);
			this->template AssembleLinear<Quadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->template AssembleLinear<ConstrainingTriangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff);
			this->template AssembleLinear<ConstrainingQuadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff);
			break;
		case 3:
			this->template AssembleLinear<Tetrahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff);
			this->template AssembleLinear<Pyramid>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff);
			this->template AssembleLinear<Prism>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff);
			this->template AssembleLinear<Hexahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff);
			this->template AssembleLinear<Octahedron>
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
				m_vConstraint[i]->adjust_linear(mat, rhs, dd, vSol->time(0));
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
				dd->template begin<TElem>(si), dd->template end<TElem>(si), si,
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
		case 1:
			this->template AssembleRhs<RegularEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->template AssembleRhs<ConstrainingEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff);
			break;
		case 2:
			this->template AssembleRhs<Triangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff);
			this->template AssembleRhs<Quadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->template AssembleRhs<ConstrainingTriangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff);
			this->template AssembleRhs<ConstrainingQuadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff);
			break;
		case 3:
			this->template AssembleRhs<Tetrahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff);
			this->template AssembleRhs<Pyramid>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff);
			this->template AssembleRhs<Prism>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff);
			this->template AssembleRhs<Hexahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff);
			this->template AssembleRhs<Octahedron>
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
				m_vConstraint[i]->adjust_rhs(rhs, rhs, dd, vSol->time(0));
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
				dd->template begin<TElem>(si), dd->template end<TElem>(si), si,
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
	std::vector<int> vType(2);
	vType[0] = CT_DIRICHLET;
	vType[1] = CT_CONSTRAINTS;

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
				m_vConstraint[i]->adjust_solution(u, dd, time);
			}
	}
	} UG_CATCH_THROW(" Cannot adjust solution.");
}


///////////////////////////////////////////////////////////////////////////////
// Error estimator (instationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
calc_error(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		   ConstSmartPtr<DoFDistribution> dd,
		   std::vector<number> vScaleMass,
		   std::vector<number> vScaleStiff,
		   vector_type* u_vtk)
{
	PROFILE_FUNC_GROUP("error_estimator");

//	get multigrid
	SmartPtr<MultiGrid> pMG = ((DoFDistribution *) dd.get())->multi_grid();

//	update the elem discs
	update_disc_items();

//	Union of Subsets
	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	create list of all subsets
	try
	{
		CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, dd->subset_handler());
	}
	UG_CATCH_THROW("'DomainDiscretization': Can not create Subset Groups and Union.");

//	get the error estimator data for all the discretizations
	std::vector<IErrEstData<TDomain>*> vErrEstData;
	for (std::size_t i = 0; i < m_vElemDisc.size(); ++i)
	{
		SmartPtr<IErrEstData<TDomain> > sp_err_est_data = m_vElemDisc[i]->err_est_data();
		IErrEstData<TDomain>* err_est_data = sp_err_est_data.get();
		if (err_est_data == NULL) continue; // no data specified
		if (std::find (vErrEstData.begin(), vErrEstData.end(), err_est_data) != vErrEstData.end())
			continue; // this one is already in the array
		if (err_est_data->consider_me()) vErrEstData.push_back(err_est_data);
	}

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

	//	Elem Disc on the subset
		std::vector<IElemDisc<TDomain>*> vSubsetElemDisc;

	//	get all element discretizations that work on the subset
		GetElemDiscOnSubset(vSubsetElemDisc, m_vElemDisc, vSSGrp, si);

	//	remove from elemDisc list those with !err_est_enabled()
		typename std::vector<IElemDisc<TDomain>*>::iterator it = vSubsetElemDisc.begin();
		while (it != vSubsetElemDisc.end())
		{
			if (!(*it)->err_est_enabled())
				it = vSubsetElemDisc.erase(it);
			else ++it;
		}

	//	assemble on suitable elements
		try
		{
			switch (dim)
			{
			case 1:
				this->template AssembleErrorEstimator<RegularEdge>
					(vSubsetElemDisc, dd, si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
				// When assembling over lower-dim manifolds that contain hanging nodes:
				this->template AssembleErrorEstimator<ConstrainingEdge>
					(vSubsetElemDisc, dd, si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
				break;
			case 2:
				this->template AssembleErrorEstimator<Triangle>
					(vSubsetElemDisc, dd, si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
				this->template AssembleErrorEstimator<Quadrilateral>
					(vSubsetElemDisc, dd, si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
				// When assembling over lower-dim manifolds that contain hanging nodes:
				this->template AssembleErrorEstimator<ConstrainingTriangle>
					(vSubsetElemDisc, dd, si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
				this->template AssembleErrorEstimator<ConstrainingQuadrilateral>
					(vSubsetElemDisc, dd, si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
				break;
			case 3:
				this->template AssembleErrorEstimator<Tetrahedron>
					(vSubsetElemDisc, dd, si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
				this->template AssembleErrorEstimator<Pyramid>
					(vSubsetElemDisc, dd, si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
				this->template AssembleErrorEstimator<Prism>
					(vSubsetElemDisc, dd, si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
				this->template AssembleErrorEstimator<Hexahedron>
					(vSubsetElemDisc, dd, si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
				this->template AssembleErrorEstimator<Octahedron>
					(vSubsetElemDisc, dd, si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
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
					m_vConstraint[i]->adjust_error(*vSol->solution(0), dd, vSol->time(0));
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
	typedef typename domain_traits<dim>::element_type elem_type;
	typedef typename SurfaceView::traits<elem_type>::const_iterator elem_iter_type;

	// default value negative in order to distinguish between newly added elements (e.g. after refinement)
	// and elements which an error indicator is known for
	pMG->template attach_to_dv<elem_type>(m_aError, -1.0);
	m_pMG = pMG;
	m_aaError = aa_type(*pMG, m_aError);

	// loop surface elements
	ConstSmartPtr<SurfaceView> sv = dd->surface_view();
	const GridLevel& gl = dd->grid_level();
	elem_iter_type elem_iter_end = sv->template end<elem_type> (gl, SurfaceView::ALL);
	for (elem_iter_type elem = sv->template begin<elem_type> (gl, SurfaceView::ALL); elem != elem_iter_end; ++elem)
	{
		// clear attachment
		m_aaError[*elem] = 0.0;

		// get corner coordinates
		std::vector<MathVector<dim> > vCornerCoords = std::vector<MathVector<dim> >(0);
		CollectCornerCoordinates(vCornerCoords, *elem, m_spApproxSpace->domain()->position_accessor(), false);

		// integrate for all estimators, then add up
		for (std::size_t ee = 0; ee < vErrEstData.size(); ++ee)
			m_aaError[*elem] += vErrEstData[ee]->scaling_factor()
								* vErrEstData[ee]->get_elem_error_indicator(*elem, &vCornerCoords[0]);
	}

//	write error estimator values to vtk
	if (u_vtk)
	{
		// local indices and solution
		LocalIndices ind; LocalVector locU;

		// cast u_vtk to grid_function
		GridFunction<TDomain,TAlgebra>* uVTK = dynamic_cast<GridFunction<TDomain,TAlgebra>*>(u_vtk);
		if (!uVTK)
		{
			UG_THROW("Argument passed as output for error function is not a GridFunction.");
		}

		// clear previous values
		uVTK->set(0.0);

		// map attachments to grid function
		ConstSmartPtr<SurfaceView> sv = uVTK->approx_space()->dof_distribution(gl)->surface_view();
		elem_iter_type elem_iter_end = sv->template end<elem_type> (gl, SurfaceView::ALL);
		for (elem_iter_type elem = sv->template begin<elem_type> (gl, SurfaceView::ALL); elem != elem_iter_end; ++elem)
		{
			// 	get global indices
			uVTK->approx_space()->dof_distribution(gl)->indices(*elem, ind, false);

			// 	adapt local algebra
			locU.resize(ind);

			// 	read local values of u
			GetLocalVector(locU, *uVTK);

			// assign error value
			locU(0,0) = m_aaError[*elem];

			// add to grid function
			AddLocalVector(*uVTK, locU);
		}
	}

	m_bErrorCalculated = true;

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
AssembleErrorEstimator(	const std::vector<IElemDisc<domain_type>*>& vElemDisc,
						ConstSmartPtr<DoFDistribution> dd,
						int si, bool bNonRegularGrid,
						std::vector<number> vScaleMass,
						std::vector<number> vScaleStiff,
						ConstSmartPtr<VectorTimeSeries<vector_type> > vSol)
{
	//	general case: assembling over all elements in subset si
	gass_type::template AssembleErrorEstimator<TElem>
		(vElemDisc, m_spApproxSpace->domain(), dd,
			dd->template begin<TElem>(si), dd->template end<TElem>(si),
				si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
}

template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
mark_for_refinement
(	IRefiner& refiner,
	number TOL,
	number refineFrac,
	int maxLevel
)
{
	// check that error indicators have been calculated
	if (!m_bErrorCalculated)
	{
		UG_THROW("Error indicators have to be calculated first by a call to 'calc_error'.");
	}

	// mark elements for refinement
	MarkElementsForRefinement<elem_type>(m_aaError, refiner,
		this->dd(GridLevel(GridLevel::TOP, GridLevel::SURFACE)),
		TOL, refineFrac, maxLevel);
}

template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
mark_for_coarsening
(	IRefiner& refiner,
	number TOL,
	number coarseFrac,
	int maxLevel
)
{
	// check that error indicators have been calculated
	if (!m_bErrorCalculated)
	{
		UG_THROW("Error indicators have to be calculated first by a call to 'calc_error'.");
	}

	// mark elements for coarsening
	MarkElementsForCoarsening<elem_type>(m_aaError, refiner,
		this->dd(GridLevel(GridLevel::TOP, GridLevel::SURFACE)),
		TOL, coarseFrac, maxLevel);
}

template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
mark_with_strategy
(	IRefiner& refiner,
	SmartPtr<IElementMarkingStrategy<TDomain> >spMarkingStrategy
)
{
	// check that error indicators have been calculated
	if (!m_bErrorCalculated)
	{
		UG_THROW("Error indicators have to be calculated first by a call to 'calc_error'.");
	}

	// mark elements for refinement
	if (spMarkingStrategy.valid())
	{
		spMarkingStrategy->mark(m_aaError, refiner, this->dd(GridLevel(GridLevel::TOP, GridLevel::SURFACE)));
	}
}

template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
void DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
invalidate_error()
{
	// check that error indicators have been calculated
	if (m_bErrorCalculated)
	{
		m_bErrorCalculated = false;
		m_pMG->template detach_from<elem_type>(m_aError);
	}
}

template <typename TDomain, typename TAlgebra, typename TGlobAssembler>
bool DomainDiscretizationBase<TDomain, TAlgebra, TGlobAssembler>::
is_error_valid()
{
	// check that error indicators have been calculated
	return m_bErrorCalculated;
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

//	call assembler's PrepareTimestep
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
		case 1:
			this->template FinishTimestepElem<RegularEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->template FinishTimestepElem<ConstrainingEdge>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			break;
		case 2:
			this->template FinishTimestepElem<Triangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			this->template FinishTimestepElem<Quadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			this->template FinishTimestepElem<ConstrainingTriangle>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			this->template FinishTimestepElem<ConstrainingQuadrilateral>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			break;
		case 3:
			this->template FinishTimestepElem<Tetrahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			this->template FinishTimestepElem<Pyramid>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			this->template FinishTimestepElem<Prism>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			this->template FinishTimestepElem<Hexahedron>
				(vSubsetElemDisc, dd, si, bNonRegularGrid, vSol);
			this->template FinishTimestepElem<Octahedron>
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
				dd->template begin<TElem>(si), dd->template end<TElem>(si), si,
					bNonRegularGrid, vSol, m_spAssTuner);
	}
}

} // end namespace ug

#endif /*__H__UG__LIB_DISC__SPATIAL_DISC__DOMAIN_DISC_IMPL__*/
