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
#include "lib_disc/spatial_disc/elem_disc/elem_disc_assemble_util.h"
#include "lib_disc/function_spaces/error_indicator_util.h"
#ifdef UG_PARALLEL
#include "lib_disc/parallelization/parallelization_util.h"
#endif

namespace ug{

template <typename TDomain, typename TAlgebra>
void DomainDiscretization<TDomain, TAlgebra>::update_elem_discs()
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

template <typename TDomain, typename TAlgebra>
void DomainDiscretization<TDomain, TAlgebra>::update_constraints()
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

template <typename TDomain, typename TAlgebra>
void DomainDiscretization<TDomain, TAlgebra>::update_disc_items()
{
	update_elem_discs();
	update_constraints();
}

///////////////////////////////////////////////////////////////////////////////
// Mass Matrix
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra>
void DomainDiscretization<TDomain, TAlgebra>::
assemble_mass_matrix(matrix_type& M, const vector_type& u,
                     ConstSmartPtr<DoFDistribution> dd)
{
	PROFILE_FUNC_GROUP("discretization");
//	update the elem discs
	update_disc_items();

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
			AssembleMassMatrix<RegularEdge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, M, u, m_spAssTuner);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			AssembleMassMatrix<ConstrainingEdge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, M, u, m_spAssTuner);
			break;
		case 2:
			AssembleMassMatrix<Triangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, M, u, m_spAssTuner);
			AssembleMassMatrix<Quadrilateral,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, M, u, m_spAssTuner);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			AssembleMassMatrix<ConstrainingTriangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, M, u, m_spAssTuner);
			AssembleMassMatrix<ConstrainingQuadrilateral,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, M, u, m_spAssTuner);
			break;
		case 3:
			AssembleMassMatrix<Tetrahedron,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, M, u, m_spAssTuner);
			AssembleMassMatrix<Pyramid,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, M, u, m_spAssTuner);
			AssembleMassMatrix<Prism,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, M, u, m_spAssTuner);
			AssembleMassMatrix<Hexahedron,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, M, u, m_spAssTuner);
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
	}UG_CATCH_THROW("DomainDiscretization::assemble_mass_matrix:"
					" Cannot execute post process.");

//	Remember parallel storage type
#ifdef UG_PARALLEL
	M.set_storage_type(PST_ADDITIVE);
	M.set_layouts(dd->layouts());
#endif
}

///////////////////////////////////////////////////////////////////////////////
// Stiffness Matrix
///////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
void DomainDiscretization<TDomain, TAlgebra>::
assemble_stiffness_matrix(matrix_type& A, const vector_type& u,
                          ConstSmartPtr<DoFDistribution> dd)
{
	PROFILE_FUNC_GROUP("discretization");
//	update the elem discs
	update_disc_items();

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
			AssembleStiffnessMatrix<RegularEdge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, A, u, m_spAssTuner);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			AssembleStiffnessMatrix<ConstrainingEdge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, A, u, m_spAssTuner);
			break;
		case 2:
			AssembleStiffnessMatrix<Triangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, A, u, m_spAssTuner);
			AssembleStiffnessMatrix<Quadrilateral,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, A, u, m_spAssTuner);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			AssembleStiffnessMatrix<ConstrainingTriangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, A, u, m_spAssTuner);
			AssembleStiffnessMatrix<ConstrainingQuadrilateral,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, A, u, m_spAssTuner);
			break;
		case 3:
			AssembleStiffnessMatrix<Tetrahedron,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, A, u, m_spAssTuner);
			AssembleStiffnessMatrix<Pyramid,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, A, u, m_spAssTuner);
			AssembleStiffnessMatrix<Prism,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, A, u, m_spAssTuner);
			AssembleStiffnessMatrix<Hexahedron,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, A, u, m_spAssTuner);
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
	}UG_CATCH_THROW("DomainDiscretization::assemble_stiffness_matrix:"
					" Cannot execute post process.");

//	Remember parallel storage type
#ifdef UG_PARALLEL
	A.set_storage_type(PST_ADDITIVE);
	A.set_layouts(dd->layouts());
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
void DomainDiscretization<TDomain, TAlgebra>::
assemble_jacobian(matrix_type& J,
                  const vector_type& u,
                  ConstSmartPtr<DoFDistribution> dd)
{
	PROFILE_FUNC_GROUP("discretization");
//	update the elem discs
	update_disc_items();

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
			AssembleJacobian<RegularEdge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, J, *pModifyU, m_spAssTuner);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			AssembleJacobian<ConstrainingEdge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, J, *pModifyU, m_spAssTuner);
			break;
		case 2:
			AssembleJacobian<Triangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, J, *pModifyU, m_spAssTuner);
			AssembleJacobian<Quadrilateral,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, J, *pModifyU, m_spAssTuner);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			AssembleJacobian<ConstrainingTriangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, J, *pModifyU, m_spAssTuner);
			AssembleJacobian<ConstrainingQuadrilateral,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, J, *pModifyU, m_spAssTuner);
			break;
		case 3:
			AssembleJacobian<Tetrahedron,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, J, *pModifyU, m_spAssTuner);
			AssembleJacobian<Pyramid,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, J, *pModifyU, m_spAssTuner);
			AssembleJacobian<Prism,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, J, *pModifyU, m_spAssTuner);
			AssembleJacobian<Hexahedron,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, J, *pModifyU, m_spAssTuner);
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
	}UG_CATCH_THROW("DomainDiscretization::assemble_jacobian:"
					" Cannot execute post process.");

//	Remember parallel storage type
#ifdef UG_PARALLEL
	J.set_storage_type(PST_ADDITIVE);
	J.set_layouts(dd->layouts());
#endif
}


///////////////////////////////////////////////////////////////////////////////
// Defect (stationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra>
void DomainDiscretization<TDomain, TAlgebra>::
assemble_defect(vector_type& d,
                const vector_type& u,
                ConstSmartPtr<DoFDistribution> dd)
{
	PROFILE_FUNC_GROUP("discretization");
//	update the elem discs
	update_disc_items();

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
			AssembleDefect<RegularEdge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, d, *pModifyU, m_spAssTuner);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			AssembleDefect<ConstrainingEdge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, d, *pModifyU, m_spAssTuner);
			break;
		case 2:
			AssembleDefect<Triangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, d, *pModifyU, m_spAssTuner);
			AssembleDefect<Quadrilateral,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, d, *pModifyU, m_spAssTuner);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			AssembleDefect<ConstrainingTriangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, d, *pModifyU, m_spAssTuner);
			AssembleDefect<ConstrainingQuadrilateral,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, d, *pModifyU, m_spAssTuner);
			break;
		case 3:
			AssembleDefect<Tetrahedron,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, d, *pModifyU, m_spAssTuner);
			AssembleDefect<Pyramid,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, d, *pModifyU, m_spAssTuner);
			AssembleDefect<Prism,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, d, *pModifyU, m_spAssTuner);
			AssembleDefect<Hexahedron,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, d, *pModifyU, m_spAssTuner);
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
void DomainDiscretization<TDomain, TAlgebra>::
assemble_linear(matrix_type& mat, vector_type& rhs,
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
			AssembleLinear<RegularEdge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, mat, rhs, m_spAssTuner);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			AssembleLinear<ConstrainingEdge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, mat, rhs, m_spAssTuner);
			break;
		case 2:
			AssembleLinear<Triangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, mat, rhs, m_spAssTuner);
			AssembleLinear<Quadrilateral,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, mat, rhs, m_spAssTuner);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			AssembleLinear<ConstrainingTriangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, mat, rhs, m_spAssTuner);
			AssembleLinear<ConstrainingQuadrilateral,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, mat, rhs, m_spAssTuner);
			break;
		case 3:
			AssembleLinear<Tetrahedron,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, mat, rhs, m_spAssTuner);
			AssembleLinear<Pyramid,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, mat, rhs, m_spAssTuner);
			AssembleLinear<Prism,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, mat, rhs, m_spAssTuner);
			AssembleLinear<Hexahedron,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, mat, rhs, m_spAssTuner);
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
	}UG_CATCH_THROW("DomainDiscretization::assemble_linear: Cannot post process.");

//	Remember parallel storage type
#ifdef UG_PARALLEL
	mat.set_storage_type(PST_ADDITIVE);
	mat.set_layouts(dd->layouts());
	rhs.set_storage_type(PST_ADDITIVE);
#endif
}

///////////////////////////////////////////////////////////////////////////////
// RHS (stationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra>
void DomainDiscretization<TDomain, TAlgebra>::
assemble_rhs(vector_type& rhs,
			const vector_type& u,
			ConstSmartPtr<DoFDistribution> dd)
{
	PROFILE_FUNC_GROUP("discretization");
//	update the elem discs
	update_disc_items();

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
			AssembleRhs<RegularEdge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, rhs, u, m_spAssTuner);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			AssembleRhs<ConstrainingEdge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, rhs, u, m_spAssTuner);
			break;
		case 2:
			AssembleRhs<Triangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, rhs, u, m_spAssTuner);
			AssembleRhs<Quadrilateral,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, rhs, u, m_spAssTuner);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			AssembleRhs<ConstrainingTriangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, rhs, u, m_spAssTuner);
			AssembleRhs<ConstrainingQuadrilateral,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, rhs, u, m_spAssTuner);
			break;
		case 3:
			AssembleRhs<Tetrahedron,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, rhs, u, m_spAssTuner);
			AssembleRhs<Pyramid,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, rhs, u, m_spAssTuner);
			AssembleRhs<Prism,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, rhs, u, m_spAssTuner);
			AssembleRhs<Hexahedron,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, rhs, u, m_spAssTuner);
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
	}UG_CATCH_THROW("DomainDiscretization::assemble_rhs:"
					" Cannot execute post process.");

//	Remember parallel storage type
#ifdef UG_PARALLEL
	rhs.set_storage_type(PST_ADDITIVE);
#endif
}

template <typename TDomain, typename TAlgebra>
void DomainDiscretization<TDomain, TAlgebra>::
assemble_rhs(vector_type& rhs,
			ConstSmartPtr<DoFDistribution> dd)
{
	assemble_rhs(rhs, rhs, dd);
}

///////////////////////////////////////////////////////////////////////////////
// set constraints (stationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra>
void DomainDiscretization<TDomain, TAlgebra>::
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
template <typename TDomain, typename TAlgebra>
void DomainDiscretization<TDomain, TAlgebra>::
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
	try{
		CreateSubsetGroups(vSSGrp, unionSubsets, m_vElemDisc, dd->subset_handler());
	}UG_CATCH_THROW("'DomainDiscretization': Can not create Subset Groups and Union.");

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
	try{
		for(size_t i = 0; i < vErrEstData.size(); ++i)
			vErrEstData[i]->alloc_err_est_data(dd->surface_view(), dd->grid_level());
	}
	UG_CATCH_THROW("DomainDiscretization::calc_error: Cannot prepare the error estimator");

//	loop subsets to assemble the estimators
	for(size_t i = 0; i < unionSubsets.size(); ++i)
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

	//	assemble on suitable elements
		try
		{
		switch(dim)
		{
		case 1:
			AssembleErrorEstimator<RegularEdge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, u);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			AssembleErrorEstimator<ConstrainingEdge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, u);
			break;
		case 2:
			AssembleErrorEstimator<Triangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, u);
			AssembleErrorEstimator<Quadrilateral,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, u);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			AssembleErrorEstimator<ConstrainingTriangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, u);
			AssembleErrorEstimator<ConstrainingQuadrilateral,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, u);
			break;
		case 3:
			AssembleErrorEstimator<Tetrahedron,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, u);
			AssembleErrorEstimator<Pyramid,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, u);
			AssembleErrorEstimator<Prism,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, u);
			AssembleErrorEstimator<Hexahedron,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, u);
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
			m_aaError[*elem] += vErrEstData[ee]->get_elem_error_indicator(*elem, &vCornerCoords[0]);
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

//	postprocess the error estimators in the discretizations
	try{
		for(std::size_t i = 0; i < vErrEstData.size(); ++i)
			vErrEstData[i]->release_err_est_data();
	}
	UG_CATCH_THROW("DomainDiscretization::calc_error: Cannot release the error estimator");
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
void DomainDiscretization<TDomain, TAlgebra>::
prepare_timestep(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
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
			PrepareTimestep<RegularEdge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vSol, m_spAssTuner);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			PrepareTimestep<ConstrainingEdge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vSol, m_spAssTuner);
			break;
		case 2:
			PrepareTimestep<Triangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vSol, m_spAssTuner);
			PrepareTimestep<Quadrilateral,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vSol, m_spAssTuner);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			PrepareTimestep<ConstrainingTriangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vSol, m_spAssTuner);
			PrepareTimestep<ConstrainingQuadrilateral,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vSol, m_spAssTuner);
			break;
		case 3:
			PrepareTimestep<Tetrahedron,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vSol, m_spAssTuner);
			PrepareTimestep<Pyramid,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vSol, m_spAssTuner);
			PrepareTimestep<Prism,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vSol, m_spAssTuner);
			PrepareTimestep<Hexahedron,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vSol, m_spAssTuner);
			break;
		default:
			UG_THROW("DomainDiscretization::prepare_timestep (instationary):"
							"Dimension "<<dim<<" (subset="<<si<<") not supported.");
		}
		}
		UG_CATCH_THROW("DomainDiscretization::prepare_timestep (instationary):"
						" Assembling of elements of Dimension " << dim << " in "
						" subset "<<si<< " failed.");
	}
}


///////////////////////////////////////////////////////////////////////////////
// Jacobian (instationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra>
void DomainDiscretization<TDomain, TAlgebra>::
assemble_jacobian(matrix_type& J,
                  ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
                  const number s_a0,
                  ConstSmartPtr<DoFDistribution> dd)
{
	PROFILE_FUNC_GROUP("discretization");
//	update the elem discs
	update_disc_items();

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
			AssembleJacobian<RegularEdge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, J, pModifyU, s_a0, m_spAssTuner);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			AssembleJacobian<ConstrainingEdge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, J, pModifyU, s_a0, m_spAssTuner);
			break;
		case 2:
			AssembleJacobian<Triangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, J, pModifyU, s_a0, m_spAssTuner);
			AssembleJacobian<Quadrilateral,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, J, pModifyU, s_a0, m_spAssTuner);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			AssembleJacobian<ConstrainingTriangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, J, pModifyU, s_a0, m_spAssTuner);
			AssembleJacobian<ConstrainingQuadrilateral,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, J, pModifyU, s_a0, m_spAssTuner);
			break;
		case 3:
			AssembleJacobian<Tetrahedron,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, J, pModifyU, s_a0, m_spAssTuner);
			AssembleJacobian<Pyramid,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, J, pModifyU, s_a0, m_spAssTuner);
			AssembleJacobian<Prism,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, J, pModifyU, s_a0, m_spAssTuner);
			AssembleJacobian<Hexahedron,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, J, pModifyU, s_a0, m_spAssTuner);
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
	}UG_CATCH_THROW("Cannot adjust jacobian.");

//	Remember parallel storage type
#ifdef UG_PARALLEL
	J.set_storage_type(PST_ADDITIVE);
	J.set_layouts(dd->layouts());
#endif
}


///////////////////////////////////////////////////////////////////////////////
// Defect (instationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra>
void DomainDiscretization<TDomain, TAlgebra>::
assemble_defect(vector_type& d,
                ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
                const std::vector<number>& vScaleMass,
                const std::vector<number>& vScaleStiff,
                ConstSmartPtr<DoFDistribution> dd)
{
	PROFILE_FUNC_GROUP("discretization");
//	update the elem discs
	update_disc_items();

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
			AssembleDefect<RegularEdge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff, m_spAssTuner);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			AssembleDefect<ConstrainingEdge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff, m_spAssTuner);
			break;
		case 2:
			AssembleDefect<Triangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff, m_spAssTuner);
			AssembleDefect<Quadrilateral,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff, m_spAssTuner);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			AssembleDefect<ConstrainingTriangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff, m_spAssTuner);
			AssembleDefect<ConstrainingQuadrilateral,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff, m_spAssTuner);
			break;
		case 3:
			AssembleDefect<Tetrahedron,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff, m_spAssTuner);
			AssembleDefect<Pyramid,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff, m_spAssTuner);
			AssembleDefect<Prism,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff, m_spAssTuner);
			AssembleDefect<Hexahedron,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff, m_spAssTuner);
			break;
		default:
			UG_THROW("DomainDiscretization::assemble_defect (instationary):"
							"Dimension "<<dim<<" (subset="<<si<<") not supported.");
		}
		}
		UG_CATCH_THROW("DomainDiscretization::assemble_defect (instationary):"
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
				m_vConstraint[i]->adjust_defect(d, *pModifyU->solution(0), dd, pModifyU->time(0), pModifyU, &vScaleMass, &vScaleStiff);
			}
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
void DomainDiscretization<TDomain, TAlgebra>::
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
			AssembleLinear<RegularEdge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff, m_spAssTuner);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			AssembleLinear<ConstrainingEdge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff, m_spAssTuner);
			break;
		case 2:
			AssembleLinear<Triangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff, m_spAssTuner);
			AssembleLinear<Quadrilateral,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff, m_spAssTuner);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			AssembleLinear<ConstrainingTriangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff, m_spAssTuner);
			AssembleLinear<ConstrainingQuadrilateral,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff, m_spAssTuner);
			break;
		case 3:
			AssembleLinear<Tetrahedron,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff, m_spAssTuner);
			AssembleLinear<Pyramid,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff, m_spAssTuner);
			AssembleLinear<Prism,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff, m_spAssTuner);
			AssembleLinear<Hexahedron,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff, m_spAssTuner);
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

///////////////////////////////////////////////////////////////////////////////
// RHS (instationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra>
void DomainDiscretization<TDomain, TAlgebra>::
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
			AssembleRhs<RegularEdge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff, m_spAssTuner);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			AssembleRhs<ConstrainingEdge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff, m_spAssTuner);
			break;
		case 2:
			AssembleRhs<Triangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff, m_spAssTuner);
			AssembleRhs<Quadrilateral,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff, m_spAssTuner);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			AssembleRhs<ConstrainingTriangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff, m_spAssTuner);
			AssembleRhs<ConstrainingQuadrilateral,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff, m_spAssTuner);
			break;
		case 3:
			AssembleRhs<Tetrahedron,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff, m_spAssTuner);
			AssembleRhs<Pyramid,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff, m_spAssTuner);
			AssembleRhs<Prism,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff, m_spAssTuner);
			AssembleRhs<Hexahedron,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff, m_spAssTuner);
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

///////////////////////////////////////////////////////////////////////////////
// set constraint values (instationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra>
void DomainDiscretization<TDomain, TAlgebra>::
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
template <typename TDomain, typename TAlgebra>
void DomainDiscretization<TDomain, TAlgebra>::
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

	//	assemble on suitable elements
		try
		{
			switch(dim)
			{
			case 1:
				AssembleErrorEstimator<RegularEdge,TDomain,TAlgebra>
					(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
				// When assembling over lower-dim manifolds that contain hanging nodes:
				AssembleErrorEstimator<ConstrainingEdge,TDomain,TAlgebra>
					(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
				break;
			case 2:
				AssembleErrorEstimator<Triangle,TDomain,TAlgebra>
					(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
				AssembleErrorEstimator<Quadrilateral,TDomain,TAlgebra>
					(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
				// When assembling over lower-dim manifolds that contain hanging nodes:
				AssembleErrorEstimator<ConstrainingTriangle,TDomain,TAlgebra>
					(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
				AssembleErrorEstimator<ConstrainingQuadrilateral,TDomain,TAlgebra>
					(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
				break;
			case 3:
				AssembleErrorEstimator<Tetrahedron,TDomain,TAlgebra>
					(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
				AssembleErrorEstimator<Pyramid,TDomain,TAlgebra>
					(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
				AssembleErrorEstimator<Prism,TDomain,TAlgebra>
					(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
				AssembleErrorEstimator<Hexahedron,TDomain,TAlgebra>
					(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vScaleMass, vScaleStiff, vSol);
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
			m_aaError[*elem] += vErrEstData[ee]->get_elem_error_indicator(*elem, &vCornerCoords[0]);
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

template <typename TDomain, typename TAlgebra>
void DomainDiscretization<TDomain, TAlgebra>::
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

template <typename TDomain, typename TAlgebra>
void DomainDiscretization<TDomain, TAlgebra>::
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

template <typename TDomain, typename TAlgebra>
void DomainDiscretization<TDomain, TAlgebra>::
invalidate_error()
{
	// check that error indicators have been calculated
	if (m_bErrorCalculated)
	{
		m_bErrorCalculated = false;
		m_pMG->template detach_from<elem_type>(m_aError);
	}
}

template <typename TDomain, typename TAlgebra>
bool DomainDiscretization<TDomain, TAlgebra>::
is_error_valid()
{
	// check that error indicators have been calculated
	return m_bErrorCalculated;
}




///////////////////////////////////////////////////////////////////////////////
// Finish Timestep (instationary)
///////////////////////////////////////////////////////////////////////////////
template <typename TDomain, typename TAlgebra>
void DomainDiscretization<TDomain, TAlgebra>::
finish_timestep(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
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
			FinishTimestep<RegularEdge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vSol, m_spAssTuner);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			FinishTimestep<ConstrainingEdge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vSol, m_spAssTuner);
			break;
		case 2:
			FinishTimestep<Triangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vSol, m_spAssTuner);
			FinishTimestep<Quadrilateral,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vSol, m_spAssTuner);
			// When assembling over lower-dim manifolds that contain hanging nodes:
			FinishTimestep<ConstrainingTriangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vSol, m_spAssTuner);
			FinishTimestep<ConstrainingQuadrilateral,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vSol, m_spAssTuner);
			break;
		case 3:
			FinishTimestep<Tetrahedron,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vSol, m_spAssTuner);
			FinishTimestep<Pyramid,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vSol, m_spAssTuner);
			FinishTimestep<Prism,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vSol, m_spAssTuner);
			FinishTimestep<Hexahedron,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vSol, m_spAssTuner);
			break;
		default:
			UG_THROW("DomainDiscretization::finish_timestep (instationary):"
							"Dimension "<<dim<<" (subset="<<si<<") not supported.");
		}
		}
		UG_CATCH_THROW("DomainDiscretization::finish_timestep (instationary):"
						" Assembling of elements of Dimension " << dim << " in "
						" subset "<<si<< " failed.");
	}

}


} // end namespace ug

#endif /*__H__UG__LIB_DISC__SPATIAL_DISC__DOMAIN_DISC_IMPL__*/
