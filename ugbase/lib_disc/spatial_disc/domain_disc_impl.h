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
			AssembleMassMatrix<Edge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, M, u, m_spAssTuner);
			break;
		case 2:
			AssembleMassMatrix<Triangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, M, u, m_spAssTuner);
			AssembleMassMatrix<Quadrilateral,TDomain,TAlgebra>
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
			AssembleStiffnessMatrix<Edge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, A, u, m_spAssTuner);
			break;
		case 2:
			AssembleStiffnessMatrix<Triangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, A, u, m_spAssTuner);
			AssembleStiffnessMatrix<Quadrilateral,TDomain,TAlgebra>
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
	SmartPtr<vector_type> pModifyMemory = NULL;
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
			AssembleJacobian<Edge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, J, *pModifyU, m_spAssTuner);
			break;
		case 2:
			AssembleJacobian<Triangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, J, *pModifyU, m_spAssTuner);
			AssembleJacobian<Quadrilateral,TDomain,TAlgebra>
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
	SmartPtr<vector_type> pModifyMemory = NULL;
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
			AssembleDefect<Edge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, d, *pModifyU, m_spAssTuner);
			break;
		case 2:
			AssembleDefect<Triangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, d, *pModifyU, m_spAssTuner);
			AssembleDefect<Quadrilateral,TDomain,TAlgebra>
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
			AssembleLinear<Edge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, mat, rhs, m_spAssTuner);
			break;
		case 2:
			AssembleLinear<Triangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, mat, rhs, m_spAssTuner);
			AssembleLinear<Quadrilateral,TDomain,TAlgebra>
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
			AssembleRhs<Edge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, rhs, u, m_spAssTuner);
			break;
		case 2:
			AssembleRhs<Triangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, rhs, u, m_spAssTuner);
			AssembleRhs<Quadrilateral,TDomain,TAlgebra>
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
	if (m_spAssTuner->single_DoFindex_assembling_enabled()) u.resize(1);

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
			PrepareTimestep<Edge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vSol, m_spAssTuner);
			break;
		case 2:
			PrepareTimestep<Triangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vSol, m_spAssTuner);
			PrepareTimestep<Quadrilateral,TDomain,TAlgebra>
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

//	pre process -  modifies the solution, used for computing the defect
	ConstSmartPtr<VectorTimeSeries<vector_type> > pModifyU = vSol;
	SmartPtr<VectorTimeSeries<vector_type> > pModifyMemory = NULL;
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
			AssembleJacobian<Edge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, J, pModifyU, s_a0, m_spAssTuner);
			break;
		case 2:
			AssembleJacobian<Triangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, J, pModifyU, s_a0, m_spAssTuner);
			AssembleJacobian<Quadrilateral,TDomain,TAlgebra>
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
	SmartPtr<VectorTimeSeries<vector_type> > pModifyMemory = NULL;
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
			AssembleDefect<Edge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff, m_spAssTuner);
			break;
		case 2:
			AssembleDefect<Triangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, d, pModifyU, vScaleMass, vScaleStiff, m_spAssTuner);
			AssembleDefect<Quadrilateral,TDomain,TAlgebra>
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
			AssembleLinear<Edge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff, m_spAssTuner);
			break;
		case 2:
			AssembleLinear<Triangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, mat, rhs, vSol, vScaleMass, vScaleStiff, m_spAssTuner);
			AssembleLinear<Quadrilateral,TDomain,TAlgebra>
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
			AssembleRhs<Edge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff, m_spAssTuner);
			break;
		case 2:
			AssembleRhs<Triangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, rhs, vSol, vScaleMass, vScaleStiff, m_spAssTuner);
			AssembleRhs<Quadrilateral,TDomain,TAlgebra>
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
	if (m_spAssTuner->single_DoFindex_assembling_enabled()) u.resize(1);

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
			FinishTimestep<Edge,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vSol, m_spAssTuner);
			break;
		case 2:
			FinishTimestep<Triangle,TDomain,TAlgebra>
				(vSubsetElemDisc, m_spApproxSpace->domain(), dd, si, bNonRegularGrid, vSol, m_spAssTuner);
			FinishTimestep<Quadrilateral,TDomain,TAlgebra>
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
