/*
 * dirichlet_post_process_util.h
 *
 *  Created on: 04.08.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__POST_PROCESS__DIRICHLET_BOUNDARY__DIRICHLET_POST_PROCESS_UTIL__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__POST_PROCESS__DIRICHLET_BOUNDARY__DIRICHLET_POST_PROCESS_UTIL__

namespace ug{

template <	typename TElem,
			typename TDiscreteFunction,
			typename TAlgebra>
bool
BNDPostProcessJacobian(	IDirichletPostProcess<TAlgebra>& elemDisc,
						typename geometry_traits<TElem>::iterator iterBegin,
						typename geometry_traits<TElem>::iterator iterEnd,
						int si,
						typename TAlgebra::matrix_type& J,
						const TDiscreteFunction& u,
						const FunctionGroup& fcts,
						number time, number s_m, number s_a)
{
	typedef typename reference_element_traits<TElem>::reference_element_type reference_element_type;
	const ReferenceObjectID refID = reference_element_type::REFERENCE_OBJECT_ID;

	// check if at least on element exist, else return
	if(iterBegin == iterEnd) return true;

	LocalIndices ind;

	// set functions
	ind.set_function_group(fcts);

	// prepare local indices for elem type
	if(!u.prepare_inner_indices(refID, si, ind))
		{UG_LOG("ERROR in AssembleJacobian: Cannot prepare indices.\n"); return false;}

	// set elem type in elem disc
	if(!elemDisc.set_geometric_object_type(refID, IDPPN_JACOBIAN))
		{UG_LOG("ERROR in AssembleJacobian: Cannot set geometric object type.\n"); return false;}

	// prepare loop
	if(!elemDisc.prepare_element_loop())
		{UG_LOG("ERROR in AssembleJacobian: Cannot prepare element loop.\n"); return false;}

	// loop over all elements
	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		// get Element
		TElem* elem = *iter;

		// get global indices
		u.update_inner_indices(elem, ind);

		// prepare element
		if(!elemDisc.prepare_element(elem, ind))
			{UG_LOG("ERROR in AssembleJacobian: Cannot prepare element.\n"); return false;}

		// Post process J
		if(!elemDisc.post_process_J(J, ind, time))
			{UG_LOG("ERROR in AssembleJacobian: Cannot assemble local Stiffness Matrix.\n"); return false;}
	}

	// finish element loop
	if(!elemDisc.finish_element_loop())
		{UG_LOG("ERROR in AssembleJacobian: Cannot finish element loop.\n"); return false;}

	return true;
}

template <	typename TElem,
			typename TDiscreteFunction,
			typename TAlgebra>
bool
BNDPostProcessDefect(	IDirichletPostProcess<TAlgebra>& elemDisc,
						typename geometry_traits<TElem>::iterator iterBegin,
						typename geometry_traits<TElem>::iterator iterEnd,
						int si,
						typename TAlgebra::vector_type& d,
						const TDiscreteFunction& u,
						const FunctionGroup& fcts,
						number time, number s_m, number s_a)
{
	typedef typename reference_element_traits<TElem>::reference_element_type reference_element_type;
	const ReferenceObjectID refID = reference_element_type::REFERENCE_OBJECT_ID;

	// check if at least on element exist, else return
	if(iterBegin == iterEnd) return true;

	LocalIndices ind;

	// set functions
	ind.set_function_group(fcts);

	// prepare local indices for elem type
	if(!u.prepare_inner_indices(refID, si, ind))
		{UG_LOG("ERROR in AssembleJacobian: Cannot prepare indices.\n"); return false;}

	// set elem type in elem disc
	if(!elemDisc.set_geometric_object_type(refID, IDPPN_JACOBIAN))
		{UG_LOG("ERROR in AssembleJacobian: Cannot set geometric object type.\n"); return false;}

	// prepare loop
	if(!elemDisc.prepare_element_loop())
		{UG_LOG("ERROR in AssembleJacobian: Cannot prepare element loop.\n"); return false;}

	// loop over all elements
	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		// get Element
		TElem* elem = *iter;

		// get global indices
		u.update_inner_indices(elem, ind);

		// prepare element
		if(!elemDisc.prepare_element(elem, ind))
			{UG_LOG("ERROR in AssembleJacobian: Cannot prepare element.\n"); return false;}

		// Post process d
		if(!elemDisc.post_process_d(d, ind, time))
			{UG_LOG("ERROR in AssembleJacobian: Cannot assemble local Stiffness Matrix.\n"); return false;}
	}

	// finish element loop
	if(!elemDisc.finish_element_loop())
		{UG_LOG("ERROR in AssembleJacobian: Cannot finish element loop.\n"); return false;}

	return true;
}

template <	typename TElem,
			typename TDiscreteFunction,
			typename TAlgebra>
bool
BNDPostProcessLinear(	IDirichletPostProcess<TAlgebra>& elemDisc,
						typename geometry_traits<TElem>::iterator iterBegin,
						typename geometry_traits<TElem>::iterator iterEnd,
						int si,
						typename TAlgebra::matrix_type& Mat,
						typename TAlgebra::vector_type& rhs,
						const TDiscreteFunction& u,
						const FunctionGroup& fcts)
{
	typedef typename reference_element_traits<TElem>::reference_element_type reference_element_type;
	const ReferenceObjectID refID = reference_element_type::REFERENCE_OBJECT_ID;

	// check if at least on element exist, else return
	if(iterBegin == iterEnd) return true;

	LocalIndices ind;

	// set functions
	ind.set_function_group(fcts);

	// prepare local indices for elem type
	if(!u.prepare_inner_indices(refID, si, ind))
		{UG_LOG("ERROR in AssembleJacobian: Cannot prepare indices.\n"); return false;}

	// set elem type in elem disc
	if(!elemDisc.set_geometric_object_type(refID, IDPPN_LINEAR))
		{UG_LOG("ERROR in AssembleJacobian: Cannot set geometric object type.\n"); return false;}

	// prepare loop
	if(!elemDisc.prepare_element_loop())
		{UG_LOG("ERROR in AssembleJacobian: Cannot prepare element loop.\n"); return false;}

	// loop over all elements
	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		// get Element
		TElem* elem = *iter;

		// get global indices
		u.update_inner_indices(elem, ind);

		// prepare element
		if(!elemDisc.prepare_element(elem))
			{UG_LOG("ERROR in AssembleJacobian: Cannot prepare element.\n"); return false;}

		// Post process d
		if(!elemDisc.post_process_J(Mat, ind))
			{UG_LOG("ERROR in AssembleJacobian: Cannot assemble local Stiffness Matrix.\n"); return false;}

		// Post process f
		if(!elemDisc.post_process_f(rhs, ind))
			{UG_LOG("ERROR in AssembleJacobian: Cannot assemble local Stiffness Matrix.\n"); return false;}
	}

	// finish element loop
	if(!elemDisc.finish_element_loop())
		{UG_LOG("ERROR in AssembleJacobian: Cannot finish element loop.\n"); return false;}

	return true;
}

template <	typename TElem,
			typename TDiscreteFunction,
			typename TAlgebra>
bool
BNDSetSolution(			IDirichletPostProcess<TAlgebra>& elemDisc,
						typename geometry_traits<TElem>::iterator iterBegin,
						typename geometry_traits<TElem>::iterator iterEnd,
						int si,
						typename TAlgebra::vector_type& x,
						const TDiscreteFunction& u,
						const FunctionGroup& fcts)
{
	typedef typename reference_element_traits<TElem>::reference_element_type reference_element_type;
	const ReferenceObjectID refID = reference_element_type::REFERENCE_OBJECT_ID;

	// check if at least on element exist, else return
	if(iterBegin == iterEnd) return true;

	LocalIndices ind;

	// set functions
	ind.set_function_group(fcts);

	// prepare local indices for elem type
	if(!u.prepare_inner_indices(refID, si, ind))
		{UG_LOG("ERROR in AssembleJacobian: Cannot prepare indices.\n"); return false;}

	// set elem type in elem disc
	if(!elemDisc.set_geometric_object_type(refID, IDPPN_LINEAR))
		{UG_LOG("ERROR in AssembleJacobian: Cannot set geometric object type.\n"); return false;}

	// prepare loop
	if(!elemDisc.prepare_element_loop())
		{UG_LOG("ERROR in AssembleJacobian: Cannot prepare element loop.\n"); return false;}

	// loop over all elements
	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		// get Element
		TElem* elem = *iter;

		// get global indices
		u.update_inner_indices(elem, ind);

		// prepare element
		if(!elemDisc.prepare_element(elem))
			{UG_LOG("ERROR in AssembleJacobian: Cannot prepare element.\n"); return false;}

		// Post process d
		if(!elemDisc.set_solution(x, ind))
			{UG_LOG("ERROR in AssembleJacobian: Cannot assemble local Stiffness Matrix.\n"); return false;}
	}

	// finish element loop
	if(!elemDisc.finish_element_loop())
		{UG_LOG("ERROR in AssembleJacobian: Cannot finish element loop.\n"); return false;}

	return true;
}

} // end namespace ug
#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__POST_PROCESS__DIRICHLET_BOUNDARY__DIRICHLET_POST_PROCESS_UTIL__ */
