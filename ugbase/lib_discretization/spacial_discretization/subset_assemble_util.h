/*
 * subset_assemble_util.h
 *
 *  Created on: 08.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__SUBSET_ASSEMBLE_UTIL__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__SUBSET_ASSEMBLE_UTIL__

// extern includes
#include <iostream>
#include <vector>

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"
#include "lib_algebra/lib_algebra.h"

namespace ug {

//////////////////////////////////
// Jacobian subset assembling
//////////////////////////////////

template <	typename TElemDisc,
			typename TDiscreteFunction,
			typename TAlgebra>
bool
AssembleJacobian(	TElemDisc& elemDisc,
					typename TAlgebra::matrix_type& J,
					const TDiscreteFunction& u,
					const std::vector<size_t>& u_comp, int si,
					number time, number s_m, number s_a)
{
	if(!AssembleJacobian<Triangle>(elemDisc, u.template begin<Triangle>(si), u.template end<Triangle>(si), J, u, u_comp, time, s_m, s_a))
		{UG_LOG("Error in 'AssembleJacobian' while calling 'assemble_jacobian<Triangle>', aborting.\n");return IAssemble_ERROR;}
	if(!AssembleJacobian<Quadrilateral>(elemDisc, u.template begin<Quadrilateral>(si), u.template end<Quadrilateral>(si), J, u, u_comp, time, s_m, s_a))
		{UG_LOG("Error in 'AssembleJacobian' while calling 'assemble_jacobian<Quadrilateral>', aborting.\n");return IAssemble_ERROR;}

	return true;
};

template <	typename TElemDisc,
			typename TDiscreteFunction,
			typename TAlgebra>
bool
AssembleJacobian(	TElemDisc& elemDisc,
					typename TAlgebra::matrix_type& J,
					const TDiscreteFunction& u,
					const std::vector<size_t>& u_comp, int si)
{
	// TODO: This is a costly quick hack, compute matrices directly (without time assembling) !
	return AssembleJacobian<TElemDisc, TDiscreteFunction, TAlgebra>(elemDisc, J, u, u_comp, si, 0.0, 0.0, 1.0);
}

//////////////////////////////////
// Defect subset assembling
//////////////////////////////////

template <	typename TElemDisc,
			typename TDiscreteFunction,
			typename TAlgebra>
bool
AssembleDefect(	TElemDisc& elemDisc,
				typename TAlgebra::vector_type& d,
				const TDiscreteFunction& u,
				const std::vector<size_t>& u_comp, int si,
				number time, number s_m, number s_a)
{
	if(!AssembleDefect<Triangle>(elemDisc, u.template begin<Triangle>(si), u.template end<Triangle>(si), d, u, u_comp, time, s_m, s_a))
		{UG_LOG("Error in AssembleDefect, aborting.\n");return false;}
	if(!AssembleDefect<Quadrilateral>(elemDisc, u.template begin<Quadrilateral>(si), u.template end<Quadrilateral>(si), d, u, u_comp, time, s_m, s_a))
		{UG_LOG("Error in AssembleDefect, aborting.\n");return false;}

	return true;
};

template <	typename TElemDisc,
			typename TDiscreteFunction,
			typename TAlgebra>
bool
AssembleDefect(	TElemDisc& elemDisc,
				typename TAlgebra::vector_type& d,
				const TDiscreteFunction& u,
				const std::vector<size_t>& u_comp, int si)
{
	// TODO: This is a costly quick hack, compute matrices directly (without time assembling) !
	return AssembleDefect<TElemDisc, TDiscreteFunction, TAlgebra>(elemDisc, d, u, u_comp, si, 0.0, 0.0, 1.0);
}


//////////////////////////////////
// Linear subset assembling
//////////////////////////////////

// TODO: Implement time dependent linear case
template <	typename TElemDisc,
			typename TDiscreteFunction,
			typename TAlgebra>
bool
AssembleLinear(	TElemDisc& elemDisc,
				typename TAlgebra::matrix_type& mat,
				typename TAlgebra::vector_type& rhs,
				const TDiscreteFunction& u,
				const std::vector<size_t>& u_comp, int si,
				number time, number s_m, number s_a)
{

/*	if(!AssembleLinear(elemDisc, u.template begin<Triangle>(si), u.template end<Triangle>(si), mat, rhs, u, u_comp, time, s_m, s_a))
		{UG_LOG("Error in AssembleLinear, aborting.\n"); return false;}
	if(!AssembleLinear(elemDisc, u.template begin<Quadrilateral>(si), u.template end<Quadrilateral>(si), mat, rhs, u, u_comp, time, s_m, s_a))
		{UG_LOG("Error in AssembleLinear, aborting.\n"); return false;}
*/
	return false;
}


template <	typename TElemDisc,
			typename TDiscreteFunction,
			typename TAlgebra>
bool
AssembleLinear(	TElemDisc& elemDisc,
				typename TAlgebra::matrix_type& mat,
				typename TAlgebra::vector_type& rhs,
				const TDiscreteFunction& u,
				const std::vector<size_t>& u_comp, int si)
{
	if(!AssembleLinear<Triangle>(elemDisc, u.template begin<Triangle>(si), u.template end<Triangle>(si), mat, rhs, u, u_comp))
		{UG_LOG("Error in AssembleLinear, aborting.\n"); return false;}
	if(!AssembleLinear<Quadrilateral>(elemDisc, u.template begin<Quadrilateral>(si), u.template end<Quadrilateral>(si), mat, rhs, u, u_comp))
		{UG_LOG("Error in AssembleLinear, aborting.\n"); return false;}

	return true;
};

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__SUBSET_ASSEMBLE_UTIL__ */
