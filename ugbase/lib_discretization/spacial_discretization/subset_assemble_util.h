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
// Stiffness Matrix subset assembling
//////////////////////////////////

template <	typename TElemDisc,
			typename TDiscreteFunction,
			typename TAlgebra>
bool
AssembleStiffnessMatrix(	TElemDisc& elemDisc,
							typename TAlgebra::matrix_type& J,
							const TDiscreteFunction& u,
							const FunctionGroup& fcts, int si, int dim)
{

	switch(dim)
	{
		case 1:
			if(!AssembleStiffnessMatrix<Edge>(elemDisc, u.template begin<Edge>(si), u.template end<Edge>(si), si, J, u, fcts))
				{UG_LOG("Error in 'AssembleStiffnessMatrix' while calling 'assemble_jacobian<Edge>', aborting.\n");return IAssemble_ERROR;}
			break;

		case 2:
			if(!AssembleStiffnessMatrix<Triangle>(elemDisc, u.template begin<Triangle>(si), u.template end<Triangle>(si), si, J, u, fcts))
				{UG_LOG("Error in 'AssembleStiffnessMatrix' while calling 'assemble_jacobian<Triangle>', aborting.\n");return IAssemble_ERROR;}
			if(!AssembleStiffnessMatrix<Quadrilateral>(elemDisc, u.template begin<Quadrilateral>(si), u.template end<Quadrilateral>(si), si, J, u, fcts))
				{UG_LOG("Error in 'AssembleStiffnessMatrix' while calling 'assemble_jacobian<Quadrilateral>', aborting.\n");return IAssemble_ERROR;}
			break;

		case 3:
			if(!AssembleStiffnessMatrix<Tetrahedron>(elemDisc, u.template begin<Tetrahedron>(si), u.template end<Tetrahedron>(si), si, J, u, fcts))
				{UG_LOG("Error in 'AssembleStiffnessMatrix' while calling 'assemble_jacobian<Tetrahedron>', aborting.\n");return IAssemble_ERROR;}
			if(!AssembleStiffnessMatrix<Prism>(elemDisc, u.template begin<Prism>(si), u.template end<Prism>(si), si, J, u, fcts))
				{UG_LOG("Error in 'AssembleStiffnessMatrix' while calling 'assemble_jacobian<Prism>', aborting.\n");return IAssemble_ERROR;}
			if(!AssembleStiffnessMatrix<Hexahedron>(elemDisc, u.template begin<Hexahedron>(si), u.template end<Hexahedron>(si), si, J, u, fcts))
				{UG_LOG("Error in 'AssembleStiffnessMatrix' while calling 'assemble_jacobian<Hexahedron>', aborting.\n");return IAssemble_ERROR;}
			break;

		default: UG_LOG("Dimension " << dim << " not supported.\n"); return false;
	}

	return true;
};


//////////////////////////////////
// Mass Matrix subset assembling
//////////////////////////////////

template <	typename TElemDisc,
			typename TDiscreteFunction,
			typename TAlgebra>
bool
AssembleMassMatrix(	TElemDisc& elemDisc,
					typename TAlgebra::matrix_type& J,
					const TDiscreteFunction& u,
					const FunctionGroup& fcts, int si, int dim)
{

	switch(dim)
	{
		case 1:
			if(!AssembleMassMatrix<Edge>(elemDisc, u.template begin<Edge>(si), u.template end<Edge>(si), si, J, u, fcts))
				{UG_LOG("Error in 'AssembleMassMatrix' while calling 'assemble_jacobian<Edge>', aborting.\n");return IAssemble_ERROR;}
			break;

		case 2:
			if(!AssembleMassMatrix<Triangle>(elemDisc, u.template begin<Triangle>(si), u.template end<Triangle>(si), si, J, u, fcts))
				{UG_LOG("Error in 'AssembleMassMatrix' while calling 'assemble_jacobian<Triangle>', aborting.\n");return IAssemble_ERROR;}
			if(!AssembleMassMatrix<Quadrilateral>(elemDisc, u.template begin<Quadrilateral>(si), u.template end<Quadrilateral>(si), si, J, u, fcts))
				{UG_LOG("Error in 'AssembleMassMatrix' while calling 'assemble_jacobian<Quadrilateral>', aborting.\n");return IAssemble_ERROR;}
			break;

		case 3:
			if(!AssembleMassMatrix<Tetrahedron>(elemDisc, u.template begin<Tetrahedron>(si), u.template end<Tetrahedron>(si), si, J, u, fcts))
				{UG_LOG("Error in 'AssembleMassMatrix' while calling 'assemble_jacobian<Tetrahedron>', aborting.\n");return IAssemble_ERROR;}
			if(!AssembleMassMatrix<Prism>(elemDisc, u.template begin<Prism>(si), u.template end<Prism>(si), si, J, u, fcts))
				{UG_LOG("Error in 'AssembleMassMatrix' while calling 'assemble_jacobian<Prism>', aborting.\n");return IAssemble_ERROR;}
			if(!AssembleMassMatrix<Hexahedron>(elemDisc, u.template begin<Hexahedron>(si), u.template end<Hexahedron>(si), si, J, u, fcts))
				{UG_LOG("Error in 'AssembleMassMatrix' while calling 'assemble_jacobian<Prism>', aborting.\n");return IAssemble_ERROR;}
			break;

		default: UG_LOG("Dimension " << dim << " not supported.\n"); return false;
	}

	return true;
};

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
					const FunctionGroup& fcts, int si, int dim,
					number time, number s_m, number s_a)
{

	switch(dim)
	{
		case 1:
			if(!AssembleJacobian<Edge>(elemDisc, u.template begin<Edge>(si), u.template end<Edge>(si), si, J, u, fcts, time, s_m, s_a))
				{UG_LOG("Error in 'AssembleJacobian' while calling 'assemble_jacobian<Edge>', aborting.\n");return IAssemble_ERROR;}
			break;

		case 2:
			if(!AssembleJacobian<Triangle>(elemDisc, u.template begin<Triangle>(si), u.template end<Triangle>(si), si, J, u, fcts, time, s_m, s_a))
				{UG_LOG("Error in 'AssembleJacobian' while calling 'assemble_jacobian<Triangle>', aborting.\n");return IAssemble_ERROR;}
			if(!AssembleJacobian<Quadrilateral>(elemDisc, u.template begin<Quadrilateral>(si), u.template end<Quadrilateral>(si), si, J, u, fcts, time, s_m, s_a))
				{UG_LOG("Error in 'AssembleJacobian' while calling 'assemble_jacobian<Quadrilateral>', aborting.\n");return IAssemble_ERROR;}
			break;

		case 3:
			if(!AssembleJacobian<Tetrahedron>(elemDisc, u.template begin<Tetrahedron>(si), u.template end<Tetrahedron>(si), si, J, u, fcts, time, s_m, s_a))
				{UG_LOG("Error in 'AssembleJacobian' while calling 'assemble_jacobian<Tetrahedron>', aborting.\n");return IAssemble_ERROR;}
			if(!AssembleJacobian<Prism>(elemDisc, u.template begin<Prism>(si), u.template end<Prism>(si), si, J, u, fcts, time, s_m, s_a))
				{UG_LOG("Error in 'AssembleJacobian' while calling 'assemble_jacobian<Prism>', aborting.\n");return IAssemble_ERROR;}
			if(!AssembleJacobian<Hexahedron>(elemDisc, u.template begin<Hexahedron>(si), u.template end<Hexahedron>(si), si, J, u, fcts, time, s_m, s_a))
				{UG_LOG("Error in 'AssembleJacobian' while calling 'assemble_jacobian<Prism>', aborting.\n");return IAssemble_ERROR;}
			break;

		default: UG_LOG("Dimension " << dim << " not supported.\n"); return false;
	}

	return true;
};

template <	typename TElemDisc,
			typename TDiscreteFunction,
			typename TAlgebra>
bool
AssembleJacobian(	TElemDisc& elemDisc,
					typename TAlgebra::matrix_type& J,
					const TDiscreteFunction& u,
					const FunctionGroup& fcts, int si, int dim)
{
	// TODO: This is a costly quick hack, compute matrices directly (without time assembling) !
	return AssembleJacobian<TElemDisc, TDiscreteFunction, TAlgebra>(elemDisc, J, u, fcts, si, dim, 0.0, 0.0, 1.0);
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
				const FunctionGroup& fcts, int si, int dim,
				number time, number s_m, number s_a)
{
	switch(dim)
	{
		case 1:
			if(!AssembleDefect<Edge>(elemDisc, u.template begin<Edge>(si), u.template end<Edge>(si), si, d, u, fcts, time, s_m, s_a))
				{UG_LOG("Error in AssembleDefect, aborting.\n");return false;}
			break;

		case 2:
			if(!AssembleDefect<Triangle>(elemDisc, u.template begin<Triangle>(si), u.template end<Triangle>(si), si, d, u, fcts, time, s_m, s_a))
				{UG_LOG("Error in AssembleDefect, aborting.\n");return false;}
			if(!AssembleDefect<Quadrilateral>(elemDisc, u.template begin<Quadrilateral>(si), u.template end<Quadrilateral>(si), si, d, u, fcts, time, s_m, s_a))
				{UG_LOG("Error in AssembleDefect, aborting.\n");return false;}
			break;

		case 3:
			if(!AssembleDefect<Tetrahedron>(elemDisc, u.template begin<Tetrahedron>(si), u.template end<Tetrahedron>(si), si, d, u, fcts, time, s_m, s_a))
				{UG_LOG("Error in AssembleDefect, aborting.\n");return false;}
			if(!AssembleDefect<Prism>(elemDisc, u.template begin<Prism>(si), u.template end<Prism>(si), si, d, u, fcts, time, s_m, s_a))
				{UG_LOG("Error in AssembleDefect, aborting.\n");return false;}
			if(!AssembleDefect<Hexahedron>(elemDisc, u.template begin<Hexahedron>(si), u.template end<Hexahedron>(si), si, d, u, fcts, time, s_m, s_a))
				{UG_LOG("Error in AssembleDefect, aborting.\n");return false;}
			break;

		default: UG_LOG("Dimension " << dim << " not supported.\n"); return false;
	}

	return true;
};

template <	typename TElemDisc,
			typename TDiscreteFunction,
			typename TAlgebra>
bool
AssembleDefect(	TElemDisc& elemDisc,
				typename TAlgebra::vector_type& d,
				const TDiscreteFunction& u,
				const FunctionGroup& fcts, int si, int dim)
{
	// TODO: This is a costly quick hack, compute matrices directly (without time assembling) !
	return AssembleDefect<TElemDisc, TDiscreteFunction, TAlgebra>(elemDisc, d, u, fcts, si, dim, 0.0, 0.0, 1.0);
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
				const FunctionGroup& fcts, int si, int dim,
				number time, number s_m, number s_a)
{

/*	if(!AssembleLinear(elemDisc, u.template begin<Triangle>(si), u.template end<Triangle>(si), mat, rhs, u, fcts, time, s_m, s_a))
		{UG_LOG("Error in AssembleLinear, aborting.\n"); return false;}
	if(!AssembleLinear(elemDisc, u.template begin<Quadrilateral>(si), u.template end<Quadrilateral>(si), mat, rhs, u, fcts, time, s_m, s_a))
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
				const FunctionGroup& fcts, int si, int dim)
{
	switch(dim)
	{
		case 1:
			UG_DLOG(LIB_DISC_ASSEMBLE, 3, "Assembling " << u.template num<Edge>(si) << " Edges on subset " << si << ".\n");
			if(!AssembleLinear<Edge>(elemDisc, u.template begin<Edge>(si), u.template end<Edge>(si), si, mat, rhs, u, fcts))
				{UG_LOG("Error in AssembleLinear, aborting.\n"); return false;}
			break;

		case 2:
			UG_DLOG(LIB_DISC_ASSEMBLE, 3, "Assembling " << u.template num<Triangle>(si) << " Triangles on subset " << si << ".\n");
			if(!AssembleLinear<Triangle>(elemDisc, u.template begin<Triangle>(si), u.template end<Triangle>(si), si, mat, rhs, u, fcts))
				{UG_LOG("Error in AssembleLinear, aborting.\n"); return false;}

			UG_DLOG(LIB_DISC_ASSEMBLE, 3, "Assembling " << u.template num<Quadrilateral>(si) << " Quadrilaterals on subset " << si << ".\n");
			if(!AssembleLinear<Quadrilateral>(elemDisc, u.template begin<Quadrilateral>(si), u.template end<Quadrilateral>(si), si, mat, rhs, u, fcts))
				{UG_LOG("Error in AssembleLinear, aborting.\n"); return false;}
			break;

		case 3:
			UG_DLOG(LIB_DISC_ASSEMBLE, 3, "Assembling " << u.template num<Tetrahedron>(si) << " Tetrahedrons on subset " << si << ".\n");
			if(!AssembleLinear<Tetrahedron>(elemDisc, u.template begin<Tetrahedron>(si), u.template end<Tetrahedron>(si), si, mat, rhs, u, fcts))
				{UG_LOG("Error in AssembleLinear, aborting.\n"); return false;}
			UG_DLOG(LIB_DISC_ASSEMBLE, 3, "Assembling " << u.template num<Prism>(si) << " Prisms on subset " << si << ".\n");
			if(!AssembleLinear<Prism>(elemDisc, u.template begin<Prism>(si), u.template end<Prism>(si), si, mat, rhs, u, fcts))
				{UG_LOG("Error in AssembleLinear, aborting.\n"); return false;}
			UG_DLOG(LIB_DISC_ASSEMBLE, 3, "Assembling " << u.template num<Prism>(si) << " Prisms on subset " << si << ".\n");
			if(!AssembleLinear<Hexahedron>(elemDisc, u.template begin<Hexahedron>(si), u.template end<Hexahedron>(si), si, mat, rhs, u, fcts))
				{UG_LOG("Error in AssembleLinear, aborting.\n"); return false;}
			break;

		default: UG_LOG("Dimension " << dim << " not supported.\n"); return false;
	}

	return true;
};

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__SUBSET_ASSEMBLE_UTIL__ */
