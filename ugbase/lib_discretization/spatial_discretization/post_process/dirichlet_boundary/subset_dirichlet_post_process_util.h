/*
 * subset_dirichlet_post_process_util.h
 *
 *  Created on: 04.08.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__POST_PROCESS__DIRICHLET_BOUNDARY__SUBSET_DIRICHLET_POST_PROCESS_UTIL__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__POST_PROCESS__DIRICHLET_BOUNDARY__SUBSET_DIRICHLET_POST_PROCESS_UTIL__

// extern includes
#include <iostream>
#include <vector>

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"
#include "lib_algebra/lib_algebra.h"

#include "./dirichlet_post_process_util.h"

namespace ug {

//////////////////////////////////
// Jacobian subset post process
//////////////////////////////////

template <	typename TDiscreteFunction,
			typename TAlgebra>
bool
BNDPostProcessJacobian(		IDirichletPostProcess<TAlgebra>& elemDisc,
							typename TAlgebra::matrix_type& J,
							const TDiscreteFunction& u,
							const FunctionGroup& fcts, int si,
							number time, number s_m, number s_a)
{
	if(!BNDPostProcessJacobian<VertexBase>(elemDisc, u.template begin<VertexBase>(si), u.template end<VertexBase>(si), si, J, u, fcts, time, s_m, s_a))
		{UG_LOG("Error in 'BNDPostProcessJacobian' while calling 'AssembleJacobian<VertexBase>', aborting.\n");return IAssemble_ERROR;}

	// TODO: Other elements
	return true;
};

template <	typename TDiscreteFunction,
			typename TAlgebra>
bool
BNDPostProcessJacobian(	IDirichletPostProcess<TAlgebra>& elemDisc,
						typename TAlgebra::matrix_type& J,
						const TDiscreteFunction& u,
						const FunctionGroup& fcts, int si)
{
	// TODO: This is a costly quick hack, compute matrices directly (without time assembling) !
	return BNDPostProcessJacobian<TDiscreteFunction, TAlgebra>(elemDisc, J, u, fcts, si, 0.0, 0.0, 1.0);
}

//////////////////////////////////
// Defect subset post process
//////////////////////////////////

template <	typename TDiscreteFunction,
			typename TAlgebra>
bool
BNDPostProcessDefect(	IDirichletPostProcess<TAlgebra>& elemDisc,
						typename TAlgebra::vector_type& d,
						const TDiscreteFunction& u,
						const FunctionGroup& fcts, int si,
						number time, number s_m, number s_a)
{
	if(!BNDPostProcessDefect<VertexBase>(elemDisc, u.template begin<VertexBase>(si), u.template end<VertexBase>(si), si, d, u, fcts, time, s_m, s_a))
		{UG_LOG("Error in BNDPostProcessDefect, aborting.\n");return false;}

	// TODO: other elements
	return true;
};

template <	typename TDiscreteFunction,
			typename TAlgebra>
bool
BNDPostProcessDefect(	IDirichletPostProcess<TAlgebra>& elemDisc,
						typename TAlgebra::vector_type& d,
						const TDiscreteFunction& u,
						const FunctionGroup& fcts, int si)
{
	// TODO: This is a costly quick hack, compute matrices directly (without time assembling) !
	return BNDPostProcessDefect<TDiscreteFunction, TAlgebra>(elemDisc, d, u, fcts, si, 0.0, 0.0, 1.0);
}


//////////////////////////////////
// Linear subset post process
//////////////////////////////////

// TODO: Implement time dependent linear case
template <	typename TDiscreteFunction,
			typename TAlgebra>
bool
BNDPostProcessLinear(	IDirichletPostProcess<TAlgebra>& elemDisc,
						typename TAlgebra::matrix_type& mat,
						typename TAlgebra::vector_type& rhs,
						const TDiscreteFunction& u,
						const FunctionGroup& fcts, int si,
						number time, number s_m, number s_a)
{
	return false;
}


template <	typename TDiscreteFunction,
			typename TAlgebra>
bool
BNDPostProcessLinear(	IDirichletPostProcess<TAlgebra>& elemDisc,
						typename TAlgebra::matrix_type& mat,
						typename TAlgebra::vector_type& rhs,
						const TDiscreteFunction& u,
						const FunctionGroup& fcts, int si)
{
	UG_DLOG(LIB_DISC_ASSEMBLE, 3, "Assembling " << u.template num<VertexBase>(si) << " Edges on subset " << si << ".\n");
	if(!BNDPostProcessLinear<VertexBase>(elemDisc, u.template begin<VertexBase>(si), u.template end<VertexBase>(si), si, mat, rhs, u, fcts))
		{UG_LOG("Error in BNDPostProcessLinear, aborting.\n"); return false;}

	// TODO: other elements
	return true;
};


//////////////////////////////////
// Solution subset post process
//////////////////////////////////

template <	typename TDiscreteFunction,
			typename TAlgebra>
bool
BNDSetSolution(			IDirichletPostProcess<TAlgebra>& elemDisc,
						typename TAlgebra::vector_type& x,
						const TDiscreteFunction& u,
						const FunctionGroup& fcts, int si, number time)
{
	UG_DLOG(LIB_DISC_ASSEMBLE, 3, "Assembling " << u.template num<VertexBase>(si) << " Edges on subset " << si << ".\n");
	if(!BNDSetSolution<VertexBase>(elemDisc, u.template begin<VertexBase>(si), u.template end<VertexBase>(si), si, x, u, fcts, time))
		{UG_LOG("Error in BNDSetSolution, aborting.\n"); return false;}

	// TODO: other elements
	return true;
};

template <	typename TDiscreteFunction,
			typename TAlgebra>
bool
BNDSetSolution(			IDirichletPostProcess<TAlgebra>& elemDisc,
						typename TAlgebra::vector_type& x,
						const TDiscreteFunction& u,
						const FunctionGroup& fcts, int si)
{
	UG_DLOG(LIB_DISC_ASSEMBLE, 3, "Assembling " << u.template num<VertexBase>(si) << " Edges on subset " << si << ".\n");
	if(!BNDSetSolution<VertexBase>(elemDisc, u.template begin<VertexBase>(si), u.template end<VertexBase>(si), si, x, u, fcts))
		{UG_LOG("Error in BNDSetSolution, aborting.\n"); return false;}

	// TODO: other elements
	return true;
};



} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__POST_PROCESS__DIRICHLET_BOUNDARY__SUBSET_DIRICHLET_POST_PROCESS_UTIL__ */
