/*
 * transfer_interface.h
 *
 *  Created on: 20.12.2011
 *      Author: andreasvogel
 */

#ifndef TRANSFER_INTERFACE_H_
#define TRANSFER_INTERFACE_H_

#include "lib_algebra/operator/operator_interface.h"
#include "lib_disc/dof_manager/grid_level.h"

namespace ug{

//predeclaration
template <typename TAlgebra>
class IConstraint;


///////////////////////////////////////////////////////////
// Prolongation Operator
///////////////////////////////////////////////////////////

template <typename TAlgebra>
class IProlongationOperator :
	public virtual ILinearOperator<	typename TAlgebra::vector_type,
									typename TAlgebra::vector_type>
{
	public:
	//	vector type
		typedef typename TAlgebra::vector_type vector_type;

	// 	Domain space
		typedef vector_type domain_function_type;

	// 	Range space
		typedef vector_type codomain_function_type;

	public:
	// 	Apply Transposed Operator u = L^T*f
		virtual void apply_transposed(vector_type& u, const vector_type& f) = 0;

	// 	Set Levels for Prolongation coarse -> fine
		virtual void set_levels(GridLevel coarseLevel, GridLevel fineLevel) = 0;

	///	clears dirichlet post processes
		virtual void clear_constraints() = 0;

	///	adds a dirichlet post process (not added if already registered)
		virtual void add_constraint(IConstraint<TAlgebra>& pp) = 0;

	///	removes a post process
		virtual void remove_constraint(IConstraint<TAlgebra>& pp) = 0;

	//	Clone
		virtual IProlongationOperator<TAlgebra>* clone() = 0;
};

///////////////////////////////////////////////////////////
// Projection Operator
///////////////////////////////////////////////////////////

template <typename X, typename Y>
class IProjectionOperator :	public virtual ILinearOperator<X,Y>
{
	public:
	// 	Domain space
		typedef X domain_function_type;

	// 	Range space
		typedef Y codomain_function_type;

	public:
	// 	Apply Transposed Operator u = L^T*f
		virtual void apply_transposed(X& u, const Y& f) = 0;

	// 	Set Levels for Prolongation coarse -> fine
		virtual void set_levels(GridLevel coarseLevel, GridLevel fineLevel) = 0;

	//	Clone
		virtual IProjectionOperator<X,Y>* clone() = 0;
};

} // end namespace ug

#endif /* TRANSFER_INTERFACE_H_ */
