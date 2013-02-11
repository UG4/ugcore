/*
 * constraint_interface.h
 *
 *  Created on: 04.08.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__CONSTRAINTS__INTERFACE__
#define __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__CONSTRAINTS__INTERFACE__

// extern headers
#include <vector>

// intern headers
#include "lib_disc/assemble_interface.h"
#include "lib_disc/common/local_algebra.h"

namespace ug{


/// Types of constraint
/**
 * This types control the order in with the constraints are performed.
 * constraints with a lower number will be performed first. The order of
 * constraints with equal number is undefined.
 */
enum ConstraintType
{
	CT_NONE = 0,
	CT_CONSTRAINTS = 1 << 0,
	CT_DIRICHLET = 1 << 1,
	CT_ALL = CT_NONE | CT_CONSTRAINTS | CT_DIRICHLET
};

/// interface for adjustment of constraints
/**
 * This class is the base class for the handling of constraints. Implementations
 * should adjust the jacobian/defect in order to guarantee the constraints. The
 * constraints can also be set in a solution vector by calling the adjust_solution
 * method.
 * Each implementation must specify the type of constraints since in some
 * situations it is important to apply the constraint modification in a certain
 * order.
 *
 * \tparam	TDomain				type of Domain
 * \tparam	TDoFDistribution	type of DoF Distribution
 * \tparam	TAlgebra			type of Algebra
 */
template <typename TAlgebra>
class IConstraint
{
	public:
	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	///	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	public:
	///	adapts jacobian to enforce constraints
		virtual void adjust_jacobian(matrix_type& J, const vector_type& u,
		                             GridLevel gl, number time = 0.0) = 0;

	///	adapts defect to enforce constraints
		virtual void adjust_defect(vector_type& d, const vector_type& u,
		                           GridLevel gl, number time = 0.0) = 0;

	///	adapts matrix and rhs (linear case) to enforce constraints
		virtual void adjust_linear(matrix_type& mat, vector_type& rhs,
		                           GridLevel gl, number time = 0.0) = 0;

	///	adapts a rhs to enforce constraints
		virtual void adjust_rhs(vector_type& rhs, const vector_type& u,
		                        GridLevel gl, number time = 0.0) = 0;

	///	sets the constraints in a solution vector
		virtual void adjust_solution(vector_type& u, GridLevel gl,
		                             number time = 0.0) = 0;

	///	sets the constraints in a solution vector
		virtual void adjust_restriction(vector_type& uCoarse, GridLevel coarseLvl,
										const vector_type& uFine, GridLevel fineLvl) {};

	///	sets the constraints in a solution vector
		virtual void adjust_prolongation(vector_type& uFine, GridLevel fineLvl,
										const vector_type& uCoarse, GridLevel coarseLvl) {};

	///	returns the type of constraints
		virtual int type() const = 0;

	///	virtual destructor
		virtual ~IConstraint() {};
};

// predeclaration
template <typename TDomain> class ApproximationSpace;

template <typename TDomain, typename TAlgebra>
class IDomainConstraint : public IConstraint<TAlgebra>
{
	public:
	///	Domain Type
		typedef TDomain domain_type;

	public:
	///	sets the approximation space
		virtual void set_approximation_space(SmartPtr<ApproximationSpace<TDomain> > approxSpace) = 0;
	///	sets the index for which the assemble operators should be build up
		virtual void set_ass_index() = 0;
		virtual void set_ass_index(size_t ind, bool index_set = true) = 0;

};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__CONSTRAINTS__INTERFACE__ */
