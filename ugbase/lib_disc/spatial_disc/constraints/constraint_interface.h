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
	CT_NONE = -1,
	CT_CONSTRAINTS = 0,
	CT_DIRICHLET = 1,
	NUM_CONSTRAINT_TYPES
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
template <	typename TDoFDistribution,
			typename TAlgebra>
class IConstraint{
	public:
	///	DoF Distribution Type
		typedef IDoFDistribution<TDoFDistribution> dof_distribution_type;

	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	///	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	public:
	///	adapts jacobian to enforce constraints
		virtual bool adjust_jacobian(matrix_type& J, const vector_type& u,
		                             const dof_distribution_type& dofDistr,
		                             number time = 0.0) = 0;

	///	adapts defect to enforce constraints
		virtual bool adjust_defect(vector_type& d, const vector_type& u,
		                           const dof_distribution_type& dofDistr,
		                           number time = 0.0) = 0;

	///	adapts matrix and rhs (linear case) to enforce constraints
		virtual bool adjust_linear(matrix_type& mat, vector_type& rhs, const vector_type& u,
		                           const dof_distribution_type& dofDistr,
		                           number time = 0.0) = 0;

	///	adapts a rhs to enforce constraints
		virtual bool adjust_rhs(vector_type& rhs, const vector_type& u,
		                        const dof_distribution_type& dofDistr,
		                        number time = 0.0) = 0;

	///	sets the constraints in a solution vector
		virtual bool adjust_solution(vector_type& u,
		                             const dof_distribution_type& dofDistr,
		                             number time = 0.0) = 0;

	///	returns the type of constraints
		virtual int type() = 0;

	///	virtual destructor
		virtual ~IConstraint() {};
};

// predeclaration
template <typename TDomain> class IApproximationSpace;

template <	typename TDomain,
			typename TDoFDistribution,
			typename TAlgebra>
class IDomainConstraint
	: public IConstraint<TDoFDistribution, TAlgebra>
{
	public:
	///	Domain Type
		typedef TDomain domain_type;

	public:
	///	sets the approximation space
		virtual void set_approximation_space(IApproximationSpace<TDomain>& approxSpace) = 0;

};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__CONSTRAINTS__INTERFACE__ */
