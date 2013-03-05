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
#include "lib_disc/dof_manager/dof_distribution.h"
#include "lib_disc/function_spaces/approximation_space.h"

namespace ug{

/// Types of constraint
/**
 * This types control the order in with the constraints are performed.
 * constraints with a lower number will be performed first. The order of
 * constraints with equal number is undefined.
 */
/*enum ConstraintType
{
	CT_NONE = 0,
	CT_CONSTRAINTS = 1 << 0,
	CT_DIRICHLET = 1 << 1,
	CT_ALL = CT_NONE | CT_CONSTRAINTS | CT_DIRICHLET
};*/

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
	/// \{
		virtual void adjust_jacobian(matrix_type& J, const vector_type& u,
		                             GridLevel gl, number time = 0.0) = 0;
		virtual void adjust_jacobian(matrix_type& J, const vector_type& u,
		                             ConstSmartPtr<DoFDistribution> dd, number time = 0.0) = 0;
	/// \}

	///	adapts defect to enforce constraints
	/// \{
		virtual void adjust_defect(vector_type& d, const vector_type& u,
		                           GridLevel gl, number time = 0.0) = 0;
		virtual void adjust_defect(vector_type& d, const vector_type& u,
		                           ConstSmartPtr<DoFDistribution> dd, number time = 0.0) = 0;
	/// \}

	///	adapts matrix and rhs (linear case) to enforce constraints
	/// \{
		virtual void adjust_linear(matrix_type& mat, vector_type& rhs,
		                           GridLevel gl, number time = 0.0) = 0;
		virtual void adjust_linear(matrix_type& mat, vector_type& rhs,
		                           ConstSmartPtr<DoFDistribution> dd, number time = 0.0)  = 0;
	/// \}

	///	adapts a rhs to enforce constraints
	/// \{
		virtual void adjust_rhs(vector_type& rhs, const vector_type& u,
		                        GridLevel gl, number time = 0.0) = 0;
		virtual void adjust_rhs(vector_type& rhs, const vector_type& u,
		                        ConstSmartPtr<DoFDistribution> dd, number time = 0.0)  = 0;
	/// \}

	///	sets the constraints in a solution vector
	/// \{
		virtual void adjust_solution(vector_type& u, GridLevel gl,
		                             number time = 0.0) = 0;
		virtual void adjust_solution(vector_type& u, ConstSmartPtr<DoFDistribution> dd,
		                             number time = 0.0)  = 0;
	/// \}


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

template <typename TDomain, typename TAlgebra>
class IDomainConstraint : public IConstraint<TAlgebra>
{
	public:
	///	Domain Type
		typedef TDomain domain_type;

	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	///	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	public:
		using IConstraint<TAlgebra>::adjust_jacobian;
		using IConstraint<TAlgebra>::adjust_defect;
		using IConstraint<TAlgebra>::adjust_linear;
		using IConstraint<TAlgebra>::adjust_rhs;
		using IConstraint<TAlgebra>::adjust_solution;

	public:
	///	sets the approximation space
		virtual void set_approximation_space(SmartPtr<ApproximationSpace<TDomain> > approxSpace)
		{
			m_spApproxSpace = approxSpace;
		}

	///	returns approximation space
		SmartPtr<ApproximationSpace<TDomain> > approximation_space()
		{
			return m_spApproxSpace;
		}

	///	returns approximation space
		ConstSmartPtr<ApproximationSpace<TDomain> > approximation_space() const
		{
			return m_spApproxSpace;
		}

	///	adapts jacobian to enforce constraints
		void adjust_jacobian(matrix_type& J, const vector_type& u,
									 GridLevel gl, number time = 0.0)
		{
			this->adjust_jacobian(J, u, dd(gl), time);
		}

	///	adapts defect to enforce constraints
		void adjust_defect(vector_type& d, const vector_type& u,
								   GridLevel gl, number time = 0.0)
		{
			this->adjust_defect(d, u, dd(gl), time);
		}

	///	adapts matrix and rhs (linear case) to enforce constraints
		void adjust_linear(matrix_type& mat, vector_type& rhs,
								   GridLevel gl, number time = 0.0)
		{
			this->adjust_linear(mat, rhs, dd(gl), time);
		}

	///	adapts a rhs to enforce constraints
		void adjust_rhs(vector_type& rhs, const vector_type& u,
								GridLevel gl, number time = 0.0)
		{
			this->adjust_rhs(rhs, u, dd(gl), time);
		}

	///	sets the constraints in a solution vector
		void adjust_solution(vector_type& u, GridLevel gl,
									 number time = 0.0)
		{
			this->adjust_solution(u, dd(gl), time);
		}

	///	returns the type of constraints
		virtual int type() const = 0;

	///	sets the assemble index for index-wise assemble routine
		void set_ass_index(){set_ass_index(0, false);}
		void set_ass_index(size_t ind, bool index_set = true)
		{
			m_AssAdapter.m_assIndex.index_set = index_set;
			m_AssAdapter.m_assIndex.index = ind;
		}

	protected:
	///	returns the level dof distribution
		ConstSmartPtr<DoFDistribution> dd(const GridLevel& gl) const
		{
			return m_spApproxSpace->dof_distribution(gl);
		}

	protected:
	///	Approximation Space
		SmartPtr<ApproximationSpace<TDomain> > m_spApproxSpace;

	///	Assemble adapter
	//	note: the member variables of the assemble adapter need
	//	to be set explicitly in this class! The vars are NOT derived from the
	//	assemble adapter used in the assemble routine!
		AssAdapter<TAlgebra> m_AssAdapter;

};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__CONSTRAINTS__INTERFACE__ */
