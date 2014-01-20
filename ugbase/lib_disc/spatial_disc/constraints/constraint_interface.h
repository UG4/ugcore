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
		                             ConstSmartPtr<DoFDistribution> dd, number time = 0.0,
		                             ConstSmartPtr<VectorTimeSeries<vector_type> > vSol = ConstSmartPtr<VectorTimeSeries<vector_type> >(),
									 const number s_a0 = 1.0) = 0;

	///	adapts defect to enforce constraints
		virtual void adjust_defect(vector_type& d, const vector_type& u,
		                           ConstSmartPtr<DoFDistribution> dd, number time = 0.0,
		                           ConstSmartPtr<VectorTimeSeries<vector_type> > vSol = ConstSmartPtr<VectorTimeSeries<vector_type> >(),
		                           const std::vector<number>* vScaleMass = NULL,
		                           const std::vector<number>* vScaleStiff = NULL) = 0;

	///	adapts matrix and rhs (linear case) to enforce constraints
		virtual void adjust_linear(matrix_type& mat, vector_type& rhs,
		                           ConstSmartPtr<DoFDistribution> dd, number time = 0.0)  = 0;

	///	adapts a rhs to enforce constraints
		virtual void adjust_rhs(vector_type& rhs, const vector_type& u,
		                        ConstSmartPtr<DoFDistribution> dd, number time = 0.0)  = 0;

	///	sets the constraints in a solution vector
		virtual void adjust_solution(vector_type& u, ConstSmartPtr<DoFDistribution> dd,
		                             number time = 0.0)  = 0;

	///	sets constraints in prolongation
		virtual void adjust_prolongation(matrix_type& P,
		                                 ConstSmartPtr<DoFDistribution> ddFine,
		                                 ConstSmartPtr<DoFDistribution> ddCoarse,
		                                 number time = 0.0) {};

	///	sets constraints in restriction
		virtual void adjust_restriction(matrix_type& R,
		                                ConstSmartPtr<DoFDistribution> ddCoarse,
		                                ConstSmartPtr<DoFDistribution> ddFine,
		                                number time = 0.0) {};

	///	sets the constraints in a solution vector
		virtual void adjust_restriction(vector_type& uCoarse, GridLevel coarseLvl,
										const vector_type& uFine, GridLevel fineLvl) {};

	///	sets the constraints in a solution vector
		virtual void adjust_prolongation(vector_type& uFine, GridLevel fineLvl,
										const vector_type& uCoarse, GridLevel coarseLvl) {};

	///	modifies solution vector before calling the assembling routine
		virtual void modify_solution(vector_type& uMod, const vector_type& u,
										ConstSmartPtr<DoFDistribution> dd) {};

	///	modify_solution for instationary case
		virtual void modify_solution(SmartPtr<VectorTimeSeries<vector_type> > vSolMod,
				ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
				ConstSmartPtr<DoFDistribution> dd) {};

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

	///	returns the type of constraints
		virtual int type() const = 0;

	///	sets the assemble adapter for the constraints
		void set_ass_tuner(ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssemblingTuner = ConstSmartPtr<AssemblingTuner<TAlgebra> >())
		{
			m_spAssTuner = spAssemblingTuner;
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
		ConstSmartPtr<AssemblingTuner<TAlgebra> > m_spAssTuner;

};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__CONSTRAINTS__INTERFACE__ */
