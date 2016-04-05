/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
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
		                             ConstSmartPtr<VectorTimeSeries<vector_type> > vSol = SPNULL,
									 const number s_a0 = 1.0) = 0;

	///	adapts defect to enforce constraints
		virtual void adjust_defect(vector_type& d, const vector_type& u,
								   ConstSmartPtr<DoFDistribution> dd, number time = 0.0,
								   ConstSmartPtr<VectorTimeSeries<vector_type> > vSol = SPNULL,
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

	///	adapts correction to enforce constraints
	/// for additive constraints (e.g. linear constraints, but NOT Dirichlet other than Dirichlet-0!),
	/// this is the same as adjust_solution, therefore default implementation
		virtual void adjust_correction(vector_type& c, ConstSmartPtr<DoFDistribution> dd,
									   number time = 0.0) {adjust_solution(c, dd, time);}

	///	adjusts the assembled error estimator values in the attachments according to the constraint
		virtual void adjust_error(const vector_type& u, ConstSmartPtr<DoFDistribution> dd,
								  number time = 0.0) {};

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
	/// constructor
		IDomainConstraint() : IConstraint<TAlgebra>(),
	  		m_bDoErrEst(false), m_spErrEstData(SPNULL) {};

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
		void set_ass_tuner(ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssemblingTuner = NULL)
		{
			m_spAssTuner = spAssemblingTuner;
		}

	// //////////////////////////
	// Error estimator
	// //////////////////////////
	public:
	///	sets the pointer to an error estimator data object (or NULL)
	/**
	 * This function sets the pointer to an error estimator data object
	 * that should be used for this discretization. Note that the ElemDisc
	 * object must use RTTI to try to convert this pointer to the type
	 * of the objects accepted by it for this purpose. If the conversion
	 * fails than an exception must be thrown since this situation is not
	 * allowed.
	 */
		void set_error_estimator(SmartPtr<IErrEstData<TDomain> > ee) {m_spErrEstData = ee; m_bDoErrEst = true;}

	/// find out whether or not a posteriori error estimation is to be performed for this disc
		bool err_est_enabled() const {return m_bDoErrEst;}

	///	returns the pointer to the error estimator data object (or NULL)
		virtual SmartPtr<IErrEstData<TDomain> > err_est_data() {return m_spErrEstData;}

	private:
	/// flag indicating whether or not a posteriori error estimation is to be performed for this disc
		bool m_bDoErrEst;

	protected:
	/// error estimation object associated to the element discretization
		SmartPtr<IErrEstData<TDomain> > m_spErrEstData;

	// ///////////////////////////


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
