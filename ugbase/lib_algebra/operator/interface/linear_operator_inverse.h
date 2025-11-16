/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__LIB_ALGEBRA__OPERATOR__INTERFACE__LINEAR_OPERATOR_INVERSE__
#define __H__LIB_ALGEBRA__OPERATOR__INTERFACE__LINEAR_OPERATOR_INVERSE__

#include "linear_operator.h"
#include "linear_iterator.h"
#include "lib_algebra/operator/convergence_check.h"
#include "common/error.h"
#include "common/util/smart_pointer.h"

namespace ug{


///////////////////////////////////////////////////////////////////////////////
// Inverse of a Linear Operator
///////////////////////////////////////////////////////////////////////////////

/// describes an inverse linear mapping X->Y
/**
 * This class is the base class for the inversion of linear operator given in
 * form of the ILinearOperator interface class. Given a operator L, the basic
 * usage of this class is to invert this operator, i.e. to compute the solution
 * u of
 *
 * 		L*u = f     i.e. u := L^{-1} f
 *
 * This application has been split up into three steps:
 *
 * 1. init():  This method initializes the inverse operator. The inverse operator
 * 			   is initialized the way that, its application will be the inverse
 * 			   application of the operator L passed in by this function. The
 * 			   prepare method can only be called, when this method has been
 * 			   called once.
 *
 * 3. apply():	This method performs the inversion. Before this method is called
 * 				the init and prepare methods have to be called.
 *
 * This splitting has been made, since initialization and preparation may be
 * computationally expansive. Thus, the user of this class has the choice
 * when to call this initialization/preparation. E.g. when the operator is
 * applied several times on the same vectors, those have only to be prepared
 * once and the init of the operator is only needed once.
 *
 * \tparam	X		domain space
 * \tparam	Y		range space
 */
template <typename X, typename Y = X>
class ILinearOperatorInverse : public ILinearIterator<X,Y>
{
	public:
	///	Domain space
		using domain_function_type = X;

	///	Range space
		using codomain_function_type = Y;

	public:
	///	constructor setting convergence check to (100, 1e-12, 1e-12, true)
		ILinearOperatorInverse()
			: m_spLinearOperator(nullptr),
			  m_spConvCheck(new StdConvCheck<X>(100, 1e-12, 1e-12, true))
		{}

	///	Default constructor
		ILinearOperatorInverse(SmartPtr<IConvergenceCheck<X> > spConvCheck)
			: m_spLinearOperator(nullptr),
			  m_spConvCheck(spConvCheck)
		{}

	/// virtual destructor
		virtual ~ILinearOperatorInverse() {};

	///	returns the name of the operator inverse
	/**
	 * This method returns the name of the inverse operator. This function is
	 * typically needed, when the inverse operator is used inside of another and
	 * some debug output should be printed
	 *
	 * \returns 	const char* 	name of inverse operator
	 */
		virtual const char* name() const = 0;


	///	returns information about configuration parameters
	/**
	 * this should return necessary information about parameters and possibly
	 * calling config_string of subcomponents.
	 *
	 * \returns std::string	necessary information about configuration parameters
	 */
		virtual std::string config_string() const { return name(); }

	///	returns if parallel solving is supported
		virtual bool supports_parallel() const = 0;

	/// initializes for the inverse for a linear operator
	/**
	 * This method passes the operator L that is inverted by this operator. In
	 * addition some preparation step can be made.
	 *
	 * \param[in]	L		linear operator to invert
	 * \returns		bool	success flag
	 */
		virtual bool init(SmartPtr<ILinearOperator<Y,X> > L)
		{
		//	remember operator
			m_spLinearOperator = L;
			return true;
		}

	/// initializes for the inverse for a linearized operator at linearization point u
	/**
	 * This method passes the linear operator J(u) that should be inverted by
	 * this operator. As second argument the linearization point is passed.
	 * This is needed e.g. for the geometric multigrid method, that inverts
	 * a linearized operator based on coarser grid operators, that have to be
	 * initialized based on the linearization point.
	 *
	 * \param[in]	J		linearized operator to invert
	 * \param[in]	u		linearization point
	 * \returns		bool	success flag
	 */
		virtual bool init(SmartPtr<ILinearOperator<Y,X> > J, const Y& u)
		{
		//	remember operator
			m_spLinearOperator = J;
			return true;
		}

	///	applies inverse operator, i.e. returns u = A^{-1} f
	/**
	 * This method applies the inverse operator, i.e. u = A^{-1} f. The
	 * domain function f remains unchanged.
	 * Note, that this method can always be implemented by creating a copy of
	 * f and calling apply_return_defect with this copy.
	 *
	 * \param[in]	f		right-hand side
	 * \param[out]	u		solution
	 * \returns		bool	success flag
	 */
		virtual bool apply(Y& u, const X& f) = 0;

	///	applies inverse operator, i.e. returns u = A^{-1} f and returns defect d := f - A*u
	/**
	 * This method applies the inverse operator, i.e. u = A^{-1} f. The
	 * domain function f is changed in the way, that the defect d := f - A*u
	 * is returned in the function. This is always useful, when the inverting
	 * algorithm can (or must) update the defect during computation (this is
	 * e.g. the case for the geometric multigrid method).
	 * Note, that this method can always be implemented by calling apply and
	 * then computing d := f - A*u.
	 *
	 * \param[in,out]	f		right-hand side
	 * \param[out]		u		solution
	 * \returns			bool	success flag
	 */
		virtual bool apply_return_defect(Y& u, X& f) = 0;

		virtual bool apply_update_defect(Y& u, X& f)
		{
			return apply_return_defect(u,f);
		}

		virtual SmartPtr<ILinearIterator<X,Y> > clone()
		{
			UG_THROW("No cloning implemented.");
			return nullptr;
		}

	///	returns the convergence check
		ConstSmartPtr<IConvergenceCheck<X> > convergence_check() const {return m_spConvCheck;}

	///	returns the convergence check
		SmartPtr<IConvergenceCheck<X> > convergence_check() {return m_spConvCheck;}

	/// returns the current defect
		number defect() const {return convergence_check()->defect();}

	/// returns the current number of steps
		int step() const {return convergence_check()->step();}

	/// returns the current relative reduction
		number reduction() const {return convergence_check()->reduction();}

	///	returns the standard offset for output
		virtual int standard_offset() const {return 3;}

	///	set the convergence check
		void set_convergence_check(SmartPtr<IConvergenceCheck<X> > spConvCheck)
		{
			m_spConvCheck = spConvCheck;
			m_spConvCheck->set_offset(standard_offset());
		};

	///	returns the current Operator this Inverse Operator is initialized for
		SmartPtr<ILinearOperator<Y,X> > linear_operator()
		{
			if(m_spLinearOperator.invalid())
				UG_THROW(name() << ": Linear Operator that should be "
				               	 "inverted has not been set.");

			return m_spLinearOperator;
		}

	protected:
	/// Operator that is inverted by this Inverse Operator
		SmartPtr<ILinearOperator<Y,X> > m_spLinearOperator;

	///	smart pointer holding the convergence check
		SmartPtr<IConvergenceCheck<X> > m_spConvCheck;
};

}
#endif /* __H__LIB_ALGEBRA__OPERATOR__INTERFACE__LINEAR_OPERATOR_INVERSE__ */
