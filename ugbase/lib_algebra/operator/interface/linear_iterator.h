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

#ifndef __H__LIB_ALGEBRA__OPERATOR__INTERFACE__OPERATOR_ITERATOR__
#define __H__LIB_ALGEBRA__OPERATOR__INTERFACE__OPERATOR_ITERATOR__

#include "lib_algebra/operator/damping.h"
#include "common/util/smart_pointer.h"

namespace ug{

///////////////////////////////////////////////////////////////////////////////
// Iterator Operator
///////////////////////////////////////////////////////////////////////////////

/// describes a linear iterator
/**
 * This class is the base class for all linear iterators. Iterators (also called
 * preconditioners) are used in iterative schemes when solving a linear system.
 * Usually, a linear problem like L*u = f is intended to be solved. This is
 * done in an iterative way by performing an iteration of
 *
 * 	start: compute d := f - L*u
 * 	iterate: 	- 	c := B*d 		(compute correction)
 * 				-	u := u + c		(update solution)
 * 				- 	d := d - L*c	(update defect)
 *
 * This iterator class describes the application of B in the scheme above.
 * The application has been split up into two parts.
 *
 * 1. init(L, u) or init(L):
 * 		These methods initialize the iterator and one of these methods has to
 * 		be called before any of the apply methods can be used. Passing the
 * 		linear operator indicates that this operator is used as underlying
 * 		for the iterator.
 *
 * 2. apply or apply_return_defect:
 * 		These methods are used to compute the correction (and to update the
 * 		defect at the same time). Note, that these methods can only be called
 * 		when the iterator has been initialized.
 *
 * This splitting has been made, since initialization may be computationally
 * expensive. Thus, the user of this class has the choice when to call this
 * initialization. E.g. when the operator is applied several times the init of
 * the iterator is only needed once.
 *
 * \tparam	X 	Domain space function
 * \tparam	Y	Range space function
 */
template <typename X, typename Y = X>
class ILinearIterator
{
	public:
	///	Domain space
		using domain_function_type = X;

	///	Range space
		using codomain_function_type = Y;

	public:
	///	returns the name of iterator
	/**
	 * This method returns the name of the iterator operator. This function is
	 * typically needed, when the iterator operator is used inside of another
	 * operator and some debug output should be printed
	 *
	 * \returns 	const char* 	name of inverse operator
	 */
		virtual const char* name() const = 0;

	///	returns if parallel solving is supported
		virtual bool supports_parallel() const = 0;

	///	initialize for operator J(u) and linearization point u
	/**
	 * This method passes the linear operator J(u) that should be used as
	 * underlying by this iterator. As second argument the linearization point
	 * is passed. This is needed e.g. for the geometric multigrid method.
	 *
	 * \param[in]	J		linearized operator to use as underlying
	 * \param[in]	u		linearization point
	 * \returns		bool	success flag
	 */
		virtual bool init(SmartPtr<ILinearOperator<Y,X> > J, const Y& u) = 0;

	///	initialize for linear operator L
	/**
	 * This method passes the operator L that used as underlying by this
	 * operator. In addition some preparation step can be made.
	 *
	 * \param[in]	L		linear operator to use as underlying
	 * \returns		bool	success flag
	 */
		virtual bool init(SmartPtr<ILinearOperator<Y,X> > L) = 0;

	///	compute new correction c = B*d
	/**
	 * This method applies the iterator operator, i.e. c = B*d. The
	 * domain function d remains unchanged.
	 * Note, that this method can always be implemented by creating a copy of
	 * d and calling apply_update_defect with this copy.
	 *
	 * \param[in]	d		defect
	 * \param[out]	c		correction
	 * \returns		bool	success flag
	 */
		virtual bool apply(Y& c, const X& d) = 0;

	///	compute new correction c = B*d and update defect d := d - A*c
	/**
	 * This method applies the inverse operator, i.e. c = B*d. The
	 * domain function d is changed in the way, that the defect d := d - A*c
	 * is returned in the function. This is always useful, when the iterating
	 * algorithm can (or must) update the defect during computation (this is
	 * e.g. the case for the geometric multigrid method).
	 * Note, that this method can always be implemented by calling apply and
	 * then computing d := d - A*c.
	 *
	 * \param[in,out]	d		defect
	 * \param[out]		u		correction
	 * \returns			bool	success flag
	 */
		virtual bool apply_update_defect(Y& c, X& d) = 0;

	///	sets a scaling for the correction
	/**
	 * Sets a scaling for the correction, i.e., once the correction has been
	 * computed, c = B*d, the correction is scaled by a factor, c := s*c, where
	 * s is provided by the passed scaling class. Note, that the scaling factor
	 * may depend on the defect and correction. The internal update of the defect
	 * as done in apply_update_defect must be performed with respect to the
	 * scaled correction.
	 */
		void set_damp(SmartPtr<IDamping<X,Y> > spScaling) {
			m_spDamping = spScaling;
		}

	///	sets the damping to a constant factor
		void set_damp(number factor) {
			m_spDamping = make_sp(new ConstantDamping<X,Y>(factor));
		}

	///	returns the scaling
		SmartPtr<IDamping<X,Y> > damping() {return m_spDamping;}

	///	clone
		virtual SmartPtr<ILinearIterator<X,Y> > clone() = 0;

	/// virtual destructor
		virtual ~ILinearIterator() {};

	///	constructor
		ILinearIterator() {set_damp(1.0);};

	///	copy constructor
		ILinearIterator(const ILinearIterator<X, Y> &parent)
		{
			set_damp(parent.m_spDamping);
		};

		virtual std::string config_string() const
		{
			std::stringstream ss; ss << name() << "( damping = " << m_spDamping->config_string() << ")"; return ss.str();
		}

	protected:
	///	the scaling
		SmartPtr<IDamping<X,Y> > m_spDamping;
};

} // end namespace ug

#endif