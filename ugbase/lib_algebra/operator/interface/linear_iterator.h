/*
 * operator_iterator.h
 *
 *  Created on: 22.02.2010
 *      Author: andreasvogel
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
		typedef X domain_function_type;

	///	Range space
		typedef Y codomain_function_type;

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

#endif /* __H__LIB_ALGEBRA__OPERATOR__INTERFACE__OPERATOR_ITERATOR__ */

