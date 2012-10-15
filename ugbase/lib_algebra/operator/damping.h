/*
 * damping.h
 *
 *  Created on: 06.08.2012
 *      Author: Andreas Vogel, Christian Wehner
 */

#ifndef __H__LIB_ALGEBRA__OPERATOR__DAMPING__
#define __H__LIB_ALGEBRA__OPERATOR__DAMPING__

#include "common/common.h"
#include "common/util/smart_pointer.h"

namespace ug{

/**
 * Base class for damping of correction in iterative schemes. An iteration for the
 * solution of a matrix problem A*x = b is given, for example, by
 * 			x = x + c,
 * where c = B*d is some proposed correction.
 *
 * The damping class now computes a damping factor \kappa, that is used to
 * damp the correction, i.e.,
 * 			x = x + \kappa c.
 *
 * In general, the damping may depend on the correction, the (old) defect and
 * the operator A itself.
 */
template <typename X, typename Y = X>
class IDamping
{
	public:
	///	returns the damping
	/**
	 * For a given correction, defect and Operator the damping is returned.
	 *
	 * @param c				the correction
	 * @param d				the defect
	 * @param spLinOp		the operator
	 * @return				the damping
	 */
		virtual number damping(const Y& c, const X& d, SmartPtr<ILinearOperator<Y,X> > spLinOp) = 0;

	///	returns if the damping is constant
	/**
	 * returns if the damping is constant, i.e. does not depend on correction,
	 * defect and linear operator.
	 *
	 * @return true if constant damping
	 */
		virtual bool constant_damping() = 0;

	///	returns the constant damping, throws exception if non-constant damping
		virtual number damping() = 0;

	///	virtual destructor
		virtual ~IDamping() {}
};

/// constant damping factor
template <typename X, typename Y = X>
class ConstantDamping : public IDamping<X,Y>
{
	public:
		ConstantDamping(number factor) : m_factor(factor) {}

	///	returns the constant damping factor
		virtual number damping(const Y& c, const X& d, SmartPtr<ILinearOperator<Y,X> > spLinOp)
		{
			return m_factor;
		}

	///	returns the constant damping factor
		virtual number damping()
		{
			return m_factor;
		}

	///	returns if damping is constant
		virtual bool constant_damping() {return true;};

	protected:
		number m_factor; ///< constant damping factor
};

/// damping computed based on the minimal residuum
template <typename X, typename Y = X>
class MinimalResiduumDamping : public IDamping<X,Y>
{
	public:
	///	returns the damping factor
		virtual number damping(const Y& c, const X& d, SmartPtr<ILinearOperator<Y,X> > spLinOp)
		{
			X Ac; Ac.create(d.size());
			spLinOp->apply(Ac, c);

		//	Compute scaling
			const number kappa = VecProd(d, Ac) / VecProd(Ac, Ac);
			
			if (kappa<0.3) return 0.3;

		//	return result
			return kappa;
		}

	///	returns if damping is constant
		virtual bool constant_damping() {return false;};

	///	returns the constant damping factor
		virtual number damping()
		{
			UG_THROW("MinimalResiduumDamping: non-constant damping.");
		}
};

/// damping computed based on the minimal energy
template <typename X, typename Y = X>
class MinimalEnergyDamping : public IDamping<X,Y>
{
	public:
	///	returns the damping factor
		virtual number damping(const Y& c, const X& d, SmartPtr<ILinearOperator<Y,X> > spLinOp)
		{
			X Ac; Ac.create(d.size());
			spLinOp->apply(Ac, c);

		//	Compute scaling
			const number kappa = VecProd(d,c) / VecProd(Ac, c);
			
			if (kappa<0.3) return 0.3;

		//	return result
			return kappa;
		}

	///	returns if damping is constant
		virtual bool constant_damping() {return false;};

	///	returns the constant damping factor
		virtual number damping()
		{
			UG_THROW("MinimalEnergyDamping: non-constant damping.");
		}
};

}

#endif /* __H__LIB_ALGEBRA__OPERATOR__DAMPING__ */
