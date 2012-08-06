/*
 * damping.h
 *
 *  Created on: 06.08.2012
 *      Author: Andreas Vogel, Christian Wehner
 */

#ifndef __H__LIB_ALGEBRA__OPERATOR__DAMPING__
#define __H__LIB_ALGEBRA__OPERATOR__DAMPING__

namespace ug{

template <typename X, typename Y = X>
class IDamping
{
	public:
		virtual number damping(const Y& c, const X& d, SmartPtr<ILinearOperator<Y,X> > spLinOp) = 0;
		virtual bool constant_damping() = 0;
		virtual ~IDamping() {}
};

template <typename X, typename Y = X>
class ConstantDamping : public IDamping<X,Y>
{
	public:
		ConstantDamping(number factor) : m_factor(factor) {}

		virtual number damping(const Y& c, const X& d, SmartPtr<ILinearOperator<Y,X> > spLinOp)
		{
			return m_factor;
		}

		virtual bool constant_damping() {return true;};

	protected:
		number m_factor;
};


template <typename X, typename Y = X>
class MinimalResiduumDamping : public IDamping<X,Y>
{
	public:
	MinimalResiduumDamping() {}

		virtual number damping(const Y& c, const X& d, SmartPtr<ILinearOperator<Y,X> > spLinOp)
		{
			X Ac; Ac.create(d.size());
			spLinOp->apply(Ac, c);

		//	Compute scaling
			const number kappa = VecProd(d, Ac) / VecProd(Ac, Ac);

		//	return result
			return kappa;
		}

		virtual bool constant_damping() {return false;};
};

}

#endif /* __H__LIB_ALGEBRA__OPERATOR__DAMPING__ */
