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
#include "interface/linear_operator.h"

namespace ug{

/**
 * Base class for damping of correction in iterative schemes. An iteration for the
 * solution of a matrix problem \f$ A*x = b \f$ is given, for example, by
 * 
 * \f[ x = x + c, \f]
 * 
 * where \f$ c = B*d \f$ is some proposed correction.
 *
 * The damping class now computes a damping factor \f$ \kappa \f$, that is used to
 * damp the correction, i.e.,
 * 
 * \f[ x = x + \kappa c. \f]
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
		virtual number damping(const Y& c, const X& d, ConstSmartPtr<ILinearOperator<Y,X> > spLinOp) const = 0;

	///	returns if the damping is constant
	/**
	 * returns if the damping is constant, i.e. does not depend on correction,
	 * defect and linear operator.
	 *
	 * @return true if constant damping
	 */
		virtual bool constant_damping() const = 0;

	///	returns the constant damping, throws exception if non-constant damping
		virtual number damping() const = 0;

	///	virtual destructor
		virtual ~IDamping() {}

	///	returns information about configuration parameters
		/**
		 * this should return necessary information about parameters and possibly
		 * calling config_string of subcomponents.
		 *
		 * \returns std::string	necessary information about configuration parameters
		 */
		virtual std::string config_string() const = 0;
};

/// constant damping factor
template <typename X, typename Y = X>
class ConstantDamping : public IDamping<X,Y>
{
	public:
		ConstantDamping(number factor) : m_factor(factor) {}

	///	returns the constant damping factor
		virtual number damping(const Y& c, const X& d, ConstSmartPtr<ILinearOperator<Y,X> > spLinOp) const
		{
			return m_factor;
		}

	///	returns the constant damping factor
		virtual number damping() const
		{
			return m_factor;
		}

	///	returns if damping is constant
		virtual bool constant_damping() const {return true;};

		virtual std::string config_string() const
		{
			std::stringstream ss; ss << "ConstantDamping(" << m_factor << ")"; return ss.str();
		}

	protected:
		number m_factor; ///< constant damping factor
};

/// damping computed based on the minimal residuum
template <typename X, typename Y = X>
class MinimalResiduumDamping : public IDamping<X,Y>
{
	public:
	///	returns the damping factor
		virtual number damping(const Y& c, const X& d, ConstSmartPtr<ILinearOperator<Y,X> > spLinOp) const
		{
			SmartPtr<X> spAc = d.clone_without_values();
			X& Ac = *spAc;

			try{
				spLinOp.cast_const()->apply(Ac, c);
			}UG_CATCH_THROW("MinimalResiduumDamping: Computing Ac failed.")

		//	Compute scaling
			try{
				const number kappa = VecProd(d, Ac) / VecProd(Ac, Ac);

				if (kappa<0.3) return 0.3;

			//	return result
				return kappa;
			}UG_CATCH_THROW("MinimalResiduumDamping: Computing (d,Ac)/(Ac,Ac) failed.")
		}

	///	returns if damping is constant
		virtual bool constant_damping() const {return false;};

	///	returns the constant damping factor
		virtual number damping() const
		{
			UG_THROW("MinimalResiduumDamping: non-constant damping.");
		}

		virtual std::string config_string() const
		{
			return "MinimalResiduumDamping";
		}
};

/// damping computed based on the minimal energy
template <typename X, typename Y = X>
class MinimalEnergyDamping : public IDamping<X,Y>
{
	public:
	///	returns the damping factor
		virtual number damping(const Y& c, const X& d, ConstSmartPtr<ILinearOperator<Y,X> > spLinOp) const
		{
			SmartPtr<X> spAc = d.clone_without_values();
			X& Ac = *spAc;
			spLinOp.cast_const()->apply(Ac, c);

		//	Compute scaling
			const number kappa = VecProd(d,c) / VecProd(Ac, c);
			
			if (kappa<0.3) return 0.3;

		//	return result
			return kappa;
		}

	///	returns if damping is constant
		virtual bool constant_damping() const {return false;};

	///	returns the constant damping factor
		virtual number damping() const
		{
			UG_THROW("MinimalEnergyDamping: non-constant damping.");
		}

		virtual std::string config_string() const
		{
			return "MinimalEnergyDamping";
		}
};

}  // end namespace ug

#endif /* __H__LIB_ALGEBRA__OPERATOR__DAMPING__ */
