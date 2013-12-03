/*
 * prolongation_operator.h
 *
 *  Created on: 04.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__TRANSFER_POST_PROCESS__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__TRANSFER_POST_PROCESS__

// extern headers
#include <iostream>

// other ug4 modules
#include "common/common.h"
#include "transfer_interface.h"

#ifdef UG_PARALLEL
#include "lib_disc/parallelization/parallelization_util.h"
#endif

namespace ug{

template <typename TDomain, typename TAlgebra>
class AverageComponent :
	virtual public ITransferPostProcess<TDomain, TAlgebra>
{
	public:
	///	Type of algebra
		typedef TAlgebra algebra_type;

	///	Type of Vector
		typedef typename TAlgebra::vector_type vector_type;

	///	Type of Domain
		typedef TDomain domain_type;

	///	GridFunction type
		typedef GridFunction<TDomain, TAlgebra> GF;

	public:
	///	Constructor setting approximation space
		AverageComponent(const std::string& fcts){m_vCmp = TokenizeTrimString(fcts);};

	///	Constructor setting approximation space
		AverageComponent(const std::vector<std::string>& vCmp){m_vCmp = vCmp;};

	public:
	/// apply Operator, interpolate function
		virtual void post_process(SmartPtr<GF> spGF);

	protected:
		template <typename TBaseElem>
		void subtract_value(SmartPtr<GridFunction<TDomain, TAlgebra> > spGF, size_t fct, number sub);

	protected:
	///	symbolic function names
		std::vector<std::string> m_vCmp;
};

} // end namespace ug

#include "transfer_post_process_impl.h"

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__TRANSFER_POST_PROCESS__ */
