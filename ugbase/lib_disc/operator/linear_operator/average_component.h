/*
 * transfer_post_process.h
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
#include "lib_disc/function_spaces/grid_function_util.h"

#ifdef UG_PARALLEL
#include "lib_disc/parallelization/parallelization_util.h"
#endif

namespace ug{

template <typename TDomain, typename TAlgebra>
class AverageComponent :
	virtual public ITransferPostProcess<TDomain, TAlgebra>
{
	public:
	///	GridFunction type
		typedef GridFunction<TDomain, TAlgebra> GF;

	public:
	///	Constructor setting approximation space
		AverageComponent(const std::string& fcts){m_vCmp = TokenizeTrimString(fcts);};

	///	Constructor setting approximation space
		AverageComponent(const std::vector<std::string>& vCmp){m_vCmp = vCmp;};

	public:
	/// apply Operator, interpolate function
		virtual void post_process(SmartPtr<GF> spGF)
		{
			AdjustMeanValue(spGF, m_vCmp);
		}

	protected:
	///	symbolic function names
		std::vector<std::string> m_vCmp;
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__TRANSFER_POST_PROCESS__ */
