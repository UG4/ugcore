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
	virtual public ITransferPostProcess<TAlgebra>
{
	public:
	///	Type of algebra
		typedef TAlgebra algebra_type;

	///	Type of Vector
		typedef typename TAlgebra::vector_type vector_type;

	///	Type of Domain
		typedef TDomain domain_type;

	public:
	///	Constructor setting approximation space
		AverageComponent(
				SmartPtr<ApproximationSpace<TDomain> > approxSpace,
				const char* fcts) :
			m_symbFct(fcts),
			m_spApproxSpace(approxSpace), m_bInit(false)
		{
			init();
		};

		virtual ~AverageComponent(){};

	public:
	///	Set levels
		virtual void set_levels(GridLevel level);

	///	initialize the operator
		virtual void init();

	/// apply Operator, interpolate function
		virtual void post_process(SmartPtr<vector_type> spU);

	///	returns new instance with same setting
		virtual SmartPtr<ITransferPostProcess<TAlgebra> > clone();

	protected:
		template <typename TBaseElem>
		void subtract_value(SmartPtr<GridFunction<TDomain, TAlgebra> > spGF, size_t fct, number sub);

	protected:
	///	symbolic function names
		std::string m_symbFct;

	///	indices of selected functions
		FunctionGroup m_fctGrp;

	///	approximation space
		SmartPtr<ApproximationSpace<TDomain> > m_spApproxSpace;

	///	fine grid level
		GridLevel m_level;

	///	initialization flag
		bool m_bInit;
};

} // end namespace ug

#include "transfer_post_process_impl.h"

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__TRANSFER_POST_PROCESS__ */
