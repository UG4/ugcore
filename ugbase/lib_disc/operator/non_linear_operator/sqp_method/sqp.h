/*
 * sqp.h
 *
 *  Created on: 11.01.2012
 *      Author: Raphael Prohl
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__SQP_METHOD__SQP__
#define __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__SQP_METHOD__SQP__

#include <cmath>

// other ug4 modules
#include "common/common.h"
#include "common/util/string_util.h"

// library intern headers
#include "lib_disc/spatial_disc/subset_assemble_util.h"
#include "lib_disc/spatial_disc/domain_disc_interface.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_assemble_util.h"
#include "lib_disc/spatial_disc/constraints/constraint_interface.h"
#include "lib_disc/spatial_disc/disc_item.h"

// for base class
#include "lib_disc/spatial_disc/domain_disc.h"

namespace ug {

template <typename TDomain, typename TAlgebra>
class SQPMethod : public DomainDiscretization<TDomain, TAlgebra>
{
	private:
	/// Base class type
		typedef DomainDiscretization<TDomain, TAlgebra>	base_type;

	public:
	//	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	///	Type of approximation space
		typedef ApproximationSpace<TDomain> approx_space_type;

	public:
		SQPMethod(SmartPtr<approx_space_type> spApproxSpace)
			: DomainDiscretization<TDomain, TAlgebra>(spApproxSpace),
			m_toleranceCheck(0.0)
			{};

		void set_tolerance_check(number TolCheck) {m_toleranceCheck = TolCheck;}

	// init Operator
		bool init();

	// prepare Operator
		bool prepare();

	// tolerance_check Operator: this operator checks,
	// if the linearized constraint is fulfilled sufficiently
		bool check_tolerance(const vector_type& u, const SurfaceDoFDistribution& dd, base_type& domDisc);

	// update Operator: this operator updates the SQP-variables
		bool update_variables(const vector_type& u, const SurfaceDoFDistribution& dd, base_type& domDisc);

		virtual ~SQPMethod(){};

	private:
	// 	Tolerance Check
		number m_toleranceCheck;
};

}

#include "sqp_impl.h"

#endif /* __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__SQP_METHOD__SQP__ */
