/*
 * sqp.h
 *
 *  Created on: 11.01.2012
 *      Author: Raphael Prohl
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__SQP_METHOD__SQP__
#define __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__SQP_METHOD__SQP__

#include <cmath>

// modul intern headers
#include "lib_disc/spatial_disc/domain_disc.h"
#include "lib_disc/common/groups_util.h"

namespace ug {

template <	typename TDomain,
			typename TDoFDistribution,
			typename TAlgebra>
class SQPMethod //:public DomainDiscretization<TDomain, TDoFDistribution, TAlgebra>
{
	public:
	///	Type of Domain
		typedef TDomain domain_type;

	//	Algebra type
		typedef TAlgebra algebra_type;

	//	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	//	DoFDistribution Type
		typedef TDoFDistribution dof_distribution_type;

	public:
		SQPMethod():
			m_tolerance_check(0.0)
			{};

		void set_tolerance_check(number TolCheck) {m_tolerance_check = TolCheck;}

	// init Operator
		virtual bool init();

	// prepare Operator
		virtual bool prepare();

	// tolerance_check Operator: this operator checks,
	// if the linearized constraint is fulfilled sufficiently
		virtual bool check_tolerance(const vector_type& u, const dof_distribution_type& dd);

	// update Operator: this operator updates the SQP-variables
		virtual bool update_variables(const vector_type& u, const dof_distribution_type& dd);

		//~SQPMethod();
		virtual ~SQPMethod(){};

	protected:
	// Tolerance Check
		number m_tolerance_check;

	///	vector holding all registered elem discs
		std::vector<IElemDisc*> m_vElemDisc;

	/// forces the assembling to regard the grid as regular
		bool m_bForceRegGrid;

};

}

#include "sqp_impl.h"

#endif /* __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__SQP_METHOD__SQP__ */
