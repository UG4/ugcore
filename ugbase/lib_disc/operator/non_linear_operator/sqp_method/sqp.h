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

template <	typename TDomain,
			typename TDoFDistribution,
			typename TAlgebra>
class SQPMethod : public DomainDiscretization<TDomain, TDoFDistribution, TAlgebra>
{
	private:
	/// Base class type
		typedef DomainDiscretization<TDomain, TDoFDistribution, TAlgebra>
				base_type;
	public:
	///	Type of Domain Disc
	/*	typedef DomainDiscretization<TDomain, TDoFDistribution, TAlgebra>
						domain_disc_type;*/

	///	Type of Domain
		//typedef TDomain domain_type;

	//	Algebra type
	//	typedef TAlgebra algebra_type;

	//	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	//	DoFDistribution Type
		typedef IDoFDistribution<TDoFDistribution> dof_distribution_type;
		//typedef TDoFDistribution dof_distribution_type;

	///	Type of approximation space
		typedef ApproximationSpace<TDomain, TDoFDistribution, TAlgebra>
					approx_space_type;

	public:
		SQPMethod(approx_space_type& pApproxSpace): DomainDiscretization<TDomain, TDoFDistribution, TAlgebra>(pApproxSpace),
			m_toleranceCheck(0.0)
			{};
		/*SQPMethod():
			m_toleranceCheck(0.0)
			{};*/

		void set_tolerance_check(number TolCheck) {m_toleranceCheck = TolCheck;}

		//void set_discretization(domain_disc_type& domDisc) {m_domDisc = domDisc;}

	// init Operator
		//virtual bool init();
		bool init();

	// prepare Operator
		//virtual bool prepare();
		bool prepare();

	// tolerance_check Operator: this operator checks,
	// if the linearized constraint is fulfilled sufficiently
		//virtual bool check_tolerance(const vector_type& u, const dof_distribution_type& dd); // const domain_disc_type& domDisc);
		bool check_tolerance(const vector_type& u, const dof_distribution_type& dd, base_type& domDisc);

	// update Operator: this operator updates the SQP-variables
		//virtual bool update_variables(const vector_type& u, const dof_distribution_type& dd);
		bool update_variables(const vector_type& u, const dof_distribution_type& dd, base_type& domDisc);

		//~SQPMethod();
		virtual ~SQPMethod(){};

	/*protected:
	///	set the approximation space in the elem discs and extract IElemDiscs
		bool update_elem_discs();
		bool update_constraints();
		bool update_disc_items();*/


	private:
	// 	Tolerance Check
		number m_toleranceCheck;

	// 	domain disc
	//	domain_disc_type m_domDisc;

	///	vector holding all registered elem discs
	//	std::vector<IElemDisc*> m_vElemDisc;

	/// forces the assembling to regard the grid as regular
	//	bool m_bForceRegGrid;

	//	vector holding all registered constraints
	//	std::vector<IDomainConstraint<TDomain, TDoFDistribution, TAlgebra>*> m_vvConstraints[NUM_CONSTRAINT_TYPES];

	///	current approximation space
	//	approx_space_type* m_pApproxSpace;

	///	vector holding all registered elem discs
	//	std::vector<IDomainElemDisc<domain_type>*> m_vDomainElemDisc;

};

}

#include "sqp_impl.h"

#endif /* __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__SQP_METHOD__SQP__ */
