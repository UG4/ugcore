/*
 * time_discretization_interface.h
 *
 *  Created on: 23.10.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__TIME_DISCRETIZATION__TIME_DISCRETIZATION_INTERFACE__
#define __H__UG__LIB_DISCRETIZATION__TIME_DISCRETIZATION__TIME_DISCRETIZATION_INTERFACE__

// extern libraries
#include <deque>

// other ug libraries
#include "common/common.h"

// module intern libraries
#include "lib_discretization/assemble_interface.h"
#include "lib_discretization/time_discretization/previous_solution.h"
#include "lib_discretization/spatial_discretization/domain_discretization_interface.h"

namespace ug{

/// \ingroup lib_disc_time_assemble
/// @{

/// Time Discretization Interface
/**
 * Defines the time discretization interface.
 *
 * This class uses a ISpatialDiscratization in order to implement the
 * IAssemble interface.
 *
 * After the method prepare step has been called, Jacobian/Defect can be computed.
 * \tparam 	TDoFDistribution	DoF Distribution Type
 * \tparam	TAlgebra			Algebra Type
 */
template <	typename TDoFDistribution,
			typename TAlgebra>
class ITimeDiscretization
	: public IAssemble<TDoFDistribution, TAlgebra>
{
	public:
	//	DoF Distribution type
		typedef IDoFDistribution<TDoFDistribution> dof_distribution_type;

	//	Algebra type
		typedef TAlgebra algebra_type;

	//	Vector type
		typedef typename algebra_type::vector_type vector_type;

	// 	Domain Discretization type
		typedef IDomainDiscretization<TDoFDistribution, algebra_type>
			domain_discretization_type;

	public:
	/// Default constructor
	/**
	 * Creates an empty Time Discretization. In order to use this class a
	 * spatial Discretization has to be set.
	 */
		ITimeDiscretization()
			: m_pDomDisc(NULL)	{}

	/// create and set domain discretization
	/**
	 * \param[in] 	dd	Domain Discretization
	 */
		ITimeDiscretization(domain_discretization_type& dd)
			: m_pDomDisc(&dd)
		{}

	///	set the domain discretization
		void set_domain_discretization(domain_discretization_type& dd)
		{
			m_pDomDisc = &dd;
		}

	///	get the domain discretization
		domain_discretization_type* get_domain_discretization()
		{
			return m_pDomDisc;
		}

	/// prepares the assembling of Defect/Jacobian for a time step
	/**
	 *	This function supplies the TimeDiscretization with previous time
	 *	steps and step size before the assembling routines can be called.
	 *
	 * \param[in] u_old 	the solution at the previous time steps
	 * \param[in] time_old	the time at the previous time steps
	 * \param[in] dt		size of time step
	 */
		virtual bool prepare_step(const PreviousSolutions<vector_type>& prevSol,
		                          number dt) = 0;

	/// returns number of previous time steps needed
		virtual size_t num_prev_steps() = 0;

	/// forces the assembling to consider the grid as regular
		virtual void force_regular_grid(bool bForce)
		{
			if(m_pDomDisc != NULL)
				m_pDomDisc->force_regular_grid(bForce);
		}

	protected:
		domain_discretization_type* m_pDomDisc; ///< Domain Discretization
};

/// @}

} // end namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__TIME_DISCRETIZATION__TIME_DISCRETIZATION_INTERFACE__ */
