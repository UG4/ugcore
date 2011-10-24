/*
 * time_disc_interface.h
 *
 *  Created on: 23.10.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__TIME_DISC__TIME_DISC_INTERFACE__
#define __H__UG__LIB_DISC__TIME_DISC__TIME_DISC_INTERFACE__

// extern libraries
#include <deque>

// other ug libraries
#include "common/common.h"

// module intern libraries
#include "lib_disc/assemble_interface.h"
#include "lib_disc/time_disc/solution_time_series.h"
#include "lib_disc/spatial_disc/domain_disc_interface.h"

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
	///	DoF Distribution type
		typedef IDoFDistribution<TDoFDistribution> dof_distribution_type;

	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Vector type
		typedef typename algebra_type::vector_type vector_type;

	///	Domain Discretization type
		typedef IDomainDiscretization<TDoFDistribution, algebra_type>
			domain_discretization_type;

	public:
	/// create and set domain discretization
	/**
	 * \param[in] 	dd	Domain Discretization
	 */
		ITimeDiscretization(domain_discretization_type& dd)
			: m_rDomDisc(dd)
		{}

	/// prepares the assembling of Defect/Jacobian for a time step
	/**
	 *	This function supplies the TimeDiscretization with previous time
	 *	steps and step size before the assembling routines can be called.
	 *
	 * \param[in] u_old 	the solution at the previous time steps
	 * \param[in] time_old	the time at the previous time steps
	 * \param[in] dt		size of time step
	 */
		virtual void prepare_step(VectorTimeSeries<vector_type>& prevSol,
		                          number dt) = 0;

	/// returns number of previous time steps needed
		virtual size_t num_prev_steps() = 0;

	/// forces the assembling to consider the grid as regular
		virtual void force_regular_grid(bool bForce)
		{
			m_rDomDisc.force_regular_grid(bForce);
		}

	///	sets a selector to exlude elements from assembling
	/**
	 * This methods sets a selector. Only elements that are selected will be
	 * assembled during assembling process. If no selector is set, this
	 * corresponds to a selector where all elements have been selected.
	 *
	 * \param[in]	sel		Selector
	 */
		virtual void set_selector(ISelector* sel = NULL)
		{
			m_pSelector = sel; forward_selector();
		}

	protected:
		void forward_selector()
		{
			m_rDomDisc.set_selector(m_pSelector);
		}

		domain_discretization_type& m_rDomDisc; ///< Domain Discretization

		ISelector* m_pSelector;
};

/// @}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__TIME_DISC__TIME_DISC_INTERFACE__ */
