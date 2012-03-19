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
 *
 * \tparam	TAlgebra			Algebra Type
 */
template <typename TAlgebra>
class ITimeDiscretization : public IAssemble<TAlgebra>
{
	public:
	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Vector type
		typedef typename algebra_type::vector_type vector_type;

	///	Domain Discretization type
		typedef IDomainDiscretization<algebra_type>	domain_discretization_type;

	public:
	/// create and set domain discretization
	/**
	 * \param[in] dd	Domain Discretization
	 */
		ITimeDiscretization(domain_discretization_type& dd)
			: m_rDomDisc(dd)
		{}

	/// prepares the assembling of Defect/Jacobian for a time step
	/**
	 *	This function supplies the TimeDiscretization with previous time
	 *	steps and step size before the assembling routines can be called.
	 *
	 * \param[in] prevSol 	the solution at the previous time steps
	 * \param[in] dt		size of time step
	 */
		virtual void prepare_step(VectorTimeSeries<vector_type>& prevSol,
		                          number dt) = 0;

	/// prepares the assembling of Defect/Jacobian for a time step
	/**
	 *	This function supplies the TimeDiscretization with previous time
	 *	steps and step size before the assembling routines can be called.
	 *	A sub-routine at element-level ("prepare_timestep_element") is called
	 *	within this function.
	 *
	 * \param[in] prevSol 	the solution at the previous time steps
	 * \param[in] dt		size of time step
	 * \param[in] dd		DoF Distribution
	 */
	/// \{
		virtual void prepare_step_elem(VectorTimeSeries<vector_type>& prevSol,
		                               number dt, GridLevel gl) = 0;
		void prepare_step_elem(VectorTimeSeries<vector_type>& prevSol,
		                       number dt)
		{prepare_step_elem(prevSol, dt, GridLevel());}
	/// \}

	/// finishes the assembling of Defect/Jacobian for a time step
	/**
	 *	This function supplies the TimeDiscretization with previous time
	 *	steps and step size after the assembling routines have been called.
	 *	A sub-routine at element-level ("finish_timestep_element") is called
	 *	within this function.
	 *
	 * \param[in] prevSol 	the solution at the previous time steps
	 * \param[in] dt		size of time step
	 * \param[in] dd		DoF Distribution
	 */
	///	\{
		virtual void finish_step_elem(VectorTimeSeries<vector_type>& prevSol,
									  number dt, GridLevel gl) = 0;
		void finish_step_elem(VectorTimeSeries<vector_type>& prevSol,
		                       number dt)
		{finish_step_elem(prevSol, dt, GridLevel());}
	///	\}

	///	returns the future time point (i.e. the one that will be computed)
		virtual number future_time() const = 0;

	/// returns number of previous time steps needed
		virtual size_t num_prev_steps() const = 0;

	///	returns the number of stages
		virtual size_t num_stages() const = 0;

	///	sets the stage
		virtual void set_stage(size_t stage) = 0;

	/// forces the assembling to consider the grid as regular
		virtual void force_regular_grid(bool bForce)
		{
			m_rDomDisc.force_regular_grid(bForce);
		}

	///	enables constraints
		virtual void enable_constraints(bool bEnable)
		{
			m_rDomDisc.enable_constraints(bEnable);
		}

	///	sets a selector to exclude elements from assembling
	/**
	 * This methods sets a selector. Only elements that are selected will be
	 * assembled during assembling process. If no selector is set, this
	 * corresponds to a selector where all elements have been selected.
	 *
	 * \param[in]	sel		Selector
	 */
		virtual void set_selector(BoolMarker* sel = NULL)
		{
			m_pBoolMarker = sel; forward_selector();
		}

	///	returns the number of post processes
		virtual size_t num_dirichlet_constraints() const
		{
			return m_rDomDisc.num_dirichlet_constraints();
		}

	///	returns the i'th post process
		virtual SmartPtr<IConstraint<TAlgebra> > dirichlet_constraint(size_t i)
		{
			return m_rDomDisc.dirichlet_constraint(i);
		}

	protected:
		void forward_selector()
		{
			m_rDomDisc.set_selector(m_pBoolMarker);
		}

		domain_discretization_type& m_rDomDisc; ///< Domain Discretization

		BoolMarker* m_pBoolMarker;
};

/// @}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__TIME_DISC__TIME_DISC_INTERFACE__ */
