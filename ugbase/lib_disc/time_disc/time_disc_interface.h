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
#include "lib_disc/spatial_disc/ass_adapter.h"

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
		typedef IDomainDiscretization<TAlgebra>	domain_discretization_type;

	public:
	/// create and set domain discretization
	/**
	 * \param[in] dd	Domain Discretization
	 */
		ITimeDiscretization(SmartPtr<IDomainDiscretization<TAlgebra> > spDD)
			: m_spDomDisc(spDD), m_pAssAdapter(NULL)
		{
			m_pAssAdapter = &get_ass_adapter();
		}

	/// prepares the assembling of Defect/Jacobian for a time step
	/**
	 *	This function supplies the TimeDiscretization with previous time
	 *	steps and step size before the assembling routines can be called.
	 *
	 * \param[in] prevSol 	the solution at the previous time steps
	 * \param[in] dt		size of time step
	 */
		virtual void prepare_step(SmartPtr<VectorTimeSeries<vector_type> > prevSol,
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
		virtual void prepare_step_elem(SmartPtr<VectorTimeSeries<vector_type> > prevSol,
		                               number dt, const GridLevel& gl) = 0;
		void prepare_step_elem(SmartPtr<VectorTimeSeries<vector_type> > prevSol,
		                       number dt)
		{prepare_step_elem(prevSol, dt, GridLevel());}
	/// \}

	/// finishes a time step and allows to adapt data depending on
	///	the current solution elementwise
	/**
	 *	This function is called after the assembling routines at the end of a
	 *	timestep.
	 *	Within this function "finish_timestep_element" is called which allows
	 *	modifying data depending on the current solution at element-level.
	 *
	 * \param[in] currSol 	the current solution
	 * \param[in] dd		DoF Distribution
	 */
	///	\{
		virtual void finish_step_elem(SmartPtr<VectorTimeSeries<vector_type> > currSol,
									  const GridLevel& gl) = 0;
		void finish_step_elem(SmartPtr<VectorTimeSeries<vector_type> > currSol)
		{finish_step_elem(currSol, GridLevel());}
	///	\}

	///	returns the future time point (i.e. the one that will be computed)
		virtual number future_time() const = 0;

	/// returns number of previous time steps needed
		virtual size_t num_prev_steps() const = 0;

	///	returns the number of stages
		virtual size_t num_stages() const = 0;

	///	sets the stage
		virtual void set_stage(size_t stage) = 0;

	///	returns the number of constraint
		virtual size_t num_constraints() const
		{
			return m_spDomDisc->num_constraints();
		}

	///	returns the i'th constraint
		virtual SmartPtr<IConstraint<TAlgebra> > constraint(size_t i)
		{
			return m_spDomDisc->constraint(i);
		}

	protected:
		AssAdapter<TAlgebra>& get_ass_adapter()
		{
			return m_spDomDisc->get_ass_adapter();
		}

		SmartPtr<IDomainDiscretization<TAlgebra> > m_spDomDisc; ///< Domain Discretization

		///	this object provides tools to adapt the assemble routine
		AssAdapter<TAlgebra>* m_pAssAdapter;
};

/// @}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__TIME_DISC__TIME_DISC_INTERFACE__ */
