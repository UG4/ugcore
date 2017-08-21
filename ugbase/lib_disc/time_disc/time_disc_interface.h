/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
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
 * This class uses a ISpatialDiscretization in order to implement the
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
			: m_spDomDisc(spDD)
		{}

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
	 *	A sub-routine at element-level ("prep_timestep_elem") is called
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
	///	the current solution elemDisc-wise
	/**
	 *	This function is called after the assembling routines
	 *	at the end of a time step.
	 *
	 * \param[in] currSol 	the solution at the previous time steps
	 * \param[in] dt		size of time step
	 */
		virtual void finish_step(SmartPtr<VectorTimeSeries<vector_type> > currSol) = 0;

	/// finishes a time step and allows to adapt data depending on
	///	the current solution element-wise
	/**
	 *	This function is called after the assembling routines at the end of a
	 *	time step.
	 *	Within this function "fsh_timestep_elem" is called which allows
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
		SmartPtr<IDomainDiscretization<TAlgebra> > m_spDomDisc; ///< Domain Discretization

	private:
	///	returns the assemble adapter
	/// \{
		SmartPtr<AssemblingTuner<TAlgebra> > ass_tuner() {return m_spDomDisc->ass_tuner();}
		ConstSmartPtr<AssemblingTuner<TAlgebra> > ass_tuner() const {return m_spDomDisc->ass_tuner();}
	/// \}
};

/// @}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__TIME_DISC__TIME_DISC_INTERFACE__ */
