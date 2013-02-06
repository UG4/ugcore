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
	 *	A sub-routine at element-level ("prepare_timestep_element") is called
	 *	within this function.
	 *
	 * \param[in] prevSol 	the solution at the previous time steps
	 * \param[in] dt		size of time step
	 * \param[in] dd		DoF Distribution
	 */
	/// \{
		virtual void prepare_step_elem(SmartPtr<VectorTimeSeries<vector_type> > prevSol,
		                               number dt, GridLevel gl) = 0;
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
									  GridLevel gl) = 0;
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

	/// forces the assembling to consider the grid as regular
		virtual void force_regular_grid(bool bForce)
		{
			m_spDomDisc->force_regular_grid(bForce);
		}

	///	returns type of constraints enabled
		virtual int constraints_enabled() const
		{
			return m_spDomDisc->constraints_enabled();
		}

	///	enables constraints
		virtual void enable_constraints(int TypesEnable)
		{
			m_spDomDisc->enable_constraints(TypesEnable);
		}

	///	returns type of boundary elem discs enabled
		virtual int elem_discs_enabled() const
		{
			return m_spDomDisc->elem_discs_enabled();
		}

	///	enables boundary elem discs
		virtual void enable_elem_discs(int TypesEnable)
		{
			m_spDomDisc->enable_elem_discs(TypesEnable);
		}

	///	sets a marker to exclude elements from assembling
	/**
	 * This methods sets a marker. Only elements that are marked will be
	 * assembled during assembling process. If no marker is set, this
	 * corresponds to a marker where all elements have been marked.
	 *
	 * \param[in]	mark	BoolMarker
	 */
		virtual void set_marker(BoolMarker* mark = NULL)
		{
			m_AssAdapter.pBoolMarker = mark; forward_marker();
		}
	///	sets a selector of elements for assembling
	/**
	 * This methods sets an element list. Only elements of this list will be
	 * assembled during assembling process. Especially the list defines the begin
	 * and end of the element-iterator in the element assembling-loop.
	 * If no element list is set, this corresponds to a assembling where the loop is
	 * carried out over all elements of a subset.
	 *
	 * \param[in]	sel		Selector
	 */
		virtual void set_selector(Selector* sel = NULL)
		{
			m_AssAdapter.pSelector = sel; forward_selector();
		}
	
		virtual void ass_index(size_t ind, bool index_set = true)
		{
			m_AssAdapter.assIndex.index = ind; m_AssAdapter.assIndex.index_set = index_set; forward_ass_index();
		}
		
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
		void forward_marker()
		{
			m_spDomDisc->set_marker(m_AssAdapter.pBoolMarker);
		}

		void forward_selector()
		{
			m_spDomDisc->set_selector(m_AssAdapter.pSelector);
		}

		void forward_ass_index()
		{
			m_spDomDisc->ass_index(m_AssAdapter.assIndex.index, m_AssAdapter.assIndex.index_set);
		}
		
		SmartPtr<IDomainDiscretization<TAlgebra> > m_spDomDisc; ///< Domain Discretization

		AssAdapter m_AssAdapter;
};

/// @}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__TIME_DISC__TIME_DISC_INTERFACE__ */
