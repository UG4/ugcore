/*
 * Copyright (c) 2010-2020:  G-CSC, Goethe University Frankfurt
 * Author: Tim Schön
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

#ifndef __H__UG__LIB_DISC__TIME_DISC__TIME_INTEGRATOR_SUBJECT
#define __H__UG__LIB_DISC__TIME_DISC__TIME_INTEGRATOR_SUBJECT


namespace ug {

/// Base class for a subject notifying observers attachment.
/** Provides the option to perform pre-/postprocessing for a (tentative step for evolving from t -> t+dt). Seven cases are distinguished
 *  1) INIT (STEP)		  : Called at t=t_i before time step is executed. Contains solution u0=u(t0).
 *  2) FINALIZE (STEP)	  : Called at t=t_i+dt, after time step has been executed and can be accepted. Provides solution u = u(t+dt).
 *  3) REWIND (STEP)	  : Called at t=t_i+dt, after time step has been executed, but must be rejected. Provides (rejected) solution u = u(t+dt).
 *  4) PREPROCESS (STEP)  : Called at t=t_i before newton solver is executed (if applicable). Contains solution u0=u(t0). This may happen multiple times per time step.
 *  5) POSTPROCESS (STEP) : Called at t=t_i+dt, after newton solver has been executed (if applicable). Provides solution u = u(t+dt). This may happen multiple times per time step.
 *  6) START  			  : Called at t=t_0, at the beginning of the time integration process (start of the simulation)
 *  7) END      		  : Called at t=T, at the end of the time integration process (end of the simulation)
 */
template<class TDomain, class TAlgebra>
class TimeIntegratorSubject
{
public:
	typedef GridFunction<TDomain, TAlgebra> grid_function_type;
	typedef ITimeIntegratorObserver<TDomain, TAlgebra> process_observer_type;
	typedef typename std::vector<SmartPtr<process_observer_type> > process_observer_container_type;

	enum observer_group_type {
		TIO_GROUP_INIT_STEP=0, 
		TIO_GROUP_REWIND_STEP, 
		TIO_GROUP_FINALIZE_STEP, 
		TIO_GROUP_PREPROCESS_STEP, 
		TIO_GROUP_POSTPROCESS_STEP, 
		TIO_GROUP_START, 
		TIO_GROUP_END, 
		TIO_GROUP_SIZE};

protected:
	// process_observer_container_type m_vProcessObservers;
	process_observer_container_type m_vProcessObservers[TIO_GROUP_SIZE];

protected:
	//! register observer (default: postprocess)
	template<int tGroup>
	void attach_to_group(SmartPtr<process_observer_type> obs)
	{
		//UG_LOG("TimeIntegratorSubject::attach_observer[" << tGroup <<"]" << this << std::endl);
		m_vProcessObservers[tGroup].push_back(obs);
	}

	//! register observer (default: postprocess)
	void attach_to_group(int tGroup, SmartPtr<process_observer_type> obs)
	{
		//UG_LOG("TimeIntegratorSubject::attach_observer[" << tGroup <<"]" << this << std::endl);
		m_vProcessObservers[tGroup].push_back(obs);
	}
public:
	//! Short-cut for
	void attach_observer(SmartPtr<process_observer_type> obs)
	{ attach_finalize_observer(obs); }

	void attach_init_observer(SmartPtr<process_observer_type> obs)
	{ attach_to_group<TIO_GROUP_INIT_STEP>(obs); }

	void attach_rewind_observer(SmartPtr<process_observer_type> obs)
	{ attach_to_group<TIO_GROUP_REWIND_STEP>(obs); }

	void attach_finalize_observer(SmartPtr<process_observer_type> obs)
	{ attach_to_group<TIO_GROUP_FINALIZE_STEP>(obs); }

	void attach_preprocess_observer(SmartPtr<process_observer_type> obs)
	{ attach_to_group<TIO_GROUP_PREPROCESS_STEP>(obs); }

	void attach_postprocess_observer(SmartPtr<process_observer_type> obs)
	{ attach_to_group<TIO_GROUP_POSTPROCESS_STEP>(obs); }

	void attach_start_observer(SmartPtr<process_observer_type> obs)
	{ attach_to_group<TIO_GROUP_START>(obs); }

	void attach_end_observer(SmartPtr<process_observer_type> obs)
	{ attach_to_group<TIO_GROUP_END>(obs); }

	void reset_observers()
	{
		m_vProcessObservers[TIO_GROUP_INIT_STEP].clear();
		m_vProcessObservers[TIO_GROUP_REWIND_STEP].clear();
		m_vProcessObservers[TIO_GROUP_FINALIZE_STEP].clear();
		m_vProcessObservers[TIO_GROUP_PREPROCESS_STEP].clear();
		m_vProcessObservers[TIO_GROUP_POSTPROCESS_STEP].clear();
		m_vProcessObservers[TIO_GROUP_START].clear();
		m_vProcessObservers[TIO_GROUP_END].clear();
	}

protected:
	/// Notify all observers for a certain group.
	template<int tGroup>
	void notify_group(SmartPtr<grid_function_type> u, int step, number time, number dt)
	{
		process_observer_container_type &observers = m_vProcessObservers[tGroup];
		for (typename process_observer_container_type::iterator it = observers.begin(); it!= observers.end(); ++it)
		{(*it)->step_process(u, step, time, dt); }
		//UG_LOG("TimeIntegratorSubject::notify_group[" << tGroup <<"]: " << observers.size() << " observers. " << this<< std::endl);
	}
public:

	/// notify all observers that time step evolution starts
	void notify_init_step(SmartPtr<grid_function_type> u, int step, number time, number dt)
	{ notify_group<TIO_GROUP_INIT_STEP>(u, step, time, dt); }

	/// Notify all observers that time step must be rewinded.
	void notify_rewind_step(SmartPtr<grid_function_type> u, int step, number time, number dt)
	{ notify_group<TIO_GROUP_REWIND_STEP>(u, step, time, dt); }

	/// notify all observers that time step has been evolved (successfully)
	void notify_finalize_step(SmartPtr<grid_function_type> u, int step, number time, number dt)
	{ notify_group<TIO_GROUP_FINALIZE_STEP>(u, step, time, dt); }

	/// notify all observers that newton solver is about to start (may happen multiple times per time step)
	void notify_preprocess_step(SmartPtr<grid_function_type> u, int step, number time, number dt)
	{ notify_group<TIO_GROUP_PREPROCESS_STEP>(u, step, time, dt); }

	/// notify all observers that newton solver has finished (may happen multiple times per time step)
	void notify_postprocess_step(SmartPtr<grid_function_type> u, int step, number time, number dt)
	{ notify_group<TIO_GROUP_POSTPROCESS_STEP>(u, step, time, dt); }

	/// notify all observers that the simulation has started
	void notify_start(SmartPtr<grid_function_type> u, int step, number time, number dt)
	{ notify_group<TIO_GROUP_START>(u, step, time, dt); }

	/// notify all observers that the simulation has ended
	void notify_end(SmartPtr<grid_function_type> u, int step, number time, number dt)
	{ notify_group<TIO_GROUP_END>(u, step, time, dt); }



};

}

#endif