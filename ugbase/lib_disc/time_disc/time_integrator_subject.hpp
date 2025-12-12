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

#include "time_integrator_observers/time_integrator_observer_interface.h"

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
template<typename TDomain, typename TAlgebra>
class TimeIntegratorSubject
{
public:
	using grid_function_type = GridFunction<TDomain, TAlgebra>;
	using process_observer_type = ITimeIntegratorObserver<TDomain, TAlgebra>;
	using init_observer_type = ITimeIntegratorStageObserver_init<TDomain, TAlgebra>;
	using rewind_observer_type = ITimeIntegratorStageObserver_rewind<TDomain, TAlgebra>;
	using finalize_observer_type = ITimeIntegratorStageObserver_finalize<TDomain, TAlgebra>;
	using preprocess_observer_type = ITimeIntegratorStageObserver_preprocess<TDomain, TAlgebra>;
	using postprocess_observer_type = ITimeIntegratorStageObserver_postprocess<TDomain, TAlgebra>;
	using start_observer_type = ITimeIntegratorStageObserver_start<TDomain, TAlgebra>;
	using end_observer_type = ITimeIntegratorStageObserver_end<TDomain, TAlgebra>;
	using process_observer_container_type = std::vector<SmartPtr<process_observer_type> >;

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

	// container for statically mapped observers
	std::vector<SmartPtr<init_observer_type> > m_vInitObservers;
	std::vector<SmartPtr<rewind_observer_type> > m_vRewindObservers;
	std::vector<SmartPtr<finalize_observer_type> > m_vFinalizeObservers;
	std::vector<SmartPtr<preprocess_observer_type> > m_vPreprocessObservers;
	std::vector<SmartPtr<postprocess_observer_type> > m_vPostprocessObservers;
	std::vector<SmartPtr<start_observer_type> > m_vStartObservers;
	std::vector<SmartPtr<end_observer_type> > m_vEndObservers;


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

#define DECLARE_CHECK_STATIC_ATTACH(stage, container) \
	bool check_attach_##stage(SmartPtr<process_observer_type> obs)\
	{\
		SmartPtr<stage##_observer_type> sp_staticObs = obs.template cast_dynamic<stage##_observer_type>();\
		if (sp_staticObs.valid())\
		{\
			(container).push_back(sp_staticObs);\
			return true;\
		}\
		return false;\
	}

	bool check_attach_init(SmartPtr<process_observer_type> obs) {
		SmartPtr<init_observer_type> sp_staticObs = obs.template cast_dynamic<init_observer_type>();
		if (sp_staticObs.valid()) {
			(m_vInitObservers).push_back(sp_staticObs);
			return true;
		}
		return false;
	}

	bool check_attach_rewind(SmartPtr<process_observer_type> obs) {
		SmartPtr<rewind_observer_type> sp_staticObs = obs.template cast_dynamic<rewind_observer_type>();
		if (sp_staticObs.valid()) {
			(m_vRewindObservers).push_back(sp_staticObs);
			return true;
		}
		return false;
	}
	bool check_attach_finalize(SmartPtr<process_observer_type> obs) {
		SmartPtr<finalize_observer_type> sp_staticObs = obs.template cast_dynamic<finalize_observer_type>();
		if (sp_staticObs.valid()) {
			(m_vFinalizeObservers).push_back(sp_staticObs);
			return true;
		}
		return false;
	}
	bool check_attach_preprocess(SmartPtr<process_observer_type> obs) {
		SmartPtr<preprocess_observer_type> sp_staticObs = obs.template cast_dynamic<preprocess_observer_type>();
		if (sp_staticObs.valid()) {
			(m_vPreprocessObservers).push_back(sp_staticObs);
			return true; } return false;
	}
	bool check_attach_postprocess(SmartPtr<process_observer_type> obs) {
		SmartPtr<postprocess_observer_type> sp_staticObs = obs.template cast_dynamic<postprocess_observer_type>();
		if (sp_staticObs.valid()) {
			(m_vPostprocessObservers).push_back(sp_staticObs);
			return true;
		}
		return false;
	}
	bool check_attach_start(SmartPtr<process_observer_type> obs) {
		SmartPtr<start_observer_type> sp_staticObs = obs.template cast_dynamic<start_observer_type>();
		if (sp_staticObs.valid()) {
			(m_vStartObservers).push_back(sp_staticObs);
			return true;
		}
		return false;
	}
	bool check_attach_end(SmartPtr<process_observer_type> obs) {
		SmartPtr<end_observer_type> sp_staticObs = obs.template cast_dynamic<end_observer_type>();
		if (sp_staticObs.valid()) {
			(m_vEndObservers).push_back(sp_staticObs);
			return true;
		} return false;
	}

	//! Attach statically mapped observers to their respective stages.
	//! Other observers are mapped to the finalize stage.
	void attach_observer(SmartPtr<process_observer_type> obs)
	{
		const bool isStaticallyAttached =
			check_attach_init(obs) ||
			check_attach_rewind(obs) ||
			check_attach_finalize(obs) ||
			check_attach_preprocess(obs) ||
			check_attach_postprocess(obs) ||
			check_attach_start(obs) ||
			check_attach_end(obs);

		if (!isStaticallyAttached)
			attach_finalize_observer(obs);
	}

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

		m_vInitObservers.clear();
		m_vRewindObservers.clear();
		m_vFinalizeObservers.clear();
		m_vPreprocessObservers.clear();
		m_vPostprocessObservers.clear();
		m_vStartObservers.clear();
		m_vEndObservers.clear();
	}

protected:
	/// Notify all observers for a certain group.
	template<int tGroup>
	bool notify_group(SmartPtr<grid_function_type> u, int step, number time, number dt)
	{
		process_observer_container_type &observers = m_vProcessObservers[tGroup];
		
		bool result = true;
		for (typename process_observer_container_type::iterator it = observers.begin(); it!= observers.end(); ++it)
		{
			result = (*it)->step_process(u, step, time, dt) && result;				 			
		}

		return result;
	}


public:

#define DECLARE_NOTIFY_STEP(functionName, stageName, stageID, container) \
	bool notify_##functionName(SmartPtr<grid_function_type> u, int step, number time, number dt)\
	{\
		bool res = notify_group<(stageID)>(u, step, time, dt);\
		const size_t numObs = (container).size();\
		for (size_t o = 0; o < numObs; ++o)\
			res &= (container)[o]->stageName##_action(u, step, time, dt);\
		return res;\
	}

	/// notify all observers that time step evolution starts
	bool notify_init_step(SmartPtr<grid_function_type> u, int step, number time, number dt) {
		bool res = notify_group<(TIO_GROUP_INIT_STEP)>(u, step, time, dt);
		const size_t numObs = (m_vInitObservers).size();
		for (size_t o = 0; o < numObs; ++o)
			res &= (m_vInitObservers)[o]->init_action(u, step, time, dt);
		return res;
	}

	/// Notify all observers that time step must be rewinded.
	bool notify_rewind_step(SmartPtr<grid_function_type> u, int step, number time, number dt) {
		bool res = notify_group<(TIO_GROUP_REWIND_STEP)>(u, step, time, dt);
		const size_t numObs = (m_vRewindObservers).size();
		for (size_t o = 0; o < numObs; ++o)
			res &= (m_vRewindObservers)[o]->rewind_action(u, step, time, dt);
		return res;
	}

	/// notify all observers that time step has been evolved (successfully)
	bool notify_finalize_step(SmartPtr<grid_function_type> u, int step, number time, number dt) {
		bool res = notify_group<(TIO_GROUP_FINALIZE_STEP)>(u, step, time, dt);
		const size_t numObs = (m_vFinalizeObservers).size();
		for (size_t o = 0; o < numObs; ++o)
			res &= (m_vFinalizeObservers)[o]->finalize_action(u, step, time, dt);
		return res;
	}

	/// notify all observers that newton solver is about to start (may happen multiple times per time step)
	bool notify_preprocess_step(SmartPtr<grid_function_type> u, int step, number time, number dt) {
		bool res = notify_group<(TIO_GROUP_PREPROCESS_STEP)>(u, step, time, dt);
		const size_t numObs = (m_vPreprocessObservers).size();
		for (size_t o = 0; o < numObs; ++o)
			res &= (m_vPreprocessObservers)[o]->preprocess_action(u, step, time, dt);
		return res;
	}

	/// notify all observers that newton solver has finished (may happen multiple times per time step)
	bool notify_postprocess_step(SmartPtr<grid_function_type> u, int step, number time, number dt) {
		bool res = notify_group<(TIO_GROUP_POSTPROCESS_STEP)>(u, step, time, dt);
		const size_t numObs = (m_vPostprocessObservers).size();
		for (size_t o = 0; o < numObs; ++o)
			res &= (m_vPostprocessObservers)[o]->postprocess_action(u, step, time, dt);
		return res;
	}

	/// notify all observers that the simulation has started
	bool notify_start(SmartPtr<grid_function_type> u, int step, number time, number dt) {
		bool res = notify_group<(TIO_GROUP_START)>(u, step, time, dt);
		const size_t numObs = (m_vStartObservers).size();
		for (size_t o = 0; o < numObs; ++o)
			res &= (m_vStartObservers)[o]->start_action(u, step, time, dt);
		return res;
	}

	/// notify all observers that the simulation has ended
	bool notify_end(SmartPtr<grid_function_type> u, int step, number time, number dt) {
		bool res = notify_group<(TIO_GROUP_END)>(u, step, time, dt);
		const size_t numObs = (m_vEndObservers).size();
		for (size_t o = 0; o < numObs; ++o)
			res &= (m_vEndObservers)[o]->end_action(u, step, time, dt); return res;
	}
};

}

#endif
