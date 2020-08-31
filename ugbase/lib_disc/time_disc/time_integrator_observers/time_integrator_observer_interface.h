/*
 * Copyright (c) 2010-2020:  G-CSC, Goethe University Frankfurt
 * Author: Arne Naegel, Andreas Kreienbuehl, Tim Schön
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

#ifndef __H__UG__LIB_DISC__TIME_DISC__TIME_INTEGRATOR_OBSERVERS__TIME_INTEGRATOR_OBSERVER_INTERFACE
#define __H__UG__LIB_DISC__TIME_DISC__TIME_INTEGRATOR_OBSERVERS__TIME_INTEGRATOR_OBSERVER_INTERFACE

#include <registry/class.h>

namespace ug {

/// Abstract base class for time integration observer
template<class TDomain, class TAlgebra>
class ITimeIntegratorObserver
{
public:
	typedef GridFunction<TDomain, TAlgebra> grid_function_type;

	virtual ~ITimeIntegratorObserver() {}
	virtual bool step_process(SmartPtr<grid_function_type> u, int step, number time, number dt)  = 0 ;
};



// template<class TDomain, class TAlgebra>
// class ProxyTimeIntegratorObserver : public ITimeIntegratorObserver<TDomain, TAlgebra>
// {
// 	typedef GridFunction<TDomain, TAlgebra> grid_function_type;
// 	typedef bool (*ProxyObserverFunction)(SmartPtr<grid_function_type> u, int step, number time, number dt);

//     public:
//         ProxyTimeIntegratorObserver(ProxyObserverFunction func, SmartPtr<void> obj) : m_method(func), m_obj(obj) {}

//         bool step_process(SmartPtr<grid_function_type> u, int step, number time, number dt)
//         {
// 			ProxyObserverFunction fp = (ProxyObserverFunction)m_method.get_raw_ptr();
//             return m_obj.get()->*fp(u, step, time, dt);
//         }

//     private:
//         bridge::MethodPtrWrapper m_method;
// 		SmartPtr<void> m_obj;
// };


}

#endif