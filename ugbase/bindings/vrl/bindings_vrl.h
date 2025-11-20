/*
 * Copyright (c) 2010-2011:  Steinbeis Forschungszentrum (STZ Ölbronn)
 * Author: Michael Hoffer
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


#ifndef BINDINGS_VRL_H
#define	BINDINGS_VRL_H

#include "registry/registry.h"
#include <jni.h>
#include "messaging.h"
#include "compile_info/compile_info.h"

#ifdef __GNUC__
#define PRETTY_FUNCTION __PRETTY_FUNCTION__
#else
#define PRETTY_FUNCTION "function name not available (not using GCC)"
#endif

namespace ug {
namespace vrl {
/**
 * Defines the registry to use for VRL.
 * @param pReg registry
 */
void SetVRLRegistry(ug::bridge::Registry* pReg);
/**
 * Defines the JVM to operate on based on the specified JNI environment.
 * The JVM can only be initialized once.
 * @param env JVM environment
 */
void initJavaVM(JNIEnv* env);

JavaVM* getJavaVM();

} // end vrl::
}// end ug::

#endif