/*
 * File:   bindings_vrl.h
 * Author: Michael Hoffer <info@michaelhoffer.de>
 *
 * Created on 5. Oktober 2010, 14:54
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

/**
 * Returns the global SVN revision of this build.
 * @return global SVN revision of this build
 */
inline std::string svnRevision() {
	return split(UGSvnRevision(), ':')[0]; //
}
} // end vrl::
}// end ug::

#endif	/*BINDINGS_VRL_H*/
