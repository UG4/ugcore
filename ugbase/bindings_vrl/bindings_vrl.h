/*
 * File:   bindings_vrl.h
 * Author: Michael Hoffer <info@michaelhoffer.de>
 *
 * Created on 5. Oktober 2010, 14:54
 */

#ifndef BINDINGS_VRL_H
#define	BINDINGS_VRL_H

#include "ug_bridge/registry.h"
#include <jni.h>
#include "messaging.h"
#include "svnrevision.h"

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
//		JNIEnv* getJNIEnv();

		JavaVM* getJavaVM();

		/**
		 * Returns the global SVN revision of this build.
		 * @return global SVN revision of this build
		 */
		inline std::string svnRevision() {
			return split(SVN_REVISION, ':')[0]; //
		}
	} // end vrl::
}// end ug::

#endif	/*BINDINGS_VRL_H*/


//************************************
//Groovycode used to test UG4 bindings
//************************************
//
//@ComponentInfo(name="UG4 Test")
//class UG4Test implements Serializable {
//
//    private static final long serialVersionUID=1;
//
//    private transient VisualCanvas mainCanvas;
//
//    @MethodInfo(noGUI = true, callOptions = "assign-to-canvas")
//    public void setMainCanvas(VisualCanvas mainCanvas){
//        this.mainCanvas = mainCanvas;
//    }
//
//    @MethodInfo(valueStyle="editor")
//    public void init(){
//        def ug4 = edu.gcsc.vrl.ug4.UG4.getUG4()
//        String[] args = [""]
//        def codes = ug4.createJavaBindings()
//
//        codes.each{it->mainCanvas.addObjectAsCode(it)}
//
//    }
//}
