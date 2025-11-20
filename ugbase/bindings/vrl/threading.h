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

#include <jni.h>
#include <stddef.h>
#include "messaging.h"
#include "bindings_vrl.h"

#ifndef THREADING_H
#define	THREADING_H
namespace ug {
	namespace vrl {
		namespace threading {

			/**
			 * Exception type
			 */
			enum ExceptionType {
				/**
				 * Attaching current thread to JVM failed.
				 */
				ATTACH_FAILED,
				/**
				 * Detaching current thread from JVM failed.
				 */
				DETACH_FAILED,
				/**
				 * Current thread not attached to JVM.
				 */
				NOT_ATTACHED
			};

			/**
			 * Exception
			 */
			class JNIThreadException {
			public:

				/**
				 * Constructor.
				 * @param type exception type
				 */
				JNIThreadException(ExceptionType type) {
					this->type = type;

					switch (type) {
						case ATTACH_FAILED:
						{
							UG_LOG("UG-VRL: Attaching thread failed in "
									<< EMPHASIZE_BEGIN
									<< PRETTY_FUNCTION
									<< EMPHASIZE_END
									<< " in line: " << __LINE__ << " !");
						}
							break;
						case DETACH_FAILED:
						{
							UG_LOG("UG-VRL: Detaching thread failed in "
									<< EMPHASIZE_BEGIN
									<< PRETTY_FUNCTION
									<< EMPHASIZE_END
									<< " in line: " << __LINE__ << " !");
						}
							break;
						case NOT_ATTACHED:
						{
							UG_LOG("UG-VRL: Thread not attached in "
									<< EMPHASIZE_BEGIN
									<< PRETTY_FUNCTION
									<< EMPHASIZE_END
									<< " in line: " << __LINE__ << " !");
						}
							break;
					}
				}

				ExceptionType type;
			};

			/**
			 * Attaches the current thread to the JVM. If the thread is already
			 * attached this is equivalent to <code>getEnv()</code>.
			 * @param javaVM Java VM to operate on
			 * @return JVM environment of the current thread
			 * @throws JNIThreadException
			 */
			inline JNIEnv* attachThread(JavaVM* javaVM) {

				// The following code raised a warning in newer GCC versions:
				// "dereferencing type-punned pointer will break strict-aliasing rules"
				// That is why we do it differently now, although this code
				// is officially used:
				//				JNIEnv* localEnv = NULL;
				//
				//				int result = javaVM->AttachCurrentThread(
				//						(void **) (&localEnv), NULL);

				JNIEnv** localEnvPtrPtr;
				JNIEnv* localEnv = NULL;
				localEnvPtrPtr = &localEnv;

				int result = javaVM->AttachCurrentThread(
						(void **) (localEnvPtrPtr), NULL);

				if (result < 0) {
					throw JNIThreadException(ATTACH_FAILED);
				}

				return localEnv;
			}

			/**
			 * Detaches the current thread from the JVM.
			 * @param javaVM Java VM to operate on
			 * @throws JNIThreadException
			 */
			inline void detachThread(JavaVM* javaVM) {

				int result = javaVM->DetachCurrentThread();

				if (result < 0) {
					throw JNIThreadException(DETACH_FAILED);
				}
			}

			/**
			 * Returns the JVM environment of the current thread.
			 * @param javaVM Java VM to operate on
			 * @return JVM environment of the current thread
			 * @throws JNIThreadException
			 */
			inline JNIEnv* getEnv(JavaVM* javaVM) {

				// The following code raised a warning in newer GCC versions:
				// "dereferencing type-punned pointer will break strict-aliasing rules"
				// That is why we do it differently now, although this code
				// is officially used:
				//				JNIEnv* localEnv = NULL;
				//
				//				jint result = javaVM->GetEnv(
				//						(void **) (&localEnv), JNI_VERSION_1_2);

				JNIEnv** localEnvPtrPtr;
				JNIEnv* localEnv = NULL;
				localEnvPtrPtr = &localEnv;

				jint result = javaVM->GetEnv(
						(void **) (localEnvPtrPtr), JNI_VERSION_1_2);

				if (result != JNI_OK) {
					throw JNIThreadException(NOT_ATTACHED);
				}

				return localEnv;
			}
		} // threading::
	} // vrl::
} // ug::

#endif
