

#include "../ug_bridge/ug_bridge.h"
#include "type_converter.h"
#include "common/common.h"
#include "../lib_discretization/spacial_discretization/user_data.h"
#include "../ug_script/ug_script.h"
#include "bindings_vrl.h"
#include "type_converter.h"
#include "threading.h"
#include <iostream>
#include <sstream>

namespace ug {
	namespace vrl {

		template <int dim>
		class VRLUserNumber : public IUserNumber<dim> {
		public:
			//	Functor Type
			typedef typename IUserNumber<dim>::functor_type functor_type;

			//	return functor

			virtual functor_type get_functor() const {
				return boost::ref(*this);
			}

		public:

			VRLUserNumber() {
				javaVM = getJavaVM();
				expression="";
			}

			void set_vrl_callback(const char* expression) {
				this->expression = expression;

//				JNIEnv* localEnv = NULL;
//				javaVM->AttachCurrentThread((void **) (&localEnv), NULL);
//
//
//				javaVM->DetachCurrentThread();
			}

			void operator() (number& c, const MathVector<dim>& x, number time = 0.0) const {

				//			for(size_t i = 0; i < /*dim*/ 2; ++i) {
				//
				//			}

				//			assert(dim==2);

				double val_x = x[0];
				double val_y = 0;

				if (dim>1) {
					val_y = x[1];
				}
				
				JNIEnv* localEnv = threading::getEnv(getJavaVM());

				jclass cls = localEnv->FindClass(
						"eu/mihosoft/vrl/types/GroovyFunction2D");

				if (localEnv->ExceptionCheck()) {
					localEnv->ExceptionDescribe();
				}

				jmethodID methodID = localEnv->GetMethodID(
						cls, "<init>", "(Ljava/lang/String;)V");

				if (localEnv->ExceptionCheck()) {
					localEnv->ExceptionDescribe();
				}

				jobject groovyFunction = localEnv->NewObject(
						cls, methodID, ug::vrl::stringC2J(localEnv, expression.c_str()));

				if (localEnv->ExceptionCheck()) {
					localEnv->ExceptionDescribe();
				}

				jmethodID runMethod = localEnv->GetMethodID(
						cls, "run",
						"(Ljava/lang/Double;Ljava/lang/Double;)Ljava/lang/Double;");

				if (localEnv->ExceptionCheck()) {
					localEnv->ExceptionDescribe();
				}

				jobject result = localEnv->CallObjectMethod(
						groovyFunction, runMethod,
						double2JObject(localEnv,val_x),
						double2JObject(localEnv,val_y));

				c = jObject2Double(localEnv,result);

				//			lua_getglobal(m_L, m_callbackName);
				//			for(size_t i = 0; i < dim; ++i)
				//				lua_pushnumber(m_L, x[i]);
				//			lua_pushnumber(m_L, time);
				//
				//			if(lua_pcall(m_L, dim + 1, 1, 0) != 0)
				//			{
				//				UG_LOG("error running diffusion callback " << m_callbackName << ": "
				//								<< lua_tostring(m_L, -1));
				//				throw(int(0));
				//			}
				//
				//			c = luaL_checknumber(m_L, -1);
				//			lua_pop(m_L, 1);
			}

		protected:
			std::string expression;
			JavaVM* javaVM;
//			jobject groovyFunction;
//			jmethodID runMethod;
		};

		template <int dim>
		void RegisterVRLUserNumber(ug::bridge::Registry& reg, const char* parentGroup) {
			std::string grp = std::string(parentGroup);

			////	Base class
			//	{
			//		stringstream ss; ss << "IUserNumberProvider" << dim << "d";
			//		reg.add_class_<IUserNumberProvider<dim> >(ss.str().c_str(), grp.c_str());
			//	}

			////	Functor
			//	{
			//	//	ConstUserNumber
			//		{
			//			typedef ConstUserNumber<dim> T;
			//			stringstream ss; ss << "ConstUserNumber" << dim << "d";
			//			reg.add_class_<T, IUserNumberProvider<dim> >(ss.str().c_str(), grp.c_str())
			//				.add_constructor()
			//				.add_method("set", &T::set)
			//				.add_method("print", &T::print);
			//		}
			//
			//	VRLUserNumber
			{
				typedef VRLUserNumber<dim> T;
				std::stringstream ss;
				ss << "VRLUserNumber" << dim << "d";
				reg.add_class_<T, IUserNumber<dim> >(ss.str().c_str(), grp.c_str())
						.add_constructor()
						.add_method("set_vrl_callback", &T::set_vrl_callback);
			}
			//	}
		}

		void RegisterVRLUserNumber(ug::bridge::Registry& reg, const char* parentGroup) {
			RegisterVRLUserNumber < 1 > (reg, parentGroup);
			RegisterVRLUserNumber < 2 > (reg, parentGroup);
			RegisterVRLUserNumber < 3 > (reg, parentGroup);
		}

	} // end vrl::
} // end ug::
