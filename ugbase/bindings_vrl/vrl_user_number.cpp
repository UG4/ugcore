

#include "../ug_bridge/ug_bridge.h"
#include "type_converter.h"
#include "common/common.h"
#include "../lib_discretization/spatial_discretization/ip_data/const_user_data.h"
#include "../ug_script/ug_script.h"
#include "bindings_vrl.h"
#include "type_converter.h"
#include "threading.h"
#include <iostream>
#include <sstream>

namespace ug {
namespace vrl {

template <int dim>
class VRLUserNumber : public IUserData<number, dim> {
public:

	/// Base class type
	typedef IUserData<number, dim> base_type;

	using base_type::num_series;
	using base_type::num_ip;
	using base_type::ip;
	using base_type::time;
	using base_type::value;


	//	Functor Type
	typedef typename IUserData<number, dim>::functor_type functor_type;

	//	return functor

	virtual functor_type get_functor() const {
		return boost::ref(*this);
	}

public:

	VRLUserNumber() {
		javaVM = getJavaVM();
		expression = "";
	}

	void set_vrl_callback(const char* expression) {
		this->expression = expression;
	}

	void operator() (number& c, const MathVector<dim>& x, number time = 0.0) const {

		double val_x = x[0];
		double val_y = 0;

		if (dim > 1) {
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
				double2JObject(localEnv, val_x),
				double2JObject(localEnv, val_y));

		c = jObject2Double(localEnv, result);
	}

	///	implement as a IPData

	virtual void compute(bool computeDeriv = false) {
		for (size_t s = 0; s < num_series(); ++s)
			for (size_t i = 0; i < num_ip(s); ++i) {
				this->operator()(value(s, i),
						ip(s, i),
						time());
			}
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
		reg.add_class_<T, IUserData<number, dim> >(ss.str().c_str(), grp.c_str())
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
