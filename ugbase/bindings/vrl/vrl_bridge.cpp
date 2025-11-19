/*
 * Copyright (c) 2014:  Steinbeis Forschungszentrum (STZ Ölbronn)
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

#include <string>
#include <vector>

#include "vrl_bridge.h"
#include "ug.h"
#include "ugbase.h"
#include "registry/registry.h"
#include "registry/class.h"
#include "common/util/path_provider.h"
#include "bridge/util.h"

#include "playground.h"
#include "basicTest.h"

#include "common/common.h"
#include "common/authors.h"
#include "common/util/string_util.h"

// PLEASE NOTE: UG_ALGEBRA is used here to exclude all algebra
//		and discretization related code.
//		This makes it possible to compile e.g. for target vrlgrid
#ifdef UG_ALGEBRA
#include "bridge/util_algebra_dependent.h"
#include "lib_algebra/operator/convergence_check.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#include "user_data.h"
#endif


namespace ug{
namespace vrl{

void Log(std::string s) {
	UG_LOG(s);
}

void Logln(std::string s) {
	UG_LOG(s << std::endl);
}

void ThrowIf(bool b, std::string s) {
	if (!b) {
		throw(ug::UGError(s.c_str()));
	}
}

void ThrowIfNot(bool b, std::string s) {
	if (!b) {
		throw(ug::UGError(s.c_str()));
	}
}

void registerMessaging(ug::bridge::Registry & reg) {
	reg.add_function("print", &Log, "UG4/Util/Messaging");
	reg.add_function("println", &Logln, "UG4/Util/Messaging");
}

void registerThrowUtil(ug::bridge::Registry & reg) {
	reg.add_function("throwIf", &ThrowIf, "UG4/Util");
	reg.add_function("throwIfNot", &ThrowIfNot, "UG4/Util");
}

class NumberArray {
private:
	std::vector<number> _vec;
public:

	NumberArray() {
	}

	NumberArray(std::vector<number> vec) {
		_vec = vec;
	}

	void setArray(std::vector<number> vec) {
		_vec = vec;
	}

	int size() const {
		return _vec.size();
	}

	number get(int i) const {

		if (i < 0 || (size_t) i >= _vec.size()) {
			throw UGError("NumberArray: index out of Bounds!");
		}

		return _vec[i];
	}
};

template<typename TVector>
SmartPtr<NumberArray> getDefects(const ug::StdConvCheck<TVector>* convCheck) {

	return SmartPtr<NumberArray>(new NumberArray(convCheck->get_defects()));
}

void registerNumberArray(ug::bridge::Registry & reg) {
	reg.add_class_<NumberArray>("NumberArray", "UG4/Util").add_constructor().add_method(
			"get", &NumberArray::get).add_method("size", &NumberArray::size);
}

void registerUGFinalize(ug::bridge::Registry & reg) {
	reg.add_function("UGFinalize", &ug::UGFinalize, "UG4/Util");
}

class VTest {
public:

	VTest() {
		UG_LOG("VTest::VTest() constructor used.\n")
	}

	VTest(const char* msg) {
		UG_LOG("VTest::VTest(const char*) constructor used.\n")
		UG_LOG("Message is: '" << msg << "'.\n");
	}

	std::string hello() {
		return "i am instantiated!";
	}
};

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality {

	/**
	 * Function called for the registration of Algebra dependent parts.
	 * All Functions and Classes depending on Algebra
	 * are to be placed here when registering. The method is called for all
	 * available Algebra types, based on the current build options.
	 *
	 * @param reg				registry
	 * @param parentGroup		group for sorting of functionality
	 */
	template<typename TAlgebra>
	static void Algebra(ug::bridge::Registry& reg, std::string parentGroup) {

		using vector_type = typename TAlgebra::vector_type;
		using matrix_type = typename TAlgebra::matrix_type;

//	suffix and tag
		std::string suffix = ug::bridge::GetAlgebraSuffix<TAlgebra>();
		std::string tag = ug::bridge::GetAlgebraTag<TAlgebra>();

		reg.add_function("GetDefects", &getDefects<vector_type>, "UG4/Util",
				"Defects");
	}

};
// end Functionality

void RegisterVRLFunctionality(ug::bridge::Registry& reg, std::string grp) {

// PLEASE NOTE: UG_ALGEBRA is used here to exclude all algebra
//		and discretization related code.
//		This makes it possible to compile e.g. for target vrlgrid
	#ifdef UG_ALGEBRA
		using Functionality = ug::vrl::Functionality ;
		ug::bridge::RegisterAlgebraDependent<Functionality>(reg, grp);
		ug::vrl::RegisterUserData(reg, "UG4/VRL");
	#endif
	
	ug::vrl::registerUGFinalize(reg);


	ug::vrl::registerBasicTest(reg);
	ug::vrl::registerMessaging(reg);
	ug::vrl::registerThrowUtil(reg);
	ug::vrl::registerNumberArray(reg);

	// ug::vrl::registerPlayground(reg);

	reg.add_class_<ug::vrl::VTest>("VTest", "UG4/VRL/Testing").add_constructor().add_constructor<
									void (*)(const char*)>().add_method("hello",
									&ug::vrl::VTest::hello);

}

} // end vrl::
} // end ug::

