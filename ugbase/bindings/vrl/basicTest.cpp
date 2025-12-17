/*
 * Copyright (c) 2012:  G-CSC, Goethe University Frankfurt
 * Author: Christian Poliwoda
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

#include "basicTest.h"
#include "bindings_vrl.h"

#include <string>

namespace ug {
namespace vrl {

	BasicTest::BasicTest() { 
		_name="basic";	
	}
	

	BasicTest::BasicTest(std::string name) {
		_name = name;
	}

	void BasicTest::setName(std::string name) {
		_name = name;
	}

	int BasicTest::size() const {
		return _name.size();
	}

	std::string BasicTest::get() const {
		return _name;
	}
	
	BasicTest::~BasicTest() {
		UG_LOG("~BasicTest\n");
	}
	

SmartPtr<BasicTest> getInstanceBySmartPtr(std::string name) {

	return SmartPtr<BasicTest > ( new BasicTest(name) );
}

void registerBasicTest(ug::bridge::Registry& reg) {
	reg.add_class_<BasicTest > ("BasicTest", "UG4/Test")
			.add_constructor()
			.add_constructor<void (*)(std::string) >()
			.add_method("get", &BasicTest::get)
			.add_method("size", &BasicTest::size);
	reg.add_function("getInstanceBySmartPtr", &getInstanceBySmartPtr, "UG4/Test", "getInstanceBySmartPtr");
}

} // end vrl::
}// end ug::