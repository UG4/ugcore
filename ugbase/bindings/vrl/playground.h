/*
 * Copyright (c) 2011-2012:  Steinbeis Forschungszentrum (STZ Ölbronn)
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

#ifndef PLAYGROUND_H
#define	PLAYGROUND_H

#include <string>
#include <vector>
#include <string>

#include "ug.h"
#include "registry/registry.h"
#include "registry/class.h"

namespace ug {
namespace vrl {

class TestClass {
public:
    TestClass();
    TestClass(std::string name);

    int performTest();

    std::string getRev();

    int add(int a, int b);

    std::string getString();

    SmartPtr<TestClass> smartTestImpl();

    ConstSmartPtr<TestClass> constSmartTestImpl();

    int print_name();

    int print();

    int print2();

    int print2() const;

    ~TestClass();
};

/**
 * Test Class for SmartPointer. Only created via SmartPtr.
 */
class SmartPtrCls {
public:
	SmartPtrCls();
	SmartPtrCls(std::string name);

    void print_name();

    void create_data(int size);

    ~SmartPtrCls();
private:
    std::string _name;
    byte* _data;
};

int SmartTestFunction(SmartPtr<TestClass> test);
////
int ConstSmartTestFunction(ConstSmartPtr<TestClass> test);
////template <typename TVector>
////TVector* CreateVector(int dim);
////
////template <typename TVector>
////SmartPtr<TVector> CreateVectorSmart(int dim);
//
void registerPlayground(ug::bridge::Registry& reg);
} // end vrl::
}// end ug::
//
#endif	/* PLAYGROUND_H */

