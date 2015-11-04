
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

