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
		UG_LOG("~BasicTest" << std::endl);
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