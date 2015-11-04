
#ifndef BASICTEST_H
#define	BASICTEST_H


#include "ug.h"
#include "ugbase.h"
#include "registry/registry.h"
#include "registry/class.h"

namespace ug {
namespace vrl {

class BasicTest {
	
private:
	std::string _name;

public:

	BasicTest();

	BasicTest(std::string name);

	void setName(std::string name);

	int size() const;

	std::string get() const;

	~BasicTest();

};

SmartPtr<BasicTest> getInstanceBySmartPtr(std::string name);

void registerBasicTest(ug::bridge::Registry& reg);

} // end vrl::
}// end ug::


#endif	/* BASICTEST_H */

