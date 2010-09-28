//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m09 d20

#include "../ugbridge/registry.h"

using namespace std;

namespace ug{
namespace interface
{

int Add(int a, int b)
{
	return a + b;
}

class Test 
{
	public:
		Test()
		{
			UG_LOG("Test created!n");
		}
		
		int add(int a, int b)
		{
			return a+b;
		}
	
		int print_name()
		{
			UG_LOG("Name is Test\n");
			return 1;
		}
};

void RegisterLibGridInterface(InterfaceRegistry& reg)
{
	reg.add_function("add", &Add, "c", "a,b");
	
	reg.add_class_<Test>("Test")
		.add_method("add", &Test::add, "c", "a,b");
}

}//	end of namespace 
}//	end of namespace 
