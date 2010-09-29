//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m09 d20

#include "ug_interface.h"
#include "registry.h"

namespace ug{
namespace interface
{




void RegisterStandardInterfaces(InterfaceRegistry& reg)
{
	RegisterLibGridInterface(reg);
	RegisterLibDiscretizationInterface(reg);
	RegisterTestInterface(reg);
}

}//	end of namespace 
}//	end of namespace 
