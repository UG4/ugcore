//	created by Sebastian Reiter, Andreas Vogel
//	s.b.reiter@googlemail.com
//	y10 m09 d20

#include "ug_bridge.h"
#include "registry.h"

namespace ug
{
namespace bridge
{

void RegisterStandardInterfaces(Registry& reg, const char* parentGroup)
{
	RegisterLibGridInterface(reg, parentGroup);
	RegisterLibAlgebraInterface(reg, parentGroup);
	RegisterLibDiscretizationInterface(reg, parentGroup);
	RegisterTestInterface(reg, parentGroup);
	RegisterInfoCommands(reg, parentGroup);
}

}//	end of namespace 
}//	end of namespace 
