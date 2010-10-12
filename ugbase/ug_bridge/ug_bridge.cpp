//	created by Sebastian Reiter, Andreas Vogel
//	s.b.reiter@googlemail.com
//	y10 m09 d20

#include "ug_bridge.h"
#include "registry.h"

namespace ug
{
namespace bridge
{

void RegisterStandardInterfaces(Registry& reg)
{
	RegisterLibGridInterface(reg);
	RegisterLibAlgebraInterface(reg);
	RegisterLibDiscretizationInterface(reg);
	RegisterTestInterface(reg);
}

}//	end of namespace 
}//	end of namespace 
