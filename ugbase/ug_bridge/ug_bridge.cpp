//	created by Sebastian Reiter, Andreas Vogel
//	s.b.reiter@googlemail.com
//	y10 m09 d20

#include "ug_bridge.h"
#include "registry.h"

namespace ug
{
namespace bridge
{

bool RegisterStandardInterfaces(Registry& reg, const char* parentGroup)
{
	bool bResult = true;

	try
	{
		bResult &= RegisterLibGridInterface(reg, parentGroup);
		bResult &= RegisterLibAlgebraInterface(reg, parentGroup);
		bResult &= RegisterLibDiscretizationInterface(reg, parentGroup);
		bResult &= RegisterTestInterface(reg, parentGroup);
		bResult &= RegisterInfoCommands(reg, parentGroup);
	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("ERROR in RegisterStandardInterfaces: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}

	return bResult;
}

}//	end of namespace 
}//	end of namespace 
