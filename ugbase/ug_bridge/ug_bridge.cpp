//	created by Sebastian Reiter, Andreas Vogel
//	s.b.reiter@googlemail.com
//	y10 m09 d20

#include "ug_bridge.h"
#include "registry.h"
#include "lib_algebra/algebra_chooser.h"

namespace ug
{
namespace bridge
{

Registry & GetUGRegistry()
{
	static Registry ugReg;
	return ugReg;
}



bool InitAlgebra(AlgebraTypeChooserInterface *algebra_type)
{
	bridge::Registry& reg = bridge::GetUGRegistry();
	bool bResult = true;
	try
	{
		bResult &= RegisterDynamicLibAlgebraInterface(reg, algebra_type->get_algebra_type());
		bResult &= RegisterDynamicLibDiscretizationInterface(reg, algebra_type->get_algebra_type());
	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterStandardInterfaces: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}
	reg.registry_changed();
	return true;
}


bool RegisterStandardInterfaces(Registry& reg, const char* parentGroup)
{
	bool bResult = true;

	try
	{
		bResult &= RegisterUserData(reg, parentGroup);
		bResult &= RegisterLibGridInterface(reg, parentGroup);
		bResult &= RegisterStaticLibAlgebraInterface(reg, parentGroup);
		bResult &= RegisterStaticLibDiscretizationInterface(reg, parentGroup);
		bResult &= RegisterTestInterface(reg, parentGroup);

		// InitAlgebra
		reg.add_function("InitAlgebra", &InitAlgebra);
	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterStandardInterfaces: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}

	return bResult;
}

}//	end of namespace 
}//	end of namespace 
