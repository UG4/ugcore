//	created by Sebastian Reiter, Andreas Vogel
//	s.b.reiter@googlemail.com
//	y10 m09 d20

#include "ug_bridge.h"
#include "registry.h"
#include "lib_algebra/algebra_selector.h"


namespace ug
{
namespace bridge
{


Registry & GetUGRegistry()
{
	static Registry ugReg;
	return ugReg;
}


#ifdef UG_ALGEBRA
	bool InitAlgebra(IAlgebraTypeSelector *algebra_type)
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
#endif

bool RegisterStandardInterfaces(Registry& reg, const char* parentGroup)
{
	bool bResult = true;
	UG_LOG("REGISTERING STANDARD INTERFACES\n");
	try
	{
		bResult &= RegisterUtilInterface(reg, parentGroup);
		bResult &= RegisterLibGridInterface(reg, parentGroup);
		bResult &= RegisterTestInterface(reg, parentGroup);
		bResult &= RegisterPCLInterface(reg, parentGroup);

		bResult &= RegisterProfileFunctions(reg, parentGroup);
		
		bResult &= RegisterMiscFunctions(reg, parentGroup);

		#ifdef UG_ALGEBRA
			bResult &= RegisterDomainInterface(reg, parentGroup);
			bResult &= RegisterUserData(reg, parentGroup);
			bResult &= RegisterStaticLibAlgebraInterface(reg, parentGroup);
			bResult &= RegisterStaticLibDiscretizationInterface(reg, parentGroup);
			// InitAlgebra
			reg.add_function("InitAlgebra", &InitAlgebra);
		#endif
	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterStandardInterfaces: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}

	return bResult;
}

bool RegisterDynamicLibDiscretizationInterface(Registry& reg, int algebra_type, const char* parentGroup)
{
	bool bResult = true;
	bResult &= RegisterDynamicLibDiscretizationInterfaceDomainIndependent(reg, algebra_type, parentGroup);
	bResult &= RegisterDynamicLibDiscretizationInterfaceDomainDependent(reg, algebra_type, parentGroup);
	bResult &= RegisterDynamicLibDiscInterfaceDiscs(reg, algebra_type, parentGroup);
	return bResult;
}

}//	end of namespace 
}//	end of namespace 
