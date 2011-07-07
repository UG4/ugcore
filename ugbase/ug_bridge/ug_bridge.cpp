//	created by Sebastian Reiter, Andreas Vogel
//	s.b.reiter@googlemail.com
//	y10 m09 d20

#include "ug_bridge.h"
#include "registry.h"

#ifdef UG_ALGEBRA
	#include "lib_algebra/algebra_selector.h"
#endif


namespace ug
{
namespace bridge
{


Registry & GetUGRegistry()
{
	static Registry ugReg;
	return ugReg;
}

///	Sets the default classes of class-groups based on a dimension suffix
/**	Only classes which end with 1d, 2d, ... are considered.
 * Note that this method only has effect onto methods which are currently
 * registered at the registry.
 */
static void SetDefaultDimension(int dim)
{
//	iterate over all groups in the registry and check whether they contain
//	a class which ends with the dim-suffix
	char suffix[8];
	sprintf(suffix, "%dd", dim);

	bridge::Registry& reg = bridge::GetUGRegistry();

	for(size_t i_grp = 0; i_grp < reg.num_class_groups(); ++i_grp){
		ClassGroupDesc* grp = reg.get_class_group(i_grp);
		for(size_t i = 0; i < grp->num_classes(); ++i){
			const IExportedClass* c = grp->get_class(i);
			size_t len = strlen(c->name());
			const char* ending = c->name() + len - 2;
			if(strcmp(ending, suffix) == 0){
				grp->set_default_class(i);
				break;
			}
		}
	}
}

#ifdef UG_ALGEBRA

bool RegisterDynamicLibDiscretizationInterface(Registry& reg, int algebra_type, const char* parentGroup)
{
	bool bResult = true;
	bResult &= RegisterDynamicLibDiscAlgebra(reg, algebra_type, parentGroup);
	bResult &= RegisterDynamicLibDiscDomain(reg, algebra_type, parentGroup);
	return bResult;
}

bool InitAlgebra(IAlgebraTypeSelector *algebra_type)
{
	bridge::Registry& reg = bridge::GetUGRegistry();
	bool bResult = true;
	try
	{
		bResult &= RegisterDynamicLibAlgebraInterface(reg, algebra_type->get_algebra_type(), "/ug4");
		bResult &= RegisterDynamicLibDiscretizationInterface(reg, algebra_type->get_algebra_type(), "/ug4");
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
	try
	{
		bResult &= RegisterUtilInterface(reg, parentGroup);
		bResult &= RegisterLibGridInterface(reg, parentGroup);
		bResult &= RegisterTestInterface(reg, parentGroup);
		bResult &= RegisterPCLInterface(reg, parentGroup);

		bResult &= RegisterProfileFunctions(reg, parentGroup);
		
		bResult &= RegisterMiscFunctions(reg, parentGroup);

		bResult &= RegisterDomainInterface(reg, parentGroup);

		bResult &= RegisterLibDiscElemDisc(reg, parentGroup);

		reg.add_function("SetDefaultDimension", &SetDefaultDimension);

		#ifdef UG_ALGEBRA
			bResult &= RegisterUserData(reg, parentGroup);
			bResult &= RegisterStaticLibAlgebraInterface(reg, parentGroup);
			bResult &= RegisterStaticLibDiscInterface(reg, parentGroup);
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

}//	end of namespace 
}//	end of namespace 
