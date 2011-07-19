//	created by Sebastian Reiter, Andreas Vogel
//	s.b.reiter@googlemail.com
//	y10 m09 d20

#include "ug_bridge.h"

#ifdef UG_ALGEBRA
	#include "lib_algebra/algebra_selector.h"
#endif

using namespace std;

namespace ug
{
namespace bridge
{


Registry & GetUGRegistry()
{
	static Registry ugReg;
	return ugReg;
}

///	Sets the default classes of class-groups based on a dimension tag
/**	If a class has a tag "dim=1", "dim=2" or "dim=3" then it will be set
 * as default - depending on the given dim.
 */
static void SetDefaultDimension(int dim)
{
	if(dim < 0 || dim > 3){
		throw(UGError("ERROR in SetDefaultDimension: Only dimensions 1, 2, 3 are supported."));
	}

//	iterate over all groups in the registry and check whether they contain
//	a classTag "dim=nd", where n is the given dimension.
	char dimTag[16];
	sprintf(dimTag, "dim=%dd", dim);

	bridge::Registry& reg = bridge::GetUGRegistry();

	for(size_t i_grp = 0; i_grp < reg.num_class_groups(); ++i_grp){
		ClassGroupDesc* grp = reg.get_class_group(i_grp);
		for(size_t i = 0; i < grp->num_classes(); ++i){
			const char* tag = grp->get_class_tag(i);
			if(strstr(tag, dimTag) != NULL){
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

	bResult &= reg.registry_changed();

	return bResult;
}

#endif

bool RegisterStandardInterfaces(Registry& reg, const char* parentGroup)
{
	bool bResult = true;
	try
	{
		bResult &= RegisterVecMathBridge(reg, parentGroup);
		bResult &= RegisterUtilInterface(reg, parentGroup);
		bResult &= RegisterLibGridInterface(reg, parentGroup);
		bResult &= RegisterTestInterface(reg, parentGroup);
		bResult &= RegisterPCLInterface(reg, parentGroup);
		bResult &= RegisterDomainInterface(reg, parentGroup);
		bResult &= RegisterRefinementBridge(reg, parentGroup);

		bResult &= RegisterProfileFunctions(reg, parentGroup);
		bResult &= RegisterMiscFunctions(reg, parentGroup);

		reg.add_function("SetDefaultDimension", &SetDefaultDimension);


		#ifdef UG_ALGEBRA
		//	depends on lib_algebra
			bResult &= RegisterStaticLibDiscInterface(reg, parentGroup);

		//	does not depend on lib_algebra
			bResult &= RegisterLibDiscElemDisc(reg, parentGroup);

		//	depends on lib_algebra
			bResult &= RegisterUserData(reg, parentGroup);
			bResult &= RegisterStaticLibAlgebraInterface(reg, parentGroup);

			// InitAlgebra
			reg.add_function("InitAlgebra", &InitAlgebra);
		#endif
	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed& ex)
	{
		UG_LOG("### ERROR in RegisterStandardInterfaces: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}

	bResult &= reg.registry_changed();

	return bResult;
}


}//	end of namespace 
}//	end of namespace 
