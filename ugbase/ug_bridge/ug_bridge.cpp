//	created by Sebastian Reiter, Andreas Vogel
//	s.b.reiter@googlemail.com
//	y10 m09 d20

#include "ug_bridge.h"

#ifdef UG_ALGEBRA
	#include "lib_algebra/algebra_selector.h"
#endif

#include "lib_discretization/dof_manager/dof_distribution_type.h"

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


///	Sets the default classes of class-groups based on a tags
/**	If a class has a tag (e.g. "dim=1d", "dim=2d" or "dim=3d") then it will be set
 * as default - depending on the given tags.
 */
void InitUG(int dim, const IAlgebraTypeSelector& algebraSel)
{
//	get algebra type
	AlgebraType algType = algebraSel.get_algebra_type();

//	get tag of algebra type
	std::string algTag = GetAlgebraTag(algType);

//	get dim tag
	std::string dimTag = GetDomainTag(dim);
	if(dim < 0 || dim > 3)
		throw(UGError("ERROR in InitUG: Only dimensions 1, 2, 3 are supported."));

//	get DoFDistribution tag
	//\todo: this is default, make others available
	std::string ddTag = GetDoFDistributionTag(DDT_P1CONFORMNONGROUPED);

	bridge::Registry& reg = bridge::GetUGRegistry();

//	iterate over all groups in the registry and check how many tags they contain
//	then find out if a class matches exactly this number of tags for the given
//	tag set.
	for(size_t i_grp = 0; i_grp < reg.num_class_groups(); ++i_grp)
	{
	//	get class group
		ClassGroupDesc* grp = reg.get_class_group(i_grp);

	//	count how many tags are contained in tag string
		int numTag = -1;
		for(size_t i = 0; i < grp->num_classes(); ++i)
		{
			const std::string& tag = grp->get_class_tag(i);
			int num = (int) count (tag.begin(), tag.end(), '=');
			if(numTag == -1) numTag = num;
			else if(numTag != num)
				throw(UGFatalError("Class Group with classes of different number"
									" of tags found."));
		}

	//	find the class with numTag matches
		for(size_t i = 0; i < grp->num_classes(); ++i)
		{
		//	get tag of class
			const std::string& tag = grp->get_class_tag(i);

		//	count matches
			int found = 0;
			if(tag.find(dimTag) != string::npos) ++found;
			if(tag.find(algTag) != string::npos) ++found;
			if(tag.find(ddTag) != string::npos) ++found;

		//	if exactly as many matches as tags, set this class
			if(found == numTag)
			{
				grp->set_default_class(i); break;
			}
		}
	}

	UG_LOG("INFO: InitUG successful. Setting is: ");
	UG_LOG(dimTag << ", " << algTag << ", " << ddTag << "\n");
}


bool RegisterStandardInterfaces(Registry& reg, string parentGroup)
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

		bResult &= RegisterLibDisc_Common(reg, parentGroup);
		bResult &= RegisterLibDisc_ElemDisc(reg, parentGroup);

		#ifdef UG_ALGEBRA
	//	depends on lib_algebra
		bResult &= RegisterLibAlgebra(reg, parentGroup);
		bResult &= RegisterLibDisc_Algebra(reg, parentGroup);
		bResult &= RegisterLibDisc_Domain(reg, parentGroup);
		bResult &= RegisterLibDisc_UserData(reg, parentGroup);
		#endif

		reg.add_function("InitUG", &InitUG);
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
