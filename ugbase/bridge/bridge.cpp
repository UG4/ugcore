//	created by Sebastian Reiter, Andreas Vogel
//	s.b.reiter@googlemail.com
//	y10 m09 d20

#include "bridge.h"
#include "lib_algebra/algebra_type.h"
#include "lib_disc/dof_manager/dof_distribution_type.h"

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
void InitUG(int dim, const AlgebraType& algType, DofDistributionType ddType)
{
//	get tag of algebra type
	std::string algTag = GetAlgebraTag(algType);

//	get dim tag
	std::string dimTag = GetDomainTag(dim);
	if(dim < 0 || dim > 3)
		UG_THROW_FATAL("ERROR in InitUG: Only dimensions 1, 2, 3 are supported.");
#ifndef UG_DIM_1
	if(dim == 1)
		UG_THROW_FATAL("ERROR in InitUG: Requested Dimension '1d' is not compiled into binary.");
#endif
#ifndef UG_DIM_2
	if(dim == 2)
		UG_THROW_FATAL("ERROR in InitUG: Requested Dimension '2d' is not compiled into binary.");
#endif
#ifndef UG_DIM_3
	if(dim == 3)
		UG_THROW_FATAL("ERROR in InitUG: Requested Dimension '3d' is not compiled into binary.");
#endif

//	get DoFDistribution tag
	std::string ddTag = GetDoFDistributionTag(ddType);
#ifndef DOF_P1
	if(ddType == DDT_P1CONFORM)
		UG_THROW_FATAL("ERROR in InitUG: Requested DoFManager 'P1' is not compiled into binary.");
#endif
#ifndef DOF_GEN
	if(ddType == DDT_CONFORM)
		UG_THROW_FATAL("ERROR in InitUG: Requested DoFManager 'GEN' is not compiled into binary.");
#endif

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
			int num = (int) count (tag.begin(), tag.end(), ';');
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
	UG_LOG(dimTag << " " << algTag << " " << ddTag << "\n");
#ifdef UG_PARALLEL
	UG_LOG("      Parallel Environment: Num Procs="<<pcl::GetNumProcesses()<<"\n");
#endif
}

///	Sets the default classes of class-groups based on a tags using default DoFManager
void InitUG(int dim, const AlgebraType& algType)
{
//	use default dof distribution
	InitUG(dim, algType, DDT_P1CONFORM);
}

///	Sets the default classes of class-groups based on a tags
void InitUG(int dim, const AlgebraType& algType, const char* ddType)
{
	std::string dd = ddType;
	if(dd == "P1") {InitUG(dim, algType, DDT_P1CONFORM); return;}
	if(dd == "GEN")  {InitUG(dim, algType, DDT_CONFORM); return;}
	UG_LOG("ERROR in 'InitUG': For DofManager choose one within ['P1', 'GEN'].\n");
	throw(UGFatalError("DoFManager not recognized."));
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

		#ifdef UG_ALGEBRA
	//	depends on lib_disc
		bResult &= RegisterLibDisc_Common(reg, parentGroup);
		bResult &= RegisterLibDisc_ElemDisc(reg, parentGroup);

	//	depends on lib_algebra
		bResult &= RegisterLibAlgebra(reg, parentGroup);
		bResult &= RegisterLibDisc_Algebra(reg, parentGroup);
		bResult &= RegisterLibDisc_Domain(reg, parentGroup);
		bResult &= RegisterLibDisc_UserData(reg, parentGroup);
		#endif


	//	build a string with all compiled dimensions
		stringstream availDims; bool first = true;
#ifdef UG_DIM_1
		if(!first) {availDims << ",";}; availDims << "1";
		first = false;
#endif
#ifdef UG_DIM_2
		if(!first) {availDims << ",";}; availDims << "2";
		first = false;
#endif
#ifdef UG_DIM_3
		if(!first) {availDims << ",";}; availDims << "3";
		first = false;
#endif

	//	build a string with all compiled DoFManager
		stringstream availDofManager; first = true;
#ifdef DOF_P1
		if(!first) {availDofManager << ",";}; availDofManager << "\"P1\"";
		first=false;
#endif
#ifdef DOF_GEN
		if(!first) {availDofManager << ",";}; availDofManager << "\"GEN\"";
		first=false;
#endif

#ifdef UG_ALGEBRA
		reg.add_function("InitUG", static_cast<void (*)(int, const AlgebraType&, const char *)>(&InitUG), "/ug4",
		                 "", string("Dimension|selection|value=[").append(availDims.str()).
		                 	 append("]#Algebra Type#DoFManager|selection|value=[").
		                 	 append(availDofManager.str()).append("]"));
		reg.add_function("InitUG", static_cast<void (*)(int, const AlgebraType&)>(&InitUG), "/ug4",
		                 "", string("Dimension|selection|value=[").append(availDims.str()).
		                 	 append("]#Algebra Type"));
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
