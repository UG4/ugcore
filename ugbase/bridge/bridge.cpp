//	created by Sebastian Reiter, Andreas Vogel
//	s.b.reiter@googlemail.com
//	y10 m09 d20

#include "bridge/bridge.h"
#include "lib_algebra/algebra_type.h"
#include "common/util/path_provider.h"
#include "common/profiler/profiler.h"
#include "bridge/util.h"

#ifdef UG_PARALLEL
#include "pcl/pcl.h"
#endif


using namespace std;

namespace ug
{
namespace bridge
{

const char* UG4_GRP = "/ug4";

Registry & GetUGRegistry()
{
	static Registry ugReg;
	return ugReg;
}

/// calls RegisterStandardInterfaces
void InitBridge()
{
	PROFILE_FUNC();
	//	initialize ug-interfaces
	RegisterStandardBridges(bridge::GetUGRegistry());
}


///	Sets the default classes of class-groups based on a tags
/**	If a class has a tag (e.g. "dim=1d", "dim=2d" or "dim=3d") then it will be set
 * as default - depending on the given tags.
 */
void InitUG(int dim, const AlgebraType& algType, bool verbose)
{
	PROFILE_FUNC();
//	get tag of algebra type
	const std::string& algTag = GetAlgebraTag(algType);
	int blocksize = algType.blocksize();
	if( (blocksize < 0 || blocksize > 4) && blocksize != AlgebraType::VariableBlockSize)
		UG_THROW("ERROR in InitUG: Only Algebra Blocksizes '1x1', '2x2', '3x3', '4x4' and 'variable' are supported.");
	
	if(algType.type() == AlgebraType::CPU)
	{
#ifndef UG_CPU_1
	if(blocksize == 1)
		UG_THROW("ERROR in InitUG: Requested Algebra CPU, Blocksize '1x1' is not compiled into binary.");
#endif
#ifndef UG_CPU_2
	if(blocksize == 2)
		UG_THROW("ERROR in InitUG: Requested Algebra CPU, Blocksize '2x2' is not compiled into binary.");
#endif
#ifndef UG_CPU_3
	if(blocksize == 3)
		UG_THROW("ERROR in InitUG: Requested Algebra CPU, Blocksize '3x3' is not compiled into binary.");
#endif
#ifndef UG_CPU_4
	if(blocksize == 4)
		UG_THROW("ERROR in InitUG: Requested Algebra CPU, Blocksize '4x4' is not compiled into binary.");
#endif
#ifndef UG_CPU_VAR
	if(blocksize == AlgebraType::VariableBlockSize)
		UG_THROW("ERROR in InitUG: Requested Algebra CPU, Blocksize 'variable' is not compiled into binary.");
#endif
	}
	else if(algType.type() == AlgebraType::CRS)
	{
#ifndef UG_CRS_1
	if(blocksize == 1)
		UG_THROW("ERROR in InitUG: Requested Algebra CRS, Blocksize '1x1' is not compiled into binary.");
#endif
#ifndef UG_CRS_2
	if(blocksize == 2)
		UG_THROW("ERROR in InitUG: Requested Algebra CRS, Blocksize '2x2' is not compiled into binary.");
#endif
#ifndef UG_CRS_3
	if(blocksize == 3)
		UG_THROW("ERROR in InitUG: Requested Algebra CRS, Blocksize '3x3' is not compiled into binary.");
#endif
#ifndef UG_CRS_4
	if(blocksize == 4)
		UG_THROW("ERROR in InitUG: Requested Algebra CRS, Blocksize '4x4' is not compiled into binary.");
#endif
#ifndef UG_CRS_VAR
	if(blocksize == AlgebraType::VariableBlockSize)
		UG_THROW("ERROR in InitUG: Requested Algebra CRS, Blocksize 'variable' is not compiled into binary.");
#endif
	}

//	get dim tag
	std::string dimTag = GetDimensionTag(dim);
	if(dim < 0 || dim > 3)
		UG_THROW("ERROR in InitUG: Only dimensions 1, 2, 3 are supported.");
#ifndef UG_DIM_1
	if(dim == 1)
		UG_THROW("ERROR in InitUG: Requested Dimension '1d' is not compiled into binary.");
#endif
#ifndef UG_DIM_2
	if(dim == 2)
		UG_THROW("ERROR in InitUG: Requested Dimension '2d' is not compiled into binary.");
#endif
#ifndef UG_DIM_3
	if(dim == 3)
		UG_THROW("ERROR in InitUG: Requested Dimension '3d' is not compiled into binary.");
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
				UG_THROW("Class Group with classes of different number"
									" of tags found.");
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

		//	if exactly as many matches as tags, set this class
			if(found == numTag)
			{
				grp->set_default_class(i); break;
			}
		}
	}

	if(verbose){
		UG_LOG("INFO: InitUG successful. Setting is: ");
		UG_LOG(dimTag << " " << algTag << "\n");
	#ifdef UG_PARALLEL
		UG_LOG("      Parallel Environment: Num Procs="<<pcl::GetNumProcesses()<<"\n");
	#endif
	}

	DefaultAlgebra::set(algType);
}

// forward for default case
void InitUG(int dim, const AlgebraType& algType)
{
	InitUG(dim, algType, false);
}


void RegisterStandardBridges(Registry& reg, string parentGroup)
{
	try
	{
		// uncomment this to register test-methods
		//RegisterBridge_Test(reg, parentGroup);

		RegisterBridge_VecMath(reg, parentGroup);
		RegisterBridge_Util(reg, parentGroup);
		RegisterBridge_Grid(reg, parentGroup);
		RegisterBridge_PCL(reg, parentGroup);
		RegisterBridge_Domain(reg, parentGroup);
		RegisterBridge_PeriodicBoundary(reg, parentGroup);
		RegisterBridge_Refinement(reg, parentGroup);
		RegisterBridge_Selection(reg, parentGroup);
		RegisterBridge_Transform(reg, parentGroup);
		RegisterBridge_LoadBalancing(reg, parentGroup);

		RegisterBridge_Profiler(reg, parentGroup);
		RegisterBridge_Misc(reg, parentGroup);

		#ifdef UG_ALGEBRA
	//	depends on lib_disc
		RegisterBridge_DiscCommon(reg, parentGroup);
		RegisterBridge_ElemDiscs(reg, parentGroup);

	//	depends on lib_algebra
		RegisterBridge_AlgebraCommon(reg, parentGroup);
		RegisterBridge_Preconditioner(reg, parentGroup);
		RegisterBridge_Solver(reg, parentGroup);
		RegisterBridge_Eigensolver(reg, parentGroup);
		RegisterBridge_DomainDependentPreconditioner(reg, parentGroup);

	//	depends on lib_disc
		RegisterBridge_DiscAlgebra(reg, parentGroup);
		RegisterBridge_DomainDisc(reg, parentGroup);
		RegisterBridge_GridFunction(reg, parentGroup);
		RegisterBridge_Interpolate(reg, parentGroup);
		RegisterBridge_MaxError(reg, parentGroup);
		RegisterBridge_Ordering(reg, parentGroup);
		RegisterBridge_UserData(reg, parentGroup);
		RegisterBridge_Constraints(reg, parentGroup);
		RegisterBridge_MultiGrid(reg, parentGroup);
		RegisterBridge_Output(reg, parentGroup);
		RegisterBridge_AdaptiveTools(reg, parentGroup);
		RegisterBridge_FiniteVolume(reg, parentGroup);
		RegisterBridge_Integrate(reg, parentGroup);
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

#ifdef UG_ALGEBRA
		reg.add_function("InitUG", static_cast<void (*)(int, const AlgebraType&, bool)>(&InitUG), "/ug4/Init",
		                 "", string("Dimension|selection|value=[").append(availDims.str()).
		                 	 append("]#AlgebraType#verbose"));
		reg.add_function("InitUG", static_cast<void (*)(int, const AlgebraType&)>(&InitUG), "/ug4/Init",
		                 "", string("Dimension|selection|value=[").append(availDims.str()).
		                 	 append("]#AlgebraType"));

	// 	AlgebraType Interface
		reg.add_class_<AlgebraType>("AlgebraType", "/ug4/Init")
			.add_constructor<void (*)(const char*, int)>("Type|selection|value=[\"CPU\"]#Blocksize|selection|value=[1,2,3,4]")
			.add_constructor<void (*)(const char*)>("Type  (Blocksize=variable)|selection|value=[\"CPU\"]")
			.set_construct_as_smart_pointer(true);
#endif

	}
	UG_REGISTRY_CATCH_THROW("RegisterStandardInterfaces")
	UG_CATCH_THROW("RegisterStandardInterfaces failed.")

	reg.registry_changed();
}

}//	end of namespace 
}//	end of namespace 
