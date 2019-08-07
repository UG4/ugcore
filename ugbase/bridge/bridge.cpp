/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Authors: Sebastian Reiter, Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include "bridge/bridge.h"
#include "lib_algebra/algebra_type.h"
#include "common/util/path_provider.h"
#include "common/profiler/profiler.h"
#include "bridge/util.h"
#include "bridge/standard_bridges.h"

#ifdef UG_PARALLEL
#include "pcl/pcl.h"
#endif


using namespace std;

namespace ug
{
namespace bridge
{

const char* UG4_GRP = "/ug4";

///	the dimension to which ug was initialized through InitUG
/** This dimension can be accessed through GetUGDim()*/
static int UG4_DIM = -1;

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
	UG4_DIM = dim;

//	get tag of algebra type
	const std::string& algTag = GetAlgebraTag(algType);
	int blocksize = algType.blocksize();
	if( (blocksize < 0 || blocksize > 6) && blocksize != AlgebraType::VariableBlockSize)
		UG_THROW("ERROR in InitUG: Only Algebra Blocksizes '1x1', '2x2', '3x3', '4x4', '5x5', 6x6 and 'variable' are supported.");
	
	#ifdef UG_ALGEBRA
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
		#ifndef UG_CPU_5
			if(blocksize == 5)
				UG_THROW("ERROR in InitUG: Requested Algebra CPU, Blocksize '5x5' is not compiled into binary.");
		#endif
		#ifndef UG_CPU_6
			if(blocksize == 6)
				UG_THROW("ERROR in InitUG: Requested Algebra CPU, Blocksize '6x6' is not compiled into binary.");
		#endif
		#ifndef UG_CPU_VAR
			if(blocksize == AlgebraType::VariableBlockSize)
				UG_THROW("ERROR in InitUG: Requested Algebra CPU, Blocksize 'variable' is not compiled into binary.");
		#endif
			}
			else if(algType.type() == AlgebraType::GPU)
			{
				if(blocksize != 1)
					UG_THROW("ERROR in InitUG: Requested Algebra GPU, Blocksize '" << blocksize << "x" << blocksize << "' is not compiled into binary.");
			}
	#endif

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
		UG_LOG("      Parallel Environment: Num Procs="<<pcl::NumProcs()<<"\n");
	#endif
	}

	#ifdef UG_ALGEBRA
		DefaultAlgebra::set(algType);
	#endif
}

// forward for default case
void InitUG(int dim, const AlgebraType& algType)
{
	InitUG(dim, algType, false);
}


int GetUGDim()
{
	return UG4_DIM;
}

void RegisterStandardBridges(Registry& reg, string parentGroup)
{
	try
	{
		// uncomment this to register test-methods
		//RegisterBridge_Test(reg, parentGroup);

		RegisterBridge_VecMath(reg, parentGroup);
		RegisterBridge_Util(reg, parentGroup);
		RegisterBridge_PCL(reg, parentGroup);

		RegisterBridge_Profiler(reg, parentGroup);
		RegisterBridge_Misc(reg, parentGroup);
		RegisterBridge_Raster(reg, parentGroup);
		RegisterBridge_OrthoPoly(reg, parentGroup);

		#ifdef UG_GRID
			RegisterBridge_Grid(reg, parentGroup);
		#endif
		
		#ifdef UG_ALGEBRA
			RegisterBridge_Selection(reg, parentGroup);
			RegisterBridge_Domain(reg, parentGroup);
			RegisterBridge_PeriodicBoundary(reg, parentGroup);
			RegisterBridge_Refinement(reg, parentGroup);
			RegisterBridge_DomainRayTracing(reg, parentGroup);
			RegisterBridge_Transform(reg, parentGroup);
			RegisterBridge_LoadBalancing(reg, parentGroup);

		//	depends on lib_disc
			RegisterBridge_DiscCommon(reg, parentGroup);
			RegisterBridge_ElemDiscs(reg, parentGroup);

		//	depends on lib_algebra
			RegisterBridge_AlgebraCommon(reg, parentGroup);
			RegisterBridge_Preconditioner(reg, parentGroup);
			RegisterBridge_Schur(reg, parentGroup);
			RegisterBridge_Obstacle(reg, parentGroup);
			RegisterBridge_PILUT(reg, parentGroup);
			RegisterBridge_Solver(reg, parentGroup);
			RegisterBridge_Eigensolver(reg, parentGroup);
			RegisterBridge_DomainDependentPreconditioner(reg, parentGroup);
			//RegisterBridge_ConstrainedLinearIterator(reg, parentGroup);

			RegisterBridge_Restart(reg, parentGroup);

		//	depends on lib_disc
			RegisterBridge_DiscAlgebra(reg, parentGroup);
			RegisterBridge_DomainDisc(reg, parentGroup);
			RegisterBridge_GridFunction(reg, parentGroup);
			RegisterBridge_Interpolate(reg, parentGroup);
			RegisterBridge_Evaluate(reg, parentGroup);
			RegisterBridge_MaxError(reg, parentGroup);
			RegisterBridge_Ordering(reg, parentGroup);
			RegisterBridge_UserData(reg, parentGroup);
			RegisterBridge_Constraints(reg, parentGroup);
			RegisterBridge_MultiGrid(reg, parentGroup);
			RegisterBridge_Output(reg, parentGroup);
			RegisterBridge_AdaptiveTools(reg, parentGroup);
			RegisterBridge_FiniteVolume(reg, parentGroup);
			RegisterBridge_Integrate(reg, parentGroup);
			RegisterBridge_ManifoldUtil(reg, parentGroup);
			RegisterBridge_ReferenceMappingTest(reg, parentGroup);
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
#endif

#ifdef UG_ALGEBRA
		reg.add_function("InitUG", static_cast<void (*)(int, const AlgebraType&, bool)>(&InitUG), "/ug4/Init",
		                 "", string("Dimension|selection|value=[").append(availDims.str()).
		                 	 append("]#AlgebraType#verbose"));
		reg.add_function("InitUG", static_cast<void (*)(int, const AlgebraType&)>(&InitUG), "/ug4/Init",
		                 "", string("Dimension|selection|value=[").append(availDims.str()).
		                 	 append("]#AlgebraType"));
		reg.add_function("GetUGDim", &GetUGDim, "/ug4", "dimension", "", "Returns the dimension to which UG was initialized.");

	// 	AlgebraType Interface
		reg.add_class_<AlgebraType>("AlgebraType", "/ug4/Init")
			.add_constructor<void (*)(const char*, int)>("Type|selection|value=[\"CPU\"]#Blocksize|selection|value=[1,2,3,4]")
			.add_constructor<void (*)(const char*)>("Type|selection|value=[\"CPU\"]", "Variable Blocksize")
			.set_construct_as_smart_pointer(true);
#endif

	}
	UG_REGISTRY_CATCH_THROW("RegisterStandardInterfaces")
	UG_CATCH_THROW("RegisterStandardInterfaces failed.")

	reg.registry_changed();
}

}//	end of namespace 
}//	end of namespace 
