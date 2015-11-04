
// extern headers
#include <iostream>
#include <sstream>
#include <string>

// include bridge
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_algebra_dependent.h"

// lib_disc includes
#include "lib_disc/domain.h"
#include "lib_disc/spatial_disc/manifold_assemble_util.h"


using namespace std;


namespace ug {
namespace bridge {

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
////////////////////////////////////////////////////////////////////////////////
//	Functionality
struct Functionality
{
	/**
	 * Function called for the registration of Domain dependent parts
	 * of the plugin. All Functions and Classes depending on the Domain
	 * are to be placed here when registering. The method is called for all
	 * available Domain types, based on the current build options.
	 *
	 * @param reg				registry
	 * @param parentGroup		group for sorting of functionality
	 */
	template <typename TAlgebra>
	static void Algebra(Registry& reg, string grp)
	{
		string suffix = GetAlgebraSuffix<TAlgebra>();
		string tag = GetAlgebraTag<TAlgebra>();

	//	registry of H slave exclusion functionality
		reg.add_function("MarkAllElemsForAssemblyButHSlaves",
				static_cast<void(*)(SmartPtr<IAssemble<TAlgebra> >, Grid&)>
		(&ug::MarkAllElemsForAssemblyButHSlaves<TAlgebra>), grp.c_str());
	}
};






////////////////////////////////////////////////////////////////////////////////
//	Register
void RegisterBridge_ManifoldUtil(Registry& reg, string grp)
{
	typedef Functionality Functionality;

	#if defined(UG_DIM_3)
		try{
			RegisterAlgebraDependent<Functionality>(reg,grp);
		}
		UG_REGISTRY_CATCH_THROW(grp);
	#endif
}

}//	end of namespace bridge
}//	end of namespace ug
