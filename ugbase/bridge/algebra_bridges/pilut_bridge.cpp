/*
 * pilut_bridge.cpp
 *
 */

// extern headers
#include <iostream>
#include <sstream>
#include <string>

// include bridge
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_algebra_dependent.h"

#if 0

// preconditioner
#include "lib_algebra/lib_algebra.h"
#include "lib_algebra/operator/preconditioner/pilut.h"
using namespace std;

namespace ug{
namespace bridge{
namespace Preconditioner{

/**
 * \defgroup precond_bridge Preconditioner Bridge
 * \ingroup algebra_bridge
 * \{
 */

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality2
{

/**
 * Function called for the registration of Algebra dependent parts.
 * All Functions and Classes depending on Algebra
 * are to be placed here when registering. The method is called for all
 * available Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TAlgebra>
static void Algebra(Registry& reg, string grp)
{
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();

//	typedefs for this algebra
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;



//	PILU Threshold
	{
		typedef PILUTPreconditioner<TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("PILUT").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Incomplete LU Decomposition with threshold")
			.add_constructor()
			.template add_constructor<void (*)(number)>("threshold parameter")
			.add_method("set_threshold", &T::set_threshold,
						"", "threshold", "sets threshold of incomplete LU factorisation")
			.add_method("set_info", &T::set_info,
						"", "info", "sets storage information output")
			//.add_method("set_group_size", &T::set_group_size,"", "groupSize", "sets group size")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "PILUT", tag);
	}
}
	

}; // end Functionality

// end group precond_bridge
/// \}

}// end Preconditioner

/// \addtogroup precond_bridge
void RegisterBridge_PILUT(Registry& reg, string grp)
{
	grp.append("/Algebra/Preconditioner");
	typedef Preconditioner::Functionality2 Functionality2;

	try{
		RegisterAlgebraDependent<Functionality2>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

} // namespace bridge
} // namespace ug

#endif

namespace ug{
namespace bridge{
/// \addtogroup precond_bridge
void RegisterBridge_PILUT(Registry& reg, std::string grp)
{

}

} // namespace bridge
} // namespace ug
