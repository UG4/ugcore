
// extern headers
#include <iostream>
#include <sstream>
#include <string>

// include bridge
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_algebra_dependent.h"

// preconditioner
#include "lib_algebra/lib_algebra.h"
#include "lib_algebra/operator/preconditioner/preconditioners.h"
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
struct Functionality
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

//	GPUJacobi
	{
		typedef Jacobi<TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("GPUJacobi").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "GPU Jacobi Preconditioner")
			.add_constructor()
			.template add_constructor<void (*)(number)>("DampingFactor")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GPUJacobi", tag);
	}
}


}; // end Functionality

// end group precond_bridge
/// \}

}// end Preconditioner

/// \addtogroup precond_bridge
void RegisterBridge_GPU(Registry& reg, string grp)
{
	grp.append("/Algebra/Preconditioner");
	typedef Preconditioner::Functionality Functionality;

	try{
		RegisterAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

} // namespace bridge
} // namespace ug
