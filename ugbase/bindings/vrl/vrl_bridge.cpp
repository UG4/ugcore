#include <string>
#include <vector>

#include "vrl_bridge.h"
#include "ug.h"
#include "ugbase.h"
#include "registry/registry.h"
#include "registry/class.h"
#include "common/util/path_provider.h"
#include "bridge/util.h"
#include "bridge/util_algebra_dependent.h"
#include "lib_algebra/operator/convergence_check.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#include "user_data.h"

#include "common/common.h"
#include "common/authors.h"
#include "common/util/string_util.h"

namespace ug {
namespace vrl {

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality {

	/**
	 * Function called for the registration of Algebra dependent parts.
	 * All Functions and Classes depending on Algebra
	 * are to be placed here when registering. The method is called for all
	 * available Algebra types, based on the current build options.
	 *
	 * @param reg				registry
	 * @param parentGroup		group for sorting of functionality
	 */
	template<typename TAlgebra>
	static void Algebra(ug::bridge::Registry& reg, std::string parentGroup) {
//	typedefs for Vector and Matrix
		typedef typename TAlgebra::vector_type vector_type;
		typedef typename TAlgebra::matrix_type matrix_type;

//	suffix and tag
		std::string suffix = ug::bridge::GetAlgebraSuffix<TAlgebra>();
		std::string tag = ug::bridge::GetAlgebraTag<TAlgebra>();

		reg.add_function("GetDefects", &getDefects<vector_type>, "UG4/Util",
				"Defects");
	}

};
// end Functionality

void RegisterVRLFunctionality(ug::bridge::Registry& reg, std::string grp) {
	typedef ug::vrl::Functionality Functionality;

	ug::bridge::RegisterAlgebraDependent<Functionality>(reg, grp);
}

} // end vrl::
} // end ug::

