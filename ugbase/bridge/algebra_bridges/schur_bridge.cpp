/*
 * schur_bridge.cpp
 *
 *  Created on: 08.01.2014
 *      Author: mrupp
 */

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
#include "lib_algebra/operator/preconditioner/schur/schur.h"

#include "../util_overloaded.h"
using namespace std;

namespace ug{

namespace bridge{
namespace Schur{

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

#ifdef UG_PARALLEL
	// 	Schur complement preconditioner
	{
		typedef SchurPrecond<TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("SchurComplement").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Schur complement preconditioner")
			.add_constructor()
			.add_method("set_dirichlet_solver", &T::set_dirichlet_solver, "","Dirichlet solver")
			.add_method("set_skeleton_solver", &T::set_skeleton_solver, "","Skeleton solver")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "SchurComplement", tag);
	}
#endif
}


}; // end Functionality

// end group precond_bridge
/// \}

}// end Preconditioner

/// \addtogroup precond_bridge
void RegisterBridge_Schur(Registry& reg, string grp)
{
	grp.append("/Algebra/Preconditioner");
	typedef Schur::Functionality Functionality;

	try{
		RegisterAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

} // namespace bridge
} // namespace ug
