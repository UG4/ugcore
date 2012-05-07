/*
 * eigensolver_bridge.cpp
 *
 *  Created on: 07.07.2011
 *      Author: mrupp
 */

// extern headers
#include <iostream>
#include <sstream>

#include "lib_algebra_bridge.h"

#include "lib_algebra/lib_algebra.h"

using namespace std;

#ifdef UG_PARALLEL
#include "lib_algebra/operator/eigensolver/pinvit.h"

namespace ug
{
namespace bridge
{

template <typename TAlgebra>
struct RegisterEigensolverClass
{
	static bool reg(Registry &reg, string parentGroup)
	{
	//	typedefs for this algebra
		typedef typename TAlgebra::vector_type vector_type;
		typedef typename TAlgebra::matrix_type matrix_type;

		//	get group string (use same as parent)
		string grp = string(parentGroup);
		stringstream grpSS2; grpSS2 << grp << "/EigenSolver";
		string grp2 = grpSS2.str();

	#ifdef LAPACK_AVAILABLE
			string subgroup = grp; // + string("/Preconditioner");

			string name = string("EigenSolver").append(GetAlgebraSuffix<TAlgebra>());
			typedef PINVIT<TAlgebra> T;
			reg.add_class_<T>(name, subgroup.c_str())
				.add_constructor()
				.add_method("add_vector", &T::add_vector,
							"", "vector")
				.add_method("set_preconditioner", &T::set_preconditioner,
							"", "Preconditioner")
				.add_method("set_linear_operator_A", &T::set_linear_operator_A,
							"", "LinearOperatorA")
				.add_method("set_linear_operator_B", &T::set_linear_operator_B,
							"", "LinearOperatorB")
				.add_method("set_max_iterations", &T::set_max_iterations,
								"", "precision")
				.add_method("set_precision", &T::set_precision,
								"", "precision")
				.add_method("set_pinvit", &T::set_pinvit, "", "iPINVIT", "1 = preconditioned inverse block iteration, 2 = preconditioned block gradient descent, 3 = LOBPCG")
				.add_method("apply", &T::apply);
			reg.add_class_to_group(name, "EigenSolver", GetAlgebraTag<TAlgebra>());
	#endif
		return true;
	}
};


bool RegisterEigensolver(Registry& reg, string parentGroup)
{
	return RegisterEigensolverClass<CPUAlgebra>::reg(reg, parentGroup);

}

}
}
#else
namespace ug
{
namespace bridge
{
bool RegisterEigensolver(Registry& reg, string parentGroup)
{
	return true;
}

}
}

#endif
