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
#include "lib_algebra/operator/operator_impl.h"

#ifdef UG_PARALLEL
#include "lib_algebra/operator/eigensolver/pinvit.h"

namespace ug
{
namespace bridge
{

template <typename TAlgebra>
struct RegisterEigensolverClass
{
	static bool reg(Registry &reg, const char *parentGroup)
	{
	//	typedefs for this algebra
		typedef CPUAlgebra algebra_type;
		typedef algebra_type::vector_type vector_type;
		typedef algebra_type::matrix_type matrix_type;

		//	get group string (use same as parent)
		std::string grp = std::string(parentGroup);
		std::stringstream grpSS2; grpSS2 << grp << "/EigenSolver";
		std::string grp2 = grpSS2.str();

	#ifdef LAPACK_AVAILABLE
			std::string subgroup = grp; // + string("/Preconditioner");

			reg.add_class_<	PINVIT<algebra_type> >("EigenSolver", subgroup.c_str())
				.add_constructor()
				.add_method("add_vector", &PINVIT<algebra_type>::add_vector,
							"", "vector")
				.add_method("set_preconditioner|interactive=false", &PINVIT<algebra_type>::set_preconditioner,
							"", "Preconditioner")
				.add_method("set_linear_operator_A|interactive=false", &PINVIT<algebra_type>::set_linear_operator_A,
							"", "LinearOperatorA")
				.add_method("set_linear_operator_B|interactive=false", &PINVIT<algebra_type>::set_linear_operator_B,
							"", "LinearOperatorB")
				.add_method("set_max_iterations|interactive=false", &PINVIT<algebra_type>::set_max_iterations,
								"", "precision")
				.add_method("set_precision|interactive=false", &PINVIT<algebra_type>::set_precision,
								"", "precision")
				.add_method("set_pinvit", &PINVIT<algebra_type>::set_pinvit, "", "iPINVIT", "1 = preconditioned inverse block iteration, 2 = preconditioned block gradient descent, 3 = LOBPCG")
				.add_method("apply", &PINVIT<algebra_type>::apply);
	#endif
		return true;
	}
};

bool RegisterEigensolver(Registry& reg, int algebra_type, const char* parentGroup)
{
	return RegisterAlgebraClass<RegisterEigensolverClass>(reg, algebra_type, parentGroup);
}
#else
bool RegisterEigensolver(Registry& reg, int algebra_type, const char* parentGroup)
{
}

}

}
}
