/*
 * preconditioner_bridge.cpp
 *
 *  Created on: 03.05.2012
 *      Author: mrupp
 */

// extern headers
#include <iostream>
#include <sstream>

#include "lib_algebra_bridge.h"

#include "lib_algebra/lib_algebra.h"

// preconditioner
#include "lib_algebra/operator/preconditioner/jacobi.h"
#include "lib_algebra/operator/preconditioner/gauss_seidel.h"
#include "lib_algebra/operator/preconditioner/ilu.h"
#include "lib_algebra/operator/preconditioner/ilut.h"
#include "lib_algebra/operator/preconditioner/iterator_product.h"

using namespace std;

namespace ug
{
namespace bridge
{

template <typename TAlgebra>
struct RegisterPreconditionerClass
{
static bool reg(Registry& reg, string parentGroup)
{
//	get group string (use same as parent)
	string grp = string(parentGroup);

//	typedefs for this algebra
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;

//	suffix and tag
	string algSuffix = GetAlgebraSuffix<TAlgebra>();
	string algTag = GetAlgebraTag<TAlgebra>();
//////////////////////
// Preconditioner
//////////////////////
//	get group string
	stringstream grpSS2; grpSS2 << grp << "/Preconditioner";
	string grp2 = grpSS2.str();

//	Jacobi
	{
		typedef Jacobi<TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("Jacobi").append(algSuffix);
		reg.add_class_<T,TBase>(name, grp2, "Jacobi Preconditioner")
			.add_constructor()
			.add_method("set_damp", &T::set_damp, "", "damp")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "Jacobi", algTag);
	}

//	GaussSeidel
	{
		typedef GaussSeidel<TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("GaussSeidel").append(algSuffix);
		reg.add_class_<T,TBase>(name, grp2, "Gauss-Seidel Preconditioner")
		.add_constructor()
		.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GaussSeidel", algTag);
	}

//	Symmetric GaussSeidel
	{
		typedef SymmetricGaussSeidel<TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("SymmetricGaussSeidel").append(algSuffix);
		reg.add_class_<T,TBase>(name, grp2)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "SymmetricGaussSeidel", algTag);
	}

//	Backward GaussSeidel
	{
		typedef BackwardGaussSeidel<TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("BackwardGaussSeidel").append(algSuffix);
		reg.add_class_<T,TBase>(name, grp2)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "BackwardGaussSeidel", algTag);
	}

//	ILU
	{
		typedef ILU<TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("ILU").append(algSuffix);
		reg.add_class_<T,TBase>(name, grp2, "Incomplete LU Decomposition")
			.add_constructor()
			.add_method("set_beta", &T::set_beta, "", "beta")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ILU", algTag);
	}

//	ILU Threshold
	{
		typedef ILUTPreconditioner<TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("ILUT").append(algSuffix);
		reg.add_class_<T,TBase>(name, grp2, "Incomplete LU Decomposition with threshold")
			.add_constructor()
			.add_method("set_threshold", &T::set_threshold,
						"", "threshold", "sets threshold of incomplete LU factorisation")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ILUT", algTag);
	}

//	LinearIteratorProduct
	{
		typedef LinearIteratorProduct<vector_type, vector_type> T;
		typedef ILinearIterator<vector_type> TBase;
		string name = string("LinearIteratorProduct").append(algSuffix);
		reg.add_class_<T,TBase>(name, grp,
						"Linear Iterator consisting of several LinearIterations")
				.add_constructor()
				.add_method("add_iteration", &T::add_iterator, "Add an iterator")
				.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "LinearIteratorProduct", algTag);
	}
	return true;
}
};

bool RegisterPreconditioner(Registry& reg, string parentGroup)
{
	return RegisterAlgebraClass<RegisterPreconditionerClass>(reg, parentGroup);
}

}
}
