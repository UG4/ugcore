/*
 * preconditioner_bridge.cpp
 *
 *  Created on: 03.05.2012
 *      Author: avogel
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
#include "lib_algebra/operator/preconditioner/preconditioners.h"
#include "lib_algebra/operator/preconditioner/ilut_scalar.h"
#include "lib_algebra/operator/linear_solver/agglomerating_solver.h"
#include "lib_algebra/operator/preconditioner/block_gauss_seidel.h"

#include "../util_overloaded.h"
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


//	Jacobi
	{
		typedef Jacobi<TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("Jacobi").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Jacobi Preconditioner")
			.add_constructor()
			.template add_constructor<void (*)(number)>("DampingFactor")
			//.add_method("set_block", &T::set_block, "", "block", "if true, use block smoothing (default), else diagonal smoothing")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "Jacobi", tag);
	}

//	GaussSeidelBase
	{
		typedef GaussSeidelBase<TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("GaussSeidelBase").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Gauss-Seidel Base")
			.add_method("set_sor_relax", &T::set_sor_relax,
					"", "sor relaxation", "sets sor relaxation parameter");
		reg.add_class_to_group(name, "GaussSeidelBase", tag);
	}

//	GaussSeidel
	{
		typedef GaussSeidel<TAlgebra> T;
		typedef GaussSeidelBase<TAlgebra> TBase;
		string name = string("GaussSeidel").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Gauss-Seidel Preconditioner")
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GaussSeidel", tag);
	}

//	Symmetric GaussSeidel
	{
		typedef SymmetricGaussSeidel<TAlgebra> T;
		typedef GaussSeidelBase<TAlgebra> TBase;
		string name = string("SymmetricGaussSeidel").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Symmetric Gauss Seidel Preconditioner")
				.add_constructor()
				.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "SymmetricGaussSeidel", tag);
	}

//	Backward GaussSeidel
	{
		typedef BackwardGaussSeidel<TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("BackwardGaussSeidel").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Backward Gauss Seidel Preconditioner")
				.add_constructor()
				.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "BackwardGaussSeidel", tag);
	}

	//	BlockGaussSeidel2
	{
		typedef BlockGaussSeidel2<TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("BlockGaussSeidel2").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Block Gauss Preconditioner2")
			.add_constructor()
			. ADD_CONSTRUCTOR( (int) )("depth")
			.add_method("set_depth", &T::set_depth)
			.set_construct_as_smart_pointer(true);

		reg.add_class_to_group(name, "BlockGaussSeidel2", tag);
	}

	//	BlockGaussSeidel
	{
		typedef BlockGaussSeidel<TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("BlockGaussSeidel").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Block Gauss Preconditioner")
			.add_constructor()
			. ADD_CONSTRUCTOR( (int) )("depth")
			.add_method("set_depth", &T::set_depth)
			.set_construct_as_smart_pointer(true);

		reg.add_class_to_group(name, "BlockGaussSeidel", tag);
	}

//	ILU
	{
		typedef ILU<TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("ILU").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Incomplete LU Decomposition")
			.add_constructor()
			.add_method("set_beta", &T::set_beta, "", "beta")
			.add_method("set_sort", &T::set_sort, "", "bSort", "if bSort=true, use a cuthill-mckey sorting to reduce fill-in. default false")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ILU", tag);
	}

//	ILU Threshold
	{
		typedef ILUTPreconditioner<TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("ILUT").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Incomplete LU Decomposition with threshold")
			.add_constructor()
			.template add_constructor<void (*)(number)>("threshold parameter")
			.add_method("set_threshold", &T::set_threshold,
						"", "threshold", "sets threshold of incomplete LU factorisation")
			.add_method("set_info", &T::set_info,
						"", "info", "sets storage information output")
			.add_method("set_sort", &T::set_sort, "", "bSort", "if bSort=true, use a cuthill-mckey sorting to reduce fill-in. default true")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ILUT", tag);
	}
//	ILU Threshold Scalar
	{
		typedef ILUTScalarPreconditioner<TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("ILUTScalar").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Scalar Incomplete LU Decomposition with threshold")
			.add_constructor()
			.template add_constructor<void (*)(number)>("threshold parameter")
			.add_method("set_threshold", &T::set_threshold,
						"", "threshold", "sets threshold of incomplete LU factorisation")
			.add_method("set_info", &T::set_info,
						"", "info", "sets storage information output")
			.add_method("set_sort", &T::set_sort, "", "bSort", "if bSort=true, use a cuthill-mckey sorting to reduce fill-in")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ILUTScalar", tag);
	}

//	LinearIteratorProduct
	{
		typedef LinearIteratorProduct<vector_type, vector_type> T;
		typedef ILinearIterator<vector_type> TBase;
		string name = string("LinearIteratorProduct").append(suffix);
		reg.add_class_<T,TBase>(name, grp,
						"Linear Iterator consisting of several LinearIterations")
				.add_constructor()
				.template add_constructor<void (*)(const std::vector<SmartPtr<ILinearIterator<vector_type,vector_type> > >&)>()
				.add_method("add_iterator",static_cast<void (T::*)(SmartPtr<TBase>)>(&T::add_iterator),
					"", "add iterator", "sets iterator")
				.add_method("add_iterator",static_cast<void (T::*)(SmartPtr<TBase>,size_t nr)>(&T::add_iterator),
					"", "add iterator", "sets iterator")	
				.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "LinearIteratorProduct", tag);
	}

//	LinearIteratorSum
	{
		typedef LinearIteratorSum<vector_type, vector_type> T;
		typedef ILinearIterator<vector_type> TBase;
		string name = string("LinearIteratorSum").append(suffix);
		reg.add_class_<T,TBase>(name, grp,
						"Linear Iterator consisting of several LinearIterations")
				.add_constructor()
				.template add_constructor<void (*)(const std::vector<SmartPtr<ILinearIterator<vector_type,vector_type> > >&)>()
				.add_method("add_iterator",static_cast<void (T::*)(SmartPtr<TBase>)>(&T::add_iterator),
					"", "add iterator", "sets iterator")
				.add_method("add_iterator",static_cast<void (T::*)(SmartPtr<TBase>,size_t nr)>(&T::add_iterator),
					"", "add iterator", "sets iterator")	
				.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "LinearIteratorSum", tag);
	}

//	Vanka
	{
		typedef Vanka<TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("Vanka").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Vanka Preconditioner")
		.add_constructor()
		.add_method("set_relax", &T::set_relax, "", "relax")
		.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "Vanka", tag);
	}

//	Diag Vanka
	{
		typedef DiagVanka<TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("DiagVanka").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Diagonal Vanka Preconditioner")
		.add_constructor()
		.add_method("set_relax", &T::set_relax, "", "relax")
		.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "DiagVanka", tag);
	}

// 	AgglomeratingPreconditioner
	{
		typedef AgglomeratingIterator<TAlgebra> T;
		typedef ILinearIterator<vector_type> TBase;
		string name = string("AgglomeratingIterator").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "AgglomeratingIterator")
			.template add_constructor<void (*)(SmartPtr<ILinearIterator<vector_type> > )>("pLinIterator")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "AgglomeratingIterator", tag);
	}
}
	

}; // end Functionality

// end group precond_bridge
/// \}

}// end Preconditioner

/// \addtogroup precond_bridge
void RegisterBridge_Preconditioner(Registry& reg, string grp)
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
