/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#include "lib_algebra/ordering_strategies/algorithms/IOrderingAlgorithm.h"

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


template <typename TDomain, typename TAlgebra>
static void DomainAlgebra(Registry& reg, string grp)
{
	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

//	AssembledTransformingSmoother
	{
		typedef AssembledTransformingSmoother<TDomain, TAlgebra> T;
		typedef ILinearIterator<typename TAlgebra::vector_type> TBase;
		typedef DebugWritingObject<TAlgebra> TBase2;
		string name = string("AssembledTransformingSmoother").append(suffix);
		reg.add_class_<T, TBase, TBase2>(name, grp)
            .ADD_CONSTRUCTOR((SmartPtr<IAssemble<TAlgebra> >,
                          SmartPtr<ILinearIterator<typename TAlgebra::vector_type> >,
                          SmartPtr<IAssemble<TAlgebra> >))
                          ("TrafoSystemAss, TrafoSystemSmoother, RightTrafoAss")
            .ADD_CONSTRUCTOR((SmartPtr<IAssemble<TAlgebra> >,
                              SmartPtr<ILinearIterator<typename TAlgebra::vector_type> >,
                              SmartPtr<IAssemble<TAlgebra> >,
                              SmartPtr<ILinearIterator<typename TAlgebra::vector_type> >))
                            ("TrafoSystemAss, TrafoSystemSmoother, RightTrafoAss, RightTrafoSmoother")
            .set_construct_as_smart_pointer(true);
        reg.add_class_to_group(name, "AssembledTransformingSmoother", tag);
	}

}

template<typename TAlgebra, typename TGSType>
static void RegisterBlockGaussSeidel(Registry& reg, string grp, string name)
{
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();

	typedef TGSType T;
	typedef IPreconditioner<TAlgebra> TBase;
	string namesuffix = name+suffix;
	reg.add_class_<T,TBase>(namesuffix, grp, name)
		.add_constructor()
		.ADD_CONSTRUCTOR( (int) )("depth")
		.add_method("set_depth", &T::set_depth)
		.set_construct_as_smart_pointer(true);

	reg.add_class_to_group(namesuffix, name, tag);
}

template<typename TAlgebra, bool forward, bool backward>
static void RegisterBlockGaussSeidelIterative(Registry& reg, string grp, string name)
{
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();

	typedef BlockGaussSeidelIterative<TAlgebra, forward, backward> T;
	typedef IPreconditioner<TAlgebra> TBase;
	string namesuffix = name+suffix;
	reg.add_class_<T,TBase>(namesuffix, grp, name)
		.add_constructor()
		. ADD_CONSTRUCTOR( (int, int) )("depth#steps")
		.add_method("set_depth", &T::set_depth)
		.add_method("set_iterative_steps", &T::set_iterative_steps)
		.set_construct_as_smart_pointer(true);
	reg.add_class_to_group(namesuffix, name, tag);
}

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
			.add_method("enable_consistent_interfaces", &T::enable_consistent_interfaces, "", "enable", "makes the matrix and defect consistent at the proc. interfaces")
			.add_method("enable_overlap", &T::enable_overlap, "", "enable", "Enables matrix overlap. This also means that interfaces are consistent.")
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
		typedef GaussSeidelBase<TAlgebra> TBase;
		string name = string("BackwardGaussSeidel").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Backward Gauss Seidel Preconditioner")
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "BackwardGaussSeidel", tag);
	}

//	BlockGaussSeidel
	{
		RegisterBlockGaussSeidel<TAlgebra, BlockGaussSeidel<TAlgebra, true, false> >(reg, grp, "BlockGaussSeidel");
		RegisterBlockGaussSeidel<TAlgebra, BlockGaussSeidel<TAlgebra, false, true> >(reg, grp, "BackwardBlockGaussSeidel");
		RegisterBlockGaussSeidel<TAlgebra, BlockGaussSeidel<TAlgebra, true, true> >(reg, grp, "SymmetricBlockGaussSeidel");

		RegisterBlockGaussSeidel<TAlgebra, SparseBlockGaussSeidel<TAlgebra, true, false> >(reg, grp, "SparseBlockGaussSeidel");
		RegisterBlockGaussSeidel<TAlgebra, SparseBlockGaussSeidel2<TAlgebra, true, false> >(reg, grp, "SparseBlockGaussSeidel2");
		RegisterBlockGaussSeidel<TAlgebra, SparseBlockGaussSeidel<TAlgebra, false, true> >(reg, grp, "SparseBackwardBlockGaussSeidel");
		RegisterBlockGaussSeidel<TAlgebra, SparseBlockGaussSeidel<TAlgebra, true, true> >(reg, grp, "SparseSymmetricBlockGaussSeidel");

		RegisterBlockGaussSeidelIterative<TAlgebra, true, false>(reg, grp, "BlockGaussSeidelIterative");
		RegisterBlockGaussSeidelIterative<TAlgebra, false, true>(reg, grp, "BackwardBlockGaussSeidelIterative");
		RegisterBlockGaussSeidelIterative<TAlgebra, true, true>(reg, grp, "SymmetricBlockGaussSeidelIterative");
	}

//	ILU
	{
		typedef ILU<TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("ILU").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Incomplete LU Decomposition")
			.add_constructor()
			.add_method("set_beta", &T::set_beta, "", "beta")
			.add_method("set_sort_eps", &T::set_sort_eps, "", "eps")
			.add_method("set_inversion_eps", &T::set_inversion_eps, "", "eps")
			.add_method("set_ordering_algorithm", &T::set_ordering_algorithm, "", "",
						"sets an ordering algorithm")
			.add_method("set_sort", &T::set_sort, "", "bSort", "if bSort=true, use a cuthill-mckey sorting to reduce fill-in. default false")
			.add_method("set_disable_preprocessing", &T::set_disable_preprocessing, "", "disable",
						"set whether preprocessing (notably, LU factorization) is to be disabled - usable when the operator has not changed; use with care")
			.add_method("enable_consistent_interfaces", &T::enable_consistent_interfaces, "", "enable", "Make Matrix consistent for connections in interfaces.")
			.add_method("enable_overlap", &T::enable_overlap, "", "enable", "Enables matrix overlap. This also means that interfaces are consistent.")
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

// 	AgglomeratingIterator
	{
		typedef AgglomeratingIterator<TAlgebra> T;
		typedef ILinearIterator<vector_type> TBase;
		string name = string("AgglomeratingIterator").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "AgglomeratingIterator")
			.template add_constructor<void (*)(SmartPtr<ILinearIterator<vector_type> > )>("pLinIterator")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "AgglomeratingIterator", tag);
	}

// 	AgglomeratingPreconditioner
	{
		typedef AgglomeratingPreconditioner<TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("AgglomeratingPreconditioner").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "AgglomeratingPreconditioner")
			.template add_constructor<void (*)(SmartPtr<ILinearIterator<vector_type> > )>("pPreconditioner")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "AgglomeratingPreconditioner", tag);
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
