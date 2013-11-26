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
#include "common/serialization.h"
using namespace std;

namespace ug{

/*
 *

		restartfilename = "uNewtonSolution_T_"..step.."_p_"..GetProcessRank()..".ug4vector"

  		restartStep = util.GetParamNumber("-restartStep", 0)
		-- apply newton solver
		if step < restartStep then
			ReadFromFile(u, restartfilename)
		else
			-- prepare newton solver
			if newtonSolver:prepare(u) == false then
				print ("Newton solver failed at step "..step.."."); exit();
			end
			if newtonSolver:apply(u) == false then
				print ("Newton solver failed at step "..step.."."); exit();
			end
			SaveToFile(u, restartfilename)
		end
 */

#ifdef UG_PARALLEL
template<typename T, class TOStream>
void Serialize(TOStream &buf, const ParallelVector<T> &v)
{
	uint t= v.get_storage_mask();
	Serialize(buf, t);
	Serialize(buf, *dynamic_cast<const T*>(&v));
}

template<typename T, class TIStream>
void Deserialize(TIStream &buf, ParallelVector<T> &v)
{
	uint t = Deserialize<uint>(buf);
	v.set_storage_type(t);
	Deserialize(buf, *dynamic_cast<T*>(&v));
}
#endif

template<typename T>
void SaveToFile(const T &v, std::string filename)
{
	fstream f(filename.c_str(), ios::out);
	Serialize(f, v);
}

template<typename T>
void ReadFromFile(T &v, std::string filename)
{
	fstream f(filename.c_str(), ios::in);
	Deserialize(f, v);
}


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


	reg.add_function("SaveToFile", static_cast<void (*)(const vector_type &, std::string)>(&SaveToFile<vector_type>), grp);
	reg.add_function("ReadFromFile", static_cast<void (*)(vector_type &, std::string)>(&ReadFromFile<vector_type>), grp);


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

//	ILU
	{
		typedef ILU<TAlgebra> T;
		typedef IPreconditioner<TAlgebra> TBase;
		string name = string("ILU").append(suffix);
		reg.add_class_<T,TBase>(name, grp, "Incomplete LU Decomposition")
			.add_constructor()
			.add_method("set_beta", &T::set_beta, "", "beta")
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

//	IProjPreconditioner
	{
		typedef IProjPreconditioner<TAlgebra> T;
		typedef ILinearIterator<vector_type> TBase;
		string name = string("IProjPreconditioner").append(suffix);
		reg.add_class_<T,TBase>(name, grp)
			.add_method("set_obstacle_constraint", &T::set_obstacle_constraint,
				"", "obstacle constraint", "sets an obstacle constraint")
			.add_method("set_sor_relax", &T::set_sor_relax,
			"", "sor relaxation", "sets sor relaxation parameter");
		reg.add_class_to_group(name, "IProjPreconditioner", tag);
	}

//	ProjGaussSeidel
	{
		typedef ProjGaussSeidel<TAlgebra> T;
		typedef IProjPreconditioner<TAlgebra> TBase;
		string name = string("ProjGaussSeidel").append(suffix);
		reg.add_class_<T,TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ProjGaussSeidel", tag);
	}

//	ProjBackwardGaussSeidel
	{
		typedef ProjBackwardGaussSeidel<TAlgebra> T;
		typedef IProjPreconditioner<TAlgebra> TBase;
		string name = string("ProjBackwardGaussSeidel").append(suffix);
		reg.add_class_<T,TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ProjBackwardGaussSeidel", tag);
	}

//	ProjSymmetricGaussSeidel
	{
		typedef ProjSymmetricGaussSeidel<TAlgebra> T;
		typedef IProjPreconditioner<TAlgebra> TBase;
		string name = string("ProjSymmetricGaussSeidel").append(suffix);
		reg.add_class_<T,TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ProjSymmetricGaussSeidel", tag);
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
