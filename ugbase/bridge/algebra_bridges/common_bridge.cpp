/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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

// operator interfaces
#include "lib_algebra/operator/fixed_convergence_check.h"
#include "lib_algebra/operator/interface/operator.h"
#include "lib_algebra/operator/interface/matrix_operator.h"
#include "lib_algebra/operator/interface/matrix_operator_inverse.h"
#include "lib_algebra/operator/interface/preconditioner.h"
#include "lib_algebra/operator/interface/preconditioned_linear_operator_inverse.h"
#include "lib_algebra/operator/interface/operator_inverse.h"
#include "lib_algebra/operator/interface/linear_iterator.h"

#include "lib_algebra/operator/debug_writer.h"

// operator util
#include "lib_algebra/operator/preconditioner/iterator_product.h"
#include "lib_algebra/operator/operator_util.h"
#include "lib_algebra/operator/vector_writer.h"
#include "../util_overloaded.h"

#include "bridge_mat_vec_operations.h"
#include "matrix_diagonal.h"

#include "lib_algebra/operator/energy_convergence_check.h"

using namespace std;

namespace ug{

namespace bridge{
namespace AlgebraCommon{

/// \addtogroup algebracommon_bridge
/// \{



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
 * @param grp				group for sorting of functionality
 */
template <typename TAlgebra>
static void Algebra(Registry& reg, string grp)
{

	using vector_type = typename TAlgebra::vector_type;
	using matrix_type = typename TAlgebra::matrix_type;

//	suffix and tag
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();

// MatrixDiagonal
	{
		using T2 = MatrixOperator<matrix_type, vector_type>;
		using T = MatrixDiagonal<matrix_type, vector_type>;
		string name = string("MatrixDiagonal").append(suffix);
		reg.add_class_<T, ILinearOperator<vector_type> >(name, grp)
				.ADD_CONSTRUCTOR( (SmartPtr<T2> mo) )()
				.set_construct_as_smart_pointer(true);

		reg.add_class_to_group(name, "MatrixDiagonal", tag);
	}

// MatrixDiagonalInverse
	{
		using T2 = MatrixOperator<matrix_type, vector_type>;
		using T = MatrixDiagonalInverse<matrix_type, vector_type>;
		string name = string("MatrixDiagonalInverse").append(suffix);
		reg.add_class_<T, ILinearOperator<vector_type> >(name, grp)
				.ADD_CONSTRUCTOR( (SmartPtr<T2> mo) )()
				.set_construct_as_smart_pointer(true);

		reg.add_class_to_group(name, "MatrixDiagonalInverse", tag);
	}
//	Vector
	{
		string name = string("Vector").append(suffix);
		reg.add_class_<vector_type>(name, grp)
		.add_constructor()
		.add_method("set|hide=true", (void (vector_type::*)(number))&vector_type::set,
								"Success", "Number")
		.add_method("size|hide=true", (size_t (vector_type::*)() const)&vector_type::size,
								"Size", "")
		.add_method("set_random|hide=true", (void (vector_type::*)(number, number))&vector_type::set_random,
								"Success", "Number")
		.add_method("print|hide=true", &vector_type::p)
#ifdef UG_PARALLEL
		.add_method("check_storage_type", &vector_type::check_storage_type)
		.add_method("enforce_consistent_type", &vector_type::enforce_consistent_type)
#endif
		.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "Vector", tag);

		reg.add_function("VecAssign",
				(void (*)(vector_type&,const vector_type &)) &VecAssign<vector_type>, grp, "", "dest#vec", "dest <- vec");
		reg.add_function("VecScaleAssign",
				(void (*)(vector_type&, number, const vector_type &)) &VecScaleAssign<vector_type>, grp, "", "dest#alpha#vec", "dest <- alpha * vec");
		reg.add_function("VecScaleAdd2", /*(void (*)(vector_type&, number, const vector_type&, number, const vector_type &)) */
				&VecScaleAdd2<vector_type>, grp, "alpha1*vec1 + alpha2*vec2",
				"dest#alpha1#vec1#alpha2#vec2");
		reg.add_function("VecScaleAdd3", /*(void (*)(vector_type&, number, const vector_type&, number, const vector_type &, number, const vector_type &))*/
				&VecScaleAdd3<vector_type>, grp, "alpha1*vec1 + alpha2*vec2 + alpha3*vec3",
				"dest#alpha1#vec1#alpha2#vec2#alpha3#vec3");
		reg.add_function("VecProd", &VecProd2<TAlgebra>);
		reg.add_function("VecProd", &VecProdOp<vector_type>);
		reg.add_function("VecProd", &VecScaleAddProd1<TAlgebra>);
		reg.add_function("VecProd", &VecScaleAddProd2<TAlgebra>);
		reg.add_function("VecNorm", &VecNorm<TAlgebra>);
		reg.add_function("VecMaxNorm", &VecMaxNorm<TAlgebra>);
		reg.add_function("VecNorm", &VecScaleAddNorm<TAlgebra>);
		reg.add_function("VecHadamardProd", (void (*)(vector_type&, const vector_type &, const vector_type &))
				&VecHadamardProd<vector_type>, grp, "", "dst#vec1#vec2", "vec1 * vec2 (elementwise)");
		reg.add_function("VecExp", (void (*)(vector_type&, const vector_type &))
				&VecExp<vector_type>, grp, "", "dst#vec", "exp(vec) (elementwise)");
		reg.add_function("VecLog", (void (*)(vector_type&, const vector_type &))
				&VecLog<vector_type>, grp, "", "dst#vec", "log(vec) (elementwise)");
	}

//	VecScaleAddClass
	{
		string name = string("VecScaleAddClass").append(suffix);
		using T = VecScaleAddClass<TAlgebra>;
		reg.add_class_<T>(name, grp)
		.add_constructor()
		.ADD_CONSTRUCTOR( (double scale1, SmartPtr<vector_type> v1) ) ()
		.ADD_CONSTRUCTOR( (double scale1, SmartPtr<vector_type> v1, double scale2, SmartPtr<vector_type> v2) ) ()
		.ADD_CONSTRUCTOR( (double scale, SmartPtr<VecScaleAddClass<TAlgebra > > vsac, double scale1, SmartPtr<vector_type> v1) ) ()
		.ADD_CONSTRUCTOR( (double scale1, SmartPtr<vector_type> v1, double scale, SmartPtr<VecScaleAddClass<TAlgebra > > vsac) ) ()
		.ADD_CONSTRUCTOR( (double scale, SmartPtr<VecScaleAddClass<TAlgebra > > vsac) ) ()
		.add_method("eval", &T::eval)
		.add_method("assign", &T::assign)
		.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "VecScaleAddClass", tag);

		reg.add_function("Eval", &Eval<TAlgebra>, grp);
		reg.add_function("Assign", &Assign<TAlgebra>, grp);
	}
//	Matrix
	{
		string name = string("Matrix").append(suffix);
		reg.add_class_<matrix_type>(name, grp)
			.add_constructor()
			.add_method("print|hide=true", &matrix_type::p)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "Matrix", tag);
	}

//	ApplyLinearSolver
	{
		reg.add_function( "ApplyLinearSolver",
						  &ApplyLinearSolver<vector_type>, grp);
	}

//  Vector Debug Writer (abstract base class)
	{
		using T = IVectorDebugWriter<vector_type>;
		string name = string("IVectorDebugWriter").append(suffix);
		reg.add_class_<T>(name, grp)
			.add_method("get_context", static_cast<SmartPtr<DebugWriterContext> (T::*) ()>(&T::get_context), "Gets the debugger writer context", "")
			.add_method("set_context", &T::set_context, "Sets the debugger writer context", "context")
		    .add_method("set_base_dir", &T::set_base_dir, "Sets the base directory for output", "dir")
		    .add_method("enter_section", &T::enter_section, "Enters a debugging section", "dirName")
		    .add_method("leave_section", &T::leave_section, "Leaves the current debugging section", "")
			;
		reg.add_class_to_group(name, "IVectorDebugWriter", tag);
	}

// Debug Writer (abstract base class)
	{
		using T = IDebugWriter<TAlgebra>;
		using TBase = IVectorDebugWriter<vector_type>;
		string name = string("IDebugWriter").append(suffix);
		reg.add_class_<T, TBase>(name, grp);
		reg.add_class_to_group(name, "IDebugWriter", tag);
	}

//  VectorDebugWritingObject
	{
		using T = VectorDebugWritingObject<vector_type>;
		string name = string("VectorDebugWritingObject").append(suffix);
		reg.add_class_<T>(name, grp)
			.template add_constructor<void (*) (SmartPtr<IVectorDebugWriter<vector_type> >)> ("VectorDebugWriter")
			.add_method("set_debug", &T::set_debug,  "", "dbgWriter", "sets a debug writer")
			.add_method("vector_debug_writer", static_cast<SmartPtr<IVectorDebugWriter<vector_type> > (T::*)()>(&T::vector_debug_writer))
			.add_method("vector_debug_writer", static_cast<ConstSmartPtr<IVectorDebugWriter<vector_type> > (T::*)() const>(&T::vector_debug_writer))
			.add_method("write_debug", static_cast<void (T::*) (const vector_type&, const char*)> (&T::write_debug))
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "VectorDebugWritingObject", tag);
	}

//  DebugWritingObject
	{
		using T = DebugWritingObject<TAlgebra>;
		using TBase = VectorDebugWritingObject<vector_type>;
		string name = string("DebugWritingObject").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_method("set_debug", static_cast<void (T::*)(SmartPtr<IDebugWriter<TAlgebra> >)>(&T::set_debug), "", "dbgWriter", "sets a debug writer")
			.add_method("debug_writer", static_cast<SmartPtr<IDebugWriter<TAlgebra> > (T::*)()>(&T::debug_writer))
			.add_method("debug_writer", static_cast<ConstSmartPtr<IDebugWriter<TAlgebra> > (T::*)() const>(&T::debug_writer));
		reg.add_class_to_group(name, "DebugWritingObject", tag);
	}

// IVectorWriter (abstract base class)
	{
		using T = IVectorWriter<vector_type>;
		string name = string("IVectorWriter").append(suffix);
		reg.add_class_<T>(name, grp)
				.add_method("update", &T::update, "", "v", "updates the vector v");
		reg.add_class_to_group(name, "IVectorWriter", tag);
	}

/////////////////////////
//	Base Classes
/////////////////////////

//	ILinearOperator
	{
		using T = ILinearOperator<vector_type>;
		using TBase = IOperator<vector_type>;
		string name = string("ILinearOperator").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_method("init", static_cast<void (T::*)()>(&T::init))
			.add_method("apply", &T::apply, "f#u", "", "calculates f = Op(u)")
			.add_method("apply_sub", &T::apply_sub, "f#u", "", "calculates f -= Op(u)");
		reg.add_class_to_group(name, "ILinearOperator", tag);
	}

// 	MatrixOperator
	{
		using TBase = ILinearOperator<vector_type>;
		using T = MatrixOperator<matrix_type, vector_type>;
		string name = string("MatrixOperator").append(suffix);
		reg.add_class_<T, TBase, matrix_type>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "MatrixOperator", tag);
	}

//	ILinearIterator
	{
		using T = ILinearIterator<vector_type>;
		string name = string("ILinearIterator").append(suffix);
		reg.add_class_<T>(name, grp)
			.add_method("set_damp", static_cast<void (T::*)(number)>(&T::set_damp), "", "damp", "set the damping to a number")
			.add_method("set_damp", static_cast<void (T::*)(SmartPtr<IDamping<vector_type> >)>(&T::set_damp), "", "damp", "set the damping to a damping object")
			.add_method("config_string", &T::config_string, "strConfiguration", "", "string to display configuration of the linear iterator")
			.add_method("clone", &T::clone, "SmartPointer to a copy of this object", "", "returns a clone of the object which can be modified independently")
			.add_method("apply", &T::apply)
			.add_method("apply_update_defect", &T::apply_update_defect)
			.add_method("init", OVERLOADED_METHOD_PTR(bool, T, init, (SmartPtr<ILinearOperator<vector_type,vector_type> > L) ))
			.add_method("init", OVERLOADED_METHOD_PTR(bool, T, init, (SmartPtr<ILinearOperator<vector_type,vector_type> > L, const vector_type &u) ))
			.add_method("name", &T::name);
		reg.add_class_to_group(name, "ILinearIterator", tag);
	}

//	IPreconditioner
	{
		using T = IPreconditioner<TAlgebra>;
		using TBase = ILinearIterator<vector_type>;
		using TBase2 = DebugWritingObject<TAlgebra>;
		string name = string("IPreconditioner").append(suffix);
		reg.add_class_<T, TBase, TBase2>(name, grp);
		reg.add_class_to_group(name, "IPreconditioner", tag);
	}

//	ILinearOperatorInverse
	{
		using T = ILinearOperatorInverse<vector_type>;
		using TBase = ILinearIterator<vector_type>;
		string name = string("ILinearOperatorInverse").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_method("init", static_cast<bool (T::*)(SmartPtr<ILinearOperator<vector_type> >)>(&T::init))
			.add_method("init", static_cast<bool (T::*)(SmartPtr<ILinearOperator<vector_type> >,const vector_type&)>(&T::init))
			.add_method("apply_return_defect", &T::apply_return_defect, "Success", "u#f",
					"Solve A*u = f, such that u = A^{-1} f by iterating u := u + B(f - A*u),  f := f - A*u becomes new defect")
			.add_method("apply", &T::apply, "Success", "u#f", "Solve A*u = f, such that u = A^{-1} f by iterating u := u + B(f - A*u), f remains constant")
			.add_method("set_convergence_check", &T::set_convergence_check)
			.add_method("convergence_check", static_cast<ConstSmartPtr<IConvergenceCheck<vector_type> > (T::*)() const>(&T::convergence_check))
			.add_method("defect", &T::defect, "the current defect")
			.add_method("step", &T::step, "the current number of steps")
			.add_method("reduction", &T::reduction, "the current relative reduction")
			.add_method("config_string", &T::config_string);
		reg.add_class_to_group(name, "ILinearOperatorInverse", tag);
	}

//	IPreconditionedLinearOperatorInverse
	{
		using T = IPreconditionedLinearOperatorInverse<vector_type>;
		using TBase = ILinearOperatorInverse<vector_type>;
		using TBase2 = VectorDebugWritingObject<vector_type>;
		string name = string("IPreconditionedLinearOperatorInverse").append(suffix);
		reg.add_class_<T, TBase, TBase2>(name, grp)
			.add_method("set_preconditioner", &T::set_preconditioner,
						"", "Preconditioner")
			.add_method("set_compute_fresh_defect_when_finished", &T::set_compute_fresh_defect_when_finished);
		reg.add_class_to_group(name, "IPreconditionedLinearOperatorInverse", tag);
	}

//	IMatrixOperatorInverse
	{
		using TBase = ILinearOperatorInverse<vector_type>;
		using T = IMatrixOperatorInverse<matrix_type, vector_type>;
		string name = string("IMatrixOperatorInverse").append(suffix);
		reg.add_class_<T, TBase>(name, grp);
		reg.add_class_to_group(name, "IMatrixOperatorInverse", tag);
	}

//	IOperator
	{
		using T = IOperator<vector_type>;
		string name = string("IOperator").append(suffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, "IOperator", tag);
	}

//	IOperatorInverse
	{
		using T = IOperatorInverse<vector_type>;
		string name = string("IOperatorInverse").append(suffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, "IOperatorInverse", tag);
	}


// 	IConvergenceCheck
	{
		using T = IConvergenceCheck<vector_type>;
		string name = string("IConvergenceCheck").append(suffix);
		reg.add_class_<T>(name, grp)
				.add_method("config_string", &T::config_string)
				.add_method("defect", &T::defect, "defect", "", "returns the current defect")
					.add_method("step", &T::step, "step", "", "returns the current number of steps")
					.add_method("reduction", &T::reduction, "reduction", "", "returns the current relative reduction")
					.add_method("iteration_ended", &T::iteration_ended)
					.add_method("avg_rate", &T::avg_rate, "", "returns the average convergence rate")
					;
		reg.add_class_to_group(name, "IConvergenceCheck", tag);
	}

// 	StandardConvCheck
	{
		using T = StdConvCheck<vector_type>;
		using TBase = IConvergenceCheck<vector_type>;
		string name = string("ConvCheck").append(suffix);
		reg.add_class_<T, TBase>(name, grp, "Convergence Check")
			.add_constructor()
			.template add_constructor<void (*)(int, number, number, bool)>
							("Maximum Steps|default|min=0;value=100#"
							 "Minimum Defect|default|min=0D;value=1e-10#"
							 "Relative Reduction|default|min=0D;value=1e-12#"
							 "Verbosity")
			.template add_constructor<void (*)(int, number, number, bool, bool)>
							("Maximum Steps|default|min=0;value=100#"
							 "Minimum Defect|default|min=0D;value=1e-10#"
							 "Relative Reduction|default|min=0D;value=1e-12#"
							 "Verbosity"
							 "supress unsuccessful return")
			.template add_constructor<void (*)(int, number, number)>
							("Maximum Steps|default|min=0;value=100#"
							 "Minimum Defect|default|min=0D;value=1e-10#"
							 "Relative Reduction|default|min=0D;value=1e-12")
			.add_method("set_maximum_steps", &T::set_maximum_steps, "", "Maximum Steps|default|min=0;value=100", "maximum number of steps to do")
			.add_method("set_minimum_defect", &T::set_minimum_defect, "", "Minimum Defect|default|min=0D;value=1e-10")
			.add_method("set_reduction", &T::set_reduction,	"", "Relative Reduction|default|min=0D;value=1e-12")
			.add_method("set_verbose", &T::set_verbose,	"", "Verbosity")
			.add_method("set_supress_unsuccessful", &T::set_supress_unsuccessful,"", "supress false return")
			.add_method("defect", &T::defect)
			.add_method("get_defect", &T::get_defect)
			.add_method("previous_defect", &T::previous_defect)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ConvCheck", tag);
	}


	{
		using T = EnergyConvCheck<vector_type>;
		using TBase = IConvergenceCheck<vector_type>;
		string name = string("EnergyConvCheck").append(suffix);
		reg.add_class_<T, TBase>(name, grp, "Energy Convergence Check")
				.add_method("set_linear_operator", &T::set_linear_operator)
				.template add_constructor<void (*)(int, number, number, bool)>
											("Maximum Steps|default|min=0;value=100#"
											 "Minimum Defect|default|min=0D;value=1e-10#"
											 "Relative Reduction|default|min=0D;value=1e-12#"
											 "Verbosity")
				.template add_constructor<void (*)(int, number, number, bool, bool)>
								("Maximum Steps|default|min=0;value=100#"
								 "Minimum Defect|default|min=0D;value=1e-10#"
								 "Relative Reduction|default|min=0D;value=1e-12#"
								 "Verbosity"
								 "supress unsuccessful return")
				.template add_constructor<void (*)(int, number, number)>
								("Maximum Steps|default|min=0;value=100#"
								 "Minimum Defect|default|min=0D;value=1e-10#"
								 "Relative Reduction|default|min=0D;value=1e-12")
				.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "EnergyConvCheck", tag);
	}

// 	FixedConvCheck
	{
		using T = FixedConvergenceCheck<vector_type>;
		using TBase = IConvergenceCheck<vector_type>;
		string name = string("FixedConvergenceCheck").append(suffix);
		reg.add_class_<T, TBase>(name, grp, "Convergence Check")
	//		.add_constructor()
			.template add_constructor<void (*)(number)>("")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "FixedConvergenceCheck", tag);
	}

}


/**
 * Function called for the registration of Domain and Algebra independent parts.
 * All Functions and Classes not depending on Domain and Algebra
 * are to be placed here when registering.
 *
 * @param reg				registry
 * @param grp				group for sorting of functionality
 */
static void Common(Registry& reg, string grp)
{

// IPositionProvider (abstract base class)
	{
		reg.add_class_<IPositionProvider<1> >("IPositionProvider1d", grp);
		reg.add_class_<IPositionProvider<2> >("IPositionProvider2d", grp);
		reg.add_class_<IPositionProvider<3> >("IPositionProvider3d", grp);
		reg.add_class_to_group("IPositionProvider1d", "IPositionProvider", GetDimensionTag<1>());
		reg.add_class_to_group("IPositionProvider2d", "IPositionProvider", GetDimensionTag<2>());
		reg.add_class_to_group("IPositionProvider3d", "IPositionProvider", GetDimensionTag<3>());
	}
	
// Debug Writer Context
	{
		using T = DebugWriterContext;
		reg.add_class_<T>("DebugWriterContext", grp)
			.add_constructor()
			.add_method("set_base_dir", &T::set_base_dir)
			.add_method("get_base_dir", &T::get_base_dir)
			.add_method("enter_section", &T::enter_section)
			.add_method("leave_section", &T::leave_section)
			.add_method("print_message", &T::leave_section)
			.add_method("compose_file_path", &T::leave_section)
			.set_construct_as_smart_pointer(true);
	}
}

}; // end Functionality

// end group algebracommon_bridge
/// \}

}// end AlgebraCommon

/// \addtogroup algebracommon_bridge
void RegisterBridge_AlgebraCommon(Registry& reg, string grp)
{
	grp.append("/Algebra");
	using Functionality = AlgebraCommon::Functionality;

	try{
		RegisterCommon<Functionality>(reg,grp);
		RegisterAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

} // end namespace bridge
} // end namespace ug
