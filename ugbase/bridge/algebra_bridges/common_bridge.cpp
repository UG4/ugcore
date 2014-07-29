/*
 * common_bridge.cpp
 *
 *  Created on: 11.10.2010
 *      Author: andreasvogel
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
 * @param parentGroup		group for sorting of functionality
 */
template <typename TAlgebra>
static void Algebra(Registry& reg, string grp)
{
//	typedefs for this algebra
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;

//	suffix and tag
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();

// MatrixDiagonal
	{
		typedef MatrixOperator<matrix_type, vector_type> T2;
		typedef MatrixDiagonal<matrix_type, vector_type> T;
		string name = string("MatrixDiagonal").append(suffix);
		reg.add_class_<T, ILinearOperator<vector_type> >(name, grp)
				.ADD_CONSTRUCTOR( (SmartPtr<T2> mo) )()
				.set_construct_as_smart_pointer(true);

		reg.add_class_to_group(name, "MatrixDiagonal", tag);
	}

// MatrixDiagonalInverse
	{
		typedef MatrixOperator<matrix_type, vector_type> T2;
		typedef MatrixDiagonalInverse<matrix_type, vector_type> T;
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

		reg.add_function("VecScaleAssign",
				(void (*)(vector_type&, number, const vector_type &)) &VecScaleAssign<vector_type>
		, grp);
		reg.add_function("VecAssign",
				(void (*)(vector_type&,const vector_type &)) &VecAssign<vector_type>, grp, "", "dest#vec", "calculates dest <- vec");
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
		reg.add_function("VecNorm", &VecScaleAddNorm<TAlgebra>);
	}

//	VecScaleAddClass
	{
		string name = string("VecScaleAddClass").append(suffix);
		typedef VecScaleAddClass<TAlgebra> T;
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
		typedef IVectorDebugWriter<vector_type> T;
		string name = string("IVectorDebugWriter").append(suffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, "IVectorDebugWriter", tag);
	}

// Debug Writer (abstract base class)
	{
		typedef IDebugWriter<TAlgebra> T;
		typedef IVectorDebugWriter<vector_type> TBase;
		string name = string("IDebugWriter").append(suffix);
		reg.add_class_<T, TBase>(name, grp);
		reg.add_class_to_group(name, "IDebugWriter", tag);
	}

//  VectorDebugWritingObject
	{
		typedef VectorDebugWritingObject<vector_type> T;
		string name = string("VectorDebugWritingObject").append(suffix);
		reg.add_class_<T>(name, grp)
			.add_method("set_debug", &T::set_debug,  "", "dbgWriter", "sets a debug writer")
			.add_method("vector_debug_writer", static_cast<SmartPtr<IVectorDebugWriter<vector_type> > (T::*)()>(&T::vector_debug_writer))
			.add_method("vector_debug_writer", static_cast<ConstSmartPtr<IVectorDebugWriter<vector_type> > (T::*)() const>(&T::vector_debug_writer));
		reg.add_class_to_group(name, "VectorDebugWritingObject", tag);
	}

//  DebugWritingObject
	{
		typedef DebugWritingObject<TAlgebra> T;
		typedef VectorDebugWritingObject<vector_type> TBase;
		string name = string("DebugWritingObject").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_method("set_debug", static_cast<void (T::*)(SmartPtr<IDebugWriter<TAlgebra> >)>(&T::set_debug), "", "dbgWriter", "sets a debug writer")
			.add_method("debug_writer", static_cast<SmartPtr<IDebugWriter<TAlgebra> > (T::*)()>(&T::debug_writer))
			.add_method("debug_writer", static_cast<ConstSmartPtr<IDebugWriter<TAlgebra> > (T::*)() const>(&T::debug_writer));
		reg.add_class_to_group(name, "DebugWritingObject", tag);
	}

// IVectorWriter (abstract base class)
	{
		typedef IVectorWriter<vector_type> T;
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
		typedef ILinearOperator<vector_type> T;
		typedef IOperator<vector_type> TBase;
		string name = string("ILinearOperator").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_method("init", static_cast<void (T::*)()>(&T::init))
			.add_method("apply", &T::apply, "f#u", "", "calculates f = Op(u)")
			.add_method("apply_sub", &T::apply_sub, "f#u", "", "calculates f -= Op(u)");
		reg.add_class_to_group(name, "ILinearOperator", tag);
	}

// 	MatrixOperator
	{
		typedef ILinearOperator<vector_type> TBase;
		typedef MatrixOperator<matrix_type, vector_type> T;
		string name = string("MatrixOperator").append(suffix);
		reg.add_class_<T, TBase, matrix_type>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "MatrixOperator", tag);
	}

//	ILinearIterator
	{
		typedef ILinearIterator<vector_type> T;
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
		typedef IPreconditioner<TAlgebra> T;
		typedef ILinearIterator<vector_type>  TBase;
		typedef DebugWritingObject<TAlgebra> TBase2;
		string name = string("IPreconditioner").append(suffix);
		reg.add_class_<T, TBase, TBase2>(name, grp);
		reg.add_class_to_group(name, "IPreconditioner", tag);
	}

//	ILinearOperatorInverse
	{
		typedef ILinearOperatorInverse<vector_type> T;
		string name = string("ILinearOperatorInverse").append(suffix);
		reg.add_class_<T>(name, grp)
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
		typedef IPreconditionedLinearOperatorInverse<vector_type> T;
		typedef ILinearOperatorInverse<vector_type> TBase;
		typedef VectorDebugWritingObject<vector_type> TBase2;
		string name = string("IPreconditionedLinearOperatorInverse").append(suffix);
		reg.add_class_<T, TBase, TBase2>(name, grp)
			.add_method("set_preconditioner", &T::set_preconditioner,
						"", "Preconditioner")
			.add_method("set_compute_fresh_defect_when_finished", &T::set_compute_fresh_defect_when_finished);
		reg.add_class_to_group(name, "IPreconditionedLinearOperatorInverse", tag);
	}

//	IMatrixOperatorInverse
	{
		typedef ILinearOperatorInverse<vector_type>  TBase;
		typedef IMatrixOperatorInverse<matrix_type, vector_type> T;
		string name = string("IMatrixOperatorInverse").append(suffix);
		reg.add_class_<T, TBase>(name, grp);
		reg.add_class_to_group(name, "IMatrixOperatorInverse", tag);
	}

//	IOperator
	{
		typedef IOperator<vector_type> T;
		string name = string("IOperator").append(suffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, "IOperator", tag);
	}

//	IOperatorInverse
	{
		typedef IOperatorInverse<vector_type> T;
		string name = string("IOperatorInverse").append(suffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, "IOperatorInverse", tag);
	}


// 	IConvergenceCheck
	{
		typedef IConvergenceCheck<vector_type> T;
		string name = string("IConvergenceCheck").append(suffix);
		reg.add_class_<T>(name, grp)
				.add_method("config_string", &T::config_string);
		reg.add_class_to_group(name, "IConvergenceCheck", tag);
	}

// 	StandardConvCheck
	{
		typedef StdConvCheck<vector_type> T;
		typedef IConvergenceCheck<vector_type> TBase;
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
			.add_method("defect", &T::defect, "defect", "", "returns the current defect")
			.add_method("step", &T::step, "step", "", "returns the current number of steps")
			.add_method("reduction", &T::reduction, "reduction", "", "returns the current relative reduction")
			.add_method("iteration_ended", &T::iteration_ended)
			.add_method("previous_defect", &T::previous_defect)
			.add_method("avg_rate", &T::avg_rate, "", "returns the average convergence rate")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ConvCheck", tag);
	}

}


/**
 * Function called for the registration of Domain and Algebra independent parts.
 * All Functions and Classes not depending on Domain and Algebra
 * are to be placed here when registering.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
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
}

}; // end Functionality

// end group algebracommon_bridge
/// \}

}// end AlgebraCommon

/// \addtogroup algebracommon_bridge
void RegisterBridge_AlgebraCommon(Registry& reg, string grp)
{
	grp.append("/Algebra");
	typedef AlgebraCommon::Functionality Functionality;

	try{
		RegisterCommon<Functionality>(reg,grp);
		RegisterAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

} // end namespace bridge
} // end namespace ug
