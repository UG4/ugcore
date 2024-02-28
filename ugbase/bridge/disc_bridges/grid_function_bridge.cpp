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
#include "bridge/util_domain_algebra_dependent.h"

// lib_disc includes
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/function_spaces/approximation_space.h"
#include "lib_disc/function_spaces/grid_function_util.h"
#include "lib_disc/function_spaces/grid_function_user_data.h"
#include "lib_disc/function_spaces/dof_position_util.h"
#include "lib_disc/function_spaces/grid_function_global_user_data.h"
#include "lib_disc/function_spaces/grid_function_user_data_explicit.h"
#include "lib_disc/function_spaces/grid_function_coordinate_util.h"
#include "lib_disc/function_spaces/metric_spaces.h"
using namespace std;

namespace ug{
namespace bridge{
namespace GridFunction{

/**
 * \defgroup gridfnct_bridge Grid Function Bridge
 * \ingroup disc_bridge
 * \{
 */

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

/**
 * Function called for the registration of Domain and Algebra dependent parts.
 * All Functions and Classes depending on both Domain and Algebra
 * are to be placed here when registering. The method is called for all
 * available Domain and Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain, typename TAlgebra, typename TRegistry>
static void DomainAlgebra(TRegistry& reg, string grp)
{
	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

//	typedef
	static const int dim = TDomain::dim;
	typedef typename TAlgebra::vector_type vector_type;
	typedef ApproximationSpace<TDomain> approximation_space_type;
	typedef ug::GridFunction<TDomain, TAlgebra> TFct;

//	group string
	grp.append("/ApproximationSpace");

//	GridFunction
	{
		string name = string("GridFunction").append(suffix);
		reg.template add_class_<TFct, vector_type>(name, grp)
			.template add_constructor<void (*)(SmartPtr<approximation_space_type>)>("ApproximationSpace")
			.template add_constructor<void (*)(SmartPtr<approximation_space_type>, int)>("ApproximationSpace#Level")
			.add_method("assign", static_cast<void (TFct::*)(const vector_type&)>(&TFct::assign),
						"Success", "Vector")
			.add_method("clone", &TFct::clone)
			.add_method("set_consistent_storage_type", &TFct::SetConsistentStorageType)
			.add_method("grid_level", &TFct::grid_level)
			.add_method("num_dofs", static_cast<size_t (TFct::*)() const>(&TFct::num_dofs))
			.add_method("approx_space", static_cast<SmartPtr<approximation_space_type> (TFct::*)()>(&TFct::approx_space))
			.add_method("redistribution_enabled", &TFct::redistribution_enabled)
			.add_method("enable_redistribution", &TFct::enable_redistribution)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GridFunction", tag);
	}

//	ExplicitGridFunctionValue
	{
		string name = string("ExplicitGridFunctionValue").append(suffix);
		typedef ExplicitGridFunctionValue<TFct> T;
		typedef CplUserData<number, dim> TBase;
		reg.template add_class_<T, TBase>(name, grp)
		   .template add_constructor<void (*)(SmartPtr<TFct>, const char*)>("GridFunction#Component")
		   .set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ExplicitGridFunctionValue", tag);
	}

//	ExplicitGridFunctionVector
	{
		string name = string("ExplicitGridFunctionVector").append(suffix);
		typedef ExplicitGridFunctionVector<TFct> T;
		typedef CplUserData<MathVector<dim>, dim> TBase;
		reg.template add_class_<T, TBase>(name, grp)
		   .template add_constructor<void (*)(SmartPtr<TFct>, const char*)>("GridFunction#Components")
		   .set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ExplicitGridFunctionVector", tag);
	}

//	ExplicitGridFunctionGradient
	{
		string name = string("ExplicitGridFunctionGradient").append(suffix);
		typedef ExplicitGridFunctionGradient<TFct> T;
		typedef CplUserData<MathVector<dim>, dim> TBase;
		reg.template add_class_<T, TBase>(name, grp)
		   .template add_constructor<void (*)(SmartPtr<TFct>, const char*)>("GridFunction#Component")
            .add_method("add_subset_coeff", &T::add_subset_coeff)
			.add_method("get_subset_coeff", &T::get_subset_coeff)
		   .set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ExplicitGridFunctionGradient", tag);
	}

//	GridFunctionNumberData
	{
		string name = string("GridFunctionNumberData").append(suffix);
		typedef GridFunctionNumberData<TFct> T;
		typedef CplUserData<number, dim> TBase;
		reg.template add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr<TFct>, const char*)>("GridFunction#Component")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GridFunctionNumberData", tag);
	}

//	GridFunctionVectorData
	{
		string name = string("GridFunctionVectorData").append(suffix);
		typedef GridFunctionVectorData<TFct> T;
		typedef CplUserData<MathVector<dim>, dim> TBase;
		reg.template add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr<TFct>, const char*)>("GridFunction#Components")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GridFunctionVectorData", tag);
	}

//	GridFunctionGradientData
	{
		string name = string("GridFunctionGradientData").append(suffix);
		typedef GridFunctionGradientData<TFct> T;
		typedef CplUserData<MathVector<dim>, dim> TBase;
		reg.template add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr<TFct>, const char*)>("GridFunction#Component")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GridFunctionGradientData", tag);
	}

//	GridFunctionGradientComponentData
	{
		string name = string("GridFunctionGradientComponentData").append(suffix);
		typedef GridFunctionGradientComponentData<TFct> T;
		typedef CplUserData<number, dim> TBase;
		reg.template add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr<TFct>, const char*, size_t)>("GridFunction#Components")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GridFunctionGradientComponentData", tag);
	}
//	GlobalGridFunctionNumberData
	{
		string name = string("GlobalGridFunctionNumberData").append(suffix);
		typedef GlobalGridFunctionNumberData<TFct> T;
		typedef CplUserData<number, dim> TBase;
		reg.template add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr<TFct>, const char*)>("GridFunction#Component")
			.add_method("evaluate", static_cast<number (T::*)(const MathVector<dim>&) const>(&T::evaluate))
			.add_method("evaluate_global", static_cast<number (T::*)(std::vector<number>)>(&T::evaluate_global))
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GlobalGridFunctionNumberData", tag);
	}

//	GlobalGridFunctionNumberData for lower-dim elem geometries (here: edges)
	if (dim > EDGE)
	{
		string name = string("GlobalEdgeGridFunctionNumberData").append(suffix);
		typedef GlobalGridFunctionNumberData<TFct, 1> T;
		typedef CplUserData<number, dim> TBase;
		reg.template add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr<TFct>, const char*)>("GridFunction#Component")
			.add_method("evaluate", static_cast<number (T::*)(const MathVector<dim>&) const>(&T::evaluate))
			.add_method("evaluate_global", static_cast<number (T::*)(std::vector<number>)>(&T::evaluate_global))
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GlobalEdgeGridFunctionNumberData", tag);
	}


//	GlobalGridFunctionGradientData
	{
		string name = string("GlobalGridFunctionGradientData").append(suffix);
		typedef GlobalGridFunctionGradientData<TFct> T;
		typedef CplUserData<MathVector<dim>, dim> TBase;
		reg.template add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr<TFct>, const char*)>("GridFunction#Component")
			.add_method("evaluate_global", static_cast<std::vector<number> (T::*)(std::vector<number>)>(&T::evaluate_global))
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GlobalGridFunctionGradientData", tag);
	}

	// IGridFunctionSpace (abstract base class)
	{
			typedef IGridFunctionSpace<TFct> T;
			string name = string("IGridFunctionSpace").append(suffix);
			reg.template add_class_<T>(name, grp)
			   .add_method("config_string", &T::config_string);
			reg.add_class_to_group(name, "IGridFunctionSpace", tag);
	}

	// IComponentSpace (abstract base class for scalar component)
	{
		typedef IComponentSpace<TFct> T;
		typedef IGridFunctionSpace<TFct> TBase;

		string name = string("IComponentSpace").append(suffix);
		reg.template add_class_<T, TBase>(name, grp)
		   .add_method("norm", static_cast<number (T::*)(TFct&) > (&T::norm))
		   .add_method("distance", static_cast<number (T::*)(TFct&, TFct&) > (&T::distance));

		reg.add_class_to_group(name, "IComponentSpace", tag);
	}

	// GridFunctionComponentSpace
	{
		typedef GridFunctionComponentSpace<TFct> T;
		typedef IComponentSpace<TFct> TBase;

		string name = string("GridFunctionComponentSpace").append(suffix);
		reg.template add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(const char *) >("function names")
			.template add_constructor<void (*)(const char *, const char *) >("function names, subset names")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GridFunctionComponentSpace", tag);
	}

	// L2ComponentSpace
	{
		typedef L2ComponentSpace<TFct> T;
		typedef IComponentSpace<TFct> TBase;
		typedef typename L2Integrand<TFct>::weight_type TWeight;

		string name = string("L2ComponentSpace").append(suffix);
		reg.template add_class_<T, TBase>(name, grp)
		   .template add_constructor<void (*)(const char *) >("fctNames")
		   .template add_constructor<void (*)(const char *, int) >("fctNames, order")
		   .template add_constructor<void (*)(const char *, int, double) >("fctNames, order, weight")
		   .template add_constructor<void (*)(const char *, int, double, const char *) >("fctNames, order, weight, ssNames")
		   .template add_constructor<void (*)(const char *, int, ConstSmartPtr<TWeight>) >("fctNames, order, weight")
		   .template add_constructor<void (*)(const char *, int, ConstSmartPtr<TWeight>, const char *) >("fctNames, order, weight, ssNames")
		 //  .template add_constructor<void (*)(const char *, const char *, int, ConstSmartPtr<TWeight>) >("fctNames, ssNames, order, weight")

		   .set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "L2ComponentSpace", tag);
	}

	// L2QuotientSpace (= L2ComponentSpace factoring out constants)
	{
		typedef L2QuotientSpace<TFct> T;
		typedef IComponentSpace<TFct> TBase;
		typedef typename L2Integrand<TFct>::weight_type TWeight;

		string name = string("L2QuotientSpace").append(suffix);
		reg.template add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(const char *) >("fctNames")
			.template add_constructor<void (*)(const char *, int) >("fctNames, order")
			.template add_constructor<void (*)(const char *, int, double) >("fctNames, order, weight")
			.template add_constructor<void (*)(const char *, int, double, const char *) >("fctNames, order, weight, ssNames")
			.template add_constructor<void (*)(const char *, int, ConstSmartPtr<TWeight>) >("fctNames, order, weight")
			.template add_constructor<void (*)(const char *, int, ConstSmartPtr<TWeight>, const char *) >("fctNames, order, weight, ssNames")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "L2QuotientSpace", tag);
	}

	// H1SemiComponentSpace
	{
		typedef H1SemiComponentSpace<TFct> T;
		typedef IComponentSpace<TFct> TBase;
		typedef typename H1SemiIntegrand<TFct>::weight_type TWeight;

		string name = string("H1SemiComponentSpace").append(suffix);
		reg.template add_class_<T, TBase>(name, grp)
		   .template add_constructor<void (*)(const char *) >("fctNames")
		   .template add_constructor<void (*)(const char *, int) >("fctNames, order")
		   .template add_constructor<void (*)(const char *, int, number) >("fctNames, order, weight")
		   .template add_constructor<void (*)(const char *, int, number, const char *) >("fctNames, order, weight, ssNames")
		   .template add_constructor<void (*)(const char *, int, ConstSmartPtr<TWeight>) >("fctNames, order, weight")
		   .template add_constructor<void (*)(const char *, int, const char *, ConstSmartPtr<TWeight>) >("fctNames, order, weight, ssNames")
		   .add_method("set_weight", &T::set_weight)
		   .add_method("get_weight", &T::get_weight)
		   .add_method("norm", static_cast<number (T::*)(TFct&) > (&T::norm))
		   .add_method("distance", static_cast<number (T::*)(TFct&, TFct&) > (&T::distance))
		   .set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "H1SemiComponentSpace", tag);
	}
/*
	// KineticEnergyComponentSpace
			{
				typedef KineticEnergyComponentSpace<TFct> T;
				typedef IComponentSpace<TFct> TBase;
				typedef typename T::weight_type TWeight;
				typedef typename T::velocity_type TVelocity;

				string name = string("KineticEnergyComponentSpace").append(suffix);
				reg.add_class_<T, TBase>(name, grp)
				   .template add_constructor<void (*)(SmartPtr<TVelocity>, const char *) >("fctNames")
				   .template add_constructor<void (*)(SmartPtr<TVelocity>,const char *, int) >("fctNames, order")
				   .template add_constructor<void (*)(SmartPtr<TVelocity>,const char *, int, number) >("fctNames, order, weight")
				   .template add_constructor<void (*)(SmartPtr<TVelocity>,const char *, int, number, const char *) >("fctNames, order, weight, ssNames")
				   .template add_constructor<void (*)(SmartPtr<TVelocity>,const char *, int, ConstSmartPtr<TWeight>) >("fctNames, order, weight")
				   .template add_constructor<void (*)(SmartPtr<TVelocity>,const char *, int, const char *, ConstSmartPtr<TWeight>) >("fctNames, order, weight, ssNames")
				   .add_method("set_weight", &T::set_weight)
				   .add_method("get_weight", &T::get_weight)
				   .set_construct_as_smart_pointer(true);
				reg.add_class_to_group(name, "KineticEnergyComponentSpace", tag);
			}
*/
	// H1EnergyComponentSpace
		{
			typedef H1EnergyComponentSpace<TFct> T;
			typedef IComponentSpace<TFct> TBase;
			typedef typename H1EnergyIntegrand<TFct>::weight_type TWeight;

			string name = string("VelEnergyComponentSpace").append(suffix);
			reg.template add_class_<T, TBase>(name, grp)
			   .template add_constructor<void (*)(const char *) >("fctNames")
			   .template add_constructor<void (*)(const char *, int) >("fctNames, order")
			   .template add_constructor<void (*)(const char *, int, number) >("fctNames, order, weight")
			   .template add_constructor<void (*)(const char *, int, number, const char *) >("fctNames, order, weight, ssNames")
			   .template add_constructor<void (*)(const char *, int, ConstSmartPtr<TWeight>) >("fctNames, order, weight")
			   //.template add_constructor<void (*)(const char *, int, const char *, ConstSmartPtr<TWeight>) >("fctNames, order, ssNames, weight")
			   .template add_constructor<void (*)(const char *, int, ConstSmartPtr<TWeight>, const char*) >("fctNames, order, weight, ssNames")
			   .add_method("set_weight", &T::set_weight)
			   .add_method("get_weight", &T::get_weight)
			   .add_method("norm", static_cast<number (T::*)(TFct&) > (&T::norm))
			   .add_method("distance", static_cast<number (T::*)(TFct&, TFct&) > (&T::distance))
			   .add_method("set_velocity", &T::set_velocity)
			   .set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "VelEnergyComponentSpace", tag);
		}

	// H1ComponentSpace
	{
		typedef H1ComponentSpace<TFct> T;
		typedef IComponentSpace<TFct> TBase;

		string name = string("H1ComponentSpace").append(suffix);
		reg.template add_class_<T, TBase>(name, grp)
		   .template add_constructor<void (*)(const char *) >("fctNames")
		   .template add_constructor<void (*)(const char *, int) >("fctNames, order")
		   .template add_constructor<void (*)(const char *, const char*, int) >("fctNames, subsetNames, order")
		   //.template add_constructor<void (*)(const char *, int, number) >("fctNames, order, scale")
		   .set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "H1ComponentSpace", tag);
	}

	// TimeDependentSpace
	/*{
		typedef TimeDependentSpace<TFct> T;
		typedef IGridFunctionSpace<TFct> TBase;

		typedef IComponentSpace<TFct> TCompSpace;

		string name = string("TimeDependentSpace").append(suffix);
		reg.template add_class_<T, TBase>(name, grp)
				.template add_constructor<void (*)(SmartPtr<TCompSpace>, number) >("component space")
				.add_method("update_time_data", &T::update_time_data)
				.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "TimeDependentSpace", tag);
	}
*/

	// CompositeSpace
	{
			typedef CompositeSpace<TFct> T;
			typedef IComponentSpace<TFct> TCompSpace;
			typedef IGridFunctionSpace<TFct> TBase;

			string name = string("CompositeSpace").append(suffix);
			reg.template add_class_<T, TBase>(name, grp)
			   .template add_constructor<void (*)() >("")
			   .add_method("add", static_cast<void (T::*)(SmartPtr<TCompSpace>) > (&T::add))
			   .add_method("add", static_cast<void (T::*)(SmartPtr<TCompSpace>, number) > (&T::add))
			   .add_method("update_time_data", &T::update_time_data)
			   .add_method("is_time_dependent", &T::is_time_dependent)
			   .set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "CompositeSpace", tag);
	}

//	AverageFunctionDifference
	{
		string name = string("AverageFunctionDifference");
		typedef ug::GridFunction<TDomain, TAlgebra> grid_function;
		typedef SmartPtr< grid_function > function_pointer;
		reg.add_function(name, static_cast<number (*)(function_pointer, std::string, std::string, std::string)>(&AverageFunctionDifference<TDomain, TAlgebra>), grp);
	}

//	CheckDoFPositions
	{
		reg.add_function("CheckDoFPositions", static_cast<bool (*)(const TFct&)>(CheckDoFPositions<TFct>), grp);
	}

//	ScaleGF
	{
		reg.add_function("ScaleGF", ScaleGF<TFct>, grp, "",
			"scaled output vector # input vector # vector of scaling factors for each function",
			"Scales the input vector using the given scaling factors for each function and writes "
			"the result to the output vector");
	}

//	AverageFunctionDifference
	{
		typedef ug::GridFunction<TDomain, TAlgebra> GF;
		reg.add_function("AdjustMeanValue", static_cast<void (*)(SmartPtr<GF>, const std::vector<std::string>&, number)>(&AdjustMeanValue<GF>), grp);
		reg.add_function("AdjustMeanValue", static_cast<void (*)(SmartPtr<GF>, const std::vector<std::string>&)>(&AdjustMeanValue<GF>), grp);
		reg.add_function("AdjustMeanValue", static_cast<void (*)(SmartPtr<GF>, const std::string&, number)>(&AdjustMeanValue<GF>), grp);
		reg.add_function("AdjustMeanValue", static_cast<void (*)(SmartPtr<GF>, const std::string&)>(&AdjustMeanValue<GF>), grp);
	}
	
//	SumGFValuesAt
	{
		typedef ug::GridFunction<TDomain, TAlgebra> GF;
		reg.add_function ("SumGFValuesAtVertices", static_cast<number (*) (GF*, const char *)> (&SumGFValuesAt<GF,Vertex>), grp);
		reg.add_function ("SumGFValuesAtVertices", static_cast<number (*) (GF*, const char *, const char *)> (&SumGFValuesAt<GF,Vertex>), grp);
	}
	
//	CheckGFforNaN
	{
		typedef ug::GridFunction<TDomain, TAlgebra> GF;
		reg.add_function ("CheckGFValuesAtVertices", static_cast<bool (*) (const GF*, const char *)> (&CheckGFforNaN<GF,Vertex>), grp);
		reg.add_function ("CheckGFValuesAtEdges", static_cast<bool (*) (const GF*, const char *)> (&CheckGFforNaN<GF,Edge>), grp);
		reg.add_function ("CheckGFValuesAtFaces", static_cast<bool (*) (const GF*, const char *)> (&CheckGFforNaN<GF,Face>), grp);
		reg.add_function ("CheckGFValuesAtVolumes", static_cast<bool (*) (const GF*, const char *)> (&CheckGFforNaN<GF,Volume>), grp);
	}

//	CheckGFValuesWithinBounds
	{
		typedef ug::GridFunction<TDomain, TAlgebra> GF;
		reg.add_function("CheckGFValuesWithinBounds", static_cast<bool (*) (ConstSmartPtr<GF>, size_t, number, number)> (&CheckGFValuesWithinBounds<GF>), grp);
		reg.add_function("CheckGFValuesWithinBounds", static_cast<bool (*) (ConstSmartPtr<GF>, const char*, number, number)> (&CheckGFValuesWithinBounds<GF>), grp);
	}

//	Move Domain by GridFunction
	{
		typedef ug::GridFunction<TDomain, TAlgebra> GF;
		reg.add_function (
			"AddFunctionValuesToGridCoordinatesP1", static_cast<void (*) (SmartPtr<GF>, const char*, size_t)>
				(&AddFunctionValuesToGridCoordinatesP1<GF>), grp);
		reg.add_function (
			"AddFunctionValuesToGridCoordinatesP1", static_cast<void (*) (SmartPtr<GF>, const char*, size_t, number)>
				(&AddFunctionValuesToGridCoordinatesP1<GF>), grp);
	}
}

/**
 * Function called for the registration of Domain dependent parts.
 * All Functions and Classes depending on the Domain
 * are to be placed here when registering. The method is called for all
 * available Domain types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain, typename TRegistry>
static void Domain(TRegistry& reg, string grp)
{
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();

//	group string
	grp.append("/ApproximationSpace");

//  ApproximationSpace
	{
		typedef ApproximationSpace<TDomain> T;
		typedef IApproximationSpace TBase;
		string name = string("ApproximationSpace").append(suffix);
		reg.template add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr<TDomain>)>("Domain")
			.template add_constructor<void (*)(SmartPtr<TDomain>, const AlgebraType&)>("Domain#AlgebraType")
			.add_method("domain", static_cast<SmartPtr<TDomain> (T::*)()>(&T::domain))
			.add_method("surface_view", static_cast<ConstSmartPtr<SurfaceView> (T::*)() const>(&T::surface_view))
			.add_method("get_dim", &T::get_dim)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ApproximationSpace", tag);
	}
}

/**
 * Function called for the registration of Dimension dependent parts.
 * All Functions and Classes depending on the Dimension
 * are to be placed here when registering. The method is called for all
 * available Dimension types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <int dim>
static void Dimension(Registry& reg, string grp)
{
	string suffix = GetDimensionSuffix<dim>();
	string tag = GetDimensionTag<dim>();

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
template <typename TAlgebra, typename TRegistry>
static void Algebra(TRegistry& reg, string grp)
{
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();

}

/**
 * Function called for the registration of Domain and Algebra independent parts.
 * All Functions and Classes not depending on Domain and Algebra
 * are to be placed here when registering.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TRegistry>
static void Common(TRegistry& reg, string grp)
{
//	GridLevel
	reg.template add_class_<GridLevel>("GridLevel", grp)
		.add_constructor()
		.template add_constructor<void (*)(int)>("Level")
		.template add_constructor<void (*)(int, std::string)>("Level#Type")
		.set_construct_as_smart_pointer(true);

//	LFEID
	{
		typedef LFEID T;
		reg.template add_class_<T>("LFEID", grp)
			.add_method("order", &T::order)
			.add_method("dim", &T::dim);
	}

//	DoFDistributionInfoProvider
	{
	typedef DoFDistributionInfoProvider T;
	reg.template add_class_<T>("DoFDistributionInfoProvider", grp)
		.add_method("print_local_dof_statistic", static_cast<void (T::*)(int) const>(&T::print_local_dof_statistic))
		.add_method("print_local_dof_statistic", static_cast<void (T::*)() const>(&T::print_local_dof_statistic))
		.add_method("num_fct", static_cast<size_t (T::*)() const>(&T::num_fct))
		.add_method("name", &T::name)
		.add_method("names", &T::names)
		.add_method("dim", &T::dim)
		.add_method("lfeid", &T::lfeid);
	}

//	IApproximationSpace
	{
	typedef IApproximationSpace T;
	typedef DoFDistributionInfoProvider TBase;
	reg.template add_class_<T, TBase>("IApproximationSpace", grp)
		.add_method("print_statistic", static_cast<void (T::*)(std::string) const>(&T::print_statistic))
		.add_method("print_statistic", static_cast<void (T::*)() const>(&T::print_statistic))
		.add_method("print_layout_statistic", static_cast<void (T::*)() const>(&T::print_layout_statistic))
		.add_method("might_contain_ghosts", static_cast<bool (T::*)() const>(&T::might_contain_ghosts))
		.add_method("num_levels", &T::num_levels)
		.add_method("init_levels", &T::init_levels)
		.add_method("init_surfaces", &T::init_surfaces)
		.add_method("init_top_surface", &T::init_top_surface)

		.add_method("clear", &T::clear)
		.add_method("add_fct", static_cast<void (T::*)(const char*, const char*, int, const char*)>(&T::add),
					"", "Name#Type|selection|value=[\"Lagrange\",\"DG\"]#Order#Subsets", "Adds a function to the Function Pattern",
					"currently no help available")
		.add_method("add_fct", static_cast<void (T::*)(const char*, const char*, int)>(&T::add),
					"", "Name#Type|selection|value=[\"Lagrange\",\"DG\"]#Order", "Adds a function to the Function Pattern",
					"currently no help available")
		.add_method("add_fct", static_cast<void (T::*)(const char*, const char*)>(&T::add),
					"", "Name#Type|selection|value=[\"crouzeix-raviart\",\"piecewise-constant\"] ", "Adds a function to the Function Pattern",
					"currently no help available")
		.add_method("add_fct", static_cast<void (T::*)(const std::vector<std::string>&, const char*, int, const std::vector<std::string>&)>(&T::add),
					"", "Name#Type|selection|value=[\"Lagrange\",\"DG\"]#Order#Subsets", "Adds a function to the Function Pattern",
					"currently no help available")
		.add_method("add_fct", static_cast<void (T::*)(const std::vector<std::string>&, const char*, int)>(&T::add),
					"", "Name#Type|selection|value=[\"Lagrange\",\"DG\"]#Order", "Adds a function to the Function Pattern",
					"currently no help available")
		.add_method("add_fct", static_cast<void (T::*)(const std::vector<std::string>&, const char*)>(&T::add),
					"", "Name#Type|selection|value=[\"crouzeix-raviart\",\"piecewise-constant\"]", "Adds a function to the Function Pattern",
					"currently no help available");

	}
}

}; // end Functionality

// end group gridfnct_bridge
/// \}

}// namespace GridFunction

/// \addtogroup gridfnct_bridge
template <typename TRegistry=Registry>
void RegisterBridge_GridFunction_(TRegistry& reg, string grp)
{
	grp.append("/Discretization");
	typedef GridFunction::Functionality TFunctionality;

	try{
		RegisterCommon<TFunctionality>(reg,grp);
		//RegisterDimensionDependent<TFunctionality, TRegistry >(reg,grp);
		RegisterDomainDependent<TFunctionality, TRegistry>(reg,grp);
		//RegisterAlgebraDependent<TFunctionality, TRegistry>(reg,grp);
		RegisterDomainAlgebraDependent<TFunctionality, TRegistry>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}//	end of namespace bridge

UG_REGISTRY_DEFINE(RegisterBridge_GridFunction);
}//	end of namespace ug
