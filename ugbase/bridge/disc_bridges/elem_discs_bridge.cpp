/*
 * elem_discs_bridge.cpp
 *
 *  Created on: 20.05.2011
 *      Author: andreasvogel
 */

// extern headers
#include <iostream>
#include <sstream>
#include <string>

// include bridge
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_domain_dependent.h"

// lib_disc includes
#include "lib_disc/domain.h"

#include "lib_disc/spatial_disc/disc_util/conv_shape_interface.h"
#include "lib_disc/spatial_disc/disc_util/conv_shape.h"

#include "lib_disc/spatial_disc/elem_disc/neumann_boundary/neumann_boundary_base.h"
#include "lib_disc/spatial_disc/elem_disc/neumann_boundary/fv1/neumann_boundary_fv1.h"
#include "lib_disc/spatial_disc/elem_disc/neumann_boundary/fv/neumann_boundary_fv.h"
#include "lib_disc/spatial_disc/elem_disc/neumann_boundary/fe/neumann_boundary_fe.h"
#include "lib_disc/spatial_disc/elem_disc/inner_boundary/inner_boundary.h"

using namespace std;

namespace ug{
namespace bridge{
namespace ElemDiscs{

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

/**
 * Function called for the registration of Domain dependent parts.
 * All Functions and Classes depending on the Domain
 * are to be placed here when registering. The method is called for all
 * available Domain types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain>
static void Domain(Registry& reg, string grp)
{
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();
	static const int dim = TDomain::dim;

	string approxGrp = grp; approxGrp.append("/ApproximationSpace");
	string elemGrp = grp; elemGrp.append("/SpatialDisc/ElemDisc");

//	DomainElemDisc base class
	{
		typedef IElemDisc<TDomain> T;
		string name = string("IElemDisc").append(suffix);
		reg.add_class_<T>(name, elemGrp)
			.add_method("set_stationary", static_cast<void (T::*)()>(&T::set_stationary));
		reg.add_class_to_group(name, "IElemDisc", tag);
	}

//	Neumann Boundary Base
	{
		typedef NeumannBoundaryBase<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("NeumannBoundaryBase").append(suffix);
		reg.add_class_<T, TBase >(name, elemGrp)
			.add_method("add", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >, const char*, const char*)>(&T::add))
			.add_method("add", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim, bool> >, const char*, const char*)>(&T::add))
			.add_method("add", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim> >, const char*, const char*)>(&T::add))
			.add_method("add", static_cast<void (T::*)(number, const char*, const char*)>(&T::add))
#ifdef UG_FOR_LUA
			.add_method("add", static_cast<void (T::*)(const char*, const char*, const char*)>(&T::add))
#endif
			.add_method("add", static_cast<void (T::*)(const vector<number>&, const char*, const char*)>(&T::add));
		reg.add_class_to_group(name, "NeumannBoundaryBase", tag);
	}

//	Neumann Boundary FV1
	{
		typedef NeumannBoundaryFV1<TDomain> T;
		typedef NeumannBoundaryBase<TDomain> TBase;
		string name = string("NeumannBoundaryFV1").append(suffix);
		reg.add_class_<T, TBase >(name, elemGrp)
			.template add_constructor<void (*)(const char*)>("Function")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NeumannBoundaryFV1", tag);
	}

//	Neumann Boundary FV
	{
		typedef NeumannBoundaryFV<TDomain> T;
		typedef NeumannBoundaryBase<TDomain> TBase;
		string name = string("NeumannBoundaryFV").append(suffix);
		reg.add_class_<T, TBase >(name, elemGrp)
			.template add_constructor<void (*)(const char*)>("Function")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NeumannBoundaryFV", tag);
	}

//	Neumann Boundary FE
	{
		typedef NeumannBoundaryFE<TDomain> T;
		typedef NeumannBoundaryBase<TDomain> TBase;
		string name = string("NeumannBoundaryFE").append(suffix);
		reg.add_class_<T, TBase >(name, elemGrp)
			.template add_constructor<void (*)(const char*)>("Function")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NeumannBoundaryFE", tag);
	}

//	Inner Boundaries
	{
		typedef FV1InnerBoundaryElemDisc<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("FV1InnerBoundary").append(suffix);
		reg.add_class_<T, TBase >(name, elemGrp);
		reg.add_class_to_group(name, "FV1InnerBoundary", tag);
	}


/////////////////////////////////////////////////////////////////////////////
// Convection Shapes
/////////////////////////////////////////////////////////////////////////////

	string upGrp = grp; upGrp.append("/SpatialDisc/Upwind");

//	IConvectionShapes
	{
		typedef IConvectionShapes<dim> T;
		string name = string("IConvectionShapes").append(suffix);
		reg.add_class_<T>(name, upGrp);
		reg.add_class_to_group(name, "IConvectionShapes", tag);
	}

//	ConvectionShapesNoUpwind
	{
		typedef ConvectionShapesNoUpwind<dim> T;
		typedef IConvectionShapes<dim> TBase;
		string name = string("NoUpwind").append(suffix);
		reg.add_class_<T, TBase>(name, upGrp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NoUpwind", tag);
	}

//	ConvectionShapesFullUpwind
	{
		typedef ConvectionShapesFullUpwind<dim> T;
		typedef IConvectionShapes<dim> TBase;
		string name = string("FullUpwind").append(suffix);
		reg.add_class_<T, TBase>(name, upGrp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "FullUpwind", tag);
	}

//	ConvectionShapesWeightedUpwind
	{
		typedef ConvectionShapesWeightedUpwind<dim> T;
		typedef IConvectionShapes<dim> TBase;
		string name = string("WeightedUpwind").append(suffix);
		reg.add_class_<T, TBase>(name, upGrp)
			.add_method("set_weight", &T::set_weight)
			.add_constructor()
			.template add_constructor<void (*)(number)>("weight")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "WeightedUpwind", tag);
	}

//	ConvectionShapesPartialUpwind
	{
		typedef ConvectionShapesPartialUpwind<dim> T;
		typedef IConvectionShapes<dim> TBase;
		string name = string("PartialUpwind").append(suffix);
		reg.add_class_<T, TBase>(name, upGrp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "PartialUpwind", tag);
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
//	get group string
	grp.append("/Discretization");
}
};

}// namespace ElemDiscs

void RegisterBridge_ElemDiscs(Registry& reg, string grp)
{
	typedef ElemDiscs::Functionality Functionality;

	try{
		RegisterCommon<Functionality>(reg,grp);
		RegisterDomainDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}//	end of namespace bridge
}//	end of namespace ug
