/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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
#include "lib_disc/spatial_disc/elem_disc/dirac_source/lagrange_dirac_source.h"

using namespace std;

namespace ug{
namespace bridge{
namespace ElemDiscs{

/**
 * \defgroup elemdisc_bridge Element Discretization Bridge
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

//	ElemDiscModifier base class
	{
		typedef IElemDiscModifier<TDomain> T;
		string name = string("IElemDiscModifier").append(suffix);
		reg.add_class_<T>(name, elemGrp);
		reg.add_class_to_group(name, "IElemDiscModifier", tag);
	}

	//	IElemError base class
	{
			typedef IElemError<TDomain> T;
			string name = string("IElemError").append(suffix);
			reg.add_class_<T>(name, elemGrp)
	 			.add_method("set_stationary", static_cast<void (T::*)()>(&T::set_stationary))
	 		//	.add_method("add_elem_modifier", &T::add_elem_modifier, "", "")
	 			.add_method("set_error_estimator", static_cast<void (T::*)(SmartPtr<IErrEstData<TDomain> >)>(&T::set_error_estimator));
	 		reg.add_class_to_group(name, "IElemError", tag);
	}


//	IElemDisc base class
	{
		typedef IElemDisc<TDomain> T;
		typedef IElemError<TDomain> TBase;
		string name = string("IElemDisc").append(suffix);
		reg.add_class_<T, TBase>(name, elemGrp)
 		//	.add_method("set_stationary", static_cast<void (T::*)()>(&T::set_stationary))
 			.add_method("add_elem_modifier", &T::add_elem_modifier, "", "");
 		//	.add_method("set_error_estimator", static_cast<void (T::*)(SmartPtr<IErrEstData<TDomain> >)>(&T::set_error_estimator));
 		reg.add_class_to_group(name, "IElemDisc", tag);
	}

//	Neumann Boundary Base
	{
		typedef NeumannBoundaryBase<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("NeumannBoundaryBase").append(suffix);
		reg.add_class_<T, TBase >(name, elemGrp)
			.add_method("add", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >, const char*, const char*)>(&T::add))
			.add_method("add", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >, const std::vector<std::string>&, const std::vector<std::string>&)>(&T::add))
			.add_method("add", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim, bool> >, const char*, const char*)>(&T::add))
			.add_method("add", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim, bool> >, const std::vector<std::string>&, const std::vector<std::string>&)>(&T::add))
			.add_method("add", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim> >, const char*, const char*)>(&T::add))
			.add_method("add", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim> >, const std::vector<std::string>&, const std::vector<std::string>&)>(&T::add))
			.add_method("add", static_cast<void (T::*)(number, const char*, const char*)>(&T::add))
			.add_method("add", static_cast<void (T::*)(number, const std::vector<std::string>&, const std::vector<std::string>&)>(&T::add))
#ifdef UG_FOR_LUA
			.add_method("add", static_cast<void (T::*)(const char*, const char*, const char*)>(&T::add))
			.add_method("add", static_cast<void (T::*)(const char*, const std::vector<std::string>&, const std::vector<std::string>&)>(&T::add))
			.add_method("add", static_cast<void (T::*)(LuaFunctionHandle, const char*, const char*)>(&T::add))
			.add_method("add", static_cast<void (T::*)(LuaFunctionHandle, const std::vector<std::string>&, const std::vector<std::string>&)>(&T::add))
#endif
			.add_method("add", static_cast<void (T::*)(const vector<number>&, const char*, const char*)>(&T::add))
			.add_method("add", static_cast<void (T::*)(const vector<number>&, const std::vector<std::string>&, const std::vector<std::string>&)>(&T::add));
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

#if 0
//	Inner Boundaries
	{
		typedef FV1InnerBoundaryElemDisc<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("FV1InnerBoundary").append(suffix);
		reg.add_class_<T, TBase >(name, elemGrp)
			.add_method("set_flux_scale", static_cast<void (T::*)(number)>(&T::set_flux_scale),
				"", "scale", "Set scale to scale (all) fluxes with.")
			.add_method("set_flux_scale", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_flux_scale),
				"", "scale", "Set scale to scale (all) fluxes with.")
#ifdef UG_FOR_LUA
			.add_method("set_flux_scale", static_cast<void (T::*)(const char*)>(&T::set_flux_scale),
				"", "scale", "Set scale to scale (all) fluxes with.")
#endif
			;
		reg.add_class_to_group(name, "FV1InnerBoundary", tag);
	}
#endif

	//	DiracSourceDisc
		{
			typedef DiracSourceDisc<TDomain> T;
			typedef IElemDisc<TDomain> TBase;
			string name = string("DiracSourceDisc").append(suffix);
			reg.add_class_<T, TBase >(name, elemGrp)
				.template add_constructor<void (*)(const char*, const char*)>("Function")
				.add_method("add_source", static_cast<void (T::*)(number, MathVector<dim> &)>(&T::add_source),
					"", "scale", "Set scale to scale (all) fluxes with.")
				.add_method("add_source", static_cast<void (T::*)(SmartPtr<UserData<number, dim> >, MathVector<dim> &)>(&T::add_source),
					"", "scale", "Set scale to scale (all) fluxes with.")
#ifdef UG_FOR_LUA
				.add_method("add_source", static_cast<void (T::*)(const char*, MathVector<dim> &)>(&T::add_source),
					"", "scale", "Set scale to scale (all) fluxes with.")
#endif
				.add_method("add_transport_sink", static_cast<void (T::*)(number)>(&T::add_transport_sink),
							"", "scale", "Set scale to scale (all) fluxes with.")
				.add_method("add_transport_sink", static_cast<void (T::*)(SmartPtr<UserData<number, dim> >)>(&T::add_transport_sink),
							"", "scale", "Set scale to scale (all) fluxes with.")
#ifdef UG_FOR_LUA
				.add_method("add_transport_sink", static_cast<void (T::*)(const char*)>(&T::add_transport_sink),
							"", "scale", "Set scale to scale (all) fluxes with.")
#endif


					.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "DiracSourceDisc", tag);
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

//	ConvectionShapesSkewedUpwind
	{
		typedef ConvectionShapesSkewedUpwind<dim> T;
		typedef IConvectionShapes<dim> TBase;
		string name = string("SkewedUpwind").append(suffix);
		reg.add_class_<T, TBase>(name, upGrp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "SkewedUpwind", tag);
	}

//	ConvectionShapesLinearProfileSkewedUpwind
	{
		typedef ConvectionShapesLinearProfileSkewedUpwind<dim> T;
		typedef IConvectionShapes<dim> TBase;
		string name = string("LinearProfileSkewedUpwind").append(suffix);
		reg.add_class_<T, TBase>(name, upGrp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "LinearProfileSkewedUpwind", tag);
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

// end group elemdisc_bridge
/// \}

}// namespace ElemDiscs

/// \addtogroup elemdisc_bridge
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
