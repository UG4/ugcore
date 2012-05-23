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
#include "../bridge.h"
#include "registry/registry.h"

// lib_disc includes
#include "lib_disc/domain.h"

#include "lib_disc/spatial_disc/elem_disc/thermohaline_flow/thermohaline_flow.h"
#include "lib_disc/spatial_disc/disc_util/conv_shape_interface.h"
#include "lib_disc/spatial_disc/disc_util/conv_shape.h"

#include "lib_disc/spatial_disc/elem_disc/density_driven_flow/density_driven_flow.h"
#include "lib_disc/spatial_disc/elem_disc/convection_diffusion/convection_diffusion.h"
#include "lib_disc/spatial_disc/elem_disc/constant_equation/constant_equation.h"

#include "lib_disc/spatial_disc/elem_disc/neumann_boundary/neumann_boundary.h"
#include "lib_disc/spatial_disc/elem_disc/inner_boundary/inner_boundary.h"
#include "lib_disc/spatial_disc/elem_disc/inner_boundary/FV1CalciumERElemDisc.h"

using namespace std;

namespace ug
{

namespace bridge
{

////////////////////////////////////////////////////////////////////////////////
//	Registering for all implementations of IElemDisc
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain>
void Register__Domain(Registry& reg, string grp)
{
//	dimension of domain
	static const int dim = TDomain::dim;

//	suffix and tag
	string dimSuffix = GetDomainSuffix<dim>();
	string dimTag = GetDomainTag<dim>();

	string approxGrp = grp; approxGrp.append("/ApproximationSpace");
	string elemGrp = grp; elemGrp.append("/SpatialDisc/ElemDisc");

//	DomainElemDisc base class
	{
		typedef IDomainElemDisc<TDomain> T;
		typedef IElemDisc TBase;
		string name = string("IDomainElemDisc").append(dimSuffix);
		reg.add_class_<T, TBase >(name, elemGrp);
		reg.add_class_to_group(name, "IDomainElemDisc", dimTag);
	}

//	Neumann Boundary
	{
		typedef NeumannBoundary<TDomain> T;
		typedef IDomainElemDisc<TDomain> TBase;
		string name = string("NeumannBoundary").append(dimSuffix);
		reg.add_class_<T, TBase >(name, elemGrp)
			.template add_constructor<void (*)(const char*)>("Function(s)")
			.add_method("add", static_cast<void (T::*)(number, const char*, const char*)>(&T::add))
			.add_method("add", static_cast<void (T::*)(SmartPtr<IPData<number, dim> >, const char*, const char*)>(&T::add))
			.add_method("add", static_cast<void (T::*)(SmartPtr<IPData<number, dim, bool> >, const char*, const char*)>(&T::add))
			.add_method("add", static_cast<void (T::*)(SmartPtr<IPData<MathVector<dim>, dim> >, const char*, const char*)>(&T::add))
#ifdef UG_FOR_LUA
			.add_method("add", static_cast<void (T::*)(const char*, const char*, const char*)>(&T::add))
#endif
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NeumannBoundary", dimTag);
	}

//	Inner Boundaries
	{
		typedef FV1InnerBoundaryElemDisc<TDomain> T;
		typedef IDomainElemDisc<TDomain> TBase;
		string name = string("FV1InnerBoundary").append(dimSuffix);
		reg.add_class_<T, TBase >(name, elemGrp);
		reg.add_class_to_group(name, "FV1InnerBoundary", dimTag);
	
		typedef FV1CalciumERElemDisc<TDomain> T1;
		typedef FV1InnerBoundaryElemDisc<TDomain> TBase1;
		name = string("FV1InnerBoundaryCalciumER").append(dimSuffix);
		reg.add_class_<T1, TBase1>(name, elemGrp)
			.template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "FV1InnerBoundaryCalciumER", dimTag);
	
	}

//	Constant Equation
	{
		typedef ConstantEquation<TDomain> T;
		typedef IDomainElemDisc<TDomain> TBase;
		string name = string("ConstantEquation").append(dimSuffix);
		reg.add_class_<T, TBase >(name, elemGrp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")
			.add_method("set_velocity", static_cast<void (T::*)(SmartPtr<IPData<MathVector<dim>, dim> >)>(&T::set_velocity), "", "Velocity Field")
			.add_method("set_velocity", static_cast<void (T::*)(number)>(&T::set_velocity), "", "Vel_x")
			.add_method("set_velocity", static_cast<void (T::*)(number,number)>(&T::set_velocity), "", "Vel_x, Vel_y")
			.add_method("set_velocity", static_cast<void (T::*)(number,number,number)>(&T::set_velocity), "", "Vel_x, Vel_y, Vel_z")
#ifdef UG_FOR_LUA
			.add_method("set_velocity", static_cast<void (T::*)(const char*)>(&T::set_velocity), "", "Velocity Field")
#endif

			.add_method("set_source", static_cast<void (T::*)(SmartPtr<IPData<number, dim> >)>(&T::set_source), "", "Source")
			.add_method("set_source", static_cast<void (T::*)(number)>(&T::set_source), "", "Source")
#ifdef UG_FOR_LUA
			.add_method("set_source", static_cast<void (T::*)(const char*)>(&T::set_source), "", "Source")
#endif

			.add_method("set_mass", static_cast<void (T::*)(SmartPtr<IPData<number, dim> >)>(&T::set_mass), "", "Mass")
			.add_method("set_mass", static_cast<void (T::*)(number)>(&T::set_mass), "", "Mass")
#ifdef UG_FOR_LUA
			.add_method("set_mass", static_cast<void (T::*)(const char*)>(&T::set_mass), "", "Mass")
#endif
			.add_method("value", &T::value)
			.add_method("gradient", &T::gradient)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ConstantEquation", dimTag);
	}

//	Convection Diffusion
	{
		typedef ConvectionDiffusion<TDomain> T;
		typedef IDomainElemDisc<TDomain> TBase;
		string name = string("ConvectionDiffusion").append(dimSuffix);
		reg.add_class_<T, TBase >(name, elemGrp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")
			.add_method("set_disc_scheme", &T::set_disc_scheme, "", "Disc Scheme|selection|value=[\"fe\",\"fv\",\"fv1\"]")
			.add_method("set_quad_order", &T::set_quad_order)
			.add_method("set_quad_order_scvf", &T::set_quad_order_scvf)
			.add_method("set_quad_order_scv", &T::set_quad_order_scv)

			.add_method("set_diffusion", static_cast<void (T::*)(SmartPtr<IPData<MathMatrix<dim, dim>, dim> >)>(&T::set_diffusion), "", "Diffusion")
			.add_method("set_diffusion", static_cast<void (T::*)(number)>(&T::set_diffusion), "", "Diagonal Diffusion")
#ifdef UG_FOR_LUA
			.add_method("set_diffusion", static_cast<void (T::*)(const char*)>(&T::set_diffusion), "", "Diffusion")
#endif

			.add_method("set_velocity", static_cast<void (T::*)(SmartPtr<IPData<MathVector<dim>, dim> >)>(&T::set_velocity), "", "Velocity Field")
			.add_method("set_velocity", static_cast<void (T::*)(number)>(&T::set_velocity), "", "Vel_x")
			.add_method("set_velocity", static_cast<void (T::*)(number,number)>(&T::set_velocity), "", "Vel_x, Vel_y")
			.add_method("set_velocity", static_cast<void (T::*)(number,number,number)>(&T::set_velocity), "", "Vel_x, Vel_y, Vel_z")
#ifdef UG_FOR_LUA
			.add_method("set_velocity", static_cast<void (T::*)(const char*)>(&T::set_velocity), "", "Velocity Field")
#endif

			.add_method("set_reaction_rate", static_cast<void (T::*)(SmartPtr<IPData<number, dim> >)>(&T::set_reaction_rate), "", "Reaction Rate")
			.add_method("set_reaction_rate", static_cast<void (T::*)(number)>(&T::set_reaction_rate), "", "Reaction Rate")
#ifdef UG_FOR_LUA
			.add_method("set_reaction_rate", static_cast<void (T::*)(const char*)>(&T::set_reaction_rate), "", "Reaction Rate")
#endif

			.add_method("set_reaction", static_cast<void (T::*)(SmartPtr<IPData<number, dim> >)>(&T::set_reaction), "", "Reaction")
			.add_method("set_reaction", static_cast<void (T::*)(number)>(&T::set_reaction), "", "Reaction")
#ifdef UG_FOR_LUA
			.add_method("set_reaction", static_cast<void (T::*)(const char*)>(&T::set_reaction), "", "Reaction")
#endif

			.add_method("set_source", static_cast<void (T::*)(SmartPtr<IPData<number, dim> >)>(&T::set_source), "", "Source")
			.add_method("set_source", static_cast<void (T::*)(number)>(&T::set_source), "", "Source")
#ifdef UG_FOR_LUA
			.add_method("set_source", static_cast<void (T::*)(const char*)>(&T::set_source), "", "Source")
#endif

			.add_method("set_mass_scale", static_cast<void (T::*)(SmartPtr<IPData<number, dim> >)>(&T::set_mass_scale), "", "Mass Scale")
			.add_method("set_mass_scale", static_cast<void (T::*)(number)>(&T::set_mass_scale), "", "Mass Scale")
#ifdef UG_FOR_LUA
			.add_method("set_mass_scale", static_cast<void (T::*)(const char*)>(&T::set_mass_scale), "", "Mass Scale")
#endif

			.add_method("set_mass", static_cast<void (T::*)(SmartPtr<IPData<number, dim> >)>(&T::set_mass), "", "Mass")
			.add_method("set_mass", static_cast<void (T::*)(number)>(&T::set_mass), "", "Mass")
#ifdef UG_FOR_LUA
			.add_method("set_mass", static_cast<void (T::*)(const char*)>(&T::set_mass), "", "Mass")
#endif

			.add_method("set_upwind", &T::set_upwind)
			.add_method("value", &T::value)
			.add_method("gradient", &T::gradient)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ConvectionDiffusion", dimTag);
	}

//	Density Driven Flow
	{
		typedef DensityDrivenFlow<TDomain> T;
		typedef IDomainElemDisc<TDomain> TBase;
		string name = string("DensityDrivenFlow").append(dimSuffix);
		reg.add_class_<T, TBase >(name, elemGrp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")
			.add_method("set_consistent_gravity", &T::set_consistent_gravity, "", "Consistent Gravity")
			.add_method("set_boussinesq_transport", &T::set_boussinesq_transport, "", "Boussinesq Transport")
			.add_method("set_boussinesq_flow", &T::set_boussinesq_flow,"", "Boussinesq Flow")
			.add_method("set_upwind", &T::set_upwind, "", "Upwind (no, part, full)")

			.add_method("set_porosity", static_cast<void (T::*)(SmartPtr<IPData<number, dim> >)>(&T::set_porosity), "", "Porosity")
			.add_method("set_porosity", static_cast<void (T::*)(number)>(&T::set_porosity), "", "Porosity")
#ifdef UG_FOR_LUA
			.add_method("set_porosity", static_cast<void (T::*)(const char*)>(&T::set_porosity), "", "Porosity")
#endif

			.add_method("set_gravity", static_cast<void (T::*)(SmartPtr<IPData<MathVector<dim>, dim> >)>(&T::set_gravity), "", "Gravity")
			.add_method("set_gravity", static_cast<void (T::*)(number)>(&T::set_gravity), "", "g_x")
			.add_method("set_gravity", static_cast<void (T::*)(number,number)>(&T::set_gravity), "", "g_x, g_y")
			.add_method("set_gravity", static_cast<void (T::*)(number,number,number)>(&T::set_gravity), "", "g_x, g_y, g_z")
#ifdef UG_FOR_LUA
			.add_method("set_gravity", static_cast<void (T::*)(const char*)>(&T::set_gravity), "", "Gravity")
#endif

			.add_method("set_molecular_diffusion", static_cast<void (T::*)(SmartPtr<IPData<MathMatrix<dim, dim>, dim> >)>(&T::set_molecular_diffusion), "", "MolecularDiffusion")
			.add_method("set_molecular_diffusion", static_cast<void (T::*)(number)>(&T::set_molecular_diffusion), "", "MolecularDiffusion")
#ifdef UG_FOR_LUA
			.add_method("set_molecular_diffusion", static_cast<void (T::*)(const char*)>(&T::set_molecular_diffusion), "", "MolecularDiffusion")
#endif

			.add_method("set_permeability", static_cast<void (T::*)(SmartPtr<IPData<MathMatrix<dim, dim>, dim> >)>(&T::set_permeability), "", "Permeability")
			.add_method("set_permeability", static_cast<void (T::*)(number)>(&T::set_permeability), "", "Permeability")
#ifdef UG_FOR_LUA
			.add_method("set_permeability", static_cast<void (T::*)(const char*)>(&T::set_permeability), "", "Permeability")
#endif

			.add_method("set_viscosity", static_cast<void (T::*)(SmartPtr<IPData<number, dim> >)>(&T::set_viscosity), "", "Viscosity")
			.add_method("set_viscosity", static_cast<void (T::*)(number)>(&T::set_viscosity), "", "Viscosity")
#ifdef UG_FOR_LUA
			.add_method("set_viscosity", static_cast<void (T::*)(const char*)>(&T::set_viscosity), "", "Viscosity")
#endif

			.add_method("set_density", &T::set_density, "", "Density")
			.add_method("darcy_velocity", &T::darcy_velocity)
			.add_method("brine", &T::brine)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "DensityDrivenFlow", dimTag);
	}

//	Thermohaline Flow
	{
		typedef ThermohalineFlow<TDomain> T;
		typedef IDomainElemDisc<TDomain> TBase;
		string name = string("ThermohalineFlow").append(dimSuffix);
		reg.add_class_<T, TBase >(name, elemGrp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")
			.add_method("set_consistent_gravity", &T::set_consistent_gravity, "", "ConsistentGravity")
			.add_method("set_boussinesq_transport", &T::set_boussinesq_transport, "", "BoussinesqTransport")
			.add_method("set_boussinesq_flow", &T::set_boussinesq_flow,"", "BoussinesqFlow")
			.add_method("set_boussinesq_density", &T::set_boussinesq_density,"", "BoussinesqDensity")
			.add_method("set_upwind", &T::set_upwind, "", "Upwind (no, part, full)")
			.add_method("set_upwind_energy", &T::set_upwind_energy, "", "Upwind (no, part, full)")

			.add_method("set_porosity", static_cast<void (T::*)(SmartPtr<IPData<number, dim> >)>(&T::set_porosity), "", "Porosity")
			.add_method("set_porosity", static_cast<void (T::*)(number)>(&T::set_porosity), "", "Porosity")
#ifdef UG_FOR_LUA
			.add_method("set_porosity", static_cast<void (T::*)(const char*)>(&T::set_porosity), "", "Porosity")
#endif

			.add_method("set_gravity", static_cast<void (T::*)(SmartPtr<IPData<MathVector<dim>, dim> >)>(&T::set_gravity), "", "Gravity")
			.add_method("set_gravity", static_cast<void (T::*)(number)>(&T::set_gravity), "", "g_x")
			.add_method("set_gravity", static_cast<void (T::*)(number,number)>(&T::set_gravity), "", "g_x, g_y")
			.add_method("set_gravity", static_cast<void (T::*)(number,number,number)>(&T::set_gravity), "", "g_x, g_y, g_z")
#ifdef UG_FOR_LUA
			.add_method("set_gravity", static_cast<void (T::*)(const char*)>(&T::set_gravity), "", "Gravity")
#endif

			.add_method("set_molecular_diffusion", static_cast<void (T::*)(SmartPtr<IPData<MathMatrix<dim, dim>, dim> >)>(&T::set_molecular_diffusion), "", "MolecularDiffusion")
			.add_method("set_molecular_diffusion", static_cast<void (T::*)(number)>(&T::set_molecular_diffusion), "", "MolecularDiffusion")
#ifdef UG_FOR_LUA
			.add_method("set_molecular_diffusion", static_cast<void (T::*)(const char*)>(&T::set_molecular_diffusion), "", "MolecularDiffusion")
#endif

			.add_method("set_thermal_conductivity", static_cast<void (T::*)(SmartPtr<IPData<MathMatrix<dim, dim>, dim> >)>(&T::set_thermal_conductivity), "", "ThermalConductivity")
			.add_method("set_thermal_conductivity", static_cast<void (T::*)(number)>(&T::set_thermal_conductivity), "", "ThermalConductivity")
#ifdef UG_FOR_LUA
			.add_method("set_thermal_conductivity", static_cast<void (T::*)(const char*)>(&T::set_thermal_conductivity), "", "ThermalConductivity")
#endif

			.add_method("set_permeability", static_cast<void (T::*)(SmartPtr<IPData<MathMatrix<dim, dim>, dim> >)>(&T::set_permeability), "", "Permeability")
			.add_method("set_permeability", static_cast<void (T::*)(number)>(&T::set_permeability), "", "Permeability")
#ifdef UG_FOR_LUA
			.add_method("set_permeability", static_cast<void (T::*)(const char*)>(&T::set_permeability), "", "Permeability")
#endif

			.add_method("set_viscosity", static_cast<void (T::*)(SmartPtr<IPData<number, dim> >)>(&T::set_viscosity), "", "Viscosity")
			.add_method("set_viscosity", static_cast<void (T::*)(number)>(&T::set_viscosity), "", "Viscosity")
#ifdef UG_FOR_LUA
			.add_method("set_viscosity", static_cast<void (T::*)(const char*)>(&T::set_viscosity), "", "Viscosity")
#endif

			.add_method("set_density", &T::set_density, "", "Density")
			.add_method("set_heat_capacity_fluid", &T::set_heat_capacity_fluid)
			.add_method("set_heat_capacity_solid", &T::set_heat_capacity_solid)
			.add_method("set_mass_density_solid", &T::set_mass_density_solid)
			.add_method("darcy_velocity", &T::darcy_velocity)
			.add_method("pressure_grad", &T::pressure_grad)
			.add_method("temperature", &T::temperature)
			.add_method("brine", &T::brine)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ThermohalineFlow", dimTag);
	}

/////////////////////////////////////////////////////////////////////////////
// Convection Shapes
/////////////////////////////////////////////////////////////////////////////

	string upGrp = grp; upGrp.append("/SpatialDisc/Upwind");

//	IConvectionShapes
	{
		typedef IConvectionShapes<dim> T;
		string name = string("IConvectionShapes").append(dimSuffix);
		reg.add_class_<T>(name, upGrp);
		reg.add_class_to_group(name, "IConvectionShapes", dimTag);
	}

//	ConvectionShapesNoUpwind
	{
		typedef ConvectionShapesNoUpwind<dim> T;
		typedef IConvectionShapes<dim> TBase;
		string name = string("NoUpwind").append(dimSuffix);
		reg.add_class_<T, TBase>(name, upGrp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "NoUpwind", dimTag);
	}

//	ConvectionShapesFullUpwind
	{
		typedef ConvectionShapesFullUpwind<dim> T;
		typedef IConvectionShapes<dim> TBase;
		string name = string("FullUpwind").append(dimSuffix);
		reg.add_class_<T, TBase>(name, upGrp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "FullUpwind", dimTag);
	}

//	ConvectionShapesWeightedUpwind
	{
		typedef ConvectionShapesWeightedUpwind<dim> T;
		typedef IConvectionShapes<dim> TBase;
		string name = string("WeightedUpwind").append(dimSuffix);
		reg.add_class_<T, TBase>(name, upGrp)
			.add_method("set_weight", &T::set_weight)
			.add_constructor()
			.template add_constructor<void (*)(number)>("weight")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "WeightedUpwind", dimTag);
	}

//	ConvectionShapesPartialUpwind
	{
		typedef ConvectionShapesPartialUpwind<dim> T;
		typedef IConvectionShapes<dim> TBase;
		string name = string("PartialUpwind").append(dimSuffix);
		reg.add_class_<T, TBase>(name, upGrp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "PartialUpwind", dimTag);
	}
}


bool RegisterElemDiscs(Registry& reg, string parentGroup)
{
//	get group string
	string grp = parentGroup; grp.append("/Discretization");

//	Elem Discs
	{
		string elemGrp = grp; elemGrp.append("/ElemDisc");
		typedef IElemDisc T;
		reg.add_class_<T>("IElemDisc", elemGrp);
	}

	try
	{
#ifdef UG_DIM_1
		Register__Domain<Domain1d>(reg, grp);
#endif
#ifdef UG_DIM_2
			Register__Domain<Domain2d>(reg, grp);
#endif
#ifdef UG_DIM_3
			Register__Domain<Domain3d>(reg, grp);
#endif
	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterElemDiscs: "
				"Registration failed (using name " << ex.name << ").\n");
		UG_THROW("Registration failed.");
	}

	return true;
}

}//	end of namespace ug
}//	end of namespace interface
