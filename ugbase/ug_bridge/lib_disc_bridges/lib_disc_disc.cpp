/*
 * thermohaline_flow_bridge.cpp
 *
 *  Created on: 20.05.2011
 *      Author: andreasvogel
 */

// extern headers
#include <iostream>
#include <sstream>
#include <string>

// include bridge
#include "../ug_bridge.h"
#include "registry/registry.h"

// lib_disc includes
#include "lib_discretization/domain.h"

#include "lib_discretization/spatial_discretization/elem_disc/thermohaline_flow/thermohaline_flow.h"
#include "lib_discretization/spatial_discretization/disc_util/conv_shape_interface.h"
#include "lib_discretization/spatial_discretization/disc_util/conv_shape.h"

#include "lib_discretization/spatial_discretization/elem_disc/navier_stokes/navier_stokes.h"
#include "lib_discretization/spatial_discretization/elem_disc/navier_stokes/upwind.h"
#include "lib_discretization/spatial_discretization/elem_disc/navier_stokes/stabilization.h"

#include "lib_discretization/spatial_discretization/elem_disc/density_driven_flow/density_driven_flow.h"
#include "lib_discretization/spatial_discretization/elem_disc/convection_diffusion/convection_diffusion.h"
#include "lib_discretization/spatial_discretization/elem_disc/constant_equation/constant_equation.h"

#include "lib_discretization/spatial_discretization/elem_disc/neumann_boundary/neumann_boundary.h"
#include "lib_discretization/spatial_discretization/elem_disc/inner_boundary/inner_boundary.h"
#include "lib_discretization/spatial_discretization/elem_disc/linear_elasticity/fe1_linear_elasticity.h"

// fe1_nonlinear_elasticity includes
#include <boost/function.hpp>
#include "lib_discretization/spatial_discretization/ip_data/ip_data.h"
#include "lib_discretization/spatial_discretization/elem_disc/nonlinear_elasticity/fe1_nonlinear_elasticity.h"

using namespace std;

namespace ug
{

namespace bridge
{

////////////////////////////////////////////////////////////////////////////////
//	Some classes for Non-Linear Elastictity
////////////////////////////////////////////////////////////////////////////////

template <int dim>
class ElasticityTensorUserData
	: public IPData<MathTensor<4,dim>, dim>,
	  public boost::function<void (MathTensor<4,dim>& value, const MathVector<dim>& x, number time)>
{
	///	Base class type
		typedef IPData<MathTensor<4,dim>, dim> base_type;

	///	Functor type
		typedef boost::function<void (MathTensor<4,dim>& value,
		                              const MathVector<dim>& x,
		                              number time)> TensorFunctor;

		using base_type::num_series;
		using base_type::num_ip;
		using base_type::ip;
		using base_type::time;
		using base_type::value;

	public:
	///	Constructor
		ElasticityTensorUserData() : TensorFunctor(boost::ref(*this)) {}

	///	virtual destructor
		virtual ~ElasticityTensorUserData()	{}

	///	evaluates the data at a given point and time
		void operator() (MathTensor<4,dim>& C, const MathVector<dim>& x, number time = 0.0)
		{
			//filling the ElasticityTensor
			C.set(0.0);

			C[0][0][0][0] = 5.;
			C[0][0][1][2] = 12.;
		}
	///	prints Elasticity Tensor C at point x and time t
		void test_evaluate(number x, number y, number z, number t)
		{
			MathTensor<4,dim> C;
			MathVector<dim> v;
			v.x = x;
			if(v.size() > 1)
				v[1] = y;
			if(v.size() > 2)
				v[2] = z;

			this->operator()(C, v, t);
			//UG_LOG("Tensor: " << C);

		}
	///	implement as a IPData
		virtual bool compute(bool bDeriv = false)
		{
			for(size_t s = 0; s < num_series(); ++s)
				for(size_t i = 0; i < num_ip(s); ++i)
				{
					this->operator() (	value(s,i),
										ip(s, i),
										time());
				}
			return true;
		}
};


////////////////////////////////////////////////////////////////////////////////
//	Registering for all implementations of IElemDisc
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain>
void RegisterIElemDiscs(Registry& reg, string grp)
{
//	dimension of domain
	static const int dim = TDomain::dim;

//	suffix and tag
	string dimSuffix = GetDomainSuffix<dim>();
	string dimTag = GetDomainTag<dim>();

//	IApproximationSpace
	{
		typedef IApproximationSpace<TDomain> T;
		typedef FunctionPattern TBase;
		string name = string("IApproximationSpace").append(dimSuffix);
		reg.add_class_<T, TBase >(name, grp)
			.add_method("assign_domain|hide=true", &T::assign_domain)
			.add_method("get_domain|hide=true", static_cast<TDomain& (T::*)()>(&T::get_domain));
		reg.add_class_to_group(name, "IApproximationSpace", dimTag);
	}

//	DomainElemDisc base class
	{
		typedef IDomainElemDisc<TDomain> T;
		typedef IElemDisc TBase;
		string name = string("IDomainElemDisc").append(dimSuffix);
		reg.add_class_<T, TBase >(name, grp)
			.add_method("set_approximation_space", &T::set_approximation_space);
		reg.add_class_to_group(name, "IDomainElemDisc", dimTag);
	}

//	Neumann Boundary
	{
		typedef boost::function<bool (number& value, const MathVector<dim>& x, number time)> BNDNumberFunctor;
		typedef FV1NeumannBoundaryElemDisc<TDomain> T;
		typedef IDomainElemDisc<TDomain> TBase;
		string name = string("FV1NeumannBoundary").append(dimSuffix);
		reg.add_class_<T, TBase >(name, grp)
			.add_constructor()
			.add_method("add", static_cast<bool (T::*)(BNDNumberFunctor&, const char*, const char*)>(&T::add));
		reg.add_class_to_group(name, "IDomainElemDisc", dimTag);
	}

//	Inner Boundary
	{
		typedef FVInnerBoundaryElemDisc<TDomain> T;
		typedef IDomainElemDisc<TDomain> TBase;
		string name = string("FV1InnerBoundary").append(dimSuffix);
		reg.add_class_<T, TBase >(name, grp)
			.add_constructor();
		reg.add_class_to_group(name, "FV1InnerBoundary", dimTag);
	}

//	Constant Equation Finite Volume
	{
		typedef FVConstantEquationElemDisc<TDomain> T;
		typedef IDomainElemDisc<TDomain> TBase;
		string name = string("FV1ConstantEquation").append(dimSuffix);
		reg.add_class_<T, TBase >(name, grp)
			.add_constructor()
			.add_method("set_velocity", &T::set_velocity)
			.add_method("set_source", &T::set_source)
			.add_method("set_mass_scale", &T::set_mass_scale)
			.add_method("get_concentration", &T::get_concentration)
			.add_method("get_concentration_grad", &T::get_concentration_grad);
		reg.add_class_to_group(name, "FV1ConstantEquation", dimTag);
	}

//	Convection Diffusion
	{
		typedef ConvectionDiffusionElemDisc<TDomain> T;
		typedef IDomainElemDisc<TDomain> TBase;
		string name = string("ConvectionDiffusion").append(dimSuffix);
		reg.add_class_<T, TBase >(name, grp)
			.add_constructor()
			.add_method("set_disc_scheme", &T::set_disc_scheme)
			.add_method("set_diffusion_tensor", &T::set_diffusion)
			.add_method("set_velocity_field", &T::set_velocity)
			.add_method("set_reaction", &T::set_reaction)
			.add_method("set_source", &T::set_source)
			.add_method("set_mass_scale", &T::set_mass_scale)
			.add_method("set_upwind", &T::set_upwind)
			.add_method("get_concentration", &T::get_concentration)
			.add_method("get_concentration_grad", &T::get_concentration_grad);
		reg.add_class_to_group(name, "ConvectionDiffusion", dimTag);
	}

//	Density Driven Flow
	{
		typedef DensityDrivenFlowElemDisc<TDomain> T2;
		typedef IDomainElemDisc<TDomain> TBase;
		string name = string("DensityDrivenFlow").append(dimSuffix);
		reg.add_class_<T2, TBase >(name, grp)
			.add_constructor()
			.add_method("set_upwind|interactive=false", &T2::set_upwind,
						"", "Upwind (no, part, full)")
			.add_method("set_boussinesq_transport|interactive=false", &T2::set_boussinesq_transport,
						"", "Boussinesq Transport")
			.add_method("set_boussinesq_flow|interactive=false", &T2::set_boussinesq_flow,
						"", "Boussinesq Flow")
			.add_method("set_porosity|interactive=false", &T2::set_porosity,
						"", "Porosity")
			.add_method("set_gravity|interactive=false", &T2::set_gravity,
						"", "Gravity")
			.add_method("set_permeability|interactive=false", &T2::set_permeability,
						"", "Permeability")
			.add_method("set_viscosity|interactive=false", &T2::set_viscosity,
						"", "Viscosity")
			.add_method("set_molecular_diffusion|interactive=false", &T2::set_molecular_diffusion,
						"", "Molecular Diffusion")
			.add_method("set_density|interactive=false", &T2::set_density,
						"", "Density")
			.add_method("set_consistent_gravity|interactive=false", &T2::set_consistent_gravity,
						"", "Consistent Gravity")
			.add_method("get_darcy_velocity", &T2::get_darcy_velocity)
			.add_method("get_brine", &T2::get_brine);
		reg.add_class_to_group(name, "DensityDrivenFlow", dimTag);
	}

//	Thermohaline Flow
	{
		typedef ThermohalineFlowElemDisc<TDomain> T2;
		typedef IDomainElemDisc<TDomain> TBase;
		string name = string("FV1ThermohalineFlow").append(dimSuffix);
		reg.add_class_<T2, TBase >(name, grp)
			.add_constructor()
			.add_method("set_upwind", &T2::set_upwind,
						"", "Upwind (no, part, full)")
			.add_method("set_upwind_energy", &T2::set_upwind_energy,
						"", "Upwind (no, part, full)")
			.add_method("set_boussinesq_transport", &T2::set_boussinesq_transport,
						"", "Boussinesq Transport")
			.add_method("set_boussinesq_flow", &T2::set_boussinesq_flow,
						"", "Boussinesq Flow")
			.add_method("set_boussinesq_density", &T2::set_boussinesq_density,
						"", "Boussinesq Density")
			.add_method("set_porosity", &T2::set_porosity,
						"", "Porosity")
			.add_method("set_gravity", &T2::set_gravity,
						"", "Gravity")
			.add_method("set_permeability", &T2::set_permeability,
						"", "Permeability")
			.add_method("set_viscosity", &T2::set_viscosity,
						"", "Viscosity")
			.add_method("set_thermal_conductivity", &T2::set_thermal_conductivity,
						"", "Molecular Diffusion")
			.add_method("set_molecular_diffusion", &T2::set_molecular_diffusion,
						"", "Molecular Diffusion")
			.add_method("set_density", &T2::set_density,
						"", "Density")
			.add_method("set_consistent_gravity", &T2::set_consistent_gravity,
						"", "Consistent Gravity")
			.add_method("set_heat_capacity_fluid", &T2::set_heat_capacity_fluid)
			.add_method("set_heat_capacity_solid", &T2::set_heat_capacity_solid)
			.add_method("set_mass_density_solid", &T2::set_mass_density_solid)
			.add_method("get_darcy_velocity", &T2::get_darcy_velocity)
			.add_method("get_pressure_grad", &T2::get_pressure_grad)
			.add_method("get_temperature", &T2::get_temperature)
			.add_method("get_brine", &T2::get_brine);
		reg.add_class_to_group(name, "FV1ThermohalineFlow", dimTag);
	}

//	Navier-Stokes
	{
		typedef FVNavierStokesElemDisc<TDomain> T;
		typedef IDomainElemDisc<TDomain> TBase;
		string name = string("FV1NavierStokes").append(dimSuffix);
		reg.add_class_<T, TBase >(name, grp)
			.add_constructor()
			.add_method("set_kinematic_viscosity", &T::set_kinematic_viscosity)
			.add_method("set_stabilization", &T::set_stabilization)
			.add_method("set_conv_upwind",  static_cast<void (T::*)(INavierStokesStabilization<dim>&)>(&T::set_conv_upwind))
			.add_method("set_conv_upwind",  static_cast<void (T::*)(INavierStokesUpwind<dim>&)>(&T::set_conv_upwind))
			.add_method("set_peclet_blend", &T::set_peclet_blend)
			.add_method("set_exact_jacobian", &T::set_exact_jacobian);
		reg.add_class_to_group(name, "FV1NavierStokes", dimTag);
	}


//	Non-Linear Elastictity
	{
		typedef FE1NonlinearElasticityElemDisc<TDomain> T;
		typedef IDomainElemDisc<TDomain> TBase;
		string name = string("FE1NonlinearElasticity").append(dimSuffix);
		reg.add_class_<T, TBase >(name, grp)
			.add_constructor()
			.add_method("set_elasticity_tensor", &T::set_elasticity_tensor);
			//methods, which are used in script-file
		reg.add_class_to_group(name, "FE1NonlinearElasticity", dimTag);
	}


/////////////////////////////////////////////////////////////////////////////
// Upwind
/////////////////////////////////////////////////////////////////////////////


//	INavierStokesUpwind
	{
		typedef INavierStokesUpwind<dim> T;
		string name = string("INavierStokesUpwind").append(dimSuffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, "INavierStokesUpwind", dimTag);
	}

//	NavierStokesNoUpwind
	{
		typedef NavierStokesNoUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		string name = string("NavierStokesNoUpwind").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor();
		reg.add_class_to_group(name, "NavierStokesNoUpwind", dimTag);
	}

//	NavierStokesFullUpwind
	{
		typedef NavierStokesFullUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		string name = string("NavierStokesFullUpwind").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor();
		reg.add_class_to_group(name, "NavierStokesFullUpwind", dimTag);
	}

//	NavierStokesSkewedUpwind
	{
		typedef NavierStokesSkewedUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		string name = string("NavierStokesSkewedUpwind").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor();
		reg.add_class_to_group(name, "NavierStokesSkewedUpwind", dimTag);
	}

//	NavierStokesLinearProfileSkewedUpwind
	{
		typedef NavierStokesLinearProfileSkewedUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		string name = string("NavierStokesLinearProfileSkewedUpwind").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor();
		reg.add_class_to_group(name, "NavierStokesLinearProfileSkewedUpwind", dimTag);
	}

//	NavierStokesPositiveUpwind
	{
		typedef NavierStokesPositiveUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		string name = string("NavierStokesPositiveUpwind").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor();
		reg.add_class_to_group(name, "NavierStokesPositiveUpwind", dimTag);
	}

/////////////////////////////////////////////////////////////////////////////
// Stabilization
/////////////////////////////////////////////////////////////////////////////


//	INavierStokesStabilization
	{
		typedef INavierStokesStabilization<dim> T;
		string name = string("INavierStokesStabilization").append(dimSuffix);
		reg.add_class_<T>(name, grp)
			.add_method("set_upwind", &T::set_upwind)
			.add_method("set_diffusion_length", &T::set_diffusion_length);
		reg.add_class_to_group(name, "INavierStokesStabilization", dimTag);
	}

//	NavierStokesFIELDSStabilization
	{
		typedef NavierStokesFIELDSStabilization<dim> T;
		typedef INavierStokesStabilization<dim> TBase;
		string name = string("NavierStokesFIELDSStabilization").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor();
		reg.add_class_to_group(name, "NavierStokesFIELDSStabilization", dimTag);
	}

/////////////////////////////////////////////////////////////////////////////
// Convection Shapes
/////////////////////////////////////////////////////////////////////////////

//	IConvectionShapes
	{
		typedef IConvectionShapes<dim> T;
		string name = string("IConvectionShapes").append(dimSuffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, "IConvectionShapes", dimTag);
	}

//	ConvectionShapesNoUpwind
	{
		typedef ConvectionShapesNoUpwind<dim> T;
		typedef IConvectionShapes<dim> TBase;
		string name = string("NoUpwind").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor();
		reg.add_class_to_group(name, "NoUpwind", dimTag);
	}

//	ConvectionShapesFullUpwind
	{
		typedef ConvectionShapesFullUpwind<dim> T;
		typedef IConvectionShapes<dim> TBase;
		string name = string("FullUpwind").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor();
		reg.add_class_to_group(name, "FullUpwind", dimTag);
	}

//	ConvectionShapesWeightedUpwind
	{
		typedef ConvectionShapesWeightedUpwind<dim> T;
		typedef IConvectionShapes<dim> TBase;
		string name = string("WeightedUpwind").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_method("set_weight", &T::set_weight)
			.add_constructor();
		reg.add_class_to_group(name, "WeightedUpwind", dimTag);
	}

//	ConvectionShapesPartialUpwind
	{
		typedef ConvectionShapesPartialUpwind<dim> T;
		typedef IConvectionShapes<dim> TBase;
		string name = string("PartialUpwind").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor();
		reg.add_class_to_group(name, "PartialUpwind", dimTag);
	}

/////////////////////////////////////////////////////////////////////////////
// ElasticityTensorUserData
/////////////////////////////////////////////////////////////////////////////
	{
		typedef ElasticityTensorUserData<dim> T;
		typedef IPData<MathTensor<4,dim>, dim> TBase;
		typedef boost::function<void (MathTensor<4,dim>& value, const MathVector<dim>& x, number time)> TBase2;
		string name = string("ElasticityTensorUserData").append(dimSuffix);
		reg.add_class_<T, TBase, TBase2>(name, grp)
			.add_constructor()
			//.add_method("operator", &T::operator ())
			.add_method("test_evaluate", &T::test_evaluate);
			//methods, which are used in script-file
		reg.add_class_to_group(name, "ElasticityTensorUserData", dimTag);
	}

}


bool RegisterLibDiscElemDisc(Registry& reg, const char* parentGroup)
{
	try
	{
	//	get group string
		string grp = parentGroup; grp.append("/Discretization");

#ifdef UG_DIM_1
	//	Domain dependend part 1D
			RegisterIElemDiscs<Domain1d>(reg, grp);
#endif

#ifdef UG_DIM_2
	//	Domain dependend part 2D
			RegisterIElemDiscs<Domain2d>(reg, grp);
#endif

#ifdef UG_DIM_3
	//	Domain dependend part 3D
			RegisterIElemDiscs<Domain3d>(reg, grp);
#endif
	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterLibDiscretizationIElemDisc: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}

	return true;
}

}//	end of namespace ug
}//	end of namespace interface
