/*
 * thermohaline_flow_bridge.cpp
 *
 *  Created on: 20.05.2011
 *      Author: andreasvogel
 */

// extern headers
#include <iostream>
#include <sstream>

// include bridge
#include "../../ug_bridge.h"
#include "../../registry.h"

// lib_disc includes
#include "lib_discretization/domain.h"

#include "lib_discretization/spatial_discretization/elem_disc/thermohaline_flow/thermohaline_flow.h"
#include "lib_discretization/spatial_discretization/disc_util/conv_shape_interface.h"
#include "lib_discretization/spatial_discretization/disc_util/conv_shape.h"

#include "lib_discretization/spatial_discretization/elem_disc/navier_stokes/navier_stokes.h"
#include "lib_discretization/spatial_discretization/elem_disc/navier_stokes/upwind.h"
#include "lib_discretization/spatial_discretization/elem_disc/navier_stokes/stabilization.h"

#include "lib_discretization/spatial_discretization/elem_disc/density_driven_flow/density_driven_flow.h"

#include "lib_discretization/spatial_discretization/elem_disc/convection_diffusion/fe1_convection_diffusion.h"
#include "lib_discretization/spatial_discretization/elem_disc/convection_diffusion/convection_diffusion.h"

#include "lib_discretization/spatial_discretization/elem_disc/constant_equation/constant_equation.h"

#include "lib_discretization/spatial_discretization/elem_disc/neumann_boundary/neumann_boundary.h"
#include "lib_discretization/spatial_discretization/elem_disc/inner_boundary/inner_boundary.h"
#include "lib_discretization/spatial_discretization/elem_disc/linear_elasticity/fe1_linear_elasticity.h"

// fe1_nonlinear_elasticity includes
#include "lib_discretization/spatial_discretization/ip_data/user_data_interface.h"
#include "lib_discretization/spatial_discretization/elem_disc/nonlinear_elasticity/fe1_nonlinear_elasticity.h"

namespace ug
{

namespace bridge
{

////////////////////////////////////////////////////////////////////////////////
//	Some classes for Non-Linear Elastictity
////////////////////////////////////////////////////////////////////////////////

template <int dim>
class ElasticityTensorUserData
	: public IUserData<MathTensor<4,dim>, dim>
{
	///	Base class type
		typedef IUserData<MathTensor<4,dim>, dim> base_type;

		using base_type::num_series;
		using base_type::num_ip;
		using base_type::ip;
		using base_type::time;
		using base_type::value;

	public:
	//	Functor Type
		typedef typename base_type::functor_type functor_type;

	//	return functor
		virtual functor_type get_functor() const {return *this;}

	public:
	///	Constructor
		ElasticityTensorUserData() {}

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
		virtual void compute(bool computeDeriv = false)
		{
			for(size_t s = 0; s < num_series(); ++s)
				for(size_t i = 0; i < num_ip(s); ++i)
				{
					this->operator() (	value(s,i),
										ip(s, i),
										time());
				}
		}
};


////////////////////////////////////////////////////////////////////////////////
//	Registering for all implementations of IElemDisc
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain>
void RegisterIElemDiscs(Registry& reg, const char* parentGroup)
{
//	typedef domain
	typedef TDomain domain_type;
	static const int dim = domain_type::dim;

	std::stringstream grpSS; grpSS << parentGroup << "/" << dim << "d";
	std::string grp = grpSS.str();

//	DomainElemDisc base class
	{
		typedef IDomainElemDisc<domain_type> T;
		typedef IElemDisc TBase;

		std::stringstream ss; ss << "IDomainElemDisc" << dim << "d";
		reg.add_class_<T, TBase >(ss.str().c_str(), grp.c_str())
			.add_method("set_approximation_space", &T::set_approximation_space);
	}

//	Neumann Boundary
	{
		typedef FV1NeumannBoundaryElemDisc<domain_type> T;
		typedef IDomainElemDisc<domain_type> TBase;
		std::stringstream ss; ss << "FV1NeumannBoundary" << dim << "d";
		reg.add_class_<T, TBase >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("add_boundary_value", (bool (T::*)(IBoundaryData<number, dim>&, const char*, const char*))&T::add_boundary_value);
	}

//	Inner Boundary
	{
		typedef FVInnerBoundaryElemDisc<domain_type> T;
		typedef IDomainElemDisc<domain_type> TBase;
		std::stringstream ss; ss << "FV1InnerBoundary" << dim << "d";
		reg.add_class_<T, TBase >(ss.str().c_str(), grp.c_str())
			.add_constructor();
	}

//	Constant Equation Finite Volume
	{
		typedef FVConstantEquationElemDisc<domain_type> T;
		typedef IDomainElemDisc<domain_type> TBase;
		std::stringstream ss; ss << "FV1ConstantEquation" << dim << "d";
		reg.add_class_<T, TBase >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set_velocity", &T::set_velocity)
			.add_method("set_source", &T::set_source)
			.add_method("set_mass_scale", &T::set_mass_scale)
			.add_method("get_concentration", &T::get_concentration)
			.add_method("get_concentration_grad", &T::get_concentration_grad);
	}

//	Convection Diffusion Finite Volume
	{
		typedef FVConvectionDiffusionElemDisc<domain_type> T;
		typedef IDomainElemDisc<domain_type> TBase;
		std::stringstream ss; ss << "FV1ConvectionDiffusion" << dim << "d";
		reg.add_class_<T, TBase >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set_diffusion_tensor", &T::set_diffusion)
			.add_method("set_velocity_field", &T::set_velocity)
			.add_method("set_reaction", &T::set_reaction)
			.add_method("set_source", &T::set_source)
			.add_method("set_mass_scale", &T::set_mass_scale)
			.add_method("set_upwind", &T::set_upwind)
			.add_method("get_concentration", &T::get_concentration)
			.add_method("get_concentration_grad", &T::get_concentration_grad);
	}

//	Convection Diffusion Finite Element
	{
		typedef FE1ConvectionDiffusionElemDisc<domain_type> T;
		typedef IDomainElemDisc<domain_type> TBase;
		std::stringstream ss; ss << "FE1ConvectionDiffusion" << dim << "d";
		reg.add_class_<T, TBase >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set_diffusion_tensor", &T::set_diffusion)
			.add_method("set_velocity_field", &T::set_velocity)
			.add_method("set_reaction", &T::set_reaction)
			.add_method("set_source", &T::set_source)
			.add_method("set_mass_scale", &T::set_mass_scale);
	}

//	Density Driven Flow
	{
		typedef DensityDrivenFlowElemDisc<domain_type> T2;
		typedef IDomainElemDisc<domain_type> TBase;
		std::stringstream ss; ss << "DensityDrivenFlow" << dim << "d";
		reg.add_class_<T2, TBase >(ss.str().c_str(), grp.c_str())
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
	}

//	Thermohaline Flow
	{
		typedef ThermohalineFlowElemDisc<domain_type> T2;
		typedef IDomainElemDisc<domain_type> TBase;
		std::stringstream ss; ss << "FV1ThermohalineFlow" << dim << "d";
		reg.add_class_<T2, TBase >(ss.str().c_str(), grp.c_str())
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
	}

//	Navier-Stokes
	{
		typedef FVNavierStokesElemDisc<domain_type> T;
		typedef IDomainElemDisc<domain_type> TBase;
		std::stringstream ss; ss << "FV1NavierStokes" << dim << "d";
		reg.add_class_<T, TBase >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set_kinematic_viscosity", &T::set_kinematic_viscosity)
			.add_method("set_stabilization", &T::set_stabilization)
			.add_method("set_conv_upwind",  (void (T::*)(INavierStokesStabilization<dim>&))&T::set_conv_upwind)
			.add_method("set_conv_upwind",  (void (T::*)(INavierStokesUpwind<dim>&))&T::set_conv_upwind)
			.add_method("set_peclet_blend", &T::set_peclet_blend)
			.add_method("set_exact_jacobian", &T::set_exact_jacobian);
	}


//	Non-Linear Elastictity
	{
		typedef FE1NonlinearElasticityElemDisc<domain_type> T;
		typedef IDomainElemDisc<domain_type> TBase;
		std::stringstream ss; ss << "FE1NonlinearElasticity" << dim << "d";
		reg.add_class_<T, TBase >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set_elasticity_tensor", &T::set_elasticity_tensor);
			//methods, which are used in script-file
	}


/////////////////////////////////////////////////////////////////////////////
// Upwind
/////////////////////////////////////////////////////////////////////////////


//	INavierStokesUpwind
	{
		typedef INavierStokesUpwind<dim> T;
		std::stringstream ss; ss << "INavierStokesUpwind" << dim << "d";
		reg.add_class_<T>(ss.str().c_str(), grp.c_str());
	}

//	NavierStokesNoUpwind
	{
		typedef NavierStokesNoUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		std::stringstream ss; ss << "NavierStokesNoUpwind" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str())
			.add_constructor();
	}

//	NavierStokesFullUpwind
	{
		typedef NavierStokesFullUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		std::stringstream ss; ss << "NavierStokesFullUpwind" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str())
			.add_constructor();
	}

//	NavierStokesSkewedUpwind
	{
		typedef NavierStokesSkewedUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		std::stringstream ss; ss << "NavierStokesSkewedUpwind" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str())
			.add_constructor();
	}

//	NavierStokesLinearProfileSkewedUpwind
	{
		typedef NavierStokesLinearProfileSkewedUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		std::stringstream ss; ss << "NavierStokesLinearProfileSkewedUpwind" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str())
			.add_constructor();
	}

//	NavierStokesPositiveUpwind
	{
		typedef NavierStokesPositiveUpwind<dim> T;
		typedef INavierStokesUpwind<dim> TBase;
		std::stringstream ss; ss << "NavierStokesPositiveUpwind" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str())
			.add_constructor();
	}

/////////////////////////////////////////////////////////////////////////////
// Stabilization
/////////////////////////////////////////////////////////////////////////////


//	INavierStokesStabilization
	{
		typedef INavierStokesStabilization<dim> T;
		std::stringstream ss; ss << "INavierStokesStabilization" << dim << "d";
		reg.add_class_<T>(ss.str().c_str(), grp.c_str())
			.add_method("set_upwind", &T::set_upwind)
			.add_method("set_diffusion_length", &T::set_diffusion_length);
	}

//	NavierStokesFIELDSStabilization
	{
		typedef NavierStokesFIELDSStabilization<dim> T;
		typedef INavierStokesStabilization<dim> TBase;
		std::stringstream ss; ss << "NavierStokesFIELDSStabilization" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str())
			.add_constructor();
	}

/////////////////////////////////////////////////////////////////////////////
// Convection Shapes
/////////////////////////////////////////////////////////////////////////////

//	IConvectionShapes
	{
		typedef IConvectionShapes<dim> T;
		std::stringstream ss; ss << "IConvectionShapes" << dim << "d";
		reg.add_class_<T>(ss.str().c_str(), grp.c_str());
	}

//	ConvectionShapesNoUpwind
	{
		typedef ConvectionShapesNoUpwind<dim> T;
		typedef IConvectionShapes<dim> TBase;
		std::stringstream ss; ss << "NoUpwind" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str())
			.add_constructor();
	}

//	ConvectionShapesFullUpwind
	{
		typedef ConvectionShapesFullUpwind<dim> T;
		typedef IConvectionShapes<dim> TBase;
		std::stringstream ss; ss << "FullUpwind" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str())
			.add_constructor();
	}

//	ConvectionShapesWeightedUpwind
	{
		typedef ConvectionShapesWeightedUpwind<dim> T;
		typedef IConvectionShapes<dim> TBase;
		std::stringstream ss; ss << "WeightedUpwind" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str())
			.add_method("set_weight", &T::set_weight)
			.add_constructor();
	}

//	ConvectionShapesPartialUpwind
	{
		typedef ConvectionShapesPartialUpwind<dim> T;
		typedef IConvectionShapes<dim> TBase;
		std::stringstream ss; ss << "PartialUpwind" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str())
			.add_constructor();
	}

/////////////////////////////////////////////////////////////////////////////
// ElasticityTensorUserData
/////////////////////////////////////////////////////////////////////////////
	{
		typedef ElasticityTensorUserData<dim> T;
		typedef IUserData<MathTensor<4,dim>, dim> TBase;
		std::stringstream ss; ss << "ElasticityTensorUserData" << dim << "d";
		reg.add_class_<T, TBase >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("operator", &T::operator ())
			.add_method("test_evaluate", &T::test_evaluate);
			//methods, which are used in script-file
	}

}


bool RegisterLibDiscElemDisc(Registry& reg, const char* parentGroup)
{
	try
	{
	//	get group string
		std::string grp = parentGroup; grp.append("/Discretization");

#ifdef UG_DIM_1
	//	Domain dependend part 1D
			RegisterIElemDiscs<Domain1d>(reg, grp.c_str());
#endif

#ifdef UG_DIM_2
	//	Domain dependend part 2D
			RegisterIElemDiscs<Domain2d>(reg, grp.c_str());
#endif

#ifdef UG_DIM_3
	//	Domain dependend part 3D
			RegisterIElemDiscs<Domain3d>(reg, grp.c_str());
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
