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
#include "../../../ug_bridge.h"
#include "../../../registry.h"

// lib_algebra includes
#include "lib_algebra/algebra_selector.h"
#include "lib_algebra/algebra_types.h"
#include "lib_algebra/operator/operator_util.h"
#include "lib_algebra/operator/operator_interface.h"
#include "lib_algebra/operator/operator_inverse_interface.h"

// lib_disc includes
#include "lib_discretization/domain.h"
#include "lib_discretization/function_spaces/grid_function.h"
#include "lib_discretization/function_spaces/approximation_space.h"
#include "lib_discretization/dof_manager/p1conform/p1conform.h"

#include "lib_discretization/spatial_discretization/elem_disc/thermohaline_flow/fv1/thermohaline_flow.h"
#include "lib_discretization/spatial_discretization/disc_helper/conv_shape_interface.h"
#include "lib_discretization/spatial_discretization/disc_helper/conv_shape.h"


namespace ug
{

namespace bridge
{

template <typename TDomain, typename TAlgebra, typename TDoFDistribution>
void RegisterThermohalineFlowObjects(Registry& reg, const char* parentGroup)
{
//	typedef domain
	typedef TDomain domain_type;
	typedef TDoFDistribution dof_distribution_type;
	typedef TAlgebra algebra_type;
	typedef typename algebra_type::vector_type vector_type;
	typedef typename algebra_type::matrix_type matrix_type;
	static const int dim = domain_type::dim;

	std::stringstream grpSS; grpSS << parentGroup << "/" << dim << "d";
	std::string grp = grpSS.str();

#ifdef UG_PARALLEL
		typedef ParallelGridFunction<GridFunction<domain_type, dof_distribution_type, algebra_type> > function_type;
#else
		typedef GridFunction<domain_type, dof_distribution_type, algebra_type> function_type;
#endif

/////////////////////////////////////////////////////////////////////////////
// Elem Disc
/////////////////////////////////////////////////////////////////////////////

//	Density Driven Flow
	{
		typedef ThermohalineFlowElemDisc<domain_type, algebra_type> T2;
		typedef IDomainElemDisc<domain_type, algebra_type> TBase;
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
}


template <typename TAlgebra, typename TDoFDistribution>
bool RegisterThermohalineFlowDisc(Registry& reg, const char* parentGroup)
{
	typedef TDoFDistribution dof_distribution_type;
	typedef TAlgebra algebra_type;
	typedef typename algebra_type::vector_type vector_type;
	typedef typename algebra_type::matrix_type matrix_type;

	try
	{
	//	get group string
		std::string grp = parentGroup; grp.append("/Discretization");

#ifdef UG_DIM_1
	//	Domain dependend part 1D
		{
			typedef Domain<1, MultiGrid, MGSubsetHandler> domain_type;
			RegisterThermohalineFlowObjects<domain_type, algebra_type, dof_distribution_type>(reg, grp.c_str());
		}
#endif

#ifdef UG_DIM_2
	//	Domain dependend part 2D
		{
			typedef Domain<2, MultiGrid, MGSubsetHandler> domain_type;
			RegisterThermohalineFlowObjects<domain_type, algebra_type, dof_distribution_type>(reg, grp.c_str());
		}
#endif

#ifdef UG_DIM_3
	//	Domain dependend part 3D
		{
			typedef Domain<3, MultiGrid, MGSubsetHandler> domain_type;
			RegisterThermohalineFlowObjects<domain_type, algebra_type, dof_distribution_type>(reg, grp.c_str());
		}
#endif
	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterThermohalineFlowDiscs: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}

	return true;
}

bool RegisterDynamicThermohalineFlowDisc(Registry& reg, int algebra_type, const char* parentGroup)
{
	bool bReturn = true;

	switch(algebra_type)
	{
	case eCPUAlgebra:		 		bReturn &= RegisterThermohalineFlowDisc<CPUAlgebra, P1ConformDoFDistribution>(reg, parentGroup); break;
//	case eCPUBlockAlgebra2x2: 		bReturn &= RegisterThermohalineFlowDisc<CPUBlockAlgebra<2>, GroupedP1ConformDoFDistribution>(reg, parentGroup); break;
	case eCPUBlockAlgebra3x3: 		bReturn &= RegisterThermohalineFlowDisc<CPUBlockAlgebra<3>, GroupedP1ConformDoFDistribution>(reg, parentGroup); break;
//	case eCPUBlockAlgebra4x4: 		bReturn &= RegisterThermohalineFlowDisc<CPUBlockAlgebra<4>, GroupedP1ConformDoFDistribution>(reg, parentGroup); break;
//	case eCPUVariableBlockAlgebra: 	bReturn &= RegisterThermohalineFlowDisc<CPUVariableBlockAlgebra, GroupedP1ConformDoFDistribution>(reg, parentGroup); break;
	default: UG_ASSERT(0, "Unsupported Algebra Type");
				UG_LOG("Unsupported Algebra Type requested.\n");
				return false;
	}

	return bReturn;
}

}//	end of namespace ug
}//	end of namespace interface
