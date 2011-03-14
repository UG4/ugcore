/*
 * lib_disc_bridge_discs.cpp
 *
 *  Created on: 22.09.2010
 *      Author: andreasvogel
 */

// extern headers
#include <iostream>
#include <sstream>

// include bridge
#include "../ug_bridge.h"
#include "../registry.h"

// lib_algebra includes
#include "lib_algebra/algebra_chooser.h"
#include "lib_algebra/operator/operator_util.h"
#include "lib_algebra/operator/operator_interface.h"
#include "lib_algebra/operator/operator_inverse_interface.h"

// lib_disc includes
#include "lib_discretization/domain.h"
#include "lib_discretization/function_spaces/grid_function.h"
#include "lib_discretization/function_spaces/grid_function_space.h"
#include "lib_discretization/dof_manager/p1conform/p1conform.h"

#include "lib_discretization/spatial_discretization/elem_disc/navier_stokes/fv/navier_stokes.h"
#include "lib_discretization/spatial_discretization/elem_disc/navier_stokes/fv/upwind.h"
#include "lib_discretization/spatial_discretization/elem_disc/navier_stokes/fv/stabilization.h"

namespace ug
{
extern enum_AlgebraType g_AlgebraType;

namespace bridge
{

template <typename TDomain, typename TAlgebra, typename TDoFDistribution>
void RegisterLibDiscDiscObjects(Registry& reg, const char* parentGroup)
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

//	Navier-Stokes
	{
		typedef FVNavierStokesElemDisc<domain_type, algebra_type> T;
		typedef IDomainElemDisc<domain_type, algebra_type> TBase;
		std::stringstream ss; ss << "FV1NavierStokes" << dim << "d";
		reg.add_class_<T, TBase >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set_kinematic_viscosity", &T::set_kinematic_viscosity)
			.add_method("set_stabilization", &T::set_stabilization)
			.add_method("set_conv_upwind",  (void (T::*)(INavierStokesStabilization<dim, algebra_type>&))&T::set_conv_upwind)
			.add_method("set_conv_upwind",  (void (T::*)(INavierStokesUpwind<dim, algebra_type>&))&T::set_conv_upwind)
			.add_method("set_peclet_blend", &T::set_peclet_blend)
			.add_method("set_exact_jacobian", &T::set_exact_jacobian);
	}

/////////////////////////////////////////////////////////////////////////////
// Upwind
/////////////////////////////////////////////////////////////////////////////

//	INavierStokesUpwind
	{
		typedef INavierStokesUpwind<dim, algebra_type> T;
		std::stringstream ss; ss << "INavierStokesUpwind" << dim << "d";
		reg.add_class_<T>(ss.str().c_str(), grp.c_str());
	}

//	NavierStokesNoUpwind
	{
		typedef NavierStokesNoUpwind<dim, algebra_type> T;
		typedef INavierStokesUpwind<dim, algebra_type> TBase;
		std::stringstream ss; ss << "NavierStokesNoUpwind" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str())
			.add_constructor();
	}

//	NavierStokesFullUpwind
	{
		typedef NavierStokesFullUpwind<dim, algebra_type> T;
		typedef INavierStokesUpwind<dim, algebra_type> TBase;
		std::stringstream ss; ss << "NavierStokesFullUpwind" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str())
			.add_constructor();
	}

//	NavierStokesSkewedUpwind
	{
		typedef NavierStokesSkewedUpwind<dim, algebra_type> T;
		typedef INavierStokesUpwind<dim, algebra_type> TBase;
		std::stringstream ss; ss << "NavierStokesSkewedUpwind" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str())
			.add_constructor();
	}

//	NavierStokesLinearProfileSkewedUpwind
	{
		typedef NavierStokesLinearProfileSkewedUpwind<dim, algebra_type> T;
		typedef INavierStokesUpwind<dim, algebra_type> TBase;
		std::stringstream ss; ss << "NavierStokesLinearProfileSkewedUpwind" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str())
			.add_constructor();
	}

//	NavierStokesPositiveUpwind
	{
		typedef NavierStokesPositiveUpwind<dim, algebra_type> T;
		typedef INavierStokesUpwind<dim, algebra_type> TBase;
		std::stringstream ss; ss << "NavierStokesPositiveUpwind" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str())
			.add_constructor();
	}

/////////////////////////////////////////////////////////////////////////////
// Stabilization
/////////////////////////////////////////////////////////////////////////////


//	INavierStokesStabilization
	{
		typedef INavierStokesStabilization<dim, algebra_type> T;
		std::stringstream ss; ss << "INavierStokesStabilization" << dim << "d";
		reg.add_class_<T>(ss.str().c_str(), grp.c_str())
			.add_method("set_upwind", &T::set_upwind)
			.add_method("set_diffusion_length", &T::set_diffusion_length);
	}

//	NavierStokesFIELDSStabilization
	{
		typedef NavierStokesFIELDSStabilization<dim, algebra_type> T;
		typedef INavierStokesStabilization<dim, algebra_type> TBase;
		std::stringstream ss; ss << "NavierStokesFIELDSStabilization" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str())
			.add_constructor();
	}

}

template <typename TDomain, typename TAlgebra, typename TDoFDistribution>
void RegisterLibDiscDiscFunctions(Registry& reg, const char* parentGroup)
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

}

template <typename TAlgebra, typename TDoFDistribution>
bool RegisterLibDiscInterfaceDiscs(Registry& reg, const char* parentGroup)
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
			RegisterLibDiscDiscObjects<domain_type, algebra_type, dof_distribution_type>(reg, grp.c_str());
			RegisterLibDiscDiscFunctions<domain_type,  algebra_type, dof_distribution_type>(reg, grp.c_str());
		}
#endif

#ifdef UG_DIM_2
	//	Domain dependend part 2D
		{
			typedef Domain<2, MultiGrid, MGSubsetHandler> domain_type;
			RegisterLibDiscDiscObjects<domain_type, algebra_type, dof_distribution_type>(reg, grp.c_str());
			RegisterLibDiscDiscFunctions<domain_type,  algebra_type, dof_distribution_type>(reg, grp.c_str());
		}
#endif

#ifdef UG_DIM_3
	//	Domain dependend part 3D
		{
			typedef Domain<3, MultiGrid, MGSubsetHandler> domain_type;
			RegisterLibDiscDiscObjects<domain_type, algebra_type, dof_distribution_type>(reg, grp.c_str());
			RegisterLibDiscDiscFunctions<domain_type,  algebra_type, dof_distribution_type>(reg, grp.c_str());
		}
#endif
	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterLibDiscretizationInterfaceForAlgebra: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}

	return true;
}

bool RegisterDynamicLibDiscInterfaceDiscs(Registry& reg, int algebra_type, const char* parentGroup)
{
	bool bReturn = true;

	switch(algebra_type)
	{
	case eCPUAlgebra:		 		bReturn &= RegisterLibDiscInterfaceDiscs<CPUAlgebra, P1ConformDoFDistribution>(reg, parentGroup); break;
//	case eCPUBlockAlgebra2x2: 		bReturn &= RegisterLibDiscInterfaceDiscs<CPUBlockAlgebra<2>, GroupedP1ConformDoFDistribution>(reg, parentGroup); break;
	case eCPUBlockAlgebra3x3: 		bReturn &= RegisterLibDiscInterfaceDiscs<CPUBlockAlgebra<3>, GroupedP1ConformDoFDistribution>(reg, parentGroup); break;
//	case eCPUBlockAlgebra4x4: 		bReturn &= RegisterLibDiscInterfaceDiscs<CPUBlockAlgebra<4>, GroupedP1ConformDoFDistribution>(reg, parentGroup); break;
//	case eCPUVariableBlockAlgebra: 	bReturn &= RegisterLibDiscInterfaceDiscs<CPUVariableBlockAlgebra, GroupedP1ConformDoFDistribution>(reg, parentGroup); break;
	default: UG_ASSERT(0, "Unsupported Algebra Type");
				UG_LOG("Unsupported Algebra Type requested.\n");
				return false;
	}

	return bReturn;
}

}//	end of namespace ug
}//	end of namespace interface
