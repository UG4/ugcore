/*
 * multigrid_bridge.cpp
 *
 *  Created on: 22.09.2010
 *      Author: andreasvogel
 */

// extern headers
#include <iostream>
#include <sstream>
#include <string>

// include bridge
#include "../bridge.h"
#include "registry/registry.h"

// lib_algebra includes
#include "lib_algebra/cpu_algebra_types.h"

// lib_disc includes
#include "lib_disc/domain.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/function_spaces/approximation_space.h"

#include "lib_disc/operator/linear_operator/projection_operator.h"
#include "lib_disc/operator/linear_operator/prolongation_operator.h"
#include "lib_disc/operator/linear_operator/multi_grid_solver/mg_solver.h"

using namespace std;

namespace ug {
namespace bridge {

template <typename TDomain, typename TAlgebra>
static void Register__Algebra_Domain(Registry& reg, string parentGroup)
{
//	typedef
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef ApproximationSpace<TDomain> approximation_space_type;

//	group string
	stringstream grpSS; grpSS << parentGroup << "/MultiGrid";
	string grp = grpSS.str();

//	suffix and tag
	string dimAlgSuffix = GetDomainSuffix<TDomain>();
	dimAlgSuffix.append(GetAlgebraSuffix<TAlgebra>());

	string dimAlgTag = GetDomainTag<TDomain>();
	dimAlgTag.append(GetAlgebraTag<TAlgebra>());


//	ProlongationOperator
	{
		typedef P1Prolongation<TDomain, TAlgebra> T;
		typedef IProlongationOperator<TAlgebra> TBase;
		string name = string("P1Prolongation").append(dimAlgSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.template add_constructor<void (*)(SmartPtr<approximation_space_type>)>("Approximation Space")
			.add_method("set_restriction_damping", &T::set_restriction_damping)
			.add_method("add_constraint", &T::add_constraint)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "P1Prolongation", dimAlgTag);
	}

//	ProjectionOperator
	{
		typedef P1Projection<TDomain, TAlgebra> T;
		typedef IProjectionOperator<vector_type, vector_type> TBase;
		string name = string("P1Projection").append(dimAlgSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.template add_constructor<void (*)(SmartPtr<approximation_space_type>)>("Approximation Space")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "P1Projection", dimAlgTag);
	}

//	AssembledMultiGridCycle
	{
		typedef AssembledMultiGridCycle<TDomain, TAlgebra> T;
		typedef ILinearIterator<vector_type> TBase;
		string name = string("GeometricMultiGrid").append(dimAlgSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >)>("Approximation Space")
			.add_method("set_discretization", &T::set_discretization, "", "Discretization")
			.add_method("set_base_level", &T::set_base_level, "", "Base Level")
			.add_method("set_parallel_base_solver", &T::set_parallel_base_solver,"", "Specifies if base solver works in parallel")
			.add_method("set_base_solver", &T::set_base_solver,"","Base Solver")
			.add_method("set_smoother", &T::set_smoother,"", "Smoother")
			.add_method("set_cycle_type", &T::set_cycle_type,"", "Cycle Type")
			.add_method("set_num_presmooth", &T::set_num_presmooth,"", "Number PreSmooth Steps")
			.add_method("set_num_postsmooth", &T::set_num_postsmooth,"", "Number PostSmooth Steps")
			.add_method("set_prolongation", &T::set_prolongation_operator,"", "Prolongation")
			.add_method("set_projection", &T::set_projection_operator,"", "Projection")
			.add_method("set_debug", &T::set_debug)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GeometricMultiGrid", dimAlgTag);
	}
}

template <typename TAlgebra>
static void Register__Algebra(Registry& reg, string parentGroup)
{
//	get group string
	string grp = parentGroup; grp.append("/Discretization");

	try
	{

#ifdef UG_DIM_1
		Register__Algebra_Domain<Domain1d, TAlgebra>(reg, grp);
#endif
#ifdef UG_DIM_2
		Register__Algebra_Domain<Domain2d, TAlgebra>(reg, grp);
#endif
#ifdef UG_DIM_3
		Register__Algebra_Domain<Domain3d, TAlgebra>(reg, grp);
#endif


		typedef typename TAlgebra::vector_type vector_type;
	//	suffix and tag
		string algSuffix = GetAlgebraSuffix<TAlgebra>();
		string algTag = GetAlgebraTag<TAlgebra>();

	//	IProlongationOperator
		{
			typedef IProlongationOperator<TAlgebra> T;
			string name = string("IProlongationOperator").append(algSuffix);
			reg.add_class_<T>(name, grp);
			reg.add_class_to_group(name, "IProlongationOperator", algTag);
		}

	//	IProjectionOperator
		{
			typedef IProjectionOperator<vector_type, vector_type> T;
			string name = string("IProjectionOperator").append(algSuffix);
			reg.add_class_<T>(name, grp);
			reg.add_class_to_group(name, "IProjectionOperator", algTag);
		}

	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in Register__Algebra_DoFDistribution: "
				"Registration failed (using name " << ex.name << ").\n");
		UG_THROW("Registration failed.");
	}
}

void RegisterBridge_MultiGrid(Registry& reg, string grp)
{
#ifdef UG_CPU_1
	Register__Algebra<CPUAlgebra>(reg, grp);
#endif
#ifdef UG_CPU_2
	Register__Algebra<CPUBlockAlgebra<2> >(reg, grp);
#endif
#ifdef UG_CPU_3
	Register__Algebra<CPUBlockAlgebra<3> >(reg, grp);
#endif
#ifdef UG_CPU_4
	Register__Algebra<CPUBlockAlgebra<4> >(reg, grp);
#endif
#ifdef UG_CPU_VAR
	Register__Algebra<CPUVariableBlockAlgebra >(reg, grp);
#endif
}

}//	end of namespace ug
}//	end of namespace interface
