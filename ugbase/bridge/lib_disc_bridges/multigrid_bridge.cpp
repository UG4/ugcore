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
#include "lib_disc/dof_manager/conform/conform.h"
#include "lib_disc/dof_manager/p1conform/p1conform.h"

#include "lib_disc/operator/linear_operator/projection_operator.h"
#include "lib_disc/operator/linear_operator/prolongation_operator.h"
#include "lib_disc/operator/linear_operator/multi_grid_solver/mg_solver.h"

using namespace std;

namespace ug {
namespace bridge {

template <typename TDomain, typename TAlgebra, typename TDoFDistribution>
static void Register__Algebra_DoFDistribution_Domain(Registry& reg, string parentGroup)
{
//	typedef
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef ApproximationSpace<TDomain, TDoFDistribution, TAlgebra> approximation_space_type;

#ifdef UG_PARALLEL
		typedef ParallelGridFunction<GridFunction<TDomain, TDoFDistribution, TAlgebra> > function_type;
#else
		typedef GridFunction<TDomain, TDoFDistribution, TAlgebra> function_type;
#endif

//	group string
	stringstream grpSS; grpSS << parentGroup << "/MultiGrid";
	string grp = grpSS.str();

//	suffix and tag
	string dimAlgDDSuffix = GetDomainSuffix<TDomain>();
	dimAlgDDSuffix.append(GetAlgebraSuffix<TAlgebra>());
	dimAlgDDSuffix.append(GetDoFDistributionSuffix<TDoFDistribution>());

	string dimAlgDDTag = GetDomainTag<TDomain>();
	dimAlgDDTag.append(GetAlgebraTag<TAlgebra>());
	dimAlgDDTag.append(GetDoFDistributionTag<TDoFDistribution>());


//	ProlongationOperator
	{
		typedef P1Prolongation<approximation_space_type, TAlgebra> T;
		typedef IProlongationOperator<vector_type, vector_type> TBase;
		string name = string("P1Prolongation").append(dimAlgDDSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.template add_constructor<void (*)(approximation_space_type&)>("Approximation Space")
			.add_method("set_restriction_damping", &T::set_restriction_damping)
			.add_method("set_dirichlet_post_process", &T::set_dirichlet_post_process);
		reg.add_class_to_group(name, "P1Prolongation", dimAlgDDTag);
	}

//	ProjectionOperator
	{
		typedef P1Projection<approximation_space_type, TAlgebra> T;
		typedef IProjectionOperator<vector_type, vector_type> TBase;
		string name = string("P1Projection").append(dimAlgDDSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.template add_constructor<void (*)(approximation_space_type&)>("Approximation Space");
		reg.add_class_to_group(name, "P1Projection", dimAlgDDTag);
	}

//	AssembledMultiGridCycle
	{
		typedef AssembledMultiGridCycle<approximation_space_type, TAlgebra> T;
		typedef ILinearIterator<vector_type, vector_type> TBase;
		string name = string("GeometricMultiGrid").append(dimAlgDDSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(approximation_space_type&)>("Approximation Space")
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
			.add_method("set_debug", &T::set_debug);
		reg.add_class_to_group(name, "GeometricMultiGrid", dimAlgDDTag);
	}
}

template <typename TAlgebra, typename TDoFDistribution>
static bool Register__Algebra_DoFDistribution(Registry& reg, string parentGroup)
{
//	get group string
	string grp = parentGroup; grp.append("/Discretization");

	try
	{

#ifdef UG_DIM_1
//	Domain dependent part 1D
	{
		typedef Domain<1, MultiGrid, MGSubsetHandler> domain_type;
		Register__Algebra_DoFDistribution_Domain<domain_type, TAlgebra, TDoFDistribution>(reg, grp);
	}
#endif

#ifdef UG_DIM_2
//	Domain dependent part 2D
	{
		typedef Domain<2, MultiGrid, MGSubsetHandler> domain_type;
		Register__Algebra_DoFDistribution_Domain<domain_type, TAlgebra, TDoFDistribution>(reg, grp);
	}
#endif

#ifdef UG_DIM_3
//	Domain dependent part 3D
	{
		typedef Domain<3, MultiGrid, MGSubsetHandler> domain_type;
		Register__Algebra_DoFDistribution_Domain<domain_type, TAlgebra, TDoFDistribution>(reg, grp);
	}
#endif

	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in Register__Algebra_DoFDistribution: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}

	return true;
}

template <typename TAlgebra>
static bool Register__Algebra(Registry& reg, string parentGroup)
{
	bool bReturn = true;
#ifdef DOF_P1
	bReturn &= Register__Algebra_DoFDistribution<TAlgebra, P1DoFDistribution>(reg, parentGroup);
#endif
#ifdef DOF_GEN
	bReturn &= Register__Algebra_DoFDistribution<TAlgebra, DoFDistribution >(reg, parentGroup);
#endif

	return bReturn;
}

bool RegisterMultiGrid(Registry& reg, string parentGroup)
{
	bool bReturn = true;
#ifdef UG_CPU_1
	bReturn &= Register__Algebra<CPUAlgebra>(*reg, grp);
#endif
#ifdef UG_CPU_2
	bReturn &= Register__Algebra<CPUBlockAlgebra<2> >(*reg, grp);
#endif
#ifdef UG_CPU_3
	bReturn &= Register__Algebra<CPUBlockAlgebra<3> >(*reg, grp);
#endif
#ifdef UG_CPU_4
	bReturn &= Register__Algebra<CPUBlockAlgebra<4> >(*reg, grp);
#endif
#ifdef UG_CPU_VAR
	bReturn &= Register__Algebra<CPUVariableBlockAlgebra >(*reg, grp);
#endif
	return bReturn;
}

}//	end of namespace ug
}//	end of namespace interface
