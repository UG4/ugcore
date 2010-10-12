/*
 * lib_discretization_bridge.cpp
 *
 *  Created on: 22.09.2010
 *      Author: andreasvogel
 */

#include "../ug_bridge.h"
#include "lib_discretization/lib_discretization.h"
#include "problems/problems.h"

namespace ug
{
namespace bridge
{

template <typename TDomain>
bool LoadDomain(TDomain& domain, const char* filename)
{
	const char * p = strstr(filename, ".ugx");
	if(p == NULL)
	{
		UG_LOG("Currently only '.ugx' format supported for domains.\n");
		return false;
	}

	return LoadGridFromUGX(domain.get_grid(), domain.get_subset_handler(), filename);
}

template <typename TDomain>
bool SaveDomain(TDomain& domain, const char* filename)
{
	const char * p = strstr(filename, ".ugx");
	if(p == NULL)
	{
		UG_LOG("Currently only '.ugx' format supported for domains.\n");
		return false;
	}

	return SaveGridToUGX(domain.get_grid(), domain.get_subset_handler(), filename);
}

bool AddP1Function(P1ConformFunctionPattern& pattern, std::string name, int dim)
{
	return pattern.add_discrete_function(name, LSFS_LAGRANGEP1, dim);
}

template <typename TGridFunction>
bool WriteGridFunctionToVTK(TGridFunction& u, const char* filename)
{
	VTKOutput<TGridFunction> out;
	return out.print(filename, u);
}

template <typename TGridFunction>
bool SaveMatrixForConnectionViewer(	TGridFunction& u,
									IMatrixOperator<typename TGridFunction::vector_type,
													typename TGridFunction::vector_type,
													typename TGridFunction::algebra_type::matrix_type>& A,
									const char* filename)
{
	const char * p = strstr(filename, ".mat");
	if(p == NULL)
	{
		UG_LOG("Currently only '.mat' format supported for domains.\n");
		return false;
	}

	static const int dim = TGridFunction::domain_type::dim;

	vector<MathVector<dim> > positions;
	ExtractPositions(u, positions);

	WriteMatrixToConnectionViewer(filename, A.get_matrix(), &positions[0], dim);
	return true;
}

template <typename TGridFunction>
bool SaveVectorForConnectionViewer(	TGridFunction& b,
									const char* filename)
{
	const char * p = strstr(filename, ".mat");
	if(p == NULL)
	{
		UG_LOG("Currently only '.mat' format supported for domains.\n");
		return false;
	}

	static const int dim = TGridFunction::domain_type::dim;

	vector<MathVector<dim> > positions;
	ExtractPositions(b, positions);

	WriteVectorToConnectionViewer(filename, b.get_vector(), &positions[0], dim);
	return true;
}

void RegisterLibDiscretizationInterface(Registry& reg)
{

//	Domain2d
	{
	typedef Domain<2, MultiGrid, MGSubsetHandler> domain_type;
	reg.add_class_<domain_type>("Domain2d")
		.add_constructor()
		.add_method("get_subset_handler", (MGSubsetHandler& (domain_type::*)()) &domain_type::get_subset_handler)
		.add_method("get_grid", (MultiGrid& (domain_type::*)()) &domain_type::get_grid);

	reg.add_function("LoadDomain2d", &LoadDomain<domain_type>);
	reg.add_function("SaveDomain2d", &SaveDomain<domain_type>);
	}

//	Domain3d
	{
	typedef Domain<3, MultiGrid, MGSubsetHandler> domain_type;
	reg.add_class_<domain_type>("Domain3d")
		.add_constructor()
		.add_method("get_subset_handler", (MGSubsetHandler& (domain_type::*)()) &domain_type::get_subset_handler)
		.add_method("get_grid", (MultiGrid& (domain_type::*)()) &domain_type::get_grid);

	reg.add_function("LoadDomain3d", &LoadDomain<domain_type>);
	reg.add_function("SaveDomain3d", &SaveDomain<domain_type>);
	}

//	FunctionPattern (Abstract Base Class)
	reg.add_class_<FunctionPattern>("FunctionPattern");

//	P1ConformFunctionPattern
	{
	typedef P1ConformFunctionPattern T;
	reg.add_class_<T, FunctionPattern>("P1ConformFunctionPattern")
		.add_constructor()
		.add_method("lock", &T::lock);
	}

//	P1ConformDoFDistribution
	{
		typedef P1ConformDoFDistribution T;
		reg.add_class_<T>("P1ConformDoFDistribution");
	}

//  Add discrete function to pattern
	reg.add_function("AddP1Function", &AddP1Function);

//  ApproximationSpace2d
	{
		typedef Domain<2, MultiGrid, MGSubsetHandler> domain_type;
		typedef ApproximationSpace<domain_type, P1ConformDoFDistribution, MartinAlgebra> T;
		reg.add_class_<T>("ApproximationSpace2d")
			.add_constructor()
			.add_method("assign_domain", &T::assign_domain)
			.add_method("assign_function_pattern", &T::assign_function_pattern)
			.add_method("get_surface_dof_distribution", &T::get_surface_dof_distribution)
			.add_method("create_surface_function", &T::create_surface_function);
	}

//  ApproximationSpace3d
	{
		typedef Domain<3, MultiGrid, MGSubsetHandler> domain_type;
		typedef ApproximationSpace<domain_type, P1ConformDoFDistribution, MartinAlgebra> T;
		reg.add_class_<T>("ApproximationSpace3d")
			.add_constructor()
			.add_method("assign_domain", &T::assign_domain)
			.add_method("assign_function_pattern", &T::assign_function_pattern)
			.add_method("get_surface_dof_distribution", &T::get_surface_dof_distribution)
			.add_method("create_surface_function", &T::create_surface_function);
	}

//	DomainDiscretization
	{
		typedef P1ConformDoFDistribution dof_distribution_type;
		typedef MartinAlgebra algebra_type;
		typedef DomainDiscretization<dof_distribution_type, algebra_type> T;

		reg.add_class_<IAssemble<dof_distribution_type, algebra_type> >("IAssemble");
		reg.add_class_<IDomainDiscretization<dof_distribution_type, algebra_type>,
						IAssemble<dof_distribution_type, algebra_type> >("IDomainDiscretization");

		reg.add_class_<T, IDomainDiscretization<dof_distribution_type, algebra_type> >("DomainDiscretization")
			.add_constructor()
			.add_method("add_dirichlet_bnd", &T::add_dirichlet_bnd)
			.add_method("add", (bool (T::*)(IElemDisc<algebra_type>&, const FunctionGroup&, int, int)) &T::add);
	}

//	Time Discretization
	{
		typedef P1ConformDoFDistribution dof_distribution_type;
		typedef MartinAlgebra algebra_type;
		typedef ThetaTimeDiscretization<dof_distribution_type, algebra_type> T;

		reg.add_class_<	ITimeDiscretization<dof_distribution_type, algebra_type>,
						IAssemble<dof_distribution_type, algebra_type> >("ITimeDiscretization");

		reg.add_class_<T, ITimeDiscretization<dof_distribution_type, algebra_type> >("ThetaTimeDiscretization")
				.add_constructor()
				.add_method("set_domain_discretization", &T::set_domain_discretization)
				.add_method("set_theta", &T::set_theta);
	}

//	DirichletBNDValues
	{
		typedef Domain<2, MultiGrid, MGSubsetHandler> domain_type;
		typedef P1ConformDoFDistribution dof_distribution_type;
		typedef MartinAlgebra algebra_type;
		typedef DirichletBNDValues<domain_type, dof_distribution_type, algebra_type> T;

	//	Base class
		reg.add_class_<IDirichletBoundaryValues<dof_distribution_type, algebra_type> >("IDirichletBoundaryValues");

	//	derived implementation
		reg.add_class_<T, IDirichletBoundaryValues<dof_distribution_type, algebra_type> >("DirichletBND2d")
			.add_constructor()
			.add_method("set_domain", &T::set_domain)
			.add_method("set_dirichlet_function", &T::set_dirichlet_function)
			.add_method("set_function", &T::set_function);
	}

//	DirichletBoundaryFunction
	{
		reg.add_class_<DirichletBoundaryFunction<1> >("DirichletBoundaryFunction1d");
		reg.add_class_<DirichletBoundaryFunction<2> >("DirichletBoundaryFunction2d");
		reg.add_class_<DirichletBoundaryFunction<3> >("DirichletBoundaryFunction3d");
	}

//	UserFunction
	{
	//	Convection - Diffusion
		reg.add_class_<IConvDiffUserFunction<2> >("IConvDiffUserFunction2d");

	//	Density - Driven - Flow
		reg.add_class_<IDensityDrivenFlowUserFunction<2> >("IDensityDrivenFlowUserFunction2d");
	}

//	Elem Discs
	{
		typedef Domain<2, MultiGrid, MGSubsetHandler> domain_type;
		typedef MartinAlgebra algebra_type;

	//	Base class
		reg.add_class_<IElemDisc<algebra_type> >("IElemDisc");

		typedef FVConvectionDiffusionElemDisc<FV1Geometry, domain_type, algebra_type> T;
		reg.add_class_<T, IElemDisc<algebra_type> >("FV1ConvectionDiffusionElemDisc")
			.add_constructor()
			.add_method("set_domain", &T::set_domain)
			.add_method("set_user_functions", &T::set_user_functions)
			.add_method("set_upwind_amount", &T::set_upwind_amount);

		typedef DensityDrivenFlowElemDisc<domain_type, algebra_type> T2;
		reg.add_class_<T2, IElemDisc<algebra_type> >("FV1DensityDrivenFlowElemDisc")
			.add_constructor()
			.add_method("set_domain", &T2::set_domain)
			.add_method("set_user_functions", &T2::set_user_functions)
			.add_method("set_upwind_amount", &T2::set_upwind_amount);
	}

//	FunctionGroup
	{
		reg.add_class_<FunctionGroup>("FunctionGroup")
			.add_constructor()
			.add_method("clear", &FunctionGroup::clear)
			.add_method("add_function", &FunctionGroup::add_function);
	}

//	AssembledLinearOperator
	{
		typedef P1ConformDoFDistribution dof_distribution_type;
		typedef MartinAlgebra algebra_type;
		typedef AssembledLinearOperator<dof_distribution_type, algebra_type> T;

		reg.add_class_<T, IMatrixOperator<algebra_type::vector_type, algebra_type::vector_type, algebra_type::matrix_type> >
						("AssembledLinearOperator")
			.add_constructor()
			.add_method("init", (bool (T::*)())&T::init)
			.add_method("set_discretization", &T::set_discretization)
			.add_method("export_rhs", &T::export_rhs)
			.add_method("set_dof_distribution", &T::set_dof_distribution)
			.add_method("set_dirichlet_values", &T::set_dirichlet_values)
			.add_method("get_rhs", &T::get_rhs);
	}

//	AssembledOperator
	{
		typedef P1ConformDoFDistribution dof_distribution_type;
		typedef MartinAlgebra algebra_type;
		typedef AssembledOperator<dof_distribution_type, algebra_type> T;

		reg.add_class_<T, IOperator<algebra_type::vector_type, algebra_type::vector_type> >
						("AssembledOperator")
			.add_constructor()
			.add_method("set_discretization", &T::set_discretization)
			.add_method("set_dof_distribution", &T::set_dof_distribution)
			.add_method("init", &T::init);
	}

//	GridFunction
	{
		typedef Domain<2, MultiGrid, MGSubsetHandler> domain_type;
		typedef P1ConformDoFDistribution dof_distribution_type;
		typedef MartinAlgebra algebra_type;

#ifdef UG_PARALLEL
		typedef ParallelGridFunction<GridFunction<domain_type, dof_distribution_type, algebra_type> > T;
#else
		typedef GridFunction<domain_type, dof_distribution_type, algebra_type> T;
#endif
		reg.add_class_<T, algebra_type::vector_type>("GridFunction2d")
			.add_constructor()
			.add_method("set", (bool (T::*)(number))&T::set)
			.add_method("assign", (bool (T::*)(const algebra_type::vector_type&))&T::assign);
	}

//	ApplyLinearSolver
	{
		typedef Domain<2, MultiGrid, MGSubsetHandler> domain_type;
		typedef P1ConformDoFDistribution dof_distribution_type;
		typedef MartinAlgebra algebra_type;

#ifdef UG_PARALLEL
		typedef ParallelGridFunction<GridFunction<domain_type, dof_distribution_type, algebra_type> > function_type;
#else
		typedef GridFunction<domain_type, dof_distribution_type, algebra_type> function_type;
#endif

		reg.add_function("ApplyLinearSolver", &ApplyLinearSolver<function_type> );

		reg.add_function("WriteGridFunctionToVTK", &WriteGridFunctionToVTK<function_type>);

		reg.add_function("SaveMatrixForConnectionViewer", &SaveMatrixForConnectionViewer<function_type>);

		reg.add_function("SaveVectorForConnectionViewer", &SaveVectorForConnectionViewer<function_type>);
	}


//	AssembledMultiGridCycle
	{
		typedef Domain<2, MultiGrid, MGSubsetHandler> domain_type;
		typedef MartinAlgebra algebra_type;
		typedef ApproximationSpace<domain_type, P1ConformDoFDistribution, algebra_type> approximation_space_type;
		typedef AssembledMultiGridCycle<approximation_space_type, algebra_type> T;

		reg.add_class_<T, ILinearIterator<algebra_type::vector_type, algebra_type::vector_type> >("GeometricMultiGridPreconditioner2d")
			.add_constructor()
			.add_method("set_discretization", &T::set_discretization)
			.add_method("set_approximation_space", &T::set_approximation_space)
			.add_method("set_surface_level", &T::set_surface_level)
			.add_method("set_base_level", &T::set_base_level)
			.add_method("set_base_solver", &T::set_base_solver)
			.add_method("set_smoother", &T::set_smoother)
			.add_method("set_cycle_type", &T::set_cycle_type)
			.add_method("set_num_presmooth", &T::set_num_presmooth)
			.add_method("set_num_postsmooth", &T::set_num_postsmooth);
	}

//	StandardLineSearch
	{
		typedef MartinAlgebra algebra_type;
		typedef StandardLineSearch<algebra_type::vector_type> T;

		reg.add_class_<ILineSearch<algebra_type::vector_type> >("ILineSearch");

		reg.add_class_<StandardLineSearch<algebra_type::vector_type>,
						ILineSearch<algebra_type::vector_type> >("StandardLineSearch")
			.add_constructor()
			.add_method("set_maximum_steps", &T::set_maximum_steps)
			.add_method("set_lambda_start", &T::set_lambda_start)
			.add_method("set_reduce_factor", &T::set_reduce_factor)
			.add_method("set_verbose_level", &T::set_verbose_level)
			.add_method("set_offset", &T::set_offset);
	}

//	NewtonSolver
	{
		typedef P1ConformDoFDistribution dof_distribution_type;
		typedef MartinAlgebra algebra_type;
		typedef NewtonSolver<dof_distribution_type, algebra_type> T;

		reg.add_class_<T, IOperatorInverse<algebra_type::vector_type, algebra_type::vector_type> >("NewtonSolver")
			.add_constructor()
			.add_method("set_linear_solver", &T::set_linear_solver)
			.add_method("set_convergence_check", &T::set_convergence_check)
			.add_method("set_line_search", &T::set_line_search)
			.add_method("init", &T::init);

	}

//	VTK Output
	{
		typedef Domain<2, MultiGrid, MGSubsetHandler> domain_type;
		typedef P1ConformDoFDistribution dof_distribution_type;
		typedef MartinAlgebra algebra_type;

#ifdef UG_PARALLEL
		typedef ParallelGridFunction<GridFunction<domain_type, dof_distribution_type, algebra_type> > T;
#else
		typedef GridFunction<domain_type, dof_distribution_type, algebra_type> T;
#endif

		reg.add_class_<VTKOutput<T> >("VTKOutput")
			.add_constructor()
			.add_method("begin_timeseries", &VTKOutput<T>::begin_timeseries)
			.add_method("end_timeseries", &VTKOutput<T>::end_timeseries)
			.add_method("print", &VTKOutput<T>::print);

//	PerformTimeStep
		reg.add_function("PerformTimeStep", &PerformTimeStep<T>);
	}

//	Register user functions
	RegisterSinusUserFunctions(reg);
	RegisterElderUserFunctions(reg);
}


}//	end of namespace ug
}//	end of namespace interface
