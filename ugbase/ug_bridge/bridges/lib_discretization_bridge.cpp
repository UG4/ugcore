/*
 * lib_discretization_bridge.cpp
 *
 *  Created on: 22.09.2010
 *      Author: andreasvogel
 */

#include "../ug_bridge.h"
#include "lib_discretization/lib_discretization.h"
#include "user_data/user_data.h"

namespace ug
{
namespace bridge
{

//////////////////////////////////
// Some global functions
//////////////////////////////////

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

bool AddP1Function(P1ConformFunctionPattern& pattern, const char* name, int dim)
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

template <typename TDomain>
void RegisterLibDiscretizationDomainDepended(Registry& reg)
{
//	typedef domain
	typedef TDomain domain_type;
	static const int dim = domain_type::dim;

//	todo: generalize
	typedef P1ConformDoFDistribution dof_distribution_type;
	typedef MartinAlgebra algebra_type;
#ifdef UG_PARALLEL
		typedef ParallelGridFunction<GridFunction<domain_type, dof_distribution_type, algebra_type> > function_type;
#else
		typedef GridFunction<domain_type, dof_distribution_type, algebra_type> function_type;
#endif

//	GridFunction
	{
		stringstream ss; ss << "GridFunction" << dim << "d";
		reg.add_class_<function_type, algebra_type::vector_type>(ss.str().c_str())
			.add_constructor()
			.add_method("set", (bool (function_type::*)(number))&function_type::set)
			.add_method("assign", (bool (function_type::*)(const algebra_type::vector_type&))&function_type::assign)
			.add_method("get_dim", &function_type::get_dim);
	}

//	Domain
	{
		stringstream ss; ss << "Domain" << dim << "d";
		reg.add_class_<domain_type>(ss.str().c_str())
			.add_constructor()
			.add_method("get_subset_handler", (MGSubsetHandler& (domain_type::*)()) &domain_type::get_subset_handler)
			.add_method("get_grid", (MultiGrid& (domain_type::*)()) &domain_type::get_grid)
			.add_method("get_dim", (int (domain_type::*)()) &domain_type::get_dim);
	}

// 	LoadDomain
	{
		stringstream ss; ss << "LoadDomain" << dim << "d";
		reg.add_function(ss.str().c_str(), &LoadDomain<domain_type>);
	}

//	SaveDomain
	{
		stringstream ss; ss << "SaveDomain" << dim << "d";
		reg.add_function(ss.str().c_str(), &SaveDomain<domain_type>);
	}

//  ApproximationSpace
	{
		typedef ApproximationSpace<domain_type, dof_distribution_type, algebra_type> T;
		stringstream ss; ss << "ApproximationSpace" << dim << "d";
		reg.add_class_<T>(ss.str().c_str())
			.add_constructor()
			.add_method("assign_domain", &T::assign_domain)
			.add_method("get_domain", (domain_type& (T::*)())&T::get_domain)
			.add_method("assign_function_pattern", &T::assign_function_pattern)
			.add_method("get_surface_dof_distribution", &T::get_surface_dof_distribution)
			.add_method("create_surface_function", &T::create_surface_function);
	}

//	DirichletBNDValues
	{
		typedef DirichletBNDValues<domain_type, dof_distribution_type, algebra_type> T;
		stringstream ss; ss << "DirichletBND" << dim << "d";
		reg.add_class_<T, IPostProcess<dof_distribution_type, algebra_type> >(ss.str().c_str())
			.add_constructor()
			.add_method("set_domain", &T::set_domain)
			.add_method("set_pattern", &T::set_pattern)
			.add_method("add_boundary_value", (bool (T::*)(IBoundaryNumberProvider<dim>&, const char*, const char*))&T::add_boundary_value);
	}

//	Neumann Boundary
	{
		typedef FVNeumannBoundaryElemDisc<FV1Geometry, domain_type, algebra_type> T;
		stringstream ss; ss << "FV1NeumannBoundaryElemDisc" << dim << "d";
		reg.add_class_<T, IElemDisc<algebra_type> >(ss.str().c_str())
			.add_constructor()
			.add_method("set_domain", &T::set_domain)
			.add_method("set_bnd_cond", &T::set_bnd_cond);
	}

//	Convection Diffusion
	{
		typedef FVConvectionDiffusionElemDisc<FV1Geometry, domain_type, algebra_type> T;
		stringstream ss; ss << "FV1ConvectionDiffusionElemDisc" << dim << "d";
		reg.add_class_<T, IElemDisc<algebra_type> >(ss.str().c_str())
			.add_constructor()
			.add_method("set_domain", &T::set_domain)
			.add_method("set_diffusion_tensor", &T::set_diffusion_tensor)
			.add_method("set_velocity_field", &T::set_velocity_field)
			.add_method("set_reaction", &T::set_reaction)
			.add_method("set_rhs", &T::set_rhs)
			.add_method("set_upwind_amount", &T::set_upwind_amount);
	}

//	Density Driven Flow
	{
		typedef DensityDrivenFlowElemDisc<domain_type, algebra_type> T2;
		stringstream ss; ss << "FV1DensityDrivenFlowElemDisc" << dim << "d";
		reg.add_class_<T2, IElemDisc<algebra_type> >(ss.str().c_str())
			.add_constructor()
			.add_method("set_domain", &T2::set_domain)
			.add_method("set_user_functions", &T2::set_user_functions)
			.add_method("set_upwind_amount", &T2::set_upwind_amount);
	}

//	ApplyLinearSolver
	{
		stringstream ss; ss << "ApplyLinearSolver" << dim << "d";
		reg.add_function(ss.str().c_str(), &ApplyLinearSolver<function_type> );
	}

//	WriteGridToVTK
	{
		stringstream ss; ss << "WriteGridFunctionToVTK" << dim << "d";
		reg.add_function(ss.str().c_str(), &WriteGridFunctionToVTK<function_type>);
	}

//	SaveMatrixForConnectionViewer
	{
		stringstream ss; ss << "SaveMatrixForConnectionViewer" << dim << "d";
		reg.add_function(ss.str().c_str(), &SaveMatrixForConnectionViewer<function_type>);
	}

//	SaveVectorForConnectionViewer
	{
		stringstream ss; ss << "SaveVectorForConnectionViewer" << dim << "d";
		reg.add_function(ss.str().c_str(), &SaveVectorForConnectionViewer<function_type>);
	}

//	ProlongationOperator
	{
		typedef ApproximationSpace<domain_type, dof_distribution_type, algebra_type> approximation_space_type;
		typedef P1ProlongationOperator<approximation_space_type, algebra_type> T;

		stringstream ss; ss << "P1ProlongationOperator" << dim << "d";
		reg.add_class_<T, IProlongationOperator<algebra_type::vector_type, algebra_type::vector_type> >(ss.str().c_str())
			.add_constructor()
			.add_method("set_approximation_space", &T::set_approximation_space)
			.add_method("set_dirichlet_post_process", &T::set_dirichlet_post_process);

	}

//	ProjectionOperator
	{
		typedef ApproximationSpace<domain_type, dof_distribution_type, algebra_type> approximation_space_type;
		typedef P1ProjectionOperator<approximation_space_type, algebra_type> T;

		stringstream ss; ss << "P1ProjectionOperator" << dim << "d";
		reg.add_class_<T, IProjectionOperator<algebra_type::vector_type, algebra_type::vector_type> >(ss.str().c_str())
			.add_constructor()
			.add_method("set_approximation_space", &T::set_approximation_space);
	}

//	AssembledMultiGridCycle
	{
		typedef ApproximationSpace<domain_type, dof_distribution_type, algebra_type> approximation_space_type;
		typedef AssembledMultiGridCycle<approximation_space_type, algebra_type> T;

		stringstream ss; ss << "GeometricMultiGridPreconditioner" << dim << "d";
		reg.add_class_<T, ILinearIterator<algebra_type::vector_type, algebra_type::vector_type> >(ss.str().c_str())
			.add_constructor()
			.add_method("set_discretization", &T::set_discretization)
			.add_method("set_approximation_space", &T::set_approximation_space)
			.add_method("set_surface_level", &T::set_surface_level)
			.add_method("set_base_level", &T::set_base_level)
			.add_method("set_base_solver", &T::set_base_solver)
			.add_method("set_smoother", &T::set_smoother)
			.add_method("set_cycle_type", &T::set_cycle_type)
			.add_method("set_num_presmooth", &T::set_num_presmooth)
			.add_method("set_num_postsmooth", &T::set_num_postsmooth)
			.add_method("set_prolongation", &T::set_prolongation_operator)
			.add_method("set_projection", &T::set_projection_operator);
	}

//	VTK Output
	{
		stringstream ss; ss << "VTKOutput" << dim << "d";
		reg.add_class_<VTKOutput<function_type> >(ss.str().c_str())
			.add_constructor()
			.add_method("begin_timeseries", &VTKOutput<function_type>::begin_timeseries)
			.add_method("end_timeseries", &VTKOutput<function_type>::end_timeseries)
			.add_method("print", &VTKOutput<function_type>::print);
	}

//	PerformTimeStep
	{
		stringstream ss; ss << "PerformTimeStep" << dim << "d";
		reg.add_function(ss.str().c_str(), &PerformTimeStep<function_type>);
	}

}


void RegisterLibDiscretizationInterface(Registry& reg)
{
	//	todo: generalize
		typedef P1ConformDoFDistribution dof_distribution_type;
		typedef MartinAlgebra algebra_type;

	//	FunctionPattern (Abstract Base Class)
		reg.add_class_<FunctionPattern>("FunctionPattern");

	//	P1ConformFunctionPattern
		{
		typedef P1ConformFunctionPattern T;
		reg.add_class_<T, FunctionPattern>("P1ConformFunctionPattern")
			.add_constructor()
			.add_method("set_subset_handler", &T::set_subset_handler)
			.add_method("lock", &T::lock);
		}

	//	P1ConformDoFDistribution
		{
			typedef P1ConformDoFDistribution T;
			reg.add_class_<T>("P1ConformDoFDistribution");
		}

	//  Add discrete function to pattern
		reg.add_function("AddP1Function", &AddP1Function);

	//	DomainDiscretization
		{
			typedef DomainDiscretization<dof_distribution_type, algebra_type> T;

			reg.add_class_<IAssemble<dof_distribution_type, algebra_type> >("IAssemble");
			reg.add_class_<IDomainDiscretization<dof_distribution_type, algebra_type>,
							IAssemble<dof_distribution_type, algebra_type> >("IDomainDiscretization");

			reg.add_class_<T, IDomainDiscretization<dof_distribution_type, algebra_type> >("DomainDiscretization")
				.add_constructor()
				.add_method("add_dirichlet_bnd", &T::add_dirichlet_bnd)
				.add_method("add", (bool (T::*)(IElemDisc<algebra_type>&, const FunctionPattern&, const char*, const char*)) &T::add);
		}

	//	Time Discretization
		{
			typedef ThetaTimeDiscretization<dof_distribution_type, algebra_type> T;

			reg.add_class_<	ITimeDiscretization<dof_distribution_type, algebra_type>,
							IAssemble<dof_distribution_type, algebra_type> >("ITimeDiscretization");

			reg.add_class_<T, ITimeDiscretization<dof_distribution_type, algebra_type> >("ThetaTimeDiscretization")
					.add_constructor()
					.add_method("set_domain_discretization", &T::set_domain_discretization)
					.add_method("set_theta", &T::set_theta);
		}

	//	Base class
		reg.add_class_<IPostProcess<dof_distribution_type, algebra_type> >("IPostProcess");

	//	UserFunction
		{
		//	Density - Driven - Flow
			reg.add_class_<IDensityDrivenFlowUserFunction<2> >("IDensityDrivenFlowUserFunction2d");
		}

	//	Elem Discs
		{
		//	Base class
			reg.add_class_<IElemDisc<algebra_type> >("IElemDisc");
		}

	//	FunctionGroup
		{
			reg.add_class_<FunctionGroup>("FunctionGroup")
				.add_constructor()
				.add_method("clear", &FunctionGroup::clear)
				.add_method("set_function_pattern", &FunctionGroup::set_function_pattern)
				.add_method("add_function", (bool (FunctionGroup::*)(const char*))&FunctionGroup::add_function);
		}

	//	AssembledLinearOperator
		{
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
			typedef AssembledOperator<dof_distribution_type, algebra_type> T;

			reg.add_class_<T, IOperator<algebra_type::vector_type, algebra_type::vector_type> >
							("AssembledOperator")
				.add_constructor()
				.add_method("set_discretization", &T::set_discretization)
				.add_method("set_dof_distribution", &T::set_dof_distribution)
				.add_method("init", &T::init);
		}

	//	StandardLineSearch
		{
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
			typedef NewtonSolver<dof_distribution_type, algebra_type> T;

			reg.add_class_<T, IOperatorInverse<algebra_type::vector_type, algebra_type::vector_type> >("NewtonSolver")
				.add_constructor()
				.add_method("set_linear_solver", &T::set_linear_solver)
				.add_method("set_convergence_check", &T::set_convergence_check)
				.add_method("set_line_search", &T::set_line_search)
				.add_method("init", &T::init);
		}


	//	Domain dependend part
		{
			typedef Domain<2, MultiGrid, MGSubsetHandler> domain_type;
			RegisterLibDiscretizationDomainDepended<domain_type>(reg);
		}



	//	Register user functions
		RegisterUserNumber(reg);
		RegisterUserVector(reg);
		RegisterUserMatrix(reg);

	//	Register Boundary functions
		RegisterBoundaryNumber(reg);

	//	todo: remove when possible
		RegisterElderUserFunctions(reg);
}

}//	end of namespace ug
}//	end of namespace interface
