/*
 * lib_disc_bridge_domain_dependent.cpp
 *
 *  Created on: 22.09.2010
 *      Author: andreasvogel
 */

// extern headers
#include <iostream>
#include <sstream>
#include <string>

// include bridge
#include "../ug_bridge.h"

// lib_algebra includes
#include "lib_algebra/algebra_selector.h"x
#include "lib_algebra/algebra_types.h"
#include "lib_algebra/operator/operator_util.h"
#include "lib_algebra/operator/operator_interface.h"
#include "lib_algebra/operator/operator_inverse_interface.h"

// lib_disc includes
#include "lib_discretization/domain.h"
#include "lib_discretization/function_spaces/grid_function.h"
#include "lib_discretization/function_spaces/approximation_space.h"
#include "lib_discretization/function_spaces/grid_function_util.h"
#include "lib_discretization/function_spaces/interpolate.h"
#include "lib_discretization/function_spaces/integrate.h"
#include "lib_discretization/function_spaces/error_indicator.h"
#include "lib_discretization/dof_manager/cuthill_mckee.h"
#include "lib_discretization/dof_manager/p1conform/p1conform.h"
#include "lib_discretization/dof_manager/cuthill_mckee.h"
#include "lib_discretization/dof_manager/lexorder.h"
#include "lib_discretization/dof_manager/conform/conform.h"
#include "lib_discretization/dof_manager/p1conform/p1conform.h"

#include "lib_discretization/io/vtkoutput.h"

#include "lib_discretization/spatial_discretization/domain_discretization.h"
#include "lib_discretization/spatial_discretization/elem_disc/elem_disc_interface.h"
#include "lib_discretization/spatial_discretization/constraints/dirichlet_boundary/lagrange_dirichlet_boundary.h"

#include "lib_discretization/operator/linear_operator/projection_operator.h"
#include "lib_discretization/operator/linear_operator/prolongation_operator.h"
#include "lib_discretization/operator/linear_operator/multi_grid_solver/mg_solver.h"

#include "lib_discretization/spatial_discretization/elem_disc/level_set/level_set.h"

using namespace std;

namespace ug
{

namespace bridge
{

/**	Calls e.g. LagrangeDirichletBoundary::assemble_dirichlet_rows.
 *
 * This method probably shouldn't be implemented here, but in some util file.
 */
template <class TMatOp, class TDirichletBnd, class TApproxSpace>
void AssembleDirichletRows(TMatOp& matOp, TDirichletBnd& dirichletBnd,
							TApproxSpace& approxSpace)
{
	dirichletBnd.assemble_dirichlet_rows(matOp.get_matrix(),
						approxSpace.get_surface_dof_distribution());
}

/**	Calls e.g. LagrangeDirichletBoundary::assemble_dirichlet_rows.
 * Also takes a time argument.
 *
 * This method probably shouldn't be implemented here, but in some util file.
 */
template <class TMatOp, class TDirichletBnd, class TApproxSpace>
void AssembleDirichletRows(TMatOp& matOp, TDirichletBnd& dirichletBnd,
							TApproxSpace& approxSpace, number time)
{
	dirichletBnd.assemble_dirichlet_rows(matOp.get_matrix(),
					approxSpace.get_surface_dof_distribution(), time);
}



/// small wrapper to write a grid function to vtk
template <typename TGridFunction>
bool WriteGridFunctionToVTK(TGridFunction& u, const char* filename)
{
	VTKOutput<TGridFunction> out;
	return out.print(filename, u);
}

template <typename TDomain, typename TAlgebra, typename TDoFDistribution>
void RegisterLibDiscDomain__Algebra_DoFDistribution_Domain(Registry& reg, string parentGroup)
{
//	typedef
	static const int dim = TDomain::dim;
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef ApproximationSpace<TDomain, TDoFDistribution, TAlgebra> approximation_space_type;

#ifdef UG_PARALLEL
		typedef ParallelGridFunction<GridFunction<TDomain, TDoFDistribution, TAlgebra> > function_type;
#else
		typedef GridFunction<TDomain, TDoFDistribution, TAlgebra> function_type;
#endif

//	group string
	stringstream grpSS; grpSS << parentGroup << "/" << dim << "d";
	string grp = grpSS.str();

//	suffix and tag
	string dimAlgDDSuffix = GetDomainSuffix<TDomain>();
	dimAlgDDSuffix.append(GetAlgebraSuffix<TAlgebra>());
	dimAlgDDSuffix.append(GetDoFDistributionSuffix<TDoFDistribution>());

	string dimAlgDDTag = GetDomainTag<TDomain>();
	dimAlgDDTag.append(GetAlgebraTag<TAlgebra>());
	dimAlgDDTag.append(GetDoFDistributionTag<TDoFDistribution>());

//	GridFunction
	{
		string name = string("GridFunction").append(dimAlgDDSuffix);
		reg.add_class_<function_type, vector_type>(name, grp)
			.add_constructor()
			.add_method("assign", static_cast<bool (function_type::*)(const vector_type&)>(&function_type::assign),
						"Success", "Vector")
			.add_method("assign_dof_distribution|hide=true", &function_type::assign_dof_distribution)
			.add_method("get_dim|hide=true", &function_type::get_dim)
			.add_method("assign_approximation_space|hide=true", &function_type::assign_approximation_space)
			.add_method("clone", &function_type::clone);
#ifdef UG_PARALLEL
		reg.get_class_<function_type>()
			.add_method("change_storage_type_by_string|hide=true", &function_type::change_storage_type_by_string)
			.add_method("set_storage_type_by_string|hide=true", &function_type::set_storage_type_by_string);
#endif
		reg.add_class_to_group(name, "GridFunction", dimAlgDDTag);
	}

//  ApproximationSpace
	{
		typedef ApproximationSpace<TDomain, TDoFDistribution, TAlgebra> T;
		typedef IApproximationSpace<TDomain> TBase;
		string name = string("ApproximationSpace").append(dimAlgDDSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("init|hide=true", &T::init)
			.add_method("set_grouping", &T::set_grouping)
			.add_method("print_statistic|hide=true", static_cast<void (T::*)(int) const>(&T::print_statistic))
			.add_method("print_statistic|hide=true", static_cast<void (T::*)() const>(&T::print_statistic))
			.add_method("print_layout_statistic|hide=true", static_cast<void (T::*)(int) const>(&T::print_layout_statistic))
			.add_method("print_layout_statistic|hide=true", static_cast<void (T::*)() const>(&T::print_layout_statistic))
			.add_method("print_local_dof_statistic|hide=true", static_cast<void (T::*)(int) const>(&T::print_local_dof_statistic))
			.add_method("print_local_dof_statistic|hide=true", static_cast<void (T::*)() const>(&T::print_local_dof_statistic))
			.add_method("defragment|hide=true", &T::defragment)
			.add_method("get_surface_view|hide=true", &T::get_surface_view)
			.add_method("get_surface_dof_distribution|hide=true",  static_cast<const typename T::dof_distribution_type& (T::*)() const>(&T::get_surface_dof_distribution))
			.add_method("create_surface_function|hide=true", &T::create_surface_function);
		reg.add_class_to_group(name, "ApproximationSpace", dimAlgDDTag);
	}

//	DomainDiscretization
	{
		typedef IDomainDiscretization<TDoFDistribution, TAlgebra> TBase;
		typedef DomainDiscretization<TDomain, TDoFDistribution, TAlgebra> T;
		typedef typename T::dof_distribution_type dof_distribution_type;
		string name = string("DomainDiscretization").append(dimAlgDDSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("set_approximation_space", &T::set_approximation_space)
			.add_method("add|interactive=false", static_cast<bool (T::*)(IConstraint<TDoFDistribution, TAlgebra>&)>(&T::add),
						"", "Post Process")
			.add_method("add|interactive=false", static_cast<bool (T::*)(IDomainElemDisc<TDomain>&)>(&T::add),
						"", "Discretization")
			.add_method("assemble_linear", static_cast<bool (T::*)(matrix_type&, vector_type&, const vector_type&)>(&T::assemble_linear))
			.add_method("assemble_solution", static_cast<bool (T::*)(vector_type&)>(&T::assemble_solution))
			.add_method("assemble_mass_matrix", static_cast<bool (T::*)(matrix_type&, const vector_type&, const dof_distribution_type&)>(&T::assemble_mass_matrix))
			.add_method("assemble_mass_matrix", static_cast<bool (T::*)(matrix_type&, const vector_type&)>(&T::assemble_mass_matrix))
			.add_method("assemble_stiffness_matrix", static_cast<bool (T::*)(matrix_type&, const vector_type&, const dof_distribution_type&)>(&T::assemble_stiffness_matrix))
			.add_method("assemble_stiffness_matrix", static_cast<bool (T::*)(matrix_type&, const vector_type&)>(&T::assemble_stiffness_matrix))
			.add_method("assemble_rhs", static_cast<bool (T::*)(vector_type&, const vector_type&, const dof_distribution_type&)>(&T::assemble_rhs))
			.add_method("assemble_rhs", static_cast<bool (T::*)(vector_type&, const vector_type&)>(&T::assemble_rhs));
		reg.add_class_to_group(name, "DomainDiscretization", dimAlgDDTag);
	}

//	Order Cuthill-McKee
	{
		reg.add_function("OrderCuthillMcKee", static_cast<bool (*)(approximation_space_type&, bool)>(&OrderCuthillMcKee));
	}

//	Order lexicographically
	{
		typedef ApproximationSpace<TDomain, TDoFDistribution, TAlgebra> T;
		reg.add_function("OrderLex", (bool (*)(T&, const char*))&OrderLex);
	}

//	DirichletBNDValues
	{
		typedef boost::function<bool (number& value, const MathVector<dim>& x, number time)> BNDNumberFunctor;
		typedef boost::function<void (number& value, const MathVector<dim>& x, number time)> NumberFunctor;
		typedef LagrangeDirichletBoundary<TDomain, TDoFDistribution, TAlgebra> T;
		typedef IConstraint<TDoFDistribution, TAlgebra> TBase;
		string name = string("DirichletBND").append(dimAlgDDSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("set_approximation_space|interactive=false", &T::set_approximation_space,
						"", "Approximation Space")
			.add_method("add", static_cast<void (T::*)(BNDNumberFunctor&, const char*, const char*)>(&T::add),
						"Success", "Value#Function#Subsets")
			.add_method("add", static_cast<void (T::*)(NumberFunctor&, const char*, const char*)>(&T::add),
						"Success", "Value#Function#Subsets")
			.add_method("add",static_cast<void (T::*)(number, const char*, const char*)>(&T::add),
						"Success", "Constant Value#Function#Subsets")
			.add_method("clear", &T::clear);
		reg.add_class_to_group(name, "DirichletBND", dimAlgDDTag);
	}

//	ProlongationOperator
	{
		typedef P1ProlongationOperator<approximation_space_type, TAlgebra> T;
		typedef IProlongationOperator<vector_type, vector_type> TBase;
		string name = string("P1ProlongationOperator").append(dimAlgDDSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("set_approximation_space", &T::set_approximation_space)
			.add_method("set_restriction_damping", &T::set_restriction_damping)
			.add_method("set_dirichlet_post_process", &T::set_dirichlet_post_process);
		reg.add_class_to_group(name, "P1ProlongationOperator", dimAlgDDTag);
	}

//	ProjectionOperator
	{
		typedef P1ProjectionOperator<approximation_space_type, TAlgebra> T;
		typedef IProjectionOperator<vector_type, vector_type> TBase;
		string name = string("P1ProjectionOperator").append(dimAlgDDSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("set_approximation_space", &T::set_approximation_space);
		reg.add_class_to_group(name, "P1ProjectionOperator", dimAlgDDTag);
	}

//	AssembledMultiGridCycle
	{
		typedef AssembledMultiGridCycle<approximation_space_type, TAlgebra> T;
		typedef ILinearIterator<vector_type, vector_type> TBase;
		string name = string("GeometricMultiGridPreconditioner").append(dimAlgDDSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("set_discretization|interactive=false", &T::set_discretization,
						"", "Discretization")
			.add_method("set_approximation_space|interactive=false", &T::set_approximation_space,
						"", "Approximation Space")
			.add_method("set_base_level|interactive=false", &T::set_base_level,
						"", "Base Level")
			.add_method("set_parallel_base_solver", &T::set_parallel_base_solver,
						"", "Specifies if base solver works in parallel")
			.add_method("set_base_solver|interactive=false", &T::set_base_solver,
						"","Base Solver")
			.add_method("set_smoother|interactive=false", &T::set_smoother,
						"", "Smoother")
			.add_method("set_cycle_type|interactive=false", &T::set_cycle_type,
						"", "Cycle Type")
			.add_method("set_num_presmooth|interactive=false", &T::set_num_presmooth,
						"", "Number PreSmooth Steps")
			.add_method("set_num_postsmooth|interactive=false", &T::set_num_postsmooth,
						"", "Number PostSmooth Steps")
			.add_method("set_prolongation|interactive=false", &T::set_prolongation_operator,
						"", "Prolongation")
			.add_method("set_projection|interactive=false", &T::set_projection_operator,
						"", "Projection")
			.add_method("set_debug", &T::set_debug);
		reg.add_class_to_group(name, "GeometricMultiGridPreconditioner", dimAlgDDTag);
	}

//	VTK Output
	{
		typedef VTKOutput<function_type> T;
		string name = string("VTKOutput").append(dimAlgDDSuffix);
		reg.add_class_<T>(name, grp)
			.add_constructor()
			.add_method("write_time_pvd", &T::write_time_pvd)
			.add_method("clear_selection", &T::clear_selection)
			.add_method("select_all", &T::select_all)
			.add_method("select_nodal_scalar", &T::select_nodal_scalar)
			.add_method("select_nodal_vector", &T::select_nodal_vector)
			.add_method("print", static_cast<bool (T::*)(const char*, function_type&, int, number)>(&T::print))
			.add_method("print", static_cast<bool (T::*)(const char*, function_type&)>(&T::print));
		reg.add_class_to_group(name, "VTKOutput", dimAlgDDTag);
	}


//	GridFunctionDebugWriter
	{
		typedef GridFunctionDebugWriter<function_type> T;
		typedef IDebugWriter<TAlgebra> TBase;
		string name = string("GridFunctionDebugWriter").append(dimAlgDDSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("set_reference_grid_function", &T::set_reference_grid_function, "", "gridFunction")
			.add_method("set_vtk_output", &T::set_vtk_output, "", "vtkOutput")
			.add_method("set_conn_viewer_output", &T::set_conn_viewer_output, "", "cvOutput");
		reg.add_class_to_group(name, "GridFunctionDebugWriter", dimAlgDDTag);
	}

//	GridFunctionPositionProvider
	{
		typedef GridFunctionPositionProvider<function_type> T;
		typedef IPositionProvider<dim> TBase;
		string name = string("GridFunctionPositionProvider").append(dimAlgDDSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("set_reference_grid_function", &T::set_reference_grid_function, "", "gridFunction");
		reg.add_class_to_group(name, "GridFunctionPositionProvider", dimAlgDDTag);
	}


	//	GridFunctionVectorWriter
	{
		typedef GridFunctionVectorWriter<function_type, vector_type> T;
		typedef IVectorWriter<vector_type> TBase;
		string name = string("GridFunctionVectorWriter").append(dimAlgDDSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("set_reference_grid_function", &T::set_reference_grid_function, "", "gridFunction")
			.add_method("set_user_data", &T::set_user_data, "", "userData");
		reg.add_class_to_group(name, "GridFunctionVectorWriter", dimAlgDDTag);
	}

	// GridFunctionVectorWriterDirichlet0
	{
		typedef GridFunctionVectorWriterDirichlet0<function_type> T;
		typedef IVectorWriter<vector_type> TBase;
		string name = string("GridFunctionVectorWriterDirichlet0").append(dimAlgDDSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("init", &T::init, "", "postProcess#approxSpace#level")
			.add_method("set_level", &T::set_level, "", "level");
		reg.add_class_to_group(name, "GridFunctionVectorWriterDirichlet0", dimAlgDDTag);
	}

// 	FV1LevelSetDisc
	{
		typedef FV1LevelSetDisc<function_type> T;
		typedef typename function_type::domain_type domain_type;
		typedef boost::function<void (number& value,
						                              const MathVector<domain_type::dim>& x,
						                              number time)> NumberFunctor;
		string name = string("FV1LevelSetDisc").append(dimAlgDDSuffix);
		reg.add_class_<T>(name, grp)
			.add_constructor()
			.add_method("set_dt", &T::set_dt)
			.add_method("set_vel_scale", &T::set_vel_scale)
			.add_method("set_reinit", &T::set_reinit)
			.add_method("set_divfree_bool", &T::set_divfree_bool)
			.add_method("set_info", &T::set_info)
			.add_method("set_timestep_nr",&T::set_timestep_nr)
			.add_method("set_nr_of_steps", &T::set_nr_of_steps)
			.add_method("advect_lsf", &T::advect_lsf)
			.add_method("add_post_process", &T::add_post_process)
			.add_method("set_delta", &T::set_delta)
			.add_method("set_limiter",&T::set_limiter)
			.add_method("set_neumann_boundary",&T::set_neumann_boundary)
			.add_method("init_function", &T::init_function)
			.add_method("set_vel_x", static_cast<void (T::*)(const NumberFunctor&)>(&T::set_vel_x))
			.add_method("set_vel_y", static_cast<void (T::*)(const NumberFunctor&)>(&T::set_vel_y))
			.add_method("set_vel_z", static_cast<void (T::*)(const NumberFunctor&)>(&T::set_vel_z))
			.add_method("set_source", static_cast<void (T::*)(const NumberFunctor&)>(&T::set_source))
			.add_method("set_vel_x", static_cast<void (T::*)(function_type&)>(&T::set_vel_x))
			.add_method("set_vel_y", static_cast<void (T::*)(function_type&)>(&T::set_vel_y))
			.add_method("set_vel_z", static_cast<void (T::*)(function_type&)>(&T::set_vel_z))
			.add_method("set_vel_x", static_cast<void (T::*)(number)>(&T::set_vel_x))
			.add_method("set_vel_y", static_cast<void (T::*)(number)>(&T::set_vel_y))
			.add_method("set_vel_z", static_cast<void (T::*)(number)>(&T::set_vel_z))
			.add_method("set_vel_x", static_cast<void (T::*)()>(&T::set_vel_x))
			.add_method("set_vel_y", static_cast<void (T::*)()>(&T::set_vel_y))
			.add_method("set_vel_z", static_cast<void (T::*)()>(&T::set_vel_z))
			.add_method("set_source", static_cast<void (T::*)(function_type&)>(&T::set_source))
			.add_method("set_source", static_cast<void (T::*)(number)>(&T::set_source))
			.add_method("set_source", static_cast<void (T::*)()>(&T::set_source))
			.add_method("set_dirichlet_data", static_cast<void (T::*)(number)>(&T::set_dirichlet_data))
			.add_method("set_dirichlet_data", static_cast<void (T::*)(const NumberFunctor&)>(&T::set_dirichlet_data))
			.add_method("set_dirichlet_data", static_cast<void (T::*)()>(&T::set_dirichlet_data))
			.add_method("fill_v_vec",&T::fill_v_vec)
			.add_method("get_time",&T::get_time)
			.add_method("runtimetest", &T::runtimetest)
			.add_method("compute_error", &T::compute_error);
		reg.add_class_to_group(name, "FV1LevelSetDisc", dimAlgDDTag);
	}

//	WriteGridToVTK
	{
		reg.add_function("WriteGridFunctionToVTK",
						 &WriteGridFunctionToVTK<function_type>, grp,
							"Success", "GridFunction#Filename|save-dialog",
							"Saves GridFunction to *.vtk file", "No help");
	}

//	SaveMatrixForConnectionViewer
	{
		reg.add_function("SaveMatrixForConnectionViewer",
						 &SaveMatrixForConnectionViewer<function_type>, grp);
	}

//	SaveVectorForConnectionViewer
	{
		reg.add_function("SaveVectorForConnectionViewer",
						 &SaveVectorForConnectionViewer<function_type>, grp);
	}

//	MarkForRefinement_GradientIndicator
	{
		reg.add_function("MarkForRefinement_GradientIndicator",
						 &MarkForRefinement_GradientIndicator<function_type>, grp);
	}

//	InterpolateFunction
	{
		typedef bool (*fct_type)(
				const boost::function<void (number& res,const MathVector<dim>& x, number time)>&,
				function_type&, const char*, number);
		reg.add_function("InterpolateFunction",
						 static_cast<fct_type>(&InterpolateFunction<function_type>),
						 grp);

		typedef bool (*fct_type_subset)(
				const boost::function<void (number& res,const MathVector<dim>& x, number time)>&,
				function_type&, const char*, number, const char*);
		reg.add_function("InterpolateFunction",
						 static_cast<fct_type_subset>(&InterpolateFunction<function_type>),
						 grp);
	}

//	L2Error
	{
		typedef number (*fct_type)(
				const boost::function<void (number& res,const MathVector<dim>& x, number time)>&,
				function_type&, const char*, number);
		reg.add_function("L2Error",
						 static_cast<fct_type>(&L2Error<function_type>),
						 grp);

		typedef number (*fct_type_subset)(
				const boost::function<void (number& res,const MathVector<dim>& x, number time)>&,
				function_type&, const char*, number, const char*);
		reg.add_function("L2Error",
						 static_cast<fct_type_subset>(&L2Error<function_type>),
						 grp);
	}

//	AssembleDirichletBoundary
	{
		typedef MatrixOperator<vector_type, vector_type, matrix_type> mat_op_type;
		typedef LagrangeDirichletBoundary<TDomain, TDoFDistribution, TAlgebra> dirichlet_type;
		typedef ApproximationSpace<TDomain, TDoFDistribution, TAlgebra> approximation_space_type;

		typedef void (*fct_type)(	mat_op_type&,
									dirichlet_type&,
									approximation_space_type&);

		reg.add_function("AssembleDirichletRows",
						static_cast<fct_type>(&AssembleDirichletRows<mat_op_type, dirichlet_type, approximation_space_type>),
						grp);

		typedef void (*fct_type2)(	mat_op_type&,
									dirichlet_type&,
									approximation_space_type&,
									number);

		reg.add_function("AssembleDirichletRows",
				static_cast<fct_type2>(&AssembleDirichletRows<mat_op_type, dirichlet_type, approximation_space_type>),
				grp);
	}
}

template <typename TAlgebra, typename TDoFDistribution>
static bool RegisterLibDiscDomain__Algebra_DoFDistribution(Registry& reg, string parentGroup)
{
//	get group string
	string grp = parentGroup; grp.append("/Discretization");

	try
	{

#ifdef UG_DIM_1
//	Domain dependent part 1D
	{
		typedef Domain<1, MultiGrid, MGSubsetHandler> domain_type;
		RegisterLibDiscDomain__Algebra_DoFDistribution_Domain<domain_type, TAlgebra, TDoFDistribution>(reg, grp);
	}
#endif

#ifdef UG_DIM_2
//	Domain dependent part 2D
	{
		typedef Domain<2, MultiGrid, MGSubsetHandler> domain_type;
		RegisterLibDiscDomain__Algebra_DoFDistribution_Domain<domain_type, TAlgebra, TDoFDistribution>(reg, grp);
	}
#endif

#ifdef UG_DIM_3
//	Domain dependent part 3D
	{
		typedef Domain<3, MultiGrid, MGSubsetHandler> domain_type;
		RegisterLibDiscDomain__Algebra_DoFDistribution_Domain<domain_type, TAlgebra, TDoFDistribution>(reg, grp);
	}
#endif

	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterLibDiscDomain__Algebra_DoFDistribution: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}

	return true;
}

template <typename TAlgebra>
static bool RegisterLibDiscDomain__Algebra(Registry& reg, string parentGroup)
{
	bool bReturn = true;
#ifdef DOF_P1
	bReturn &= RegisterLibDiscDomain__Algebra_DoFDistribution<TAlgebra, P1DoFDistribution>(reg, parentGroup);
#endif
#ifdef DOF_GEN
	bReturn &= RegisterLibDiscDomain__Algebra_DoFDistribution<TAlgebra, DoFDistribution >(reg, parentGroup);
#endif

	return bReturn;
}

bool RegisterLibDisc_Domain(Registry& reg, string parentGroup)
{
	bool bReturn = true;
	bReturn &= RegisterLibDiscDomain__Algebra<CPUAlgebra>(reg, parentGroup);
//	bReturn &= RegisterLibDiscDomain__Algebra<CPUBlockAlgebra<2> >(reg, parentGroup);
	bReturn &= RegisterLibDiscDomain__Algebra<CPUBlockAlgebra<3> >(reg, parentGroup);
//	bReturn &= RegisterLibDiscDomain__Algebra<CPUBlockAlgebra<4> >(reg, parentGroup);
//	bReturn &= RegisterLibDiscDomain__Algebra<CPUVariableBlockAlgebra >(reg, parentGroup);
	return bReturn;
}

}//	end of namespace ug
}//	end of namespace interface
