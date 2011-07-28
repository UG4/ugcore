/*
 * lib_disc_bridge_domain_dependent.cpp
 *
 *  Created on: 22.09.2010
 *      Author: andreasvogel
 */

// extern headers
#include <iostream>
#include <sstream>

// include bridge
#include "../ug_bridge.h"

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
#include "lib_discretization/function_spaces/grid_function_util.h"
#include "lib_discretization/function_spaces/interpolate.h"
#include "lib_discretization/function_spaces/integrate.h"
#include "lib_discretization/function_spaces/error_indicator.h"
#include "lib_discretization/dof_manager/cuthill_mckee.h"
#include "lib_discretization/dof_manager/p1conform/p1conform.h"
#include "lib_discretization/dof_manager/conform/conform.h"

#include "lib_discretization/io/vtkoutput.h"

#include "lib_discretization/spatial_discretization/elem_disc/elem_disc_interface.h"
#include "lib_discretization/spatial_discretization/constraints/dirichlet_boundary/lagrange_dirichlet_boundary.h"

#include "lib_discretization/operator/linear_operator/projection_operator.h"
#include "lib_discretization/operator/linear_operator/prolongation_operator.h"
#include "lib_discretization/operator/linear_operator/multi_grid_solver/mg_solver.h"

#include "lib_discretization/spatial_discretization/elem_disc/level_set/level_set.h"

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
void RegisterLibDiscretizationDomainObjects(Registry& reg, const char* parentGroup)
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

//	GridFunction
	{
		std::stringstream ss; ss << "GridFunction" << dim << "d";
		reg.add_class_<function_type, vector_type>(ss.str().c_str(), grp.c_str())
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
	}

//  ApproximationSpace
	{
		typedef ApproximationSpace<domain_type, dof_distribution_type, algebra_type> T;
		std::stringstream ss; ss << "ApproximationSpace" << dim << "d";
		reg.add_class_<T,  IApproximationSpace<domain_type> >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("init|hide=true", &T::init)
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
	}

//	Order Cuthill-McKee
	{
		typedef ApproximationSpace<domain_type, dof_distribution_type, algebra_type> T;
		reg.add_function("OrderCuthillMcKee", static_cast<bool (*)(T&, bool)>(&OrderCuthillMcKee));
	}

//	DirichletBNDValues
	{
		typedef boost::function<bool (number& value, const MathVector<dim>& x, number time)> BNDNumberFunctor;
		typedef boost::function<void (number& value, const MathVector<dim>& x, number time)> NumberFunctor;
		typedef LagrangeDirichletBoundary<domain_type, dof_distribution_type, algebra_type> T;
		std::stringstream ss; ss << "DirichletBND" << dim << "d";
		reg.add_class_<T, IConstraint<dof_distribution_type, algebra_type> >(ss.str().c_str(), grp.c_str())
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
	}

//	ProlongationOperator
	{
		typedef ApproximationSpace<domain_type, dof_distribution_type, algebra_type> approximation_space_type;
		typedef P1ProlongationOperator<approximation_space_type, algebra_type> T;

		std::stringstream ss; ss << "P1ProlongationOperator" << dim << "d";
		reg.add_class_<T, IProlongationOperator<vector_type, vector_type> >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set_approximation_space", &T::set_approximation_space)
			.add_method("set_restriction_damping", &T::set_restriction_damping)
			.add_method("set_dirichlet_post_process", &T::set_dirichlet_post_process);

	}

//	ProjectionOperator
	{
		typedef ApproximationSpace<domain_type, dof_distribution_type, algebra_type> approximation_space_type;
		typedef P1ProjectionOperator<approximation_space_type, algebra_type> T;

		std::stringstream ss; ss << "P1ProjectionOperator" << dim << "d";
		reg.add_class_<T, IProjectionOperator<vector_type, vector_type> >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set_approximation_space", &T::set_approximation_space);
	}

//	AssembledMultiGridCycle
	{
		typedef ApproximationSpace<domain_type, dof_distribution_type, algebra_type> approximation_space_type;
		typedef AssembledMultiGridCycle<approximation_space_type, algebra_type> T;

		std::stringstream ss; ss << "GeometricMultiGridPreconditioner" << dim << "d";
		reg.add_class_<T, ILinearIterator<vector_type, vector_type> >(ss.str().c_str(), grp.c_str())
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
	}

//	VTK Output
	{
		typedef VTKOutput<function_type> T;
		std::stringstream ss; ss << "VTKOutput" << dim << "d";
		reg.add_class_<T>(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("write_time_pvd", &T::write_time_pvd)
			.add_method("clear_selection", &T::clear_selection)
			.add_method("select_all", &T::select_all)
			.add_method("select_nodal_scalar", &T::select_nodal_scalar)
			.add_method("select_nodal_vector", &T::select_nodal_vector)
			.add_method("print", static_cast<bool (T::*)(const char*, function_type&, int, number)>(&T::print))
			.add_method("print", static_cast<bool (T::*)(const char*, function_type&)>(&T::print));
	}


	//	GridFunctionDebugWriter
		{
			typedef GridFunctionDebugWriter<function_type> T;
			typedef IDebugWriter<typename function_type::algebra_type> TBase;
			std::stringstream ss; ss << "GridFunctionDebugWriter" << dim << "d";
			reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str())
				.add_constructor()
				.add_method("set_reference_grid_function", &T::set_reference_grid_function, "", "gridFunction")
				.add_method("set_vtk_output", &T::set_vtk_output, "", "vtkOutput")
				.add_method("set_conn_viewer_output", &T::set_conn_viewer_output, "", "cvOutput");

		}

	//	GridFunctionPositionProvider
	{
		typedef GridFunctionPositionProvider<function_type> T;
		typedef IPositionProvider<dim> TBase;
		std::stringstream ss; ss << "GridFunctionPositionProvider" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set_reference_grid_function", &T::set_reference_grid_function, "", "gridFunction");
	}


	//	GridFunctionVectorWriter
	{
		typedef GridFunctionVectorWriter<function_type, vector_type> T;
		typedef IVectorWriter<vector_type> TBase;
		std::stringstream ss; ss << "GridFunctionVectorWriter" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set_reference_grid_function", &T::set_reference_grid_function, "", "gridFunction")
			.add_method("set_user_data", &T::set_user_data, "", "userData");
	}

	// GridFunctionVectorWriterDirichlet0
	{
		typedef GridFunctionVectorWriterDirichlet0<function_type> T;
		std::stringstream ss; ss << "GridFunctionVectorWriterDirichlet0" << dim << "d";
		typedef IVectorWriter<vector_type> TBase;
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("init", &T::init, "", "postProcess#approxSpace#level")
			.add_method("set_level", &T::set_level, "", "level");
	}

	// FV1LevelSetDisc
	{
		typedef FV1LevelSetDisc<function_type> T;
		std::stringstream ss; ss << "FV1LevelSetDisc" << dim << "d";
		reg.add_class_<T>(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set_dt", &T::set_dt)
			.add_method("set_vel_scale", &T::set_vel_scale)
			.add_method("set_reinit", &T::set_reinit)
			.add_method("set_div_bool", &T::set_div_bool)
			.add_method("set_solution_nr", &T::set_solution_nr)
			.add_method("set_velocity_nr", &T::set_velocity_nr)
			.add_method("set_source", &T::set_source)
			.add_method("set_bool_source", &T::set_bool_source)
			.add_method("set_static_values_type", &T::set_static_values_type)
			.add_method("set_info", &T::set_info)
			.add_method("set_nr_of_steps", &T::set_nr_of_steps)
			.add_method("advect_lsf", &T::advect_lsf)
			.add_method("add_post_process", &T::add_post_process)
			.add_method("compute_error", &T::compute_error);
	}
}

template <typename TDomain, typename TAlgebra, typename TDoFDistribution>
void RegisterLibDiscretizationDomainFunctions(Registry& reg, const char* parentGroup)
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

//	WriteGridToVTK
	{
		reg.add_function("WriteGridFunctionToVTK",
						 &WriteGridFunctionToVTK<function_type>, grp.c_str(),
							"Success", "GridFunction#Filename|save-dialog",
							"Saves GridFunction to *.vtk file", "No help");
	}

//	SaveMatrixForConnectionViewer
	{
		reg.add_function("SaveMatrixForConnectionViewer",
						 &SaveMatrixForConnectionViewer<function_type>, grp.c_str());
	}

//	SaveVectorForConnectionViewer
	{
		reg.add_function("SaveVectorForConnectionViewer",
						 &SaveVectorForConnectionViewer<function_type>, grp.c_str());
	}

//	MarkForRefinement_GradientIndicator
	{
		reg.add_function("MarkForRefinement_GradientIndicator",
						 &MarkForRefinement_GradientIndicator<function_type>, grp.c_str());
	}

//	InterpolateFunction
	{
		typedef bool (*fct_type)(
				const boost::function<void (number& res,const MathVector<function_type::domain_type::dim>& x, number time)>&,
				function_type&, const char*, number);
		reg.add_function("InterpolateFunction",
						 static_cast<fct_type>(&InterpolateFunction<function_type>),
						 grp.c_str());

		typedef bool (*fct_type_subset)(
				const boost::function<void (number& res,const MathVector<function_type::domain_type::dim>& x, number time)>&,
				function_type&, const char*, number, const char*);
		reg.add_function("InterpolateFunction",
						 static_cast<fct_type_subset>(&InterpolateFunction<function_type>),
						 grp.c_str());
	}

//	L2Error
	{
		typedef number (*fct_type)(
				const boost::function<void (number& res,const MathVector<function_type::domain_type::dim>& x, number time)>&,
				function_type&, const char*, number);
		reg.add_function("L2Error",
						 static_cast<fct_type>(&L2Error<function_type>),
						 grp.c_str());

		typedef number (*fct_type_subset)(
				const boost::function<void (number& res,const MathVector<function_type::domain_type::dim>& x, number time)>&,
				function_type&, const char*, number, const char*);
		reg.add_function("L2Error",
						 static_cast<fct_type_subset>(&L2Error<function_type>),
						 grp.c_str());
	}

//	AssembleDirichletBoundary
	{
	//todo: This should work for all IDirichletPostProcess
		typedef MatrixOperator<vector_type, vector_type, matrix_type> mat_op_type;
		typedef LagrangeDirichletBoundary<domain_type, dof_distribution_type, algebra_type> dirichlet_type;
		typedef ApproximationSpace<domain_type, dof_distribution_type, algebra_type> approximation_space_type;

		typedef void (*fct_type)(	mat_op_type&,
									dirichlet_type&,
									approximation_space_type&);

		reg.add_function("AssembleDirichletRows",
						static_cast<fct_type>(&AssembleDirichletRows<mat_op_type, dirichlet_type, approximation_space_type>),
						grp.c_str());

		typedef void (*fct_type2)(	mat_op_type&,
									dirichlet_type&,
									approximation_space_type&,
									number);

		reg.add_function("AssembleDirichletRows",
				static_cast<fct_type2>(&AssembleDirichletRows<mat_op_type, dirichlet_type, approximation_space_type>),
				grp.c_str());
	}
}

template <typename TAlgebra, typename TDoFDistribution>
bool RegisterLibDiscForDomain(Registry& reg, const char* parentGroup)
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
	//	Domain dependent part 1D
		{
			typedef Domain<1, MultiGrid, MGSubsetHandler> domain_type;
			RegisterLibDiscretizationDomainObjects<domain_type, algebra_type, dof_distribution_type>(reg, grp.c_str());
			RegisterLibDiscretizationDomainFunctions<domain_type,  algebra_type, dof_distribution_type>(reg, grp.c_str());
		}
#endif

#ifdef UG_DIM_2
	//	Domain dependent part 2D
		{
			typedef Domain<2, MultiGrid, MGSubsetHandler> domain_type;
			RegisterLibDiscretizationDomainObjects<domain_type, algebra_type, dof_distribution_type>(reg, grp.c_str());
			RegisterLibDiscretizationDomainFunctions<domain_type,  algebra_type, dof_distribution_type>(reg, grp.c_str());
		}
#endif

#ifdef UG_DIM_3
	//	Domain dependent part 3D
		{
			typedef Domain<3, MultiGrid, MGSubsetHandler> domain_type;
			RegisterLibDiscretizationDomainObjects<domain_type, algebra_type, dof_distribution_type>(reg, grp.c_str());
			RegisterLibDiscretizationDomainFunctions<domain_type,  algebra_type, dof_distribution_type>(reg, grp.c_str());
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

bool RegisterDynamicLibDiscDomain(Registry& reg, int algebra_type, const char* parentGroup)
{
	bool bReturn = true;

	typedef P1DoFDistribution<false> dd_type_single;
	typedef P1DoFDistribution<true> dd_type_grouped;

//	typedef DoFDistribution dd_type_single;
//	typedef DoFDistribution dd_type_grouped;

	switch(algebra_type)
	{
	case eCPUAlgebra:		 		bReturn &= RegisterLibDiscForDomain<CPUAlgebra, dd_type_single>(reg, parentGroup); break;
//	case eCPUBlockAlgebra2x2: 		bReturn &= RegisterLibDiscForDomain<CPUBlockAlgebra<2>, dd_type_groupded>(reg, parentGroup); break;
	case eCPUBlockAlgebra3x3: 		bReturn &= RegisterLibDiscForDomain<CPUBlockAlgebra<3>, dd_type_grouped>(reg, parentGroup); break;
//	case eCPUBlockAlgebra4x4: 		bReturn &= RegisterLibDiscForDomain<CPUBlockAlgebra<4>, dd_type_grouped>(reg, parentGroup); break;
//	case eCPUVariableBlockAlgebra: 	bReturn &= RegisterLibDiscForDomain<CPUVariableBlockAlgebra, dd_type_grouped>(reg, parentGroup); break;
	default: UG_ASSERT(0, "Unsupported Algebra Type");
				UG_LOG("Unsupported Algebra Type requested.\n");
				return false;
	}

	return bReturn;
}

}//	end of namespace ug
}//	end of namespace interface
