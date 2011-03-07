/*
 * lib_disc_bridge_domain_dependent.cpp
 *
 *  Created on: 22.09.2010
 *      Author: andreasvogel
 */

#include "../ug_bridge.h"
#include "common/math/ugmath.h"
#include "common/math/misc/math_util.h"
#include "lib_algebra/lib_algebra.h"
#include "lib_discretization/lib_discretization.h"
#ifdef UG_PARALLEL
#include "pcl/pcl.h"
#endif
#include <iostream>
#include <sstream>

namespace ug
{
extern enum_AlgebraType g_AlgebraType;

namespace bridge
{

#ifdef UG_PARALLEL
template <class grid_function_type>
void OneToManyTests(grid_function_type& func)
{
	IndexLayout oneToManyMasterLayout;
	IndexLayout oneToManySlaveLayout;
	UG_LOG("executing OneToManyTests...\n");
	LogIndexLayout(func.get_master_layout());
	BuildOneToManyLayout(oneToManyMasterLayout,
						 oneToManySlaveLayout,
						 0,
						 func.get_master_layout(),
						 func.get_slave_layout(),
						 pcl::ProcessCommunicator(pcl::PCD_WORLD));
}

template <class grid_function_type>
void BuildDomainDecompositionLayoutsTest(grid_function_type& func,
										pcl::IDomainDecompositionInfo& ddinfo)
{
	IndexLayout subdomMasters;
	IndexLayout subdomSlaves;
	IndexLayout processMasters;
	IndexLayout processSlaves;
	IndexLayout deltaNbrMasters;
	IndexLayout deltaNbrSlaves;
	IndexLayout crossPointMasters;
	IndexLayout crossPointSlaves;

	UG_LOG("executing BuildDomainDecompositionLayoutsTest...\n");
	UG_LOG("standard master layout: ");
	LogIndexLayout(func.get_master_layout());
	UG_LOG("standard slave layout: ");
	LogIndexLayout(func.get_slave_layout());

	BuildDomainDecompositionLayouts(subdomMasters, subdomSlaves,
					processMasters, processSlaves, deltaNbrMasters,
					deltaNbrSlaves, crossPointMasters, crossPointSlaves,
					func.get_master_layout(), func.get_slave_layout(),
					(int)func.num_dofs() - 1, ddinfo);

	UG_LOG("done\n");

	UG_LOG("subdomMasters: ");
	LogIndexLayout(subdomMasters);
	UG_LOG("subdomSlaves: ");
	LogIndexLayout(subdomSlaves);

	UG_LOG("processMasters: ");
	LogIndexLayout(processMasters);
	UG_LOG("processSlaves: ");
	LogIndexLayout(processSlaves);

	UG_LOG("deltaNbrMasters: ");
	LogIndexLayout(deltaNbrMasters);
	UG_LOG("deltaNbrSlaves: ");
	LogIndexLayout(deltaNbrSlaves);

	UG_LOG("crossPointMasters: ");
	LogIndexLayout(crossPointMasters);
	UG_LOG("crossPointSlaves: ");
	LogIndexLayout(crossPointSlaves);

}

template <class grid_function_type>
void TestGridFunctionLayout(grid_function_type& func)
{
	pcl::ParallelCommunicator<IndexLayout> com;
	pcl::TestLayout(com, func.get_master_layout(),
					func.get_slave_layout(), true);
}

#endif

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
//	check that extension '.mat' is chosen
	const char * p = strstr(filename, ".mat");
	if(p == NULL)
	{
		UG_LOG("Currently only '.mat' format supported for domains.\n");
		return false;
	}

	std::string name(filename);

#ifdef UG_PARALLEL
//	search for ending
	size_t found = name.find_first_of(".");

//	remove endings
	name.resize(found);

//	add new ending, containing process number
	int rank = pcl::GetProcRank();
	char ext[20];
	sprintf(ext, "_p%04d.mat", rank);
	name.append(ext);
#endif

	static const int dim = TGridFunction::domain_type::dim;

	std::vector<MathVector<dim> > positions;
	ExtractPositions(u, positions);

	WriteMatrixToConnectionViewer(name.c_str(), A.get_matrix(), &positions[0], dim);
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

	std::string name(filename);

#ifdef UG_PARALLEL
	size_t found = name.find_first_of(".");
	name.resize(found);

	int rank = pcl::GetProcRank();
	char ext[20];
	sprintf(ext, "_p%04d.mat", rank);

	name.append(ext);
#endif


	std::vector<MathVector<dim> > positions;
	ExtractPositions(b, positions);

	WriteVectorToConnectionViewer(name.c_str(), b, &positions[0], dim);
	return true;
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

//	IApproximationSpace
	{
		typedef IApproximationSpace<domain_type> T;
		std::stringstream ss; ss << "IApproximationSpace" << dim << "d";
		reg.add_class_<T, FunctionPattern>(ss.str().c_str(), grp.c_str())
			.add_method("assign_domain|hide=true", &T::assign_domain)
			.add_method("get_domain|hide=true", (domain_type& (T::*)())&T::get_domain);
	}

//	GridFunction
	{
		std::stringstream ss; ss << "GridFunction" << dim << "d";
		reg.add_class_<function_type, vector_type>(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("assign|hide=true", (bool (function_type::*)(const vector_type&))&function_type::assign,
						"Success", "Vector")
			.add_method("assign_dof_distribution|hide=true", &function_type::assign_dof_distribution)
			.add_method("get_dim|hide=true", &function_type::get_dim)
			.add_method("assign_approximation_space|hide=true", &function_type::assign_approximation_space)
			.add_method("set_space|interactive=false", &function_type::assign_surface_approximation_space,
								"", "Approximation Space")
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
			.add_method("print_statistic|hide=true", &T::print_statistic)
			.add_method("print_layout_statistic|hide=true", &T::print_layout_statistic)
			.add_method("get_surface_dof_distribution|hide=true",  (const typename T::dof_distribution_type& (T::*)() const) &T::get_surface_dof_distribution)
			.add_method("create_surface_function|hide=true", &T::create_surface_function);
	}

//	DirichletBNDValues
	{
		typedef P1DirichletBoundary<domain_type, dof_distribution_type, algebra_type> T;
		std::stringstream ss; ss << "DirichletBND" << dim << "d";
		reg.add_class_<T, IPostProcess<dof_distribution_type, algebra_type> >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set_domain|hide=true", &T::set_domain)
			.add_method("set_pattern|hide=true", &T::set_pattern)
			.add_method("set_approximation_space|interactive=false", &T::set_approximation_space,
						"", "Approximation Space")
			.add_method("add_boundary_value", (bool (T::*)(IBoundaryData<number, dim>&, const char*, const char*))&T::add_boundary_value,
						"Success", "Value#Function#Subsets")
			.add_method("add_constant_boundary_value", &T::add_constant_boundary_value,
						"Success", "Constant Value#Function#Subsets");
	}

//	DomainElemDisc base class
	{
		typedef IDomainElemDisc<domain_type, algebra_type> T;
		typedef IElemDisc<algebra_type> TBase;

		std::stringstream ss; ss << "IDomainElemDisc" << dim << "d";
		reg.add_class_<T, TBase >(ss.str().c_str(), grp.c_str())
			.add_method("set_approximation_space", &T::set_approximation_space);
	}

//	Neumann Boundary
	{
		typedef FVNeumannBoundaryElemDisc<domain_type, algebra_type> T;
		typedef IDomainElemDisc<domain_type, algebra_type> TBase;
		std::stringstream ss; ss << "FV1NeumannBoundary" << dim << "d";
		reg.add_class_<T, TBase >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("add_boundary_value", (bool (T::*)(IBoundaryData<number, dim>&, const char*, const char*))&T::add_boundary_value);
	}

//	Inner Boundary
	{
		typedef FVInnerBoundaryElemDisc<domain_type, algebra_type> T;
		typedef IDomainElemDisc<domain_type, algebra_type> TBase;
		std::stringstream ss; ss << "FV1InnerBoundary" << dim << "d";
		reg.add_class_<T, TBase >(ss.str().c_str(), grp.c_str())
			.add_constructor();
	}

//	Convection Diffusion Finite Volume
	{
		typedef FVConvectionDiffusionElemDisc<domain_type, algebra_type> T;
		typedef IDomainElemDisc<domain_type, algebra_type> TBase;
		std::stringstream ss; ss << "FV1ConvectionDiffusion" << dim << "d";
		reg.add_class_<T, TBase >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set_diffusion_tensor", &T::set_diffusion)
			.add_method("set_velocity_field", &T::set_velocity)
			.add_method("set_reaction", &T::set_reaction)
			.add_method("set_rhs", &T::set_rhs)
			.add_method("set_mass_scale", &T::set_mass_scale)
			.add_method("set_upwind_amount", &T::set_upwind_amount);
	}

//	Convection Diffusion Finite Element
	{
		typedef FE1ConvectionDiffusionElemDisc<domain_type, algebra_type> T;
		typedef IDomainElemDisc<domain_type, algebra_type> TBase;
		std::stringstream ss; ss << "FE1ConvectionDiffusion" << dim << "d";
		reg.add_class_<T, TBase >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set_diffusion_tensor", &T::set_diffusion)
			.add_method("set_velocity_field", &T::set_velocity)
			.add_method("set_reaction", &T::set_reaction)
			.add_method("set_rhs", &T::set_rhs)
			.add_method("set_mass_scale", &T::set_mass_scale);
	}


//	Density Driven Flow
	{
		typedef DensityDrivenFlowElemDisc<FV1Geometry, domain_type, algebra_type> T2;
		typedef IDomainElemDisc<domain_type, algebra_type> TBase;
		std::stringstream ss; ss << "DensityDrivenFlow" << dim << "d";
		reg.add_class_<T2, TBase >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set_approximation_space", &T2::set_approximation_space,
						"", "Approximation Space")
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

//	ProlongationOperator
	{
		typedef ApproximationSpace<domain_type, dof_distribution_type, algebra_type> approximation_space_type;
		typedef P1ProlongationOperator<approximation_space_type, algebra_type> T;

		std::stringstream ss; ss << "P1ProlongationOperator" << dim << "d";
		reg.add_class_<T, IProlongationOperator<vector_type, vector_type> >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set_approximation_space", &T::set_approximation_space)
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
			.add_method("set_surface_level", &T::set_surface_level)
			.add_method("set_base_level|interactive=false", &T::set_base_level,
						"", "Base Level")
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
		std::stringstream ss; ss << "VTKOutput" << dim << "d";
		reg.add_class_<VTKOutput<function_type> >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("begin_timeseries", &VTKOutput<function_type>::begin_timeseries)
			.add_method("end_timeseries", &VTKOutput<function_type>::end_timeseries)
			.add_method("print", &VTKOutput<function_type>::print);
	}

//	GridFunctionDebugWriter
	{
		typedef GridFunctionDebugWriter<function_type> T;
		typedef IDebugWriter<typename function_type::algebra_type> TBase;
		std::stringstream ss; ss << "GridFunctionDebugWriter" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set_reference_grid_function", &T::set_reference_grid_function)
			.add_method("set_vtk_output", &T::set_vtk_output)
			.add_method("set_conn_viewer_output", &T::set_conn_viewer_output);

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

#ifdef UG_PARALLEL
	//	todo: only temporary
		{
			std::stringstream ss; ss << "OneToManyTests" << dim << "d";

			reg.add_function(ss.str().c_str(),
							&OneToManyTests<function_type>);
		}

	//	todo: only temporary
		{
			std::stringstream ss; ss << "BuildDomainDecompositionLayoutsTest"
								<< dim << "d";

			reg.add_function(ss.str().c_str(),
							&BuildDomainDecompositionLayoutsTest<function_type>);
		}

	//	todo: only temporary
		{
			std::stringstream ss; ss << "TestGridFunctionLayout"
								<< dim << "d";

			reg.add_function(ss.str().c_str(),
							&TestGridFunctionLayout<function_type>);
		}
#endif

	//	PerformTimeStep
		{
			std::stringstream ss; ss << "PerformTimeStep" << dim << "d";
			reg.add_function(ss.str().c_str(), &PerformTimeStep<function_type>, grp.c_str(),
							"Success", "Solver#GridFunction#Discretization#Timesteps#StartTimestep#Time#Timestep Size#Visualization Output#FileName#DoOutput");
		}

	//	ApplyLinearSolver
		{
			std::stringstream ss; ss << "ApplyLinearSolver" << dim << "d";
			reg.add_function(ss.str().c_str(), &ApplyLinearSolver<vector_type>, grp.c_str());
		}

	//	WriteGridToVTK
		{
			std::stringstream ss; ss << "WriteGridFunctionToVTK" << dim << "d";
			reg.add_function(ss.str().c_str(), &WriteGridFunctionToVTK<function_type>, grp.c_str(),
								"Success", "GridFunction#Filename|save-dialog",
								"Saves GridFunction to *.tvk file", "No help");
		}

	//	SaveMatrixForConnectionViewer
		{
			std::stringstream ss; ss << "SaveMatrixForConnectionViewer" << dim << "d";
			reg.add_function(ss.str().c_str(), &SaveMatrixForConnectionViewer<function_type>, grp.c_str());
		}

	//	SaveVectorForConnectionViewer
		{
			std::stringstream ss; ss << "SaveVectorForConnectionViewer" << dim << "d";
			reg.add_function(ss.str().c_str(), &SaveVectorForConnectionViewer<function_type>, grp.c_str());
		}

	//	InterpolateFunction
		{
			std::stringstream ss; ss << "InterpolateFunction" << dim << "d";
			reg.add_function(ss.str().c_str(),
								&InterpolateFunction<function_type>,
								grp.c_str());
		}
}

template <typename TAlgebra, typename TDoFDistribution>
bool RegisterLibDiscretizationInterfaceForAlgebraDomainDedependent(Registry& reg, const char* parentGroup)
{
	typedef TDoFDistribution dof_distribution_type;
	typedef TAlgebra algebra_type;
	typedef typename algebra_type::vector_type vector_type;
	typedef typename algebra_type::matrix_type matrix_type;

	try
	{
	//	get group string
		std::string grp = parentGroup; grp.append("/Discretization");

	//	Domain dependend part 1D
/*		{
			typedef Domain<1, MultiGrid, MGSubsetHandler> domain_type;
			RegisterLibDiscretizationDomainObjects<domain_type, algebra_type, dof_distribution_type>(reg, grp.c_str());
			RegisterLibDiscretizationDomainFunctions<domain_type,  algebra_type, dof_distribution_type>(reg, grp.c_str());
		}
*/
	//	Domain dependend part 2D
		{
			typedef Domain<2, MultiGrid, MGSubsetHandler> domain_type;
			RegisterLibDiscretizationDomainObjects<domain_type, algebra_type, dof_distribution_type>(reg, grp.c_str());
			RegisterLibDiscretizationDomainFunctions<domain_type,  algebra_type, dof_distribution_type>(reg, grp.c_str());
		}

	//	Domain dependend part 3D
/*		{
			typedef Domain<3, MultiGrid, MGSubsetHandler> domain_type;
			RegisterLibDiscretizationDomainObjects<domain_type, algebra_type, dof_distribution_type>(reg, grp.c_str());
			RegisterLibDiscretizationDomainFunctions<domain_type,  algebra_type, dof_distribution_type>(reg, grp.c_str());
		}
*/
	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterLibDiscretizationInterfaceForAlgebra: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}

	return true;
}

bool RegisterDynamicLibDiscretizationInterfaceDomainDependent(Registry& reg, int algebra_type, const char* parentGroup)
{
	bool bReturn = true;

	switch(algebra_type)
	{
	case eCPUAlgebra:		 		bReturn &= RegisterLibDiscretizationInterfaceForAlgebraDomainDedependent<CPUAlgebra, P1ConformDoFDistribution>(reg, parentGroup); break;
//	case eCPUBlockAlgebra2x2: 		bReturn &= RegisterLibDiscretizationInterfaceForAlgebraDomainDedependent<CPUBlockAlgebra<2>, GroupedP1ConformDoFDistribution>(reg, parentGroup); break;
	case eCPUBlockAlgebra3x3: 		bReturn &= RegisterLibDiscretizationInterfaceForAlgebraDomainDedependent<CPUBlockAlgebra<3>, GroupedP1ConformDoFDistribution>(reg, parentGroup); break;
//	case eCPUBlockAlgebra4x4: 		bReturn &= RegisterLibDiscretizationInterfaceForAlgebraDomainDedependent<CPUBlockAlgebra<4>, GroupedP1ConformDoFDistribution>(reg, parentGroup); break;
//	case eCPUVariableBlockAlgebra: 	bReturn &= RegisterLibDiscretizationInterfaceForAlgebraDomainDedependent<CPUVariableBlockAlgebra, GroupedP1ConformDoFDistribution>(reg, parentGroup); break;
	default: UG_ASSERT(0, "Unsupported Algebra Type");
				UG_LOG("Unsupported Algebra Type requested.\n");
				return false;
	}

	return bReturn;
}

}//	end of namespace ug
}//	end of namespace interface
