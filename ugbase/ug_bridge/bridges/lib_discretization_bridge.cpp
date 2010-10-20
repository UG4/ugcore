/*
 * lib_discretization_bridge.cpp
 *
 *  Created on: 22.09.2010
 *      Author: andreasvogel
 */

#include "../ug_bridge.h"
#include "lib_discretization/lib_discretization.h"
#include "user_data/user_data.h"
#include <iostream>
#include <sstream>

namespace ug
{
namespace bridge
{

//////////////////////////////////
// Some global functions
//////////////////////////////////

template <typename TDomain>
bool DistributeDomain(TDomain& domainOut)
{
#ifdef UG_PARALLEL
//	typedefs
	typedef typename TDomain::subset_handler_type subset_handler_type;
	typedef typename TDomain::distributed_grid_manager_type distributed_grid_manager_type;

//	get distributed grid manager
	distributed_grid_manager_type* pDistGridMgr = domainOut.get_distributed_grid_manager();

//	check that manager exists
	if(!pDistGridMgr)
	{
		UG_LOG("DistibuteDomain: Cannot find Distributed Grid Manager.\n");
		return false;
	}
	distributed_grid_manager_type& distGridMgrOut = *pDistGridMgr;

//	get subset handler
	subset_handler_type& sh = domainOut.get_subset_handler();

//	get number of processes
	const int numProcs = pcl::GetNumProcesses();
	if(numProcs == 1) return true;

//	check, that grid is a multigrid
	MultiGrid* pMG = dynamic_cast<MultiGrid*>(distGridMgrOut.get_assigned_grid());
	if(pMG == NULL)
	{
		UG_LOG("DistibuteDomain: MultiGrid-Domain required in current implementation.\n");
		return false;
	}
	MultiGrid& mg = *pMG;

//	get Grid Layout
	GridLayoutMap& glmOut = distGridMgrOut.grid_layout_map();

//	make sure that each grid has a position attachment - even if no data
//	will be received.
	typedef typename TDomain::position_attachment_type position_attachment_type;
	position_attachment_type& domPosition = domainOut.get_position_attachment();
	bool tmpPosAttachment = false;
	if(!mg.has_vertex_attachment(aPosition))
	{
	// convert to 3d positions (FVGeometry depends on PositionCoordinates)
       mg.attach_to_vertices(aPosition);
       ConvertMathVectorAttachmentValues<VertexBase>(mg, domPosition, aPosition);

       UG_LOG("DistributeDomain: temporarily adding Position Attachment.\n");
       tmpPosAttachment = true;
	}

//	AdjustFunctions
	FuncAdjustGrid funcAdjustGrid = DefaultAdjustGrid;
	FuncPartitionGrid funcPartitionGrid = PartitionGrid_Bisection;

//	process 0 loads and distributes the grid. The others receive it.
	if(pcl::GetProcRank() == 0)
	{
	//	adjust the grid
		funcAdjustGrid(mg, sh);

		LOG("  Performing load balancing ... \n");

	//	perform load-balancing
	//TODO: if grid partitioning fails, the whole process should be aborted.
	//		this has to be communicated to the other processes.
		SubsetHandler shPartition(mg);
		if(!funcPartitionGrid(shPartition, mg, sh, numProcs)){
			UG_LOG("  grid partitioning failed. proceeding anyway...\n");
		}

	//	get min and max num elements
		int maxElems = 0;
		int minElems = 0;
		if(mg.num<Volume>() > 0){
			minElems = mg.num<Volume>();
			for(int i = 0; i < shPartition.num_subsets(); ++i){
				minElems = min(minElems, (int)shPartition.num<Volume>(i));
				maxElems = max(maxElems, (int)shPartition.num<Volume>(i));
			}
		}
		else if(mg.num<Face>() > 0){
			minElems = mg.num<Face>();
			for(int i = 0; i < shPartition.num_subsets(); ++i){
				minElems = min(minElems, (int)shPartition.num<Face>(i));
				maxElems = max(maxElems, (int)shPartition.num<Face>(i));
			}
		}

		LOG("  Element Distribution - min: " << minElems << ", max: " << maxElems << endl);

		const char* partitionMapFileName = "partitionMap.obj";
		LOG("saving partition map to " << partitionMapFileName << endl);
		SaveGridToFile(mg, partitionMapFileName, shPartition);

		//	distribute the grid.
		LOG("  Distributing grid... ");
		if(DistributeGrid_KeepSrcGrid(mg, sh, glmOut, shPartition, 0))
		{
			UG_LOG("done!\n");
		}
		else{
			UG_LOG("failed\n");
		}
	}
	else
	{
	//	a grid will only be received, if the process-rank is smaller than numProcs
		if(pcl::GetProcRank() < numProcs){
			if(!ReceiveGrid(mg, sh, glmOut, 0, true)){
				UG_LOG("  ReceiveGrid failed on process " << pcl::GetProcRank() <<
				". Aborting...\n");
				return false;
			}
		}
	}

	if(tmpPosAttachment)
	{
	// convert to 3d positions (FVGeometry depends on PositionCoordinates)
       ConvertMathVectorAttachmentValues<VertexBase>(mg, aPosition, domPosition);
       mg.detach_from_vertices(aPosition);

       UG_LOG("DistributeDomain: removing temporary Position Attachment.\n");
 	}

//	tell the distGridMgr that the associated layout changed.
	distGridMgrOut.grid_layouts_changed(true);
#endif

//	in serial case: do nothing
	return true;
}

template <typename TDomain>
void GlobalRefineParallelDomain(TDomain& domain)
{
#ifdef UG_PARALLEL
//	get distributed grid manager
	typedef typename TDomain::distributed_grid_manager_type distributed_grid_manager_type;
	distributed_grid_manager_type* pDistGridMgr = domain.get_distributed_grid_manager();

//	check that manager exists
	if(!pDistGridMgr)
	{
		UG_LOG("GlobalRefineParallelDomain: Cannot find Distributed Grid Manager.\n");
		throw(int(1));
	}
	distributed_grid_manager_type& distGridMgr = *pDistGridMgr;

//	create Refiner
	ParallelGlobalMultiGridRefiner refiner(distGridMgr);
#else
	GlobalMultiGridRefiner refiner();
	refiner.assign_grid(domain.get_grid());
#endif

//	perform refinement.
	refiner.refine();
}

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
	typedef CPUAlgebra algebra_type;
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
			.add_method("get_function_pattern", &T::get_function_pattern)
			.add_method("assign_function_pattern", &T::assign_function_pattern)
			.add_method("get_surface_dof_distribution", &T::get_surface_dof_distribution)
			.add_method("create_surface_function", &T::create_surface_function);
	}

//	DirichletBNDValues
	{
		typedef P1DirichletBoundary<domain_type, dof_distribution_type, algebra_type> T;
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
			.add_method("add_boundary_value", (bool (T::*)(IBoundaryNumberProvider<dim>&, const char*, const char*))&T::add_boundary_value);
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

//	DistributeDomain
	{
		stringstream ss; ss << "DistributeDomain" << dim << "d";
		reg.add_function(ss.str().c_str(), &DistributeDomain<domain_type>);
	}

//	GlobalRefineParallelDomain
	{
		stringstream ss; ss << "GlobalRefineParallelDomain" << dim << "d";
		reg.add_function(ss.str().c_str(), &GlobalRefineParallelDomain<domain_type>);
	}

}


void RegisterLibDiscretizationInterface(Registry& reg, const char* parentGroup)
{
	//	get group string
		std::stringstream groupString; groupString << parentGroup << "/Discretization";
		const char* grp = groupString.str().c_str();

	//	todo: generalize
		typedef P1ConformDoFDistribution dof_distribution_type;
		typedef CPUAlgebra algebra_type;

	//	FunctionPattern (Abstract Base Class)
		reg.add_class_<FunctionPattern>("FunctionPattern", grp);

	//	P1ConformFunctionPattern
		{
		typedef P1ConformFunctionPattern T;
		reg.add_class_<T, FunctionPattern>("P1ConformFunctionPattern", grp)
			.add_constructor()
			.add_method("set_subset_handler", &T::set_subset_handler)
			.add_method("lock", &T::lock);
		}

	//	P1ConformDoFDistribution
		{
			typedef P1ConformDoFDistribution T;
			reg.add_class_<T>("P1ConformDoFDistribution", grp);
		}

	//  Add discrete function to pattern
		reg.add_function("AddP1Function", &AddP1Function, grp);

	//	DomainDiscretization
		{
			typedef DomainDiscretization<dof_distribution_type, algebra_type> T;

			reg.add_class_<IAssemble<dof_distribution_type, algebra_type> >("IAssemble", grp);
			reg.add_class_<IDomainDiscretization<dof_distribution_type, algebra_type>,
							IAssemble<dof_distribution_type, algebra_type> >("IDomainDiscretization", grp);

			reg.add_class_<T, IDomainDiscretization<dof_distribution_type, algebra_type> >("DomainDiscretization", grp)
				.add_constructor()
				.add_method("add_post_process", &T::add_post_process)
				.add_method("add_elem_disc", (bool (T::*)(IElemDisc<algebra_type>&)) &T::add_elem_disc);
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
		reg.add_class_<IPostProcess<dof_distribution_type, algebra_type> >("IPostProcess", grp);

	//	UserFunction
		{
		//	Density - Driven - Flow
			reg.add_class_<IDensityDrivenFlowUserFunction<2> >("IDensityDrivenFlowUserFunction2d", grp);
		}

	//	Elem Discs
		{
		//	Base class
			typedef IElemDisc<algebra_type> T;
			reg.add_class_<T>("IElemDisc", grp)
				.add_method("set_pattern", &T::set_pattern)
				.add_method("set_functions", (bool (T::*)(const char*))&T::set_functions)
				.add_method("set_subsets",  (bool (T::*)(const char*))&T::set_subsets);
		}

	//	FunctionGroup
		{
			reg.add_class_<FunctionGroup>("FunctionGroup", grp)
				.add_constructor()
				.add_method("clear", &FunctionGroup::clear)
				.add_method("set_function_pattern", &FunctionGroup::set_function_pattern)
				.add_method("add_function", (bool (FunctionGroup::*)(const char*))&FunctionGroup::add_function);
		}

	//	AssembledLinearOperator
		{
			typedef AssembledLinearOperator<dof_distribution_type, algebra_type> T;

			reg.add_class_<T, IMatrixOperator<algebra_type::vector_type, algebra_type::vector_type, algebra_type::matrix_type> >
							("AssembledLinearOperator", grp)
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
							("AssembledOperator", grp)
				.add_constructor()
				.add_method("set_discretization", &T::set_discretization)
				.add_method("set_dof_distribution", &T::set_dof_distribution)
				.add_method("init", &T::init);
		}

	//	StandardLineSearch
		{
			typedef StandardLineSearch<algebra_type::vector_type> T;

			reg.add_class_<ILineSearch<algebra_type::vector_type> >("ILineSearch", grp);

			reg.add_class_<StandardLineSearch<algebra_type::vector_type>,
							ILineSearch<algebra_type::vector_type> >("StandardLineSearch", grp)
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

			reg.add_class_<T, IOperatorInverse<algebra_type::vector_type, algebra_type::vector_type> >("NewtonSolver", grp)
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
