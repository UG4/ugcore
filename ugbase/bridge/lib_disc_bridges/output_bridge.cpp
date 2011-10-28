/*
 * output_bridge.cpp
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
#include "lib_disc/function_spaces/grid_function_util.h"
#include "lib_disc/function_spaces/approximation_space.h"
#include "lib_disc/dof_manager/conform/conform.h"
#include "lib_disc/dof_manager/p1conform/p1conform.h"

#include "lib_disc/io/vtkoutput.h"

using namespace std;

namespace ug {
namespace bridge {

/// small wrapper to write a grid function to vtk
template <typename TGridFunction>
void WriteGridFunctionToVTK(TGridFunction& u, const char* filename)
{
	VTKOutput<TGridFunction> out;
	out.print(filename, u, true); // TODO: setting of last argument (intended to skip "make consistent", for writing of "raw" data; see (grid_function_util.h')
}

template <typename TDomain, typename TAlgebra, typename TDoFDistribution>
static void Register__Algebra_DoFDistribution_Domain(Registry& reg, string parentGroup)
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
	stringstream grpSS; grpSS << parentGroup << "/Output";
	string grp = grpSS.str();

//	suffix and tag
	string dimAlgDDSuffix = GetDomainSuffix<TDomain>();
	dimAlgDDSuffix.append(GetAlgebraSuffix<TAlgebra>());
	dimAlgDDSuffix.append(GetDoFDistributionSuffix<TDoFDistribution>());

	string dimAlgDDTag = GetDomainTag<TDomain>();
	dimAlgDDTag.append(GetAlgebraTag<TAlgebra>());
	dimAlgDDTag.append(GetDoFDistributionTag<TDoFDistribution>());

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
			.add_method("print", static_cast<bool (T::*)(const char*, function_type&, int, number, bool)>(&T::print))
			.add_method("print", static_cast<bool (T::*)(const char*, function_type&, bool)>(&T::print));
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
			.add_method("set_conn_viewer_output", &T::set_conn_viewer_output, "", "cvOutput")
			.add_method("set_print_raw_data",  &T::set_print_raw_data, "", "printRawData");
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

//	WriteGridToVTK
	{
		reg.add_function("WriteGridFunctionToVTK",
						 &WriteGridFunctionToVTK<function_type>, grp,
							"", "GridFunction#Filename|save-dialog",
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

	//	SaveVectorCSV
	/*{
		reg.add_function("SaveVectorCSV",
						 &SaveVectorCSV<function_type>, grp);
	}*/

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

bool RegisterOutput(Registry& reg, string grp)
{
	bool bReturn = true;
#ifdef UG_CPU_1
	bReturn &= Register__Algebra<CPUAlgebra>(reg, grp);
#endif
#ifdef UG_CPU_2
	bReturn &= Register__Algebra<CPUBlockAlgebra<2> >(reg, grp);
#endif
#ifdef UG_CPU_3
	bReturn &= Register__Algebra<CPUBlockAlgebra<3> >(reg, grp);
#endif
#ifdef UG_CPU_4
	bReturn &= Register__Algebra<CPUBlockAlgebra<4> >(reg, grp);
#endif
#ifdef UG_CPU_VAR
	bReturn &= Register__Algebra<CPUVariableBlockAlgebra >(reg, grp);
#endif
	return bReturn;
}

}//	end of namespace ug
}//	end of namespace interface
