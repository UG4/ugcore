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

#include "lib_disc/io/vtkoutput.h"
#include "common/profiler/profiler.h"

using namespace std;

namespace ug {
namespace bridge {

/// small wrapper to write a grid function to vtk
template <typename TGridFunction>
void WriteGridFunctionToVTK(TGridFunction& u, const char* filename)
{
	PROFILE_FUNC();
	VTKOutput<TGridFunction> out;
	out.print(filename, u, true); // TODO: setting of last argument (intended to skip "make consistent", for writing of "raw" data; see (grid_function_util.h')
}

template <typename TDomain, typename TAlgebra>
static void Register__Algebra_Domain(Registry& reg, string parentGroup)
{
//	typedef
	static const int dim = TDomain::dim;
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef ApproximationSpace<TDomain> approximation_space_type;

	typedef GridFunction<TDomain, SurfaceDoFDistribution, TAlgebra> function_type;

//	group string
	stringstream grpSS; grpSS << parentGroup << "/Output";
	string grp = grpSS.str();

//	suffix and tag
	string dimAlgSuffix = GetDomainSuffix<TDomain>();
	dimAlgSuffix.append(GetAlgebraSuffix<TAlgebra>());

	string dimAlgTag = GetDomainTag<TDomain>();
	dimAlgTag.append(GetAlgebraTag<TAlgebra>());

//	VTK Output
	{
		typedef VTKOutput<function_type> T;
		string name = string("VTKOutput").append(dimAlgSuffix);
		reg.add_class_<T>(name, grp)
			.add_constructor()
			.add_method("write_time_pvd", &T::write_time_pvd)
			.add_method("clear_selection", &T::clear_selection)
			.add_method("select_all", &T::select_all)
			.add_method("select_nodal_scalar", &T::select_nodal_scalar)
			.add_method("select_nodal_vector", &T::select_nodal_vector)
			.add_method("print", static_cast<void (T::*)(const char*, function_type&, int, number, bool)>(&T::print))
			.add_method("print", static_cast<void (T::*)(const char*, function_type&, int, number)>(&T::print))
			.add_method("print", static_cast<void (T::*)(const char*, function_type&, bool)>(&T::print))
			.add_method("print", static_cast<void (T::*)(const char*, function_type&)>(&T::print))
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "VTKOutput", dimAlgTag);
	}


//	GridFunctionDebugWriter
	{
		typedef GridFunctionDebugWriter<TDomain, TAlgebra> T;
		typedef IDebugWriter<TAlgebra> TBase;
		string name = string("GridFunctionDebugWriter").append(dimAlgSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >)>("")
			.add_method("reset", &T::reset, "", "")
			.add_method("set_vtk_output", &T::set_vtk_output, "", "vtkOutput")
			.add_method("set_conn_viewer_output", &T::set_conn_viewer_output, "", "cvOutput")
			.add_method("set_print_consistent",  &T::set_print_consistent, "", "printConsistent")
		    .add_method("set_base_dir", &T::set_base_dir, "", "")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GridFunctionDebugWriter", dimAlgTag);
	}

//	GridFunctionPositionProvider
	{
		typedef GridFunctionPositionProvider<function_type> T;
		typedef IPositionProvider<dim> TBase;
		string name = string("GridFunctionPositionProvider").append(dimAlgSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("set_reference_grid_function", &T::set_reference_grid_function, "", "gridFunction")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GridFunctionPositionProvider", dimAlgTag);
	}


	//	GridFunctionVectorWriter
	{
		typedef GridFunctionVectorWriter<function_type, vector_type> T;
		typedef IVectorWriter<vector_type> TBase;
		string name = string("GridFunctionVectorWriter").append(dimAlgSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("set_reference_grid_function", &T::set_reference_grid_function, "", "gridFunction")
			.add_method("set_user_data", &T::set_user_data, "", "userData")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GridFunctionVectorWriter", dimAlgTag);
	}

	// GridFunctionVectorWriterDirichlet0
	{
		typedef GridFunctionVectorWriterDirichlet0<function_type> T;
		typedef IVectorWriter<vector_type> TBase;
		string name = string("GridFunctionVectorWriterDirichlet0").append(dimAlgSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("init", &T::init, "", "postProcess#approxSpace#level")
			.add_method("set_level", &T::set_level, "", "level")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GridFunctionVectorWriterDirichlet0", dimAlgTag);
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
		typedef MatrixOperator<matrix_type,	vector_type> matOp;

		reg.add_function("SaveVectorForConnectionViewer", static_cast<void (*)(function_type& ,const char*)>(&SaveVectorForConnectionViewer<function_type>), grp);
		reg.add_function("SaveVectorForConnectionViewer", static_cast<bool (*)(function_type& , matOp&, const char*)>(&SaveVectorForConnectionViewer<function_type>), grp);
		reg.add_function("SaveVectorForConnectionViewer", static_cast<bool (*)(function_type& , function_type& , matOp&, const char*)>(&SaveVectorForConnectionViewer<function_type>), grp);
	}

//	SaveVectorCSV
	{
		reg.add_function("SaveVectorCSV",
						 &SaveVectorCSV<function_type>, grp);
	}
}

template <typename TAlgebra>
static void Register__Algebra(Registry& reg, string parentGroup)
{
//	get group string
	string grp = parentGroup; grp.append("/Discretization");

	try{
#ifdef UG_DIM_1
		Register__Algebra_Domain<Domain1d, TAlgebra>(reg, grp);
#endif
#ifdef UG_DIM_2
		Register__Algebra_Domain<Domain2d, TAlgebra>(reg, grp);
#endif
#ifdef UG_DIM_3
		Register__Algebra_Domain<Domain3d, TAlgebra>(reg, grp);
#endif
	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterOutput: "
				"Registration failed (using name " << ex.name << ").\n");
		UG_THROW_FATAL("Registration failed.");
	}
}

bool RegisterOutput(Registry& reg, string grp)
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
	return true;
}

}//	end of namespace ug
}//	end of namespace interface
