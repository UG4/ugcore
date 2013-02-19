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
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"

// lib_disc includes
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/function_spaces/grid_function_util.h"
#include "lib_disc/function_spaces/approximation_space.h"

#include "lib_disc/io/vtkoutput.h"
#include "common/profiler/profiler.h"

using namespace std;

namespace ug{
namespace bridge{
namespace Output{

/// small wrapper to write a grid function to vtk
template <typename TGridFunction>
void WriteGridFunctionToVTK(TGridFunction& u, const char* filename)
{
	PROFILE_FUNC();
	VTKOutput<TGridFunction::dim> out;
	out.print(filename, u, true); // TODO: setting of last argument (intended to skip "make consistent", for writing of "raw" data; see (grid_function_util.h')
}

template <int dim>
void SaveDomainToVTK(const char* filename, Domain<dim>& domain)
{
	VTKOutput<dim>::print(filename, domain);
}

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

/**
 * Function called for the registration of Domain and Algebra dependent parts.
 * All Functions and Classes depending on both Domain and Algebra
 * are to be placed here when registering. The method is called for all
 * available Domain and Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain, typename TAlgebra>
static void DomainAlgebra(Registry& reg, string grp)
{
	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

//	typedef
	static const int dim = TDomain::dim;
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef ApproximationSpace<TDomain> approximation_space_type;
	typedef GridFunction<TDomain, TAlgebra> function_type;

//	VTK Output
	{
		typedef VTKOutput<dim> T;
		reg.get_class_<T>()
			.add_method("write_time_pvd", static_cast<void (T::*)(const char*, function_type&)>(&T::write_time_pvd))
			.add_method("print", static_cast<void (T::*)(const char*, function_type&, int, number, bool)>(&T::print))
			.add_method("print", static_cast<void (T::*)(const char*, function_type&, int, number)>(&T::print))
			.add_method("print", static_cast<void (T::*)(const char*, function_type&, bool)>(&T::print))
			.add_method("print", static_cast<void (T::*)(const char*, function_type&)>(&T::print));
	}


//	GridFunctionDebugWriter
	{
		typedef GridFunctionDebugWriter<TDomain, TAlgebra> T;
		typedef IDebugWriter<TAlgebra> TBase;
		string name = string("GridFunctionDebugWriter").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >)>("")
			.add_method("reset", &T::reset, "", "")
			.add_method("set_vtk_output", &T::set_vtk_output, "", "vtkOutput")
			.add_method("set_conn_viewer_output", &T::set_conn_viewer_output, "", "cvOutput")
			.add_method("set_print_consistent",  &T::set_print_consistent, "", "printConsistent")
		    .add_method("set_base_dir", &T::set_base_dir, "", "")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GridFunctionDebugWriter", tag);
	}

//	GridFunctionPositionProvider
	{
		typedef GridFunctionPositionProvider<function_type> T;
		typedef IPositionProvider<dim> TBase;
		string name = string("GridFunctionPositionProvider").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("set_reference_grid_function", &T::set_reference_grid_function, "", "gridFunction")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GridFunctionPositionProvider", tag);
	}


	//	GridFunctionVectorWriter
	{
		typedef GridFunctionVectorWriter<function_type, vector_type> T;
		typedef IVectorWriter<vector_type> TBase;
		string name = string("GridFunctionVectorWriter").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("set_reference_grid_function", &T::set_reference_grid_function, "", "gridFunction")
			.add_method("set_user_data", &T::set_user_data, "", "userData")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GridFunctionVectorWriter", tag);
	}

	// GridFunctionVectorWriterDirichlet0
	{
		typedef GridFunctionVectorWriterDirichlet0<function_type> T;
		typedef IVectorWriter<vector_type> TBase;
		string name = string("GridFunctionVectorWriterDirichlet0").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("init", &T::init, "", "postProcess#approxSpace#level")
			.add_method("set_level", &T::set_level, "", "level")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GridFunctionVectorWriterDirichlet0", tag);
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
		reg.add_function("SaveVectorForConnectionViewer", static_cast<void (*)(function_type& , matOp&, const char*)>(&SaveVectorForConnectionViewer<function_type>), grp);
		reg.add_function("SaveVectorForConnectionViewer", static_cast<void (*)(function_type& , function_type& , matOp&, const char*)>(&SaveVectorForConnectionViewer<function_type>), grp);
	}

//	SaveVectorCSV
	{
		reg.add_function("SaveVectorCSV",
						 &SaveVectorCSV<function_type>, grp);
	}
}


/**
 * Function called for the registration of Dimension dependent parts.
 * All Functions and Classes depending on the Dimension
 * are to be placed here when registering. The method is called for all
 * available Dimension types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <int dim>
static void Dimension(Registry& reg, string grp)
{
	string suffix = GetDimensionSuffix<dim>();
	string tag = GetDimensionTag<dim>();

	reg.add_function("SaveDomainToVTK", &SaveDomainToVTK<dim>);

//	VTK Output
	{
		typedef VTKOutput<dim> T;
		string name = string("VTKOutput").append(suffix);
		reg.add_class_<T>(name, grp)
			.add_constructor()
			.add_method("clear_selection", &T::clear_selection)
			.add_method("select_all", &T::select_all)
			.add_method("select_nodal", static_cast<void (T::*)(const char*, const char*)>(&T::select_nodal))
			.add_method("select_nodal", static_cast<void (T::*)(const vector<string>&, const char*)>(&T::select_nodal))
			.add_method("select_nodal", static_cast<void (T::*)(SmartPtr<UserData<number, dim> >, const char*)>(&T::select_nodal))
			.add_method("select_nodal", static_cast<void (T::*)(SmartPtr<UserData<MathVector<dim>, dim> >, const char*)>(&T::select_nodal))
			.add_method("select_element", static_cast<void (T::*)(const char*, const char*)>(&T::select_element))
			.add_method("select_element", static_cast<void (T::*)(const vector<string>&, const char*)>(&T::select_element))
			.add_method("select_element", static_cast<void (T::*)(SmartPtr<UserData<number, dim> >, const char*)>(&T::select_element))
			.add_method("select_element", static_cast<void (T::*)(SmartPtr<UserData<MathVector<dim>, dim> >, const char*)>(&T::select_element))
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "VTKOutput", tag);
	}
}

/**
 * Function called for the registration of Domain and Algebra independent parts.
 * All Functions and Classes not depending on Domain and Algebra
 * are to be placed here when registering.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
static void Common(Registry& reg, string grp)
{
}

}; // end Functionality
}// end Output

void RegisterBridge_Output(Registry& reg, string grp)
{
	grp.append("/Discretization/Output");
	typedef Output::Functionality Functionality;

	try{
		RegisterCommon<Functionality>(reg,grp);
		RegisterDimensionDependent<Functionality>(reg,grp);
		RegisterDomainAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}


}//	end of namespace bridge
}//	end of namespace ug
