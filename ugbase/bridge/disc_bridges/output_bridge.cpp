/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
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

#include "../util_overloaded.h"

using namespace std;

namespace ug{
namespace bridge{
namespace Output{

/**
 * \defgroup output_bridge Output Bridge
 * \ingroup disc_bridge
 * \{
 */

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
	VTKOutput<dim> out;
	out.print(filename, domain);
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
	typedef GridFunction<TDomain, TAlgebra> function_type;

//	VTK Output
	{
		typedef VTKOutput<dim> T;
		reg.get_class_<T>()
			.add_method("write_time_pvd", static_cast<void (T::*)(const char*, function_type&)>(&T::write_time_pvd))
			.add_method("write_time_pvd_subset", static_cast<void (T::*)(const char*, function_type&, int)>(&T::write_time_pvd_subset), "", "name # printed grid function # subset index", "", "")
			.add_method("print", static_cast<void (T::*)(const char*, function_type&, int, number, bool)>(&T::print))
			.add_method("print", static_cast<void (T::*)(const char*, function_type&, int, number)>(&T::print))
			.add_method("print", static_cast<void (T::*)(const char*, function_type&, bool)>(&T::print))
			.add_method("print", static_cast<void (T::*)(const char*, function_type&)>(&T::print))
			.add_method("print_subset", static_cast<void (T::*)(const char*, function_type&, int, int, number, bool)>(&T::print_subset))
			.add_method("print_subset", static_cast<void (T::*)(const char*, function_type&, int, int, number)>(&T::print_subset))
			.add_method("print_subsets", static_cast<void (T::*)(const char*, function_type&, const char*, int, number, bool)>(&T::print_subsets))
			.add_method("print_subsets", static_cast<void (T::*)(const char*, function_type&, const char*, int, number)>(&T::print_subsets))
			.add_method("print_subsets", static_cast<void (T::*)(const char*, function_type&, const char*, bool)>(&T::print_subsets))
			.add_method("print_subsets", static_cast<void (T::*)(const char*, function_type&, const char*)>(&T::print_subsets))
			;
			
	}


//	GridFunctionDebugWriter
	{
		typedef GridFunctionDebugWriter<TDomain, TAlgebra> T;
		typedef IDebugWriter<TAlgebra> TBase;
		string name = string("GridFunctionDebugWriter").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >)>("")
			.add_method("reset", &T::reset, "", "")
			.add_method("set_vtk_output", &T::set_vtk_output, "", "bVtkOutput")
			.add_method("set_conn_viewer_output", &T::set_conn_viewer_output, "", "bCVOutput")
			.add_method("set_conn_viewer_indices", &T::set_conn_viewer_indices, "", "bIndicesOutput")
			.add_method("set_print_consistent",  &T::set_print_consistent, "", "printConsistent")
		    .add_method("set_base_dir", &T::set_base_dir, "", "")
		    .add_method("set_grid_level", &T::set_grid_level, "Sets the grid level", "GridLevel")
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
			. ADD_CONSTRUCTOR( (const function_type&) )("gridFunction")
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
							"", "GridFunction#Filename|save-dialog|endings=[\"vtk\"];description=\"VTK-Files\"",
							"Saves GridFunction to *.vtk file", "No help");
	}

//	SaveMatrixForConnectionViewer
	{
		reg.add_function("SaveMatrixForConnectionViewer",
						 &SaveMatrixForConnectionViewer<function_type>, grp, "", "u#A#Filename|save-dialog|endings=[\"mat\"]");
	}

//	SaveVectorForConnectionViewer
	{
		typedef MatrixOperator<matrix_type,	vector_type> matOp;

		reg.add_function("SaveVectorForConnectionViewer",OVERLOADED_FUNCTION_PTR
				(void, SaveVectorForConnectionViewer<function_type>, (function_type& ,const char*) ) ,
				grp, "", "vec#Filename|save-dialog|endings=[\"vec\"]", "save vector as .vec for ConnectionViewer");
		reg.add_function("SaveVectorDiffForConnectionViewer", OVERLOADED_FUNCTION_PTR
				(void, SaveVectorDiffForConnectionViewer<function_type>, (function_type& ,function_type&, const char*) ),
				grp, "", "vecA#vecB#Filename|save-dialog|endings=[\"vec\"]", "compare two vectors a and b and save difference in .vec for ConnectionViewer");
		reg.add_function("SaveVectorForConnectionViewer", OVERLOADED_FUNCTION_PTR
				(void, SaveVectorForConnectionViewer<function_type>, (function_type& , matOp&, const char*) ),
				grp, "", "vec#matrix#Filename|save-dialog|endings=[\"vec\"]", "save vector as .vec for ConnectionViewer, use matrix connections as a \"grid\"");
		reg.add_function("SaveVectorDiffForConnectionViewer", OVERLOADED_FUNCTION_PTR
				(void, SaveVectorDiffForConnectionViewer<function_type>, (function_type& , function_type& , matOp&, const char*) ),
				grp, "", "vecA#veccB#matrix#Filename|save-dialog|endings=[\"vec\"]", "compare two vectors a and b and save difference in .vec for ConnectionViewer, use matrix connections as a \"grid\"");
		reg.add_function("LoadVector", OVERLOADED_FUNCTION_PTR
				(void, LoadVector<function_type>, (function_type& ,const char*) ),
				grp, "", "vec#Filename|save-dialog|endings=[\"vec\"]", "save vector as .vec for ConnectionViewer");
	}

//	SaveVectorCSV
	{
		reg.add_function("SaveVectorCSV",
						 &SaveVectorCSV<function_type>, grp, "", "b#filename|save-dialog");
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

	reg.add_function("SaveDomainToVTK", &SaveDomainToVTK<dim>, grp, "", "filename|save-dialog#Domain");

//	VTK Output
	{
		typedef VTKOutput<dim> T;
		string name = string("VTKOutput").append(suffix);
		reg.add_class_<T>(name, grp)
			.add_constructor()
			.add_method("clear_selection", &T::clear_selection, "", "", "clears the selected output")
			.add_method("select_all", &T::select_all, "", "", "schedules that all components of the passed discrete functions will be written to file")
			.add_method("select", static_cast<void (T::*)(const char*, const char*)>(&T::select), "", "fctName#name",
					"selects a value of a grid function value to be written. example:\nfctName = \"p\"; name = \"pressure\"\nfctNames = \"u,v,w\"; name = \"velocity\"")
			.add_method("select", static_cast<void (T::*)(const vector<string>&, const char*)>(&T::select), "", "fctName#name")
			.add_method("select", static_cast<void (T::*)(SmartPtr<UserData<number, dim> >, const char*)>(&T::select))
			.add_method("select", static_cast<void (T::*)(SmartPtr<UserData<MathVector<dim>, dim> >, const char*)>(&T::select))
			.add_method("select_nodal", static_cast<void (T::*)(const char*, const char*)>(&T::select_nodal), "", "", "selects a nodal value of a grid function value to be written")
			.add_method("select_nodal", static_cast<void (T::*)(const vector<string>&, const char*)>(&T::select_nodal))
			.add_method("select_nodal", static_cast<void (T::*)(SmartPtr<UserData<number, dim> >, const char*)>(&T::select_nodal))
			.add_method("select_nodal", static_cast<void (T::*)(SmartPtr<UserData<MathVector<dim>, dim> >, const char*)>(&T::select_nodal))
			.add_method("select_element", static_cast<void (T::*)(const char*, const char*)>(&T::select_element))
			.add_method("select_element", static_cast<void (T::*)(const vector<string>&, const char*)>(&T::select_element), "", "", "selects a element value of a grid function value to be written")
			.add_method("select_element", static_cast<void (T::*)(SmartPtr<UserData<number, dim> >, const char*)>(&T::select_element))
			.add_method("select_element", static_cast<void (T::*)(SmartPtr<UserData<MathVector<dim>, dim> >, const char*)>(&T::select_element))
			.add_method("set_binary", &T::set_binary, "", "bBinary", "should values be printed in binary (base64 encoded way ) or plain ascii")
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
#ifdef UG_CPU_1
// SaveMatrixToMTX
	{
		typedef MatrixOperator<CPUAlgebra::matrix_type, CPUAlgebra::vector_type> matOp;
// 		reg.add_function( "SaveMatrixToMTX", static_cast<void (*)(const char*, matOp&)>(&SaveMatrixToMTX), grp );
		reg.add_function( "SaveMatrixToMTX", static_cast<void (*)(const char*, matOp&, std::string)>(&SaveMatrixToMTX), grp,
				"", "filename.mtx|save-dialog|endings=[\"mtx\"];description=\"MatrixMarket Files\"#mat#comment", "Save the assembled matrix of a matrix operator to MatrixMarket format");
	}
#endif
}

}; // end Functionality

// end group output_bridge
/// \}

}// end Output

/// \addtogroup output_bridge
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
