/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__IO__VTKOUTPUT__
#define __H__UG__LIB_DISC__IO__VTKOUTPUT__

// extern libraries
#include <string>
#include <vector>
#include <map>

// other ug modules
#include "common/util/string_util.h"
#include "common/util/base64_file_writer.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/domain.h"
#include "lib_disc/spatial_disc/user_data/user_data.h"

#include "lib_algebra/ordering_strategies/algorithms/lua_ordering.h"

namespace ug{
// todo this should avoid refactoring all the signatures of VTKOutput, remove it later
typedef Base64FileWriter VTKFileWriter;

template <typename T>
struct IteratorProvider
{
	template <typename TElem>
	struct traits
	{
		typedef typename T::template traits<TElem>::iterator iterator;
		typedef typename T::template traits<TElem>::const_iterator const_iterator;
	};

	template <typename TElem>
	static typename traits<TElem>::const_iterator begin(const T& provider, int si)
	{
		if(si < 0) return provider.template begin<TElem>();
		else return provider.template begin<TElem>(si);
	}

	template <typename TElem>
	static typename traits<TElem>::const_iterator end(const T& provider, int si)
	{
		if(si < 0) return provider.template end<TElem>();
		else return provider.template end<TElem>(si);
	}
};

template <>
struct IteratorProvider<MGSubsetHandler>
{
	private:
	typedef MGSubsetHandler T;

	public:
	template <typename TElem>
	struct traits
	{
		typedef typename T::template traits<TElem>::iterator iterator;
		typedef typename T::template traits<TElem>::const_iterator const_iterator;
	};

	template <typename TElem>
	static typename traits<TElem>::const_iterator begin(const T& provider, int si)
	{
		const int lev = provider.num_levels() - 1;
		if(si < 0) return provider.multi_grid()->template begin<TElem>(lev);
		else return provider.template begin<TElem>(si, lev);
	}

	template <typename TElem>
	static typename traits<TElem>::const_iterator end(const T& provider, int si)
	{
		const int lev = provider.num_levels() - 1;
		if(si < 0) return provider.multi_grid()->template end<TElem>(lev);
		else return provider.template end<TElem>(si, lev);
	}
};


/// output writer to the VTK file format
/**
 * This class can be used to write grid and data associated with the grid to
 * the vtk file format.
 *
 * The produced files are:
 *
 * 1)) TIME DEPENDENT
 * 1a) In case of functions that are only defined for some subsets
 * - filename_pXXXX_sXXXX_tXXXX.vtu  (piece of grid for proc, subset, timestep)
 * - filename_sXXXX_tXXXX.pvtu 	     (piece of grid for subset, timestep)
 * - filename_tXXXX.pvd              (piece of grid for timestep)
 * - filename.pvd                    (group of timeseries)
 *
 * 1b) In case that all functions are defined globally
 * - filename_p0000_t0000.vtu		(piece of grid for proc, timestep)
 * - filename_t0000.pvtu			(piece of grid for timestep)
 *
 * 2)) TIME INDEPENDENT
 * 2a) In case of functions that are only defined for some subsets
 * - filename_pXXXX_sXXXX.vtu  (piece of grid for proc, subset)
 * - filename_sXXXX.pvtu 	   (piece of grid for subset)
 * - filename.pvd              (group of subsets)
 *
 * 2b) In case that all functions are defined globally
 * - filename_p0000.vtu		(piece of grid for proc)
 * - filename.pvtu			(group of procs)
 *
 * In case that the simulation is run in serial (or in parallel environment but
 * with only one process), the *.pvtu are not written and the *.vtu file
 * contains all information.
 *
 * ATTENTION: This class uses heavily the mark-function of the grid.
 * 			  Do not use any member function while having called begin_mark()
 *
 * \tparam		dim 	world dimension
 */
template <int TDim>
class VTKOutput
{
		public:
	///	clears the selected output
		void clear_selection() {m_vSymbFct.clear(); m_vSymbFctNodal.clear(); m_vSymbFctElem.clear();}

	///	clears the selected output
		void clear_data_selection()
		{
			m_vScalarNodalData.clear();
			m_vScalarElemData.clear();

			m_vVectorNodalData.clear();
			m_vVectorElemData.clear();
		}

	///	schedules that all components of the passed discrete functions will be written to file
		void select_all(bool bAll) {m_bSelectAll = bAll;}

	///	selects a value of a grid function value to be written
	/**
	 * This function schedules the component(s) passed by symbolic name(s) of an
	 * approximation space to be written to the vtk file under a specified name.
	 * The type of the data (nodal/element) will be determined based on the
	 * trial space of the components (i.e if continuous available or not)
	 *
	 * If more than one component is passed, the data will be interpreted as a
	 *  vector and #dim arguments must be passed.
	 *
	 * example: fctName = "p"; name = "pressure"
	 * example: fctNames = "u,v,w"; name = "velocity"
	 *
	 * \param[in]	fctName		symbolic name of component
	 * \param[in]	name		name that will appear in the vtk file for the data
	 */
	/// \{
		void select(const char* fctName, const char* name);
		void select(const std::vector<std::string>& vFct, const char* name);
	/// \}

	///	selects a nodal value of a grid function value to be written
	/**
	 * This function schedules the component(s) passed by symbolic name(s) of an
	 * approximation space to be
	 * written to the vtk file under a specified name. Note, that for the
	 * trial space of the component an evaluation of the data at the nodes
	 * must be available (continuous). If more than one component is passed, the
	 * data will be interpreted as a vector and #dim arguments must be passed.
	 *
	 * example: fctName = "p"; name = "pressure"
	 * example: fctNames = "u,v,w"; name = "velocity"
	 *
	 * \param[in]	fctName		symbolic name of component
	 * \param[in]	name		name that will appear in the vtk file for the data
	 */
	/// \{
		void select_nodal(const char* fctName, const char* name);
		void select_nodal(const std::vector<std::string>& vFct, const char* name);
	/// \}

	///	selects a element value of a grid function value to be written
	/**
	 * This function schedules the component(s) passed by symbolic name(s) of an
	 * approximation space to be written to the vtk file under a specified name.
	 * If more than one component is passed, the
	 * data will be interpreted as a vector and #dim arguments must be passed.
	 *
	 * example: fctName = "p"; name = "pressure"
	 * example: fctNames = "u,v,w"; name = "velocity"
	 *
	 * \param[in]	fctName		symbolic name of component
	 * \param[in]	name		name that will appear in the vtk file for the data
	 */
	/// \{
		void select_element(const char* fctName, const char* name);
		void select_element(const std::vector<std::string>& vFct, const char* name);
	/// \}

	///	selects a data value to be written
	/**
	 * This function schedules a user data to be written to the vtk file under
	 * a specified name.
	 * The type of the data (nodal/element) will be determined based on the
	 * trial space of the components (i.e if continuous available or not)
	 *
	 * \param[in]	spData		data to be written
	 * \param[in]	name		name that will appear in the vtk file for the data
	 */
	/// \{
		void select(SmartPtr<UserData<number, TDim> > spData, const char* name);
		void select(SmartPtr<UserData<MathVector<TDim>, TDim> > spData, const char* name);
	/// \}

	///	selects a nodal data value to be written
	/**
	 * This function schedules a user data to be written to the vtk file under
	 * a specified name. Note, that for the the data at the nodes must be
	 * available (continuous).
	 *
	 * \param[in]	spData		data to be written
	 * \param[in]	name		name that will appear in the vtk file for the data
	 */
	/// \{
		void select_nodal(SmartPtr<UserData<number, TDim> > spData, const char* name);
		void select_nodal(SmartPtr<UserData<MathVector<TDim>, TDim> > spData, const char* name);
	/// \}

	///	selects a element data value to be written
	/**
	 * This function schedules a user data to be written to the vtk file under
	 * a specified name.
	 *
	 * \param[in]	spData		data to be written
	 * \param[in]	name		name that will appear in the vtk file for the data
	 */
	/// \{
		void select_element(SmartPtr<UserData<number, TDim> > spData, const char* name);
		void select_element(SmartPtr<UserData<MathVector<TDim>, TDim> > spData, const char* name);
	/// \}

	/**
	 * This function writes the grid to a *.vtu file. It calls the function
	 * print_subset for each subset if necessary and produces a grouping
	 * *.pvd file for paraview. If all subsets have to be written, then only
	 * a single file is produced containing the whole grid.
	 *
	 * \param[in]	filename		filename for produced files
	 * \param[in]	u				grid function
	 * \param[in]	step 			time step counter (-1 indicates stationary case)
	 * \param[in]	time			time point corresponding to timestep
	 * \param[in]	makeConsistent	flag if data shall be made consistent before printing
	 */
		template <typename TFunction>
		void print(const char*  filename, TFunction& u,
		           int step, number time,
				   bool makeConsistent);

	/**	Calls print with makeConsistent enabled.*/
		template <typename TFunction>
		void print(const char*  filename, TFunction& u,
		           int step, number time)
		{
			return print(filename, u, step, time, true);
		}

	/**
	 * This function simply calls 'print' using step = -1 to indicate the
	 * stationary case. It is intended to write time independent data.
	 *
	 * \param[in]	filename		filename for produced files
	 * \param[in]	u				grid function
	 * \param[in]	makeConsistent	flag if data shall be made consistent before printing
	 */
		template <typename TFunction>
		void print(const char*  filename, TFunction& u,
				   bool makeConsistent)
		{
			print(filename, u, -1, 0.0, makeConsistent);
		}

	/**	Calls print with makeConsistent enabled.*/
		template <typename TFunction>
		void print(const char*  filename, TFunction& u)
		{
			return print(filename, u, true);
		}

	///	prints the domain to file
		void print(const char* filename, Domain<TDim>& domain);

	/**
	 * This function writes the subset si of the grid (or the whole grid if
	 * si < 0) to the file "filename.vtu".
	 *
	 * If step >= 0 is passed, this indicates that a time step is written and
	 * to the filename a "*_pXXXX*.vtu" is added, indicating the timestep.
	 *
	 * If the computation is done in parallel and the number of Processes is
	 * greater than one, a "*_pXXXX*.vtu" is added, indicating the process. Then,
	 * in addition a "filename*.pvtu" file is written, grouping all *vtu files
	 * from different processes.
	 *
	 * If only a part of the grid, i.e. a subset, is written, the to the filename
	 * a "*_sXXXX*.vtu" is added, to indicate the subset.
	 *
	 * \param[in, out]	filename 		base name for output file(s)
	 * \param[in]		u				grid function
	 * \param[in]		si				Subset (si < 0 indicates whole grid)
	 * \param[in]		step			counter for timestep (-1 means stationary)
	 * \param[in]		time			time point of timestep
	 * \param[in]	makeConsistent	flag if data shall be made consistent before printing
	 */
		template <typename TFunction>
		void print_subset(const char* filename, TFunction& u,
		                  int si, int step, number time,
						  bool makeConsistent);

	/**	Calls print_subset with makeConsistent enabled.*/
		template <typename TFunction>
		void print_subset(const char* filename, TFunction& u,
		                  int si, int step = -1, number time = 0.0)
		{
			return print_subset(filename, u, si, step, time, true);
		}

	/**
	 * This function writes the given subsets ssGrp of the grid "filename.vtu".
	 *
	 * If step >= 0 is passed, this indicates that a time step is written and
	 * to the filename a "*_pXXXX*.vtu" is added, indicating the timestep.
	 *
	 * If the computation is done in parallel and the number of Processes is
	 * greater than one, a "*_pXXXX*.vtu" is added, indicating the process. Then,
	 * in addition a "filename*.pvtu" file is written, grouping all *vtu files
	 * from different processes.
	 *
	 * \param[in, out]	filename 		base name for output file(s)
	 * \param[in]		u				grid function
	 * \param[in]		ssGrp			subset group
	 * \param[in]		step			counter for timestep (-1 means stationary)
	 * \param[in]		time			time point of timestep
	 * \param[in]	makeConsistent	flag if data shall be made consistent before printing
	 */
		template <typename TFunction>
		void print_subsets(const char* filename, TFunction& u,
		                  SubsetGroup& ssGrp, int step = -1, number time = 0.0,
						  bool makeConsistent = true);

	/**
	 * This function writes the given subsets ssGrp of the grid "filename.vtu".
	 *
	 * If step >= 0 is passed, this indicates that a time step is written and
	 * to the filename a "*_pXXXX*.vtu" is added, indicating the timestep.
	 *
	 * If the computation is done in parallel and the number of Processes is
	 * greater than one, a "*_pXXXX*.vtu" is added, indicating the process. Then,
	 * in addition a "filename*.pvtu" file is written, grouping all *vtu files
	 * from different processes.
	 *
	 * \param[in, out]	filename 		base name for output file(s)
	 * \param[in]		u				grid function
	 * \param[in]		ssNames			subset names
	 * \param[in]		step			counter for timestep (-1 means stationary)
	 * \param[in]		time			time point of timestep
	 * \param[in]	makeConsistent	flag if data shall be made consistent before printing
	 */
		template <typename TFunction>
		void print_subsets(const char* filename, TFunction& u,
		                  const char* ssNames, int step = -1, number time = 0.0,
						  bool makeConsistent = true);

	/**	Calls print_subsets with makeConsistent enabled.*/
		template <typename TFunction>
		void print_subsets(const char*  filename, TFunction& u, const char* ssNames,
		           int step, number time)
		{
			return print_subsets(filename, u, ssNames, step, time, true);
		}

	/**
	 * This function simply calls 'print_subsets' using step = -1 to indicate the
	 * stationary case. It is intended to write time independent data.
	 *
	 * \param[in]	filename		filename for produced files
	 * \param[in]	u				grid function
	 * \param[in]	ssNames			subset names
	 * \param[in]	makeConsistent	flag if data shall be made consistent before printing
	 */
		template <typename TFunction>
		void print_subsets(const char*  filename, TFunction& u, const char* ssNames,
				   bool makeConsistent)
		{
			print_subsets(filename, u, ssNames, -1, 0.0, makeConsistent);
		}

	/**	Calls print_subsets with makeConsistent enabled.*/
		template <typename TFunction>
		void print_subsets(const char*  filename, TFunction& u, const char* ssNames)
		{
			return print_subsets(filename, u, ssNames, true);
		}

	/**
	 * When a time series has been computed, this function can be used to procduce
	 * a grouping *.pvd file for paraview visualization.
	 *
	 * \param[in]		filename		filename used in time series
	 * \param[in]		u				grid function
	 */
		template <typename TFunction>
		void write_time_pvd(const char* filename, TFunction& u);

	/**
	 * When a time series has been computed, this function can be used to procduce
	 * a grouping *_processwise.pvd file for paraview visualization.
	 *
	 * \param[in]		filename		filename used in time series
	 * \param[in]		u				grid function
	 */
		template <typename TFunction>
		void write_time_processwise_pvd(const char* filename, TFunction& u);

	/**
	 * When a time series has been computed, this function can be used to procduce
	 * a grouping *.pvd file for paraview visualization.
	 *
	 * \param[in]		filename		filename used in time series
	 * \param[in]		u				grid function
	 * \param[in]		si				subset index
	 */
		template <typename TFunction>
		void write_time_pvd_subset(const char* filename, TFunction& u, int si);
protected:
	/**
	 * This function counts the number of vertices, elements and connections for
	 * a given subset (or the whole grid if si < 0).
	 *
	 * \param[in]		u			discrete function
	 * \param[in]		si			subset
	 * \param[in]		dim			dimension of subset
	 * \param[out]		numVert		number of vertices
	 * \param[out]		numElem		number of elements
	 * \param[out]		numConn		number of connections
	 */
		template <typename T>
		void
		count_piece_sizes(Grid& grid, const T& iterContainer, int si, int dim,
		                  int& numVert, int& numElem, int& numConn);

	/**
	 * This function counts the number of vertices, elements and connections that
	 * are in a subset of the grid (or whole grid, if and only if si < 0).
	 *
	 * NOTE: that the number of vertices that are needed to describe the subset
	 * may be more than the vertices of the subset. This happens, when elements are
	 * part of the subset, but the vertices of these elements are part of another
	 * subset. This function counts also those vertices, i.e. it counts the closure
	 * of all vertices needed to describe the subset.
	 *
	 * \param[in]		u			discrete function
	 * \param[in]		si			subset index (si < 0 indicates whole grid)
	 * \param[out]		numVert		number of vertices needed to compose subset
	 * \param[out]		numElem		number of elements in the subset
	 * \param[out]		numConn		number of connections
	 */
		template <typename TElem, typename T>
		void
		count_sizes(Grid& grid, const T& iterContainer, int si,
		            int& numVert, int& numElem, int& numConn);

		template <typename T>
		void
		write_points_cells_piece(VTKFileWriter& File,
		                         Grid::VertexAttachmentAccessor<Attachment<int> >& aaVrtIndex,
		                         const Grid::VertexAttachmentAccessor<Attachment<MathVector<TDim> > >& aaPos,
		                         Grid& grid, const T& iterContainer, int si, int dim,
		                         int numVert, int numElem, int numConn, MGSubsetHandler& sh);

		template <typename T>
		void
		write_points_cells_piece(VTKFileWriter& File,
		                         Grid::VertexAttachmentAccessor<Attachment<int> >& aaVrtIndex,
		                         const Grid::VertexAttachmentAccessor<Attachment<MathVector<TDim> > >& aaPos,
		                         Grid& grid, const T& iterContainer, SubsetGroup& ssGrp, int dim,
		                         int numVert, int numElem, int numConn);

		template <typename T>
		void
		write_grid_piece(VTKFileWriter& File,
		                 Grid::VertexAttachmentAccessor<Attachment<int> >& aaVrtIndex,
		                 const Grid::VertexAttachmentAccessor<Attachment<MathVector<TDim> > >& aaPos,
		                 Grid& grid, const T& iterContainer, int si, int dim, MGSubsetHandler& sh);

public:
		// maybe somebody wants to do this from outside
		static void write_empty_grid_piece(VTKFileWriter& File,
				bool binary = true);

		void
		set_user_defined_comment(const char* comment);

protected:

	/**
	 * This function writes the vertices of a piece of the grid to a vtk file. If
	 * si >= 0 only all vertices needed to form the subset are written; if si < 0
	 * the whole grid is plotted.
	 *
	 * \param[in,out]	File		file to write the points
	 * \param[in]		u			discrete function
	 * \param[in]		si			subset
	 * \param[in]		dim			dimension of subset
	 * \param[in]		numVert		number of vertices
	 */

		void
		write_comment(VTKFileWriter& File);

		void
		write_comment_printf(FILE* File);

		template <typename T>
		void
		write_points(VTKFileWriter& File,
		             Grid::VertexAttachmentAccessor<Attachment<int> >& aaVrtIndex,
		             const Grid::VertexAttachmentAccessor<Attachment<MathVector<TDim> > >& aaPos,
		             Grid& grid, const T& iterContainer, int si, int dim,
		             int numVert);

	/**
	 * This function writes the vertices of a piece (the specified subsets) of
	 * the grid to a vtk file.
	 *
	 * \param[in,out]	File		file to write the points
	 * \param[in]		u			discrete function
	 * \param[in]		ssGrp		specified subsets
	 * \param[in]		dim			dimension of subsets
	 * \param[in]		numVert		number of vertices
	 */
		template <typename T>
		void
		write_points(VTKFileWriter& File,
		             Grid::VertexAttachmentAccessor<Attachment<int> >& aaVrtIndex,
		             const Grid::VertexAttachmentAccessor<Attachment<MathVector<TDim> > >& aaPos,
		             Grid& grid, const T& iterContainer, SubsetGroup& ssGrp, int dim,
		             int numVert);

	/**
	 * This method loops all elements of a subset and writes the vertex
	 * coordinates to a file. All vertices of the element are written reguardless
	 * if the vertex itself is contained in the subset or another. Written
	 * vertices are marked and not written again.
	 *
	 * \param[in]	File		file stream
	 * \param[in]	u			grid function
	 * \param[in]	si			subset index
	 * \param[in]	n			counter for vertices
	 */
		template <typename TElem, typename T>
		void
		write_points_elementwise(VTKFileWriter& File,
		                         Grid::VertexAttachmentAccessor<Attachment<int> >& aaVrtIndex,
		                         const Grid::VertexAttachmentAccessor<Attachment<MathVector<TDim> > >& aaPos,
		                         Grid& grid, const T& iterContainer, int si, int& n);

	/**
	 * This function writes the elements that are part of a given subset. If si < 0
	 * is passed, the whole grid is written.
	 *
	 * \param[in,out]	File		file to write the points
	 * \param[in]		u			discrete function
	 * \param[in]		si			subset
	 * \param[in]		dim			dimension of subset
	 * \param[out]		numElem		number of elements
	 * \param[out]		numConn		number of connections
	 */
		template <typename T>
		void
		write_cells(VTKFileWriter& File,
		            Grid::VertexAttachmentAccessor<Attachment<int> >& aaVrtIndex,
		            Grid& grid, const T& iterContainer, int si,
		            int dim, int numElem, int numConn, MGSubsetHandler& sh);

	/**
	 * This function writes the elements that are part of the specified subsets.
	 *
	 * \param[in,out]	File		file to write the points
	 * \param[in]		u			discrete function
	 * \param[in]		ssGrp		specified subsets
	 * \param[in]		dim			dimension of subsets
	 * \param[out]		numElem		number of elements
	 * \param[out]		numConn		number of connections
	 */
		template <typename T>
		void
		write_cells(VTKFileWriter& File,
		            Grid::VertexAttachmentAccessor<Attachment<int> >& aaVrtIndex,
		            Grid& grid, const T& iterContainer, SubsetGroup& ssGrp,
		            int dim, int numElem, int numConn);


	/**
	 * This method writes the 'connectivity' for each element of a subset. The
	 * connectivity are the indices of all vertex the element is formed of.
	 *
	 * \param[in]	File		file stream
	 * \param[in]	u			grid function
	 * \param[in]	si			subset index
	 */
		template <typename TElem, typename T>
		void
		write_cell_connectivity(VTKFileWriter& File,
		                        Grid::VertexAttachmentAccessor<Attachment<int> >& aaVrtIndex,
		                        Grid& grid, const T& iterContainer, int si);

		template <typename T>
		void
		write_cell_connectivity(VTKFileWriter& File,
		                        Grid::VertexAttachmentAccessor<Attachment<int> >& aaVrtIndex,
		                        Grid& grid, const T& iterContainer, int si, int dim,
		                        int numConn);
		                        
		template <typename T>
		void
		write_cell_connectivity(VTKFileWriter& File,
		                        Grid::VertexAttachmentAccessor<Attachment<int> >& aaVrtIndex,
		                        Grid& grid, const T& iterContainer, SubsetGroup& ssGrp, int dim,
		                        int numConn);

	/**
	 * This method writes the 'offset' for each element of a subset.
	 *
	 * \param[in]	File		file stream
	 * \param[in]	u			grid function
	 * \param[in]	si			subset index
	 * \param[in]	n			counter for vertices
	 */
		template <typename TElem, typename T>
		void
		write_cell_offsets(VTKFileWriter& File, const T& iterContainer, int si, int& n);

		template <typename T>
		void
		write_cell_offsets(VTKFileWriter& File, const T& iterContainer, int si, int dim,
		                   int numElem);

		template <typename T>
		void
		write_cell_offsets(VTKFileWriter& File, const T& iterContainer, SubsetGroup& ssGrp, int dim,
		                   int numElem);

	/**
	 * This method writes the 'type' for each element of a subset. The type is
	 * a vtk-format specified number for a given element type.
	 *
	 * \param[in]	File		file stream
	 * \param[in]	u			grid function
	 * \param[in]	si			subset index
	 */
		template <typename TElem, typename T>
		void
		write_cell_types(VTKFileWriter& File, const T& iterContainer, int si);

		template <typename T>
		void
		write_cell_types(VTKFileWriter& File, const T& iterContainer, int si, int dim,
		                 int numElem);

		template <typename T>
		void
		write_cell_types(VTKFileWriter& File, const T& iterContainer, SubsetGroup& ssGrp, int dim,
		                 int numElem);

	/**
	 * This method writes the 'region' for each element of a subset. The region is
	 * the ug4 subset index for a given element.
	 *
	 * \param[in]	File		file stream
	 * \param[in]	u			grid function
	 * \param[in]	si			subset index
	 * \param[in]	sh			multigrid subset handler
	 */
		template <typename TElem, typename T>
		void
		write_cell_subsets(VTKFileWriter& File, const T& iterContainer, int si, MGSubsetHandler& sh);

		template <typename T>
		void
		write_cell_subsets(VTKFileWriter& File, const T& iterContainer, int si, int dim,
		                 int numElem, MGSubsetHandler& sh);

		template <typename T>
		void
		write_cell_subsets(VTKFileWriter& File, const T& iterContainer, SubsetGroup& ssGrp, int dim,
		                 int numElem, MGSubsetHandler& sh);

	/**
	 * This method writes the subset names.
	 *
	 * \param[in]	File		file stream
	 * \param[in]	sh			multigrid subset handler
	 */
		void
		write_cell_subset_names(VTKFileWriter& File, MGSubsetHandler& sh);

	/**
	 * This method writes the proc ranks for each cell.
	 *
	 * \param[in]	File		file stream
	 * \param[in]	u			grid function
	 * \param[in]	si			subset index
	 * \param[in]	sh			multigrid subset handler
	 */
		template <typename TElem, typename T>
		void
		write_cell_proc_ranks(VTKFileWriter& File, const T& iterContainer, int si, MGSubsetHandler& sh);

		template <typename T>
		void
		write_cell_proc_ranks(VTKFileWriter& File, const T& iterContainer, int si, int dim,
		                 int numElem, MGSubsetHandler& sh);

		template <typename T>
		void
		write_cell_proc_ranks(VTKFileWriter& File, const T& iterContainer, SubsetGroup& ssGrp, int dim,
		                 int numElem, MGSubsetHandler& sh);

	protected:
	/**
	 * This function writes a piece of the grid to the vtk file. First the
	 * geometric data is written (points, cells, connectivity). Then the data
	 * associated to the points or cells is written.
	 *
	 * \param[in,out]	File		file to write the points
	 * \param[in]		u			discrete function
	 * \param[in]		si			subset
	 * \param[in]		dim			dimension of subset
	 */
		template <typename TFunction>
		void
		write_grid_solution_piece(VTKFileWriter& File,
								  Grid::VertexAttachmentAccessor<Attachment<int> >& aaVrtIndex,
								  Grid& grid,
								  TFunction& u, number time, int si, int dim);

	/**
	 * This function writes a piece of the grid to the vtk file. First the
	 * geometric data is written (points, cells, connectivity). Then the data
	 * associated to the points or cells is written.
	 *
	 * \param[in,out]	File		file to write the points
	 * \param[in]		u			discrete function
	 * \param[in]		ssGrp		subsets
	 * \param[in]		dim			dimension of subset
	 */
		template <typename TFunction>
		void
		write_grid_solution_piece(VTKFileWriter& File,
								  Grid::VertexAttachmentAccessor<Attachment<int> >& aaVrtIndex,
								  Grid& grid,
								  TFunction& u, number time, SubsetGroup& ssGrp, int dim);

	///////////////////////////////////////////////////////////////////////////
	// nodal data

	/**
	 * This method writes the nodal values of a function to file. If the function
	 * group consists of more than on component, than a vector of components is
	 * expected. For each function, access to the (unique) nodal value is
	 * required.
	 *
	 * \param[in]	File		file stream
	 * \param[in]	u			grid function
	 * \param[in]	vFct		components to be written
	 * \param[in]	grid		Grid
	 * \param[in]	si			subset index
	 */
		template <typename TElem, typename TFunction>
		void write_nodal_values_elementwise(VTKFileWriter& File, TFunction& u,
		                                    const std::vector<size_t>& vFct,
		                                    Grid& grid, int si);

	/**
	 * This function writes the values of a function as a \<DataArray\> field to
	 * the vtu file.
	 *
	 * \param[in,out]	File		file to write the points
	 * \param[in]		u			discrete function
	 * \param[in]		vFct		components to be written
	 * \param[in]		name		name to appear for the component
	 * \param[in]		grid		Grid
	 * \param[in]		si			subset
	 * \param[in]		dim			dimension of subset
	 * \param[in]		numVert		number of vertices
	 */
		template <typename TFunction>
		void write_nodal_values(VTKFileWriter& File, TFunction& u,
								const std::vector<size_t>& vFct,
								const std::string& name,
								Grid& grid, int si, int dim, int numVert);

	/**
	 * This function writes the values of a function as a \<DataArray\> field to
	 * the vtu file.
	 *
	 * \param[in,out]	File		file to write the points
	 * \param[in]		u			discrete function
	 * \param[in]		vFct		components to be written
	 * \param[in]		name		name to appear for the component
	 * \param[in]		grid		Grid
	 * \param[in]		ssGrp		subsets
	 * \param[in]		dim			dimension of subsets
	 * \param[in]		numVert		number of vertices
	 */
		template <typename TFunction>
		void write_nodal_values(VTKFileWriter& File, TFunction& u,
								const std::vector<size_t>& vFct,
								const std::string& name,
								Grid& grid, SubsetGroup& ssGrp, int dim, int numVert);

	///	writes the nodal scalar data
	/// \{
		template <typename TElem, typename TFunction, typename TData>
		void write_nodal_data_elementwise(VTKFileWriter& File, TFunction& u,
		                                  number time,
		                                  SmartPtr<UserData<TData, TDim> > spData,
		                                  Grid& grid, int si);
		template <typename TFunction, typename TData>
		void write_nodal_data(VTKFileWriter& File, TFunction& u, number time,
		                      SmartPtr<UserData<TData, TDim> > spData,
		                      const int numCmp,
		                      const std::string& name,
		                      Grid& grid, int si, int dim, int numVert);
		template <typename TFunction, typename TData>
		void write_nodal_data(VTKFileWriter& File, TFunction& u, number time,
		                      SmartPtr<UserData<TData, TDim> > spData,
		                      const int numCmp,
		                      const std::string& name,
		                      Grid& grid, SubsetGroup& ssGrp, int dim, int numVert);
	/// \}

		template <typename TFunction>
		void write_nodal_values_piece(VTKFileWriter& File, TFunction& u, number time,
		                              Grid& grid, int si, int dim, int numVert, Grid::VertexAttachmentAccessor<Attachment<int> >&aaVrtIndex);

		template <typename TFunction>
		void write_nodal_values_piece(VTKFileWriter& File, TFunction& u, number time,
		                              Grid& grid, SubsetGroup& ssGrp, int dim, int numVert, Grid::VertexAttachmentAccessor<Attachment<int> >&aaVrtIndex);

	///////////////////////////////////////////////////////////////////////////
	// cell data

	/**
	 * This method writes the cell values of a function to file. If the function
	 * group consists of more than on component, than a vector of components is
	 * expected.
	 *
	 * \param[in]	File		file stream
	 * \param[in]	u			grid function
	 * \param[in]	vFct		components to be written
	 * \param[in]	grid		Grid
	 * \param[in]	si			subset index
	 */
		template <typename TElem, typename TFunction>
		void write_cell_values_elementwise(VTKFileWriter& File, TFunction& u,
										   const std::vector<size_t>& vFct,
										   Grid& grid, int si);

	/**
	 * This function writes the values of a function as a \<DataArray\> field to
	 * the vtu file.
	 *
	 * \param[in,out]	File		file to write the points
	 * \param[in]		u			discrete function
	 * \param[in]		vFct		components to be written
	 * \param[in]		name		name to appear for the component
	 * \param[in]		grid		Grid
	 * \param[in]		si			subset
	 * \param[in]		dim			dimension of subset
	 * \param[in]		numVert		number of vertices
	 */
		template <typename TFunction>
		void write_cell_values(VTKFileWriter& File, TFunction& u,
								const std::vector<size_t>& vFct,
								const std::string& name,
								Grid& grid, int si, int dim, int numElem);

	/**
	 * This function writes the values of a function as a \<DataArray\> field to
	 * the vtu file.
	 *
	 * \param[in,out]	File		file to write the points
	 * \param[in]		u			discrete function
	 * \param[in]		vFct		components to be written
	 * \param[in]		name		name to appear for the component
	 * \param[in]		grid		Grid
	 * \param[in]		ssGrp		subset
	 * \param[in]		dim			dimension of subsets
	 * \param[in]		numVert		number of vertices
	 */
		template <typename TFunction>
		void write_cell_values(VTKFileWriter& File, TFunction& u,
								const std::vector<size_t>& vFct,
								const std::string& name,
								Grid& grid, SubsetGroup& ssGrp, int dim, int numElem);

	///	writes the nodal cell data
	/// \{
		template <typename TElem, typename TFunction, typename TData>
		void write_cell_data_elementwise(VTKFileWriter& File, TFunction& u, number time,
										  SmartPtr<UserData<TData, TDim> > spData,
										  Grid& grid, int si);
		template <typename TFunction, typename TData>
		void write_cell_data(VTKFileWriter& File, TFunction& u, number time,
							  SmartPtr<UserData<TData, TDim> > spData,
							  const int numCmp,
							  const std::string& name,
							  Grid& grid, int si, int dim, int numElem);
		template <typename TFunction, typename TData>
		void write_cell_data(VTKFileWriter& File, TFunction& u, number time,
							  SmartPtr<UserData<TData, TDim> > spData,
							  const int numCmp,
							  const std::string& name,
							  Grid& grid, SubsetGroup& ssGrp, int dim, int numElem);
	/// \}

		template <typename TFunction>
		void write_cell_values_piece(VTKFileWriter& File, TFunction& u, number time,
									  Grid& grid, int si, int dim, int numElem);

		template <typename TFunction>
		void write_cell_values_piece(VTKFileWriter& File, TFunction& u, number time,
									  Grid& grid, SubsetGroup& ssGrp, int dim, int numElem);

	///////////////////////////////////////////////////////////////////////////
	// file names

	///	writes a grouping *.pvtu file
		template <typename TFunction>
		void write_pvtu(TFunction& u, const std::string&  filename,
		                int si, int step, number time);

	public:
	///	writes a grouping *.pvd file, grouping all data from different subsets
		static void write_subset_pvd(int numSubset, const std::string&  filename,
		                             int step = -1, number time = 0.0);

	///	creates the needed vtu file name
		static void vtu_filename(std::string& nameOut, std::string nameIn,
		                         int rank, int si, int maxSi, int step);

	///	create the needed pvtu file name
		static void pvtu_filename(std::string& nameOut, std::string nameIn,
		                          int si, int maxSi, int step);

	///	creates the needed pvd file name
		static void pvd_filename(std::string& nameOut, std::string nameIn);

	///	creates the needed pvd file name to group the time steps
		static void pvd_time_filename(std::string& nameOut, std::string nameIn,
		                              int step);

	public:
	///	default constructor
		VTKOutput()	: m_bSelectAll(true), m_bBinary(true), m_bWriteGrid(true), m_bWriteSubsetIndices(false), m_bWriteProcRanks(false) {} //TODO: maybe true?

	/// should values be printed in binary (base64 encoded way ) or plain ascii
		void set_binary(bool b);

		void set_write_grid(bool b);

		void set_write_subset_indices(bool b);

		void set_write_proc_ranks(bool b);

		void set_write_grid_points_ordering(bool b);

		void set_write_grid_function_ordering(bool b);

		void set_write_sorted_grid_function_ordering(bool b);

		void set_lua_ordering(SmartPtr<LuaOrdering> lua_ordering){
			m_spLuaOrdering = lua_ordering;
		}

	protected:
	///	returns true if name for vtk-component is already used
		bool vtk_name_used(const char* name) const;

	///	writes data to stream
	/**
	 * The purpose of the function is to convert a double data to binary float
	 * (the Float32-format as used in vtk) or ASCII representation. Note that
	 * the ASCII output should contain only numbers in the range of the normalized
	 * Float32 representation, i.e. with the minimum absolute value
	 * 1.1755e-38 = 2^-126. As 'float'-variables can potentially store denormalized
	 * numbers (and this is really the case in the practice) with the minimum
	 * absolute value 1.4013e-45 = 2^-149 that are not correctly read by some
	 * visualization tools (in some operating systems), these function replace
	 * those values (lying near to 0) with 0.
	 */
	/// \{
		inline void write_item_to_file(VTKFileWriter& File, float data);
		inline void write_item_to_file(VTKFileWriter& File, double data);
		inline void write_item_to_file(VTKFileWriter& File, const ug::MathVector<1>& data);
		inline void write_item_to_file(VTKFileWriter& File, const ug::MathVector<2>& data);
		inline void write_item_to_file(VTKFileWriter& File, const ug::MathVector<3>& data);
	/// \}


	protected:
	///	scheduled components to be printed
		bool m_bSelectAll;
	/// print values in binary (base64 encoded way) or plain ascii
		bool m_bBinary;
		std::map<std::string, std::vector<std::string> > m_vSymbFct;
		std::map<std::string, std::vector<std::string> > m_vSymbFctNodal;
		std::map<std::string, std::vector<std::string> > m_vSymbFctElem;
		typedef typename std::map<std::string,
				std::vector<std::string> >::iterator ComponentsIterator;

	///	scheduled scalar data to be printed
		std::map<std::string, SmartPtr<UserData<number, TDim> > > m_vScalarNodalData;
		std::map<std::string, SmartPtr<UserData<number, TDim> > > m_vScalarElemData;
		typedef typename std::map<std::string,
				SmartPtr<UserData<number, TDim> > >::iterator ScalarDataIterator;

	///	scheduled vector data to be printed
		std::map<std::string, SmartPtr<UserData<MathVector<TDim>, TDim> > > m_vVectorNodalData;
		std::map<std::string, SmartPtr<UserData<MathVector<TDim>, TDim> > > m_vVectorElemData;
		typedef typename std::map<std::string,
				SmartPtr<UserData<MathVector<TDim>, TDim> > >::iterator VectorDataIterator;

	///	map storing the time points
		std::map<std::string, std::vector<number> > m_mTimestep;

		bool m_bWriteGrid;
		std::string m_sComment;

		bool m_bWriteSubsetIndices;
		bool m_bWriteProcRanks;
		bool m_bWriteGridPointsOrdering;
		bool m_bWriteGridFunctionOrdering;
		bool m_bWriteSortedGridFunctionOrdering;

		SmartPtr<LuaOrdering> m_spLuaOrdering;
};

} // namespace ug

#include "vtkoutput_impl.h"

#endif /* __H__UG__LIB_DISC__IO__VTKOUTPUT__ */
