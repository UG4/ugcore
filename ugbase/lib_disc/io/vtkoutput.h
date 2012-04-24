/*
 * vtkoutput.h
 *
 *  Created on: 06.07.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__IO__VTKOUTPUT__
#define __H__UG__LIB_DISC__IO__VTKOUTPUT__

// extern libraries
#include <vector>

// other ug modules
#include "lib_grid/lg_base.h"
#include "common/util/string_util.h"
#include "lib_disc/common/function_group.h"

namespace ug{

/// output writer to the VTK file format
/**
 * This class can be used to write grid and data associated with the grid to
 * then vtk file format.
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
 * \tparam		TFunction 	discrete grid function used
 */
template <typename TFunction>
class VTKOutput{
	public:
	///	type of grid function
		typedef TFunction function_type;

	public:
	///	clears the selected output
		void clear_selection() {m_vSymbFct.clear();}

	///	schedules that all components of the passed discrete functions will be written to file
		void select_all(bool bAll) {m_bSelectAll = bAll;}

	///	selects a nodal scalar value to be written
	/**
	 * This function schedules the component passed by symbolic name to be
	 * written to the vtk file under a specified name. Note, that for the
	 * ansatz space of the component an evaluation of the data at the nodes
	 * must be available (continuous).
	 *
	 * example: fctName = "p"; name = "pressure"
	 *
	 * \param[in]	fctName		symbolic name of component
	 * \param[in]	name		name that will appear in the vtk file for the data
	 */
		void select_nodal_scalar(const char* fctName, const char* name);

	///	selects a nodal vector value to be written
	/**
	 * This function schedules a vector to be written to the vtk file. The
	 * symbolic names of the vector components have to be specified as components
	 * of the passed grid function. Evaluation at nodes must be available for
	 * all components.
	 *
	 * example: fctNames = "u,v,w"; name = "velocity"
	 *
	 * \param[in]	fctName		symbolic name of components
	 * \param[in]	name		name that will appear in the vtk file for the data
	 */
		void select_nodal_vector(const char* fctNames, const char* name);

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
		void print(const char*  filename, function_type& u,
		           int step, number time,
				   bool makeConsistent);

	/**	Calls print with makeConsistent enabled.*/
		void print(const char*  filename, function_type& u,
		           int step, number time)
		{
			return print(filename, u, step, time, true);
		}

	/**
	 * This function simply calles 'print' using step = -1 to indicate the
	 * stationary case. It is intended to write time independent data.
	 *
	 * \param[in]	filename		filename for produced files
	 * \param[in]	u				grid function
	 * \param[in]	makeConsistent	flag if data shall be made consistent before printing
	 */
		void print(const char*  filename, function_type& u,
				   bool makeConsistent)
		{
			print(filename, u, -1, 0.0, makeConsistent);
		}

	/**	Calls print with makeConsistent enabled.*/
		void print(const char*  filename, function_type& u)
		{
			return print(filename, u, true);
		}

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
		void print_subset(const char* filename, function_type& u,
		                  int si, int step = -1, number time = 0.0,
						  bool makeConsistent = true);

	/**	Calls print_subset with makeConsistent enabled.*/
		void print_subset(const char* filename, function_type& u,
		                  int si, int step = -1, number time = 0.0)
		{
			return print_subset(filename, u, si, step, time, true);
		}

	/**
	 * When a time series has been computed, this function can be used to procduce
	 * a grouping *.pvd file for paraview visualization.
	 *
	 * \param[in]		filename		filename used in time series
	 * \param[in]		u				grid function
	 */
		void write_time_pvd(const char* filename, function_type& u);

	///	sets the verbosity, if true some output status info is written to log
		void set_verbose(bool bVerb) {m_bVerb = bVerb;}

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
		void count_piece_sizes(function_type& u, int si, int dim,
		                       int& numVert, int& numElem, int& numConn);

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
		void write_piece(FILE* File, function_type& u, int si, int dim);

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
		void write_points(FILE* File, function_type& u, int si,
		                  int dim, int numVert);

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
		void write_cells(FILE* File, function_type& u, int si,
		                 int dim, int numElem, int numConn);


	/**
	 * This function writes the values of a function as a <DataArray> field to
	 * the vtu file.
	 *
	 * \param[in,out]	File		file to write the points
	 * \param[in]		u			discrete function
	 * \param[in]		vFct		components to be written
	 * \param[in]		name		name to appear for the component
	 * \param[in]		si			subset
	 * \param[in]		dim			dimension of subset
	 * \param[in]		numVert		number of vertices
	 */
		void write_nodal_values(FILE* File, function_type& u,
		                        const FunctionGroup& vFct,
		                        const std::string& name,
		                        int si, int dim, int numVert);

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
		template <typename TElem>
		void count_sizes(function_type& u, int si,
		                     int& numVert, int& numElem, int& numConn);

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
		template <typename TElem>
		void write_points_elementwise(FILE* File, function_type& u, int si, int& n);

	/**
	 * This method writes the 'connectivity' for each element of a subset. The
	 * connectivity are the indices of all vertex the element is formed of.
	 *
	 * \param[in]	File		file stream
	 * \param[in]	u			grid function
	 * \param[in]	si			subset index
	 */
		template <typename TElem>
		void write_cell_connectivity(FILE* File, function_type& u, int si);

	/**
	 * This method writes the 'offset' for each element of a subset.
	 *
	 * \param[in]	File		file stream
	 * \param[in]	u			grid function
	 * \param[in]	si			subset index
	 * \param[in]	n			counter for vertices
	 */
		template <typename TElem>
		void write_cell_offsets(FILE* File, function_type& u, int si, int& n);

	/**
	 * This method writes the 'type' for each element of a subset. The type is
	 * a vtk-format specified number for a given element type.
	 *
	 * \param[in]	File		file stream
	 * \param[in]	u			grid function
	 * \param[in]	si			subset index
	 */
		template <typename TElem>
		void write_cell_types(FILE* File, function_type& u, int si);

	/**
	 * This method writes the nodal values of a function to file. If the function
	 * group consists of more than on component, than a vector of components is
	 * expected. For each function, access to the (unique) nodal value is
	 * requiered.
	 *
	 * \param[in]	File		file stream
	 * \param[in]	u			grid function
	 * \param[in]	vFct		components to be written
	 * \param[in]	si			subset index
	 */
		template <typename TElem>
		void write_nodal_values_elementwise(FILE* File, function_type& u,
		                                   const FunctionGroup& vFct, int si);

	private:
	///	writes a grouping *.pvtu file
		void write_pvtu(function_type& u, const std::string&  filename,
		                int si, int step, number time);

	///	writes a grouping *.pvd file, grouping all data from different subsets
		void write_subset_pvd(function_type& u, const std::string&  filename,
		                      int step = -1, number time = 0.0);

	///	creates the needed vtu file name
		void vtu_filename(std::string& nameOut, std::string nameIn,
		                  int rank, int si, int maxSi, int step);

	///	create the needed pvtu file name
		void pvtu_filename(std::string& nameOut, std::string nameIn,
		                   int si, int maxSi, int step);

	///	creates the needed pvd file name
		void pvd_filename(std::string& nameOut, std::string nameIn);

	///	creates the needed pvd file name to group the time steps
		void pvd_time_filename(std::string& nameOut, std::string nameIn,
		                       int step);

	///	checks if chosen filename is admissible
		bool is_valid_filename(std::string& nameIn);

	public:
	///	default constructor
		VTKOutput()	: m_bVerb(false), m_bSelectAll(true), m_pGrid(NULL) {}

	protected:
	///	verbosity flag
		bool m_bVerb;

	///	scheduled components to be printed
		bool m_bSelectAll;
		std::vector<std::pair<std::string, std::string> > m_vSymbFct;

	///	Grid, working on
		Grid* m_pGrid;

	///	attachment for dofs
		typedef ug::Attachment<int> ADOFIndex;
		ADOFIndex m_aDOFIndex;
		Grid::VertexAttachmentAccessor<ADOFIndex> m_aaDOFIndexVRT;

	protected:
	///	map storing the time points
		std::map<std::string, std::vector<number> > m_mTimestep;
};



} // namespace ug

#include "vtkoutput_impl.h"

#endif /* __H__UG__LIB_DISC__IO__VTKOUTPUT__ */
