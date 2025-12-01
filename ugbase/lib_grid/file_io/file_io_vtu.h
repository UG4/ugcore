/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#ifndef __H__LIB_GRID__FILE_IO_VTU__
#define __H__LIB_GRID__FILE_IO_VTU__

#include <iostream>
#include <vector>
#include <utility>
#include "common/parser/rapidxml/rapidxml.hpp"
#include "lib_grid/grid/grid.h"
#include "lib_grid/multi_grid.h"
#include "lib_grid/tools/subset_handler_interface.h"
#include "lib_grid/tools/selector_interface.h"
#include "lib_grid/common_attachments.h"
#include "lib_grid/grid_objects/grid_objects.h"
#include "lib_grid/callbacks/basic_callbacks.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
///	Reads a grid to an vtu (vtk) file. internally uses GridReaderVTK.
/**	The position attachment can be specified. Since the type of the
 *	position attachment is a template parameter, MathVector attachments
 * 	of any dimension are supported. Especially ug::aPosition, ug::aPostion2
 *	and ug::aPosition1.
 */
template <typename TAPosition>
bool LoadGridFromVTU(Grid& grid, ISubsetHandler& sh,
					const char* filename, APosition& aPos);

///	Reads a grid to a vtu (vtk unstructured mesh) file.
/**	Before reading a grid from file, this method searches for the
 *	attached standard position attachment with the highest dimension.
 *	This will be used as position-attachment in a call to the overloaded
 *	version of LoadGridFromVTU.
 *
 *	If no standard attachment is found, aPosition will be attached and used.
 */
bool LoadGridFromVTU(Grid& grid, ISubsetHandler& sh,
					const char* filename);
					

///	Writes a grid to a vtu (vtk unstructured mesh) file
template <typename TAPosition>
bool SaveGridToVTU(Grid& grid, ISubsetHandler* psh, const char* filename,
				   TAPosition& aPos);


////////////////////////////////////////////////////////////////////////					
////////////////////////////////////////////////////////////////////////
///	Grants write access to vtu files.
/**	Make sure that all elements added via one of the add_* methods
 *	exist until the FileAccessor is destroyed.
 */
class GridWriterVTU
{
	public:
		GridWriterVTU();
	///	Pass a pointer to an ostream to which the data shall be written.
	/**	Make sure, that the stream is open and stays valid while write operations
	 * are performed.*/
		explicit GridWriterVTU(std::ostream* out);

		virtual ~GridWriterVTU() = default;

	///	Pass a pointer to an ostream to which the data shall be written.
	/**	Make sure, that the stream is open and stays valid while write operations
	 * are performed.*/
		void set_stream(std::ostream* out);

		void finish();

	/**	TPositionAttachments value type has to be compatible with MathVector.
	 * Make sure that aPos is attached to the vertices of the grid.
	 * If a SubsetHandler is specified through 'psh' (nullptr is valid too), it will
	 * automatically be passed to 'add_subset_handler' with the name "regions".
	 * If you pass a subset handler to this method, furthermore only those elements
	 * which are assigned to subsets will be written to the file as cells.
	 * If you don't pass a subset-handler, only elements of highest dimension are
	 * written to the file.*/
		template <typename TPositionAttachment>
		bool new_piece(Grid& grid, ISubsetHandler* psh,
					   TPositionAttachment& aPos);

	///	You may add subset-handlers which will be written as regions to the vtk-file.
	/**	Make sure to add subset-handlers before 'end_cell_data' is called.
	 * Note that if you pass a subset-handler to 'new_piece', then it will be
	 * automatically registered with the name "regions".*/
		void add_subset_handler(ISubsetHandler& sh, const std::string& name);

		void begin_point_data();
		//...
		void end_point_data();

		void begin_cell_data();
		void end_cell_data();

	protected:
		using AAVrtIndex = Grid::VertexAttachmentAccessor<AInt>;
		using AAEdgeIndex = Grid::EdgeAttachmentAccessor<AInt>;
		using AAFaceIndex = Grid::FaceAttachmentAccessor<AInt>;
		using AAVolIndex = Grid::VolumeAttachmentAccessor<AInt>;

		inline std::ostream& out_stream();

	/** \param name:			  Set to "" to omit this parameter in the output.
	 * \param numberOfComponents: Set to 0 to omit this parameter in the output.*/
		void write_data_array_header(const char* type, const char* name,
									 int numberOfComponents);

		void write_data_array_footer();

		template <typename TElem, typename TAttachment>
		void write_vector_data(Grid& grid,
							   TAttachment aData,
							   const char* name = "",
							   typename Grid::traits<TElem>::callback consider_elem = ConsiderAll());

		template <typename TElem>
		void collect_cells(std::vector<GridObject*>& cellsOut,
						   Grid& grid,
						   typename Grid::traits<TElem>::callback consider_elem = ConsiderAll());

		void write_cells(const std::vector<GridObject*>& cells,
						 Grid& grid,
						 AAVrtIndex aaInd);

		void end_piece();
		

		enum Mode{
			NONE,
			OPEN,
			CLOSED
		};

		std::ostream*	m_pout;
		Mode			m_pieceMode;
		Mode			m_pointDataMode;
		Mode			m_cellDataMode;

		Grid*		 	m_curGrid;
		ISubsetHandler*	m_curSH;
				
		std::vector<GridObject*>		m_cells;
		std::vector<std::pair<ISubsetHandler*, std::string> >	m_pieceSubsetHandlers;
};



////////////////////////////////////////////////////////////////////////
///	Grants read access to vtu (vtk) files.
/**	Before any data can be retrieved using the get_* methods, a file
 *	has to be successfully loaded using load_file.
 *
 *	\todo: Improve performance by using in-situ stringstreams during element creation.
 */
class GridReaderVTU
{
	public:
		GridReaderVTU() = default;
		virtual ~GridReaderVTU() = default;

	///	parses an xml file
		bool parse_file(const char* filename);

	///	returns the number of grids in the given file
		[[nodiscard]] size_t num_grids() const {return m_entries.size();}

	///	returns the i-th grid.
	/**	TPositionAttachments value type has to be compatible with MathVector.
	 *	Make sure that a file has already been loaded.*/
		template <typename TPositionAttachment>
		bool grid(Grid& gridOut, size_t index, TPositionAttachment& aPos);

	///	returns the name of the i-th grid
		[[nodiscard]] const char* get_grid_name(size_t index) const;

	///	returns the number of subset handlers for the given grid
		[[nodiscard]] size_t num_subset_handlers(size_t refGridIndex) const;

	///	returns the name of the given subset handler
		[[nodiscard]] const char* get_subset_handler_name(size_t refGridIndex,
											size_t subsetHandlerIndex) const;

	///	fills the given subset-handler
		bool subset_handler(ISubsetHandler& shOut,
							size_t refGridIndex,
							size_t subsetHandlerIndex);

		static std::string getRegionOfInterestIdentifyer()
		{ return m_regionOfInterest; }

		static void setRegionOfInterestIdentifier( std::string const & regOfInt )
		{ m_regionOfInterest = regOfInt; }

	protected:
		struct SubsetHandlerEntry
		{
			explicit SubsetHandlerEntry(rapidxml::xml_node<>* n) : node(n), sh(nullptr) {}

			rapidxml::xml_node<>* 	node;
			ISubsetHandler*			sh;
		};

		struct GridEntry
		{
			explicit GridEntry(rapidxml::xml_node<>* n) : node(n), grid(nullptr), mg(nullptr)	{}

			rapidxml::xml_node<>* 	node;
			Grid* 					grid;
			MultiGrid* 				mg;
			std::vector<SubsetHandlerEntry>	subsetHandlerEntries;
			std::vector<Vertex*> 			vertices;
			std::vector<GridObject*>		cells;
		};

	protected:
	///	initializes internal arrays
	/**	searches for all grid-nodes and stores, resizes m_entries and stores
	 *	the node for each entry.
	 *
	 *	If you create your own version of this method, don't forget to call the
	 *	base-class implementation!*/
		virtual bool new_document_parsed();

	///	creates vertices from a vertex-node.
	/**	if aaPos has more coordinates per vertex than the vrtNode,
	 *	0's will be appended. If it has less, unused coordinates will
	 *	be ignored.*/
		template <typename TAAPos>
		bool create_vertices(std::vector<Vertex*>& vrtsOut, Grid& grid,
							rapidxml::xml_node<>* vrtNode, TAAPos aaPos);

		bool create_cells(std::vector<GridObject*>& cellsOut,
						  Grid& grid,
						  rapidxml::xml_node<>* node,
						  std::vector<Vertex*> vertices,
						  size_t pieceVrtOffset);

		template <typename T>
		void read_scalar_data(std::vector<T>& dataOut,
							  rapidxml::xml_node<>* dataNode,
							  bool clearData = true);

		static void trafoDblVec2Int( std::vector<double> const & dblVec, std::vector<int> & intVec );

		template <typename T>
		void check_indices(std::vector<T>& inds, size_t first, size_t num, size_t max);

		static rapidxml::xml_node<>* find_child_node_by_argument_value(
													rapidxml::xml_node<>* parent,
													const char* nodeName,
													const char* argName,
													const char* argValue);
	protected:
	///	stores the file-name
		std::string					m_filename;

	///	the xml_document which stores the data
		rapidxml::xml_document<>	m_doc;

	///	holds grids which already have been created
		std::vector<GridEntry>		m_entries;



		static std::string m_regionOfInterest; // ProMesh standard = "regions", in Braunschweig case	often "Material Id", but not always

};

}//	end of namespace

////////////////////////////////
//	include implementation
#include "file_io_vtu_impl.h"

#endif
