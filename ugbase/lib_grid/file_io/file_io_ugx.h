/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__LIB_GRID__FILE_IO_UGX__
#define __H__LIB_GRID__FILE_IO_UGX__


#include <errno.h>
#include <cstdlib>
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
#include "lib_grid/refinement/projectors/refinement_projector.h"
#include "common/math/misc/shapes.h"	// AABox

namespace ug
{

class ProjectionHandler;

////////////////////////////////////////////////////////////////////////
///	Writes a grid to an ugx file. internally uses GridWriterUGX.
/**	The position attachment can be specified. Since the type of the
 *	position attachment is a template parameter, MathVector attachments
 * 	of any dimension are supported. Especially ug::aPosition, ug::aPostion2
 *	and ug::aPosition1.
 */
template <class TAPosition>
bool SaveGridToUGX(Grid& grid, ISubsetHandler& sh,
				   const char* filename, TAPosition& aPos);

///	Writes a grid to a ugx file.
/**	Before writing a grid to file, this method searches for the
 *	attached standard position attachment with the highest dimension.
 *	This will be used as position-attachment in a call to the overloaded
 *	version of SaveGridToUGX.
 */
bool SaveGridToUGX(Grid& grid, ISubsetHandler& sh,
				   const char* filename);

////////////////////////////////////////////////////////////////////////
///	Reads a grid to an ugx file. internally uses GridReaderUGX.
/**	The position attachment can be specified. Since the type of the
 *	position attachment is a template parameter, MathVector attachments
 * 	of any dimension are supported. Especially ug::aPosition, ug::aPostion2
 *	and ug::aPosition1.
 */
template <class TAPosition>
bool LoadGridFromUGX(Grid& grid, ISubsetHandler& sh,
					const char* filename, APosition& aPos);

///	Reads a grid to an ugx file.
/**	Before reading a grid from file, this method searches for the
 *	attached standard position attachment with the highest dimension.
 *	This will be used as position-attachment in a call to the overloaded
 *	version of LoadGridFromUGX.
 *
 *	If no standard attachment is found, aPosition will be attached and used.
 */
bool LoadGridFromUGX(Grid& grid, ISubsetHandler& sh,
					const char* filename);
					
					

////////////////////////////////////////////////////////////////////////					
////////////////////////////////////////////////////////////////////////
///	Grants write access to ugx files.
/**	Make sure that all elements added via one of the add_* methods
 *	exist until the FileAccessor is destroyed.
 */
class GridWriterUGX
{
	public:
		GridWriterUGX();
		virtual ~GridWriterUGX();

	/**	TPositionAttachments value type has to be compatible with MathVector.
	 *	Make sure that aPos is attached to the vertices of the grid.*/
		template <class TPositionAttachment>
		bool add_grid(Grid& grid, const char* name,
					  TPositionAttachment& aPos);

//		template <class TPositionAttachment>
//		void add_grid(MultiGrid& mg, const char* name,
//					  TPositionAttachment& aPos);

		void add_subset_handler(ISubsetHandler& sh, const char* name,
								size_t refGridIndex);

		void add_selector(ISelector& sel, const char* name,
						  size_t refGridIndex);

		void add_projection_handler(ProjectionHandler& ph, const char* name,
									size_t refGridIndex);

		template <class TElem, class TAttachment>
		void add_attachment(TAttachment attachment,
							const char* name,
							size_t refGridIndex);

		virtual bool write_to_stream(std::ostream& out);

		bool write_to_file(const char* filename);

	protected:
		typedef Grid::VertexAttachmentAccessor<AInt> AAVrtIndex;
		typedef Grid::EdgeAttachmentAccessor<AInt> AAEdgeIndex;
		typedef Grid::FaceAttachmentAccessor<AInt> AAFaceIndex;
		typedef Grid::VolumeAttachmentAccessor<AInt> AAVolIndex;

	protected:
		void init_grid_attachments(Grid& grid);
		
	//	VERTICES
		template <class TAAPos>
		rapidxml::xml_node<>*
		create_vertex_node(RegularVertexIterator vrtsBegin,
						  RegularVertexIterator vrtsEnd,
						  TAAPos& aaPos);

		template <class TAAPos>
		rapidxml::xml_node<>*
		create_constrained_vertex_node(ConstrainedVertexIterator vrtsBegin,
										ConstrainedVertexIterator vrtsEnd,
										TAAPos& aaPos,
										AAEdgeIndex aaIndEDGE,
										AAFaceIndex aaIndFACE);
						
														
	///	adds grid elements (edges, faces, volumes) to the given node.
		void add_elements_to_node(rapidxml::xml_node<>* node,
								  Grid& grid);

	//	EDGES
		rapidxml::xml_node<>*
		create_edge_node(RegularEdgeIterator edgesBegin,
						 RegularEdgeIterator edgesEnd,
						 AAVrtIndex aaIndVRT);

		rapidxml::xml_node<>*
		create_constraining_edge_node(ConstrainingEdgeIterator edgesBegin,
									  ConstrainingEdgeIterator edgesEnd,
									  AAVrtIndex aaIndVRT);
							  
		rapidxml::xml_node<>*
		create_constrained_edge_node(ConstrainedEdgeIterator edgesBegin,
									 ConstrainedEdgeIterator edgesEnd,
									 AAVrtIndex aaIndVRT,
									 AAEdgeIndex aaIndEDGE,
									 AAFaceIndex aaIndFACE);
									 
	//	FACES
		rapidxml::xml_node<>*
		create_triangle_node(TriangleIterator trisBegin,
							 TriangleIterator trisEnd,
							 AAVrtIndex aaIndVRT);

		rapidxml::xml_node<>*
		create_constraining_triangle_node(ConstrainingTriangleIterator trisBegin,
										  ConstrainingTriangleIterator trisEnd,
										  AAVrtIndex aaIndVRT);

		rapidxml::xml_node<>*
		create_constrained_triangle_node(ConstrainedTriangleIterator trisBegin,
							 			 ConstrainedTriangleIterator trisEnd,
										 AAVrtIndex aaIndVRT,
									 	 AAFaceIndex aaIndFACE);
							 
		rapidxml::xml_node<>*
		create_quadrilateral_node(QuadrilateralIterator quadsBegin,
								  QuadrilateralIterator quadsEnd,
								  AAVrtIndex aaIndVRT);

		rapidxml::xml_node<>*
		create_constraining_quadrilateral_node(ConstrainingQuadrilateralIterator quadsBegin,
											   ConstrainingQuadrilateralIterator quadsEnd,
								  			   AAVrtIndex aaIndVRT);

		rapidxml::xml_node<>*
		create_constrained_quadrilateral_node(ConstrainedQuadrilateralIterator quadsBegin,
											  ConstrainedQuadrilateralIterator quadsEnd,
											  AAVrtIndex aaIndVRT,
											  AAFaceIndex aaIndFACE);
	//	VOLUMES
		rapidxml::xml_node<>*
		create_tetrahedron_node(TetrahedronIterator tetsBegin,
								  TetrahedronIterator tetsEnd,
								  AAVrtIndex aaIndVRT);

		rapidxml::xml_node<>*
		create_hexahedron_node(HexahedronIterator hexasBegin,
								  HexahedronIterator hexasEnd,
								  AAVrtIndex aaIndVRT);

		rapidxml::xml_node<>*
		create_prism_node(PrismIterator prismsBegin,
							PrismIterator prismsEnd,
							AAVrtIndex aaIndVRT);

		rapidxml::xml_node<>*
		create_pyramid_node(PyramidIterator pyrasBegin,
							PyramidIterator pyrasEnd,
							AAVrtIndex aaIndVRT);

		rapidxml::xml_node<>*
		create_octahedron_node(OctahedronIterator octsBegin,
								OctahedronIterator octsEnd,
								AAVrtIndex aaIndVRT);

		void add_subset_attributes(rapidxml::xml_node<>* targetNode,
								   ISubsetHandler& sh, size_t subsetIndex);

		template <class TGeomObj>
		rapidxml::xml_node<>*
		create_subset_element_node(const char* name,
								   const ISubsetHandler& sh,
								   size_t si);

		template <class TGeomObj>
		rapidxml::xml_node<>*
		create_selector_element_node(const char* name, const ISelector& sel);


		rapidxml::xml_node<>*
		create_projector_node(RefinementProjector& proj, const char* nodeName);

		template <class TElem>
		void process_global_attachments(Grid& grid, rapidxml::xml_node<>* gridNode);

		template <class TElem>
		const char* attachment_node_name();

	protected:
	///	entries are stored for each grid.
	/**	an entry holds a pointer to a grid together with its xml_node.*/
		struct Entry{
			Entry()	{}
			Entry(Grid* g, rapidxml::xml_node<>* n) :
				grid(g), node(n)	{}

			Grid* grid;
			std::vector<const ISubsetHandler*> subsetHandlers;
			rapidxml::xml_node<>* node;
		};

	///	the xml_document which stores the data
		rapidxml::xml_document<> m_doc;

	///	List of accessible grids
		std::vector<Entry> m_vEntries;

	///	attached to vertices of each grid during add_grid.
		AInt	m_aInt;
};


////////////////////////////////////////////////////////////////////////
///	Grants read access to ugx files.
/**	Before any data can be retrieved using the get_* methods, a file
 *	has to be successfully loaded using load_file.
 *
 *	\todo: Improve performance by using in-situ stringstreams during element creation.
 */
class GridReaderUGX
{
	public:
		GridReaderUGX();
		virtual ~GridReaderUGX();

	///	parses an xml file
		bool parse_file(const char* filename);

	///	returns the number of grids
		inline size_t num_grids() const	{return m_entries.size();}

	///	returns the i-th grid.
	/**	TPositionAttachments value type has to be compatible with MathVector.
	 *	Make sure that a file has already been loaded.*/
		template <class TPositionAttachment>
		bool grid(Grid& gridOut, size_t index,
					  TPositionAttachment& aPos);

	///	returns the name of the i-th grid
		const char* get_grid_name(size_t index) const;

	///	returns the number of subset handlers for the given grid
		size_t num_subset_handlers(size_t refGridIndex) const;

	///	returns the name of the given subset handler
		const char* get_subset_handler_name(size_t refGridIndex,
											size_t subsetHandlerIndex) const;

	///	fills the given subset-handler
		bool subset_handler(ISubsetHandler& shOut,
								size_t subsetHandlerIndex,
								size_t refGridIndex);

	///	returns the number of selectors for the given grid
		size_t num_selectors(size_t refGridIndex) const;

	///	returns the name of the given selector
		const char* get_selector_name(size_t refGridIndex, size_t selectorIndex) const;

	///	fills the given selector
		bool selector(ISelector& selOut, size_t selectorIndex, size_t refGridIndex);

	///	returns the number of projection-handlers for the given grid
		size_t num_projection_handlers(size_t refGridIndex) const;

	///	returns the name of the given projection-handler
		const char* get_projection_handler_name(size_t refGridIndex, size_t phIndex) const;

	///	returns the subset handler index for a projection handler
		size_t get_projection_handler_subset_handler_index(size_t phIndex, size_t refGridIndex);

	///	fills the given projection-handler
		bool projection_handler(ProjectionHandler& phOut, size_t phIndex, size_t refGridIndex);
		
	protected:
		struct SubsetHandlerEntry
		{
			SubsetHandlerEntry(rapidxml::xml_node<>* n) : node(n), sh(NULL) {}

			rapidxml::xml_node<>* 	node;
			ISubsetHandler*			sh;
		};

		struct SelectorEntry
		{
			SelectorEntry(rapidxml::xml_node<>* n) : node(n), sel(NULL) {}

			rapidxml::xml_node<>* 	node;
			ISelector*				sel;
		};

		struct GridEntry
		{
			GridEntry(rapidxml::xml_node<>* n) : node(n), grid(NULL), mg(NULL)	{}

			rapidxml::xml_node<>* node;
			Grid* 		grid;
			MultiGrid* 	mg;
			std::vector<SubsetHandlerEntry>	subsetHandlerEntries;
			std::vector<SelectorEntry>		selectorEntries;
			std::vector<rapidxml::xml_node<>*>	projectionHandlerEntries;
			std::vector<Vertex*> 		vertices;
			std::vector<Edge*> 			edges;
			std::vector<Face*>				faces;
			std::vector<Volume*>			volumes;
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
		template <class TAAPos>
		bool create_vertices(std::vector<Vertex*>& vrtsOut, Grid& grid,
							rapidxml::xml_node<>* vrtNode, TAAPos aaPos);

		template <class TAAPos>
		bool create_constrained_vertices(std::vector<Vertex*>& vrtsOut,
							std::vector<std::pair<int, int> >& constrainingObjsOut,
							Grid& grid, rapidxml::xml_node<>* vrtNode, TAAPos aaPos);
							
		bool create_edges(std::vector<Edge*>& edgesOut,
						  Grid& grid, rapidxml::xml_node<>* node,
			 			  std::vector<Vertex*>& vrts);

		bool create_constraining_edges(std::vector<Edge*>& edgesOut,
						  Grid& grid, rapidxml::xml_node<>* node,
			 			  std::vector<Vertex*>& vrts);

		bool create_constrained_edges(std::vector<Edge*>& edgesOut,
						  std::vector<std::pair<int, int> >& constrainingObjsOut,
						  Grid& grid, rapidxml::xml_node<>* node,
			 			  std::vector<Vertex*>& vrts);

		bool create_triangles(std::vector<Face*>& facesOut,
							  Grid& grid, rapidxml::xml_node<>* node,
							  std::vector<Vertex*>& vrts);

		bool create_constraining_triangles(std::vector<Face*>& facesOut,
							  Grid& grid, rapidxml::xml_node<>* node,
							  std::vector<Vertex*>& vrts);

		bool create_constrained_triangles(std::vector<Face*>& facesOut,
							  std::vector<std::pair<int, int> >& constrainingObjsOut,
							  Grid& grid, rapidxml::xml_node<>* node,
							  std::vector<Vertex*>& vrts);
							  
		bool create_quadrilaterals(std::vector<Face*>& facesOut,
								   Grid& grid, rapidxml::xml_node<>* node,
								   std::vector<Vertex*>& vrts);

		bool create_constraining_quadrilaterals(std::vector<Face*>& facesOut,
							  Grid& grid, rapidxml::xml_node<>* node,
							  std::vector<Vertex*>& vrts);

		bool create_constrained_quadrilaterals(std::vector<Face*>& facesOut,
							  std::vector<std::pair<int, int> >& constrainingObjsOut,
							  Grid& grid, rapidxml::xml_node<>* node,
							  std::vector<Vertex*>& vrts);
							  
		bool create_tetrahedrons(std::vector<Volume*>& volsOut,
								 Grid& grid, rapidxml::xml_node<>* node,
								 std::vector<Vertex*>& vrts);

		bool create_hexahedrons(std::vector<Volume*>& volsOut,
								Grid& grid, rapidxml::xml_node<>* node,
								std::vector<Vertex*>& vrts);

		bool create_prisms(std::vector<Volume*>& volsOut,
							Grid& grid, rapidxml::xml_node<>* node,
							std::vector<Vertex*>& vrts);

		bool create_pyramids(std::vector<Volume*>& volsOut,
							Grid& grid, rapidxml::xml_node<>* node,
							std::vector<Vertex*>& vrts);

		bool create_octahedrons(std::vector<Volume*>& volsOut,
								Grid& grid, rapidxml::xml_node<>* node,
								std::vector<Vertex*>& vrts);

		template <class TGeomObj>
		bool read_subset_handler_elements(ISubsetHandler& shOut,
										 const char* elemNodeName,
										 rapidxml::xml_node<>* subsetNode,
										 int subsetIndex,
										 std::vector<TGeomObj*>& vElems);

		template <class TGeomObj>
		bool read_selector_elements(ISelector& selOut,
									 const char* elemNodeName,
									 rapidxml::xml_node<>* selNode,
									 std::vector<TGeomObj*>& vElems);

		template <class TElem>
		bool read_attachment(Grid& grid, rapidxml::xml_node<>* node);

		SPRefinementProjector
		read_projector(rapidxml::xml_node<>* projNode);

	protected:
	///	the xml_document which stores the data
		rapidxml::xml_document<> m_doc;

	///	holds grids which already have been created
		std::vector<GridEntry>	m_entries;
};


class UGXFileInfo{
	public:
		UGXFileInfo();

		bool parse_file(const char* filename);

		size_t num_grids() const;
		size_t num_subset_handlers(size_t gridInd) const;
		size_t num_subsets(size_t gridInd, size_t shInd) const;

		std::string grid_name(size_t gridInd) const;
		std::string subset_handler_name(size_t gridInd, size_t shInd) const;
		std::string subset_name(size_t gridInd, size_t shInd, size_t subsetInd) const;

		bool grid_has_vertices(size_t gridInd) const;
		bool grid_has_edges(size_t gridInd) const;
		bool grid_has_faces(size_t gridInd) const;
		bool grid_has_volumes(size_t gridInd) const;

	///	Returns the physical dimension of the given grid.
	/** We define the 'maximal range' as the maximum of the ranges of the
	 * 	the particular coordinates. Then the result is
	 * 	3 - iff the z-coordinate are in a range that is larger than
	 * 		SMALL times the maximal range;
	 * 	2 - iff it is not 3 and the y-coordinate are in a range that is
	 * 		larger than SMALL times the maximal range;
	 * 	1 - iff it is not 0 or 1 and the x-coordinate are in a range that is
	 * 		larger than SMALL times the maximal range;
	 * 	0 - iff it is not 3 or 2 or 1 (i.e. iff the geometry resides in one point).
	 */
		size_t physical_grid_dimension(size_t gridInd) const;

	///	Returns the topological dimension of the given grid.
	/** That is the dimension of the element of highest dimension in the given grid.
	 */
		size_t topological_grid_dimension(size_t gridInd) const;

	///	Returns the dimension of the world-coordinates required for the given grid.
	/** We define the 'maximal range' as the maximum of the ranges of the
	 * 	the particular coordinates. Then the result is
	 * 	3 - iff the z-coordinate are in a range that is larger than
	 * 		SMALL times the maximal range;
	 * 	2 - iff it is not 3 and the y-coordinate are in a range that is
	 * 		larger than SMALL times the maximal range;
	 * 	1 - iff it is not 0 or 1 and the x-coordinate are in a range that is
	 * 		larger than SMALL times the maximal range;
	 * 	0 - iff it is not 3 or 2 or 1 (i.e. iff the geometry resides in one point).
	 *
	 *	@note		The functionality of this method has been changed slightly as of
	 *				2015-03-13 and is the same as in physical_grid_dimension(size_t gridInd);
	 *				the previous functionality is still available in
	 *				topological_grid_dimension(size_t gridInd).
	 * 	@deprecated	This method is marked deprecated and might be removed
	 * 				in a future update.
	 *				Please use physical_grid_dimension(size_t gridInd) and
	 *				topological_grid_dimension(size_t gridInd) instead.
	 */
		size_t grid_world_dimension(size_t gridInd) const;


	private:
		struct SubsetInfo{
			std::string	m_name;
		};

		struct SubsetHandlerInfo{
			std::string	m_name;
			std::vector<SubsetInfo>	m_subsets;
		};

		struct GridInfo{
			std::string	m_name;
			bool	m_hasVertices;
			bool	m_hasEdges;
			bool	m_hasFaces;
			bool	m_hasVolumes;
			std::vector<SubsetHandlerInfo>	m_subsetHandlers;
			vector3 m_extension;
		};

		std::vector<GridInfo>	m_grids;
		bool					m_fileParsed;


	private:
	///	returns the name-attribute of the node or "" if none exists
		std::string node_name(rapidxml::xml_node<>* n) const;

	///	throws an error if no file has been parsed yet
		void check_file_parsed() const;

	///	return the queried grid info and throws an error if the grid index is out of range
	/**	Also calls check_file_parsed().*/
		const GridInfo& grid_info(size_t index) const;

	///	throws an error if the subset handler index is out of range
	/**	Also calls check_grid_index.*/
		const SubsetHandlerInfo& subset_handler_info(size_t gridInd, size_t shInd) const;

	///	throws an error if the subset index is out of range.
	/**	Also calls check_subset_handler_index.*/
		const SubsetInfo& subset_info(size_t gridInd, size_t shInd, size_t subsetInd) const;

	/// calculates the bounding box of a group of vertices
	/**
	 *
	 * @param[in] vrtNode	node in the xml file (containing vertex information)
	 * @param[out] bb		output bounding box
	 *
	 * @return true iff at least one valid (coordinate dimension in {0,1,2,3}) vertex is contained
	 */
		bool calculate_vertex_node_bbox(rapidxml::xml_node<>* vrtNode, AABox<vector3>& bb) const;
};

}//	end of namespace

////////////////////////////////
//	include implementation
#include "file_io_ugx_impl.hpp"

#endif
