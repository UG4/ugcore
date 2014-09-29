//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com

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

namespace ug
{

////////////////////////////////////////////////////////////////////////
///	Reads a grid to an vtu (vtk) file. internally uses GridReaderVTK.
/**	The position attachment can be specified. Since the type of the
 *	position attachment is a template parameter, MathVector attachments
 * 	of any dimension are supported. Especially ug::aPosition, ug::aPostion2
 *	and ug::aPosition1.
 */
template <class TAPosition>
bool LoadGridFromVTU(Grid& grid, ISubsetHandler& sh,
					const char* filename, APosition& aPos);

///	Reads a grid to an vtu (vtk) file.
/**	Before reading a grid from file, this method searches for the
 *	attached standard position attachment with the highest dimension.
 *	This will be used as position-attachment in a call to the overloaded
 *	version of LoadGridFromVTU.
 *
 *	If no standard attachment is found, aPosition will be attached and used.
 */
bool LoadGridFromVTU(Grid& grid, ISubsetHandler& sh,
					const char* filename);
					


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
		GridReaderVTU();
		virtual ~GridReaderVTU();

	///	parses an xml file
		bool parse_file(const char* filename);

	///	returns the number of grids in the given file
		size_t num_grids() const	{return m_entries.size();}

	///	returns the i-th grid.
	/**	TPositionAttachments value type has to be compatible with MathVector.
	 *	Make sure that a file has already been loaded.*/
		template <class TPositionAttachment>
		bool grid(Grid& gridOut, size_t index, TPositionAttachment& aPos);

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

	protected:
		struct SubsetHandlerEntry
		{
			SubsetHandlerEntry(rapidxml::xml_node<>* n) : node(n), sh(NULL) {}

			rapidxml::xml_node<>* 	node;
			ISubsetHandler*			sh;
		};

		struct GridEntry
		{
			GridEntry(rapidxml::xml_node<>* n) : node(n), grid(NULL), mg(NULL)	{}

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
		template <class TAAPos>
		bool create_vertices(std::vector<Vertex*>& vrtsOut, Grid& grid,
							rapidxml::xml_node<>* vrtNode, TAAPos aaPos);

		bool create_cells(std::vector<GridObject*>& cellsOut,
						  Grid& grid,
						  rapidxml::xml_node<>* node,
						  std::vector<Vertex*> vertices,
						  size_t pieceVrtOffset);

		template <class T>
		void read_scalar_data(std::vector<T>& dataOut,
							  rapidxml::xml_node<>* dataNode,
							  bool clearData = true);

		template <class T>
		void check_indices(std::vector<T>& inds, size_t first, size_t num, size_t max);
		
		template <class TGeomObj>
		bool read_subset_handler_elements(ISubsetHandler& shOut,
										 const char* elemNodeName,
										 rapidxml::xml_node<>* subsetNode,
										 int subsetIndex,
										 std::vector<TGeomObj*>& vElems);

	protected:
	///	stores the file-name
		std::string					m_filename;

	///	the xml_document which stores the data
		rapidxml::xml_document<>	m_doc;

	///	holds grids which already have been created
		std::vector<GridEntry>		m_entries;
};

}//	end of namespace

////////////////////////////////
//	include implementation
#include "file_io_vtu_impl.h"

#endif
