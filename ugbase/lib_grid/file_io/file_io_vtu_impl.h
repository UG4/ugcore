// created by Sebastian Reiter
// s.b.reiter@gmail.com

#ifndef __H__UG_file_io_vtu_impl
#define __H__UG_file_io_vtu_impl

#include <sstream>
#include <cstring>
#include "lib_grid/algorithms/debug_util.h"

namespace ug{


template <class TAPosition>
bool LoadGridFromVTU(Grid& grid, ISubsetHandler& sh, const char* filename,
					 TAPosition& aPos)
{
	GridReaderVTU vtuReader;
	if(!vtuReader.parse_file(filename)){
		UG_LOG("ERROR in LoadGridFromVTU: File not found: " << filename << std::endl);
		return false;
	}

	if(vtuReader.num_grids() < 1){
		UG_LOG("ERROR in LoadGridFromVTU: File contains no grid.\n");
		return false;
	}

	vtuReader.grid(grid, 0, aPos);

	//if(vtuReader.num_subset_handlers(0) > 0)
		vtuReader.subset_handler(sh, 0, 0);

	return true;
}



template <class TPositionAttachment>
bool GridReaderVTU::
grid(Grid& gridOut, size_t index, TPositionAttachment& aPos)
{
	using namespace rapidxml;
	using namespace std;

//	make sure that a node at the given index exists
	if(num_grids() <= index){
		UG_LOG("  GridReaderVTU::read: bad grid index!\n");
		return false;
	}

	Grid& grid = gridOut;

//	Since we have to create all elements in the correct order and
//	since we have to make sure that no elements are created in between,
//	we'll first disable all grid-options and reenable them later on
	uint gridopts = grid.get_options();
	grid.set_options(GRIDOPT_NONE);

//	access node data
	if(!grid.has_vertex_attachment(aPos)){
		grid.attach_to_vertices(aPos);
	}

	Grid::VertexAttachmentAccessor<TPositionAttachment> aaPos(grid, aPos);

//	store the grid in the grid-vector and assign indices to the vertices
	m_entries[index].grid = &grid;

//	get the grid-node and the vertex-vector
	xml_node<>* gridNode = m_entries[index].node;
	vector<Vertex*>& vertices = m_entries[index].vertices;
	vector<GridObject*>& cells = m_entries[index].cells;

//	iterate over all pieces
	xml_node<>* pieceNode = gridNode->first_node("Piece");
	for(; pieceNode; pieceNode = pieceNode->next_sibling()){
	//	first we'll create all points and cells, then we'll parse point- and cell-data
		xml_node<>* pointsNode = pieceNode->first_node("Points");
		UG_COND_THROW(pointsNode == NULL, "Missing Points node in UnstructuredGrid node!")

		size_t vrtOffset = vertices.size();
		create_vertices(vertices, grid, pointsNode, aaPos);

		xml_node<>* cellsNode = pieceNode->first_node("Cells");
		UG_COND_THROW(cellsNode == NULL, "Missing Cells node in UnstructuredGrid node!")		
		create_cells(cells, grid, cellsNode, vertices, vrtOffset);
	}



//	iterate through the nodes in the grid and create the entries
	// xml_node<>* curNode = gridNode->first_node();
	// for(;curNode; curNode = curNode->next_sibling()){
	// 	bool bSuccess = true;
	// 	const char* name = curNode->name();
	// 	if(strcmp(name, "Points") == 0)
	// 	 	bSuccess = create_vertices(vertices, grid, curNode, aaPos);
	// 	// else if(strcmp(name, "constrained_vertices") == 0)
	// 	// 	bSuccess = create_constrained_vertices(vertices, constrainingObjsVRT,
	// 	// 										   grid, curNode, aaPos);
	// 	// else if(strcmp(name, "edges") == 0)
	// 	// 	bSuccess = create_edges(edges, grid, curNode, vertices);
	// 	// else if(strcmp(name, "constraining_edges") == 0)
	// 	// 	bSuccess = create_constraining_edges(edges, grid, curNode, vertices);
	// 	// else if(strcmp(name, "constrained_edges") == 0)
	// 	// 	bSuccess = create_constrained_edges(edges, constrainingObjsEDGE,
	// 	// 										grid, curNode, vertices);
	// 	// else if(strcmp(name, "triangles") == 0)
	// 	// 	bSuccess = create_triangles(faces, grid, curNode, vertices);
	// 	// else if(strcmp(name, "constraining_triangles") == 0)
	// 	// 	bSuccess = create_constraining_triangles(faces, grid, curNode, vertices);
	// 	// else if(strcmp(name, "constrained_triangles") == 0)
	// 	// 	bSuccess = create_constrained_triangles(faces, constrainingObjsTRI,
	// 	// 										grid, curNode, vertices);
	// 	// else if(strcmp(name, "quadrilaterals") == 0)
	// 	// 	bSuccess = create_quadrilaterals(faces, grid, curNode, vertices);
	// 	// else if(strcmp(name, "constraining_quadrilaterals") == 0)
	// 	// 	bSuccess = create_constraining_quadrilaterals(faces, grid, curNode, vertices);
	// 	// else if(strcmp(name, "constrained_quadrilaterals") == 0)
	// 	// 	bSuccess = create_constrained_quadrilaterals(faces, constrainingObjsQUAD,
	// 	// 										grid, curNode, vertices);
	// 	// else if(strcmp(name, "tetrahedrons") == 0)
	// 	// 	bSuccess = create_tetrahedrons(volumes, grid, curNode, vertices);
	// 	// else if(strcmp(name, "hexahedrons") == 0)
	// 	// 	bSuccess = create_hexahedrons(volumes, grid, curNode, vertices);
	// 	// else if(strcmp(name, "prisms") == 0)
	// 	// 	bSuccess = create_prisms(volumes, grid, curNode, vertices);
	// 	// else if(strcmp(name, "pyramids") == 0)
	// 	// 	bSuccess = create_pyramids(volumes, grid, curNode, vertices);
	// 	// else if(strcmp(name, "octahedrons") == 0)
	// 	// 	bSuccess = create_octahedrons(volumes, grid, curNode, vertices);


	// 	if(!bSuccess){
	// 		grid.set_options(gridopts);
	// 		return false;
	// 	}
	// }

//	reenable the grids options.
	grid.set_options(gridopts);

	return true;
}

template <class TAAPos>
bool GridReaderVTU::
create_vertices(std::vector<Vertex*>& vrtsOut, Grid& grid,
				rapidxml::xml_node<>* vrtNode, TAAPos aaPos)
{
	using namespace rapidxml;
	using namespace std;

	xml_node<>* dataNode = vrtNode->first_node("DataArray");
	UG_COND_THROW(!dataNode, "Missing DataArray node in Points node");

	int numSrcCoords = -1;
	xml_attribute<>* attrib = dataNode->first_attribute("NumberOfComponents");
	if(attrib)
		numSrcCoords = atoi(attrib->value());

	int numDestCoords = (int)TAAPos::ValueType::Size;

	assert(numDestCoords > 0 && "bad position attachment type");

	if(numSrcCoords < 1 || numDestCoords < 1)
		return false;

//	create a buffer with which we can access the data
	string str(dataNode->value(), dataNode->value_size());
	stringstream ss(str, ios_base::in);

//	if numDestCoords == numSrcCoords parsing will be faster
	if(numSrcCoords == numDestCoords){
		while(!ss.eof()){
		//	read the data
			typename TAAPos::ValueType v;

			for(int i = 0; i < numSrcCoords; ++i)
				ss >> v[i];

		//	make sure that everything went right
			if(ss.fail())
				break;

		//	create a new vertex
			RegularVertex* vrt = *grid.create<RegularVertex>();
			vrtsOut.push_back(vrt);

		//	set the coordinates
			aaPos[vrt] = v;
		}
	}
	else{
	//	we have to be careful with reading.
	//	if numDestCoords < numSrcCoords we'll ignore some coords,
	//	in the other case we'll add some 0's.
		int minNumCoords = min(numSrcCoords, numDestCoords);
		typename TAAPos::ValueType::value_type dummy = 0;

		while(!ss.eof()){
		//	read the data
			typename TAAPos::ValueType v;

			int iMin;
			for(iMin = 0; iMin < minNumCoords; ++iMin)
				ss >> v[iMin];

		//	ignore unused entries in the input buffer
			for(int i = iMin; i < numSrcCoords; ++i)
				ss >> dummy;

		//	add 0's to the vector
			for(int i = iMin; i < numDestCoords; ++i)
				v[i] = 0;

		//	make sure that everything went right
			if(ss.fail()){
				UG_LOG("GridReaderVTU::create_vertices: Failed to read vertex.\n");
				return false;
			}

		//	create a new vertex
			RegularVertex* vrt = *grid.create<RegularVertex>();
			vrtsOut.push_back(vrt);

		//	set the coordinates
			aaPos[vrt] = v;
		}
	}

	return true;
}


template <class T>
void GridReaderVTU::
read_scalar_data(std::vector<T>& dataOut,
				 rapidxml::xml_node<>* dataNode,
				 bool clearData)
{
	using namespace std;

	if(clearData)
		dataOut.clear();

//	create a buffer with which we can access the data
	string str(dataNode->value(), dataNode->value_size());
	stringstream ss(str, ios_base::in);

	while(!ss.eof()){
	//	read the data
		T d;
		ss >> d;

	//	make sure that everything went right
		if(ss.fail())
			break;

		dataOut.push_back(d);
	}
}

template <class T>
void GridReaderVTU::
check_indices(std::vector<T>& inds, size_t first, size_t num, size_t validSize)
{
	UG_COND_THROW(first + num > inds.size(),
				  "Bad index range encountered during parsing of " << m_filename);

	for(size_t i = first; i < num; ++i){
		UG_COND_THROW(inds[i] >= validSize,
					  "Bad index encountered during parsing of " << m_filename);
	}
}

}//	end of namespace

#endif	//__H__file_io_vtu_impl
