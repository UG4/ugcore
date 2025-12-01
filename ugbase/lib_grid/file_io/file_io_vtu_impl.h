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

#ifndef __H__UG_file_io_vtu_impl
#define __H__UG_file_io_vtu_impl

#include <sstream>
#include <cstring>
#include "lib_grid/algorithms/debug_util.h"
#include "lib_grid/callbacks/subset_callbacks.h"

namespace ug{


template <typename TAPosition>
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

	vtuReader.subset_handler(sh, 0, 0);

	return true;
}

template <typename TAPosition>
bool SaveGridToVTU(Grid& grid, ISubsetHandler* psh, const char* filename,
				   TAPosition& aPos)
{
	GridWriterVTU vtuWriter;
	std::ofstream out(filename);

	vtuWriter.set_stream(&out);

	vtuWriter.new_piece(grid, psh, aPos);

	vtuWriter.finish();

	return true;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//	GridWriterVTU
////////////////////////////////////////////////////////////////////////////////

std::ostream& GridWriterVTU::
out_stream()
{
	UG_COND_THROW(!m_pout, "Invalid ostream specified fo GridWriterVTU!");
	return *m_pout;
}

template <typename TPositionAttachment>
bool GridWriterVTU::
new_piece(Grid& grid, ISubsetHandler* psh, TPositionAttachment& aPos)
{
	using namespace std;

	if(m_pieceMode == Mode::OPEN)
		end_piece();

	m_pieceMode = Mode::OPEN;
	m_pointDataMode = Mode::NONE;
	m_cellDataMode = Mode::NONE;

	ostream& out = out_stream();

	m_curGrid = &grid;
	m_curSH = psh;
	
	m_pieceSubsetHandlers.clear();
	m_cells.clear();

//	if a subset-handler is specified, all grid elements which are assigned to
//	subsets are considered as cells.
//	If not, only elements of highest dimension are considered as cells.

	if(psh){
		add_subset_handler(*psh, string("regions"));
		collect_cells<Vertex>(m_cells, grid, IsNotInSubset(*psh, -1));
		collect_cells<Edge>(m_cells, grid, IsNotInSubset(*psh, -1));
		collect_cells<Face>(m_cells, grid, IsNotInSubset(*psh, -1));
		collect_cells<Volume>(m_cells, grid, IsNotInSubset(*psh, -1));
	}
	else{
		if(grid.num<Volume>() > 0)
			collect_cells<Volume>(m_cells, grid, ConsiderAll());
		else if(grid.num<Face>() > 0)
			collect_cells<Face>(m_cells, grid, ConsiderAll());
		else if(grid.num<Edge>() > 0)
			collect_cells<Edge>(m_cells, grid, ConsiderAll());
		else
			collect_cells<Vertex>(m_cells, grid, ConsiderAll());
	}

	out << "    <Piece NumberOfPoints=\"" << grid.num<Vertex>() << "\""
		  << " NumberOfCells=\"" << m_cells.size() << "\">" << endl;


//	write points
	out << "      <Points>" << endl;
	write_vector_data<Vertex>(grid, aPos, "");
	out << "      </Points>" << endl;

//	write cells
//	first we'll assign indices to the vertices, which can then be used as
//	references by the cells.
	AInt aInd;
	grid.attach_to_vertices(aInd);
	Grid::VertexAttachmentAccessor<AInt> aaInd(grid, aInd);
	AssignIndices(grid.begin<Vertex>(), grid.end<Vertex>(), aaInd, 0);
	
	write_cells(m_cells, grid, aaInd);

	grid.detach_from_vertices(aInd);

	return true;
}


template <typename TElem, typename TAttachment>
void GridWriterVTU::
write_vector_data(Grid& grid,
				  TAttachment aData,
				  const char* name,
				  typename Grid::traits<TElem>::callback consider_elem)
{
	using namespace std;
	ostream& out = out_stream();

	using vector_t = typename TAttachment::ValueType;
	using iter_t = typename Grid::traits<TElem>::iterator;

	Grid::AttachmentAccessor<TElem, TAttachment> aaData(grid, aData);

	write_data_array_header("Float32", "", vector_t::Size);
	out << "         ";

	for(iter_t i = grid.begin<TElem>(); i != grid.end<TElem>(); ++i){
		vector_t& v = aaData[*i];
		for(size_t c = 0; c < vector_t::Size; ++c){
			out << " " << v[c];
		}
	}
	out << endl;
	write_data_array_footer();
}

template <typename TElem>
void GridWriterVTU::
collect_cells(std::vector<GridObject*>& cellsOut, Grid& grid,
			  typename Grid::traits<TElem>::callback consider_elem)
{
	using iter_t = typename Grid::traits<TElem>::iterator;
	for(iter_t i = grid.begin<TElem>(); i != grid.end<TElem>(); ++i){
		if(consider_elem(*i))
			cellsOut.push_back(*i);
	}
}




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//	GridReaderVTU
////////////////////////////////////////////////////////////////////////////////
template <typename TPositionAttachment>
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
	grid.set_options(GridOptions::GRIDOPT_NONE);

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
	xml_node<>* pieceNode = gridNode;
//	first we'll create all points and cells, then we'll parse point- and cell-data
	xml_node<>* pointsNode = pieceNode->first_node("Points");
	UG_COND_THROW(pointsNode == nullptr, "Missing Points node in UnstructuredGrid node!")

	size_t vrtOffset = vertices.size();
	create_vertices(vertices, grid, pointsNode, aaPos);

	xml_node<>* cellsNode = pieceNode->first_node("Cells");
	UG_COND_THROW(cellsNode == nullptr, "Missing Cells node in UnstructuredGrid node!")
	create_cells(cells, grid, cellsNode, vertices, vrtOffset);

//	reenable the grids options.
	grid.set_options(gridopts);

	return true;
}

template <typename TAAPos>
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

	int numDestCoords = static_cast<int>(TAAPos::ValueType::Size);

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


template <typename T>
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

template <typename T>
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

#endif