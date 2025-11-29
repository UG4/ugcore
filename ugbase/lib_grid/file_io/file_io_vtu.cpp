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

#include "file_io_vtu.h"
#include <string>
#include <iostream>

using namespace std;
using namespace rapidxml;

namespace ug{

enum VTKCellTypes{
	VTK_VERTEX					= 1,
	VTK_POLY_VERTEX				= 2,
	VTK_LINE					= 3,
	VTK_POLY_LINE				= 4,
	VTK_TRIANGLE				= 5,
	VTK_TRIANGLE_STRIP			= 6,
	VTK_POLYGON					= 7,
	VTK_PIXEL					= 8,
	VTK_QUAD					= 9,
	VTK_TETRA					= 10,
	VTK_VOXEL					= 11,
	VTK_HEXAHEDRON				= 12,
	VTK_WEDGE					= 13,
	VTK_PYRAMID					= 14,

	VTK_QUADRATIC_EDGE			= 21,
	VTK_QUADRATIC_TRIANGLE		= 22,
	VTK_QUADRATIC_QUAD			= 23,
	VTK_QUADRATIC_TETRA			= 24,
	VTK_QUADRATIC_HEXAHEDRON	= 25,

	VTK_NUM_TYPES
};

const char* VTKCellNames[] = {	"UNDEFINED",
								"VERTEX",
								"POLY_VERTEX",
								"LINE",
								"POLY_LINE",
								"TRIANGLE",
								"TRIANGLE_STRIP",
								"POLYGON",
								"PIXEL",
								"QUAD",
								"TETRA",
								"VOXEL",
								"HEXAHEDRON",
								"WEDGE",
								"PYRAMID",
								"UNDEFINED",
								"UNDEFINED",
								"UNDEFINED",
								"UNDEFINED",
								"UNDEFINED",
								"UNDEFINED",
								"QUADRATIC_EDGE",
								"QUADRATIC_TRIANGLE",
								"QUADRATIC_QUAD",
								"QUADRATIC_TETRA",
								"QUADRATIC_HEXAHEDRON"};

const int ugRefObjIdToVTKCellType[] = {
										VTKCellTypes::VTK_VERTEX,
										VTKCellTypes::VTK_LINE,
										VTKCellTypes::VTK_TRIANGLE,
										VTKCellTypes::VTK_QUAD,
										VTKCellTypes::VTK_TETRA,
										VTKCellTypes::VTK_HEXAHEDRON,
										VTKCellTypes::VTK_WEDGE,
										VTKCellTypes::VTK_PYRAMID};

bool LoadGridFromVTU(Grid& grid, ISubsetHandler& sh,
					const char* filename)
{
	if(grid.has_vertex_attachment(aPosition))
		return LoadGridFromVTU(grid, sh, filename, aPosition);
	else if(grid.has_vertex_attachment(aPosition2))
		return LoadGridFromVTU(grid, sh, filename, aPosition2);
	else if(grid.has_vertex_attachment(aPosition1))
		return LoadGridFromVTU(grid, sh, filename, aPosition1);

//	no standard position attachments are available.
//	Attach aPosition and use it.
	grid.attach_to_vertices(aPosition);
	return LoadGridFromVTU(grid, sh, filename, aPosition);
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//	GridWriterVTU
////////////////////////////////////////////////////////////////////////////////
GridWriterVTU::
GridWriterVTU() :
	m_pout(nullptr),
	m_pieceMode(Mode::NONE),
	m_pointDataMode(Mode::NONE),
	m_cellDataMode(Mode::NONE),
	m_curGrid(nullptr)
{

}

GridWriterVTU::
GridWriterVTU(std::ostream* pout) :
	m_pout(nullptr),
	m_pieceMode(Mode::NONE),
	m_pointDataMode(Mode::NONE),
	m_cellDataMode(Mode::NONE),
	m_curGrid(nullptr)
{
	set_stream(pout);
}



void GridWriterVTU::
set_stream(std::ostream* pout)
{
	if(m_pout)
		finish();
	m_pout = pout;

	if(!m_pout)
		return;

	std::ostream& out = *m_pout;
	UG_COND_THROW(!out, "Invalid ostream specified fo GridWriterVTU!");

//	todo: check endianiess for binary writes!
	out << "<?xml version=\"1.0\"?>" << endl;
	out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
	out << "  <UnstructuredGrid>" << endl;
}

void GridWriterVTU::
add_subset_handler(ISubsetHandler& sh, const std::string& name)
{
	UG_COND_THROW(m_cellDataMode == Mode::CLOSED,
				  "A subset handler has to be added before a call to 'end_cell_data'");

	m_pieceSubsetHandlers.push_back(make_pair(&sh, name));
}

void GridWriterVTU::
write_data_array_header(const char* type, const char* name,
						int numberOfComponents)
{
	ostream& out = out_stream();

	out << "        <DataArray type=\"" << type << "\"";
	if(strlen(name) > 0)
		out << " Name=\"" << name << "\"";
	if(numberOfComponents > 0)
		out << " NumberOfComponents=\"" << numberOfComponents << "\"";
	out << " format=\"ascii\">" << endl;
}

void GridWriterVTU::
write_data_array_footer()
{
	out_stream() << "        </DataArray>" << endl;
}

void GridWriterVTU::
begin_point_data()
{
	UG_COND_THROW(m_pointDataMode != Mode::NONE,
				  "begin_point_data can only be called once per piece!");
	
	UG_COND_THROW(m_cellDataMode == Mode::OPEN,
				  "begin_point_data may not be called between "
				  "begin_cell_data and end_cell_data!");

	out_stream() << "      <PointData>" << endl;

	m_pointDataMode = Mode::OPEN;
}

void GridWriterVTU::
end_point_data()
{
	UG_COND_THROW(m_pointDataMode != Mode::OPEN,
				  "end_point_data has to be called after begin_point_data!");
	
	out_stream() << "      </PointData>" << endl;

	m_pointDataMode = Mode::CLOSED;
}


void GridWriterVTU::
begin_cell_data()
{
	UG_COND_THROW(m_cellDataMode != Mode::NONE,
				  "begin_cell_data can only be called once per piece!");
	
	UG_COND_THROW(m_pointDataMode == Mode::OPEN,
				  "begin_cell_data may not be called between "
				  "begin_point_data and end_point_data!");

	out_stream() << "      <CellData>" << endl;

	m_cellDataMode = Mode::OPEN;
}


void GridWriterVTU::
end_cell_data()
{
	UG_COND_THROW(m_cellDataMode != Mode::OPEN,
				  "end_cell_data has to be called after begin_cell_data!");

	ostream& out = out_stream();

	for(size_t i = 0; i < m_pieceSubsetHandlers.size(); ++i){
		ISubsetHandler* sh = m_pieceSubsetHandlers[i].first;
		const string& name = m_pieceSubsetHandlers[i].second;
		write_data_array_header("Int32", name.c_str(), 1);

		out << "         ";
		for(size_t icell = 0; icell < m_cells.size(); ++icell){
			out << " " << sh->get_subset_index(m_cells[icell]);
		}
		out << endl;

		write_data_array_footer();
	}

	out << "      </CellData>" << endl;

	m_cellDataMode = Mode::CLOSED;

//	finally write regionos-info
	for(size_t i = 0; i < m_pieceSubsetHandlers.size(); ++i){
		ISubsetHandler& sh = *m_pieceSubsetHandlers[i].first;
		const string& name = m_pieceSubsetHandlers[i].second;
		out << "      <RegionInfo Name=\"" << name << "\">" << endl;
		
		for(int si = 0; si < sh.num_subsets(); ++si){
			out << "        <Region Name=\"" << sh.subset_info(si).name << "\">";// << endl;
			out << "</Region>" << endl;
			//out << "        </Region>" << endl;
		}

		out << "      </RegionInfo>" << endl;
	}
}

void GridWriterVTU::
write_cells(std::vector<GridObject*>& cells, Grid & grid,
			AAVrtIndex aaInd)
{
	ostream& out = out_stream();

	out << "      <Cells>" << endl;

	write_data_array_header("Int32", "connectivity", 0);
	const char* dataPrefix = "         ";
	out << dataPrefix;

	Grid::vertex_traits::secure_container vrts;
	for(size_t icell = 0; icell < cells.size(); ++icell){
		GridObject* cell = cells[icell];
		grid.associated_elements(vrts, cell);
	//	Prisms have a different vertex order in ug than in vtk.
		if(cell->reference_object_id() == ReferenceObjectID::ROID_PRISM){
			out << " " << aaInd[vrts[1]];
			out << " " << aaInd[vrts[0]];
			out << " " << aaInd[vrts[2]];
			out << " " << aaInd[vrts[4]];
			out << " " << aaInd[vrts[3]];
			out << " " << aaInd[vrts[5]];
		}
		else{
			for(size_t ivrt = 0; ivrt < vrts.size(); ++ivrt){
				out << " " << aaInd[vrts[ivrt]];
			}
		}
	}
	
	out << endl;
	write_data_array_footer();


	write_data_array_header("Int32", "offsets", 0);
	out << dataPrefix << " ";

	size_t offset = 0;
	for(size_t icell = 0; icell < cells.size(); ++icell){
		GridObject* cell = cells[icell];
		grid.associated_elements(vrts, cell);
		offset += vrts.size();
		out << " " << offset;
	}
	
	out << endl;
	write_data_array_footer();


	write_data_array_header("Int8", "types", 0);
	out << dataPrefix << " ";

	for(size_t icell = 0; icell < cells.size(); ++icell){
		GridObject* cell = cells[icell];
		int roid = cell->reference_object_id();
		if(roid >= 0 && roid <= ReferenceObjectID::ROID_PYRAMID){
			out << " " << ugRefObjIdToVTKCellType[roid];
		}
		else{
			UG_THROW("Can't map grid-object with ROID " << roid << " to vtk cell type!");
		}
	}
	
	out << endl;
	write_data_array_footer();

	out << "      </Cells>" << endl;
}


void GridWriterVTU::
end_piece()
{
	if(m_pieceMode == Mode::OPEN){
		if(m_pointDataMode == Mode::NONE)
			begin_point_data();
		if(m_pointDataMode == Mode::OPEN)
			end_point_data();

		if(m_cellDataMode == Mode::NONE)
			begin_cell_data();
		if(m_cellDataMode == Mode::OPEN)
			end_cell_data();

		m_pieceMode = Mode::CLOSED;
		out_stream() << "    </Piece>" << endl;
	}

	m_pieceSubsetHandlers.clear();
	m_cells.clear();
}

void GridWriterVTU::
finish()
{
//	make sure that point and cell data are present in the file
	if(m_pieceMode == Mode::NONE){
		m_pieceMode = Mode::OPEN;
		m_pieceSubsetHandlers.clear();
		m_cells.clear();
		out_stream() << "    <Piece NumberOfPoints=\"0\" NumberOfCells=\"0\">" << endl;
		out_stream() << "		<Points></Points>" << endl;
		out_stream() << "		<Cells></Cells>" << endl;
	}

	if(m_pieceMode == Mode::OPEN)
		end_piece();

	out_stream() << "  </UnstructuredGrid>" << endl;
	out_stream() << "</VTKFile>" << endl;

//	invalidate out-stream
	m_pout = nullptr;
}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//	GridReaderVTU
////////////////////////////////////////////////////////////////////////////////


#if 0
// original function
bool GridReaderVTU::
parse_file(const char* filename)
{
	ifstream in(filename, ios::binary);
	if(!in)
		return false;

	m_filename = filename;

//	get the length of the file
	streampos posStart = in.tellg();
	in.seekg(0, ios_base::end);
	streampos posEnd = in.tellg();
	streamsize size = posEnd - posStart;

//	go back to the start of the file
	in.seekg(posStart);

//	read the whole file en-block and terminate it with 0
	char* fileContent = m_doc.allocate_string(0, size + 1);
	in.read(fileContent, size);
	fileContent[size] = 0;
	in.close();

//	parse the xml-data
	m_doc.parse<0>(fileContent);

//	notify derived classes that a new document has been parsed.
	return new_document_parsed();
}
#else

bool GridReaderVTU::
parse_file(const char* filename)
{
	ifstream in(filename, ios::binary);
	if(!in)
		return false;

	m_filename = filename;

//	get the length of the file
	streampos posStart = in.tellg();
	in.seekg(0, ios_base::end);
	streampos posEnd = in.tellg();
	streamsize size = posEnd - posStart;

//	go back to the start of the file
	in.seekg(posStart);

	char * fileContentOriginal = new char[size + 1];

	in.read(fileContentOriginal,size);
	fileContentOriginal[size] = 0;
	in.close();

	std::string fiCo2Str(fileContentOriginal);

	delete [] fileContentOriginal;
	fileContentOriginal = nullptr;

	std::string regInf("RegionInfo");

	if ( fiCo2Str.find(regInf) == std::string::npos && fiCo2Str.find(m_regionOfInterest) != std::string::npos  )
	{
		// we need to insert the additional string

		std::string regInfLines;

		regInfLines.append( "\n<RegionInfo Name=\"" );

		regInfLines.append( m_regionOfInterest );

		regInfLines.append("\">\n");

		regInfLines.append( "</RegionInfo>" );

		std::string insAft( "</CellData>" );

		size_t cedava = fiCo2Str.find( insAft );

		if( cedava == std::string::npos )
			return false;

		size_t insVal = cedava + insAft.size();

		fiCo2Str.insert( insVal, regInfLines );

	}

	char* fileContent = m_doc.allocate_string(0, fiCo2Str.size() );

	strcpy(fileContent, fiCo2Str.c_str());


//	parse the xml-data
	m_doc.parse<0>(fileContent);

//	notify derived classes that a new document has been parsed.
	return new_document_parsed();
}

#endif



bool GridReaderVTU::
new_document_parsed()
{
//	update entries
	m_entries.clear();

	xml_node<>* vtkNode = m_doc.first_node("VTKFile");
	UG_COND_THROW(!vtkNode, "Specified file is not a valid VTKFile!");

	xml_node<>* ugridNode = vtkNode->first_node("UnstructuredGrid");
	UG_COND_THROW(!ugridNode, "Specified file does not contain an unstructured grid!");

//	iterate through all grids
	xml_node<>* curNode = ugridNode->first_node("Piece");
	while(curNode){
		m_entries.push_back(GridEntry(curNode));
		GridEntry& gridEntry = m_entries.back();


	//	collect associated subset handlers
		xml_node<>* curSHNode = curNode->first_node("RegionInfo");
		while(curSHNode){
			gridEntry.subsetHandlerEntries.push_back(SubsetHandlerEntry(curSHNode));
			curSHNode = curSHNode->next_sibling("RegionInfo");
		}

		curNode = curNode->next_sibling("Piece");
	}

	return true;
}


const char* GridReaderVTU::
get_grid_name(size_t index) const
{
	assert(index < num_grids() && "Bad index!");
	xml_attribute<>* attrib = m_entries[index].node->first_attribute("name");
	if(attrib)
		return attrib->value();
	return "";
}

size_t GridReaderVTU::num_subset_handlers(size_t refGridIndex) const
{
//	access the referred grid-entry
	if(refGridIndex >= m_entries.size()){
		UG_LOG("GridReaderVTU::num_subset_handlers: bad refGridIndex. Aborting.\n");
		return 0;
	}

	return m_entries[refGridIndex].subsetHandlerEntries.size();
}

const char* GridReaderVTU::
get_subset_handler_name(size_t refGridIndex, size_t subsetHandlerIndex) const
{
	assert(refGridIndex < num_grids() && "Bad refGridIndex!");
	const GridEntry& ge = m_entries[refGridIndex];
	assert(subsetHandlerIndex < ge.subsetHandlerEntries.size() && "Bad subsetHandlerIndex!");

	xml_attribute<>* attrib = ge.subsetHandlerEntries[subsetHandlerIndex].node->first_attribute("Name");
	if(attrib)
		return attrib->value();
	return "";
}

bool GridReaderVTU::
subset_handler(ISubsetHandler& shOut,
				size_t refGridIndex,
				size_t subsetHandlerIndex)
{
	if(refGridIndex >= m_entries.size()){
		UG_THROW("GridReaderVTU::subset_handler: bad refGridIndex " << refGridIndex
				 << ". Only " << m_entries.size() << " grids available in the file "
				 << m_filename);
	}

 	GridEntry& gridEntry = m_entries[refGridIndex];

//	get the referenced subset-handler entry
	if(subsetHandlerIndex >= gridEntry.subsetHandlerEntries.size()){
		UG_LOG("GridReaderVTU::subset_handler: bad subsetHandlerIndex."
			   << " Assigning all cells to subset 0.\n");
		vector<GridObject*>& cells = gridEntry.cells;
	 	for(size_t i = 0; i < cells.size(); ++i){
	 		shOut.assign_subset(cells[i], 0);
	 	}
		return false;
	}

	SubsetHandlerEntry& shEntry = gridEntry.subsetHandlerEntries[subsetHandlerIndex];
	shEntry.sh = &shOut;
	string shName = get_subset_handler_name(refGridIndex, subsetHandlerIndex);

//	read region-infos
	xml_node<>* regionNode = shEntry.node->first_node("Region");
	int subsetInd = 0;
	while(regionNode)
	{
	//	set subset info
	//	retrieve an initial subset-info from shOut, so that initialised values are kept.
		SubsetInfo& si = shOut.subset_info(subsetInd);

		xml_attribute<>* attrib = regionNode->first_attribute("Name");
		if(attrib)
			si.name = attrib->value();

		regionNode = regionNode->next_sibling("Region");
		++subsetInd;
	}

//	read subset-index for each cell
	xml_node<>* cellDataNode = gridEntry.node->first_node("CellData");
	UG_COND_THROW(cellDataNode == nullptr, "CellData has to be available if regions are defined!");

	xml_node<>* regionDataNode = find_child_node_by_argument_value(cellDataNode,
																   "DataArray",
																   "Name",
																   shName.c_str());
	UG_COND_THROW(regionDataNode == nullptr, "No Cell-Data-Array with name " << shName
				  << " has been defined!");

	vector<GridObject*>& cells = gridEntry.cells;
	vector<int> subsetIndices;
	read_scalar_data(subsetIndices, regionDataNode, true);

	if( subsetIndices.size() != cells.size() )
	{
		vector<double> subsetIndicesDbl;
		read_scalar_data(subsetIndicesDbl, regionDataNode, true);

		trafoDblVec2Int( subsetIndicesDbl, subsetIndices );
	}

	UG_COND_THROW(subsetIndices.size() != cells.size(),
				  "Mismatch regarding number of cells and number of region-indices!");

	for(size_t i = 0; i < cells.size(); ++i){
		shOut.assign_subset(cells[i], subsetIndices[i]);
	}


// 		attrib = subsetNode->first_attribute("color");
// 		if(attrib){
// 			stringstream ss(attrib->value(), ios_base::in);
// 			for(size_t i = 0; i < 4; ++i)
// 				ss >> si.color[i];
// 		}

// 		attrib = subsetNode->first_attribute("state");
// 		if(attrib){
// 			stringstream ss(attrib->value(), ios_base::in);
// 			size_t state;
// 			ss >> state;
// 			si.subsetState = (uint)state;
// 		}

// 		shOut.set_subset_info(subsetInd, si);

// 	//	read elements of this subset
// 		if(shOut.elements_are_supported(SHE_VERTEX))
// 			read_subset_handler_elements<Vertex>(shOut, "vertices",
// 													 subsetNode, subsetInd,
// 													 gridEntry.vertices);
// 		if(shOut.elements_are_supported(SHE_EDGE))
// 			read_subset_handler_elements<Edge>(shOut, "edges",
// 													 subsetNode, subsetInd,
// 													 gridEntry.edges);
// 		if(shOut.elements_are_supported(SHE_FACE))
// 			read_subset_handler_elements<Face>(shOut, "faces",
// 												 subsetNode, subsetInd,
// 												 gridEntry.faces);
// 		if(shOut.elements_are_supported(SHE_VOLUME))
// 			read_subset_handler_elements<Volume>(shOut, "volumes",
// 												 subsetNode, subsetInd,
// 												 gridEntry.volumes);
// 	//	next subset
// 		subsetNode = subsetNode->next_sibling("subset");
// 		++subsetInd;
// 	}

	return true;
}



// This macro should only be used during cell creation in 'create_cells'
// locInd specified the vertex relative to the current offset
#define VRT(locInd)	vertices[connectivity[curOffset + (locInd)] + pieceVrtOffset]

bool GridReaderVTU::
create_cells(std::vector<GridObject*>& cellsOut,
			 Grid& grid,
			 rapidxml::xml_node<>* cellNode,
			 std::vector<Vertex*> vertices,
			 size_t pieceVrtOffset)
{
	const size_t numPieceVrts = vertices.size() - pieceVrtOffset;

	std::vector<size_t> connectivity;
	std::vector<size_t> offsets;
	std::vector<size_t> types;

	xml_node<>* dataNode = cellNode->first_node("DataArray");
	for(; dataNode; dataNode = dataNode->next_sibling()){
		xml_attribute<>* attrib = dataNode->first_attribute("Name");
		if(!attrib)
			attrib = dataNode->first_attribute("name");
		if(!attrib)
			continue;

		if(strcmp(attrib->value(), "connectivity") == 0){
			read_scalar_data(connectivity, dataNode);
		}
		else if(strcmp(attrib->value(), "offsets") == 0){
			read_scalar_data(offsets, dataNode);
		}
		else if(strcmp(attrib->value(), "types") == 0){
			read_scalar_data(types, dataNode);
		}
	}

	UG_COND_THROW(offsets.size() != types.size(),
				  "VTU parsing error in file " << m_filename
				  << ": There have to be as many cell-offsets "
				  "as there are cell-types in a 'cells' node.");

	// UG_LOG("connectivity:");
	// for(size_t i = 0; i < connectivity.size(); ++i){
	// 	UG_LOG(" " << connectivity[i]);
	// }
	// UG_LOG("\n");

	// UG_LOG("offsets:");
	// for(size_t i = 0; i < offsets.size(); ++i){
	// 	UG_LOG(" " << offsets[i]);
	// }
	// UG_LOG("\n");

	// UG_LOG("types:");
	// for(size_t i = 0; i < types.size(); ++i){
	// 	UG_LOG(" " << types[i]);
	// }
	// UG_LOG("\n");

	size_t curOffset = 0;
	for(size_t icell = 0; icell < types.size(); ++icell){
		try{
			size_t nextOffset = offsets[icell];
			switch(types[icell]){
				case VTKCellTypes::VTK_VERTEX:{
					check_indices(connectivity, curOffset, 1, numPieceVrts);
					cellsOut.push_back(VRT(0));
				}break;

				case VTKCellTypes::VTK_LINE:{
					check_indices(connectivity, curOffset, 2, numPieceVrts);
					RegularEdge* e = *grid.create<RegularEdge>(EdgeDescriptor(VRT(0), VRT(1)));
					cellsOut.push_back(e);
				}break;
				
				case VTKCellTypes::VTK_TRIANGLE:{
					check_indices(connectivity, curOffset, 3, numPieceVrts);
					Triangle* f =
						*grid.create<Triangle>(
							TriangleDescriptor(VRT(0), VRT(1), VRT(2)));
					cellsOut.push_back(f);
				}break;
				
				// case VTK_TRIANGLE_STRIP:{

				// }break;
				
				case VTKCellTypes::VTK_QUAD:{
					check_indices(connectivity, curOffset, 4, numPieceVrts);
					Quadrilateral* f =
						*grid.create<Quadrilateral>(
							QuadrilateralDescriptor(VRT(0), VRT(1), VRT(2), VRT(3)));
					cellsOut.push_back(f);
				}break;
				
				case VTKCellTypes::VTK_TETRA:{
					check_indices(connectivity, curOffset, 4, numPieceVrts);
					Volume* v = *grid.create<Tetrahedron>(
									TetrahedronDescriptor(VRT(0), VRT(1),
														  VRT(2), VRT(3)));
					cellsOut.push_back(v);
				}break;
				
				case VTKCellTypes::VTK_HEXAHEDRON:{
					check_indices(connectivity, curOffset, 8, numPieceVrts);
					Volume* v = *grid.create<Hexahedron>(
									HexahedronDescriptor(VRT(0), VRT(1), VRT(2), VRT(3),
														 VRT(4), VRT(5), VRT(6), VRT(7)));
					cellsOut.push_back(v);
				}break;
				
				case VTKCellTypes::VTK_WEDGE:{
					check_indices(connectivity, curOffset, 6, numPieceVrts);
					Volume* v = *grid.create<Prism>(
									PrismDescriptor(VRT(1), VRT(0), VRT(2),
													VRT(4), VRT(3), VRT(5)));
					cellsOut.push_back(v);
				}break;
				
				case VTKCellTypes::VTK_PYRAMID:{
					check_indices(connectivity, curOffset, 5, numPieceVrts);
					Volume* v = *grid.create<Pyramid>(
									PyramidDescriptor(VRT(0), VRT(1), VRT(2),
													VRT(3), VRT(4)));
					cellsOut.push_back(v);
				}break;
				
				default:{
					UG_THROW("VTU parsing error in file " << m_filename
					  << ": Unsupported cell type encountered. ");
				}break;
			}
			curOffset = nextOffset;
		}
		catch(UGError& err){
			UG_LOG("icell: " << icell << ", types[icell]: " << types[icell] << endl);

			string cellName = "UNDEFINED";
			if(types[icell] > 0 && types[icell] < VTKCellTypes::VTK_NUM_TYPES)
				cellName = VTKCellNames[types[icell]];

			std::stringstream ss;
			ss << "Parsed cell type: " << cellName;
			err.push_msg(ss.str(),__FILE__,__LINE__);
			throw(err);
		}
	}

	return true;
}

void GridReaderVTU::
trafoDblVec2Int( std::vector<double> const & dblVec, std::vector<int> & intVec )
{
	intVec = std::vector<int>();

	for( auto d : dblVec )
	{
		intVec.push_back( static_cast<int>( d ) );
	}
}


xml_node<>* GridReaderVTU::
find_child_node_by_argument_value(rapidxml::xml_node<>* parent,
								  const char* nodeName,
								  const char* argName,
								  const char* argValue)
{
	xml_node<>* curNode = parent->first_node(nodeName);
	while(curNode){
		xml_attribute<>* attrib = curNode->first_attribute("Name");
		if(attrib){
			if(strcmp(attrib->value(), argValue) == 0)
				return curNode;
		}
		curNode = curNode->next_sibling(nodeName);
	}
	return nullptr;
}

std::string GridReaderVTU::m_regionOfInterest = "regions";

}//	end of namespace
