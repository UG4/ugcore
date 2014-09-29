// created by Sebastian Reiter
// s.b.reiter@gmail.com

#include "file_io_vtu.h"

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


GridReaderVTU::GridReaderVTU()
{
}

GridReaderVTU::~GridReaderVTU()
{
}


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

bool GridReaderVTU::
new_document_parsed()
{
//	update entries
	m_entries.clear();

	xml_node<>* vtkNode = m_doc.first_node("VTKFile");
	UG_COND_THROW(!vtkNode, "Specified file is not a valid VTKFile!");

//	iterate through all grids
	xml_node<>* curNode = vtkNode->first_node("UnstructuredGrid");
	while(curNode){
		m_entries.push_back(GridEntry(curNode));
		GridEntry& gridEntry = m_entries.back();

	//	collect associated subset handlers
		xml_node<>* curSHNode = curNode->first_node("Regions");
		while(curSHNode){
			gridEntry.subsetHandlerEntries.push_back(SubsetHandlerEntry(curSHNode));
			curSHNode = curSHNode->next_sibling("Regions");
		}

		curNode = curNode->next_sibling("UnstructuredGrid");
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

	xml_attribute<>* attrib = ge.subsetHandlerEntries[subsetHandlerIndex].node->first_attribute("name");
	if(attrib)
		return attrib->value();
	return "";
}

bool GridReaderVTU::
subset_handler(ISubsetHandler& shOut,
				size_t subsetHandlerIndex,
				size_t refGridIndex)
{
	if(refGridIndex >= m_entries.size()){
		UG_THROW("GridReaderVTU::subset_handler: bad refGridIndex " << refGridIndex
				 << ". Only " << m_entries.size() << " grids available in the file "
				 << m_filename);
	}

 	GridEntry& gridEntry = m_entries[refGridIndex];


//	todo: use the "Regions" CellData to define regions
 	vector<GridObject*>& cells = gridEntry.cells;
 	for(size_t i = 0; i < cells.size(); ++i){
 		shOut.assign_subset(cells[i], 0);
 	}

// //	get the referenced subset-handler entry
// 	if(subsetHandlerIndex >= gridEntry.subsetHandlerEntries.size()){
// 		UG_LOG("GridReaderVTU::subset_handler: bad subsetHandlerIndex. Aborting.\n");
// 		return false;
// 	}

// 	SubsetHandlerEntry& shEntry = gridEntry.subsetHandlerEntries[subsetHandlerIndex];
// 	shEntry.sh = &shOut;

// 	xml_node<>* subsetNode = shEntry.node->first_node("subset");
// 	size_t subsetInd = 0;
// 	while(subsetNode)
// 	{
// 	//	set subset info
// 	//	retrieve an initial subset-info from shOut, so that initialised values are kept.
// 		SubsetInfo si = shOut.subset_info(subsetInd);

// 		xml_attribute<>* attrib = subsetNode->first_attribute("name");
// 		if(attrib)
// 			si.name = attrib->value();

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

template <class TGeomObj>
bool GridReaderVTU::
read_subset_handler_elements(ISubsetHandler& shOut,
							 const char* elemNodeName,
							 rapidxml::xml_node<>* subsetNode,
							 int subsetIndex,
							 std::vector<TGeomObj*>& vElems)
{
	xml_node<>* elemNode = subsetNode->first_node(elemNodeName);

	while(elemNode)
	{
	//	read the indices
		stringstream ss(elemNode->value(), ios_base::in);

		size_t index;
		while(!ss.eof()){
			ss >> index;
			if(ss.fail())
				continue;

			if(index < vElems.size()){
				shOut.assign_subset(vElems[index], subsetIndex);
			}
			else{
				UG_LOG("Bad element index in subset-node " << elemNodeName <<
						": " << index << ". Ignoring element.\n");
				return false;
			}
		}

	//	get next element node
		elemNode = elemNode->next_sibling(elemNodeName);
	}

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
				case VTK_VERTEX:{
					check_indices(connectivity, curOffset, 1, numPieceVrts);
					cellsOut.push_back(VRT(0));
				}break;

				case VTK_LINE:{
					check_indices(connectivity, curOffset, 2, numPieceVrts);
					RegularEdge* e = *grid.create<RegularEdge>(EdgeDescriptor(VRT(0), VRT(1)));
					cellsOut.push_back(e);
				}break;
				
				case VTK_TRIANGLE:{
					check_indices(connectivity, curOffset, 3, numPieceVrts);
					Triangle* f =
						*grid.create<Triangle>(
							TriangleDescriptor(VRT(0), VRT(1), VRT(2)));
					cellsOut.push_back(f);
				}break;
				
				// case VTK_TRIANGLE_STRIP:{

				// }break;
				
				case VTK_QUAD:{
					check_indices(connectivity, curOffset, 4, numPieceVrts);
					Quadrilateral* f =
						*grid.create<Quadrilateral>(
							QuadrilateralDescriptor(VRT(0), VRT(1), VRT(2), VRT(3)));
					cellsOut.push_back(f);
				}break;
				
				case VTK_TETRA:{
					check_indices(connectivity, curOffset, 4, numPieceVrts);
					Volume* v = *grid.create<Tetrahedron>(
									TetrahedronDescriptor(VRT(0), VRT(1),
														  VRT(2), VRT(3)));
					cellsOut.push_back(v);
				}break;
				
				case VTK_HEXAHEDRON:{
					check_indices(connectivity, curOffset, 8, numPieceVrts);
					Volume* v = *grid.create<Hexahedron>(
									HexahedronDescriptor(VRT(0), VRT(1), VRT(2), VRT(3),
														 VRT(4), VRT(5), VRT(6), VRT(7)));
					cellsOut.push_back(v);
				}break;
				
				case VTK_WEDGE:{
					check_indices(connectivity, curOffset, 6, numPieceVrts);
					Volume* v = *grid.create<Prism>(
									PrismDescriptor(VRT(1), VRT(0), VRT(2),
													VRT(4), VRT(3), VRT(5)));
					cellsOut.push_back(v);
				}break;
				
				case VTK_PYRAMID:{
					check_indices(connectivity, curOffset, 5, numPieceVrts);
					Volume* v = *grid.create<Pyramid>(
									PyramidDescriptor(VRT(1), VRT(0), VRT(2),
													VRT(4), VRT(3)));
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
			if(types[icell] > 0 && types[icell] < VTK_NUM_TYPES)
				cellName = VTKCellNames[types[icell]];

			std::stringstream ss;
			ss << "Parsed cell type: " << cellName;
			err.push_msg(ss.str(),__FILE__,__LINE__);
			throw(err);
		}
	}

	return true;
}

}//	end of namespace
