//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m06 d16

#include <sstream>
#include "file_io_ugx.h"
#include "common/parser/rapidxml/rapidxml_print.hpp"
#include "lib_grid/algorithms/attachment_util.h"

using namespace std;
using namespace rapidxml;

namespace ug
{




GridWriterUGX::GridWriterUGX()
{
	xml_node<>* decl = m_doc.allocate_node(node_declaration);
	decl->append_attribute(m_doc.allocate_attribute("version", "1.0"));
	decl->append_attribute(m_doc.allocate_attribute("encoding", "utf-8"));
	m_doc.append_node(decl);
}

GridWriterUGX::~GridWriterUGX()
{
//	detach aInt from the vertices of the grid
	for(size_t i = 0; i < m_vEntries.size(); ++i)
		m_vEntries[i].grid->detach_from_vertices(m_aInt);
}

bool GridWriterUGX::
write_to_stream(std::ostream& out)
{
	out << m_doc;
	return true;
}

bool GridWriterUGX::
write_to_file(const char* filename)
{
	ofstream out(filename);
	if(out){
		return write_to_stream(out);
	}
	return false;
}

void GridWriterUGX::
add_subset_attributes(rapidxml::xml_node<>* targetNode,
					  ISubsetHandler& sh, size_t subsetIndex)
{
	stringstream ss;
	SubsetInfo& si = sh.subset_info(subsetIndex);
//	write color
	for(size_t i = 0; i < 4; ++i){
		ss << si.color[i];
		if(i < 3)
			ss << " ";
	}
		
	targetNode->append_attribute(m_doc.allocate_attribute("name", si.name.c_str()));
	targetNode->append_attribute(m_doc.allocate_attribute("color", ss.str().c_str()));
}

void GridWriterUGX::
add_subset_handler(SubsetHandler& sh, const char* name,
					size_t refGridIndex)
{
//	get the node of the referenced grid
	if(refGridIndex >= m_vEntries.size()){
		UG_LOG("GridWriterUGX::add_subset_handler: bad refGridIndex. Aborting.\n");
		return;
	}
	
	xml_node<>* parentNode = m_vEntries[refGridIndex].node;
	
//	create the subset-handler node
	xml_node<>* ndSH = m_doc.allocate_node(node_element, "subset_handler");
	ndSH->append_attribute(m_doc.allocate_attribute("name", name));
	
//	add the subset-handler-node to the grid-node.
	parentNode->append_node(ndSH);
	
//	add the subsets
	for(size_t i = 0; i < sh.num_subsets(); ++i){
		xml_node<>* ndSubset = m_doc.allocate_node(node_element, "subset");
		add_subset_attributes(ndSubset, sh, i);
		ndSH->append_node(ndSubset);
		
	//	add elements
	}
}
						
void GridWriterUGX::
add_subset_handler(MGSubsetHandler& mgsh, const char* name,
					size_t refGridIndex)
{

}
						
void GridWriterUGX::
add_elements_to_node(rapidxml::xml_node<>* node,
					  Grid& grid)
{
//	assign indices to the vertices
	grid.attach_to_vertices(m_aInt);
	Grid::VertexAttachmentAccessor<AInt> aaIndVRT(grid, m_aInt);	
	AssignIndices(grid.begin<Vertex>(), grid.end<Vertex>(), aaIndVRT, 0);
	
//	write edges
	if(grid.num<Edge>() > 0)
		node->append_node(create_edge_node(grid.begin<Edge>(),
										grid.end<Edge>(), aaIndVRT));
										
//TODO: write constrained / constraining edges

//	write triangles
	if(grid.num<Triangle>() > 0)
		node->append_node(create_triangle_node(grid.begin<Triangle>(),
												grid.end<Triangle>(), aaIndVRT));

//	write quadrilaterals
	if(grid.num<Quadrilateral>() > 0)
		node->append_node(create_quadrilateral_node(grid.begin<Quadrilateral>(),
													grid.end<Quadrilateral>(), aaIndVRT));

//	write tetrahedrons
	if(grid.num<Tetrahedron>() > 0)
		node->append_node(create_tetrahedron_node(grid.begin<Tetrahedron>(),
													grid.end<Tetrahedron>(), aaIndVRT));

//	write hexahedrons
	if(grid.num<Hexahedron>() > 0)
		node->append_node(create_hexahedron_node(grid.begin<Hexahedron>(),
													grid.end<Hexahedron>(), aaIndVRT));

}

rapidxml::xml_node<>* GridWriterUGX::
create_edge_node(EdgeIterator edgesBegin,
				 EdgeIterator edgesEnd,
				 AAVrtIndex aaIndVRT)
{
//	write the elements to a temporary stream
	stringstream ss;
	for(EdgeIterator iter = edgesBegin; iter != edgesEnd; ++iter)
	{
		ss << aaIndVRT[(*iter)->vertex(0)] << " " << aaIndVRT[(*iter)->vertex(1)] << " ";
	}
	
	if(ss.str().size() > 0){
	//	allocate a string and erase last character(' ')
		char* nodeData = m_doc.allocate_string(ss.str().c_str(), ss.str().size());
		nodeData[ss.str().size()-1] = 0;
	//	create and return the node
		return m_doc.allocate_node(node_element, "edges", nodeData);
	}
	else{
	//	return an emtpy node
		return m_doc.allocate_node(node_element, "edges");
	}
}

rapidxml::xml_node<>* GridWriterUGX::
create_triangle_node(TriangleIterator trisBegin,
				 	 TriangleIterator trisEnd,
				 	 AAVrtIndex aaIndVRT)
{
//	write the elements to a temporary stream
	stringstream ss;
	for(TriangleIterator iter = trisBegin; iter != trisEnd; ++iter)
	{
		ss << aaIndVRT[(*iter)->vertex(0)] << " " << aaIndVRT[(*iter)->vertex(1)]
			<< " " << aaIndVRT[(*iter)->vertex(2)] << " " ;
	}
	
	if(ss.str().size() > 0){
	//	allocate a string and erase last character(' ')
		char* nodeData = m_doc.allocate_string(ss.str().c_str(), ss.str().size());
		nodeData[ss.str().size()-1] = 0;
	//	create and return the node
		return m_doc.allocate_node(node_element, "triangles", nodeData);
	}
	else{
	//	return an emtpy node
		return m_doc.allocate_node(node_element, "triangles");
	}
}

rapidxml::xml_node<>* GridWriterUGX::
create_quadrilateral_node(QuadrilateralIterator quadsBegin,
						  QuadrilateralIterator quadsEnd,
						  AAVrtIndex aaIndVRT)
{
//	write the elements to a temporary stream
	stringstream ss;
	for(QuadrilateralIterator iter = quadsBegin; iter != quadsEnd; ++iter)
	{
		ss << aaIndVRT[(*iter)->vertex(0)] << " " << aaIndVRT[(*iter)->vertex(1)] << " " 
			<< aaIndVRT[(*iter)->vertex(2)] << " " << aaIndVRT[(*iter)->vertex(3)] << " " ;
	}
	
	if(ss.str().size() > 0){
	//	allocate a string and erase last character(' ')
		char* nodeData = m_doc.allocate_string(ss.str().c_str(), ss.str().size());
		nodeData[ss.str().size()-1] = 0;
	//	create and return the node
		return m_doc.allocate_node(node_element, "quadrilaterals", nodeData);
	}
	else{
	//	return an emtpy node
		return m_doc.allocate_node(node_element, "quadrilaterals");
	}
}

rapidxml::xml_node<>* GridWriterUGX::
create_tetrahedron_node(TetrahedronIterator tetsBegin,
						  TetrahedronIterator tetsEnd,
						  AAVrtIndex aaIndVRT)
{
//	write the elements to a temporary stream
	stringstream ss;
	for(TetrahedronIterator iter = tetsBegin; iter != tetsEnd; ++iter)
	{
		ss << aaIndVRT[(*iter)->vertex(0)] << " " << aaIndVRT[(*iter)->vertex(1)] << " " 
			<< aaIndVRT[(*iter)->vertex(2)] << " " << aaIndVRT[(*iter)->vertex(3)] << " " ;
	}
	
	if(ss.str().size() > 0){
	//	allocate a string and erase last character(' ')
		char* nodeData = m_doc.allocate_string(ss.str().c_str(), ss.str().size());
		nodeData[ss.str().size()-1] = 0;
	//	create and return the node
		return m_doc.allocate_node(node_element, "tetrahedrons", nodeData);
	}
	else{
	//	return an emtpy node
		return m_doc.allocate_node(node_element, "tetrahedrons");
	}
}

rapidxml::xml_node<>* GridWriterUGX::
create_hexahedron_node(HexahedronIterator hexasBegin,
						  HexahedronIterator hexasEnd,
						  AAVrtIndex aaIndVRT)
{
//	write the elements to a temporary stream
	stringstream ss;
	for(HexahedronIterator iter = hexasBegin; iter != hexasEnd; ++iter)
	{
		ss << aaIndVRT[(*iter)->vertex(0)] << " " << aaIndVRT[(*iter)->vertex(1)] << " " 
			<< aaIndVRT[(*iter)->vertex(2)] << " " << aaIndVRT[(*iter)->vertex(3)] << " "
			<< aaIndVRT[(*iter)->vertex(4)] << " " << aaIndVRT[(*iter)->vertex(5)] << " "
			<< aaIndVRT[(*iter)->vertex(6)] << " " << aaIndVRT[(*iter)->vertex(7)] << " ";
	}
	
	if(ss.str().size() > 0){
	//	allocate a string and erase last character(' ')
		char* nodeData = m_doc.allocate_string(ss.str().c_str(), ss.str().size());
		nodeData[ss.str().size()-1] = 0;
	//	create and return the node
		return m_doc.allocate_node(node_element, "hexahedrons", nodeData);
	}
	else{
	//	return an emtpy node
		return m_doc.allocate_node(node_element, "hexahedrons");
	}
}

rapidxml::xml_node<>* GridWriterUGX::
create_prism_node(PrismIterator prismsBegin,
					PrismIterator prismsEnd,
					AAVrtIndex aaIndVRT)
{
//	write the elements to a temporary stream
	stringstream ss;
	for(PrismIterator iter = prismsBegin; iter != prismsEnd; ++iter)
	{
		ss << aaIndVRT[(*iter)->vertex(0)] << " " << aaIndVRT[(*iter)->vertex(1)] << " " 
			<< aaIndVRT[(*iter)->vertex(2)] << " " << aaIndVRT[(*iter)->vertex(3)] << " "
			<< aaIndVRT[(*iter)->vertex(4)] << " " << aaIndVRT[(*iter)->vertex(5)] << " ";
	}
	
	if(ss.str().size() > 0){
	//	allocate a string and erase last character(' ')
		char* nodeData = m_doc.allocate_string(ss.str().c_str(), ss.str().size());
		nodeData[ss.str().size()-1] = 0;
	//	create and return the node
		return m_doc.allocate_node(node_element, "prisms", nodeData);
	}
	else{
	//	return an emtpy node
		return m_doc.allocate_node(node_element, "prisms");
	}
}

rapidxml::xml_node<>* GridWriterUGX::
create_pyramid_node(PyramidIterator pyrasBegin,
					PyramidIterator pyrasEnd,
					AAVrtIndex aaIndVRT)
{
//	write the elements to a temporary stream
	stringstream ss;
	for(PyramidIterator iter = pyrasBegin; iter != pyrasEnd; ++iter)
	{
		ss << aaIndVRT[(*iter)->vertex(0)] << " " << aaIndVRT[(*iter)->vertex(1)] << " " 
			<< aaIndVRT[(*iter)->vertex(2)] << " " << aaIndVRT[(*iter)->vertex(3)] << " "
			<< aaIndVRT[(*iter)->vertex(4)] << " ";
	}
	
	if(ss.str().size() > 0){
	//	allocate a string and erase last character(' ')
		char* nodeData = m_doc.allocate_string(ss.str().c_str(), ss.str().size());
		nodeData[ss.str().size()-1] = 0;
	//	create and return the node
		return m_doc.allocate_node(node_element, "pyramids", nodeData);
	}
	else{
	//	return an emtpy node
		return m_doc.allocate_node(node_element, "pyramids");
	}
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	implementation of GridReaderUGX
GridReaderUGX::GridReaderUGX()
{
}

GridReaderUGX::~GridReaderUGX()
{
}

const char* GridReaderUGX::
get_grid_name(size_t index) const
{
	assert(index <= num_grids() && "Bad index!");
	xml_attribute<>* attrib = m_entries[index].node->first_attribute("name");
	if(attrib)
		return attrib->value();
	return NULL;
}

bool GridReaderUGX::
parse_file(const char* filename)
{
	ifstream in(filename, ios::binary);
	if(!in)
		return false;

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

bool GridReaderUGX::
new_document_parsed()
{
//	update entries
	m_entries.clear();
	
	xml_node<>* curNode = m_doc.first_node("grid");
	while(curNode){
		m_entries.push_back(GridEntry(curNode));
		curNode = curNode->next_sibling("grid");
	}
	
	return true;
}

bool GridReaderUGX::
create_edges(Grid& grid, rapidxml::xml_node<>* node,
			 std::vector<VertexBase*>& vrts)
{
//	create a buffer with which we can access the data
	string str(node->value(), node->value_size());
	stringstream ss(str, ios_base::in);
	
//	read the edges
	int i1, i2;
	while(!ss.eof()){
	//	read the indices
		ss >> i1 >> i2;
		
	//	make sure that the indices are valid
		int maxInd = (int)vrts.size() - 1;
		if(i1 < 0 || i1 > maxInd ||
		   i2 < 0 || i2 > maxInd)
		{
			UG_LOG("  ERROR in GridReaderUGX::create_edges: invalid vertex index.\n");
			return false;
		}
		
	//	create the edge
		grid.create<Edge>(EdgeDescriptor(vrts[i1], vrts[i2]));
	}
	
	return true;
}

bool GridReaderUGX::
create_triangles(Grid& grid, rapidxml::xml_node<>* node,
			 std::vector<VertexBase*>& vrts)
{
//	create a buffer with which we can access the data
	string str(node->value(), node->value_size());
	stringstream ss(str, ios_base::in);
	
//	read the triangles
	int i1, i2, i3;
	while(!ss.eof()){
	//	read the indices
		ss >> i1 >> i2 >> i3;
		
	//	make sure that the indices are valid
		int maxInd = (int)vrts.size() - 1;
		if(i1 < 0 || i1 > maxInd ||
		   i2 < 0 || i2 > maxInd ||
		   i3 < 0 || i3 > maxInd)
		{
			UG_LOG("  ERROR in GridReaderUGX::create_triangles: invalid vertex index.\n");
			return false;
		}
		
	//	create the triangle
		grid.create<Triangle>(TriangleDescriptor(vrts[i1], vrts[i2], vrts[i3]));
	}
	
	return true;
}

bool GridReaderUGX::
create_quadrilaterals(Grid& grid, rapidxml::xml_node<>* node,
			 		  std::vector<VertexBase*>& vrts)
{
//	create a buffer with which we can access the data
	string str(node->value(), node->value_size());
	stringstream ss(str, ios_base::in);
	
//	read the quadrilaterals
	int i1, i2, i3, i4;
	while(!ss.eof()){
	//	read the indices
		ss >> i1 >> i2 >> i3 >> i4;
		
	//	make sure that the indices are valid
		int maxInd = (int)vrts.size() - 1;
		if(i1 < 0 || i1 > maxInd ||
		   i2 < 0 || i2 > maxInd ||
		   i3 < 0 || i3 > maxInd ||
		   i4 < 0 || i4 > maxInd)
		{
			UG_LOG("  ERROR in GridReaderUGX::create_quadrilaterals: invalid vertex index.\n");
			return false;
		}
		
	//	create the quad
		grid.create<Quadrilateral>(QuadrilateralDescriptor(vrts[i1], vrts[i2], vrts[i3], vrts[i4]));
	}
	
	return true;
}

bool GridReaderUGX::
create_tetrahedrons(Grid& grid, rapidxml::xml_node<>* node,
			 		  std::vector<VertexBase*>& vrts)
{
//	create a buffer with which we can access the data
	string str(node->value(), node->value_size());
	stringstream ss(str, ios_base::in);
	
//	read the tetrahedrons
	int i1, i2, i3, i4;
	while(!ss.eof()){
	//	read the indices
		ss >> i1 >> i2 >> i3 >> i4;
		
	//	make sure that the indices are valid
		int maxInd = (int)vrts.size() - 1;
		if(i1 < 0 || i1 > maxInd ||
		   i2 < 0 || i2 > maxInd ||
		   i3 < 0 || i3 > maxInd ||
		   i4 < 0 || i4 > maxInd)
		{
			UG_LOG("  ERROR in GridReaderUGX::create_tetrahedrons: invalid vertex index.\n");
			return false;
		}
		
	//	create the element
		grid.create<Tetrahedron>(TetrahedronDescriptor(vrts[i1], vrts[i2], vrts[i3], vrts[i4]));
	}
	
	return true;
}

bool GridReaderUGX::
create_hexahedrons(Grid& grid, rapidxml::xml_node<>* node,
					std::vector<VertexBase*>& vrts)
{
//	create a buffer with which we can access the data
	string str(node->value(), node->value_size());
	stringstream ss(str, ios_base::in);
	
//	read the hexahedrons
	int i1, i2, i3, i4, i5, i6, i7, i8;
	while(!ss.eof()){
	//	read the indices
		ss >> i1 >> i2 >> i3 >> i4 >> i5 >> i6 >> i7 >> i8;
		
	//	make sure that the indices are valid
		int maxInd = (int)vrts.size() - 1;
		if(i1 < 0 || i1 > maxInd ||
		   i2 < 0 || i2 > maxInd ||
		   i3 < 0 || i3 > maxInd ||
		   i4 < 0 || i4 > maxInd ||
		   i5 < 0 || i5 > maxInd ||
		   i6 < 0 || i6 > maxInd ||
		   i7 < 0 || i7 > maxInd ||
		   i8 < 0 || i8 > maxInd)
		{
			UG_LOG("  ERROR in GridReaderUGX::create_hexahedrons: invalid vertex index.\n");
			return false;
		}
		
	//	create the element
		grid.create<Hexahedron>(HexahedronDescriptor(vrts[i1], vrts[i2], vrts[i3], vrts[i4],
													 vrts[i5], vrts[i6], vrts[i7], vrts[i8]));
	}
	
	return true;
}

bool GridReaderUGX::
create_prisms(Grid& grid, rapidxml::xml_node<>* node,
				std::vector<VertexBase*>& vrts)
{
//	create a buffer with which we can access the data
	string str(node->value(), node->value_size());
	stringstream ss(str, ios_base::in);
	
//	read the hexahedrons
	int i1, i2, i3, i4, i5, i6;
	while(!ss.eof()){
	//	read the indices
		ss >> i1 >> i2 >> i3 >> i4 >> i5 >> i6;
		
	//	make sure that the indices are valid
		int maxInd = (int)vrts.size() - 1;
		if(i1 < 0 || i1 > maxInd ||
		   i2 < 0 || i2 > maxInd ||
		   i3 < 0 || i3 > maxInd ||
		   i4 < 0 || i4 > maxInd ||
		   i5 < 0 || i5 > maxInd ||
		   i6 < 0 || i6 > maxInd)
		{
			UG_LOG("  ERROR in GridReaderUGX::create_prisms: invalid vertex index.\n");
			return false;
		}
		
	//	create the element
		grid.create<Prism>(PrismDescriptor(vrts[i1], vrts[i2], vrts[i3], vrts[i4],
										   vrts[i5], vrts[i6]));
	}
	
	return true;
}

bool GridReaderUGX::
create_pyramids(Grid& grid, rapidxml::xml_node<>* node,
				std::vector<VertexBase*>& vrts)
{
//	create a buffer with which we can access the data
	string str(node->value(), node->value_size());
	stringstream ss(str, ios_base::in);
	
//	read the hexahedrons
	int i1, i2, i3, i4, i5;
	while(!ss.eof()){
	//	read the indices
		ss >> i1 >> i2 >> i3 >> i4 >> i5;
		
	//	make sure that the indices are valid
		int maxInd = (int)vrts.size() - 1;
		if(i1 < 0 || i1 > maxInd ||
		   i2 < 0 || i2 > maxInd ||
		   i3 < 0 || i3 > maxInd ||
		   i4 < 0 || i4 > maxInd ||
		   i5 < 0 || i5 > maxInd)
		{
			UG_LOG("  ERROR in GridReaderUGX::create_pyramids: invalid vertex index.\n");
			return false;
		}
		
	//	create the element
		grid.create<Pyramid>(PyramidDescriptor(vrts[i1], vrts[i2], vrts[i3],
												vrts[i4], vrts[i5]));
	}
	
	return true;
}

}//	end of namespace
