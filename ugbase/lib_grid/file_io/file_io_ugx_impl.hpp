//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m06 d16

#ifndef __H__LIB_GRID__FILE_IO_UGX_IMPL__
#define __H__LIB_GRID__FILE_IO_UGX_IMPL__

#include <sstream>

namespace ug
{

////////////////////////////////////////////////////////////////////////
template <class TPositionAttachment>
bool GridWriterUGX::
add_grid(Grid& grid, const char* name,
		 TPositionAttachment& aPos)
{
	using namespace rapidxml;
	using namespace std;
//	access node data
	if(!grid.has_vertex_attachment(aPos)){
		UG_LOG("  position attachment missing in grid " << name << endl);
		return false;
	}
	
	Grid::VertexAttachmentAccessor<TPositionAttachment> aaPos(grid, aPos);
	
//	create a new grid-node
	xml_node<>* gridNode = m_doc.allocate_node(node_element, "grid");
	gridNode->append_attribute(m_doc.allocate_attribute("name", name));
	
//	store the grid and the node in an entry
	m_vEntries.push_back(Entry(&grid, gridNode));
	
//	append it to the document
	m_doc.append_node(gridNode);
	
//	write vertices
	if(grid.num<Vertex>() > 0)
		gridNode->append_node(create_vertex_node(grid.begin<Vertex>(),
										grid.end<Vertex>(), aaPos));
//TODO: write hanging-vertices

//	add the remaining grid elements to the nodes
	add_elements_to_node(gridNode, grid);
	
	return true;
}
			  
template <class TPositionAttachment>
void GridWriterUGX::
add_grid(MultiGrid& mg, const char* name,
		 TPositionAttachment& aPos)
{
}

template <class TAttachment>
void GridWriterUGX::
add_vertex_attachment(const TAttachment& attachment,
						const char* name,
						size_t refGridIndex)
{
}

template <class TAAPos>
rapidxml::xml_node<>*
GridWriterUGX::
create_vertex_node(VertexIterator vrtsBegin,
				  VertexIterator vrtsEnd,
				  TAAPos& aaPos)
{
	using namespace rapidxml;
	using namespace std;
//	the number of coordinates
	const int numCoords = (int)TAAPos::ValueType::Size;
	
//	write the vertices to a temporary stream
	stringstream ss;
	for(VertexIterator iter = vrtsBegin; iter != vrtsEnd; ++iter)
	{
		for(int i = 0; i < numCoords; ++i)
			ss << aaPos[*iter][i] << " ";
	}
	
//	create the node
	xml_node<>* node = NULL;
	
	if(ss.str().size() > 0){
	//	allocate a string and erase last character(' ')
		char* nodeData = m_doc.allocate_string(ss.str().c_str(), ss.str().size());
		nodeData[ss.str().size()-1] = 0;
	//	create a node with some data
		node = m_doc.allocate_node(node_element, "vertices", nodeData);
	}
	else{
	//	create an emtpy node
		node = m_doc.allocate_node(node_element, "vertices");
	}

	char* buff = m_doc.allocate_string(NULL, 10);
	sprintf(buff, "%d", numCoords);
	node->append_attribute(m_doc.allocate_attribute("coords", buff));
	
//	return the node
	return node;
}



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	implementation of GridReaderUGX
template <class TPositionAttachment>
bool GridReaderUGX::
get_grid(Grid& gridOut, size_t index,
		  TPositionAttachment& aPos)
{
	using namespace rapidxml;
	using namespace std;

//	make sure that a node at the given index exists
	if(num_grids() <= index){
		UG_LOG("  GridReaderUGX::read: bad grid index!\n");
		return false;
	}
	
	Grid& grid = gridOut;
	
//	access node data
	if(!grid.has_vertex_attachment(aPos)){
		grid.attach_to_vertices(aPos);
	}
	
	Grid::VertexAttachmentAccessor<TPositionAttachment> aaPos(grid, aPos);
	
//	store the grid in the grid-vector and assign indices to the vertices
	m_entries[index].grid = &grid;
	
//	get the grid-node and the vertex-vector
	xml_node<>* gridNode = m_entries[index].node;
	vector<VertexBase*>& vertices = m_entries[index].vertices;
	
//	iterate through the nodes in the grid and create the entries
	xml_node<>* curNode = gridNode->first_node();
	for(;curNode; curNode = curNode->next_sibling()){
		bool bSuccess = true;
		const char* name = curNode->name();
		if(strcmp(name, "vertices") == 0)
			bSuccess = create_vertices(grid, curNode, vertices, aaPos);
		else if(strcmp(name, "edges") == 0)
			bSuccess = create_edges(grid, curNode, vertices);
		else if(strcmp(name, "triangles") == 0)
			bSuccess = create_triangles(grid, curNode, vertices);
		else if(strcmp(name, "quadrilaterals") == 0)
			bSuccess = create_quadrilaterals(grid, curNode, vertices);
		else if(strcmp(name, "tetrahedrons") == 0)
			bSuccess = create_tetrahedrons(grid, curNode, vertices);
		else if(strcmp(name, "hexahedrons") == 0)
			bSuccess = create_hexahedrons(grid, curNode, vertices);
		else if(strcmp(name, "prisms") == 0)
			bSuccess = create_prisms(grid, curNode, vertices);
		else if(strcmp(name, "pyramids") == 0)
			bSuccess = create_pyramids(grid, curNode, vertices);
								
		if(!bSuccess)
			return false;
	}
	
//TODO: read hanging-vertices

	
	return true;
}

template <class TAAPos>
bool GridReaderUGX::
create_vertices(Grid& grid, rapidxml::xml_node<>* vrtNode,
				std::vector<VertexBase*>& vrts, TAAPos aaPos)
{
	using namespace rapidxml;
	using namespace std;
	
	int numSrcCoords = -1;
	xml_attribute<>* attrib = vrtNode->first_attribute("coords");
	if(attrib)
		numSrcCoords = atoi(attrib->value());
	
	int numDestCoords = (int)TAAPos::ValueType::Size;
	
	assert(numDestCoords > 0 && "bad position attachment type");
	
	if(numSrcCoords < 1 || numDestCoords < 1)
		return false;
	
//	create a buffer with which we can access the data
	string str(vrtNode->value(), vrtNode->value_size());
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
			Vertex* vrt = *grid.create<Vertex>();
			vrts.push_back(vrt);
			
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
			if(ss.fail())
				break;
			
		//	create a new vertex
			Vertex* vrt = *grid.create<Vertex>();
			vrts.push_back(vrt);
			
		//	set the coordinates
			aaPos[vrt] = v;
		}
	}
	
	return true;
}

}//	end of namespace

#endif
