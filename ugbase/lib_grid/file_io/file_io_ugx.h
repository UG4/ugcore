//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m06 d16

#include <sstream>
#include <iostream>
#include <vector>
#include "../lg_base.h"
#include "../common_attachments.h"
#include "common/parser/rapidxml/rapidxml.hpp"

#ifndef __H__LIB_GRID__FILE_IO_UGX__
#define __H__LIB_GRID__FILE_IO_UGX__

namespace ug
{

////////////////////////////////////////////////////////////////////////
///	Grants access to ugx files.
class FileAccessorUGX
{
	public:
		FileAccessorUGX();
		virtual ~FileAccessorUGX();
		
	/**	TPositionAttachments value type has to be compatible with MathVector.
	 *	Make sure that aPos is attached to the vertices of the grid.*/
		template <class TPositionAttachment>
		bool add_grid(Grid& grid, const char* name,
					  TPositionAttachment& aPos);
					  
		template <class TPositionAttachment>
		void add_grid(MultiGrid* mg, const char* name,
					  TPositionAttachment& aPos);
					  		
		void add_subset_handler(SubsetHandler* sh, const char* name,
								size_t refGridIndex);
								
		void add_subset_handler(MGSubsetHandler* mgsh, const char* name,
								size_t refGridIndex);
		
		template <class TAttachment>
		void add_vertex_attachment(const TAttachment& attachment,
									const char* name,
									size_t refGridIndex);
									
		
		virtual bool write_to_stream(std::ostream& out);
		
	protected:
		template <class TAAPos>
		rapidxml::xml_node<>*
		create_vertex_node(VertexBaseIterator vrtsBegin,
						  VertexBaseIterator vrtsEnd,
						  TAAPos& aaPos);

	protected:
		rapidxml::xml_document<> m_doc;
};


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
template <class TAAPos>
rapidxml::xml_node<>*
FileAccessorUGX::
create_vertex_node(VertexBaseIterator vrtsBegin,
				  VertexBaseIterator vrtsEnd,
				  TAAPos& aaPos)
{
	using namespace rapidxml;
	using namespace std;
//	the number of coordinates
	const size_t numCoords = TAAPos::ValueType::Size;
	
//	write the vertices to a temporary stream
	stringstream ss;
	for(VertexBaseIterator iter = vrtsBegin; iter != vrtsEnd; ++iter)
	{
		for(size_t i = 0; i < numCoords; ++i)
			ss << aaPos[*iter][i] << " ";
	}
	
//	allocate a string
	char* nodeData = m_doc.allocate_string(ss.str().c_str(), ss.str().size());
	
//	create a vertex-node
	xml_node<>* node = m_doc.allocate_node(node_element, "vertex", nodeData);

	char buff[10];
	sprintf(buff, "%d", numCoords);
	node->append_attribute(m_doc.allocate_attribute("coords", buff));
	
//	return the node
	return node;
}

template <class TPositionAttachment>
bool FileAccessorUGX::
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
	
//	append it to the document
	m_doc.append_node(gridNode);
	
//	write vertices
	gridNode->append_node(create_vertex_node(grid.vertices_begin(),
										grid.vertices_end(), aaPos));
	return true;
}
			  
template <class TPositionAttachment>
void FileAccessorUGX::
add_grid(MultiGrid* mg, const char* name,
		 TPositionAttachment& aPos)
{
}

template <class TAttachment>
void FileAccessorUGX::
add_vertex_attachment(const TAttachment& attachment,
						const char* name,
						size_t refGridIndex)
{
}
}//	end of namespace

#endif
