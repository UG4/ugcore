//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m06 d16

#include "file_io_ugx.h"
#include "common/parser/rapidxml/rapidxml_print.hpp"

using namespace std;
using namespace rapidxml;

namespace ug
{




FileAccessorUGX::FileAccessorUGX()
{
	xml_node<>* decl = m_doc.allocate_node(node_declaration);
	decl->append_attribute(m_doc.allocate_attribute("version", "1.0"));
	decl->append_attribute(m_doc.allocate_attribute("encoding", "utf-8"));
	m_doc.append_node(decl);
}

FileAccessorUGX::~FileAccessorUGX()
{
}
					
void FileAccessorUGX::
add_subset_handler(SubsetHandler* sh, const char* name,
				   size_t refGridIndex)
{
}
						
void FileAccessorUGX::
add_subset_handler(MGSubsetHandler* mgsh, const char* name,
				   size_t refGridIndex)
{
}			


bool FileAccessorUGX::
write_to_stream(std::ostream& out)
{
	out << m_doc;
	return true;
}

}//	end of namespace
