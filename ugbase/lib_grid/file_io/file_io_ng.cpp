//	created by Sebastian Reiter, Nicolas Tessore
//	s.b.reiter@googlemail.com
//	y08 m12 d02

#include "file_io_ng.h"

#include <vector>

#include "lib_grid/lg_base.h"

extern "C" {
#include <lib_grid/externals/include/ng/ng.h>
}

using namespace std;

enum ng_volume_type
{
	ng_tetrahedra = 4,
	ng_pyramid = 5,
	ng_prism = 6,
	ng_hexahedra = 8
};

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	ImportGridFromNG
bool ImportGridFromNG(Grid& grid,
                      const char* filename,
                      AVector3& aPos,
                      ISubsetHandler* pSubdomainHandler)
{
	// create ng object
	ng* n = ng_new();

	// read ng file
	ng_info* ninfo = ng_info_new();
	if(ng_read(filename, n, ninfo))
	{
		LOG("WARNING in ImportGridFromNG: " << ninfo->err_msg << endl);
		ng_info_delete(ninfo);
		ng_delete(n);
		return false;
	}
	ng_info_delete(ninfo);

	//	set up vertex attachment
	if(!grid.has_vertex_attachment(aPos))
		grid.attach_to_vertices(aPos);

	//	set up vertex attachment accessor
	Grid::VertexAttachmentAccessor<AVector3> aaPosition(grid, aPos);

	//	read nodes and store them in an array for index access
	vector<Vertex*>	vVertices;
	vVertices.reserve(n->num_bnodes + n->num_inodes);

	// read boundary nodes
	for(int i = 0; i < n->num_bnodes; ++i)
	{
		// get boundary node
		ng_bnode* node = &n->bnodes[i];

		// create and store vertex
		Vertex* vert = *grid.create<Vertex>();
		vVertices.push_back(vert);

		// set vertex coordinates
		aaPosition[vert] = vector3(
			(number)node->coords[0],
			(number)node->coords[1],
			(number)node->coords[2]
		);
	}

	// read interior nodes
	for(int i = 0; i < n->num_inodes; ++i)
	{
		// get interior node
		ng_inode* node = &n->inodes[i];

		// create and store vertex
		Vertex* vert = *grid.create<Vertex>();
		vVertices.push_back(vert);

		// set vertex coordinates
		aaPosition[vert] = vector3(
			(number)node->coords[0],
			(number)node->coords[1],
			(number)node->coords[2]
		);
	}

	//	read volumes
	for(int i = 0; i < n->num_elements; ++i)
	{
		ng_element* elem = &n->elements[i];

		Volume* vol = NULL;

		// create volume
		switch(elem->num_nodes)
		{
			case ng_tetrahedra:
				vol = *grid.create<Tetrahedron>(TetrahedronDescriptor(
						vVertices[elem->nodes[0]],
						vVertices[elem->nodes[1]],
						vVertices[elem->nodes[2]],
						vVertices[elem->nodes[3]]
				));
				break;

			case ng_pyramid:
				vol = *grid.create<Pyramid>(PyramidDescriptor(
						vVertices[elem->nodes[0]],
						vVertices[elem->nodes[1]],
						vVertices[elem->nodes[2]],
						vVertices[elem->nodes[3]],
						vVertices[elem->nodes[4]]));
				break;
			case ng_prism:
				vol = *grid.create<Prism>(PrismDescriptor(
						vVertices[elem->nodes[0]],
						vVertices[elem->nodes[1]],
						vVertices[elem->nodes[2]],
						vVertices[elem->nodes[3]],
						vVertices[elem->nodes[4]],
						vVertices[elem->nodes[5]]));
				break;
			case ng_hexahedra:
				vol = *grid.create<Hexahedron>(HexahedronDescriptor(
						vVertices[elem->nodes[0]],
						vVertices[elem->nodes[1]],
						vVertices[elem->nodes[2]],
						vVertices[elem->nodes[3]],
						vVertices[elem->nodes[4]],
						vVertices[elem->nodes[5]],
						vVertices[elem->nodes[6]],
						vVertices[elem->nodes[7]]));
				break;
			default:
				LOG("WARNING in ImportGridFromNG: Volume type not implemented!" << endl);
				break;
		}

		// add volume to subset
		if(vol != NULL && pSubdomainHandler != NULL)
			pSubdomainHandler->assign_subset(vol, elem->subdomain - 1);
	}

	// done importing!

	// delete ng object
	ng_delete(n);

	return true;
}

}//	end of namespace
