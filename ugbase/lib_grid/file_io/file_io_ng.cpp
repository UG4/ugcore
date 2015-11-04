#include "file_io_ng.h"

#include <vector>

#include "lib_grid/lg_base.h"

extern "C" {
#include <lib_grid/file_io/externals/include/ng/ng.h>
}

using namespace std;

enum ng_face_type
{
	ng_triangle = 3,
	ng_quadrilateral = 4,
};

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
		ng_info_delete(ninfo);
		ng_delete(n);
		n = ng_new();
		n->dim = 2;
		ninfo = ng_info_new();
		
	//TODO: 2d parsing fails, since it is not fully supported.
	//		There are problems in the element description (awaits F).
		if(ng_read(filename, n, ninfo)){
			LOG("WARNING in ImportGridFromNG: " << ninfo->err_msg << endl);
			ng_info_delete(ninfo);
			ng_delete(n);
			return false;
		}
	}
	ng_info_delete(ninfo);

	//	set up vertex attachment
	if(!grid.has_vertex_attachment(aPos))
		grid.attach_to_vertices(aPos);

	//	set up vertex attachment accessor
	Grid::VertexAttachmentAccessor<AVector3> aaPosition(grid, aPos);

	//	read nodes and store them in an array for index access
	vector<RegularVertex*>	vVertices;
	vVertices.reserve(n->num_bnodes + n->num_inodes);

	// read boundary nodes
	for(int i = 0; i < n->num_bnodes; ++i)
	{
		// get boundary node
		ng_bnode* node = &n->bnodes[i];

		// create and store vertex
		RegularVertex* vert = *grid.create<RegularVertex>();
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
		RegularVertex* vert = *grid.create<RegularVertex>();
		vVertices.push_back(vert);

		// set vertex coordinates
		aaPosition[vert] = vector3(
			(number)node->coords[0],
			(number)node->coords[1],
			(number)node->coords[2]
		);
	}

//	if we're in 2d, set all z-coords to 0
	if(n->dim == 2){
		for(size_t i = 0; i < vVertices.size(); ++i)
			aaPosition[vVertices[i]].z() = 0;
	}

//	create the elements
	if(n->dim == 2){
		//	read faces
		for(int i = 0; i < n->num_elements; ++i)
		{
			ng_element* elem = &n->elements[i];

			Face* face = NULL;

			// create face
			switch(elem->num_nodes)
			{
				case ng_triangle:
					face = *grid.create<Triangle>(TriangleDescriptor(
							vVertices[elem->nodes[0]],
							vVertices[elem->nodes[1]],
							vVertices[elem->nodes[2]]));
					break;

				case ng_quadrilateral:
					face = *grid.create<Quadrilateral>(QuadrilateralDescriptor(
							vVertices[elem->nodes[0]],
							vVertices[elem->nodes[1]],
							vVertices[elem->nodes[2]],
							vVertices[elem->nodes[3]]));
					break;
				default:
					LOG("WARNING in ImportGridFromNG: Face type not implemented!" << endl);
					break;
			}

			// add face to subset
			if(face != NULL && pSubdomainHandler != NULL)
				pSubdomainHandler->assign_subset(face, elem->subdomain - 1);
		}
	}
	else{
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
	}
	// done importing!

	// delete ng object
	ng_delete(n);

	return true;
}

}//	end of namespace
