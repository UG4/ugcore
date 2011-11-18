/*
 * tetrakaidekaeder_generator.cpp
 *
 *  Created on: 18.11.2011
 *      Author: marscher
 */

#include "tetrakaidecahedron_generator.h"

using namespace std;

namespace ug{

void GenerateTetrakaidecahedron(Grid& grid, number& height, number& baseEdgeLength, number& diameter)
{
	std::vector<vector3> positions;
	std::vector<int> indices;


//	fill the arrays
	GenerateTetrakaidecahedron(positions, indices, height, baseEdgeLength, diameter);

//	access vertex position data
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);

//	generate vertices in the grid and store them in an array so that we can index them
	std::vector<VertexBase*> vertices(positions.size());

	for(size_t i = 0; i < positions.size(); ++i){
		VertexBase* v = *grid.create<Vertex>();
		aaPos[v] = positions[i];
		vertices[i] = v;
	}

//	the VolumeDescriptor will be used to create new volumes
	VolumeDescriptor vd;

//	create the elements from the given indices (there should be 2 tets)
	for(size_t i = 0; i < indices.size();){
		int num = indices[i++];
		vd.set_num_vertices(num);
		for(int j = 0; j < num; ++j)
			vd.set_vertex(j, vertices[indices[i++]]);

		switch(num){
			case 4:	grid.create<Tetrahedron>(vd); break;
			case 5:	grid.create<Pyramid>(vd); break;
			case 6:	grid.create<Prism>(vd); break;
			case 8:	grid.create<Hexahedron>(vd); break;
		}
	}
}

void GenerateTetrakaidecahedron(std::vector<vector3>& posOut, std::vector<int>& indsOut,
							  number& height, number& baseEdgeLength, number& diameter)
{

}

}//	end of namespace
