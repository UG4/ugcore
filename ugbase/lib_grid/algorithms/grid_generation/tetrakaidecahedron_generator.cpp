/*
 * tetrakaidekaeder_generator.cpp
 *
 *  Created on: 18.11.2011
 *      Author: marscher
 */

#include "tetrakaidecahedron_generator.h"

using namespace std;

namespace ug{

namespace tkdGenerator{

// global index for nodes
static unsigned int index = 0;

void GenerateTetrakaidecahedron(Grid& grid, number& height, number& baseEdgeLength, number& diameter)
{
	CoordsArray positions;
	IndexArray indices;


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

//	create the elements from the given indices
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

/**
 *
 */
void GenerateTetrakaidecahedron(CoordsArray& posOut, IndexArray& indsOut,
							  number& height, number& baseEdgeLength, number& diameter)
{
	// quantity s, overlap of two aligned tkd's
	// number overlap = 1 / pow(3, 0.5) * (diameter - 2 * baseEdgeLength);
	number a = baseEdgeLength;
	number h = height; // / 3;

	// height of base triangle of top inner prism
	number h_Ap = sqrt(3) * a / 2;

	//works
//	createPrism(vector3(0, 0, 0),
//				vector3(a, 0, 0),
//				vector3(a/2, 0, h_Ap),
//				vector3(0, h, 0),
//				vector3(a, h, 0),
//				vector3(a/2, h, h_Ap),
//				posOut, indsOut);

	// rotate prism around (0, y, 0) with angle theta
	// TODO check why this does not work
	for(number theta = 0; theta < 360; theta += 60) {
		UG_LOG("creating prism with angle " << theta << endl);
		createPrism(vector3(0, 0, 0),
				vector3(a*cos(theta), 0, a*sin(theta)),
				vector3((-2 * h_Ap * sin(theta) + a*cos(theta)) /2, 0, (a*sin(theta) + 2*h_Ap*cos(theta)) /2),
				vector3(0, h, 0),
				vector3(a*cos(theta), h, a*sin(theta)),
				vector3((-2*h_Ap *sin(theta) + a*cos(theta) ) /2, h, (a*sin(theta) + 2*h_Ap *cos(theta)) /2),
				posOut, indsOut);
	}
}

void createPrism(vec3Ref v1, vec3Ref v2, vec3Ref v3,
				 vec3Ref v4, vec3Ref v5, vec3Ref v6,
				 CoordsArray& posOut, IndexArray& indsOut) {
	posOut.push_back(v1);
	posOut.push_back(v2);
	posOut.push_back(v3);

	posOut.push_back(v4);
	posOut.push_back(v5);
	posOut.push_back(v6);

	// 6 nodes
	indsOut.push_back(6);
	// enum each node of prism to assign unique node index
	for(unsigned int current = index; index < current + 6; index++) {
		indsOut.push_back(index);
	}
}

}// end of namespace tkdGenerator
}//	end of namespace ug
