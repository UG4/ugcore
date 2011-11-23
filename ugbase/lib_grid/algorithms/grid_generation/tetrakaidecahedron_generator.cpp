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
	// baseEdgeLength
	number a = baseEdgeLength;
	// height of one level of decomposition which is 1/3 of overall height of the tkd
	number h = height/ 3;
	// quantity s, overlap of two aligned tkd's
	number s = 1 / pow(3, 0.5) * (diameter - 2 * baseEdgeLength);
	// height of base triangle of top inner prism
	number hAp = sqrt(3) * a / 2;


	// create G(Ki -> ObenInnen) = 6 prism with equilateral sites
	// rotate prism around (0, y, 0) with step angle 60° = PI/3 (rad)
	/**
	 * rotation matrix:
	 *  +cos(t)  0  - sin(t)+
        |                   |
        |  0     1     0    |
        |                   |
        +sin(t)  0   cos(t) +
	 */
	for(number t = 0; t < 2*PI; t += PI/3) {
		number x = (-2*hAp * sin(t) + a*cos(t)) /2;
		number z = (a*sin(t) + 2*hAp * cos(t)) /2;

		createPrism(vector3(0, 0, 0),				// invariant
					vector3(a*cos(t), 0, a*sin(t)), // (a, 0, 0)
					vector3(x, 0, z), 				// (a/2, 0, hAp)
					vector3(0, h, 0), 				// invariant
					vector3(a*cos(t), h, a*sin(t)), // (a, h, 0)
					vector3(x, h, z), 				// (a/2, h, hAp)
					posOut, indsOut);
	}

	// TODO create G(Ki -> ObenAussenPr2T)
	for(number t = 0; t < 2*PI; t += (2./3 * PI)) {

	}


	// create G(Ki -> ObenAussenPr):
	// create 3 prisms every 120° or 2/3 PI (rad)
	for(number t = 0; t < 2*PI; t += (2./3 * PI)) {
		// for vertex 1 and v2
		number x = (-2*hAp*sin(t) + a*cos(t))/2;
		number z = (a*sin(t) + 2*hAp*cos(t))/2;

		// for v3
		number xs = ((-2*s - 2*hAp)*sin(t) + a* cos(t))/2;
		number zs = (a*sin(t) + (2*s + 2*hAp)*cos(t))/2;

		// for v4: equal to x, z but with different sign
		number xn = (-2*hAp*sin(t) - a*cos(t))/2;
		number zn = (-a*sin(t) + 2*hAp*cos(t))/2;

		// for v5:
		number xng = (-2*hAp*sin(t) - a*cos(t))/2;
		number zng = (-a*sin(t) + 2*hAp*cos(t))/2;

		// for v6:
		number xns = ((-2*s - 2*hAp)*sin(t) - a*cos(t))/2;
		number zns = (-a*sin(t) + (2*s + 2*hAp)*cos(t))/2;

		createPrism(vector3(x, 0, z), // a/2, 0, hAp
					vector3(x, h, z), // a/2, h, hAp
					vector3(xs, 0, zs), // a/2, 0, hAp + s

					vector3(xn, 0, zn), // -a/2, 0, hAp
					vector3(xng, h, zng), // -a/2, h, hAp
					vector3(xns, 0, zns), // -a/2, 0, hAp + s
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
