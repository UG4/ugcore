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

/*
*  +cos(t)  0  - sin(t)+
   |                   |
   |  0     1     0    |
   |                   |
   +sin(t)  0   cos(t) +
*/
class rotationMatrix {
public:
	void setAngle(const number& theta) {
		R[0][0] = cos(theta);
		R[0][2] = -sin(theta);
		R[2][0] = sin(theta);
		R[2][2] = cos(theta);
	}

	rotationMatrix(const number& theta) {
			R[0][0] = cos(theta);
			R[0][1] = 0;
			R[0][2] = -sin(theta);

			R[1][0] = 0;
			R[1][1] = 1;
			R[1][2] = 0;

			R[2][0] = sin(theta);
			R[2][1] = 0;
			R[2][2] = cos(theta);
	}

	const vector3& operator*(vector3& v) {
		// temp values
		number x, y, z;
		x = R[0][0] * v[0] + R[0][1] * v[1] + R[0][2] * v[2];
		y = R[1][0] * v[0] + R[1][1] * v[1] + R[1][2] * v[2];
		z = R[2][0] * v[0] + R[2][1] * v[1] + R[2][2] * v[2];

		v[0] = x;
		v[1] = y;
		v[2] = z;
		return v;
	}

protected:
	number R[3][3];
};

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
	number h = height / 3;
	// quantity s, overlap of two aligned tkd's
	number s = 1 / sqrt(3) * (diameter - 2 * baseEdgeLength);
	// height of base triangle of top inner prism
	number g = sqrt(3) * a / 2;
	// TODO whats this again?
	number b = sqrt(3) * s / 2;

	UG_LOG("g: " << g << endl);
	UG_LOG("s: " << s << endl);

	// origin of construction of tkd
	const vector3 origin(30, 0, 30);
	// Rotation matrix with initial rotation angle of 0 (rad)
	rotationMatrix R(0);

	// create G(Ki -> ObenInnen) = 6 prism with equilateral sites
	// rotate prism around (0, y, 0) with step angle 60° = PI/3 (rad)
	for(number t = 0; t < 2*PI; t += PI/3) {
		// switch angle of rotation matrix
		R.setAngle(t);

		// begin in origin
		// copy origin, because operator works in place
		vector3 copy = vector3(origin);
		vector3 v1 = R * copy;
		// create surrounding vertices relative to origin
		vector3 v2 = R * (vector3(a, 0, 0) += origin);
		vector3 v3 = R * (vector3(a/2, 0, g) += origin);
		vector3 v4 = R * (vector3(0, h, 0) += origin);
		vector3 v5 = R * (vector3(a, h, 0) += origin);
		vector3 v6 = R * (vector3(a/2, h, g) += origin);

		createPrism(v1,	v2, v3, v4, v5, v6,	posOut, indsOut);
	}

	// create G(Ki -> ObenAussenPr2T)
	for(number t = 2./3 * PI; t < 2*PI; t += (2./3 * PI)) {
		R.setAngle(t);
		vector3 v1_a = R * (vector3(a/2, 0, g) += origin);
		vector3 v2_a = R * (vector3(a/2, 0, g + s) += origin);
		vector3 v3_a = R * (vector3(a/2 + b, 0, g + s/2) += origin);
		vector3 v4_a = R * (vector3(a/2, h, g) += origin);
		// left tetrahedron right of prism of ObenAussenPr
		createTetrahedron(v1_a, v2_a, v3_a, v4_a, posOut, indsOut);

		//TODO calculate correct offset
//		number offset = a+2*b;
//		vector3 v1_b = R * (vector3(a/2, 0, g) += origin);
//		vector3 v2_b = R * (vector3(a/2, 0, g + offset) += origin);
//		vector3 v3_b = R * (vector3(a/2 + b, 0, g + offset) += origin);
//		vector3 v4_b = R * (vector3(a/2, h, g) += origin);
//		// right
//		createTetrahedron(v1_b, v2_b, v3_b, v4_b, posOut, indsOut);

		// Prism of ObenAussenPr2T
		vector3 v1p = R * (vector3(a/2, 0, g) += origin);
		vector3 v2p = R * (vector3(a/2, h, g) += origin);
		vector3 v3p = R * (vector3(a/2 + b, 0, g + s/2) += origin);
		// OK
		// vorzeichen wechsel x
		vector3 v4p = R * (vector3(a/2 + b + a, 0, g) += origin);
		vector3 v5p = R * (vector3(a/2, h, g) += origin);
		vector3 v6p = R * (vector3(a/2 + b + a, 0, g + s/2) += origin);

		createPrism(v1p, v2p, v3p, v4p, v5p, v6p, posOut, indsOut);

		// create G(Ki -> ObenAussenPr):
		// create 3 prisms every 120° or 2/3 PI (rad)
		vector3 v1 = R * (vector3(a/2, 0, g) += origin);
		vector3 v2 = R * (vector3(a/2, h, g) += origin);
		vector3 v3 = R * (vector3(a/2, 0, g + s) += origin);
		vector3 v4 = R * (vector3(-a/2, 0, g) += origin);
		vector3 v5 = R * (vector3(-a/2, h, g) += origin);
		vector3 v6 = R * (vector3(-a/2, 0, g + s) += origin);

		createPrism(v1, v2, v3, v4, v5, v6, posOut, indsOut);
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


/**
 * please be sure to pass the vertices in the correct order:
 * v1, v2, v3: bottom-vertices in counterclockwise order (if viewed from the top).
 * v4: top
 */
void createTetrahedron(vec3Ref v1, vec3Ref v2, vec3Ref v3, vec3Ref v4,
		CoordsArray& posOut, IndexArray& indsOut) {
		posOut.push_back(v1);
		posOut.push_back(v2);
		posOut.push_back(v3);
		// top
		posOut.push_back(v4);

		// 4 nodes
		indsOut.push_back(4);
		// enum each node of prism to assign unique node index
		for(unsigned int current = index; index < current + 4; index++) {
			indsOut.push_back(index);
		}
}

}// end of namespace tkdGenerator
}//	end of namespace ug
