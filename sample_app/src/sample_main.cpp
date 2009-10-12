//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m10 d06

#include <iostream>
#include "lib_grid/lib_grid.h"

using namespace std;
using namespace ug;

////////////////////////////////////////////////////////////////////////
//	main
int main(int argc, char* argv[])
{
//	creates a triangle and stores it in a .obj file.

//	create a grid.
	Grid grid;
	
//	create vertices
	cout << "creating vertices...\n";
	Vertex* v1 = *grid.create<Vertex>();
	Vertex* v2 = *grid.create<Vertex>();
	Vertex* v3 = *grid.create<Vertex>();

//	attach position data
	cout << "attaching position data...\n";
	grid.attach_to_vertices(aPosition);
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);

//	assign coordinates
	cout << "assigning coordinates...\n";
	aaPos[v1] = vector3(-1, -1, -1);
	aaPos[v2] = vector3(1, -1, -1);
	aaPos[v3] = vector3(0, 1, 1);

//	create the triangle
	cout << "creating triangle...\n";
	grid.create<Triangle>(TriangleDescriptor(v1, v2, v3));

//	save the grid
	cout << "saving to sample_triangle.obj...\n";
	SaveGridToFile(grid, "sample_triangle.obj");

//	done
	cout << "done\n";

	return 0;
}
