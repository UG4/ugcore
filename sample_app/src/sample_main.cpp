//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m10 d06

////////////////////////////////////////////////////////////////////////
//	This sample project is primarily intended to demonstrate how a
//	project that uses ug can be set up.
//	It also shows how the profiler that comes with ug (Shiny Profiler)
//	can be used.
//	Used profiler instructions are:
//		PROFILE_FUNC();
//		PROFILE_BEGIN(sectionName);
//		PROFILE_END();
//		PROFILE_CODE(SaveGridToFile(grid, "sample_triangle.obj"));
//		PROFILER_UPDATE();
//		PROFILER_OUTPUT();
//
//	Code that is executed during PROFILE_CODE will be executed even
//	if the profiler is disabled via SHINY_PROFILER = FALSE.
//	For more information about the profiler take a look at
//	ugbase/common/profiler/
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "lib_grid/lib_grid.h"
#include "common/profiler/profiler.h"

using namespace std;
using namespace ug;

////////////////////////////////////////////////////////////////////////
//	build_geometry
void build_geometry(Grid& grid)
{
//	enables profiling for the whole method.
	PROFILE_FUNC();

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

}

////////////////////////////////////////////////////////////////////////
//	main
int main(int argc, char* argv[])
{
//	creates a triangle and stores it in a .obj file.

//	begin profiling for the main section.
//	main_section is an arbitrary name.
	PROFILE_BEGIN(main_section);

//	create a grid.
	Grid grid;
	
//	build geometry
	build_geometry(grid);

//	save the grid
//	we want to measure the time that SaveGridToFile takes.
	cout << "saving to sample_triangle.obj...\n";
	PROFILE_CODE(SaveGridToFile(grid, "sample_triangle.obj"));

//	done
	cout << "done\n\n";

	PROFILE_END();
	
//	call this for output.
	PROFILER_UPDATE();
	PROFILER_OUTPUT();

	return 0;
}
