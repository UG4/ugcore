// created by Sebastian Reiter
// s.b.reiter@gmail.com

#include "grid_bridges.h"
#include "lib_grid/grid/grid.h"
#include "lib_grid/grid_objects/grid_objects.h"
#include "lib_grid/multi_grid.h"
using namespace std;

namespace ug{
namespace bridge{

template <class T>
static
bool IsValidPtr(T* o){
	return o != NULL;
}


void RegisterGridBridge_Grid(Registry& reg, string parentGroup)
{
	string grp = parentGroup;
//	Geometric Objects
	reg.add_class_<GridObject>("GridObject", grp);
	reg.add_class_<Vertex, GridObject>("Vertex", grp);
	reg.add_class_<Edge, GridObject>("Edge", grp)
		.add_method("num_vertices", &Edge::num_vertices, grp)
		.add_method("vertex", &Edge::vertex, grp);
	reg.add_class_<Face, GridObject>("Face", grp)
		.add_method("num_vertices", &Face::num_vertices, grp)
		.add_method("vertex", &Face::vertex, grp);
	reg.add_class_<Volume, GridObject>("Volume", grp)
		.add_method("num_vertices", &Volume::num_vertices, grp)
		.add_method("vertex", &Volume::vertex, grp);

	reg.add_function("IsValid", &IsValidPtr<Vertex>, grp);
	reg.add_function("IsValid", &IsValidPtr<Edge>, grp);
	reg.add_function("IsValid", &IsValidPtr<Face>, grp);
	reg.add_function("IsValid", &IsValidPtr<Volume>, grp);

//	Grid
	reg.add_class_<Grid>("Grid", grp)
		.add_constructor()
		.add_method("clear", static_cast<void (Grid::*)()>(&Grid::clear))
		.add_method("clear_geometry", &Grid::clear_geometry)
		.add_method("num_vertices", &Grid::num_vertices)
		.add_method("num_edges", &Grid::num_edges)
		.add_method("num_faces", &Grid::num_faces)
		.add_method("num_triangles", &Grid::num<Triangle>)
		.add_method("num_quadrilaterals", &Grid::num<Quadrilateral>)
		.add_method("num_volumes", &Grid::num_volumes)
		.add_method("num_tetrahedrons", &Grid::num<Tetrahedron>)
		.add_method("num_pyramids", &Grid::num<Pyramid>)
		.add_method("num_prisms", &Grid::num<Prism>)
		.add_method("num_hexahedrons", &Grid::num<Hexahedron>)
		.add_method("reserve_vertices", &Grid::reserve<Vertex>, "", "num")
		.add_method("reserve_edges", &Grid::reserve<Edge>, "", "num")
		.add_method("reserve_faces", &Grid::reserve<Face>, "", "num")
		.add_method("reserve_volumes", &Grid::reserve<Volume>, "", "num")
		.set_construct_as_smart_pointer(true);

//	MultiGrid
	reg.add_class_<MultiGrid, Grid>("MultiGrid", grp)
		.add_constructor()
		.add_method("num_levels", &MultiGrid::num_levels)

		.add_method("num_vertices", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<Vertex>)
		.add_method("num_edges", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<Edge>)
		.add_method("num_faces", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<Face>)
		.add_method("num_triangles", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<Triangle>)
		.add_method("num_quadrilaterals", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<Quadrilateral>)
		.add_method("num_volumes", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<Volume>)
		.add_method("num_tetrahedrons", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<Tetrahedron>)
		.add_method("num_pyramids", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<Pyramid>)
		.add_method("num_prisms", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<Prism>)
		.add_method("num_hexahedrons", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<Hexahedron>)
		.set_construct_as_smart_pointer(true);

//	standard attachments
	reg.add_class_<APosition1>("APosition1");
	reg.add_class_<APosition2>("APosition2");
	reg.add_class_<APosition3>("APosition3");
}

}//	end of namespace
}//	end of namespace
