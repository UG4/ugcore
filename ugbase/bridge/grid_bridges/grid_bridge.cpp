/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

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
