// created by Sebastian Reiter
// s.b.reiter@gmail.com

#include "grid_bridges.h"
#include "lib_grid/selector.h"
#include "lib_grid/grid_objects/grid_objects.h"

using namespace std;

namespace ug{
namespace bridge{

void RegisterGridBridge_Selector(Registry& reg, string parentGroup)
{
	string grp = parentGroup;

	reg.add_class_<ISelector>("ISelector", grp)
		.add_method("select", static_cast<void (ISelector::*)(Vertex*)>(&ISelector::select<Vertex>))
		.add_method("select", static_cast<void (ISelector::*)(Edge*)>(&ISelector::select<Edge>))
		.add_method("select", static_cast<void (ISelector::*)(Face*)>(&ISelector::select<Face>))
		.add_method("select", static_cast<void (ISelector::*)(Volume*)>(&ISelector::select<Volume>))
		.add_method("deselect", static_cast<void (ISelector::*)(Vertex*)>(&ISelector::deselect<Vertex>))
		.add_method("deselect", static_cast<void (ISelector::*)(Edge*)>(&ISelector::deselect<Edge>))
		.add_method("deselect", static_cast<void (ISelector::*)(Face*)>(&ISelector::deselect<Face>))
		.add_method("deselect", static_cast<void (ISelector::*)(Volume*)>(&ISelector::deselect<Volume>))
		.add_method("is_selected", static_cast<bool (ISelector::*)(Vertex*) const>(&ISelector::is_selected<Vertex>))
		.add_method("is_selected", static_cast<bool (ISelector::*)(Edge*) const>(&ISelector::is_selected<Edge>))
		.add_method("is_selected", static_cast<bool (ISelector::*)(Face*) const>(&ISelector::is_selected<Face>))
		.add_method("is_selected", static_cast<bool (ISelector::*)(Volume*) const>(&ISelector::is_selected<Volume>));

	reg.add_class_<Selector, ISelector>("Selector", grp)
		.add_constructor<void (*)(Grid&)>()
		.add_method("num_vertices", static_cast<size_t (Selector::*)() const>(&Selector::num<Vertex>))
		.add_method("num_edges", static_cast<size_t (Selector::*)() const>(&Selector::num<Edge>))
		.add_method("num_faces", static_cast<size_t (Selector::*)() const>(&Selector::num<Face>))
		.add_method("num_triangles", static_cast<size_t (Selector::*)() const>(&Selector::num<Triangle>))
		.add_method("num_quadrilaterals", static_cast<size_t (Selector::*)() const>(&Selector::num<Quadrilateral>))
		.add_method("num_volumes", static_cast<size_t (Selector::*)() const>(&Selector::num<Volume>))
		.add_method("num_tetrahedrons", static_cast<size_t (Selector::*)() const>(&Selector::num<Tetrahedron>))
		.add_method("num_pyramids", static_cast<size_t (Selector::*)() const>(&Selector::num<Pyramid>))
		.add_method("num_prisms", static_cast<size_t (Selector::*)() const>(&Selector::num<Prism>))
		.add_method("num_hexahedrons", static_cast<size_t (Selector::*)() const>(&Selector::num<Hexahedron>))
		.set_construct_as_smart_pointer(true);
}

}//	end of namespace
}//	end of namespace
