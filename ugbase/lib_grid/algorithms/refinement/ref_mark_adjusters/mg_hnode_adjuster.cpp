// created by Sebastian Reiter
// s.b.reiter@gmail.com
// 16.11.2012 (d,m,y)

#include "mg_hnode_adjuster.h"
#include "lib_grid/algorithms/debug_util.h"

namespace ug{

MGHNodeAdjuster::AdjustRetVal MGHNodeAdjuster::
ref_marks_changed(IRefiner& ref,
			   	  const std::vector<VertexBase*>& vrts,
			   	  const std::vector<EdgeBase*>& edges,
			   	  const std::vector<Face*>& faces,
			   	  const std::vector<Volume*>& vols)
{
	UG_ASSERT(ref.grid(), "A refiner has to operate on a grid, before marks can be adjusted!");
	if(!ref.grid()){
		return CONTINUE_IF_MARKED_NEW;
	}
	
	MultiGrid* pmg = dynamic_cast<MultiGrid*>(ref.grid());
	UG_ASSERT(pmg, "MGHNodeAdjuster can only operate on multi-grids, not on standard grids.");
	if(!pmg)
		return CONTINUE_IF_MARKED_NEW;
	MultiGrid& mg = *pmg;

//	select all associated vertices of marked objects,
//	since we have to create new vertices in the next levels of the hierarchies.
//	only vertices which do not already have child vertices are selected.

	for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge){
		EdgeBase* e = edges[i_edge];
		for(size_t i = 0; i < e->num_vertices(); ++i){
			if(!mg.has_children(e->vertex(i)))
				ref.mark(e->vertex(i));
		}
	}

	for(size_t i_face = 0; i_face < faces.size(); ++i_face){
		Face* f = faces[i_face];
		for(size_t i = 0; i < f->num_vertices(); ++i){
			if(!mg.has_children(f->vertex(i)))
				ref.mark(f->vertex(i));
		}
	}

	for(size_t i_vol = 0; i_vol < vols.size(); ++i_vol){
		Volume* v = vols[i_vol];
		for(size_t i = 0; i < v->num_vertices(); ++i){
			if(!mg.has_children(v->vertex(i)))
				ref.mark(v->vertex(i));
		}
	}

	return CONTINUE_IF_MARKED_NEW;
}
}// end of namespace
