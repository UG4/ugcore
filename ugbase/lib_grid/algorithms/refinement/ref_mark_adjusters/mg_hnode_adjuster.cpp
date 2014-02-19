// created by Sebastian Reiter
// s.b.reiter@gmail.com
// 16.11.2012 (d,m,y)

#include "mg_hnode_adjuster.h"
#include "lib_grid/algorithms/debug_util.h"

namespace ug{

void MGHNodeAdjuster::
ref_marks_changed(IRefiner& ref,
			   	  const std::vector<Vertex*>& vrts,
			   	  const std::vector<EdgeBase*>& edges,
			   	  const std::vector<Face*>& faces,
			   	  const std::vector<Volume*>& vols)
{
	UG_ASSERT(ref.grid(), "A refiner has to operate on a grid, before marks can be adjusted!");
	if(!ref.grid()){
		return;
	}
	
	MultiGrid* pmg = dynamic_cast<MultiGrid*>(ref.grid());
	UG_ASSERT(pmg, "MGHNodeAdjuster can only operate on multi-grids, not on standard grids.");
	if(!pmg)
		return;
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


//	since we have to make sure that surface elements which meet at a vertex have
//	a level-distance of at most 1, we now mark parent vertices of marked vertices
//	which are connected to a constrained edge.
//	Those marked parents are then used to mark associated edges, faces and volumes
	Grid::edge_traits::secure_container assEdges;
	Grid::face_traits::secure_container assFaces;
	Grid::volume_traits::secure_container assVols;
	for(size_t i_vrt = 0; i_vrt < vrts.size(); ++i_vrt){
		Vertex* vrt = vrts[i_vrt];
		if(vrt->is_constrained())
			continue;

		if(mg.num_child_vertices(vrt) > 0){
		//	we have to mark all associated edges, faces and volumes
			mg.associated_elements(assEdges, vrt);
			for(size_t i = 0; i < assEdges.size(); ++i)
				ref.mark(assEdges[i]);
			mg.associated_elements(assFaces, vrt);
			for(size_t i = 0; i < assFaces.size(); ++i)
				ref.mark(assFaces[i]);
			mg.associated_elements(assVols, vrt);
			for(size_t i = 0; i < assVols.size(); ++i)
				ref.mark(assVols[i]);
		}
		else if(ref.get_mark(vrt) != RM_DUMMY){
		//	we don't have to select parents of dummy vertices, since we assume
		//	that the maximum level-distance is 1
			Vertex* parent = dynamic_cast<Vertex*>(mg.get_parent(vrt));
			if(parent)
				ref.mark(parent, RM_DUMMY);
		}
	}
}
}// end of namespace
