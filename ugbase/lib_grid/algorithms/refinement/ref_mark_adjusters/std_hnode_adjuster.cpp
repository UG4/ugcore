// created by Sebastian Reiter
// s.b.reiter@gmail.com
// 16.11.2012 (d,m,y)

#include "./std_hnode_adjuster.h"
#include "lib_grid/tools/periodic_boundary_manager.h"
#include "lib_grid/algorithms/debug_util.h"

namespace ug{

// marks geometric object e for refinement if it is periodic
template <class TElem>
static void mark_if_periodic(IRefiner& ref, TElem* e) {
	if(!ref.grid())
		return;
	if(!ref.grid()->has_periodic_boundaries())
		return;

	PeriodicBoundaryManager& pbm = *ref.grid()->periodic_boundary_manager();

	// ensure element is periodic
	if(!pbm.is_periodic(e))
		return;

	RefinementMark refMark = ref.get_mark(e);
	if(pbm.is_master(e))
	{
		typedef typename PeriodicBoundaryManager::Group<TElem>::SlaveContainer SlaveContainer;
		typedef typename PeriodicBoundaryManager::Group<TElem>::SlaveIterator SlaveIterator;
		SlaveContainer& slaves = *pbm.slaves(e);
		for (SlaveIterator iter = slaves.begin(); iter != slaves.end(); ++iter)
			ref.mark(*iter, refMark);
	}
	else { // is slave
		ref.mark(pbm.master(e), refMark);
	}
}

void StdHNodeAdjuster::
ref_marks_changed(IRefiner& ref,
			   	  const std::vector<Vertex*>& vrts,
			   	  const std::vector<Edge*>& edges,
			   	  const std::vector<Face*>& faces,
			   	  const std::vector<Volume*>& vols)
{
	if(!ref.grid())
		return;
	Grid& grid = *ref.grid();

	Grid::edge_traits::secure_container		assEdges;
	Grid::face_traits::secure_container		assFaces;
	Grid::volume_traits::secure_container 	assVols;

////////////////////////////////
//	VERTICES
	if(node_dependency_order_1_enabled()){
	//	make sure that a hanging node is never constrained by an element which
	//	is has constrained vertices as corners
		for(size_t i_vrt = 0; i_vrt < vrts.size(); ++i_vrt){
			Vertex* vrt = vrts[i_vrt];
			if(!vrt->is_constrained())
				continue;
			ConstrainedVertex* hv = dynamic_cast<ConstrainedVertex*>(vrt);
			if(!hv)
				continue;

		//	make sure that all parents are marked
			GridObject* co = hv->get_constraining_object();
			if(co){
				if(Edge* e = dynamic_cast<Edge*>(co)){
					ref.mark(e);
				}
				else if(Face* f = dynamic_cast<Face*>(co)){
					ref.mark(f);
				}
			}
		}
	}

////////////////////////////////
//	EDGES
	for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge){
		Edge* e = edges[i_edge];

	//	check whether hangingNodeOrder1 is enabled. If so, we have to check
	//	for associated hanging vertices and push them to qHVrts.
		if(node_dependency_order_1_enabled()){
			for(size_t i = 0; i < 2; ++i){
				if(e->vertex(i)->is_constrained()){
					Vertex* v = e->vertex(i);
					ref.mark(v);
				}
			}
		}

	//	depending on the type of the edge, we have to perform different operations
		if(ConstrainedEdge* cde = dynamic_cast<ConstrainedEdge*>(e))
		{
		//	the edge is a constrained edge. Make sure that its constraining edge
		//	or face will be refined.
			if(ConstrainingEdge* cge = dynamic_cast<ConstrainingEdge*>(
										cde->get_constraining_object()))
			{
				ref.mark(cge);
			}
			else if(ConstrainingFace* cgf = dynamic_cast<ConstrainingFace*>(
										cde->get_constraining_object()))
			{
				ref.mark(cgf);
			}
			else{
				UG_ASSERT(grid.is_parallel(),
						  "Constrained edge doesn't have a constraining edge. In "
						  "a serial environment this should always be the case!"
						  << " At: " << GetGridObjectCenter(grid, e));
			}
		}
		else if(ConstrainingEdge* cge = dynamic_cast<ConstrainingEdge*>(e))
		{
		//	associated faces and volumes have to be marked
			if(grid.num_faces() > 0){
				grid.associated_elements(assFaces, cge);
				for(size_t i = 0; i < assFaces.size(); ++i){
					ref.mark(assFaces[i]);
				}
			}

			if(grid.num_volumes() > 0){
				grid.associated_elements(assVols, cge);
				for(size_t i = 0; i < assVols.size(); ++i){
					ref.mark(assVols[i]);
				}
			}
		}
	}

////////////////////////////////
//	FACES
	for(size_t i_face = 0; i_face < faces.size(); ++i_face){
		Face* f = faces[i_face];
		RefinementMark refMark = ref.get_mark(f);
	//	check whether hangingNodeOrder1 is enabled. If so, we have to check
	//	for associated hanging vertices and push them to qHVrts.
		if(node_dependency_order_1_enabled()){
			for(size_t i = 0; i < f->num_vertices(); ++i){
				if(f->vertex(i)->is_constrained())
					ref.mark(f->vertex(i));
			}
		}

	//	we have to make sure that all associated edges are marked.
		grid.associated_elements(assEdges, f);
		for(size_t i = 0; i < assEdges.size(); ++i){
			if(refMark > ref.get_mark(assEdges[i]))
				ref.mark(assEdges[i], refMark);
		}

	//	constrained and constraining faces require special treatment
		if(ConstrainedFace* cdf = dynamic_cast<ConstrainedFace*>(f)){
		//	make sure that its constraining face will be refined
			if(ConstrainingFace* cgf = dynamic_cast<ConstrainingFace*>(
										cdf->get_constraining_object()))
			{
				ref.mark(cgf);
			}
			else{
				UG_ASSERT(grid.is_parallel(),
						  "Constrained face doesn't have a constraining face. In "
						  "a serial environment this should always be the case!"
						  << " At: " << GetGridObjectCenter(grid, f));
			}
		}
		else if(ConstrainingFace* cgf = dynamic_cast<ConstrainingFace*>(f)){
		//	associated volumes have to be marked
			if(grid.num_volumes() > 0){
				grid.associated_elements(assVols, cgf);
				for(size_t i = 0; i < assVols.size(); ++i){
					ref.mark(assVols[i]);
				}
			}
		}
	}

////////////////////////////////
//	VOLUMES
	for(size_t i_vol = 0; i_vol < vols.size(); ++i_vol){
		Volume* v = vols[i_vol];
		RefinementMark refMark = ref.get_mark(v);
	//	we have to make sure that all associated edges and faces are marked.
		grid.associated_elements(assEdges, v);
		for(size_t i = 0; i < assEdges.size(); ++i){
			if(refMark > ref.get_mark(assEdges[i]))
				ref.mark(assEdges[i], refMark);
		}

		grid.associated_elements(assFaces, v);
		for(size_t i = 0; i < assFaces.size(); ++i){
			if(refMark > ref.get_mark(assFaces[i]))
				ref.mark(assFaces[i], refMark);
		}
	}

////////////////////////////////
// Periodic boundaries
	if(grid.has_periodic_boundaries()){
		for(size_t i_vrt = 0; i_vrt < vrts.size(); ++i_vrt)
			mark_if_periodic(ref, vrts[i_vrt]);

		for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge)
			mark_if_periodic(ref, edges[i_edge]);

		for(size_t i_face = 0; i_face < faces.size(); ++i_face)
			mark_if_periodic(ref, faces[i_face]);

		// omit volumes, as these are not meant to be periodic
	}
}
}// end of namespace
