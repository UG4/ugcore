// created by Sebastian Reiter
// s.b.reiter@gmail.com
// 16.11.2012 (d,m,y)

#include "std_hnode_adjuster.h"
#include "lib_grid/algorithms/debug_util.h"

namespace ug{

StdHNodeAdjuster::AdjustRetVal StdHNodeAdjuster::
ref_marks_changed(IRefiner& ref,
			   	  const std::vector<VertexBase*>& vrts,
			   	  const std::vector<EdgeBase*>& edges,
			   	  const std::vector<Face*>& faces,
			   	  const std::vector<Volume*>& vols)
{
	if(!ref.grid())
		return CONTINUE_IF_MARKED_NEW;
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
			VertexBase* vrt = vrts[i_vrt];
			if(!vrt->is_constrained())
				continue;
			ConstrainedVertex* hv = dynamic_cast<ConstrainedVertex*>(vrt);
			if(!hv)
				continue;

		//	make sure that all parents are marked
			GeometricObject* co = hv->get_constraining_object();
			if(co){
				if(EdgeBase* e = dynamic_cast<EdgeBase*>(co)){
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
		EdgeBase* e = edges[i_edge];

	//	check whether hangingNodeOrder1 is enabled. If so, we have to check
	//	for associated hanging vertices and push them to qHVrts.
		if(node_dependency_order_1_enabled()){
			for(size_t i = 0; i < 2; ++i){
				if(e->vertex(i)->is_constrained()){
					ref.mark(e->vertex(i));
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
						  << " At: " << GetGeometricObjectCenter(grid, e));
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

	//	check whether hangingNodeOrder1 is enabled. If so, we have to check
	//	for associated hanging vertices and push them to qHVrts.
		if(node_dependency_order_1_enabled()){
			for(size_t i = 0; i < f->num_vertices(); ++i){
				if(f->vertex(i)->is_constrained())
					ref.mark(f->vertex(i));
			}
		}

	//	if the face is not marked anisotropic, then
	//	we have to make sure that all associated edges are marked.
		if(ref.get_mark(f) != RM_ANISOTROPIC){
			grid.associated_elements(assEdges, f);
			for(size_t i = 0; i < assEdges.size(); ++i){
				ref.mark(assEdges[i]);
			}
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
						  << " At: " << GetGeometricObjectCenter(grid, f));
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

		if(ref.get_mark(v) != RM_ANISOTROPIC){
		//	we have to make sure that all associated edges and faces are marked.
			grid.associated_elements(assEdges, v);
			for(size_t i = 0; i < assEdges.size(); ++i){
				ref.mark(assEdges[i]);
			}

			grid.associated_elements(assFaces, v);
			for(size_t i = 0; i < assFaces.size(); ++i){
				ref.mark(assFaces[i]);
			}
		}
	}

	return CONTINUE_IF_MARKED_NEW;
}
}// end of namespace
