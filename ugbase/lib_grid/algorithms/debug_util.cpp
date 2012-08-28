// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 13.01.2011 (m,d,y)
 
#include "debug_util.h"
using namespace std;

namespace ug{

void PrintElementNumbers(const GeometricObjectCollection& goc)
{
	UG_LOG("grid element numbers:\n");
	for(size_t i = 0; i < goc.num_levels(); ++i)
	{
		if(goc.num_levels() > 1){
			UG_LOG("level " << i << endl);
		}
		UG_LOG("  vertices total:\t" << goc.num<VertexBase>(i) << endl);
		if(goc.num<VertexBase>() > 0){
			UG_LOG("    normal vrts:\t" << goc.num<Vertex>(i) << endl);
			UG_LOG("    hanging vrts:\t" << goc.num<ConstrainedVertex>(i) << endl);
		}

		UG_LOG("  edges total:\t\t" << goc.num<EdgeBase>(i) << endl);
		if(goc.num<EdgeBase>() > 0){
			UG_LOG("    normal edges:\t" << goc.num<Edge>(i) << endl);
			UG_LOG("    constraining edges:\t" << goc.num<ConstrainingEdge>(i) << endl);
			UG_LOG("    constrained edges:\t" << goc.num<ConstrainedEdge>(i) << endl);
		}

		UG_LOG("  faces total:\t\t" << goc.num<Face>(i) << endl);
		if(goc.num<Face>() > 0){
			UG_LOG("    normal triangles:\t" << goc.num<Triangle>(i) << endl);
			UG_LOG("    constraining tris:\t" << goc.num<ConstrainingTriangle>(i) << endl);
			UG_LOG("    constrained tris:\t" << goc.num<ConstrainedTriangle>(i) << endl);

			UG_LOG("    normal quads:\t" << goc.num<Quadrilateral>(i) << endl);
			UG_LOG("    constraining quads:\t" << goc.num<ConstrainingQuadrilateral>(i) << endl);
			UG_LOG("    constrained quads:\t" << goc.num<ConstrainedQuadrilateral>(i) << endl);
		}

		UG_LOG("  volumes total:\t" << goc.num<Volume>(i) << endl);
		if(goc.num<Volume>() > 0){
			UG_LOG("    tetrahedrons:\t" << goc.num<Tetrahedron>(i) << endl);
			UG_LOG("    pyramids:\t" << goc.num<Pyramid>(i) << endl);
			UG_LOG("    prisms:\t" << goc.num<Prism>(i) << endl);
			UG_LOG("    hexahedrons:\t" << goc.num<Hexahedron>(i) << endl);
		}

		UG_LOG(endl);
	}
}

void PrintGridElementNumbers(Grid& grid)
{
	PrintElementNumbers(grid.get_geometric_objects());
}

void PrintGridElementNumbers(MultiGrid& mg)
{
	PrintElementNumbers(mg.get_geometric_objects());
}

void PrintGridElementNumbers(GridSubsetHandler& sh)
{
	PrintElementNumbers(sh.get_geometric_objects());
}

template <class TGeomObj>
void PrintAttachmentInfo(Grid& grid)
{
	typedef typename Grid::traits<TGeomObj>::AttachmentPipe	AttachmentPipe;
	typedef typename AttachmentPipe::ConstAttachmentEntryIterator AttIter;

//	iterate over all attachments of the grid
	AttachmentPipe& pipe = grid.get_attachment_pipe<TGeomObj>();

	int counter = 1;
	size_t totalSize = 0;
	for(AttIter iter = pipe.attachments_begin();
		iter != pipe.attachments_end(); ++iter, ++counter)
	{
	//	name
		IAttachment* att = iter->m_pAttachment;
		UG_LOG("Attachment " << counter << " (" << att->get_name() << "): ");

	//	size
		IAttachmentDataContainer* con = iter->m_pContainer;
		UG_LOG(con->occupied_memory() << " bytes\n");
		totalSize += con->occupied_memory();
	}

	UG_LOG(counter - 1 << " attachments with a total size of "
			<< totalSize << " bytes.\n");
}

void PrintAttachmentInfo(Grid& grid)
{
	UG_LOG("Vertex Attachments:\n");
	PrintAttachmentInfo<VertexBase>(grid);

	UG_LOG("\nEdge Attachments:\n");
	PrintAttachmentInfo<EdgeBase>(grid);

	UG_LOG("\nFace Attachments:\n");
	PrintAttachmentInfo<Face>(grid);

	UG_LOG("\nVolume Attachments:\n");
	PrintAttachmentInfo<Volume>(grid);
}

bool CheckHangingNodeConsistency(Grid& g)
{
	bool isConsistent = true;

//	iterate over all hanging nodes and check whether the associated parent
//	contains the node in its list of constraiend objects
	for(Grid::traits<ConstrainedVertex>::iterator iter = g.begin<ConstrainedVertex>();
		iter != g.end<ConstrainedVertex>(); ++iter)
	{
		ConstrainedVertex* hnode = *iter;
		GeometricObject* constrObj = hnode->get_constraining_object();
		if(!constrObj){
			UG_ERR_LOG("CheckHangingNodeConsistency: Hanging Vertex has no constraining object!\n");
			UG_ERR_LOG("  at: " << GetGeometricObjectCenter(g, hnode) << endl);
			isConsistent = false;
			continue;
		}

		if(ConstrainingEdge* ce = dynamic_cast<ConstrainingEdge*>(constrObj)){
		//	check whether hnode is a constrained object of ce
			if(!ce->is_constrained_object(hnode)){
				UG_ERR_LOG("CheckHangingNodeConsistency: Hanging Vertex is not constrained by parent edge!\n");
				UG_ERR_LOG("  at: " << GetGeometricObjectCenter(g, hnode) << endl);
				isConsistent = false;
				continue;
			}
		}
		else if(ConstrainingFace* cf = dynamic_cast<ConstrainingFace*>(constrObj)){
		//	check whether hnode is a constraiend object of cf
			if(!cf->is_constrained_object(hnode)){
				UG_ERR_LOG("CheckHangingNodeConsistency: Hanging Vertex is not constrained by parent face!\n");
				UG_ERR_LOG("  at: " << GetGeometricObjectCenter(g, hnode) << endl);
				isConsistent = false;
				continue;
			}
		}
		else{
			UG_ERR_LOG("CheckHangingNodeConsistency: Parent of Hanging Vertex is not a constraining object!\n");
			UG_ERR_LOG("  at: " << GetGeometricObjectCenter(g, hnode) << endl);
			isConsistent = false;
			continue;
		}
	}

//	iterate over all constrained edges and check whether the associated constraining
//	object contains the edge in its list of constraiend objects
	for(Grid::traits<ConstrainedEdge>::iterator iter = g.begin<ConstrainedEdge>();
		iter != g.end<ConstrainedEdge>(); ++iter)
	{
		ConstrainedEdge* e = *iter;
		GeometricObject* parent = e->get_constraining_object();

		if(!parent){
			UG_ERR_LOG("CheckHangingNodeConsistency: Constrained Edge has no parent!\n");
			UG_ERR_LOG("  at: " << GetGeometricObjectCenter(g, e) << endl);
			isConsistent = false;
			continue;
		}

		if(!parent->is_constraining()){
			UG_ERR_LOG("CheckHangingNodeConsistency: Parent of Constrained Edge is not a constraining object!\n");
			UG_ERR_LOG("  at: " << GetGeometricObjectCenter(g, e) << endl);
			isConsistent = false;
			continue;
		}

		if(ConstrainingEdge* ce = dynamic_cast<ConstrainingEdge*>(parent)){
		//	check whether e is a constraiend object of ce
			if(!ce->is_constrained_object(e)){
				UG_ERR_LOG("CheckHangingNodeConsistency: Constrained Edge is not constrained by parent edge!\n");
				UG_ERR_LOG("  at: " << GetGeometricObjectCenter(g, e) << endl);
				isConsistent = false;
				continue;
			}
		}
		else if(ConstrainingFace* cf = dynamic_cast<ConstrainingFace*>(parent)){
		//	check whether hnode is a constraiend object of cf
			if(!cf->is_constrained_object(e)){
				UG_ERR_LOG("CheckHangingNodeConsistency: Constrained Edge is not constrained by parent face!\n");
				UG_ERR_LOG("  at: " << GetGeometricObjectCenter(g, e) << endl);
				isConsistent = false;
				continue;
			}
		}
		else{
			UG_ERR_LOG("CheckHangingNodeConsistency: Unknown parent type of Constrained Edge!\n");
			UG_ERR_LOG("  at: " << GetGeometricObjectCenter(g, e) << endl);
			isConsistent = false;
			continue;
		}
	}

//	CHECK FACES
//	for convenience, we'll add all constraining faces to a vector, first
	vector<ConstrainingFace*> constrainingFaces;
	constrainingFaces.assign(g.begin<ConstrainingTriangle>(), g.end<ConstrainingTriangle>());
	constrainingFaces.assign(g.begin<ConstrainingQuadrilateral>(), g.end<ConstrainingQuadrilateral>());
	
	Grid::edge_traits::secure_container	edges;

	for(vector<ConstrainingFace*>::iterator iter = constrainingFaces.begin();
		iter != constrainingFaces.end(); ++iter)
	{
		ConstrainingFace* cf = *iter;
		
	//	make sure that all associated edges are constraining, too
		g.associated_elements(edges, cf);
		for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge){
			if(!edges[i_edge]->is_constraining()){
				UG_ERR_LOG("CheckHangingNodeConsistency: Non constraining side of constraining face detected!\n"
						   "  at " << GetGeometricObjectCenter(g, edges[i_edge]) << endl);
				isConsistent = false;
			}
		}
	}

	return isConsistent;
}

bool CheckHangingNodeConsistency(MultiGrid& mg)
{
	Grid& g = mg;
	bool isConsistent = CheckHangingNodeConsistency(g);

//	iterate over all edges if an edge has children, then collect its ass. faces.
//	If not all such faces have children, then the edge has to be a constraining
//	edge and its children have to be constrained.
	vector<Face*> faces;
	for(EdgeBaseIterator iter = mg.begin<EdgeBase>();
		iter != mg.end<EdgeBase>(); ++iter)
	{
		EdgeBase* e = *iter;
		if(mg.has_children<EdgeBase>(e)){
		//	e may not be a constrained edge
			if(e->is_constrained()){
				UG_ERR_LOG("CheckHangingNodeConsistency: Constrained Edge may not have children!\n");
				UG_ERR_LOG("  at: " << GetGeometricObjectCenter(g, e) << endl);
				isConsistent = false;
				continue;
			}

		//	check whether all ass. faces have children.
			bool allNbrsAreParents= true;
			bool hasConstrainingNbrFaces = false;
			CollectAssociated(faces, mg, e);
			for(size_t i = 0; i < faces.size(); ++i){
				if(!mg.has_children<Face>(faces[i]))
					allNbrsAreParents = false;
				if(faces[i]->is_constraining())
					hasConstrainingNbrFaces = true;
			}
			

		//	depending whether e is constrained or normal, either all nbrs
		//	or only some have to have children.
			if(e->is_constraining()){
				if(allNbrsAreParents && (!hasConstrainingNbrFaces)){
				//	e should be a normal edge
					UG_ERR_LOG("CheckHangingNodeConsistency: At least one neighbor face of a"
							" constraining edge should not have children!\n");
					UG_ERR_LOG("  at: " << GetGeometricObjectCenter(g, e) << endl);
					isConsistent = false;
					continue;
				}
				else{
				//	make sure that e has two constrained edge-children and a
				//	constrained vertex child.
					ConstrainingEdge* ce = dynamic_cast<ConstrainingEdge*>(e);
					UG_ASSERT(ce, "Only ConstrainingEdges should return true in EdgeBase::is_constrained()!");
					if(ce){
						if(ce->num_constrained_edges() != 2){
							UG_ERR_LOG("CheckHangingNodeConsistency: ConstrainingEdge "
									" has to constrain 2 edges not "
									<< ce->num_constrained_edges() << "!\n");
							UG_ERR_LOG("  at: " << GetGeometricObjectCenter(g, e) << endl);
							isConsistent = false;
							continue;
						}

						if(!(ce->constrained_edge(0)->is_constrained())
							 && ce->constrained_edge(1)->is_constrained())
						{
							UG_ERR_LOG("CheckHangingNodeConsistency: Both edges constrained by "
									" a constraining edge have to be ConstrainedEdges!\n");
							UG_ERR_LOG("  at: " << GetGeometricObjectCenter(g, e) << endl);
							isConsistent = false;
							continue;
						}

						if(ce->num_constrained_vertices() != 1){
							UG_ERR_LOG("CheckHangingNodeConsistency: ConstrainingEdge "
									" has to constrain 1 vertex, not "
									<< ce->num_constrained_vertices() << "!\n");
							UG_ERR_LOG("  at: " << GetGeometricObjectCenter(g, e) << endl);
							isConsistent = false;
							continue;
						}

						if(!ce->constrained_vertex(0)->is_constrained()){
							UG_ERR_LOG("CheckHangingNodeConsistency: a vertex constrained by "
									" a constraining edge has to be a hangingVertex!\n");
							UG_ERR_LOG("  at: " << GetGeometricObjectCenter(g, e) << endl);
							isConsistent = false;
							continue;
						}
					}
				}
			}
			else{
			//	e is not constraining
			//	all neighbors should be parents, too!
				if(!allNbrsAreParents){
				//	e should be a normal edge
					UG_ERR_LOG("CheckHangingNodeConsistency: All neighbor faces of a"
							" normal edge should have children!\n");
					UG_ERR_LOG("  at: " << GetGeometricObjectCenter(g, e) << endl);
					isConsistent = false;
					continue;
				}
			}
		}
	}

//	todo: add check for faces

	return isConsistent;
}

}// end of namespace
