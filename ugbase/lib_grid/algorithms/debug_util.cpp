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



//	the following define is only used by the CheckHangingNodeConsistency methods
//	and is highly specialized for them.
#define FAILED_CHECK(elem, msg)	{UG_ERR_LOG("CheckHangingNodeConsistency: " << msg << endl);\
								UG_ERR_LOG("  at: " << GetGeometricObjectCenter(g, elem) << endl);\
								isConsistent = false;\
								continue;}


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
			FAILED_CHECK(hnode, "Hanging Vertex has no constraining object!");
		}

		if(ConstrainingEdge* ce = dynamic_cast<ConstrainingEdge*>(constrObj)){
		//	check whether hnode is a constrained object of ce
			if(!ce->is_constrained_object(hnode)){
				FAILED_CHECK(hnode, "Hanging Vertex is not constrained by parent edge!");
			}
		}
		else if(ConstrainingFace* cf = dynamic_cast<ConstrainingFace*>(constrObj)){
		//	check whether hnode is a constraiend object of cf
			if(!cf->is_constrained_object(hnode)){
				FAILED_CHECK(hnode, "Hanging Vertex is not constrained by parent face!");
			}
		}
		else{
			FAILED_CHECK(hnode, "Parent of Hanging Vertex is not a constraining object!");
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
			FAILED_CHECK(e, "Constrained Edge has no parent!");
		}

		if(!parent->is_constraining()){
			FAILED_CHECK(e, "Parent of Constrained Edge is not a constraining object!");
		}

		if(ConstrainingEdge* ce = dynamic_cast<ConstrainingEdge*>(parent)){
		//	check whether e is a constraiend object of ce
			if(!ce->is_constrained_object(e)){
				FAILED_CHECK(e, "Constrained Edge is not constrained by parent edge!");
			}
			
		//	since the constraining object is an edge, we'll make sure, that the constrained edge
		//	is connected to exactly one constrained vertex.
			if(e->vertex(0)->is_constrained() && e->vertex(1)->is_constrained()){
				FAILED_CHECK(e, "Constrained Edge (constrained by edge) "
							 " may not be connected to two constrained vertices!");
			}
		}
		else if(ConstrainingFace* cf = dynamic_cast<ConstrainingFace*>(parent)){
		//	check whether hnode is a constraiend object of cf
			if(!cf->is_constrained_object(e)){
				FAILED_CHECK(e, "Constrained Edge is not constrained by parent face!");
			}
		}
		else{
			FAILED_CHECK(e, "Unknown parent type of Constrained Edge!");
		}
		
	//	make sure that the edge connected to at least one constrained vertex
	//	we check above, that an edge which is constrained by an edge is not
	//	connected to two constrained vertices.
		if(!(e->vertex(0)->is_constrained() || e->vertex(1)->is_constrained())){
			FAILED_CHECK(e, "Constrained Edge has to be connected to a constrained vertex!");
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
				FAILED_CHECK(edges[i_edge], "Non constraining side of constraining face detected!");
			}
		}
	}

	return isConsistent;
}

bool CheckHangingNodeConsistency(MultiGrid& mg)
{
	Grid& g = mg;
	bool isConsistent = CheckHangingNodeConsistency(g);
UG_LOG("1\n");
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
				FAILED_CHECK(e, "Constrained Edge may not have children!");
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
					FAILED_CHECK(e, "At least one neighbor face of a"
								 " constraining edge should not have children!");
				}
				else{
				//	make sure that e has two constrained edge-children and a
				//	constrained vertex child.
					ConstrainingEdge* ce = dynamic_cast<ConstrainingEdge*>(e);
					UG_ASSERT(ce, "Only ConstrainingEdges should return true in EdgeBase::is_constrained()!");
					if(ce){
						if(ce->num_constrained_edges() != 2){
							FAILED_CHECK(e, "ConstrainingEdge has to constrain 2 edges not "
										 << ce->num_constrained_edges() << "!");
						}

						if(!(ce->constrained_edge(0)->is_constrained())
							 && ce->constrained_edge(1)->is_constrained())
						{
							FAILED_CHECK(e, "Both edges constrained by "
										 " a constraining edge have to be ConstrainedEdges!");
						}

						if(ce->num_constrained_vertices() != 1){
							FAILED_CHECK(e, "ConstrainingEdge has to constrain 1 vertex, not "
										<< ce->num_constrained_vertices() << "!");
						}

						if(!ce->constrained_vertex(0)->is_constrained()){
							FAILED_CHECK(e, "a vertex constrained by a constraining "
										 "edge has to be a hangingVertex!");
						}
					}
				}
			}
			else{
			//	e is not constraining
			//	all neighbors should be parents, too!
				/*this check is only valid if we don't create a closure
				if(!allNbrsAreParents){
				//	e should be a normal edge
					FAILED_CHECK(e, "All neighbor faces of a normal edge should have children!");
				}*/
			}
		}
	}
UG_LOG("2\n");
//	check faces
	for(FaceIterator iter = g.begin<Face>(); iter != g.end<Face>(); ++iter){
		Face* f = *iter;
		
		if(f->is_constraining()){
			UG_LOG("2.1/n");
			ConstrainingFace* cf = dynamic_cast<ConstrainingFace*>(f);
			UG_ASSERT(cf, "All constraining faces should derive from ConstrainingFace");
			
		//	make sure that the face has the right number of children
			if(mg.num_children<Face>(f) != 4){
				FAILED_CHECK(f, "Face has bad number of child faces. "
						     "4 required, " << mg.num_children<Face>(f) << " found.");
			}
			
			if(mg.num_children<EdgeBase>(f) != f->num_vertices()){
				FAILED_CHECK(f, "Face has bad number of child edges. "
						     << f->num_vertices() << " required, "
						     << mg.num_children<EdgeBase>(f) << " found.");
			}
			
			if((f->num_vertices() > 3) && mg.num_children<VertexBase>(f) != 1){
				FAILED_CHECK(f, "Face has bad number of child vertices. "
							 "1 required, " << mg.num_children<VertexBase>(f) << " found.");
			}
			
		//	make sure that number of children and number of constrained elements match
			if(mg.num_children<Face>(f) != cf->num_constrained_faces())
				FAILED_CHECK(f, "Number of child faces of constraining face does not match number of constrained faces.");

			if(mg.num_children<EdgeBase>(f) != cf->num_constrained_edges())
				FAILED_CHECK(f, "Number of child edges of constraining face does not match number of constrained edges.");

			if(mg.num_children<VertexBase>(f) != cf->num_constrained_vertices())
				FAILED_CHECK(f, "Number of child vertices of constraining face does not match number of constrained vertices.");
			
		//	make sure that all children are constrained and that they are contained
		//	in the list of constrained elements
			for(size_t i = 0; i < mg.num_children<Face>(f); ++i){
				Face* child = mg.get_child<Face>(f, i);
				if(!child->is_constrained()){
					FAILED_CHECK(f, "All child faces of a constraining face have to be constrained faces.");
				}
				else{
				//	make sure that the child is contained in the list of constrained objects of cf
					if(!cf->is_constrained_object(child))
						FAILED_CHECK(f, "Child face of constraining face is not in list of constrained faces.");
				}
			}
			for(size_t i = 0; i < mg.num_children<EdgeBase>(f); ++i){
				EdgeBase* child = mg.get_child<EdgeBase>(f, i);
				if(!child->is_constrained()){
					FAILED_CHECK(f, "All child edges of a constraining face have to be constrained edges.");
				}
				else{
				//	make sure that the child is contained in the list of constrained objects of cf
					if(!cf->is_constrained_object(child))
						FAILED_CHECK(f, "Child edge of constraining face is not in list of constrained edges.");
				}
			}
			for(size_t i = 0; i < mg.num_children<VertexBase>(f); ++i){
				VertexBase* child = mg.get_child<VertexBase>(f, i);
				if(!child->is_constrained()){
					FAILED_CHECK(f, "All child vertices of a constraining face have to be constrained vertices.");
				}
				else{
				//	make sure that the child is contained in the list of constrained objects of cf
					if(!cf->is_constrained_object(child))
						FAILED_CHECK(f, "Child vertex of constraining face is not in list of constrained vertices.");
				}
			}
		}
		
		else if(f->is_constrained()){
			UG_LOG("2.2.0/n");
			ConstrainedFace* cdf = dynamic_cast<ConstrainedFace*>(f);
			UG_ASSERT(cdf, "All constrained faces should derive from ConstrainedFace");
			UG_LOG("2.2.1/n");
		//	we don't have to check all interconnections, since we already checked a lot of
		//	stuff for constraining faces. So just do the rest now.
			ConstrainingFace* cf = dynamic_cast<ConstrainingFace*>(cdf->get_constraining_object());
			UG_LOG("2.2.2/n");
			if(!cf){
				FAILED_CHECK(cdf, "No constraining face found for given constrained face."); 
			}
			
			if(cf != mg.get_parent(cdf)){
				FAILED_CHECK(cdf, "The parent of a constrained face should always be its constraining face!");
			}
		}
	}
UG_LOG("3\n");
	return isConsistent;
}

}// end of namespace

