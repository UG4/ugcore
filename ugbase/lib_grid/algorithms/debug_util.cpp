// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 13.01.2011 (m,d,y)
 
#include "debug_util.h"
#include "attachment_util.h"

#ifdef UG_PARALLEL
#include "lib_grid/parallelization/distributed_grid.h"
#include "lib_grid/parallelization/util/compol_copy_attachment.h"
#include "lib_grid/parallelization/util/compol_binary_or_attachment.h"
#include "pcl/pcl.h"
#endif

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


vector3 GetGeometricObjectCenter(Grid& g, GeometricObject* elem)
{
	switch(elem->base_object_id()){
		case VERTEX:	return GetGeometricObjectCenter(g, static_cast<VertexBase*>(elem));
		case EDGE:		return GetGeometricObjectCenter(g, static_cast<EdgeBase*>(elem));
		case FACE:		return GetGeometricObjectCenter(g, static_cast<Face*>(elem));
		case VOLUME:	return GetGeometricObjectCenter(g, static_cast<Volume*>(elem));
		default:		UG_THROW("Unknown base object type."); break;
	}
	return vector3(0, 0, 0);
}


template <class TElem>
static void CheckMultiGridConsistencyImpl(MultiGrid& mg)
{
	DistributedGridManager* dgm = mg.distributed_grid_manager();

	for(size_t lvl = 0; lvl < mg.num_levels(); ++lvl){
		for(typename MultiGrid::traits<TElem>::iterator iter = mg.begin<TElem>(lvl);
			iter != mg.end<TElem>(lvl); ++iter)
		{
			TElem* e = *iter;
		//	make sure that all children have the local element as parent
			for(size_t i_child = 0; i_child < mg.num_children<VertexBase>(e); ++i_child)
			{
				if(e != mg.get_parent(mg.get_child<VertexBase>(e, i_child))){
					UG_THROW("parent is not referenced by children!");
				}
			}

			for(size_t i_child = 0; i_child < mg.num_children<EdgeBase>(e); ++i_child)
			{
				if(e != mg.get_parent(mg.get_child<EdgeBase>(e, i_child))){
					UG_THROW("parent is not referenced by children!");
				}
			}

			for(size_t i_child = 0; i_child < mg.num_children<Face>(e); ++i_child)
			{
				if(e != mg.get_parent(mg.get_child<Face>(e, i_child))){
					UG_THROW("parent is not referenced by children!");
				}
			}

			for(size_t i_child = 0; i_child < mg.num_children<Volume>(e); ++i_child)
			{
				if(e != mg.get_parent(mg.get_child<Volume>(e, i_child))){
					UG_THROW("parent is not referenced by children!");
				}
			}

		//	also make sure that each element with a parent is contained in the
		//	children list of its parent
			GeometricObject* parent = mg.get_parent(e);
			if(parent){
				bool gotIt = false;
				for(size_t i = 0; i < mg.num_children<TElem>(parent); ++i){
					if(mg.get_child<TElem>(parent, i) == e){
						gotIt = true;
						break;
					}
				}

				if(!gotIt){
					UG_THROW("Element not contained in child list of its parent");
				}
			}
			else{
				if(lvl > 0){
					#ifdef UG_PARALLEL
						if(!dgm->contains_status(e, INT_V_SLAVE))
						{
							UG_THROW("Each element in a higher level has to have a parent"
									 " unless it is a v-slave!");
						}
					#else
						UG_THROW("Each element in a higher level has to have a parent"
								 " unless it is a v-slave!");
					#endif
				}
			}
		}
	}
}

void CheckMultiGridConsistency(MultiGrid& mg)
{
	CheckMultiGridConsistencyImpl<VertexBase>(mg);
	CheckMultiGridConsistencyImpl<EdgeBase>(mg);
	CheckMultiGridConsistencyImpl<Face>(mg);
	CheckMultiGridConsistencyImpl<Volume>(mg);
}

//	the following define is only used by the CheckHangingNodeConsistency methods
//	and is highly specialized for them.
#define FAILED_CHECK(elem, msg)\
	{UG_ERR_LOG("CheckHangingNodeConsistency: " << msg << endl);\
	 UG_ERR_LOG("  at: " << GetGeometricObjectCenter(g, elem));\
	 if(MultiGrid* mg = dynamic_cast<MultiGrid*>(&g))\
		{UG_ERR_LOG(" on level " << mg->get_level(elem));}\
	 UG_ERR_LOG(endl);\
	 isConsistent = false;\
	}//continue;}


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

//	iterate over all hanging nodes and check whether the associated parent
//	matches the constraining object. Other checks have already been performed!
	for(Grid::traits<ConstrainedVertex>::iterator iter = g.begin<ConstrainedVertex>();
		iter != g.end<ConstrainedVertex>(); ++iter)
	{
		ConstrainedVertex* hnode = *iter;
		GeometricObject* co = hnode->get_constraining_object();
		GeometricObject* parent = mg.get_parent(hnode);

		if(co != parent){
			FAILED_CHECK(hnode, "Hanging Vertex: Constraining object and parent do not match!");
		}

		if(co && ((mg.get_level(co) + 1) != mg.get_level(hnode))){
			FAILED_CHECK(hnode, "Hanging Vertex has to be one level higher than its constraining object!");
		}
	}

//	iterate over all edges if an edge has children, then collect its ass. faces.
//	If not all such faces have children, then the edge has to be a constraining
//	edge and its children have to be constrained.
	vector<Face*> faces;
	for(EdgeBaseIterator iter = mg.begin<EdgeBase>();
		iter != mg.end<EdgeBase>(); ++iter)
	{
		EdgeBase* e = *iter;

		if(mg.get_parent(e) && (mg.get_level(mg.get_parent(e)) + 1 != mg.get_level(e))){
			FAILED_CHECK(e, "Edge and parent are not in consecutive levels!");
		}

		for(size_t i = 0; i < 2; ++i){
			if(mg.get_level(e) != mg.get_level(e->vertex(i))){
				FAILED_CHECK(e, "Edge-Vertex is not on the same level as the edge itself!");
				UG_ERR_LOG("  Vertex at " << GetGeometricObjectCenter(mg, e->vertex(i))
							<< " on level " << mg.get_level(e->vertex(i)) << endl);
			}
		}

		if(mg.has_children<EdgeBase>(e)){
		//	e may not be a constrained edge
			if(e->is_constrained()){
				FAILED_CHECK(e, "Constrained Edge may not have children!");
			}
/*
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
*/

		//	depending whether e is constrained or normal, either all nbrs
		//	or only some have to have children.
			if(e->is_constraining()){
				/* This causes false positives for distributed grids!
				if(allNbrsAreParents && (!hasConstrainingNbrFaces)){
				//	e should be a normal edge
					FAILED_CHECK(e, "At least one neighbor face of a"
								 " constraining edge should not have children!");
				}
				else{*/
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
				//}
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


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
enum ConstraintTypes{
	CT_NONE = 0,
	CT_CONSTRAINING = 1,
	CT_CONSTRAINED = 1 << 1
};

template <class TElem>
static bool CheckDistributedObjectConstraintTypes(MultiGrid& mg)
{
	typedef typename Grid::traits<TElem>::iterator ElemIter;
	bool retVal = true;

	AInt aState;
	mg.attach_to<TElem>(aState);
	Grid::AttachmentAccessor<TElem, AInt> aaState(mg, aState);

//	set up local states
	for(ElemIter iter = mg.begin<TElem>(); iter != mg.end<TElem>(); ++iter){
		TElem* e = *iter;
		aaState[e] = CT_NONE;
		if(e->is_constraining())
			aaState[e] |= CT_CONSTRAINING;
		if(e->is_constrained())
			aaState[e] |= CT_CONSTRAINED;
	}

//	communicate states
	#ifdef UG_PARALLEL
		typedef typename GridLayoutMap::Types<TElem>::Layout Layout;
		pcl::InterfaceCommunicator<Layout> com;
		DistributedGridManager& distGridMgr = *mg.distributed_grid_manager();
		GridLayoutMap& glm = distGridMgr.grid_layout_map();

		ComPol_BinaryOrAttachment<Layout, AInt> compolOr(mg, aState);
		com.exchange_data(glm, INT_H_SLAVE, INT_H_MASTER, compolOr);
		com.communicate();

		com.exchange_data(glm, INT_H_MASTER, INT_H_SLAVE, compolOr);
		com.communicate();

		com.exchange_data(glm, INT_V_SLAVE, INT_V_MASTER, compolOr);
		com.communicate();

		com.exchange_data(glm, INT_V_MASTER, INT_V_SLAVE, compolOr);
		com.communicate();

		com.exchange_data(glm, INT_V_SLAVE, INT_V_MASTER, compolOr);
		com.communicate();
	#endif

//	check whether communicated states match the actual element states on this proc.
	for(ElemIter iter = mg.begin<TElem>(); iter != mg.end<TElem>(); ++iter){
		TElem* e = *iter;
		int state = CT_NONE;
		if(e->is_constraining())
			state |= CT_CONSTRAINING;
		if(e->is_constrained())
			state |= CT_CONSTRAINED;

		if(state != aaState[e]){
			UG_LOG("ERROR: Distributed object has different constraint states on different procs. "
					<< "At: " << GetGeometricObjectCenter(mg, e)
					<< " on level " << mg.get_level(e) << endl);
			retVal = false;
		}
	}

	mg.detach_from<TElem>(aState);
	return retVal;
}

bool CheckDistributedObjectConstraintTypes(MultiGrid& mg)
{
	bool retVal = true;

//	assign constraint states
	UG_LOG("Checking constraint types of VERTICES\n");
	retVal &= CheckDistributedObjectConstraintTypes<VertexBase>(mg);
	UG_LOG("Checking constraint types of EDGES\n");
	retVal &= CheckDistributedObjectConstraintTypes<EdgeBase>(mg);
	UG_LOG("Checking constraint types of FACES\n");
	retVal &= CheckDistributedObjectConstraintTypes<Face>(mg);
	UG_LOG("Checking constraint types DONE\n");

//	make sure that we return the same value globally
	#ifdef UG_PARALLEL
		retVal = pcl::AllProcsTrue(retVal);
	#endif

	return retVal;
}

#ifdef UG_PARALLEL
template <class TLayout>
class ComPol_CheckDistributedParentStates : public pcl::ICommunicationPolicy<TLayout>
{
	public:
		typedef TLayout								Layout;
		typedef typename Layout::Type				GeomObj;
		typedef typename Layout::Element			Element;
		typedef typename Layout::Interface			Interface;
		typedef typename Interface::const_iterator	InterfaceIter;

		ComPol_CheckDistributedParentStates(MultiGrid& mg) :
			m_mg(mg),
			m_dgm(mg.distributed_grid_manager()),
			m_comparisionFailed(false),
			m_performMasterCheck(false)
		{
		}

		virtual ~ComPol_CheckDistributedParentStates()	{}

		virtual int
		get_required_buffer_size(const Interface& interface)
		{return interface.size() * sizeof(int);}

		virtual bool
		collect(ug::BinaryBuffer& buff, const Interface& interface)
		{
			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				Element elem = interface.get_element(iter);
				int parentType = m_mg.parent_type(elem);
				Serialize(buff, parentType);
			}
			return true;
		}

		virtual bool
		extract(ug::BinaryBuffer& buff, const Interface& interface)
		{
			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				Element elem = interface.get_element(iter);
				int parentType;
				Deserialize(buff, parentType);

				if(m_performMasterCheck){
					GeometricObject* parent = m_mg.get_parent(elem);
					if(parent && (parentType != parent->base_object_id())){
						UG_LOG("  PARENT-TYPE MISMATCH AT CHILD ELEMENT WITH OBJECT ID " << elem->base_object_id()
								<< " at " << GetGeometricObjectCenter(m_mg, elem) << " on level " << m_mg.get_level(elem) << "\n");
						UG_LOG("    Parent object id: " << parent->base_object_id() << ", received id: " << parentType << "\n");
						m_comparisionFailed = true;
					}
					else if((!parent) && (parentType != -1)){
						UG_LOG("  PARENT-TYPE MISMATCH AT CHILD ELEMENT WITH OBJECT ID " << elem->base_object_id()
								<< " at " << GetGeometricObjectCenter(m_mg, elem) << " on level " << m_mg.get_level(elem) << "\n");
						UG_LOG("  The element hasn't got a parent but received parent id is != -1. Received id: " << parentType << "\n");
						m_comparisionFailed = true;
					}
				}
				else if(parentType != m_mg.parent_type(elem)){
					UG_LOG("  PARENT-TYPE MISMATCH AT ELEMENT WITH OBJECT ID " << elem->base_object_id()
							<< " at " << GetGeometricObjectCenter(m_mg, elem) << " on level " << m_mg.get_level(elem) << "\n");
					UG_LOG("    Parent object id: " <<  (int)m_mg.parent_type(elem) << ", received id: " << parentType << "\n");
					m_comparisionFailed = true;
				}
			}
			return true;
		}

		bool exchange_data()
		{
			pcl::InterfaceCommunicator<TLayout> com;

			m_comparisionFailed = false;
			GridLayoutMap& glm = m_dgm->grid_layout_map();
			m_performMasterCheck = false;
			if(glm.has_layout<GeomObj>(INT_H_SLAVE))
				com.send_data(glm.get_layout<GeomObj>(INT_H_SLAVE), *this);
			if(glm.has_layout<GeomObj>(INT_H_MASTER))
				com.receive_data(glm.get_layout<GeomObj>(INT_H_MASTER), *this);
			com.communicate();

			m_performMasterCheck = true;
			if(glm.has_layout<GeomObj>(INT_V_SLAVE))
				com.send_data(glm.get_layout<GeomObj>(INT_V_SLAVE), *this);
			if(glm.has_layout<GeomObj>(INT_V_MASTER))
				com.receive_data(glm.get_layout<GeomObj>(INT_V_MASTER), *this);
			com.communicate();

			return !m_comparisionFailed;
		}

	private:
		MultiGrid& m_mg;
		DistributedGridManager* m_dgm;
		bool m_comparisionFailed;
		bool m_performMasterCheck;
};
#endif

template <class TElem>
bool CheckLocalParentTypes(MultiGrid& mg)
{
	bool success = true;
	typedef typename MultiGrid::traits<TElem>::iterator iter_t;
	for(iter_t iter = mg.begin<TElem>(); iter != mg.end<TElem>(); ++iter){
		TElem* e = *iter;
		GeometricObject* parent = mg.get_parent(e);
		if(parent){
			if(mg.parent_type(e) != parent->base_object_id()){
				UG_LOG("  LOCAL PARENT-TYPE MISMATCH AT ELEMENT WITH OBJECT ID " << e->base_object_id()
						<< " at " << GetGeometricObjectCenter(mg, e) << " on level " << mg.get_level(e) << "\n");
				UG_LOG("    Stored parent id: " <<  (int)mg.parent_type(e) << ", actual parent id: " << parent->base_object_id() << "\n");
				success = false;
			}
		}
	}
	return success;
}

bool CheckDistributedParentTypes(MultiGrid& mg)
{
	UG_LOG("DEBUG: Checking distributed parent types...\n");
#ifdef UG_PARALLEL
	ComPol_CheckDistributedParentStates<VertexLayout>	vrtChecker(mg);
	ComPol_CheckDistributedParentStates<EdgeLayout>		edgeChecker(mg);
	ComPol_CheckDistributedParentStates<FaceLayout> 	faceChecker(mg);
	ComPol_CheckDistributedParentStates<VolumeLayout>	volChecker(mg);

	bool success = true;
	UG_LOG(" checking vertices...\n");
	success &= CheckLocalParentTypes<VertexBase>(mg);
	success &= vrtChecker.exchange_data();
	UG_LOG(" checking edges...\n");
	success &= CheckLocalParentTypes<EdgeBase>(mg);
	success &= edgeChecker.exchange_data();
	UG_LOG(" checking faces...\n");
	success &= CheckLocalParentTypes<Face>(mg);
	success &= faceChecker.exchange_data();
	UG_LOG(" checking volumes...\n");
	success &= CheckLocalParentTypes<Volume>(mg);
	success &= volChecker.exchange_data();
	UG_LOG(" checking done with status ");

	success = pcl::AllProcsTrue(success);

	if(success){
		UG_LOG("SUCCESS\n");
	}
	else{
		UG_LOG("FAILURE\n");
	}
	return success;
#else
	return true;
#endif
}


bool CheckElementConsistency(MultiGrid& mg, VertexBase* v)
{
	bool success = true;
	UG_LOG("DEBUG: Checking vertex at " << GetGeometricObjectCenter(mg, v) << endl);
	UG_LOG("  vertex type:");
	if(v->is_constrained()){
		UG_LOG(" constrained");
	}
	else{
		UG_LOG(" normal");
	}
	UG_LOG("\n");

	if(v->is_constrained()){
		ConstrainedVertex* cdv = dynamic_cast<ConstrainedVertex*>(v);
		(void) cdv; // removes unused warning
		UG_ASSERT(cdv, "Bad type!");
		UG_ASSERT(cdv->get_constraining_object() == mg.get_parent(cdv),
				  "ConstrainingObject / Parent mismatch");
	}
	else{
	}

	return success;
}

bool CheckElementConsistency(MultiGrid& mg, EdgeBase* e)
{
	bool success = true;
	UG_LOG("DEBUG: Checking edge at " << GetGeometricObjectCenter(mg, e) << endl);
	UG_LOG("  edge type:");
	if(e->is_constrained()){
		UG_LOG(" constrained");
	}
	else if(e->is_constraining()){
		UG_LOG(" constraining");
	}
	else{
		UG_LOG(" normal");
	}
	UG_LOG("\n");

//	constrained/constraining checks
	if(e->is_constrained()){
		ConstrainedEdge* cde = dynamic_cast<ConstrainedEdge*>(e);
		(void) cde; // removes unused warning
		UG_ASSERT(cde, "Bad type!");
		UG_ASSERT(cde->get_constraining_object() == mg.get_parent(cde),
				  "ConstrainingObject / Parent mismatch");

		UG_ASSERT(cde->vertex(0)->is_constrained() || cde->vertex(1)->is_constrained(),
				  "At least one corner of a constrained edge has to be a constrained vertex");
	}
	else if(e->is_constraining()){
		ConstrainingEdge* cge = dynamic_cast<ConstrainingEdge*>(e);
		UG_ASSERT(cge, "Bad type!");
		UG_ASSERT(cge->constrained_vertex(0), "A constrained vertex has to exist");
		UG_ASSERT(mg.get_child_vertex(cge) == cge->constrained_vertex(0),
				  "Mismatch between vertex child and constrained vertex.");
		UG_ASSERT(cge->num_constrained_edges() == 2, "2 constrained edges have to exist!");

		for(size_t i1 = 0; i1 < cge->num_constrained_edges(); ++i1){
			bool constrainedEdgeMatch = false;
			for(size_t i2 = 0; i2 < mg.num_children<EdgeBase>(cge); ++i2){
				if(mg.get_child<EdgeBase>(cge, i2) == cge->constrained_edge(i1)){
					constrainedEdgeMatch = true;
					break;
				}
			}
			UG_ASSERT(constrainedEdgeMatch, "no matching child found to constrained edge");
		}

		CheckElementConsistency(mg, cge->constrained_vertex(0));
	}

//	check vertices
	for(size_t i = 0; i < e->num_vertices(); ++i){
		success &= CheckElementConsistency(mg, e->vertex(i));
	}

	return success;
}

bool CheckElementConsistency(MultiGrid& mg, Face* f)
{
	bool success = true;
	UG_LOG("DEBUG: Checking face at " << GetGeometricObjectCenter(mg, f) << endl);

	UG_LOG("  face type:");
	if(f->is_constrained()){
		UG_LOG(" constrained");
	}
	else if(f->is_constraining()){
		UG_LOG(" constraining");
	}
	else{
		UG_LOG(" normal");
	}
	UG_LOG("\n");

//	check sides
	Grid::edge_traits::secure_container edges;
	mg.associated_elements(edges, f);

	for(size_t i = 0; i < edges.size(); ++i){
		success &= CheckElementConsistency(mg, edges[i]);
	}

	return success;
}

}// end of namespace
