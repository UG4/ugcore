// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y11 m01 d11

#include <vector>
#include "hanging_node_refiner_base.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"

using namespace std;

namespace ug{

HangingNodeRefinerBase::
HangingNodeRefinerBase(IRefinementCallback* refCallback) :
	m_refCallback(refCallback),
	m_pGrid(NULL)
{
}

HangingNodeRefinerBase::
HangingNodeRefinerBase(const HangingNodeRefinerBase&)
{
	throw(UGError("no copy construction allowed."));
}

HangingNodeRefinerBase::
~HangingNodeRefinerBase()
{
	if(m_pGrid)
	{
		m_pGrid->unregister_observer(this);
	}
}

void HangingNodeRefinerBase::
set_refinement_callback(IRefinementCallback* refCallback)
{
	m_refCallback = refCallback;
}

void HangingNodeRefinerBase::
set_grid(Grid* grid)
{
	if(m_pGrid)
	{
		m_pGrid->unregister_observer(this);
		m_selMarkedElements.assign_grid(NULL);
		m_pGrid = NULL;
	}

	if(grid){
		m_pGrid = grid;
		grid->register_observer(this, OT_GRID_OBSERVER);
		m_selMarkedElements.assign_grid(*grid);
		m_selMarkedElements.enable_autoselection(false);
		m_selMarkedElements.enable_selection_inheritance(false);
	}
}

void HangingNodeRefinerBase::grid_to_be_destroyed(Grid* grid)
{
	m_pGrid = NULL;
}

void HangingNodeRefinerBase::clear_marks()
{
	m_selMarkedElements.clear();
}

void HangingNodeRefinerBase::mark_for_refinement(EdgeBase* e)
{
	assert(m_pGrid && "ERROR in HangingNodeRefinerBase::mark_for_refinement(...): No grid assigned.");

	mark(e);
}

void HangingNodeRefinerBase::mark_for_refinement(Face* f)
{
	assert(m_pGrid && "ERROR in HangingNodeRefinerBase::mark_for_refinement(...): No grid assigned.");

	mark(f);
}

void HangingNodeRefinerBase::mark_for_refinement(Volume* v)
{
	assert(m_pGrid && "ERROR in HangingNodeRefinerBase::mark_for_refinement(...): No grid assigned.");

	mark(v);
}

void HangingNodeRefinerBase::refine()
{
	if(!m_pGrid)
		throw(UGError("ERROR in HangingNodeRefinerBase::refine(...): No grid assigned."));

	if(m_selMarkedElements.get_assigned_grid() != m_pGrid)
		throw(UGError("selector not initialized properly. Use HangingNodeRefinerBase::set_grid."));

	Grid& grid = *m_pGrid;

//	check if a refinement-callback is set.
//	if not, we'll automatically set one, if a position attachment is available
	bool localRefCallbackSet = false;
	if(!m_refCallback){
		if(grid.has_vertex_attachment(aPosition)){
			localRefCallbackSet = true;
			m_refCallback = new RefinementCallbackLinear<APosition>(grid, aPosition);
		}
		else if(grid.has_vertex_attachment(aPosition2)){
			localRefCallbackSet = true;
			m_refCallback = new RefinementCallbackLinear<APosition2>(grid, aPosition2);
		}
	}

//	check grid options.
	if(!grid.option_is_enabled(GRIDOPT_FULL_INTERCONNECTION))
	{
		LOG("WARNING in HangingNodeRefiner_IR1::refine(): grid option GRIDOPT_FULL_INTERCONNECTION auto-enabled." << endl);
		grid.enable_options(GRIDOPT_FULL_INTERCONNECTION);
	}

//	containers used for temporary results
	vector<EdgeBase*> 	vEdges;
	vector<Face*>	 	vFaces;
	vector<Volume*>		vVols;

	UG_DLOG(LIB_GRID, 1, "performing hanging-node-refine:\n");
//	fills the queues with the elements that have to be refined.
	collect_objects_for_refine();

//	call pre_refine to allow derived classes to perform some actions
	pre_refine();

////////////////////////////////
//	ConstrainedEdges
//	constrained edges can be ignored, since they will be refined together with
//	their constraining edges.

	UG_DLOG(LIB_GRID, 1, "  constraining edges.\n");
////////////////////////////////
//	ConstrainingEdges
//	iterate through scheduled cg-edges
	{
		ConstrainingEdgeIterator iter = m_selMarkedElements.begin<ConstrainingEdge>();
		while(iter != m_selMarkedElements.end<ConstrainingEdge>()){
			ConstrainingEdge* cge = *iter;
			++iter;

			refine_constraining_edge(cge);
		}
	}

	UG_DLOG(LIB_GRID, 1, "  normal edges.\n");
////////////////////////////////
//	Normal Edges
//	iterate through scheduled edges
	{
		EdgeIterator iter = m_selMarkedElements.begin<Edge>();
		while(iter != m_selMarkedElements.end<Edge>())
		{
			Edge* e = *iter;
			++iter;

		//	check whether all connected elements are marked.
		//	If so, we can perform a normal refine.
		//	If not, we'll create a new hanging vertex and new constrained edges.
		//	The edge itself will be transformed to a constraining edge.
			bool bAllMarked = true;

			if(grid.num_faces() > 0){
				CollectFaces(vFaces, grid, e);
				for(size_t i = 0; i < vFaces.size(); ++i){
					if(!is_marked(vFaces[i])){
						bAllMarked = false;
						break;
					}
				}
			}

		//	volumes only have to be checked if all faces were marked.
			if(bAllMarked && (grid.num_volumes() > 0)){
				CollectVolumes(vVols, grid, e);
				for(size_t i = 0; i < vVols.size(); ++i){
					if(!is_marked(vVols[i])){
						bAllMarked = false;
						break;
					}
				}
			}

			if(bAllMarked)
				refine_edge_with_normal_vertex(e);
			else
				refine_edge_with_hanging_vertex(e);
		}
	}


////////////////////////////////
//	Faces
//	constrained faces are ignored, since they are automatically refined during
//	constraining-face-refine.

	UG_DLOG(LIB_GRID, 1, "  constraining triangles.\n");
//	constraining triangles
	{
		ConstrainingTriangleIterator iter = m_selMarkedElements.begin<ConstrainingTriangle>();
		while(iter != m_selMarkedElements.end<ConstrainingTriangle>()){
			ConstrainingTriangle* cgf = *iter;
			++iter;
			refine_constraining_face(cgf);
		}
	}

	UG_DLOG(LIB_GRID, 1, "  constraining quadrilaterals.\n");
//	constraining quadrilaterals
	{
		ConstrainingQuadrilateralIterator iter = m_selMarkedElements.begin<ConstrainingQuadrilateral>();
		while(iter != m_selMarkedElements.end<ConstrainingQuadrilateral>()){
			ConstrainingQuadrilateral* cgf = *iter;
			++iter;
			refine_constraining_face(cgf);
		}
	}

	UG_DLOG(LIB_GRID, 1, "  normal triangles.\n");
//	normal triangles
	{
		TriangleIterator iter = m_selMarkedElements.begin<Triangle>();
		while(iter != m_selMarkedElements.end<Triangle>()){
			Face* f = *iter;
			++iter;
		//	check whether all associated volumes are marked
			bool bAllMarked = true;
			if(grid.num_volumes() > 0){
				CollectVolumes(vVols, grid, f);
				for(size_t i = 0; i < vVols.size(); ++i){
					if(!is_marked(vVols[i])){
						bAllMarked = false;
						break;
					}
				}
			}

		//	if all are marked, the face can be normally refined.
		//	If not, we'll create constrained faces.
			if(bAllMarked)
				refine_face_with_normal_vertex(f);
			else{
				refine_face_with_hanging_vertex(f);
			}
		}
	}

	UG_DLOG(LIB_GRID, 1, "  normal quadrilaterals.\n");
//	normal quadrilaterals
	{
		QuadrilateralIterator iter = m_selMarkedElements.begin<Quadrilateral>();
		while(iter != m_selMarkedElements.end<Quadrilateral>()){
			Face* f = *iter;
			++iter;
		//	check whether all associated volumes are marked
			bool bAllMarked = true;
			if(grid.num_volumes() > 0){
				CollectVolumes(vVols, grid, f);
				for(size_t i = 0; i < vVols.size(); ++i){
					if(!is_marked(vVols[i])){
						bAllMarked = false;
						break;
					}
				}
			}

		//	if all are marked, the face can be normally refined.
		//	If not, we'll create constrained faces.
			if(bAllMarked)
				refine_face_with_normal_vertex(f);
			else{
				refine_face_with_hanging_vertex(f);
			}
		}
	}


////////////////////////////////
//	Volumes
	UG_DLOG(LIB_GRID, 1, "  volumes.\n");
	{
		VolumeIterator iter = m_selMarkedElements.begin<Volume>();
		while(iter != m_selMarkedElements.end<Volume>()){
			Volume* vol = *iter;
			++iter;
			refine_volume_with_normal_vertex(vol);
		}
	}

	UG_DLOG(LIB_GRID, 1, "  refinement done.\n");


////////////////////////////////
//	call post_refine to allow derived classes to perform some actions
	post_refine();


////////////////////////////////
//	Clean up
//	clear the refinement-callback if necessary
	if(localRefCallbackSet){
		delete m_refCallback;
		m_refCallback = NULL;
	}

	UG_DLOG(LIB_GRID, 1, "  done.\n");
}

void HangingNodeRefinerBase::collect_objects_for_refine()
{
//	Algorithm Layout:
//		The marks have to be adjusted. All edges on which a vertex has to be generated
//		should be marked.

//	This variable determines whether a marked edge leads to the refinement of
//	associated faces and volumes.
	static const bool automarkHigherDimensionalObjects = false;

//	comfortable grid access.
	Grid& grid = *m_pGrid;

//	containers used for temporary results
	vector<EdgeBase*> vEdges;
	vector<Face*> vFaces;
	vector<Volume*> vVols;

//	queues will be used to avoid recursion during element selection
	queue<EdgeBase*>	qEdges;
	queue<Face*> 		qFaces;
	queue<Volume*> 		qVols;

//	assert correct selection. see HangingVertexRefiner description.
//	collect all unmarked edges, faces and volumes adjacent to marked elements.
//	note that the queues may possibly contain several elements multiple times.
	collect_associated_unmarked_edges(qEdges, grid,
						m_selMarkedElements.begin<Face>(),
						m_selMarkedElements.end<Face>());

	collect_associated_unmarked_edges(qEdges, grid,
						m_selMarkedElements.begin<Volume>(),
						m_selMarkedElements.end<Volume>());

	if(grid.num_faces() > 0){
		if(automarkHigherDimensionalObjects){
		//	automarking of higher dimensional elements
		//	if a edge is selected, all associated faces will be refined.
			collect_associated_unmarked_faces(qFaces, grid,
								m_selMarkedElements.begin<EdgeBase>(),
								m_selMarkedElements.end<EdgeBase>());
		}

		collect_associated_unmarked_faces(qFaces, grid,
							m_selMarkedElements.begin<Volume>(),
							m_selMarkedElements.end<Volume>());
	}

	if(grid.num_volumes() > 0){
		if(automarkHigherDimensionalObjects){
		//	automarking of higher dimensional elements
		//	if an edge or face is selected, associated volumes will be refined, too.
			collect_associated_unmarked_volumes(qVols, grid,
								m_selMarkedElements.begin<EdgeBase>(),
								m_selMarkedElements.end<EdgeBase>());

			collect_associated_unmarked_volumes(qVols, grid,
								m_selMarkedElements.begin<Face>(),
								m_selMarkedElements.end<Face>());
		}
	}

//	add unmarked constraining edges and faces of constrained ones to the queues.
//	NOTE: This step is only required if automarking of higher dimensional
//		  objects is disabled.
	if(!automarkHigherDimensionalObjects){
		for(ConstrainedEdgeIterator iter = m_selMarkedElements.begin<ConstrainedEdge>();
			iter != m_selMarkedElements.end<ConstrainedEdge>(); ++iter)
		{
			if(ConstrainingEdge* cge = dynamic_cast<ConstrainingEdge*>((*iter)->get_constraining_object()))
			{
				if(!is_marked(cge))
					qEdges.push(cge);
			}
			else if(ConstrainingFace* cgf = dynamic_cast<ConstrainingFace*>((*iter)->get_constraining_object()))
			{
				if(!is_marked(cgf))
					qFaces.push(cgf);
			}
		}
		for(ConstrainedTriangleIterator iter = m_selMarkedElements.begin<ConstrainedTriangle>();
			iter != m_selMarkedElements.end<ConstrainedTriangle>(); ++iter)
		{
			if(ConstrainingFace* cgf = dynamic_cast<ConstrainingFace*>((*iter)->get_constraining_object()))
			{
				if(!is_marked(cgf))
					qFaces.push(cgf);
			}
		}
		for(ConstrainedQuadrilateralIterator iter = m_selMarkedElements.begin<ConstrainedQuadrilateral>();
			iter != m_selMarkedElements.end<ConstrainedQuadrilateral>(); ++iter)
		{
			if(ConstrainingFace* cgf = dynamic_cast<ConstrainingFace*>((*iter)->get_constraining_object()))
			{
				if(!is_marked(cgf))
					qFaces.push(cgf);
			}
		}
	}



//	we'll now iterate over the queues and adjust the marks
//	as long as at least one queue contains something, we'll continue looping.
	while(!(qEdges.empty() && qFaces.empty() && qVols.empty()))
	{
	////////////////////////////////
	//	process edges in qEdges
		while(!qEdges.empty()){
		//	get and remove the first edge from the queue
			EdgeBase* e = qEdges.front();
			qEdges.pop();

		//	if the edge is already marked, we'll continue with the next one.
			if(is_marked(e))
				continue;

		//	mark the edge
			mark(e);

		//	depending on the type of the edge, we have to perform different operations
			if(ConstrainedEdge* cde = dynamic_cast<ConstrainedEdge*>(e))
			{
			//	the edge is a constrained edge. Make sure that its constraining edge
			//	or face will be refined.
				if(ConstrainingEdge* cge = dynamic_cast<ConstrainingEdge*>(
											cde->get_constraining_object()))
				{
					if(!is_marked(cge))
						qEdges.push(cge);
				}
				else if(ConstrainingFace* cgf = dynamic_cast<ConstrainingFace*>(
											cde->get_constraining_object()))
				{
					if(!is_marked(cgf))
						qFaces.push(cgf);
				}
				else{
					assert(!"This point shouldn't be reached. A constrained edge should always be constrained by a constraining object.");
				}
			}
			else if(ConstrainingEdge* cge = dynamic_cast<ConstrainingEdge*>(e))
			{
			//	associated faces and volumes have to be marked
				if(grid.num_faces() > 0){
					CollectFaces(vFaces, grid, cge);
					for(size_t i = 0; i < vFaces.size(); ++i){
						if(!is_marked(vFaces[i]))
							qFaces.push(vFaces[i]);
					}
				}

				if(grid.num_volumes() > 0){
					CollectVolumes(vVols, grid, cge);
					for(size_t i = 0; i < vVols.size(); ++i){
						if(!is_marked(vVols[i]))
							qVols.push(vVols[i]);
					}
				}
			}
			// normal edges don't require a special treatment at this point.
		}

	////////////////////////////////
	//	process faces in qFaces
		while(!qFaces.empty()){
		//	get and remove the first face from the queue
			Face* f = qFaces.front();
			qFaces.pop();

		//	if the face is already marked, we'll continue with the next one.
			if(is_marked(f))
				continue;

		//	mark the face
			mark(f);

		//	we have to make sure that all associated edges are marked.
			CollectEdges(vEdges, grid, f);
			for(size_t i = 0; i < vEdges.size(); ++i){
				if(!is_marked(vEdges[i]))
					qEdges.push(vEdges[i]);
			}

		//	constrained and constraining faces require special treatment
			if(ConstrainedFace* cdf = dynamic_cast<ConstrainedFace*>(f)){
			//	make sure that its constraining face will be refined
				if(ConstrainingFace* cgf = dynamic_cast<ConstrainingFace*>(
											cdf->get_constraining_object()))
				{
					if(!is_marked(cgf))
						qFaces.push(cgf);
				}
				else{
					assert(!"This point shouldn't be reached. A constrained face should always be constrained by a constraining face.");
				}
			}
			else if(ConstrainingFace* cgf = dynamic_cast<ConstrainingFace*>(f)){
			//	associated volumes have to be marked
				if(grid.num_volumes() > 0){
					CollectVolumes(vVols, grid, cgf);
					for(size_t i = 0; i < vVols.size(); ++i){
						if(!is_marked(vVols[i]))
							qVols.push(vVols[i]);
					}
				}
			}
		}

	////////////////////////////////
	//	process volumes in qVols
		while(!qVols.empty()){
		//	get and remove the first volume from the queue
			Volume* v = qVols.front();
			qVols.pop();

		//	if the volume is already marked, we'll continue with the next one.
			if(is_marked(v))
				continue;

		//	mark the volume
			mark(v);

		//	we have to make sure that all associated edges and faces are marked.
			CollectEdges(vEdges, grid, v);
			for(size_t i = 0; i < vEdges.size(); ++i){
				if(!is_marked(vEdges[i]))
					qEdges.push(vEdges[i]);
			}

			CollectFaces(vFaces, grid, v);
			for(size_t i = 0; i < vFaces.size(); ++i){
				if(!is_marked(vFaces[i]))
					qFaces.push(vFaces[i]);
			}
		}
	}
}

template <class TIterator>
void HangingNodeRefinerBase::
collect_associated_unmarked_edges(std::queue<EdgeBase*>& qEdgesOut, Grid& grid,
						 		  TIterator elemsBegin, TIterator elemsEnd)
{
	vector<EdgeBase*> vEdges;
	for(TIterator iter = elemsBegin; iter != elemsEnd; ++iter){
		CollectEdges(vEdges, grid, *iter);
		for(size_t i = 0; i < vEdges.size(); ++i){
			if(!is_marked(vEdges[i]))
				qEdgesOut.push(vEdges[i]);
		}
	}
}

template <class TIterator>
void HangingNodeRefinerBase::
collect_associated_unmarked_faces(std::queue<Face*>& qFacesOut, Grid& grid,
						 		  TIterator elemsBegin, TIterator elemsEnd)
{
	vector<Face*> vFaces;
	for(TIterator iter = elemsBegin; iter != elemsEnd; ++iter){
		CollectFaces(vFaces, grid, *iter);
		for(size_t i = 0; i < vFaces.size(); ++i){
			if(!is_marked(vFaces[i]))
				qFacesOut.push(vFaces[i]);
		}
	}
}

template <class TIterator>
void HangingNodeRefinerBase::
collect_associated_unmarked_volumes(std::queue<Volume*>& qVolsOut, Grid& grid,
									TIterator elemsBegin, TIterator elemsEnd)
{
	vector<Volume*> vVols;
	for(TIterator iter = elemsBegin; iter != elemsEnd; ++iter){
		CollectVolumes(vVols, grid, *iter);
		for(size_t i = 0; i < vVols.size(); ++i){
			if(!is_marked(vVols[i]))
				qVolsOut.push(vVols[i]);
		}
	}
}

////////////////////////////////////////////////////////////////////////
//	implementation of refine-methods.
void HangingNodeRefinerBase::
refine_constraining_edge(ConstrainingEdge* cge)
{
//	make sure that there is one hanging vertex and two constrained edges.
	assert(cge->num_constrained_vertices() == 1 && "bad number of constrained vertices. There has to be exactly 1.");
	assert(cge->num_constrained_edges() == 2 && "bad number of constrained edges. There have to be exactly 2.");

//	the grid
	Grid& grid = *m_pGrid;

//	the central hanging vertex has to be transformed into a normal vertex
	HangingVertex* centralHV = dynamic_cast<HangingVertex*>(*cge->constrained_vertices_begin());

	if(!centralHV){
		UG_LOG("The central hanging vertex of a constraining edge is missing. ignoring edge.\n");
		return;
	}

//	replace the central vertex with a normal vertex
	Vertex* centerVrt = *grid.create_and_replace<Vertex>(centralHV);

//	Iterate over the constrained edges.
//	Unmarked constrained edges will be replaced by a normal edge.
//	Marked ones will be replaced by a ConstrainingEdge. Additionally
//	associated constrained edges will be created together with the
//	new central vertex.
	for(EdgeBaseIterator iter = cge->constrained_edges_begin();
		iter != cge->constrained_edges_end(); ++iter)
	{
		if(is_marked(*iter)){
			refine_edge_with_hanging_vertex(*iter);
		}
		else{
		//	the constrained-edge can be transformed to a normal edge
			grid.create_and_replace<Edge>(*iter);
		}
	}

	cge->clear_constrained_objects();
	set_center_vertex(cge, centerVrt);
}

void HangingNodeRefinerBase::
refine_edge_with_normal_vertex(EdgeBase* e, VertexBase** newCornerVrts)
{
//	the grid
	Grid& grid = *m_pGrid;

	Vertex* nVrt = *grid.create<Vertex>(e);
	set_center_vertex(e, nVrt);

//	allow refCallback to calculate a new position
	if(m_refCallback)
		m_refCallback->new_vertex(nVrt, e);

//	split the edge
	vector<EdgeBase*> vEdges(2);
	e->refine(vEdges, nVrt, newCornerVrts);
	assert((vEdges.size() == 2) && "EdgeBase::refine - produced wrong number of edges.");
	grid.register_element(vEdges[0], e);
	grid.register_element(vEdges[1], e);
}

void HangingNodeRefinerBase::
refine_edge_with_hanging_vertex(EdgeBase* e, VertexBase** newCornerVrts)
{
	Grid& grid = *m_pGrid;
//	we have to insert a hanging node.
//	e has to be transformed to a constraining edge at the same time.
	assert(!ConstrainingEdge::type_match(e) && "invalid operation. e is a constraining edge.");
	//assert(!ConstrainedEdge::type_match(e) && "invalid operation. e is a constrained edge.");

	ConstrainingEdge* ce = *grid.create_and_replace<ConstrainingEdge>(e);

	HangingVertex* hv = *grid.create<HangingVertex>(ce);

//	allow refCallback to calculate a new position
	if(m_refCallback)
		m_refCallback->new_vertex(hv, ce);

	set_center_vertex(ce, hv);
	hv->set_parent(ce);
	ce->add_constrained_object(hv);

	hv->set_local_coordinate_1(0.5);

//	two constrained edges have to be created.
	ConstrainedEdge* nEdge[2];
	if(newCornerVrts){
		nEdge[0] = *grid.create<ConstrainedEdge>(EdgeDescriptor(newCornerVrts[0], hv), ce);
		nEdge[1] = *grid.create<ConstrainedEdge>(EdgeDescriptor(hv, newCornerVrts[1]), ce);
	}
	else{
		nEdge[0] = *grid.create<ConstrainedEdge>(EdgeDescriptor(ce->vertex(0), hv), ce);
		nEdge[1] = *grid.create<ConstrainedEdge>(EdgeDescriptor(hv, ce->vertex(1)), ce);
	}

	for(uint i = 0; i < 2; ++i)
	{
		ce->add_constrained_object(nEdge[i]);
		nEdge[i]->set_constraining_object(ce);
	}
}

void HangingNodeRefinerBase::
refine_face_with_normal_vertex(Face* f, VertexBase** newCornerVrts)
{
//UG_LOG("refine_face_with_normal_vertex\n");
	Grid& grid = *m_pGrid;

	vector<EdgeBase*> 	vEdges(f->num_edges());
	vector<VertexBase*> vNewEdgeVertices(f->num_edges());
	vector<Face*>		vFaces(f->num_vertices());// heuristic

//	collect all associated edges.
	CollectEdges(vEdges, grid, f);
	size_t numEdges = vEdges.size();

	assert(numEdges == f->num_edges() && "ERROR in RefineFaceWithNormalVertex(...): associated edges missing.");

//	each should have an associated vertex. sort them into vNewEdgeVertices.
	for(size_t i = 0; i < numEdges; ++i)
	{
		EdgeBase* e = vEdges[i];
		int edgeIndex = GetEdgeIndex(f, e);

		assert((edgeIndex >= 0) && (edgeIndex < (int)vEdges.size()) && "ERROR in RefineFaceWithNormalVertex(...): unknown problem in CollectEdges / GetEdgeIndex.");
		assert((get_center_vertex(e) != NULL) && "ERROR in RefineFaceWithNormalVertex(...): no new vertex on refined edge.");
		vNewEdgeVertices[edgeIndex] = get_center_vertex(e);
	}

//	we'll perform a regular refine
	VertexBase* nVrt = NULL;
	f->refine_regular(vFaces, &nVrt, vNewEdgeVertices, NULL,
					  Vertex(), newCornerVrts);

//	if a new vertex has been created during refine, then register it at the grid.
	if(nVrt)
	{
		grid.register_element(nVrt, f);

	//	allow refCallback to calculate a new position
		if(m_refCallback)
			m_refCallback->new_vertex(nVrt, f);
	}

	for(uint i = 0; i < vFaces.size(); ++i)
	{
		grid.register_element(vFaces[i], f);
	}

	set_center_vertex(f, nVrt);
}

void HangingNodeRefinerBase::
refine_face_with_hanging_vertex(Face* f, VertexBase** newCornerVrts)
{
//UG_LOG("refine_face_with_hanging_vertex\n");
	Grid& grid = *m_pGrid;
	size_t numVrts = f->num_vertices();

	vector<EdgeBase*> 	vEdges(f->num_edges());
	vector<VertexBase*> vNewEdgeVertices(f->num_edges());
	vector<Face*>		vFaces(numVrts);// heuristic
//	collect all associated edges.
	CollectEdges(vEdges, grid, f);
	size_t numEdges = vEdges.size();

	assert(numEdges == f->num_edges() && "ERROR in RefineFaceWithNormalVertex(...): associated edges missing.");

//	each should have an associated vertex. sort them into vNewEdgeVertices.
	for(size_t i = 0; i < numEdges; ++i)
	{
		EdgeBase* e = vEdges[i];
		int edgeIndex = GetEdgeIndex(f, e);

		assert((edgeIndex >= 0) && (edgeIndex < (int)vEdges.size()) && "ERROR in RefineFaceWithNormalVertex(...): unknown problem in CollectEdges / GetEdgeIndex.");
		assert((get_center_vertex(e) != NULL) && "ERROR in RefineFaceWithNormalVertex(...): no new vertex on refined edge.");
		vNewEdgeVertices[edgeIndex] = get_center_vertex(e);
	}


	ConstrainingFace* cgf = NULL;
	HangingVertex* hv = NULL;

//	the face has to be replaced by a constraining face.
//	we'll perform a switch here depending on the number of vertices
	switch(numVrts)
	{
		case 3:
			{
			//	create the constraining triangle and replace the old face.
				cgf = *grid.create_and_replace<ConstrainingTriangle>(f);

			//	create the constrained faces.
			//	the following triangle will not be registered at the grid. Just a temporary one.
				ConstrainedTriangle constrainedTri(cgf->vertex(0),
													cgf->vertex(1),
													cgf->vertex(2));

			//	refine the constrained tri
				VertexBase* tmpVrt;
				constrainedTri.refine(vFaces, &tmpVrt, &vNewEdgeVertices.front(),
									  NULL, newCornerVrts);
			}
			break;
		case 4:
			{
				cgf = *grid.create_and_replace<ConstrainingQuadrilateral>(f);

			//	a central hanging vertex is required
				hv = *grid.create<HangingVertex>(cgf);

			//	allow refCallback to calculate a new position
				if(m_refCallback)
					m_refCallback->new_vertex(hv, cgf);

				set_center_vertex(cgf, hv);
				hv->set_parent(cgf);
				cgf->add_constrained_object(hv);
				hv->set_local_coordinates(0.5, 0.5);

			//	create the constrained faces.
			//	the following quadrilateral will not be registered at the grid. Just a temporary one.
				ConstrainedQuadrilateral cdf(cgf->vertex(0), cgf->vertex(1),
											 cgf->vertex(2), cgf->vertex(3));

			//	refine the constrained quad
				VertexBase* tmpVrt;
				cdf.refine(vFaces, &tmpVrt, &vNewEdgeVertices.front(), hv, newCornerVrts);
			}
			break;
		default:
			assert(!"unsupported element type.");
			break;
	}

	if(cgf)
	{
	//	register the new faces
		for(size_t i = 0; i < vFaces.size(); ++i)
		{
			ConstrainedFace* cdf = dynamic_cast<ConstrainedFace*>(vFaces[i]);
			assert(cdf && "constrained face refine did produce faces which are not constrained.");
			if(cdf)
			{
				grid.register_element(cdf, cgf);
				cdf->set_constraining_object(cgf);
				cgf->add_constrained_object(cdf);
			}
		}

	//	we have to link the new constrained edges which have been auto-generated between the constrained faces.
	//	Since only edges that lie inside of the constraining face are newly created, and since only those
	//	have to be linked with the constraining face, the following algorithm will be ok for
	//	triangles and quadrilaterals.
	//	Check for each new edge-vertex, if an edge exists with the new center vertex or with it's next neighbor.
		uint numNewEdgeVertices = vNewEdgeVertices.size();
		for(uint i = 0; i < numNewEdgeVertices; ++i)
		{
			if(hv)
			{
				ConstrainedEdge* e = dynamic_cast<ConstrainedEdge*>(grid.get_edge(vNewEdgeVertices[i], hv));
				if(e)
				{
				//	link e with the constraining face
					e->set_constraining_object(cgf);
					cgf->add_constrained_object(e);
				}
			}
			else{
			//	check if a constrained edge exists between the vertex and its next neighbor
				VertexBase* vNext = vNewEdgeVertices[(i + 1) % numNewEdgeVertices];
				ConstrainedEdge* e = dynamic_cast<ConstrainedEdge*>(grid.get_edge(vNewEdgeVertices[i], vNext));
				if(e)
				{
				//	link e with the constraining face
					e->set_constraining_object(cgf);
					cgf->add_constrained_object(e);
				}
			}
		}
	}
}

void HangingNodeRefinerBase::
refine_constraining_face(ConstrainingFace* cgf)
{
	size_t numVrts = cgf->num_vertices();

//	make sure that there is one hanging vertex and two constrained edges.
	assert(cgf->num_constrained_edges() == numVrts && "bad number of constrained edges. There have to be as many as vertices.");
	assert(cgf->num_constrained_faces() == 4 && "bad number of constrained faces. There have to be exactly 4.");

//	the grid
	Grid& grid = *m_pGrid;

	HangingVertex* centralHV = NULL;
	Vertex* centerVrt = NULL;

	if(numVrts == 4){
	//	the central hanging vertex has to be transformed into a normal vertex
		centralHV = dynamic_cast<HangingVertex*>(*cgf->constrained_vertices_begin());

		if(!centralHV){
			UG_LOG("The central hanging vertex of a constraining face is missing. ignoring face.\n");
			return;
		}

	//	replace the central vertex with a normal vertex
		centerVrt = *grid.create_and_replace<Vertex>(centralHV);
	}

//	Iterate over the constrained edges.
//	Unmarked constrained edges will be replaced by a normal edge.
//	Marked ones will be replaced by a ConstrainingEdge. Additionally
//	associated constrained edges will be created together with the
//	new central vertex.
	for(EdgeBaseIterator iter = cgf->constrained_edges_begin();
		iter != cgf->constrained_edges_end(); ++iter)
	{
		if(is_marked(*iter)){
			refine_edge_with_hanging_vertex(*iter);
		}
		else{
		//	the constrained-edge can be transformed to a normal edge
			grid.create_and_replace<Edge>(*iter);
		}
	}

//	iterate over the constrained faces.
//	If it is marked, we'll replace it by a constraining face and create
//	associated constrained faces.
//	if not, it will simply be transformed to a normal face.
//	To ease implementation we will transform it anyway and if required we will
//	call refine_face_with_hanging_vertex(...).
	for(FaceIterator iter = cgf->constrained_faces_begin();
		iter != cgf->constrained_faces_end(); ++iter)
	{
		Face* f = *iter;
		if(is_marked(f)){
		//	refine it using hanging_node_refinement.
			refine_face_with_hanging_vertex(f);
		}
		else{
		//	replace it by a normal face
			if(f->num_vertices() == 3)
				f = *grid.create_and_replace<Triangle>(f);
			else
				f = *grid.create_and_replace<Quadrilateral>(f);
		}
	}

	cgf->clear_constrained_objects();
	set_center_vertex(cgf, centerVrt);
}

void HangingNodeRefinerBase::
refine_volume_with_normal_vertex(Volume* v, VertexBase** newCornerVrts)
{
	Grid& grid = *m_pGrid;

	vector<EdgeBase*> 	vEdges(v->num_edges());
	vector<VertexBase*> vNewEdgeVertices(v->num_edges());
	vector<Face*>		vFaces(v->num_faces());
	vector<VertexBase*>	vNewFaceVertices(v->num_faces());
	vector<Volume*>		vVolumes(8);// heuristic
//	collect all associated edges.
	CollectEdges(vEdges, grid, v);
	uint numEdges = vEdges.size();

	assert(numEdges == v->num_edges() && "ERROR in RefineVolumeWithNormalVertex(...): associated edges missing.");

//	each should have an associated vertex. sort them into vNewEdgeVertices.
	for(uint i = 0; i < numEdges; ++i)
	{
		EdgeBase* e = vEdges[i];
		int edgeIndex = GetEdgeIndex(v, e);

		assert((edgeIndex >= 0) && (edgeIndex < (int)vEdges.size()) && "ERROR in RefineVolumeWithNormalVertex(...): unknown problem in CollectEdges / GetEdgeIndex.");
		assert((get_center_vertex(e) != NULL) && "ERROR in RefineVolumeWithNormalVertex(...): no new vertex on refined edge.");
		vNewEdgeVertices[edgeIndex] = get_center_vertex(e);
	}

//	get the vertices contained by adjacent faces
	CollectFaces(vFaces, grid, v);
	assert(vFaces.size() == v->num_faces() && "Bad number of adjacent faces.");

	for(uint i = 0; i < vFaces.size(); ++i)
	{
		Face* f = vFaces[i];
		int faceIndex = GetFaceIndex(v, f);

		if(f->num_vertices() == 3)
			vNewFaceVertices[faceIndex] = NULL;
		else{
			assert(get_center_vertex(f) && "center-vertex of associated face is missing.");
			vNewFaceVertices[faceIndex] = get_center_vertex(f);
		}
	}

//	refine the volume and register new volumes at the grid.
	VertexBase* createdVrt = NULL;
	v->refine(vVolumes, &createdVrt, &vNewEdgeVertices.front(),
			  &vNewFaceVertices.front(), NULL, Vertex(), newCornerVrts);

	if(createdVrt){
	//	register the new vertex
		grid.register_element(createdVrt, v);

	//	allow refCallback to calculate a new position
		if(m_refCallback)
			m_refCallback->new_vertex(createdVrt, v);
	}

	for(uint i = 0; i < vVolumes.size(); ++i)
		grid.register_element(vVolumes[i], v);
}

}// end of namespace
