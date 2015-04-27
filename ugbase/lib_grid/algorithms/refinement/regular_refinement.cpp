// created by Sebastian Reiter
// s.b.reiter@googlemail.com
//	y09 m03 d16

#include <vector>
#include "regular_refinement.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "lib_grid/algorithms/selection_util.h"
#include "lib_grid/algorithms/debug_util.h"

using namespace std;

namespace ug
{

bool Refine(Grid& grid, Selector& sel,
			IRefinementCallback* refCallback)
{
	AInt aInt;
	if(grid.num<Face>() > 0)
		grid.attach_to_edges(aInt);
	if(grid.num<Volume>() > 0)
		grid.attach_to_faces(aInt);
		
	bool bSuccess = Refine(grid, sel, aInt, refCallback);

	if(grid.num<Face>() > 0)
		grid.detach_from_edges(aInt);
	if(grid.num<Volume>() > 0)
		grid.detach_from_faces(aInt);

	return bSuccess;
}

static void AdjustSelection(Grid& grid, Selector& sel)
{
//	select all edges of selected faces
	vector<Edge*> vEdges;
	for(FaceIterator iter = sel.begin<Face>(); iter != sel.end<Face>(); ++iter){
		CollectEdges(vEdges, grid, *iter);
		for(size_t i = 0; i < vEdges.size(); ++i)
			sel.select(vEdges[i]);
	}
	
//	select all edges of selected volumes
	for(VolumeIterator iter = sel.begin<Volume>();
		iter != sel.end<Volume>(); ++iter)
	{
		CollectEdges(vEdges, grid, *iter);
		for(size_t i = 0; i < vEdges.size(); ++i)
			sel.select(vEdges[i]);
	}
	
//	select all faces and volumes which are adjacent to selected edges
	vector<Face*> vFaces;
	vector<Volume*> vVols;
	for(EdgeIterator iter = sel.begin<Edge>();
		iter != sel.end<Edge>(); ++iter)
	{
		CollectFaces(vFaces, grid, *iter);
		for(size_t i = 0; i < vFaces.size(); ++i)
			sel.select(vFaces[i]);
			
		CollectVolumes(vVols, grid, *iter);
		for(size_t i = 0; i < vVols.size(); ++i)
			sel.select(vVols[i]);
	}
}


bool Refine(Grid& grid, Selector& sel, AInt& aInt,
			IRefinementCallback* refCallback)
{
//	position data is required
	if(!grid.has_vertex_attachment(aPosition)){
		LOG("  WARNING in Refine: aPosition is not attached to the vertices of the grid. Aborting.\n");
		return false;
	}

//	if the refinement-callback is empty, use a linear one.
	RefinementCallbackLinear<APosition> LinRefCallback(grid, aPosition);
	if(!refCallback)
		refCallback = &LinRefCallback;
		
//	make sure that GRIDOPT_VERTEXCENTRIC_INTERCONNECTION is enabled
	if(grid.num_edges() && (!grid.option_is_enabled(VRTOPT_STORE_ASSOCIATED_EDGES))){
		LOG("  INFO in Refine: autoenabling VRTOPT_STORE_ASSOCIATED_EDGES\n");
		grid.enable_options(VRTOPT_STORE_ASSOCIATED_EDGES);
	}
	if(grid.num_faces() && (!grid.option_is_enabled(VRTOPT_STORE_ASSOCIATED_FACES))){
		LOG("  INFO in Refine: autoenabling VRTOPT_STORE_ASSOCIATED_FACES\n");
		grid.enable_options(VRTOPT_STORE_ASSOCIATED_FACES);
	}
	if(grid.num_volumes() && (!grid.option_is_enabled(VRTOPT_STORE_ASSOCIATED_VOLUMES))){
		LOG("  INFO in Refine: autoenabling VRTOPT_STORE_ASSOCIATED_VOLUMES\n");
		grid.enable_options(VRTOPT_STORE_ASSOCIATED_VOLUMES);
	}

//	make sure that FACEOPT_AUTOGENERATE_EDGES is enabled
	if(grid.num<Face>() > 0 && (!grid.option_is_enabled(FACEOPT_AUTOGENERATE_EDGES))){
		LOG("  INFO in Refine: autoenabling FACEOPT_AUTOGENERATE_EDGES\n");
		grid.enable_options(FACEOPT_AUTOGENERATE_EDGES);
	}
	
//	if there are volumes, make sure that VOLOPT_AUTOGENERATE_FACES is enabled.
	if(grid.num<Volume>() > 0 && (!grid.option_is_enabled(VOLOPT_AUTOGENERATE_FACES)))
	{
		LOG("  INFO in Refine: autoenabling VOLOPT_AUTOGENERATE_FACES\n");
		grid.enable_options(VOLOPT_AUTOGENERATE_FACES);
	}
	
//	adjust selection
	AdjustSelection(grid, sel);

//	we will select associated vertices, too, since we have to
//	notify the refinement-callback, that they are involved in refinement.
	sel.clear<Vertex>();
	SelectAssociatedVertices(sel, sel.begin<Edge>(), sel.end<Edge>());
	SelectAssociatedVertices(sel, sel.begin<Face>(), sel.end<Face>());
	SelectAssociatedVertices(sel, sel.begin<Volume>(), sel.end<Volume>());
	

//	aInt has to be attached to the edges of the grid
	if(sel.num<Face>() > 0 && (!grid.has_edge_attachment(aInt))){
		LOG("  WARNING in Refine: aInt is not attached to the edges of the grid. Aborting.\n");
		return false;
	}

//	if there are selected volumes,
//	aInt has to be attached to the faces of the grid
	if(sel.num<Volume>() && (!grid.has_face_attachment(aInt))){
		LOG("  WARNING in Refine: aInt is not attached to the faces of the grid. Aborting.\n");
		return false;
	}

//	number of edges, faces and volumes that will be refined
	const size_t numRefEdges = sel.num<Edge>();
	const size_t numRefFaces = sel.num<Face>();
	const size_t numRefVols = sel.num<Volume>();
	
//	we need several arrays.
//	this one stores pointers to all edges that shall be refined
	vector<Edge*> edges(numRefEdges);
//	one that stores the new vertex for each edge that shall be refined.
	vector<RegularVertex*>	edgeVrts(numRefEdges);
//	one that stores the selected faces
	vector<Face*> faces(numRefFaces);
//	one that stores vertices which are created on faces
	vector<Vertex*> faceVrts;
	if(numRefVols > 0)
		faceVrts.resize(numRefFaces);
//	one that stores selected volumes
	vector<Volume*> vols(sel.num<Volume>());
	
//	acces the int-attachment
	Grid::EdgeAttachmentAccessor<AInt> aaIntEDGE;
	Grid::FaceAttachmentAccessor<AInt> aaIntFACE;
	
	if(numRefFaces > 0)
		aaIntEDGE.access(grid, aInt);
	if(numRefVols > 0)
		aaIntFACE.access(grid, aInt);
	
//	notify refinement callbacks about encountered vertices
	for(VertexIterator iter = sel.begin<Vertex>();
		iter != sel.end<Vertex>(); ++iter)
	{
		refCallback->flat_grid_vertex_encountered(*iter);
	}

////////////////////////////////
//	fill the edges- and edgeVrts-array and assign indices to selected edges
	{
		EdgeIterator edgesEnd = sel.end<Edge>();
		int i = 0;
		for(EdgeIterator iter = sel.begin<Edge>();
			iter != edgesEnd; ++iter, ++i)
		{
		//	store the edge
			edges[i] = *iter;
			if(numRefFaces > 0)		
				aaIntEDGE[*iter] = i;
				
		//	create the new vertex
			edgeVrts[i] = *grid.create<RegularVertex>(*iter);
			sel.select(edgeVrts[i]);
		//	calculate new position
			refCallback->new_vertex(edgeVrts[i], *iter);
		}
	}

////////////////////////////////
//	refine the selected edges
	vector<Edge*> newEdges;
	newEdges.reserve(2);
	for(size_t i = 0; i < edges.size(); ++i){
		Edge* e = edges[i];
		if(e->refine(newEdges, edgeVrts[i])){
			for(size_t j = 0; j < newEdges.size(); ++j)
				grid.register_element(newEdges[j], e);
		}
		else{
			LOG("  WARNING in Refine: could not refine edge.\n");
		}
	}
	
////////////////////////////////
//	set up face arrays
	{		
		FaceIterator facesEnd = sel.end<Face>();
		int i = 0;
		for(FaceIterator iter = sel.begin<Face>();
			iter != facesEnd; ++iter, ++i)
		{
			Face* f = *iter;
			faces[i] = f;
			if(numRefVols > 0)
				aaIntFACE[f] = i;
		}
	}

////////////////////////////////
//	refine the selected faces
	vector<Face*> newFaces;
	newFaces.reserve(4);
//	we need a container that stores the vertex for each edge of a face
//	entries will be set to NULL if the associated edge will not be refined
	vector<Vertex*> faceEdgeVrts;
	faceEdgeVrts.reserve(4);
	
	for(size_t i = 0; i < faces.size(); ++i){
		Face* f = faces[i];
		Vertex* newVrt;
		
	//	collect vertices of associated edges
		faceEdgeVrts.clear();
		for(uint j = 0; j < f->num_edges(); ++j){
			Edge* e = grid.get_edge(f, j);
			if(sel.is_selected(e))
				faceEdgeVrts.push_back(edgeVrts[aaIntEDGE[e]]);
			else
				faceEdgeVrts.push_back(NULL);
		}

		if(f->refine(newFaces, &newVrt, &faceEdgeVrts.front())){
		//	if a new vertex was generated, we have to register it
			if(newVrt){
				grid.register_element(newVrt, f);
				refCallback->new_vertex(newVrt, f);
				sel.select(newVrt);
			}
			
		//	if volumes are refined too, we have to store the vertex
			if(numRefVols > 0)
				faceVrts[i] = newVrt;
				
		//	register the new faces
			for(size_t j = 0; j < newFaces.size(); ++j)
				grid.register_element(newFaces[j], f);
		}
		else{
			LOG("  WARNING in Refine: could not refine face.\n");
		}
	}

////////////////////////////////
//	set up volume arrays
	{							
		VolumeIterator volsEnd = sel.end<Volume>();
		int i = 0;
		for(VolumeIterator iter = sel.begin<Volume>();
			iter != volsEnd; ++iter, ++i)
		{
			Volume* v = *iter;
			vols[i] = v;
		}
	}
	
////////////////////////////////
//	refine the selected volumes
	vector<Volume*> newVols;
	newVols.reserve(8);
//	we need a container that stores the vertex for each edge of a volume
//	entries will be set to NULL if the associated edge will not be refined
	vector<Vertex*> volEdgeVrts;
	volEdgeVrts.reserve(12);
//	we need a container that stores the vertex for each face of a volume
//	entries will be set to NULL if the associated face will not be refined
	vector<Vertex*> volFaceVrts;
	volFaceVrts.reserve(6);
	
//	only used for tetrahedron refinement
	vector<vector3> corners(4, vector3(0, 0, 0));

// //	DEBUG
// 	UG_LOG("> VOL-REF-BEGIN\n");
// 	UG_LOG("> DEBUG-ACCESSOR...\n");
// 	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);


	for(size_t i = 0; i < vols.size(); ++i){
		Volume* v = vols[i];
		Vertex* newVrt;
		
	//	collect vertices of associated edges
		volEdgeVrts.clear();
		for(uint j = 0; j < v->num_edges(); ++j){
			Edge* e = grid.get_edge(v, j);
			if(sel.is_selected(e))
				volEdgeVrts.push_back(edgeVrts[aaIntEDGE[e]]);
			else
				volEdgeVrts.push_back(NULL);
		}

	//	collect vertices of associated faces
		volFaceVrts.clear();
		for(uint j = 0; j < v->num_faces(); ++j){
			Face* f = grid.get_face(v, j);
			if(sel.is_selected(f))
				volFaceVrts.push_back(faceVrts[aaIntFACE[f]]);
			else
				volFaceVrts.push_back(NULL);
		}
		
	//	if we're performing tetrahedral refinement, we have to collect
	//	the corner coordinates, so that the refinement algorithm may choose
	//	the best interior diagonal.
		vector3* pCorners = NULL;
		if((v->num_vertices() == 4) && refCallback){
			for(size_t i = 0; i < 4; ++i){
				refCallback->current_pos(&corners[i].x(), v->vertex(i), 3);
			}
			pCorners = &corners.front();
		}

		// UG_LOG("v " << CalculateCenter(v, aaPos) << "\n");

		if(v->refine(newVols, &newVrt, &volEdgeVrts.front(),
					&volFaceVrts.front(), NULL, RegularVertex(), NULL, pCorners))
		{
		//	if a new vertex was generated, we have to register it
			if(newVrt){
				grid.register_element(newVrt, v);
				refCallback->new_vertex(newVrt, v);
				sel.select(newVrt);
			}
				
		//	register the new volumes
			for(size_t j = 0; j < newVols.size(); ++j)
				grid.register_element(newVols[j], v);
		}
		else{
			LOG("  WARNING in Refine: could not refine volume.\n");
		}
	}
// UG_LOG("> VOL-REF-END\n");

//	erase old volumes
	grid.erase(vols.begin(), vols.end());
//	erase old faces
	grid.erase(faces.begin(), faces.end());
//	erase old edges
	grid.erase(edges.begin(), edges.end());

	return true;
}

}//	end of namespace
