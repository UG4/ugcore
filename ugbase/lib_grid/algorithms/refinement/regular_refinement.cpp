// created by Sebastian Reiter
// s.b.reiter@googlemail.com
//	y09 m03 d16

#include <vector>
#include "regular_refinement.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"

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
	vector<EdgeBase*> vEdges;
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
	for(EdgeBaseIterator iter = sel.begin<EdgeBase>();
		iter != sel.end<EdgeBase>(); ++iter)
	{
		CollectFaces(vFaces, grid, *iter);
		for(size_t i = 0; i < vFaces.size(); ++i)
			sel.select(vFaces[i]);
			
		CollectVolumes(vVols, grid, *iter);
		for(size_t i = 0; i < vVols.size(); ++i)
			sel.select(vVols[i]);
	}
}
/*
bool Refine(Grid& grid, Selector& sel, AInt& aInt)
{
//	aInt has to be attached to the edges of the grid
	if(grid.num<Face>() > 0 &! grid.has_edge_attachment(aInt)){
		LOG("  WARNING in Refine: aInt is not attached to the edges of the grid. Aborting.\n");
		return false;
	}

//	if there are volumes in the grid, 
//	aInt has to be attached to the faces of the grid
	if(grid.num<Volume>() &! grid.has_face_attachment(aInt)){
		LOG("  WARNING in Refine: aInt is not attached to the faces of the grid. Aborting.\n");
		return false;
	}
		
//	position data is required
	if(!grid.has_vertex_attachment(aPosition)){
		LOG("  WARNING in Refine: aPosition is not attached to the vertices of the grid. Aborting.\n");
		return false;
	}

//	make sure that GRIDOPT_VERTEXCENTRIC_INTERCONNECTION is enabled
	if(!grid.option_is_enabled(GRIDOPT_VERTEXCENTRIC_INTERCONNECTION)){
		LOG("  INFO in Refine: autoenabling GRIDOPT_VERTEXCENTRIC_INTERCONNECTION\n");
		grid.enable_options(GRIDOPT_VERTEXCENTRIC_INTERCONNECTION);
	}

//	make sure that FACEOPT_AUTOGENERATE_EDGES is enabled
	if(grid.num<Face>() > 0 &! grid.option_is_enabled(FACEOPT_AUTOGENERATE_EDGES)){
		LOG("  INFO in Refine: autoenabling FACEOPT_AUTOGENERATE_EDGES\n");
		grid.enable_options(FACEOPT_AUTOGENERATE_EDGES);
	}
	
//	if there are volumes, make sure that VOLOPT_AUTOGENERATE_FACES is enabled.
	if(grid.num<Volume>() > 0 &! grid.option_is_enabled(VOLOPT_AUTOGENERATE_FACES))
	{
		LOG("  INFO in Refine: autoenabling VOLOPT_AUTOGENERATE_FACES\n");
		grid.enable_options(VOLOPT_AUTOGENERATE_FACES);
	}
	
//	adjust selection
	AdjustSelection(grid, sel);
	
//	we need several arrays.
//	this one stores pointers to all edges that shall be refined
	vector<EdgeBase*> edges(sel.num<EdgeBase>());
//	one that stores the new vertex for each edge that shall be refined.
	vector<Vertex*>	edgeVrts(sel.num<EdgeBase>());
//	one that stores the selected faces
	vector<Face*> faces(sel.num<Face>());
//	one that stores vertices which are created on faces
	vector<VertexBase*> faceVrts;
//	one that stores selected volumes
	vector<Volume*> vols(sel.num<Volume>());
	
//	one that stores the vertex for each edge of each face
//	entries will be set to NULL if the associated edge will not be refined
	vector<VertexBase*> faceEdgeVrts;
//	one that stores the offset to faceEdgeVrts for each face.
//	size is num_faces + 1 (the last entry marks the end)
	vector<int> faceOffsets(sel.num<Face>() + 1);

//	one that stores the vertex for each edge of each volume
//	entries will be set to NULL if the associated edge will not be refined
	vector<VertexBase*> volEdgeVrts;
//	one that stores the offset to volEdgeVrts for each volume
	vector<VertexBase*> volOffsetsEV;
//	one that stores the vertex for each face of each volume
//	entries will be set to NULL if the associated face will not be refined
	vector<VertexBase*> volFaceVrts;
//	one that stores the offset to volFaceVrts for each volume
	vector<VertexBase*> volOffsetsFV;
	
//	number of faces and volumes that will be refined
	const size_t numRefFaces = sel.num<Face>();
	const size_t numRefVols = sel.num<Volume>();
	
//	acces the int-attachment
	Grid::EdgeAttachmentAccessor<AInt> aaIntEDGE;
	Grid::FaceAttachmentAccessor<AInt> aaIntFACE;
	
	if(numRefFaces > 0)
		aaIntEDGE.access(grid, aInt);
	if(numRefVols > 0)
		aaIntFACE.access(grid, aInt);
	
//	access the position-attachment
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);
	
//	fill the edges- and edgeVrts-array and assign indices to selected edges
	{
		EdgeBaseIterator edgesEnd = sel.end<EdgeBase>();
		int i = 0;
		for(EdgeBaseIterator iter = sel.begin<EdgeBase>();
			iter != edgesEnd; ++iter, ++i)
		{
		//	store the edge
			edges[i] = *iter;
			if(numRefFaces > 0)		
				aaIntEDGE[*iter] = i;
				
		//	create the new vertex
			edgeVrts[i] = *grid.create<Vertex>(*iter);
			sel.select(edgeVrts[i]);
		//	calculate new position
			aaPos[edgeVrts[i]] = CalculateCenter(*iter, aaPos);
		}
	}

//	set up face arrays
	{
	//	this estimate will be exact in most cases
		faceEdgeVrts.reserve(sel.num<Triangle>() * 3 + sel.num<Quadrilateral>() * 4 + 1);
		
		FaceIterator facesEnd = sel.end<Face>();
		int i = 0;
		for(FaceIterator iter = sel.begin<Face>();
			iter != facesEnd; ++iter, ++i)
		{
			Face* f = *iter;
			faces[i] = f;
			if(numRefVols > 0)
				aaIntFACE[f] = i;
				
		//	assign the offset
			faceOffsets[i] = faceEdgeVrts.size();
			
		//	assign the edge indices
			for(uint j = 0; j < f->num_edges(); ++j){
				EdgeBase* e = grid.get_edge(f, j);
				if(sel.is_selected(e))
					faceEdgeVrts.push_back(edgeVrts[aaIntEDGE[e]]);
				else
					faceEdgeVrts.push_back(NULL);
			}
		}
		
	//	assign the last entry in faceOffset
	//	(required to retrieve the number of edges in that face later on)
		faceOffsets[i] = faceEdgeVrts.size();
	}
	

//	refine the selected edges
	vector<EdgeBase*> newEdges;
	newEdges.reserve(2);
	for(size_t i = 0; i < edges.size(); ++i){
		EdgeBase* e = edges[i];
		if(e->refine(newEdges, edgeVrts[i])){
			for(size_t j = 0; j < newEdges.size(); ++j)
				grid.register_element(newEdges[j], e);
		}
		else{
			LOG("  WARNING in Refine: could not refine edge.\n");
		}
	}

//	refine the selected faces
	vector<Face*> newFaces;
	newFaces.reserve(4);
	
	for(size_t i = 0; i < faces.size(); ++i){
		Face* f = faces[i];
		VertexBase* newVrt;
		if(f->refine(newFaces, &newVrt, &faceEdgeVrts[faceOffsets[i]])){
		//	if a new vertex was generated, we have to register it
			if(newVrt){
				grid.register_element(newVrt, f);
				aaPos[newVrt] = CalculateCenter(f, aaPos);
				sel.select(newVrt);
			}

		//	register the new faces
			for(size_t j = 0; j < newFaces.size(); ++j)
				grid.register_element(newFaces[j], f);
		}
		else{
			LOG("  WARINING in Refine: could not refine face.\n");
		}
	}

//	set up volume arrays
	{
	//	the estimates will be exact in most cases
		volEdgeVrts.reserve(sel.num<Tetrahedron>() * 6
							+ sel.num<Hexahedron>() * 12
							+ sel.num<Pyramid>() * 8
							+ sel.num<Prism>() * 9 + 1);
		volFaceVrts.reserve(sel.num<Tetrahedron>() * 4
							+ sel.num<Hexahedron>() * 6
							+ sel.num<Pyramid>() * 5
							+ sel.num<Prism>() * 5 + 1);
							
		VolumeIterator volsEnd = sel.end<Volume>();
		int i = 0;
		for(VolumeIterator iter = sel.begin<Volume>();
			iter != volsEnd; ++iter, ++i)
		{
			Volume* v = *iter;
			vols[i] = v;
			
		//	assign vertices offsets 
			volOffsetsEV[i] = volEdgeVrts.size();
			
		}
	}
	
//	erase old faces
	grid.erase(faces.begin(), faces.end());
//	erase old edges
	grid.erase(edges.begin(), edges.end());
	
	return true;
}
*/
bool Refine(Grid& grid, Selector& sel, AInt& aInt,
			IRefinementCallback* refCallback)
{
//	aInt has to be attached to the edges of the grid
	if(grid.num<Face>() > 0 &! grid.has_edge_attachment(aInt)){
		LOG("  WARNING in Refine: aInt is not attached to the edges of the grid. Aborting.\n");
		return false;
	}

//	if there are volumes in the grid, 
//	aInt has to be attached to the faces of the grid
	if(grid.num<Volume>() &! grid.has_face_attachment(aInt)){
		LOG("  WARNING in Refine: aInt is not attached to the faces of the grid. Aborting.\n");
		return false;
	}
		
//	position data is required
	if(!grid.has_vertex_attachment(aPosition)){
		LOG("  WARNING in Refine: aPosition is not attached to the vertices of the grid. Aborting.\n");
		return false;
	}

//	if the refinement-callback is empty, use a linear one.
	RefinementCallbackLinear LinRefCallback(grid, aPosition);
	if(!refCallback)
		refCallback = &LinRefCallback;
		
//	make sure that GRIDOPT_VERTEXCENTRIC_INTERCONNECTION is enabled
	if(!grid.option_is_enabled(GRIDOPT_VERTEXCENTRIC_INTERCONNECTION)){
		LOG("  INFO in Refine: autoenabling GRIDOPT_VERTEXCENTRIC_INTERCONNECTION\n");
		grid.enable_options(GRIDOPT_VERTEXCENTRIC_INTERCONNECTION);
	}

//	make sure that FACEOPT_AUTOGENERATE_EDGES is enabled
	if(grid.num<Face>() > 0 &! grid.option_is_enabled(FACEOPT_AUTOGENERATE_EDGES)){
		LOG("  INFO in Refine: autoenabling FACEOPT_AUTOGENERATE_EDGES\n");
		grid.enable_options(FACEOPT_AUTOGENERATE_EDGES);
	}
	
//	if there are volumes, make sure that VOLOPT_AUTOGENERATE_FACES is enabled.
	if(grid.num<Volume>() > 0 &! grid.option_is_enabled(VOLOPT_AUTOGENERATE_FACES))
	{
		LOG("  INFO in Refine: autoenabling VOLOPT_AUTOGENERATE_FACES\n");
		grid.enable_options(VOLOPT_AUTOGENERATE_FACES);
	}
	
//	adjust selection
	AdjustSelection(grid, sel);
	
//	number of edges, faces and volumes that will be refined
	const size_t numRefEdges = sel.num<EdgeBase>();
	const size_t numRefFaces = sel.num<Face>();
	const size_t numRefVols = sel.num<Volume>();
	
//	we need several arrays.
//	this one stores pointers to all edges that shall be refined
	vector<EdgeBase*> edges(numRefEdges);
//	one that stores the new vertex for each edge that shall be refined.
	vector<Vertex*>	edgeVrts(numRefEdges);
//	one that stores the selected faces
	vector<Face*> faces(numRefFaces);
//	one that stores vertices which are created on faces
	vector<VertexBase*> faceVrts;
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
	
//	access the position-attachment
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);
	
////////////////////////////////
//	fill the edges- and edgeVrts-array and assign indices to selected edges
	{
		EdgeBaseIterator edgesEnd = sel.end<EdgeBase>();
		int i = 0;
		for(EdgeBaseIterator iter = sel.begin<EdgeBase>();
			iter != edgesEnd; ++iter, ++i)
		{
		//	store the edge
			edges[i] = *iter;
			if(numRefFaces > 0)		
				aaIntEDGE[*iter] = i;
				
		//	create the new vertex
			edgeVrts[i] = *grid.create<Vertex>(*iter);
			sel.select(edgeVrts[i]);
		//	calculate new position
			refCallback->new_vertex(edgeVrts[i], *iter);
		}
	}

////////////////////////////////
//	refine the selected edges
	vector<EdgeBase*> newEdges;
	newEdges.reserve(2);
	for(size_t i = 0; i < edges.size(); ++i){
		EdgeBase* e = edges[i];
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
	vector<VertexBase*> faceEdgeVrts;
	faceEdgeVrts.reserve(4);
	
	for(size_t i = 0; i < faces.size(); ++i){
		Face* f = faces[i];
		VertexBase* newVrt;
		
	//	collect vertices of associated edges
		faceEdgeVrts.clear();
		for(uint j = 0; j < f->num_edges(); ++j){
			EdgeBase* e = grid.get_edge(f, j);
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
			LOG("  WARINING in Refine: could not refine face.\n");
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
	vector<VertexBase*> volEdgeVrts;
	volEdgeVrts.reserve(4);
//	we need a container that stores the vertex for each face of a volume
//	entries will be set to NULL if the associated face will not be refined
	vector<VertexBase*> volFaceVrts;
	volFaceVrts.reserve(12);
	
	for(size_t i = 0; i < vols.size(); ++i){
		Volume* v = vols[i];
		VertexBase* newVrt;
		
	//	collect vertices of associated edges
		volEdgeVrts.clear();
		for(uint j = 0; j < v->num_edges(); ++j){
			EdgeBase* e = grid.get_edge(v, j);
			if(sel.is_selected(e))
				volEdgeVrts.push_back(edgeVrts[aaIntEDGE[e]]);
			else
				volEdgeVrts.push_back(NULL);
		}

	//	collect vertices of associated faces
		volFaceVrts.clear();
		for(uint j = 0; j < v->num_edges(); ++j){
			Face* f = grid.get_face(v, j);
			if(sel.is_selected(f))
				volFaceVrts.push_back(faceVrts[aaIntFACE[f]]);
			else
				volFaceVrts.push_back(NULL);
		}
		
		if(v->refine(newVols, &newVrt, &volEdgeVrts.front(),
					&volFaceVrts.front(), NULL, Vertex()))
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
			LOG("  WARINING in Refine: could not refine volume.\n");
		}
	}

//	erase old volumes
	grid.erase(vols.begin(), vols.end());
//	erase old faces
	grid.erase(faces.begin(), faces.end());
//	erase old edges
	grid.erase(edges.begin(), edges.end());

	return true;
}

}//	end of namespace
