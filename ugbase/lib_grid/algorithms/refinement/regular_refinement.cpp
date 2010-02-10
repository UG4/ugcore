// created by Sebastian Reiter
// s.b.reiter@googlemail.com
//	y09 m03 d16

#include <vector>
#include "regular_refinement.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"

using namespace std;

namespace ug
{

static void AdjustSelection(Grid& grid, Selector& sel)
{
//	select all edges of selected faces
	vector<EdgeBase*> vEdges;
	for(FaceIterator iter = sel.begin<Face>(); iter != sel.end<Face>(); ++iter){
		CollectEdges(vEdges, grid, *iter);
		for(size_t i = 0; i < vEdges.size(); ++i)
			sel.select(vEdges[i]);
	}
	
//	select all faces which are adjacent to selected edges
	vector<Face*> vFaces;
	for(EdgeBaseIterator iter = sel.begin<EdgeBase>();
		iter != sel.end<EdgeBase>(); ++iter)
	{
		CollectFaces(vFaces, grid, *iter);
		for(size_t i = 0; i < vFaces.size(); ++i)
			sel.select(vFaces[i]);
	}
}

bool Refine(Grid& grid, Selector& sel, AInt& aInt)
{
//	aInt has to be attached to the edges of the grid
	if(!grid.has_edge_attachment(aInt)){
		LOG("  WARNING in Refine: aInt is not attached to the edges of the grid. Aborting.\n");
		return false;
	}
	
//	position data is required
	if(!grid.has_vertex_attachment(aPosition)){
		LOG("  WARNING in Refine: aPosition is not attached to the vertices of the grid. Aborting.\n");
		return false;
	}
	
//	make sure that FACEOPT_AUTOGENERATE_EDGES is enabled
	if(!grid.option_is_enabled(FACEOPT_AUTOGENERATE_EDGES)){
		LOG("  INFO in Refine: autoenabling FACEOPT_AUTOGENERATE_EDGES\n");
		grid.enable_options(FACEOPT_AUTOGENERATE_EDGES);
	}
	
	if(!grid.option_is_enabled(GRIDOPT_VERTEXCENTRIC_INTERCONNECTION)){
		LOG("  INFO in Refine: autoenabling GRIDOPT_VERTEXCENTRIC_INTERCONNECTION\n");
		grid.enable_options(GRIDOPT_VERTEXCENTRIC_INTERCONNECTION);
	}
	
LOG("1");
//	adjust selection
	AdjustSelection(grid, sel);
	
//	we need several arrays.
//	this one stores pointers to all edges that shall be refined
	vector<EdgeBase*> edges(sel.num<EdgeBase>());
//	one that stores the new vertex for each edge that shall be refined.
	vector<Vertex*>	edgeVrts(sel.num<EdgeBase>());
//	one that stores the selected faces
	vector<Face*> faces(sel.num<Face>());
//	one that stores the vertex for each edge of each face
//	entries will be set to NULL if the associated edge will not be refined
	vector<VertexBase*> faceEdgeVrts;
//	one that stores the offset to faceEdgeVrtInds for each face.
//	size is num_faces + 1 (the last entry marks the end)
	vector<int> faceOffsets(sel.num<Face>() + 1);

//	acces the int-attachment
	Grid::EdgeAttachmentAccessor<AInt> aaInt(grid, aInt);
	
//	access the position-attachment
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);
	
LOG("2");
//	fill the edges- and edgeVrts-array and assign indices to selected edges
	{
		EdgeBaseIterator edgesEnd = sel.end<EdgeBase>();
		int i = 0;
		for(EdgeBaseIterator iter = sel.begin<EdgeBase>();
			iter != edgesEnd; ++iter, ++i)
		{
		//	store the edge
			edges[i] = *iter;
			aaInt[*iter] = i;
		//	create the new vertex
			edgeVrts[i] = *grid.create<Vertex>(*iter);
		//	caculate new position
			aaPos[edgeVrts[i]] = CalculateCenter(*iter, aaPos);
		}
	}
	
LOG("3");
//	set up face arrays
	{
	//	this estimate will be exact in most cases
		faceEdgeVrts.reserve(sel.num<Triangle>() * 3 + sel.num<Quadrilateral>() * 4);
		
		FaceIterator facesEnd = sel.end<Face>();
		int i = 0;
		for(FaceIterator iter = sel.begin<Face>();
			iter != facesEnd; ++iter, ++i)
		{
			Face* f = *iter;
			faces[i] = f;
			
		//	assign the offset
			faceOffsets[i] = faceEdgeVrts.size();
			
		//	assign the edge indices
			for(uint j = 0; j < f->num_edges(); ++j){
				EdgeBase* e = grid.get_edge(f, j);
				if(sel.is_selected(e))
					faceEdgeVrts.push_back(edgeVrts[aaInt[e]]);
				else
					faceEdgeVrts.push_back(NULL);
			}
		}
		
	//	assign the last entry in faceOffset
	//	(required to retrieve the number of edges in that face later on)
		faceOffsets[i] = faceEdgeVrts.size();
	}
	
LOG("4");
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
	
LOG("5");
//	refine the selected faces
	vector<Face*> newFaces;
	newFaces.reserve(4);
	
LOG(endl);
	for(size_t i = 0; i < faces.size(); ++i){
		Face* f = faces[i];
		VertexBase* newVrt;
		if(f->refine(newFaces, &newVrt, &faceEdgeVrts[faceOffsets[i]])){
		//	if a new vertex was generated, we have to register it
			if(newVrt)
				grid.register_element(newVrt, f);

		//	register the new faces
			for(size_t j = 0; j < newFaces.size(); ++j)
				grid.register_element(newFaces[j], f);
		}
		else{
			LOG("  WARINING in Refine: could not refine face.\n");
		}
	}
	
LOG("6");
//	erase old faces
	grid.erase(faces.begin(), faces.end());
//	erase old edges
	grid.erase(edges.begin(), edges.end());
LOG("7");
}

}//	end of namespace
