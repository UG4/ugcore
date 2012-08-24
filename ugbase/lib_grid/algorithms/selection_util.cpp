// created by Sebastian Reiter
// y09 m11 d16
// s.b.reiter@googlemail.com

#include <vector>
#include <stack>
#include <queue>
#include "lib_grid/selector.h"
#include "selection_util.h"
#include "geom_obj_util/geom_obj_util.h"
#include "graph/graph.h"
#include "lib_grid/grid/grid_util.h"

using namespace std;

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	instatiate template specializations
template void SelectInnerSelectionVertices<Selector>(Selector&);
template void SelectInnerSelectionVertices<MGSelector>(MGSelector&);

template void SelectInnerSelectionEdges<Selector>(Selector&);
template void SelectInnerSelectionEdges<MGSelector>(MGSelector&);

template void SelectInnerSelectionFaces<Selector>(Selector&);
template void SelectInnerSelectionFaces<MGSelector>(MGSelector&);

////////////////////////////////////////////////////////////////////////
template void DeselectBoundarySelectionVertices<Selector>(Selector&);
template void DeselectBoundarySelectionVertices<MGSelector>(MGSelector&);

template void DeselectBoundarySelectionEdges<Selector>(Selector&);
template void DeselectBoundarySelectionEdges<MGSelector>(MGSelector&);

template void DeselectBoundarySelectionFaces<Selector>(Selector&);
template void DeselectBoundarySelectionFaces<MGSelector>(MGSelector&);

////////////////////////////////////////////////////////////////////////
template void EraseSelectedObjects<Selector>(Selector&);
template void EraseSelectedObjects<MGSelector>(MGSelector&);

////////////////////////////////////////////////////////////////////////
template void InvertSelection<Selector>(Selector&);
template void InvertSelection<MGSelector>(MGSelector&);


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
size_t CollectVerticesTouchingSelection(std::vector<VertexBase*>& vrtsOut,
										ISelector& sel)
{
	vrtsOut.clear();
	Grid* pGrid = sel.grid();
	if(!pGrid)
		return 0;

	Grid& grid = *pGrid;

	grid.begin_marking();

//	get the goc and iterate over all elements
	GeometricObjectCollection goc = sel.get_geometric_objects();
	for(size_t lvl = 0; lvl < goc.num_levels(); ++lvl){
		for(VertexBaseIterator iter = goc.begin<VertexBase>(lvl);
			iter != goc.end<VertexBase>(lvl); ++iter)
		{
			if(!grid.is_marked(*iter)){
				grid.mark(*iter);
				vrtsOut.push_back(*iter);
			}
		}

		for(EdgeBaseIterator iter = goc.begin<EdgeBase>(lvl);
			iter != goc.end<EdgeBase>(lvl); ++iter)
		{
			for(size_t i = 0; i < (*iter)->num_vertices(); ++i){
				VertexBase* vrt = (*iter)->vertex(i);
				if(!grid.is_marked(vrt)){
					grid.mark(vrt);
					vrtsOut.push_back(vrt);
				}
			}
		}

		for(FaceIterator iter = goc.begin<Face>(lvl);
			iter != goc.end<Face>(lvl); ++iter)
		{
			for(size_t i = 0; i < (*iter)->num_vertices(); ++i){
				VertexBase* vrt = (*iter)->vertex(i);
				if(!grid.is_marked(vrt)){
					grid.mark(vrt);
					vrtsOut.push_back(vrt);
				}
			}
		}

		for(VolumeIterator iter = goc.begin<Volume>(lvl);
			iter != goc.end<Volume>(lvl); ++iter)
		{
			for(size_t i = 0; i < (*iter)->num_vertices(); ++i){
				VertexBase* vrt = (*iter)->vertex(i);
				if(!grid.is_marked(vrt)){
					grid.mark(vrt);
					vrtsOut.push_back(vrt);
				}
			}
		}
	}

	grid.end_marking();

	return vrtsOut.size();
}

////////////////////////////////////////////////////////////////////////
template <class TSelector>
void EraseSelectedObjects(TSelector& sel)
{
	if(!sel.grid())
		return;
	
	Grid& grid = *sel.grid();
	
	for(size_t i = 0; i < sel.num_levels(); ++i)
	{
		EraseElements<VertexBase>(grid, sel.template begin<VertexBase>(i),
								  sel.template end<VertexBase>(i));
		EraseElements<EdgeBase>(grid, sel.template begin<EdgeBase>(i),
					  			sel.template end<EdgeBase>(i));
		EraseElements<Face>(grid, sel.template begin<Face>(i),
					  		sel.template end<Face>(i));
		EraseElements<Volume>(grid, sel.template begin<Volume>(i),
					  		  sel.template end<Volume>(i));
	}
}

////////////////////////////////////////////////////////////////////////
template <class TSelector>
void InvertSelection(TSelector& sel)
{
	if(!sel.grid())
		return;
	
	Grid& grid = *sel.grid();
	
	InvertSelection(sel, grid.begin<VertexBase>(),
					grid.end<VertexBase>());
	InvertSelection(sel, grid.begin<EdgeBase>(),
					grid.end<EdgeBase>());
	InvertSelection(sel, grid.begin<Face>(),
					grid.end<Face>());
	InvertSelection(sel, grid.begin<Volume>(),
					grid.end<Volume>());
}

////////////////////////////////////////////////////////////////////////
//	SelectAssociatedGeometricObjects
void SelectAssociatedGeometricObjects(Selector& sel, ISelector::status_t status)
{
	if(!sel.grid()){
		UG_LOG("ERROR in SelectAssociatedGeometricObjects: Selector has to be assigned to a grid.\n");
		return;
	}
	
	Grid& grid = *sel.grid();
	
//	select associated elements of selected elements
	SelectAssociatedFaces(sel, sel.begin<Volume>(),
						  sel.end<Volume>(), status);
	if(!grid.option_is_enabled(VOLOPT_AUTOGENERATE_FACES)){
	//	if faces are not automatically generated, there may be edges connected
	//	to the volume, which are not connected to a face. Same goes for vertices.
		SelectAssociatedEdges(sel, sel.begin<Volume>(),
							  sel.end<Volume>(), status);
		if(!grid.option_is_enabled(VOLOPT_AUTOGENERATE_EDGES))
			SelectAssociatedVertices(sel, sel.begin<Volume>(),
									 sel.end<Volume>(), status);
	}
	
	SelectAssociatedEdges(sel, sel.begin<Face>(), sel.end<Face>(), status);
	if(!grid.option_is_enabled(FACEOPT_AUTOGENERATE_EDGES))
		SelectAssociatedVertices(sel, sel.begin<Face>(),
								 sel.end<Face>(), status);
		
	SelectAssociatedVertices(sel, sel.begin<EdgeBase>(),
							 sel.end<EdgeBase>(), status);
}


////////////////////////////////////////////////////////////////////////
//	SelectAssociatedGeometricObjects
void SelectAssociatedGeometricObjects(MGSelector& msel, ISelector::status_t status)
{
	if(!msel.multi_grid()){
		UG_LOG("ERROR in SelectAssociatedGeometricObjects: Selector has to be assigned to a grid.\n");
		return;
	}
	
	Grid& grid = *msel.multi_grid();
	
//	select associated elements of selected elements on each level
	for(size_t i = 0; i < msel.num_levels(); ++i)
	{
		SelectAssociatedFaces(msel, msel.begin<Volume>(i),
							  msel.end<Volume>(i), status);
		if(!grid.option_is_enabled(VOLOPT_AUTOGENERATE_FACES)){
			SelectAssociatedEdges(msel, msel.begin<Volume>(i),
								  msel.end<Volume>(i), status);
			if(!grid.option_is_enabled(VOLOPT_AUTOGENERATE_EDGES))
				SelectAssociatedVertices(msel, msel.begin<Volume>(i),
										 msel.end<Volume>(i), status);
		}
		
		SelectAssociatedEdges(msel, msel.begin<Face>(i),
							  msel.end<Face>(i), status);
		if(!grid.option_is_enabled(FACEOPT_AUTOGENERATE_EDGES))
			SelectAssociatedVertices(msel, msel.begin<Face>(i),
									 msel.end<Face>(i), status);
			
		SelectAssociatedVertices(msel, msel.begin<EdgeBase>(i),
								 msel.end<EdgeBase>(i), status);
	}
}

////////////////////////////////////////////////////////////////////////
//	SelectParents
///	helper for SelectAssociatedGenealogy.
template <class TIterator>
static void SelectParents(MultiGrid& mg, MGSelector& msel,
						  TIterator iterBegin, TIterator iterEnd)
{
	while(iterBegin != iterEnd)
	{
	//	if the object has a parent, then select it.
		GeometricObject* parent = mg.get_parent(*iterBegin);
		if(parent)
			msel.select(parent);

		iterBegin++;
	}
}

////////////////////////////////////////////////////////////////////////
//	ExtendSelection
void ExtendSelection(Selector& sel, size_t extSize)
{
	if(!sel.grid()){
		UG_LOG("ERROR in ExtendSelection: Selector has to be assigned to a grid.\n");
		return;
	}
	
	Grid& grid = *sel.grid();
	
//	first select associated elements of volumes, faces and edges.
//	then select associated elements of selected vertices.
//	do this extSize times.
//	elements that have already been processed are marked.
	
	grid.begin_marking();
	
//	perform iteration
	for(size_t extIters = 0; extIters < extSize; ++extIters)
	{
//TODO: speed-up by only calling SelectAssociatedGeometricObjects once before the loop.
//		During the loop only newly selected elements should be checked for associated elements.

	//	select associated elements
		SelectAssociatedGeometricObjects(sel);
	
	//	iterate over all selected vertices.
		for(VertexBaseIterator iter = sel.begin<VertexBase>();
			iter != sel.end<VertexBase>(); ++iter)
		{
			VertexBase* vrt = *iter;
		//	all marked vertices have already been processed.
			if(!grid.is_marked(vrt)){
				grid.mark(vrt);
				
			//	select associated volumes, faces and edges.
				for(Grid::AssociatedEdgeIterator asIter = grid.associated_edges_begin(vrt);
					asIter != grid.associated_edges_end(vrt); ++asIter)
				{
					sel.select(*asIter);
				}
				
				for(Grid::AssociatedFaceIterator asIter = grid.associated_faces_begin(vrt);
					asIter != grid.associated_faces_end(vrt); ++asIter)
				{
					sel.select(*asIter);
				}
				
				for(Grid::AssociatedVolumeIterator asIter = grid.associated_volumes_begin(vrt);
					asIter != grid.associated_volumes_end(vrt); ++asIter)
				{
					sel.select(*asIter);
				}
			}
		}
	}
	
	grid.end_marking();
}

template <class TGeomObj>
void SelectionFill(Selector& sel)
{
	typedef typename geometry_traits<TGeomObj>::iterator GeomObjIter;
	typedef typename TGeomObj::lower_dim_base_object Side;

	if(sel.grid() == 0){
		UG_LOG("WARNING in SelectionFill: A grid has to be assigned! Aborting.\n");
		return;
	}

	Grid& grid = *sel.grid();

	vector<Side*> sides;
	vector<TGeomObj*> objs;
	queue<TGeomObj*> qCandidates;

//	all initially selected objects are candidates
	for(GeomObjIter iter = sel.begin<TGeomObj>();
		iter != sel.end<TGeomObj>(); ++iter)
	{
		qCandidates.push(*iter);
	}

//	while there are candidates left
	while(!qCandidates.empty())
	{
	//	get the candidate
		TGeomObj* o = qCandidates.front();
		qCandidates.pop();

	//	collect all sides in a vector
		CollectAssociated(sides, grid, o);
		for(size_t i = 0; i < sides.size(); ++i){
		//	if the side is selected, we don't have to process it
			if(sel.is_selected(sides[i]))
				continue;

		//	all associated unselected geom-objs have to be selected
		//	and are new candidates
			CollectAssociated(objs, grid, sides[i]);
			for(size_t j = 0; j < objs.size(); ++j){
				if(!sel.is_selected(objs[j])){
					sel.select(objs[j]);
					qCandidates.push(objs[j]);
				}
			}
		}
	}
}

//	Only those template specializations make sense.
template void SelectionFill<EdgeBase>(Selector&);
template void SelectionFill<Face>(Selector&);
template void SelectionFill<Volume>(Selector&);



////////////////////////////////////////////////////////////////////////
//	SelectAssociatedGenealogy
void SelectAssociatedGenealogy(MGSelector& msel, bool selectAssociatedElements)
{
	MultiGrid* mg = msel.multi_grid();
	if(!mg)
		return;

//	we'll iterate through the levels from top to bottom.
//	in each level we'll select the parents of all selected elements.
	for(int i = (int)msel.num_levels() - 1; i >= 0; --i)
	{
		if(selectAssociatedElements)
		{
			SelectAssociatedVertices(msel, msel.begin<EdgeBase>(i), msel.end<EdgeBase>(i));
			SelectAssociatedVertices(msel, msel.begin<Face>(i), msel.end<Face>(i));
			SelectAssociatedVertices(msel, msel.begin<Volume>(i), msel.end<Volume>(i));
			
			SelectAssociatedEdges(msel, msel.begin<Face>(i), msel.end<Face>(i));
			SelectAssociatedEdges(msel, msel.begin<Volume>(i), msel.end<Volume>(i));
			
			SelectAssociatedFaces(msel, msel.begin<Volume>(i), msel.end<Volume>(i));
		}
		if(i > 0)
		{
			SelectParents(*mg, msel, msel.vertices_begin(i), msel.vertices_end(i));
			SelectParents(*mg, msel, msel.edges_begin(i), msel.edges_end(i));
			SelectParents(*mg, msel, msel.faces_begin(i), msel.faces_end(i));
			SelectParents(*mg, msel, msel.volumes_begin(i), msel.volumes_end(i));
		}
	}

//	thats it. done!
}

////////////////////////////////////////////////////////////////////////
//	SelectSmoothEdgePath
void SelectSmoothEdgePath(Selector& sel, number thresholdDegree,
							bool stopAtSelVrts, APosition& aPos)
{
	bool bMinimalNormalDeviation = true;
	
	if(!sel.grid())
		return;
	
	Grid& grid = *sel.grid();
	
//	access the position attachment
	assert(grid.has_vertex_attachment(aPos) &&
			"INFO in SelectSmoothEdgePath: missing position attachment.");
	
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPos);
	
//	make sure that associated edges can be easily accessed.
	if(!grid.option_is_enabled(VRTOPT_STORE_ASSOCIATED_EDGES)){
		UG_LOG("  INFO in SelectSmoothEdgePath: auto-enabling VRTOPT_STORE_ASSOCIATED_EDGES.\n");
		grid.enable_options(VRTOPT_STORE_ASSOCIATED_EDGES);
	}
	
//	convert the thresholdDegree to thresholdDot
	number thresholdDot = cos(deg_to_rad(thresholdDegree));
	
//	here we'll store candidates.
	stack<VertexBase*>	m_candidates;
	
//	initially mark all vertices of selected edges as candidates
	for(EdgeBaseIterator iter = sel.begin<EdgeBase>();
		iter != sel.end<EdgeBase>(); ++iter)
	{
	//	we don't care if vertices are pushed twice.
	//	if a vertex is selected and stopAtSelVrts is true,
	//	then we don't need to push them on the stack.
		for(size_t i = 0; i < 2; ++i){
			if(!(stopAtSelVrts && sel.is_selected((*iter)->vertex(i))))
				m_candidates.push((*iter)->vertex(i));
		}
	}
		
//	while there are candidates left
	while(!m_candidates.empty())
	{
		VertexBase* srcVrt = m_candidates.top();
		m_candidates.pop();

	//	search for associated selected edges (there has to be at last one!)
		EdgeBase* lastEdge = NULL;
		int counter = 0;
		
		for(Grid::AssociatedEdgeIterator iter = grid.associated_edges_begin(srcVrt);
			iter != grid.associated_edges_end(srcVrt); ++iter)
		{
			EdgeBase* e = *iter;
			if(sel.is_selected(e)){
				lastEdge = e;
				++counter;
			}
		}
		
		assert(lastEdge && "there has to be at least one selected associated edge!");
		
	//	if more than one associated selected edge have been found,
	//	then the vertex has already been completly handled.
		if(counter > 1)
			continue;
			
	//	the direction of the last edge
		vector3 lastDir;
		VecSubtract(lastDir, aaPos[GetConnectedVertex(lastEdge, srcVrt)],
					aaPos[srcVrt]);
		VecNormalize(lastDir, lastDir);
		
		vector3 lastNormal;
		int numInvolvedFaces = CalculateNormal(lastNormal, grid, lastEdge, aaPos);
		bool bLastNormalValid = (numInvolvedFaces > 0 && numInvolvedFaces < 3);
		
	//	follow the smooth path
		while(srcVrt)
		{
		//	check smoothness for each connected unselected edge
			EdgeBase* bestEdge = NULL;
			number bestDot = -1.1;
			number bestNormalDot = -1.1;
			vector3 bestDir(0, 0, 0);
			vector3 bestNormal(0, 0, 0);
			bool bBestNormalValid = false;
		
		//	if invalid normals are involved we have to skip further normal
		//	tests for this edge-set.
			bool ignoreNormalChecks = false;
			
			int counter = 0;
			Grid::AssociatedEdgeIterator iterEnd = grid.associated_edges_end(srcVrt);
			for(Grid::AssociatedEdgeIterator iter = grid.associated_edges_begin(srcVrt);
				iter != iterEnd; ++iter)
			{
				EdgeBase* e = *iter;
				if(!sel.is_selected(e)){
				//	check smoothness
					vector3 dir;
					VecSubtract(dir, aaPos[srcVrt],
								aaPos[GetConnectedVertex(e, srcVrt)]);
					VecNormalize(dir, dir);
					
					number d = VecDot(lastDir, dir);
					if(d > thresholdDot){
						bool moreChecks = true;
						bool bNormalValid = false;
						vector3 n(0, 0, 0);
					//	if minimal normal deviation is activated, then first try do perform
					//	the following checks. Take into account, that edges can be connected
					//	to an arbitrary amount of faces.
						if(bMinimalNormalDeviation){
							int numAdjacentFaces = CalculateNormal(n, grid, e, aaPos);
							bNormalValid = (numAdjacentFaces > 0 && numAdjacentFaces < 3);
							
							if(bLastNormalValid && bNormalValid && (!ignoreNormalChecks)){
								moreChecks = false;
							//	check whether the normal dot is better than the last one.
								number nd = VecDot(lastNormal, n);
							//	weight the dots to ensure that the better edge is found even
							//	for equal normal-dots.
								if((0.9*nd + 0.1*d) > (0.9*bestNormalDot + 0.1*bestDot)){
									bestEdge = e;
									bestDot = d;
									bestDir = dir;
									bestNormal = n;
									bestNormalDot = nd;
									bBestNormalValid = true;
								}
							}
						}
						
					//	either bMinimalNormalDeviation was false, or a normal of one
					//	of the edges was not valid.
						if(moreChecks){
							if(d > bestDot){
								bestEdge = e;
								bestDot = d;
								bestDir = dir;

								if(bMinimalNormalDeviation){
									bestNormal = n;
									bBestNormalValid = bNormalValid;
									ignoreNormalChecks = true;
								}
							}
						}
					}
				}
				else {
					counter++;
				}

			}
			
			if((bestEdge != NULL) && (counter < 2)){
				sel.select(bestEdge);
			//	the next vertex has to be checked
				srcVrt = GetConnectedVertex(bestEdge, srcVrt);
			//	make sure that we stop at selected vertices - if desired
				if(stopAtSelVrts && sel.is_selected(srcVrt))
					srcVrt = NULL;
					
				lastEdge = bestEdge;
				lastDir = bestDir;
				bLastNormalValid = bBestNormalValid;
				lastNormal = bestNormal;
			}
			else{
				srcVrt = NULL;
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////
//	SelectInnerSelectionVertices
template <class TSelector>
void SelectInnerSelectionVertices(TSelector& sel)
{
	if(!sel.grid())
		return;
	
	Grid& grid = *sel.grid();
	
	grid.begin_marking();
	
//	we'll first collect all vertices that we want to check
	vector<VertexBase*> vrts;
	
//	iterate over all levels
	for(size_t lvl = 0; lvl < sel.num_levels(); ++lvl)
	{
		for(VolumeIterator iter = sel.template begin<Volume>(lvl);
			iter != sel.template end<Volume>(lvl); ++iter)
		{
			Volume* vol = *iter;
			for(size_t i = 0; i < vol->num_vertices(); ++i){
				VertexBase* v = vol->vertex(i);
				if(!grid.is_marked(v)){
					grid.mark(v);
					vrts.push_back(v);
				}
			}
		}

		for(FaceIterator iter = sel.template begin<Face>(lvl);
			iter != sel.template end<Face>(lvl); ++iter)
		{
			Face* f = *iter;
			for(size_t i = 0; i < f->num_vertices(); ++i){
				VertexBase* v = f->vertex(i);
				if(!grid.is_marked(v)){
					grid.mark(v);
					vrts.push_back(v);
				}
			}
		}

		for(EdgeBaseIterator iter = sel.template begin<EdgeBase>(lvl);
			iter != sel.template end<EdgeBase>(lvl); ++iter)
		{
			EdgeBase* e = *iter;
			for(size_t i = 0; i < 2; ++i){
				VertexBase* v = e->vertex(i);
				if(!grid.is_marked(v)){
					grid.mark(v);
					vrts.push_back(v);
				}
			}
		}
	}
	grid.end_marking();


//	now check for each vertex if an unselected element is associated
	for(size_t i = 0; i < vrts.size(); ++i)
	{
		VertexBase* v = vrts[i];
	
	//	check whether all associated elements are selected
		bool foundUnselected = false;
		
	//	volumes
		for(Grid::AssociatedVolumeIterator aIter = grid.associated_volumes_begin(v);
			aIter != grid.associated_volumes_end(v); ++ aIter)
		{
			if(!sel.is_selected(*aIter)){
				foundUnselected = true;
				break;
			}
		}

	//	face		
		if(!foundUnselected){
			for(Grid::AssociatedFaceIterator aIter = grid.associated_faces_begin(v);
				aIter != grid.associated_faces_end(v); ++ aIter)
			{
				if(!sel.is_selected(*aIter)){
					foundUnselected = true;
					break;
				}
			}
		}
	
	//	edge		
		if(!foundUnselected){
			for(Grid::AssociatedEdgeIterator aIter = grid.associated_edges_begin(v);
				aIter != grid.associated_edges_end(v); ++ aIter)
			{
				if(!sel.is_selected(*aIter)){
					foundUnselected = true;
					break;
				}
			}
		}

		if(!foundUnselected)
			sel.select(v);
	}
}

////////////////////////////////////////////////////////////////////////
//	SelectInnerSelectionEdges
template <class TSelector>
void SelectInnerSelectionEdges(TSelector& sel)
{
	if(!sel.grid())
		return;
	
	Grid& grid = *sel.grid();
	
	grid.begin_marking();
	
//	we'll first collect all edges that we want to check
	vector<EdgeBase*> edges;
	vector<EdgeBase*> vAssEdges;
	
//	iterate over all levels
	for(size_t lvl = 0; lvl < sel.num_levels(); ++lvl)
	{
		for(VolumeIterator iter = sel.template begin<Volume>(lvl);
			iter != sel.template end<Volume>(lvl); ++iter)
		{
			CollectEdges(vAssEdges, grid, *iter);
			for(size_t i = 0; i < vAssEdges.size(); ++i){
				EdgeBase* e = vAssEdges[i];
				if(!grid.is_marked(e)){
					grid.mark(e);
					edges.push_back(e);
				}
			}
		}

		for(FaceIterator iter = sel.template begin<Face>(lvl);
			iter != sel.template end<Face>(lvl); ++iter)
		{
			CollectEdges(vAssEdges, grid, *iter);
			for(size_t i = 0; i < vAssEdges.size(); ++i){
				EdgeBase* e = vAssEdges[i];
				if(!grid.is_marked(e)){
					grid.mark(e);
					edges.push_back(e);
				}
			}
		}
	}
	grid.end_marking();


//	now check for each edge if an unselected element is associated
	vector<Face*> vAssFaces;
	vector<Volume*> vAssVols;
	
	for(size_t i = 0; i < edges.size(); ++i)
	{
		EdgeBase* e = edges[i];
	
	//	check whether all associated elements are selected
		bool foundUnselected = false;
		
	//	volumes
		CollectVolumes(vAssVols, grid, e);
		for(size_t j = 0; j < vAssVols.size(); ++j)
		{
			if(!sel.is_selected(vAssVols[j])){
				foundUnselected = true;
				break;
			}
		}

	//	face		
		if(!foundUnselected){
			CollectFaces(vAssFaces, grid, e);
			for(size_t j = 0; j < vAssFaces.size(); ++j)
			{
				if(!sel.is_selected(vAssFaces[j])){
					foundUnselected = true;
					break;
				}
			}
		}

		if(!foundUnselected)
			sel.select(e);
	}
}


////////////////////////////////////////////////////////////////////////
//	SelectInnerSelectionFaces
template <class TSelector>
void SelectInnerSelectionFaces(TSelector& sel)
{
	if(!sel.grid())
		return;
	
	Grid& grid = *sel.grid();
	
	grid.begin_marking();
	
//	iterate through selected volumes and check for each side
//	whether it is connected to any unselected volumes.
	vector<Face*> vAssFaces;
	vector<Volume*> vAssVols;
	
//	iterate over all levels
	for(size_t lvl = 0; lvl < sel.num_levels(); ++lvl)
	{
		for(VolumeIterator iter = sel.template begin<Volume>(lvl);
			iter != sel.template end<Volume>(lvl); ++iter)
		{
			CollectFaces(vAssFaces, grid, *iter);
			for(size_t i = 0; i < vAssFaces.size(); ++i){
				Face* f = vAssFaces[i];
				if(!grid.is_marked(f)){
					grid.mark(f);
					CollectVolumes(vAssVols, grid, f);
					bool foundUnselected = false;
					for(size_t j = 0; j < vAssVols.size(); ++j){
						if(!sel.is_selected(vAssVols[j])){
							foundUnselected = true;
							break;
						}
					}
					
					if(!foundUnselected)
						sel.select(f);
				}
			}
		}
	}
	
	grid.end_marking();
}




////////////////////////////////////////////////////////////////////////
//	DeselectBoundarySelectionVertices
template <class TSelector>
void DeselectBoundarySelectionVertices(TSelector& sel)
{
	if(!sel.grid())
		return;
	
	Grid& grid = *sel.grid();

//	check each selected vertex of each level
	for(size_t lvl = 0; lvl < sel.num_levels(); ++lvl)
	{
		for(VertexBaseIterator iter = sel.template begin<VertexBase>(lvl);
			iter != sel.template end<VertexBase>(lvl);)
		{
			VertexBase* v = *iter;
		//	increase iterator here, since v may get deselected.
			++iter;
			
		//	check whether there is an unselected associated element
			bool foundUnselected = false;
			
		//	volumes
			for(Grid::AssociatedVolumeIterator aIter = grid.associated_volumes_begin(v);
				aIter != grid.associated_volumes_end(v); ++ aIter)
			{
				if(!sel.is_selected(*aIter)){
					foundUnselected = true;
					break;
				}
			}

		//	face		
			if(!foundUnselected){
				for(Grid::AssociatedFaceIterator aIter = grid.associated_faces_begin(v);
					aIter != grid.associated_faces_end(v); ++ aIter)
				{
					if(!sel.is_selected(*aIter)){
						foundUnselected = true;
						break;
					}
				}
			}
		
		//	edge		
			if(!foundUnselected){
				for(Grid::AssociatedEdgeIterator aIter = grid.associated_edges_begin(v);
					aIter != grid.associated_edges_end(v); ++ aIter)
				{
					if(!sel.is_selected(*aIter)){
						foundUnselected = true;
						break;
					}
				}
			}

			if(foundUnselected)
				sel.deselect(v);
		}
	}
}

////////////////////////////////////////////////////////////////////////
//	DeselectBoundarySelectionEdges
template <class TSelector>
void DeselectBoundarySelectionEdges(TSelector& sel)
{
	if(!sel.grid())
		return;
	
	Grid& grid = *sel.grid();

	vector<Face*> vAssFaces;
	vector<Volume*> vAssVols;

//	check each selected vertex of each level
	for(size_t lvl = 0; lvl < sel.num_levels(); ++lvl)
	{
		for(EdgeBaseIterator iter = sel.template begin<EdgeBase>(lvl);
			iter != sel.template end<EdgeBase>(lvl);)
		{
			EdgeBase* e = *iter;
		//	increase iterator here, since e may get deselected.
			++iter;
			
		//	check whether there is an unselected associated element
			bool foundUnselected = false;
			
		//	volumes
			CollectVolumes(vAssVols, grid, e);
			for(size_t j = 0; j < vAssVols.size(); ++j)
			{
				if(!sel.is_selected(vAssVols[j])){
					foundUnselected = true;
					break;
				}
			}

		//	face		
			if(!foundUnselected){
				CollectFaces(vAssFaces, grid, e);
				for(size_t j = 0; j < vAssFaces.size(); ++j)
				{
					if(!sel.is_selected(vAssFaces[j])){
						foundUnselected = true;
						break;
					}
				}
			}
			
			if(foundUnselected)
				sel.deselect(e);
		}
	}
}


////////////////////////////////////////////////////////////////////////
//	DeselectBoundarySelectionFaces
template <class TSelector>
void DeselectBoundarySelectionFaces(TSelector& sel)
{
	if(!sel.grid())
		return;
	
	Grid& grid = *sel.grid();
	
	vector<Volume*> vAssVols;
	
//	iterate over all levels
	for(size_t lvl = 0; lvl < sel.num_levels(); ++lvl)
	{
		for(FaceIterator iter = sel.template begin<Face>(lvl);
			iter != sel.template end<Face>(lvl);)
		{
			Face* f = *iter;
		//	increase iterator here, since f may get deselected.
			++iter;
			
			CollectVolumes(vAssVols, grid, f);
			for(size_t i = 0; i < vAssVols.size(); ++i){
				if(!sel.is_selected(vAssVols[i])){
					sel.deselect(f);
					break;
				}
			}
		}
	}
}


////////////////////////////////////////////////////////////////////////
// SelectLinkedFlatFaces
void SelectLinkedFlatFaces(Selector& sel, number maxDeviationAngle,
						   bool traverseFlipped, bool stopAtSelectedEdges,
						   APosition& aPos)
{
	if(!sel.grid())
		return;
	
	Grid& grid = *sel.grid();
	
	if(!grid.has_vertex_attachment(aPosition))
		return;
		
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);
	
//	convert the thresholdDegree to thresholdDot
	number thresholdDot = cos(deg_to_rad(maxDeviationAngle));
	
//	all initially selected faces are candidates
	queue<Face*> qCandidates;
	for(FaceIterator iter = sel.begin<Face>(); iter != sel.end<Face>(); ++iter)
		qCandidates.push(*iter);
	
//	we'll collect neighbours in this vector
	vector<Face*> vNeighbours;
	vector<EdgeBase*> edges;
	
//	while there are candidates left
	while(!qCandidates.empty())
	{
		Face* f = qCandidates.front();
		qCandidates.pop();
		
	//	calculate the normal
		vector3 n;
		CalculateNormal(n, f, aaPos);
		
		//CollectNeighbors(vNeighbours, f, grid);
	//	collect associated edges
		CollectAssociated(edges, grid, f);
		
	//	iterate through all neighbours
		for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge)
		{
			EdgeBase* e = edges[i_edge];
			if(stopAtSelectedEdges && sel.is_selected(e))
				continue;

			CollectAssociated(vNeighbours, grid, e);
			for(size_t i = 0; i < vNeighbours.size(); ++i){
				Face* nbr = vNeighbours[i];
				if(nbr == f)
					continue;
				
				if(!sel.is_selected(nbr)){
				//	compare normals
					vector3 nNbr;
					CalculateNormal(nNbr, nbr, aaPos);

				//	check dots
					number d = VecDot(n, nNbr);
					if(traverseFlipped)
						d = fabs(d);

					if(d >= thresholdDot){
					//	nbr is a linked flat face
						sel.select(nbr);
						qCandidates.push(nbr);
					}
				}
			}
		}
	}
}


void SelectLinkedFlatAndDegeneratedFaces(Selector& sel,
										 number maxDeviationAngle,
										 bool traverseFlipped,
										 bool stopAtSelectedEdges,
										 number degThreshold,
						   	   	   	     APosition& aPos)
{
	if(!sel.grid())
		return;

	Grid& grid = *sel.grid();

	if(!grid.has_vertex_attachment(aPosition))
		return;

	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);

//	convert the thresholdDegree to thresholdDot
	number thresholdDot = cos(deg_to_rad(maxDeviationAngle));
	number degThresholdSq = degThreshold * degThreshold;

//	all initially selected faces are candidates
//	with each candidate we'll also store its normal (this is
//	required since we'll traverse degenerated faces)
	queue<Face*> 	qCandidates;
	queue<vector3>	qNormals;
	for(FaceIterator iter = sel.begin<Face>(); iter != sel.end<Face>(); ++iter){
		qCandidates.push(*iter);
	//	calculate the normal
		vector3 n;
		CalculateNormal(n, *iter, aaPos);
		qNormals.push(n);
	}

//	temporary vectors for edges and faces
	vector<EdgeBase*> edges;
	vector<Face*> faces;

//	while there are candidates left
	while(!qCandidates.empty())
	{
		Face* f = qCandidates.front();
		qCandidates.pop();

	//	calculate the normal
		vector3 n = qNormals.front();
		qNormals.pop();

	//	get associated edges
		CollectAssociated(edges, grid, f);

	//	iterate through all edges
		for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge)
		{
			EdgeBase* e = edges[i_edge];

		//	only proceed if the edge is not degenerated
			if(EdgeLengthSq(e, aaPos) <= degThresholdSq)
				continue;

		//	if the edge is selected we might have to ignore it
			if(stopAtSelectedEdges && sel.is_selected(e))
				continue;

		//	check associated faces
			CollectAssociated(faces, grid, e);
			for(size_t i_face = 0; i_face < faces.size(); ++i_face){
				Face* nbr = faces[i_face];
				if(!sel.is_selected(nbr)){
				//	if the neighbor is degenerated, we can immediately select it.
					if(IsDegenerated(nbr, aaPos, degThreshold)){
						sel.select(nbr);
						qCandidates.push(nbr);
					//	use the normal from the face from which the degenerated
					//	face was encountered
						qNormals.push(n);
					}
					else{
					//	compare normals
						vector3 nNbr;
						CalculateNormal(nNbr, nbr, aaPos);

					//	check dots
						number d = VecDot(n, nNbr);
						if(traverseFlipped)
							d = fabs(d);

						if(d >= thresholdDot){
						//	nbr is a linked flat face
							sel.select(nbr);
							qCandidates.push(nbr);
							qNormals.push(nNbr);
						}
					}
				}
			}
		}
	}
}

}//	end of namespace
