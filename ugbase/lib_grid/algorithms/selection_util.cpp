// created by Sebastian Reiter
// y09 m11 d16
// s.b.reiter@googlemail.com

#include <vector>
#include <stack>
#include "selection_util.h"
#include "geom_obj_util/vertex_util.h"
#include "geom_obj_util/edge_util.h"
using namespace std;

namespace ug
{
////////////////////////////////////////////////////////////////////////
void SelectAreaBoundaryEdges(ISelector& sel, FaceIterator facesBegin,
							  FaceIterator facesEnd)
{
	if(!sel.get_assigned_grid())
		return;
	
	Grid& grid = *sel.get_assigned_grid();

//	iterate over associated edges of faces.
//	mark edges if they have not yet been examined, unmark them if they
//	have been (-> inner edge). Use Grid::mark to avoid reselection
//	and to keep existing selections.
	grid.begin_marking();

	vector<EdgeBase*> vEdges;
	while(facesBegin != facesEnd){
		Face* f = *facesBegin;
		++facesBegin;
		CollectEdges(vEdges, grid, f);
		for(size_t i = 0; i < vEdges.size(); ++i){
			EdgeBase* e = vEdges[i];
			if(!grid.is_marked(e)){
			//	if the edge was initially selected, it should stay that way
				if(!sel.is_selected(e)){
					grid.mark(e);
					sel.select(e);
				}
			}
			else{
			//	if the edge is not selected, then it already is an inner edge
				sel.deselect(e);
			}
		}
	}

	grid.end_marking();
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
//	SelectAssociatedGenealogy
void SelectAssociatedGenealogy(MGSelector& msel, bool selectAssociatedElements)
{
	MultiGrid* mg = msel.get_assigned_multi_grid();
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
	
	if(!sel.get_assigned_grid())
		return;
	
	Grid& grid = *sel.get_assigned_grid();
	
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
		
		for(EdgeBaseIterator iter = grid.associated_edges_begin(srcVrt);
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
			EdgeBaseIterator iterEnd = grid.associated_edges_end(srcVrt);
			for(EdgeBaseIterator iter = grid.associated_edges_begin(srcVrt);
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
							
							if(bLastNormalValid && bNormalValid &! ignoreNormalChecks){
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
void SelectInnerSelectionVertices(Selector& sel)
{
	if(!sel.get_assigned_grid())
		return;
	
	Grid& grid = *sel.get_assigned_grid();
	
	grid.begin_marking();
	
//	we'll first collect all vertices that we want to check
	vector<VertexBase*> vrts;
	
	for(VolumeIterator iter = sel.begin<Volume>();
		iter != sel.end<Volume>(); ++iter)
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

	for(FaceIterator iter = sel.begin<Face>();
		iter != sel.end<Face>(); ++iter)
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

	for(EdgeBaseIterator iter = sel.begin<EdgeBase>();
		iter != sel.end<EdgeBase>(); ++iter)
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
		
	grid.end_marking();


//	now check for each vertex if an unselected element is associated
	for(size_t i = 0; i < vrts.size(); ++i)
	{
		VertexBase* v = vrts[i];
	
	//	check whether all associated elements are selected
		bool foundUnselected = false;
		
	//	volumes
		for(VolumeIterator aIter = grid.associated_volumes_begin(v);
			aIter != grid.associated_volumes_end(v); ++ aIter)
		{
			if(!sel.is_selected(*aIter)){
				foundUnselected = true;
				break;
			}
		}

	//	face		
		if(!foundUnselected){
			for(FaceIterator aIter = grid.associated_faces_begin(v);
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
			for(EdgeBaseIterator aIter = grid.associated_edges_begin(v);
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
void SelectInnerSelectionEdges(Selector& sel)
{
	if(!sel.get_assigned_grid())
		return;
	
	Grid& grid = *sel.get_assigned_grid();
	
	grid.begin_marking();
	
//	we'll first collect all edges that we want to check
	vector<EdgeBase*> edges;
	vector<EdgeBase*> vAssEdges;
	
	for(VolumeIterator iter = sel.begin<Volume>();
		iter != sel.end<Volume>(); ++iter)
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

	for(FaceIterator iter = sel.begin<Face>();
		iter != sel.end<Face>(); ++iter)
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
void SelectInnerSelectionFaces(Selector& sel)
{
	if(!sel.get_assigned_grid())
		return;
	
	Grid& grid = *sel.get_assigned_grid();
	
	grid.begin_marking();
	
//	iterate through selected volumes and check for each side
//	whether it is connected to any unselected volumes.
	vector<Face*> vAssFaces;
	vector<Volume*> vAssVols;
	
	for(VolumeIterator iter = sel.begin<Volume>();
		iter != sel.end<Volume>(); ++iter)
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

	grid.end_marking();
}

}//	end of namespace
