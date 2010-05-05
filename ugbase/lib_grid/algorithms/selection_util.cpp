// created by Sebastian Reiter
// y09 m11 d16
// s.b.reiter@googlemail.com

#include <vector>
#include <stack>
#include "selection_util.h"
#include "geom_obj_util/vertex_util.h"

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
							APosition& aPos)
{
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
		m_candidates.push((*iter)->vertex(0));
		m_candidates.push((*iter)->vertex(1));
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
		
	//	follow the smooth path
		while(srcVrt)
		{
		//	check smoothness for each connected unselected edge
			EdgeBase* bestEdge = NULL;
			number bestDot = -1.1;
			vector3 bestDir(0, 0, 0);
			
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
					if(d > bestDot && d > thresholdDot){
						bestEdge = e;
						bestDot = d;
						bestDir = dir;
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
				lastEdge = bestEdge;
				lastDir = bestDir;
			}
			else{
				srcVrt = NULL;
			}
		}
	}
}

}//	end of namespace
