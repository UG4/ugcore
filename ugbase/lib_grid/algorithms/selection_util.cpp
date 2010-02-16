// created by Sebastian Reiter
// y09 m11 d16
// s.b.reiter@googlemail.com

#include <vector>
#include "selection_util.h"

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

}//	end of namespace
