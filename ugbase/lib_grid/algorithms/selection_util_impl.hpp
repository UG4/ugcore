// created by Sebastian Reiter
// y10 m05 d18
// s.b.reiter@googlemail.com

#ifndef __H__LIB_GRID__SELECTION_UTIL_IMPL__
#define __H__LIB_GRID__SELECTION_UTIL_IMPL__

namespace ug
{

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	selection util methods

////////////////////////////////////////////////////////////////////////
//	SelectAssociatedVertices
template <class TSelector, class TElemIterator>
void SelectAssociatedVertices(TSelector& sel, TElemIterator elemsBegin,
								TElemIterator elemsEnd)
{
	while(elemsBegin != elemsEnd)
	{
		uint numVrts = (*elemsBegin)->num_vertices();
		for(uint i = 0; i < numVrts; ++i)
			sel.select((*elemsBegin)->vertex(i));
		elemsBegin++;
	}
}

////////////////////////////////////////////////////////////////////////
//	SelectAssociatedEdges
template <class TSelector, class TElemIterator>
void SelectAssociatedEdges(TSelector& sel,
								TElemIterator elemsBegin,
								TElemIterator elemsEnd)
{
	Grid* pGrid = sel.get_assigned_grid();
	if(pGrid)
	{
		Grid& grid = *pGrid;
		std::vector<EdgeBase*> vEdges;
		while(elemsBegin != elemsEnd)
		{
			CollectEdges(vEdges, grid, *elemsBegin);
			for(uint i = 0; i < vEdges.size(); ++i)
				sel.select(vEdges[i]);
			elemsBegin++;
		}
	}
}

////////////////////////////////////////////////////////////////////////
//	SelectAssociatedFaces
template <class TSelector, class TElemIterator>
void SelectAssociatedFaces(TSelector& sel,
								TElemIterator elemsBegin,
								TElemIterator elemsEnd)
{
	Grid* pGrid = sel.get_assigned_grid();
	if(pGrid)
	{
		Grid& grid = *pGrid;
		std::vector<Face*> vFaces;
		while(elemsBegin != elemsEnd)
		{
			CollectFaces(vFaces, grid, *elemsBegin);
			for(uint i = 0; i < vFaces.size(); ++i)
				sel.select(vFaces[i]);
			elemsBegin++;
		}
	}
}

}// end of namespace

#endif
