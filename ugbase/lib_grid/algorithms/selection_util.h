// created by Sebastian Reiter
// y09 m11 d06
// s.b.reiter@googlemail.com

#ifndef __H__LIB_GRID__SELECTION_UTIL__
#define __H__LIB_GRID__SELECTION_UTIL__

namespace ug
{

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	selection util methods

////////////////////////////////////////////////////////////////////////
//	SelectAssociatedVertices
///	selects all associated vertices of the elements between elemsBegin and elemsEnd
/**
 * TElemIterator has to be a stl-compatible iterator.
 * The underlying element-type has to be a pointer to a class that
 * features the following methods:
 * 
 * VertexBase* vertex(int i);//returns the i-th vertex of the element.
 * uint num_vertices();//returns the number of vertices that the element holds.
 *
 * Valid classes are for example EdgeBase, Face and Volume.
 *
 * Make sure that the elements only reference vertices that belong to the grid
 * at which the selector is registered.
 */
template <class TElemIterator>
void SelectAssociatedVertices(Selector& sel, TElemIterator elemsBegin,
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
///	selects all associated edges of the elements between elemsBegin and elemsEnd
/**
 * TElemIterator has to be a stl-compatible iterator.
 * The underlying element-type has to be a pointer to a class that
 * is supported by libGrid::CollectEdges(...)
 *
 * Valid classes are for example Face and Volume.
 *
 * Make sure that the elements only reference edges that belong to the grid
 * at which the selector is registered.
 */
template <class TElemIterator>
void SelectAssociatedEdges(Selector& sel,
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
///	selects all associated faces of the elements between elemsBegin and elemsEnd
/**
 * TElemIterator has to be a stl-compatible iterator.
 * The underlying element-type has to be a pointer to a class that
 * is supported by libGrid::CollectFaces(...)
 *
 * A valid classe is for example Volume.
 *
 * Make sure that the elements only reference faces that belong to the grid
 * at which the selector is registered.
 */
template <class TElemIterator>
void SelectAssociatedFaces(Selector& sel,
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
