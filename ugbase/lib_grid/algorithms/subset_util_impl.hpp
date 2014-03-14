//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m07 d21

#ifndef __H__LIB_GRID__SUBSET_UTIL_IMPL__
#define __H__LIB_GRID__SUBSET_UTIL_IMPL__

#include "subset_util.h"
#include "selection_util.h"
#include "callback_util.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	FindFirstFreeSubset
template <class TElem>
int GetMaxSubsetIndex(SubsetHandler& sh)
{
//	go from back to front
	for(int i = (int)sh.num_subsets() - 1; i >= 0; --i)
	{
		if(sh.num_elements<TElem>(i) > 0)
		{
		//	this is the highest subset that contains elements of type TElem
			return i;
		}
	}

//	no subset contains elements of type TElem.
	return -1;
}

////////////////////////////////////////////////////////////////////////
//	MakeSubsetsConsecutive
template <class TElem>
void MakeSubsetsConsecutive(SubsetHandler& sh)
{
//	TODO: this algo could be slightly improved regarding runtime.

//	iterate through all subsets.
	for(int i = 0; i < sh.num_subsets(); ++i)
	{
	//	check whether the subset is empty
		if(sh.num_elements<TElem>(i) == 0)
		{
		//	it is. find the next filled one.
			for(int j = i + 1; j < sh.num_subsets(); ++j)
			{
				if(sh.num_elements<TElem>(j) > 0)
				{
				//	this subset is filled. move it to position i.
					sh.move_subset(j, i);
					break;
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////
//	EraseEmptySubsets
///	Erases all subsets which do not contain any geometric objects
void EraseEmptySubsets(ISubsetHandler& sh);


template <class TElem>
void SeparateSubsetsByLowerDimSubsets(Grid& grid, SubsetHandler& sh,
									  bool appendAtEnd)
{
	SeparateSubsetsByLowerDimSeparators<TElem>(grid, sh, appendAtEnd,
												IsNotInSubset(sh, -1));
}

template <class TElem>
void SeparateSubsetsByLowerDimSelection(Grid& grid, SubsetHandler& sh,
										Selector& sel, bool appendAtEnd)
{
	SeparateSubsetsByLowerDimSeparators<TElem>(grid, sh, appendAtEnd,
												IsSelected(sel));
}

template <class TElem>
void SeparateSubsetsByLowerDimSeparators(Grid& grid, SubsetHandler& sh,
					bool appendAtEnd,
					boost::function<bool (typename TElem::lower_dim_base_object*)>
						cbIsSeparator)

{
	using namespace std;

//	the element type of separating elements
	typedef typename TElem::lower_dim_base_object	TSide;

//	assign all elements to subset -1
	sh.assign_subset(grid.begin<TElem>(), grid.end<TElem>(), -1);

//	we'll keep all unassigned volumes in a selector.
	Selector sel(grid);
	sel.select(grid.begin<TElem>(), grid.end<TElem>());

//	those vectors will be used to gather element neighbours.
	vector<TSide*> vSides;
	vector<TElem*> vElems;

//	this stack contains all volumes that we still have to check for neighbours.
	stack<TElem*> stkElems;

//	now - while there are unassigned elements.
	int subsetIndex = 0;
	if(appendAtEnd)
		subsetIndex = sh.num_subsets();
	
	while(!sel.empty())
	{
	//	choose the element with which we want to start
	//	TODO: if material-points are supplied, this should be the
	//		the element that contains the i-th material point.
		stkElems.push(*sel.begin<TElem>());
		while(!stkElems.empty())
		{
			TElem* elem = stkElems.top();
			stkElems.pop();
		//	if the volume is unselected it has already been processed.
			if(!sel.is_selected(elem))
				continue;
			sel.deselect(elem);

		//	assign elem to its new subset
			sh.assign_subset(elem, subsetIndex);

		//	check neighbour-elements, whether they belong to the same subset.
		//	iterate through the sides of the element
			for(uint i = 0; i < elem->num_sides(); ++i)
			{
			//	get the i-th side
				TSide* side = grid.get_side(elem, i);

			//	check whether the side is regarded as a separator.
			//	If not, we'll add all associated elements.
				if(!cbIsSeparator(side))
				{
					CollectAssociated(vElems, grid, side);

				//	add all elements that are still selected (elem is not selected anymore).
					for(uint j = 0; j < vElems.size(); ++j)
					{
						if(sel.is_selected(vElems[j]))
							stkElems.push(vElems[j]);
					}
				}
			}
		}
	//	the stack is empty. increase subset index.
		subsetIndex++;
	}
}


////////////////////////////////////////////////////////////////////////
template <class TIterator>
void CopySubsetIndicesToSides(ISubsetHandler& sh, TIterator elemsBegin,
							TIterator elemsEnd, bool toUnassignedOnly)
{
	typedef typename PtrTypeToGeomObjBaseType<typename TIterator::value_type>::base_type TElem;
	typedef typename TElem::side TSide;

	if(!TElem::HAS_SIDES)
		return;

	UG_ASSERT(sh.grid(), "No grid assigned to subset-handler");

	Grid& grid = *sh.grid();

	typename Grid::traits<TSide>::secure_container	sides;
	for(TIterator iter = elemsBegin; iter != elemsEnd; ++iter){
		TElem* e = *iter;

		int si = sh.get_subset_index(e);

		grid.associated_elements(sides, e);
		for(size_t i = 0; i < sides.size(); ++i){
			TSide* s = sides[i];
			if(toUnassignedOnly){
				if(sh.get_subset_index(s) == -1)
					sh.assign_subset(s, si);
			}
			else
				sh.assign_subset(s, si);
		}
	}
}


////////////////////////////////////////////////////////////////////////
template <class TElem, class TSubsetHandler>
void AssignUnassignedElemsToSubset(TSubsetHandler& sh, int si)
{
	typedef typename geometry_traits<TElem>::iterator 	ElemIter;

//	access the grid on which sh operates.
	if(!sh.grid())
		return;

	Grid& grid = *sh.grid();

//	first make sure, that all elems are assigned to a subset, since
//	those won't be processed later on.

//	num is not part of ISubsetHandler and thus causes problems, if sh has type ISubsetHandler
	//if(sh.template num<TElem>() != grid.num<TElem>()){
		for(ElemIter iter = grid.begin<TElem>();
			iter != grid.end<TElem>(); ++iter)
		{
			if(sh.get_subset_index(*iter) == -1)
				sh.assign_subset(*iter, si);
		}
	//}
}

////////////////////////////////////////////////////////////////////////
template <class TSubsetHandler>
void AdjustSubsetsForSimulation(TSubsetHandler& sh,
								bool preserveExistingSubsets)
{
//	access the grid on which sh operates.
	if(!sh.grid())
		return;

	Grid& grid = *sh.grid();
	Selector sel(grid);

	if(grid.num_volumes() > 0){
	//	deselect all elements of lower dimension, if existing subsets are
	//	not to be preserved.
		if(!preserveExistingSubsets){
			sh.assign_subset(grid.begin<Face>(), grid.end<Face>(), -1);
			sh.assign_subset(grid.begin<Edge>(), grid.end<Edge>(), -1);
			sh.assign_subset(grid.begin<Vertex>(), grid.end<Vertex>(), -1);
		}

		AssignUnassignedElemsToSubset<Volume>(sh, GetFirstFreeSubset(sh));
		AssignSidesToSubsets<Volume>(sh);

		SelectInterfaceElements(sel, sh, grid.begin<Face>(), grid.end<Face>());
		SelectBoundaryElements(sel, grid.begin<Face>(), grid.end<Face>());
		SelectAssociatedEdges(sel, sel.begin<Face>(), sel.end<Face>());
		SelectAssociatedVertices(sel, sel.begin<Face>(), sel.end<Face>());

		AssignSidesToSubsets<Face>(sh, &sel);
		AssignSidesToSubsets<Edge>(sh, &sel);

	//	finally assign vertices on edge interfaces
		sel.clear<Edge>();
		sel.clear<Vertex>();
		SelectInterfaceElements(sel, sh, grid.begin<Edge>(),
								grid.end<Edge>(), true);
		sel.clear<Face>();
		SelectAssociatedVertices(sel, sel.begin<Edge>(), sel.end<Edge>());

		AssignSidesToSubsets<Edge>(sh, &sel);
	}
	else if(grid.num_faces() > 0){
	//	deselect all elements of lower dimension, if existing subsets are
	//	not to be preserved.
		if(!preserveExistingSubsets){
			sh.assign_subset(grid.begin<Edge>(), grid.end<Edge>(), -1);
			sh.assign_subset(grid.begin<Vertex>(), grid.end<Vertex>(), -1);
		}

		AssignUnassignedElemsToSubset<Face>(sh, GetFirstFreeSubset(sh));
		AssignSidesToSubsets<Face>(sh);

		SelectInterfaceElements(sel, sh, grid.begin<Edge>(), grid.end<Edge>());
		SelectBoundaryElements(sel, grid.begin<Edge>(), grid.end<Edge>());
		SelectAssociatedVertices(sel, sel.begin<Edge>(), sel.end<Edge>());

		AssignSidesToSubsets<Edge>(sh, &sel);
	}
	else if(grid.num_edges() > 0){
	//	deselect all elements of lower dimension, if existing subsets are
	//	not to be preserved.
		if(!preserveExistingSubsets){
			sh.assign_subset(grid.begin<Vertex>(), grid.end<Vertex>(), -1);
		}

		AssignUnassignedElemsToSubset<Edge>(sh, GetFirstFreeSubset(sh));
		AssignSidesToSubsets<Edge>(sh);
	}
	else{
		AssignUnassignedElemsToSubset<Vertex>(sh, GetFirstFreeSubset(sh));
	}

//	erase empty subsets again
	EraseEmptySubsets(sh);
}

////////////////////////////////////////////////////////////////////////
template <class TAAPosVRT>
number FaceArea(ISubsetHandler& sh, int si, size_t lvl, TAAPosVRT& aaPos)
{
	number sum = 0.;
	GridObjectCollection goc = sh.get_grid_objects_in_subset(si);

	if (goc.num<Face>(lvl) == 0) {
		UG_WARNING("WARNING: Given subset doesn't contain any faces on given level.");
	} else {
		typedef geometry_traits<Face>::const_iterator CIT;
		for (CIT cit = goc.faces_begin(lvl); cit != goc.faces_end(lvl); cit++)
			sum += FaceArea(*cit, aaPos);
	}

	return sum;
}


////////////////////////////////////////////////////////////////////////
//	AssignAssociatedVerticesToSubset
template <class TIterator>
void AssignAssociatedVerticesToSubset(ISubsetHandler& sh, TIterator elemsBegin,
										TIterator elemsEnd, int subsetIndex)
{
//	iterate through the elements
	for(;elemsBegin != elemsEnd; elemsBegin++)
	{
		typename TIterator::value_type elem = *elemsBegin;
		uint numVrts = elem->num_vertices();
	//	iterate through the vertices of elem and assign them
		for(uint i = 0; i < numVrts; ++i)
			sh.assign_subset(elem->vertex(i), subsetIndex);
	}
}

////////////////////////////////////////////////////////////////////////
template <class TElem, class TSubsetHandler>
void AssignAssociatedVerticesToSubsets(TSubsetHandler& sh,
									const ISubsetHandler& srcIndHandler)
{
	typedef typename geometry_traits<TElem>::const_iterator iterator;
	for(size_t l  = 0; l < sh.num_levels(); ++l){
		for(int si = 0; si < sh.num_subsets(); ++si){
			for(iterator iter = sh.template begin<TElem>(si, l);
				iter != sh.template end<TElem>(si, l); ++iter)
			{
				TElem* e = *iter;
				for(size_t i = 0; i < e->num_vertices(); ++i)
				{
					Vertex* vrt = e->vertex(i);
					sh.assign_subset(vrt, srcIndHandler.get_subset_index(vrt));
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////
template <class TElem, class TSubsetHandler>
void AssignAssociatedEdgesToSubsets(TSubsetHandler& sh,
									const ISubsetHandler& srcIndHandler)
{
	typedef typename geometry_traits<TElem>::const_iterator iterator;
	std::vector<Edge*> vEdges;

	for(size_t l  = 0; l < sh.num_levels(); ++l){
		for(int si = 0; si < sh.num_subsets(); ++si){
			for(iterator iter = sh.template begin<TElem>(si, l);
				iter != sh.template end<TElem>(si, l); ++iter)
			{
				TElem* e = *iter;
				CollectEdges(vEdges, *sh.grid(), e);

				for(size_t i = 0; i < vEdges.size(); ++i)
				{
					Edge* edge = vEdges[i];
					sh.assign_subset(edge, srcIndHandler.get_subset_index(edge));
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////
template <class TElem, class TSubsetHandler>
void AssignAssociatedFacesToSubsets(TSubsetHandler& sh,
									const ISubsetHandler& srcIndHandler)
{
	typedef typename geometry_traits<TElem>::const_iterator iterator;
	std::vector<Face*> vFaces;

	for(size_t l  = 0; l < sh.num_levels(); ++l){
		for(int si = 0; si < sh.num_subsets(); ++si){
			for(iterator iter = sh.template begin<TElem>(si, l);
				iter != sh.template end<TElem>(si, l); ++iter)
			{
				TElem* e = *iter;
				CollectFaces(vFaces, *sh.grid(), e);

				for(size_t i = 0; i < vFaces.size(); ++i)
				{
					Face* f = vFaces[i];
					sh.assign_subset(f, srcIndHandler.get_subset_index(f));
				}
			}
		}
	}
}

template <class TElem, class TSubsetHandlerDest, class TSubsetHandlerSrc>
void AssignAssociatedSidesToSubsets(TSubsetHandlerDest& sh,
									const TSubsetHandlerSrc& srcIndHandler)
{
	typedef typename geometry_traits<TElem>::const_iterator iterator;
	typedef typename TElem::lower_dim_base_object Side;
	std::vector<Side*> vSides;
	Grid& grid = *sh.grid();

	for(size_t l  = 0; l < sh.num_levels(); ++l){
		for(int si = 0; si < sh.num_subsets(); ++si){
			for(iterator iter = sh.template begin<TElem>(si, l);
				iter != sh.template end<TElem>(si, l); ++iter)
			{
				TElem* e = *iter;
				CollectAssociated(vSides, grid, e);

				for(size_t i = 0; i < vSides.size(); ++i)
				{
					Side* s = vSides[i];
					sh.assign_subset(s, srcIndHandler.get_subset_index(s));
				}
			}
		}
	}
}

///	helper with with dummy-param for compile-time function selection.
template <class TElem, class TSubsetHandlerDest, class TSubsetHandlerSrc>
void AssignAssociatedLowerDimElemsToSubsets(TSubsetHandlerDest& sh,
									const TSubsetHandlerSrc& srcIndHandler,
									const Volume&)
{
//	we have to find all associated elements of lower dimension.
	if(srcIndHandler.template num<Face>() > 0)
		AssignAssociatedFacesToSubsets<TElem>(sh, srcIndHandler);
	if(srcIndHandler.template num<Edge>() > 0)
		AssignAssociatedEdgesToSubsets<TElem>(sh, srcIndHandler);
	if(srcIndHandler.template num<Vertex>() > 0)
		AssignAssociatedVerticesToSubsets<TElem>(sh, srcIndHandler);
}

///	helper with with dummy-param for compile-time function selection.
template <class TElem, class TSubsetHandlerDest, class TSubsetHandlerSrc>
void AssignAssociatedLowerDimElemsToSubsets(TSubsetHandlerDest& sh,
									const TSubsetHandlerSrc& srcIndHandler,
									const Face&)
{
//	we have to find all associated elements of lower dimension.
	if(srcIndHandler.template num<Edge>() > 0)
		AssignAssociatedEdgesToSubsets<TElem>(sh, srcIndHandler);
	if(srcIndHandler.template num<Vertex>() > 0)
		AssignAssociatedVerticesToSubsets<TElem>(sh, srcIndHandler);
}

///	helper with with dummy-param for compile-time function selection.
template <class TElem, class TSubsetHandlerDest, class TSubsetHandlerSrc>
void AssignAssociatedLowerDimElemsToSubsets(TSubsetHandlerDest& sh,
									const TSubsetHandlerSrc& srcIndHandler,
									const Edge&)
{
//	we have to find all associated elements of lower dimension.
	if(srcIndHandler.template num<Vertex>() > 0)
		AssignAssociatedVerticesToSubsets<TElem>(sh, srcIndHandler);
}

////////////////////////////////////////////////////////////////////////
template <class TElem, class TSubsetHandlerDest, class TSubsetHandlerSrc>
void AssignAssociatedLowerDimElemsToSubsets(TSubsetHandlerDest& sh,
									const TSubsetHandlerSrc& srcIndHandler)
{
	AssignAssociatedLowerDimElemsToSubsets<TElem>(sh,
											srcIndHandler, TElem());
}

////////////////////////////////////////////////////////////////////////
template <typename TBaseObj>
void FindSubsetGroups
(
	std::vector<int> & minCondInd,
	const std::vector<bool> & isMarked,
	const ISubsetHandler & sh,
	const NeighborhoodType nbhType
)
{
	typedef typename geometry_traits<TBaseObj>::const_iterator elem_iterator;
	
	UG_ASSERT (((int) isMarked.size ()) == sh.num_subsets (), "FindSubsetGroups: array size mismatch");
	
	std::vector<TBaseObj*> neighbours;
	
//	Prepare minCondInd
	minCondInd.resize (sh.num_subsets ());
	for (size_t si = 0; si < minCondInd.size (); si++)
		minCondInd [si] = (isMarked [si])? si : -1;
	
//	Loop over the subsets:
	for (size_t si = 0; si < minCondInd.size (); si++)
	{
		int min_si;
		
	//	Marked subset?
		if ((min_si = minCondInd [si]) < 0)
			continue; // no, we do not treat this subset
		
	//	Yes, loop over the elements in the subdomain (in the grid level 0):
		GridObjectCollection goc = sh.get_grid_objects_in_subset (si);
		elem_iterator e_end = goc.end<TBaseObj> (0);
		bool is_empty = true;
		for (elem_iterator e_iter = goc.begin<TBaseObj> (0); e_iter != e_end; ++e_iter)
		{
			is_empty = false;
		//	Loop over the neighbours:
			CollectNeighbors (neighbours, *e_iter, *sh.grid(), nbhType);
			for (size_t k = 0; k < neighbours.size (); k++)
			{
				int min_nbr_si;
				int nbr_si = sh.get_subset_index (neighbours [k]);
				
				if (nbr_si < 0 || nbr_si >= (int) minCondInd.size ())
					UG_THROW ("FindSubsetGroups: Illegal neighbour subset index.");
				if ((min_nbr_si = minCondInd [nbr_si]) < 0)
					continue; // we do not treat this subset
				
			//	Set the same smallest index to both groups of the subsets:
				if (min_nbr_si < min_si)
				{
					for (size_t l = 0; l < minCondInd.size (); l++)
						if (minCondInd [l] == min_si)
							minCondInd [l] = min_nbr_si;
				}
				else if (min_nbr_si > min_si)
				{
					for (size_t l = 0; l < minCondInd.size (); l++)
						if (minCondInd [l] == min_nbr_si)
							minCondInd [l] = min_si;
				}
			}
		}
		if (is_empty) minCondInd [si] = -2;
	}
}

}//	end of namespace

#endif
