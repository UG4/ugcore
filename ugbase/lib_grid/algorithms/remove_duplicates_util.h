// created by Sebastian Reiter
// s.b.reiter@gmail.com

#ifndef __H__UG_remove_duplicates_util
#define __H__UG_remove_duplicates_util

#include "lib_grid/grid/grid.h"
#include "lib_grid/grid/grid_util.h"

namespace ug{

////////////////////////////////////////////////////////////////////////
///	Removes elements which share the same set of vertices with other elements in the given range
template <class TElemIter>
void RemoveDuplicates(Grid& grid, TElemIter elemsBegin, TElemIter elemsEnd)
{
	typedef typename TElemIter::value_type					TElemPtr;
	typedef typename PtrToValueType<TElemPtr>::base_type	TElem;

	grid.begin_marking();

//	the first is the double, the second the original
	std::vector<std::pair<TElem*, TElem*> > doubles;
	typename Grid::traits<TElem>::secure_container	elems;

	for(TElemIter iter = elemsBegin; iter != elemsEnd; ++iter){
		TElem* e = *iter;
		if(!grid.is_marked(e)){
			typename TElem::ConstVertexArray vrts = e->vertices();
			const size_t numVrts = e->num_vertices();

			bool allMarked = true;
			for(size_t i = 0; i < numVrts; ++i){
				if(!grid.is_marked(vrts[i])){
					allMarked = false;
					break;
				}
			}

			bool isDouble = false;

			if(allMarked){
			//	a necessary condition is met. However not yet sufficient.
			//	find marked elems wich connect the same vertices.
				grid.associated_elements(elems, vrts[0]);
				for(size_t i = 0; i < elems.size(); ++i){
					TElem* te = elems[i];
					if(grid.is_marked(te) && CompareVertices(e, te)){
					//	e is a double
						isDouble = true;
						doubles.push_back(std::make_pair(e, te));
						break;
					}
				}
			}

		//	finally mark e and its vertices (every processed non-double element is marked).
			if(!isDouble){
				grid.mark(e);
				for(size_t i = 0; i < numVrts; ++i)
					grid.mark(vrts[i]);
			}
		}
	}

	grid.end_marking();

//	now erase all doubles
	for(size_t i = 0; i < doubles.size(); ++i){
	//	this allows listeners to take actions
		grid.objects_will_be_merged(doubles[i].second, doubles[i].second,
									doubles[i].first);
		grid.erase(doubles[i].first);
	}
}

}//	end of namespace

#endif	//__H__UG_remove_duplicates_util
