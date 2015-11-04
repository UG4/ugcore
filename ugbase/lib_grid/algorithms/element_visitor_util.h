#ifndef __H__UG__element_visitor_util__
#define __H__UG__element_visitor_util__

#include <vector>
#include "grid/grid.h"
#include "common/util/metaprogramming_util.h"

namespace ug
{

///	Visits all elements between begin and end and executes the visitorCallback on them
template <class TIter>
void VisitAll(const TIter begin, const TIter end,
			  boost::function<void (typename TIter::value_type)> visitorCallback)
{
	TIter iter = begin;
	while(iter != end){
		typename TIter::value_type val = *iter;
		++iter;
		visitorCallback(val);
	}
}

///	Visits all boundary elements of the area specified through the iterators.
/** WARNING: This method uses Grid::mark
 * TIter::value_type has to be compatible to Edge*, Face* or Volume*.
 * Lets say that TIter::value_type equals TElem*.
 *
 * You have to specify a callback cbBelongsToArea with the signature
 * 'void (TElem*)'. which tells for each element of TIter::value_type, whether
 * it lies in the area or not.
 * Note that this callback should optimally return true for all elements between
 * begin and end, and false for all other elements.
 *
 * You furthermore have to specify a callback cbVisitSide, with the signature
 * 'void (TElem::side*)'. This callback is executed for all sides of the
 * given area.
 */
template <class TIter>
void VisitAreaBoundary(Grid& g, const TIter begin, const TIter end,
		      boost::function<bool (typename TIter::value_type)> cbBelongsToArea,
			  boost::function<void (typename Pointer2Value<typename TIter::value_type>::type::side)> cbVisitSide)
{
	typedef typename Pointer2Value<typename TIter::value_type>::type	TElem;
	typedef typename TElem::side										TSide;

	g.begin_marking();

	std::vector<TSide*> sides;
	std::vector<TElem*> elems;
	TIter iter = begin;
	while(iter != facesEnd){
		TElem* elem = *iter;
		++iter;
		CollectAssociated(sides, g, elem);
		for(size_t i = 0; i < sides.size(); ++i){
			TSide* side = sides[i];
			if(!g.is_marked(side)){
				g.mark(side);
			//	collect associated elems and check whether all belong to the area
				CollectAssociated(elems, g, side);
				int numInArea = 0;
				for(size_t i_elem = 0; i_elem < elems.size(); ++i_elem){
					if(cbBelongsToArea(elems[i_elem]))
						++numInArea;
				}

				if(numInArea == 1){
					cbVisitSide(side);
				}
			}
		}
	}

	grid.end_marking();
}

}//	end of namespace

#endif
