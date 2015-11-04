#ifndef __H__UG__hanging_node_refiner_t__
#define __H__UG__hanging_node_refiner_t__

#include "hanging_node_refiner_grid.h"
#include "hanging_node_refiner_multi_grid.h"

namespace ug
{

///	Gives access to a hanging node refiner, depending on the grid-type
template <class TGrid>
class THangingNodeRefiner;

template <>
class THangingNodeRefiner<Grid> : public HangingNodeRefiner_Grid
{
	public:
		THangingNodeRefiner(IRefinementCallback* refCallback = NULL) :
			HangingNodeRefiner_Grid(refCallback)	{}

		THangingNodeRefiner(Grid& grid,
							IRefinementCallback* refCallback = NULL) :
			HangingNodeRefiner_Grid(grid, refCallback)	{}
};

template <>
class THangingNodeRefiner<MultiGrid> : public HangingNodeRefiner_MultiGrid
{
	public:
		THangingNodeRefiner(IRefinementCallback* refCallback = NULL) :
			HangingNodeRefiner_MultiGrid(refCallback)	{}

		THangingNodeRefiner(MultiGrid& mg,
							IRefinementCallback* refCallback = NULL) :
			HangingNodeRefiner_MultiGrid(mg, refCallback)	{}
};

}//	end of namespace

#endif
