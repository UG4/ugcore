// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 21.04.2011 (m,d,y)

#ifndef __H__UG__load_balancing__
#define __H__UG__load_balancing__

#include <vector>
#include "lib_grid/lg_base.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "util/parallel_callbacks.h"
#include "lib_grid/algorithms/graph/dual_graph.h"

//	if we're using metis, then include the header now
#ifdef UG_METIS
	extern "C" {
		#include "metis.h"
	}
#endif


namespace ug
{

////////////////////////////////////////////////////////////////////////
//	PartitionElementsByRepeatedIntersection
///	partitions a grid into subsets that contain the same number of elements (as far as possible).
/**
 * assigns elements to subsets.
 * Tries to construct subsets that all have the same size.
 * Runs with something like O(n*log(n)).
 * It can not be guaranteed, that subsets are connected.
 */
template <class TElem, int IDimension, class TAPosition>
bool PartitionElementsByRepeatedIntersection(SubsetHandler& shOut,
										Grid& grid,
										int numSubsets,
										TAPosition& aVrtPos,
										int startDim = 0);

////////////////////////////////////////////////////////////////////////
//	PartitionElementsByRepeatedIntersection
///	partitions a grid into subsets that contain the same number of elements (as far as possible).
/**
 * assigns elements to subsets.
 * Tries to construct subsets that all have the same size.
 * Runs with something like O(n*log(n)).
 * It can not be guaranteed, that subsets are connected.
 */
template <class TElem, int IDimension, class TAPosition>
bool PartitionElementsByRepeatedIntersection(SubsetHandler& shOut,
										MultiGrid& mg,
										int level,
										int numSubsets,
										TAPosition& aVrtPos,
										int startDim = 0);

////////////////////////////////////////////////////////////////////////////////
///	Partitions the elements in the grid by sorting them into a regular grid.
/**	Let xInd, yInd be the indices of the cell in which an element lies.
 * the associated subset index is then calculated by
 * \code
 * 		subsetIndex = yInd * numCellsX + xInd;
 * \endcode
 *
 * \param bucketSubset	All elements which shall not be considered are assigned
 * 						to this subset (default -1).
 */
template <class TElem, class TIterator, class TAAPos>
bool PartitionElements_RegularGrid(SubsetHandler& shOut,
								TIterator begin, TIterator end,
								int numCellsX, int numCellsY,
								TAAPos& aaPos,
								boost::function<bool (TElem*)> cbConsiderElem
									= (bool(*)(TElem*))ConsiderAll,
								int bucketSubset = -1);

////////////////////////////////////////////////////////////////////////////////
///	Partitions the elements in the grid using the METIS library
/**	This method calls METIS_PartGraphKway. Note that METIS is an external library
 * developed at Karypis Labs (http://glaros.dtc.umn.edu/gkhome/)
 *
 * Note that this method is best suited for partitions with more than 8 procs.
 * For less than 8 procs metis features other, better suited methods.
 */
template <class TGeomBaseObj>
bool PartitionGrid_MetisKway(SubsetHandler& shPartitionOut,
							 Grid& grid, int numParts)
{
	#ifdef UG_METIS
		using namespace std;

		if(numParts == 1){
			shPartitionOut.assign_subset(grid.begin<TGeomBaseObj>(),
										grid.end<TGeomBaseObj>(), 0);
			return true;
		}

	//	construct the dual graph to the grid
		vector<idxtype> adjacencyMapStructure;
		vector<idxtype> adjacencyMap;
		ConstructDualGraph<TGeomBaseObj, idxtype>(adjacencyMapStructure, adjacencyMap, grid);

	//	partition the graph using metis
		int n = (int)adjacencyMapStructure.size() - 1;
		int wgtflag = 0;
		int numflag = 0;
		int options[5]; options[0] = 0;
		int edgeCut;
		vector<idxtype> partitionMap(n);

		METIS_PartGraphKway(&n, &adjacencyMapStructure.front(), &adjacencyMap.front(),
							NULL, NULL, &wgtflag, &numflag, &numParts, options, &edgeCut, &partitionMap.front());

	//	assign the subsets to the subset-handler
		int counter = 0;
		typedef typename geometry_traits<TGeomBaseObj>::iterator iterator;
		for(iterator iter = grid.begin<TGeomBaseObj>();
			iter != grid.end<TGeomBaseObj>(); ++iter)
		{
			shPartitionOut.assign_subset(*iter, partitionMap[counter++]);
		}

		return true;
	#else
		UG_LOG("WARNING in PartitionGrid_MetisKway: METIS is not available in "
				"the current build. Please consider recompiling with METIS "
				"support enabled.\n");
		return false;
	#endif
}
}//	end of namespace


////////////////////////////////
//	include implementation
#include "load_balancing_impl.hpp"

#endif
