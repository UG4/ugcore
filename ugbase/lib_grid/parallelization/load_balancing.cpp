// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 09.06.2011 (m,d,y)
 
#include "load_balancing.h"
#include "lib_grid/algorithms/graph/dual_graph.h"

//	if we're using metis, then include the header now
#ifdef UG_METIS
	extern "C" {
		#include "metis.h"
	}
#endif


using namespace std;


namespace ug{

////////////////////////////////////////////////////////////////////////////////
template <class TGeomBaseObj>
bool PartitionGrid_MetisKway(SubsetHandler& shPartitionOut,
							 Grid& grid, int numParts)
{
	#ifdef UG_METIS
		if(numParts == 1){
			shPartitionOut.assign_subset(grid.begin<TGeomBaseObj>(),
										grid.end<TGeomBaseObj>(), 0);
			return true;
		}

	//	construct the dual graph to the grid
		vector<idx_t> adjacencyMapStructure;
		vector<idx_t> adjacencyMap;
		ConstructDualGraph<TGeomBaseObj, idx_t>(adjacencyMapStructure,
												  adjacencyMap, grid);

	//	partition the graph using metis
		int n = (int)adjacencyMapStructure.size() - 1;
		int wgtflag = 0;
		int numflag = 0;
		int options[5]; options[0] = 0;
		int edgeCut;
		vector<idx_t> partitionMap(n);

		UG_LOG("CALLING METIS\n");
		METIS_PartGraphKway(&n, &adjacencyMapStructure.front(),
							&adjacencyMap.front(), NULL, NULL, &wgtflag,
							&numflag, &numParts, options, &edgeCut,
							&partitionMap.front());
		UG_LOG("METIS IS DONE\n");

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

//////////////////////////////
//	explicit instantiation
template bool PartitionGrid_MetisKway<EdgeBase>(SubsetHandler&, Grid&, int);
template bool PartitionGrid_MetisKway<Face>(SubsetHandler&, Grid&, int);
template bool PartitionGrid_MetisKway<Volume>(SubsetHandler&, Grid&, int);


////////////////////////////////////////////////////////////////////////////////
template <class TGeomBaseObj>
bool PartitionMultiGrid_MetisKway(SubsetHandler& shPartitionOut,
							 	  MultiGrid& mg, int numParts,
							 	  size_t baseLevel,
							 	  int hWeight, int vWeight)
{
#ifdef UG_METIS
	typedef TGeomBaseObj	TElem;
	typedef typename geometry_traits<TGeomBaseObj>::iterator	ElemIter;

//	only call metis if more than 1 part is required
	if(numParts > 1){
	//	here we'll store the dual graph
		vector<idx_t> adjacencyMapStructure;
		vector<idx_t> adjacencyMap;
		vector<idx_t> edgeWeightMap;

	//todo	add baseLevel to ConstructDualGraphMG.
		ConstructDualGraphMG<TGeomBaseObj, idx_t>(adjacencyMapStructure,
													adjacencyMap, &edgeWeightMap,
													mg, baseLevel, hWeight, vWeight);

	//	partition the graph using metis
		int n = (int)adjacencyMapStructure.size() - 1;
		int wgtflag = 1;//weights on edges only
		int numflag = 0;
		int options[5]; options[0] = 0;
		int edgeCut;
		vector<idx_t> partitionMap(n);

		UG_LOG("CALLING METIS\n");
		METIS_PartGraphKway(&n, &adjacencyMapStructure.front(),
							&adjacencyMap.front(), NULL,
							&edgeWeightMap.front(), &wgtflag,
							&numflag, &numParts, options, &edgeCut,
							&partitionMap.front());
		UG_LOG("METIS DONE\n");

	//	assign the subsets to the subset-handler
		int counter = 0;
		for(size_t lvl = baseLevel; lvl < mg.num_levels(); ++lvl){
			typedef typename geometry_traits<TGeomBaseObj>::iterator iterator;
			for(iterator iter = mg.begin<TGeomBaseObj>(lvl);
				iter != mg.end<TGeomBaseObj>(lvl); ++iter)
			{
				shPartitionOut.assign_subset(*iter, partitionMap[counter++]);
			}
		}
	}
	else{
	//	assign all elements to subset 0.
		for(size_t lvl = baseLevel; lvl < mg.num_levels(); ++lvl){
			shPartitionOut.assign_subset(mg.begin<TGeomBaseObj>(lvl),
										 mg.end<TGeomBaseObj>(lvl), 0);
		}
	}
	return true;
#else
	UG_LOG("WARNING in PartitionMultiGrid_MetisKway: METIS is not available in "
			"the current build. Please consider recompiling with METIS "
			"support enabled.\n");
	return false;
#endif
}

//////////////////////////////
//	explicit instantiation
template bool PartitionMultiGrid_MetisKway<EdgeBase>(SubsetHandler&, MultiGrid&,
													 int, size_t, int, int);
template bool PartitionMultiGrid_MetisKway<Face>(SubsetHandler&, MultiGrid&,
												 int, size_t, int, int);
template bool PartitionMultiGrid_MetisKway<Volume>(SubsetHandler&, MultiGrid&,
												   int, size_t, int, int);

}// end of namespace
