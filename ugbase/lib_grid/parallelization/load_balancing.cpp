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
		vector<idxtype> adjacencyMapStructure;
		vector<idxtype> adjacencyMap;
		ConstructDualGraph<TGeomBaseObj, idxtype>(adjacencyMapStructure,
												  adjacencyMap, grid);

	//	partition the graph using metis
		int n = (int)adjacencyMapStructure.size() - 1;
		int wgtflag = 0;
		int numflag = 0;
		int options[5]; options[0] = 0;
		int edgeCut;
		vector<idxtype> partitionMap(n);

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
							 	  int hWeight = 1, int vWeight = 1,
							 	  int baseLevel = 0)
{
#ifdef UG_METIS
	typedef TGeomBaseObj	TElem;
	typedef typename geometry_traits<TGeomBaseObj>::iterator	ElemIter;

//	here we'll store the dual graph
	vector<idxtype> adjacencyMapStructure;
	vector<idxtype> adjacencyMap;
	vector<idxtype> edgeWeightMap;

//todo	add baseLevel to ConstructDualGraphMG.
	ConstructDualGraphMG<TGeomBaseObj, idxtype>(adjacencyMapStructure,
												adjacencyMap, &edgeWeightMap,
												mg, hWeight, vWeight);

//	partition the graph using metis
	int n = (int)adjacencyMapStructure.size() - 1;
	int wgtflag = 1;//weights on edges only
	int numflag = 0;
	int options[5]; options[0] = 0;
	int edgeCut;
	vector<idxtype> partitionMap(n);

	UG_LOG("CALLING METIS\n");
	METIS_PartGraphKway(&n, &adjacencyMapStructure.front(),
						&adjacencyMap.front(), NULL,
						&edgeWeightMap.front(), &wgtflag,
						&numflag, &numParts, options, &edgeCut,
						&partitionMap.front());
	UG_LOG("METIS DONE\n");

//	assign the subsets to the subset-handler
	int counter = 0;
	typedef typename geometry_traits<TGeomBaseObj>::iterator iterator;
	for(iterator iter = mg.begin<TGeomBaseObj>();
		iter != mg.end<TGeomBaseObj>(); ++iter)
	{
		shPartitionOut.assign_subset(*iter, partitionMap[counter++]);
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
													 int, int, int, int);
template bool PartitionMultiGrid_MetisKway<Face>(SubsetHandler&, MultiGrid&,
												 int, int, int, int);
template bool PartitionMultiGrid_MetisKway<Volume>(SubsetHandler&, MultiGrid&,
												   int, int, int, int);

}// end of namespace
