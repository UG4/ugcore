/**
 * \file agglomeration.h
 *
 * \author Martin Rupp
 *
 * \date 03.08.2011
 *
 * Goethe-Center for Scientific Computing 2011.
 */

#include <vector>
#include <map>

namespace ug
{

// Agglomeration
//---------------------------------------------------------------------------
/**
 * serial load distribution algorithm
 *
 * @param sizes 		the sizes of the processing units
 * @param connections 	connections[i][j] = connection size between unit i and unit j
 * @param mergeWith		(returned) : if mergeWith[i].size() == 1, the processor to send matrix to (if == pcl::GetProcRank, no mergin),
 * 						otherwise, list of nodes which are merged, with mergeWith[i][0] = i.
 * @param minimalSize	force all merged nodes to be greater than minimalSize
 * @param preferredSize
 *
 * Example:
 * minimalSize = 400
 * preferredSize = 4000
 * this setting ensures that if we have to agglomerate on one level, we try to do
 * as much agglomeration as possible, so that we have only 1 or 2 levels with many agglomerations and all other levels
 * are free of agglomeration (we want to prevent levels where only 2 processors agglomerate)
 *
 * mergeWith[0] = {0} : no merging, processor 1 stays
 * mergeWith[1] = {1} : no merging, processor 1 stays
 * mergeWith[2] = {5} : processor 2 send matrix to processor 5
 * mergeWith[3] = {5} : processor 3 send matrix to processor 5
 * mergeWith[4] = {5} : processor 4 send matrix to processor 5
 * mergeWith[5] = {5, 3, 4, 2}. : processor 5 receives matrix from processors 3, 4 and 2.
 *
 * @note this is a very easy serial agglomeration algorithm
 * @todo add agglomeration interface and other implementations like METIS
 */
void EasyAgglomeration(const std::vector<size_t> sizes,
		const std::vector<std::map<int, size_t> > connections,
		std::vector<std::vector<int> > &mergeWith,
		size_t minimalSize, size_t preferredSize);

}
