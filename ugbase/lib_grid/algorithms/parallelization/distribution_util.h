//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m11 d17

#ifndef __H__LIB_GRID__DISTRIBUTION_UTIL__
#define __H__LIB_GRID__DISTRIBUTION_UTIL_

#include "parallel_node_layout.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	typedefs
typedef pcl::ParallelNodeLayout<VertexBase*>	ParallelVertexLayout;
typedef pcl::ParallelNodeLayout<EdgeBase*>		ParallelEdgeLayout;
typedef pcl::ParallelNodeLayout<Face*>			ParallelFaceLayout;
typedef pcl::ParallelNodeLayout<Volume*>		ParallelVolumeLayout;

struct ParallelGridLayout
{
	ParallelVertexLayout	vertexLayout;
	ParallelEdgeLayout		edgeLayout;
	ParallelFaceLayout		faceLayout;
	ParallelVolumeLayout	volumeLayout;
};

////////////////////////////////////////////////////////////////////////
//	DistributeGrid
///	distributes a grid to several processes.
/**
 * Partitions the grid into parts specified by the \sa SubsetHandler sh.
 * For partitioning the grid the algorithm \sa CreateGridLayouts
 * is used internally. Have a look at its documentation to see how you
 * have to specify your parts in the \sa SubsetHandler.
 * The part for the local (the calling) process (specified by localProcID)
 * won't be send through the network, instead it will be directly written
 * to pLocalGridOut and pLocalGridCommSetOut - if those are specified.
 * Both should be empty before calling this method. Passing NULL only makes
 * sense if you specified a processMap that has no entry for the calling
 * process.
 * You may optionally specify a process-map. This map is represented by
 * a std::vector<int> which should have as many entries as there are
 * subsets in the SubsetHandler. Each entry specifies the target process
 * for the associated subset. You may pass NULL as parameter. The grids
 * are distributed to processes 0 to sh.num_subsets()-1 in this case.
 *
 * Grids distributed through this method may be received by \sa ReveiveGrid.
 */
void DistributeGrid(MultiGrid& mg, SubsetHandler& sh, int localProcID,
					MultiGrid* pLocalGridOut = NULL,
					ParallelGridLayout* pLocalGridLayoutOut = NULL,
					std::vector<int>* pProcessMap = NULL);

////////////////////////////////////////////////////////////////////////
//	ReceiveGrid
///	receives a part of a grid that was distributed through \sa DistributeGrid.
/**
 * gridOut and gridLayoutOut should be empty when passed to this method.
 * srcProcID has to specify the process that distributes the grids.
 */
void ReceiveGrid(MultiGrid& mgOut, ParallelGridLayout& gridLayoutOut,
					int srcProcID);

////////////////////////////////////////////////////////////////////////
//	CreateGridLayouts
///	Creates parallel layouts for vertices, edges, faces and volumes
/**
 * Given a MultiGrid and a SubsetHandler, this method creates parallel
 * layouts for vertices, edges, ...
 * For each subset a new layout is created. That means that i.e.
 * vertexLayoutsOut[k] holds the vertex-layout for the k-th subset.
 *
 * If you pass a pointer to a valid selector (which is registered at mg),
 * the selector will be used for internal calculations. The only reason
 * for this parameter is a speed increase.
 * You shouldn't assume anything about the content of pSel after the
 * method finished.
 */
void CreateGridLayouts(	std::vector<ParallelVertexLayout>& vertexLayoutsOut,
						std::vector<ParallelEdgeLayout>& edgeLayoutsOut,
						std::vector<ParallelFaceLayout>& faceLayoutsOut,
						std::vector<ParallelVolumeLayout>& volumeLayoutsOut,
						MultiGrid& mg, SubsetHandler& sh,
						MGSelector* pSel = NULL);
						
////////////////////////////////////////////////////////////////////////
//	SerializeGridAndLayouts
/**
 * Writes the elements of a grid, which are referenced by the given
 * gridLayout to a binary stream.
 * You may pass a selector via pSel, which increases performance of this
 * method. After this method finished pSel will contain all elements
 * that have been written to the stream (in the same order as they were
 * written).
 */
void SerializeGridAndLayouts(std::ostream& out, MultiGrid& mg,
							ParallelVertexLayout& vrtLayout,
							ParallelEdgeLayout& edgeLayout,
							ParallelFaceLayout& faceLayout,
							ParallelVolumeLayout& volLayout,
							AInt& aLocalIndVRT, AInt& aLocalIndEDGE,
							AInt& aLocalIndFACE, AInt& aLocalIndVOL,
							MGSelector* pSel = NULL,
							std::vector<int>* pProcessMap = NULL);

////////////////////////////////////////////////////////////////////////
//	SerializeLayoutInterfaces
template <class TLayout>
void SerializeLayoutInterfaces(std::ostream& out, TLayout& layout,
								std::vector<int>* pProcessMap = NULL);

////////////////////////////////////////////////////////////////////////
//	DeserializeGridAndLayouts
void DeserializeGridAndLayouts(MultiGrid& mgOut,
							ParallelGridLayout& gridLayoutOut,
							std::istream& in);

////////////////////////////////////////////////////////////////////////
//	DeserializeLayoutInterfaces
template <class TLayout>
void DeserializeLayoutInterfaces(TLayout& layoutOut,
								std::istream& in);
}//	end of namespace

////////////////////////////////
//	include implementation
#include "distribution_util_impl.hpp"

#endif
