//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m11 d17

#ifndef __H__LIB_GRID__DISTRIBUTION_UTIL__
#define __H__LIB_GRID__DISTRIBUTION_UTIL__

#include "lib_grid/lg_base.h"
#include "distribution_node_layout.h"
#include "grid_distribution.h"
#include "parallel_grid_layout.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	typedefs
typedef DistributionNodeLayout<VertexBase*>	DistributionVertexLayout;
typedef DistributionNodeLayout<EdgeBase*>	DistributionEdgeLayout;
typedef DistributionNodeLayout<Face*>		DistributionFaceLayout;
typedef DistributionNodeLayout<Volume*>		DistributionVolumeLayout;

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
void CreateGridLayouts(	std::vector<DistributionVertexLayout>& vertexLayoutsOut,
						std::vector<DistributionEdgeLayout>& edgeLayoutsOut,
						std::vector<DistributionFaceLayout>& faceLayoutsOut,
						std::vector<DistributionVolumeLayout>& volumeLayoutsOut,
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
							DistributionVertexLayout& vrtLayout,
							DistributionEdgeLayout& edgeLayout,
							DistributionFaceLayout& faceLayout,
							DistributionVolumeLayout& volLayout,
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
							GridLayoutMap& gridLayoutOut,
							std::istream& in);

////////////////////////////////////////////////////////////////////////
//	DeserializeLayoutInterfaces
/**
 * TLayoutMap has to be compatible with an
 * std::map<int, ParallelELEMENTLayout>, where ELEMENT can be either
 * Vertex, Edge, Face or Volume.
 */
template <class TGeomObj, class TLayoutMap>
void DeserializeLayoutInterfaces(TLayoutMap& layoutMapOut,
								std::vector<TGeomObj*> vGeomObjs,
								std::istream& in);
}//	end of namespace

////////////////////////////////
//	include implementation
#include "distribution_util_impl.hpp"

#endif
