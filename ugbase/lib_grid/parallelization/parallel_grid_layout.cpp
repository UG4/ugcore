// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 08.03.2011 (m,d,y)
 
#include "parallel_grid_layout.h"
#include "pcl/pcl_util.h"

namespace ug{
AGeomObjID 	aGeomObjID("globalID", false);

template <>
unsigned long hash_key<GeomObjID>(const GeomObjID& key)
{
//	of course this hash does not completly avoid collisions.
//	One should check whether the chosen key is fine.
	return (unsigned long)(99971 * key.first + key.second * key.second);
}

std::ostream& operator<<(std::ostream& out, const GeomObjID& goId)
{
	out << "(" << goId.first << ", " << goId.second << ")";
	return out;
}

bool operator<(const GeomObjID& gid1, const GeomObjID& gid2)
{
	if(gid1.first < gid2.first)
		return true;
	if(gid1.first > gid2.first)
		return false;
	return gid1.second < gid2.second;
}

///	A helper method for GridLayoutMap::remove_empty_interfaces()
template <class TGeomObj>
static void RemoveEmptyInterfaces(
		typename GridLayoutMap::Types<TGeomObj>::Map& map)
{
	typedef typename GridLayoutMap::Types<TGeomObj>::Map TMap;
	typedef typename TMap::iterator TIterator;

	typedef typename GridLayoutMap::Types<TGeomObj>::Layout TLayout;
	typedef typename TLayout::iterator TInterfaceIter;
	typedef typename TLayout::Interface TInterface;

	for(TIterator layoutIter = map.begin(); layoutIter != map.end(); ++layoutIter)
	{
		TLayout& layout = layoutIter->second;
		RemoveEmptyInterfaces(layout);
	}
}

void GridLayoutMap::remove_empty_interfaces()
{
	RemoveEmptyInterfaces<VertexBase>(m_vertexLayoutMap);
	RemoveEmptyInterfaces<EdgeBase>(m_edgeLayoutMap);
	RemoveEmptyInterfaces<Face>(m_faceLayoutMap);
	RemoveEmptyInterfaces<Volume>(m_volumeLayoutMap);
}

}// end of namespace
