// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 08.03.2011 (m,d,y)

#ifndef __H__UG__LIB_GRID__parallelization_util_impl__
#define __H__UG__LIB_GRID__parallelization_util_impl__

#include "copy_policy.h"
#include "pcl/pcl_interface_communicator.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
template <class TGeomObj>
void CreateAndDistributeGlobalIDs(Grid& g, GridLayoutMap& glm,
								  AGeomObjID& aID)
{
//	make sure that aID is attached to the objects of the grid
	if(!g.has_attachment<TGeomObj>(aID)){
		g.attach_to<TGeomObj>(aID);
	}

	Grid::AttachmentAccessor<TGeomObj, AGeomObjID> aaID(g, aID);

//	set up local ids
	int rank = pcl::GetProcRank();

	typedef typename geometry_traits<TGeomObj>::iterator GeomObjIter;

	size_t count = 0;
	for(GeomObjIter iter = g.begin<TGeomObj>();
		iter != g.end<TGeomObj>(); ++iter, ++count)
	{
		aaID[*iter] = MakeGeomObjID(rank, count);
	}

//	distribute the ids master->slave
	typedef typename GridLayoutMap::Types<TGeomObj>::Layout Layout;
	pcl::InterfaceCommunicator<Layout> com;
	CopyPolicy<Layout, AGeomObjID> compolCopy(g, aID);

	com.exchange_data(glm, INT_H_MASTER, INT_H_SLAVE, compolCopy);
	com.exchange_data(glm, INT_V_MASTER, INT_V_SLAVE, compolCopy);
	com.communicate();
}

}//	end of namespace

#endif
