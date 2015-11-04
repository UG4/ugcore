#ifndef __H__UG__LIB_GRID__parallelization_util_impl__
#define __H__UG__LIB_GRID__parallelization_util_impl__

#include "util/compol_copy_attachment.h"
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
	int rank = pcl::ProcRank();

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
	ComPol_CopyAttachment<Layout, AGeomObjID> compolCopy(g, aID);

	com.exchange_data(glm, INT_H_MASTER, INT_H_SLAVE, compolCopy);
	com.communicate();

//	we copy data from vertical slaves to vertical masters and not vice versa for
//	the following reason:
//	Multiple copies of low dimensional elements may reside in vertical master
//	interfaces and they do not necessarily have the same global id yet.
//	Elements in horizontal interfaces, however, already have the same ids yet.
//	Since all vertical slaves are either unique or in a horizontal interface,
//	we can simply copy their global ids to vertical master elements.
	com.exchange_data(glm, INT_V_SLAVE, INT_V_MASTER, compolCopy);
	com.communicate();
}

}//	end of namespace

#endif
