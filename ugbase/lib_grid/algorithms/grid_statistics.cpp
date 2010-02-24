//	created by Martin Stepniewski
//	mastep@gmx.de
//	y09 m11 d11

#include "grid_statistics.h"
#include "geom_obj_util/volume_util.h"

using namespace std;

namespace ug
{
//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************
//								definitions
//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	AssignTetrahedronAttributesByAspectRatio
/// assigns attributes to tetrahedral elements of a grid respecting their aspect ratio
bool AssignTetrahedronAttributesByAspectRatio(	Grid& grid,
												SubsetHandler& shVolume,
												AInt& aTetrahedronAspectRatioClass,
												vector<double>& offsets)
{
	VolumeSelector VolumeSel(grid);
	Grid::VolumeAttachmentAccessor<AInt> aaTARC(grid, aTetrahedronAspectRatioClass);

//	order the volume subsets correctly, regarding user-specified aElementAttribute
	for(int i = 1; i < (int)offsets.size(); ++i)
	{
		int moveTo = i;
		for(int j = i-1; j >= 0; --j)
		{
			if(offsets[i] > offsets[j])
				moveTo = j;
			else
				break;
		}

		if(moveTo != i)
		{
			swap(offsets[i], offsets[moveTo]);
		}
	}

//	add zero as last offset
	offsets.push_back(0.0);

//	iterate through all tetrahedrons and assign attributes by their aspect ratio
	for(VolumeIterator vIter = grid.volumes_begin(); vIter != grid.volumes_end(); ++vIter)
	{
		Volume* v = *vIter;
		number AspectRatio = CalculateTetrahedronAspectRatio(grid, *v);

		for(size_t i = 0; i<offsets.size(); ++i)
		{
			if(AspectRatio >= offsets[i])
			{
				shVolume.assign_subset(v, i);
				aaTARC[v] = i;
				break;
			}
		}
	}

	return true;
}

}//	end of namespace
