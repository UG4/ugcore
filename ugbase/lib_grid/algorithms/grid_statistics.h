//	created by Martin Stepniewski
//	mastep@gmx.de
//	y09 m11 d11

#ifndef __H__UG__GRID_STATISTICS__
#define __H__UG__GRID_STATISTICS__

#include <vector>
#include "lib_grid/lg_base.h"

namespace ug
{
	
//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************
//								declarations
//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************

////////////////////////////////////////////////////////////////////////
//	AssignTetrahedronAttributesByAspectRatio - mstepnie
/// assigns tetrahedral elements of a grid to subsets respecting their aspect ratio
bool AssignTetrahedronAttributesByAspectRatio(Grid& grid,
											  SubsetHandler& shVolume,
											  AInt& aTetrahedronAspectRatioClass,
											  std::vector<double>& offsets);


}//	end of namespace

#endif
