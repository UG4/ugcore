#include <sstream>
#include "marker_points.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "lib_grid/file_io/file_io.h"

using namespace std;

namespace ug
{

MarkerPoint::MarkerPoint() :
	localCoord(0, 0, 0),
	associatedObj(NULL)
{
}

MarkerPoint::MarkerPoint(const char* strName, const vector3& position,
						 const vector3& normal) :
	name(strName),
	pos(position),
	norm(normal),
	localCoord(0, 0, 0),
	associatedObj(NULL)
{
}

////////////////////////////////////////////////////////////////////////
///	Loads marker points from a file
bool LoadMarkerPointsFromFile(MarkerPointManager& manager,
							  const char* filename)
{
//	load a grid and copy the points
	Grid grid;
	if(LoadGridFromFile(grid, filename)){
	//	copy the points to the marker-file
		Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);
		Grid::VertexAttachmentAccessor<ANormal> aaNorm;
		if(grid.has_vertex_attachment(aNormal))
			aaNorm.access(grid, aNormal);

		int i = 0;
		for(VertexIterator iter = grid.begin<Vertex>();
			iter != grid.end<Vertex>(); ++iter, ++i)
		{
			stringstream ss;
			ss << "marker_" << i;
			if(aaNorm.valid())
				manager.add_marker(MarkerPoint(ss.str().c_str(), aaPos[*iter],
									aaNorm[*iter]));
			else
				manager.add_marker(MarkerPoint(ss.str().c_str(), aaPos[*iter],
									vector3(0, 0, 0)));
		}

		return true;
	}

	return false;
}

////////////////////////////////////////////////////////////////////////
void SnapMarkerPointToGridVertex(MarkerPoint& markerInOut, Grid& grid,
								 number normalOffset,
								 Grid::VertexAttachmentAccessor<APosition>& aaPos,
								 Grid::VertexAttachmentAccessor<ANormal>* paaNorm)
{
	MarkerPoint& marker = markerInOut;
//	find the closest vertex
//	the specified point
	vector3& p = marker.pos;

//	find the closest point in the grid
	Vertex* vrt;
	{
		VertexIterator iter = grid.begin<Vertex>();
		vrt = *iter;
		number distSq = VecDistanceSq(aaPos[vrt], p);
		iter++;
		for(; iter != grid.end<Vertex>(); ++iter){
			number d = VecDistanceSq(aaPos[*iter], p);
			if(d < distSq){
				distSq = d;
				vrt = *iter;
			}
		}
	}

//	vrt now holds the closest vertex
//	get the normal offset
	vector3 nOff(0, 0, 0);
	if(normalOffset != 0){
		if(paaNorm)
			VecScale(nOff, (*paaNorm)[vrt], normalOffset);
		else{
			CalculateVertexNormal(nOff, grid, vrt, aaPos);
			VecScale(nOff, nOff, normalOffset);
		}
	}

//	the new position
	VecAdd(marker.pos, aaPos[vrt], nOff);

//	the new local coordinate
	marker.localCoord = vector3(0, 0, 0);

//	the new associated geometric object
	marker.associatedObj = vrt;
}

}//	end of namespace
