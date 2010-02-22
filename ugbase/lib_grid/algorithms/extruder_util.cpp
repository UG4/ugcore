//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m02 d22

#include <vector>
#include "extruder_util.h"
#include "lib_grid/algorithms/remeshing/grid_adaption.h"
#include "lib_grid/algorithms/selection_util.h"

using namespace std;

namespace ug
{

//	this method actually performs the extrusion. It is only callable from
//	inside this file.
static bool ExtrudeCylinder(Grid& grid, SubsetHandler* sh, VertexBase* vrt,
					const vector3& direction, number height, number radius,
					Grid::VertexAttachmentAccessor<APosition>& aaPos,
					int bottomSubInd, int cylSubInd,
					Selector* pSel, float minDot)
{
//	we'll use this selector if no selector was specified
	Selector locSel;
	if(!pSel){
		locSel.assign_grid(grid);
		pSel = &locSel;
	}
	Selector& sel = *pSel;
	
	if(!AdaptSurfaceGridToCylinder(sel, grid, vrt, direction, radius, minDot)){
		LOG("  WARNING: AdaptSurfaceGridToCylinder failed during ExtrudeCylinder.\n");
		return false;
	}

//	select boundary edges
	sel.clear<EdgeBase>();
	SelectAreaBoundaryEdges(sel, sel.begin<Face>(), sel.end<Face>());

//	gather faces and edges for extrusion
	vector<EdgeBase*> vEdges;
	vector<Face*> vFaces;

	for(EdgeBaseIterator iter = sel.begin<EdgeBase>();
		iter != sel.end<EdgeBase>(); ++iter)
		vEdges.push_back(*iter);

	for(FaceIterator iter = sel.begin<Face>();
		iter != sel.end<Face>(); ++iter)
		vFaces.push_back(*iter);

//	everything that is extruded shall go into a separate subset
	bool bSInh = false;
	int defSubInd = -1;
	if(sh){
	//	assign faces that will be extruded into a separate subset
		sh->assign_subset(vFaces.begin(), vFaces.end(), bottomSubInd);
	//	the rest gets into the next subset
		bSInh = sh->subset_inheritance_enabled();
		defSubInd = sh->get_default_subset_index();
		sh->enable_subset_inheritance(false);
		sh->set_default_subset_index(cylSubInd);		
	}

	vector3 scaledDir;
	VecScale(scaledDir, direction, height);
	
	Extrude(grid, NULL, &vEdges, &vFaces, scaledDir, EO_CREATE_FACES, aPosition, (height < 0));

	if(sh){
	//	restore subset handler
		sh->enable_subset_inheritance(bSInh);
		sh->set_default_subset_index(defSubInd);
	}

	return true;
}

bool ExtrudeCylinder(Grid& grid, SubsetHandler& sh, VertexBase* vrt,
					const vector3& direction, number height, number radius,
					Grid::VertexAttachmentAccessor<APosition>& aaPos,
					int bottomSubInd, int cylSubInd,
					Selector* pSel, float minDot)
{
	return ExtrudeCylinder(grid, &sh, vrt, direction, height, radius, aaPos,
							bottomSubInd, cylSubInd, pSel, minDot);

}

bool ExtrudeCylinder(Grid& grid, VertexBase* vrt,
					const vector3& direction, number height, number radius,
					Grid::VertexAttachmentAccessor<APosition>& aaPos,
					Selector* pSel, float minDot)
{
	return ExtrudeCylinder(grid, NULL, vrt, direction, height, radius, aaPos,
							-1, -1, pSel, minDot);
}
					
}//	end of namespace
