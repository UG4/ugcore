#include <vector>
#include "extrude.h"
#include "cylinder_extrusion.h"
#include "lib_grid/algorithms/remeshing/grid_adaption.h"
#include "lib_grid/algorithms/selection_util.h"

using namespace std;

namespace ug
{

//	this method actually performs the extrusion. It is only callable from
//	inside this file.
static bool ExtrudeCylinder(Grid& grid, SubsetHandler* sh, Vertex* vrt,
					const vector3& direction, number height, number radius,
					number rimSnapThreshold,
					Grid::VertexAttachmentAccessor<APosition>& aaPos,
					int bottomSubInd, int cylSubInd,
					Selector* pSel)
{
//	we'll use this selector if no selector was specified
	Selector locSel;
	if(!pSel){
		locSel.assign_grid(grid);
		pSel = &locSel;
	}
	Selector& sel = *pSel;
	
	if(!AdaptSurfaceGridToCylinder(sel, grid, vrt, direction, radius, rimSnapThreshold)){
		LOG("  WARNING: AdaptSurfaceGridToCylinder failed during ExtrudeCylinder.\n");
		return false;
	}

//	select boundary edges
	sel.clear<Edge>();
	SelectAreaBoundary(sel, sel.begin<Face>(), sel.end<Face>());

//	gather faces and edges for extrusion
	vector<Edge*> vEdges;
	vector<Face*> vFaces;

	for(EdgeIterator iter = sel.begin<Edge>();
		iter != sel.end<Edge>(); ++iter)
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
	
	Extrude(grid, NULL, &vEdges, &vFaces, scaledDir, EO_CREATE_FACES, aPosition);

	if(sh){
	//	restore subset handler
		sh->enable_subset_inheritance(bSInh);
		sh->set_default_subset_index(defSubInd);
	}

	return true;
}

bool ExtrudeCylinder(Grid& grid, SubsetHandler& sh, Vertex* vrt,
					const vector3& direction, number height, number radius,
					number rimSnapThreshold,
					Grid::VertexAttachmentAccessor<APosition>& aaPos,
					int bottomSubInd, int cylSubInd,
					Selector* pSel)
{
	return ExtrudeCylinder(grid, &sh, vrt, direction, height, radius,
						 rimSnapThreshold, aaPos, bottomSubInd, cylSubInd, pSel);

}

bool ExtrudeCylinder(Grid& grid, Vertex* vrt,
					const vector3& direction, number height, number radius,
					number rimSnapThreshold,
					Grid::VertexAttachmentAccessor<APosition>& aaPos,
					Selector* pSel)
{
	return ExtrudeCylinder(grid, NULL, vrt, direction, height, radius,
							rimSnapThreshold, aaPos, -1, -1, pSel);
}
					
}//	end of namespace
