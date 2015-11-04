#ifndef __H__UG__LIB_GRID__FRACTALS__
#define __H__UG__LIB_GRID__FRACTALS__

#include "lib_grid/algorithms/refinement/hanging_node_refiner_grid.h"
#include "lib_grid/algorithms/refinement/refinement_projectors/fractal_projector.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
///	Repeatedly refines a grid and moves new vertices along the normal by a given factor.
template <class TAPosition>
bool CreateFractal_NormalScale(Grid& grid, HangingNodeRefiner_Grid& href,
							   number scaleFac, size_t numIterations,
							   TAPosition& aPosVRT)
{
	if(!grid.has_vertex_attachment(aPosVRT)){
		UG_LOG("WARNING in CreateFractal_NormalScale: given position attachment is ");
		UG_LOG("not attached to the grids vertices. Aborting...\n");
		return false;
	}
		
//	store the old refinement callback of href
	IRefinementCallback* oldCallback = href.get_refinement_callback();
	
//	create the new one.
	FractalProjector refCallback(grid, scaleFac);
	href.set_refinement_callback(&refCallback);
	
//	iterate for the specified number of times
	for(size_t i = 0; i < numIterations; ++i){

		if(grid.num_volumes() > 0){
		//	iterate over all faces and mark them for refinement, if they are boundary faces.
			for(FaceIterator iter = grid.faces_begin();
				iter != grid.faces_end(); ++iter)
			{
				if(IsVolumeBoundaryFace(grid, *iter))
					href.mark(*iter);
			}
		}
		else if(grid.num_faces() > 0){
		//	markall faces
			href.mark(grid.faces_begin(), grid.faces_end());
		}
		else{
		//	mark all edges
			href.mark(grid.edges_begin(), grid.edges_end());
		}

	//	refine them
		href.refine();

	//	change the scalefac
		refCallback.set_scale_fac(-0.5 * refCallback.get_scale_fac());
		//refCallback.set_scale_fac(-0.5 * refCallback.get_scale_fac());
		//refCallback.set_scale_fac(refCallback.get_scale_fac() * refCallback.get_scale_fac());

	}

//	done. restore href
	href.set_refinement_callback(oldCallback);
	return true;
}

inline bool CreateFractal_NormalScale(Grid& grid,
									   HangingNodeRefiner_Grid& href,
									   number scaleFac, size_t numIterations)
{
	return CreateFractal_NormalScale(grid, href, scaleFac,
									 numIterations, aPosition);
}
							   
}

#endif
