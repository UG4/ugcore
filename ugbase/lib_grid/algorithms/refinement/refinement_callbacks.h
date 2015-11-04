#ifndef __H__LIB_GRID__REFINEMENT_CALLBACKS__
#define __H__LIB_GRID__REFINEMENT_CALLBACKS__

#include <vector>
#include "lib_grid/lg_base.h"

namespace ug
{

///	\addtogroup lib_grid_algorithms_refinement
///	@{

////////////////////////////////////////////////////////////////////////
///	can be used to calculate new positions of vertices created during refinement.
class IRefinementCallback
{
	public:
		virtual ~IRefinementCallback()	{}
	///	called when a new vertex was created from an old vertex.
		virtual void new_vertex(Vertex* vrt, Vertex* parent) = 0;
	///	called when a new vertex was created from an old edge.
		virtual void new_vertex(Vertex* vrt, Edge* parent) = 0;
	///	called when a new vertex was created from an old face.
		virtual void new_vertex(Vertex* vrt, Face* parent) = 0;
	///	called when a new vertex was created from an old volume.
		virtual void new_vertex(Vertex* vrt, Volume* parent) = 0;

	///	callback for vertices in flat grids.
	/**	called for old vertices in flat grids which are used during refinement,
	 *	since no new vertices will be created here.
	 *	Calls new_vertex(vrt, vrt) by default.
	 *
	 *	Please note that this method won't be called during multigrid refinement.*/
		virtual void flat_grid_vertex_encountered(Vertex* vrt)	{new_vertex(vrt, vrt);}

	///	Gives access to the position of the current vertex
	/**	The method returns the number of coordinates associated with the vertex
	 * and writes those coordinates to coordsOut. No more than maxCoords will
	 * be written to coordsOut. coordsOut thus has to be at least as big as maxCoords.*/
		virtual int current_pos(number* coordsOut, Vertex* vrt, int maxCoords) = 0;

	protected:
	///	A helper implementation of current_pos available for derived classes.
		template <class TAttachmentAccessor>
		int current_pos_helper(number* coordsOut, Vertex* vrt, int maxCoords,
							   TAttachmentAccessor& aaPos);
};



////////////////////////////////////////////////////////////////////////
///	calculates new positions of vertices through linear interpolation.
/**	Make sure to initialize the callback correctly. Use the same grid
 *	on which the refinement-operations will be performed. Make sure
 *	that aPos (given in the constructor) is attached to the vertices
 *	of the grid.
 *
 *	An uninitialized refinement-callback may not be used during refinement.
 */
template <class TAPosition>
class RefinementCallbackLinear : public IRefinementCallback
{
	public:
		RefinementCallbackLinear();
		
	///	make sure that aPos is attached to the vertices of the grid.
		RefinementCallbackLinear(Grid& grid, TAPosition& aPos);
	
		virtual ~RefinementCallbackLinear();
		
		virtual void new_vertex(Vertex* vrt, Vertex* parent);
		virtual void new_vertex(Vertex* vrt, Edge* parent);
		virtual void new_vertex(Vertex* vrt, Face* parent);
		virtual void new_vertex(Vertex* vrt, Volume* parent);
		
		virtual int current_pos(number* coordsOut, Vertex* vrt, int maxCoords);

	protected:
		typedef typename TAPosition::ValueType		pos_type;
		
		Grid* 										m_pGrid;
		Grid::VertexAttachmentAccessor<TAPosition>	m_aaPos;
};

/// @}

}// end of namespace

////////////////////////////////
//	include implementation
#include "refinement_callbacks_impl.hpp"

#endif
