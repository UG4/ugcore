#ifndef __H__UG__fractal_projector__
#define __H__UG__fractal_projector__

#include "../refinement_callbacks.h"

namespace ug{

///	\addtogroup lib_grid_algorithms_refinement
///	@{

////////////////////////////////////////////////////////////////////////
///	calculates new positions by offsetting new vertices to along its parents normal.
/**
 *	Make sure to initialise the callback correctly. Use the same grid
 *	on which the refinement-operations will be performed. Make sure
 *	that aPos (given in the constructor) is attached to the vertices
 *	of the grid.
 *
 *	An uninitialized refinement-callback may not be used during refinement.
 */
class FractalProjector : public RefinementCallbackLinear<APosition>
{
	public:
		using RefinementCallbackLinear<APosition>::new_vertex;

	public:
		FractalProjector();

	///	make sure that aPos is attached to the vertices of the grid.
		FractalProjector(Grid& grid, number scaleFac,
								  APosition& aPos = aPosition);

		virtual ~FractalProjector();

		virtual void new_vertex(Vertex* vrt, Edge* parent);
		virtual void new_vertex(Vertex* vrt, Face* parent);

		inline void set_scale_fac(number scaleFac)	{m_scaleFac = scaleFac;}
		inline number get_scale_fac()				{return m_scaleFac;}

	protected:
		number m_scaleFac;
};

/// @}

}// end of namespace

#endif
