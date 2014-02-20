// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Feb 20, 2014

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
class RefinementCallbackFractal : public RefinementCallbackLinear<APosition>
{
	public:
		using RefinementCallbackLinear<APosition>::new_vertex;

	public:
		RefinementCallbackFractal();

	///	make sure that aPos is attached to the vertices of the grid.
		RefinementCallbackFractal(Grid& grid, number scaleFac,
								  APosition& aPos = aPosition);

		virtual ~RefinementCallbackFractal();

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
