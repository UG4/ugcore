/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG__misc_refinement_projectors__
#define __H__UG__misc_refinement_projectors__

#include "refinement_callbacks.h"

namespace ug{

///	\addtogroup lib_grid_algorithms_refinement
///	@{

///	If a refined edge intersects the given cylinder, the new vertex will be placed at the intersection.
class RefinementCallback_IntersectCylinder : public RefinementCallbackLinear<APosition>
{
	public:
		using RefinementCallbackLinear<APosition>::new_vertex;

	public:
		RefinementCallback_IntersectCylinder();

	///	make sure that aPos is attached to the vertices of the grid.
		RefinementCallback_IntersectCylinder(Grid& grid, const vector3& center,
											 const vector3& axis, number radius,
											 APosition& aPos = aPosition);

		virtual ~RefinementCallback_IntersectCylinder();

		virtual void new_vertex(Vertex* vrt, Edge* parent);

	protected:
		vector3 m_center;
		vector3 m_axis;
		number m_radius;
};

////////////////////////////////////////////////////////////////////////
///	calculates new positions by cutting parent edges with a plane
/**	For each edge the intersection of the edge with the initially
 *	given plane is calculated and used as new point. Vertices
 *	created on other geometric objects are treated as in the linear case.
 *
 *	Make sure to initialise the callback correctly. Use the same grid
 *	on which the refinement-operations will be performed. Make sure
 *	that aPos (given in the constructor) is attached to the vertices
 *	of the grid.
 *
 *	An uninitialized refinement-callback may not be used during refinement.
 */
class RefinementCallbackEdgePlaneCut : public RefinementCallbackLinear<APosition>
{
	public:
		using RefinementCallbackLinear<APosition>::new_vertex;

	public:
		RefinementCallbackEdgePlaneCut();

	///	make sure that aPos is attached to the vertices of the grid.
		RefinementCallbackEdgePlaneCut(Grid& grid, const vector3& p,
										const vector3& n,
										APosition& aPos = aPosition);

		virtual ~RefinementCallbackEdgePlaneCut();

		virtual void new_vertex(Vertex* vrt, Edge* parent);

	protected:
		vector3 m_p;
		vector3 m_n;
};




/// @}

}// end of namespace

#endif
