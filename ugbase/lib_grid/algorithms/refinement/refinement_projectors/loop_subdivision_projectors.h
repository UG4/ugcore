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

#ifndef __H__UG__loop_subdivision_projectors__
#define __H__UG__loop_subdivision_projectors__

#include "../refinement_callbacks.h"

namespace ug{

///	\addtogroup lib_grid_algorithms_refinement
///	@{

///	Applies smooth subdivision rules on boundary elements.
/**	Inner vertices are positioned using the standard linear scheme,
 *  boundary vertices are positioned using b-spline weights.
 *
 *  Please note that the refiner works on different position attachments:
 *  A source and a target position attachment. If you use the refiner in the
 *  context of multigrid refinement, you can use the same attachment for both.
 *  However, if the callback is used to refine a flat grid, then it's crucial
 *  to use different attachments..
 */
template <class TAPosition>
class SubdivisionLoopBoundaryProjector : public RefinementCallbackLinear<TAPosition>
{
	private:
		typedef RefinementCallbackLinear<TAPosition> BaseClass;
		typedef typename TAPosition::ValueType		pos_type;

	public:
		using BaseClass::new_vertex;

	public:
		SubdivisionLoopBoundaryProjector();

	///	make sure that aPos and aTargetPos are attached to the vertices of the grid.
	/**	For flat grids it is crucial that aPos and aTargetPos differ.
	 *  If g is a multigrid, that aPos and aTargetPos can safely be the same.*/
		SubdivisionLoopBoundaryProjector(Grid& g,
										 TAPosition& aPos,
										 TAPosition& aTargetPos);

		virtual ~SubdivisionLoopBoundaryProjector();

		virtual void new_vertex(Vertex* vrt, Vertex* parent);
		virtual void new_vertex(Vertex* vrt, Edge* parent);
		virtual void new_vertex(Vertex* vrt, Face* parent);
		virtual void new_vertex(Vertex* vrt, Volume* parent);

	protected:
		virtual bool is_crease_vertex(Vertex* vrt);
		virtual bool is_crease_edge(Edge* edge);

	protected:
		Grid::VertexAttachmentAccessor<TAPosition>	m_aaTargetPos;
};


///	Refines grids with an extended loop scheme.
/**	Inner vertices are positioned using the loop scheme, boundary vertices
 *  are positioned using b-spline weights.
 *
 *  Please note that the refiner works on different position attachments:
 *  A source and a target position attachment. If you use the refiner in the
 *  context of multigrid refinement, you can use the same attachment for both.
 *  However, if the callback is used to refine a flat grid, then it's crucial
 *  to use different attachments.
 */
template <class TAPosition>
class SubdivisionLoopProjector : public RefinementCallbackLinear<TAPosition>
{
	private:
		typedef RefinementCallbackLinear<TAPosition> BaseClass;
		typedef typename TAPosition::ValueType		pos_type;

	public:
		using BaseClass::new_vertex;

	public:
		SubdivisionLoopProjector();

	///	make sure that aPos and aTargetPos are attached to the vertices of the grid.
	/**	For flat grids it is crucial that aPos and aTargetPos differ.
	 *  If g is a multigrid, that aPos and aTargetPos can safely be the same.*/
		SubdivisionLoopProjector(Grid& g,
										  TAPosition& aPos,
										  TAPosition& aTargetPos);

		virtual ~SubdivisionLoopProjector();

		virtual void new_vertex(Vertex* vrt, Vertex* parent);
		virtual void new_vertex(Vertex* vrt, Edge* parent);
		virtual void new_vertex(Vertex* vrt, Face* parent);
		virtual void new_vertex(Vertex* vrt, Volume* parent);

		virtual void consider_as_crease_edge(Grid::edge_traits::callback cbIsCrease);

	protected:
		virtual bool is_crease_vertex(Vertex* vrt);
		virtual bool is_crease_edge(Edge* edge);

	protected:
		Grid::VertexAttachmentAccessor<TAPosition>	m_aaTargetPos;
		Grid::edge_traits::callback					m_cbIsCrease;
};
/// @}

}// end of namespace

#include "loop_subdivision_projectors_impl.h"

#endif
