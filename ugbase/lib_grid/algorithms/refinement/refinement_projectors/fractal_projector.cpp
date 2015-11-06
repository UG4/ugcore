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

#include "fractal_projector.h"

namespace ug{


////////////////////////////////////////////////////////////////////////
FractalProjector::
FractalProjector()
{
}

////////////////////////////////////////////////////////////////////////
FractalProjector::
FractalProjector(Grid& grid, number scaleFac, APosition& aPos) :
	RefinementCallbackLinear<APosition>(grid, aPos),
	m_scaleFac(scaleFac)
{
}

////////////////////////////////////////////////////////////////////////
FractalProjector::
~FractalProjector()
{
}

////////////////////////////////////////////////////////////////////////
void FractalProjector::
new_vertex(Vertex* vrt, Edge* parent)
{
//	set the vertex to the center by calling the parents method
	RefinementCallbackLinear<APosition>::new_vertex(vrt, parent);

//	calculate the normal of the edge
	vector3 n;
	CalculateNormal(n, *m_pGrid, parent, m_aaPos);
//	first scale it to the same length as the edge, then scale it with the
//	given scaleFac
	number len = VecDistance(m_aaPos[parent->vertex(0)], m_aaPos[parent->vertex(1)]);
	VecScale (n, n, m_scaleFac * len);

//	offset the vertex by the calculated normal
	VecAdd(m_aaPos[vrt], m_aaPos[vrt], n);
}

////////////////////////////////////////////////////////////////////////
void FractalProjector::
new_vertex(Vertex* vrt, Face* parent)
{
//	set the vertex to the center by calling the parents method
	RefinementCallbackLinear<APosition>::new_vertex(vrt, parent);

//	calculate the normal of the face
	vector3 n;
	CalculateNormal(n, parent, m_aaPos);

//	we have to alter the length.
//	first scale it to the average length of the faces edges,
//	then scale it with the given scaleFac.
	number avLen = 0;

	size_t numVrts = parent->num_vertices();
	for(size_t i = 0; i < numVrts; ++i){
		avLen += VecDistance(m_aaPos[parent->vertex(i)],
							 m_aaPos[parent->vertex((i+1)% numVrts)]);
	}
	avLen /= (number)numVrts;
	VecScale (n, n, m_scaleFac * avLen);

//	offset the vertex by the calculated normal
	VecAdd(m_aaPos[vrt], m_aaPos[vrt], n);
}

}// end of namespace
