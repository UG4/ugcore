/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#include "grid_objects_1d.h"

namespace ug{

////////////////////////////////////////////////////////////////////////
//	RegularEdge
bool RegularEdge::refine(std::vector<Edge*>& vNewEdgesOut, Vertex* newVertex,
				Vertex** pSubstituteVrts)
{
	vNewEdgesOut.clear();
	if(pSubstituteVrts)
	{
		vNewEdgesOut.push_back(new RegularEdge(pSubstituteVrts[0], newVertex));
		vNewEdgesOut.push_back(new RegularEdge(newVertex, pSubstituteVrts[1]));
	}
	else
	{
		vNewEdgesOut.push_back(new RegularEdge(vertex(0), newVertex));
		vNewEdgesOut.push_back(new RegularEdge(newVertex, vertex(1)));
	}
	return true;
}

bool RegularEdge::refine(std::vector<RegularEdge*>& vNewEdgesOut, Vertex* newVertex,
				  Vertex** pSubstituteVrts)
{
	return refine(reinterpret_cast<std::vector<Edge*>&>(vNewEdgesOut),
					newVertex, pSubstituteVrts);
}

////////////////////////////////////////////////////////////////////////
//	ConstrainedEdge
bool ConstrainedEdge::refine(std::vector<Edge*>& vNewEdgesOut, Vertex* newVertex,
				  Vertex** pSubstituteVrts)
{
	vNewEdgesOut.clear();
	if(pSubstituteVrts)
	{
		vNewEdgesOut.push_back(new ConstrainedEdge(pSubstituteVrts[0], newVertex));
		vNewEdgesOut.push_back(new ConstrainedEdge(newVertex, pSubstituteVrts[1]));
	}
	else
	{
		vNewEdgesOut.push_back(new ConstrainedEdge(vertex(0), newVertex));
		vNewEdgesOut.push_back(new ConstrainedEdge(newVertex, vertex(1)));
	}
	return true;
}

bool ConstrainedEdge::refine(std::vector<ConstrainedEdge*>& vNewEdgesOut, Vertex* newVertex,
				  Vertex** pSubstituteVrts)
{
	return refine(reinterpret_cast<std::vector<Edge*>&>(vNewEdgesOut),
				newVertex, pSubstituteVrts);
}

////////////////////////////////////////////////////////////////////////
//	ConstrainingEdge
bool ConstrainingEdge::refine(std::vector<Edge*>& vNewEdgesOut, Vertex* newVertex,
				  Vertex** pSubstituteVrts)
{
	vNewEdgesOut.clear();
	if(pSubstituteVrts)
	{
		vNewEdgesOut.push_back(new ConstrainingEdge(pSubstituteVrts[0], newVertex));
		vNewEdgesOut.push_back(new ConstrainingEdge(newVertex, pSubstituteVrts[1]));
	}
	else
	{
		vNewEdgesOut.push_back(new ConstrainingEdge(vertex(0), newVertex));
		vNewEdgesOut.push_back(new ConstrainingEdge(newVertex, vertex(1)));
	}
	return true;
}

bool ConstrainingEdge::refine(std::vector<ConstrainingEdge*>& vNewEdgesOut, Vertex* newVertex,
				  Vertex** pSubstituteVrts)
{
	return refine(reinterpret_cast<std::vector<Edge*>&>(vNewEdgesOut),
					newVertex, pSubstituteVrts);
}

template <> size_t
ConstrainingEdge::
num_constrained<Vertex>() const
{
	return num_constrained_vertices();
}

template <> size_t
ConstrainingEdge::
num_constrained<Edge>() const
{
	return num_constrained_edges();
}

template <> Vertex*
ConstrainingEdge::
constrained<Vertex>(size_t ind) const
{
	return constrained_vertex(ind);
}

template <> Edge*
ConstrainingEdge::
constrained<Edge>(size_t ind) const
{
	return constrained_edge(ind);
}

}// end of namespace
