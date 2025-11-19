/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Stepniewski
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

#include "common/math/misc/math_util.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "subdivision_loop.h"

using namespace std;

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	ProjectToLimitLoop
/// projects surface vertices to their limit subdivision surface position
bool ProjectToLimitLoop(Grid& grid, APosition& aProjPos)
{
//	grid management
	Grid::VertexAttachmentAccessor aaPos(grid,aPosition);
	Grid::VertexAttachmentAccessor aaProjPos(grid, aProjPos);

//	needed variables

	double x = 0;
	double y = 0;
	double z = 0;
	int valence = 0;
	constexpr int numPrecalculated = 10;
	double beta[numPrecalculated];
	double b = 0;
	double chi = 0;

//	calculate weights for subdivision mask
	for(int i = 1; i < numPrecalculated; ++i)
	{
		double tmp = 0.375 + 0.25 * cos( (2.0 * M_PI) / (float)i );
		beta[i] = ( 0.625 - tmp * tmp ) / (float)i ;
	}

	beta[0] = 0;
	beta[6] = 0.0625;

//	iterate through all vertices, evaluate their limit positions and save them in their projection attachment
	for(VertexIterator vIter = grid.vertices_begin(); vIter != grid.vertices_end(); ++vIter)
	{
		Vertex* v = *vIter;
		valence = 0;
		x = 0;
		y = 0;
		z = 0;

		for(auto eIter = grid.associated_edges_begin(v); eIter != grid.associated_edges_end(v); ++eIter)
		{
			Edge* e = *eIter;
			valence++;

			if(valence >= numPrecalculated)
			{
				double tmp = 0.375 + 0.25 * cos( (2.0*M_PI) / (float)valence );
				b = (0.625 - tmp*tmp) / (float)valence;
			}

			else
				b = beta[valence];

			chi = 1.0 / (0.375 / b + valence);

			if(aaPos[v].x() == aaPos[e->vertex(0)].x() && aaPos[v].y() == aaPos[e->vertex(0)].y() && aaPos[v].z() == aaPos[e->vertex(0)].z())
			{
				x += aaPos[e->vertex(1)].x();
				y += aaPos[e->vertex(1)].y();
				z += aaPos[e->vertex(1)].z();
			}

			else
			{
				x += aaPos[e->vertex(0)].x();
				y += aaPos[e->vertex(0)].y();
				z += aaPos[e->vertex(0)].z();
			}
		}

		x*=chi;
		y*=chi;
		z*=chi;

		aaProjPos[v].x() = aaPos[v].x() * (1.0 - (float)valence * chi);
		aaProjPos[v].y() = aaPos[v].y() * (1.0 - (float)valence * chi);
		aaProjPos[v].z() = aaPos[v].z() * (1.0 - (float)valence * chi);

		aaProjPos[v].x() += x;
		aaProjPos[v].y() += y;
		aaProjPos[v].z() += z;
	}

	return true;
}

}//	end of namespace

