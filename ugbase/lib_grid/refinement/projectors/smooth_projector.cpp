/*
 * Copyright (c) 2016:  G-CSC, Goethe University Frankfurt
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

#include "smooth_projector.h"

namespace ug{

void SmoothProjector::
refinement_ends ()	
{
	Grid::edge_traits::secure_container		edges;
	Grid::face_traits::secure_container		faces;
	Grid::volume_traits::secure_container	volumes;

	Grid& grid = geom().grid();

	std::vector<Vertex*>	localNbrs;

	for(int iteration = 0; iteration < m_iterations; ++iteration){
		for(size_t inewVrt = 0; inewVrt < m_newVrts.size(); ++inewVrt){
			Vertex* vrt = m_newVrts[inewVrt];
			localNbrs.clear();
			grid.associated_elements(edges, vrt);
			grid.associated_elements(faces, vrt);
			grid.associated_elements(volumes, vrt);

			grid.begin_marking();
			grid.mark(vrt);
			for(size_t iass = 0; iass < edges.size(); ++iass){
				Edge* e = edges[iass];
				for(size_t ivrt = 0; ivrt < e->num_vertices(); ++ivrt){
					Vertex* nbr = e->vertex(ivrt);
					if(!grid.is_marked(nbr)){
						grid.mark(nbr);
						localNbrs.push_back(nbr);
					}
				}
			}
			for(size_t iass = 0; iass < faces.size(); ++iass){
				Face* e = faces[iass];
				for(size_t ivrt = 0; ivrt < e->num_vertices(); ++ivrt){
					Vertex* nbr = e->vertex(ivrt);
					if(!grid.is_marked(nbr)){
						grid.mark(nbr);
						localNbrs.push_back(nbr);
					}
				}
			}
			for(size_t iass = 0; iass < volumes.size(); ++iass){
				Volume* e = volumes[iass];
				for(size_t ivrt = 0; ivrt < e->num_vertices(); ++ivrt){
					Vertex* nbr = e->vertex(ivrt);
					if(!grid.is_marked(nbr)){
						grid.mark(nbr);
						localNbrs.push_back(nbr);
					}
				}
			}
			grid.end_marking();

			const size_t numNbrs = localNbrs.size();
			if(numNbrs > 0){
				number wgt = m_changeRate / (number)numNbrs;
				vector3 weigtedCenter(0, 0, 0);
				for(size_t _vfeI = 0; _vfeI < localNbrs.size(); ++_vfeI){ Vertex* nbr = localNbrs[_vfeI];{
					vector3 nbrPos = pos(nbr);
					nbrPos *= wgt;
					weigtedCenter += nbrPos;
				}};

				vector3 newPos = pos(vrt);
				newPos *= (1. - m_changeRate);
				newPos += weigtedCenter;
				set_pos(vrt, newPos);
			}
		}
	}
}

}//	end of namespace
