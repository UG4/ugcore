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

#include "subdivision_projector.h"
#include "../../subdivision/subdivision_rules_piecewise_loop.h"

using namespace std;

namespace ug{

void SubdivisionProjector::
refinement_begins(const ISubGrid& sg)
{
	m_newPositions.clear();

//	calculate new positions of all selected vertices and store them in m_newPositions
	SubdivRules_PLoop& subdiv = SubdivRules_PLoop::inst();
	Grid& g = geom().grid();
	Grid::edge_traits::secure_container edges;

	bool volumesExist = g.num<Volume>() > 0;

	lg_for_each_const (Vertex, vrt, sg.goc())
	{
		if(volumesExist) {
			if(!IsBoundaryVertex3D(g, vrt)){
				m_newPositions.push_back(make_pair(vrt, pos(vrt)));
				continue;
			}
		}

		g.associated_elements (edges, vrt);

		if(is_crease_vertex(vrt)){
		//	get the neighboured crease edges
			Edge* nbrs[2];
			size_t numNbrs = 0;
			for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge){
				Edge* e = edges[i_edge];
				if(is_crease_edge(e)){
					nbrs[numNbrs] = e;
					++numNbrs;
					if(numNbrs == 2)
						break;
				}
			}

			if(numNbrs == 2){
				vector3 p0 = pos(GetConnectedVertex(nbrs[0], vrt));
				vector3 p1 = pos(GetConnectedVertex(nbrs[1], vrt));
				vector3 w = subdiv.ref_even_crease_weights();

				vector3 p;
				VecScaleAdd(p, w.x(), pos(vrt),
							w.y(), p0, w.z(), p1);
				m_newPositions.push_back(make_pair(vrt, p));
			}
			else{
				m_newPositions.push_back(make_pair(vrt, pos(vrt)));
			}
		}
		else{
		//	perform loop subdivision on even vertices
		//	first get neighboured vertices
			size_t valence = 0;
			vector3 p;
			VecSet(p, 0);

			for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge){
				Edge* e = edges[i_edge];
				if((!volumesExist) || IsBoundaryEdge3D(g, e)){
					VecAdd(p, p, pos(GetConnectedVertex(e, vrt)));
					++valence;
				}
			}

			number centerWgt = subdiv.ref_even_inner_center_weight(valence);
			number nbrWgt = subdiv.ref_even_inner_nbr_weight(valence);

			vector3 np;
			VecScaleAdd(np, centerWgt, pos(vrt), nbrWgt, p);
			m_newPositions.push_back(make_pair(vrt, np));
		}
	} lg_end_for;
}

void SubdivisionProjector::
refinement_ends(const ISubGrid& sg)
{
//	we have to adjust positions of old vertices (no multigrid) or vertex children
//	(multigrid) here.
	Grid& g = geom().grid();
	MultiGrid* pmg = dynamic_cast<MultiGrid*>(&g);

	if(pmg){
		MultiGrid& mg = *pmg;
		for(new_pos_vec_t::iterator i = m_newPositions.begin();
			i != m_newPositions.end(); ++i)
		{
			Vertex* child = mg.get_child_vertex(i->first);
			if(child)
				set_pos(child, i->second);
			else
				set_pos(i->first, i->second);
		}
	}
	else{
		for(new_pos_vec_t::iterator i = m_newPositions.begin();
			i != m_newPositions.end(); ++i)
		{
			set_pos(i->first, i->second);
		}
	}
}

number SubdivisionProjector::
new_vertex(Vertex* vrt, Edge* parent)
{
	SubdivRules_PLoop& subdiv = SubdivRules_PLoop::inst();

//	check whether the parent edge lies inside a volume geometry. If it does,
//	perform linear refinement.
	Grid& g = geom().grid();

	if(g.num<Volume>() > 0){
		if(!IsBoundaryEdge3D(g, parent)){
			set_pos(vrt, geom().element_center(parent));
			return 1;
		}
	}

	if(is_crease_edge(parent)){
		vector2 wghts = subdiv.ref_odd_crease_weights();
		vector3 p;
		VecScaleAdd(p, wghts.x(), pos(parent->vertex(0)),
					wghts.y(), pos(parent->vertex(1)));
		set_pos(vrt, p);
	}
	else{
	//	apply loop-subdivision on inner elements
	//	get the neighboured triangles
		Face* f[2];
		int numAssociatedBndFaces = 0;
		if(g.num<Volume>() > 0){
			vector<Face*> faces;
			CollectAssociated(faces, g, parent);
			for(size_t i = 0; i < faces.size(); ++i){
				if(IsBoundaryFace3D(g, faces[i])){
					if(numAssociatedBndFaces < 2){
						f[numAssociatedBndFaces] = faces[i];
					}
					++numAssociatedBndFaces;
				}
			}
		}
		else{
			numAssociatedBndFaces = GetAssociatedFaces(f, g, parent, 2);
		}
		if(numAssociatedBndFaces == 2){
			if(f[0]->num_vertices() == 3 && f[1]->num_vertices() == 3){
			//	the 4 vertices that are important for the calculation
				Vertex* v[4];
				v[0] = parent->vertex(0); v[1] = parent->vertex(1);
				v[2] = GetConnectedVertex(parent, f[0]);
				v[3] = GetConnectedVertex(parent, f[1]);

				vector4 wghts;

				bool isCrease0 = is_crease_vertex(v[0]);
				bool isCrease1 = is_crease_vertex(v[1]);
			//	if exactly one of the two is a crease vertex, special
			//	weighting has to be performed
			//todo: this does not yet work for inner creases.
				if((isCrease0 && !isCrease1) || (!isCrease0 && isCrease1))
				{
				//	the crease vertex has to be in v[0]
					if(isCrease1)
						swap(v[0], v[1]);

				//	get the number of edges that are connected to v[0]
				//	todo: only count edges that are on the correct side of the crease.
					Grid::edge_traits::secure_container edges;
					g.associated_elements(edges, v[0]);

					wghts = subdiv.ref_odd_inner_weights(edges.size());
				}
				else{
					wghts = subdiv.ref_odd_inner_weights();
				}

				vector3 p;
				VecScaleAdd(p, wghts.x(), pos(v[0]), wghts.y(), pos(v[1]),
							   wghts.z(), pos(v[2]), wghts.w(), pos(v[3]));
				set_pos(vrt, p);
			}
			else
				RefinementProjector::new_vertex(vrt, parent);
		}
		else
			RefinementProjector::new_vertex(vrt, parent);
	}

	return 1;
}

bool SubdivisionProjector::
is_crease_vertex(Vertex* vrt)
{
	Grid& grid = geom().grid();
	if(grid.num<Volume>() > 0)
		return false;

	if(!IsRegularSurfaceVertex(grid, vrt))
		return true;

	Grid::edge_traits::secure_container edges;
	grid.associated_elements(edges, vrt);

	for(size_t i = 0; i < edges.size(); ++i){
		if(is_crease_edge(edges[i]))
			return true;
	}

	return false;
}

bool SubdivisionProjector::
is_crease_edge(Edge* edge)
{
	Grid& grid = geom().grid();
	if(grid.num<Volume>() > 0)
		return false;
	return (NumAssociatedFaces(grid, edge) != 2 ||
			m_cbIsCrease(edge));
}

}//	end of namespace
