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
#include "../../algorithms/subdivision/subdivision_rules_piecewise_loop.h"

using namespace std;

namespace ug{

void SubdivisionProjector::
refinement_begins(const ISubGrid* psg)
{
//	PLEASE NOTE:	The implementation below is able to perform subdivision on
//					manifolds which consist of sides of volume-elements, even
//					if volume-elements are on both sides of those manifolds.
	RefinementProjector::refinement_begins(psg);
	const ISubGrid& sg = *psg;

	m_newPositions.clear();

//	calculate new positions of all selected vertices and store them in m_newPositions
	SubdivRules_PLoop& subdiv = SubdivRules_PLoop::inst();
	Grid& g = geom().grid();
	Grid::edge_traits::secure_container edges;

	for(auto _feI = sg.goc().begin<Vertex>(); _feI != sg.goc().end<Vertex>(); ++_feI){
		Vertex* vrt = *_feI;
		g.associated_elements (edges, vrt);

		Edge* creaseNbrs[2];
		const size_t numCreaseNbrs = nbr_crease_edges(vrt, &edges, creaseNbrs);
		if(numCreaseNbrs == 2){
			vector3 p0 = pos(GetConnectedVertex(creaseNbrs[0], vrt));
			vector3 p1 = pos(GetConnectedVertex(creaseNbrs[1], vrt));
			vector3 w = subdiv.ref_even_crease_weights();

			vector3 p;
			VecScaleAdd(p, w.x(), pos(vrt),
						w.y(), p0, w.z(), p1);
			m_newPositions.push_back(make_pair(vrt, p));
		}
		else if(numCreaseNbrs > 0){
			m_newPositions.emplace_back(vrt, pos(vrt));
		}
		else{
		//	perform loop subdivision on even vertices
			size_t valence = 0;
			vector3 p;
			VecSet(p, 0);

			Grid::face_traits::secure_container faces;
			for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge){
				Edge* e = edges[i_edge];
				g.associated_elements(faces, e);
				if(concerned_nbr_faces(e, &faces)){
					VecAdd(p, p, pos(GetConnectedVertex(e, vrt)));
					++valence;
				}
			}

			number centerWgt = subdiv.ref_even_inner_center_weight(valence);
			number nbrWgt = subdiv.ref_even_inner_nbr_weight(valence);

			vector3 np;
			VecScaleAdd(np, centerWgt, pos(vrt), nbrWgt, p);
			m_newPositions.emplace_back(vrt, np);
		}
	}
}

void SubdivisionProjector::
refinement_ends()
{
//	we have to adjust positions of old vertices (no multigrid) or vertex children
//	(multigrid) here.
	Grid& g = geom().grid();
	auto* pmg = dynamic_cast<MultiGrid*>(&g);

	if(pmg){
		MultiGrid& mg = *pmg;
		for(auto i = m_newPositions.begin(); i != m_newPositions.end(); ++i)
		{
			Vertex* child = mg.get_child_vertex(i->first);
			if(child)
				set_pos(child, i->second);
			else
				set_pos(i->first, i->second);
		}
	}
	else{
		for(auto i = m_newPositions.begin(); i != m_newPositions.end(); ++i)
		{
			set_pos(i->first, i->second);
		}
	}
}

number SubdivisionProjector::
new_vertex(Vertex* vrt, Edge* parent)
{
//todo: Adjust to new helper methods.
	
	const SubdivRules_PLoop& subdiv = SubdivRules_PLoop::inst();

//	check whether the parent edge lies inside a volume geometry. If it does,
//	perform linear refinement.
	Grid& g = geom().grid();

	const bool parentIsCrease = m_cbIsCrease(parent);
	Face* concernedNbrs[2];
	size_t numConcernedFaces = 0;

	if(!parentIsCrease)
		numConcernedFaces = concerned_nbr_faces (parent, nullptr, concernedNbrs);

	if(parentIsCrease || numConcernedFaces != 2){
		vector2 wghts = subdiv.ref_odd_crease_weights();
		vector3 p;
		VecScaleAdd(p, wghts.x(), pos(parent->vertex(0)),
					wghts.y(), pos(parent->vertex(1)));
		set_pos(vrt, p);
		return 1;
	}
	
	if(concernedNbrs[0]->num_vertices() == 3 && concernedNbrs[1]->num_vertices() == 3){
	//	the 4 vertices that are important for the calculation
		Vertex* v[4];
		v[0] = parent->vertex(0); v[1] = parent->vertex(1);
		v[2] = GetConnectedVertex(parent, concernedNbrs[0]);
		v[3] = GetConnectedVertex(parent, concernedNbrs[1]);

		vector4 wghts;

		bool isCrease0 = nbr_crease_edges(v[0]) > 0;
		bool isCrease1 = nbr_crease_edges(v[1]) > 0;
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

	return 1;
}

size_t SubdivisionProjector::
nbr_crease_edges (Vertex* vrt,
				  Grid::edge_traits::secure_container* assEdges,
				  Edge* creaseEdgesOut[2])
{
	Grid& grid = geom().grid();

	Grid::edge_traits::secure_container	localAssEdges;
	if(!assEdges){
		grid.associated_elements(localAssEdges, vrt);
		assEdges = &localAssEdges;
	}

	size_t creaseEdges = 0;
	Grid::face_traits::secure_container faces;
	if(creaseEdgesOut){
		for(size_t i = 0; i < assEdges->size(); ++i){
			Edge* edge = (*assEdges)[i];
			if(m_cbIsCrease(edge)){
				if(creaseEdges < 2)
					creaseEdgesOut[creaseEdges] = edge;
				++creaseEdges;
			}
			else{
				grid.associated_elements(faces, edge);
				size_t num = 0;
				if(!faces.empty())
					num = concerned_nbr_faces(edge, &faces);
				if(faces.empty() || ((num > 0) && (num!= 2))){
					if(creaseEdges < 2)
						creaseEdgesOut[creaseEdges] = edge;
					++creaseEdges;
				}
			}
		}
	}
	else{
		for(size_t i = 0; i < assEdges->size(); ++i){
			Edge* edge = (*assEdges)[i];
			if(m_cbIsCrease(edge))
				++creaseEdges;
			else{
				grid.associated_elements(faces, edge);
				if(faces.empty())
					++creaseEdges;
				else{
					size_t num = concerned_nbr_faces(edge, &faces);
					creaseEdges += static_cast<int>((num > 0) && (num != 2));
				}
			}
		}
	}

	return creaseEdges;
}

size_t SubdivisionProjector::
concerned_nbr_faces (Edge* edge,
					 Grid::face_traits::secure_container* assFaces,
					 Face* facesOut[2])
{
	Grid::face_traits::secure_container	localAssFaces;
	if(!assFaces){
		Grid& grid = geom().grid();
		grid.associated_elements(localAssFaces, edge);
		assFaces = &localAssFaces;
	}

	size_t numConcernedFaces = 0;
	if(facesOut){
		for(size_t i = 0; i < assFaces->size(); ++i){
			if(is_concerned((*assFaces)[i])){
				if(numConcernedFaces < 2)
					facesOut[numConcernedFaces] = (*assFaces)[i];
				++numConcernedFaces;
			}
		}
	}
	else{
		for(size_t i = 0; i < assFaces->size(); ++i)
			numConcernedFaces += static_cast<size_t>(is_concerned((*assFaces)[i]));
	}

	return numConcernedFaces;
}

}//	end of namespace
