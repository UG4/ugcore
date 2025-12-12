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

#ifndef __H__UG_plane_cut_projector
#define __H__UG_plane_cut_projector

#include "refinement_projector.h"

namespace ug {

///	calculates new positions by cutting parent edges with a plane
/**	For each edge the intersection of the edge with the initially
 *	given plane is calculated and used as new point. Vertices
 *	created on other geometric objects are treated as in the linear case.
 */
class PlaneCutProjector : public RefinementProjector {
public:
	PlaneCutProjector () :
		m_p (0, 0, 0),
		m_n (0, 0, 1)
	{}
	
	PlaneCutProjector (const vector3& position,
						const vector3& normal) :
		m_p (position),
		m_n (normal)
	{}

/**	\sa ug::RefinementProjector::RefinementProjector*/
	PlaneCutProjector (SPIGeometry3d geometry,
					   const vector3& position,
					   const vector3& normal) :
		RefinementProjector (geometry),
		m_p (position),
		m_n (normal)
	{}

	~PlaneCutProjector () override = default;

	void set_position (const vector3& position) {m_p = position;}
	const vector3& position () const {return m_p;}

	void set_normal (const vector3& normal) {m_n = normal;}
	[[nodiscard]] const vector3& normal () const {return m_n;}

///	called when a new vertex was created from an old edge.
	number new_vertex(Vertex* vrt, Edge* parent) override {
		number t;
		vector3 v;
		vector3 dir;
		VecSubtract(dir, pos(parent->vertex(1)), pos(parent->vertex(0)));
		
		if(RayPlaneIntersection(v, t, pos(parent->vertex(0)), dir, m_p, m_n)) {
			set_pos(vrt, v);
		}
		else{
			RefinementProjector::new_vertex(vrt, parent);
		}
		return 1;
	}

private:
	friend class boost::serialization::access;

	template <typename Archive>
	void serialize( Archive& ar, const unsigned int version)
	{
		ar & make_nvp("position", m_p);
		ar & make_nvp("normal", m_n);
		UG_EMPTY_BASE_CLASS_SERIALIZATION(PlaneCutProjector, RefinementProjector);
	}

	vector3	m_p;
	vector3	m_n;
};

}//	end of namespace

#endif