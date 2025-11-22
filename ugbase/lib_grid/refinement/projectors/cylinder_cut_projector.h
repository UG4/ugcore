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

#ifndef __H__UG_cylinder_cut_projector
#define __H__UG_cylinder_cut_projector

#include "refinement_projector.h"

namespace ug{

/// Refines linearly except for when a refined edge intersects the given cylinder
/**	Refines linearly except for when a refined edge intersects the given cylinder.
 * the new vertex will be placed at the intersection of the edge with the cylinder
 * in this case.
 *
 * \note	This projector is not commonly used during grid adaption. It serves
 *			very special purposes for specialized algorithms. Used without
 *			specialized meshing operations it may lead to bad elements.
 */
class CylinderCutProjector : public RefinementProjector {
public:
	CylinderCutProjector () :
		m_center (0, 0, 0),
		m_axis (0, 0, 1),
		m_radius (1)
	{}
	
	CylinderCutProjector (const vector3& center,
						  const vector3& axis,
					 	  number radius) :
		m_center (center),
		m_axis (axis),
		m_radius (radius)
	{}

/**	\sa ug::RefinementProjector::RefinementProjector*/
	CylinderCutProjector (SPIGeometry3d geometry,
						  const vector3& center,
						  const vector3& axis,
						  number radius) :
		RefinementProjector (geometry),
		m_center (center),
		m_axis (axis),
		m_radius (radius)
	{}



	virtual ~CylinderCutProjector () = default;

	void set_center (const vector3& center)		{m_center = center;}
	const vector3& center () const				{return m_center;}

	void set_axis (const vector3& axis)			{m_axis = axis;}
	const vector3& axis () const				{return m_axis;}

	void set_radius (number radius)				{m_radius = radius;}
	number radius () const						{return m_radius;}

///	called when a new vertex was created from an old edge.
	virtual number new_vertex(Vertex* vrt, Edge* parent)
	{
		number t0, t1;
		vector3 from = pos(parent->vertex(0));
		vector3 dir;
		VecSubtract(dir, pos(parent->vertex(1)), from);

		if(RayCylinderIntersection(t0, t1, from, dir, m_center, m_axis, m_radius))
		{
		//	if there are two intersections with parameters between 0 and 1,
		//	we'll return their median.
			bool t0IsFine = (t0 >= 0) && (t0 <= 1);
			bool t1IsFine = (t1 >= 0) && (t1 <= 1);
			if(t0IsFine){
				if(t1IsFine) {
					vector3 v0, v1, p;
					VecScaleAdd(v0, 1., from, t0, dir);
					VecScaleAdd(v1, 1., from, t1, dir);
					VecScaleAdd(p, 0.5, v0, 0.5, v1);
					set_pos(vrt, p);
				}
				else {
					vector3 p;
					VecScaleAdd(p, 1., from, t0, dir);
					set_pos(vrt, p);
				}
				return 1;
			}
			else if(t1IsFine) {
				vector3 p;
				VecScaleAdd(p, 1., from, t1, dir);
				set_pos(vrt, p);
				return 1;
			}
		}

		RefinementProjector::new_vertex(vrt, parent);
		return 1;
	}

private:
	friend class boost::serialization::access;

	template <typename Archive>
	void serialize( Archive& ar, const unsigned int version)
	{
		ar & make_nvp("center", m_center);
		ar & make_nvp("axis", m_axis);
		ar & make_nvp("radius", m_radius);
		UG_EMPTY_BASE_CLASS_SERIALIZATION(CylinderCutProjector, RefinementProjector);
	}

	vector3	m_center;
	vector3	m_axis;
	number	m_radius;
};

}//	end of namespace

#endif