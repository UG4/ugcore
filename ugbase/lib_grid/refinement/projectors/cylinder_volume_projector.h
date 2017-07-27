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

#ifndef __H__UG_cylinder_volume_projector_new
#define __H__UG_cylinder_volume_projector_new

#include "common/math/misc/math_util.h"
#include "refinement_projector.h"

#include <algorithm>

namespace ug{
///	Projects new vertices to concentric cylinders.
/** For projection during refinement, the distance to the center of a newly inserted vertex
 * is calculated as follows:
 * (a) If the parent vertex center is near the axis, use it as projected point.
 * (b) Otherwise, if all parent vertices have the same distance to the axis,
 *     then project their center to the same distance.
 * (c) Otherwise, if the parent vertices are at two distinct distances from the axis, then
 *     use their center and project to the average of these two distances.
 * Other cases cannot (?) happen if one of the cases is true for all the parent elements.
 */
class CylinderVolumeProjector : public RefinementProjector {
public:
	CylinderVolumeProjector ()
	: m_center(0, 0, 0), m_axis(0, 0, 1)
	{}
	
	CylinderVolumeProjector (const vector3& center, const vector3& axis)
	: m_center(center), m_axis(axis)
	{}

	CylinderVolumeProjector(SPIGeometry3d geometry, const vector3& center, const vector3& axis)
	: RefinementProjector(geometry), m_center(center), m_axis(axis)
	{}

	virtual ~CylinderVolumeProjector() {}

	void set_center (const vector3& center) {m_center = center;}
	const vector3& center () const {return m_center;}

	void set_axis (const vector3& axis) {m_axis = axis;}
	const vector3& axis () const {return m_axis;}

	///	called when a new vertex was created from an old edge.
	virtual number new_vertex(Vertex* vrt, Edge* parent)
	{return perform_projection(vrt, parent);}

	///	called when a new vertex was created from an old face.
	virtual number new_vertex(Vertex* vrt, Face* parent)
	{return perform_projection(vrt, parent);}

	///	called when a new vertex was created from an old volume.
	virtual number new_vertex(Vertex* vrt, Volume* parent)
	{return perform_projection(vrt, parent);}

private:
	struct CmpAlmostEqual
	{
		CmpAlmostEqual(number thresh) : t(thresh) {};

		bool operator()(const number& a, const number& b)
		{return fabs(a-b) < t;}

		number t;
	};

	template <class TElem>
	number perform_projection(Vertex* vrt, TElem* parent)
	{
		// calculate the new position by linear interpolation
		// and project that point onto the cylinder.
		typename TElem::ConstVertexArray vrts = parent->vertices();
		size_t numVrts = parent->num_vertices();

		if (numVrts == 0)
		{
			set_pos(vrt, vector3(0, 0, 0));
			return 1;
		}

		UG_COND_THROW(numVrts == 1, "How can that be?");

		std::vector<number> vDist(numVrts);
		vector3 parentCenter(0, 0, 0);
		number avDist = 0.0;
		for (size_t i = 0; i < numVrts; ++i)
		{
			const vector3& p = pos(vrts[i]);
			avDist += vDist[i] = DistancePointToRay(p, m_center, m_axis);
			parentCenter += p;
		}

		// sort distances
		std::sort(vDist.begin(), vDist.end());
		std::vector<number>::iterator newEnd = std::unique(vDist.begin(), vDist.end(),
		    CmpAlmostEqual(0.01 * (vDist[vDist.size()-1] - vDist[0])));
		vDist.erase(newEnd, vDist.end());

		// implement (a)-(c)
		avDist /= (number)numVrts;
		VecScale(parentCenter, parentCenter, 1.0 / (number)numVrts);
		vector3 proj, v;
		ProjectPointToRay(proj, parentCenter, m_center, m_axis);
		VecSubtract(v, parentCenter, proj);
		number len = VecLength(v);

		// (a)
		if (len < 0.01*avDist)
		{
            set_pos(vrt, parentCenter);
            return 1;
		}

		// (b)
		if (vDist.size() == 1)
		{
		    VecScale(v, v, vDist[0]/len);
            proj += v;
            set_pos(vrt, proj);
            return 1;
		}

		// (c)
		if (vDist.size() == 2)
		{
		    VecScale(v, v, 0.5*(vDist[0] + vDist[1]) / len);
            proj += v;
            set_pos(vrt, proj);
            return 1;
		}

		UG_THROW("This should not be reached.");

		return 1;
	}


	friend class boost::serialization::access;

	template <class Archive>
	void serialize( Archive& ar, const unsigned int version)
	{
		ar & make_nvp("center", m_center);
		ar & make_nvp("axis", m_axis);
		UG_EMPTY_BASE_CLASS_SERIALIZATION(CylinderVolumeProjector, RefinementProjector);
	}

	vector3	m_center;
	vector3	m_axis;
};

}//	end of namespace

#endif	//__H__UG_cylinder_volume_projector_new
