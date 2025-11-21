/*
 * Copyright (c) 2016:  G-CSC, Goethe University Frankfurt
 * Author: Markus Breit
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

#include "elliptic_cylinder_projector.h"

#include "common/math/misc/math_util.h"
#include "refinement_projector.h"


namespace ug {


EllipticCylinderProjector::EllipticCylinderProjector()
: m_center(0, 0, 0),
  m_cylinder_axis(0, 0, 1),
  m_ellipse_axis1(1, 0, 0),
  m_ellipse_axis2(0, 1, 0),
  m_radius(-1),
  m_influenceRadius(-1),
  m_ellipseNormal(0, 0, 1),
  m_a(1.0),
  m_b(1.0)
{}


EllipticCylinderProjector::EllipticCylinderProjector
(
	const vector3& center,
	const vector3& cylAxis,
	const vector3& ellipseAxis1,
	const vector3& ellipseAxis2
)
: m_center(center),
  m_cylinder_axis(cylAxis),
  m_ellipse_axis1(ellipseAxis1),
  m_ellipse_axis2(ellipseAxis2),
  m_radius(-1),
  m_influenceRadius(-1),
  m_a(VecLength(ellipseAxis1)),
  m_b(VecLength(ellipseAxis2))
{
	VecCross(m_ellipseNormal, m_ellipse_axis1, m_ellipse_axis2);
}


EllipticCylinderProjector::EllipticCylinderProjector
(
	const vector3& center,
	const vector3& cylAxis,
	const vector3& ellipseAxis1,
	const vector3& ellipseAxis2,
	number radius
)
: m_center(center),
  m_cylinder_axis(cylAxis),
  m_ellipse_axis1(ellipseAxis1),
  m_ellipse_axis2(ellipseAxis2),
  m_radius(radius),
  m_influenceRadius(-1),
  m_a(VecLength(ellipseAxis1)),
  m_b(VecLength(ellipseAxis2))
{
	UG_LOGN("Constructor5");
	VecCross(m_ellipseNormal, m_ellipse_axis1, m_ellipse_axis2);
}


EllipticCylinderProjector::EllipticCylinderProjector
(
	const vector3& center,
	const vector3& cylAxis,
	const vector3& ellipseAxis1,
	const vector3& ellipseAxis2,
	number radius,
	number influenceRadius
)
: m_center(center),
  m_cylinder_axis(cylAxis),
  m_ellipse_axis1(ellipseAxis1),
  m_ellipse_axis2(ellipseAxis2),
  m_radius(radius),
  m_influenceRadius(influenceRadius),
  m_a(VecLength(ellipseAxis1)),
  m_b(VecLength(ellipseAxis2))
{
	VecCross(m_ellipseNormal, m_ellipse_axis1, m_ellipse_axis2);
}


EllipticCylinderProjector::EllipticCylinderProjector
(
	SPIGeometry3d geometry,
	const vector3& center,
	const vector3& cylAxis,
	const vector3& ellipseAxis1,
	const vector3& ellipseAxis2,
	number radius,
	number influenceRadius
)
: RefinementProjector(geometry),
  m_center(center),
  m_cylinder_axis(cylAxis),
  m_ellipse_axis1(ellipseAxis1),
  m_ellipse_axis2(ellipseAxis2),
  m_radius(radius),
  m_influenceRadius(influenceRadius),
  m_a(VecLength(ellipseAxis1)),
  m_b(VecLength(ellipseAxis2))
{
	VecCross(m_ellipseNormal, m_ellipse_axis1, m_ellipse_axis2);
}



void EllipticCylinderProjector::set_center(const vector3& center)
{
	m_center = center;
}

const vector3& EllipticCylinderProjector::center() const
{
	return m_center;
}


void EllipticCylinderProjector::set_cylinder_axis(const vector3& axis)
{
	m_cylinder_axis = axis;
}

const vector3& EllipticCylinderProjector::cylinder_axis() const
{
	return m_cylinder_axis;
}


void EllipticCylinderProjector::set_ellipse_axis1(const vector3& axis)
{
	m_ellipse_axis1 = axis;
	VecCross(m_ellipseNormal, m_ellipse_axis1, m_ellipse_axis2);
	m_a = VecLength(axis);
}

const vector3& EllipticCylinderProjector::ellipse_axis1() const
{
	return m_ellipse_axis1;
}


void EllipticCylinderProjector::set_ellipse_axis2(const vector3& axis)
{
	m_ellipse_axis2 = axis;
	VecCross(m_ellipseNormal, m_ellipse_axis1, m_ellipse_axis2);
	m_b = VecLength(axis);
}

const vector3& EllipticCylinderProjector::ellipse_axis2() const
{
	return m_ellipse_axis2;
}


void EllipticCylinderProjector::set_radius(number radius)
{
	m_radius = radius;
}

number EllipticCylinderProjector::radius() const
{
	return m_radius;
}


void EllipticCylinderProjector::set_influence_radius(number influenceRadius)
{
	m_influenceRadius = influenceRadius;
}

number EllipticCylinderProjector::influence_radius() const
{
	return m_influenceRadius;
}



number EllipticCylinderProjector::new_vertex(Vertex* vrt, Edge* parent)
{
	return perform_projection(vrt, parent);
}

number EllipticCylinderProjector::new_vertex(Vertex* vrt, Face* parent)
{
	return perform_projection(vrt, parent);
}

number EllipticCylinderProjector::new_vertex(Vertex* vrt, Volume* parent)
{
	return perform_projection(vrt, parent);
}



number EllipticCylinderProjector::radial_ellipse_coord(const vector3& v)
{
	// project coordinates along cylinder axis to ellipse plane
	vector3 proj2Ellipse;
	number dummy;
	RayPlaneIntersection(proj2Ellipse, dummy, v, m_cylinder_axis, m_center, m_ellipseNormal);

	// get x and y coordinates
	const number x = VecDot(m_ellipse_axis1, proj2Ellipse) / m_a;
	const number y = VecDot(m_ellipse_axis2, proj2Ellipse) / m_b;

	return sqrt(x*x/(m_a*m_a) + y*y/(m_b*m_b));
}


number EllipticCylinderProjector::scale_point_to_radius(vector3& vIO, number r)
{
	// project coordinates along cylinder axis to ellipse plane
	vector3 proj2Ellipse;
	number dummy;
	RayPlaneIntersection(proj2Ellipse, dummy, vIO, m_cylinder_axis, m_center, m_ellipseNormal);

	VecSubtract(vIO, vIO, proj2Ellipse);

	// current radius
	const number x = VecDot(m_ellipse_axis1, proj2Ellipse) / m_a;
	const number y = VecDot(m_ellipse_axis2, proj2Ellipse) / m_b;
	const number rcur = sqrt(x*x/(m_a*m_a) + y*y/(m_b*m_b));

	// scale and undo projection
	if (rcur > SMALL * r)
		VecScaleAdd(vIO, 1.0, vIO, r / rcur, proj2Ellipse);
	else
		// if current position is in the center of the cylinder, leave it there
		VecAdd(vIO, vIO, m_center);

	return rcur;
}


template <typename TElem>
number EllipticCylinderProjector::perform_projection(Vertex* vrt, TElem* parent)
{
	// General method:
	// The "radius" of all parent vertices is averaged
	// and used for the radius of the projected child.
	// The (Cartesian) coordinates of all parents are averaged.
	// The resulting coordinates are scaled around the axis
	// within the ellipse plane to reach the previously determined radius.
	// TODO: It might be better to also average arc lengths and axial coords of parents
	//       and project the child to averaged radius, arc length and axial coords.

	// average parent vertex positions and radii
	typename TElem::ConstVertexArray vrts = parent->vertices();
	const size_t numVrts = parent->num_vertices();

	if (numVrts == 0)
	{
		set_pos(vrt, vector3(0, 0, 0));
		return 1;
	}

	number avgR = 0;
	vector3 proj(0, 0, 0);
	for (size_t i = 0; i < numVrts; ++i)
	{
		const vector3& p = pos(vrts[i]);
		avgR += radial_ellipse_coord(p);
		proj += p;
	}
	avgR /= numVrts;
	VecScale(proj, proj, 1.0 / numVrts);

	// move averaged position to new radius
	const number curR = scale_point_to_radius(proj, avgR);
	set_pos(vrt, proj);

	if (m_influenceRadius > 0)
	{
		if (m_radius > m_influenceRadius)
		{
			const number dist = m_radius - m_influenceRadius;
			return clip<number>((curR - m_influenceRadius) / dist, 0, 1);
		}
		else if (m_radius >= 0)
		{
			const number dist = m_influenceRadius - m_radius;
			if (dist > 0)
				return clip<number>(1 - (curR - m_radius) / dist, 0, 1);
			return curR < m_radius ? 1 : 0;
		}
		else
			return clip<number>(1 - curR / m_influenceRadius, 0, 1);
	}

	return 1;
}


} // namespace ug
