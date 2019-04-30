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

#ifndef UG__LIB_GRID__REFINEMENT__PROJECTORS__ELLIPTIC_CYLINDER_PROJECTOR_H
#define UG__LIB_GRID__REFINEMENT__PROJECTORS__ELLIPTIC_CYLINDER_PROJECTOR_H

#include "common/math/math_vector_matrix/math_vector_functions.h"  // VecLength etc.
#include "refinement_projector.h"


namespace ug {
///	Projects new vertices onto a cylinder with an elliptic base.
/** For projection during refinement the radius property is ignored. Instead
 * the distance to the center of a newly inserted vertex is calculated
 * as the average distance of the vertices of the parent element to the center.
 * The radius property thus defaults to -1.
 *
 * You may still specify a radius. This radius can be used for auto-fitting of
 * the center and for reprojecting a set of vertices onto the sphere.
 */
class EllipticCylinderProjector
: public RefinementProjector
{
	public:
		EllipticCylinderProjector();

		EllipticCylinderProjector
		(
			const vector3& center,
			const vector3& cylAxis,
			const vector3& ellipseAxis1,
			const vector3& ellipseAxis2
		);

		EllipticCylinderProjector
		(
			const vector3& center,
			const vector3& cylAxis,
			const vector3& ellipseAxis1,
			const vector3& ellipseAxis2,
			number radius
		);

		EllipticCylinderProjector
		(
			const vector3& center,
			const vector3& cylAxis,
			const vector3& ellipseAxis1,
			const vector3& ellipseAxis2,
			number radius,
			number influenceRadius
		);

		EllipticCylinderProjector
		(
			SPIGeometry3d geometry,
			const vector3& center,
			const vector3& cylAxis,
			const vector3& ellipseAxis1,
			const vector3& ellipseAxis2,
			number radius,
			number influenceRadius
		);

		virtual ~EllipticCylinderProjector();

		void set_center(const vector3& center);
		const vector3& center() const;

		void set_cylinder_axis(const vector3& axis);
		const vector3& cylinder_axis() const;

		void set_ellipse_axis1(const vector3& axis);
		const vector3& ellipse_axis1() const;

		void set_ellipse_axis2(const vector3& axis);
		const vector3& ellipse_axis2() const;

		void set_radius(number radius);
		number radius() const;

		void set_influence_radius(number influenceRadius);
		number influence_radius() const;
	

		///	called when a new vertex was created from an old edge
		virtual number new_vertex(Vertex* vrt, Edge* parent);

		///	called when a new vertex was created from an old face
		virtual number new_vertex(Vertex* vrt, Face* parent);

		///	called when a new vertex was created from an old volume
		virtual number new_vertex(Vertex* vrt, Volume* parent);


	private:
		number radial_ellipse_coord(const vector3& v);
		number scale_point_to_radius(vector3& vIO, number r);

		template <class TElem>
		number perform_projection(Vertex* vrt, TElem* parent);

		friend class boost::serialization::access;

		template <class Archive>
		void serialize(Archive& ar, const unsigned int version)
		{
			ar & make_nvp("center", m_center);
			ar & make_nvp("cylAxis", m_cylinder_axis);
			ar & make_nvp("ellipseAxis1", m_ellipse_axis1);
			ar & make_nvp("ellipseAxis2", m_ellipse_axis2);
			ar & make_nvp("radius", m_radius);
			ar & make_nvp("influence radius", m_influenceRadius);
			UG_EMPTY_BASE_CLASS_SERIALIZATION(EllipticCylinderProjector, RefinementProjector);

			// this is only needed during load, but does not hurt during save either
			VecCross(m_ellipseNormal, m_ellipse_axis1, m_ellipse_axis2);
			m_a = VecLength(m_ellipse_axis1);
			m_b = VecLength(m_ellipse_axis2);
		}


	private:
		vector3	m_center;
		vector3	m_cylinder_axis;
		vector3	m_ellipse_axis1;
		vector3	m_ellipse_axis2;
		number	m_radius;
		number	m_influenceRadius;

		vector3	m_ellipseNormal;
		number m_a;
		number m_b;
};

} // namespace ug

#endif  // UG__LIB_GRID__REFINEMENT__PROJECTORS__ELLIPTIC_CYLINDER_PROJECTOR_H
