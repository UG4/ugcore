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

#ifndef __H__UG_sphere_projector_new
#define __H__UG_sphere_projector_new

#include "refinement_projector.h"

namespace ug{

///	Projects new vertices onto a sphere during refinement
/** For projection during refinement the radius property is ignored. Instead
 * the distance to the center of a newly inserted vertex is calculated
 * as the average distance of the vertices of the parent element to the center.
 * The radius property thus defaults to -1.
 *
 * You may still specify a radius. This radius can be used for auto-fitting of
 * the center and for reprojecting a set of vertices onto the sphere.
 */
class SphereProjector : public RefinementProjector {
public:
	SphereProjector () :
		m_center (0, 0, 0),
		m_radius (-1),
		m_influenceRadius (-1)
	{}
	
	SphereProjector (const vector3& center) :
		m_center (center),
		m_radius (-1),
		m_influenceRadius (-1)
	{}

	SphereProjector (const vector3& center,
						number radius) :
		m_center (center),
		m_radius (radius),
		m_influenceRadius (-1)
	{}

	SphereProjector (const vector3& center,
						number radius,
						number influenceRadius) :
		m_center (center),
		m_radius (radius),
		m_influenceRadius (influenceRadius)
	{}

/**	\sa ug::RefinementProjector::RefinementProjector*/
	SphereProjector (SPIGeometry3d geometry,
						const vector3& center) :
		RefinementProjector (geometry),
		m_center (center),
		m_radius (-1),
		m_influenceRadius (-1)
	{}

	/**	\sa ug::RefinementProjector::RefinementProjector*/
	SphereProjector (SPIGeometry3d geometry,
						const vector3& center,
						number radius) :
		RefinementProjector (geometry),
		m_center (center),
		m_radius (radius),
		m_influenceRadius (-1)
	{}

/**	\sa ug::RefinementProjector::RefinementProjector*/
	template <typename TGeomProvider>
	SphereProjector (const TGeomProvider& geometry,
						const vector3& center,
						number radius,
					 	number influenceRadius) :
		RefinementProjector (geometry),
		m_center (center),
		m_radius (radius),
		m_influenceRadius (influenceRadius)
	{}

	virtual ~SphereProjector ()					{}

	void set_center (const vector3& center)		{m_center = center;}
	const vector3& center () const				{return m_center;}

	void set_radius (number radius)				{m_radius = radius;}
	number radius () const						{return m_radius;}

	void set_influence_radius (number influenceRadius)	{m_influenceRadius = influenceRadius;}
	number influence_radius () const					{return m_influenceRadius;}

///	called when a new vertex was created from an old edge.
	virtual number new_vertex(Vertex* vrt, Edge* parent)
	{
		return perform_projection(vrt, parent);
	}

///	called when a new vertex was created from an old face.
	virtual number new_vertex(Vertex* vrt, Face* parent)
	{
		return perform_projection(vrt, parent);
	}

///	called when a new vertex was created from an old volume.
	virtual number new_vertex(Vertex* vrt, Volume* parent)
	{
		return perform_projection(vrt, parent);
	}

private:

	template <typename TElem>
	number perform_projection(Vertex* vrt, TElem* parent)
	{
	//	first calculate the average distance of corners of the parent to the
	//	parents center
		typename TElem::ConstVertexArray vrts = parent->vertices();
		size_t numVrts = parent->num_vertices();

		if(numVrts == 0){
			set_pos(vrt, vector3(0, 0, 0));
			return 1;
		}

		number avDist = 0;
		vector3 parentCenter(0, 0, 0);

		for(size_t i = 0; i < numVrts; ++i){
			vector3 p = pos(vrts[i]);
			avDist += VecDistance(p, m_center);
			parentCenter += p;
		}

		avDist /= (number)numVrts;
		VecScale(parentCenter, parentCenter, 1. / (number)numVrts);

	//	calculate projection
		vector3 proj;
		VecSubtract(proj, parentCenter, m_center);
		number len = VecLength(proj);
		if(len > SMALL * avDist){	// if avDist is very small, len may be small, too
			VecScale(proj, proj, avDist / len);
			proj += m_center;
			set_pos(vrt, proj);
		}
		else
			set_pos(vrt, parentCenter);

		if(m_influenceRadius > 0) {
			if(m_radius > m_influenceRadius){
				const number dist = m_radius - m_influenceRadius;
				if(dist > 0)
					return clip<number>((len - m_influenceRadius) / dist, 0, 1);
				return len > m_radius ? 1 : 0;
			}
			else if(m_radius >= 0){
				const number dist = m_influenceRadius - m_radius;
				if(dist > 0)
					return clip<number>(1 - (len - m_radius) / dist, 0, 1);
				return len < m_radius ? 1 : 0;
			}
			else
				return clip<number>(1 - len / m_influenceRadius, 0, 1);
		}

		return 1;
	}


	friend class boost::serialization::access;

	template <typename Archive>
	void serialize( Archive& ar, const unsigned int version)
	{
		ar & make_nvp("center", m_center);
		ar & make_nvp("radius", m_radius);
		ar & make_nvp("influence radius", m_influenceRadius);
		UG_EMPTY_BASE_CLASS_SERIALIZATION(SphereProjector, RefinementProjector);
	}

	vector3	m_center;
	number	m_radius;
	number	m_influenceRadius;
};


}//	end of namespace

#endif