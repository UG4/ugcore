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

#ifndef __H__UG_refinement_projector
#define __H__UG_refinement_projector

#include "common/boost_serialization.h"
#include "common/error.h"
#include "lib_grid/grid/geometry.h"
#include "lib_grid/grid/sub_grid.h"

namespace ug{

///	Adjusts vertex coordinates during refinement
/** The refinement projector serves as a base class for other refinement projectors
 * and serves as a simple linear projector at the same time. I.E. new vertices are
 * always placed at the center of their parent element.
 *
 * Internally, refinement projectors are implemented for 3d coordinates. Geometries
 * of lower dimension can still be treated since they are wrapped in a
 * IGeometry3d class. While this adds a small overhead it simplifies implementation
 * and usability considerably, since no additional template dependencies have to
 * be considered.*/
class RefinementProjector {
public:
	RefinementProjector ()			{}
	
	RefinementProjector (SPIGeometry3d geometry) :
		m_geometry (geometry)
	{}

	virtual ~RefinementProjector ()	{}

	virtual void set_geometry (SPIGeometry3d geometry)
	{
		m_geometry = geometry;
	}

	virtual SPIGeometry3d geometry () const				{return m_geometry;}

///	called before refinement begins
/**	The specified sub-grid contains all elements that will be refined.*/
	virtual void refinement_begins(const ISubGrid& sg)
	{
		UG_COND_THROW(m_geometry.invalid(),
					  "Please set a valid geometry to RefinementProjectors"
					  "before using them during refinement.");
	}

///	called when refinement is done
/**	The specified sub-grid contains all elements which resulted from refinement.*/
	virtual void refinement_ends(const ISubGrid& sg)	{}

///	called when a new vertex was created from an old vertex.
	virtual void new_vertex(Vertex* vrt, Vertex* parent)
	{
		set_pos(vrt, pos(parent));
	}

///	called when a new vertex was created from an old edge.
	virtual void new_vertex(Vertex* vrt, Edge* parent)
	{
		set_pos(vrt, geom().element_center(parent));
	}

///	called when a new vertex was created from an old face.
	virtual void new_vertex(Vertex* vrt, Face* parent)
	{
		set_pos(vrt, geom().element_center(parent));
	}

///	called when a new vertex was created from an old volume.
	virtual void new_vertex(Vertex* vrt, Volume* parent)
	{
		set_pos(vrt, geom().element_center(parent));
	}

protected:
	vector3 pos (Vertex* v) const
	{
		// UG_COND_THROW(m_geometry.invalid(),
		// 			  "Invalid Geometry in IRefinementProjector. Please make "
		// 			  "sure to set a valid adaptor through 'set_geometry_adaptor'.");
		return m_geometry->pos(v);
	}

	void set_pos(Vertex* v, const vector3& p)
	{
		// UG_COND_THROW(m_geometry.invalid(),
		// 			  "Invalid Geometry in IRefinementProjector. Please make "
		// 			  "sure to set a valid adaptor through 'set_geometry_adaptor'.");
		m_geometry->set_pos(v, p);
	}

	IGeometry3d& geom ()				{return *m_geometry;}
	const IGeometry3d& geom () const	{return *m_geometry;}

private:
	friend class boost::serialization::access;

	template <class Archive>
	void serialize( Archive& ar, const unsigned int version)
	{
	}

	SPIGeometry3d	m_geometry;
};


typedef SmartPtr<RefinementProjector>	SPRefinementProjector;


}//	end of namespace

#endif	//__H__UG_refinement_projector
