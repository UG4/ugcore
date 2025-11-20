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

#ifndef __H__UG_smooth_projector
#define __H__UG_smooth_projector

#include "refinement_projector.h"

namespace ug{


///	Smoothes vertices during refinement
class SmoothProjector : public RefinementProjector {
public:
	SmoothProjector () :
		m_iterations (10),
		m_changeRate (0.25)
	{}
	
	SmoothProjector (int iterations,
					number changeRate) :
		m_iterations (iterations),
		m_changeRate (changeRate)
	{}

/**	\sa ug::RefinementProjector::RefinementProjector*/
	SmoothProjector (SPIGeometry3d geometry,
					int iterations,
					number changeRate) :
		RefinementProjector (geometry),
		m_iterations (iterations),
		m_changeRate (changeRate)
	{}

	int iterations () const	{return m_iterations;}
	void set_iterations (int iterations)	{m_iterations = iterations;}

	number change_rate () const	{return m_changeRate;}
	void set_change_rate (number changeRate)	{m_changeRate = changeRate;}

///	called before refinement begins
	virtual void refinement_begins(const ISubGrid* sg)
	{
		RefinementProjector::refinement_begins(sg);
		m_newVrts.clear();
	}

///	called when refinement is done.
/**	The actual smoothing is performed here*/
	virtual void refinement_ends();

///	called when a new vertex was created from an old edge.
	virtual number new_vertex(Vertex* vrt, Edge* parent)
	{
		m_newVrts.push_back(vrt);
		return RefinementProjector::new_vertex(vrt, parent);
	}

///	called when a new vertex was created from an old face.
	virtual number new_vertex(Vertex* vrt, Face* parent)
	{
		m_newVrts.push_back(vrt);
		return RefinementProjector::new_vertex(vrt, parent);
	}

///	called when a new vertex was created from an old volume.
	virtual number new_vertex(Vertex* vrt, Volume* parent)
	{
		m_newVrts.push_back(vrt);
		return RefinementProjector::new_vertex(vrt, parent);
	}

private:
	friend class boost::serialization::access;

	template <typename Archive>
	void serialize( Archive& ar, const unsigned int version)
	{
		ar & make_nvp("iterations", m_iterations);
		ar & make_nvp("change rate", m_changeRate);
		UG_EMPTY_BASE_CLASS_SERIALIZATION(SmoothProjector, RefinementProjector);
	}

	std::vector<Vertex*> m_newVrts;
	int		m_iterations;
	number	m_changeRate;
};

}//	end of namespace

#endif