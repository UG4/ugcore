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

#ifndef __H__UG_projection_handler_new
#define __H__UG_projection_handler_new

#include <vector>
#include <algorithm>
#include "common/assert.h"
#include "refinement_projector.h"
#include "lib_grid/tools/subset_handler_interface.h"
#include "lib_grid/callbacks/selection_callbacks.h"
#include "lib_grid/callbacks/subset_callbacks.h"

namespace ug{

///	Associates different projectors with individual subsets
/**	In many cases it is useful to use different RefinementProjectors on different
 * parts of a grid. To do so one may use a SubsetHandler to partition the grid into
 * different parts. This SubsetHandler can then be passed to the ProjectionHandler
 * and different RefinementProjectors can be associated with each subset through
 * 'set_projector'.
 */
UG_API class ProjectionHandler : public RefinementProjector {
public:
	ProjectionHandler ();

/**	Please makre sure that the given subset-handler outlives the ProjectionHandler.
 * Please note that an alternative constructor taking a smart-pointer to a
 * SubsetHandler exists.
 * \sa ug::RefinementProjector::RefinementProjector*/
	ProjectionHandler (ISubsetHandler* psh);

/**	\sa ug::RefinementProjector::RefinementProjector*/
	ProjectionHandler (SmartPtr<ISubsetHandler> psh);

/**	Please makre sure that the given subset-handler outlives the ProjectionHandler.
 * Please note that an alternative constructor taking a smart-pointer to a
 * SubsetHandler exists.
 * \sa ug::RefinementProjector::RefinementProjector*/
	ProjectionHandler (SPIGeometry3d geometry,
					   ISubsetHandler* psh);

/**	\sa ug::RefinementProjector::RefinementProjector*/
	template <typename TGeomProvider>
	ProjectionHandler (const TGeomProvider& geometry,
					   SmartPtr<ISubsetHandler> psh) :
		RefinementProjector (geometry),
		m_sh (psh.get()),
		m_spSH (psh)
	{
		m_defaultProjector = make_sp(new RefinementProjector);
	}

	virtual ~ProjectionHandler () = default;;

	void clear();

	virtual void set_geometry (SPIGeometry3d geometry);

///	Sets the geometry of the ProjectionHandler and of all associated projectors
/** \note	If you'd like to change the geometry of the ProjectionHandler only,
 *			please use 'set_geometry' instead.*/
	virtual void set_geometry_all (SPIGeometry3d geometry);

/**	Please make sure that the given subset-handler outlives the ProjectionHandler.
 * Please note that an alternative 'set_subset_handler' exists, which takes a
 * smart-pointer to a SubsetHandler.*/
	void set_subset_handler (ISubsetHandler* psh);

	void set_subset_handler (SmartPtr<ISubsetHandler> psh);


///	return the subset handler that the projection handler is based on
	inline const ISubsetHandler* subset_handler() const			{return m_sh;}

	void set_default_projector (SPRefinementProjector projector);
	
///	associate a projector with a given subsetIndex. Note that '-1' is a valid index, too.
/**	\note	if no geometry-object was set to the ProjectionHandler and if the given
 *			projector has an associated geometry-object, the handler will copy
 *			the geometry-object of the projector for its own use.*/
	void set_projector (int subsetIndex, SPRefinementProjector projector);

///	associate a projector with a subset of the given name.
/** \note	a valid subset-handler has to be set to ProjectionHandler before a
 *			subset can be specified by a name. Note that it is possible to
 *			specify a subset by index before setting a SubsetHandler.
 *
 * \note	if no geometry-object was set to the ProjectionHandler and if the given
 *			projector has an associated geometry-object, the handler will copy
 *			the geometry-object of the projector for its own use.*/
	void set_projector (const char* subsetName, SPRefinementProjector projector);

	inline size_t num_projectors () const		{return m_projectors.size() > 0 ? m_projectors.size() - 1 : 0;}
	
	inline SPRefinementProjector
	projector (size_t i)						{projector_required((int)i); return m_projectors.at(i + 1);}

	inline ConstSmartPtr<RefinementProjector>
	projector (size_t i)	const				{return m_projectors.at(i + 1);}

	inline SPRefinementProjector
	default_projector ()						{return m_defaultProjector;}

	inline ConstSmartPtr<RefinementProjector>
	default_projector ()	const				{return m_defaultProjector;}

////////////////////////////////////////
//	IMPLEMENTATION OF RefinementProjector
	virtual bool refinement_begins_requires_subgrid () const;

///	prepares associated projectors for refinement
/**	If an associated projector hasn't got an associated geometry, the geometry
 * of the ProjectionHandler will automatically be assigned.*/
	virtual void refinement_begins (const ISubGrid* psg);

	virtual void refinement_ends();

///	called when a new vertex was created from an old vertex.
	virtual number new_vertex (Vertex* vrt, Vertex* parent);

///	called when a new vertex was created from an old edge.
	virtual number new_vertex (Vertex* vrt, Edge* parent);

///	called when a new vertex was created from an old face.
	virtual number new_vertex (Vertex* vrt, Face* parent);

///	called when a new vertex was created from an old volume.
	virtual number new_vertex (Vertex* vrt, Volume* parent);


private:
	friend class boost::serialization::access;

	void projector_required(int index);

	template <typename TParent>
	number handle_new_vertex (Vertex* vrt, TParent* parent);

	ISubsetHandler*								m_sh;
	SmartPtr<ISubsetHandler>					m_spSH;
	std::vector<SmartPtr<RefinementProjector> >	m_projectors;
	SmartPtr<RefinementProjector>				m_defaultProjector;
};

using SPProjectionHandler = SmartPtr<ProjectionHandler>;

}//	end of namespace

#endif