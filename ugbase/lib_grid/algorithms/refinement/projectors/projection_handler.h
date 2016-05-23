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
#include "common/assert.h"
#include "../refinement_projector.h"
#include "lib_grid/tools/subset_handler_interface.h"
#include "lib_grid/callbacks/selection_callbacks.h"

namespace ug{

///	Associates different projectors with individual subsets
/**	In many cases it is useful to use different RefinementProjectors on different
 * parts of a grid. To do so one may use a SubsetHandler to partition the grid into
 * different parts. This SubsetHandler can then be passed to the ProjectionHandler
 * and different RefinementProjectors can be associated with each subset through
 * 'set_projector'.
 */
class ProjectionHandler : public RefinementProjector {
public:
	ProjectionHandler () :
		m_sh (NULL)
	{
		m_defaultProjector = make_sp(new RefinementProjector);
	}

/**	Please makre sure that the given subset-handler outlives the ProjectionHandler.
 * Please note that an alternative constructor taking a smart-pointer to a
 * SubsetHandler exists.*/
	ProjectionHandler (SPIGeometry3d geometry, ISubsetHandler* psh) :
		RefinementProjector (geometry),
		m_sh (psh),
		m_linearProjector (geometry)
	{
		m_defaultProjector = make_sp(new RefinementProjector);
	}

	ProjectionHandler (SPIGeometry3d geometry, SmartPtr<ISubsetHandler> psh) :
		RefinementProjector (geometry),
		m_sh (psh.get()),
		m_spSH (psh),
		m_linearProjector (geometry)
	{
		m_defaultProjector = make_sp(new RefinementProjector);
	}

	virtual ~ProjectionHandler ()	{}

	virtual void set_geometry (SPIGeometry3d geometry)
	{
		RefinementProjector::set_geometry (geometry);
		m_linearProjector.set_geometry (geometry);
	}

///	Sets the geometry of the ProjectionHandler and of all associated projectors
/** \note	If you'd like to change the geometry of the ProjectionHandler only,
 *			please use 'set_geometry' instead.*/
	virtual void set_geometry_all (SPIGeometry3d geometry)
	{
		set_geometry(geometry);
		m_defaultProjector->set_geometry(geometry);
		for(size_t i = 0; i < m_projectors.size(); ++i){
			m_projectors[i]->set_geometry(geometry);
		}
	}

/**	Please makre sure that the given subset-handler outlives the ProjectionHandler.
 * Please note that an alternative 'set_subset_handler' exists, which takes a
 * smart-pointer to a SubsetHandler.*/
	void set_subset_handler (ISubsetHandler* psh)
	{
		if(m_spSH.valid())
			m_spSH = SmartPtr<ISubsetHandler>();
		m_sh = psh;
	}

	void set_subset_handler (SmartPtr<ISubsetHandler> psh)
	{
		m_spSH = psh;
		m_sh = psh.get();
	}

	void set_default_projector (SPRefinementProjector projector)
	{
		m_defaultProjector = projector;
	}
	
///	associate a projector with a given subsetIndex. Note that '-1' is a valid index, too.
/**	\note	if no geometry-object was set to the ProjectionHandler and if the given
 *			projector has an associated geometry-object, the handler will copy
 *			the geometry-object of the projector for its own use.*/
	void set_projector (int subsetIndex, SPRefinementProjector projector)
	{
		UG_COND_THROW (subsetIndex < -1,
					   "Bad subset-index in ProjectionHandler::set_projector: "
					   << subsetIndex << ". Indices have to be >= -1");

		projector_required(subsetIndex);

		m_projectors [subsetIndex + 1] = projector;

	//	we want to make sure that a default geometry exists
		if(geometry().invalid() && projector->geometry().valid())
			set_geometry (projector->geometry());
	}

///	associate a projector with a subset of the given name.
/** \note	a valid subset-handler has to be set to ProjectionHandler before a
 *			subset can be specified by a name. Note that it is possible to
 *			specify a subset by index before setting a SubsetHandler.
 *
 * \note	if no geometry-object was set to the ProjectionHandler and if the given
 *			projector has an associated geometry-object, the handler will copy
 *			the geometry-object of the projector for its own use.*/
	void set_projector (const char* subsetName, SPRefinementProjector projector)
	{
		UG_COND_THROW (!m_sh, "Please set a valid SubsetHandler to 'ProjectionHandler' "
					   "before calling 'set_projector' with a subset-name.");
		int si = m_sh->get_subset_index(subsetName);
		set_projector (si, projector);
	}

	size_t num_projectors () const				{return m_projectors.size() > 0 ? m_projectors.size() - 1 : 0;}
	
	SPRefinementProjector
	projector (size_t i)						{projector_required(i); return m_projectors.at(i + 1);}

	ConstSmartPtr<RefinementProjector>
	projector (size_t i)	const				{return m_projectors.at(i + 1);}

	SPRefinementProjector
	default_projector ()						{return m_defaultProjector;}

	ConstSmartPtr<RefinementProjector>
	default_projector ()	const				{return m_defaultProjector;}

////////////////////////////////////////
//	IMPLEMENTATION OF RefinementProjector
	virtual bool refinement_begins_requires_subgrid () const
	{
		for(size_t i = 0; i < m_projectors.size(); ++i){
			if(m_projectors[i]->refinement_begins_requires_subgrid())
				return true;
		}
		return false;
	}

///	prepares associated projectors for refinement
/**	If an associated projector hasn't got an associated geometry, the geometry
 * of the ProjectionHandler will automatically be assigned.*/
	virtual void refinement_begins (const ISubGrid* psg)
	{
		RefinementProjector::refinement_begins(psg);

		UG_COND_THROW(!m_sh, "Please set a valid SubsetHandler to "
					  "'ProjectionHandler' before using it during refinement.");

	//	make sure that the internal vector of projectors is at least as large
	//	as there are subsets in the associated subset handler
		projector_required((int)m_sh->num_subsets()-1);

	//	this selector will be used if we have to pass sub-grids to individual projectors
		Selector sel;

		for(size_t i = 0; i < m_projectors.size(); ++i){
			RefinementProjector& proj = *m_projectors[i];
			if(!proj.geometry().valid())
				proj.set_geometry(geometry());

			if(!proj.refinement_begins_requires_subgrid()) {
				UG_LOG("DEBUG: REF-PATH-0\n");
				proj.refinement_begins(NULL);
			}
			else{
				const ISubGrid& sg = *psg;

				if(!sel.grid())
					sel.assign_grid(geom().grid());
				
				sel.clear();
				bool selectionPerformed = false;

			//	this is a small optimization. We can either iterate over the elements
			//	of sg or of subset (i-1), i>0. Choose the one with less elements.
				GridObjectCollection subsetGoc;
				if(i > 0){
					subsetGoc = m_sh->get_grid_objects_in_subset(i - 1);
				//	todo:	this check would be more precise if only those element-types
				//			were considered which exist in both goc's
					if(subsetGoc.num_vertices() + subsetGoc.num_edges()
							+ subsetGoc.num_faces() + subsetGoc.num_volumes()
						< sg.goc().num_vertices() + sg.goc().num_edges()
							+ sg.goc().num_faces() + sg.goc().num_volumes())
					{
						if(sg.goc().num_vertices() > 0){
							lg_for_each_const(Vertex, e, subsetGoc){
								if(sg.is_contained(e))
									sel.select(e);
							}lg_end_for;
						}
						if(sg.goc().num_edges() > 0){
							lg_for_each_const(Edge, e, subsetGoc){
								if(sg.is_contained(e))
									sel.select(e);
							}lg_end_for;
						}
						if(sg.goc().num_faces() > 0){
							lg_for_each_const(Face, e, subsetGoc){
								if(sg.is_contained(e))
									sel.select(e);
							}lg_end_for;
						}
						if(sg.goc().num_volumes() > 0){
							lg_for_each_const(Volume, e, subsetGoc){
								if(sg.is_contained(e))
									sel.select(e);
							}lg_end_for;
						}

						selectionPerformed = true;
						UG_LOG("DEBUG: REF-PATH-1-1\n");
					}
				}
				
				if(!selectionPerformed){
					if(i == 0 || subsetGoc.num_vertices() > 0){
						lg_for_each_const(Vertex, e, sg.goc()){
							if(m_sh->get_subset_index(e) + 1 == (int)i)
								sel.select(e);
						}lg_end_for;
					}
					if(i == 0 || subsetGoc.num_edges() > 0){
						lg_for_each_const(Edge, e, sg.goc()){
							if(m_sh->get_subset_index(e) + 1 == (int)i)
								sel.select(e);
						}lg_end_for;
					}
					if(i == 0 || subsetGoc.num_faces() > 0){
						lg_for_each_const(Face, e, sg.goc()){
							if(m_sh->get_subset_index(e) + 1 == (int)i)
								sel.select(e);
						}lg_end_for;
					}
					if(i == 0 || subsetGoc.num_volumes() > 0){
						lg_for_each_const(Volume, e, sg.goc()){
							if(m_sh->get_subset_index(e) + 1 == (int)i)
								sel.select(e);
						}lg_end_for;
					}
					selectionPerformed = true;
					UG_LOG("DEBUG: REF-PATH-1-2\n");
				}
				
				SubGrid<IsSelected> tsg(sel.get_grid_objects(), IsSelected(sel));
				proj.refinement_begins(&tsg);
			}
		}
	}

	virtual void refinement_ends()
	{
		RefinementProjector::refinement_ends();
		for(size_t i = 0; i < m_projectors.size(); ++i)
			m_projectors[i]->refinement_ends();
	}

///	called when a new vertex was created from an old vertex.
	virtual number new_vertex (Vertex* vrt, Vertex* parent)
	{
		return handle_new_vertex (vrt, parent);
	}

///	called when a new vertex was created from an old edge.
	virtual number new_vertex (Vertex* vrt, Edge* parent)
	{
		return handle_new_vertex (vrt, parent);
	}

///	called when a new vertex was created from an old face.
	virtual number new_vertex (Vertex* vrt, Face* parent)
	{
		return handle_new_vertex (vrt, parent);
	}

///	called when a new vertex was created from an old volume.
	virtual number new_vertex (Vertex* vrt, Volume* parent)
	{
		return handle_new_vertex (vrt, parent);
	}


private:
	friend class boost::serialization::access;

	void projector_required(int index) {
	//	in order to automatically treat subset '-1' correctly, internally we add 1 
	//	to each subset-index
		if(index + 1 >= (int)m_projectors.size())
			m_projectors.resize(index + 2, m_defaultProjector);
	}

	template <class TParent>
	number handle_new_vertex (Vertex* vrt, TParent* parent)
	{
		int si = m_sh->get_subset_index(parent);

		UG_ASSERT (si + 1 < (int) m_projectors.size(),
				   "Please make sure to call 'refinement_begins' before calling "
				   "'new_vertex' for the first time for a given subset-assignment.");
		
		number ia = m_projectors [si + 1]->new_vertex (vrt, parent);
		if(ia < 1){
			vector3 p = pos (vrt);
			m_linearProjector.new_vertex (vrt, parent);
			vector3 lp = pos (vrt);
			p *= ia;
			lp *= (1. - ia);
			p += lp;
			set_pos(vrt, p);
		}
		return 1;
	}

	ISubsetHandler*								m_sh;
	SmartPtr<ISubsetHandler>					m_spSH;
	std::vector<SmartPtr<RefinementProjector> >	m_projectors;
	SmartPtr<RefinementProjector>				m_defaultProjector;
	RefinementProjector							m_linearProjector;
};

}//	end of namespace

#endif	//__H__UG_projection_handler_new
