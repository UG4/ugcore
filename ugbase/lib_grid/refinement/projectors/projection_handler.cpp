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

#include "projection_handler.h"

namespace ug{

ProjectionHandler::
ProjectionHandler () :
	m_sh (nullptr)
{
	m_defaultProjector = make_sp(new RefinementProjector);
}

ProjectionHandler::
ProjectionHandler (ISubsetHandler* psh) :
	m_sh (psh)
{
	m_defaultProjector = make_sp(new RefinementProjector);
}

ProjectionHandler::
ProjectionHandler (SmartPtr<ISubsetHandler> psh) :
	m_sh (psh.get()),
	m_spSH (psh)
{
	m_defaultProjector = make_sp(new RefinementProjector);
}

ProjectionHandler::
ProjectionHandler (SPIGeometry3d geometry,
				   ISubsetHandler* psh) :
	RefinementProjector (geometry),
	m_sh (psh)
{
	m_defaultProjector = make_sp(new RefinementProjector);
}

ProjectionHandler::
~ProjectionHandler ()	{}

void ProjectionHandler::
clear()
{
	m_projectors.clear();
	m_defaultProjector = make_sp(new RefinementProjector);
}

void ProjectionHandler::
set_geometry (SPIGeometry3d geometry)
{
	RefinementProjector::set_geometry (geometry);
	m_defaultProjector->set_geometry(geometry);
}

void ProjectionHandler::
set_geometry_all (SPIGeometry3d geometry)
{
	set_geometry(geometry);
	m_defaultProjector->set_geometry(geometry);
	for(size_t i = 0; i < m_projectors.size(); ++i){
		if(m_projectors[i].valid())
			m_projectors[i]->set_geometry(geometry);
	}
}

void ProjectionHandler::
set_subset_handler (ISubsetHandler* psh)
{
	if(m_spSH.valid())
		m_spSH = SmartPtr<ISubsetHandler>();
	m_sh = psh;
}

void ProjectionHandler::
set_subset_handler (SmartPtr<ISubsetHandler> psh)
{
	m_spSH = psh;
	m_sh = psh.get();
}


void ProjectionHandler::
set_default_projector (SPRefinementProjector projector)
{
	m_defaultProjector = projector;
}

void ProjectionHandler::
set_projector (int subsetIndex, SPRefinementProjector projector)
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

void ProjectionHandler::
set_projector (const char* subsetName, SPRefinementProjector projector)
{
	UG_COND_THROW (!m_sh, "Please set a valid SubsetHandler to 'ProjectionHandler' "
				   "before calling 'set_projector' with a subset-name.");
	int si = m_sh->get_subset_index(subsetName);
	set_projector (si, projector);
}

////////////////////////////////////////
//	IMPLEMENTATION OF RefinementProjector
bool ProjectionHandler::
refinement_begins_requires_subgrid () const
{
	for(size_t i = 0; i < m_projectors.size(); ++i){
		if(m_projectors[i].valid()
			&& m_projectors[i]->refinement_begins_requires_subgrid())
		{
			return true;
		}
	}
	return false;
}

void ProjectionHandler::
refinement_begins (const ISubGrid* psg)
{
	RefinementProjector::refinement_begins(psg);

	UG_COND_THROW(!m_sh, "Please set a valid SubsetHandler to "
				  "'ProjectionHandler' before using it during refinement.");


	if(!m_defaultProjector->geometry().valid())
		m_defaultProjector->set_geometry(geometry());

//	make sure that the internal vector of projectors is at least as large
//	as there are subsets in the associated subset handler
	projector_required((int)m_sh->num_subsets()-1);

//	this selector will be used if we have to pass sub-grids to individual projectors
	Selector sel;

	for(size_t i = 0;
		i < std::min<size_t>(m_projectors.size(), m_sh->num_subsets() + 1);
		++i)
	{
		if(!m_projectors[i].valid())
			continue;

		RefinementProjector& proj = *m_projectors[i];
		if(!proj.geometry().valid())
			proj.set_geometry(geometry());

	//todo:	move this somewhere where it is only triggered if required.
	//		e.g. to 'set_subset_handler' and 'add_projector'
		proj.set_concerned_elements(make_sp(new IsInSubset(*m_sh, i-1)));

		if(!proj.refinement_begins_requires_subgrid()) {
			proj.refinement_begins(nullptr);
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
						for(Grid::traits<Vertex>::const_iterator _feI = subsetGoc.begin<Vertex>(); _feI != subsetGoc.end<Vertex>(); ++_feI){ Vertex* e = *_feI;{
							if(sg.is_contained(e))
								sel.select(e);
						}};
					}
					if(sg.goc().num_edges() > 0){
						for(Grid::traits<Edge>::const_iterator _feI = subsetGoc.begin<Edge>(); _feI != subsetGoc.end<Edge>(); ++_feI){ Edge* e = *_feI;{
							if(sg.is_contained(e))
								sel.select(e);
						}};
					}
					if(sg.goc().num_faces() > 0){
						for(Grid::traits<Face>::const_iterator _feI = subsetGoc.begin<Face>(); _feI != subsetGoc.end<Face>(); ++_feI){ Face* e = *_feI;{
							if(sg.is_contained(e))
								sel.select(e);
						}};
					}
					if(sg.goc().num_volumes() > 0){
						for(Grid::traits<Volume>::const_iterator _feI = subsetGoc.begin<Volume>(); _feI != subsetGoc.end<Volume>(); ++_feI){ Volume* e = *_feI;{
							if(sg.is_contained(e))
								sel.select(e);
						}};
					}

					selectionPerformed = true;
				}
			}
			
			if(!selectionPerformed){
				if(i == 0 || subsetGoc.num_vertices() > 0){
					for(Grid::traits<Vertex>::const_iterator _feI = sg.goc().begin<Vertex>(); _feI != sg.goc().end<Vertex>(); ++_feI){ Vertex* e = *_feI;{
						if(m_sh->get_subset_index(e) + 1 == (int)i)
							sel.select(e);
					}};
				}
				if(i == 0 || subsetGoc.num_edges() > 0){
					for(Grid::traits<Edge>::const_iterator _feI = sg.goc().begin<Edge>(); _feI != sg.goc().end<Edge>(); ++_feI){ Edge* e = *_feI;{
						if(m_sh->get_subset_index(e) + 1 == (int)i)
							sel.select(e);
					}};
				}
				if(i == 0 || subsetGoc.num_faces() > 0){
					for(Grid::traits<Face>::const_iterator _feI = sg.goc().begin<Face>(); _feI != sg.goc().end<Face>(); ++_feI){ Face* e = *_feI;{
						if(m_sh->get_subset_index(e) + 1 == (int)i)
							sel.select(e);
					}};
				}
				if(i == 0 || subsetGoc.num_volumes() > 0){
					for(Grid::traits<Volume>::const_iterator _feI = sg.goc().begin<Volume>(); _feI != sg.goc().end<Volume>(); ++_feI){ Volume* e = *_feI;{
						if(m_sh->get_subset_index(e) + 1 == (int)i)
							sel.select(e);
					}};
				}
			}
			
			SubGrid<IsSelected> tsg(sel.get_grid_objects(), IsSelected(sel));
			proj.refinement_begins(&tsg);
		}
	}
}

void ProjectionHandler::
refinement_ends()
{
	RefinementProjector::refinement_ends();
	for(size_t i = 0; i < m_projectors.size(); ++i){
		if(m_projectors[i].valid())
			m_projectors[i]->refinement_ends();
	}
}

number ProjectionHandler::
new_vertex (Vertex* vrt, Vertex* parent)
{
	return handle_new_vertex (vrt, parent);
}

number ProjectionHandler::
new_vertex (Vertex* vrt, Edge* parent)
{
	return handle_new_vertex (vrt, parent);
}

number ProjectionHandler::
new_vertex (Vertex* vrt, Face* parent)
{
	return handle_new_vertex (vrt, parent);
}

number ProjectionHandler::
new_vertex (Vertex* vrt, Volume* parent)
{
	return handle_new_vertex (vrt, parent);
}



void ProjectionHandler::
projector_required(int index) {
//	in order to automatically treat subset '-1' correctly, internally we add 1 
//	to each subset-index
	if(index + 1 >= (int)m_projectors.size())
		m_projectors.resize(index + 2);
}

template <class TParent>
number ProjectionHandler::
handle_new_vertex (Vertex* vrt, TParent* parent)
{
	int si = m_sh->get_subset_index(parent);

	UG_ASSERT (si + 1 < (int) m_projectors.size(),
			   "Please make sure to call 'refinement_begins' before calling "
			   "'new_vertex' for the first time for a given subset-assignment.");
	

	number ia = 0;
	if(m_projectors [si + 1].valid())
	{
	    ia = m_projectors [si + 1]->new_vertex (vrt, parent);

	    if(ia < 1){
	        vector3 p = pos (vrt);
	        m_defaultProjector->new_vertex (vrt, parent);
	        vector3 lp = pos (vrt);
	        p *= ia;
	        lp *= (1. - ia);
	        p += lp;
	        set_pos(vrt, p);
	    }
	}
	else m_defaultProjector->new_vertex (vrt, parent);

	return 1;
}

}//	end of namespace
