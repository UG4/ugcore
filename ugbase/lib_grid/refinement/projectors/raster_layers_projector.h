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

#ifndef __H__UG_raster_layers_projector
#define __H__UG_raster_layers_projector

#include "refinement_projector.h"
#include "lib_grid/algorithms/raster_layer_util.h"

namespace ug{

class RasterLayersProjector : public RefinementProjector {
public:
	typedef ANumber					rel_z_attachment_t;
	typedef Grid::VertexAttachmentAccessor<ANumber>	rel_z_attachment_accessor_t;

	RasterLayersProjector ()	{}

	RasterLayersProjector (SPIGeometry3d geometry) :
		RefinementProjector (geometry)
	{add_attachments();}

	RasterLayersProjector (SPIGeometry3d geometry, SPRasterLayers layers) :
		RefinementProjector (geometry)
	{add_attachments();
	 set_layers(layers);}

	void set_geometry(SPIGeometry3d geometry)
	{
		remove_attachments();
		RefinementProjector::set_geometry(geometry);
		add_attachments();
	}

	void set_layers(SPRasterLayers layers)
	{m_layers = layers;}

	rel_z_attachment_t			rel_z_attachment () const			{return m_aRelZ;}
	rel_z_attachment_accessor_t	rel_z_attachment_accessor () const	{return m_aaRelZ;}

	number average_rel_z(Vertex* e) const	{return m_aaRelZ[e];}

	template <class TElem>
	number average_rel_z(TElem* e) const
	{
		typename TElem::ConstVertexArray	vrts = e->vertices();
		const size_t numVrts = e->num_vertices();
		if(numVrts == 0)
			return 0;

		number relZ = 0;
		for(size_t i = 0; i < numVrts; ++i){
			relZ += m_aaRelZ[vrts[i]];
		}

		return relZ / (number)numVrts;
	}

	virtual number new_vertex(Vertex* vrt, Vertex* parent)
	{
		m_aaRelZ[vrt] = m_aaRelZ[parent];
		return 1;
	}

	virtual number new_vertex(Vertex* vrt, Edge* parent)
	{
		return new_vertex_impl(vrt, parent);
	}

	virtual number new_vertex(Vertex* vrt, Face* parent)
	{
		return new_vertex_impl(vrt, parent);
	}

	virtual number new_vertex(Vertex* vrt, Volume* parent)
	{
		return new_vertex_impl(vrt, parent);
	}

private:
	template <class TParent>
	number new_vertex_impl(Vertex* vrt, TParent* parent)
	{
		const number relZ = average_rel_z(parent);
		m_aaRelZ[vrt] = relZ;
		vector3 p = geom().element_center(parent);
		p.z() = m_layers->relative_to_global_height(vector2(p.x(), p.y()), relZ);
		set_pos(vrt, p);
		return 1;
	}

	void add_attachments ()
	{
		IGeometry3d& g = geom();
		g.grid().attach_to_vertices(m_aRelZ);
		m_aaRelZ.access(g.grid(), m_aRelZ);
	}

	void remove_attachments ()
	{
		SPIGeometry3d geom = geometry();
		if(geom.valid() && geom->grid().has_vertex_attachment(m_aRelZ)){
			geom->grid().detach_from_vertices(m_aRelZ);
			m_aaRelZ.invalidate();
		}
	}

	friend class boost::serialization::access;

	template <class Archive>
	void serialize( Archive& ar, const unsigned int version)
	{
		using namespace ug;
		if(ArchiveInfo<Archive>::TYPE == AT_DATA){
			UG_LOG("AT_DATA\n");
		}
		else if (ArchiveInfo<Archive>::TYPE == AT_GUI){
			UG_LOG("AT_GUI\n");
		}

		UG_EMPTY_BASE_CLASS_SERIALIZATION(RasterLayersProjector, RefinementProjector);
	}
	
	SPRasterLayers				m_layers;
	rel_z_attachment_t			m_aRelZ;
	rel_z_attachment_accessor_t	m_aaRelZ;
};

typedef SmartPtr<RasterLayersProjector>	SPRasterLayersProjector;

}//	end of namespace

#endif	//__H__UG_raster_layers_projector
