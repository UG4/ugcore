/*
 * Copyright (c) 2018:  G-CSC, Goethe University Frankfurt
 * Author: Stephan Grein
 *
 * This file is part of UG4.
 *
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3):
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

#ifndef __H__UG_cylinder_soma_projector
#define __H__UG_cylinder_soma_projector

#include "common/math/misc/math_util.h"
#include "refinement_projector.h"
#include "lib_grid/tools/copy_attachment_handler.h"
#include "lib_grid/global_attachments.h"

namespace ug {
///	Projects new vertices onto a sphere during refinement.
/** For projection during refinement the radius property is ignored. Instead
 * the distance to the center of a newly inserted vertex is calculated
 * as the average distance of the vertices of the parent element to the center.
 * The radius property thus defaults to -1.
 *
 * You may still specify a radius. This radius can be used for auto-fitting of
 * the center and for reprojecting a set of vertices onto the sphere.
 *
 * Only vertices which are at the soma (axial parameter: -1) are projected
 * on the surface. This works well for good natured connecting neurites.
 */
class SomaProjector : public RefinementProjector {
private:
	Attachment<NeuriteProjector::SurfaceParams> m_aSurfParams;
    Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> > m_aaSurfParams;
    CopyAttachmentHandler<Vertex, Attachment<NeuriteProjector::SurfaceParams> > m_cah;
protected:
    virtual void set_geometry(SPIGeometry3d geometry)
    {
        // call base class method
        RefinementProjector::set_geometry(geometry);
		attach_surf_params();
    }
public:
		SomaProjector () :
			m_center (0, 0, 0),
			m_axis (0, 0, 1),
			m_soma (0, 0, 0),
			m_somaRad(0)
		{
			using NPSurfParam = Attachment<NeuriteProjector::SurfaceParams>;
				if (!GlobalAttachments::is_declared("npSurfParams"))
					GlobalAttachments::declare_attachment<NPSurfParam>("npSurfParams", true);

		}

		SomaProjector (const vector3& center,
					  const vector3& axis,
					  const vector3& soma,
					  const number& rad) :

		m_center (center),
		m_axis (axis),
		m_soma (soma),
		m_somaRad(rad)
	{
		using NPSurfParam = Attachment<NeuriteProjector::SurfaceParams>;
			if (!GlobalAttachments::is_declared("npSurfParams"))
				GlobalAttachments::declare_attachment<NPSurfParam>("npSurfParams", true);
	}

public:
/**	\sa ug::RefinementProjector::RefinementProjector*/
	SomaProjector (SPIGeometry3d geometry,
				  const vector3& center,
				  const vector3& axis,
				  const vector3& soma,
				  const number& rad) :
		RefinementProjector (geometry),
		m_center (center),
		m_axis (axis),
		m_soma (soma),
		m_somaRad(rad)
	{
		attach_surf_params();
	}


public:
	virtual ~SomaProjector ()					{}

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
/// check if global attachment was declared
	void check_attachment() {
		using NPSurfParam = Attachment<NeuriteProjector::SurfaceParams>;
		if (!GlobalAttachments::is_declared("npSurfParams"))
			GlobalAttachments::declare_attachment<NPSurfParam>("npSurfParams", true);
	}


/// helper method to attach surface parameters and do a consistency check
	void attach_surf_params() {
		check_attachment();

		Grid& grid = this->geometry()->grid();

		// make sure attachment was declared
		UG_COND_THROW(!GlobalAttachments::is_declared("npSurfParams"), "GlobalAttachment 'npSurfParams' not declared.");
		m_aSurfParams = GlobalAttachments::attachment<Attachment<NeuriteProjector::SurfaceParams> >("npSurfParams");

		// make sure surfaceParams attachment is attached
	    UG_COND_THROW(!grid.has_vertex_attachment(m_aSurfParams),
        "No surface parameter attachment for neurite projector attached to grid.");
		 m_aaSurfParams.access(grid, m_aSurfParams);

		// handle attachment values also on higher grid levels (if required)
		MultiGrid* mg = dynamic_cast<MultiGrid*>(&grid);
		if (mg) {
			SmartPtr<MultiGrid> spMG(mg);

			// never destroy the grid from here - we did not create it
		    ++(*spMG.refcount_ptr());

		    // attach copy attachment handler to propagate attachment to higher levels
		    m_cah.set_attachment(m_aSurfParams);
		    m_cah.set_grid(spMG);
	   }
	}


	template <typename TElem>
	number perform_projection(Vertex* vrt, TElem* parent)
	{
	//	calculate the new position by linear interpolation and project that point
	//	onto the cylinder.
		typename TElem::ConstVertexArray vrts = parent->vertices();
		size_t numVrts = parent->num_vertices();

		if(numVrts == 0){
			set_pos(vrt, vector3(0, 0, 0));
			return 1;
		}

		number avDist = 0;
		vector3 parentCenter (0, 0, 0);

		for(size_t i = 0; i < numVrts; ++i){
			vector3 p = pos(vrts[i]);
			avDist += DistancePointToRay(p, m_center, m_axis);
			parentCenter += p;
		}

		avDist /= (number)numVrts;
		VecScale(parentCenter, parentCenter, 1. / (number)numVrts);

		vector3 proj, v;
		ProjectPointToRay(proj, parentCenter, m_center, m_axis);
		VecSubtract(v, parentCenter, proj);
		number len = VecLength(v);
		if(len > SMALL * avDist){	// if avDist is very small, len may be small, too
			VecScale(v, v, avDist/len);
			proj += v;
			set_pos(vrt, proj);
		}
		else
			set_pos(vrt, parentCenter);

		if (vertex_at_soma_surf(vrt, parent))
			project_to_soma_surface(vrt, parent);

		return 1;
	}

	template <typename TElem>
	bool vertex_at_soma_surf(Vertex* vrt, TElem* parent) {
		size_t numSomaVerts = 0;
	    size_t nVrt = parent->num_vertices();
	    for (size_t i = 0; i < nVrt; ++i) {
	    	 if(m_aaSurfParams[parent->vertex(i)].axial == -1.0) {
	    		 numSomaVerts++;
	    	 }
	    }
	    return numSomaVerts >= 2;
	}

	template <typename TElem>
	void project_to_soma_surface(Vertex* vrt, TElem* parent) {
		/// old vertex position
		vector3 v0 = this->pos(vrt);
		/// 1. P = (x',y',z') = (x - x0, y - y0, z - z0)
		vector3 P;
		VecSubtract(P, v0, m_soma);
		/// 2. |P| = sqrt(x'^2 + y'^2 + z'^2)
		number Plength = VecLength(P);
		/// 3. Q = (radius/|P|)*P
		vector3 Q = P;
		VecScale(Q, P, m_somaRad / Plength);
		/// 4. R = Q + (x0,y0,z0)
		vector3 R;
		VecAdd(R, Q, m_soma);
		/// 5. Set new position
		set_pos(vrt, R);
		/// 6. set the surface parameters for the new vertex
		m_aaSurfParams[vrt].axial = -1;
	}


	friend class boost::serialization::access;

	template <typename Archive>
	void serialize( Archive& ar, const unsigned int version)
	{
		ar & make_nvp("center", m_center);
		ar & make_nvp("axis", m_axis);
		ar & make_nvp("soma", m_soma);;
		ar & make_nvp("somaRad", m_somaRad);;
		UG_EMPTY_BASE_CLASS_SERIALIZATION(SomaProjector, RefinementProjector);
	}

	vector3	m_center;
	vector3	m_axis;
	vector3 m_soma;
	number m_somaRad;
};

}//	end of namespace

#endif