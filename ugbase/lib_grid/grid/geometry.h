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

#ifndef __H__UG_geometry
#define __H__UG_geometry

#include "common/util/smart_pointer.h"
#include "common/math/ugmath_types.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "lib_grid/grid/grid.h"

namespace ug{

///	provides a grid and access to the coordinates of the vertices
/**	The class serves as a base class for wrappers around grid coordinates. It
 * allows for the implementation of fixed-dim classes which can operate on
 * coordinates arbitrary dimension.
 *
 * Vertex coordinates of the fixed dimension 'dim' can be accessed through the
 * method 'pos' and can be set through the method 'set_pos'.
 *
 * For code that has to operate on the original position attachment, the following
 * methods are provided: The dimension of the original 'attached'
 * coordinates can be accessed through 'position_attachment_dim'. The attachement
 * for those original coordinates can be accessed through 'position_attachment'.
 *
 * If 'attachmentDim == position_attachment_dim()', then the return value of
 * 'position_attachment' must always be castable to
 * 'Attachment<MathVector<attachmentDim> >'.
 *
 * \sa Geometry
 */
template <int dim>
class IGeometry {
public:
	using vector_t = MathVector<dim>;

	explicit IGeometry (Grid& g) :
		m_grid(g)
	{}

	virtual ~IGeometry () = default;

	virtual vector_t pos (Vertex* vrt) const = 0;
	virtual void set_pos (Vertex* vrt, const vector_t& p) = 0;

	[[nodiscard]] Grid& grid() const {return m_grid;}

	[[nodiscard]] virtual int position_attachment_dim() const = 0;
	[[nodiscard]] virtual const IAttachment& position_attachment() const = 0;

	virtual vector_t element_center(Vertex* e) const	= 0;
	virtual vector_t element_center(Edge* e) const		= 0;
	virtual vector_t element_center(Face* e) const		= 0;
	virtual vector_t element_center(Volume* e) const	= 0;

private:
	Grid&	m_grid;
};


using IGeometry1d = IGeometry<1>;
using IGeometry2d = IGeometry<2>;
using IGeometry3d = IGeometry<3>;

using SPIGeometry1d = SmartPtr<IGeometry1d>;
using SPIGeometry2d = SmartPtr<IGeometry2d>;
using SPIGeometry3d = SmartPtr<IGeometry3d>;


///	provides a grid and access to the coordinates of the vertices
/**	The dimension of the output coordinates ('dim') may differ from the
 * dimension of the coordinates which are actually attached to the grid ('attachedDim').
 */
template <int dim, int attachmentDim>
class Geometry : public IGeometry<dim> {
public:
	using base_t = IGeometry<dim>;
	using vector_t = typename base_t::vector_t;
	using attached_vector_t = MathVector<attachmentDim>;
	using position_attachment_t = Attachment<attached_vector_t>;

	Geometry (Grid& g, position_attachment_t a) :
		base_t (g),
		m_aPos(a)
	{
		UG_COND_THROW (!g.has_vertex_attachment(a),
					   "the specified attachment is not attached to the "
					   "vertices of the specified grid.");
		m_aaPos.access(g, a);
	}

	~Geometry () override = default;

	vector_t pos (Vertex* vrt) const override {
		return vector_t::from(m_aaPos[vrt]);
	}

	void set_pos (Vertex* vrt, const vector_t& p) override {
		m_aaPos[vrt] = attached_vector_t::from(p);
	}

	[[nodiscard]] int position_attachment_dim () const override {return attachmentDim;}
	[[nodiscard]] const IAttachment& position_attachment () const override {return m_aPos;}


	vector_t element_center (Vertex* e) const override {
		return vector_t::from(m_aaPos[e]);
	}

	vector_t element_center (Edge* e) const override {
		return vector_t::from(CalculateCenter(e, m_aaPos));
	}

	vector_t element_center (Face* e) const override {
		return vector_t::from(CalculateCenter(e, m_aaPos));
	}

	vector_t element_center (Volume* e) const override {
		return vector_t::from(CalculateCenter(e, m_aaPos));
	}


private:
	Grid::VertexAttachmentAccessor<position_attachment_t> m_aaPos;
	position_attachment_t m_aPos;
};



///	Utility method to construct an IGeometry3d for a given grid and position attachment
template <class TAPos>
SPIGeometry3d
MakeGeometry3d (Grid& grid, TAPos aPos)
{
	return make_sp(new Geometry<3, TAPos::ValueType::Size>(grid, aPos));
}
}//	end of namespace

#endif