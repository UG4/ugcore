/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
 * Author: Michael Lampe
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

/*
 * Singular sources and sinks for the FV discretizations: Implementation of the functions
 */

// ug4 headers
#include "lib_grid/algorithms/ray_element_intersection_util.h"
#include "lib_disc/reference_element/reference_element_util.h"

namespace ug{

/**
 * Test whether a source/sink point corresponds to a given corner of the element.
 * Note that \c TElem must specify the dimensionality of the element. Thus,
 * \c GridObject is not a proper type for \c TElem. It should be at least
 * \c Edge, \c Face, \c Volume or derived classes.
 */
template <int dim, typename TData>
template <typename TElem, typename TAAPos, typename TFVGeom>
bool FVPointSourceOrSink<dim, TData>::corresponds_to
(
	TElem* elem, ///< the element
	Grid& grid, ///< the grid
	TAAPos& aaPos, ///< position of the vertices
	const TFVGeom& geo, ///< FV geometry (initialized for 'elem')
	size_t co ///< corner to get the contribution for
)
{
// restrict point to element
	if (!ContainsPoint(elem, point, aaPos))
		return false;

// restrict point to scv
	for (size_t i = 0; i < geo.num_scvf(); i++)
	{
		const typename TFVGeom::SCVF& scvf = geo.scvf(i);
		MathVector<dim> d;
		VecSubtract(d, point, scvf.global_ip());
		if (scvf.from() == co)
		{
			if (VecDot(d, scvf.normal()) >= 0.0)
				return false;
		}
		else if (scvf.to() == co)
		{
			if (VecDot(d, scvf.normal()) < 0.0)
				return false;
		}
	}
	return true;
};

/**
 * Test whether a line source/sink corresponds to a given corner of an element.
 * If yes, the function returns true and sets ls and le to the beginning and the
 * end of the subsegment of the source/sink.
 *
 * Note that \c TElem must specify the dimensionality of the element. Thus,
 * \c GridObject is not a proper type for \c TElem. It should be at least
 * \c Edge, \c Face, \c Volume or derived classes.
 *
 * \todo This function works only for there is only one scv per corner. This is not true for pyramids!
 */
template <int dim, typename TData>
template <typename TElem, typename TAAPos, typename TFVGeom>
bool FVLineSourceOrSink<dim, TData>::corresponds_to
(
	TElem* elem, ///< [in] the element
	Grid& grid, ///< [in] the grid
	TAAPos& aaPos, ///< [in] position of the vertices
	const TFVGeom& geo, ///< [in] FV geometry (initialized for 'elem')
	size_t co, ///< [in] corner to get the contribution for
	MathVector<dim>& ls, ///< [out] beginning of the subsegment
	MathVector<dim>& le ///< [out] end of the subsegment
)
{
//	get the dimensionality of the element
	int elem_dim = ReferenceElementDimension(elem->reference_object_id());

//	restrict line segment to element (so that ls and le are the intersection points with its sides)
	ls = point1;
	le = point2;
	MathVector<dim> dir;
	number lambda_min, lambda_max;
	VecSubtract(dir, le, ls);
	if (!RayElementIntersection(lambda_min, lambda_max, ls, dir, elem, grid, aaPos))
		return false;
	
	if (elem_dim == dim) // full-dimensional element, two different points
	{
	// extend the ends of the segment to the sides of the element
		if (ContainsPoint(elem, ls, aaPos))
			lambda_min = 0.0;
		if (ContainsPoint(elem, le, aaPos))
			lambda_max = 1.0;
	
	// skip the segment if its line cuts the element before or after the segment
		if (lambda_min < 0.0 || lambda_max > 1.0)
			return false;
	
	// compute the corrected ends
		VecScaleAdd(le, 1.0, ls, lambda_max, dir);
		VecScaleAdd(ls, 1.0, ls, lambda_min, dir);

	//	restrict line segment to scv (so that ls and le are the intersection points with the sides of the scv)
		MathVector<dim> p;
		for (size_t i = 0; i < geo.num_scvf(); i++)
		{
			const typename TFVGeom::SCVF& scvf = geo.scvf(i);
			number d_min, d_max, lambda;
			if (scvf.from() == co)
			{
				VecSubtract(p, ls, scvf.global_ip());
				d_min = VecDot(p, scvf.normal());
				VecSubtract(p, le, scvf.global_ip());
				d_max = VecDot(p, scvf.normal());
				if (d_min*d_max < 0.0)
				{
					lambda = fabs(d_min)/(fabs(d_min)+fabs(d_max));
					VecScaleAdd(p, 1.0-lambda, ls, lambda, le);
					if (d_max > 0.0)
						le = p;
					else
						ls = p;
				}
				else if (d_min > 0.0) // if so than d_max >= 0, too
					return false;
			}
			else if (scvf.to() == co)
			{
				VecSubtract(p, ls, scvf.global_ip());
				d_min = VecDot(p, scvf.normal());
				VecSubtract(p, le, scvf.global_ip());
				d_max = VecDot(p, scvf.normal());
				if (d_min*d_max < 0.0)
				{
					lambda = fabs(d_min)/(fabs(d_min)+fabs(d_max));
					VecScaleAdd(p, 1.0-lambda, ls, lambda, le);
					if (d_max < 0.0)
						le = p;
					else
						ls = p;
				}
				else if (d_min < 0.0) // if so than d_max <= 0, too
					return false;
			}
		}
	}
	else if (elem_dim + 1 == dim) // low-dimensional element, we accept only one point
	{
	// skip the segment if its line cuts the element before or after the segment
		if (lambda_min < 0.0 || lambda_max > 1.0)
			return false;
		
	//	exclude the case the the line in the plane
		if (std::fabs (lambda_max - lambda_min) > SMALL)
			UG_THROW ("FVLineSourceOrSink: Line source or sink lyes in a low-dimensional element!");

	//	check if the point is in the scv
		VecScaleAdd(ls, 1.0, point1, lambda_min, dir);
		for (size_t i = 0; i < geo.num_scvf(); i++)
		{
			const typename TFVGeom::SCVF& scvf = geo.scvf(i);
			MathVector<dim> d;
			VecSubtract(d, ls, scvf.global_ip());
			if (scvf.from() == co)
			{
				if (VecDot(d, scvf.normal()) > 0.0)
					return false;
			}
			else if (scvf.to() == co)
			{
				if (VecDot(d, scvf.normal()) < 0.0)
					return false;
			}
		}
		le = ls;
	}
	else
		UG_THROW ("FVLineSourceOrSink: Dimension mismatch with a grid element!");
	
	return true;
};

/**
 * Finds the proper point source or sink starting with a given index
 */
template <int dim, typename TPointData, typename TLineData>
template <typename TElem, typename TAAPos, typename TFVGeom>
void FVSingularSourcesAndSinks<dim, TPointData, TLineData>::point_iterator<TElem, TAAPos, TFVGeom>::next_sss
(
	size_t index ///< index to start the search with
)
{
	for (m_index = index; m_index < m_sss->num_points (); m_index++)
	{
		point_sss_type * pnt_sss = m_sss->ListP[m_index].get();
		if (! m_elem_bbox.contains_point (pnt_sss->position ())) // a quick test
			continue;
		if (pnt_sss->corresponds_to (m_elem, m_grid, m_aaPos, m_geo, m_co))
			break;
	}
};

/**
 * Finds the proper point source or sink starting with a given index
 */
template <int dim, typename TPointData, typename TLineData>
template <typename TElem, typename TAAPos, typename TFVGeom>
void FVSingularSourcesAndSinks<dim, TPointData, TLineData>::line_iterator<TElem, TAAPos, TFVGeom>::next_sss
(
	size_t index, ///< index to start the search with
	MathVector<dim>& ls, ///< the 1st of the intersection points with the bnd of the scv
	MathVector<dim>& le ///< the 2nd of the intersection points with the bnd of the scv
)
{
	for (m_index = index; m_index < m_sss->num_lines (); m_index++)
	{
		line_sss_type * line_sss = m_sss->ListL[m_index].get();
		if (! m_elem_bbox.overlaps_line (line_sss->from_position (), line_sss->to_position ())) // a quick test
			continue;
		if (line_sss->corresponds_to (m_elem, m_grid, m_aaPos, m_geo, m_co, ls, le))
			break;
	}
};

} // end of namespace ug

