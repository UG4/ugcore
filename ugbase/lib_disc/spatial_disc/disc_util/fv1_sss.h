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
 * Singular sources and sinks for the FV discretizations:
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__FV_SINGULAR_SOURCES_AND_SINKS__
#define __H__UG__LIB_DISC__SPATIAL_DISC__FV_SINGULAR_SOURCES_AND_SINKS__

#include <vector>

// ug4 headers
#include "lib_grid/algorithms/bounding_box_util.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"

namespace ug {

/// Base class for the point sources and sinks
template <int dim, typename TData>
class FVPointSourceOrSink : public TData
{
private:
	MathVector<dim> point; ///< the point of the source/sink
	
public:

///	class constructor (with the default constructor for data)
	FVPointSourceOrSink
	(
		const MathVector<dim>& _point ///< coordinates of the source/sink
	)
	: point(_point) {}

///	class constructor (with the copy constructor for data)
	FVPointSourceOrSink
	(
		const MathVector<dim>& _point, ///< coordinates of the source/sink
		const TData& _data ///< the data to copy
	)
	: TData(_data), point(_point) {}

///	class constructor (with the default constructor for data)
	FVPointSourceOrSink
	(
		const std::vector<number>& _point ///< coordinates of the source/sink
	)
	{
		if (_point.size () != dim)
			UG_THROW ("FVPointSourceOrSink: Number of coordinates does not match the dimension.");
		for (size_t i = 0; i < dim; i++) point[i] = _point[i];
	}

///	returns the point of the source/sink
	const MathVector<dim>& position () {return point;}
	
///	test whether a source/sink point corresponds to a given corner of the element
    template <typename TElem, typename TAAPos, typename TFVGeom>
    bool corresponds_to
    (
		TElem* elem, ///< the element
		Grid& grid, ///< the grid
		TAAPos& aaPos, ///< position of the vertices
		const TFVGeom& geo, ///< FV geometry (initialized for 'elem')
		size_t co ///< corner to get the contribution for
	);
};

/// Base class for line sources and sinks
template<int dim, typename TData>
class FVLineSourceOrSink : public TData
{
private:
    MathVector<dim> point1; ///< beginning of the line segment
    MathVector<dim> point2; ///< end of the line segment

public:
	
///	class constructor (with the default constructor for data)
    FVLineSourceOrSink
    (
    	const MathVector<dim>& _point1, ///< beginning of the line segment
		const MathVector<dim>& _point2 ///< end of the line segment
	)
	: point1(_point1), point2(_point2) {}

///	class constructor (with the copy constructor for data)
    FVLineSourceOrSink
    (
    	const MathVector<dim>& _point1, ///< beginning of the line segment
		const MathVector<dim>& _point2, ///< end of the line segment
		const TData& _data ///< the data to copy
	)
	: TData(_data), point1(_point1), point2(_point2) {}

///	class constructor (with the default constructor for data)
    FVLineSourceOrSink
    (
    	const std::vector<number>& _point1, ///< beginning of the line segment
		const std::vector<number>& _point2 ///< end of the line segment
	)
	{
		if (_point1.size () != dim || _point2.size () != dim)
			UG_THROW ("FVLineSourceOrSink: Number of coordinates does not match the dimension.");
		for (size_t i = 0; i < dim; i++)
		{
			point1[i] = _point1[i]; point2[i] = _point2[i];
		}
	}

///	returns the beginning of the line segment
	const MathVector<dim>& from_position () {return point1;}
	
///	returns the end of the line segment
	const MathVector<dim>& to_position () {return point2;}
	
///	test whether a source/sink point corresponds to a given corner of the element
    template <typename TElem, typename TAAPos, typename TFVGeom>
    bool corresponds_to
    (
		TElem* elem, ///< [in] the element
		Grid& grid, ///< [in] the grid
		TAAPos& aaPos, ///< [in] position of the vertices
		const TFVGeom& geo, ///< [in] FV geometry (initialized for 'elem')
		size_t co, ///< [in] corner to get the contribution for
		MathVector<dim>& ls, ///< [out] beginning of the subsegment
		MathVector<dim>& le ///< [out] end of the subsegment
	);
};

/// Partial specialization of for 1d: No line sources and sinks in 1d
/// Base class for line sources and sinks
template<typename TData>
class FVLineSourceOrSink<1, TData> : public TData
{
private:
	static constexpr int dim = 1;
    MathVector<dim> point1; ///< beginning of the line segment
    MathVector<dim> point2; ///< end of the line segment

public:
	
///	class constructor (with the default constructor for data)
    FVLineSourceOrSink
    (
    	const MathVector<dim>& _point1, ///< beginning of the line segment
		const MathVector<dim>& _point2 ///< end of the line segment
	)
	: point1(_point1), point2(_point2) {}

///	class constructor (with the copy constructor for data)
    FVLineSourceOrSink
    (
    	const MathVector<dim>& _point1, ///< beginning of the line segment
		const MathVector<dim>& _point2, ///< end of the line segment
		const TData& _data ///< the data to copy
	)
	: TData(_data), point1(_point1), point2(_point2) {}

///	class constructor (with the default constructor for data)
    FVLineSourceOrSink
    (
    	const std::vector<number>& _point1, ///< beginning of the line segment
		const std::vector<number>& _point2 ///< end of the line segment
	)
	{
		UG_THROW ("FVLineSourceOrSink: Line sources and sinks are not supported in 1d.");
	}

///	returns the beginning of the line segment
	const MathVector<dim>& from_position () {return point1;}
	
///	returns the end of the line segment
	const MathVector<dim>& to_position () {return point2;}
	
///	test whether a source/sink point corresponds to a given corner of the element
    template <typename TElem, typename TAAPos, typename TFVGeom>
    bool corresponds_to
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
		UG_THROW ("FVLineSourceOrSink: Line sources and sinks are not supported in 1d.");
		return false;
	}
};

/// Manager class for point and line sources and sinks
template <int dim, typename TPointData, typename TLineData = TPointData>
class FVSingularSourcesAndSinks
{
public:
	using point_sss_type = FVPointSourceOrSink<dim, TPointData>;
	using line_sss_type = FVLineSourceOrSink<dim, TLineData>;

private:
	typename std::vector<SmartPtr<point_sss_type> > ListP;
	typename std::vector<SmartPtr<line_sss_type> > ListL;
	
public:

///	class constructor
	FVSingularSourcesAndSinks () {}

///	add a point source or sink
	void add_point
	(
		SmartPtr<point_sss_type> point_sss ///< the object to add
	)
	{
		ListP.push_back (point_sss);
	}
	
///	add a line source or sink
	void add_line
	(
		SmartPtr<line_sss_type> line_sss ///< the object to add
	)
	{
		ListL.push_back (line_sss);
	}
	
///	returns the number of the point sources and sinks
	size_t num_points () {return ListP.size ();}
	
///	return a point sink with a given index
	point_sss_type * point (size_t i) {return ListP[i].get();}
	
///	calls the init-functions of all point sinks (if there are any in TPointData)
	void init_all_point_sss ()
	{
		for (size_t i = 0; i < ListP.size (); i++)
			((TPointData *) ListP[i].get())->init ();
	}
	
///	class of point source and sink iterator for a given scv
	template <typename TElem, typename TAAPos, typename TFVGeom>
	class point_iterator
	{
		using this_type = point_iterator<TElem, TAAPos, TFVGeom>;
		using master_type = FVSingularSourcesAndSinks<dim, TPointData, TLineData>;
		
		master_type * m_sss;
		TElem * m_elem;
		Grid & m_grid;
		TAAPos & m_aaPos;
		const TFVGeom & m_geo;
		size_t m_co;
		AABox<MathVector<dim> > m_elem_bbox;
		size_t m_index; ///< index of the current source/sink
		
	public:
	
	///	class constructor (sets the iterator to the beginnging of the list)
		point_iterator
		(
			master_type* sss, /// the manager
			TElem* elem, ///< the element to find the sources and sinks for
			Grid& grid, ///< the grid
			TAAPos& aaPos, ///< position attachment for the grid
			const TFVGeom& geo, ///< FV geometry for the element
			size_t co ///< corner to check the sources and sinks for
		)
		:	m_sss (sss), m_elem (elem), m_grid (grid), m_aaPos (aaPos), m_geo (geo), m_co (co)
		{
			m_elem_bbox = CalculateBoundingBox(m_elem, m_aaPos);
			next_sss (0);
		}
		
	///	copy constructor
		point_iterator (this_type& op)
		:	m_elem (op.m_elem), m_grid (op.m_grid), m_aaPos (op.m_aaPos),
			m_geo (op.m_geo), m_co (op.m_co), m_elem_bbox (op.m_elem_bbox),
			m_index (op.m_index)
		{}
		
	///	moves the iterator to the beginning of the list
		this_type& rewind () {m_index = 0; return * this;}
	
	///	access operator
		point_sss_type * operator * () {return m_sss->ListP[m_index].get();}
		
	///	checks whether we are at the end
		bool is_over () {return m_index >= m_sss->num_points ();}
		
	///	moves to the next valid source or sink
		this_type & operator++ ()
		{
			next_sss (m_index + 1);
			return * this;
		}
		
	private:
		void next_sss (size_t index);
	};
	
///	returns the number of the line sources and sinks
	size_t num_lines () {return ListL.size ();}
	
///	return a line sink with a given index
	line_sss_type * line (size_t i) {return ListL[i].get();}
	
///	calls the init-functions of all line sinks (if there are any in TLineData)
	void init_all_line_sss ()
	{
		for (size_t i = 0; i < ListL.size (); i++)
			((TLineData *) ListL[i].get())->init ();
	}
	
///	class of line source and sink iterator for a given scv
	template <typename TElem, typename TAAPos, typename TFVGeom>
	class line_iterator
	{
		using this_type = line_iterator<TElem, TAAPos, TFVGeom>;
		using master_type = FVSingularSourcesAndSinks<dim, TPointData, TLineData>;
		
		master_type * m_sss;
		TElem * m_elem;
		Grid & m_grid;
		TAAPos & m_aaPos;
		const TFVGeom & m_geo;
		size_t m_co;
		AABox<MathVector<dim> > m_elem_bbox;
		
		size_t m_index; ///< index of the current source/sink
		MathVector<dim> m_ls; ///< the 1st of the intersection points with the bnd of the scv
		MathVector<dim> m_le; ///< the 2nd of the intersection points with the bnd of the scv
		
	public:
	
	///	class constructor (sets the iterator to the beginnging of the list)
		line_iterator
		(
			master_type* sss, /// the manager
			TElem* elem, ///< the element to find the sources and sinks for
			Grid& grid, ///< the grid
			TAAPos& aaPos, ///< position attachment for the grid
			const TFVGeom& geo, ///< FV geometry for the element
			size_t co ///< corner to check the sources and sinks for
		)
		:	m_sss (sss), m_elem (elem), m_grid (grid), m_aaPos (aaPos), m_geo (geo), m_co (co)
		{
			m_elem_bbox = CalculateBoundingBox(m_elem, m_aaPos);
			next_sss (0, m_ls, m_le);
		}
		
	///	copy constructor
		line_iterator (this_type& op)
		:	m_elem (op.m_elem), m_grid (op.m_grid), m_aaPos (op.m_aaPos),
			m_geo (op.m_geo), m_co (op.m_co), m_elem_bbox (op.m_elem_bbox),
			m_index (op.m_index)
		{}
		
	///	moves the iterator to the beginning of the list
		this_type& rewind () {m_index = 0; return * this;}
	
	///	access operator
		line_sss_type * operator * () {return m_sss->ListL[m_index].get();}
		
	///	the 1st of the intersection points with the bnd of the scv
		const MathVector<dim>& seg_start () {return m_ls;}
		
	///	the 2nd of the intersection points with the bnd of the scv
		const MathVector<dim>& seg_end () {return m_le;}
		
	///	checks whether we are at the end
		bool is_over () {return m_index >= m_sss->num_lines ();}
		
	///	moves to the next valid source or sink
		this_type & operator++ ()
		{
			next_sss (m_index + 1, m_ls, m_le);
			return * this;
		}
		
	private:
		void next_sss (size_t index, MathVector<dim>& ls, MathVector<dim>& le);
	};
};

} // namespace ug

#include "fv1_sss_impl.h"

#endif
