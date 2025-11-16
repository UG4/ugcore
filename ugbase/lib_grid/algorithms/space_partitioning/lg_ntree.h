/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__lg_ntree__
#define __H__UG__lg_ntree__

#include "common/space_partitioning/ntree.h"
#include "common/space_partitioning/ntree_traverser.h"
#include "common/math/misc/shapes.h"
#include "lib_grid/algorithms/ray_element_intersection_util.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"

namespace ug{


template <int world_dim>
class NTreeGridData
{
	public:
		using position_t = MathVector<world_dim >;
		using position_attachment_t = Attachment<position_t>;
		using position_accessor_t = Grid::VertexAttachmentAccessor<position_attachment_t>;

		NTreeGridData() : m_pGrid(nullptr)	{}

		NTreeGridData(Grid& grid, position_attachment_t aPos)
		{
			m_pGrid = &grid;
			if(!grid.has_vertex_attachment(aPos))
				grid.attach_to_vertices(aPos);
			m_aaPos.access(grid, aPos);
		}

		const position_t& position(Vertex* v) const
		{
			UG_ASSERT(m_aaPos.valid(),
					  "Make sure to pass an instance of NTreeGridData to lg_ntree::set_common_data");
			return m_aaPos[v];
		}

		position_accessor_t position_accessor() const	{return m_aaPos;}

		Grid* grid_ptr() const	{return m_pGrid;}

	private:
		position_accessor_t	m_aaPos;
		Grid*				m_pGrid;

};


template <int tree_dim, int world_dim, class elem_t_, class common_data_t_>
struct lg_ntree_traits_base
{
	using real_t = number;
	using vector_t = MathVector<world_dim>;
	using box_t = AABox<vector_t>;
	using common_data_t = common_data_t_;
	using elem_t = elem_t_;

	static void calculate_center(vector_t& centerOut, const elem_t& e,
								 const common_data_t& commonData)
	{
		Grid::vertex_traits::secure_container vrts;
		commonData.grid_ptr()->associated_elements(vrts, e);

		assert(vrts.size() > 0);
		centerOut = commonData.position(vrts[0]);
		for(size_t i = 1; i < vrts.size(); ++i)
			VecAdd(centerOut, centerOut, commonData.position(vrts[i]));
		VecScale(centerOut, centerOut, 1. / (real_t)vrts.size());
	}

	static void calculate_bounding_box(box_t& boxOut, const elem_t& e,
									   const common_data_t& commonData)
	{
		Grid::vertex_traits::secure_container vrts;
		commonData.grid_ptr()->associated_elements(vrts, e);

		assert(vrts.size() > 0);
		boxOut.min = boxOut.max = commonData.position(vrts[0]);
		for(size_t i = 1; i < vrts.size(); ++i)
			boxOut = box_t(boxOut, commonData.position(vrts[i]));
	}

	static void grow_box(box_t& boxOut, const box_t& box,
						 const vector_t& offset)
	{
		VecSubtract(boxOut.min, box.min, offset);
		VecAdd(boxOut.max, box.max, offset);
	}

	static vector_t box_diagonal(const box_t& box)
	{
		vector_t d;
		VecSubtract(d, box.max, box.min);
		return d;
	}

	static bool box_contains_point(const box_t& box, const vector_t& point)
	{
		return box.contains_point(point);
	}

///	returns true if the given boxes intersect
	static bool box_box_intersection(const box_t& box1, const box_t& box2)
	{
		for(int i = 0; i < world_dim; ++i){
			if(box1.min[i] > box2.max[i] || box1.max[i] < box2.min[i])
				return false;
		}
		return true;
	}

	static bool ray_box_intersection(const vector_t& from,
									 const vector_t& dir,
									 const box_t& box)
	{
		return RayBoxIntersection(from, dir, box.min, box.max);
	}

///	returns the smallest box that contains both box1 and box2
	static void merge_boxes(box_t& boxOut, const box_t& box1, const box_t& box2)
	{
		boxOut = box_t(box1, box2);
	}

	static bool contains_point(const elem_t& e, const vector_t& point,
					 	 	   const common_data_t& commonData)
	{
		return ContainsPoint(e, point, commonData.position_accessor());

	//	todo: think about some fast-checks like the following (but faster!!!)
	//	first check whether the bounding box contains the point, then perform
	//	the exact check (slow e.g. for tetrahedrons...)
//		typename common_data_t::position_accessor_t aaPos = commonData.position_accessor();
//		box_t box(aaPos[e->vertex(0)], aaPos[e->vertex(1)]);
//		for(size_t i = 2; i < e->num_vertices(); ++i)
//			box = box_t(box, aaPos[e->vertex(i)]);
//
//		if(box.contains_point(point))
//			return ContainsPoint(e, point, aaPos);
//		return false;
	}


	static bool intersects_ray( Vertex* e,
								const vector_t& rayFrom,
								const vector_t& rayDir,
								const common_data_t& cd,
								number& s0out,
								number& s1out)
	{
		UG_THROW("intersects_ray not yet implemented for vertices");
		return false;
	}

	static bool intersects_ray( Edge* e,
								const vector_t& rayFrom,
								const vector_t& rayDir,
								const common_data_t& cd,
								number& s0out,
								number& s1out,
								number sml = SMALL)
	{
		UG_COND_THROW(!cd.grid_ptr(), "No grid assigned to ntree::common_data.");
		return RayElementIntersection(s0out, s1out, rayFrom, rayDir,
									  e, *cd.grid_ptr(), cd.position_accessor(),
									  sml);
	}

	static bool intersects_ray( Face* e,
								const vector_t& rayFrom,
								const vector_t& rayDir,
								const common_data_t& cd,
								number& s0out,
								number& s1out,
								number sml = SMALL)
	{
		UG_COND_THROW(!cd.grid_ptr(), "No grid assigned to ntree::common_data.");
		return RayElementIntersection(s0out, s1out, rayFrom, rayDir,
									  e, *cd.grid_ptr(), cd.position_accessor(),
									  sml);
	}

	static bool intersects_ray( Volume* e,
								const vector_t& rayFrom,
								const vector_t& rayDir,
								const common_data_t& cd,
								number& s0out,
								number& s1out,
								number sml = SMALL)
	{
		UG_COND_THROW(!cd.grid_ptr(), "No grid assigned to ntree::common_data.");
		return RayElementIntersection(s0out, s1out, rayFrom, rayDir,
									  e, *cd.grid_ptr(), cd.position_accessor(),
									  sml);
	}
};


template <class elem_t>
struct ntree_traits<1, 1, elem_t, NTreeGridData<1> > :
	public lg_ntree_traits_base<1, 1, elem_t, NTreeGridData<1> >
{
	using vector_t = MathVector<1>;
	using box_t = AABox<vector_t>;

	static void split_box(box_t* boxesOut, const box_t& box, const vector_t& splitPoint)
	{
		boxesOut[0] = box_t(box.min, splitPoint);
		boxesOut[1] = box_t(splitPoint, box.max);
	}
};


template <class elem_t>
struct ntree_traits<1, 2, elem_t, NTreeGridData<2> > :
	public lg_ntree_traits_base<1, 2, elem_t, NTreeGridData<2> >
{
	using vector_t = MathVector<2>;
	using box_t = AABox<vector_t>;

	static void split_box(box_t* boxesOut, const box_t& box, const vector_t& splitPoint)
	{
		vector_t splitPointYMin = splitPoint;
		splitPointYMin.y() = box.min.y();
		vector_t splitPointYMax = splitPoint;
		splitPointYMax.y() = box.max.y();
		boxesOut[0] = box_t(box.min, splitPointYMax);
		boxesOut[1] = box_t(splitPointYMin, box.max);
	}
};


template <class elem_t>
struct ntree_traits<2, 2, elem_t, NTreeGridData<2> > :
	public lg_ntree_traits_base<2, 2, elem_t, NTreeGridData<2> >
{
	using vector_t = MathVector<2>;
	using box_t = AABox<vector_t>;

	static void split_box(box_t* boxesOut, const box_t& box, const vector_t& splitPoint)
	{
		boxesOut[0] = box_t(box.min, splitPoint);
		boxesOut[1] = box_t(vector_t(splitPoint.x(), box.min.y()),
							vector_t(box.max.x(), splitPoint.y()));
		boxesOut[2] = box_t(vector_t(box.min.x(), splitPoint.y()),
							vector_t(splitPoint.x(), box.max.y()));
		boxesOut[3] = box_t(splitPoint, box.max);
	}
};


template <class elem_t>
struct ntree_traits<2, 3, elem_t, NTreeGridData<3> > :
	public lg_ntree_traits_base<2, 3, elem_t, NTreeGridData<3> >
{
	using vector_t = MathVector<3>;
	using box_t = AABox<vector_t>;

	static void split_box(box_t* boxesOut, const box_t& box, const vector_t& splitPoint)
	{
		vector_t splitPointZMin = splitPoint;
		splitPointZMin.z() = box.min.z();
		vector_t splitPointZMax = splitPoint;
		splitPointZMax.z() = box.max.z();

		boxesOut[0] = box_t(box.min, splitPointZMax);
		boxesOut[1] = box_t(vector_t(splitPoint.x(), box.min.y(), box.min.z()),
							vector_t(box.max.x(), splitPoint.y(), box.max.z()));
		boxesOut[2] = box_t(vector_t(box.min.x(), splitPoint.y(), box.min.z()),
							vector_t(splitPoint.x(), box.max.y(), box.max.z()));
		boxesOut[3] = box_t(splitPointZMin, box.max);
	}
};

template <class elem_t>
struct ntree_traits<3, 3, elem_t, NTreeGridData<3> > :
	public lg_ntree_traits_base<3, 3, elem_t, NTreeGridData<3> >
{
	using vector_t = MathVector<3>;
	using box_t = AABox<vector_t>;

	static void split_box(box_t* boxesOut, const box_t& box, const vector_t& splitPoint)
	{
		boxesOut[0] = box_t(box.min, splitPoint);
		boxesOut[1] = box_t(vector_t(splitPoint.x(), box.min.y(), box.min.z()),
							vector_t(box.max.x(), splitPoint.y(), splitPoint.z()));
		boxesOut[2] = box_t(vector_t(box.min.x(), splitPoint.y(), box.min.z()),
							vector_t(splitPoint.x(), box.max.y(), splitPoint.z()));
		boxesOut[3] = box_t(vector_t(splitPoint.x(), splitPoint.y(), box.min.z()),
							vector_t(box.max.x(), box.max.y(), splitPoint.z()));

		boxesOut[4] = box_t(vector_t(box.min.x(), box.min.y(), splitPoint.z()),
							vector_t(splitPoint.x(), splitPoint.y(), box.max.z()));
		boxesOut[5] = box_t(vector_t(splitPoint.x(), box.min.y(), splitPoint.z()),
							vector_t(box.max.x(), splitPoint.y(), box.max.z()));
		boxesOut[6] = box_t(vector_t(box.min.x(), splitPoint.y(), splitPoint.z()),
							vector_t(splitPoint.x(), box.max.y(), box.max.z()));
		boxesOut[7] = box_t(splitPoint, box.max);
	}
};



template <int tree_dim, int world_dim, class grid_elem_t>
class lg_ntree : public ntree<tree_dim, world_dim, grid_elem_t*, NTreeGridData<world_dim> >
{
	public:
		using base_t = ntree<tree_dim, world_dim, grid_elem_t*, NTreeGridData<world_dim> >;
		using position_attachment_t = typename NTreeGridData<world_dim>::position_attachment_t;

		lg_ntree()
		{}

		lg_ntree(Grid& grid, position_attachment_t aPos) :
			m_gridData(grid, aPos)
		{}

		void set_grid(Grid& grid, position_attachment_t aPos)
		{
			m_gridData = NTreeGridData<world_dim>(grid, aPos);
		}

		template <class TIterator>
		void create_tree(TIterator elemsBegin, TIterator elemsEnd)
		{
			base_t::set_common_data(m_gridData);

			base_t::clear();

			while(elemsBegin != elemsEnd){
				base_t::add_element(*elemsBegin);
				++elemsBegin;
			}

			base_t::rebalance();
		}

	private:
		NTreeGridData<world_dim>	m_gridData;
};


}// end of namespace

#endif
