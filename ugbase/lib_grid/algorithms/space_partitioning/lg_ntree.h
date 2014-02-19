// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Sep 4, 2013 (d,m,y)

#ifndef __H__UG__lg_ntree__
#define __H__UG__lg_ntree__

#include "common/space_partitioning/ntree.h"
#include "common/space_partitioning/ntree_traverser.h"
#include "common/math/misc/shapes.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"

namespace ug{


template <int world_dim>
class NTreeGridData
{
	public:
		typedef MathVector<world_dim >									position_t;
		typedef Attachment<position_t>									position_attachment_t;
		typedef Grid::VertexAttachmentAccessor<position_attachment_t> 	position_accessor_t;

		NTreeGridData() : m_pGrid(NULL)	{}

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
	typedef number					real_t;
	typedef MathVector<world_dim>	vector_t;
	typedef AABox<vector_t>			box_t;
	typedef common_data_t_			common_data_t;
	typedef elem_t_					elem_t;

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

	static bool box_contains_point(const box_t& box, const vector_t& point)
	{
		return box.contains_point(point);
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
};


template <class elem_t>
struct ntree_traits<1, 1, elem_t, NTreeGridData<1> > :
	public lg_ntree_traits_base<1, 1, elem_t, NTreeGridData<1> >
{
	typedef MathVector<1>			vector_t;
	typedef AABox<vector_t>			box_t;

	static void split_box(box_t* boxesOut, const box_t& box, const vector_t& splitPoint)
	{
		boxesOut[0] = box_t(box.min, splitPoint);
		boxesOut[1] = box_t(splitPoint, box.max);
	}
};


template <class elem_t>
struct ntree_traits<2, 2, elem_t, NTreeGridData<2> > :
	public lg_ntree_traits_base<2, 2, elem_t, NTreeGridData<2> >
{
	typedef MathVector<2>			vector_t;
	typedef AABox<vector_t>			box_t;

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
	typedef MathVector<3>			vector_t;
	typedef AABox<vector_t>			box_t;

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
	typedef MathVector<3>			vector_t;
	typedef AABox<vector_t>			box_t;

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
		typedef ntree<tree_dim, world_dim, grid_elem_t*, NTreeGridData<world_dim> >	base_t;
		typedef typename NTreeGridData<world_dim>::position_attachment_t	position_attachment_t;

		lg_ntree(Grid& grid, position_attachment_t aPos) :
			m_gridData(grid, aPos)
		{}

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
