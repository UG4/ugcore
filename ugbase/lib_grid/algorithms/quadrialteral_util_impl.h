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

#ifndef __H__UG_quadrialteral_util_impl
#define __H__UG_quadrialteral_util_impl

#include <set>
#include <boost/iterator/transform_iterator.hpp>
#include "quadrilateral_util.h"
#include "lib_grid/algorithms/geom_obj_util/vertex_util.h"
#include "common/util/hash.h"

namespace ug{

template <class face_iter_t, class TAAPos>
void ReplaceByQuadrilaterals_FaceBased(
		Grid& g,
		face_iter_t facesBegin,
		face_iter_t facesEnd,
		TAAPos aaPos)
{
	std::vector<Edge*> edges;
	GetInnerEdgesOfFaceSoup(edges, g, facesBegin, facesEnd);

	if(!edges.empty())
	{
		ReplaceByQuadrilaterals_EdgeBased (
				g,
				edges.begin(),
				edges.end(),
				aaPos);
	}
}

namespace detail{
	namespace quadUtil{
		struct EdgeToQuadInfo {
			EdgeToQuadInfo () : e(NULL), quality(-1) {};
			EdgeToQuadInfo (Edge* _e, number _q) : e(_e), quality(_q) {};

			typedef Edge* result_type;

			Edge* e;
			number quality;

			bool operator < (const EdgeToQuadInfo& other) const
			{return quality > other.quality;}

			inline Edge* operator() (const EdgeToQuadInfo& v) const
			{return v.e;}
		};
	}//	end of namespace
}//	end of namespace detail

template <class edge_iter_t, class TAAPos>
void ReplaceByQuadrilaterals_EdgeBased(
		Grid& g,
		edge_iter_t edgesBegin,
		edge_iter_t edgesEnd,
		TAAPos aaPos)
{
	using namespace detail::quadUtil;

	Grid::face_traits::secure_container faces;
	FaceDescriptor fd(4);
	std::multiset<EdgeToQuadInfo> ms;

	for(edge_iter_t eiter = edgesBegin; eiter != edgesEnd; ++eiter){
		Edge* e = *eiter;
		g.associated_elements(faces, e);
		if(		faces.size() == 2
			&&	faces[0]->num_vertices() == 3
			&& faces[1]->num_vertices() == 3)
		{
			int i0 = GetConnectedVertexIndex(faces[0], e);
			if(i0 != -1){
				fd.set_vertex(0, faces[0]->vertex(i0));
				fd.set_vertex(1, faces[0]->vertex((i0+1)%3));
				fd.set_vertex(2, GetConnectedVertex(e, faces[1]));
				fd.set_vertex(3, faces[0]->vertex((i0+2)%3));

				ms.insert(EdgeToQuadInfo(e, FaceQuality(&fd, aaPos)));
			}
		}
	}

	typedef boost::transform_iterator<
						EdgeToQuadInfo,
						std::multiset<EdgeToQuadInfo>::iterator>
			trans_iter_t;

	ReplaceByQuadrilaterals_EdgeBasedNoSort(
			g,
			trans_iter_t(ms.begin()),
			trans_iter_t(ms.end()));
}


template <class face_iter_t>
void ReplaceByQuadrilaterals_FaceBasedNoSort(
		Grid& g,
		face_iter_t facesBegin,
		face_iter_t facesEnd)
{
	std::vector<Edge*> edges;
	GetInnerEdgesOfFaceSoup(edges, g, facesBegin, facesEnd);

	if(!edges.empty())
	{
		ReplaceByQuadrilaterals_EdgeBasedNoSort (
				g,
				edges.begin(),
				edges.end());
	}
}


template <class edge_iter_t>
void ReplaceByQuadrilaterals_EdgeBasedNoSort(
		Grid& g,
		edge_iter_t edgesBegin,
		edge_iter_t edgesEnd)
{
	Grid::face_traits::secure_container faces;
	for(edge_iter_t eiter = edgesBegin; eiter != edgesEnd;)
	{
		Edge* e = *eiter;
		++eiter;
		g.associated_elements(faces, e);
		if(faces.size() == 2) {
			try{
			//NOTE: This may fail. faces[0] and faces[1] don't have to be
			//		triangles or may be identical etc. This would, however,
			//		raise an UGError.
				ReplaceByQuadrilateral(g, faces[0], faces[1]);
			}
			catch(UGError&) {}
		}
	}
}

}//	end of namespace

#endif	//__H__UG_quadrialteral_util_impl
