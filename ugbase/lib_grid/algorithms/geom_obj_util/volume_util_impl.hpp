/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__LIB_GRID__VOLUME_UTIL_IMPL__
#define __H__LIB_GRID__VOLUME_UTIL_IMPL__

#include <algorithm>
#include <vector>
#include "lib_grid/lg_base.h"
#include "common/static_assert.h"
#include "common/util/vec_for_each.h"
#include "lib_grid/grid_objects/hexahedron_rules.h"
#include "lib_grid/grid_objects/prism_rules.h"
#include "lib_grid/grid_objects/pyramid_rules.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
template <typename TAAPos>
bool
ContainsPoint(Volume* vol, const vector3& p, TAAPos aaPos)
{
//	iterate over face descriptors of the sides and check whether the point
//	lies inside or outside
	FaceDescriptor fd;
	vector3 n, dir;

//	to minimize rounding errors we'll compare against a small-constant which is
//	relative to the first edge of the examined volume.
	EdgeDescriptor ed;
	vol->edge_desc(0, ed);
	number len = EdgeLength(&ed, aaPos);

	// the constant should be relative to the same geometric measure as what it is
	// compared against later on, i.e. length*area, since otherwise problems arise
	// with geometries scaled to very small extensions;
	// which is why I changed sqrt(lenSq) to lenSq^1.5 (mbreit, 2015-05-11)
	const number locSmall = len * len * len * SMALL;

	for(size_t i = 0; i < vol->num_faces(); ++i){
		vol->face_desc(i, fd);
		CalculateNormalNoNormalize(n, &fd, aaPos);
		VecSubtract(dir, aaPos[fd.vertex(0)], p);

		if(VecDot(dir, n) < -locSmall)
			return false;
	}
	return true;
}

////////////////////////////////////////////////////////////////////////
//	PointIsInsideTetrahedron
inline bool
PointIsInsideTetrahedron(const vector3& v, Tetrahedron* tet,
						 Grid::VertexAttachmentAccessor<APosition>& aaPos)
{
	return PointIsInsideTetrahedron(v, aaPos[tet->vertex(0)], aaPos[tet->vertex(1)],
									aaPos[tet->vertex(2)], aaPos[tet->vertex(3)]);
}

////////////////////////////////////////////////////////////////////////
//	CalculateCenter
template<typename TVertexPositionAttachmentAccessor>
typename TVertexPositionAttachmentAccessor::ValueType
CalculateCenter(const VolumeVertices* vol, TVertexPositionAttachmentAccessor& aaPosVRT)
{
	typename TVertexPositionAttachmentAccessor::ValueType v;
//	init v with 0.
	VecSet(v, 0);

	size_t numVrts = vol->num_vertices();
	VolumeVertices::ConstVertexArray vrts = vol->vertices();

//	sum up
	for(size_t i = 0; i < numVrts; ++i)
	{
		VecAdd(v, v, aaPosVRT[vrts[i]]);
	}

//	average
	if(numVrts > 0)
		VecScale(v, v, 1./(number)numVrts);

	return v;
}

////////////////////////////////////////////////////////////////////////
template<typename TAAPosVRT, typename TAAWeightVRT>
UG_API
typename TAAPosVRT::ValueType
CalculateCenter(const VolumeVertices* vol, TAAPosVRT& aaPos, TAAWeightVRT& aaWeight)
{
	typename TAAPosVRT::ValueType v;
	using weight_t = typename TAAWeightVRT::ValueType;

//	init v with 0.
	VecSet(v, 0);

	size_t numVrts = vol->num_vertices();
	VolumeVertices::ConstVertexArray vrts = vol->vertices();

//	sum up
	weight_t totalWeight = 0;
	for(size_t i = 0; i < numVrts; ++i)
	{
		weight_t w = aaWeight[vrts[i]];
		VecScaleAppend(v, w, aaPos[vrts[i]]);
		totalWeight += w;
	}

//	average
	if(totalWeight != 0)
		VecScale(v, v, 1./(number)totalWeight);

	return v;
}


///	Can be used to compare vertices of their grids through their hash-value.
template <typename TElem>
class CmpVrtsByHash{
	public:
		CmpVrtsByHash(TElem* e) : m_e(e) {};
		bool operator () (int i0, int i1) {
			return m_e->vertex(i0)->get_hash_value() < m_e->vertex(i1)->get_hash_value();
		}
	private:
		TElem* m_e;
};


template <typename TVolIter>
void ConvertToTetrahedra (
		Grid& grid,
		TVolIter volsBegin,
		TVolIter volsEnd)
{
	using namespace std;

//	first we'll collect all quadrilaterals that are connected to selected
//	volumes. Avoid 'grid.mark' here, since it may be used by 'grid.associated_elements'
	Grid::face_traits::secure_container	faces;
	vector<Face*>	quads;

	for(TVolIter iv = volsBegin; iv != volsEnd; ++iv){
		grid.associated_elements(faces, *iv);
		for(size_t _vfeI = 0; _vfeI < faces.size(); ++_vfeI){ Face* f = faces[_vfeI];{
			if(f->num_vertices() == 4)
				quads.push_back(f);
		}};
	}

//	remove double entries
	grid.begin_marking();
	size_t offset = 0;

	for(size_t i = 0; i + offset < quads.size();){
		if(offset > 0)
			quads[i] = quads[i + offset];
		
		if(!grid.is_marked(quads[i])){
			grid.mark(quads[i]);
			++i;
		}
		else{
			++offset;
		}
	}

	if(offset > 0)
		quads.resize(quads.size() - offset);

	grid.end_marking();

	for(size_t _vfeI = 0; _vfeI < quads.size(); ++_vfeI){ Face* f = quads[_vfeI];{
//todo	in a parallel environment, global id's should be compared here
		CmpVrtsByHash<Face> cmp(f);
	//	get the smallest vertex of the face
		int smallest = 0;
		for(int i = 1; i < 4; ++i){
			if(cmp(i, smallest))
				smallest = i;
		}

		int i0 = smallest;
		int i1 = (smallest + 1) % 4;
		int i2 = (smallest + 2) % 4;
		int i3 = (smallest + 3) % 4;
		grid.create<Triangle>(TriangleDescriptor(f->vertex(i0), f->vertex(i1), f->vertex(i2)), f);
		grid.create<Triangle>(TriangleDescriptor(f->vertex(i2), f->vertex(i3), f->vertex(i0)), f);
	}};


//	now convert the given volume-elements
	UG_STATIC_ASSERT((prism_rules::MAX_NUM_CONVERT_TO_TETS_INDS_OUT < 
							hex_rules::MAX_NUM_CONVERT_TO_TETS_INDS_OUT) &&
					 	(pyra_rules::MAX_NUM_CONVERT_TO_TETS_INDS_OUT <
							hex_rules::MAX_NUM_CONVERT_TO_TETS_INDS_OUT),
					 HEX_RULES_MAX_NUM_CONVERT_TO_TETS_INDS_OUT__considered_to_be_highest_among_prism_pyra_and_hex);
	static constexpr int arrayLen = hex_rules::MAX_NUM_CONVERT_TO_TETS_INDS_OUT;
	int inds[arrayLen];
	
	vector<Volume*> volsToErase;
	for(TVolIter iv = volsBegin; iv != volsEnd; ++iv){
		Volume* vol = *iv;
		const ReferenceObjectID roid = vol->reference_object_id();
		CmpVrtsByHash<Volume> cmp(vol);
		size_t numEntries = 0;

		switch(roid){
			case ReferenceObjectID::ROID_PYRAMID:
				numEntries = pyra_rules::ConvertToTetrahedra(inds, cmp);
				break;

			case ReferenceObjectID::ROID_PRISM:
				numEntries = prism_rules::ConvertToTetrahedra(inds, cmp);
				break;

			case ReferenceObjectID::ROID_HEXAHEDRON:
				UG_THROW("ConvertToTetrahedra for hexahedra not yet implemented!");
				break;
			default:
				break;
		}

		if(numEntries > 0){
			volsToErase.push_back(vol);
			size_t i = 0;
			Volume::ConstVertexArray vrts = vol->vertices();
			while(i < numEntries){
				int goid = inds[i++];
				UG_COND_THROW(goid != GridObjectID::GOID_TETRAHEDRON,
							  "Only tetrahedra may result from ConvertToTetrahedra");
				int i0 = inds[i++];
				int i1 = inds[i++];
				int i2 = inds[i++];
				int i3 = inds[i++];
				grid.create<Tetrahedron>(
					TetrahedronDescriptor(vrts[i0], vrts[i1], vrts[i2], vrts[i3]),
					vol);
			}
		}
	}

//	finally erase all unused volumes and quadrilaterals
	for(size_t _vfeI = 0; _vfeI < volsToErase.size(); ++_vfeI){ Volume* v = volsToErase[_vfeI];{
		grid.erase(v);
	}};

	Grid::volume_traits::secure_container	assVols;
	for(size_t _vfeI = 0; _vfeI < quads.size(); ++_vfeI){ Face* f = quads[_vfeI];{
		grid.associated_elements(assVols, f);
		if(assVols.empty())
			grid.erase(f);
	}};
}

}//	end of namespace

#endif
