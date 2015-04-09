//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m02 d19

#ifndef __H__LIB_GRID__VOLUME_UTIL_IMPL__
#define __H__LIB_GRID__VOLUME_UTIL_IMPL__

#include <algorithm>
#include <vector>
#include "lib_grid/lg_base.h"
#include "common/util/vec_for_each.h"
#include "lib_grid/grid_objects/hexahedron_rules.h"
#include "lib_grid/grid_objects/prism_rules.h"
#include "lib_grid/grid_objects/pyramid_rules.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
template <class TAAPos>
bool
ContainsPoint(Volume* vol, const vector3& p, TAAPos aaPos)
{
	using std::max;
//	iterate over face descriptors of the sides and check whether the point
//	lies inside or outside
	FaceDescriptor fd;
	vector3 n, dir;

//	to minimize rouding errors we'll compare against a small-constant which is
//	relative to the longest edge of the examined volume.
	number lenSq = 0;
	EdgeDescriptor ed;
	for(size_t i = 0; i < vol->num_edges(); ++i){
		vol->edge_desc(i, ed);
		lenSq = max(lenSq, EdgeLengthSq(&ed, aaPos));
	}
	const number locSmall = sqrt(lenSq) * SMALL;

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
template<class TVertexPositionAttachmentAccessor>
typename TVertexPositionAttachmentAccessor::ValueType
CalculateCenter(const VolumeVertices* vol, TVertexPositionAttachmentAccessor& aaPosVRT)
{
	typename TVertexPositionAttachmentAccessor::ValueType v;
//	init v with 0.
	VecSet(v, 0);

	uint numVrts = vol->num_vertices();
	VolumeVertices::ConstVertexArray vrts = vol->vertices();

//	sum up
	for(uint i = 0; i < numVrts; ++i)
	{
		VecAdd(v, v, aaPosVRT[vrts[i]]);
	}

//	average
	if(numVrts > 0)
		VecScale(v, v, 1./(number)numVrts);

	return v;
}

////////////////////////////////////////////////////////////////////////
template<class TAAPosVRT, class TAAWeightVRT>
UG_API
typename TAAPosVRT::ValueType
CalculateCenter(const VolumeVertices* vol, TAAPosVRT& aaPos, TAAWeightVRT& aaWeight)
{
	typename TAAPosVRT::ValueType v;
	typedef typename TAAWeightVRT::ValueType weight_t;

//	init v with 0.
	VecSet(v, 0);

	uint numVrts = vol->num_vertices();
	VolumeVertices::ConstVertexArray vrts = vol->vertices();

//	sum up
	weight_t totalWeight = 0;
	for(uint i = 0; i < numVrts; ++i)
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

////////////////////////////////////////////////////////////////////////
//	CheckOrientation
template<class TAAPosVRT>
bool
CheckOrientation(Volume* vol, TAAPosVRT& aaPosVRT)
{
//	some typedefs
	typedef typename TAAPosVRT::ValueType vector_t;
	
//	First calculate the center of the volume
	vector_t volCenter = CalculateCenter(vol, aaPosVRT);
	
//	now check for each side whether it points away from the center.
	size_t numFaces = vol->num_faces();
	FaceDescriptor fd;
	vector_t normal;
	for(size_t i = 0; i < numFaces; ++i){
		vol->face_desc(i, fd);
		CalculateNormal(normal, &fd, aaPosVRT);
		
	//	in order to best approximate quadrilateral faces, we'll calculate the
	//	center of the face and compare that to the volCenter.
		vector_t faceCenter = CalculateCenter(&fd, aaPosVRT);
		
	//	now compare normal and center
		vector_t dir;
		VecSubtract(dir, faceCenter, volCenter);
		if(VecDot(dir, normal) < 0)
			return false;
	}
	
//	all center / normal checks succeeded. Orientation is fine.
	return true;
}

template<class TVolIterator, class TAAPosVRT>
int
FixOrientation(Grid& grid, TVolIterator volsBegin, TVolIterator volsEnd,
			   TAAPosVRT& aaPosVRT)
{
	int numFlips = 0;
//	iterate through all volumes
	for(VolumeIterator iter = volsBegin; iter != volsEnd; ++iter){
	//	check whether the orientation is fine
		if(!CheckOrientation(*iter, aaPosVRT)){
			grid.flip_orientation(*iter);
			++numFlips;
		}
	}
	
	return numFlips;
}


///	Can be used to compare vertices of their grids through their hash-value.
template <class TElem>
class CmpVrtsByHash{
	public:
		CmpVrtsByHash(TElem* e) : m_e(e) {};
		bool operator () (int i0, int i1) {
			return m_e->vertex(i0)->get_hash_value() < m_e->vertex(i1)->get_hash_value();
		}
	private:
		TElem* m_e;
};


template <class TVolIter>
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
		for_each_in_vec(Face* f, faces){
			if(f->num_vertices() == 4)
				quads.push_back(f);
		}end_for;
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

	for_each_in_vec(Face* f, quads){
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
	}end_for;


//	now convert the given volume-elements
	static const int arrayLen = max(pyra_rules::MAX_NUM_CONVERT_TO_TETS_INDS_OUT,
									max(prism_rules::MAX_NUM_CONVERT_TO_TETS_INDS_OUT,
										hex_rules::MAX_NUM_CONVERT_TO_TETS_INDS_OUT));
	int inds[arrayLen];
	
	vector<Volume*> volsToErase;
	for(TVolIter iv = volsBegin; iv != volsEnd; ++iv){
		Volume* vol = *iv;
		const ReferenceObjectID roid = vol->reference_object_id();
		CmpVrtsByHash<Volume> cmp(vol);
		size_t numEntries = 0;

		switch(roid){
			case ROID_PYRAMID:
				numEntries = pyra_rules::ConvertToTetrahedra(inds, cmp);
				break;

			case ROID_PRISM:
				numEntries = prism_rules::ConvertToTetrahedra(inds, cmp);
				break;

			case ROID_HEXAHEDRON:
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
				UG_COND_THROW(goid != GOID_TETRAHEDRON,
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
	for_each_in_vec(Volume* v, volsToErase){
		grid.erase(v);
	}end_for;

	Grid::volume_traits::secure_container	assVols;
	for_each_in_vec(Face* f, quads){
		grid.associated_elements(assVols, f);
		if(assVols.empty())
			grid.erase(f);
	}end_for;
}

}//	end of namespace

#endif
