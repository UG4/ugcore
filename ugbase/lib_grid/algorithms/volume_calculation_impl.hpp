// created by Sebastian Reiter, Martin Stepniewski, Martin Scherer
// s.b.reiter@gmail.com
// Nov 28, 2013

#ifndef __H__UG__volume_calculation_impl__
#define __H__UG__volume_calculation_impl__

#include "common/math/misc/math_util.h"
#include "geom_obj_util/face_util.h"

namespace ug{


template <class TAAPos>
number CalculateVolume(Volume* elem, TAAPos aaPos)
{
	switch (elem->reference_object_id()) {
	case ROID_TETRAHEDRON:
		return CalculateVolume(static_cast<Tetrahedron*>(elem), aaPos);
	case ROID_PRISM:
		return CalculateVolume(static_cast<Prism*>(elem), aaPos);
	case ROID_PYRAMID:
		return CalculateVolume(static_cast<Pyramid*>(elem), aaPos);
	case ROID_HEXAHEDRON:
		return CalculateVolume(static_cast<Hexahedron*>(elem), aaPos);
	default:
		UG_THROW("Unknown volume type");
		break;
	}

	return NAN;
}

template <class TAAPos>
number CalculateVolume(Tetrahedron* elem, TAAPos aaPos)
{
	return CalculateTetrahedronVolume(aaPos[elem->vertex(0)],
									aaPos[elem->vertex(1)],
									aaPos[elem->vertex(2)],
									aaPos[elem->vertex(3)]);
}

template <class TAAPos>
number CalculateVolume(Pyramid* elem, TAAPos aaPos)
{
	return CalculatePyramidVolume(aaPos[elem->vertex(0)],
								aaPos[elem->vertex(1)],
								aaPos[elem->vertex(2)],
								aaPos[elem->vertex(3)],
								aaPos[elem->vertex(4)]);
}

template <class TAAPos>
number CalculateVolume(Prism* elem, TAAPos aaPos)
{
	return CalculatePrismVolume(aaPos[elem->vertex(0)],
								aaPos[elem->vertex(1)],
								aaPos[elem->vertex(2)],
								aaPos[elem->vertex(3)],
								aaPos[elem->vertex(4)],
								aaPos[elem->vertex(5)]);
}

template <class TAAPos>
number CalculateVolume(Hexahedron* elem, TAAPos aaPos)
{
	return CalculateHexahedronVolume(aaPos[elem->vertex(0)],
									aaPos[elem->vertex(1)],
									aaPos[elem->vertex(2)],
									aaPos[elem->vertex(3)],
									aaPos[elem->vertex(4)],
									aaPos[elem->vertex(5)],
									aaPos[elem->vertex(6)],
									aaPos[elem->vertex(7)]);
}


template <class TAAPos>
number CalculateVolume(Face* elem, TAAPos aaPos)
{
	return FaceArea(elem, aaPos);
}


template <class TAAPos>
number CalculateVolume(Edge* elem, TAAPos aaPos)
{
	return EdgeLength(elem, aaPos);
}


template <class TIterator, class TAAPos>
number CalculateVolume(TIterator begin, TIterator end, TAAPos aaPos)
{
	number totalVolume = 0;
	for(TIterator iter = begin; iter != end; ++iter){
		totalVolume += CalculateVolume(*iter, aaPos);
	}
	return totalVolume;
}

}// end of namespace

#endif
