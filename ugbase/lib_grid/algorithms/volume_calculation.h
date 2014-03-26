// created by Sebastian Reiter, Martin Stepniewski, Martin Scherer
// s.b.reiter@gmail.com
// Nov 28, 2013

#ifndef __H__UG__volume_calculation__
#define __H__UG__volume_calculation__

namespace ug{

/**
 * \brief	contains methods for the calculation of the volumes of
 * 			edge, face and volume elements
 *
 * \defgroup lib_grid_algorithms_volume_calculation volume_calculation
 * \ingroup lib_grid_algorithms
 * \{
 */

///	Calculates the volume of the given element.
/**	Make sure that aaPos has at least the same dimensionality as the specified element!
 * \note	If Volume-elements feature bended sides, the returned value is only an
 * 			approximation to the actual volume.
 * \{ */
//	VOLUMES
template <class TAAPos>
inline
number CalculateVolume(Volume* elem, TAAPos aaPos);

template <class TAAPos>
inline
number CalculateVolume(Tetrahedron* elem, TAAPos aaPos);

template <class TAAPos>
inline
number CalculateVolume(Pyramid* elem, TAAPos aaPos);

template <class TAAPos>
inline
number CalculateVolume(Prism* elem, TAAPos aaPos);

template <class TAAPos>
inline
number CalculateVolume(Hexahedron* elem, TAAPos aaPos);

//	FACES
template <class TAAPos>
inline
number CalculateVolume(FaceVertices* elem, TAAPos aaPos);

//	EDGES
template <class TAAPos>
inline
number CalculateVolume(EdgeVertices* elem, TAAPos aaPos);

template <class TAAPos>
inline
number CalculateVolume(Vertex* elem, TAAPos aaPos);

template <class TIterator, class TAAPos>
inline
number CalculateVolume(TIterator begin, TIterator end, TAAPos aaPos);
/** \} */


/** \} */	//end of ingroup

}// end of namespace


////////////////////////////////////////
//	include implementation
#include "volume_calculation_impl.hpp"

#endif
