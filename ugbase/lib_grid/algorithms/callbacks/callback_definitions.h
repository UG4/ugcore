// created by Sebastian Reiter
// y10 m12 d13
// s.b.reiter@googlemail.com

#ifndef __H__LIB_GRID__CALLBACK_DEFINITIONS__
#define __H__LIB_GRID__CALLBACK_DEFINITIONS__

#include <boost/function.hpp>
#include "lib_grid/lg_base.h"

namespace ug
{

/**
 * Callbacks that allow algorithms to query whether they should consider
 * an element in their computations.
 *
 * \defgroup lib_grid_algorithms_callbacks callbacks
 * \ingroup lib_grid_algorithms
 * \{
 */
 
/**
 *\{
 *	\brief Callback definition used by several algorithms.
 *
 *	Allows to apply algorithms on arbitrary parts of a mesh.
 *	You may implement your own callbacks, which have to return true
 *	if the given geometric object should be considered in the algorithm.
 */
typedef boost::function<bool (VertexBase*)> Callback_ConsiderVertex;
typedef boost::function<bool (EdgeBase*)>	Callback_ConsiderEdge;
typedef boost::function<bool (Face*)>		Callback_ConsiderFace;
typedef boost::function<bool (Volume*)>		Callback_ConsiderVolume;

/** \} */

/**
 *\{
 *	\brief A callback that returns true for all given objects.
 */
inline bool ConsiderAllVertices(VertexBase*)	{return true;}
inline bool ConsiderAllEdges(EdgeBase*)			{return true;}
inline bool ConsiderAllFaces(Face*)				{return true;}
inline bool ConsiderAllVolumes(Volume*)			{return true;}
/** \} */


/** \} */	// end of group definition

}//	end of namespace

#endif
