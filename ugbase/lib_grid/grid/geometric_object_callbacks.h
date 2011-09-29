// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 28.09.2011 (m,d,y) (originally created on y10 m12 d13)

#ifndef __H__UG__geometric_object_callbacks__
#define __H__UG__geometric_object_callbacks__

#include <boost/function.hpp>
#include "geometric_base_objects.h"

namespace ug
{

/**
 * Callbacks that allow algorithms to query whether they should consider
 * an element in their computations.
 *
 * \defgroup lib_grid_callbacks callbacks
 * \ingroup lib_grid
 * \{
 */

////////////////////////////////////////////////////////////////////////
/**
 *\{
 *	\brief Callback definition associating a boolean with a geometric object
 *
 *	Allows to apply algorithms on arbitrary parts of a mesh.
 *	You may implement your own callbacks, which have to return true
 *	if the given geometric object should be considered in the algorithm.
 */
typedef boost::function<bool (VertexBase*)> CB_ConsiderVertex;
typedef boost::function<bool (EdgeBase*)>	CB_ConsiderEdge;
typedef boost::function<bool (Face*)>		CB_ConsiderFace;
typedef boost::function<bool (Volume*)>		CB_ConsiderVolume;
/** \} */

////////////////////////////////////////////////////////////////////////
/**
 *\{
 *	\brief A callback that returns true for all given objects.
 */
inline bool ConsiderAllVertices(VertexBase*)	{return true;}
inline bool ConsiderAllEdges(EdgeBase*)			{return true;}
inline bool ConsiderAllFaces(Face*)				{return true;}
inline bool ConsiderAllVolumes(Volume*)			{return true;}

inline bool ConsiderAll(VertexBase*)	{return true;}
inline bool ConsiderAll(EdgeBase*)		{return true;}
inline bool ConsiderAll(Face*)			{return true;}
inline bool ConsiderAll(Volume*)		{return true;}

/** \} */

/** \} */	// end of group definition

}//	end of namespace

#endif
