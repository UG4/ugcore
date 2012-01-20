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
 *	Check out boost::function for more information on how such callbacks should
 *	look like.
 */
typedef boost::function<bool (VertexBase*)> CB_ConsiderVertex;
typedef boost::function<bool (EdgeBase*)>	CB_ConsiderEdge;
typedef boost::function<bool (Face*)>		CB_ConsiderFace;
typedef boost::function<bool (Volume*)>		CB_ConsiderVolume;
/** \} */

////////////////////////////////////////////////////////////////////////
/**
 *\{
 *	\brief Callback definition used to execute arbitrary code on an object
 *
 *	Allows to apply algorithms on arbitrary parts of a mesh.
 *	You may implement your own callbacks, which can execute arbitrary code
 *	on the given element. Check out boost::function for more information
 *	on how such callbacks should look like.
 */
typedef boost::function<void (VertexBase*)> CB_VisitVertex;
typedef boost::function<void (EdgeBase*)>	CB_VisitEdge;
typedef boost::function<void (Face*)>		CB_VisitFace;
typedef boost::function<void (Volume*)>		CB_VisitVolume;
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
