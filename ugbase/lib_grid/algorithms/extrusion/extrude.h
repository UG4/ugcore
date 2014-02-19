// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y09 m05 d28

#ifndef __H__LIBGRID__EXTRUDE__
#define __H__LIBGRID__EXTRUDE__

#include <vector>
#include "lib_grid/lg_base.h"

namespace ug
{
/// \addtogroup lib_grid_algorithms_extrusion
///	@{

///	used to specify the behavior of ug::Extrude.
enum ExtrusionOptions
{
	EO_CREATE_FACES 	= 1,
	EO_CREATE_VOLUMES	= 1 << 1,

	EO_DEFAULT 			= 	EO_CREATE_FACES |
							EO_CREATE_VOLUMES
};

////////////////////////////////////////////////////////////////////////
//	Extrude
///	extrudes geometry and creates new edges, faces and volumes.
/**
 * For each vertex in pvVerticesInOut a new edge will be created.
 * NULL is a valid parameter for pvVerticesInOut.
 *
 * For each edge in pvEdgesInOut a new edge.
 * if EO_CREATE_FACES is enabled, a new face that connects the two
 * edges will be created too.
 * NULL is a valid parameter for pvEdgesInOut.
 *
 * For each face in pvFacesInOut a new face will be created.
 * if EO_CREATE_VOLUMES is enabled a new volume that connects the
 * two faces will be created too.
 * NULL is a valid parameter for pvFacesInOut.
 *
 * The element from which new elements are extruded is passed to
 * each new element as the parent element.
 *
 * Pass any or-combination of const in ExtrusionOptions as extrusionOptions.
 *
 * After the algorithm has finished, the in-out-vectors hold the elements
 * which have been created from the input-elements of those vectors.
 * The in-out-vectors thus can be directly used in a new call to Extrude.
 *
 * If you need to have access to all newly created elements you could use
 * a ug::Selector with enabled autoselection.
 */
void Extrude(Grid& grid,
			std::vector<Vertex*>* pvVerticesInOut,
			std::vector<EdgeBase*>* pvEdgesInOut,
			std::vector<Face*>* pvFacesInOut,
			const vector3& direction,
			uint extrusionOptions = EO_DEFAULT,
			APosition& aPos = aPosition);


/*
void Extrude(Grid& grid,
			Selector& sel,
			const vector3& direction,
			uint extrusionOptions = EO_DEFAULT,
			APosition& aPos = aPosition);
*/

/// @}

}//	end of namespace


#endif
