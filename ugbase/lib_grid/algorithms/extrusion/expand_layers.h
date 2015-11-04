#ifndef __H__UG__expand_layers__
#define __H__UG__expand_layers__

#include <vector>
#include "lib_grid/lg_base.h"

namespace ug
{

/// Used to tell ExpandLayers_... which subsets should be regarded as layers.
/**
 * 	- subsetIndex defines the subset of the source (low dimensional) layer.
 *	- newSubsetIndex defines the subset into which the newly generated elements will go.
 *	- width describes the width to which a layer shall be expanded.
 */
struct FractureInfo{
	FractureInfo(int subsetInd, int newSubsetInd, double w) :
		subsetIndex(subsetInd), newSubsetIndex(newSubsetInd), width(w)	{}

	int subsetIndex;
	int newSubsetIndex;
	double width;
};

/**
 * This algorithm indirectly uses Grid::mark.
 *
 * 1 dimensional fractures specified in fracInfos are expanded to 2 dimensional subsets.
 * the resulting fractures will then consist of 2 layers of quadrilaterals. On the
 * boundaries triangles are inserted.
 *
 * Through expandFracBoundaries you can tell the algorithm whether inner fracture
 * boundaries shall be expanded. Note that this means that an additional node is
 * introduced at each inner fracture boundary vertex and that the associated
 * fracture elements are connected at two sides.
 * Note that fractures are always expanded at boundaries which lie on the geometries
 * boundary.
 *
 *	This algorithm requires the option FACEOPT_AUTOGENERATE_EDGES.
 *	The option is automatically enabled if required.
 */
bool ExpandFractures2d(Grid& grid, SubsetHandler& sh,
						const std::vector<FractureInfo>& fracInfos,
						bool expandInnerFracBnds, bool expandOuterFracBnds);


/**
 * This algorithm indirectly uses Grid::mark.
 *
 * 2 dimensional fractures specified in fracInfos are expanded to 3 dimensional subsets.
 * the resulting fractures will then consist of 2 layers of hexahedrons. On the
 * boundaries tetrahedrons, prisms and pyramids are inserted.
 *
 * Through expandFracBoundaries you can tell the algorithm whether inner fracture
 * boundaries shall be expanded. Note that this means that an additional node is
 * introduced at each inner fracture boundary vertex and that the associated
 * fracture elements are connected at two sides.
 * Note that fractures are always expanded at boundaries which lie on the geometries
 * boundary.
 *
 *	This algorithm requires the option FACEOPT_AUTOGENERATE_EDGES.
 *	The option is automatically enabled if required.
 *
 *	This algorithm requires the option VOLOPT_AUTOGENERATE_FACES.
 *	The option is automatically enabled if required.
 */
bool ExpandFractures3d(Grid& grid, SubsetHandler& sh,
						const std::vector<FractureInfo>& fracInfos,
						bool expandInnerFracBnds, bool expandOuterFracBnds);

}//	end of namespace

#endif
