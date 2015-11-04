#ifndef __H__LIB_GRID__GRID_CONSTANTS__
#define __H__LIB_GRID__GRID_CONSTANTS__

namespace ug
{

/// \addtogroup lib_grid
/// @{

////////////////////////////////////////////////////////////////////////
//	VertexOptions
///	Used to specify the way in which Grid manages vertex-specific data.
enum VertexOptions
{
	VRTOPT_NONE = 0x00000000,
	VRTOPT_STORE_ASSOCIATED_EDGES = 0x00000001,
	VRTOPT_STORE_ASSOCIATED_FACES = 0x00000002,
	VRTOPT_STORE_ASSOCIATED_VOLUMES = 0x00000004,
};

////////////////////////////////////////////////////////////////////////
//	EdgeOptions
///	Used to specify the way in which Grid manages edge-specific data.
enum EdgeOptions
{
	EDGEOPT_NONE = 0x00000000,
	EDGEOPT_STORE_ASSOCIATED_FACES = 0x00000100,
	EDGEOPT_STORE_ASSOCIATED_VOLUMES = 0x00000200,
};

////////////////////////////////////////////////////////////////////////
//	FaceOptions
///	Used to specify the way in which Grid manages face-specific data.
enum FaceOptions
{
	FACEOPT_NONE = 0x00000000,
	FACEOPT_STORE_ASSOCIATED_EDGES = 0x00010000,	///< minor speed-improvement for grid.get_edge(Face*, int)
	FACEOPT_STORE_ASSOCIATED_VOLUMES = 0x00020000,
	FACEOPT_AUTOGENERATE_EDGES = 0x00080000
};

////////////////////////////////////////////////////////////////////////
//	VolumeOptions
///	Used to specify the way in which Grid manages volume-specific data.
enum VolumeOptions
{
	VOLOPT_NONE = 0x00000000,
	VOLOPT_STORE_ASSOCIATED_EDGES = 0x01000000,		///< minor speed-improvement for grid.get_edge(Volume*, int)
	VOLOPT_STORE_ASSOCIATED_FACES = 0x02000000,		///< speed-improvement for grid.get_face(Face*, int) ~15%
	VOLOPT_AUTOGENERATE_EDGES = 0x08000000,
	VOLOPT_AUTOGENERATE_FACES = 0x10000000
};

////////////////////////////////////////////////////////////////////////
//	GridOptions
///	Specify how references between associated objects are stored in a grid.
enum GridOptions
{
	GRIDOPT_NONE = 0x00000000,
	GRIDOPT_NO_INTERCONNECTION = 0x00000000,
///	vertices store lists of associated geometric objects.
/**	Note that this is the minimal required interconnection for many dynamic
 * algorithms (i.e. deleting an object from a grid automatically enables this
 * option.*/
	GRIDOPT_VERTEXCENTRIC_INTERCONNECTION =	  VRTOPT_STORE_ASSOCIATED_EDGES
											| VRTOPT_STORE_ASSOCIATED_FACES
											| VRTOPT_STORE_ASSOCIATED_VOLUMES,

///	sides are automatically created
	GRIDOPT_AUTOGENERATE_SIDES =	  FACEOPT_AUTOGENERATE_EDGES
									| VOLOPT_AUTOGENERATE_FACES,

///	All elements store references to associated lower dimensional geometric objects
/**	Additionally GRIDOPT_VERTEXCENTRIC_INTERCONNECTION is used.*/
	GRIDOPT_STANDARD_INTERCONNECTION =	  GRIDOPT_VERTEXCENTRIC_INTERCONNECTION
										| GRIDOPT_AUTOGENERATE_SIDES
										| FACEOPT_STORE_ASSOCIATED_EDGES
										| VOLOPT_STORE_ASSOCIATED_EDGES
										| VOLOPT_STORE_ASSOCIATED_FACES,

///	All elements store references to all associated elements
/**	This includes GRIDOPT_VERTEXCENTRIC_INTERCONNECTION and
 * GRIDOPT_STANDARD_INTERCONNECTION.*/
	GRIDOPT_FULL_INTERCONNECTION =	  GRIDOPT_STANDARD_INTERCONNECTION
									| EDGEOPT_STORE_ASSOCIATED_FACES
									| EDGEOPT_STORE_ASSOCIATED_VOLUMES
									| FACEOPT_STORE_ASSOCIATED_VOLUMES,

	GRIDOPT_DEFAULT =	GRIDOPT_VERTEXCENTRIC_INTERCONNECTION
					  | GRIDOPT_AUTOGENERATE_SIDES
};

/// @}
}

#endif
