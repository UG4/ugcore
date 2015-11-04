// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 22.12.2011 (m,d,y)

#ifndef __H__UG__grid_base_object_traits__
#define __H__UG__grid_base_object_traits__

#include "element_storage.h"
#include "generic_grid_object_iterator.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	The geometry_traits. This class can be specialized by each element-type.
/**
 * In order to use a custom geometric object with libGrid, you have to
 * supply a specialization for geometry_traits.
 * Specializations have to specify the following types and methods:
 *
 * MANDATORY:
 * Types:
 * - grid_base_object:		the geometric object from which TElem derives.
 * 					has to be either Vertex, Edge, Face or Volume.
 * - iterator:		An iterator that iterates over ElementContainer<BaseClass>
 * 					and which has a constructor that takes
 * 					ElementContainer<BaseClass>::iterator as an argument.
 * 					casts should be checked when in DEBUG mode!
 *
 * constants:
 * - CONTAINER_SECTION: This constant should hold the pipes section in which your objects should be placed after creation, starting from 0. See the Grid-documentation for more information.
 * - BASE_OBJECT_ID: Has to hold one of the GridBaseObjectId constants, or -1.
 *
 * OPTIONAL:
 * Types:
 * - Descriptor:	a class which can be passed to the constructor of the element.
 */
template <class TElem>
class geometry_traits
{};


////////////////////////////////////////////////////////////////////////////////
///	This Iterator will be used as base-class for iterators of specialized geometric objects.
typedef ElementStorage<Vertex>::SectionContainer::iterator			VertexIterator;
typedef ElementStorage<Vertex>::SectionContainer::const_iterator	ConstVertexIterator;

typedef ElementStorage<Edge>::SectionContainer::iterator			EdgeIterator;
typedef ElementStorage<Edge>::SectionContainer::const_iterator		ConstEdgeIterator;

typedef ElementStorage<Face>::SectionContainer::iterator				FaceIterator;
typedef ElementStorage<Face>::SectionContainer::const_iterator			ConstFaceIterator;

typedef ElementStorage<Volume>::SectionContainer::iterator				VolumeIterator;
typedef ElementStorage<Volume>::SectionContainer::const_iterator		ConstVolumeIterator;



////////////////////////////////////////////////////////////////////////////////
template <>
class geometry_traits<GridObject>
{
	public:

		enum
		{
			CONTAINER_SECTION = -1,
			BASE_OBJECT_ID = -1
		};
};

template <>
class geometry_traits<Vertex>
{
	public:
		typedef VertexIterator		iterator;
		typedef ConstVertexIterator	const_iterator;

		typedef Vertex	grid_base_object;

		enum
		{
			CONTAINER_SECTION = -1,
			BASE_OBJECT_ID = VERTEX
		};
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_VERTEX;
};

template <>
class geometry_traits<Edge>
{
	public:
		typedef EdgeIterator		iterator;
		typedef ConstEdgeIterator	const_iterator;

		typedef Edge	grid_base_object;
		typedef EdgeDescriptor GeneralDescriptor;

		enum
		{
			CONTAINER_SECTION = -1,
			BASE_OBJECT_ID = EDGE
		};
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_EDGE;
};

template <>
class geometry_traits<Face>
{
	public:
		typedef FaceIterator		iterator;
		typedef ConstFaceIterator	const_iterator;

		typedef Face	grid_base_object;
		typedef FaceDescriptor GeneralDescriptor;

		enum
		{
			CONTAINER_SECTION = -1,
			BASE_OBJECT_ID = FACE
		};
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_UNKNOWN;
};

template <>
class geometry_traits<Volume>
{
	public:
		typedef VolumeIterator			iterator;
		typedef ConstVolumeIterator		const_iterator;

		typedef Volume		grid_base_object;
		typedef VolumeDescriptor GeneralDescriptor;

		enum
		{
			CONTAINER_SECTION = -1,
			BASE_OBJECT_ID = VOLUME
		};
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_UNKNOWN;
};

}//	end of namespace

#endif
