//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m12 d08

#ifndef __H__LIB_GRID__PARALLEL_GRID_LAYOUT__
#define __H__LIB_GRID__PARALLEL_GRID_LAYOUT__

#include <vector>
#include <list>
#include <map>
#include "pcl/pcl.h"
#include "lib_grid/lg_base.h"

//	specialize pcl::type_traits for VertexBase, EdgeBase, Face and Volume
namespace pcl
{
///	VertexBase interfaces and layouts store elements of type VertexBase*
template <>
struct type_traits<ug::VertexBase>
{
	typedef ug::VertexBase* Elem;
};

///	EdgeBase interfaces and layouts store elements of type VertexBase*
template <>
struct type_traits<ug::EdgeBase>
{
	typedef ug::EdgeBase* Elem;
};

///	Face interfaces and layouts store elements of type VertexBase*
template <>
struct type_traits<ug::Face>
{
	typedef ug::Face* Elem;
};

///	Volume interfaces and layouts store elements of type VertexBase*
template <>
struct type_traits<ug::Volume>
{
	typedef ug::Volume* Elem;
};

}//	end of namespace pcl

namespace ug
{

////////////////////////////////////////////////////////////////////////
///	The types of interface-entries.
enum InterfaceNodeTypes
{
	INT_UNKNOWN =	0,
	INT_MASTER =	1,
	INT_SLAVE =		1<<1,
	INT_LINK =		1<<2
};

//	declare vertex-, edge-, face- and volume-layouts
//	we're using std::list as interface-element container, since we
//	require interface-element-iterators that stay valid even if the
//	interface is altered.
typedef pcl::MultiLevelLayout<
		pcl::OrderedInterface<VertexBase, std::list> >	VertexLayout;
typedef pcl::MultiLevelLayout<
		pcl::OrderedInterface<EdgeBase, std::list> >	EdgeLayout;
typedef pcl::MultiLevelLayout<
		pcl::OrderedInterface<Face, std::list> >		FaceLayout;
typedef pcl::MultiLevelLayout<
		pcl::OrderedInterface<Volume, std::list> >		VolumeLayout;

//	declare GridLayoutMap
typedef pcl::LayoutMap<pcl::MultiLevelLayout,
						pcl::OrderedInterface,
						int,
						std::list>						GridLayoutMap;


}//	end of namespace

#endif
