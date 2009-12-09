//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m12 d08

#ifndef __H__LIB_GRID__PARALLEL_GRID_LAYOUT__
#define __H__LIB_GRID__PARALLEL_GRID_LAYOUT__

#include <vector>
#include <map>
#include "pcl/pcl.h"
#include "lib_grid/lg_base.h"

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
	
////////////////////////////////////////////////////////////////////////
//	ParallelNodeLayout
///	The parallel node layout holds a pcl::Layout on each of its levels.
/**
 * The parallel node layout may be used to store Layouts for the same
 * local-id-type in multiple layers.
 * Please note that the ParallelNodeLayout will not be used as a
 * plc compliant layout - instead grants access to pcl comliant layouts.
 * \sa ParallelVertexLayout, ParallelEdgeLayout, ParallelFaceLayout, ParallelVolumeLayout
 */
template <class TLocalID>
class ParallelNodeLayout
{	
	public:
		typedef TLocalID				LocalID;
		typedef std::vector<TLocalID>	Interface;
		typedef pcl::Layout<Interface>	Layout;
		
	public:
	///	returns the number of levels.
		inline int num_levels()					{return m_vLayouts.size();}

	///	returns the layout at the given level.
	/**	If level >= num_levels() then the layouts in between
		will be automatically created.*/
		inline Layout& layout(int level = 0)	{require_level(level); return *m_vLayouts[level];}
		
	protected:
	///	adds num new levels.
		inline void new_levels(int num)			{for(int i = 0; i < num; ++i) m_vLayouts.push_back(new Layout);}
		
	///	if the required level doesn't exist yet, it will created.
		inline void require_level(int level)	{if(level >= num_levels()) new_levels(level - num_levels() + 1);}
		
	protected:
		std::vector<Layout*>	m_vLayouts;
};

////////////////////////////////////////////////////////////////////////
///	The vertex layout uses pointers to vertices as local id.
typedef ParallelNodeLayout<VertexBase*>	ParallelVertexLayout;

////////////////////////////////////////////////////////////////////////
///	The edge layout uses pointers to edges as local id.
typedef ParallelNodeLayout<EdgeBase*>	ParallelEdgeLayout;

////////////////////////////////////////////////////////////////////////
///	The face layout uses pointers to faces as local id.
typedef ParallelNodeLayout<Face*>	ParallelFaceLayout;

////////////////////////////////////////////////////////////////////////
///	The volume layout uses pointers to volumes as local id.
typedef ParallelNodeLayout<Volume*>	ParallelVolumeLayout;


////////////////////////////////////////////////////////////////////////
//	The ElementGroups for vertices, edges, faces and volumes don't do
//	much in the current implementation. Indeed they only define some
//	types, that can be used by the default implementation of pcl::group_traits.
template <class TElement>
class ParallelNodeGroup
{
	public:
		typedef TElement	Element;
		typedef TElement*	LocalID;
		typedef typename ParallelNodeLayout<LocalID>::Interface Interface;
		typedef typename ParallelNodeLayout<LocalID>::Layout	Layout;
};

typedef ParallelNodeGroup<VertexBase>	ParallelVertexGroup;
typedef ParallelNodeGroup<EdgeBase>		ParallelEdgeGroup;
typedef ParallelNodeGroup<Face>			ParallelFaceGroup;
typedef ParallelNodeGroup<Volume>		ParallelVolumeGroup;

////////////////////////////////////////////////////////////////////////
//	ParallelGridLayout
///	the ParallelGridLayout is an aggregation of element-layouts.
/**
 * Element layouts can be accessed using an integer key.
 * This allows to store layouts for different node-types.
 * E.g. master-layouts or slave-layouts.
 */
class ParallelGridLayout
{
	public:
		typedef std::map<int, ParallelVertexLayout>	VertexLayoutMap;
		typedef std::map<int, ParallelEdgeLayout>	EdgeLayoutMap;
		typedef std::map<int, ParallelFaceLayout>	FaceLayoutMap;
		typedef std::map<int, ParallelVolumeLayout>	VolumeLayoutMap;
		
	public:
	///	returns true if the layout for the given key exists.
		inline bool has_vertex_layout(int key)		{return m_vertexLayoutMap.find(key) != m_vertexLayoutMap.end();}
	///	returns true if the layout for the given key exists.
		inline bool has_edge_layout(int key)		{return m_edgeLayoutMap.find(key) != m_edgeLayoutMap.end();}
	///	returns true if the layout for the given key exists.
		inline bool has_face_layout(int key)		{return m_faceLayoutMap.find(key) != m_faceLayoutMap.end();}
	///	returns true if the layout for the given key exists.
		inline bool has_volume_layout(int key)		{return m_volumeLayoutMap.find(key) != m_volumeLayoutMap.end();}

	///	returns the ParallelVertexLayout that is associated with key.
	/**	if the layout doesn't already exist then it will be created.*/
		inline ParallelVertexLayout& vertex_layout(int key)		{return m_vertexLayoutMap[key];}
	///	returns the ParallelEdgeLayout that is associated with key.
	/**	if the layout doesn't already exist then it will be created.*/
		inline ParallelEdgeLayout& edge_layout(int key)			{return m_edgeLayoutMap[key];}
	///	returns the ParallelFaceLayout that is associated with key.
	/**	if the layout doesn't already exist then it will be created.*/
		inline ParallelFaceLayout& face_layout(int key)			{return m_faceLayoutMap[key];}
	///	returns the ParallelVolumeLayout that is associated with key.
	/**	if the layout doesn't already exist then it will be created.*/
		inline ParallelVolumeLayout& volume_layout(int key)		{return m_volumeLayoutMap[key];}
		
	///	the layout map should only be accessed if it can't be avoided.
	/**	returns a std::map that holds pairs of int and VertexLayoutMaps.*/
		inline VertexLayoutMap& vertex_layout_map()				{return m_vertexLayoutMap;}
	///	the layout map should only be accessed if it can't be avoided.
	/**	returns a std::map that holds pairs of int and EdgeLayoutMaps.*/
		inline EdgeLayoutMap& edge_layout_map()					{return m_edgeLayoutMap;}
	///	the layout map should only be accessed if it can't be avoided.
	/**	returns a std::map that holds pairs of int and FaceLayoutMaps.*/
		inline FaceLayoutMap& face_layout_map()					{return m_faceLayoutMap;}
	///	the layout map should only be accessed if it can't be avoided.
	/**	returns a std::map that holds pairs of int and VolumeLayoutMaps.*/
		inline VolumeLayoutMap& volume_layout_map()				{return m_volumeLayoutMap;}
		
	protected:
		VertexLayoutMap		m_vertexLayoutMap;
		EdgeLayoutMap		m_edgeLayoutMap;
		FaceLayoutMap		m_faceLayoutMap;
		VolumeLayoutMap		m_volumeLayoutMap;
};


}//	end of namespace

#endif
