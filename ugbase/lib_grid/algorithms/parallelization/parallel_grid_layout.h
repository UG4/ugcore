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
//	NodeLayout
///	specialization-scheme of pcl::Layout for geometric objects
template <class TGeomObj>
class NodeLayout : public pcl::Layout<pcl::BasicInterface<TGeomObj*> >
{
	public:
		typedef TGeomObj	GeomObj;
};

////////////////////////////////////////////////////////////////////////
///	The vertex layout stores pointers to vertices in the interfaces.
typedef NodeLayout<VertexBase>	VertexLayout;

////////////////////////////////////////////////////////////////////////
///	The edge layout stores pointers to edges in the interfaces.
typedef NodeLayout<EdgeBase>	EdgeLayout;

////////////////////////////////////////////////////////////////////////
///	The face layout stores pointers to faces in the interfaces.
typedef NodeLayout<Face>		FaceLayout;

////////////////////////////////////////////////////////////////////////
///	The volume layout stores pointers to volumes in the interfaces.
typedef NodeLayout<Volume>		VolumeLayout;


////////////////////////////////////////////////////////////////////////
//	NodeLayoutHierarchy
///	Holds a pcl::Layout for each level.
/**
 * May be used to store Layouts for the same element-type in multiple
 * layers.
 * 
 * TNodeLayout has to be a specialization of NodeLayout for either
 * VertexBase, EdgeBase, Face or Volume.
 *
 * \sa VertexLayout, EdgeLayout, FaceLayout, VolumeLayout.
 * \sa VertexLayoutHierarchy, EdgeLayoutHierarchy,
 * \sa FaceLayoutHierarchy, VolumeLayoutHierarchy.
 */
template <class TNodeLayout>
class NodeLayoutHierarchy
{	
	public:
		typedef TNodeLayout	Layout;
		typedef typename Layout::GeomObj	GeomObj;
		typedef typename Layout::Interface	Interface;
		
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
///	Stores VertexLayouts in different levels.
typedef NodeLayoutHierarchy<VertexLayout>	VertexLayoutHierarchy;

////////////////////////////////////////////////////////////////////////
///	Stores EdgeLayouts in different levels.
typedef NodeLayoutHierarchy<EdgeLayout>	EdgeLayoutHierarchy;

////////////////////////////////////////////////////////////////////////
///	Stores FaceLayouts in different levels
typedef NodeLayoutHierarchy<FaceLayout>		FaceLayoutHierarchy;

////////////////////////////////////////////////////////////////////////
///	Stores VolumeLayouts in different levels
typedef NodeLayoutHierarchy<VolumeLayout>		VolumeLayoutHierarchy;


////////////////////////////////////////////////////////////////////////
//	The ElementGroups for vertices, edges, faces and volumes don't do
//	much in the current implementation. Indeed they only define some
//	types, that can be used by the default implementation of pcl::group_traits.
/*
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
*/
////////////////////////////////////////////////////////////////////////
//	GridLayoutMap
///	the GridLayoutMap is an aggregation of element-layouts.
/**
 * Element layouts can be accessed using an integer key.
 * This allows to store layouts for different node-types.
 * E.g. master-layouts or slave-layouts.
 */
class GridLayoutMap
{
	public:
		typedef std::map<int, VertexLayoutHierarchy>	VertexLayoutHierarchyMap;
		typedef std::map<int, EdgeLayoutHierarchy>		EdgeLayoutHierarchyMap;
		typedef std::map<int, FaceLayoutHierarchy>		FaceLayoutHierarchyMap;
		typedef std::map<int, VolumeLayoutHierarchy>	VolumeLayoutHierarchyMap;
		
	public:
	///	returns true if the layout for the given key exists.
		inline bool has_vertex_layouts(int key)		{return m_vertexLayoutMap.find(key) != m_vertexLayoutMap.end();}
	///	returns true if the layout for the given key exists.
		inline bool has_edge_layouts(int key)		{return m_edgeLayoutMap.find(key) != m_edgeLayoutMap.end();}
	///	returns true if the layout for the given key exists.
		inline bool has_face_layouts(int key)		{return m_faceLayoutMap.find(key) != m_faceLayoutMap.end();}
	///	returns true if the layout for the given key exists.
		inline bool has_volume_layouts(int key)		{return m_volumeLayoutMap.find(key) != m_volumeLayoutMap.end();}

	///	returns the ParallelVertexLayout that is associated with key.
	/**	if the layout doesn't already exist then it will be created.*/
		inline VertexLayoutHierarchy& vertex_layout_hierarchy(int key)	{return m_vertexLayoutMap[key];}
	///	returns the ParallelEdgeLayout that is associated with key.
	/**	if the layout doesn't already exist then it will be created.*/
		inline EdgeLayoutHierarchy& edge_layout_hierarchy(int key)		{return m_edgeLayoutMap[key];}
	///	returns the ParallelFaceLayout that is associated with key.
	/**	if the layout doesn't already exist then it will be created.*/
		inline FaceLayoutHierarchy& face_layout_hierarchy(int key)		{return m_faceLayoutMap[key];}
	///	returns the ParallelVolumeLayout that is associated with key.
	/**	if the layout doesn't already exist then it will be created.*/
		inline VolumeLayoutHierarchy& volume_layout_hierarchy(int key)	{return m_volumeLayoutMap[key];}
		
	///	the layout map should only be accessed if it can't be avoided.
	/**	returns a std::map that holds pairs of int and VertexLayoutMaps.*/
		inline VertexLayoutHierarchyMap& vertex_layout_hierarchy_map()	{return m_vertexLayoutMap;}
	///	the layout map should only be accessed if it can't be avoided.
	/**	returns a std::map that holds pairs of int and EdgeLayoutMaps.*/
		inline EdgeLayoutHierarchyMap& edge_layout_hierarchy_map()		{return m_edgeLayoutMap;}
	///	the layout map should only be accessed if it can't be avoided.
	/**	returns a std::map that holds pairs of int and FaceLayoutMaps.*/
		inline FaceLayoutHierarchyMap& face_layout_hierarchy_map()		{return m_faceLayoutMap;}
	///	the layout map should only be accessed if it can't be avoided.
	/**	returns a std::map that holds pairs of int and VolumeLayoutMaps.*/
		inline VolumeLayoutHierarchyMap& volume_layout_hierarchy_map()	{return m_volumeLayoutMap;}
		
	protected:
		VertexLayoutHierarchyMap	m_vertexLayoutMap;
		EdgeLayoutHierarchyMap		m_edgeLayoutMap;
		FaceLayoutHierarchyMap		m_faceLayoutMap;
		VolumeLayoutHierarchyMap	m_volumeLayoutMap;
};


}//	end of namespace

#endif
