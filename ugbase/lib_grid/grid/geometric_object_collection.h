// created by Sebastian Reiter
// y09 m07 d31
// s.b.reiter@googlemail.com

#ifndef __H__LIB_GRID__GEOMETRIC_OBJECT_COLLECTION__
#define __H__LIB_GRID__GEOMETRIC_OBJECT_COLLECTION__

#include <list>
#include "geometric_base_objects.h"
#include "common/util/section_container.h"

namespace ug
{

///	a section-container that holds geometric-objects in a std::list
typedef SectionContainer<GeometricObject*, std::list<GeometricObject*> >
		GeometricObjectSectionContainer;

////////////////////////////////////////////////////////////////////////
///	a helper class that holds a collection of possibly unconnected geometric-objects.
/**
 * This class is a simple helper class..
 * Its purpose is to make it easy to pass a collection of geometric-objects
 * to a function while maintaining the possibility to iterate over different
 * sub-types of geometric-objects separatly.
 *
 * Please note that a GeometricObjectCollection is only valid as long as
 * the object from which you received the collection still exists.
 *
 * A GeometricObjectCollection can only be queried for iterators and 
 * element-counts. You may not insert new elements or remove
 * old ones (at least not directly).
 *
 * Classes that you can query for their GeometricObjectCollections are for example
 * \sa Grid, \sa Selector, \sa SubsetHandler.
 *
 * As long as the object that provides the GeometricObjectCollection is still
 * valid, the GeometricObjectCollection will always hold the current
 * geometric objects of the source-object (i.e. a grid, a selector or a subset-handler).
 *
 * Please note that a GeometricObjectCollection does not represent a topological
 * closed part of the grid. A Collection can for example hold faces without their
 * associated vertices.
 *
 * How to use GeometricObjectCollection:
 * Once you retrieved an instance (let's call it goc) you can query it for
 * iterators like this:
 * VertexBaseIterator iter = goc.vertices_begin();
 * or if you want to iterate over triangles type the following:
 * TriangleIterator iter = goc.begin<Triangle>();
 *
 * if you want to get the number of hexahedrons you would go like this:
 * uint numHexas = goc.num<Hexahedron>();
 */
 /*
class GeometricObjectCollection
{
	public:
		GeometricObjectCollection();
									
		GeometricObjectCollection(	GeometricObjectSectionContainer* pVrtSection,
									GeometricObjectSectionContainer* pEdgeSection,
									GeometricObjectSectionContainer* pFaceSection,
									GeometricObjectSectionContainer* pVolSection);
		
	//	copy constructor.
		GeometricObjectCollection(const GeometricObjectCollection& goc);
		
		GeometricObjectCollection& operator =(const GeometricObjectCollection& goc);
		
	//	Iterators
	//	begin
		template <class TGeomObj>
		typename geometry_traits<TGeomObj>::iterator
		begin();

	//	end
		template <class TGeomObj>
		typename geometry_traits<TGeomObj>::iterator
		end();

		inline VertexBaseIterator	vertices_begin()	{return begin<VertexBase>();}
		inline VertexBaseIterator	vertices_end()		{return end<VertexBase>();}
		inline EdgeBaseIterator		edges_begin()		{return begin<EdgeBase>();}
		inline EdgeBaseIterator		edges_end()			{return end<EdgeBase>();}
		inline FaceIterator			faces_begin()		{return begin<Face>();}
		inline FaceIterator			faces_end()			{return end<Face>();}
		inline VolumeIterator		volumes_begin()		{return begin<Volume>();}
		inline VolumeIterator		volumes_end()		{return end<Volume>();}
		
	//	element numbers
		template <class TGeomObj>
		uint num();
		
		inline uint num_vertices()	{return num<VertexBase>();}
		inline uint num_edges()		{return num<EdgeBase>();}
		inline uint num_faces()		{return num<Face>();}
		inline uint num_volumes()	{return num<Volume>();}
		
	protected:
		GeometricObjectSectionContainer*	m_pSectionContainers[4];
};
*/
////////////////////////////////////////////////////////////////////////
//	GeometricObjectCollection
///	a helper class that holds a collection of possibly unconnected geometric-objects.
/**
 * This class is a simple helper class..
 * Its purpose is to make it easy to pass a collection of geometric-objects
 * to a function while maintaining the possibility to iterate over different
 * sub-types of geometric-objects seperatly.
 *
 * In contrary to \sa GeometricObjectCollection, the
 * GeometricObjectCollection allows access to the elements through
 * different levels.
 *
 * Please note that a GeometricObjectCollection is only valid as long as
 * the object from which you received the collection still exists.
 *
 * A GeometricObjectCollection can only be queried for iterators and 
 * element-counts. You may not insert new elements or remove
 * old ones (at least not directly).
 *
 * Classes that you can query for their GeometricObjectCollection
 * are for example \sa MultiGrid, \sa MGSelector.
 *
 * As long as the object that provides the GeometricObjectCollection
 * is still valid, the GeometricObjectCollection will always hold the current
 * geometric objects of the source-object (i.e. a grid, a selector or a subset-handler).
 *
 * Please note that a GeometricObjectCollection does not represent
 * a topological closed part of the grid.
 * A Collection can for example hold faces without their
 * associated vertices.
 *
 * How to use GeometricObjectCollection:
 * Once you retrieved an instance (let's call it mlgoc) you can query it for
 * iterators like this:
 * VertexBaseIterator iter = mlgoc.vertices_begin(1);
 * or if you want to iterate over triangles of level 1 type the following:
 * TriangleIterator iter = goc.begin<Triangle>(1);
 *
 * if you want to get the number of hexahedrons in level 0 you would go like this:
 * uint numHexas = goc.num<Hexahedron>(0);
 */
class GeometricObjectCollection
{
	public:
	///	initializes the instance with an estimate of the number of levels.
	/**	The estimate does not have to match exactly. However, if it does
	 *  it makes things faster.*/
		GeometricObjectCollection(size_t levelEstimate = 1);
		
	///	initializes level 0 with the given sections.
		GeometricObjectCollection(GeometricObjectSectionContainer* pVrtSection,
											GeometricObjectSectionContainer* pEdgeSection,
											GeometricObjectSectionContainer* pFaceSection,
											GeometricObjectSectionContainer* pVolSection);

	//	copy constructor.
		GeometricObjectCollection(const GeometricObjectCollection& mgoc);
		
		GeometricObjectCollection& operator =(const GeometricObjectCollection& mgoc);
		
	///	only used during creation by the methods that create the collection
		void add_level(	GeometricObjectSectionContainer* pVrtSection,
						GeometricObjectSectionContainer* pEdgeSection,
						GeometricObjectSectionContainer* pFaceSection,
						GeometricObjectSectionContainer* pVolSection);

	///	returns the number of levels
		inline size_t num_levels() const		{return m_levels.size();}
//		inline const GeometricObjectCollection& get_geometric_object_collection(int level)	{return m_levels[level];}
		
	//	Iterators
	//	begin
		template <class TGeomObj>
		inline
		typename geometry_traits<TGeomObj>::iterator
		begin(size_t level = 0);

	//	end
		template <class TGeomObj>
		inline
		typename geometry_traits<TGeomObj>::iterator
		end(size_t level = 0);

		inline VertexBaseIterator	vertices_begin(size_t level = 0)	{return begin<VertexBase>(level);}
		inline VertexBaseIterator	vertices_end(size_t level = 0)		{return end<VertexBase>(level);}
		inline EdgeBaseIterator		edges_begin(size_t level = 0)		{return begin<EdgeBase>(level);}
		inline EdgeBaseIterator		edges_end(size_t level = 0)			{return end<EdgeBase>(level);}
		inline FaceIterator			faces_begin(size_t level = 0)		{return begin<Face>(level);}
		inline FaceIterator			faces_end(size_t level = 0)			{return end<Face>(level);}
		inline VolumeIterator		volumes_begin(size_t level = 0)		{return begin<Volume>(level);}
		inline VolumeIterator		volumes_end(size_t level = 0)		{return end<Volume>(level);}
		
	//	element numbers
		template <class TGeomObj>
		size_t num();
		
		inline size_t num_vertices()	{return num<VertexBase>();}
		inline size_t num_edges()		{return num<EdgeBase>();}
		inline size_t num_faces()		{return num<Face>();}
		inline size_t num_volumes()		{return num<Volume>();}
		
		template <class TGeomObj>
		inline
		size_t num(size_t level);
		
		inline size_t num_vertices(size_t level)	{return num<VertexBase>(level);}
		inline size_t num_edges(size_t level)		{return num<EdgeBase>(level);}
		inline size_t num_faces(size_t level)		{return num<Face>(level);}
		inline size_t num_volumes(size_t level)		{return num<Volume>(level);}
		
	protected:
		void assign(const GeometricObjectCollection& goc);

		template <class TGeomObj> inline
		GeometricObjectSectionContainer* get_container(size_t level);
		
	protected:
		struct ContainerCollection{
			ContainerCollection()	{}
			ContainerCollection(GeometricObjectSectionContainer* vrtCon,
								GeometricObjectSectionContainer* edgeCon,
								GeometricObjectSectionContainer* faceCon,
								GeometricObjectSectionContainer* volCon);

			GeometricObjectSectionContainer*	pSectionContainers[4];
		};
		
		typedef std::vector<ContainerCollection> ContainerVec;
		//typedef std::vector<GeometricObjectCollection> GOCVec;
		
	protected:
		ContainerVec	m_levels;
};

}//	end of namespace

////////////////////////////////////////////////
//	include implementation
#include "geometric_object_collection_impl.hpp"

#endif
