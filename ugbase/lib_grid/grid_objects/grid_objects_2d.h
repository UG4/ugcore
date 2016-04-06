/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG__grid_objects_2d__
#define __H__UG__grid_objects_2d__

#include "../grid/grid.h"
#include "common/math/ugmath.h"
#include "common/assert.h"
#include "grid_objects_0d.h"
#include "grid_objects_1d.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
///	These numbers define where in the face-section-container a face will be stored.
/**	The order of the constants must not be changed! Algorithms may exist that rely on it.*/
enum FaceContainerSections
{
	CSFACE_NONE = -1,
	CSFACE_TRIANGLE = 0,
	CSFACE_QUADRILATERAL = 1,
	CSFACE_CONSTRAINED_TRIANGLE = 2,
	CSFACE_CONSTRAINED_QUADRILATERAL = 3,
	CSFACE_CONSTRAINING_TRIANGLE = 4,
	CSFACE_CONSTRAINING_QUADRILATERAL = 5,

	CSFACE_USER	// always last
};


////////////////////////////////////////////////////////////////////////////////
//	NORMAL FACES
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//	TriangleDescriptor
///	only used to initialize a triangle. for all other tasks you should use FaceDescriptor.
class UG_API TriangleDescriptor
{
	public:
		TriangleDescriptor()	{}
		TriangleDescriptor(const TriangleDescriptor& td);
		TriangleDescriptor(Vertex* v1, Vertex* v2, Vertex* v3);

		inline uint num_vertices() const					{return 3;}
		inline void set_vertex(uint index, Vertex* v)	{m_vertex[index] = v;}
		inline Vertex* vertex(size_t index) const			{return m_vertex[index];}

	protected:
		Vertex*	m_vertex[3];
};



////////////////////////////////////////////////////////////////////////
//	CustomTriangle
///	Concrete types share this base-type. It is not intended for direct use.
/**
 * BaseClass has to be derived from Face (or simply should be Face).
 * The ConcreteTriangleType is used in methods like refine, etc. as the type
 * of newly created objects.
 */
template <class ConcreteTriangleType, class BaseClass>
class UG_API CustomTriangle : public BaseClass
{
	public:
		static const size_t NUM_VERTICES = 3;

	public:
		CustomTriangle()	{}
		CustomTriangle(const TriangleDescriptor& td);
		CustomTriangle(Vertex* v1, Vertex* v2, Vertex* v3);

		virtual GridObject* create_empty_instance() const	{return new ConcreteTriangleType;}
		virtual ReferenceObjectID reference_object_id() const {return ROID_TRIANGLE;}

		virtual Vertex* vertex(size_t index) const	{return m_vertices[index];}
		virtual Face::ConstVertexArray vertices() const		{return m_vertices;}
		virtual size_t num_vertices() const	{return 3;}

		virtual EdgeDescriptor edge_desc(int index) const
			{return EdgeDescriptor(m_vertices[index], m_vertices[(index+1) % 3]);}

		virtual void edge_desc(int index, EdgeDescriptor& edOut) const
			{edOut.set_vertices(m_vertices[index], m_vertices[(index+1) % 3]);}

		virtual std::pair<GridBaseObjectId, int> get_opposing_object(Vertex* vrt) const;

	///	Refines a Triangle by inserting new vertices. \sa Face::refine.
		virtual bool refine(std::vector<Face*>& vNewFacesOut,
							Vertex** newFaceVertexOut,
							Vertex** newEdgeVertices,
							Vertex* newFaceVertex = NULL,
							Vertex** pSubstituteVertices = NULL,
							int snapPointIndex = -1);

		virtual bool collapse_edge(std::vector<Face*>& vNewFacesOut,
								int edgeIndex, Vertex* newVertex,
								Vertex** pSubstituteVertices = NULL);

		virtual bool collapse_edges(std::vector<Face*>& vNewFacesOut,
								std::vector<Vertex*>& vNewEdgeVertices,
								Vertex** pSubstituteVertices = NULL);

//	BEGIN Depreciated
		virtual void create_faces_by_edge_split(int splitEdgeIndex,
							Vertex* newVertex,
							std::vector<Face*>& vNewFacesOut,
							Vertex** pSubstituteVertices = NULL);

	protected:
		virtual void set_vertex(uint index, Vertex* pVrt)	{m_vertices[index] = pVrt;}

	protected:
		Vertex* m_vertices[3];
};



////////////////////////////////////////////////////////////////////////
//	Triangle
///	the most simple form of a face
/**
 *
 * \ingroup lib_grid_grid_objects
 */
class UG_API Triangle : public CustomTriangle<Triangle, Face>
{
	typedef CustomTriangle<Triangle, Face> BaseClass;
	public:
		inline static bool type_match(GridObject* pObj)	{return dynamic_cast<Triangle*>(pObj) != NULL;}

		Triangle() : BaseClass()	{}
		Triangle(const TriangleDescriptor& td) : BaseClass(td)	{}
		Triangle(Vertex* v1, Vertex* v2, Vertex* v3) : BaseClass(v1, v2, v3)	{}

		virtual int container_section() const	{return CSFACE_TRIANGLE;}

	protected:
		virtual Edge* create_edge(int index)
			{
				return new RegularEdge(m_vertices[index], m_vertices[(index+1) % 3]);
			}
};

template <>
class geometry_traits<Triangle>
{
	public:
		typedef GenericGridObjectIterator<Triangle*, FaceIterator>			iterator;
		typedef ConstGenericGridObjectIterator<Triangle*, FaceIterator,
															ConstFaceIterator>	const_iterator;

		typedef TriangleDescriptor Descriptor;	///< Faces can't be created directly
		typedef Face	grid_base_object;

		enum
		{
			CONTAINER_SECTION = CSFACE_TRIANGLE,
			BASE_OBJECT_ID = FACE
		};
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_TRIANGLE;
};

typedef geometry_traits<Triangle>::iterator			TriangleIterator;
typedef geometry_traits<Triangle>::const_iterator	ConstTriangleIterator;



////////////////////////////////////////////////////////////////////////
//	QuadrilateralDescriptor
///	only used to initialize a quadrilateral. for all other tasks you should use FaceDescriptor.
class UG_API QuadrilateralDescriptor
{
	public:
		QuadrilateralDescriptor()	{}
		QuadrilateralDescriptor(const QuadrilateralDescriptor& qd);
		QuadrilateralDescriptor(Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4);

		inline uint num_vertices() const					{return 4;}
		inline void set_vertex(uint index, Vertex* v)	{m_vertex[index] = v;}
		inline Vertex* vertex(size_t index) const			{return m_vertex[index];}

	protected:
		Vertex*	m_vertex[4];
};



////////////////////////////////////////////////////////////////////////
//	CustomQuadrilateral
///	Concrete types share this base-type. It is not intended for direct use.
/**
 * BaseClass has to be derived from Face (or simply should be Face).
 * The ConcreteQuadrilateralType is used in methods like refine, etc. as the type
 * of newly created objects.
 */
template <class ConcreteQuadrilateralType, class BaseClass>
class UG_API CustomQuadrilateral : public BaseClass
{
	public:
		static const size_t NUM_VERTICES = 4;

	public:
		using typename BaseClass::ConstVertexArray;

		CustomQuadrilateral()	{}
		CustomQuadrilateral(const QuadrilateralDescriptor& qd);
		CustomQuadrilateral(Vertex* v1, Vertex* v2,
							Vertex* v3, Vertex* v4);

		virtual GridObject* create_empty_instance() const	{return new ConcreteQuadrilateralType;}
		virtual ReferenceObjectID reference_object_id() const {return ROID_QUADRILATERAL;}

		virtual Vertex* vertex(size_t index) const	{return m_vertices[index];}
		virtual Face::ConstVertexArray vertices() const		{return m_vertices;}
		virtual size_t num_vertices() const	{return 4;}

		virtual EdgeDescriptor edge_desc(int index) const
			{return EdgeDescriptor(m_vertices[index], m_vertices[(index+1) % 4]);}

		virtual void edge_desc(int index, EdgeDescriptor& edOut) const
			{edOut.set_vertices(m_vertices[index], m_vertices[(index+1) % 4]);}


	///	fills the edge-descriptor with the edge that lies opposed to the specified one
	/**	If the specified edge is not part of the face, false is returned.*/
		virtual bool get_opposing_side(EdgeVertices* e, EdgeDescriptor& edOut) const;

		virtual std::pair<GridBaseObjectId, int> get_opposing_object(Vertex* vrt) const;

	///	Refines a Quadrilateral by inserting new vertices. \sa Face::refine.
		virtual bool refine(std::vector<Face*>& vNewFacesOut,
							Vertex** newFaceVertexOut,
							Vertex** newEdgeVertices,
							Vertex* newFaceVertex = NULL,
							Vertex** pSubstituteVertices = NULL,
							int snapPointIndex = -1);

		virtual bool collapse_edge(std::vector<Face*>& vNewFacesOut,
								int edgeIndex, Vertex* newVertex,
								Vertex** pSubstituteVertices = NULL);

		virtual bool collapse_edges(std::vector<Face*>& vNewFacesOut,
								std::vector<Vertex*>& vNewEdgeVertices,
								Vertex** pSubstituteVertices = NULL);

//	BEGIN Depreciated
		virtual void create_faces_by_edge_split(int splitEdgeIndex,
							Vertex* newVertex,
							std::vector<Face*>& vNewFacesOut,
							Vertex** pSubstituteVertices = NULL);

	protected:
		virtual void set_vertex(uint index, Vertex* pVrt)	{m_vertices[index] = pVrt;}

	protected:
		Vertex* m_vertices[4];
};



////////////////////////////////////////////////////////////////////////
//	Quadrilateral
///	a face with four points.
/**
 * \ingroup lib_grid_grid_objects
 */
class UG_API Quadrilateral : public CustomQuadrilateral<Quadrilateral, Face>
{
	public:
		typedef CustomQuadrilateral<Quadrilateral, Face> BaseClass;

		inline static bool type_match(GridObject* pObj)	{return dynamic_cast<Quadrilateral*>(pObj) != NULL;}

		Quadrilateral()	{}
		Quadrilateral(const QuadrilateralDescriptor& td) : BaseClass(td)	{}
		Quadrilateral(Vertex* v1, Vertex* v2,
					  Vertex* v3, Vertex* v4) : BaseClass(v1, v2, v3, v4)	{}

		virtual int container_section() const	{return CSFACE_QUADRILATERAL;}

	protected:
		virtual Edge* create_edge(int index)
		{
			return new RegularEdge(m_vertices[index], m_vertices[(index+1) % 4]);
		}
};

template <>
class geometry_traits<Quadrilateral>
{
	public:
		typedef GenericGridObjectIterator<Quadrilateral*, FaceIterator>		iterator;
		typedef ConstGenericGridObjectIterator<Quadrilateral*, FaceIterator,
																ConstFaceIterator>	const_iterator;

		typedef QuadrilateralDescriptor Descriptor;	///< Faces can't be created directly
		typedef Face	grid_base_object;

		enum
		{
			CONTAINER_SECTION = CSFACE_QUADRILATERAL,
			BASE_OBJECT_ID = FACE
		};
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_QUADRILATERAL;
};

typedef geometry_traits<Quadrilateral>::iterator		QuadrilateralIterator;
typedef geometry_traits<Quadrilateral>::const_iterator	ConstQuadrilateralIterator;



////////////////////////////////////////////////////////////////////////////////
//	CONSTRAINED FACES
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//	ConstrainedFace
///	This class stores the constraining object.
/**
 * Please note, that the user is has to link and unlink constraining
 * objects manually.
 */
class UG_API ConstrainedFace : public Face
{
	public:
		inline static bool type_match(GridObject* pObj)	{return dynamic_cast<ConstrainedFace*>(pObj) != NULL;}

		ConstrainedFace() : m_pConstrainingObject(NULL), m_parentBaseObjectId(-1)	{}
		virtual ~ConstrainedFace()	{}

		inline void set_constraining_object(GridObject* pObj)
		{
			m_pConstrainingObject = pObj;
			if(pObj)
				m_parentBaseObjectId = pObj->base_object_id();
		}

		inline GridObject* get_constraining_object()			{return m_pConstrainingObject;}

		inline int get_parent_base_object_id()				{return m_parentBaseObjectId;}
		inline void set_parent_base_object_id(int id)
		{
			if((m_parentBaseObjectId != -1) && (m_parentBaseObjectId != id)){
				UG_THROW("Bad parent base object id specified! The given id"
						" has to match the id of the constraining object if that"
						" is present. Call this method only, if no constraining"
						" object has been set!");
			}
			m_parentBaseObjectId = id;
		}

		virtual bool is_constrained() const							{return true;}

	protected:
		GridObject*	m_pConstrainingObject;
		int					m_parentBaseObjectId;
};


////////////////////////////////////////////////////////////////////////
//	ConstrainedTriangle
///	a triangle constrained by another object.
/**
 * \ingroup lib_grid_grid_objects
 */
class UG_API ConstrainedTriangle : public CustomTriangle<ConstrainedTriangle, ConstrainedFace>
{
	typedef CustomTriangle<ConstrainedTriangle, ConstrainedFace> BaseTriangle;

	public:
		inline static bool type_match(GridObject* pObj)	{return dynamic_cast<ConstrainedTriangle*>(pObj) != NULL;}

		ConstrainedTriangle() :
			BaseTriangle()	{}

		ConstrainedTriangle(const TriangleDescriptor& td) :
			BaseTriangle(td)	{}

		ConstrainedTriangle(Vertex* v1, Vertex* v2, Vertex* v3) :
			BaseTriangle(v1, v2, v3)	{}

		virtual int container_section() const	{return CSFACE_CONSTRAINED_TRIANGLE;}

	protected:
		virtual Edge* create_edge(int index)
			{
				return new ConstrainedEdge(m_vertices[index], m_vertices[(index+1) % 3]);
			}
};

template <>
class geometry_traits<ConstrainedTriangle>
{
	public:
		typedef GenericGridObjectIterator<ConstrainedTriangle*, FaceIterator>		iterator;
		typedef ConstGenericGridObjectIterator<ConstrainedTriangle*, FaceIterator,
																	ConstFaceIterator>	const_iterator;

		typedef TriangleDescriptor Descriptor;	///< Faces can't be created directly
		typedef Face	grid_base_object;

		enum
		{
			CONTAINER_SECTION = CSFACE_CONSTRAINED_TRIANGLE,
			BASE_OBJECT_ID = FACE
		};
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_TRIANGLE;
};

typedef geometry_traits<ConstrainedTriangle>::iterator			ConstrainedTriangleIterator;
typedef geometry_traits<ConstrainedTriangle>::const_iterator	ConstConstrainedTriangleIterator;



////////////////////////////////////////////////////////////////////////
//	ConstrainedQuadrilateral
///	a quadrilateral constrained by another object.
/**
 * \ingroup lib_grid_grid_objects
 */
class UG_API ConstrainedQuadrilateral : public CustomQuadrilateral<ConstrainedQuadrilateral, ConstrainedFace>
{
	typedef CustomQuadrilateral<ConstrainedQuadrilateral, ConstrainedFace> BaseClass;

	public:
		inline static bool type_match(GridObject* pObj)	{return dynamic_cast<ConstrainedQuadrilateral*>(pObj) != NULL;}

		ConstrainedQuadrilateral() : BaseClass()	{}
		ConstrainedQuadrilateral(const QuadrilateralDescriptor& qd) : BaseClass(qd)	{}
		ConstrainedQuadrilateral(Vertex* v1, Vertex* v2,
								 Vertex* v3, Vertex* v4) : BaseClass(v1, v2, v3, v4)	{}

		virtual int container_section() const	{return CSFACE_CONSTRAINED_QUADRILATERAL;}

	protected:
		virtual Edge* create_edge(int index)
			{
				return new ConstrainedEdge(m_vertices[index], m_vertices[(index+1) % 4]);
			}
};


template <>
class geometry_traits<ConstrainedQuadrilateral>
{
	public:
		typedef GenericGridObjectIterator<ConstrainedQuadrilateral*, FaceIterator>	iterator;
		typedef ConstGenericGridObjectIterator<ConstrainedQuadrilateral*,
													FaceIterator, ConstFaceIterator>	const_iterator;

		typedef QuadrilateralDescriptor Descriptor;	///< Faces can't be created directly
		typedef Face	grid_base_object;

		enum
		{
			CONTAINER_SECTION = CSFACE_CONSTRAINED_QUADRILATERAL,
			BASE_OBJECT_ID = FACE
		};
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_QUADRILATERAL;
};

typedef geometry_traits<ConstrainedQuadrilateral>::iterator			ConstrainedQuadrilateralIterator;
typedef geometry_traits<ConstrainedQuadrilateral>::const_iterator	ConstConstrainedQuadrilateralIterator;



////////////////////////////////////////////////////////////////////////////////
//	CONSTRAINING FACES
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//	ConstrainingFace
///	This class is used to store constrained geometric objects.
/**
 * Please note, that the user is has to link and unlink constrained
 * objects manually.
 */
class UG_API ConstrainingFace : public Face
{

	public:
		inline static bool type_match(GridObject* pObj)	{return dynamic_cast<ConstrainingFace*>(pObj) != NULL;}

		virtual ~ConstrainingFace()	{}

		virtual bool is_constraining() const					{return true;}

		inline void add_constrained_object(Vertex* pObj)
			{
				UG_ASSERT(!is_constrained_object(pObj), "vertex is already registered at constraining face");
					m_constrainedVertices.push_back(pObj);
			}

		inline void add_constrained_object(Edge* pObj)
			{
				UG_ASSERT(!is_constrained_object(pObj), "edge is already registered at constraining face");
					m_constrainedEdges.push_back(pObj);
			}

		inline void add_constrained_object(Face* pObj)
			{
				UG_ASSERT(!is_constrained_object(pObj), "face is already registered at constraining face");
					m_constrainedFaces.push_back(pObj);
			}

		inline bool is_constrained_object(Vertex* vrt)
			{
				std::vector<Vertex*>::iterator iter = find(m_constrainedVertices.begin(),
															m_constrainedVertices.end(), vrt);
				return iter != m_constrainedVertices.end();
			}

		inline bool is_constrained_object(Edge* edge)
			{
				std::vector<Edge*>::iterator iter = find(m_constrainedEdges.begin(),
															m_constrainedEdges.end(), edge);
				return iter != m_constrainedEdges.end();
			}

		inline bool is_constrained_object(Face* face)
			{
				std::vector<Face*>::iterator iter = find(m_constrainedFaces.begin(),
														m_constrainedFaces.end(), face);
				return iter != m_constrainedFaces.end();
			}

		inline void unconstrain_object(const Vertex* vrt)
			{
				std::vector<Vertex*>::iterator iter = find(m_constrainedVertices.begin(),
															 m_constrainedVertices.end(), vrt);
				if(iter != m_constrainedVertices.end())
					m_constrainedVertices.erase(iter);
			}

		inline void unconstrain_object(const Edge* edge)
			{
				std::vector<Edge*>::iterator iter = find(m_constrainedEdges.begin(),
															m_constrainedEdges.end(), edge);
				if(iter != m_constrainedEdges.end())
					m_constrainedEdges.erase(iter);
			}

		inline void unconstrain_object(const Face* face)
			{
				std::vector<Face*>::iterator iter = find(m_constrainedFaces.begin(),
														m_constrainedFaces.end(), face);
				if(iter != m_constrainedFaces.end())
					m_constrainedFaces.erase(iter);
			}

		inline void clear_constrained_vertices()	{m_constrainedVertices.clear();}
		inline void clear_constrained_edges()		{m_constrainedEdges.clear();}
		inline void clear_constrained_faces()		{m_constrainedFaces.clear();}
		inline void clear_constrained_objects()
			{
				clear_constrained_vertices();
				clear_constrained_edges();
				clear_constrained_faces();
			}

		inline size_t num_constrained_vertices() const	{return m_constrainedVertices.size();}
		inline size_t num_constrained_edges() const		{return m_constrainedEdges.size();}
		inline size_t num_constrained_faces() const		{return m_constrainedFaces.size();}

		template <class TElem> size_t num_constrained() const;


		inline Vertex* constrained_vertex(size_t ind) const
			{
				UG_ASSERT(ind < m_constrainedVertices.size(), "bad index.");
				return m_constrainedVertices[ind];
			}

		inline Edge* constrained_edge(size_t ind) const
			{
				UG_ASSERT(ind < m_constrainedEdges.size(), "bad index.");
				return m_constrainedEdges[ind];
			}

		inline Face* constrained_face(size_t ind) const
			{
				UG_ASSERT(ind < m_constrainedFaces.size(), "bad index.");
				return m_constrainedFaces[ind];
			}

		template <class TElem> TElem* constrained(size_t ind) const;

	protected:
		std::vector<Vertex*>	m_constrainedVertices;
		std::vector<Edge*>		m_constrainedEdges;
		std::vector<Face*>			m_constrainedFaces;
};



////////////////////////////////////////////////////////////////////////
//	ConstrainingTriangle
///	a triangle constraining other objects.
/**
 * \ingroup lib_grid_grid_objects
 */
class UG_API ConstrainingTriangle : public CustomTriangle<ConstrainingTriangle, ConstrainingFace>
{
	typedef CustomTriangle<ConstrainingTriangle, ConstrainingFace> BaseTriangle;

	public:
		inline static bool type_match(GridObject* pObj)	{return dynamic_cast<ConstrainingTriangle*>(pObj) != NULL;}

		ConstrainingTriangle() :
			BaseTriangle()	{reserve_memory();}
		ConstrainingTriangle(const TriangleDescriptor& td) :
			BaseTriangle(td)	{reserve_memory();}
		ConstrainingTriangle(Vertex* v1, Vertex* v2, Vertex* v3) :
			BaseTriangle(v1, v2, v3)	{reserve_memory();}

		virtual int container_section() const	{return CSFACE_CONSTRAINING_TRIANGLE;}

	protected:
		void reserve_memory()
			{
				m_constrainedEdges.reserve(3);
				m_constrainedFaces.reserve(4);
			}

		virtual Edge* create_edge(int index)
			{
				return new RegularEdge(m_vertices[index], m_vertices[(index+1) % 3]);
			}
};

template <>
class geometry_traits<ConstrainingTriangle>
{
	public:
		typedef GenericGridObjectIterator<ConstrainingTriangle*, FaceIterator>		iterator;
		typedef ConstGenericGridObjectIterator<ConstrainingTriangle*, FaceIterator,
																	ConstFaceIterator>	const_iterator;

		typedef TriangleDescriptor Descriptor;	///< Faces can't be created directly
		typedef Face	grid_base_object;

		enum
		{
			CONTAINER_SECTION = CSFACE_CONSTRAINING_TRIANGLE,
			BASE_OBJECT_ID = FACE
		};
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_TRIANGLE;
};

typedef geometry_traits<ConstrainingTriangle>::iterator			ConstrainingTriangleIterator;
typedef geometry_traits<ConstrainingTriangle>::const_iterator	ConstConstrainingTriangleIterator;



////////////////////////////////////////////////////////////////////////
//	ConstrainingQuadrilateral
///	a quadrilateral constraining other objects.
/**
 * \ingroup lib_grid_grid_objects
 */
class UG_API ConstrainingQuadrilateral : public CustomQuadrilateral<ConstrainingQuadrilateral, ConstrainingFace>
{
	typedef CustomQuadrilateral<ConstrainingQuadrilateral, ConstrainingFace> BaseClass;

	public:
		inline static bool type_match(GridObject* pObj)	{return dynamic_cast<ConstrainingQuadrilateral*>(pObj) != NULL;}

		ConstrainingQuadrilateral() :
			BaseClass()	{reserve_memory();}
		ConstrainingQuadrilateral(const QuadrilateralDescriptor& qd) :
			BaseClass(qd)	{reserve_memory();}
		ConstrainingQuadrilateral(Vertex* v1, Vertex* v2,
								  Vertex* v3, Vertex* v4) :
			BaseClass(v1, v2, v3, v4)	{reserve_memory();}

		virtual int container_section() const	{return CSFACE_CONSTRAINING_QUADRILATERAL;}

	protected:
		void reserve_memory()
			{
				m_constrainedVertices.reserve(1);
				m_constrainedEdges.reserve(4);
				m_constrainedFaces.reserve(4);
			}

		virtual Edge* create_edge(int index)
			{
				return new RegularEdge(m_vertices[index], m_vertices[(index+1) % 4]);
			}
};

template <>
class geometry_traits<ConstrainingQuadrilateral>
{
	public:
		typedef GenericGridObjectIterator<ConstrainingQuadrilateral*, FaceIterator>	iterator;
		typedef ConstGenericGridObjectIterator<ConstrainingQuadrilateral*,
														FaceIterator, ConstFaceIterator>	const_iterator;

		typedef QuadrilateralDescriptor Descriptor;	///< Faces can't be created directly
		typedef Face	grid_base_object;

		enum
		{
			CONTAINER_SECTION = CSFACE_CONSTRAINING_QUADRILATERAL,
			BASE_OBJECT_ID = FACE
		};
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_QUADRILATERAL;
};

typedef geometry_traits<ConstrainingQuadrilateral>::iterator		ConstrainingQuadrilateralIterator;
typedef geometry_traits<ConstrainingQuadrilateral>::const_iterator	ConstConstrainingQuadrilateralIterator;


}//	end of namespace

#endif
