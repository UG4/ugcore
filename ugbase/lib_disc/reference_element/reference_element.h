/*
 * reference_element.h
 *
 *  Created on: 13.05.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__REFERENCE_ELEMENT__REFERENCE_ELEMENT__
#define __H__UG__LIB_DISC__REFERENCE_ELEMENT__REFERENCE_ELEMENT__

#include <cassert>
#include <iostream>
#include <sstream>
#include "common/common.h"
#include "common/math/ugmath.h"
#include "lib_grid/grid/geometric_base_objects.h"

namespace ug{

///////////////////////////////////////////////////////////////////////////////
// Reference Element Common Base Class
///////////////////////////////////////////////////////////////////////////////

/// base class for reference elements
/**
 * Reference element interface. A reference element describes in local
 * coordinates the structure an element type. Physical elements of a grid are
 * thought to be constructed by a mapping from a reference element into the
 * real world space.
 *
 * Each ReferenceElement may be constructed from other (lower dimensional)
 * geometric objects, that are themselves a mapping from a (lower dimensional)
 * reference element. (E.g. a triangle is constructed by three edges and
 * three vertices) Thus, these relationships are also specified by the reference
 * element and methods of this function provide the number of constructing
 * sub-geometric objects and the relationship between those.
 *
 * Note, that we use one base class in order to implement all reference elements
 * providing enough space to store all data for each derivation. This enlarges
 * the memory consumption slightly but allows fast and inlined code. Since
 * usually only one reference element (singleton) per program is created memory
 * consumption is not an issue.
 */
class ReferenceElement
{
	public:
	///	Constructor filling the arrays
		ReferenceElement();

	/// returns the reference object id
		ReferenceObjectID roid() const
			{return m_vRefElemType[m_dim][0];}

	/// returns the dimension where reference element lives
		int dimension() const {return m_dim;}

	/// returns the size (e.g. area or volume) of the reference element
		number size() const	{return m_size;}

	/// returns the number of geometric objects of dim
	/**
	 * A reference element is constructed by several geometric objects, that
	 * are mapped by a reference element by themselves. This method returns how
	 * many (sub-)geometric objects of a given dimension are contained in this
	 * reference element.
	 *
	 * \param[in]	dim		dimension
	 * \returns		number of objects of the dimension contained in the ref elem
	 */
		size_t num(int dim) const	{return m_vNum[dim];}

	/// returns the number of object of dim for a sub-geometric object
	/**
	 * A reference element is constructed by several geometric objects, that
	 * are mapped by a reference element by themselves. This method returns how
	 * many (sub-)geometric objects of a given dimension are contained in a
	 * (sub-)geometric object of this reference element.
	 *
	 * \param[in]	dim_i		dimension of sub geometric object
	 * \param[in]	i			number of sub geometric object
	 * \param[in]	dim_j		dimension for elems contained in the sub-object
	 * \returns		number of objects of the dimension dim_j that are
	 * 				contained in the i*th (sub-)geom object of dimension dim_i
	 */
		size_t num(int dim_i, size_t i, int dim_j) const
			{return m_vSubNum[dim_i][i][dim_j];}

	/// id of object j in dimension dim_j of obj i in dimension dim_i
	/**
	 * A reference element is constructed by several geometric objects, that
	 * are mapped by a reference element by themselves. This method returns the
	 * id (w.r.t. this reference element) of a sub-geometric object that is
	 * part of a sub-geometric object of the reference element
	 *
	 * \param[in]	dim_i		dimension of sub geometric object
	 * \param[in]	i			id of sub geometric object
	 * \param[in]	dim_j		dimension for obj contained in the sub-object
	 * \param[in]	j			number of obj contained in the sub-object
	 * \returns		id of the j'th object of the dimension dim_j that are
	 * 				contained in the i*th (sub-)geom object of dimension dim_i
	 */
		int id(int dim_i, size_t i, int dim_j, size_t j) const
			{return m_id[dim_i][i][dim_j][j];}

	/// number of reference elements this element contains
		size_t num(ReferenceObjectID type) const
			{return m_vNumRefElem[type];}

	/// reference element type of obj nr i in dimension dim_i
		ReferenceObjectID roid(int dim_i, size_t i) const
			{return m_vRefElemType[dim_i][i];}

	/// print informations about the reference element
		void print_info() const;

	protected:
	/// to make it more readable
		enum{POINT = 0, EDGE = 1, FACE = 2, VOLUME= 3};

	///	maximum number of Objects in all dimensions
		enum{MAXOBJECTS = 12};

	///	maximum dimension
		enum{MAXDIM = 3};

	///	dimension of the reference world
		int m_dim;

	///	size of reference element
		number m_size;

	/// number of Geometric Objects of a dimension
		size_t m_vNum[MAXDIM+1];

	/// number of Geometric Objects contained in a (Sub-)Geometric Object of the Element
		size_t m_vSubNum[MAXDIM+1][MAXOBJECTS][MAXDIM+1];

	/// indices of GeomObjects
		int m_id[MAXDIM+1][MAXOBJECTS][MAXDIM+1][MAXOBJECTS];

	///	number of reference elements
		size_t m_vNumRefElem[NUM_REFERENCE_OBJECTS];

	///	type of reference elements
		ReferenceObjectID m_vRefElemType[MAXDIM+1][MAXOBJECTS];
};

/// dimension dependent base class for reference elements
/**
 * This is the base class for reference elements with their dimension. It
 * simply adds to the ReferenceElement base class the corner position of the
 * reference element vertices in local coordinates.
 *
 * \tparam 		d		dimension, where reference element lives
 */
template <int d>
class DimReferenceElement : public ReferenceElement
{
	public:
	///	dimension, where the reference element is defined
		static const int dim = d;

	/// coordinates of reference corner in a vector
		const MathVector<dim>* corners() const {return &m_vCorner[0];}

	/// coordinates of reference corner (i = 0 ... num(0))
		const MathVector<dim>& corner(size_t i) const {return m_vCorner[i];}

	/// coordinates of reference corner as integer
		const MathVector<dim,int>* corner() const {return m_vCoInt;}

	/// print informations about the reference element
		void print_info() const;

	protected:
	///	maximum number of corners for fixed reference elements
		enum{MAXCORNERS = 8};

	///	coordinates of Reference Corner
		MathVector<dim> m_vCorner[MAXCORNERS];
		MathVector<dim, int> m_vCoInt[MAXCORNERS];
};

///////////////////////////////////////////////////////////////////////////////
// Reference Element Providers
///////////////////////////////////////////////////////////////////////////////

/// Exception thrown when reference element not found
struct UGError_ReferenceElementMissing : public UGError
{
	UGError_ReferenceElementMissing(int dim_, ReferenceObjectID roid_)
	: UGError(""), dim(dim_), roid(roid_)
	{
		std::stringstream ss; ss << "Reference Element not found for "
							<<roid<<" (dim="<<dim<<")";
		UGError::push_msg(ss.str());
	}
	int dim;
	ReferenceObjectID roid;
};

/// Provider for Reference Elements
class ReferenceElementProvider
{
	private:
	///	constructor
		ReferenceElementProvider();

	//	intentionally left unimplemented
		ReferenceElementProvider(const ReferenceElementProvider&){};
		ReferenceElementProvider& operator=(const ReferenceElementProvider&);

	///	provide instance of singleton
		static ReferenceElementProvider& instance()
		{
			static ReferenceElementProvider inst;
			return inst;
		}

	///	adds a Reference Element
		static bool add_elem(const ReferenceElement& elem);

	///	returns a Reference Element
		static const ReferenceElement& get_elem(ReferenceObjectID roid);

	///	vector storing all ReferenceElement
		static const ReferenceElement* m_vElem[NUM_REFERENCE_OBJECTS];

	///	adds a Reference Element
		template <int dim>
		static bool add_dim_elem(const DimReferenceElement<dim>& elem);

	///	returns a Reference Element
		template <int dim>
		static const DimReferenceElement<dim>& get_dim_elem(ReferenceObjectID roid)
		{
			UG_ASSERT(roid >= 0, "roid ="<<roid<<" wrong")
			UG_ASSERT(roid < NUM_REFERENCE_OBJECTS, "roid ="<<roid<<" wrong")
			static const DimReferenceElement<dim>** vDimElem = get_vector<dim>();
			UG_ASSERT(vDimElem[roid] != NULL, "Null pointer for roid ="<<roid);
			return *vDimElem[roid];
		}

	///	returns vector of DimReferenceElement
		template <int dim>
		static const DimReferenceElement<dim>** get_vector()
		{
			static const DimReferenceElement<dim>* sVec[NUM_REFERENCE_OBJECTS];
			return sVec;
		}

	public:
	///	returns a dimension dependent Reference Element
		template <int dim>
		inline static const DimReferenceElement<dim>& get(ReferenceObjectID roid)
		{
			return instance().get_dim_elem<dim>(roid);
		}

	///	returns a Reference Element
		inline static const ReferenceElement& get(ReferenceObjectID roid)
		{
			return instance().get_elem(roid);
		}
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Concrete Reference Elements
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Reference Vertex
///////////////////////////////////////////////////////////////////////////////

class ReferenceVertex : public DimReferenceElement<1>
{
	public:
	///	type of reference element
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_VERTEX;

	///	dimension of reference element
		static const int dim = 0;

	///	number of corners
		static const int numCorners = 1;

	///	number of eges
		static const int numEdges = 0;

	///	number of faces
		static const int numFaces = 0;

	///	number of volumes
		static const int numVolumes = 0;

	///	number of sides
		static const int numSides = 0;

	public:
	///	Constructor filling the arrays
		ReferenceVertex();

	/// \copydoc ug::ReferenceElement::reference_object_id()
		ReferenceObjectID reference_object_id() const {return REFERENCE_OBJECT_ID;}

	/// \copydoc ug::ReferenceElement::dimension()
		int dimension() const {return dim;}

	/// \copydoc ug::ReferenceElement::size()
		number size() const	{return 1.0;}
};

///////////////////////////////////////////////////////////////////////////////
// Reference Edge
///////////////////////////////////////////////////////////////////////////////

class ReferenceEdge : public DimReferenceElement<1>
{
	public:
	///	type of reference element
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_EDGE;

	///	dimension of reference element
		static const int dim = 1;

	///	number of corners
		static const int numCorners = 2;

	///	number of eges
		static const int numEdges = 1;

	///	number of faces
		static const int numFaces = 0;

	///	number of volumes
		static const int numVolumes = 0;

	///	number of sides
		static const int numSides = numCorners;

	public:
	///	Constructor
		ReferenceEdge();

	/// \copydoc ug::ReferenceElement::reference_object_id()
		ReferenceObjectID reference_object_id() const {return REFERENCE_OBJECT_ID;}

	/// \copydoc ug::ReferenceElement::dimension()
		int dimension() const {return dim;}

	/// \copydoc ug::ReferenceElement::size()
		number size() const	{return 0.5;}

	///	\copydoc ug::DimReferenceElement::check_position()
		inline static void check_position(const MathVector<dim>& pos)
		{
			UG_ASSERT(pos[0] >= 0.0 && pos[0] <= 1.0,
			          "Local position "<<pos<<" outside Reference Element");
		}
};

///////////////////////////////////////////////////////////////////////////////
// Reference Triangle
///////////////////////////////////////////////////////////////////////////////

class ReferenceTriangle : public DimReferenceElement<2>
{
	public:
	///	type of reference element
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_TRIANGLE;

	///	dimension of reference element
		static const int dim = 2;

	///	number of corners
		static const int numCorners = 3;

	///	number of eges
		static const int numEdges = 3;

	///	number of faces
		static const int numFaces = 1;

	///	number of volumes
		static const int numVolumes = 0;

	///	number of sides
		static const int numSides = numEdges;

	public:
	///	Constructor filling the arrays
		ReferenceTriangle();

	/// \copydoc ug::ReferenceElement::reference_object_id()
		ReferenceObjectID reference_object_id() const {return REFERENCE_OBJECT_ID;}

	/// \copydoc ug::ReferenceElement::dimension()
		int dimension() const {return dim;}

	/// \copydoc ug::ReferenceElement::size()
		number size() const	{return 0.5;}

	///	\copydoc ug::DimReferenceElement::check_position()
		inline static void check_position(const MathVector<dim>& pos)
		{
			UG_ASSERT(pos[0] >= 0.0 && pos[0] <= 1.0 &&
			          pos[1] >= 0.0 && pos[1] <= 1.0 &&
			          pos[0]+pos[1] <= 1.0,
					  "Local position "<<pos<<" outside Reference Element");
		}
};

///////////////////////////////////////////////////////////////////////////////
// Reference Quadrilateral
///////////////////////////////////////////////////////////////////////////////

class ReferenceQuadrilateral : public DimReferenceElement<2>
{
	public:
	///	type of reference element
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_QUADRILATERAL;

	///	dimension of reference element
		static const int dim = 2;

	///	number of corners
		static const int numCorners = 4;

	///	number of eges
		static const int numEdges = 4;

	///	number of faces
		static const int numFaces = 1;

	///	number of volumes
		static const int numVolumes = 0;

	///	number of sides
		static const int numSides = numEdges;

	public:
		ReferenceQuadrilateral();

	/// \copydoc ug::ReferenceElement::reference_object_id()
		ReferenceObjectID reference_object_id() const {return REFERENCE_OBJECT_ID;}

	/// \copydoc ug::ReferenceElement::dimension()
		int dimension() const {return dim;}

	/// \copydoc ug::ReferenceElement::size()
		number size() const	{return 1.0;}

	///	\copydoc ug::DimReferenceElement::check_position()
		inline static void check_position(const MathVector<dim>& pos)
		{
			UG_ASSERT(pos[0] >= 0.0 && pos[0] <= 1.0 &&
					  pos[1] >= 0.0 && pos[1] <= 1.0,
					  "Local position "<<pos<<" outside Reference Element");
		}
};

///////////////////////////////////////////////////////////////////////////////
// Reference Tetrahedron
///////////////////////////////////////////////////////////////////////////////

class ReferenceTetrahedron : public DimReferenceElement<3>
{
	public:
	///	type of reference element
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_TETRAHEDRON;

	///	dimension of reference element
		static const int dim = 3;

	///	number of corners
		static const int numCorners = 4;

	///	number of eges
		static const int numEdges = 6;

	///	number of faces
		static const int numFaces = 4;

	///	number of volumes
		static const int numVolumes = 1;

	///	number of sides
		static const int numSides = numFaces;

	public:
	///	Constructor
		ReferenceTetrahedron();

	/// \copydoc ug::ReferenceElement::reference_object_id()
		ReferenceObjectID reference_object_id() const {return REFERENCE_OBJECT_ID;}

	/// \copydoc ug::ReferenceElement::dimension()
		int dimension() const {return dim;}

	/// \copydoc ug::ReferenceElement::size()
		number size() const	{return 1.0/6.0;}

	///	\copydoc ug::DimReferenceElement::check_position()
		inline static void check_position(const MathVector<dim>& pos)
		{
			UG_ASSERT(pos[0] >= 0.0 && pos[0] <= 1.0 &&
					  pos[1] >= 0.0 && pos[1] <= 1.0 &&
					  pos[2] >= 0.0 && pos[2] <= 1.0 &&
					  pos[0]+pos[1]+pos[2] <= 1.0,
					  "Local position "<<pos<<" outside Reference Element");
		}
};

///////////////////////////////////////////////////////////////////////////////
// Reference Pyramid
///////////////////////////////////////////////////////////////////////////////

class ReferencePyramid : public DimReferenceElement<3>
{
	public:
	///	type of reference element
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_PYRAMID;

	///	dimension of reference element
		static const int dim = 3;

	///	number of corners
		static const int numCorners = 5;

	///	number of eges
		static const int numEdges = 8;

	///	number of faces
		static const int numFaces = 5;

	///	number of volumes
		static const int numVolumes = 1;

	///	number of sides
		static const int numSides = numFaces;

	public:
	///	Constructor
		ReferencePyramid();

	/// \copydoc ug::ReferenceElement::reference_object_id()
		ReferenceObjectID reference_object_id() const {return REFERENCE_OBJECT_ID;}

	/// \copydoc ug::ReferenceElement::dimension()
		int dimension() const {return dim;}

	/// \copydoc ug::ReferenceElement::size()
		number size() const	{return 1.0/3.0;}

	///	\copydoc ug::DimReferenceElement::check_position()
		inline static void check_position(const MathVector<dim>& pos)
		{
			//\todo: add check
		}
};


///////////////////////////////////////////////////////////////////////////////
// Reference Prism
///////////////////////////////////////////////////////////////////////////////

class ReferencePrism : public DimReferenceElement<3>
{
	public:
	///	type of reference element
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_PRISM;

	///	dimension of reference element
		static const int dim = 3;

	///	number of corners
		static const int numCorners = 6;

	///	number of eges
		static const int numEdges = 9;

	///	number of faces
		static const int numFaces = 5;

	///	number of volumes
		static const int numVolumes = 1;

	///	number of sides
		static const int numSides = numFaces;

	public:
	///	Constructor
		ReferencePrism();

	/// \copydoc ug::ReferenceElement::reference_object_id()
		ReferenceObjectID reference_object_id() const {return REFERENCE_OBJECT_ID;}

	/// \copydoc ug::ReferenceElement::dimension()
		int dimension() const {return dim;}

	/// \copydoc ug::ReferenceElement::size()
		number size() const	{return 0.5;}

	///	\copydoc ug::DimReferenceElement::check_position()
		inline static void check_position(const MathVector<dim>& pos)
		{
			UG_ASSERT(pos[0] >= 0.0 && pos[0] <= 1.0 &&
					  pos[1] >= 0.0 && pos[1] <= 1.0 &&
					  pos[0]+pos[1] <= 1.0 &&
					  pos[2] >= 0.0 && pos[2] <= 1.0,
					  "Local position "<<pos<<" outside Reference Element");
		}
};

///////////////////////////////////////////////////////////////////////////////
// Reference Hexahedron
///////////////////////////////////////////////////////////////////////////////

///	reference element for a hexahedron
class ReferenceHexahedron : public DimReferenceElement<3>
{
	public:
	///	type of reference element
		static const ReferenceObjectID REFERENCE_OBJECT_ID = ROID_HEXAHEDRON;

	///	dimension of reference element
		static const int dim = 3;

	///	number of corners
		static const int numCorners = 8;

	///	number of eges
		static const int numEdges = 12;

	///	number of faces
		static const int numFaces = 6;

	///	number of volumes
		static const int numVolumes = 1;

	///	number of sides
		static const int numSides = numFaces;

	public:
	///	Constructor filling the arrays
		ReferenceHexahedron();

	/// \copydoc ug::ReferenceElement::reference_object_id()
		ReferenceObjectID reference_object_id() const {return REFERENCE_OBJECT_ID;}

	/// \copydoc ug::ReferenceElement::dimension()
		int dimension() const {return dim;}

	/// \copydoc ug::ReferenceElement::size()
		number size() const	{return 1.0;}

	///	\copydoc ug::DimReferenceElement::check_position()
		inline static void check_position(const MathVector<dim>& pos)
		{
			UG_ASSERT(pos[0] >= 0.0 && pos[0] <= 1.0 &&
					  pos[1] >= 0.0 && pos[1] <= 1.0 &&
					  pos[2] >= 0.0 && pos[2] <= 1.0,
					  "Local position "<<pos<<" outside Reference Element");
		}
};

} // end namespace ug

#include "reference_element_traits.h"

#endif /* __H__UG__LIB_DISC__REFERENCE_ELEMENT__REFERENCE_ELEMENT__ */
