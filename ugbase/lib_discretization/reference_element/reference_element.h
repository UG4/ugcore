/*
 * reference_element.h
 *
 *  Created on: 13.05.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_ELEMENT__
#define __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_ELEMENT__

#include <cassert>
#include <iostream>
#include <sstream>
#include "common/common.h"
#include "common/math/ugmath.h"
#include "lib_grid/lg_base.h"

namespace ug{

///////////////////////////////////////////////////////////////////////////////
// Reference Element
///////////////////////////////////////////////////////////////////////////////

/// virtual reference element interface
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
 */
class ReferenceElement
{
	public:
	/// returns the reference object id
		virtual ReferenceObjectID reference_object_id() const = 0;

	/// returns the dimension where reference element lives
		virtual int dimension() const = 0;

	/// returns the size (e.g. area or volume) of the reference element
		virtual number size() const = 0;

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
		virtual size_t num(int dim) const = 0;

	/// returns the number of object of dim for a sub-geometric object
	/**
	 * A reference element is constructed by several geometric objects, that
	 * are mapped by a reference element by themselves. This method returns how
	 * many (sub-)geometric objects of a given dimension are contained in a
	 * (sub-)geometric object of this reference element.
	 *
	 * \param[in]	dim_i		dimension of sub geometric object
	 * \param[in]	i			number of sub geomertic object
	 * \param[in]	dim_j		dimension for elems contained in the sub-object
	 * \returns		number of objects of the dimension dim_j that are
	 * 				contained in the i*th (sub-)geom object of dimension dim_i
	 */
		virtual size_t num(int dim_i, size_t i, int dim_j) const = 0;

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
		virtual int id(int dim_i, size_t i, int dim_j, size_t j) const = 0;

	/// number of reference elements this element contains
		virtual size_t num_ref_elem(ReferenceObjectID type) const = 0;

	/// reference element type of obj nr i in dimension dim_i
		virtual ReferenceObjectID ref_elem_type(int dim_i, size_t i) const = 0;

	/// virtual destructor
		virtual ~ReferenceElement() {}

	/// print informations about the reference element
		virtual void print_info() const;
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

	public:
	/// coordinates of reference corner (i = 0 ... num(0))
		virtual const MathVector<dim>& corner(size_t i) const = 0;

	/// coordinates of reference corner as integer
		virtual const MathVector<dim,int>* corner() const = 0;

	/// print informations about the reference element
		virtual void print_info() const;

	///	checks in debug mode if a local position is inside of the reference element
		virtual void check_position(const MathVector<dim>& pos) const = 0;
};

/// wrapper class for reference elements
/**
 * This class wraps a reference element into the ReferenceElement-interface to
 * make it available through the interface on the price of virtual functions.
 *
 * \tparam		TRefElem		the reference element to wrap
 */
template <typename TRefElem>
class ReferenceElementWrapper
	: public ReferenceElement, protected TRefElem
{
	public:
	/// \copydoc ug::ReferenceElement::reference_object_id()
		ReferenceObjectID reference_object_id() const
			{return TRefElem::reference_object_id();}

	/// \copydoc ug::ReferenceElement::dimension()
		int dimension() const {return TRefElem::dimension();}

	/// \copydoc ug::ReferenceElement::size()
		number size() const {return TRefElem::size();}

	/// \copydoc ug::ReferenceElement::num(int)
		size_t num(int dim) const {return TRefElem::num(dim);}

	/// \copydoc ug::ReferenceElement::num(int, size_t, int)
		size_t num(int dim_i, size_t i, int dim_j) const
			{return TRefElem::num(dim_i, i, dim_j);}

	/// \copydoc ug::ReferenceElement::id()
		int id(int dim_i, size_t i, int dim_j, size_t j) const
			{return TRefElem::id(dim_i, i, dim_j, j);}

	/// \copydoc ug::ReferenceElement::num_ref_elem()
		size_t num_ref_elem(ReferenceObjectID type) const
			{return TRefElem::num_ref_elem(type);}

	/// \copydoc ug::ReferenceElement::ref_elem_type()
		ReferenceObjectID ref_elem_type(int dim_i, size_t i) const
			{return TRefElem::ref_elem_type(dim_i, i);}
};

/// wrapper class for dimension dependent reference elements
/**
 * This class wraps a reference element into the DimReferenceElement-interface to
 * make it available through the interface on the price of virtual functions.
 *
 * \tparam		TRefElem		the reference element to wrap
 */
template <typename TRefElem>
class DimReferenceElementWrapper
	: public DimReferenceElement<TRefElem::dim>, protected TRefElem
{
	public:
	///	\copydoc ug::DimReferenceElement<d>::dim
		static const int dim = TRefElem::dim;

	public:
	///	\copydoc ug::DimReferenceElement<d>::reference_object_id()
		ReferenceObjectID reference_object_id() const
			{return TRefElem::reference_object_id();}

	///	\copydoc ug::DimReferenceElement<d>::dimension()
		int dimension() const {return TRefElem::dimension();}

	///	\copydoc ug::DimReferenceElement<d>::size()
		number size() const {return TRefElem::size();}

	///	\copydoc ug::DimReferenceElement<d>::num()
		size_t num(int dim) const {return TRefElem::num(dim);}

	///	\copydoc ug::DimReferenceElement<d>::num()
		size_t num(int dim_i, size_t i, int dim_j) const
			{return TRefElem::num(dim_i, i, dim_j);}

	///	\copydoc ug::DimReferenceElement<d>::id()
		int id(int dim_i, size_t i, int dim_j, size_t j) const
			{return TRefElem::id(dim_i, i, dim_j, j);}

	///	\copydoc ug::DimReferenceElement<d>::num_ref_elem()
		size_t num_ref_elem(ReferenceObjectID type) const
			{return TRefElem::num_ref_elem(type);}

	///	\copydoc ug::DimReferenceElement<d>::ref_elem_type()
		ReferenceObjectID ref_elem_type(int dim_i, size_t i) const
			{return TRefElem::ref_elem_type(dim_i, i);}

	///	\copydoc ug::DimReferenceElement<d>::corner()
		const MathVector<dim>& corner(size_t i) const {return TRefElem::corner(i);}

	/// \copydoc ug::DimReferenceElement::corner()
		const MathVector<dim,int>* corner() const {return TRefElem::corner();}

	///	\copydoc ug::DimReferenceElement::check_position()
		void check_position(const MathVector<dim>& pos) const
			{TRefElem::check_position(pos);}
};


///////////////////////////////////////////////////////////////////////////////
// Reference Element Providers
///////////////////////////////////////////////////////////////////////////////

/// Exception thrown when reference element not found
struct UG_ERROR_ReferenceElementMissing
	: public UGFatalError
{
	UG_ERROR_ReferenceElementMissing(int dim_, ReferenceObjectID roid_)
	: UGFatalError(""), dim(dim_), roid(roid_)
	{
		std::stringstream ss; ss << "Refernce Element not found for "
							<<roid<<" (dim="<<dim<<")";
		UGFatalError::set_msg(ss.str());
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

} // end namespace ug

#include "reference_element_mapping.h"

#include "reference_vertex.h"
#include "reference_edge.h"
#include "reference_triangle.h"
#include "reference_quadrilateral.h"
#include "reference_tetrahedron.h"
#include "reference_pyramid.h"
#include "reference_prism.h"
#include "reference_hexahedron.h"

#include "reference_element_traits.h"

#endif /* __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_ELEMENT__ */
