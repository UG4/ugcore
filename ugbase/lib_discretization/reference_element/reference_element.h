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
 */
class ReferenceElement
{
	public:
	/// returns the reference object id
		virtual ReferenceObjectID reference_object_id() const = 0;

	/// returns the dimension where reference element lives
		virtual int dimension() const = 0;

	/// returns the size of reference element
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
		virtual size_t num_obj(int dim) const = 0;

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
		virtual size_t num_obj_of_obj(int dim_i, size_t i, int dim_j) const = 0;

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
	/// coordinates of reference corner (i = 0 ... num_obj(0))
		virtual const MathVector<dim>& corner(int i) const = 0;

	/// print informations about the reference element
		virtual void print_info() const;
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

	/// \copydoc ug::ReferenceElement::num_obj()
		size_t num_obj(int dim) const {return TRefElem::num_obj(dim);}

	/// \copydoc ug::ReferenceElement::num_obj_of_obj()
		size_t num_obj_of_obj(int dim_i, size_t i, int dim_j) const
			{return TRefElem::num_obj_of_obj(dim_i, i, dim_j);}

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
template <typename TRefElem, int d>
class DimReferenceElementWrapper
	: public DimReferenceElement<d>, protected TRefElem
{
	public:
	///	\copydoc ug::DimReferenceElement<d>::dim
		static const int dim = d;

	public:
	///	\copydoc ug::DimReferenceElement<d>::reference_object_id()
		ReferenceObjectID reference_object_id() const
			{return TRefElem::reference_object_id();}

	///	\copydoc ug::DimReferenceElement<d>::dimension()
		int dimension() const {return TRefElem::dimension();}

	///	\copydoc ug::DimReferenceElement<d>::size()
		number size() const {return TRefElem::size();}

	///	\copydoc ug::DimReferenceElement<d>::num_obj()
		size_t num_obj(int dim) const {return TRefElem::num_obj(dim);}

	///	\copydoc ug::DimReferenceElement<d>::num_obj_of_obj()
		size_t num_obj_of_obj(int dim_i, size_t i, int dim_j) const
			{return TRefElem::num_obj_of_obj(dim_i, i, dim_j);}

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
		const MathVector<dim>& corner(int i) const {return TRefElem::corner(i);}
};


///////////////////////////////////////////////////////////////////////////////
// Reference Element Factories
///////////////////////////////////////////////////////////////////////////////

// predeclaration
bool RegisterStandardDimReferenceElements();

template <int d>
class DimReferenceElementFactory{
	private:
		DimReferenceElementFactory(){init();};
		DimReferenceElementFactory(const DimReferenceElementFactory&){};
		DimReferenceElementFactory& operator=(const DimReferenceElementFactory&);

		static const DimReferenceElement<d>& get_elem(int refID)
		{
			UG_ASSERT((size_t)refID < m_vElem.size(),
			          "ReferenceObjectID does not exist.");
			UG_ASSERT(m_vElem[refID] != 0, "Object does not exist.");

			return *m_vElem[refID];
		}

		inline static DimReferenceElementFactory<d>& instance()
		{
			static DimReferenceElementFactory<d> inst;
			return inst;
		}

		bool init()
		{
			static bool isInit = false;
			if(!isInit)
			{
				return RegisterStandardDimReferenceElements();
			}
			else return true;
		}


		static std::vector<const DimReferenceElement<d>*> m_vElem;

	public:
		static bool register_reference_element(const DimReferenceElement<d>& elem)
		{
			const int refID = elem.reference_object_id();
			if((int) m_vElem.size() <= refID) m_vElem.resize(refID+1, 0);
			if(m_vElem[refID] != 0)
			{
				UG_LOG("Reference element for this ID already registered.\n");
				return false;
			}

			m_vElem[refID] = &elem;
			return true;
		}

		inline static const DimReferenceElement<d>& get_reference_element(int refID)
		{
			return instance().get_elem(refID);
		}
};

template <int d>
std::vector<const DimReferenceElement<d>* > DimReferenceElementFactory<d>::m_vElem =
	std::vector<const DimReferenceElement<d>* >();


class ReferenceElementFactory{
	private:
		ReferenceElementFactory(){init();};
		ReferenceElementFactory(const ReferenceElementFactory&){};
		ReferenceElementFactory& operator=(const ReferenceElementFactory&);

		static const ReferenceElement& get_elem(int refID)
		{
			UG_ASSERT((size_t)refID < m_vElem.size(),
			          "ReferenceObjectID does not exist.");
			UG_ASSERT(m_vElem[refID] != 0, "Object does not exist.");

			return *m_vElem[refID];
		}

		static ReferenceElementFactory& instance()
		{
			static ReferenceElementFactory inst;
			return inst;
		}

		bool init()
		{
			static bool isInit = false;
			if(!isInit)
			{
				return RegisterStandardDimReferenceElements();
			}
			else return true;
		}

		static std::vector<const ReferenceElement*> m_vElem;

	public:
		static bool register_reference_element(const ReferenceElement& elem)
		{
			const int refID = elem.reference_object_id();
			if((int) m_vElem.size() <= refID) m_vElem.resize(refID+1, 0);
			if(m_vElem[refID] != 0)
			{
				UG_LOG("Reference element for this ID already registered.\n");
				return false;
			}

			m_vElem[refID] = &elem;
			return true;
		}

		inline static const ReferenceElement& get_reference_element(int refID)
		{
			return instance().get_elem(refID);
		}
};

// Singleton, holding a Reference Elements
class ReferenceElementProvider {

	// 	private constructor
		ReferenceElementProvider();

	// 	disallow copy and assignment (intentionally left unimplemented)
		ReferenceElementProvider(const ReferenceElementProvider&);
		ReferenceElementProvider& operator=(const ReferenceElementProvider&);

	// 	private destructor
		~ReferenceElementProvider(){};

	// 	geometry provider
		template <typename TRefElem>
		inline static const TRefElem& inst()
		{
			static TRefElem myInst;
			return myInst;
		};

	public:
	//	get the reference element by reference element type
		template <typename TRefElem>
		inline static const TRefElem& get()
		{
			return inst<TRefElem>();
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
