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


class ReferenceElement {
	public:
		/// reference object id
		virtual ReferenceObjectID reference_object_id() const = 0;

		/// Dimension where reference element lives
		virtual int dimension() const = 0;

		/// size of reference triangle
		virtual number size() const = 0;

		/// number of objects of dim
		virtual size_t num_obj(int dim) const = 0;

		/// number of object of dim
		virtual size_t num_obj_of_obj(int dim_i, size_t i, int dim_j) const = 0;

		/// id of object j in dimension dim_j of obj i in dimension dim_i
		virtual int id(int dim_i, size_t i, int dim_j, size_t j) const = 0;

		/// number of reference elements this element is contained of
		virtual size_t num_ref_elem(ReferenceObjectID type) const = 0;

		/// reference element type of obj nr i in dimension dim_i */
		virtual ReferenceObjectID ref_elem_type(int dim_i, size_t i) const = 0;

		/// virtual destructor
		virtual ~ReferenceElement()
		{}

		/// print informations about the reference element
		virtual void print_info() const
		{
			using namespace std;

			string GeomObjects[4] ={"Corner", "Edge", "Face", "Volume"};

			cout << "Reference Element Info: " << endl;
			cout << "----------------------- " << endl;

			cout << "Size: " << size() << endl;
			cout << "Dimension where Reference Element lives: " << dimension() << endl;

			for(int i = dimension(); i>=0 ;i--)
				cout << "Number of " << GeomObjects[i] << "s: " << num_obj(i) << endl;

			for(int dim_i = dimension(); dim_i>=0 ;dim_i--)
			{
				for(size_t i=0; i < num_obj(dim_i); i++)
				{
					cout << GeomObjects[dim_i] << " with id '" << i << "' contains the following GeomObjects:" << endl;
					for(int dim_j=dim_i; dim_j>= 0; dim_j--)
					{
						cout << num_obj_of_obj(dim_i,i,dim_j) << " " << GeomObjects[dim_j] << "s with id: ";
						for(size_t j=0; j< num_obj_of_obj(dim_i,i,dim_j); j++)
						{
							cout << id(dim_i,i,dim_j,j) << " ";
						}
						cout << endl;
					}
				}
			}
		}
};

template <int d>
class DimReferenceElement : public ReferenceElement
{
	public:
		static const int dim = d;

	public:
		/// coordinates of reference corner (i = 0 ... num_obj(0))
		virtual const MathVector<dim>& corner(int i) const = 0;

		/// print informations about the reference element
		virtual void print_info() const
		{
			using namespace std;

			ReferenceElement::print_info();

			cout << "corners:\n";
			for(size_t i = 0; i< num_obj(0); i++)
			{
				cout << i << ":" << corner(i) << "\n";
			}
		}
};


template <typename TRefElem>
class ReferenceElementWrapper : public ReferenceElement, protected TRefElem
{
	public:
		/// reference object id
		ReferenceObjectID reference_object_id() const {return TRefElem::reference_object_id();}

		/// Dimension where reference element lives
		int dimension() const {return TRefElem::dimension();}

		/// size of reference triangle
		number size() const {return TRefElem::size();}

		/// number of objects of dim
		size_t num_obj(int dim) const {return TRefElem::num_obj(dim);}

		/// number of object of dim
		size_t num_obj_of_obj(int dim_i, size_t i, int dim_j) const {return TRefElem::num_obj_of_obj(dim_i, i, dim_j);}

		/// id of object j in dimension dim_j of obj i in dimension dim_i
		int id(int dim_i, size_t i, int dim_j, size_t j) const {return TRefElem::id(dim_i, i, dim_j, j);}

		/// number of reference elements this element is contained of
		size_t num_ref_elem(ReferenceObjectID type) const {return TRefElem::num_ref_elem(type);}

		/// reference element type of obj nr i in dimension dim_i */
		ReferenceObjectID ref_elem_type(int dim_i, size_t i) const {return TRefElem::ref_elem_type(dim_i, i);}
};

template <typename TRefElem, int d>
class DimReferenceElementWrapper : public DimReferenceElement<d>, protected TRefElem
{
		//UG_STATIC_ASSERT(TRefElem::dim == d, Dimension_does_not_match);
	public:
		static const int dim = d;

	public:
		/// reference object id
		ReferenceObjectID reference_object_id() const {return TRefElem::reference_object_id();}

		/// Dimension where reference element lives
		int dimension() const {return TRefElem::dimension();}

		/// size of reference triangle
		number size() const {return TRefElem::size();}

		/// number of objects of dim
		size_t num_obj(int dim) const {return TRefElem::num_obj(dim);}

		/// number of object of dim
		size_t num_obj_of_obj(int dim_i, size_t i, int dim_j) const {return TRefElem::num_obj_of_obj(dim_i, i, dim_j);}

		/// id of object j in dimension dim_j of obj i in dimension dim_i
		int id(int dim_i, size_t i, int dim_j, size_t j) const {return TRefElem::id(dim_i, i, dim_j, j);}

		/// number of reference elements this element is contained of
		size_t num_ref_elem(ReferenceObjectID type) const {return TRefElem::num_ref_elem(type);}

		/// reference element type of obj nr i in dimension dim_i */
		ReferenceObjectID ref_elem_type(int dim_i, size_t i) const {return TRefElem::ref_elem_type(dim_i, i);}

		/// coordinates of reference corner (i = 0 ... num_obj(0))
		const MathVector<dim>& corner(int i) const {return TRefElem::corner(i);}
};

template <int d>
class DimReferenceElementFactory{
	private:
		DimReferenceElementFactory(){};
		DimReferenceElementFactory(const DimReferenceElementFactory&){};
		DimReferenceElementFactory& operator=(const DimReferenceElementFactory&);

		static const DimReferenceElement<d>& get_elem(int refID)
		{
			UG_ASSERT((size_t)refID < m_vElem.size(), "ReferenceObjectID does not exist.");
			UG_ASSERT(m_vElem[refID] != 0, "Object does not exist.");

			return *m_vElem[refID];
		}

		static DimReferenceElementFactory<d>& instance()
		{
			static DimReferenceElementFactory<d> inst;
			return inst;
		}


		static std::vector<const DimReferenceElement<d>*> m_vElem;

	public:
		static bool register_reference_element(const DimReferenceElement<d>& elem)
		{
			const int refID = elem.reference_object_id();
			if((int) m_vElem.size() <= refID) m_vElem.resize(refID+1, 0);
			if(m_vElem[refID] != 0)
				{UG_LOG("Reference element for this ID already registered.\n"); return false;}

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
		ReferenceElementFactory(){};
		ReferenceElementFactory(const ReferenceElementFactory&){};
		ReferenceElementFactory& operator=(const ReferenceElementFactory&);

		static const ReferenceElement& get_elem(int refID)
		{
			UG_ASSERT((size_t)refID < m_vElem.size(), "ReferenceObjectID does not exist.");
			UG_ASSERT(m_vElem[refID] != 0, "Object does not exist.");

			return *m_vElem[refID];
		}

		static ReferenceElementFactory& instance()
		{
			static ReferenceElementFactory inst;
			return inst;
		}


		static std::vector<const ReferenceElement*> m_vElem;

	public:
		static bool register_reference_element(const ReferenceElement& elem)
		{
			const int refID = elem.reference_object_id();
			if((int) m_vElem.size() <= refID) m_vElem.resize(refID+1, 0);
			if(m_vElem[refID] != 0)
				{UG_LOG("Reference element for this ID already registered.\n"); return false;}

			m_vElem[refID] = &elem;
			return true;
		}

		inline static const ReferenceElement& get_reference_element(int refID)
		{
			return instance().get_elem(refID);
		}
};


template <class TElem>
class reference_element_traits
{};



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

#endif /* __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_ELEMENT__ */
