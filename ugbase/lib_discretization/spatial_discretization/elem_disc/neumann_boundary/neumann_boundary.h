/*
 * neumann_boundary.h
 *
 *  Created on: 26.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__NEUMANN_BOUNDARY__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__NEUMANN_BOUNDARY__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_discretization/spatial_discretization/elem_disc/elem_disc_interface.h"

#include "lib_discretization/spatial_discretization/ip_data/const_user_data.h"

namespace ug{

template<typename TDomain>
class FV1NeumannBoundaryElemDisc
	: public IDomainElemDisc<TDomain>
{
	private:
	///	Base class type
		typedef IDomainElemDisc<TDomain> base_type;

	///	Base class type
		typedef FV1NeumannBoundaryElemDisc<TDomain> this_type;

	///	explicitly forward function
		using base_type::time;

	public:
	///	Domain type
		typedef typename base_type::domain_type domain_type;

	///	World dimension
		static const int dim = base_type::dim;

	///	Position type
		typedef typename base_type::position_type position_type;

	///	Local matrix type
		typedef typename base_type::local_matrix_type local_matrix_type;

	///	Local vector type
		typedef typename base_type::local_vector_type local_vector_type;

	///	Local index type
		typedef typename base_type::local_index_type local_index_type;

	protected:
	///	type of bnd number
		typedef typename IBoundaryData<number, dim>::functor_type BNDNumberFunctor;

	public:
	///	default constructor
		FV1NeumannBoundaryElemDisc();

	///	add a boundary value
		bool add_boundary_value(IBoundaryData<number, dim>& user, const char* function, const char* subsets);

	///	add a boundary value
		bool add_boundary_value(IBoundaryData<number, dim>& user, size_t fct, SubsetGroup bndSubsetGroup);

	private:
	//	Functor, function grouping
		struct UserDataFunction
		{
			UserDataFunction(size_t fct_, BNDNumberFunctor functor_)
				: loc_fct(fct_), functor(functor_) {}

			size_t loc_fct;
			BNDNumberFunctor functor;
		};

		std::map<int, std::vector<UserDataFunction> > m_mBoundarySegment;
		size_t m_numFct;

	public:
	///	number of functons required
		virtual size_t num_fct(){return m_numFct;}

	///	type of function required
		virtual LocalShapeFunctionSetID local_shape_function_set_id(size_t loc_fct)
		{
			return LocalShapeFunctionSetID(LocalShapeFunctionSetID::LAGRANGE, 1);
		}

	///	switches between non-regular and regular grids
		virtual bool treat_non_regular_grid(bool bNonRegular)
		{
		//	switch, which assemble functions to use.
			if(bNonRegular)
			{
				UG_LOG("ERROR in 'FVNeumannBoundaryElemDisc::treat_non_regular_grid':"
						" Non-regular grid not implemented.\n");
				return false;
			}

		//	this disc supports regular grids
			return true;
		}


	private:
		template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
		bool prepare_element_loop();

		template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
		bool prepare_element(TElem* elem, const local_vector_type& u, const local_index_type& glob_ind);

		template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
		bool finish_element_loop();

		template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
		bool assemble_JA(local_matrix_type& J, const local_vector_type& u);

		template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
		bool assemble_JM(local_matrix_type& J, const local_vector_type& u);

		template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
		bool assemble_A(local_vector_type& d, const local_vector_type& u);

		template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
		bool assemble_M(local_vector_type& d, const local_vector_type& u);

		template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
		bool assemble_f(local_vector_type& d);

	private:
	// 	position access
		const position_type* m_vCornerCoords;

	private:
		void register_all_fv1_funcs(bool bHang);

		template <template <class Elem, int WorldDim> class TFVGeom>
		struct RegisterFV1 {
				RegisterFV1(this_type* pThis) : m_pThis(pThis){}
				this_type* m_pThis;
				template< typename TElem > void operator()(TElem&)
				{m_pThis->register_fv1_func<TElem, TFVGeom>();}
		};

		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		void register_fv1_func();

};

} // end namespac ug

#endif /*__H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__NEUMANN_BOUNDARY__*/
