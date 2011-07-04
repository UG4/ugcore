/*
 * neumann_boundary.h
 *
 *  Created on: 26.02.2010
 *      Author: markusbreit
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__INNER_BOUNDARY__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__INNER_BOUNDARY__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_discretization/spatial_discretization/elem_disc/elem_disc_interface.h"

#include "lib_discretization/spatial_discretization/ip_data/const_user_data.h"

namespace ug{

/// Finite Volume Element Discretization for an inner BndCond that depends on the unknowns (on the bnd)
/**
 * This class implements the IElemDisc interface to provide element local
 * assemblings for the unknown-dependent Neumann-flux over an inner boundary.
 * The equation of this flux should be given on the script level.
 * 
 * \tparam	TDomain		Domain
 * \tparam	TAlgebra	Algebra
 */


template<typename TDomain>
class FVInnerBoundaryElemDisc
: public IDomainElemDisc<TDomain>
{
	private:
	///	Base class type
		typedef IDomainElemDisc<TDomain> base_type;

	///	own type
		typedef FVInnerBoundaryElemDisc<TDomain> this_type;

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
		typedef typename IBoundaryData<number, dim>::functor_type BNDNumberFunctor;

	public:
		FVInnerBoundaryElemDisc()
		{
			register_all_fv1_funcs();
		}
		
	private:
	//	number of functions required on manifold
		// \todo: handle flexible, using flexible DataLinker UserFunction
		static const size_t m_numFct = 3;
	
	public:	// inherited from IElemDisc
	///	number of functions used
		virtual size_t num_fct(){return m_numFct;}
		
	///	type of trial space for each function used
		virtual LSFSID local_shape_function_set_id(size_t loc_fct)
		{
			return LSFSID(LSFSID::LAGRANGE, 1);
		}

	///	switches between non-regular and regular grids
		virtual bool treat_non_regular_grid(bool bNonRegular)
		{
		//	switch, which assemble functions to use.
			if(bNonRegular)
			{
				UG_LOG("ERROR in 'DensityDrivenFlowElemDisc::treat_non_regular_grid':"
						" Non-regular grid not implemented.\n");
				return false;
			}

		//	this disc supports regular grids
			return true;
		}

	///	returns if hanging nodes are used
		virtual bool use_hanging() const {return false;}

	private:
	
	///	prepares the loop over all elements
	/**
	 * This method prepares the loop over all elements. It resizes the Position
	 * array for the corner coordinates and schedules the local ip positions
	 * at the data imports.
	 */
		template<typename TElem, template <class Elem, int Dim> class TFVGeom>
		bool prepare_element_loop();

	///	prepares the element for assembling
	/**
	 * This methods prepares an element for the assembling. The Positions of
	 * the Element Corners are read and the Finite Volume Geometry is updated.
	 * The global ip positions are scheduled at the data imports.
	 */
		template<typename TElem, template <class Elem, int Dim> class TFVGeom>
		bool prepare_element(TElem* elem, const local_vector_type& u, const local_index_type& glob_ind);

	///	finishes the loop over all elements
		template<typename TElem, template <class Elem, int Dim> class TFVGeom>
		bool finish_element_loop();

	///	assembles the local stiffness matrix using a finite volume scheme
		template<typename TElem, template <class Elem, int Dim> class TFVGeom>
		bool assemble_JA(local_matrix_type& J, const local_vector_type& u);

	///	assembles the local mass matrix using a finite volume scheme
		template<typename TElem, template <class Elem, int Dim> class TFVGeom>
		bool assemble_JM(local_matrix_type& J, const local_vector_type& u);

	///	assembles the stiffness part of the local defect
		template<typename TElem, template <class Elem, int Dim> class TFVGeom>
		bool assemble_A(local_vector_type& d, const local_vector_type& u);

	///	assembles the mass part of the local defect
		template<typename TElem, template <class Elem, int Dim> class TFVGeom>
		bool assemble_M(local_vector_type& d, const local_vector_type& u);

	///	assembles the local right hand side
		template<typename TElem, template <class Elem, int Dim> class TFVGeom>
		bool assemble_f(local_vector_type& d);

	private:
		// position access
		const position_type* m_vCornerCoords;

	private:
		void register_all_fv1_funcs();

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

}

#endif /*__H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__INNER_BOUNDARY__*/
