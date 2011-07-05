/*
 * convection_diffusion.h
 *
 *  Created on: 02.08.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__CONVECTION_DIFFUSION__FE1__CONVECTION_DIFFUSION__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__CONVECTION_DIFFUSION__FE1__CONVECTION_DIFFUSION__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_discretization/spatial_discretization/elem_disc/elem_disc_interface.h"
#include "lib_discretization/spatial_discretization/ip_data/data_import_export.h"

namespace ug{


template<typename TDomain>
class FE1ConvectionDiffusionElemDisc
	: public IDomainElemDisc<TDomain>
{
	private:
	///	Base class type
		typedef IDomainElemDisc<TDomain> base_type;

	/// own type
		typedef FE1ConvectionDiffusionElemDisc<TDomain> this_type;

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

	public:
	///	Constructor
		FE1ConvectionDiffusionElemDisc();

	///	sets the diffusion tensor
	/**
	 * This method sets the Diffusion tensor used in computations. If no
	 * Tensor is set, a zero value is assumed.
	 */
		void set_diffusion(IPData<MathMatrix<dim, dim>, dim>& user) {m_imDiffusion.set_data(user);}

	///	sets the velocity field
	/**
	 * This method sets the Velocity field. If no field is provided a zero
	 * value is assumed.
	 */
		void set_velocity(IPData<MathVector<dim>, dim>& user) {m_imVelocity.set_data(user);}

	///	sets the reaction
	/**
	 * This method sets the Reaction. A zero value is assumed as default.
	 */
		void set_reaction(IPData<number, dim>& user) {m_imReaction.set_data(user);}

	///	sets the right-hand side
	/**
	 * This method sets the right hand side value. A zero value is assumed as
	 * default.
	 */
		void set_source(IPData<number, dim>& user)	{m_imSource.set_data(user);}


	///	sets mass scale
	/**
	 * This method sets the mass scale value. A value of 1.0 is assumed as
	 * default.
	 */
		void set_mass_scale(IPData<number, dim>& user)	{m_imMassScale.set_data(user);}

	public:
	///	number of functions used
		virtual size_t num_fct(){return 1;}

	///	type of trial space for each function used
		virtual LSFSID local_shape_function_set_id(size_t loc_fct)
		{
			return LSFSID(LSFSID::LAGRANGE, 1);
		}

	///	switches between non-regular and regular grids
		virtual bool treat_non_regular_grid(bool bNonRegular)
		{
		//	this disc does not need to take special care for non-regular grids
			return true;
		}

	///	does not use hanging at all for assembling
		virtual bool use_hanging() const {return false;}

	private:
	///	prepares the discretization for time dependent discretization
	/**
	 * This function prepares the discretization for time-dependent problems.
	 * It sets the time in the imports.
	 *
	 * \param[in]	time	new time point
	 * \returns 	true	indicates, that old values are needed
	 */
		virtual bool time_point_changed(number time);

		template <typename TElem>
		bool elem_loop_prepare_fe();

		template <typename TElem>
		bool elem_prepare_fe(TElem* elem, const local_vector_type& u);

		template <typename TElem>
		bool elem_loop_finish_fe();

		template <typename TElem>
		bool elem_JA_fe(local_matrix_type& J, const local_vector_type& u);

		template <typename TElem>
		bool elem_JM_fe(local_matrix_type& J, const local_vector_type& u);

		template <typename TElem>
		bool elem_dA_fe(local_vector_type& d, const local_vector_type& u);

		template <typename TElem>
		bool elem_dM_fe(local_vector_type& d, const local_vector_type& u);

		template <typename TElem>
		bool elem_rhs_fe(local_vector_type& d);

	private:
	///	Corner Coordinates
		const position_type* m_vCornerCoords;

	///	abbreviation for the local solution
		static const size_t _C_ = 0;

	///	Data import for Diffusion
		DataImport<MathMatrix<dim,dim>, dim> m_imDiffusion;

	///	Data import for the Velocity field
		DataImport<MathVector<dim>, dim > m_imVelocity;

	///	Data import for the reaction term
		DataImport<number, dim> m_imReaction;

	///	Data import for the right-hand side
		DataImport<number, dim> m_imSource;

	///	Data import for the right-hand side
		DataImport<number, dim> m_imMassScale;

	private:
		void register_all_fe1_funcs();

		struct RegisterFE1 {
				RegisterFE1(this_type* pThis) : m_pThis(pThis){}
				this_type* m_pThis;
				template< typename TElem > void operator()(TElem&)
				{m_pThis->register_fe1_func<TElem>();}
		};

		template <typename TElem>
		void register_fe1_func();

};

}

#endif /*__H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__CONVECTION_DIFFUSION__FE1__CONVECTION_DIFFUSION__*/
