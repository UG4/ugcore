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
#include "lib_discretization/spatial_discretization/disc_helper/fvgeom.h"
#include "lib_discretization/spatial_discretization/elem_disc/elem_disc_interface.h"
#include "lib_discretization/common/local_algebra.h"
#include "lib_discretization/spatial_discretization/ip_data/ip_data.h"
#include "lib_discretization/spatial_discretization/ip_data/data_export.h"
#include "lib_discretization/spatial_discretization/ip_data/data_import.h"

namespace ug{


template<typename TDomain, typename TAlgebra>
class FE1ConvectionDiffusionElemDisc
: public IDomainElemDisc<TDomain, TAlgebra>
{
	private:
	///	Base class type
		typedef IDomainElemDisc<TDomain, TAlgebra> base_type;

	public:
	///	Domain type
		typedef typename base_type::domain_type domain_type;

	///	World dimension
		static const int dim = base_type::dim;

	///	Position type
		typedef typename base_type::position_type position_type;

	///	Algebra type
		typedef typename base_type::algebra_type algebra_type;

	///	Local matrix type
		typedef typename base_type::local_matrix_type local_matrix_type;

	///	Local vector type
		typedef typename base_type::local_vector_type local_vector_type;

	///	Local index type
		typedef typename base_type::local_index_type local_index_type;

	public:
	///	Constructor
		FE1ConvectionDiffusionElemDisc()
		{
		//	register assemling functions
			register_ass_funcs(Int2Type<dim>());

		//	register imports
			register_import(m_Diff);
			register_import(m_ConvVel);
			register_import(m_Reaction);
			register_import(m_Rhs);
			register_import(m_MassScale);
		}

	///	sets the diffusion tensor
	/**
	 * This method sets the Diffusion tensor used in computations. If no
	 * Tensor is set, a zero value is assumed.
	 */
		void set_diffusion(IPData<MathMatrix<dim, dim>, dim>& user) {m_Diff.set_data(user);}

	///	sets the velocity field
	/**
	 * This method sets the Velocity field. If no field is provided a zero
	 * value is assumed.
	 */
		void set_velocity(IPData<MathVector<dim>, dim>& user) {m_ConvVel.set_data(user);}

	///	sets the reaction
	/**
	 * This method sets the Reaction. A zero value is assumed as default.
	 */
		void set_reaction(IPData<number, dim>& user) {m_Reaction.set_data(user);}

	///	sets the right-hand side
	/**
	 * This method sets the right hand side value. A zero value is assumed as
	 * default.
	 */
		void set_source(IPData<number, dim>& user)	{m_Rhs.set_data(user);}


	///	sets mass scale
	/**
	 * This method sets the mass scale value. A value of 1.0 is assumed as
	 * default.
	 */
		void set_mass_scale(IPData<number, dim>& user)	{m_MassScale.set_data(user);}

	public:
	///	number of functions used
		virtual size_t num_fct(){return 1;}

	///	type of trial space for each function used
		virtual LocalShapeFunctionSetID local_shape_function_set_id(size_t loc_fct)
		{
			return LocalShapeFunctionSetID(LocalShapeFunctionSetID::LAGRANGE, 1);
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
		template <typename TElem>
		bool prepare_element_loop();

		template <typename TElem>
		bool prepare_element(TElem* elem, const local_vector_type& u, const local_index_type& glob_ind);

		template <typename TElem>
		bool finish_element_loop();

		template <typename TElem>
		bool assemble_JA(local_matrix_type& J, const local_vector_type& u);

		template <typename TElem>
		bool assemble_JM(local_matrix_type& J, const local_vector_type& u);

		template <typename TElem>
		bool assemble_A(local_vector_type& d, const local_vector_type& u);

		template <typename TElem>
		bool assemble_M(local_vector_type& d, const local_vector_type& u);

		template <typename TElem>
		bool assemble_f(local_vector_type& d);

	private:
	///	Corner Coordinates
		const position_type* m_vCornerCoords;

	///	abbreviation for the local solution
		static const size_t _C_ = 0;

	///	Data import for Diffusion
		DataImport<MathMatrix<dim,dim>, dim, algebra_type> m_Diff;

	///	Data import for the Velocity field
		DataImport<MathVector<dim>, dim, algebra_type > m_ConvVel;

	///	Data import for the reaction term
		DataImport<number, dim, algebra_type> m_Reaction;

	///	Data import for the right-hand side
		DataImport<number, dim, algebra_type> m_Rhs;

	///	Data import for the right-hand side
		DataImport<number, dim, algebra_type> m_MassScale;

	private:
		///////////////////////////////////////
		// registering for reference elements
		///////////////////////////////////////

	/// register for 1D
		void register_ass_funcs(Int2Type<1>)
		{
			register_all_ass_funcs<Edge>(ROID_EDGE);
		}

	/// register for 2D
		void register_ass_funcs(Int2Type<2>)
		{
			register_ass_funcs(Int2Type<1>());
			register_all_ass_funcs<Triangle>(ROID_TRIANGLE);
			register_all_ass_funcs<Quadrilateral>(ROID_QUADRILATERAL);
		}

	/// register for 3D
		void register_ass_funcs(Int2Type<3>)
		{
			register_ass_funcs(Int2Type<2>());
			register_all_ass_funcs<Tetrahedron>(ROID_TETRAHEDRON);
			register_all_ass_funcs<Pyramid>(ROID_PYRAMID);
			register_all_ass_funcs<Prism>(ROID_PRISM);
			register_all_ass_funcs<Hexahedron>(ROID_HEXAHEDRON);
		}

	///	register all functions for on element type
		template <typename TElem>
		void register_all_ass_funcs(int id)
		{
			typedef FE1ConvectionDiffusionElemDisc T;

			register_prepare_element_loop_function(	id, &T::template prepare_element_loop<TElem>);
			register_prepare_element_function(		id, &T::template prepare_element<TElem>);
			register_finish_element_loop_function(	id, &T::template finish_element_loop<TElem>);
			register_assemble_JA_function(			id, &T::template assemble_JA<TElem>);
			register_assemble_JM_function(			id, &T::template assemble_JM<TElem>);
			register_assemble_A_function(			id, &T::template assemble_A<TElem>);
			register_assemble_M_function(			id, &T::template assemble_M<TElem>);
			register_assemble_f_function(			id, &T::template assemble_f<TElem>);
		}

};

}

#include "fe1_convection_diffusion_impl.h"

#endif /*__H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__CONVECTION_DIFFUSION__FE1__CONVECTION_DIFFUSION__*/
