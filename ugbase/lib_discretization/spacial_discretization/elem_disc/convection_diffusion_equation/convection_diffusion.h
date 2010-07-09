/*
 * convection_diffusion.h
 *
 *  Created on: 26.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__CONVECTION_DIFFUSION_EQUATION__CONVECTION_DIFFUSION__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__CONVECTION_DIFFUSION_EQUATION__CONVECTION_DIFFUSION__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"
#include "lib_algebra/lib_algebra.h"

// library intern headers
#include "lib_discretization/spacial_discretization/disc_helper/fvgeom.h"
#include "lib_discretization/spacial_discretization/elem_disc/elem_disc_interface.h"


namespace ug{


template<typename TDomain, typename TAlgebra>
class ConvectionDiffusionElemDisc : public IElemDisc<TAlgebra>
{
	public:
		// domain type
		typedef TDomain domain_type;

		// world dimension
		static const int dim = TDomain::dim;

		// position type
		typedef typename TDomain::position_type position_type;

		// algebra type
		typedef TAlgebra algebra_type;

		// local matrix type
		typedef typename algebra_type::matrix_type::local_matrix_type local_matrix_type;

		// local vector type
		typedef typename algebra_type::vector_type::local_vector_type local_vector_type;

		// local index type
		typedef typename algebra_type::vector_type::local_index_type local_index_type;

	protected:
		typedef void (*Diff_Tensor_fct)(MathMatrix<dim,dim>&, const position_type&, number);
		typedef void (*Conv_Vel_fct)(position_type&, const position_type&, number);
		typedef void (*Reaction_fct)(number&, const position_type&, number);
		typedef void (*Rhs_fct)(number&, const position_type&, number);

	public:
		ConvectionDiffusionElemDisc(TDomain& domain, number upwind_amount,
									Diff_Tensor_fct diff, Conv_Vel_fct vel, Reaction_fct reac, Rhs_fct rhs);

		virtual size_t num_fct(){return 1;}

		virtual LocalShapeFunctionSetID local_shape_function_set_id(size_t loc_fct) {return LSFS_LAGRANGEP1;}

	private:
		template <typename TElem>
		inline size_t num_total_sh(){ 	typedef typename reference_element_traits<TElem>::reference_element_type reference_element_type;
										return reference_element_type::num_corners;};

		template <typename TElem>
		inline size_t num_sh(size_t loc_fct){ 	typedef typename reference_element_traits<TElem>::reference_element_type reference_element_type;
												return reference_element_type::num_corners;};

		template <typename TElem>
		inline bool prepare_element_loop();

		template <typename TElem>
		inline bool prepare_element(TElem* elem, const local_vector_type& u, const local_index_type& glob_ind);

		template <typename TElem>
		inline bool finish_element_loop();

		template <typename TElem>
		inline bool assemble_JA(local_matrix_type& J, const local_vector_type& u, number time=0.0);

		template <typename TElem>
		inline bool assemble_JM(local_matrix_type& J, const local_vector_type& u, number time=0.0);

		template <typename TElem>
		inline bool assemble_A(local_vector_type& d, const local_vector_type& u, number time=0.0);

		template <typename TElem>
		inline bool assemble_M(local_vector_type& d, const local_vector_type& u, number time=0.0);

		template <typename TElem>
		inline bool assemble_f(local_vector_type& d, number time=0.0);

	private:
		// domain
		TDomain& m_domain;

		// position access
		position_type* m_corners;
		typename TDomain::position_accessor_type m_aaPos;

		// Finite Volume Element Geometry
		template <typename TElem>
		inline FVElementGeometry<TElem, dim>& get_fvgeom()
		{
			static FVElementGeometry<TElem, dim> geo;
			return geo;
		}

		// amount of upwind (1.0 == full upwind, 0.0 == no upwind)
		number m_upwind_amount;

		// User functions
		Diff_Tensor_fct m_Diff_Tensor;
		Conv_Vel_fct m_Conv_Vel;
		Reaction_fct m_Reaction;
		Rhs_fct m_Rhs;

	private:
		// help function
		template <typename TElem>
		void register_all_assemble_functions(int id)
		{
			register_num_total_sh_function(			id, &ConvectionDiffusionElemDisc::template num_total_sh<TElem>);
			register_num_sh_function(				id, &ConvectionDiffusionElemDisc::template num_sh<TElem>);
			register_prepare_element_loop_function(	id, &ConvectionDiffusionElemDisc::template prepare_element_loop<TElem>);
			register_prepare_element_function(		id, &ConvectionDiffusionElemDisc::template prepare_element<TElem>);
			register_finish_element_loop_function(	id, &ConvectionDiffusionElemDisc::template finish_element_loop<TElem>);
			register_assemble_JA_function(			id, &ConvectionDiffusionElemDisc::template assemble_JA<TElem>);
			register_assemble_JM_function(			id, &ConvectionDiffusionElemDisc::template assemble_JM<TElem>);
			register_assemble_A_function(			id, &ConvectionDiffusionElemDisc::template assemble_A<TElem>);
			register_assemble_M_function(			id, &ConvectionDiffusionElemDisc::template assemble_M<TElem>);
			register_assemble_f_function(			id, &ConvectionDiffusionElemDisc::template assemble_f<TElem>);
		}

};

}

#include "convection_diffusion_impl.h"

#endif /*__H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__CONVECTION_DIFFUSION_EQUATION__CONVECTION_DIFFUSION__*/
