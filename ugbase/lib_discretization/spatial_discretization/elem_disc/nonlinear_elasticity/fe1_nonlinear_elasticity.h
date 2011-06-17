/*
 * fe1_nonlinear_elasticity.h
 *
 *  Created on: 18.05.2011
 *      Author: raphaelprohl
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NONLINEAR_ELASTICITY__FE1_NONLINEAR_ELASTICITY__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NONLINEAR_ELASTICITY__FE1_NONLINEAR_ELASTICITY__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_discretization/spatial_discretization/elem_disc/elem_disc_interface.h"
#include "lib_discretization/common/local_algebra.h"

namespace ug{


template<typename TDomain, typename TAlgebra>
class FE1NonlinearElasticityElemDisc
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

	protected:
		typedef void (*Elasticity_Tensor_fct)(MathTensor<4,dim>&);
		typedef void (*Stress_Tensor_fct)(MathTensor<2,dim>&);

	public:
		FE1NonlinearElasticityElemDisc();

	///	set the elasticity tensor
		void set_elasticity_tensor(IUserData<MathTensor<4,dim>, dim>& elast)
		{
			m_ElasticityTensorFunctor = elast.get_functor();
		}

	///	number of functions used
		virtual size_t num_fct(){return dim;}

	///	type of trial space for each function used
		virtual LocalShapeFunctionSetID local_shape_function_set_id(size_t loc_fct)
		{
			return LocalShapeFunctionSetID(LocalShapeFunctionSetID::LAGRANGE, 1);
		}

	///	switches between non-regular and regular grids
		virtual bool treat_non_regular_grid(bool bNonRegular)
		{
		//	no special care for non-regular grids
			return true;
		}

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
		// position access
		const position_type* m_corners;

		// User functions
		typename IUserData<MathTensor<4,dim>, dim>::functor_type m_ElasticityTensorFunctor;

		Elasticity_Tensor_fct m_ElasticityTensorFct;
		MathTensor<4, dim> m_ElasticityTensor;

		Stress_Tensor_fct m_StressTensorFct;
		MathTensor<2, dim> m_StressTensor;

	private:
		///////////////////////////////////////
		// registering for reference elements
		///////////////////////////////////////
		template <int dim> class numType{};

		void register_assemble_functions()
		{
			numType<TDomain::dim> dummy;
			register_assemble_functions(dummy);
		}

		// register for 1D
		void register_assemble_functions(numType<1> dummy)
		{
			register_all_assemble_functions<Edge>(ROID_EDGE);
		}

		// register for 2D
		void register_assemble_functions(numType<2> dummy)
		{
			register_all_assemble_functions<Edge>(ROID_EDGE);
			register_all_assemble_functions<Triangle>(ROID_TRIANGLE);
			register_all_assemble_functions<Quadrilateral>(ROID_QUADRILATERAL);
		}

		// register for 3D
		void register_assemble_functions(numType<3> dummy)
		{
			register_all_assemble_functions<Edge>(ROID_EDGE);
			register_all_assemble_functions<Triangle>(ROID_TRIANGLE);
			register_all_assemble_functions<Quadrilateral>(ROID_QUADRILATERAL);
			register_all_assemble_functions<Tetrahedron>(ROID_TETRAHEDRON);
			register_all_assemble_functions<Prism>(ROID_PYRAMID);
			register_all_assemble_functions<Prism>(ROID_PRISM);
			register_all_assemble_functions<Prism>(ROID_HEXAHEDRON);
		}

		// help function
		template <typename TElem>
		void register_all_assemble_functions(int id)
		{
			register_prepare_element_loop_function(	id, &FE1NonlinearElasticityElemDisc::template prepare_element_loop<TElem>);
			register_prepare_element_function(		id, &FE1NonlinearElasticityElemDisc::template prepare_element<TElem>);
			register_finish_element_loop_function(	id, &FE1NonlinearElasticityElemDisc::template finish_element_loop<TElem>);
			register_assemble_JA_function(			id, &FE1NonlinearElasticityElemDisc::template assemble_JA<TElem>);
			register_assemble_JM_function(			id, &FE1NonlinearElasticityElemDisc::template assemble_JM<TElem>);
			register_assemble_A_function(			id, &FE1NonlinearElasticityElemDisc::template assemble_A<TElem>);
			register_assemble_M_function(			id, &FE1NonlinearElasticityElemDisc::template assemble_M<TElem>);
			register_assemble_f_function(			id, &FE1NonlinearElasticityElemDisc::template assemble_f<TElem>);
		}

};

}

#include "fe1_nonlinear_elasticity_impl.h"

#endif /*__H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NONLINEAR_ELASTICITY__FE1_NONLINEAR_ELASTICITY__*/
