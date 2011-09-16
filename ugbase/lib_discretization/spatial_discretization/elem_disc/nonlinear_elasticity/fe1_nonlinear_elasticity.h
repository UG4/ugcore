/*
 * fe1_nonlinear_elasticity.h
 *
 *  Created on: 18.05.2011
 *      Author: raphaelprohl
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NONLINEAR_ELASTICITY__FE1_NONLINEAR_ELASTICITY__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NONLINEAR_ELASTICITY__FE1_NONLINEAR_ELASTICITY__

#include <boost/function.hpp>

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_discretization/spatial_discretization/elem_disc/elem_disc_interface.h"

namespace ug{


template<typename TDomain>
class FE1NonlinearElasticityElemDisc
	: public IDomainElemDisc<TDomain>
{
	private:
	///	Base class type
		typedef IDomainElemDisc<TDomain> base_type;

	///	own type
		typedef FE1NonlinearElasticityElemDisc<TDomain> this_type;

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
		typedef void (*Elasticity_Tensor_fct)(MathTensor<4,dim>&);
		typedef void (*Stress_Tensor_fct)(MathTensor<2,dim>&);
		typedef boost::function<void (MathTensor<4,dim>& value,
		                              const MathVector<dim>& x,
		                              number time)> TensorFunctor;

	public:
		FE1NonlinearElasticityElemDisc();

	///	set the elasticity tensor
		void set_elasticity_tensor(const TensorFunctor& elast)
		{
			m_ElasticityTensorFunctor = elast;
		}

	///	number of functions used
		virtual size_t num_fct(){return dim;}

	///	type of trial space for each function used
		virtual bool request_finite_element_id(const std::vector<LFEID>& vLfeID)
		{
		//	check number
			if(vLfeID.size() != num_fct()) return false;

		//	check that Lagrange 1st order
			for(size_t i = 0; i < vLfeID.size(); ++i)
				if(vLfeID[i] != LFEID(LFEID::LAGRANGE, 1)) return false;
			return true;
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
		bool prepare_element(TElem* elem, const local_vector_type& u);

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
		TensorFunctor m_ElasticityTensorFunctor;

		Elasticity_Tensor_fct m_ElasticityTensorFct;
		MathTensor<4, dim> m_ElasticityTensor;

		Stress_Tensor_fct m_StressTensorFct;
		MathTensor<2, dim> m_StressTensor;

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

#endif /*__H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NONLINEAR_ELASTICITY__FE1_NONLINEAR_ELASTICITY__*/
