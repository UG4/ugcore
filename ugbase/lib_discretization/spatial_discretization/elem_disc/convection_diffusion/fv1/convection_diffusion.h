/*
 * convection_diffusion.h
 *
 *  Created on: 26.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__CONVECTION_DIFFUSION__FV1__CONVECTION_DIFFUSION__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__CONVECTION_DIFFUSION__FV1__CONVECTION_DIFFUSION__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"
#include "lib_algebra/lib_algebra.h"

// library intern headers
#include "lib_discretization/spatial_discretization/disc_helper/disc_helper.h"
#include "lib_discretization/spatial_discretization/elem_disc/elem_disc_interface.h"
#include "lib_discretization/common/local_algebra.h"
#include "lib_discretization/spatial_discretization/ip_data/user_data.h"

namespace ug{

template<	template <class TElem, int TWorldDim> class TFVGeom,
			typename TDomain,
			typename TAlgebra>
class FVConvectionDiffusionElemDisc
	: public IElemDisc<TAlgebra>
{
	public:
	// 	Domain type
		typedef TDomain domain_type;

	// 	World dimension
		static const int dim = TDomain::dim;

	// 	Position type
		typedef typename TDomain::position_type position_type;

	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Local matrix type
		typedef LocalMatrix<typename TAlgebra::matrix_type::value_type> local_matrix_type;

	// 	Local vector type
		typedef LocalVector<typename TAlgebra::vector_type::value_type> local_vector_type;

	// 	Local index type
		typedef LocalIndices local_index_type;

	protected:
	//	Functor type for Diffusion
	//	typedef typename  IUserData<MathMatrix<dim, dim>, dim>::functor_type DiffusionFunctor;

	//	Functor type for Convection Velocity
	//	typedef typename  IUserData<MathVector<dim>, dim>::functor_type VelocityFunctor;

	public:
	//	Constructor
		FVConvectionDiffusionElemDisc()
		 : m_pDomain(NULL), m_upwindAmount(0.0)
			{
			//	register assemling functions
				register_assemble_functions(Int2Type<dim>());

			//	register imports
				register_import(m_Diff);
				register_import(m_ConvVel);
				register_import(m_Reaction);
				register_import(m_Rhs);
			}

		void set_upwind_amount(number amount) {m_upwindAmount = amount;}
		void set_domain(domain_type& domain) {m_pDomain = &domain;}

		void set_diffusion(IPData<MathMatrix<dim, dim>, dim>& user) {m_Diff.set_data(user);}
		void set_velocity(IPData<MathVector<dim>, dim>& user) {m_ConvVel.set_data(user);}
		void set_reaction(IPData<number, dim>& user) {m_Reaction.set_data(user);}
		void set_rhs(IPData<number, dim>& user)	{m_Rhs.set_data(user);}

		virtual size_t num_fct(){return 1;}

		virtual LocalShapeFunctionSetID local_shape_function_set_id(size_t loc_fct)
		{
			return LocalShapeFunctionSetID(LocalShapeFunctionSetID::LAGRANGE, 1);
		}

		virtual bool use_hanging() const {return TFVGeom<Edge, dim>::usesHangingNodes;}

	private:
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

		template <typename TElem>
		bool lin_defect_velocity(const local_vector_type& u)
		{
		// get finite volume geometry
			TFVGeom<TElem, dim>& geo = FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

		// loop Sub Control Volume Faces (SCVF)
			size_t ip = 0;
			for(size_t i = 0; i < geo.num_scvf(); ++i)
			{
			// get current SCVF
				const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(i);

			// loop integration point of SCVF
				for(size_t i = 0; i < scvf.num_ip(); ++i, ++ip)
				{
				// loop shape functions
					for(size_t j = 0; j < scvf.num_sh(); ++j)
					{
						m_ConvVel.lin_defect(ip, _C_, j) = 0.0;
					}

				// central part convection
					number scale = (1.- m_upwindAmount);

					number shape_u = 0.0;
					for(size_t j = 0; j < scvf.num_sh(); ++j)
						shape_u += u(_C_,j) * scvf.shape(j, i);

					scale *= shape_u;

				// upwind part convection
					const number flux = m_upwindAmount
										* VecDot(m_ConvVel[ip], scvf.normal());

					if(flux >= 0.0) scale += m_upwindAmount* u(_C_, scvf.from());
					else scale += m_upwindAmount * u(_C_, scvf.to());

					MathVector<dim> linDefect = scvf.normal();
					linDefect *= scale;

					m_ConvVel.lin_defect(ip, _C_, scvf.from()) += linDefect;
					m_ConvVel.lin_defect(ip, _C_, scvf.to()) -= linDefect;
				}
			}

		//	we're done
			return true;
		}

	private:
		// domain
		TDomain* m_pDomain;

		// position access
		std::vector<position_type> m_vCornerCoords;
		typename TDomain::position_accessor_type m_aaPos;

		// to make it more readable
		static const size_t _C_ = 0;

		// amount of upwind (1.0 == full upwind, 0.0 == no upwind)
		number m_upwindAmount;

		// User functions
		DataImport<MathMatrix<dim,dim>, dim, algebra_type> m_Diff;
		DataImport<MathVector<dim>, dim, algebra_type > m_ConvVel;
		DataImport<number, dim, algebra_type> m_Reaction;
		DataImport<number, dim, algebra_type> m_Rhs;

	private:
		///////////////////////////////////////
		// registering for reference elements
		///////////////////////////////////////
		// register for 1D
		void register_assemble_functions(Int2Type<1>)
		{
			register_all_assemble_functions<Edge>(ROID_EDGE);
		}

		// register for 2D
		void register_assemble_functions(Int2Type<2>)
		{
			register_assemble_functions(Int2Type<1>());
			register_all_assemble_functions<Triangle>(ROID_TRIANGLE);
			register_all_assemble_functions<Quadrilateral>(ROID_QUADRILATERAL);
		}

		// register for 3D
		void register_assemble_functions(Int2Type<3>)
		{
			register_assemble_functions(Int2Type<2>());
			register_all_assemble_functions<Tetrahedron>(ROID_TETRAHEDRON);
			register_all_assemble_functions<Pyramid>(ROID_PYRAMID);
			register_all_assemble_functions<Prism>(ROID_PRISM);
			register_all_assemble_functions<Hexahedron>(ROID_HEXAHEDRON);
		}

		// help function
		template <typename TElem>
		void register_all_assemble_functions(int id)
		{
			typedef FVConvectionDiffusionElemDisc T;

			register_prepare_element_loop_function(	id, &T::template prepare_element_loop<TElem>);
			register_prepare_element_function(		id, &T::template prepare_element<TElem>);
			register_finish_element_loop_function(	id, &T::template finish_element_loop<TElem>);
			register_assemble_JA_function(			id, &T::template assemble_JA<TElem>);
			register_assemble_JM_function(			id, &T::template assemble_JM<TElem>);
			register_assemble_A_function(			id, &T::template assemble_A<TElem>);
			register_assemble_M_function(			id, &T::template assemble_M<TElem>);
			register_assemble_f_function(			id, &T::template assemble_f<TElem>);

			m_ConvVel.register_lin_defect_func(id, this, &T::template lin_defect_velocity<TElem>);
		}

};

}

#include "convection_diffusion_impl.h"

#endif /*__H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__CONVECTION_DIFFUSION__FV1__CONVECTION_DIFFUSION__*/
