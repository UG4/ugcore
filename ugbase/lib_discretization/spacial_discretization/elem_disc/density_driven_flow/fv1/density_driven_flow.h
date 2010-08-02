/*
 * density_driven_flow.h
 *
 *  Created on: 26.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__DENSITY_DRIVEN_FLOW__DENSITY_DRIVEN_FLOW__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__DENSITY_DRIVEN_FLOW__DENSITY_DRIVEN_FLOW__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"
#include "lib_algebra/lib_algebra.h"

// library intern headers
#include "lib_discretization/spacial_discretization/disc_helper/fvgeom.h"
#include "lib_discretization/spacial_discretization/elem_disc/elem_disc_interface.h"
#include "lib_discretization/common/local_algebra.h"

namespace ug{

template<typename TDomain, typename TAlgebra>
class DensityDrivenFlowElemDisc  : public IElemDisc<TAlgebra> {
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
		typedef LocalMatrix<typename TAlgebra::matrix_type::entry_type> local_matrix_type;

		// local vector type
		typedef LocalVector<typename TAlgebra::vector_type::entry_type> local_vector_type;

		// local index type
		typedef LocalIndices local_index_type;

	protected:
		typedef void (*Pososity_fct)(number&);
		typedef void (*Viscosity_fct)(number&, number);
		typedef void (*Density_fct)(number&, number);
		typedef void (*D_Density_fct)(number&, number);
		typedef void (*Mol_Diff_Tensor_fct)(MathMatrix<dim,dim>&);
		typedef void (*Permeability_Tensor_fct)(MathMatrix<dim,dim>&);
		typedef void (*Gravity_fct)(MathVector<dim>&);

	public:
		DensityDrivenFlowElemDisc(	TDomain& domain, number upwind_amount,
									Pososity_fct Porosity, Viscosity_fct Viscosity, Density_fct Density, D_Density_fct D_Density,
									Mol_Diff_Tensor_fct Mol_Diff, Permeability_Tensor_fct Permeability_Tensor, Gravity_fct Gravity);

		virtual size_t num_fct() {return 2;}

		virtual LocalShapeFunctionSetID local_shape_function_set_id(size_t loc_fct)	{return LSFS_LAGRANGEP1;}

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

	private:
		// domain
		TDomain& m_domain;

		// position access
		typename TDomain::position_type* m_corners;
		typename TDomain::position_accessor_type m_aaPos;

		template <typename TElem>
		ReferenceMapping<typename reference_element_traits<TElem>::reference_element_type, dim>& get_mapping()
		{
			static ReferenceMapping<typename reference_element_traits<TElem>::reference_element_type, dim> mapping;
			return mapping;
		}

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
		Pososity_fct m_Porosity;
		Viscosity_fct m_Viscosity;
		Density_fct m_Density;
		D_Density_fct m_D_Density;
		Mol_Diff_Tensor_fct m_Mol_Diff_Tensor;
		Permeability_Tensor_fct m_Permeability_Tensor;
		Gravity_fct m_Gravity;

		// constant values, read in once
		number m_porosity;

	private:
		// local constants for readability (local function 0 == _C_, local function 1 == _P_)
		static const size_t _C_ = 0;
		static const size_t _P_ = 1;

	protected:
		bool compute_ip_Darcy_velocity(MathVector<dim>& Darcy_vel, number c_ip, const MathVector<dim>& grad_p_ip);

		template <typename TElem>
		bool compute_D_ip_Darcy_velocity(	const SubControlVolumeFace<TElem, dim>& scvf,
											MathVector<dim>& Darcy_vel, MathVector<dim> D_Darcy_vel_c[], MathVector<dim> D_Darcy_vel_p[],
											number c_ip, const MathVector<dim>& grad_p_ip);

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
			// TODO: Register 3D Ref-Elems
		}

		// help function
		template <typename TElem>
		void register_all_assemble_functions(int id)
		{
			register_prepare_element_loop_function(	id, &DensityDrivenFlowElemDisc::template prepare_element_loop<TElem>);
			register_prepare_element_function(		id, &DensityDrivenFlowElemDisc::template prepare_element<TElem>);
			register_finish_element_loop_function(	id, &DensityDrivenFlowElemDisc::template finish_element_loop<TElem>);
			register_assemble_JA_function(			id, &DensityDrivenFlowElemDisc::template assemble_JA<TElem>);
			register_assemble_JM_function(			id, &DensityDrivenFlowElemDisc::template assemble_JM<TElem>);
			register_assemble_A_function(			id, &DensityDrivenFlowElemDisc::template assemble_A<TElem>);
			register_assemble_M_function(			id, &DensityDrivenFlowElemDisc::template assemble_M<TElem>);
			register_assemble_f_function(			id, &DensityDrivenFlowElemDisc::template assemble_f<TElem>);
		}
};


} // end namespace ug

#include "density_driven_flow_impl.h"

#endif /*__H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__DENSITY_DRIVEN_FLOW__DENSITY_DRIVEN_FLOW__*/
