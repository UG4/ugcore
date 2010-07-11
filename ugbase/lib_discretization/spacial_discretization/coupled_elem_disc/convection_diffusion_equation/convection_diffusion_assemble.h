
#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__CPL_CONVECTION_DIFFUSION_EQUATION__CPL_CONVECTION_DIFFUSION_ASSEMBLE__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__CPL_CONVECTION_DIFFUSION_EQUATION__CPL_CONVECTION_DIFFUSION_ASSEMBLE__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"
#include "lib_algebra/lib_algebra.h"

// library intern headers
#include "lib_discretization/spacial_discretization/disc_helper/fvgeom.h"

#include "../coupled_elem_disc_interface.h"
#include "../disc_coupling/element_data.h"


namespace ug{


template<typename TDomain, int ref_dim, typename TAlgebra>
class CplConvectionDiffusionElemDisc : public ICoupledElemDisc<TAlgebra> {
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
		typedef void (*Conv_Scale_fct)(number&, const position_type&, number);
		typedef void (*Mass_Scale_fct)(number&, const position_type&, number);
		typedef void (*Mass_Const_fct)(number&, const position_type&, number);
		typedef void (*Reaction_fct)(number&, const position_type&, number);
		typedef void (*Rhs_fct)(number&, const position_type&, number);

	public:
		CplConvectionDiffusionElemDisc(TDomain& domain, number upwind_amount,
										Diff_Tensor_fct diff, Conv_Scale_fct conv_scale,
										Mass_Scale_fct mass_scale, Mass_Const_fct mass_const,
										Reaction_fct reac, Rhs_fct rhs);


	public:
		virtual size_t num_fct(){return 1;}

		virtual LocalShapeFunctionSetID local_shape_function_set_id(size_t loc_fct) {return LSFS_LAGRANGEP1;}

	public:
		virtual size_t num_imports() {return 1;}

		virtual DataImportItem* import(size_t i)
		{
			if(i != 0) return NULL;
			return &m_Velocity;
		}

		virtual bool register_exports(DataContainer& Cont){ return true;}

		virtual bool unregister_exports(DataContainer& Cont) {return true;}

		virtual bool register_imports(DataContainer& Cont)
		{
			if(Cont.register_item(m_Velocity) != true) return false;
			return true;
		}

		virtual bool unregister_imports(DataContainer& Cont)
		{
			if(Cont.unregister_item(m_Velocity) != true) return false;
			return true;
		}

		virtual bool set_sys_id(size_t sys_id) {return true;}

	protected:
		DataImport<MathVector<dim>, MathVector<ref_dim> > m_Velocity;

	public:
		template <typename TElem>
		inline size_t num_total_sh(){	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
										return ref_elem_type::num_corners;};

		template <typename TElem>
		inline size_t num_sh(size_t loc_fct){		typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
													return ref_elem_type::num_corners;};

		template <typename TElem>
		inline bool prepare_element_loop();

		template <typename TElem>
		inline bool finish_element_loop();

		template <typename TElem>
		inline bool prepare_element(TElem* elem, const local_vector_type& u, const local_index_type& glob_ind);

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
		Conv_Scale_fct m_Conv_Scale;
		Mass_Scale_fct m_Mass_Scale;
		Mass_Const_fct m_Mass_Const;
		Reaction_fct m_Reaction;
		Rhs_fct m_Rhs;

	private:
		// help function
		template <typename TElem>
		void register_all_assemble_functions(int id)
		{
			register_num_total_sh_function(			id, &CplConvectionDiffusionElemDisc::template num_total_sh<TElem>);
			register_num_sh_function(				id, &CplConvectionDiffusionElemDisc::template num_sh<TElem>);
			register_prepare_element_loop_function(	id, &CplConvectionDiffusionElemDisc::template prepare_element_loop<TElem>);
			register_prepare_element_function(		id, &CplConvectionDiffusionElemDisc::template prepare_element<TElem>);
			register_finish_element_loop_function(	id, &CplConvectionDiffusionElemDisc::template finish_element_loop<TElem>);
			register_assemble_JA_function(			id, &CplConvectionDiffusionElemDisc::template assemble_JA<TElem>);
			register_assemble_JM_function(			id, &CplConvectionDiffusionElemDisc::template assemble_JM<TElem>);
			register_assemble_A_function(			id, &CplConvectionDiffusionElemDisc::template assemble_A<TElem>);
			register_assemble_M_function(			id, &CplConvectionDiffusionElemDisc::template assemble_M<TElem>);
			register_assemble_f_function(			id, &CplConvectionDiffusionElemDisc::template assemble_f<TElem>);
		}
};


} // end namespace ug

#include "convection_diffusion_assemble_impl.h"

#endif /*__H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__CPL_CONVECTION_DIFFUSION_EQUATION__CPL_CONVECTION_DIFFUSION_ASSEMBLE__*/
