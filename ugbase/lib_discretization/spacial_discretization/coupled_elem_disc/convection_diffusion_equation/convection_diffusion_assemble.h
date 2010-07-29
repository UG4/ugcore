
#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__CPL_CONVECTION_DIFFUSION_EQUATION__CPL_CONVECTION_DIFFUSION_ASSEMBLE__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__CPL_CONVECTION_DIFFUSION_EQUATION__CPL_CONVECTION_DIFFUSION_ASSEMBLE__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"
#include "lib_algebra/lib_algebra.h"

// library intern headers
#include "lib_discretization/spacial_discretization/disc_helper/fvgeom.h"

#include "../coupled_elem_disc_interface.h"
#include "../elem_data/element_data.h"
#include "lib_discretization/common/local_algebra.h"

namespace ug{


template<typename TDomain, typename TAlgebra>
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
		typedef LocalMatrix<typename TAlgebra::matrix_type::entry_type> local_matrix_type;

		// local vector type
		typedef LocalVector<typename TAlgebra::vector_type::entry_type> local_vector_type;

		// local index type
		//typedef typename algebra_type::vector_type::local_index_type local_index_type;
		typedef LocalIndices local_index_type;

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
		DataImport<MathVector<dim> > m_Velocity;

	public:
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

		// to make it more readable
		const static size_t _C_ = 0;

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
