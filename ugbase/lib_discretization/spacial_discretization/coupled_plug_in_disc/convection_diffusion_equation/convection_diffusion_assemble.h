
#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__PLUG_IN_DISC__CPL_CONVECTION_DIFFUSION_EQUATION__CPL_CONVECTION_DIFFUSION_ASSEMBLE__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__PLUG_IN_DISC__CPL_CONVECTION_DIFFUSION_EQUATION__CPL_CONVECTION_DIFFUSION_ASSEMBLE__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"
#include "lib_algebra/lib_algebra.h"

// library intern headers
#include "lib_discretization/spacial_discretization/disc_helper/fvgeom.h"

#include "lib_discretization/spacial_discretization/plug_in_disc/plug_in_element_disc_interface.h"
#include "lib_discretization/spacial_discretization/disc_coupling/element_data.h"


namespace ug{


template<typename TDomain, typename TAlgebra, typename TElem >
class CplConvectionDiffusionEquation {
	public:
		typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

		// domain type
		typedef TDomain domain_type;

		// world dimension
		static const int dim = TDomain::dim;
		static const int ref_dim = ref_elem_type::dim;

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
		CplConvectionDiffusionEquation(TDomain& domain, number upwind_amount,
										Diff_Tensor_fct diff, Conv_Vel_fct vel, Reaction_fct reac, Rhs_fct rhs,
										DataImport<MathVector<dim>, MathVector<dim> >& Velocity)
		: m_domain(domain), m_upwind_amount(upwind_amount),
			m_Diff_Tensor(diff), m_Conv_Vel(vel), m_Reaction(reac), m_Rhs(rhs),
			m_Velocity(Velocity)
		{};

	public:

		// total number of shape functions on elements of type 'TElem'
		inline size_t num_sh(){return ref_elem_type::num_corners;};

		// number of shape functions on elements of type 'TElem' for the 'i'-th fundamental function
		inline size_t num_sh(size_t loc_fct){return ref_elem_type::num_corners;};

		// prepares the loop. Must be called, before prepare_element can be used
		inline IPlugInReturn prepare_element_loop();

		// prepares the element. Must be called before assemble_element_XXX can be used. Must be called after prepare_element_loop
		inline IPlugInReturn prepare_element(TElem* elem, const local_vector_type& u, const local_index_type& glob_ind);

		inline IPlugInReturn assemble_element_JA(local_matrix_type& J, const local_vector_type& u, number time=0.0);

		inline IPlugInReturn assemble_element_JM(local_matrix_type& J, const local_vector_type& u, number time=0.0);

		inline IPlugInReturn assemble_element_A(local_vector_type& d, const local_vector_type& u, number time=0.0);

		inline IPlugInReturn assemble_element_M(local_vector_type& d, const local_vector_type& u, number time=0.0);

		inline IPlugInReturn assemble_element_f(local_vector_type& d, number time=0.0);

		inline IPlugInReturn finish_element_loop();

	private:
		// domain
		TDomain& m_domain;

		// position access
		typename TDomain::position_type m_corners[ref_elem_type::num_corners];
		typename TDomain::position_accessor_type m_aaPos;

		// Finite Volume Element Geometry
		FVElementGeometry<TElem, dim>* m_geo;

		// amount of upwind (1.0 == full upwind, 0.0 == no upwind)
		number m_upwind_amount;

		// User functions
		Diff_Tensor_fct m_Diff_Tensor;
		Conv_Vel_fct m_Conv_Vel;
		Reaction_fct m_Reaction;
		Rhs_fct m_Rhs;

		DataImport<MathVector<dim>, MathVector<ref_elem_type::dim> >& m_Velocity;
};



template <	typename TDomain,
			int ref_dim,
			typename TAlgebra>
class CplConvectionDiffusionEquationPlugIn : public IPlugInElementDiscretization<TAlgebra>{

	public:
		// domain type
		typedef TDomain domain_type;

		// world dimension
		static const int dim = domain_type::dim;

		// position type
		typedef typename domain_type::position_type position_type;

		// algebra type
		typedef TAlgebra algebra_type;

		// local matrix type
		typedef typename algebra_type::matrix_type::local_matrix_type local_matrix_type;

		// local vector tyoe
		typedef typename algebra_type::vector_type::local_vector_type local_vector_type;

		// local index type
		typedef typename algebra_type::vector_type::local_index_type local_index_type;

	protected:
		typedef void (*Diff_Tensor_fct)(MathMatrix<dim,dim>&, const position_type&, number);
		typedef void (*Conv_Vel_fct)(position_type&, const position_type&, number);
		typedef void (*Reaction_fct)(number&, const position_type&, number);
		typedef void (*Rhs_fct)(number&, const position_type&, number);

	public:
		CplConvectionDiffusionEquationPlugIn(	size_t fct, domain_type& domain, number upwind_amount,
											Diff_Tensor_fct diff, Conv_Vel_fct vel, Reaction_fct reac, Rhs_fct rhs) :
			m_Velocity("Velocity"),
			m_fct(fct),
			m_domain(domain),
			m_upwind_amount(upwind_amount),
			m_diff(diff),
			m_vel(vel),
			m_reac(reac),
			m_rhs(rhs)
			{};


	public:
		size_t num_imports() {return 1;}

		DataImportItem* import(size_t i)
		{
			if(i != 0) return NULL;
			return &m_Velocity;
		}

		bool register_exports(DataContainer& Cont){ return true;}

		bool unregister_exports(DataContainer& Cont) {return true;}

		bool register_imports(DataContainer& Cont)
		{
			if(Cont.register_item(m_Velocity) != true)
			{
				UG_ASSERT(0, "Must work.");
				return false;
			}
			return true;
		}

		bool unregister_imports(DataContainer& Cont)
		{
			if(Cont.unregister_item(m_Velocity) != true)
			{
				UG_ASSERT(0, "Must work.");
				return false;
			}
			return true;
		}

	protected:
		DataImport<MathVector<dim>, MathVector<ref_dim> > m_Velocity;

	public:
		// number of fundamental functions required for this assembling
		inline size_t num_fct(){return 1;}

		// local shape function set required for the 'i'-th fundamental function
		inline LocalShapeFunctionSetID local_shape_function_set(size_t loc_fct)
		{
			UG_ASSERT(loc_fct < num_fct(), "Accessing fundamental function, that is not contained in this assembling.");
			return LSFS_LAGRANGEP1;
		}

		// global number of fundamental function of local fundamental function
		inline size_t fct(size_t loc_fct)
		{
			UG_ASSERT(loc_fct == 0, "Convection Diffusion Assembling has only one component.");
			return m_fct;
		}

	protected:
		// number of fundamental function, where this assembling works
		size_t m_fct;

		template<typename TElem>
		CplConvectionDiffusionEquation<domain_type, algebra_type, TElem>&
		get_inst(	domain_type& domain, number upwind_amount, Diff_Tensor_fct diff, Conv_Vel_fct vel,
					Reaction_fct reac, Rhs_fct rhs, DataImport<MathVector<dim>, MathVector<ref_dim> >& Velocity)
		{
			static CplConvectionDiffusionEquation<domain_type, algebra_type, TElem>
						inst(domain, upwind_amount, diff, vel, reac, rhs, Velocity);
			return inst;
		}

	public:
		// support assembling on triangles
		template <typename TElem>
		inline size_t num_sh(TElem* elem)
		{ return get_inst<TElem>(m_domain, m_upwind_amount, m_diff, m_vel, m_reac, m_rhs, m_Velocity).num_sh();};

		template <typename TElem>
		inline size_t num_sh(TElem* elem, size_t fct)
		{ return get_inst<TElem>(m_domain, m_upwind_amount, m_diff, m_vel, m_reac, m_rhs, m_Velocity).num_sh(fct);};

		template <typename TElem>
		inline IPlugInReturn prepare_element_loop(TElem* elem)
		{ return get_inst<TElem>(m_domain, m_upwind_amount, m_diff, m_vel, m_reac, m_rhs, m_Velocity).prepare_element_loop(); };

		template <typename TElem>
		inline IPlugInReturn prepare_element(TElem* elem, const local_vector_type& u, const local_index_type& glob_ind)
		{ return get_inst<TElem>(m_domain, m_upwind_amount, m_diff, m_vel, m_reac, m_rhs, m_Velocity).prepare_element(elem, u, glob_ind); };

		template <typename TElem>
		inline IPlugInReturn assemble_element_JA(TElem* elem, local_matrix_type& J, const local_vector_type& u, number time=0.0)
		{ return get_inst<TElem>(m_domain, m_upwind_amount, m_diff, m_vel, m_reac, m_rhs, m_Velocity).assemble_element_JA(J, u, time); };

		template <typename TElem>
		inline IPlugInReturn assemble_element_JM(TElem* elem, local_matrix_type& J, const local_vector_type& u, number time=0.0)
		{ return get_inst<TElem>(m_domain, m_upwind_amount, m_diff, m_vel, m_reac, m_rhs, m_Velocity).assemble_element_JM(J, u, time); };

		template <typename TElem>
		inline IPlugInReturn assemble_element_A(TElem* elem, local_vector_type& d, const local_vector_type& u, number time=0.0)
		{ return get_inst<TElem>(m_domain, m_upwind_amount, m_diff, m_vel, m_reac, m_rhs, m_Velocity).assemble_element_A(d, u, time); };

		template <typename TElem>
		inline IPlugInReturn assemble_element_M(TElem* elem, local_vector_type& d, const local_vector_type& u, number time=0.0)
		{ return get_inst<TElem>(m_domain, m_upwind_amount, m_diff, m_vel, m_reac, m_rhs, m_Velocity).assemble_element_M(d, u, time); };

		template <typename TElem>
		inline IPlugInReturn assemble_element_f(TElem* elem, local_vector_type& d, number time=0.0)
		{ return get_inst<TElem>(m_domain, m_upwind_amount, m_diff, m_vel, m_reac, m_rhs, m_Velocity).assemble_element_f(d, time); };

		template <typename TElem>
		inline IPlugInReturn finish_element_loop(TElem* elem)
		{ return get_inst<TElem>(m_domain, m_upwind_amount, m_diff, m_vel, m_reac, m_rhs, m_Velocity).finish_element_loop(); };

	protected:
		domain_type& m_domain;
		number m_upwind_amount;
		Diff_Tensor_fct m_diff;
		Conv_Vel_fct m_vel;
		Reaction_fct m_reac;
		Rhs_fct m_rhs;

};

}

#include "convection_diffusion_assemble_impl.h"

#endif /*__H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__PLUG_IN_DISC__CPL_CONVECTION_DIFFUSION_EQUATION__CPL_CONVECTION_DIFFUSION_ASSEMBLE__*/
