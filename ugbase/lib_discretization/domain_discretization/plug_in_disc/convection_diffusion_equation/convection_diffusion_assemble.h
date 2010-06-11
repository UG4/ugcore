
#ifndef __H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__PLUG_IN_DISC__CONVECTION_DIFFUSION_EQUATION__CONVECTION_DIFFUSION_ASSEMBLE__
#define __H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__PLUG_IN_DISC__CONVECTION_DIFFUSION_EQUATION__CONVECTION_DIFFUSION_ASSEMBLE__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lib_grid.h"
#include "lib_algebra/lib_algebra.h"

// library intern headers
#include "lib_discretization/domain_discretization/disc_helper/fvgeom.h"

#include "lib_discretization/domain_discretization/plug_in_disc/plug_in_element_disc_interface.h"
#include "lib_discretization/domain_discretization/disc_coupling/element_data.h"


namespace ug{


template<typename TDomain, typename TAlgebra, typename TElem >
class ConvectionDiffusionEquation {
	public:
		// forward constants and types

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
		ConvectionDiffusionEquation(TDomain& domain, number upwind_amount, Diff_Tensor_fct diff, Conv_Vel_fct vel, Reaction_fct reac, Rhs_fct rhs, DataImport<MathVector<dim>, MathVector<dim> >& Velocity)
		: m_domain(domain), m_upwind_amount(upwind_amount), m_Diff_Tensor(diff), m_Conv_Vel(vel), m_Reaction(reac), m_Rhs(rhs), m_Velocity(Velocity)
		{};

	public:

		// total number of shape functions on elements of type 'TElem'
		inline size_t num_sh(){return reference_element_traits<TElem>::num_corners;};

		// number of shape functions on elements of type 'TElem' for the 'i'-th fundamental function
		inline size_t num_sh(size_t loc_fct){return reference_element_traits<TElem>::num_corners;};

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
		typename TDomain::position_type m_corners[reference_element_traits<TElem>::num_corners];
		typename TDomain::position_accessor_type m_aaPos;

		// Finite Volume Element Geometry
		FVElementGeometry<TElem>* m_geo;

		// amount of upwind (1.0 == full upwind, 0.0 == no upwind)
		number m_upwind_amount;

		// User functions
		Diff_Tensor_fct m_Diff_Tensor;
		Conv_Vel_fct m_Conv_Vel;
		Reaction_fct m_Reaction;
		Rhs_fct m_Rhs;

		DataImport<MathVector<dim>, MathVector<dim> >& m_Velocity;
};



template <typename TDomain, typename TAlgebra>
class ConvectionDiffusionEquationPlugIn : public IPlugInElementDiscretization<TAlgebra>{

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
		ConvectionDiffusionEquationPlugIn(size_t fct, domain_type& domain, number upwind_amount, Diff_Tensor_fct diff, Conv_Vel_fct vel, Reaction_fct reac, Rhs_fct rhs) :
			m_fct(fct),
			m_Velocity("Velocity"),
			m_ImpTriangle(domain, upwind_amount, diff, vel, reac, rhs, m_Velocity),
			m_ImpQuadrilateral(domain, upwind_amount, diff, vel, reac, rhs, m_Velocity)
			{};

	public:
		/* GENERAL INFORMATIONS */
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

	public:
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
		DataImport<MathVector<dim>, MathVector<dim> > m_Velocity;

	/* ELEMENT WISE ASSEMBLNG */
	public:
		// support assembling on triangles
		inline size_t num_sh(Triangle* elem)
		{ return m_ImpTriangle.num_sh();};

		inline size_t num_sh(Triangle* elem, size_t fct)
		{ return m_ImpTriangle.num_sh(fct);};

		inline IPlugInReturn prepare_element_loop(Triangle* elem)
		{ return m_ImpTriangle.prepare_element_loop(); };

		inline IPlugInReturn prepare_element(Triangle* elem, const local_vector_type& u, const local_index_type& glob_ind)
		{ return m_ImpTriangle.prepare_element(elem, u, glob_ind); };

		inline IPlugInReturn assemble_element_JA(Triangle* elem, local_matrix_type& J, const local_vector_type& u, number time=0.0)
		{ return m_ImpTriangle.assemble_element_JA(J, u, time); };

		inline IPlugInReturn assemble_element_JM(Triangle* elem, local_matrix_type& J, const local_vector_type& u, number time=0.0)
		{ return m_ImpTriangle.assemble_element_JM(J, u, time); };

		inline IPlugInReturn assemble_element_A(Triangle* elem, local_vector_type& d, const local_vector_type& u, number time=0.0)
		{ return m_ImpTriangle.assemble_element_A(d, u, time); };

		inline IPlugInReturn assemble_element_M(Triangle* elem, local_vector_type& d, const local_vector_type& u, number time=0.0)
		{ return m_ImpTriangle.assemble_element_M(d, u, time); };

		inline IPlugInReturn assemble_element_f(Triangle* elem, local_vector_type& d, number time=0.0)
		{ return m_ImpTriangle.assemble_element_f(d, time); };

		inline IPlugInReturn finish_element_loop(Triangle* elem)
		{ return m_ImpTriangle.finish_element_loop(); };

	protected:
		ConvectionDiffusionEquation<domain_type, algebra_type, Triangle> m_ImpTriangle;


	public:
		// support assembling on triangles
		inline size_t num_sh(Quadrilateral* elem)
		{ return m_ImpQuadrilateral.num_sh();};

		inline size_t num_sh(Quadrilateral* elem, size_t fct)
		{ return m_ImpQuadrilateral.num_sh(fct);};

		inline IPlugInReturn prepare_element_loop(Quadrilateral* elem)
		{ return m_ImpQuadrilateral.prepare_element_loop(); };

		inline IPlugInReturn prepare_element(Quadrilateral* elem, const local_vector_type& u, const local_index_type& glob_ind)
		{ return m_ImpQuadrilateral.prepare_element(elem, u, glob_ind); };

		inline IPlugInReturn assemble_element_JA(Quadrilateral* elem, local_matrix_type& J, const local_vector_type& u, number time=0.0)
		{ return m_ImpQuadrilateral.assemble_element_JA(J, u, time); };

		inline IPlugInReturn assemble_element_JM(Quadrilateral* elem, local_matrix_type& J, const local_vector_type& u, number time=0.0)
		{ return m_ImpQuadrilateral.assemble_element_JM(J, u, time); };

		inline IPlugInReturn assemble_element_A(Quadrilateral* elem, local_vector_type& d, const local_vector_type& u, number time=0.0)
		{ return m_ImpQuadrilateral.assemble_element_A(d, u, time); };

		inline IPlugInReturn assemble_element_M(Quadrilateral* elem, local_vector_type& d, const local_vector_type& u, number time=0.0)
		{ return m_ImpQuadrilateral.assemble_element_M(d, u, time); };

		inline IPlugInReturn assemble_element_f(Quadrilateral* elem, local_vector_type& d, number time=0.0)
		{ return m_ImpQuadrilateral.assemble_element_f(d, time); };

		inline IPlugInReturn finish_element_loop(Quadrilateral* elem)
		{ return m_ImpQuadrilateral.finish_element_loop(); };

	protected:
		ConvectionDiffusionEquation<domain_type, algebra_type, Quadrilateral> m_ImpQuadrilateral;

};

}

#include "convection_diffusion_assemble_impl.h"

#endif /*__H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__PLUG_IN_DISC__CONVECTION_DIFFUSION_EQUATION__CONVECTION_DIFFUSION_ASSEMBLE__*/
