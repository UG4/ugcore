
#ifndef CONV_DIFF_H
#define CONV_DIFF_H

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lib_grid.h"
#include "lib_algebra/lib_algebra.h"

// library intern headers
#include "lib_discretization/domain_discretization/disc_helper/fvgeom.h"

namespace ug{


template<typename TDomain, typename TElem >
class ConvectionDiffusionEquation {
	public:
		// forward constants and types

		// domain type
		typedef TDomain domain_type;

		// world dimension
		static const int dim = TDomain::dim;

		// position type
		typedef typename TDomain::position_type position_type;

		// local matrix type
		typedef FlexLocalMatrix local_matrix_type;

		// local vector tyoe
		typedef FlexLocalVector local_vector_type;

	protected:
		typedef void (*Diff_Tensor_fct)(MathMatrix<dim,dim>&, const position_type&, number);
		typedef void (*Conv_Vel_fct)(position_type&, const position_type&, number);
		typedef void (*Reaction_fct)(number&, const position_type&, number);
		typedef void (*Rhs_fct)(number&, const position_type&, number);

	public:
		ConvectionDiffusionEquation(TDomain& domain, number upwind_amount, Diff_Tensor_fct diff, Conv_Vel_fct vel, Reaction_fct reac, Rhs_fct rhs)
		: m_domain(domain), m_upwind_amount(upwind_amount), m_Diff_Tensor(diff), m_Conv_Vel(vel), m_Reaction(reac), m_Rhs(rhs)
		{};

	public:

		// number of fundamental functions required for this assembling
		inline uint num_fct(){return 2;}

		// total number of shape functions on elements of type 'TElem'
		inline uint num_sh(){return 2*reference_element_traits<TElem>::num_corners;};

		// local shape function set required for the 'i'-th fundamental function
		inline LocalShapeFunctionSetID trial_space_type(uint i){return LSFS_LAGRANGEP1;}

		// number of shape functions on elements of type 'TElem' for the 'i'-th fundamental function
		inline uint num_sh(uint i){return reference_element_traits<TElem>::num_corners;};

		// prepares the loop. Must be called, before prepare_element can be used
		inline void prepare_element_loop();

		// prepares the element. Must be called before assemble_element_XXX can be used. Must be called after prepare_element_loop
		inline void prepare_element(TElem* elem);

		/**
		 * Assembling routines for local stiffness matrices. Can only be called after prepare_element has been set.
		 *
		 * totol number of degrees of freedom are: total_num_sh = sum_{i=0}^num_fct num_sh(i)
		 * It is required:
		 * total_num_sh == local_vector_type.row_size() == local_matrix_type.size() == local_vector_type.size()
		 *
		 */
		inline void assemble_element_JA(local_matrix_type& J, const local_vector_type& u, number time=0.0);
		inline void assemble_element_JA_bnd(local_matrix_type& J, const local_vector_type& u, number time=0.0){};

		inline void assemble_element_JM(local_matrix_type& J, const local_vector_type& u, number time=0.0);

		inline void assemble_element_A(local_vector_type& d, const local_vector_type& u, number time=0.0);
		inline void assemble_element_A_bnd(local_vector_type& d, const local_vector_type& u, number time=0.0){};

		inline void assemble_element_M(local_vector_type& d, const local_vector_type& u, number time=0.0);

		inline void assemble_element_f(local_vector_type& d, number time=0.0);
		inline void assemble_element_f_bnd(local_vector_type& d, number time=0.0){};

		inline void finish_element_loop();

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
};

/*
template <typename TDomain, typename TElem, template <typename TDomain, typename TElem> class ElemDisc>
class ElemAssembleFct{
	public:
		inline void prepare_element_discretization();
		inline void prepare_element(TElem* elem);
		inline void assemble_element_JA(TElem* elem, number mat_values[], number u_values[], const uint num_dofs, number time=0.0);
		inline void assemble_element_JM(TElem* elem, number mat_values[], number u_values[], const uint num_dofs, number time=0.0);
		inline void assemble_element_A(TElem* elem, number def_values[], number u_values[], const uint num_dofs, number time=0.0);
		inline void assemble_element_M(TElem* elem, number def_values[], number u_values[], const uint num_dofs, number time=0.0);
		inline void assemble_element_f(TElem* elem, number def_values[], const uint num_dofs, number time=0.0);
		inline void prepare_element_discretization();

	private:
		ElemDisc<TDomain, TElem> m_Disc;
};

template <typename TDomain>
class BndAssembleFct{

};

template <typename TDomain, template <typename TDomain, typename TElem> class ElemDisc>
class PlugIn{

	ElemAssembleFct<TDomain, Triangle> m_ElemDiscTriangle;


};
*/


/*
template <typename TDomain>
class ConvectionDiffusionPlugIn{

	ConvectionDiffusionPlugIn(TDomain& domain, number upwind_amout, Diff_Tensor_fct diff, Conv_Vel_fct vel, Reaction_fct reac, Rhs_fct rhs) :
		m_ImpTriangle(domain, upwind_amount, diff, vel, reac, rhs),
		m_ImpQuadrilateral(domain, upwind_amount, diff, vel, reac, rhs) {};


	public:
		// support assembling on triangles
		inline void prepare_element_discretization();
		inline void prepare_element(Triangle* elem);
		inline void assemble_element_JA(Triangle* elem, number mat_values[], number u_values[], const uint num_dofs, number time=0.0);
		inline void assemble_element_JM(Triangle* elem, number mat_values[], number u_values[], const uint num_dofs, number time=0.0);
		inline void assemble_element_A(Triangle* elem, number def_values[], number u_values[], const uint num_dofs, number time=0.0);
		inline void assemble_element_M(Triangle* elem, number def_values[], number u_values[], const uint num_dofs, number time=0.0);
		inline void assemble_element_f(Triangle* elem, number def_values[], const uint num_dofs, number time=0.0);
		inline void prepare_element_discretization();
	protected:
		ConvectionDiffusionEquation<TDomain, Triangle> m_ImpTriangle;


	public:
		// support assembling on quadrilaterals
		inline void prepare_element(Quadrilateral* elem);
		inline void assemble_element_JA(Quadrilateral* elem, number mat_values[], number u_values[], const uint num_dofs, number time=0.0);
		inline void assemble_element_JM(Quadrilateral* elem, number mat_values[], number u_values[], const uint num_dofs, number time=0.0);
		inline void assemble_element_A(Quadrilateral* elem, number def_values[], number u_values[], const uint num_dofs, number time=0.0);
		inline void assemble_element_M(Quadrilateral* elem, number def_values[], number u_values[], const uint num_dofs, number time=0.0);
		inline void assemble_element_f(Quadrilateral* elem, number def_values[], const uint num_dofs, number time=0.0);
	protected:
		ConvectionDiffusionEquation<TDomain, Quadrilateral> m_ImpQuadrilateral;
};*/

}

#include "convection_diffusion_assemble_impl.h"

#endif
