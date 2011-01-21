/*
 * convection_diffusion.h
 *
 *  Created on: 26.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__CONVECTION_DIFFUSION__FV1__CONVECTION_DIFFUSION__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__CONVECTION_DIFFUSION__FV1__CONVECTION_DIFFUSION__

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

/// \ingroup lib_disc_elem_disc
/// @{

/// Finite Volume Element Discretization for the Convection-Diffusion Equation
/**
 * This class implements the IElemDisc interface to provide element local
 * assemblings for the convection diffusion equation.
 * The Equation has the form
 * \f[
 * 	\partial_t c - \nabla \left( D \nabla c - \vec{v} c \right) + r \cdot c
 * 		= f
 * \f]
 * with
 * <ul>
 * <li>	\f$ c \f$ is the unknown solution
 * <li>	\f$ D \equiv D(\vec{x},t) \f$ is the Diffusion Tensor
 * <li>	\f$ v \equiv \vec{v}(\vec{x},t) \f$ is the Velocity Field
 * <li>	\f$ r \equiv r(\vec{x},t) \f$ is the Reaction Term
 * <li>	\f$ f \equiv f(\vec{x},t) \f$ is a Source Term
 * </ul>
 *
 * \tparam	TFVGeom		Finite Volume Geometry used
 * \tparam	TDomain		Domain
 * \tparam	TAlgebra	Algebra
 */
template<	template <class TElem, int TWorldDim> class TFVGeom,
			typename TDomain,
			typename TAlgebra>
class FVConvectionDiffusionElemDisc
	: public IElemDisc<TAlgebra>
{
	public:
	///	Domain type
		typedef TDomain domain_type;

	///	World dimension
		static const int dim = TDomain::dim;

	///	Position type
		typedef typename TDomain::position_type position_type;

	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Local matrix type
		typedef typename IElemDisc<TAlgebra>::local_matrix_type local_matrix_type;

	///	Local vector type
		typedef typename IElemDisc<TAlgebra>::local_vector_type local_vector_type;

	///	Local index type
		typedef typename IElemDisc<TAlgebra>::local_index_type local_index_type;

	public:
	///	Constructor
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
				register_import(m_MassScale);
			}

	///	set the amount of upwind
	/**
	 * This method sets the amount of upwind to use for the velocity
	 * discretization. A value of 0.0 gives a central scheme, while
	 * 1.0 corresponds to full upwinding.
	 *
	 * \param	amount		Amount of upwind
	 */
		void set_upwind_amount(number amount) {m_upwindAmount = amount;}

	///	sets the domain
		void set_domain(domain_type& domain) {m_pDomain = &domain;}

	///	sets the diffusion tensor
	/**
	 * This method sets the Diffusion tensor used in computations. If no
	 * Tensor is set, a zero value is assumed.
	 */
		void set_diffusion(IPData<MathMatrix<dim, dim>, dim>& user) {m_Diff.set_data(user);}

	///	sets the velocity field
	/**
	 * This method sets the Velocity field. If no field is provided a zero
	 * value is assumed.
	 */
		void set_velocity(IPData<MathVector<dim>, dim>& user) {m_ConvVel.set_data(user);}

	///	sets the reaction
	/**
	 * This method sets the Reaction. A zero value is assumed as default.
	 */
		void set_reaction(IPData<number, dim>& user) {m_Reaction.set_data(user);}

	///	sets the right-hand side
	/**
	 * This method sets the right hand side value. A zero value is assumed as
	 * default.
	 */
		void set_rhs(IPData<number, dim>& user)	{m_Rhs.set_data(user);}


	///	sets mass scale
	/**
	 * This method sets the mass scale value. A value of 1.0 is assumed as
	 * default.
	 */
		void set_mass_scale(IPData<number, dim>& user)	{m_MassScale.set_data(user);}

	public:
	///	number of functions used
		virtual size_t num_fct(){return 1;}

	///	type of trial space for each function used
		virtual LocalShapeFunctionSetID local_shape_function_set_id(size_t loc_fct)
		{
			return LocalShapeFunctionSetID(LocalShapeFunctionSetID::LAGRANGE, 1);
		}

	///	returns if hanging nodes are used
		virtual bool use_hanging() const {return TFVGeom<Edge, dim>::usesHangingNodes;}

	private:
	///	prepares the loop over all elements
	/**
	 * This method prepares the loop over all elements. It resizes the Position
	 * array for the corner coordinates and schedules the local ip positions
	 * at the data imports.
	 */
		template <typename TElem>
		inline bool prepare_element_loop();

	///	prepares the element for assembling
	/**
	 * This methods prepares an element for the assembling. The Positions of
	 * the Element Corners are read and the Finite Volume Geometry is updated.
	 * The global ip positions are scheduled at the data imports.
	 */
		template <typename TElem>
		inline bool prepare_element(TElem* elem,
		                            const local_vector_type& u,
		                            const local_index_type& glob_ind);

	///	finishes the loop over all elements
		template <typename TElem>
		inline bool finish_element_loop();

	///	assembles the local stiffness matrix using a finite volume scheme
		template <typename TElem>
		inline bool assemble_JA(local_matrix_type& J,
		                        const local_vector_type& u, number time=0.0);

	///	assembles the local mass matrix using a finite volume scheme
		template <typename TElem>
		inline bool assemble_JM(local_matrix_type& J,
		                        const local_vector_type& u, number time=0.0);

	///	assembles the stiffness part of the local defect
		template <typename TElem>
		inline bool assemble_A(local_vector_type& d,
		                       const local_vector_type& u, number time=0.0);

	///	assembles the mass part of the local defect
		template <typename TElem>
		inline bool assemble_M(local_vector_type& d,
		                       const local_vector_type& u, number time=0.0);

	///	assembles the local right hand side
		template <typename TElem>
		inline bool assemble_f(local_vector_type& d, number time=0.0);

	protected:
	///	computes the linearized defect w.r.t to the velocity
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
	///	Domain
		TDomain* m_pDomain;

	///	Corner Coordinates
		std::vector<position_type> m_vCornerCoords;

	///	position accessor
		typename TDomain::position_accessor_type m_aaPos;

	///	abbreviation for the local solution
		static const size_t _C_ = 0;

	/// Amount of upwind (1.0 == full upwind, 0.0 == no upwind)
		number m_upwindAmount;

	///	Data import for Diffusion
		DataImport<MathMatrix<dim,dim>, dim, algebra_type> m_Diff;

	///	Data import for the Velocity field
		DataImport<MathVector<dim>, dim, algebra_type > m_ConvVel;

	///	Data import for the reaction term
		DataImport<number, dim, algebra_type> m_Reaction;

	///	Data import for the right-hand side
		DataImport<number, dim, algebra_type> m_Rhs;

	///	Data import for the right-hand side
		DataImport<number, dim, algebra_type> m_MassScale;

	private:
	/// register for 1D
		void register_assemble_functions(Int2Type<1>)
		{
			register_all_assemble_functions<Edge>(ROID_EDGE);
		}

	/// register for 2D
		void register_assemble_functions(Int2Type<2>)
		{
			register_assemble_functions(Int2Type<1>());
			register_all_assemble_functions<Triangle>(ROID_TRIANGLE);
			register_all_assemble_functions<Quadrilateral>(ROID_QUADRILATERAL);
		}

	/// register for 3D
		void register_assemble_functions(Int2Type<3>)
		{
			register_assemble_functions(Int2Type<2>());
			register_all_assemble_functions<Tetrahedron>(ROID_TETRAHEDRON);
			register_all_assemble_functions<Pyramid>(ROID_PYRAMID);
			register_all_assemble_functions<Prism>(ROID_PRISM);
			register_all_assemble_functions<Hexahedron>(ROID_HEXAHEDRON);
		}

	///	register all functions for on element type
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

		//	set computation of linearized defect w.r.t velocity
			m_ConvVel.register_lin_defect_func(id, this, &T::template lin_defect_velocity<TElem>);
		}
};

/// @}

} // end namespace ug

#include "convection_diffusion_impl.h"

#endif /*__H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__CONVECTION_DIFFUSION__FV1__CONVECTION_DIFFUSION__*/
