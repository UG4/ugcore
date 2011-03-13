/*
 * navier_stokes.h
 *
 *  Created on: 20.09.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__NAVIER_STOKES__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__NAVIER_STOKES__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_discretization/spatial_discretization/elem_disc/elem_disc_interface.h"

#include "lib_discretization/spatial_discretization/disc_helper/geometry_provider.h"
#include "lib_discretization/spatial_discretization/disc_helper/finite_volume_geometry.h"
#include "lib_discretization/common/local_algebra.h"
#include "lib_discretization/spatial_discretization/ip_data/user_data.h"

#include "new_upwind.h"
#include "new_stabilization.h"

namespace ug{


/// \ingroup lib_disc_elem_disc
/// @{

/// Finite Volume Element Discretization for the incompressible Navier-Stokes Equation
/**
 * This class implements the IElemDisc interface to provide element local
 * assemblings for the incompressible Navier-Stokes equation.
 *
 * The unknowns of the equation are:
 * <ul>
 * <li> \f$ \vec{u} \f$ velocity
 * <li> \f$ p \f$ pressure.
 * </ul>
 *
 * The equation takes the form
 * \f{align*}
 * 	\frac{\partial \rho \vec{u}}{\partial t}
 * 	- \nabla \left( \rho \nu (\nabla \vec{u} + (\nabla \vec{u})^T) \right)
 * 	+ \nabla \cdot \left( \rho \vec{u} \vec{u}^T \right)
 *  + \nabla p
 *  &= \vec{f}\\
 *
 *  \nabla \cdot (\rho \vec{u}) &= 0
 * \f}
 *
 * with
 * <ul>
 * <li>	\f$ \rho \f$ is the constant density
 * <li>	\f$ \nu \f$ is the kinematic viscosity (temporarily constant, implement)
 * <li>	\f$ \vec{f} \equiv f(\vec{x},t) \f$ is a Source Term (not implemented yet)
 * </ul>
 *
 * The first equation models the conservation of momentum and is therefore
 * referred to as the <b>Momentum equation</b>. The second equation models the
 * conservation of mass and is known as the <b>Mass continuity equation</b> or
 * simply <b>Continuity equation</b>.
 *
 * In component-wise notation, the equation reads
 *
 * \f{align*}
 * \frac{\partial \rho u_{d_1}}{\partial t}
 * - \sum_{d_2 = 1}^{\text{dim}} \frac{\partial}{\partial x_{d_2}}
 * 		\left( \rho \nu \left( \frac{\partial u_{d_1}}{\partial x_{d_2}}
 * 			+ \frac{\partial u_{d_2}}{\partial x_{d_1}} \right)\right)
 * + \sum_{d_2 = 1}^{\text{dim}}  \frac{\partial}{\partial x_{d_2}}
 * 		\left( \rho u_{d_2} u_{d_1} \right)
 * + \frac{\partial p}{\partial x_{d_1}}
 * &= f_{d_1} \qquad \forall \quad d_1 = 1, \dots, \text{dim}\\
 *
 * \sum_{d = 1}^{\text{dim}}  \frac{\partial \rho u_d}{\partial x_{d}} &= 0
 * \f}
 *
 * Since the density is assumed to be constant, it set to \f$\rho \equiv 1 \f$.
 *
 * The finite volume element discretization uses a dual mesh consisting of
 * Control Volumes \f$\mathcal{B}_h\f$, that cover the domain. For each
 * control volume \f$B \in \mathcal{B}_h\f$ the following equation is solved
 *
 * \f{align*}
 * 	\frac{\partial}{\partial t} \int_{B} \begin{pmatrix} \vec{u} \\ 0 \end{pmatrix} dV
 *  + \int_{\partial B}
 *  	\begin{pmatrix}
 *  		- \rho \nu \left(\nabla \vec{u} + (\nabla \vec{u})^T \right) \vec{n}
 *  		+ \vec{u} \vec{u}^T \vec{n} + p \vec{n} \\
 * 		 	\vec{u} \vec{n}
 * 		\end{pmatrix} dS
 * = \int_B \begin{pmatrix} \vec{f} \\ 0 \end{pmatrix} dV
 * \f}
 *
 * where \f$ \vec{n}\f$ is the outer normal on the boundary of the control volume
 * and \f$\int_B \cdot \; dV \f$ and \f$\int_{\partial B} \cdot \; dS \f$ are the
 * integration over the control volume and the integration over the boundary of
 * the control volume. By \$f B_{sh} \$f is denoted the control volume associated
 * to a shape function (i.e. vertices, since we use P1 Lagrange ansatz functions),
 * T usually denotes the Element.
 *
 * In order to number the local unknowns, we use the following notation.
 * <ul>
 * <li> \f$ sh = 0, \dots, n_{sh} - 1\f$ loop the \f$n_{sh}\f$ shape functions. For
 * 		this implementation these are the corners of the element
 * <li> \f$ d = 0, \dots, \text{dim} - 1 \f$ are the velocity component of each
 * 		corner
 * <li> _P_ := dim indicates the pressure component
 * </ul>
 * The access to local unknowns of local vectors and matrices is now given by
 * \f{align*}
 * 	&\mathbf{d}(d_1, sh), \mathbf{d}(\text{\_P\_}, sh),\\
 *  &\mathcal{J}(d_1, sh_1, d_2, sh_2), \mathcal{J}(d_1, sh_1, \text{\_P\_}, sh_2),
 * \f}
 *
 * \tparam	TDomain		Domain
 * \tparam	TAlgebra	Algebra
 */
template<	typename TDomain,
			typename TAlgebra>
class FVNavierStokesElemDisc
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

	public:
	///	Constructor (setting default values)
		FVNavierStokesElemDisc()
			: m_bExactJacobian(true), m_pStab(NULL)
		{
		//	set default options
			set_PAC(false);
			set_peclet_blend(false);

		//	register imports
			register_import(m_imSource);
			register_import(m_imKinViscosity);

		//	register assemble functions
			register_assemble_functions(Int2Type<dim>());
		}

	///	sets the kinematic viscosity
	/**
	 * This method sets the kinematic viscosity value.
	 *
	 * \param[in]	nu		kinematic Viscosity
	 */
		void set_kinematic_viscosity(IPData<number, dim>& data)
			{m_imKinViscosity.set_data(data);}

	///	sets the source function
	/**
	 * This method sets the source value. A zero value is assumed as default.
	 * \param[in]	data		source data
	 */
		void set_source(IPData<MathVector<dim>, dim>& data)
			{m_imSource.set_data(data);}

	///	sets the stabilization used to compute the stabilized velocities
	/**
	 * \param[in]	stab		Stabilization
	 */
		void set_stabilization(INavierStokesStabilization<dim, algebra_type>& stab)
			{m_pStab = & stab;}

	///	sets if Physical Advection Correction is used
        void set_PAC(bool bPAC){m_bPAC = bPAC;}

    ///	sets if peclet blending is used in momentum equation
        void set_peclet_blend(bool pecletBlend) {m_bPecletBlend = pecletBlend;}

    ///	sets if the exact jacobian is computed (fixpoint approximation else)
        void set_exact_jacobian(bool bExactJacobian) {m_bExactJacobian = bExactJacobian;}

	public:
	///	returns number of functions needed for this elem disc
		virtual size_t num_fct(){return dim+1;}

	///	returns shape function set extected for each function
		virtual LocalShapeFunctionSetID local_shape_function_set_id(size_t loc_fct)
		{
			return LocalShapeFunctionSetID(LocalShapeFunctionSetID::LAGRANGE, 1);
		}

	///	switches between non-regular and regular grids
	/**
	 * \param[in]	bNonRegular		flag if non-regular grid needed.
	 */
		virtual bool treat_non_regular_grid(bool bNonRegular)
		{
		//	switch, which assemble functions to use.
			if(bNonRegular)
			{
				UG_LOG("ERROR in 'FVNavierStokesElemDisc::treat_non_regular_grid':"
						" Non-regular grid not implemented.\n");
				return false;
			}

		//	this disc supports regular grids
			return true;
		}

	private:
	///	prepares the element loop
	/**
	 * This function is used to prepare the element loop. It gets called, before
	 * the loop over all elements starts. Here, general checks are performed:
	 * - Is the correct finite volume geometry used
	 * - Has the stabilization been specified
	 *
	 * In addition, local coordinates of evaluation points are requested by
	 * the DataImports in case of element-fixed points.
	 */
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		inline bool prepare_element_loop();

	///	prepares the element for evaluation
	/**
	 * This function prepare a specific element for assembling. This function
	 * is always called before any of the assemble_... functions is performed.
	 * Here, the Finite Volume Geometry is prepared and global positions of
	 * the evaluation points of the DataImports are requested (and local
	 * positions in case of hanging node assembling)
	 *
	 * \param[in]	elem		element to prepare the data for
	 * \param[in]	u			current local solution vector
	 * \param[in]	glob_ind	global indices of the local vector components
	 */
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		inline bool prepare_element(TElem* elem, const local_vector_type& u, const local_index_type& glob_ind);

	///	finishes the element loop
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		inline bool finish_element_loop();

	///	adds the stiffness part to the local jacobian
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		inline bool assemble_JA(local_matrix_type& J, const local_vector_type& u, number time=0.0);

	///	adds the stiffness part to the local defect
	/**
	 * This function adds the contribution of one element of the stiffness part to
	 * the local defect.
	 *
	 * The local solution \f$ \mathbf{u} \f$ passes the p1 conform lagrange
	 * dofs of the element, i.e. the unknowns in the corners of the element.
	 *
	 * Define the functional matrix of the velocity at some integration point by
	 * \f{align*}
	 * 	\nabla \vec{u}|_{ip} := \left(\left.\frac{\partial u_i}{\partial x_j}
	 * 							\right|_{ip}\right)_{i,j = 1,\dots,dim}
	 * \f}
	 *
	 * and the diffusive flux to
	 * \f{align*}
	 * 	\vec{f}^{\text{diff}}|_{ip} := - \nu \left( \nabla \vec{u}|_{ip} +
	 * 			(\nabla \vec{u}|_{ip})^T) \right) \cdot \vec{n}_{ip}
	 * \f}
	 *
	 * Suppose that a procedure producing some upwind velocity \f$ \vec{u}^{up}\f$
	 * is available for each ip. Then, we define the convective flux to
	 * \f{align*}
	 * 	\vec{f}^{\text{conv}}|_{ip} := \left( \vec{u}^{up}|_{ip}  \vec{u}^{up}|_{ip}^T
	 * 									\right) \cdot \vec{n}_{ip}
	 * \f}
	 *
	 * Define the pressure at the ip to be interpolated as
	 * \f{align*}
	 * 	p|_{ip} := \sum_{sh=0}^{n_{sh}-1} \mathbf{u}(\text{\_P\_}, sh) \phi_{sh}
	 * \f}
	 *
	 * The contribution added for the <b>momentum equation</b> is:
	 * \f{align*}
	 *  \forall d=0,\dots, \text{dim}-1, \quad
	 *  \forall sh=0,\dots, n_{sh}-1: \qquad
	 *	\mathbf{d}(d, sh) &\text{ +=}
	 *		\int_{\partial B_{sh} \cap T} \vec{f}^{\text{diff}}_d
	 *			+ \vec{f}^{\text{conv}}_d
	 *			+ p \vec{n}_d \;dS
	 *
	 *		\approx \sum_{i=1}^{n_{\partial B_{sh}}}
	 *			 	|\partial B^i_{sh} \cap T|
	 *			 	\left( \left.\vec{f}^{\text{diff}}_d\right|_{ip}
	 *			 	+ \left.\vec{f}^{\text{conv}}_d\right|_{ip}
	 *			 	+ p|_{ip} \vec{n}|_{ip} \right)\\
	 * \f}
	 *
	 * Assume, that some stabilized velocity \f$ \vec{u}^{stab}\f$ is available
	 * for each ip.
	 * Then, the contribution added for the <b>continuity equation</b> is:
	 * \f{align*}
	 *  \forall sh=0,\dots, n_{sh}-1: \qquad
	 *	\mathbf{d}(\text{\_P\_}, sh) &\text{ +=}
	 *		\int_{\partial B_{sh} \cap T} \vec{u}^{stab} \cdot \vec{n} \;dS
	 *
	 *		\approx \sum_{i=1}^{n_{\partial B_{sh}}}
	 *			 	|\partial B^i_{sh} \cap T| \; \vec{u}^{stab}|_{ip} \cdot \vec{n}|_{ip}\\
	 * \f}
	 *
	 * \param[in,out]	d		local defect with values added on output
	 * \param[in]		u		local vector of current solution
	 * \param[in]		time	current time step
	 *
	 * \tparam	TElem 	Element type
	 * \tparam	TFVGeom	Finite Volume Geometry
	 */
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		inline bool assemble_A(local_vector_type& d, const local_vector_type& u, number time=0.0);

	///	adds the mass part to the local jacobian
	/**
	 * This function adds the contribution of one element of the mass part to
	 * the local jacobian.
	 *
	 * The local solution \f$ \mathbf{u} \f$ passes the p1 conform lagrange
	 * dofs of the element, i.e. the unknowns in the corners of the element.
	 *
	 * The contribution added is:
	 * \f{align*}
	 *  \forall d=0,\dots, \text{dim}-1, \quad
	 *  \forall sh=0,\dots, n_{sh}-1: \qquad
	 *	\mathcal{J}(d, sh, d, sh) &\text{ +=} \int_{B_{sh} \cap T} dV
	 *							  = |B_{sh} \cap T|\\
	 * \f}
	 *
	 * \param[in,out]	J		local jacobian with values added on output
	 * \param[in]		u		local vector of current solution
	 * \param[in]		time	current time step
	 *
	 * \tparam	TElem 	Element type
	 * \tparam	TFVGeom Finite Volume Geometry
	 */
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		inline bool assemble_JM(local_matrix_type& J, const local_vector_type& u, number time=0.0);

	///	adds the mass part to the local defect
	/**
	 * This function adds the contribution of one element of the mass part to
	 * the local defect.
	 *
	 * The local solution \f$ \mathbf{u} \f$ passes the p1 conform lagrange
	 * dofs of the element, i.e. the unknowns in the corners of the element.
	 *
	 * The contribution added is:
	 * \f{align*}
	 *  \forall d=0,\dots, \text{dim}-1, \quad
	 *  \forall sh=0,\dots, n_{sh}-1: \qquad
	 *	\mathbf{d}(d, sh) &\text{ +=} \int_{B_{sh} \cap T} \vec{u}_d dV
	 *					  \approx |B_{sh} \cap T| \mathbf{u}(d, sh)\\
	 *
	 *	\mathbf{d}(\text{\_P\_, sh}) &\text{ +=} 0
	 * \f}
	 *
	 * \param[in,out]	d		local defect with values added on output
	 * \param[in]		u		local vector of current solution
	 * \param[in]		time	current time step
	 *
	 * \tparam	TElem 	Element type
	 * \tparam	TFVGeom	Finite Volume Geometry
	 */
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		inline bool assemble_M(local_vector_type& d, const local_vector_type& u, number time=0.0);

	///	adds the source part to the local defect
	/**
	 * This function adds the contribution of one element of the source part to
	 * the local defect, if a source function \f$ \vec{f} \f$ has been specified.
	 * Otherwise a zero source function is assumed.
	 *
	 * The contribution added is:
	 * \f{align*}
	 *  \forall d=0,\dots, \text{dim}-1, \quad
	 *  \forall sh=0,\dots, n_{sh}-1: \qquad
	 *	\mathbf{d}(d, sh) &\text{ +=} \int_{B_{sh} \cap T} \vec{f}_d \; dV
	 *					  \approx |B_{sh} \cap T| \; \vec{f}_d|_{ip}
	 * \f}
	 *
	 * \param[in,out]	d		local defect with values added on output
	 * \param[in]		time	current time step
	 *
	 * \tparam	TElem 	Element type
	 * \tparam	TFVGeom	Finite Volume Geometry
	 */
		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		inline bool assemble_f(local_vector_type& d, number time=0.0);

	///	computes the pecled blended Upwind veloctity
	/**
	 * This function compute the peclet blended velocity. First, the
	 * standard interpolated (no upwind) velocity is computed at the integration
	 * point of the scvf:
	 * \f{align*}
	 * 	\vec{u}^{\text{std}}_d|_{ip}
	 * 			:= \sum_{sh=0}^{n_{sh}-1} \mathbf{u}(d, sh) \; \phi_{sh}
	 * \f}
	 *
	 * Defining the Peclet Number
	 * \f{align*}
	 * 		Pe := \frac{\vec{u}^{\text{std}}|_{ip} \cdot \vec{n}|_{ip} \; L}{\nu}
	 * \f}
	 * with the Viscosity \f$ \nu \f$ and the diffusion length \f$ L \f$, a
	 * weighting factor is computed:
	 * \f{align*}
	 * 	\omega := \frac{Pe^2}{5 + Pe^2}
	 * \f}
	 *
	 * and the blended velocity is computed to
	 * \f{align*}
	 * 	\vec{u}^{\text{blend}}|_{ip} := \omega \vec{u}^{\text{up}}
	 * 									+ (1-\omega) \vec{u}^{\text{std}}
	 * \f}
	 *
	 * \param[in,out]	UpwindVel		upwind vel on entry, blended vel on exit
	 * \param[in]		scvf			SubControlVolumeFace of the Velocity
	 * \param[in]		u				local solution vector
	 * \param[in]		kinVisco		kinematic Viscosity at scvf
	 * \return			\f$\omega\f$ 	weighting factor
	 */
		template <typename SCVF>
		inline number peclet_blend(MathVector<dim>& UpwindVel, const SCVF& scvf,
		                        const local_vector_type& u, number kinVisco);

	private:
	///	flag if using PAC
		int m_bPAC;

	///	flag if using Peclet Blending
		int m_bPecletBlend;

	///	flag if computing exact jacobian
		bool m_bExactJacobian;

	///	Data import for source
		DataImport<MathVector<dim>, dim, algebra_type> m_imSource;

	///	Data import for kinematic viscosity
		DataImport<number, dim, algebra_type> m_imKinViscosity;

	///	Stabilization for velocity in continuity equation
		INavierStokesStabilization<dim, algebra_type>* m_pStab;

	/// position access
		const position_type* m_vCornerCoords;

	/// abbreviation for pressure
		static const size_t _P_ = dim;

	private:
		///////////////////////////////////////
		// registering for reference elements
		///////////////////////////////////////
	/// register for 1D
		void register_assemble_functions(Int2Type<1>)
		{
			register_all_assemble_functions<Edge>(ROID_EDGE);
		}

	/// register for 2D
		void register_assemble_functions(Int2Type<2>)
		{
			register_all_assemble_functions<Triangle>(ROID_TRIANGLE);
			register_all_assemble_functions<Quadrilateral>(ROID_QUADRILATERAL);
		}

	/// register for 3D
		void register_assemble_functions(Int2Type<3>)
		{
			register_all_assemble_functions<Tetrahedron>(ROID_TETRAHEDRON);
			register_all_assemble_functions<Pyramid>(ROID_PYRAMID);
			register_all_assemble_functions<Prism>(ROID_PRISM);
			register_all_assemble_functions<Hexahedron>(ROID_HEXAHEDRON);
		}

	/// help function to register assemble functions for one element
		template <typename TElem>
		void register_all_assemble_functions(int id)
		{
			register_prepare_element_loop_function(	id, &FVNavierStokesElemDisc::template prepare_element_loop<TElem, FV1Geometry>);
			register_prepare_element_function(		id, &FVNavierStokesElemDisc::template prepare_element<TElem, FV1Geometry>);
			register_finish_element_loop_function(	id, &FVNavierStokesElemDisc::template finish_element_loop<TElem, FV1Geometry>);
			register_assemble_JA_function(			id, &FVNavierStokesElemDisc::template assemble_JA<TElem, FV1Geometry>);
			register_assemble_JM_function(			id, &FVNavierStokesElemDisc::template assemble_JM<TElem, FV1Geometry>);
			register_assemble_A_function(			id, &FVNavierStokesElemDisc::template assemble_A<TElem, FV1Geometry>);
			register_assemble_M_function(			id, &FVNavierStokesElemDisc::template assemble_M<TElem, FV1Geometry>);
			register_assemble_f_function(			id, &FVNavierStokesElemDisc::template assemble_f<TElem, FV1Geometry>);
		}
};

/// @}

} // end namespace ug

// include implementation
#include "navier_stokes_impl.h"

#endif /*__H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__NAVIER_STOKES__*/
