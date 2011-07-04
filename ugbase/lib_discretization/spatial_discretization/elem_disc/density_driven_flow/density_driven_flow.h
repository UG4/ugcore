/*
 * density_driven_flow.h
 *
 *  Created on: 26.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__DENSITY_DRIVEN_FLOW__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__DENSITY_DRIVEN_FLOW__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_discretization/spatial_discretization/elem_disc/elem_disc_interface.h"

#include "lib_discretization/spatial_discretization/ip_data/ip_data.h"
#include "lib_discretization/spatial_discretization/ip_data/data_export.h"
#include "lib_discretization/spatial_discretization/ip_data/data_import.h"

#include "lib_discretization/spatial_discretization/disc_util/conv_shape_interface.h"
#include "lib_discretization/spatial_discretization/disc_util/consistent_gravity.h"

namespace ug{

/// \ingroup lib_disc_elem_disc
/// @{

/// Finite Volume Element Discretization for Density Driven Flow
/**
 * This class implements the IElemDisc interface for the density driven
 * flow equations. It is a system of two coupled equations of the form.
 * \f{align*}
 * 	\partial_t (\phi \rho) + \nabla(\rho \vec{q}) &= 0\\
 * 	\partial_t (\phi \rho c) - \nabla(\rho D \nabla c -\rho c \vec{q}) &= 0\\
 * 	\vec{q} &:= - \frac{K}{\mu}(\nabla p - \rho \vec{g})
 * \f}
 * with
 * <ul>
 * <li> \f$ c \f$		unknown brine mass fraction
 * <li>	\f$ p \f$ 		unknown pressure
 * </ul>
 * and given data
 * <ul>
 * <li>	\f$ \phi(\vec{x}, t) \f$	Porosity
 * <li>	\f$ K(\vec{x}, t) \f$		Permeability tensor
 * <li>	\f$ \vec{g} \f$	constant Gravity field
 * <li>	\f$ D(\vec{q}, \vec{x}, t) := \phi(\vec{x}, t) D_d(\vec{x}, t) + D_m(\vec{q})  \f$
 * 			Diffusion-Dispersion tensor
 * <li>	\f$ D_m(\vec{q})  \f$	mechanical Dispersion tensor
 * <li>	\f$ D_d(\vec{x}, t) \f$	molecular Diffusion tensor (prop. to Tortuosity)
 * <li>	\f$ \mu(c) \f$		Viscosity
 * <li>	\f$ \rho(c)\f$	Density
 * </ul>
 * as well as the abbreviation
 * <ul>
 * <li>	\f$ \vec{q} \f$	Darcy velocity
 * </ul>
 *
 * The first equation is usually referred to as the flow equation. The second
 * equation is known as the transport equation.
 *
 * \tparam	TFVGeom		Finite Volume Geometry
 * \tparam	TDomain		Domain
 * \tparam	TAlgebra	Algebra
 */
template<typename TDomain>
class DensityDrivenFlowElemDisc
	: public IDomainElemDisc<TDomain>
{
	private:
	///	Base class type
		typedef IDomainElemDisc<TDomain> base_type;

	///	own type
		typedef DensityDrivenFlowElemDisc<TDomain> this_type;

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

	public:
	///	Constructor
		DensityDrivenFlowElemDisc();

	///	sets usage of consistent gravity
		void set_consistent_gravity(bool bUse)
		{
			m_bConsGravity = bUse;
			register_all_fv1_funcs();
		}

	///	sets usage of boussinesq approximation for transport equation
		void set_boussinesq_transport(bool bUse) {m_BoussinesqTransport = bUse;}

	///	sets usage of boussinesq approcimation for flow equation
		void set_boussinesq_flow(bool bUse) {m_BoussinesqFlow = bUse;}

	///	sets the type of upwind
	/**
	 * This method sets the procedure that compute the upwinded flux.
	 */
		void set_upwind(IConvectionShapes<dim>& shape)
		{
			m_pUpwind = &shape;
		}

	///	sets the porosity
	/**
	 * This method sets the Porosity. (Dimensionless)
	 */
		void set_porosity(IPData<number, dim>& user)
		{
			m_imPorosityScv.set_data(user);
			m_imPorosityScvf.set_data(user);
		}

	///	sets the gravity vector
	/**
	 * This method sets the Gravity. (Unit is \f$ \frac{m}{s^2} \f$)
	 */
		void set_gravity(IPData<MathVector<dim>, dim>& user)
		{
			m_imConstGravity.set_data(user);
		}

	///	sets the molecular diffusion tensor
	/**
	 * This method sets the molecular Diffusion tensor.
	 */
		void set_molecular_diffusion(IPData<MathMatrix<dim, dim>, dim>& user)
		{
			m_imMolDiffusionScvf.set_data(user);
		}

	///	sets the permeability tensor
	/**
	 * This method sets the Permeability tensor.
	 */
		void set_permeability(IPData<MathMatrix<dim, dim>, dim>& user)
		{
			m_imPermeabilityScvf.set_data(user);
		}

	///	sets the viscosity tensor
	/**
	 * This method sets the Viscosity.
	 */
		void set_viscosity(IPData<number, dim>& user)
		{
			m_imViscosityScvf.set_data(user);
		}

	///	set density
		void set_density(IPData<number,dim>& data)
		{
		//	remove old data
			IIPData* oldData = m_imDensityScv.get_data();
			if (oldData != NULL)
				m_exDarcyVel.remove_needed_data(*oldData);
			oldData = m_imDensityScvf.get_data();
			if (oldData != NULL)
				m_exDarcyVel.remove_needed_data(*oldData);

		//	connect to import
			m_imDensityScv.set_data(data);
			m_imDensityScvf.set_data(data);

		//	darcy velocity depends on density
			m_exDarcyVel.add_needed_data(data);
		}

	public:
	///	number of functions used
		virtual size_t num_fct() {return 2;}

	///	trial space for functions
		virtual LocalShapeFunctionSetID local_shape_function_set_id(size_t loc_fct)
		{
			return LocalShapeFunctionSetID(LocalShapeFunctionSetID::LAGRANGE, 1);
		}

	///	switches between non-regular and regular grids
		virtual bool treat_non_regular_grid(bool bNonRegular)
		{
		//	switch, which assemble functions to use.
			if(bNonRegular)
			{
				UG_LOG("ERROR in 'DensityDrivenFlowElemDisc::treat_non_regular_grid':"
						" Non-regular grid not implemented.\n");
				return false;
			}

		//	this disc supports regular grids
			return true;
		}

	private:
	///	prepares the discretization for time dependent discretization
	/**
	 * This function prepares the discretization for time-dependent problems.
	 * It sets the time in the imports.
	 *
	 * \param[in]	time	new time point
	 * \returns 	true	indicates, that old values are needed
	 */
		virtual bool time_point_changed(number time);

		template <typename TElem>
		bool prepare_element_loop();

		template <typename TElem>
		bool prepare_element(TElem* elem, const local_vector_type& u,
		                            const local_index_type& glob_ind);

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
	///	strategy to compute the upwind shapes
		IConvectionShapes<dim>* m_pUpwind;

	///	flag if using Consistent Gravity
		bool m_bConsGravity;

	///	flag if using boussinesq transport
		bool m_BoussinesqTransport;

	///	flag if using boussinesq flow
		bool m_BoussinesqFlow;

	///	Corner Coordinates
		const position_type* m_vCornerCoords;

	/// constant Gravity, read in once
		MathVector<dim> m_Gravity;

	///	abbreviation for local function: brine mass fraction
		static const size_t _C_ = 0;

	///	abbreviation for local function: pressure
		static const size_t _P_ = 1;

	private:
	///	Data import for brine mass fraction at scvf ips
		DataImport<number, dim> m_imBrineScvf;
		DataImport<MathVector<dim>, dim> m_imBrineGradScvf;

	///	Data import for pressure gradient at scvf ips
		DataImport<MathVector<dim>, dim> m_imPressureGradScvf;

	///	Data import for the reaction term
		DataImport<number, dim> m_imPorosityScv;
		DataImport<number, dim> m_imPorosityScvf;

	///	Data import for gravity (must be constant in current implementation)
		DataImport<MathVector<dim>, dim> m_imConstGravity;

	///	Data import for permeability tensor
		DataImport<MathMatrix<dim, dim>, dim> m_imPermeabilityScvf;

	///	Data import for molecular diffusion tensor
		DataImport<MathMatrix<dim, dim>, dim> m_imMolDiffusionScvf;

	///	Data import for viscosity
		DataImport<number, dim> m_imViscosityScvf;

	///	Data import for density
		DataImport<number, dim> m_imDensityScv;
		DataImport<number, dim> m_imDensityScvf;

	///	Data import for Darcy Velocity
		DataImport<MathVector<dim>, dim> m_imDarcyVelScvf;

	public:
	///	returns the export of the darcy velocity
		IPData<MathVector<dim>, dim>& get_darcy_velocity() {return m_exDarcyVel;}

	///	returns the export of brine mass fracture
		IPData<number, dim>& get_brine() {return m_exBrine;}

	///	returns the export of brine mass fracture
		IPData<MathVector<dim>, dim>& get_brine_grad() {return m_exBrineGrad;}

	///	returns the export of brine mass fracture
		IPData<MathVector<dim>, dim>& get_pressure_grad() {return m_exPressureGrad;}

	protected:
	///	compute darcy velocity at one ip
		void compute_darcy_velocity_ip_std(MathVector<dim>& DarcyVel,
		                                   const MathMatrix<dim, dim>& Permeability,
		                                   number Viscosity,
		                                   const MathVector<dim>& DensityTimesGravity,
		                                   const MathVector<dim>& PressureGrad,
		                                   bool compDeriv,
		                                   MathVector<dim>* DarcyVel_c,
		                                   MathVector<dim>* DarcyVel_p,
		                                   const MathVector<dim>* DensityTimesGravity_c,
		                                   const number* Viscosity_c,
		                                   const MathVector<dim>* PressureGrad_p,
		                                   size_t numSh);

	///	computes the darcy velocity using consistent gravity
		template <typename TElem>
		bool compute_darcy_export_std(const local_vector_type& u, bool compDeriv);

	///	computes the darcy velocity using consistent gravity
		template <typename TElem>
		bool compute_darcy_export_cons_grav(const local_vector_type& u, bool compDeriv);

	///	computes the value of the brine mass fraction
		template <typename TElem>
		bool compute_brine_export(const local_vector_type& u, bool compDeriv);

	///	computes the value of the gradient of the brine mass fraction
		template <typename TElem>
		bool compute_brine_grad_export(const local_vector_type& u, bool compDeriv);

	///	computes the value of the gradient of the pressure
		template <typename TElem>
		bool compute_pressure_grad_export(const local_vector_type& u, bool compDeriv);

	///	Export for the Darcy velocity
		DataExport<MathVector<dim>, dim> m_exDarcyVel;

	///	Export for the brine mass fraction
		DataExport<number, dim> m_exBrine;

	///	Export for the gradient of brine mass fraction
		DataExport<MathVector<dim>, dim> m_exBrineGrad;

	///	Export for the gradient of brine mass fraction
		DataExport<MathVector<dim>, dim> m_exPressureGrad;

	private:
		void register_all_fv1_funcs();

		struct RegisterFV1 {
				RegisterFV1(this_type* pThis) : m_pThis(pThis){}
				this_type* m_pThis;
				template< typename TElem > void operator()(TElem&)
				{m_pThis->register_fv1_func<TElem>();}
		};

		template <typename TElem>
		void register_fv1_func();

};

///@}

} // end namespace ug

#endif /*__H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__DENSITY_DRIVEN_FLOW__*/
