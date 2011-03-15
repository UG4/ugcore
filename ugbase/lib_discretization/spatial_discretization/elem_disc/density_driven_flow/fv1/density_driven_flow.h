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
#include "lib_discretization/spatial_discretization/disc_helper/geometry_provider.h"
#include "lib_discretization/spatial_discretization/disc_helper/fvgeom.h"
#include "lib_discretization/common/local_algebra.h"
#include "lib_discretization/spatial_discretization/ip_data/ip_data.h"
#include "lib_discretization/spatial_discretization/ip_data/data_export.h"
#include "lib_discretization/spatial_discretization/ip_data/data_import.h"
#include "consistent_gravity.h"
#include "convection_shape.h"

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
template<template <	class TElem, int TWorldDim> class TFVGeom,
					typename TDomain,
					typename TAlgebra>
class DensityDrivenFlowElemDisc
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
	///	Constructor
		DensityDrivenFlowElemDisc() :
			m_Upwind(NO_UPWIND), m_bConsGravity(true),
			m_BoussinesqTransport(true), m_BoussinesqFlow(true),
			m_imBrineScvf(false), m_imBrineGradScvf(false),
			m_imPressureGradScvf(false),
			m_imDensityScv(false), m_imDensityScvf(false),
			m_imDarcyVelScvf(false)
		{
		//	register assemble functions
			register_assemble_functions(Int2Type<dim>());
			register_export_functions(Int2Type<dim>());

		//	register export
			m_exDarcyVel.add_needed_data(m_exBrine);
			register_export(m_exDarcyVel);
			register_export(m_exBrine);
			register_export(m_exBrineGrad);
			register_export(m_exPressureGrad);

		//	register import
			register_import(m_imBrineScvf);
			register_import(m_imBrineGradScvf);
			register_import(m_imPressureGradScvf);
			register_import(m_imPorosityScvf);
			register_import(m_imPorosityScv);
			register_import(m_imPermeabilityScvf);
			register_import(m_imMolDiffusionScvf);
			register_import(m_imViscosityScvf);
			register_import(m_imDensityScv);
			register_import(m_imDensityScvf);
			register_import(m_imDarcyVelScvf);

		//	connect to own export
			m_imBrineScvf.set_data(m_exBrine);
			m_imBrineGradScvf.set_data(m_exBrineGrad);
			m_imPressureGradScvf.set_data(m_exPressureGrad);
			m_imDarcyVelScvf.set_data(m_exDarcyVel);
		}

	///	sets usage of consistent gravity
		void set_consistent_gravity(bool bUse)
		{
			m_bConsGravity = bUse;
			register_export_functions(Int2Type<dim>());
		}

	///	sets usage of boussinesq approximation for transport equation
		void set_boussinesq_transport(bool bUse) {m_BoussinesqTransport = bUse;}

	///	sets usage of boussinesq approcimation for flow equation
		void set_boussinesq_flow(bool bUse) {m_BoussinesqFlow = bUse;}

	///	sets the type of upwind
	/**
	 * This methods lets the user choose the type of upwind
	 * Currently, there are three implementations:
	 * <ul>
	 * <li> no
	 * <li> full
	 * <li> part
	 * </ul>
	 */
		bool set_upwind(std::string upwind)
		{
			if(upwind == "no")		   m_Upwind = NO_UPWIND;
			else if(upwind == "full")  m_Upwind = FULL_UPWIND;
			else if(upwind == "part")  m_Upwind = PART_UPWIND;
			else
			{
#ifndef FOR_VRL
				UG_LOG("Upwind type not found.\n");
				return false;
#endif
			}
			return true;
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
		template <typename TElem>
		inline bool prepare_element_loop();

		template <typename TElem>
		inline bool prepare_element(TElem* elem, const local_vector_type& u,
		                            const local_index_type& glob_ind);

		template <typename TElem>
		inline bool finish_element_loop();

		template <typename TElem>
		inline bool assemble_JA(local_matrix_type& J,
		                        const local_vector_type& u, number time=0.0);

		template <typename TElem>
		inline bool assemble_JM(local_matrix_type& J,
		                        const local_vector_type& u, number time=0.0);

		template <typename TElem>
		inline bool assemble_A(local_vector_type& d,
		                       const local_vector_type& u, number time=0.0);

		template <typename TElem>
		inline bool assemble_M(local_vector_type& d,
		                       const local_vector_type& u, number time=0.0);

		template <typename TElem>
		inline bool assemble_f(local_vector_type& d, number time=0.0);

	private:
	/// Upwind (1.0 == full upwind, 0.0 == no upwind)
		enum UpwindType{ NO_UPWIND = 0, FULL_UPWIND, PART_UPWIND};
		int m_Upwind;

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
		DataImport<number, dim, algebra_type> m_imBrineScvf;
		DataImport<MathVector<dim>, dim, algebra_type> m_imBrineGradScvf;

	///	Data import for pressure gradient at scvf ips
		DataImport<MathVector<dim>, dim, algebra_type> m_imPressureGradScvf;

	///	Data import for the reaction term
		DataImport<number, dim, algebra_type> m_imPorosityScv;
		DataImport<number, dim, algebra_type> m_imPorosityScvf;

	///	Data import for gravity (must be constant in current implementation)
		DataImport<MathVector<dim>, dim, algebra_type> m_imConstGravity;

	///	Data import for permeability tensor
		DataImport<MathMatrix<dim, dim>, dim, algebra_type> m_imPermeabilityScvf;

	///	Data import for molecular diffusion tensor
		DataImport<MathMatrix<dim, dim>, dim, algebra_type> m_imMolDiffusionScvf;

	///	Data import for viscosity
		DataImport<number, dim, algebra_type> m_imViscosityScvf;

	///	Data import for density
		DataImport<number, dim, algebra_type> m_imDensityScv;
		DataImport<number, dim, algebra_type> m_imDensityScvf;

	///	Data import for Darcy Velocity
		DataImport<MathVector<dim>, dim, algebra_type> m_imDarcyVelScvf;

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
		DataExport<MathVector<dim>, dim, algebra_type> m_exDarcyVel;

	///	Export for the brine mass fraction
		DataExport<number, dim, algebra_type> m_exBrine;

	///	Export for the gradient of brine mass fraction
		DataExport<MathVector<dim>, dim, algebra_type> m_exBrineGrad;

	///	Export for the gradient of brine mass fraction
		DataExport<MathVector<dim>, dim, algebra_type> m_exPressureGrad;

	private:
		// register for 1D
		void register_assemble_functions(Int2Type<1>)
		{
			register_all_assemble_functions<Edge>(ROID_EDGE);
		}

		// register for 2D
		void register_assemble_functions(Int2Type<2>)
		{
			//register_assemble_functions(Int2Type<1>());
			register_all_assemble_functions<Triangle>(ROID_TRIANGLE);
			register_all_assemble_functions<Quadrilateral>(ROID_QUADRILATERAL);
		}

		// register for 3D
		void register_assemble_functions(Int2Type<3>)
		{
			//register_assemble_functions(Int2Type<2>());
			register_all_assemble_functions<Tetrahedron>(ROID_TETRAHEDRON);
			register_all_assemble_functions<Pyramid>(ROID_PYRAMID);
			register_all_assemble_functions<Prism>(ROID_PRISM);
			register_all_assemble_functions<Hexahedron>(ROID_HEXAHEDRON);
		}

		// help function
		template <typename TElem>
		void register_all_assemble_functions(int id)
		{
			typedef DensityDrivenFlowElemDisc T;

			register_prepare_element_loop_function(	id, &T::template prepare_element_loop<TElem>);
			register_prepare_element_function(		id, &T::template prepare_element<TElem>);
			register_finish_element_loop_function(	id, &T::template finish_element_loop<TElem>);
			register_assemble_JA_function(			id, &T::template assemble_JA<TElem>);
			register_assemble_JM_function(			id, &T::template assemble_JM<TElem>);
			register_assemble_A_function(			id, &T::template assemble_A<TElem>);
			register_assemble_M_function(			id, &T::template assemble_M<TElem>);
			register_assemble_f_function(			id, &T::template assemble_f<TElem>);
		}

	// register for 1D
		void register_export_functions(Int2Type<1>)
		{
			register_all_export_functions<Edge>(ROID_EDGE);
		}

		// register for 2D
		void register_export_functions(Int2Type<2>)
		{
			//register_export_functions(Int2Type<1>());
			register_all_export_functions<Triangle>(ROID_TRIANGLE);
			register_all_export_functions<Quadrilateral>(ROID_QUADRILATERAL);
		}

		// register for 3D
		void register_export_functions(Int2Type<3>)
		{
			//register_export_functions(Int2Type<2>());
			register_all_export_functions<Tetrahedron>(ROID_TETRAHEDRON);
			register_all_export_functions<Pyramid>(ROID_PYRAMID);
			register_all_export_functions<Prism>(ROID_PRISM);
			register_all_export_functions<Hexahedron>(ROID_HEXAHEDRON);
		}

		// help function
		template <typename TElem>
		void register_all_export_functions(int id)
		{
			typedef DensityDrivenFlowElemDisc T;

			if(m_bConsGravity)
				m_exDarcyVel.register_export_func(id, this, &T::template compute_darcy_export_cons_grav<TElem>);
			else
				m_exDarcyVel.register_export_func(id, this, &T::template compute_darcy_export_std<TElem>);

			m_exBrine.register_export_func(id, this, &T::template compute_brine_export<TElem>);
			m_exBrineGrad.register_export_func(id, this, &T::template compute_brine_grad_export<TElem>);
			m_exPressureGrad.register_export_func(id, this, &T::template compute_pressure_grad_export<TElem>);
		}

};

///@}

} // end namespace ug

#include "density_driven_flow_impl.h"

#endif /*__H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__DENSITY_DRIVEN_FLOW__*/
