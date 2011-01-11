/*
 * density_driven_flow.h
 *
 *  Created on: 26.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__DENSITY_DRIVEN_FLOW__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__DENSITY_DRIVEN_FLOW__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"
#include "lib_algebra/lib_algebra.h"

// library intern headers
#include "lib_discretization/spatial_discretization/disc_helper/disc_helper.h"
#include "lib_discretization/spatial_discretization/disc_helper/fvgeom.h"
#include "lib_discretization/spatial_discretization/elem_disc/elem_disc_interface.h"
#include "lib_discretization/common/local_algebra.h"
#include "consistent_gravity.h"
#include "convection_shape.h"

namespace ug{

template <int dim>
class IDensityDrivenFlowUserFunction
{
	public:
	//	Function Types
		typedef void (*Viscosity_fct)(number&, number);
		typedef void (*Density_fct)(number&, number);
		typedef void (*D_Density_fct)(number&, number);

	public:
		virtual Viscosity_fct get_viscosity_function() const = 0;
		virtual Density_fct get_density_function() const = 0;
		virtual D_Density_fct get_d_density_function() const = 0;

		virtual ~IDensityDrivenFlowUserFunction(){}
};

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
class DensityDrivenFlowElemDisc  : public IElemDisc<TAlgebra> {
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

	protected:
		typedef void (*Viscosity_fct)(number&, number);
		typedef void (*Density_fct)(number&, number);
		typedef void (*D_Density_fct)(number&, number);

	public:
	///	Constructor
		DensityDrivenFlowElemDisc() :
			m_pDomain(NULL), m_Upwind(NO_UPWIND), m_UseConsistentGravity(true),
			m_BoussinesqTransport(true), m_BoussinesqFlow(true),
			m_ViscosityFct(NULL),
			m_DensityFct(NULL), m_DDensityFct(NULL)
		{
		//	register assemble functions
			register_assemble_functions(Int2Type<dim>());

		//	register export
			register_export(m_DarcyVelExport);

		//	register import
			register_import(m_scvfPorosity);
			register_import(m_scvPorosity);
			register_import(m_scvfPermeability);
			register_import(m_scvfMolDiffusion);
		}

	///	sets usage of consistent gravity
		void set_consistent_gravity(bool bUse) {m_UseConsistentGravity = bUse;}

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

	///	sets the approximation space
		void set_approximation_space(IApproximationSpace<domain_type>& approxSpace)
		{
			m_pDomain = &approxSpace.get_domain();
			set_pattern(approxSpace.get_function_pattern());
		}

	///	sets the domain
		void set_domain(domain_type& domain) {m_pDomain = &domain;}

	///	sets the user functions
		void set_user_functions(IDensityDrivenFlowUserFunction<dim>& user) {
			m_ViscosityFct = user.get_viscosity_function();
			m_DensityFct = user.get_density_function();
			m_DDensityFct = user.get_d_density_function();
		}

	///	sets the porosity
	/**
	 * This method sets the Porosity. (Dimensionless)
	 */
		void set_porosity(IPData<number, dim>& user)
		{
			m_scvPorosity.set_data(user);
			m_scvfPorosity.set_data(user);
		}

	///	sets the gravity vector
	/**
	 * This method sets the Gravity. (Unit is \f$ \frac{m}{s^2} \f$)
	 */
		void set_gravity(IPData<MathVector<dim>, dim>& user)
		{
			m_constGravity.set_data(user);
		}

	///	sets the molecular diffusion tensor
	/**
	 * This method sets the molecular Diffusion tensor.
	 */
		void set_molecular_diffusion(IPData<MathMatrix<dim, dim>, dim>& user)
		{
			m_scvfMolDiffusion.set_data(user);
		}

	///	sets the permeability tensor
	/**
	 * This method sets the Permeability tensor.
	 */
		void set_permeability(IPData<MathMatrix<dim, dim>, dim>& user)
		{
			m_scvfPermeability.set_data(user);
		}

	public:
	///	number of functions used
		virtual size_t num_fct() {return 2;}

	///	trial space for functions
		virtual LocalShapeFunctionSetID local_shape_function_set_id(size_t loc_fct)
		{
			return LocalShapeFunctionSetID(LocalShapeFunctionSetID::LAGRANGE, 1);
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
	///	Domain
		TDomain* m_pDomain;

	///	Position access
		typename TDomain::position_accessor_type m_aaPos;

	/// Upwind (1.0 == full upwind, 0.0 == no upwind)
		enum UpwindType{ NO_UPWIND = 0, FULL_UPWIND, PART_UPWIND};
		int m_Upwind;

	///	flag if using Consistent Gravity
		bool m_UseConsistentGravity;

	///	flag if using boussinesq transport
		bool m_BoussinesqTransport;

	///	flag if using boussinesq flow
		bool m_BoussinesqFlow;

	///	max num of element Corners (this is to much in 3d)
		static const size_t m_MaxNumCorners = dim*dim;

	///	Corner Coordinates
		std::vector<position_type> m_vCornerCoords;

	/// constant Gravity, read in once
		MathVector<dim> m_Gravity;

	///	help Matrix for Diffusion
		MathMatrix<dim, dim> m_Diffusion;

		// density at corners
		number m_vDensity[m_MaxNumCorners];
		number m_vDDensity[m_MaxNumCorners];

		// consistent gravity
		MathVector<dim> m_vConsGravity[m_MaxNumCorners]; // Consistent Gravity at corners
		MathVector<dim> m_vvDConsGravity[m_MaxNumCorners][m_MaxNumCorners]; //Derivative at corners

	///	abbreviation for local function: brine mass fraction
		static const size_t _C_ = 0;

	///	abbreviation for local function: pressure
		static const size_t _P_ = 1;

	private:
	///	Data import for the reaction term
		DataImport<number, dim, algebra_type> m_scvPorosity;
		DataImport<number, dim, algebra_type> m_scvfPorosity;

	///	Data import for gravity (must be constant in current implementation)
		DataImport<MathVector<dim>, dim, algebra_type> m_constGravity;

	///	Data import for permeability tensor
		DataImport<MathMatrix<dim, dim>, dim, algebra_type> m_scvfPermeability;

	///	Data import for molecular diffusion tensor
		DataImport<MathMatrix<dim, dim>, dim, algebra_type> m_scvfMolDiffusion;

	private:
		// User functions
		Viscosity_fct m_ViscosityFct;
		Density_fct m_DensityFct;
		D_Density_fct m_DDensityFct;

	protected:
	///	computes the darcy velocity
		bool compute_ip_Darcy_velocity(
				MathVector<dim>& DarcyVel,
				number c,
				const MathVector<dim>& grad_p,
				const MathMatrix<dim, dim>& JTInv,
				const std::vector<MathVector<dim> >& vLocalGrad,
				size_t scvfIP);

	///	computes the darcy velocity and its derivative w.r.t to c,p
		bool compute_D_ip_Darcy_velocity(
				MathVector<dim>& DarcyVel,
				MathVector<dim> D_DarcyVel_c[],
				MathVector<dim> D_DarcyVel_p[],
				number c,
				const MathVector<dim>& grad_p,
				const MathMatrix<dim, dim>& JTInv,
				const std::vector<number>& vShape,
				const std::vector<MathVector<dim> >& vLocalGrad,
				const std::vector<MathVector<dim> >& vLGlobalGrad,
				size_t scvfIP);

	protected:
	///	computes the export
		template <typename TElem>
		bool compute_darcy_export(const local_vector_type& u, bool compDeriv)
		{
		// 	Get finite volume geometry
			static const TFVGeom<TElem, dim>& geo =
						FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

		//	Number of Corners
			typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
			static const size_t num_co = ref_elem_type::num_corners;

		//	Prepare Consistent Gravity if needed
			if(m_UseConsistentGravity)
			{
			// 	Prepare Density in Corners
				if(!PrepareConsistentGravity<dim>(	&m_vConsGravity[0],
													num_co,
													&m_vCornerCoords[0],
													&m_vDensity[0],
													m_Gravity))
				{
					UG_LOG("ERROR in assemble_JA: Cannot "
							"prepare Consistent Gravity.\n");
					return false;
				}

			// 	Prepare DensityDerivative in Corners
				number DCoVal[num_co];
				memset(DCoVal, 0, sizeof(number)*num_co);

				for(size_t sh = 0; sh < num_co; sh++)
				{
					DCoVal[sh] = m_vDDensity[sh];
					if(!PrepareConsistentGravity<dim>(	&m_vvDConsGravity[sh][0],
														num_co,
														&m_vCornerCoords[0],
														&DCoVal[0],
														m_Gravity))
					{
						UG_LOG("ERROR in assemble_JA: Cannot "
								"prepare Consistent Gravity.\n");
						return false;
					}
					DCoVal[sh] = 0.0;
				}
			}

			static const int refDim =
					reference_element_traits<TElem>::reference_element_type::dim;

			number c_ip;
			MathVector<dim> grad_p_ip, grad_c_ip;	// Gradients

			for(size_t s = 0; s < m_DarcyVelExport.num_series(); ++s)
			{
			//	currently only for FV1 elem dofs implemented
				if(m_DarcyVelExport.template local_ips<refDim>(s) != geo.scvf_local_ips())
				{
					UG_LOG("ERROR in 'compute_darcy_export': Currently export"
							" of Darcy Velocity only implemented for FV1 and"
							"SCVF integration points.\n");
					return false;
				}

		//	Loop Sub Control Volume Faces (SCVF)
			size_t ip = 0;
			for(size_t i = 0; i < geo.num_scvf(); ++i)
			{
			// 	Get current SCVF
				const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(i);

			// 	Loop integration point of SCVF
				for(size_t j = 0; j < scvf.num_ip(); ++j, ++ip)
				{
					///////////////////////////////////////////
					// IP Values
					///////////////////////////////////////////

				//	Compute Gradients and concentration at ip
					VecSet(grad_p_ip, 0.0); VecSet(grad_c_ip, 0.0); c_ip = 0.0;
					for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					{
						VecScaleAppend(grad_p_ip, u(_P_,sh), scvf.global_grad(sh, j));
						VecScaleAppend(grad_c_ip, u(_C_,sh), scvf.global_grad(sh, j));
						c_ip += u(_C_, sh) * scvf.shape(sh, j);
					}

				//	Compute density and its derivative at ip
					number rho_ip, DRho_ip;
					if(!m_BoussinesqTransport || !m_BoussinesqFlow)
					{
						m_DensityFct(rho_ip, c_ip);
						m_DDensityFct(DRho_ip, c_ip);
					}

					///////////////////////////////////////////
					// Darcy Velocity at ip
					///////////////////////////////////////////

				//	Compute Darcy Velocity
					compute_D_ip_Darcy_velocity(	m_DarcyVelExport.value(s, ip),
					                            	m_DarcyVelExport.deriv(s, ip, _C_),
					                            	m_DarcyVelExport.deriv(s, ip, _P_),
													c_ip, grad_p_ip,
													scvf.JTInv(j),
													scvf.shape_vector(j),
													scvf.local_grad_vector(j),
													scvf.global_grad_vector(j),
													ip);
				}
			}

			}
		//	we're done
			return true;
		}

	///	computes the export
		template <typename TElem>
		bool compute_brine_export(const local_vector_type& u, bool compDeriv)
		{
		// 	Get finite volume geometry
			static const TFVGeom<TElem, dim>& geo =
						FVGeometryProvider::get_geom<TFVGeom, TElem,dim>();

			typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
			static const size_t refDim = ref_elem_type::dim;
			static const size_t num_co = ref_elem_type::num_corners;

			for(size_t s = 0; s < m_BrineExport.num_series(); ++s)
			{
			//	FV1 SCVF ip
				if(m_BrineExport.template local_ips<refDim>(s)
						== geo.scvf_local_ips())
				{
				//	Loop Sub Control Volume Faces (SCVF)
					size_t ip = 0;
					for(size_t i = 0; i < geo.num_scvf(); ++i)
					{
					// 	Get current SCVF
						const typename TFVGeom<TElem, dim>::SCVF& scvf = geo.scvf(i);

					// 	Loop integration point of SCVF
						for(size_t j = 0; j < scvf.num_ip(); ++j, ++ip)
						{
						//	Compute Gradients and concentration at ip
							number& cIP = m_BrineExport.value(s, ip);
							for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
								cIP += u(_C_, sh) * scvf.shape(sh, j);

							if(compDeriv)
							{
								number* cIP_c = m_BrineExport.deriv(s, ip, _C_);
								number* cIP_p = m_BrineExport.deriv(s, ip, _P_);

								for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
								{
									cIP_c[sh] = scvf.shape(sh, j);
									cIP_p[sh] = 0.0;
								}
							}
						}
					}
				}

			//	FV1 SCV ip
				if(m_BrineExport.template local_ips<refDim>(s)
						== geo.scv_local_ips())
				{
				//	Loop Corners
					for(size_t ip = 0; ip < num_co; ip++)
					{
					//	Compute Gradients and concentration at ip
						m_BrineExport.value(s, ip) = u(_C_, ip);

						if(compDeriv)
						{
							number* cIP_c = m_BrineExport.deriv(s, ip, _C_);
							number* cIP_p = m_BrineExport.deriv(s, ip, _P_);

							for(size_t sh = 0; sh < num_co; ++sh)
							{
								cIP_c[sh] = (ip==sh) ? 1.0 : 0.0;
								cIP_p[sh] = 0.0;
							}
						}
					}
				}

				// others not implemented
				UG_LOG("Evaluation not implemented.");
				return false;
			}

		//	we're done
			return true;
		}

	///	Export for the Darcy velocity
		DataExport<MathVector<dim>, dim, algebra_type> m_DarcyVelExport;

	///	Export for the brine mass fraction
		DataExport<number, dim, algebra_type> m_BrineExport;

	public:
	///	returns the export of the darcy velocity
		IPData<MathVector<dim>, dim>& get_darcy_velocity() {return m_DarcyVelExport;}

	///	returns the export of brine mass fracture
		IPData<number, dim>& get_brine() {return m_BrineExport;}

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

			m_DarcyVelExport.register_export_func(id, this, &T::template compute_darcy_export<TElem>);
			m_BrineExport.register_export_func(id, this, &T::template compute_brine_export<TElem>);

		}
};

///@}

} // end namespace ug

#include "density_driven_flow_impl.h"

#endif /*__H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__DENSITY_DRIVEN_FLOW__*/
