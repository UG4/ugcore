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
#include "lib_discretization/spacial_discretization/disc_helper/disc_helper.h"
#include "lib_discretization/spacial_discretization/disc_helper/fvgeom.h"
#include "lib_discretization/spacial_discretization/elem_disc/elem_disc_interface.h"
#include "lib_discretization/common/local_algebra.h"
#include "consistent_gravity.h"
#include "convection_shape.h"

namespace ug{

template <int dim>
class IDensityDrivenFlowUserFunction
{
	public:
	//	Function Types
		typedef void (*Porosity_fct)(number&);
		typedef void (*Viscosity_fct)(number&, number);
		typedef void (*Density_fct)(number&, number);
		typedef void (*D_Density_fct)(number&, number);
		typedef void (*Mol_Diff_Tensor_fct)(MathMatrix<dim,dim>&);
		typedef void (*Permeability_Tensor_fct)(MathMatrix<dim,dim>&);
		typedef void (*Gravity_fct)(MathVector<dim>&);

	public:
		virtual Porosity_fct get_porosity_function() const = 0;
		virtual Viscosity_fct get_viscosity_function() const = 0;
		virtual Density_fct get_density_function() const = 0;
		virtual D_Density_fct get_d_density_function() const = 0;
		virtual Mol_Diff_Tensor_fct get_mol_diff_tensor_function() const = 0;
		virtual Permeability_Tensor_fct get_perm_tensor_function() const = 0;
		virtual Gravity_fct get_gravity_function() const = 0;

		virtual ~IDensityDrivenFlowUserFunction(){}
};

template<template <class TElem, int TWorldDim> class TFVGeom, typename TDomain, typename TAlgebra>
class DensityDrivenFlowElemDisc  : public IElemDisc<TAlgebra> {
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
		typedef LocalMatrix<typename TAlgebra::matrix_type::value_type> local_matrix_type;

		// local vector type
		typedef LocalVector<typename TAlgebra::vector_type::value_type> local_vector_type;

		// local index type
		typedef LocalIndices local_index_type;

	protected:
		typedef void (*Porosity_fct)(number&);
		typedef void (*Viscosity_fct)(number&, number);
		typedef void (*Density_fct)(number&, number);
		typedef void (*D_Density_fct)(number&, number);
		typedef void (*Mol_Diff_Tensor_fct)(MathMatrix<dim,dim>&);
		typedef void (*Permeability_Tensor_fct)(MathMatrix<dim,dim>&);
		typedef void (*Gravity_fct)(MathVector<dim>&);

	public:
		DensityDrivenFlowElemDisc(	TDomain& domain, number upwind_amount,
									Porosity_fct Porosity, Viscosity_fct Viscosity, Density_fct Density, D_Density_fct D_Density,
									Mol_Diff_Tensor_fct Mol_Diff, Permeability_Tensor_fct Permeability_Tensor, Gravity_fct Gravity);

		DensityDrivenFlowElemDisc() :
			m_pDomain(NULL), m_Upwind(NO_UPWIND), m_UseConsistentGravity(true),
			m_BoussinesqTransport(true), m_BoussinesqFlow(true),
			m_PorosityFct(NULL), m_ViscosityFct(NULL), m_DensityFct(NULL), m_DDensityFct(NULL),
			m_MolDiffTensorFct(NULL), m_PermeabilityTensorFct(NULL), m_GravityFct(NULL)
		{
			register_assemble_functions(Int2Type<dim>());
		}

		void set_consistent_gravity(bool bUse) {m_UseConsistentGravity = bUse;}
		void set_boussinesq_transport(bool bUse) {m_BoussinesqTransport = bUse;}
		void set_boussinesq_flow(bool bUse) {m_BoussinesqFlow = bUse;}
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
		void set_approximation_space(IApproximationSpace<domain_type>& approxSpace)
		{
			m_pDomain = &approxSpace.get_domain();
			set_pattern(approxSpace.get_function_pattern());
		}
		void set_domain(domain_type& domain) {m_pDomain = &domain;}
		void set_user_functions(IDensityDrivenFlowUserFunction<dim>& user) {
			m_PorosityFct = user.get_porosity_function();
			m_ViscosityFct = user.get_viscosity_function();
			m_DensityFct = user.get_density_function();
			m_DDensityFct = user.get_d_density_function();
			m_MolDiffTensorFct = user.get_mol_diff_tensor_function();
			m_PermeabilityTensorFct = user.get_perm_tensor_function();
			m_GravityFct = user.get_gravity_function();
		}

		virtual size_t num_fct() {return 2;}

		virtual LocalShapeFunctionSetID local_shape_function_set_id(size_t loc_fct)
		{
			return LocalShapeFunctionSetID(LocalShapeFunctionSetID::LAGRANGE, 1);
		}

	private:
		template <typename TElem>
		inline bool prepare_element_loop();

		template <typename TElem>
		inline bool prepare_element(TElem* elem, const local_vector_type& u, const local_index_type& glob_ind);

		template <typename TElem>
		inline bool finish_element_loop();

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
		TDomain* m_pDomain;

		// position access
		typename TDomain::position_accessor_type m_aaPos;

		// amount of upwind (1.0 == full upwind, 0.0 == no upwind)
		enum UpwindType{ NO_UPWIND = 0, FULL_UPWIND, PART_UPWIND};
		int m_Upwind;
		bool m_UseConsistentGravity;
		bool m_BoussinesqTransport;
		bool m_BoussinesqFlow;

		// max num of element Corners (this is to much in 3d)
		static const size_t m_MaxNumCorners = dim*dim;

		// positions
		position_type m_vCornerCoords[m_MaxNumCorners];

		// constant values, read in once
		number m_Porosity;
		MathMatrix<dim, dim> m_Diffusion;
		MathVector<dim> m_Gravity;
		MathMatrix<dim, dim> m_Permeability;

		// density at corners
		number m_vDensity[m_MaxNumCorners];
		number m_vDDensity[m_MaxNumCorners];

		// consistent gravity
		MathVector<dim> m_vConsGravity[m_MaxNumCorners]; // Consistent Gravity at corners
		MathVector<dim> m_vvDConsGravity[m_MaxNumCorners][m_MaxNumCorners]; //Derivative at corners

		// local constants for readability
		// _C_ == local function 0 (concentration)
		// _P_ == local function 1 (pressure)
		static const size_t _C_ = 0;
		static const size_t _P_ = 1;

	private:
		// User functions
		Porosity_fct m_PorosityFct;
		Viscosity_fct m_ViscosityFct;
		Density_fct m_DensityFct;
		D_Density_fct m_DDensityFct;
		Mol_Diff_Tensor_fct m_MolDiffTensorFct;
		Permeability_Tensor_fct m_PermeabilityTensorFct;
		Gravity_fct m_GravityFct;

	protected:
		bool compute_ip_Darcy_velocity(
				MathVector<dim>& DarcyVel,
				number c,
				const MathVector<dim>& grad_p,
				const MathMatrix<dim, dim>& JTInv,
				const std::vector<MathVector<dim> >& vLocalGrad);
		bool compute_D_ip_Darcy_velocity(
				MathVector<dim>& DarcyVel,
				MathVector<dim> D_DarcyVel_c[],
				MathVector<dim> D_DarcyVel_p[],
				number c,
				const MathVector<dim>& grad_p,
				const MathMatrix<dim, dim>& JTInv,
				const std::vector<number>& vShape,
				const std::vector<MathVector<dim> >& vLocalGrad,
				const std::vector<MathVector<dim> >& vLGlobalGrad);

	private:
		///////////////////////////////////////
		// registering for reference elements
		///////////////////////////////////////
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
			register_prepare_element_loop_function(	id, &DensityDrivenFlowElemDisc::template prepare_element_loop<TElem>);
			register_prepare_element_function(		id, &DensityDrivenFlowElemDisc::template prepare_element<TElem>);
			register_finish_element_loop_function(	id, &DensityDrivenFlowElemDisc::template finish_element_loop<TElem>);
			register_assemble_JA_function(			id, &DensityDrivenFlowElemDisc::template assemble_JA<TElem>);
			register_assemble_JM_function(			id, &DensityDrivenFlowElemDisc::template assemble_JM<TElem>);
			register_assemble_A_function(			id, &DensityDrivenFlowElemDisc::template assemble_A<TElem>);
			register_assemble_M_function(			id, &DensityDrivenFlowElemDisc::template assemble_M<TElem>);
			register_assemble_f_function(			id, &DensityDrivenFlowElemDisc::template assemble_f<TElem>);
		}
};


} // end namespace ug

#include "density_driven_flow_impl.h"

#endif /*__H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__DENSITY_DRIVEN_FLOW__FV1__DENSITY_DRIVEN_FLOW__*/
