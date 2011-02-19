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
#include "lib_algebra/lib_algebra.h"

// library intern headers
#include "lib_discretization/spatial_discretization/disc_helper/disc_helper.h"
#include "lib_discretization/spatial_discretization/disc_helper/fvgeom.h"
#include "lib_discretization/spatial_discretization/elem_disc/elem_disc_interface.h"
#include "lib_discretization/common/local_algebra.h"
#include "lib_discretization/spatial_discretization/ip_data/user_data.h"
#include "upwind.h"
#include "stabilization.h"

namespace ug{

template<template <class TElem, int TWorldDim> class TFVGeom, typename TDomain, typename TAlgebra>
class FVNavierStokesElemDisc : public IElemDisc<TAlgebra>
{
	public:
		// domain type
		typedef TDomain domain_type;

		// world dimension
		static const size_t dim = TDomain::dim;

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

	public:
	//	Constructor (setting default values)
		FVNavierStokesElemDisc()
		 : m_pDomain(NULL),m_Viscosity(1.0)
		   {
                set_upwind("FULL_UPWIND");
                set_stabilization("FIELDS");
                set_diffusionlength("RAW");
                set_physicalAdvectionCorrection("NOPAC");
                set_pecletBlend("NOPEBLEND");
				register_assemble_functions(Int2Type<dim>());
			}


	public:
	//	Setup
		void set_domain(domain_type& domain) {m_pDomain = &domain;}
		void set_kinematicViscosity(number nu) {m_Viscosity = nu;}

	private:
		int m_UpwindMethod;
        int m_StabOptions[STAB_OPTIONS];

	public:
		bool set_upwind(const std::string& upwind)
		{
			if      (upwind == "NoUpwind")  m_UpwindMethod = NO_UPWIND;
            else if (upwind == "Full")      m_UpwindMethod = FULL_UPWIND;
            else if (upwind == "LPS")       m_UpwindMethod = LPS_UPWIND;
            else if (upwind == "NUM")       m_UpwindMethod = NUM_UPWIND;
            else {UG_LOG("Upwind Type not recognized.\n"); return false;}
			return true;
		}

		bool set_stabilization(const std::string& stabilization)
		{
			if      (stabilization == "FIELDS")  m_StabOptions[STAB_METHOD] = FIELDS;
            else if (stabilization == "FLOW")    m_StabOptions[STAB_METHOD] = FLOW;
            else {UG_LOG("Stabilization Type not recognized.\n"); return false;}
			return true;
		}

        bool set_diffusionlength(const std::string& diffusionlength)
		{
			if      (diffusionlength == "RAW")        m_StabOptions[DIFF_LENGTH_METHOD] = RAW;
            else if (diffusionlength == "FIVEPOINT")  m_StabOptions[DIFF_LENGTH_METHOD] = FIVEPOINT;
            else if (diffusionlength == "COR")        m_StabOptions[DIFF_LENGTH_METHOD] = COR;
            else {UG_LOG("Diffusion Length calculation method not recognized.\n"); return false;}
			return true;
		}

        bool set_physicalAdvectionCorrection(const std::string& physicalAdvectionCorrection)
		{
			if      (physicalAdvectionCorrection == "PAC")    m_StabOptions[PAC_SWITCH] = PAC;
            else if (physicalAdvectionCorrection == "NOPAC")  m_StabOptions[PAC_SWITCH] = NOPAC;
            else {UG_LOG("PAC term usage definition not recognized.\n"); return false;}
			return true;
		}

        bool set_pecletBlend(const std::string& pecletBlend)
		{
			if      (pecletBlend == "PEBLEND")    m_StabOptions[PECLET_BLENDING_SWITCH] = PEBLEND;
            else if (pecletBlend == "NOPEBLEND")  m_StabOptions[PECLET_BLENDING_SWITCH] = NOPEBLEND;
            else {UG_LOG("Peclet Blending setup not recognized.\n"); return false;}
			return true;
		}

	public:
	//	Implement Interface
		virtual size_t num_fct(){return dim+1;}

		virtual LocalShapeFunctionSetID local_shape_function_set_id(size_t loc_fct)
		{
			return LocalShapeFunctionSetID(LocalShapeFunctionSetID::LAGRANGE, 1);
		}

		virtual bool use_hanging() const {return TFVGeom<Edge, dim>::usesHangingNodes;}

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
		std::vector<position_type> m_vCornerCoords;
		typename TDomain::position_accessor_type m_aaPos;

		// abbreviation for pressure
        static const size_t _U_ = 0;
        static const size_t _V_ = 1;
        static const size_t _W_ = 2;
		static const size_t _P_ = dim;

		// User functions
		number m_Viscosity;

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
			register_assemble_functions(Int2Type<1>());
			register_all_assemble_functions<Triangle>(ROID_TRIANGLE);
			register_all_assemble_functions<Quadrilateral>(ROID_QUADRILATERAL);
		}

		// register for 3D
		void register_assemble_functions(Int2Type<3>)
		{
			register_assemble_functions(Int2Type<2>());
			register_all_assemble_functions<Tetrahedron>(ROID_TETRAHEDRON);
			register_all_assemble_functions<Pyramid>(ROID_PYRAMID);
			register_all_assemble_functions<Prism>(ROID_PRISM);
			register_all_assemble_functions<Hexahedron>(ROID_HEXAHEDRON);
		}

		// help function
		template <typename TElem>
		void register_all_assemble_functions(int id)
		{
			register_prepare_element_loop_function(	id, &FVNavierStokesElemDisc::template prepare_element_loop<TElem>);
			register_prepare_element_function(		id, &FVNavierStokesElemDisc::template prepare_element<TElem>);
			register_finish_element_loop_function(	id, &FVNavierStokesElemDisc::template finish_element_loop<TElem>);
			register_assemble_JA_function(			id, &FVNavierStokesElemDisc::template assemble_JA<TElem>);
			register_assemble_JM_function(			id, &FVNavierStokesElemDisc::template assemble_JM<TElem>);
			register_assemble_A_function(			id, &FVNavierStokesElemDisc::template assemble_A<TElem>);
			register_assemble_M_function(			id, &FVNavierStokesElemDisc::template assemble_M<TElem>);
			register_assemble_f_function(			id, &FVNavierStokesElemDisc::template assemble_f<TElem>);
		}

};

}

#include "navier_stokes_impl.h"

#endif /*__H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__NAVIER_STOKES__FV__NAVIER_STOKES__*/
