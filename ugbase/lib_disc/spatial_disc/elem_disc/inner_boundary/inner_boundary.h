/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Markus Breit
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

/*
 * Finite Volume Element Discretization for an inner BndCond that depends on the unknowns (on the bnd)
 *
 * This class implements the IElemDisc interface to provide element local
 * assemblings for the unknown-dependent Neumann-flux over an inner boundary.
 * The equation of this flux and its derivative should be given
 * in a concretization of this class.
 */

#ifndef __H__UG__LIB_DISC__SPACIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__INNER_BOUNDARY__
#define __H__UG__LIB_DISC__SPACIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__INNER_BOUNDARY__

#include <boost/function.hpp>
#include <vector>
#include <string>
#include <utility>      // std::pair

// other ug4 modules
#include "common/common.h"

// library intern headers
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/hfv1_geom.h"



namespace ug
{

struct InnerBoundaryConstants
{
	public:
		/// index value for which a flux is ignored
		static const size_t _IGNORE_ = -1;
};


/// struct that holds information about the flux densities and from where to where the flux occurs
struct InnerBoundaryFluxCond
{
	// vector of fluxFctValues
	std::vector<number> flux;
	std::vector<size_t> from;
	std::vector<size_t> to;
};

/// struct that holds information about the derivatives of the flux densities
/// and from where to where the flux occurs
struct InnerBoundaryFluxDerivCond
{
	// vector of fluxFctDerivValues
	// rows: fluxIndex, cols: dependency;
	// in the pair: first: with respect to which local function index, second: value
	std::vector<std::vector<std::pair<size_t,number> > > fluxDeriv;
	std::vector<size_t> from;
	std::vector<size_t> to;
};


template <typename TImpl, typename TDomain>
class FV1InnerBoundaryElemDisc
: public IElemDisc<TDomain>
{
	public:
		typedef InnerBoundaryFluxCond FluxCond;
		typedef InnerBoundaryFluxDerivCond FluxDerivCond;

	public:
	///	Domain type
		typedef TDomain domain_type;

	private:
	///	Base class type
		typedef IElemDisc<domain_type> base_type;

	///	own type
		typedef FV1InnerBoundaryElemDisc<TImpl, TDomain> this_type;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	///	Position type
		typedef typename base_type::position_type position_type;

	/// error estimator type
		typedef MultipleSideAndElemErrEstData<domain_type> err_est_type;


	public:

	/// Constructor with c-strings
		FV1InnerBoundaryElemDisc(const char* functions = "", const char* subsets = "")
			: IElemDisc<domain_type>(functions, subsets), m_bNonRegularGrid(false), m_bCurrElemIsHSlave(false), m_si(0)
		{
			register_all_fv1_funcs();
		}

	/// Constructor with functions
		FV1InnerBoundaryElemDisc(const std::vector<std::string>& functions, const std::vector<std::string>& subsets)
			: IElemDisc<domain_type>(functions, subsets), m_bNonRegularGrid(false), m_bCurrElemIsHSlave(false), m_si(0)
		{
			register_all_fv1_funcs();
		}

	/// Destructor
		virtual ~FV1InnerBoundaryElemDisc() {};

    /// Setting the flux function
        //void set_fluxFunction(UserData<number, dim>& fluxFct) {m_fluxFct.set_data(fluxFct);}
	
	/// Setting a scaling factor for the flux
		void set_flux_scale(number scale);
		void set_flux_scale(SmartPtr<CplUserData<number, dim> > scaleFct);
#ifdef UG_FOR_LUA
		void set_flux_scale(const char* luaScaleFctName);
#endif

	public:	// inherited from IElemDisc
	///	type of trial space for each function used
		virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid) override;

	///	returns if hanging nodes are used
		virtual bool use_hanging() const override;

	private:
	///	prepares the loop over all elements
	/**
	 * This method prepares the loop over all elements. It resizes the Position
	 * array for the corner coordinates and schedules the local ip positions
	 * at the data imports.
	 */
		template<typename TElem, typename TFVGeom>
		void prep_elem_loop(const ReferenceObjectID roid, const int si);

	///	prepares the element for assembling
	/**
	 * This methods prepares an element for the assembling. The Positions of
	 * the Element Corners are read and the Finite Volume Geometry is updated.
	 * The global ip positions are scheduled at the data imports.
	 */
		template<typename TElem, typename TFVGeom>
		void prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[]);

	///	finishes the loop over all elements
		template<typename TElem, typename TFVGeom>
		void fsh_elem_loop();

	///	assembles the local stiffness matrix using a finite volume scheme
		template<typename TElem, typename TFVGeom>
		void add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	assembles the local mass matrix using a finite volume scheme
		template<typename TElem, typename TFVGeom>
		void add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	assembles the stiffness part of the local defect
		template<typename TElem, typename TFVGeom>
		void add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	assembles the mass part of the local defect
		template<typename TElem, typename TFVGeom>
		void add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	assembles the local right hand side
		template<typename TElem, typename TFVGeom>
		void add_rhs_elem(LocalVector& rhs, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	prepares the loop over all elements of one type for the computation of the error estimator
		template <typename TElem, typename TFVGeom>
		void prep_err_est_elem_loop(const ReferenceObjectID roid, const int si);

	///	prepares the element for assembling the error estimator
		template <typename TElem, typename TFVGeom>
		void prep_err_est_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	computes the error estimator contribution for one element
		template <typename TElem, typename TFVGeom>
		void compute_err_est_A_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale);

	///	postprocesses the loop over all elements of one type in the computation of the error estimator
		template <typename TElem, typename TFVGeom>
		void fsh_err_est_elem_loop();

	private:
		DataImport<number, dim> m_fluxScale;  // data import for scaling of fluxes

	private:
		void register_all_fv1_funcs();

		struct RegisterFV1
		{
				// friend class this_type;
				RegisterFV1(this_type* pThis) : m_pThis(pThis){}
				this_type* m_pThis;
				template< typename TElem > void operator()(TElem&)
				{
					if (m_pThis->m_bNonRegularGrid)
						m_pThis->register_fv1_func<TElem, HFV1ManifoldGeometry<TElem, dim> >();
					else
						m_pThis->register_fv1_func<TElem, FV1ManifoldGeometry<TElem, dim> >();
				}
		};

		template <typename TElem, typename TFVGeom>
		void register_fv1_func();


		/// struct holding values of shape functions in IPs
		struct ShapeValues
		{
			public:
				void resize(size_t nSip, size_t _nSh)
				{
					nSh = _nSh;
					sideVals.resize(nSip);
					for (size_t i = 0; i < nSip; i++) sideVals[i].resize(nSh);
				}
				number& shapeAtSideIP(size_t sh, size_t ip)
				{
					UG_ASSERT(ip < sideVals.size(), "Requested data for IP " << ip << ", but only " << sideVals.size() << " IPs present.");
					UG_ASSERT(sh < sideVals[ip].size(), "Requested data for shape fct " << sh << ", but only " << sideVals[ip].size() << " shape fcts present.");
					return sideVals.at(ip).at(sh);
				}
				number* shapesAtSideIP(size_t ip) {return &sideVals[ip][0];}
				size_t num_sh() {return nSh;}
			private:
				size_t nSh;
				std::vector<std::vector<number> > sideVals;
		} m_shapeValues;

	private:
		bool m_bNonRegularGrid;
		bool m_bCurrElemIsHSlave;

		int m_si;
};

}

#endif /*__H__UG__LIB_DISC__SPACIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__INNER_BOUNDARY__*/
