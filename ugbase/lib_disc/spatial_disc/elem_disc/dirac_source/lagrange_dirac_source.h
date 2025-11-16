/*
 * Copyright (c) 2019-2020:  G-CSC, Goethe University Frankfurt
 * Author: Arne Naegel
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
 * 
 * \tparam	TDomain		Domain
 * \tparam	TAlgebra	Algebra
 */

#ifndef __H__UG__LIB_DISC__SPACIAL_DISCRETIZATION__ELEM_DISC__DIRAC_SOURCE__LAGRANGE_DIRAC_SOURCE_H__
#define __H__UG__LIB_DISC__SPACIAL_DISCRETIZATION__ELEM_DISC__DIRAC_SOURCE__LAGRANGE_DIRAC_SOURCE_H__

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




namespace ug
{



template<typename TDomain>
class DiracSourceDisc
: public IElemDisc<TDomain>
{

	public:
	///	Base class type
		using base_type = IElemDisc<TDomain>;

	///	Domain type
		using domain_type = typename base_type::domain_type;

	///	World dimension
		static constexpr int dim = base_type::dim;

		static constexpr int _C_ = 0;

	///	Position type
		using position_type = typename base_type::position_type;

	private:
	///	own type
		using this_type = DiracSourceDisc<TDomain>;


	public:

	/// Constructor with c-strings
		DiracSourceDisc(const char* functions = "", const char* subsets = "")
			: IElemDisc<TDomain>(functions, subsets), m_bNonRegularGrid(false) /*, m_bCurrElemIsHSlave(false)*/
		{
			register_all_funcs();
		}

	/// Constructor with functions
		DiracSourceDisc(const std::vector<std::string>& functions, const std::vector<std::string>& subsets)
			: IElemDisc<TDomain>(functions, subsets), m_bNonRegularGrid(false) /*, m_bCurrElemIsHSlave(false)*/
		{
			register_all_funcs();
		}

	/// Destructor
		virtual ~DiracSourceDisc() {};


	
	/// Setting a scaling factor for the flux
		void add_source(number scale, MathVector<dim> &srcCoord);
		void add_source(SmartPtr<UserData<number, dim> > srcData, MathVector<dim> &srcCoord);
#ifdef UG_FOR_LUA
		void add_source(const char* luaScaleFctName, MathVector<dim> &srcCoord);
#endif

		/// Setting a scaling factor for the flux
		void add_transport_sink(SmartPtr<UserData<number, dim> > snkData);
		void add_transport_sink(number snk);

#ifdef UG_FOR_LUA
		void add_transport_sink(const char* luaScaleFctName);
#endif

	public:	// inherited from IElemDisc
	///	type of trial space for each function used
		virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);

	///	returns if hanging nodes are used
		virtual bool use_hanging() const;

	protected:
	
	
	///	prepares the loop over all elements
	/**
	 * This method prepares the loop over all elements. It resizes the Position
	 * array for the corner coordinates and schedules the local ip positions
	 * at the data imports.
	 */
		template<typename TElem, typename TFVGeom>
		void prep_elem_loop(const ReferenceObjectID roid, const int si)
		{}

	///	prepares the element for assembling
	/**
	 * This methods prepares an element for the assembling. The Positions of
	 * the Element Corners are read and the Finite Volume Geometry is updated.
	 * The global ip positions are scheduled at the data imports.
	 */
		template<typename TElem, typename TFVGeom>
		void prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
		{}

	///	finishes the loop over all elements
		template<typename TElem, typename TFVGeom>
		void fsh_elem_loop()
		{}

	///	assembles the local stiffness matrix using a finite volume scheme
		template<typename TElem, typename TFVGeom>
		void add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	assembles the local mass matrix using a finite volume scheme
		template<typename TElem, typename TFVGeom>
		void add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{}

	///	assembles the stiffness part of the local defect
		template<typename TElem, typename TFVGeom>
		void add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	assembles the mass part of the local defect
		template<typename TElem, typename TFVGeom>
		void add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{}

	///	assembles the local right hand side
		template<typename TElem, typename TFVGeom>
		void add_rhs_elem(LocalVector& rhs, GridObject* elem, const MathVector<dim> vCornerCoords[]);



	///	prepares the loop over all elements of one type for the computation of the error estimator
		template <typename TElem, typename TFVGeom>
		void prep_err_est_elem_loop(const ReferenceObjectID roid, const int si)
		{ UG_THROW("Not Implemented!"); }

	///	prepares the element for assembling the error estimator
		template <typename TElem, typename TFVGeom>
		void prep_err_est_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{ UG_THROW("Not Implemented!"); }

	///	computes the error estimator contribution for one element
		template <typename TElem, typename TFVGeom>
		void compute_err_est_A_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale)
		{ UG_THROW("Not Implemented!"); }

	///	postprocesses the loop over all elements of one type in the computation of the error estimator
		template <typename TElem, typename TFVGeom>
		void fsh_err_est_elem_loop()
		{ UG_THROW("Not Implemented!"); }

	private:
		SmartPtr<UserData<number, dim> > m_srcData;  // source/sink: adding 'source' to rhs
		SmartPtr<UserData<number, dim> > m_snkTransportData;  // transport sink: subtracting u*sink
		MathVector<dim> m_srcCoord;

	private:
		void register_all_funcs();

		template <typename TElem, typename TFVGeom>
		void register_func();
	public:
		using NumberExport = SmartPtr<CplUserData<number, dim> >;

	protected:
		///	Export for the concentration
		SmartPtr<DataExport<number, dim> > m_exRate;

	private:
		bool m_bNonRegularGrid;
	//	bool m_bCurrElemIsHSlave;
};

} // namespace ug

#include "lagrange_dirac_source_impl.h"


#endif /*__H__UG__LIB_DISC__SPACIAL_DISCRETIZATION__ELEM_DISC__DIRAC_SOURCE__LAGRANGE_DIRAC_SOURCE__*/
