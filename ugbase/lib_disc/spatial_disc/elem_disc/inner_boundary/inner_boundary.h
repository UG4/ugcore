/*
 * inner_boundary.h
 *
 *  Created on: 26.02.2010
 *      Author: markusbreit
 */

#ifndef __H__UG__LIB_DISC__SPACIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__INNER_BOUNDARY__
#define __H__UG__LIB_DISC__SPACIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__INNER_BOUNDARY__

#include <boost/function.hpp>
#include <string>

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/ip_data/data_import_export.h"


/// Finite Volume Element Discretization for an inner BndCond that depends on the unknowns (on the bnd)
/**
 * This class implements the IElemDisc interface to provide element local
 * assemblings for the unknown-dependent Neumann-flux over an inner boundary.
 * The equation of this flux should be given on the script level.
 * 
 * \tparam	TDomain		Domain
 * \tparam	TAlgebra	Algebra
 */

namespace ug{

template<typename TDomain>
class FV1InnerBoundaryElemDisc
: public IDomainElemDisc<TDomain>
{
	private:
	///	Base class type
		typedef IDomainElemDisc<TDomain> base_type;

	///	own type
		typedef FV1InnerBoundaryElemDisc<TDomain> this_type;

	public:
	///	Domain type
		typedef typename base_type::domain_type domain_type;

	///	World dimension
		static const int dim = base_type::dim;

	///	Position type
		typedef typename base_type::position_type position_type;

	public:
	
	/// Constructor
		FV1InnerBoundaryElemDisc(size_t numFct, const char* functions, const char* subsets)
			: IDomainElemDisc<TDomain>(numFct, functions, subsets), _numFct(numFct)
		{
			register_all_fv1_funcs();
		}
	
	/// Setting the flux function
		void set_fluxFunction(IPData<number, dim>& fluxFct) {m_fluxFct.set_data(fluxFct);}
	
	public:	// inherited from IElemDisc
	///	type of trial space for each function used
		virtual bool request_finite_element_id(const std::vector<LFEID>& vLfeID)
		{
		//	check number
			if(vLfeID.size() != _numFct) return false;

		//	check that Lagrange 1st order
			for(size_t i = 0; i < vLfeID.size(); ++i)
				if(vLfeID[i] != LFEID(LFEID::LAGRANGE, 1)) return false;
			return true;
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

	///	returns if hanging nodes are used
		virtual bool use_hanging() const {return false;}

	private:
	
	///	number of unknowns involved
		size_t _numFct;
		
	///	Data import for flux
		DataImport<number, dim> m_fluxFct;

	
	///	prepares the loop over all elements
	/**
	 * This method prepares the loop over all elements. It resizes the Position
	 * array for the corner coordinates and schedules the local ip positions
	 * at the data imports.
	 */
		template<typename TElem, template <class Elem, int Dim> class TFVGeom>
		bool prepare_element_loop();

	///	prepares the element for assembling
	/**
	 * This methods prepares an element for the assembling. The Positions of
	 * the Element Corners are read and the Finite Volume Geometry is updated.
	 * The global ip positions are scheduled at the data imports.
	 */
		template<typename TElem, template <class Elem, int Dim> class TFVGeom>
		bool prepare_element(TElem* elem, const LocalVector& u);

	///	finishes the loop over all elements
		template<typename TElem, template <class Elem, int Dim> class TFVGeom>
		bool finish_element_loop();

	///	assembles the local stiffness matrix using a finite volume scheme
		template<typename TElem, template <class Elem, int Dim> class TFVGeom>
		bool assemble_JA(LocalMatrix& J, const LocalVector& u);

	///	assembles the local mass matrix using a finite volume scheme
		template<typename TElem, template <class Elem, int Dim> class TFVGeom>
		bool assemble_JM(LocalMatrix& J, const LocalVector& u);

	///	assembles the stiffness part of the local defect
		template<typename TElem, template <class Elem, int Dim> class TFVGeom>
		bool assemble_A(LocalVector& d, const LocalVector& u);

	///	assembles the mass part of the local defect
		template<typename TElem, template <class Elem, int Dim> class TFVGeom>
		bool assemble_M(LocalVector& d, const LocalVector& u);

	///	assembles the local right hand side
		template<typename TElem, template <class Elem, int Dim> class TFVGeom>
		bool assemble_f(LocalVector& d);

	private:
		// position access
		const position_type* m_vCornerCoords;

	private:
		void register_all_fv1_funcs();

		//template <template <class Elem, int WorldDim> class TFVGeom>
		template <template <class Elem, int WorldDim> class TFVGeom>
		struct RegisterFV1 {
				RegisterFV1(this_type* pThis) : m_pThis(pThis){}
				this_type* m_pThis;
				template< typename TElem > void operator()(TElem&)
				{m_pThis->register_fv1_func<TElem, TFVGeom>();}
		};

		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		void register_fv1_func();

};

}

#endif /*__H__UG__LIB_DISC__SPACIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__INNER_BOUNDARY__*/
