/*
 * inner_boundary.h
 *
 * Finite Volume Element Discretization for an inner BndCond that depends on the unknowns (on the bnd)
 *
 * This class implements the IElemDisc interface to provide element local
 * assemblings for the unknown-dependent Neumann-flux over an inner boundary.
 * The equation of this flux should be given in a concretization of this class.
 * 
 * \tparam	TDomain		Domain
 * \tparam	TAlgebra	Algebra
 * 
 *  Created on: 26.02.2010
 *      Author: markusbreit
 */

#ifndef __H__UG__LIB_DISC__SPACIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__INNER_BOUNDARY__
#define __H__UG__LIB_DISC__SPACIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__INNER_BOUNDARY__

#include <boost/function.hpp>
#include <vector>
#include <string>

// other ug4 modules
#include "common/common.h"

// library intern headers
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"



namespace ug
{

/// struct that holds information about the flux densities and from where to where the flux occurs
struct FluxCond
{
	// vector of fluxFctValues
	std::vector<number> flux;
	std::vector<size_t> from;
	std::vector<size_t> to;
};

struct FluxDerivCond
{
	// vector of fluxFctDerivValues (wrt fct, flux number)
	std::vector<std::vector<number> > fluxDeriv;
	std::vector<size_t> from;
	std::vector<size_t> to;
};


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
        FV1InnerBoundaryElemDisc(const char* functions, const char* subsets)
        	: IDomainElemDisc<TDomain>(functions, subsets)
        {
        	register_all_fv1_funcs();
        }

    /// Setting the flux function
        //void set_fluxFunction(UserData<number, dim>& fluxFct) {m_fluxFct.set_data(fluxFct);}
	
	public:	// inherited from IElemDisc
	///	type of trial space for each function used
		virtual bool request_finite_element_id(const std::vector<LFEID>& vLfeID)
		{
		//	check number
			if(vLfeID.size() != this->num_fct()) return false;

		//	check that Lagrange 1st order
			for(size_t i = 0; i < vLfeID.size(); ++i)
				if(vLfeID[i] != LFEID(LFEID::LAGRANGE, 1)) return false;
			return true;
		}

	///	switches between non-regular and regular grids
		virtual bool request_non_regular_grid(bool bNonRegular)
		{
		//	switch, which assemble functions to use.
			if(bNonRegular)
			{
				UG_LOG("ERROR in 'FV1InnerBoundaryElemDisc::request_non_regular_grid':"
						" Non-regular grid not implemented.\n");
				return false;
			}

		//	this disc supports regular grids
			return true;
		}

	///	returns if hanging nodes are used
		virtual bool use_hanging() const {return false;}

	private:
	
	/// the flux function
	/**	This is the actual flux function defining the flux density over the boundary
	 *	depending on the unknowns on the boundary;
	 *	shall be defined in a specialized class that is derived from FV1InnerBoundaryElemDisc.
	 */
		virtual bool fluxDensityFct(const LocalVector& u, size_t node_id, FluxCond& fc) = 0;	/// the flux function
	
	/**	This is the flux derivative function defining the flux density derivatives over the boundary
	 *	depending on the unknowns on the boundary;
	 *	shall be defined in a specialized class that is derived from FV1InnerBoundaryElemDisc.
	 */
		virtual bool fluxDensityDerivFct(const LocalVector& u, size_t node_id, FluxDerivCond& fdc) = 0;
	
	///	prepares the loop over all elements
	/**
	 * This method prepares the loop over all elements. It resizes the Position
	 * array for the corner coordinates and schedules the local ip positions
	 * at the data imports.
	 */
		template<typename TElem, template <class Elem, int Dim> class TFVGeom>
		void prepare_element_loop();

	///	prepares the element for assembling
	/**
	 * This methods prepares an element for the assembling. The Positions of
	 * the Element Corners are read and the Finite Volume Geometry is updated.
	 * The global ip positions are scheduled at the data imports.
	 */
		template<typename TElem, template <class Elem, int Dim> class TFVGeom>
		void prepare_element(TElem* elem, const LocalVector& u);

	///	finishes the loop over all elements
		template<typename TElem, template <class Elem, int Dim> class TFVGeom>
		void finish_element_loop();

	///	assembles the local stiffness matrix using a finite volume scheme
		template<typename TElem, template <class Elem, int Dim> class TFVGeom>
		void ass_JA_elem(LocalMatrix& J, const LocalVector& u);

	///	assembles the local mass matrix using a finite volume scheme
		template<typename TElem, template <class Elem, int Dim> class TFVGeom>
		void ass_JM_elem(LocalMatrix& J, const LocalVector& u);

	///	assembles the stiffness part of the local defect
		template<typename TElem, template <class Elem, int Dim> class TFVGeom>
		void ass_dA_elem(LocalVector& d, const LocalVector& u);

	///	assembles the mass part of the local defect
		template<typename TElem, template <class Elem, int Dim> class TFVGeom>
		void ass_dM_elem(LocalVector& d, const LocalVector& u);

	///	assembles the local right hand side
		template<typename TElem, template <class Elem, int Dim> class TFVGeom>
		void ass_rhs_elem(LocalVector& d);

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
