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

template<typename TDomain>
class FV1InnerBoundaryElemDisc
: public IElemDisc<TDomain>
{
	public:
		/// struct that holds information about the flux densities and from where to where the flux occurs
		struct FluxCond
		{
			// vector of fluxFctValues
			std::vector<number> flux;
			std::vector<size_t> from;
			std::vector<size_t> to;
		};

		/// struct that holds information about the derivatives of the flux densities
		/// and from where to where the flux occurs
		struct FluxDerivCond
		{
			// vector of fluxFctDerivValues (wrt fct, flux number)
			std::vector<std::vector<number> > fluxDeriv;
			std::vector<size_t> from;
			std::vector<size_t> to;
		};

	private:
	///	Base class type
		typedef IElemDisc<TDomain> base_type;

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
        	: IElemDisc<TDomain>(functions, subsets)
        {
        	register_all_fv1_funcs();
        }

    /// Setting the flux function
        //void set_fluxFunction(UserData<number, dim>& fluxFct) {m_fluxFct.set_data(fluxFct);}
	
	public:	// inherited from IElemDisc
	///	type of trial space for each function used
		virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
		{
			if(bNonRegularGrid)
				UG_THROW("FV1InnerBoundary: only regular grid implemented.");

		//	check that Lagrange 1st order
			for(size_t i = 0; i < vLfeID.size(); ++i)
				if(vLfeID[i] != LFEID(LFEID::LAGRANGE, 1))
					UG_THROW("FV1InnerBoundary: 1st order lagrange expected.");
		}

	///	returns if hanging nodes are used
		virtual bool use_hanging() const {return false;}

	private:
	
	/// the flux function
	/**	This is the actual flux function defining the flux density over the boundary
	 *	depending on the unknowns on the boundary;
	 *	shall be defined in a specialized class that is derived from FV1InnerBoundaryElemDisc.
	 */
		virtual bool fluxDensityFct(const LocalVector& u, size_t node_id, int si, FluxCond& fc) = 0;
	
	/**	This is the flux derivative function defining the flux density derivatives over the boundary
	 *	depending on the unknowns on the boundary;
	 *	shall be defined in a specialized class that is derived from FV1InnerBoundaryElemDisc.
	 */
		virtual bool fluxDensityDerivFct(const LocalVector& u, size_t node_id, int si, FluxDerivCond& fdc) = 0;
	
	///	prepares the loop over all elements
	/**
	 * This method prepares the loop over all elements. It resizes the Position
	 * array for the corner coordinates and schedules the local ip positions
	 * at the data imports.
	 */
		template<typename TElem, template <class Elem, int Dim> class TFVGeom>
		void prep_elem_loop(const ReferenceObjectID roid, const int si);

	///	prepares the element for assembling
	/**
	 * This methods prepares an element for the assembling. The Positions of
	 * the Element Corners are read and the Finite Volume Geometry is updated.
	 * The global ip positions are scheduled at the data imports.
	 */
		template<typename TElem, template <class Elem, int Dim> class TFVGeom>
		void prep_elem(const LocalVector& u, GeometricObject* elem, const MathVector<dim> vCornerCoords[]);

	///	finishes the loop over all elements
		template<typename TElem, template <class Elem, int Dim> class TFVGeom>
		void fsh_elem_loop();

	///	assembles the local stiffness matrix using a finite volume scheme
		template<typename TElem, template <class Elem, int Dim> class TFVGeom>
		void add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GeometricObject* elem, const MathVector<dim> vCornerCoords[]);

	///	assembles the local mass matrix using a finite volume scheme
		template<typename TElem, template <class Elem, int Dim> class TFVGeom>
		void add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GeometricObject* elem, const MathVector<dim> vCornerCoords[]);

	///	assembles the stiffness part of the local defect
		template<typename TElem, template <class Elem, int Dim> class TFVGeom>
		void add_def_A_elem(LocalVector& d, const LocalVector& u, GeometricObject* elem, const MathVector<dim> vCornerCoords[]);

	///	assembles the mass part of the local defect
		template<typename TElem, template <class Elem, int Dim> class TFVGeom>
		void add_def_M_elem(LocalVector& d, const LocalVector& u, GeometricObject* elem, const MathVector<dim> vCornerCoords[]);

	///	assembles the local right hand side
		template<typename TElem, template <class Elem, int Dim> class TFVGeom>
		void add_rhs_elem(LocalVector& rhs, GeometricObject* elem, const MathVector<dim> vCornerCoords[]);

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
