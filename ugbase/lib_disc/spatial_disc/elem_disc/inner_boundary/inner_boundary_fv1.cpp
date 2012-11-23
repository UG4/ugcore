/*
 * inner_boundary_impl.h
 *
 *  Created on: 14.10.2010
 *      Author: markusbreit
 */

#ifndef __H__UG__LIB_DISC__SPACIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__INNER_BOUNDARY_IMPL__
#define __H__UG__LIB_DISC__SPACIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__INNER_BOUNDARY_IMPL__

#include "inner_boundary.h"
#include "common/util/provider.h"
#include "lib_disc/spatial_disc/disc_util/finite_volume_geometry.h"



namespace ug
{


template<typename TDomain>
template<typename TElem, template <class Elem, int Dim> class TFVGeom>
void FV1InnerBoundaryElemDisc<TDomain>::
prepare_element_loop()
{}


template<typename TDomain>
template<typename TElem, template <class Elem, int Dim> class TFVGeom>
void FV1InnerBoundaryElemDisc<TDomain>::
finish_element_loop()
{}


template<typename TDomain>
template<typename TElem, template <class Elem, int Dim> class TFVGeom>
void FV1InnerBoundaryElemDisc<TDomain>::
prepare_element(TElem* elem, const LocalVector& u)
{
//	get corners
	m_vCornerCoords = this->template element_corners<TElem>(elem);

	// update Geometry for this element
	TFVGeom<TElem, dim>& geo = Provider<TFVGeom<TElem,dim> >::get();
	if(!geo.update(elem, &m_vCornerCoords[0], &(this->subset_handler())))
		UG_THROW("FV1InnerBoundaryElemDisc::prepare_element: "
						"Cannot update Finite Volume Geometry.\n");
}

// assemble stiffness part of Jacobian
template<typename TDomain>
template<typename TElem, template <class Elem, int Dim> class TFVGeom>
void FV1InnerBoundaryElemDisc<TDomain>::
ass_JA_elem(LocalMatrix& J, const LocalVector& u)
{
	// get finite volume geometry
	const static TFVGeom<TElem, dim>& fvgeom = Provider<TFVGeom<TElem,dim> >::get();

	for (size_t i = 0; i < fvgeom.num_bf(); ++i)
	{
		// get current BF
		const typename TFVGeom<TElem, dim>::BF& bf = fvgeom.bf(i);

		// get associated node
		const int co = bf.node_id();
		
		FluxDerivCond fdc;
		if (!fluxDensityDerivFct(u, co, fdc))
			UG_THROW("FV1InnerBoundaryElemDisc::ass_JA_elem:"
							" Call to fluxDensityDerivFct resulted did not succeed.");
		
		// scale with volume of BF
		for (size_t j=0; j<fdc.fluxDeriv.size(); j++)
			for (size_t k=0; k<fdc.fluxDeriv[j].size(); k++)
				fdc.fluxDeriv[j][k] *= bf.volume();
		
		// add to Jacobian
		for (size_t j=0; j<fdc.fluxDeriv.size(); j++)
		{
			for (size_t k=0; k<fdc.fluxDeriv[j].size(); k++)
			{
				J(fdc.from[j],co,k,co)	+= fdc.fluxDeriv[j][k];
				J(fdc.to[j],co,k,co)	-= fdc.fluxDeriv[j][k];
			}
		}	
	}
}

template<typename TDomain>
template<typename TElem, template <class Elem, int Dim> class TFVGeom>
void FV1InnerBoundaryElemDisc<TDomain>::
ass_JM_elem(LocalMatrix& J, const LocalVector& u)
{}


template<typename TDomain>
template<typename TElem, template <class Elem, int Dim> class TFVGeom>
void FV1InnerBoundaryElemDisc<TDomain>::
ass_dA_elem(LocalVector& d, const LocalVector& u)
{
	// get finite volume geometry
	static TFVGeom<TElem, dim>& fvgeom = Provider<TFVGeom<TElem,dim> >::get();

	// loop Boundary Faces
	for (size_t i = 0; i < fvgeom.num_bf(); ++i)
	{
		// get current BF
		const typename TFVGeom<TElem, dim>::BF& bf = fvgeom.bf(i);
		
		// get associated node
		const int co = bf.node_id();
		
		// get flux densities in that node
		FluxCond fc;
		if (!fluxDensityFct(u, co, fc))
			UG_THROW("FV1InnerBoundaryElemDisc::ass_dA_elem:"
						" Call to fluxDensityFct did not succeed.");

		/*if (fc.flux[0] < -1.0e-7)
		{
			//const MathVector<dim>& corners = bf.global_corner(co);
			//std::cout << "flux density: " << fc.flux[0]/bf.volume() << ";   " << corners << std::endl;
			number caCyt = u(0, co);	// cytosolic Ca2+ concentration
			number caER = u(1, co);		// ER Ca2+ concentration
			std::cout << "flux density: " << fc.flux[0] << "   cyt: " << caCyt << "   er: " << caER << std::endl;
		}*/

		// scale with volume of BF
		for (size_t j=0; j<fc.flux.size(); j++) fc.flux[j] *= bf.volume();
		
		// add to defect
		for (size_t j=0; j<fc.flux.size(); j++)
		{
			d(fc.from[j], co) += fc.flux[j];
			d(fc.to[j], co) -= fc.flux[j];
		}
	}
}

template<typename TDomain>
template<typename TElem, template <class Elem, int Dim> class TFVGeom>
void FV1InnerBoundaryElemDisc<TDomain>::
ass_dM_elem(LocalVector& d, const LocalVector& u)
{}

template<typename TDomain>
template<typename TElem, template <class Elem, int Dim> class TFVGeom>
void FV1InnerBoundaryElemDisc<TDomain>::
ass_rhs_elem(LocalVector& d)
{}


////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

// register for 1D
template<typename TDomain>
void
FV1InnerBoundaryElemDisc<TDomain>::
register_all_fv1_funcs()
{
//	get all grid element types in this dimension and below
	typedef typename domain_traits<dim>::ManifoldElemList ElemList;

//	switch assemble functions
	boost::mpl::for_each<ElemList>(RegisterFV1<FV1ManifoldBoundary>(this));
}

template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void
FV1InnerBoundaryElemDisc<TDomain>::
register_fv1_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	this->enable_fast_ass_elem(true);
	this->set_prep_elem_loop_fct(id, &T::template prepare_element_loop<TElem, TFVGeom>);
	this->set_prep_elem_fct(	 id, &T::template prepare_element<TElem, TFVGeom>);
	this->set_fsh_elem_loop_fct( id, &T::template finish_element_loop<TElem, TFVGeom>);
	this->set_ass_JA_elem_fct(	 id, &T::template ass_JA_elem<TElem, TFVGeom>);
	this->set_ass_JM_elem_fct(	 id, &T::template ass_JM_elem<TElem, TFVGeom>);
	this->set_ass_dA_elem_fct(	 id, &T::template ass_dA_elem<TElem, TFVGeom>);
	this->set_ass_dM_elem_fct(	 id, &T::template ass_dM_elem<TElem, TFVGeom>);
	this->set_ass_rhs_elem_fct(	 id, &T::template ass_rhs_elem<TElem, TFVGeom>);
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

template class FV1InnerBoundaryElemDisc<Domain1d>;
template class FV1InnerBoundaryElemDisc<Domain2d>;
template class FV1InnerBoundaryElemDisc<Domain3d>;



} // namespace ug


#endif /*__H__UG__LIB_DISC__SPACIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__INNER_BOUNDARY_IMPL__*/
