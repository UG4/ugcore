/*
 * neumann_boundary_impl.h
 *
 *  Created on: 14.10.2010
 *      Author: markusbreit
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__INNER_BOUNDARY_IMPL__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__INNER_BOUNDARY_IMPL__

#include "inner_boundary.h"


// constants (all according to De Young & Keizer)
#define C1 0.185		// (ER vol) / (cytosolic vol)

						// receptor dissociation
#define D1 1.3e-7		// IP3
#define D2 1.049e-5	//1.049e-4	//1.049e-6		// Ca2+ inhibition
#define D3 9.434e-7		// IP3
#define D5 8.234e-7	//6.234e-6	//8.234e-8		// Ca2+ activation

#define NY1 6.0			// max Ca2+ channel flux factor (when all channels are open)
#define NY2 0.0	//0.0011 //0.011	//0.11		// leakage flux
#define NY3 9.0e-7		// max Ca2+ uptake (SERCA)

#define K3 1.0e-7		// activation constant for ATP-Ca2+ pump

// RyR
#define NY4 0.15
#define KA4 1.92e-26
#define KB3 2.573e-19
#define KC  7.1e-3	//5.71e-2 

#define KA4I 1.92e-2
#define KB3I 2.573e-1

#define V_SERCA 1.2e-10
#define K_SERCA 1.8e-7

#define DAMPENING 1.0



namespace ug{




template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int Dim> class TFVGeom>
bool
FVInnerBoundaryElemDisc<TDomain, TAlgebra>::
prepare_element_loop()
{
//	we're done
	return true;
}


template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int Dim> class TFVGeom>
inline
bool
FVInnerBoundaryElemDisc<TDomain, TAlgebra>::
finish_element_loop()
{
//	we're done
	return true;
}


template<typename TDomain, typename TAlgebra>

template<typename TElem, template <class Elem, int Dim> class TFVGeom>
inline
bool
FVInnerBoundaryElemDisc<TDomain, TAlgebra>::
prepare_element(TElem* elem, const local_vector_type& u, const local_index_type& glob_ind)
{
//	get corners
	m_vCornerCoords = this->template get_element_corners<TElem>(elem);

	// update Geometry for this element
	TFVGeom<TElem, dim>& geo = GeomProvider::get<TFVGeom<TElem, dim> >();
	if(!geo.update(elem, this->get_subset_handler(), &m_vCornerCoords[0]))
	{
		UG_LOG("ERROR in 'FVInnerBoundaryElemDisc::prepare_element: "
				"Cannot update Finite Volume Geometry.\n"); return false;}

	return true;
}



// assemble stiffness part of Jacobian
template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int Dim> class TFVGeom>
inline
bool
FVInnerBoundaryElemDisc<TDomain, TAlgebra>::
assemble_JA(local_matrix_type& J, const local_vector_type& u)
{
	// get finite volume geometry
	const static TFVGeom<TElem, dim>& fvgeom = GeomProvider::get<TFVGeom<TElem,dim> >();

	for (size_t i = 0; i < fvgeom.num_bf(); ++i)
	{
		// get current BF
		const typename TFVGeom<TElem, dim>::BF& bf = fvgeom.bf(i);

		// get associated node
		const int co = bf.node_id();
		
		// get values of the unknowns in associated node
		number caCyt = u(0, co);
		number caER = u(1, co);
		number ip3 = u(2, co);
		//INT n_surf;
		//surf_id = BNDP_SurfaceId(V_BNDP(vert_to), &n_surf, 0);
		
		
		number d_dCyt = 0.0;
		number d_dER = 0.0;
		number d_dIP3 = 0.0;

		
		// IP3 parts
		number schlonz1 = caCyt*ip3 + ip3*D2 + D1*D2 + caCyt*D3;
		number schlonz = schlonz1 * (caCyt+D5);
		number x110 = (caCyt*ip3*D2) / schlonz;
		
		number blubb1 = NY1*3.0*x110*x110*ip3*D2 * (1.0 - caCyt/schlonz * ( (ip3+D3)*(caCyt+D5) + schlonz1 )) / schlonz;
		number blubb1a = NY1*3.0*x110*x110*caCyt*D2 * (1.0 - ip3/schlonz * (caCyt+D2)*(caCyt+D5)) / schlonz; 
		number blubb2 = NY1*x110*x110*x110 + NY2;
		
		number release_factor = 1.0; //ip3_channel_density(surf_id);
		
		d_dCyt += release_factor * C1 * (blubb1 * (caER-caCyt) - blubb2);
		d_dER += release_factor * C1 * blubb2;
		d_dIP3 += release_factor * C1 * blubb1a * (caER-caCyt);
		
	
		// RyR parts
		number inner = 1.0 + pow(caCyt,3) / KB3;
		number dInner_dCyt = 3.0 * caCyt * caCyt / KB3;
		number dKlShice_dCyt = dInner_dCyt - 4.0 * KA4 / pow(caCyt,5);
		number dKlShice_dCyt_a = 3.0 * inner * inner / KB3I - 4.0 * KA4I / pow(inner,5);
		number klShice = 1.0 + KA4 / pow(caCyt,4) + pow(caCyt,3) / KB3;
		number klShice_a = 1.0 + KA4I / pow(inner,4) + pow(inner,3) / KB3I;
		number grShice_a = 1.0 + KC * klShice_a;
		number w = KC * klShice_a / grShice_a;
		number open_prob =  w / klShice;
		number dOpenProb_dCyt = (dInner_dCyt * dKlShice_dCyt_a / grShice_a - klShice_a / klShice * dKlShice_dCyt) 
								* KC / (grShice_a * klShice);
		
		release_factor = 1.0; //ryr_channel_density(surf_id);
		
		d_dCyt += release_factor * NY4 * (dOpenProb_dCyt * (caER-caCyt) - open_prob);
		d_dER  += release_factor * NY4 * open_prob;
	
	
		// SERCA parts
		d_dCyt -= V_SERCA*K_SERCA / (caER*(K_SERCA+caCyt)*(K_SERCA+caCyt));	//	2.0*NY3*K3*K3 * caCyt / ((caCyt*caCyt + K3*K3)*(caCyt*caCyt + K3*K3));
		d_dER  -= - V_SERCA*caCyt / ((K_SERCA+caCyt)*caER*caER);
	
	
		// scale with volume of BF
		d_dCyt *= bf.volume() * DAMPENING;
		d_dER *= bf.volume() * DAMPENING;
		d_dIP3 *= bf.volume() * DAMPENING;
		
		// add to Jacobian					<-- hier entsteht der segfault
		J(1,co,1,co)	+= d_dER;
		J(1,co,0,co)	+= d_dCyt;
		J(0,co,1,co)	-= d_dER;
		J(0,co,0,co)	-= d_dCyt;
		J(1,co,2,co)	+= d_dIP3;
		J(0,co,2,co)	-= d_dIP3;
	}
	
	return true;
}



template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int Dim> class TFVGeom>
inline
bool
FVInnerBoundaryElemDisc<TDomain, TAlgebra>::
assemble_JM(local_matrix_type& J, const local_vector_type& u)
{
	// nothing to be done
	return true;
}



template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int Dim> class TFVGeom>
inline
bool
FVInnerBoundaryElemDisc<TDomain, TAlgebra>::
assemble_A(local_vector_type& d, const local_vector_type& u)
{
	// get finite volume geometry
	static TFVGeom<TElem, dim>& fvgeom = GeomProvider::get<TFVGeom<TElem,dim> >();

	// loop Boundary Faces
	for (size_t i = 0; i < fvgeom.num_bf(); ++i)
	{
		// get current BF
		const typename TFVGeom<TElem, dim>::BF& bf = fvgeom.bf(i);
		
		// get associated node
		const int co = bf.node_id();
		
		
		number flux = 0.0;
		
		number caCyt = u(0, co);	// cytosolic Ca2+ concentration
		number caER = u(1, co);		// ER Ca2+ concentration
		number ip3 = u(2, co);		// IP3 concentration
		//INT n_surf;
		//INT surf_id = BNDP_SurfaceId(V_BNDP(vert_to), &n_surf, 0);
		
		
		
		// IP3 parts
		number x110 = (caCyt*ip3*D2) / ((caCyt*ip3 + ip3*D2 + D1*D2 + caCyt*D3) * (caCyt+D5));
		number IP3flux = C1*(NY1*x110*x110*x110 + NY2) * (caER-caCyt);
		
		number release_factor = 1.0; //ip3_channel_density(surf_id);
			
		flux += release_factor * IP3flux;
		
		// RyR parts
		number inner = 1.0 + pow(caCyt,3)/KB3;
		number w = (1.0 + KA4I/pow(inner,4) + pow(inner,3)/KB3I) / (1.0 + 1.0/KC + KA4I/pow(inner,4) + pow(inner,3)/KB3I);
		number open_prob =  w / (1.0 + KA4/pow(caCyt,4) + pow(caCyt,3)/KB3);	//printf("open_prob: %g\n", open_prob);
		
		release_factor = 1.0; //ryr_channel_density(surf_id);

		flux += release_factor * NY4 * open_prob * (caER-caCyt);

		// SERCA parts
		flux -=	V_SERCA*caCyt / ((K_SERCA + caCyt) * caER);
		
		
		
		// scale with volume of BF
		flux *= bf.volume() * DAMPENING;

		// add to defect					<-- hier entsteht der segfault
		d(0, co) -= flux;	// cytosol
		d(1, co) += flux;	// ER
	}

	return true;
}



template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int Dim> class TFVGeom>
inline
bool
FVInnerBoundaryElemDisc<TDomain, TAlgebra>::
assemble_M(local_vector_type& d, const local_vector_type& u)
{
	// nothing to be done
	return true;
}



template<typename TDomain, typename TAlgebra>
template<typename TElem, template <class Elem, int Dim> class TFVGeom>
inline
bool
FVInnerBoundaryElemDisc<TDomain, TAlgebra>::
assemble_f(local_vector_type& d)
{
	// nothing to be done
	return true;
}


} // namespace ug


#endif /*__H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__INNER_BOUNDARY_IMPL__*/
