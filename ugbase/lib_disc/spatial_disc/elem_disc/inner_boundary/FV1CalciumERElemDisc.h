/*
 *  FV1CalciumERElemDisc.h
 *
 *  Created on: 20.12.2011
 *      Author: markusbreit
 */

#ifndef __H__UG__LIB_DISC__SPACIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__INNER_BOUNDARY__CALCIUM_ER__
#define __H__UG__LIB_DISC__SPACIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__INNER_BOUNDARY__CALCIUM_ER__


#include "inner_boundary.h"


/// Finite Volume Element Discretization for the inner BndCond on an ER membrane
/**
 * This class implements the InnerBoundary interface to provide element local
 * assemblings for the unknown-dependent Calcium-Neumann-flux over the ER.
 */


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



namespace ug
{

template<typename TDomain>
class FV1CalciumERElemDisc
: public FV1InnerBoundaryElemDisc<TDomain>
{
	public:
		
		FV1CalciumERElemDisc(const char* functions, const char* subsets)
			: FV1InnerBoundaryElemDisc<TDomain>(functions, subsets){};
	
	private:
		virtual bool fluxDensityFct(const LocalVector& u, size_t node_id, FluxCond& fc)
		{
			if (u.num_fct() < 3)
			{
				UG_LOG( "ERROR in 'FV1InnerBoundaryElemDisc:fluxFct':"
						" LocalVector u does not have enough functions"
						" (has " << u.num_fct() << ", but needs at least 3).");
				return false;
			}
			
			number flux = 0.0;
		
			number caCyt = u(0, node_id);	// cytosolic Ca2+ concentration
			number caER = u(1, node_id);		// ER Ca2+ concentration
			number ip3 = u(2, node_id);		// IP3 concentration
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
			
			fc.flux.resize(1);	fc.flux[0] = flux;
			fc.from.resize(1);	fc.from[0] = 1;		// ER
			fc.to.resize(1);	fc.to[0] = 0;		// cytosol
			
			return true;
		}
		
		virtual bool fluxDensityDerivFct(const LocalVector& u, size_t node_id, FluxDerivCond& fdc)
		{
			// get values of the unknowns in associated node
			number caCyt = u(0, node_id);
			number caER = u(1, node_id);
			number ip3 = u(2, node_id);
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
			
			// add to Jacobian
			fdc.fluxDeriv.resize(3);
			for (size_t i=0; i<3; i++) fdc.fluxDeriv[i].resize(1);
			
			fdc.fluxDeriv[0][0] = d_dCyt;
			fdc.fluxDeriv[1][0] = d_dER;
			fdc.fluxDeriv[2][0] = d_dIP3;
			
			fdc.from.resize(1);	fdc.from[0] = 1;		// ER
			fdc.to.resize(1);	fdc.to[0] = 0;		// cytosol
			
			return true;
		}

};


////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

template class FV1CalciumERElemDisc<Domain1d>;
template class FV1CalciumERElemDisc<Domain2d>;
template class FV1CalciumERElemDisc<Domain3d>;

}


#endif /*__H__UG__LIB_DISC__SPACIAL_DISCRETIZATION__ELEM_DISC__NEUMANN_BOUNDARY__FV1__INNER_BOUNDARY__CALCIUM_ER__*/
