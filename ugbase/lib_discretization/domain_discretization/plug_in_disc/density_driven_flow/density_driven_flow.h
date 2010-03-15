
#ifndef DENSITY_DRIVEN_FLOW_H
#define DENSITY_DRIVEN_FLOW_H

#include "lib_grid/lib_grid.h"
#include "lib_algebra/lib_algebra.h"
#include "lib_discretization/lib_discretization.h"
#include "common/common.h"
#include "fvgeom.h"

namespace ug{


static number _porosity;

template
<	int dim,
	void Pososity_fct(number&),
	void Viscosity_fct(number&, number),
	void Density_fct(number&, number),
	void D_Density_fct(number&, number),
	void Mol_Diff_Tensor_fct(MathMatrix<dim,dim>&),
	void Permeability_Tensor_fct(MathMatrix<dim,dim>&),
	void Gravity_fct(MathVector<dim>&)>
class DensityDrivenFlowSystem {
	public:
		template <typename TElem>
		static inline void prepare_element(TElem* elem, FVElementGeometry<TElem>& geo, Domain<dim>& domain);

		template <typename TElem>
		static inline void assemble_element_JA(FVElementGeometry<TElem>& geo, number mat_values[], number u_values[], const uint num_dofs, number time=0.0);

		template <typename TElem>
		static inline void assemble_element_JM(FVElementGeometry<TElem>& geo, number mat_values[], number u_values[], const uint num_dofs, number time=0.0);

		template <typename TElem>
		static inline void assemble_element_A(FVElementGeometry<TElem>& geo, number def_values[], number u_values[], const uint num_dofs, number time=0.0);

		template <typename TElem>
		static inline void assemble_element_M(FVElementGeometry<TElem>& geo, number def_values[], number u_values[], const uint num_dofs, number time=0.0);

		template <typename TElem>
		static inline void assemble_element_f(FVElementGeometry<TElem>& geo, number def_values[], const uint num_dofs, number time=0.0);

	public:
		static number _upwind_amount;

	private:
		static inline void compute_ip_Darcy_velocity(MathVector<dim>& Darcy_vel, number c_ip, MathVector<dim>& grad_p_ip)
		{
			number s, viscosity_ip;
			MathVector<dim> vel;
			MathMatrix<dim, dim> K;

			Density_fct(s, c_ip);
			Viscosity_fct(viscosity_ip, c_ip);
			Gravity_fct(vel);
			Permeability_Tensor_fct(K);
			VecScale(vel, vel, s);
			VecSubtract(vel, vel, grad_p_ip);
			MatVecMult(Darcy_vel, K, vel);
			VecScale(Darcy_vel, Darcy_vel, 1./viscosity_ip);
		};

		template <typename TElem>
		static inline void compute_D_ip_Darcy_velocity(	const SubControlVolumeFace<TElem>& scvf,
														MathVector<dim>& Darcy_vel, MathVector<dim> D_Darcy_vel_c[], MathVector<dim> D_Darcy_vel_p[],
														number c_ip, MathVector<dim>& grad_p_ip)
		{
			const int num_co = reference_element_traits<TElem>::num_corners;
			number s, mu_ip;
			MathVector<dim> vel, gravity;
			MathVector<dim> D_vel_c[num_co], D_vel_p[num_co];
			MathMatrix<dim, dim> K;
			const SD_Values<TElem>& sdv = scvf.sdv();

			Density_fct(s, c_ip);
			Gravity_fct(gravity);
			Viscosity_fct(mu_ip, c_ip);
			Permeability_Tensor_fct(K);
			VecScale(vel, gravity, s);
			VecSubtract(vel, vel, grad_p_ip);
			MatVecMult(Darcy_vel, K, vel);
			VecScale(Darcy_vel, Darcy_vel, 1./mu_ip);

			D_Density_fct(s, c_ip);
			for(int co = 0; co < num_co; ++co)
			{
				VecScale(D_vel_c[co], gravity, s*sdv.shape(co));
				VecScale(D_vel_p[co], sdv.grad_global(co), -1);
				MatVecMult(D_Darcy_vel_c[co], K, D_vel_c[co]);
				MatVecMult(D_Darcy_vel_p[co], K, D_vel_p[co]);

				VecScale(D_Darcy_vel_c[co],D_Darcy_vel_c[co],1./mu_ip);
				VecScale(D_Darcy_vel_p[co],D_Darcy_vel_p[co],1./mu_ip);
			}

			// D_Viscosity == 0 !!!!
		};


};


template
<	int dim,
	void Pososity_fct(number&),
	void Viscosity_fct(number&, number),
	void Density_fct(number&, number),
	void D_Density_fct(number&, number),
	void Mol_Diff_Tensor_fct(MathMatrix<dim,dim>&),
	void Permeability_Tensor_fct(MathMatrix<dim,dim>&),
	void Gravity_fct(MathVector<dim>&)>
number
DensityDrivenFlowSystem< dim, Pososity_fct, Viscosity_fct, Density_fct, D_Density_fct, Mol_Diff_Tensor_fct, Permeability_Tensor_fct, Gravity_fct>::_upwind_amount = 0.0;



template
<	int dim,
	void Pososity_fct(number&),
	void Viscosity_fct(number&, number),
	void Density_fct(number&, number),
	void D_Density_fct(number&, number),
	void Mol_Diff_Tensor_fct(MathMatrix<dim,dim>&),
	void Permeability_Tensor_fct(MathMatrix<dim,dim>&),
	void Gravity_fct(MathVector<dim>&)>
template <typename TElem>
inline
void
DensityDrivenFlowSystem< dim, Pososity_fct, Viscosity_fct, Density_fct, D_Density_fct, Mol_Diff_Tensor_fct, Permeability_Tensor_fct, Gravity_fct>::
prepare_element(TElem* elem, FVElementGeometry<TElem>& geo, Domain<dim>& domain)
{
	// update Geometry for this element
	typename Domain<dim>::position_type corners[reference_element_traits<TElem>::num_corners];
	typename Domain<dim>::position_accessor_type aaPos = domain.get_position_accessor();
	for(int i = 0; i < reference_element_traits<TElem>::num_corners; ++i)
	{
		VertexBase* vert = elem->vertex(i);
		corners[i] = aaPos[vert];
	}

	geo.update(corners);

	// user function
	Pososity_fct(_porosity);
}

#define _C_ 0
#define _P_ 1
#define J(fct1, fct2, i, j) (J_values[num_dofs *((num_dofs/2)*(fct1) + i) +((num_dofs/2)*(fct2) + j)])
#define d(fct, i)    (d_values[num_dofs/2*(fct) + (i)])
#define u(fct, i)    (u_values[num_dofs/2*(fct) + (i)])

template
<	int dim,
	void Pososity_fct(number&),
	void Viscosity_fct(number&, number),
	void Density_fct(number&, number),
	void D_Density_fct(number&, number),
	void Mol_Diff_Tensor_fct(MathMatrix<dim,dim>&),
	void Permeability_Tensor_fct(MathMatrix<dim,dim>&),
	void Gravity_fct(MathVector<dim>&)>
template <typename TElem>
inline
void
DensityDrivenFlowSystem< dim, Pososity_fct, Viscosity_fct, Density_fct, D_Density_fct, Mol_Diff_Tensor_fct, Permeability_Tensor_fct, Gravity_fct>::
assemble_element_JA(FVElementGeometry<TElem>& geo, number J_values[], number u_values[], const uint num_dofs, number time)
{
	const int num_co = reference_element_traits<TElem>::num_corners;
	number flux, flux_c, flux_p;
	MathMatrix<dim,dim> D;
	MathVector<dim> grad_p_ip, grad_c_ip;
	MathVector<dim> Darcy_vel, D_Darcy_vel_c[num_co], D_Darcy_vel_p[num_co];
	number c_ip;
	MathVector<dim> Dgrad;
	for(uint i = 0; i < geo.num_scvf(); ++i)
	{
		const SubControlVolumeFace<TElem>& scvf = geo.scvf(i);

		for(uint ip = 0; ip < scvf.num_ip(); ++ip)
		{
			const SD_Values<TElem>& sdv = scvf.sdv();

			VecSet(grad_p_ip, 0.0); VecSet(grad_c_ip, 0.0);
			c_ip = 0.0;
			for(uint j = 0; j < sdv.num_sh(); ++j)
			{
				VecScaleAppend(grad_p_ip, u(_P_,j), sdv.grad_global(j));
				VecScaleAppend(grad_c_ip, u(_C_,j), sdv.grad_global(j));
				c_ip += u(_C_, j) * sdv.shape(j);
			}

			compute_D_ip_Darcy_velocity(scvf, Darcy_vel, D_Darcy_vel_c, D_Darcy_vel_p, c_ip, grad_p_ip);
			Mol_Diff_Tensor_fct(D);

			for(uint j = 0; j < sdv.num_sh(); ++j)
			{
				////////////////////////////////////
				// diffusiv term (central discretization)
				MatVecMult(Dgrad, D, sdv.grad_global(j));
				flux = VecDot(Dgrad, scvf.normal());

				J(_C_, _C_, scvf.from(), j) -= flux;
				J(_C_, _C_, scvf.to(), j) += flux;

				//J(_C_, _P_, scvf.from(), j) -= 0.0;
				//J(_C_, _P_, scvf.to(), j) += 0.0;;
				////////////////////////////////////
				// convective term
				// (upwinding_amount == 1.0 -> full upwind;
				//  upwinding_amount == 0.0 -> full central disc)

				// central part convection
				if(_upwind_amount != 1.0)
				{
					flux_c = (1.- _upwind_amount) * (sdv.shape(j) * VecDot(Darcy_vel, scvf.normal()) + c_ip * VecDot(D_Darcy_vel_c[j], scvf.normal()));
					flux_p = (1.- _upwind_amount) * (											0.0	 + c_ip * VecDot(D_Darcy_vel_p[j], scvf.normal()));

					// coupling 'from' with j  (i.e. A[from][j]) and 'to' with j (i.e. A[to][j])
					J(_C_, _C_, scvf.from(), j) += flux_c;
					J(_C_, _C_, scvf.to(),j) -= flux_c;

					J(_C_, _P_, scvf.from(), j) += flux_p;
					J(_C_, _P_, scvf.to(),j) -= flux_p;
				}
			}
			// upwind part convection
			if(_upwind_amount != 0.0)
			{
				uint up;
				flux_c = _upwind_amount * VecDot(Darcy_vel, scvf.normal());
				if(flux_c >= 0.0) up = scvf.from(); else up = scvf.to();

				for(uint j = 0; j < sdv.num_sh(); ++j)
				{
					if(j == up) flux_c = _upwind_amount * ( sdv.shape(up) * VecDot(Darcy_vel, scvf.normal()) );
					else flux_c = 0.0;
					flux_c += _upwind_amount * u(_C_, up) * VecDot(D_Darcy_vel_c[j], scvf.normal());

					J(_C_, _C_,scvf.from(), j) += flux_c;
					J(_C_, _C_,scvf.to(), j) -= flux_c;

					flux_p =  _upwind_amount * ( u(_C_, up) * VecDot(D_Darcy_vel_p[j], scvf.normal()));
					J(_P_, _P_,scvf.from(), j) += flux;
					J(_P_, _P_,scvf.to(), j) -= flux;
				}
			}


			///////////////////
			///// flow equation
			for(uint j = 0; j < sdv.num_sh(); ++j)
			{
				flux_c = VecDot(D_Darcy_vel_c[j], scvf.normal());
				flux_p = VecDot(D_Darcy_vel_p[j], scvf.normal());


				J(_P_, _C_, scvf.from(), j) += flux_c;
				J(_P_, _C_, scvf.to(), j) -= flux_c;

				J(_P_, _P_, scvf.from(), j) += flux_p;
				J(_P_, _P_, scvf.to(), j) -= flux_p;
			}
		}
	}

}


template
<	int dim,
	void Pososity_fct(number&),
	void Viscosity_fct(number&, number),
	void Density_fct(number&, number),
	void D_Density_fct(number&, number),
	void Mol_Diff_Tensor_fct(MathMatrix<dim,dim>&),
	void Permeability_Tensor_fct(MathMatrix<dim,dim>&),
	void Gravity_fct(MathVector<dim>&)>
template <typename TElem>
inline
void
DensityDrivenFlowSystem< dim, Pososity_fct, Viscosity_fct, Density_fct, D_Density_fct, Mol_Diff_Tensor_fct, Permeability_Tensor_fct, Gravity_fct>::
assemble_element_JM(FVElementGeometry<TElem>& geo, number J_values[], number u_values[], const uint num_dofs, number time)
{
	int co;
	for(uint i = 0; i < geo.num_scv(); ++i)
	{
		const SubControlVolume<TElem>& scv = geo.scv(i);

		co = scv.local_corner_id();

		J(_C_, _C_, co, co) += _porosity * scv.volume();
		//J(_C_, _P_, co, co) += 0;
		//J(_P_, _C_, co, co) += 0;
		//J(_P_, _P_, co, co) += 0;

	}
}


template
<	int dim,
	void Pososity_fct(number&),
	void Viscosity_fct(number&, number),
	void Density_fct(number&, number),
	void D_Density_fct(number&, number),
	void Mol_Diff_Tensor_fct(MathMatrix<dim,dim>&),
	void Permeability_Tensor_fct(MathMatrix<dim,dim>&),
	void Gravity_fct(MathVector<dim>&)>
template <typename TElem>
inline
void
DensityDrivenFlowSystem< dim, Pososity_fct, Viscosity_fct, Density_fct, D_Density_fct, Mol_Diff_Tensor_fct, Permeability_Tensor_fct, Gravity_fct>::
assemble_element_A(FVElementGeometry<TElem>& geo, number d_values[], number u_values[], const uint num_dofs, number time)
{
	number flux;
	MathVector<dim> grad_p_ip, grad_c_ip;
	number c_ip;
	MathMatrix<dim,dim> D;
	MathVector<dim> Dgrad_c_ip;
	MathVector<dim> Darcy_vel;
	for(uint i = 0; i < geo.num_scvf(); ++i)
	{
		const SubControlVolumeFace<TElem>& scvf = geo.scvf(i);

		for(uint ip = 0; ip < scvf.num_ip(); ++ip)
		{
			const SD_Values<TElem>& sdv = scvf.sdv();

			VecSet(grad_p_ip, 0.0); VecSet(grad_c_ip, 0.0);
			c_ip = 0.0;
			for(uint j = 0; j < sdv.num_sh(); ++j)
			{
				VecScaleAppend(grad_p_ip, u(_P_,j), sdv.grad_global(j));
				VecScaleAppend(grad_c_ip, u(_C_,j), sdv.grad_global(j));
				c_ip += u(_C_, j) * sdv.shape(j);
			}

			compute_ip_Darcy_velocity(Darcy_vel, c_ip, grad_p_ip);

			Mol_Diff_Tensor_fct(D);


			////////////////////////////////////
			// diffusiv term (central discretization)
			MatVecMult(Dgrad_c_ip, D, grad_c_ip);
			flux = VecDot(Dgrad_c_ip, scvf.normal());

			d(_C_, scvf.from()) -= flux;
			d(_C_,scvf.to()) += flux;

			////////////////////////////////////
			// convective term
			// (upwinding_amount == 1.0 -> full upwind;
			//  upwinding_amount == 0.0 -> full central disc)

			// central part convection
			if(_upwind_amount != 1.0)
			{
				flux = (1.- _upwind_amount) * c_ip * VecDot(Darcy_vel, scvf.normal());

				d(_C_,scvf.from()) += flux;
				d(_C_,scvf.to()) -= flux;

			}

			// upwind part convection
			if(_upwind_amount != 0.0)
			{
				flux = _upwind_amount * VecDot(Darcy_vel, scvf.normal());
				if(flux >= 0.0) flux *= u(_C_,scvf.from()); else flux *= u(_C_,scvf.to());
				d(_C_,scvf.from()) += flux;
				d(_C_,scvf.to()) -= flux;
			}

			flux = VecDot(Darcy_vel, scvf.normal());
			d(_P_,scvf.from()) += flux;
			d(_P_,scvf.from()) -= flux;
		}
	}
}


template
<	int dim,
	void Pososity_fct(number&),
	void Viscosity_fct(number&, number),
	void Density_fct(number&, number),
	void D_Density_fct(number&, number),
	void Mol_Diff_Tensor_fct(MathMatrix<dim,dim>&),
	void Permeability_Tensor_fct(MathMatrix<dim,dim>&),
	void Gravity_fct(MathVector<dim>&)>
template <typename TElem>
inline
void
DensityDrivenFlowSystem< dim, Pososity_fct, Viscosity_fct, Density_fct, D_Density_fct, Mol_Diff_Tensor_fct, Permeability_Tensor_fct, Gravity_fct>::
assemble_element_M(FVElementGeometry<TElem>& geo, number d_values[], number u_values[], const uint num_dofs, number time)
{
	int co;
	for(uint i = 0; i < geo.num_scv(); ++i)
	{
		const SubControlVolume<TElem>& scv = geo.scv(i);

		co = scv.local_corner_id();

		d(_C_,co) += _porosity * u_values[co] * scv.volume();
		d(_P_,co) += _porosity * scv.volume();
	}
}


template
<	int dim,
	void Pososity_fct(number&),
	void Viscosity_fct(number&, number),
	void Density_fct(number&, number),
	void D_Density_fct(number&, number),
	void Mol_Diff_Tensor_fct(MathMatrix<dim,dim>&),
	void Permeability_Tensor_fct(MathMatrix<dim,dim>&),
	void Gravity_fct(MathVector<dim>&)>
template <typename TElem>
inline
void
DensityDrivenFlowSystem< dim, Pososity_fct, Viscosity_fct, Density_fct, D_Density_fct, Mol_Diff_Tensor_fct, Permeability_Tensor_fct, Gravity_fct>::
assemble_element_f(FVElementGeometry<TElem>& geo, number d_values[], const uint num_dofs, number time)
{
}

}

#endif
