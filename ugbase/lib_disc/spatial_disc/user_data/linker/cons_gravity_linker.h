/*
 * Copyright (c) 2026:  CEMSE, KAUST
 * Author: Dmitry Logashenko
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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__CONSISTENT_GRAVITY_LINKER__
#define __H__UG__LIB_DISC__SPATIAL_DISC__CONSISTENT_GRAVITY_LINKER__

#include <vector>

#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/spatial_disc/disc_util/consistent_gravity.h"
#include "linker.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{


////////////////////////////////////////////////////////////////////////////////
// Consistent Gravity linker
////////////////////////////////////////////////////////////////////////////////

/// Linker for the consistent gravity (according to P. Frolkovic)
/**
 * This linker computes the s.c. consistent gravity \f$[ \rho \mathbf{g}  ]_{consistent} \f$, according to Frolkovic,
 * where
 * <ul>
 * <li> \f$ \rho \f$    the density
 * <li> \f$ \mathbf{g} \f$  constant gravity
 * </ul>
 * is the input parameter. ( \f$ \mathbf{g} \f$ should be constant.)
 *
 * References:
 * <ul>
 *  <li> P. Frolkovic, P. Knabner, Consistent Velocity Approximations in Finite
 *       Element or Volume Discretizations of Density Driven Flow, In: Computational
 *       Methods in WaterResources XI, Vol. 1 (A.A. Aldama et al., eds.),
 *       Computational Mechanics Publication, Southhampten, 1996, p. 93-100
 *  </li>
 * </ul>
 */
template <int dim>
class ConsistentGravityLinker
	: public StdDataLinker<ConsistentGravityLinker<dim>, MathVector<dim>, dim>
{
	///	Base class type
		typedef StdDataLinker<ConsistentGravityLinker<dim>, MathVector<dim>, dim> base_type;
	
	/// The world dimension
		static const int world_dim = dim;
		
	public:
		ConsistentGravityLinker() :
			m_spDensity(NULL), m_spDDensity(NULL),
			m_gravity(0), m_noDerivatives(false)
		{
			this->set_num_input(1);
		}


		inline void evaluate (MathVector<dim>& value,
		                      const MathVector<dim>& globIP,
		                      number time, int si) const
		{
			UG_THROW("ConsistentGravityLinker: Element required for evaluation.")
		}

		template <int refDim>
		inline void evaluate(MathVector<dim> vValue[],
		                     const MathVector<dim> vGlobIP[],
		                     number time, int si,
		                     GridObject* elem,
		                     const MathVector<dim> vCornerCoords[],
		                     const MathVector<refDim> vLocalIP[],
		                     const size_t nip,
		                     LocalVector* u,
		                     const MathMatrix<refDim, dim>* vJT = NULL) const
		{
			const ReferenceObjectID roid = elem->reference_object_id();
			const DimReferenceElement<refDim>& refElem = ReferenceElementProvider::get<refDim>(roid);
			const size_t nco = refElem.num(ROID_VERTEX);
						
			const LocalShapeFunctionSet<refDim>& trialSpace =
				LocalFiniteElementProvider::get<refDim>(roid, LFEID(LFEID::LAGRANGE, refDim, 1));
	
		// consistent gravity
			StdLinConsistentGravity<refDim> ConsGravityMethod;
			std::vector<MathVector<refDim> > vConsGravity(nco);
			std::vector<MathVector<refDim> > vLocalGrad(nco);
			
		// coefficients at the integration points
			MathVector<dim> gravity; // we assume that the gravity is constant in the elem
			std::vector<number> vDensity(nco);

			(*m_spDensity)(&(vDensity[0]), vCornerCoords, time, si,
							elem, vCornerCoords, refElem.corners(), nco, u, NULL);
			
			try
			{
				ConsGravityMethod.template prepare<dim>
					(&(vConsGravity[0]), nco, vCornerCoords, &(vDensity[0]), m_gravity);
			}
			UG_CATCH_THROW ("ConsistentDarcyVelLinker: Cannot prepare the consistent gravity.");
			
		// get jacobians of the transformation mapping if not passed
			std::vector<MathMatrix<refDim, dim> > vJT_;
			if(vJT == NULL)
			{
				DimReferenceMapping<refDim, dim>& mapping
					= ReferenceMappingProvider::get<refDim, dim>(roid, vCornerCoords);
				vJT_.resize(nip);
				mapping.jacobian_transposed(&(vJT_[0]), vLocalIP, nip);
				vJT = &(vJT_[0]);
			}
			
			for(size_t ip = 0; ip < nip; ++ip)
			{
			//	get the local gradient (assuming the Lagrange-1 basis functions)
				trialSpace.grads(&(vLocalGrad[0]), vLocalIP[ip]);
				
			//	get the inverse Jacobian
				MathMatrix<dim,refDim> JTInv;
				RightInverse(JTInv, vJT[ip]);

			//	compute [rho*g]_consistent
				ConsGravityMethod.template compute<dim>
					(vValue[ip], vLocalIP[ip], JTInv, &(vLocalGrad[0]), &(vConsGravity[0]));
			}
		}
		
		template <int refDim>
		void prepare_dim_elem(GridObject* elem,
							const ReferenceObjectID roid,
							const MathVector<dim> vCornerCoords[])
		{
			const DimReferenceElement<refDim>& refElem = ReferenceElementProvider::get<refDim>(roid);
			
			m_densityCornerS = m_spDensity->template register_local_ip_series<refDim>
				(refElem.corners(), refElem.num(0), this->time_point(), false);
			m_spDensity->set_global_ips(m_densityCornerS, vCornerCoords, refElem.num(0));
		}

		virtual void prepare_element(GridObject* elem, const MathVector<dim> vCornerCoords[])
		{
			const ReferenceObjectID roid = elem->reference_object_id();
			const ReferenceElement& refElem = ReferenceElementProvider::get(roid);
			const int ref_dim = refElem.dimension();
			switch(ref_dim)
			{
				case 1: this->template prepare_dim_elem<1> (elem, roid, vCornerCoords); break;
				case 2: this->template prepare_dim_elem<2> (elem, roid, vCornerCoords); break;
				case 3: this->template prepare_dim_elem<3> (elem, roid, vCornerCoords); break;
				default: UG_THROW("ConsistentDarcyVelLinker: Ref. dimension " << ref_dim << " not supported.");
			}
		}

		template <int refDim>
		void eval_and_deriv(MathVector<dim> vValue[],
		                    const MathVector<dim> vGlobIP[],
		                    number time, int si,
		                    GridObject* elem,
		                    const MathVector<dim> vCornerCoords[],
		                    const MathVector<refDim> vLocalIP[],
		                    const size_t nip,
		                    LocalVector* u,
		                    bool bDeriv,
		                    int s,
		                    std::vector<std::vector<MathVector<dim> > > vvvDeriv[],
		                    const MathMatrix<refDim, dim>* vJT = NULL) const
		{
			const ReferenceObjectID roid = elem->reference_object_id();
			const DimReferenceElement<refDim>& refElem = ReferenceElementProvider::get<refDim>(roid);
			const size_t nco = refElem.num(ROID_VERTEX);
						
			const LocalShapeFunctionSet<refDim>& trialSpace =
				LocalFiniteElementProvider::get<refDim>(roid, LFEID(LFEID::LAGRANGE, refDim, 1));
	
		//	consistent gravity
			StdLinConsistentGravity<refDim> ConsGravityMethod;
			std::vector<MathVector<refDim> > vConsGravity(nco);
			
		//	inverse jacobians of the transformation at the integration points
			std::vector<MathMatrix<dim,refDim> > vJTInv(nip);
			
		//	local gradients at the integration points
			std::vector<MathVector<refDim> > vLocalGrad(nco);
			
		//	get the data of the ip series
			const number* vDensity = m_spDensity->values(m_densityCornerS);

		//	prepare the consistent gravity term
			try
			{
				ConsGravityMethod.template prepare<dim>
					(&(vConsGravity[0]), nco, vCornerCoords, vDensity, m_gravity);
			}
			UG_CATCH_THROW ("ConsistentGravityLinker: Cannot prepare the consistent gravity.");

			// get jacobians of the transformation mapping if not passed
			std::vector<MathMatrix<refDim, dim> > vJT_;
			if(vJT == NULL)
			{
				DimReferenceMapping<refDim, dim>& mapping
					= ReferenceMappingProvider::get<refDim, dim>(roid, vCornerCoords);
				vJT_.resize(nip);
				mapping.jacobian_transposed(&(vJT_[0]), vLocalIP, nip);
				vJT = &(vJT_[0]);
			}
			
			for(size_t ip = 0; ip < nip; ++ip)
			{
			//	get the local gradient (assuming the Lagrange-1 basis functions)
				trialSpace.grads(&(vLocalGrad[0]), vLocalIP[ip]);
				
			//	get the inverse Jacobian
				RightInverse(vJTInv[ip], vJT[ip]);

			//	compute [rho*g]_consistent
				ConsGravityMethod.template compute<dim>
					(vValue[ip], vLocalIP[ip], vJTInv[ip], &(vLocalGrad[0]), &(vConsGravity[0]));
			}

			if(!bDeriv)
				return;
		
        //	Compute the derivatives at all ips

			this->set_zero(vvvDeriv, nip);
			
			if(m_noDerivatives || this->zero_derivative() || (!m_spDDensity.valid()) || m_spDDensity->zero_derivative())
				return;

        //  prepare derivatives of the primary function at the corners
            std::vector<std::vector<MathVector<refDim> > > vvDConsGravity(nco);
            try
            {
                std::vector<number> DCoVal(nco);
                DCoVal.assign(nco, 0);
                for (size_t co = 0; co < nco; co++)
                {
                    DCoVal[co] = 1;
                    vvDConsGravity[co].resize(nco);
                    ConsGravityMethod.template prepare<dim>
                        (&(vvDConsGravity[co][0]), nco, vCornerCoords, &(DCoVal[0]), m_gravity);
                    DCoVal[co] = 0;
                }
            }
            UG_CATCH_THROW ("ConsistentGravityLinker: Cannot prepare Consistent Gravity for its derivatives.");
        
		//	compute the derivatives of the density
			if(m_noDerivatives || m_spDDensity->zero_derivative())
				return;
			
			for(size_t fct = 0; fct < m_spDDensity->num_fct(); ++fct) // fct = concentration, pressure, ...
			{
			//	get common fct id for this function
				const size_t commonFct = this->input_common_fct(_RHO_, fct);
				if(this->num_sh(commonFct) != nco)
					UG_THROW ("ConsistentGravityLinker: Number of shapes mismatch.");
				
				for(size_t ip = 0; ip < nip; ++ip) // we derive the cons. grav. at this ip
				{
					trialSpace.grads(&(vLocalGrad[0]), vLocalIP[ip]);
					
					for(size_t co = 0; co < nco; ++co) // w.r.t to the DoF at this corner (shape idx.)
					{
						const number* vDDensity = m_spDDensity->deriv(m_densityCornerS, co, fct);
						MathVector<dim>& deriv = vvvDeriv[ip][commonFct][co];
						
						ConsGravityMethod.template compute<dim>(deriv,
							vLocalIP[ip], vJTInv[ip], &(vLocalGrad[0]), &(vvDConsGravity[co][0]));
						VecScale(deriv, deriv, vDDensity[co]);
						
						/* Note that here we assume that the density has the local dependence
						 * on the arguments (concentration, pressure, ...): The density at
						 * corner co depends only on the DoFs at that corner. However, this
						 * excludes any differential operators, interpolations etc.
						 */
					}
				}
			}
		}

	public:
	
	///	set density import
		void set_density(SmartPtr<CplUserData<number, dim> > data)
		{
			m_spDensity = data;
			m_spDDensity = data.template cast_dynamic<DependentUserData<number, dim> >();
			base_type::set_input(_RHO_, data, data);
		}

		void set_density(number val)
		{
			set_density(make_sp(new ConstUserNumber<dim>(val)));
		}

	///	set gravity vector
		void set_gravity(MathVector<dim> g)
		{
			m_gravity = g;
		}
	
	/// set gravity from an array
		void set_gravity(const std::vector<number>& vGravity)
		{
			if(vGravity.size() != dim)
				UG_THROW("ConsistentGravityLinker: Illegal dimension of the specified gravity vector.");
			for(size_t i = 0; i < dim; i++) m_gravity[i] = vGravity[i];
		}
	
	///	set gravity in the z direction
		void set_gravity(number g)
		{
			m_gravity = 0;
			m_gravity[dim-1] = g;
		}
	
	protected:
	///	import density
		static const size_t _RHO_ = 0;
		SmartPtr<CplUserData<number, dim> > m_spDensity;
		SmartPtr<DependentUserData<number, dim> > m_spDDensity;
		size_t m_densityCornerS;

	///	constant gravity
		MathVector<dim> m_gravity;


	public:
		void set_no_derivatives(bool v) { m_noDerivatives = v; }

	protected:
		// disable the derivatives
		bool m_noDerivatives;
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__CONSISTENT_GRAVITY_LINKER__ */
