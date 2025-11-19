/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Ivo Muha, Andreas Vogel
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
 *      andreasvogel scale_add_linker_impl.h used as template
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__INVERSE_LINKER_IMPL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__INVERSE_LINKER_IMPL__

#include "inverse_linker.h"
#include "linker_traits.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#include "common/util/number_util.h"
namespace ug{

////////////////////////////////////////////////////////////////////////////////
//	InverseLinker
////////////////////////////////////////////////////////////////////////////////

template <int dim>
InverseLinker<dim>::
InverseLinker(const InverseLinker& linker)
{
	if(linker.m_vpDividendData.size() != linker.m_vpDivisorData.size())
		UG_THROW("InverseLinker: number of inputs mismatch.");

	for(size_t i = 0; i < linker.m_vpDivisorData.size(); ++i)
	{
		this->divide(linker.m_vpDividendData[i], linker.m_vpDivisorData[i]);
	}
}


template <int dim>
void InverseLinker<dim>::
divide(SmartPtr<CplUserData<number, dim> > dividend, SmartPtr<CplUserData<number, dim> > divisor)
{

//	current number of inputs
	const size_t numInput = base_type::num_input() / 2;

//	resize scaling
	m_vpDivisorData.resize(numInput+1);
	m_vpDependData.resize(numInput+1);
	m_vpDividendData.resize(numInput+1);
	m_vpDividendDependData.resize(numInput+1);

//	remember userdata
	UG_ASSERT(divisor.valid(), "Null Pointer as Input set.");
	m_vpDivisorData[numInput] = divisor;
	m_vpDependData[numInput] = divisor.template cast_dynamic<DependentUserData<number, dim> >();

//	remember userdata
	UG_ASSERT(dividend.valid(), "Null Pointer as Scale set.");
	m_vpDividendData[numInput] = dividend;
	m_vpDividendDependData[numInput] = dividend.template cast_dynamic<DependentUserData<number, dim> >();

//	increase number of inputs by one and set inputs at base class
	base_type::set_num_input(2*numInput+2);
	base_type::set_input(2*numInput, divisor, divisor);
	base_type::set_input(2*numInput+1, dividend, dividend);
}

template <int dim>
void InverseLinker<dim>::
divide(number dividend, SmartPtr<CplUserData<number, dim> > divisor)
{
	divide(CreateConstUserData<dim>(dividend, number()), divisor);
}

template <int dim>
void InverseLinker<dim>::
divide(SmartPtr<CplUserData<number, dim> > dividend, number divisor)
{
	divide(dividend, CreateConstUserData<dim>(divisor, number()));
}

template <int dim>
void InverseLinker<dim>::
divide(number dividend, number divisor)
{
	divide(CreateConstUserData<dim>(dividend, number()),
	    CreateConstUserData<dim>(divisor, number()));
}

//Scale ist Dividend! UserData ist Divisor
template <int dim>
void InverseLinker<dim>::
evaluate (number& value,
          const MathVector<dim>& globIP,
          number time, int si) const
{
	//	reset value
	value = 1.0;

	number valDivisor=0;
	number valDividend=0;

	for(size_t c = 0; c < m_vpDivisorData.size(); ++c)
	{
		(*m_vpDivisorData[c])(valDivisor, globIP, time, si);
		(*m_vpDividendData[c])(valDividend, globIP, time, si);
		UG_ASSERT(valDivisor!=0, "DIVISOR IS 0");
		UG_COND_THROW(valDivisor==0, "DIVISOR IS 0");

		value *= valDividend/valDivisor;
	}
}

template <int dim>
template <int refDim>
void InverseLinker<dim>::
evaluate(number vValue[],
         const MathVector<dim> vGlobIP[],
         number time, int si,
         GridObject* elem,
         const MathVector<dim> vCornerCoords[],
         const MathVector<refDim> vLocIP[],
         const size_t nip,
         LocalVector* u,
         const MathMatrix<refDim, dim>* vJT) const
{
	//	reset value
	for(size_t ip = 0; ip < nip; ++ip)
		vValue[ip] = 1.0;

	std::vector<number> vValData(nip);
	std::vector<number> vValScale(nip);

//	add contribution of each summand
	for(size_t c = 0; c < m_vpDivisorData.size(); ++c)
	{
		(*m_vpDivisorData[c])(&vValData[0], vGlobIP, time, si,
						elem, vCornerCoords, vLocIP, nip, u, vJT);
		(*m_vpDividendData[c])(&vValScale[0], vGlobIP, time, si,
							elem, vCornerCoords, vLocIP, nip, u, vJT);

		for(size_t ip = 0; ip < nip; ++ip)
			vValue[ip] *=  vValScale[ip]/vValData[ip];
	}
}

template <int dim>
template <int refDim>
void InverseLinker<dim>::
eval_and_deriv(number vValue[],
                    const MathVector<dim> vGlobIP[],
                    number time, int si,
                    GridObject* elem,
                    const MathVector<dim> vCornerCoords[],
                    const MathVector<refDim> vLocIP[],
                    const size_t nip,
                    LocalVector* u,
                    bool bDeriv,
                    int s,
                    std::vector<std::vector<number> > vvvDeriv[],
                    const MathMatrix<refDim, dim>* vJT) const
{
//	check that size of Scalings and inputs is equal
	UG_ASSERT(m_vpDivisorData.size() == m_vpDividendData.size(), "Wrong num Scales.");

//	compute value
	for(size_t ip = 0; ip < nip; ++ip)
	{
	//	reset value
		vValue[ip] = 1.0;

	//	add contribution of each summand
		for(size_t c2 = 0; c2 < m_vpDivisorData.size(); ++c2)
		{
			UG_COND_THROW(CloseToZero(divisor_value(c2,s,ip)), "DIVISOR IS 0");
			vValue[ip] *= dividend_value(c2,s,ip)/divisor_value(c2,s,ip);
		}
	}

//	compute value
//	check if derivative is required
	if(!bDeriv || this->zero_derivative()) return;

//	check sizes
	UG_ASSERT(m_vpDependData.size() == m_vpDividendDependData.size(),
	          	  	  	  	  	  	  	  	  	  "Wrong num Scales.");

//	clear all derivative values
	this->set_zero(vvvDeriv, nip);

// (dividend/divisor)' = (u/v)' = (u'v - uv')/v^2  = u'/v - uv'/v^2
//
// now u = prod_i u_i and v = prod_i v_i
// set nu_i = prod_{j != i} u_j  and  nv_i = prod_{j != i} v_j
// note that prod_i u_i = u_j nu_j  for all j.
// then u'_i = sum_j  nu_j u_j'
//
// (u/v)' = ((sum_j  nu_j u_j') prod_i v_i - (sum_j  nv_j v_j') prod_i u_i ) / prod_i v_i^2
// 		  = ((sum_j  nu_j u_j' v_j nv_j - sum_j  nv_j v_j' u_j nu_j ) / v_j^2 nv_j^2
//        =  sum_j (  nu_j / nv_j  * (u_j' v_j - u_j v_j') / v_j^2 )


//	loop all inputs
	for(size_t c = 0; c < m_vpDivisorData.size(); ++c)
	{
		for(size_t ip = 0; ip < nip; ++ip)
		{
			double vValue = 1.0; // vValue = nu_j / nv_j = prod_{k ! = j} u_k/v_k
			for(size_t c2 = 0; c2 < m_vpDivisorData.size(); ++c2)
			{
				if(c == c2) continue;
				vValue *= dividend_value(c2,s,ip)/divisor_value(c2,s,ip);
			}

		//	check if Divisor has derivative
			if(!m_vpDivisorData[c]->zero_derivative())
			{

			//	loop functions
				for(size_t fct = 0; fct < divisor_num_fct(c); ++fct)
				{
				//	get common fct id for this function
					const size_t commonFct = divisor_common_fct(c, fct);

				//	loop dofs
					for(size_t sh = 0; sh < this->num_sh(fct); ++sh)
					{
						if(dividend_value(c,s,ip)!=0)
						{

							vvvDeriv[ip][commonFct][sh] +=
									vValue *
									(-1.0)*
									(dividend_value(c,s,ip))/divisor_value(c,s,ip)/divisor_value(c,s,ip)
									*divisor_deriv(c, s, ip, fct, sh);
							UG_COND_THROW(IsFiniteAndNotTooBig(vvvDeriv[ip][commonFct][sh])==false, "");
						}
						else
							vvvDeriv[ip][commonFct][sh] += 0;
						 //        input_deriv(c, s, ip, fct, sh),
						 //        scale_value(c, s, ip));
					}
				}
			}

		//	check if Dividend has derivative
			if(!m_vpDividendData[c]->zero_derivative())
			{

			//	loop functions
				for(size_t fct = 0; fct < dividend_num_fct(c); ++fct)
				{
				//	get common fct id for this function
					const size_t commonFct = dividend_common_fct(c, fct);

				//	loop dofs
					for(size_t sh = 0; sh < this->num_sh(fct); ++sh)
					{
						if(dividend_value(c,s,ip)!=0)
						{
							vvvDeriv[ip][commonFct][sh] +=
									vValue
									* dividend_deriv(c, s, ip, fct, sh) / divisor_value(c,s,ip);
							UG_COND_THROW(IsFiniteAndNotTooBig(vvvDeriv[ip][commonFct][sh])==false, "");
						}
						else
							vvvDeriv[ip][commonFct][sh] += 0;

							//linker_traits<TData, TDataScale>::
						//mult_add(this->deriv(s, ip, commonFct, sh),
						//		 input_value(c, s, ip),
						//		 scale_deriv(c, s, ip, fct, sh));
					}
				}
			}

		}


	}

}


} // end namespace ug

#endif