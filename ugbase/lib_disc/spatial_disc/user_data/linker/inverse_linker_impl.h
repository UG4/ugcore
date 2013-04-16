/*
 * inverse_linker_impl.h
 *
 *  Created on: 12.02.2013
 *      Author: ivomuha, andreasvogel scale_add_linker_impl.h used as template
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__INVERSE_LINKER_IMPL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__INVERSE_LINKER_IMPL__

#include "inverse_linker.h"
#include "linker_traits.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"

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
	m_vpDivisorData.resize(numInput+1, NULL);
	m_vpDependData.resize(numInput+1, NULL);
	m_vpDividendData.resize(numInput+1, NULL);
	m_vpDividendDependData.resize(numInput+1, NULL);

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

		value *= valDividend/valDivisor;
	}
}

template <int dim>
template <int refDim>
void InverseLinker<dim>::
evaluate(number vValue[],
         const MathVector<dim> vGlobIP[],
         number time, int si,
         GeometricObject* elem,
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
                    GeometricObject* elem,
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
		for(size_t c = 0; c < m_vpDivisorData.size(); ++c)
		{
			UG_ASSERT(divisor_value(c,s,ip)!=0, "DIVISOR IS 0");
			vValue[ip] *= dividend_value(c,s,ip)/divisor_value(c,s,ip);
		}
	}

//	check if derivative is required
	if(!bDeriv || this->zero_derivative()) return;

//	check sizes
	UG_ASSERT(m_vpDependData.size() == m_vpDividendDependData.size(),
	          	  	  	  	  	  	  	  	  	  "Wrong num Scales.");

//	clear all derivative values
	this->set_zero(vvvDeriv, nip);

//	loop all inputs
	for(size_t c = 0; c < m_vpDivisorData.size(); ++c)
	{
	//	check if Divisor has derivative
		if(!m_vpDivisorData[c]->zero_derivative())
		{
			for(size_t ip = 0; ip < nip; ++ip)
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
							UG_ASSERT(divisor_value(c,s,ip)!=0, "DIVISOR IS 0");
							vvvDeriv[ip][commonFct][sh] += vValue[ip] / (dividend_value(c,s,ip)/divisor_value(c,s,ip))*(-1.0)*(dividend_value(c,s,ip))/divisor_value(c,s,ip)/divisor_value(c,s,ip)*divisor_deriv(c, s, ip, fct, sh);
						}
						else
							vvvDeriv[ip][commonFct][sh] += 0;
						 //        input_deriv(c, s, ip, fct, sh),
						 //        scale_value(c, s, ip));
					}
				}
			}
		}

	//	check if Dividend has derivative
		if(!m_vpDividendData[c]->zero_derivative())
		{
			for(size_t ip = 0; ip < nip; ++ip)
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
							UG_ASSERT(divisor_value(c,s,ip)!=0, "DIVISOR IS 0");
							vvvDeriv[ip][commonFct][sh] += vValue[ip] / (dividend_value(c,s,ip)/divisor_value(c,s,ip))*dividend_deriv(c, s, ip, fct, sh)/divisor_value(c,s,ip);
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

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__INVERSE_LINKER_IMPL__ */
