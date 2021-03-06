/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__SCALE_ADD_LINKER_IMPL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__SCALE_ADD_LINKER_IMPL__

#include "scale_add_linker.h"
#include "linker_traits.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
//	ScaleAddLinker
////////////////////////////////////////////////////////////////////////////////

template <typename TData, int dim, typename TDataScale, typename TRet>
ScaleAddLinker<TData,dim,TDataScale,TRet>::
ScaleAddLinker(const ScaleAddLinker& linker)
{
	if(linker.m_vpUserData.size() != linker.m_vpScaleData.size())
		UG_THROW("ScaleAddLinker: number of scaling factors and data mismatch.");

	for(size_t i = 0; i < linker.m_vpUserData.size(); ++i)
	{
		this->add(linker.m_vpScaleData[i], linker.m_vpUserData[i]);
	}
}


template <typename TData, int dim, typename TDataScale, typename TRet>
void ScaleAddLinker<TData,dim,TDataScale,TRet>::
add(SmartPtr<CplUserData<TDataScale, dim> > scale, SmartPtr<CplUserData<TData, dim> > data)
{
//	current number of inputs
	const size_t numInput = base_type::num_input() / 2;

//	resize scaling
	m_vpUserData.resize(numInput+1);
	m_vpDependData.resize(numInput+1);
	m_vpScaleData.resize(numInput+1);
	m_vpScaleDependData.resize(numInput+1);

//	remember userdata
	UG_ASSERT(data.valid(), "Null Pointer as Input set.");
	m_vpUserData[numInput] = data;
	m_vpDependData[numInput] = data.template cast_dynamic<DependentUserData<TData, dim> >();

//	remember userdata
	UG_ASSERT(scale.valid(), "Null Pointer as Scale set.");
	m_vpScaleData[numInput] = scale;
	m_vpScaleDependData[numInput] = scale.template cast_dynamic<DependentUserData<TDataScale, dim> >();

//	increase number of inputs by one and set inputs at base class
	base_type::set_num_input(2*numInput+2);
	base_type::set_input(2*numInput, data, data);
	base_type::set_input(2*numInput+1, scale, scale);
}

template <typename TData, int dim, typename TDataScale, typename TRet>
void ScaleAddLinker<TData,dim,TDataScale,TRet>::
add(number scale, SmartPtr<CplUserData<TData, dim> > data)
{
	add(CreateConstUserData<dim>(scale, TDataScale()), data);
}

template <typename TData, int dim, typename TDataScale, typename TRet>
void ScaleAddLinker<TData,dim,TDataScale,TRet>::
add(SmartPtr<CplUserData<TDataScale, dim> > scale, number data)
{
	add(scale, CreateConstUserData<dim>(data, TData()));
}

template <typename TData, int dim, typename TDataScale, typename TRet>
void ScaleAddLinker<TData,dim,TDataScale,TRet>::
add(number scale, number data)
{
	add(CreateConstUserData<dim>(scale, TDataScale()),
	    CreateConstUserData<dim>(data, TData()));
}


template <typename TData, int dim, typename TDataScale, typename TRet>
void ScaleAddLinker<TData,dim,TDataScale,TRet>::
evaluate (TRet& value,
          const MathVector<dim>& globIP,
          number time, int si) const
{
	//	reset value
	value = 0.0;

	TData valData;
	TDataScale valScale;

//	add contribution of each summand
	for(size_t c = 0; c < m_vpUserData.size(); ++c)
	{
		(*m_vpUserData[c])(valData, globIP, time, si);
		(*m_vpScaleData[c])(valScale, globIP, time, si);

		linker_traits<TData, TDataScale,TRet>::
		mult_add(value, valData, valScale);
	}
}

template <typename TData, int dim, typename TDataScale, typename TRet>
template <int refDim>
void ScaleAddLinker<TData,dim,TDataScale,TRet>::
evaluate(TRet vValue[],
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
		vValue[ip] = 0.0;

	std::vector<TData> vValData(nip);
	std::vector<TDataScale> vValScale(nip);

//	add contribution of each summand
	for(size_t c = 0; c < m_vpUserData.size(); ++c)
	{
		(*m_vpUserData[c])(&vValData[0], vGlobIP, time, si,
						elem, vCornerCoords, vLocIP, nip, u, vJT);
		(*m_vpScaleData[c])(&vValScale[0], vGlobIP, time, si,
							elem, vCornerCoords, vLocIP, nip, u, vJT);

		for(size_t ip = 0; ip < nip; ++ip)
			linker_traits<TData, TDataScale,TRet>::
			mult_add(vValue[ip], vValData[ip], vValScale[ip]);
	}
}

template <typename TData, int dim, typename TDataScale, typename TRet>
template <int refDim>
void ScaleAddLinker<TData,dim,TDataScale,TRet>::
eval_and_deriv(TRet vValue[],
                    const MathVector<dim> vGlobIP[],
                    number time, int si,
                    GridObject* elem,
                    const MathVector<dim> vCornerCoords[],
                    const MathVector<refDim> vLocIP[],
                    const size_t nip,
                    LocalVector* u,
                    bool bDeriv,
                    int s,
                    std::vector<std::vector<TRet> > vvvDeriv[],
                    const MathMatrix<refDim, dim>* vJT) const
{
//	check that size of Scalings and inputs is equal
	UG_ASSERT(m_vpUserData.size() == m_vpScaleData.size(), "Wrong num Scales.");

//	compute value
	for(size_t ip = 0; ip < nip; ++ip)
	{
	//	reset value
		vValue[ip] = 0.0;

	//	add contribution of each summand
		for(size_t c = 0; c < m_vpUserData.size(); ++c)
		{
			linker_traits<TData, TDataScale,TRet>::
			mult_add(vValue[ip],
					 input_value(c, s, ip),
					 scale_value(c, s, ip));
		}
	}

//	check if derivative is required
	if(!bDeriv || this->zero_derivative()) return;

//	check sizes
	UG_ASSERT(m_vpDependData.size() == m_vpScaleDependData.size(),
	          	  	  	  	  	  	  	  	  	  "Wrong num Scales.");

//	clear all derivative values
	this->set_zero(vvvDeriv, nip);

//	loop all inputs
	for(size_t c = 0; c < m_vpUserData.size(); ++c)
	{
	//	check if input has derivative
		if(!m_vpUserData[c]->zero_derivative())
		{
			for(size_t ip = 0; ip < nip; ++ip)
			{
			//	loop functions
				for(size_t fct = 0; fct < input_num_fct(c); ++fct)
				{
				//	get common fct id for this function
					const size_t commonFct = input_common_fct(c, fct);

				//	loop dofs
					for(size_t sh = 0; sh < this->num_sh(fct); ++sh)
					{
						linker_traits<TData, TDataScale,TRet>::
						mult_add(vvvDeriv[ip][commonFct][sh],
								 input_deriv(c, s, ip, fct, sh),
								 scale_value(c, s, ip));
					}
				}
			}
		}

	//	check if scaling has derivative
		if(!m_vpScaleData[c]->zero_derivative())
		{
			for(size_t ip = 0; ip < nip; ++ip)
			{
			//	loop functions
				for(size_t fct = 0; fct < scale_num_fct(c); ++fct)
				{
				//	get common fct id for this function
					const size_t commonFct = scale_common_fct(c, fct);

				//	loop dofs
					for(size_t sh = 0; sh < this->num_sh(fct); ++sh)
					{
						linker_traits<TData, TDataScale,TRet>::
						mult_add(vvvDeriv[ip][commonFct][sh],
								 input_value(c, s, ip),
								 scale_deriv(c, s, ip, fct, sh));
					}
				}
			}
		}
	}
}


} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__SCALE_ADD_LINKER_IMPL__ */
