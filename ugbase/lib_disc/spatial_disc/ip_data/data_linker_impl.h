/*
 * data_linker_impl.h
 *
 *  Created on: 04.07.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISCRETIZATION__DATA_LINKER_IMPL__
#define __H__UG__LIB_DISC__SPATIAL_DISCRETIZATION__DATA_LINKER_IMPL__

#include "data_linker.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
//	Data Linker
////////////////////////////////////////////////////////////////////////////////

template <typename TData, int dim>
bool DataLinker<TData,dim>::zero_derivative() const
{
	bool bRet = true;

//	loop inputs
	for(size_t i = 0; i < m_vpIIPData.size(); ++i)
	{
	//	skip unset data ( null as default )
		if(m_vpIIPData[i] == NULL) continue;

	//	flag iff ipdata is dependent
		bRet = bRet && m_vpIIPData[i]->zero_derivative();
	}

	return bRet;
}

template <typename TData, int dim>
bool DataLinker<TData,dim>::is_ready() const
{
//	check, that all inputs are set
	for(size_t i = 0; i < num_input(); ++i)
		if(m_vpIIPData[i] == NULL)
		{
			UG_LOG("ERROR in 'DataLinker::is_ready': Input number "<<
					i << " missing.\n");
			return false;
		}

//	everything ok
	return true;
}

template <typename TData, int dim>
bool DataLinker<TData,dim>::update_function_group()
{
//	collect all function groups
	std::vector<const FunctionGroup*> vFctGrp(num_input(), NULL);
	for(size_t i = 0; i < m_vpIDependData.size(); ++i)
		if(m_vpIDependData[i] != NULL)
			vFctGrp[i] = &(m_vpIDependData[i]->get_function_group());

//	create union of all function groups
	if(!CreateUnionOfFunctionGroups(m_commonFctGroup, vFctGrp, true))
	{
		UG_LOG("ERROR in 'DataLinker::update_function_group': Cannot create"
				" common function group.\n");
		return false;
	}

//	create FunctionIndexMapping for each Disc
	m_vMap.resize(vFctGrp.size());
	for(size_t i = 0; i < vFctGrp.size(); ++i)
	{
		if(vFctGrp[i] != NULL)
			if(!CreateFunctionIndexMapping(m_vMap[i], *vFctGrp[i],
										   m_commonFctGroup))
			{
				UG_LOG("ERROR in 'DataLinker::update_function_group':"
						"Cannot create Function Index Mapping for input "<<i<<".\n");
				return false;
			}
	}

//	set common function group as the function group the data depends on
	this->set_function_group(m_commonFctGroup);

//	we're done
	return true;
}

template <typename TData, int dim>
void DataLinker<TData,dim>::
local_ips_added()
{
//	 we need a series id for all inputs
	m_vvSeriesID.resize(m_vpIIPData.size());

//	loop inputs
	for(size_t i = 0; i < m_vpIIPData.size(); ++i)
	{
	//	resize series ids
		m_vvSeriesID[i].resize(this->num_series());

	//	skip unset data
		UG_ASSERT(m_vpIIPData[i] != NULL, "No Input set, but requested.");

	//	request local ips for all series at input data
		for(size_t s = 0; s < m_vvSeriesID[i].size(); ++s)
		{
			switch(this->dim_local_ips())
			{
				case 1:
					m_vvSeriesID[i][s] =
							m_vpIIPData[i]->template register_local_ip_series<1>
									(this->template local_ips<1>(s), num_ip(s));
					break;
				case 2:
					m_vvSeriesID[i][s] =
							m_vpIIPData[i]->template register_local_ip_series<2>
									(this->template local_ips<2>(s), num_ip(s));
					break;
				case 3:
					m_vvSeriesID[i][s] =
							m_vpIIPData[i]->template register_local_ip_series<3>
									(this->template local_ips<3>(s), num_ip(s));
					break;
				default: throw(UGFatalError("Dimension not supported."));
			}
		}
	}

//	resize data fields
	DependentIPData<TData, dim>::local_ips_added();
}

template <typename TData, int dim>
void DataLinker<TData,dim>::
global_ips_changed(size_t s, const MathVector<dim>* vPos, size_t numIP)
{
//	loop inputs
	for(size_t i = 0; i < m_vpIIPData.size(); ++i)
	{
	//	skip unset data
		UG_ASSERT(m_vpIIPData[i] != NULL, "No Input set, but requested.");

	//	adjust global ids of imported data
		m_vpIIPData[i]->set_global_ips(m_vvSeriesID[i][s], vPos, numIP);
	}
}

////////////////////////////////////////////////////////////////////////////////
//	DataLinkerEqualData
////////////////////////////////////////////////////////////////////////////////

template <typename TData, int dim, typename TDataIn>
void DataLinkerEqualData<TData,dim,TDataIn>::set_num_input(size_t num)
{
//	resize arrays
	m_vpIPData.resize(num, NULL);
	m_vpDependData.resize(num, NULL);

//	forward size to base class
	base_type::set_num_input(num);
}

template <typename TData, int dim, typename TDataIn>
bool DataLinkerEqualData<TData,dim,TDataIn>::
set_input(size_t i, IPData<TDataIn, dim>& data)
{
	UG_ASSERT(i < m_vpIPData.size(), "Input not needed");
	UG_ASSERT(i < m_vpDependData.size(), "Input not needed");

//	check input number
	if(i >= num_input())
	{
		UG_LOG("ERROR in 'DataLinker::set_input': Only " << num_input()
		       << " inputs can be set. Use 'set_num_input' to increase"
		       " the number of needed inputs.\n");
		return false;
	}

//	remember ipdata
	m_vpIPData[i] = &data;

//	cast to dependent data
	m_vpDependData[i] = dynamic_cast<DependentIPData<TDataIn, dim>*>(&data);

//	forward to base class
	base_type::set_input(i, &data);

//	we're done
	return true;
}

////////////////////////////////////////////////////////////////////////////////
//	ScaleAddLinker
////////////////////////////////////////////////////////////////////////////////

template <typename TData, int dim, typename TDataScale>
bool ScaleAddLinker<TData,dim,TDataScale>::
add(IPData<TDataScale, dim>& scale, IPData<TData, dim>& data)
{
//	current number of inputs
	const size_t numInput = base_type::num_input() / 2;

//	resize scaling
	resize_scaling(numInput+1);

//	remember ipdata
	m_vpIPData[numInput] = &data;
	UG_ASSERT(m_vpIPData[numInput] != NULL, "Null Pointer as Input set.");
	m_vpDependData[numInput] = dynamic_cast<DependentIPData<TData, dim>*>(&data);

//	remember ipdata
	m_vpScaleData[numInput] = &scale;
	UG_ASSERT(m_vpScaleData[numInput] != NULL, "Null Pointer as Scale set.");
	m_vpScaleDependData[numInput]
	              = dynamic_cast<DependentIPData<TDataScale, dim>*>(&scale);

//	increase number of inputs by one
	base_type::set_num_input(2*numInput+2);

//	add this input
	base_type::set_input(2*numInput, &data);
	base_type::set_input(2*numInput+1, &scale);

//	done
	return true;
}

template <typename TData, int dim, typename TDataScale>
bool ScaleAddLinker<TData,dim,TDataScale>::compute(bool bDeriv)
{
//	check that size of Scalings and inputs is equal
	UG_ASSERT(m_vpIPData.size() == m_vpScaleData.size(), "Wrong num Scales.");

//	compute value
	for(size_t s = 0; s < num_series(); ++s)
		for(size_t ip = 0; ip < num_ip(s); ++ip)
		{
		//	reset value
			value(s,ip) = 0.0;

		//	add contribution of each summand
			for(size_t c = 0; c < m_vpIPData.size(); ++c)
			{
				linker_traits<TData, TDataScale>::
				mult_add(value(s, ip),
				         input_value(c, s, ip),
				         scale_value(c, s, ip));
			}
		}

//	check if derivative is required
	if(!bDeriv || this->zero_derivative()) return true;

//	check sizes
	UG_ASSERT(m_vpDependData.size() == m_vpScaleDependData.size(),
	          	  	  	  	  	  	  	  	  	  "Wrong num Scales.");

//	clear all derivative values
	this->clear_derivative_values();

//	loop all inputs
	for(size_t c = 0; c < m_vpIPData.size(); ++c)
	{
	//	check if input has derivative
		if(!m_vpIPData[c]->zero_derivative())
		{
			for(size_t s = 0; s < num_series(); ++s)
				for(size_t ip = 0; ip < num_ip(s); ++ip)
				{
				//	loop functions
					for(size_t fct = 0; fct < input_num_fct(c); ++fct)
					{
					//	get common fct id for this function
						const size_t commonFct = input_common_fct(c, fct);

					//	loop dofs
						for(size_t sh = 0; sh < num_sh(s, fct); ++sh)
						{
							linker_traits<TData, TDataScale>::
							mult_add(deriv(s, ip, commonFct, sh),
							         input_deriv(c, s, ip, fct, sh),
							         scale_value(c, s, ip));
						}
					}
				}
		}

	//	check if scaling has derivative
		if(!m_vpScaleData[c]->zero_derivative())
		{
			for(size_t s = 0; s < num_series(); ++s)
				for(size_t ip = 0; ip < num_ip(s); ++ip)
				{
				//	loop functions
					for(size_t fct = 0; fct < scale_num_fct(c); ++fct)
					{
					//	get common fct id for this function
						const size_t commonFct = scale_common_fct(c, fct);

					//	loop dofs
						for(size_t sh = 0; sh < num_sh(s, fct); ++sh)
						{
							linker_traits<TData, TDataScale>::
							mult_add(deriv(s, ip, commonFct, sh),
									 input_value(c, s, ip),
									 scale_deriv(c, s, ip, fct, sh));
						}
					}
				}
		}
	}

	return true;
}


} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISCRETIZATION__DATA_LINKER_IMPL__ */
