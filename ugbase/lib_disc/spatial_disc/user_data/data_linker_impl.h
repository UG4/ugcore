/*
 * data_linker_impl.h
 *
 *  Created on: 04.07.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DATA_LINKER_IMPL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DATA_LINKER_IMPL__

#include "data_linker.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
//	Data Linker
////////////////////////////////////////////////////////////////////////////////

template <typename TData, int dim>
bool DataLinker<TData,dim>::zero_derivative() const
{
	bool bRet = true;
	for(size_t i = 0; i < m_vspICplUserData.size(); ++i)
		bRet &= m_vspICplUserData[i]->zero_derivative();
	return bRet;
}

template <typename TData, int dim>
void DataLinker<TData,dim>::check_setup() const
{
//	check, that all inputs are set
	for(size_t i = 0; i < num_input(); ++i)
		if(!m_vspICplUserData[i].valid())
			UG_THROW("DataLinker::check_setup: Input number "<<i<<" missing.");
}

template <typename TData, int dim>
void DataLinker<TData,dim>::update_function_group()
{
//	collect all function groups
	std::vector<const FunctionGroup*> vFctGrp(num_input(), NULL);
	for(size_t i = 0; i < m_vspICplUserData.size(); ++i)
		if(m_vspICplUserData[i].valid())
			if(m_vspICplUserData[i]->num_fct() != 0)
				vFctGrp[i] = &(m_vspICplUserData[i]->function_group());

//	create union of all function groups
	try{
		CreateUnionOfFunctionGroups(m_commonFctGroup, vFctGrp, true);
	}UG_CATCH_THROW("'DataLinker::update_function_group': Cannot create"
					" common function group.");

//	create FunctionIndexMapping for each Disc
	m_vMap.resize(vFctGrp.size());
	for(size_t i = 0; i < vFctGrp.size(); ++i)
	{
		if(vFctGrp[i] != NULL)
		{
			try{
				CreateFunctionIndexMapping(m_vMap[i], *vFctGrp[i], m_commonFctGroup);
			}UG_CATCH_THROW("'DataLinker::update_function_group':"
							"Cannot create Function Index Mapping for input "<<i<<".");
		}
	}

//	set common function group as the function group the data depends on
	this->set_function_group(m_commonFctGroup);
}

template <typename TData, int dim>
void DataLinker<TData,dim>::
local_ip_series_added(const size_t seriesID)
{
	const size_t s = seriesID;

//	 we need a series id for all inputs
	m_vvSeriesID.resize(m_vspICplUserData.size());

//	loop inputs
	for(size_t i = 0; i < m_vspICplUserData.size(); ++i)
	{
	//	check unset data
		UG_ASSERT(m_vspICplUserData[i].valid(), "No Input set, but requested.");

	//	resize series ids
		m_vvSeriesID[i].resize(s+1);

	//	request local ips for series at input data
		switch(this->dim_local_ips())
		{
			case 1:
				m_vvSeriesID[i][s] =
						m_vspICplUserData[i]->template register_local_ip_series<1>
								(this->template local_ips<1>(s), this->num_ip(s),
								 this->m_vMayChange[s]);
				break;
			case 2:
				m_vvSeriesID[i][s] =
						m_vspICplUserData[i]->template register_local_ip_series<2>
								(this->template local_ips<2>(s), this->num_ip(s),
								 this->m_vMayChange[s]);
				break;
			case 3:
				m_vvSeriesID[i][s] =
						m_vspICplUserData[i]->template register_local_ip_series<3>
								(this->template local_ips<3>(s), this->num_ip(s),
								 this->m_vMayChange[s]);
				break;
			default: UG_THROW("Dimension not supported."); break;
		}
	}

//	resize data fields
	DependentUserData<TData, dim>::local_ip_series_added(seriesID);
}


template <typename TData, int dim>
void DataLinker<TData,dim>::
local_ips_changed(const size_t seriesID, const size_t newNumIP)
{
	const size_t s = seriesID;

//	loop inputs
	for(size_t i = 0; i < m_vspICplUserData.size(); ++i)
	{
	//	skip unset data
		UG_ASSERT(m_vspICplUserData[i].valid(), "No Input set, but requested.");

		switch(this->dim_local_ips())
		{
			case 1: m_vspICplUserData[i]->template set_local_ips<1>
					(m_vvSeriesID[i][s], this->template local_ips<1>(s), this->num_ip(s));
				break;
			case 2: m_vspICplUserData[i]->template set_local_ips<2>
					(m_vvSeriesID[i][s], this->template local_ips<2>(s), this->num_ip(s));
				break;
			case 3: m_vspICplUserData[i]->template set_local_ips<3>
					(m_vvSeriesID[i][s], this->template local_ips<3>(s), this->num_ip(s));
				break;
			default: UG_THROW("Dimension not supported."); break;
		}
	}

//	resize data fields
	DependentUserData<TData, dim>::local_ips_changed(seriesID, newNumIP);
}

template <typename TData, int dim>
void DataLinker<TData,dim>::
global_ips_changed(const size_t seriesID, const MathVector<dim>* vPos, const size_t numIP)
{
//	loop inputs
	for(size_t i = 0; i < m_vspICplUserData.size(); ++i)
	{
	//	skip unset data
		UG_ASSERT(m_vspICplUserData[i].valid(), "No Input set, but requested.");

	//	adjust global ids of imported data
		m_vspICplUserData[i]->set_global_ips(m_vvSeriesID[i][seriesID], vPos, numIP);
	}
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DATA_LINKER_IMPL__ */
