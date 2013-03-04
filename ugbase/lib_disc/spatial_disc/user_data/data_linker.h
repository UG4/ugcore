/*
 * data_linker.h
 *
 *  Created on: 12.11.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DATA_LINKER__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DATA_LINKER__

#include "user_data.h"
#include "lib_disc/common/groups_util.h"

namespace ug{

/// combines several UserDatas to a new UserData of a specified type
/**
 * This class provides data at integration points and implements the
 * DependentUserData interface.
 *
 * \tparam 	TData		output Data type
 * \tparam 	dim			World dimension
 */
template <typename TData, int dim>
class DataLinker
	: public DependentUserData<TData, dim>
{
	public:
	///	constructor
		DataLinker() {m_vpIUserData.clear(); m_vpIDependData.clear();}

	///	returns if derivative is zero
		virtual bool zero_derivative() const;

	///	returns if the derivative of the i'th input is zero
		bool zero_derivative(size_t i) const
		{
			if(!m_vpIUserData[i].valid()) return true;
			return m_vpIUserData[i]->zero_derivative();
		}

	///	sets the number of inputs
		void set_num_input(size_t num)
		{
			m_vpIUserData.resize(num, NULL);
			m_vpIDependData.resize(num, NULL);
		}

	///	sets an input
		virtual void set_input(size_t i, SmartPtr<IDimUserData<dim> > input)
		{
			UG_ASSERT(i < m_vpIUserData.size(), "invalid index");
			m_vpIUserData[i] = input;
			m_vpIDependData[i] = input;
		}

	///	number of inputs
		virtual size_t num_input() const {return num_needed_data();}

	///	number of other Data this data depends on
		virtual size_t num_needed_data() const {return m_vpIUserData.size();}

	///	return needed data
		virtual SmartPtr<IUserData> needed_data(size_t i)
		{
			UG_ASSERT(i < m_vpIUserData.size(), "Input not needed");
			UG_ASSERT(m_vpIUserData[i].valid(), "Data input not valid");
			return m_vpIUserData[i];
		}

	///	returns if data is ok
		virtual void check_setup() const;

	///	updates the function group
		void update_function_group();

	protected:
	///	returns number of functions the input depends on
		size_t input_num_fct(size_t i) const
		{
			UG_ASSERT(i < m_vpIDependData.size(), "Input invalid");
			if(!m_vpIDependData[i].valid()) return 0;
			return m_vpIDependData[i]->num_fct();
		}

	///	returns the number in the common FctGrp for a fct of an input
		size_t input_common_fct(size_t i, size_t fct) const
		{
			UG_ASSERT(i < m_vMap.size(), "Input Map invalid");
			UG_ASSERT(fct < m_vMap[i].num_fct(), "Input Map invalid for fct");
			return m_vMap[i][fct];
		}

	///	returns the series id set for the i'th input
		size_t series_id(size_t i, size_t s) const
		{
			UG_ASSERT(i < m_vvSeriesID.size(), "invalid index");
			UG_ASSERT(s < m_vvSeriesID[i].size(), "invalid index");
			return m_vvSeriesID[i][s];
		}

	///	requests series id's from input data
		virtual void local_ip_series_added(const size_t seriesID);

	///	forwards the local positions to the data inputs
		virtual void local_ips_changed(const size_t seriesID, const size_t newNumIP);

	///	forwards the global positions to the data inputs
		virtual void global_ips_changed(const size_t seriesID, const MathVector<dim>* vPos, const size_t numIP);

	protected:
	///	data input
		std::vector<SmartPtr<IDimUserData<dim> > > m_vpIUserData;

	///	data input casted to IDependend data
		std::vector<SmartPtr<IUserData> > m_vpIDependData;

	///	common functions the data depends on
		FunctionGroup m_commonFctGroup;

	///	Function mapping for each input relative to common FunctionGroup
		std::vector<FunctionIndexMapping> m_vMap;

	///	series id the linker uses to get data from input
		std::vector<std::vector<size_t> > m_vvSeriesID;
};

} // end namespace ug

//	include implementation
#include "data_linker_impl.h"

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DATA_LINKER__ */
