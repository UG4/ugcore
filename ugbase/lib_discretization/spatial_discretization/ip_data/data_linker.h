/*
 * data_linker.h
 *
 *  Created on: 12.11.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_LINKER__
#define __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_LINKER__

#include "ip_data.h"
#include "lib_discretization/common/groups_util.h"

namespace ug{

/// combines several IPDatas of the one data type to a new IPData of a second type
/**
 * This class provides data at integration points and implements the
 * DependentIPData interface. It gets passed several data of one C++-type
 * and combines them to a new value of a second type. In addition also the
 * derivative w.r.t. the primal unknowns are passed.
 *
 * \tparam 	TData		output Data type
 * \tparam 	dim			World dimension
 * \tparam	TDataIn		Input data type
 */
template <typename TData, int dim, typename TDataIn>
class DataLinker
	: public DependentIPData<TData, dim>
{
	public:
	///	Base class type
		typedef DependentIPData<TData, dim> base_type;

	//	explicitly forward methods of IIPData
		using base_type::num_series;
		using base_type::num_ip;
		using base_type::time;

	//	explicitly forward methods of IPData
		using base_type::ip;
		using base_type::value;

	public:
	///	constructor
		DataLinker()
		{
			m_vpIPData.clear();
			m_vpDependData.clear();
		}

	///	compute method
		virtual void compute(bool compDeriv) = 0;

	///	returns if derivative is zero
		virtual bool zero_derivative() const
		{
			bool bRet = true;

		//	loop inputs
			for(size_t i = 0; i < m_vpIPData.size(); ++i)
			{
			//	skip unset data ( null as default )
				if(m_vpIPData[i] == NULL) continue;

			//	flag iff ipdata is dependent
				bRet = bRet && m_vpIPData[i]->zero_derivative();
			}

			return bRet;
 		}

	///	returns if the derivative of the i'th input is zero
		bool zero_derivative(size_t i) const
		{
			if(m_vpIPData[i] == NULL) return true;
			return m_vpIPData[i]->zero_derivative();
		}

	///	set number of needed inputs
		void set_num_input(size_t num)
		{
			m_vpIPData.clear(); m_vpIPData.resize(num, NULL);
			m_vpDependData.clear(); m_vpDependData.resize(num, NULL);
		}

	///	set input i
		bool set_input(size_t i, IPData<TDataIn, dim>& data)
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

		//	we're done
			return true;
		}

	///	returns if data is ok
		virtual bool make_ready()
		{
		//	check, that all inputs are set
			for(size_t i = 0; i < num_input(); ++i)
			{
				if(m_vpIPData[i] == NULL)
				{
					UG_LOG("ERROR in 'DataLinker::make_ready': Input number "<<
							i << " missing.\n");
					return false;
				}
			}

		//	if data is dependent, forward function group
			if(!update_function_group())
			{
				UG_LOG("ERROR in 'DataLinker::set_input': Cannot build function"
						" group of linker.\n");
				return false;
			}

			return true;
		}

	///	number of inputs
		size_t num_input() const {return num_needed_data();}

	///	number of other Data this data depends on
		virtual size_t num_needed_data() const {return m_vpIPData.size();}

	///	return needed data
		virtual IIPData* needed_data(size_t i)
		{
			UG_ASSERT(i < m_vpIPData.size(), "Input not needed");
			UG_ASSERT(m_vpIPData[i] != NULL, "Data input not valid");
			return m_vpIPData[i];
		}

	protected:
	///	updates the function group
		bool update_function_group()
		{
		//	collect all function groups
			std::vector<const FunctionGroup*> vFctGrp(num_input(), NULL);
			for(size_t i = 0; i < m_vpDependData.size(); ++i)
				if(m_vpDependData[i] != NULL)
					vFctGrp[i] = &(m_vpDependData[i]->get_function_group());

		//	create union of all function groups
			if(!CreateUnionOfFunctionGroups(m_commonFctGroup, vFctGrp, true))
			{
				UG_LOG("ERROR in 'DataLinker::update_function_group': Cannot create"
						"common function group.\n");
				return false;
			}

		//	create FunctionIndexMapping for each Disc
			m_vMap.resize(m_vpDependData.size());
			for(size_t i = 0; i < m_vpDependData.size(); ++i)
			{
				if(m_vpDependData[i] != NULL)
				if(!CreateFunctionIndexMapping(m_vMap[i],
				                               (m_vpDependData[i]->get_function_group()),
				                               m_commonFctGroup))
				{
					UG_LOG("ERROR in 'DataLinker::update_function_group':"
							"Cannot create Function Index Mapping.\n");
					return false;
				}
			}

		//	set common function group as the function group the data depends on
			this->set_function_group(m_commonFctGroup);

		//	we're done
			return true;
		}

	///	returns number of functions the input depends on
		size_t input_num_fct(size_t i) const
		{
			UG_ASSERT(i < m_vpDependData.size(), "Input invalid");
			if(m_vpDependData[i] == NULL) return 0;
			return m_vpDependData[i]->num_fct();
		}

	///	returns the number in the common FctGrp for a fct of an input
		size_t input_common_fct(size_t i, size_t fct) const
		{
			UG_ASSERT(i < m_vMap.size(), "Input Map invalid");
			UG_ASSERT(fct < m_vMap[i].num_fct(), "Input Map invalid for fct");
			return m_vMap[i][fct];
		}

	///	data at ip of input
		const TData& input_value(size_t i, size_t s, size_t ip) const
		{
			UG_ASSERT(i < m_vpIPData.size(), "Input not needed");
			UG_ASSERT(m_vpIPData[i] != NULL, "Input invalid");
			UG_ASSERT(i < m_vvSeriesID.size(), "Input Series missing");
			UG_ASSERT(s < m_vvSeriesID[i].size(), "Input Series invalid");
			return m_vpIPData[i]->value(m_vvSeriesID[i][s], ip);
		}

	///	data at ip of input
		TData& input_value(size_t i, size_t s, size_t ip)
		{
			UG_ASSERT(i < m_vpIPData.size(), "Input not needed");
			UG_ASSERT(m_vpIPData[i] != NULL, "Input invalid");
			UG_ASSERT(i < m_vvSeriesID.size(), "Input Series missing");
			UG_ASSERT(s < m_vvSeriesID[i].size(), "Input Series invalid");
			return m_vpIPData[i]->value(m_vvSeriesID[i][s], ip);
		}

	///	derivative of data at input at ip
		const TData& input_deriv(size_t i, size_t s, size_t ip, size_t fct, size_t dof) const
		{
			UG_ASSERT(i < m_vpDependData.size(), "Input not needed");
			UG_ASSERT(m_vpDependData[i] != NULL, "Input invalid");
			UG_ASSERT(i < m_vvSeriesID.size(), "Input Series missing");
			UG_ASSERT(s < m_vvSeriesID[i].size(), "Input Series invalid");
			return m_vpDependData[i]->deriv(m_vvSeriesID[i][s], ip, fct, dof);
		}

	///	derivative of data at input at ip
		TData& input_deriv(size_t i, size_t s, size_t ip, size_t fct, size_t dof)
		{
			UG_ASSERT(i < m_vpDependData.size(), "Input not needed");
			UG_ASSERT(m_vpDependData[i] != NULL, "Input invalid");
			UG_ASSERT(i < m_vvSeriesID.size(), "Input Series missing");
			UG_ASSERT(s < m_vvSeriesID[i].size(), "Input Series invalid");
			return m_vpDependData[i]->deriv(m_vvSeriesID[i][s], ip, fct, dof);
		}

	protected:
	///	get series id's from input data and resize data
		virtual void adjust_global_ips_and_data(const std::vector<size_t>& vNumIP)
		{
		//	 we need a series id for all inputs
			m_vvSeriesID.resize(m_vpIPData.size());

		//	loop inputs
			for(size_t i = 0; i < m_vpIPData.size(); ++i)
			{
			//	resize series ids
				m_vvSeriesID[i].resize(vNumIP.size());

			//	skip unset data
				UG_ASSERT(m_vpIPData[i] != NULL, "No Input set, but requested.");

			//	request local ips for all series at input data
				for(size_t s = 0; s < vNumIP.size(); ++s)
				{
					switch(this->dim_local_ips())
					{
						case 1:
							m_vvSeriesID[i][s] =
								m_vpIPData[i]->template register_local_ip_series<1>
											(this->template local_ips<1>(s), num_ip(s));
							break;
						case 2:
							m_vvSeriesID[i][s] =
								m_vpIPData[i]->template register_local_ip_series<2>
											(this->template local_ips<2>(s), num_ip(s));
							break;
						case 3:
							m_vvSeriesID[i][s] =
								m_vpIPData[i]->template register_local_ip_series<3>
											(this->template local_ips<3>(s), num_ip(s));
							break;
						default: throw(UGFatalError("Dimension not supported."));
					}
				}
			}

		//	resize data fields
			DependentIPData<TData, dim>::adjust_global_ips_and_data(vNumIP);
		}

	protected:
	///	data input
		std::vector<IPData<TDataIn, dim>*> m_vpIPData;

	///	data input casted to dependend data
		std::vector<DependentIPData<TDataIn, dim>*> m_vpDependData;

	///	series id the linker uses to get data from input
		std::vector<std::vector<size_t> > m_vvSeriesID;

	///	common functions the data depends on
		FunctionGroup m_commonFctGroup;

	///	Function mapping for each input relative to common FunctionGroup
		std::vector<FunctionIndexMapping> m_vMap;
};

} // end namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_LINKER__ */
