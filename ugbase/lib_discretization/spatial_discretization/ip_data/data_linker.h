/*
 * data_linker.h
 *
 *  Created on: 12.11.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_LINKER__
#define __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_LINKER__

#include "ip_data.h"

namespace ug{

template <typename TData, int dim>
class ScalarLinker : public DependentIPData<TData, dim>
{
	public:
		ScalarLinker() : m_pIPData(NULL){}

		virtual void compute(bool compDeriv)
		{
			for(size_t s = 0; s < this->num_series(); ++s)
			{
				for(size_t ip = 0; ip < this->num_ip(s); ++ip)
				{
					this->value(s, ip) =
							1e3 + 0.2e3 * m_pIPData->value(m_vSeriesID[s], ip);
				}
			}

			if(!compDeriv || m_pIPData->zero_derivative()) return;

			for(size_t s = 0; s < this->num_series(); ++s)
				for(size_t ip = 0; ip < this->num_ip(s); ++ip)
					for(size_t fct = 0; fct < this->num_fct(); ++fct)
						for(size_t dof = 0; dof < this->num_sh(s, fct); ++dof)
						{
							this->deriv(s, ip, fct, dof) =
									0.2e3 *
									m_pDependData->deriv(m_vSeriesID[s], ip, fct, dof);
						}

		}

		virtual bool zero_derivative() const
		{
			if(m_pIPData == NULL) return true;
			else return m_pIPData->zero_derivative();
 		}

		void set_input(IPData<TData, dim>& data)
		{
			m_pIPData = &data;
			m_pDependData = dynamic_cast<DependentIPData<TData, dim>*>(&data);

			if(m_pDependData != NULL)
				this->set_function_group(m_pDependData->get_function_group());
		}

	///	number of other Data this data depends on
		virtual size_t num_needed_data() const {return 1;}

	///	return needed data
		virtual IIPData* needed_data(size_t i)
		{
			UG_ASSERT(i == 0, "Only one element needed");
			return m_pIPData;
		}

	protected:
		virtual void adjust_global_ips_and_data(const std::vector<size_t>& vNumIP)
		{
			if(m_pIPData != NULL)
			{
			//	resize series ids
				m_vSeriesID.resize(vNumIP.size());

			//	request local ips for all series at imported data
				for(size_t s = 0; s < vNumIP.size(); ++s)
				{
					switch(this->dim_local_ips())
					{
						case 1:
							m_vSeriesID[s] =
								m_pIPData->template register_local_ip_series<1>
											(this->template local_ips<1>(s), this->num_ip(s));
							break;
						case 2:
							m_vSeriesID[s] =
								m_pIPData->template register_local_ip_series<2>
											(this->template local_ips<2>(s), this->num_ip(s));
							break;
						case 3:
							m_vSeriesID[s] =
								m_pIPData->template register_local_ip_series<3>
											(this->template local_ips<3>(s), this->num_ip(s));
							break;
						default: throw(UGFatalError("Dimension not supported."));
					}
				}
			}

		//	resize data fields
			DependentIPData<TData, dim>::adjust_global_ips_and_data(vNumIP);
		}

	protected:
		IPData<TData, dim>* m_pIPData;

		DependentIPData<TData, dim>* m_pDependData;

		std::vector<size_t> m_vSeriesID;
};


} // end namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_LINKER__ */
