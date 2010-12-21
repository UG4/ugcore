/*
 * ip_data.h
 *
 *  Created on: 17.12.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__IP_DATA__IP_DATA__
#define __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__IP_DATA__IP_DATA__

namespace ug{

/// Base class for IP Data
/**
 * Base class for all ip data
 */
class IIPData
{
	public:
		IIPData() :
			m_time(0.0)
		{
			m_vNumIP.clear();
			m_locPosDim = -1;
			m_pvLocIP1d.clear(); m_pvLocIP2d.clear(); m_pvLocIP3d.clear();
		}

	///	set evaluation time
		void set_time(number time) {m_time = time;}

	///	get evaluation time
		number time() const {return m_time;}

	///	returns the number of ip series
		size_t num_series() const {return m_vNumIP.size();}

	/// returns the number of integration points
		size_t num_ip(size_t s) const
		{
			UG_ASSERT(s < num_series(), "Invalid series");
			return m_vNumIP[s];
		}

	///	clear all data
		void clear_ips()
		{
			m_vNumIP.clear();
			m_locPosDim = -1;
			m_pvLocIP1d.clear(); m_pvLocIP2d.clear(); m_pvLocIP3d.clear();
			adjust_global_ips_and_data(m_vNumIP);
		}

	///	set local positions
	/**
	 *
	 * \return number of series
	 */
		template <int ldim>
		size_t register_local_ip_series(const MathVector<ldim>* vPos, size_t numIP)
		{
		//	check, that dimension is ok.
			if(m_locPosDim == -1) m_locPosDim = ldim;
			else if(m_locPosDim != ldim)
				throw(UGFatalError("Local IP dimension conflict"));

		//	get local positions
			std::vector<const MathVector<ldim>*>& vvIP = get_local_ips(Int2Type<ldim>());

		//	search for ips
			for(size_t s = 0; s < vvIP.size(); ++s)
			{
			//	return series number iff exists
				if(vvIP[s] == vPos) return s;
			}

		//	if series not yet registered, add it
			vvIP.push_back(vPos);
			m_vNumIP.push_back(numIP);

		//	resize global ips and data
			adjust_global_ips_and_data(m_vNumIP);

		//	return new seris id
			return m_vNumIP.size() - 1;
		}

	///	returns local ips
		template <int ldim>
		const MathVector<ldim>* local_ips(size_t s) const
		{
		//	check, that dimension is ok.
			if(m_locPosDim != ldim)
				throw(UGFatalError("Local IP dimension conflict"));

			UG_ASSERT(s < num_series(), "Wrong series id");
			UG_ASSERT( (const_cast<IIPData*>(this)->get_local_ips(Int2Type<ldim>()))[s] != NULL, "Zero local ip pointer.");

			return (const_cast<IIPData*>(this)->get_local_ips(Int2Type<ldim>()))[s];
		}

	/// returns local ip
		template <int ldim>
		const MathVector<ldim>& local_ip(size_t s, size_t ip) const
		{
		//	check, that dimension is ok.
			if(m_locPosDim != ldim)
				throw(UGFatalError("Local IP dimension conflict"));

			UG_ASSERT(s < num_series(), "Wrong series id");
			UG_ASSERT(ip < num_ip(s), "Invalid index.");

			return get_local_ips(Int2Type<ldim>())[s][ip];
		}

	///	returns if data is constant
		virtual bool constant_data() const {return false;}

	/// compute values (and derivatives iff compDeriv == true)
		virtual void compute(bool compDeriv = false) = 0;

	///	virtual desctructor
		~IIPData() {};

	protected:
	///	callback to adjust global ips and data
		virtual void adjust_global_ips_and_data(const std::vector<size_t>& vNumIP) = 0;

	///	help function to get local ips
		std::vector<const MathVector<1>*>& get_local_ips(Int2Type<1>) {return m_pvLocIP1d;}
		std::vector<const MathVector<2>*>& get_local_ips(Int2Type<2>) {return m_pvLocIP2d;}
		std::vector<const MathVector<3>*>& get_local_ips(Int2Type<3>) {return m_pvLocIP3d;}

	protected:
	/// number of evaluation points (-1 indicates no ips set)
		std::vector<size_t> m_vNumIP;

	/// dimension of local position (-1 indicates no dim set)
		int m_locPosDim;

	/// local ips of dimension 1d-3d
		std::vector<const MathVector<1>*> m_pvLocIP1d;
		std::vector<const MathVector<2>*> m_pvLocIP2d;
		std::vector<const MathVector<3>*> m_pvLocIP3d;

	///	time for evaluation
		number m_time;
};

/**
 * \tparam	TData	Data
 * \tparam	dim		world dimension
 */
template <typename TData, int dim>
class IPData : public IIPData
{
	public:
	///	returns the value at ip
		const TData& value(size_t s, size_t ip) const
		{
			UG_ASSERT(s < num_series(), "Wrong series id "<<s);
			UG_ASSERT(s < m_vvValue.size(), "Invalid index "<<s);
			UG_ASSERT(ip < num_ip(s), "Invalid index "<<ip);
			UG_ASSERT(ip < m_vvValue[s].size(), "Invalid index "<<ip);
			return m_vvValue[s][ip];
		}

	///	returns all values for a series
		const TData* values(size_t s) const
		{
			UG_ASSERT(s < num_series(), "Wrong series id");
			UG_ASSERT(s < m_vvValue.size(), "Invalid index");
			return &(m_vvValue[s][0]);
		}

	///	returns the value at ip
		TData& value(size_t s, size_t ip)
		{
			UG_ASSERT(s < num_series(), "Wrong series id"<<s);
			UG_ASSERT(s < m_vvValue.size(), "Invalid index "<<s);
			UG_ASSERT(ip < num_ip(s), "Invalid index "<<ip);
			UG_ASSERT(ip < m_vvValue[s].size(), "Invalid index "<<ip);
			return m_vvValue[s][ip];
		}

	///	returns all values for a series
		TData* values(size_t s)
		{
			UG_ASSERT(s < num_series(), "Wrong series id");
			UG_ASSERT(s < m_vvValue.size(), "Invalid index");
			return &(m_vvValue[s][0]);
		}

	///	set global positions
		void set_global_ips(size_t s, const MathVector<dim>* vPos, size_t numIP)
		{
			UG_ASSERT(s < num_series(), "Wrong series id");

		//	check number of ips (must match local ip number)
			if(numIP != num_ip(s))
				throw(UGFatalError("Num ip does not match."));

		//	remember global positions
			m_vvGlobPos[s] = vPos;
		}

	///	returns global ips
		const MathVector<dim>* ips(size_t s) const
		{
			UG_ASSERT(s < num_series(), "Wrong series id");
			return m_vvGlobPos[s];
		}

	/// returns global ip
		const MathVector<dim>& ip(size_t s, size_t ip) const
		{
			UG_ASSERT(s < num_series(), "Wrong series id.");
			UG_ASSERT(s < m_vvGlobPos.size(), "Invalid index.");
			UG_ASSERT(ip < num_ip(s), "Invalid index.");
			UG_ASSERT(m_vvGlobPos[s] != NULL, "Local IP not set.");

			return m_vvGlobPos[s][ip];
		}

	protected:
		virtual void adjust_global_ips_and_data(const std::vector<size_t>& vNumIP)
		{
		//	adjust Global positions pointer
			m_vvGlobPos.resize(vNumIP.size());

		//	adjust data arrays
			m_vvValue.resize(vNumIP.size());
			for(size_t s = 0; s < vNumIP.size(); ++s)
				m_vvValue[s].resize(vNumIP[s]);
		}

	protected:
		/// global ips
		std::vector<const MathVector<dim>*> m_vvGlobPos;

		/// data at ip
		std::vector<std::vector<TData> > m_vvValue;
};

} // end namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__IP_DATA__IP_DATA__ */
