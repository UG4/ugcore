/*
 * ip_data.h
 *
 *  Created on: 17.12.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__IP_DATA__IP_DATA__
#define __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__IP_DATA__IP_DATA__

#include "lib_discretization/common/local_algebra.h"

namespace ug{

/// Base class for IP Data
/**
 * This is the base class for all coupled data at integration point. It handles
 * the set of local integration points and stores the current time.
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
	 * \return returns the series id, that is needed for access
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

		//	return new series id
			return m_vNumIP.size() - 1;
		}

	///	returns current local ip dimension
		int dim_local_ips() const {return m_locPosDim;}

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

	///	returns if data depends on solution
		virtual bool zero_derivative() const {return true;}

	///	number of other Data this data depends on
		virtual size_t num_needed_data() const {return 0;}

	///	return needed data
		virtual IIPData* needed_data(size_t i) {return NULL;}

	/// compute values (and derivatives iff compDeriv == true)
		virtual void compute(bool compDeriv = false) = 0;

	///	virtual desctructor
		virtual ~IIPData() {};

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

/// Type based IP Data
/**
 * This class is the base class for all integration point data for a templated
 * type. It provides the access to the data and handles the global integration
 * points.
 *
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

		/// data at ip (size: (0,...num_series-1) x (0,...,num_ip-1))
		std::vector<std::vector<TData> > m_vvValue;
};



/// Base class for solution dependent Data
/**
 * All classes that derive from this class may depend on the unknwon numerical
 * solution, when computing the data.
 */
class IDependentIPData
{
	public:
	///	resize arrays
		virtual void resize(const LocalIndices& ind) = 0;

	/// set	Function Group of functions (by copy)
		void set_function_group(const FunctionGroup& fctGrp)
			{m_fctGrp = fctGrp;}

	///	Function Group of functions
		const FunctionGroup& get_function_group() const {return m_fctGrp;}

	///	number of fuctions this export depends on
		size_t num_fct() const {return m_fctGrp.num_fct();}

	///	virtual destructor
		virtual ~IDependentIPData() {}

	protected:
	/// functions the data depends on
		FunctionGroup m_fctGrp;

};

/// Dependent IP Data
/**
 * This class extends the IPData by the derivatives of the data w.r.t. to
 * unknown solutions.
 */
template <typename TData, int dim>
class DependentIPData : public IPData<TData, dim>,
						virtual public IDependentIPData
{
	public:
		using IPData<TData, dim>::num_series;
		using IPData<TData, dim>::num_ip;
		using IPData<TData, dim>::local_ips;

	public:
	/// number of shapes for local function
		size_t num_sh(size_t s, size_t fct) const
		{
			const size_t ip = 0;
			UG_ASSERT(s < num_series(), "Wrong series id"<<s);
			UG_ASSERT(s < m_vvvDeriv.size(), "Invalid index "<<s);
			UG_ASSERT(ip < num_ip(s), "Invalid index "<<ip);
			UG_ASSERT(ip < m_vvvDeriv[s].size(), "Invalid index "<<ip);
			UG_ASSERT(fct < m_vvvDeriv[s][ip].size(), "Invalid index.");
			return m_vvvDeriv[s][ip][fct].size();
		}

	///	returns the derivative of the local function, at ip and for a dof
		const TData& deriv(size_t s, size_t ip, size_t fct, size_t dof) const
		{
			UG_ASSERT(s < num_series(), "Wrong series id"<<s);
			UG_ASSERT(s < m_vvvDeriv.size(), "Invalid index "<<s);
			UG_ASSERT(ip < num_ip(s), "Invalid index "<<ip);
			UG_ASSERT(ip < m_vvvDeriv[s].size(), "Invalid index "<<ip);
			UG_ASSERT(fct < m_vvvDeriv[s][ip].size(), "Invalid index.");
			UG_ASSERT(dof < m_vvvDeriv[s][ip][fct].size(), "Invalid index.");
			return m_vvvDeriv[s][ip][fct][dof];
		}

	///	returns the derivative of the local function, at ip and for a dof
		TData& deriv(size_t s, size_t ip, size_t fct, size_t dof)
		{
			UG_ASSERT(s < num_series(), "Wrong series id"<<s);
			UG_ASSERT(s < m_vvvDeriv.size(), "Invalid index "<<s);
			UG_ASSERT(ip < num_ip(s), "Invalid index "<<ip);
			UG_ASSERT(ip < m_vvvDeriv[s].size(), "Invalid index "<<ip);
			UG_ASSERT(fct < m_vvvDeriv[s][ip].size(), "Invalid index.");
			UG_ASSERT(dof < m_vvvDeriv[s][ip][fct].size(), "Invalid index.");
			return m_vvvDeriv[s][ip][fct][dof];
		}

	///	returns the derivatives of the local function, at ip
		TData* deriv(size_t s, size_t ip, size_t fct)
		{
			UG_ASSERT(s < num_series(), "Wrong series id"<<s);
			UG_ASSERT(s < m_vvvDeriv.size(), "Invalid index "<<s);
			UG_ASSERT(ip < num_ip(s), "Invalid index "<<ip);
			UG_ASSERT(ip < m_vvvDeriv[s].size(), "Invalid index "<<ip);
			UG_ASSERT(fct < m_vvvDeriv[s][ip].size(), "Invalid index.");
			return &(m_vvvDeriv[s][ip][fct][0]);
		}

	///	returns the derivatives of the local function, at ip
		const TData* deriv(size_t s, size_t ip, size_t fct) const
		{
			UG_ASSERT(s < num_series(), "Wrong series id"<<s);
			UG_ASSERT(s < m_vvvDeriv.size(), "Invalid index "<<s);
			UG_ASSERT(ip < num_ip(s), "Invalid index "<<ip);
			UG_ASSERT(ip < m_vvvDeriv[s].size(), "Invalid index "<<ip);
			UG_ASSERT(fct < m_vvvDeriv[s][ip].size(), "Invalid index.");
			return &(m_vvvDeriv[s][ip][fct][0]);
		}

	///	resize lin defect arrays
		virtual void resize(const LocalIndices& ind)
		{
		//	resize num fct
			for(size_t s = 0; s < m_vvvDeriv.size(); ++s)
				for(size_t ip = 0; ip < m_vvvDeriv[s].size(); ++ip)
				{
				//	resize num fct
					m_vvvDeriv[s][ip].resize(ind.num_fct());

				//	resize dofs
					for(size_t fct = 0; fct < ind.num_fct(); ++fct)
						m_vvvDeriv[s][ip][fct].resize(ind.num_dofs(fct));
				}
		}

	protected:
		virtual void adjust_global_ips_and_data(const std::vector<size_t>& vNumIP)
		{
		//	adjust data arrays
			m_vvvDeriv.resize(vNumIP.size());
			for(size_t s = 0; s < vNumIP.size(); ++s)
				m_vvvDeriv[s].resize(vNumIP[s]);

		//	resize values
			IPData<TData, dim>::adjust_global_ips_and_data(vNumIP);
		}

	protected:

	///	Derivatives
	// Data (size: (0,...,num_series-1) x (0,...,num_ip-1) x (0,...,num_fct-1) x (0,...,num_sh(fct) )
		std::vector<std::vector<std::vector<std::vector<TData> > > > m_vvvDeriv;
};



} // end namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__IP_DATA__IP_DATA__ */
