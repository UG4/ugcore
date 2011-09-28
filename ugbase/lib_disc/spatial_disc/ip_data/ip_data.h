/*
 * ip_data.h
 *
 *  Created on: 17.12.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__IP_DATA__IP_DATA__
#define __H__UG__LIB_DISC__SPATIAL_DISC__IP_DATA__IP_DATA__

#include <vector>
#include "lib_disc/common/local_algebra.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
//	IP Data
////////////////////////////////////////////////////////////////////////////////

/// Base class for IP Data
/**
 * This is the base class for all coupled data at integration point. It handles
 * the set of local integration points and stores the current time.
 */
class IIPData
{
	public:
	///	default constructor
		IIPData();

	///	set evaluation time
		void set_time(number time) {m_time = time;}

	///	get evaluation time
		number time() const {return m_time;}

	///	returns the number of ip series
		size_t num_series() const {return m_vNumIP.size();}

	/// returns the number of integration points
		size_t num_ip(size_t s) const {UG_ASSERT(s < num_series(), "Invalid series"); return m_vNumIP[s];}

	///	clear all data
		void clear();

	///	set local positions, returns series id
		template <int ldim>
		size_t register_local_ip_series(const MathVector<ldim>* vPos, size_t numIP);

	///	returns current local ip dimension
		int dim_local_ips() const {return m_locPosDim;}

	///	returns local ips
		template <int ldim>
		const MathVector<ldim>* local_ips(size_t s) const;

	/// returns local ip
		template <int ldim>
		const MathVector<ldim>& local_ip(size_t s, size_t ip) const;

	///	returns if data is constant
		virtual bool constant_data() const {return false;}

	///	returns if data depends on solution
		virtual bool zero_derivative() const {return true;}

	///	number of other Data this data depends on
		virtual size_t num_needed_data() const {return 0;}

	///	return needed data
		virtual IIPData* needed_data(size_t i) {return NULL;}

	/// compute values (and derivatives iff compDeriv == true)
		virtual bool compute(bool bDeriv = false) = 0;

	///	virtual desctructor
		virtual ~IIPData() {};

	protected:
	///	callback invoked after local ips have been added to the series
	/**
	 * This callback is invoked when local ips have been added. It can
	 * be used by derived classes to react on this fact, e.g. to forward the
	 * local_ips or to adapt data field sizes.
	 * Note: The number of series can only be increased and the number of ips
	 * 		 for a series can not be changed once set. This is important to
	 * 		 allow derived classes to only increase their data fields as well,
	 * 		 leaving accessing pointers invariant. If the local ip series must
	 * 		 be changed, this can only be done through a clear(), that will
	 * 		 involke the local_ips_to_be_cleared() callback, and adding all local
	 * 		 series again.
	 */
		virtual void local_ips_added() = 0;

	///	callback invoked, when local ips are cleared
		virtual void local_ips_to_be_cleared() = 0;

	///	help function to get local ips
		std::vector<const MathVector<1>*>& get_local_ips(Int2Type<1>) {return m_pvLocIP1d;}
		std::vector<const MathVector<2>*>& get_local_ips(Int2Type<2>) {return m_pvLocIP2d;}
		std::vector<const MathVector<3>*>& get_local_ips(Int2Type<3>) {return m_pvLocIP3d;}
		const std::vector<const MathVector<1>*>& get_local_ips(Int2Type<1>) const {return m_pvLocIP1d;}
		const std::vector<const MathVector<2>*>& get_local_ips(Int2Type<2>) const {return m_pvLocIP2d;}
		const std::vector<const MathVector<3>*>& get_local_ips(Int2Type<3>) const {return m_pvLocIP3d;}

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

/// World dimension based IP Data
/**
 * This class is the dimension dependent base class for all integration point data.
 * It provides the access to the data and handles the global integration
 * points.
 *
 * \tparam	dim		world dimension
 */
template <int dim>
class IIPDimData : virtual public IIPData
{
	public:
	///	set global positions
		void set_global_ips(size_t s, const MathVector<dim>* vPos, size_t numIP);

	///	returns global ips
		const MathVector<dim>* ips(size_t s) const {check_s(s); return m_vvGlobPos[s];}

	/// returns global ip
		const MathVector<dim>& ip(size_t s, size_t ip) const{check_s_ip(s,ip); return m_vvGlobPos[s][ip];}

	protected:
	///	callback invoked after global ips have been changed
	/**
	 * This callback is invoked when the global ips have been changed. It can
	 * be used by derived classes to react on this fact, e.g. to forward the
	 * global_ips.
	 */
		virtual void global_ips_changed(size_t s, const MathVector<dim>* vPos, size_t numIP) {}

	///	implement callback, called when num of local IPs changes
		virtual void local_ips_added()
		{
		//	adjust Global positions pointer to current number of series
			m_vvGlobPos.resize(num_series());
		}

	///	implement callback, invoked when local ips are cleared
		virtual void local_ips_to_be_cleared() {m_vvGlobPos.clear();}

	///	checks in debug mode the correct usage of indices
		inline void check_s(size_t s) const;

	///	checks in debug mode the correct usage of indices
		inline void check_s_ip(size_t s, size_t ip) const;

	/// global ips
		std::vector<const MathVector<dim>*> m_vvGlobPos;
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
class IPData : public IIPDimData<dim>
{
	public:
	///	type of base class
		typedef IIPDimData<dim> base_type;

	///	explicitly forward some functions
		using base_type::num_series;
		using base_type::num_ip;

	public:
	///	returns the value at ip
		const TData& value(size_t s, size_t ip) const
			{check_series_ip(s,ip); return m_vvValue[s][ip];}

	///	returns all values for a series
		const TData* values(size_t s) const
			{check_series(s); return &(m_vvValue[s][0]);}

	///	returns the value at ip
		TData& value(size_t s, size_t ip)
			{check_series_ip(s,ip);return m_vvValue[s][ip];}

	///	returns all values for a series
		TData* values(size_t s)
			{check_series(s); return &(m_vvValue[s][0]);}

	///	destructor
		~IPData() {local_ips_to_be_cleared();}

	protected:
	///	checks in debug mode the correct index
		inline void check_series(size_t s) const;

	///	checks in debug mode the correct index
		inline void check_series_ip(size_t s, size_t ip) const;

	///	resizes the data field, when local ip changed signaled
		virtual void local_ips_added();

	///	free the data field memory and set series to zero
		virtual void local_ips_to_be_cleared();

	protected:
	/// data at ip (size: (0,...num_series-1) x (0,...,num_ip-1))
		std::vector<TData*> m_vvValue;
};

////////////////////////////////////////////////////////////////////////////////
//	Dependent IP Data
////////////////////////////////////////////////////////////////////////////////

/// Base class for solution dependent Data
/**
 * All classes that derive from this class may depend on the unknwon numerical
 * solution, when computing the data.
 */
class IDependentIPData : virtual public IIPData
{
	protected:
	///	type of local vector
		typedef IElemDisc::local_vector_type local_vector_type;

	public:
	///	default constructor
		IDependentIPData() : m_bCompNeedsSol(false) {};

	///	resize arrays
		virtual void resize(const LocalIndices& ind,
		                    const FunctionIndexMapping& map) = 0;

	///	computation of data depending on current solution
		virtual bool compute(bool bDeriv)
		{
			UG_LOG("ERROR in 'IDependentIPData::compute': No implementation found.\n");
			return false;
		}

	///	computation of data depending on current solution
		virtual bool compute(const local_vector_type& u, bool bDeriv)
		{
			UG_LOG("ERROR in 'IDependentIPData::compute': No implementation found.\n");
			return false;
		}

	///	returns if the computation needs the current solution
		bool comp_needs_sol() const {return m_bCompNeedsSol;}

	///	updates the function group (needed for Linker)
		virtual bool update_function_group() {return true;}

	/// set	Function Group of functions (by copy)
		void set_function_group(const FunctionGroup& fctGrp) {m_fctGrp = fctGrp;}

	///	Function Group of functions
		const FunctionGroup& get_function_group() const {return m_fctGrp;}

	///	number of fuctions this export depends on
		size_t num_fct() const {return m_fctGrp.num_fct();}

	///	returns if the dependent data is ready for evaluation
		virtual bool is_ready() const = 0;

	///	virtual destructor
		virtual ~IDependentIPData() {}

	protected:
	/// functions the data depends on
		FunctionGroup m_fctGrp;

	///	flag if computation needs the solution
		bool m_bCompNeedsSol;
};

/// Dependent IP Data
/**
 * This class extends the IPData by the derivatives of the data w.r.t. to
 * unknown solutions.
 */
template <typename TData, int dim>
class DependentIPData : public IPData<TData, dim>,
						public IDependentIPData
{
	public:
	///	Base class type
		typedef IPData<TData, dim> base_type;

	//	explicitly forward methods of IIPData
		using base_type::num_series;
		using base_type::num_ip;
		using base_type::local_ips;

	public:
	/// number of shapes for local function
		size_t num_sh(size_t s, size_t fct) const
		{
			const size_t ip = 0; check_s_ip_fct(s,ip,fct);
			return m_vvvvDeriv[s][ip][fct].size();
		}

	///	returns the derivative of the local function, at ip and for a dof
		const TData& deriv(size_t s, size_t ip, size_t fct, size_t dof) const
			{check_s_ip_fct_dof(s,ip,fct,dof);return m_vvvvDeriv[s][ip][fct][dof];}

	///	returns the derivative of the local function, at ip and for a dof
		TData& deriv(size_t s, size_t ip, size_t fct, size_t dof)
			{check_s_ip_fct_dof(s,ip,fct,dof);return m_vvvvDeriv[s][ip][fct][dof];}

	///	returns the derivatives of the local function, at ip
		TData* deriv(size_t s, size_t ip, size_t fct)
			{check_s_ip_fct(s,ip,fct);return &(m_vvvvDeriv[s][ip][fct][0]);}

	///	returns the derivatives of the local function, at ip
		const TData* deriv(size_t s, size_t ip, size_t fct) const
			{check_s_ip_fct(s,ip,fct);return &(m_vvvvDeriv[s][ip][fct][0]);}

	///	resize lin defect arrays
		virtual void resize(const LocalIndices& ind, const FunctionIndexMapping& map);

	///	sets all derivative values to zero
		void clear_derivative_values();

	///	returns if the dependent data is ready for evaluation
		virtual bool is_ready() const = 0;

	protected:
	///	checks in debug mode the correct usage of indices
		inline void check_s_ip(size_t s, size_t ip) const;

	///	checks in debug mode the correct usage of indices
		inline void check_s_ip_fct(size_t s, size_t ip, size_t fct) const;

	///	checks in debug mode the correct usage of indices
		inline void check_s_ip_fct_dof(size_t s, size_t ip, size_t fct, size_t dof) const;

	///	resizes the derivative field when local ip change is signaled
		virtual void local_ips_added();

	protected:
	// 	Data (size: (0,...,num_series-1) x (0,...,num_ip-1) x (0,...,num_fct-1) x (0,...,num_sh(fct) )
	///	Derivatives
		std::vector<std::vector<std::vector<std::vector<TData> > > > m_vvvvDeriv;
};

} // end namespace ug

//include implementation
#include "ip_data_impl.h"

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__IP_DATA__IP_DATA__ */
