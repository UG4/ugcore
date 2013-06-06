/*
 * user_data.h
 *
 *  Created on: 17.12.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__USER_DATA__USER_DATA__
#define __H__UG__LIB_DISC__SPATIAL_DISC__USER_DATA__USER_DATA__

#include <vector>
#include "lib_disc/common/local_algebra.h"
#include "lib_disc/common/function_group.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
//	UserData Info
////////////////////////////////////////////////////////////////////////////////

/// base class providing runtime-info on dimension and type
class UserDataInfo {
	public:
	///	returns dimension
		virtual int get_dim() const = 0;

	///	returns type of data as string (e.g. "Number", "Vector", "Matrix")
		virtual std::string type() const = 0;

	///	returns if provided data is continuous over geometric object boundaries
		virtual bool continuous() const = 0;

	/// virtual destructor
		virtual ~UserDataInfo() {}

	public:
	///	returns if grid function is needed for evaluation
		virtual bool requires_grid_fct() const = 0;

	///	sets the function pattern for a possibly needed grid function
		virtual void set_function_pattern(ConstSmartPtr<FunctionPattern> fctPatt) {
			m_fctGrp.set_function_pattern(fctPatt);
		}

	///	Function Group of functions
		const FunctionGroup& function_group() const {return m_fctGrp;}

	///	get function mapping
		const FunctionIndexMapping& map() const{return m_map;}

	///	number of functions this export depends on
		size_t num_fct() const {return m_map.num_fct();}

	protected:
	/// functions the data depends on
		FunctionGroup m_fctGrp;

	///	Mapping for import fct
		FunctionIndexMapping m_map;
};

////////////////////////////////////////////////////////////////////////////////
//	UserData
////////////////////////////////////////////////////////////////////////////////

// Traits
template <typename TData>
struct user_data_traits{static std::string name() 							{return "(unknown)";}};
template <>
struct user_data_traits<number>{static std::string name() 					{return "Number";}};
template <std::size_t dim>
struct user_data_traits< MathVector<dim> >{static std::string name() 		{return "Vector";}};
template <std::size_t dim>
struct user_data_traits< MathMatrix<dim,dim> >{static std::string name() 	{return "Matrix";}};
template <std::size_t dim>
struct user_data_traits< MathTensor<4,dim> >{static std::string name() 		{return "Tensor4";}};

/// Type based UserData
/**
 * This class is the base class for all integration point data for a templated
 * type. It provides the access to the data and handles the global integration
 * points.
 *
 * \tparam	TData	Data
 * \tparam	dim		world dimension
 * \tparam	TRet	Type of return flag (bool or void)
 */
template <typename TData, int dim, typename TRet = void>
class UserData : virtual public UserDataInfo
{
	public:
	///	returns dimension
		int get_dim() const {return dim;}

	///	returns type of data as string (e.g. "Number", "Vector", "Matrix")
		std::string type() const {return user_data_traits<TData>::name();}

	///	returns if provided data is continuous over geometric object boundaries
		virtual bool continuous() const = 0;

	///	returns if grid function is needed for evaluation
		virtual bool requires_grid_fct() const = 0;

	public:
	///	returns value for a global position
		virtual TRet operator() (TData& value,
								 const MathVector<dim>& globIP,
								 number time, int si) const = 0;

	///	returns value for global positions
		virtual void operator()(TData vValue[],
								const MathVector<dim> vGlobIP[],
								number time, int si, const size_t nip) const = 0;

	///	returns value for local and global position
	///	\{
		TRet operator() (TData& value,
						 const MathVector<dim>& globIP,
						 number time, int si,
						 GeometricObject* elem,
						 const MathVector<dim> vCornerCoords[],
						 const MathVector<1>& locIP,
						 LocalVector* u) const {
			operator()(&value, &globIP, time, si, elem, vCornerCoords, &locIP, 1, u);
		}

		TRet operator() (TData& value,
						 const MathVector<dim>& globIP,
						 number time, int si,
						 GeometricObject* elem,
						 const MathVector<dim> vCornerCoords[],
						 const MathVector<2>& locIP,
						 LocalVector* u) const {
			operator()(&value, &globIP, time, si, elem, vCornerCoords, &locIP, 1, u);
		}

		TRet operator() (TData& value,
						 const MathVector<dim>& globIP,
						 number time, int si,
						 GeometricObject* elem,
						 const MathVector<dim> vCornerCoords[],
						 const MathVector<3>& locIP,
	                     LocalVector* u) const {
			operator()(&value, &globIP, time, si, elem, vCornerCoords, &locIP, 1, u);
		}
	///	\}

	///	returns values for local and global positions
	///	\{
		virtual void operator()(TData vValue[],
		                        const MathVector<dim> vGlobIP[],
		                        number time, int si,
		                        GeometricObject* elem,
		                        const MathVector<dim> vCornerCoords[],
		                        const MathVector<1> vLocIP[],
		                        const size_t nip,
		                        LocalVector* u,
		                        const MathMatrix<1, dim>* vJT = NULL) const = 0;

		virtual void operator()(TData vValue[],
		                        const MathVector<dim> vGlobIP[],
		                        number time, int si,
		                        GeometricObject* elem,
		                        const MathVector<dim> vCornerCoords[],
		                        const MathVector<2> vLocIP[],
		                        const size_t nip,
		                        LocalVector* u,
		                        const MathMatrix<2, dim>* vJT = NULL) const = 0;

		virtual void operator()(TData vValue[],
		                        const MathVector<dim> vGlobIP[],
		                        number time, int si,
		                        GeometricObject* elem,
		                        const MathVector<dim> vCornerCoords[],
		                        const MathVector<3> vLocIP[],
		                        const size_t nip,
		                        LocalVector* u,
		                        const MathMatrix<3, dim>* vJT = NULL) const = 0;
	///	\}
};
////////////////////////////////////////////////////////////////////////////////
//	UserData Interface
////////////////////////////////////////////////////////////////////////////////

/// Base class for UserData
/**
 * This is the base class for all coupled data at integration point. It handles
 * the set of local integration points and stores the current time.
 *
 * \tparam	dim		world dimension
 */
template <int dim>
class ICplUserData : virtual public UserDataInfo
{
	public:
	///	default constructor
		ICplUserData();

	///	clear all data
		void clear();

	public:
	///	set the subset of evaluation
		void set_subset(int si) {m_si = si;}

	///	returns the subset of evaluation
		int subset() const {return m_si;}

	///	set evaluation time
		void set_times(const std::vector<number>& vTime) {m_vTime = vTime;}

	/// sets the current time point
		void set_time_point(size_t timePoint) {m_timePoint = timePoint;}

	///	get evaluation time
		number time() const {return m_vTime[m_timePoint];}

	public:
	///	returns if data is constant
		virtual bool constant() const {return false;}

	///	number of other Data this data depends on
		virtual size_t num_needed_data() const {return 0;}

	///	return needed data
		virtual SmartPtr<ICplUserData> needed_data(size_t i) {return NULL;}

	/// compute values (and derivatives iff compDeriv == true)
		virtual void compute(LocalVector* u,
		                     GeometricObject* elem,
		                     const MathVector<dim> vCornerCoords[],
		                     bool bDeriv = false) = 0;

	///	returns if the dependent data is ready for evaluation
		virtual void check_setup() const {}

	///	virtual desctructor
		virtual ~ICplUserData() {};

	public:
	///	returns if data depends on solution
		virtual bool zero_derivative() const {return true;}

	///	resize arrays
		virtual void update_dof_sizes(const LocalIndices& ind) {}

	public:
	///	returns the number of ip series
		size_t num_series() const {return m_vNumIP.size();}

	/// returns the number of integration points
		size_t num_ip(size_t s) const {UG_ASSERT(s < num_series(), "Invalid series"); return m_vNumIP[s];}

	///	set local positions, returns series id
	/**
	 * This method registers a local ip series. If the position of points may
	 * change during computations, this can be specified.
	 * IMPORTANT: the memory of the local ip values must remain valid until the
	 * UserData is deleted.
	 *
	 * \returns size_t 		series id
	 */
		template <int ldim>
		size_t register_local_ip_series(const MathVector<ldim>* vPos,
										const size_t numIP,
										bool bMayChange = true);

	///	sets new local ip positions for a local ip series
	/**
	 * This method set new local positions for an already registered ip series.
	 * Of coarse this is only possible for a ip series, that has the bMayChange
	 * flag set to true.
	 */
		template <int ldim>
		void set_local_ips(const size_t seriesId, const MathVector<ldim>* vPos,
						   const size_t numIP);

	///	returns current local ip dimension
		int dim_local_ips() const {return m_locPosDim;}

	///	returns local ips
		template <int ldim>
		const MathVector<ldim>* local_ips(size_t s) const;

	/// returns local ip
		template <int ldim>
		const MathVector<ldim>& local_ip(size_t s, size_t ip) const;

	///	set global positions
		void set_global_ips(size_t s, const MathVector<dim>* vPos, size_t numIP);

	///	returns global ips
		const MathVector<dim>* ips(size_t s) const {check_s(s); return m_vvGlobPos[s];}

	/// returns global ip
		const MathVector<dim>& ip(size_t s, size_t ip) const{check_s_ip(s,ip); return m_vvGlobPos[s][ip];}

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
	 * 		 invoke the local_ip_series_to_be_cleared() callback, and adding all local
	 * 		 series again.
	 */
		virtual void local_ip_series_added(const size_t seriesID){m_vvGlobPos.resize(seriesID+1);}

	///	callback invoked, if a local ip series has been changed
		virtual void local_ips_changed(const size_t seriesID, const size_t newNumIP) = 0;

	///	callback invoked, when local ips are cleared
		virtual void local_ip_series_to_be_cleared() {m_vvGlobPos.clear();}

	///	callback invoked after global ips have been changed
	/**
	 * This callback is invoked when the global ips have been changed. It can
	 * be used by derived classes to react on this fact, e.g. to forward the
	 * global_ips.
	 */
		virtual void global_ips_changed(const size_t seriesID, const MathVector<dim>* vPos, const size_t numIP) {};

	///	checks in debug mode the correct usage of indices
		inline void check_s(size_t s) const;

	///	checks in debug mode the correct usage of indices
		inline void check_s_ip(size_t s, size_t ip) const;

	protected:
	///	help function to get local ips
		std::vector<const MathVector<1>*>& get_local_ips(Int2Type<1>) {return m_pvLocIP1d;}
		std::vector<const MathVector<2>*>& get_local_ips(Int2Type<2>) {return m_pvLocIP2d;}
		std::vector<const MathVector<3>*>& get_local_ips(Int2Type<3>) {return m_pvLocIP3d;}
		const std::vector<const MathVector<1>*>& get_local_ips(Int2Type<1>) const {return m_pvLocIP1d;}
		const std::vector<const MathVector<2>*>& get_local_ips(Int2Type<2>) const {return m_pvLocIP2d;}
		const std::vector<const MathVector<3>*>& get_local_ips(Int2Type<3>) const {return m_pvLocIP3d;}

	protected:
	///	flags if local ips may change
		std::vector<bool> m_vMayChange;

	/// number of evaluation points (-1 indicates no ips set)
		std::vector<size_t> m_vNumIP;

	/// dimension of local position (-1 indicates no dim set)
		int m_locPosDim;

	/// local ips of dimension 1d-3d
		std::vector<const MathVector<1>*> m_pvLocIP1d;
		std::vector<const MathVector<2>*> m_pvLocIP2d;
		std::vector<const MathVector<3>*> m_pvLocIP3d;

	/// global ips
		std::vector<const MathVector<dim>*> m_vvGlobPos;

	///	time for evaluation
		std::vector<number> m_vTime;

	///	current time point
		size_t m_timePoint;

	///	subset for evaluation
		int m_si;
};

////////////////////////////////////////////////////////////////////////////////
//	UserData
////////////////////////////////////////////////////////////////////////////////

// predeclaration
template <typename TData, int dim> class DataImport;

/// Type based UserData
/**
 * This class is the base class for all integration point data for a templated
 * type. It provides the access to the data and handles the global integration
 * points.
 *
 * \tparam	TData	Data
 * \tparam	dim		world dimension
 * \tparam	TRet	Type of return flag (bool or void)
 */
template <typename TData, int dim, typename TRet = void>
class CplUserData : public ICplUserData<dim>, public UserData<TData,dim,TRet>
{
	public:
	///	type of base class
		typedef ICplUserData<dim> base_type;

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

	///	returns flag, if data is evaluated (for conditional data)
		bool defined(size_t s, size_t ip) const
			{check_series_ip(s,ip); return m_vvBoolFlag[s][ip];}

	///	destructor
		~CplUserData() {local_ip_series_to_be_cleared();}

	///	register external callback, invoked when data storage changed
		void register_storage_callback(DataImport<TData,dim>* obj, void (DataImport<TData,dim>::*func)());

	///	register all callbacks registered by class
		void unregister_storage_callback(DataImport<TData,dim>* obj);

	protected:
	///	checks in debug mode the correct index
		inline void check_series(size_t s) const;

	///	checks in debug mode the correct index
		inline void check_series_ip(size_t s, size_t ip) const;

	///	resizes the data field, when local ip changed signaled
		virtual void local_ip_series_added(const size_t seriesID);

	///	free the data field memory and set series to zero
		virtual void local_ip_series_to_be_cleared();

	///	implement callback, called when local IPs changed
		virtual void local_ips_changed(const size_t seriesID, const size_t newNumIP);

	///	callback, invoked when storage of data has changed for a series
		virtual void value_storage_changed(const size_t seriesID) {}

	///	calls are registered external storage callbacks
		void call_storage_callback() const;

	private:
	/// data at ip (size: (0,...num_series-1) x (0,...,num_ip-1))
		std::vector<std::vector<TData> > m_vvValue;

	/// bool flag at ip (size: (0,...num_series-1) x (0,...,num_ip-1))
		std::vector<std::vector<bool> > m_vvBoolFlag;

	///	registered callbacks
//		typedef void (DataImport<TData,dim>::*CallbackFct)();
		typedef boost::function<void ()> CallbackFct;
		std::vector<std::pair<DataImport<TData,dim>*, CallbackFct> > m_vCallback;

};

////////////////////////////////////////////////////////////////////////////////
//	Dependent UserData
////////////////////////////////////////////////////////////////////////////////

/// Dependent UserData
/**
 * This class extends the UserData by the derivatives of the data w.r.t. to
 * unknown solutions.
 */
template <typename TData, int dim>
class DependentUserData : public CplUserData<TData, dim>
{
	public:
	///	Base class type
		typedef CplUserData<TData, dim> base_type;

	//	explicitly forward methods of ICplUserData
		using base_type::num_series;
		using base_type::num_ip;
		using base_type::local_ips;

	public:
	///	default constructor
		DependentUserData() {}

	///	sets the associated symbolic functions
	/// \{
		DependentUserData(const char* symbFct) {set_functions(symbFct);}
		DependentUserData(const std::string& symbFct) {set_functions(symbFct);}
		DependentUserData(const std::vector<std::string>& symbFct) {set_functions(symbFct);}
	/// \}

	public:
	/// number of shapes for local function
		size_t num_sh(size_t fct) const
		{
			UG_ASSERT(fct < m_vvNumDoFPerFct.size(), "Wring index");
			return m_vvNumDoFPerFct[fct];
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

	///	sets all derivative values to zero
		static void set_zero(std::vector<std::vector<TData> > vvvDeriv[], const size_t nip);

	public:
	///	returns that data depends on solution
		virtual bool zero_derivative() const {return false;}

	///	returns if grid function is needed for evaluation
		virtual bool requires_grid_fct() const {return true;}

	///	resize lin defect arrays
		virtual void update_dof_sizes(const LocalIndices& ind);

	///	sets the associated function pattern
		virtual void set_function_pattern(ConstSmartPtr<FunctionPattern> fctPatt);

	///	sets the associated symbolic functions
	/// \{
		void set_functions(const char* symbFct);
		void set_functions(const std::string& symbFct);
		void set_functions(const std::vector<std::string>& symbFct);
	/// \}

	protected:
	///	extracts the function group
		void extract_fct_grp();

	protected:
	///	string of symbolic functions required
		std::vector<std::string> m_SymbFct;

	protected:
	///	checks in debug mode the correct usage of indices
		inline void check_s_ip(size_t s, size_t ip) const;

	///	checks in debug mode the correct usage of indices
		inline void check_s_ip_fct(size_t s, size_t ip, size_t fct) const;

	///	checks in debug mode the correct usage of indices
		inline void check_s_ip_fct_dof(size_t s, size_t ip, size_t fct, size_t dof) const;

	///	resizes the derivative field when local ip change is signaled
		virtual void local_ip_series_added(const size_t seriesID);

	///	implement callback, called when local IPs changed
		virtual void local_ips_changed(const size_t seriesID, const size_t newNumIP);

	///	implement callback, invoked when local ips are cleared
		virtual void local_ip_series_to_be_cleared();

	///	resizes the derivative arrays for current number of ips.
		void resize_deriv_array();

	///	resizes the derivative arrays for current number of ips of a single series
		void resize_deriv_array(const size_t seriesID);

	protected:
	///	number of functions and their dofs
		std::vector<size_t> m_vvNumDoFPerFct;

	// 	Data (size: (0,...,num_series-1) x (0,...,num_ip-1) x (0,...,num_fct-1) x (0,...,num_sh(fct) )
	///	Derivatives
		std::vector<std::vector<std::vector<std::vector<TData> > > > m_vvvvDeriv;
};

} // end namespace ug

//include implementation
#include "user_data_impl.h"

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__USER_DATA__USER_DATA__ */
