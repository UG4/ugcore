/*
 * data_import_export.h
 *
 *  Created on: 17.12.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_IMPORT_EXPORT__
#define __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_IMPORT_EXPORT__

#include "ip_data.h"
#include "lib_discretization/spatial_discretization/elem_disc/elem_disc_interface.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// Data Import
////////////////////////////////////////////////////////////////////////////////

/// Base class for data import
/**
 * An IDataImport is the base class for importing data to ElemDiscs
 */
class IDataImport
{
	typedef IElemDisc::local_matrix_type local_matrix_type;

	public:
	/// Constructor
		IDataImport(bool compLinDefect = true)
			: m_pIDependentIPData(NULL), m_bInMassPart(false),
			  m_id(-1), m_bCompLinDefect(compLinDefect)
		{}
		virtual ~IDataImport()	{}

	///	sets if import is located in mass part (for time dependent problems)
		void set_mass_part(bool bInMassPart) {m_bInMassPart = bInMassPart;}

	///	returns if import is located in mass part (for time dependent problems)
		bool in_mass_part() const {return m_bInMassPart;}

	/// returns if data is set
		virtual bool data_given() const = 0;

	/// returns if data is constant
	/**
	 * This method, returns if the connected data is constant.
	 */
		virtual bool constant_data() const = 0;

	///	returns if data depends on unknown functions
	/**
	 * This method returns if the data depends on the unknown functions.
	 */
		bool zero_derivative() const
		{
			if(m_pIDependentIPData == NULL) return true;
			else return !m_bCompLinDefect;
		}

	/// returns the connected ip data
		virtual IIPData* get_data() = 0;

	///	set function group for linearization of defect
		void set_function_group(const FunctionGroup& fctGrp){m_fctGrp = fctGrp;}

	///	get funtion group
		const FunctionGroup& get_function_group() const{return m_fctGrp;}

	/// number of functions
		size_t num_fct() const {return m_fctGrp.num_fct();}

	///	resize arrays
		virtual void resize(const LocalIndices& ind,
		                    const FunctionIndexMapping& map) = 0;

	///	add jacobian entries introduced by this import
		virtual void assemble_jacobian(local_matrix_type& J) = 0;

	protected:
	/// connected iexport
		IDependentIPData* m_pIDependentIPData;

	///	function group for linear defect
		FunctionGroup m_fctGrp;

	protected:
	///	type of local vector
		typedef IElemDisc::local_vector_type local_vector_type;

	///	type of evaluation function
		typedef bool (IElemDisc::*LinDefectFunc)(const local_vector_type& u);

	///	function pointers for all elem types
		std::vector<LinDefectFunc>	m_vLinDefectFunc;

	///	flag to indicate if import is located in mass part
		bool m_bInMassPart;

	/// current Geom Object
		int m_id;

	///	elem disc
		IElemDisc* m_pObj;

	///	indicates iff lin defect should be computed
		bool m_bCompLinDefect;

	public:
	///	sets the geometric object type
		bool set_geometric_object_type(int id);

	///	register evaluation of linear defect for a element
		template <typename TFunc>
		void reg_lin_defect_fct(int id, IElemDisc* obj, TFunc func);

	///	compute lin defect
		bool compute_lin_defect(const local_vector_type& u)
		{
			return (m_pObj->*(m_vLinDefectFunc[m_id]))(u);
		}
};

/// Data import
/**
 * A DataImport is used to import data into an ElemDisc.
 *
 * \todo some data could be cached to allow faster access than using virtual fct
 */
template <typename TData, int dim>
class DataImport : public IDataImport
{
	typedef IElemDisc::local_matrix_type local_matrix_type;

	public:
	/// Constructor
		DataImport(bool bLinDefect = true) : IDataImport(bLinDefect),
			m_seriesID(-1),	m_pIPData(NULL), m_pDependentIPData(NULL)
		{}

	///	set the user data
		void set_data(IPData<TData, dim>& data);

	/// returns the connected IIPData
		IIPData* get_data(){return m_pIPData;}

	///	returns true if data given
		virtual bool data_given() const {return !(m_pIPData == NULL);}

	///	sets the evaluation time point
		void set_time(number time) {if(m_pIPData != NULL) m_pIPData->set_time(time);}


	/////////////////////////////////////////
	// Data
	/////////////////////////////////////////

	/// \copydoc IDataImport::constant_data()
		virtual bool constant_data() const
		{
			if(m_pIPData == NULL) return true;
			else return m_pIPData->constant_data();
		}

	///	returns the data value at ip
		const TData& operator[](size_t ip) const
		{
			UG_ASSERT(m_pIPData != NULL, "No Data set");
			UG_ASSERT(m_seriesID >= 0, "No series ticket set");
			UG_ASSERT(ip < num_ip(), "Invalid index");
			return const_cast<const IPData<TData, dim>*>(m_pIPData)->value(m_seriesID, ip);
		}

	///	returns the data value at ip
		const TData* values() const
		{
			UG_ASSERT(m_pIPData != NULL, "No Data set");
			UG_ASSERT(m_seriesID >= 0, "No series ticket set");
			return const_cast<const IPData<TData, dim>*>(m_pIPData)->values(m_seriesID);
		}

	///	return the derivative w.r.t to local function at ip
		const TData* deriv(size_t ip, size_t fct) const
		{
			UG_ASSERT(m_pIPData != NULL, "No Data set");
			UG_ASSERT(m_seriesID >= 0, "No series ticket set");
			return dynamic_cast<const DependentIPData<TData, dim>*>(m_pIPData)->deriv(m_seriesID, ip, fct);
		}

	///	return the derivative w.r.t to local function and dof at ip
		const TData& deriv(size_t ip, size_t fct, size_t dof) const
		{
			UG_ASSERT(m_pIPData != NULL, "No Data set");
			UG_ASSERT(m_seriesID >= 0, "No series ticket set");
			return dynamic_cast<const DependentIPData<TData, dim>*>(m_pIPData)->deriv(m_seriesID, ip, fct, dof);
		}

	/////////////////////////////////////////
	// Positions
	/////////////////////////////////////////

	/// number of integration points
		size_t num_ip() const
		{
			if(m_pIPData == NULL) return 0;
			else return m_pIPData->num_ip(m_seriesID);
		}

	///	set the local integration points
		template <int ldim>
		void set_local_ips(const MathVector<ldim>* vPos, size_t numIP);

	///	sets the global positions
		void set_global_ips(const MathVector<dim>* vPos, size_t numIP);

	///	position of ip
		const MathVector<dim>& position(size_t i) const
		{
			if(m_pIPData == NULL) throw(UGFatalError("No Data set"));
			return m_pIPData->ip(m_seriesID, i);
		}

	/////////////////////////////////////////
	// Linearization of Defect
	/////////////////////////////////////////

	/// number of shapes for local function
		size_t num_sh(size_t fct) const
		{
			const size_t ip = 0; check_ip_fct(ip,fct);
			return m_vvvLinDefect[ip][fct].size();
		}

	///	returns the pointer to all  linearized defects at one ip
		TData* lin_defect(size_t ip, size_t fct)
			{check_ip_fct(ip,fct);return &(m_vvvLinDefect[ip][fct][0]);}

	///	returns the pointer to all  linearized defects at one ip
		const TData* lin_defect(size_t ip, size_t fct) const
			{check_ip_fct(ip,fct);return &(m_vvvLinDefect[ip][fct][0]);}

	///	returns the linearized defect
		TData& lin_defect(size_t ip, size_t fct, size_t sh)
			{check_ip_fct_sh(ip,fct,sh);return m_vvvLinDefect[ip][fct][sh];}

	/// const access to lin defect
		const TData& lin_defect(size_t ip, size_t fct, size_t sh) const
			{check_ip_fct_sh(ip,fct,sh);return m_vvvLinDefect[ip][fct][sh];}

	///	resets all values of the linearized defect to zero
		void clear_lin_defect();

	/// compute jacobian for derivative w.r.t. non-system owned unknowns
		void assemble_jacobian(local_matrix_type& J);

	///	resize lin defect arrays
		virtual void resize(const LocalIndices& ind, const FunctionIndexMapping& map);

	protected:
	///	checks in debug mode the correct index
		inline void check_ip_fct(size_t ip, size_t fct) const;

	///	checks in debug mode the correct index
		inline void check_ip_fct_sh(size_t ip, size_t fct, size_t sh) const;

	///	series number provided by export
		int m_seriesID;

	/// connected IP Data
		IPData<TData, dim>* m_pIPData;

	/// connected export (if depended data)
		DependentIPData<TData, dim>* m_pDependentIPData;

	/// linearized defect (num_ip) x (num_fct) x (num_dofs(i))
		std::vector<std::vector<std::vector<TData> > > m_vvvLinDefect;
};

////////////////////////////////////////////////////////////////////////////////
// Data Export
////////////////////////////////////////////////////////////////////////////////


/// Base class for Data Export
/**
 * An base class for all data exports
 */
class IDataExport
{
	protected:
	///	type of local vector
		typedef IElemDisc::local_vector_type local_vector_type;

	///	type of evaluation function
		typedef bool (IElemDisc::*ExportFunc)(const local_vector_type& u, bool bDeriv);

	public:
	///	Constructor
		IDataExport() : m_id(-1), m_pObj(NULL) {}

	///	sets the geometric object type
		bool set_geometric_object_type(int id);

	///	sets the function group
		virtual void set_function_group(const FunctionGroup& fctGrp) = 0;

	///	register evaluation of linear defect for a element
		template <typename TFunc>
		void reg_export_fct(int id, IElemDisc* obj, TFunc func);

	///	get corresponding element disc
		IElemDisc* get_elem_disc() {return m_pObj;}

	///	returns if the dependent data is ready for evaluation
		virtual bool is_ready() const;

	///	virtual destructor
		virtual ~IDataExport() {}

	protected:
	///	function pointers for all elem types
		std::vector<ExportFunc>	m_vExportFunc;

	/// current Geom Object
		int m_id;

	///	elem disc
		IElemDisc* m_pObj;
};

/// Data export
/**
 * A DataExport is user data produced by an element discretization.
 */
template <typename TData, int dim>
class DataExport : 	public DependentIPData<TData, dim>,
					public IDataExport
{
	public:
	///	default constructor
		DataExport() {this->m_bCompNeedsSol = true;}

	//	implement compute() method: Not available
		virtual bool compute(bool bDeriv)
		{
			UG_LOG("ERROR in 'DataExport::compute()': Computation of Export "
				 	 "without current solution called. Cannot evaluate.\n");
			return false;
		}

	///	compute export
		virtual bool compute(const local_vector_type& u, bool bDeriv)
		{
			UG_ASSERT(m_id >=0, "ElemType id is not set correctly.");
			UG_ASSERT((size_t)m_id < m_vExportFunc.size(), "id "<<m_id<<" not registered");
			UG_ASSERT(m_vExportFunc[m_id] != NULL, "Func pointer is NULL");
			return (m_pObj->*(m_vExportFunc[m_id]))(u, bDeriv);
		}

	///	returns if data depends on solution
		virtual bool zero_derivative() const {return false;}

	///	clear dependent data
		void clear() {m_vDependData.clear();}

	///	add data dependency
		void add_needed_data(IIPData& data) {m_vDependData.push_back(&data);}

	///	remove needed data
		void remove_needed_data(IIPData& data) {remove(m_vDependData.begin(),
		                                               m_vDependData.end(),
		                                               &data);}

	///	number of other Data this data depends on
		virtual size_t num_needed_data() const {return m_vDependData.size();}

	///	return needed data
		virtual IIPData* needed_data(size_t i) {return m_vDependData.at(i);}

	///	returns if the dependent data is ready for evaluation
		virtual bool is_ready() const {return IDataExport::is_ready();}

	///	sets the function group
		virtual void set_function_group(const FunctionGroup& fctGrp)
			{return IDependentIPData::set_function_group(fctGrp);}

	protected:
		std::vector<IIPData*> m_vDependData;
};

} // end namespace ug

// include implementation
#include "data_import_export_impl.h"

#endif /* __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_IMPORT_EXPORT__ */
