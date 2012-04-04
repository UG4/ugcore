/*
 * data_import_export.h
 *
 *  Created on: 17.12.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DATA_IMPORT_EXPORT__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DATA_IMPORT_EXPORT__

#include "ip_data.h"
//#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"

namespace ug{

// predeclaration
class IElemDisc;

////////////////////////////////////////////////////////////////////////////////
// Data Import
////////////////////////////////////////////////////////////////////////////////

/// Base class for data import
/**
 * An IDataImport is the base class for importing data to ElemDiscs
 */
class IDataImport
{
	public:
	/// Constructor
		IDataImport(bool compLinDefect = true)
			: m_spIDependentIPData(NULL),
			 m_bInMassPart(false), m_bInRhsPart(false),
			 m_bCompLinDefect(compLinDefect)
		{}

		virtual ~IDataImport()	{}

	///	sets if import is located in mass part (for time dependent problems)
		void set_mass_part(bool bInMassPart) {m_bInMassPart = bInMassPart;}

	///	returns if import is located in mass part (for time dependent problems)
		bool in_mass_part() const {return m_bInMassPart;}

	///	sets if import is located in rhs part
		void set_rhs_part(bool bInRhsPart) {m_bInRhsPart = bInRhsPart;}

	///	returns if import is located in rhs part
		bool in_rhs_part() const {return m_bInRhsPart;}

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
			if(!m_spIDependentIPData.valid()) return true;
			else return !m_bCompLinDefect;
		}

	/// returns the connected ip data
		virtual SmartPtr<IIPData> data() = 0;

	///	set function group for linearization of defect
		void set_function_group(const FunctionGroup& fctGrp){m_fctGrp = fctGrp;}

	///	get funtion group
		const FunctionGroup& function_group() const{return m_fctGrp;}

	/// number of functions
		size_t num_fct() const {return m_fctGrp.num_fct();}

	///	sets the geometric object type
		virtual bool set_roid(ReferenceObjectID id) = 0;

	///	compute lin defect
		virtual void compute_lin_defect(const LocalVector& u) = 0;

	///	resize arrays
		virtual void resize(const LocalIndices& ind,
		                    const FunctionIndexMapping& map) = 0;

	///	add jacobian entries introduced by this import
		virtual void assemble_jacobian(LocalMatrix& J) = 0;

	///	removes the positions
		virtual void clear_ips() = 0;

	protected:
	/// connected iexport
		SmartPtr<IDependentIPData> m_spIDependentIPData;

	///	function group for linear defect
		FunctionGroup m_fctGrp;

	protected:
	///	flag to indicate if import is located in mass part
		bool m_bInMassPart;

	///	flag to indicate if import is located in rhs part
		bool m_bInRhsPart;

	///	indicates iff lin defect should be computed
		bool m_bCompLinDefect;
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
	public:
	/// Constructor
		DataImport(bool bLinDefect = true) : IDataImport(bLinDefect),
			m_id(ROID_UNKNOWN),
			m_seriesID(-1),	m_spIPData(NULL), m_vValue(NULL),
			m_numIP(0), m_spDependentIPData(NULL)
		{clear_fct();}

	///	set the user data
		void set_data(SmartPtr<IPData<TData, dim> > spData);

	/// returns the connected IIPData
		SmartPtr<IIPData> data() {return m_spIPData.template cast_dynamic<IIPData>();}

	/// returns the connected IIPData
		SmartPtr<IPData<TData, dim> > ip_data(){return m_spIPData;}

	///	returns true if data given
		virtual bool data_given() const {return m_spIPData.valid();}

	///	sets the evaluation time point
		void set_time(number time) {if(data_given()) m_spIPData->set_time(time);}


	/////////////////////////////////////////
	// Data
	/////////////////////////////////////////

	/// \copydoc IDataImport::constant_data()
		virtual bool constant_data() const
		{
			if(data_given()) return m_spIPData->constant_data();
			else return false;
		}

	///	returns the data value at ip
		const TData& operator[](size_t ip) const{check_ip(ip); return m_vValue[ip];}

	///	returns the data value at ip
		const TData* values() const {check_values(); return m_vValue;}

	///	return the derivative w.r.t to local function at ip
		const TData* deriv(size_t ip, size_t fct) const
		{
			UG_ASSERT(m_spDependentIPData.valid(), "No Dependent Data set");
			UG_ASSERT(m_seriesID >= 0, "No series ticket set");
			return m_spDependentIPData->deriv(m_seriesID, ip, fct);
		}

	///	return the derivative w.r.t to local function and dof at ip
		const TData& deriv(size_t ip, size_t fct, size_t dof) const
		{
			UG_ASSERT(m_spDependentIPData.valid(), "No Dependent Data set");
			UG_ASSERT(m_seriesID >= 0, "No series ticket set");
			return m_spDependentIPData->deriv(m_seriesID, ip, fct, dof);
		}

	/////////////////////////////////////////
	// Positions
	/////////////////////////////////////////

	/// number of integration points
		size_t num_ip() const {return m_numIP;}

	///	set the local integration points
		template <int ldim>
		void set_local_ips(const MathVector<ldim>* vPos, size_t numIP,
		                   bool bMayChange = true);

	///	set the local integration points
		void set_local_ips(const MathVector<dim>* vPos, size_t numIP,
		                   bool bMayChange = true);

	///	sets the global positions
		void set_global_ips(const MathVector<dim>* vPos, size_t numIP);

	///	position of ip
		const MathVector<dim>& position(size_t i) const
		{
			if(data_given()) return m_spIPData->ip(m_seriesID, i);
			 UG_THROW_FATAL("DataLinker::position: "
					 	 	 "No Data set, but positions requested.");
		}

	///	removes the positions
		virtual void clear_ips();

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

	/// compute jacobian for derivative w.r.t. non-system owned unknowns
		void assemble_jacobian(LocalMatrix& J);

	///	resize lin defect arrays
		virtual void resize(const LocalIndices& ind, const FunctionIndexMapping& map);

	public:
	///	type of evaluation function
		typedef boost::function<void (const LocalVector& u,
		                              std::vector<std::vector<TData> > vvvLinDefect[],
		                              const size_t nip)> LinDefectFunc;

	///	sets the geometric object type
		virtual bool set_roid(ReferenceObjectID id);

	///	register evaluation of linear defect for a element
		template <typename TClass>
		void set_fct(ReferenceObjectID id, TClass* obj,
		             void (TClass::*func)(
		            		 const LocalVector& u,
		            		 std::vector<std::vector<TData> > vvvLinDefect[],
		            		 const size_t nip));

	///	register evaluation of linear defect for a element
		void set_fct(ReferenceObjectID id,
					 void (*func)(
							 const LocalVector& u,
							 std::vector<std::vector<TData> > vvvLinDefect[],
							 const size_t nip));

	///	clear all evaluation functions
		void clear_fct();

	///	compute lin defect
		virtual void compute_lin_defect(const LocalVector& u)
		{
			UG_ASSERT(m_vLinDefectFunc[m_id] != NULL, "No evaluation function.");
			(m_vLinDefectFunc[m_id])(u, &m_vvvLinDefect[0], m_numIP);
		}

	protected:
	///	checks in debug mode the correct index
		inline void check_ip_fct(size_t ip, size_t fct) const;

	///	checks in debug mode the correct index
		inline void check_ip_fct_sh(size_t ip, size_t fct, size_t sh) const;

	///	checks in debug mode the correct index
		inline void check_ip(size_t ip) const;

	///	checks in debug mode the correct index
		inline void check_values() const;

	///	resizes the lin defect arrays for current number of ips.
		void resize_defect_array();

	/// current Geom Object
		ReferenceObjectID m_id;

	///	function pointers for all elem types
		LinDefectFunc m_vLinDefectFunc[NUM_REFERENCE_OBJECTS];

	///	series number provided by export
		int m_seriesID;

	/// connected IP Data
		SmartPtr<IPData<TData, dim> > m_spIPData;

	///	cached access to the IPData field
		const TData* m_vValue;

	///	number of ips
		size_t m_numIP;

	/// connected export (if depended data)
		SmartPtr<DependentIPData<TData, dim> > m_spDependentIPData;

	///	number of functions and their dofs
		std::vector<size_t> m_vvNumDoFPerFct;

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
	public:
	///	Constructor
		IDataExport() {}

	///	sets the geometric object type
		virtual bool set_roid(ReferenceObjectID id) = 0;

	///	sets the function group
		virtual void set_function_group(const FunctionGroup& fctGrp) = 0;

	///	virtual destructor
		virtual ~IDataExport() {}
};

/// Data export
/**
 * A DataExport is user data produced by an element discretization.
 */
template <typename TData, int dim>
class DataExport : 	public DependentIPData<TData, dim>,
					public IDataExport
{
  using DependentIPData<TData, dim>::compute;

	public:
	///	default constructor
		DataExport();

	//	implement compute() method of IIPData: Not available
		virtual void compute(bool bDeriv = false);

	///	compute export (implements IDependendIPData::compute)
		virtual void compute_with_sol(const LocalVector& u, bool bDeriv);

	///	sets the geometric object type
		virtual bool set_roid(ReferenceObjectID id);

	///	register evaluation of export function
		template <typename T, int refDim>
		void set_fct(ReferenceObjectID id, IElemDisc* obj,
		             void (T::*func)(const LocalVector& u,
		            		 	 	 const MathVector<dim> vGlobIP[],
		            		 	 	 const MathVector<refDim> vLocIP[],
		            		 	 	 const size_t nip,
		            		 	 	 TData vValue[],
		            		 	 	 bool bDeriv,
		            		 	 	 std::vector<std::vector<TData> > vvvDeriv[]));

	///	clears all export functions
		void clear_fct();

	///	returns if data depends on solution
		virtual bool zero_derivative() const {return false;}

	///	clear dependent data
		void clear() {m_vDependData.clear();}

	///	add data dependency
		void add_needed_data(SmartPtr<IIPData> data);

	///	remove needed data
		void remove_needed_data(SmartPtr<IIPData> data);

	///	number of other Data this data depends on
		virtual size_t num_needed_data() const {return m_vDependData.size();}

	///	return needed data
		virtual SmartPtr<IIPData> needed_data(size_t i) {return m_vDependData.at(i);}

	///	returns if the dependent data is ready for evaluation
		virtual bool is_ready() const;

	///	sets the function group
		virtual void set_function_group(const FunctionGroup& fctGrp)
			{return IDependentIPData::set_function_group(fctGrp);}

	protected:
		template <typename T, int refDim>
		inline void comp(const LocalVector& u, bool bDeriv);

	/// current Geom Object
		ReferenceObjectID m_id;

	///	corresponding elem disc
		IElemDisc* m_pObj;

	///	function pointers for all elem types
		typedef void (IElemDisc::*DummyMethod)();
		DummyMethod	m_vExportFunc[NUM_REFERENCE_OBJECTS];

		typedef void (DataExport::*CompFct)(const LocalVector&, bool);
		CompFct m_vCompFct[NUM_REFERENCE_OBJECTS];

	///	data the export depends on
		std::vector<SmartPtr<IIPData> > m_vDependData;
};

} // end namespace ug

// include implementation
#include "data_import_export_impl.h"

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DATA_IMPORT_EXPORT__ */
