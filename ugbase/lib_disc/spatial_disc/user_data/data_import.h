/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DATA_IMPORT__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DATA_IMPORT__

#include "user_data.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// Data Import
////////////////////////////////////////////////////////////////////////////////

enum DiscPart {	NONE = 0,
				MASS = 1 << 0,
				STIFF = 1 << 1,
				RHS = 1 << 2,
				EXPL = 1 << 3,
				MAX_PART};

inline
std::ostream& operator<< (std::ostream& outStream, DiscPart part)
{
	switch(part)
	{
		case NONE: outStream << "(none)"; break;
		case MASS: outStream << "mass"; break;
		case STIFF: outStream << "stiff"; break;
		case RHS: outStream << "rhs"; break;
		case EXPL: outStream << "expl"; break;
		default: UG_THROW("Unknown DiscPart in operator<<");
	}
	return outStream;
};

/// Base class for data import
/**
 * An IDataImport is the base class for importing data to ElemDiscs
 */
template <int dim>
class IDataImport
{
	public:
	/// Constructor
		IDataImport(bool compLinDefect = true)
			: m_spICplUserData(nullptr), m_part(STIFF),
			 m_bCompLinDefect(compLinDefect)
		{}

		virtual ~IDataImport()	{}

	///	sets part of import
		void set_part(DiscPart part) {m_part = part;}

	///	set to no part
		void set_none_part() {m_part = NONE;}

	///	set to explicit part
		void set_expl_part() {m_part = EXPL;}

	///	sets if import is located in mass part (for time dependent problems)
		void set_mass_part() {m_part = MASS;}

	///	sets if import is located in rhs part
		void set_rhs_part() {m_part = RHS;}

	///	sets if import is located in stiff part (default)
		void set_stiff_part() {m_part = STIFF;}

	///	sets if lin defect is to be computed
		void set_comp_lin_defect(bool b) {m_bCompLinDefect = b;}

	///	returns if import is located in mass part (for time dependent problems)
		DiscPart part() const {return m_part;}

	/// returns if data is set
		virtual bool data_given() const = 0;

	/// returns if data is constant
	/**
	 * This method, returns if the connected data is constant.
	 */
		virtual bool constant() const = 0;

	///	returns if data depends on unknown functions
	/**
	 * This method returns if the data depends on the unknown functions.
	 */
		bool zero_derivative() const
		{
			if(!m_spICplUserData.valid()) return true;
			else if (m_spICplUserData->zero_derivative()) return true;
			else return !m_bCompLinDefect;
		}

	/// returns the connected user data
		virtual SmartPtr<ICplUserData<dim> > data() = 0;

	///	set function group for linearization of defect
		void set_map(const FunctionIndexMapping& map){m_map = map;}

	///	get function mapping
		const FunctionIndexMapping& map() const{return m_map;}

	///	get function mapping of dependent data
		const FunctionIndexMapping& conn_map() const{return m_spICplUserData->map();}

	/// number of functions
		size_t num_fct() const {return m_map.num_fct();}

	///	sets the geometric object type
		virtual void set_roid(ReferenceObjectID id) = 0;

	///	checks if ready for evaluation
		virtual void check_setup() = 0;

	///	compute lin defect
		virtual void compute_lin_defect(LocalVector& u) = 0;

	///	resize arrays
		virtual void update_dof_sizes(const LocalIndices& ind) = 0;

	///	add jacobian entries introduced by this import
		virtual void add_jacobian(LocalMatrix& J, const number scale) = 0;

	///	removes the positions
		virtual void clear_ips() = 0;

	protected:
	/// connected iexport
		SmartPtr<ICplUserData<dim> > m_spICplUserData;

	///	Mapping for import fct
		FunctionIndexMapping m_map;

	protected:
	///	flag to indicate where import is located
		DiscPart m_part;

	///	indicates iff lin defect should be computed
		bool m_bCompLinDefect;
};

/// Data import
/**
 * A DataImport is used to import data into an ElemDisc.
 */
template <typename TData, int dim>
class DataImport : public IDataImport<dim>
{
	public:
	/// Constructor
		DataImport(bool bLinDefect = true) : IDataImport<dim>(bLinDefect),
			m_id(ROID_UNKNOWN),
			m_seriesID(-1),	m_spUserData(nullptr), m_vValue(nullptr),
			m_numIP(0), m_spDependentUserData(nullptr)
		{clear_fct();
		}

	///	Destructor
		~DataImport();

	///	set the user data
		void set_data(SmartPtr<CplUserData<TData, dim> > spData);

	/// returns the connected ICplUserData
		SmartPtr<ICplUserData<dim> > data() {return m_spUserData.template cast_dynamic<ICplUserData<dim> >();}

	/// returns the connected ICplUserData
		SmartPtr<CplUserData<TData, dim> > user_data(){return m_spUserData;}

	///	returns true if data given
		virtual bool data_given() const {return m_spUserData.valid();}


	/////////////////////////////////////////
	// Data
	/////////////////////////////////////////

	/// \copydoc IDataImport::constant()
		virtual bool constant() const
		{
			if(data_given()) return m_spUserData->constant();
			else return false;
		}

	///	returns the data value at ip
		const TData& operator[](size_t ip) const{check_ip(ip); return m_vValue[ip];}

	///	returns the data value at ip
		const TData* values() const {check_values(); return m_vValue;}

	///	return the derivative w.r.t to local function at ip
		const TData* deriv(size_t ip, size_t fct) const
		{
			UG_ASSERT(m_spDependentUserData.valid(), "No Dependent Data set");
			UG_ASSERT(m_seriesID >= 0, "No series ticket set");
			return m_spDependentUserData->deriv(m_seriesID, ip, fct);
		}

	///	return the derivative w.r.t to local function and dof at ip
		const TData& deriv(size_t ip, size_t fct, size_t dof) const
		{
			UG_ASSERT(m_spDependentUserData.valid(), "No Dependent Data set");
			UG_ASSERT(m_seriesID >= 0, "No series ticket set");
			return m_spDependentUserData->deriv(m_seriesID, ip, fct, dof);
		}

	/////////////////////////////////////////
	// Positions
	/////////////////////////////////////////

	/// number of integration points
		size_t num_ip() const {return m_numIP;}

	///	set the local integration points
		template <int ldim>
		void set_local_ips(const MathVector<ldim>* vPos, size_t numIP, int timePointSpec,
		                   bool bMayChange = true);

	///	set the local integration points
		void set_local_ips(const MathVector<dim>* vPos, size_t numIP, int timePointSpec,
		                   bool bMayChange = true);

	///	set the local integration points
		template <int ldim>
		void set_local_ips(const MathVector<ldim>* vPos, size_t numIP,
		                   bool bMayChange = true);

	///	set the local integration points
		void set_local_ips(const MathVector<dim>* vPos, size_t numIP,
		                   bool bMayChange = true);

	///	set the time point specification
		void set_time_point(int timePointSpec);

	///	sets the global positions
		void set_global_ips(const MathVector<dim>* vPos, size_t numIP);

	///	position of ip
		const MathVector<dim>& position(size_t i) const
		{
			if(data_given()) return m_spUserData->ip(m_seriesID, i);
			 UG_THROW("DataImport::position: "
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
			UG_ASSERT(fct <  m_vvNumDoFPerFct[fct], "Invalid index");
			return m_vvNumDoFPerFct[fct];
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
		void add_jacobian(LocalMatrix& J, const number scale);

	///	resize lin defect arrays
		virtual void update_dof_sizes(const LocalIndices& ind);

	public:
	///	type of evaluation function
		using LinDefectFunc = boost::function<void (const LocalVector& u,
	                                            std::vector<std::vector<TData> > vvvLinDefect[],
	                                            const size_t nip)>;

	///	sets the geometric object type
		virtual void set_roid(ReferenceObjectID id);

	///	checks if ready for evaluation
		virtual void check_setup();

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
		virtual void compute_lin_defect(LocalVector& u)
		{
		///	compute the linearization only if the export parameter is 'at current time'
			if (! m_spUserData->at_current_time (m_seriesID))
				return;
		///	compute the linearization
			UG_ASSERT(m_vLinDefectFunc[m_id] != nullptr, "No evaluation function.");
			UG_ASSERT(num_ip() == 0 || m_vvvLinDefect.size() >= num_ip(),
			          "DataImport: Num ip "<<num_ip()<<", but memory: "<<m_vvvLinDefect.size());
			u.access_by_map(this->map());
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

	///	caches data access
		void cache_data_access();

	///	resizes the lin defect arrays for current number of ips.
		void resize_defect_array();

	/// current Geom Object
		ReferenceObjectID m_id;

	///	function pointers for all elem types
		LinDefectFunc m_vLinDefectFunc[NUM_REFERENCE_OBJECTS];

	///	series number provided by export
		int m_seriesID;

	/// connected UserData
		SmartPtr<CplUserData<TData, dim> > m_spUserData;

	///	cached access to the UserData field
		const TData* m_vValue;

	///	number of ips
		size_t m_numIP;

	/// connected export (if depended data)
		SmartPtr<DependentUserData<TData, dim> > m_spDependentUserData;

	///	number of functions and their dofs
		std::vector<size_t> m_vvNumDoFPerFct;

	/// linearized defect (num_ip) x (num_fct) x (num_dofs(i))
		std::vector<std::vector<std::vector<TData> > > m_vvvLinDefect;
};

} // end namespace ug

// include implementation
#include "data_import_impl.h"

#endif