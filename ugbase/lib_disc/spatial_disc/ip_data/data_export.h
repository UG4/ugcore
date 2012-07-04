/*
 * data_export.h
 *
 *  Created on: 04.07.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EXPORT__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EXPORT__

#include "data_import.h"

namespace ug{
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
class DataExport : 	public DependentUserData<TData, dim>,
					public IDataExport
{
  using DependentUserData<TData, dim>::compute;

	public:
	///	default constructor
		DataExport();

	///	implement compute() method of IUserData
		virtual void compute(LocalVector* u, GeometricObject* elem, bool bDeriv = false);

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
		void add_needed_data(SmartPtr<IUserData> data);

	///	remove needed data
		void remove_needed_data(SmartPtr<IUserData> data);

	///	number of other Data this data depends on
		virtual size_t num_needed_data() const {return m_vDependData.size();}

	///	return needed data
		virtual SmartPtr<IUserData> needed_data(size_t i) {return m_vDependData.at(i);}

	///	returns if the dependent data is ready for evaluation
		virtual void check_setup() const;

	///	sets the function group
		virtual void set_function_group(const FunctionGroup& fctGrp)
			{return IUserData::set_function_group(fctGrp);}

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
		std::vector<SmartPtr<IUserData> > m_vDependData;
};

} // end namespace ug

#include "data_export_impl.h"

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EXPORT__ */
