/*
 * data_export.h
 *
 *  Created on: 17.12.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_EXPORT__
#define __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_EXPORT__

#include "ip_data.h"
#include "../elem_disc/elem_disc_interface.h"

namespace ug{

/// Base class for Data Export
/**
 * An base class for all data exports
 */
template <typename TAlgebra>
class IDataExport : virtual public IDependentIPData
{
	protected:
	///	type of local vector
		typedef typename IElemDisc<TAlgebra>::local_vector_type local_vector_type;

	///	type of evaluation function
		typedef bool (IElemDisc<TAlgebra>::*ExportFunc)(const local_vector_type& u,
														bool compDeriv);

	public:
	///	Constructor
		IDataExport() : m_id(-1), m_pObj(NULL) {}
		~IDataExport()	{}

	/// clear ips
		virtual void clear_export_ips() = 0;

	///	sets the geometric object type
		bool set_geometric_object_type(int id)
		{
			if(id < 0 || (size_t)id >= m_vExportFunc.size() || m_vExportFunc[id] == NULL)
			{
				UG_LOG("No or not all functions registered "
						"for object with reference object id " << id << ".\n");
				m_id = -1; return false;
			}
			else{m_id = id;	return true;}
		}

	///	register evaluation of linear defect for a element
		template <typename TFunc>
		void register_export_func(int id, IElemDisc<TAlgebra>* obj, TFunc func)
		{
		//	make sure that there is enough space
			if((size_t)id >= m_vExportFunc.size())
				m_vExportFunc.resize(id+1, NULL);

			m_vExportFunc[id] = (ExportFunc)func;

			if(m_pObj == NULL) m_pObj = obj;
			else if(m_pObj != obj)
				throw(UGFatalError("Exports assume to be used by on object for all functions."));
		}

	///	compute export
		bool compute_export(const local_vector_type& u, bool compDeriv)
		{
			UG_ASSERT(m_id >=0, "ElemType id is not set correctly.");
			UG_ASSERT((size_t)m_id < m_vExportFunc.size(), "id "<<m_id<<" not registered");
			UG_ASSERT(m_vExportFunc[m_id] != NULL, "Func pointer is NULL");
			return (m_pObj->*(m_vExportFunc[m_id]))(u, compDeriv);
		}

	protected:
	///	function pointers for all elem types
		std::vector<ExportFunc>	m_vExportFunc;

	/// current Geom Object
		int m_id;

	///	elem disc
		IElemDisc<TAlgebra>* m_pObj;
};

/// Data export
/**
 * A DataExport is user data produced by an element discretization.
 */
template <typename TData, int dim, typename TAlgebra>
class DataExport : 	public DependentIPData<TData, dim>,
					public IDataExport<TAlgebra>
{
	public:
		virtual void compute(bool computeDeriv)
			{throw(UGFatalError("Should not be called."));}

	///	returns if data depends on solution
		virtual bool zero_derivative() const {return false;}

	/// clear ips
		virtual void clear_export_ips()
		{
			this->clear_ips();
		}

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

	protected:
		std::vector<IIPData*> m_vDependData;
};

} // end namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_EXPORT__ */
