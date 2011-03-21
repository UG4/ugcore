/*
 * data_import.h
 *
 *  Created on: 17.12.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_IMPORT__
#define __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_IMPORT__

#include "user_data.h"
#include "../elem_disc/elem_disc_interface.h"
#include "data_export.h"

namespace ug{

// predeclaration
template <typename TAlgebra> class IElemDisc;

/// Base class for data import
/**
 * An IDataImport is the base class for importing data to ElemDiscs
 */
template <typename TAlgebra>
class IDataImport
{
	typedef typename IElemDisc<TAlgebra>::local_matrix_type local_matrix_type;

	public:
	/// Constructor
		IDataImport(bool compLinDefect = true)
			: m_pIDataExport(NULL), m_pIDependentIPData(NULL),
			  m_id(-1), m_bCompLinDefect(compLinDefect)
		{}
		virtual ~IDataImport()	{}

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
		void set_function_group(const FunctionGroup& fctGrp){m_pFctGrp = &fctGrp;}

	///	get funtion group
		const FunctionGroup& get_function_group() const
		{
			UG_ASSERT(m_pFctGrp != NULL, "No func group set.");
			return *m_pFctGrp;
		}

	/// number of functions
		size_t num_fct() const
		{
			UG_ASSERT(m_pFctGrp != NULL, "No func group set.");
			return m_pFctGrp->num_fct();
		}

	///	resize arrays
		virtual void resize(const LocalIndices& ind) = 0;

	///	add jacobian entries introduced by this import
		virtual void assemble_jacobian(local_matrix_type& J) = 0;

	protected:
	/// connected iexport
		IDataExport<TAlgebra>* m_pIDataExport;
		IDependentIPData* m_pIDependentIPData;

	///	function group for linear defect
		const FunctionGroup* m_pFctGrp;

	protected:
	///	type of local vector
		typedef typename IElemDisc<TAlgebra>::local_vector_type local_vector_type;

	///	type of evaluation function
		typedef bool (IElemDisc<TAlgebra>::*LinDefectFunc)(const local_vector_type& u);

	///	function pointers for all elem types
		std::vector<LinDefectFunc>	m_vLinDefectFunc;

	/// current Geom Object
		int m_id;

	///	elem disc
		IElemDisc<TAlgebra>* m_pObj;

	///	indicates iff lin defect should be computed
		bool m_bCompLinDefect;

	public:
	///	sets the geometric object type
		bool set_geometric_object_type(int id)
		{
		//	if lin defect is not supposed to be computed, we're done
			if(!m_bCompLinDefect) return true;

		//	Check for evaluation function and choose it if present
			if(id < (int)m_vLinDefectFunc.size() && m_vLinDefectFunc[id] != NULL)
			{
				m_id = id;
				return true;
			}
		//	return error else
			else
			{
				UG_LOG("No or not all lin defect functions registered "
						"for object with reference object id " << id << ".\n");
				m_id = -1; return false;
			}
		}

	///	register evaluation of linear defect for a element
		template <typename TFunc>
		void register_lin_defect_func(int id, IElemDisc<TAlgebra>* obj, TFunc func)
		{
		//	make sure that there is enough space
			if((size_t)id >= m_vLinDefectFunc.size())
				m_vLinDefectFunc.resize(id+1, NULL);

			m_vLinDefectFunc[id] = (LinDefectFunc)func;
			m_pObj = obj;
		}

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
template <typename TData, int dim, typename TAlgebra>
class DataImport : public IDataImport<TAlgebra>
{
	typedef typename IElemDisc<TAlgebra>::local_matrix_type local_matrix_type;
	using IDataImport<TAlgebra>::num_fct;

	public:
	/// Constructor
		DataImport(bool compLinDefect = true) :
			IDataImport<TAlgebra>(compLinDefect),
			m_seriesID(-1), m_pIPData(NULL), m_pDependentIPData(NULL)
		{}

	///	set the user data
		void set_data(IPData<TData, dim>& data)
		{
		//	remember IPData
			m_pIPData = &data;

		//	remember iexport
			this->m_pIDataExport = dynamic_cast<IDataExport<TAlgebra>*>(&data);
			this->m_pIDependentIPData = dynamic_cast<IDependentIPData*>(&data);

		//	remember dependent data (i.e. is NULL iff no dependent data given)
			m_pDependentIPData = dynamic_cast<DependentIPData<TData, dim>*>(&data);
		}

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
			UG_ASSERT((dynamic_cast<const DependentIPData<TData, dim>*>(m_pIPData)) != NULL, "No Data set");
			UG_ASSERT(m_seriesID >= 0, "No series ticket set");
			return dynamic_cast<const DependentIPData<TData, dim>*>(m_pIPData)->deriv(m_seriesID, ip, fct);
		}

	///	return the derivative w.r.t to local function and dof at ip
		const TData& deriv(size_t ip, size_t fct, size_t dof) const
		{
			UG_ASSERT((dynamic_cast<const DependentIPData<TData, dim>*>(m_pIPData)) != NULL, "No Data set");
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
			else
			{
				UG_ASSERT(m_seriesID >= 0, "No series ticket set");
				return m_pIPData->num_ip(m_seriesID);
			}
		}

	///	set the local integration points
		template <int ldim>
		void set_local_ips(const MathVector<ldim>* vPos, size_t numIP)
		{
		//	if no data set, skip
			if(m_pIPData == NULL) return;

		//	request series
			m_seriesID = m_pIPData->template
						register_local_ip_series<ldim>(vPos,numIP);
		}

	///	sets the global positions
		void set_global_ips(const MathVector<dim>* vPos, size_t numIP)
		{
		//  if no data set, skip
			if(m_pIPData == NULL) return;

		//	set global ips for series ID
			UG_ASSERT(m_seriesID >= 0, "Wrong series id.");
			m_pIPData->set_global_ips(m_seriesID,vPos,numIP);
		}

	///	position of ip
		const MathVector<dim>& position(size_t i) const
		{
			if(m_pIPData == NULL)
				throw(UGFatalError("No Data set"));

			return m_pIPData->ip(m_seriesID, i);
		}

	/////////////////////////////////////////
	// Linearization of Defect
	/////////////////////////////////////////

	/// number of shapes for local function
		size_t num_sh(size_t fct) const
		{
			const size_t ip = 0;
			UG_ASSERT(ip < m_vvvLinDefect.size(), "Invalid index.");
			UG_ASSERT(fct < m_vvvLinDefect[ip].size(), "Invalid index.");
			return m_vvvLinDefect[ip][fct].size();
		}

	///	returns the linearized defect
		TData& lin_defect(size_t ip, size_t fct, size_t dof)
		{
			UG_ASSERT(ip  < m_vvvLinDefect.size(), "Invalid index.");
			UG_ASSERT(fct < m_vvvLinDefect[ip].size(), "Invalid index.");
			UG_ASSERT(dof < m_vvvLinDefect[ip][fct].size(), "Invalid index.");
			return m_vvvLinDefect[ip][fct][dof];
		}

	/// const access to lin defect
		const TData& lin_defect(size_t ip, size_t fct, size_t dof) const
		{
			UG_ASSERT(ip  < m_vvvLinDefect.size(), "Invalid index.");
			UG_ASSERT(fct < m_vvvLinDefect[ip].size(), "Invalid index.");
			UG_ASSERT(dof < m_vvvLinDefect[ip][fct].size(), "Invalid index.");
			return m_vvvLinDefect[ip][fct][dof];
		}

	/// compute jacobian for derivative w.r.t. non-system owned unknowns
		void assemble_jacobian(local_matrix_type& J)
		{
			UG_ASSERT(m_pDependentIPData != NULL, "No Export set.");

			for(size_t fct1 = 0; fct1 < num_fct(); ++fct1)
				for(size_t fct2 = 0; fct2 < m_pDependentIPData->num_fct(); ++fct2)
					for(size_t dof1 = 0; dof1 < num_sh(fct1); ++dof1)
						for(size_t dof2 = 0; dof2 < m_pDependentIPData->num_sh(m_seriesID, fct2); ++dof2)
							for(size_t ip = 0; ip < num_ip(); ++ip)
							{
								number prod = lin_defect(ip, fct1, dof1)
												* m_pDependentIPData->deriv(m_seriesID, ip, fct2, dof2);
								J(fct1, dof1, fct2, dof2) += prod;
							}
		}

	///	resize lin defect arrays
		virtual void resize(const LocalIndices& ind)
		{
		//	resize ips
			m_vvvLinDefect.resize(num_ip());

		//	resize num fct
			for(size_t ip = 0; ip < num_ip(); ++ip)
			{
			//	resize num fct
				m_vvvLinDefect[ip].resize(ind.num_fct());

			//	resize dofs
				for(size_t fct = 0; fct < ind.num_fct(); ++fct)
					m_vvvLinDefect[ip][fct].resize(ind.num_dofs(fct));
			}
		}

	protected:
	///	series number provided by export
		int m_seriesID;

	/// connected IP Data
		IPData<TData, dim>* m_pIPData;

	/// connected export (if depended data)
		DependentIPData<TData, dim>* m_pDependentIPData;

	/// linearized defect (num_ip) x (num_fct) x (num_dofs(i))
		std::vector<std::vector<std::vector<TData> > > m_vvvLinDefect;
};

} // end namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_IMPORT__ */
