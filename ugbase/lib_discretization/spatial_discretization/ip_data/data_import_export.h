/*
 * data_import_export.h
 *
 *  Created on: 17.12.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_IMPORT_EXPORT__
#define __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_IMPORT_EXPORT__

#include "user_data.h"
#include "../elem_disc/elem_disc_interface.h"

namespace ug{

/// Base class for Data Export
/**
 * An base class for all data exports
 */
template <typename TAlgebra>
class IDataExport
{
	public:
	///	Constructor
		IDataExport() : m_id(-1), m_pObj(NULL) {}

	/// clear ips
		virtual void clear_ips_virt() = 0;

	/// set	Function Group of functions
		void set_function_group(const FunctionGroup& fctGrp)
			{m_pFctGrp = &fctGrp;}

	///	Function Group of functions
		const FunctionGroup& get_function_group() const
		{
			UG_ASSERT(m_pFctGrp != NULL, "No func group set.");
			return *m_pFctGrp;
		}

	///	number of fuctions this export depends on
		size_t num_fct() const
		{
			UG_ASSERT(m_pFctGrp != NULL, "No func group set.");
			return m_pFctGrp->num_fct();
		}

	///	resize arrays
		virtual void resize(const LocalIndices& ind) = 0;

	protected:
	/// number of functions the data depends on
		const FunctionGroup* m_pFctGrp;

	protected:
	///	type of local vector
		typedef typename IElemDisc<TAlgebra>::local_vector_type local_vector_type;

	///	type of evaluation function
		typedef bool (IElemDisc<TAlgebra>::*ExportFunc)(const local_vector_type& u);

	///	function pointers for all elem types
		std::vector<ExportFunc>	m_vExportFunc;

	/// current Geom Object
		int m_id;

	///	elem disc
		IElemDisc<TAlgebra>* m_pObj;

	public:
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

	///	compute lin defect
		bool compute_export(const local_vector_type& u)
		{
			UG_ASSERT(m_id >=0, "ElemType id is not set correctly.");
			UG_ASSERT((size_t)m_id < m_vExportFunc.size(), "id "<<m_id<<" not registered");
			UG_ASSERT(m_vExportFunc[m_id] != NULL, "Func pointer is NULL");
			return (m_pObj->*(m_vExportFunc[m_id]))(u);
		}
};

/// Data export
/**
 * A DataExport is user data produced by an element discretization.
 */
template <typename TData, int dim, typename TAlgebra>
class DataExport : 	public IPData<TData, dim>,
					public IDataExport<TAlgebra>
{
	public:
	using IPData<TData, dim>::num_series;
	using IPData<TData, dim>::num_ip;
	using IPData<TData, dim>::local_ips;

	public:
	/// Constructor
		DataExport()
		{}

	///	clear ips
		virtual void clear_ips_virt(){this->clear_ips();}

	/// compute values (and derivatives iff compDeriv == true)
		virtual void compute(bool compDeriv = false)
			{throw(UGFatalError("Not Implemented."));}

	////////////////////////////////
	// Derivatives
	////////////////////////////////

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
		IDataImport() : m_pIDataExport(NULL), m_id(-1) {}

	/// returns if data is set
		virtual bool data_set() const = 0;

	/// returns if data is constant
	/**
	 * This method, returns if the connected data is constant.
	 */
		virtual bool constant_data() const = 0;

	///	returns if data depends on unknown functions
	/**
	 * This method returns if the data depends on the unknown functions.
	 */
		bool non_zero_derivative() const
		{
			if(m_pIDataExport == NULL) return false;
			else return true;
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

	public:
	///	sets the geometric object type
		bool set_geometric_object_type(int id)
		{
			if(m_vLinDefectFunc[id] != NULL){m_id = id;	return true;}
			else
			{
				UG_LOG("No or not all functions registered "
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
		DataImport() :
			m_seriesID(-1), m_pIPData(NULL), m_pDataExport(NULL)
		{}

	///	set the user data
		void set_data(IPData<TData, dim>& data)
		{
		//	remember IPData
			m_pIPData = &data;

		//	remember iexport
			this->m_pIDataExport = dynamic_cast<IDataExport<TAlgebra>*>(&data);

		//	remember export (i.e. is NULL iff no export given)
			m_pDataExport = dynamic_cast<DataExport<TData, dim, TAlgebra>*>(&data);
		}

	/// returns the connected IIPData
		IIPData* get_data()
		{
			return m_pIPData;
		}

	/////////////////////////////////////////
	// Data
	/////////////////////////////////////////

	///	returns if zero data given
		virtual bool data_set() const {return !zero_data();}

	///	returns if zero data given
		bool zero_data() const
			{return m_pIPData == NULL;}

	/// \copydoc IDataImport::non_constant_data()
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

	/// returns true iff depends on global position
		bool depends_on_global_position() const
		{
			if(m_pDataExport == NULL)
				return false;
			else
				return true;
		}

	/// returns true iff depends on global position
		bool depends_on_local_position() const
		{
			if(m_pDataExport == NULL)
				return false; // only Exports depend on local positions
			else
				return true;
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

	///	returns if import depends on other solutions
		bool non_zero_derivative() const
		{
			return m_pDataExport != NULL;
		}

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
			UG_ASSERT(m_pDataExport != NULL, "No Export set.");

			for(size_t fct1 = 0; fct1 < num_fct(); ++fct1)
				for(size_t fct2 = 0; fct2 < m_pDataExport->num_fct(); ++fct2)
					for(size_t dof1 = 0; dof1 < num_sh(fct1); ++dof1)
						for(size_t dof2 = 0; dof2 < m_pDataExport->num_sh(m_seriesID, fct2); ++dof2)
							for(size_t ip = 0; ip < num_ip(); ++ip)
							{
								number prod = lin_defect(ip, fct1, dof1)
												* m_pDataExport->deriv(m_seriesID, ip, fct2, dof2);
/*								if(prod < -1e-10 ||prod > 1e-10)
								{
									UG_LOG("Adding "<<prod<<"for lin="<<lin_defect(ip, fct1, dof1)<<
									       "and deriv="<<m_pDataExport->deriv(m_seriesID, ip, fct2, dof2)<<"\n");
								}
*/								J(fct1, dof1, fct2, dof2) += prod;
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
		DataExport<TData, dim, TAlgebra>* m_pDataExport;

	/// linearized defect (num_ip) x (num_fct) x (num_dofs(i))
		std::vector<std::vector<std::vector<TData> > > m_vvvLinDefect;
};


template <typename TAlgebra>
class DataEvaluator
{
	typedef typename IElemDisc<TAlgebra>::local_vector_type local_vector_type;
	typedef typename IElemDisc<TAlgebra>::local_matrix_type local_matrix_type;

	public:
		DataEvaluator(const std::vector<IElemDisc<TAlgebra>*>& vElemDisc)
			: m_vElemDisc(vElemDisc)
		{
		//	create list of all needed functions
			if(!CreateUnionOfFunctions(m_commonFctGroup, m_vElemDisc))
			{
				UG_LOG("Cannot create union of all needed functions.\n");
				throw(UGFatalError("List of function not created."));
			}

		//	create FunctionIndexMapping for each Disc
			m_vMap.resize(m_vElemDisc.size());
			for(size_t i = 0; i < m_vElemDisc.size(); ++i)
			{
				if(!CreateFunctionIndexMapping(m_vMap[i],
				                               m_vElemDisc[i]->get_function_group(),
				                               m_commonFctGroup))
				{
					UG_LOG("Cannot create Function Index Mapping for disc.\n");
					throw(UGFatalError("Function Mapping not created."));
				}
			}

			extract_imports();
			extract_exports();
		}

		const FunctionIndexMapping& map(size_t i) const
		{
			UG_ASSERT(i < m_vMap.size(), "Wrong index");
			return m_vMap[i];
		}

	///	Function group of all needed functions
		const FunctionGroup& fct_group() const {return m_commonFctGroup;}

	///	Function group for a disc
		const FunctionGroup& fct_group(size_t i) const
		{
			UG_ASSERT(i < m_vElemDisc.size(), "Invalid index");
			return m_vElemDisc[i]->get_function_group();
		}

		void clear_imports()
		{
				m_vConstImport.clear();
				m_vPosImport.clear();
				m_vIDataImport.clear();
		}

		void clear_exports()
		{
				m_vConstData.clear();
				m_vPosData.clear();
				m_vIDataExport.clear();
		}

		void extract_imports()
		{
		//	clear imports
			clear_imports();

		//	loop elem discs
			for(size_t d = 0; d < m_vElemDisc.size(); ++d)
			{
			//	loop imports
				for(size_t i = 0; i < m_vElemDisc[d]->num_imports(); ++i)
				{
				//	get import
					IDataImport<TAlgebra>* iimp = m_vElemDisc[d]->import(i);

				//	skip zero data
					if(!iimp->data_set()) continue;

				//	detect constant import
					if(iimp->constant_data())
					{
						m_vConstImport.push_back(iimp);
						continue;
					}

				//	detect position dependent import
					if(!iimp->non_zero_derivative())
					{
						m_vPosImport.push_back(iimp);
						continue;
					}

				//	detect data export dependent import
					m_vIDataImport.push_back(iimp);

				//	remember func map (the same map as for elem disc)
					m_vMapImp.push_back(m_vMap[d]);
				}
			}
		}

		void extract_exports()
		{
		//	clear exports
			clear_exports();

		//	get connected const ipdata
			for(size_t i = 0; i < m_vConstImport.size(); ++i)
				m_vConstData.push_back(m_vConstImport[i]->get_data());

		//	get connected position dependend ipdata
			for(size_t i = 0; i < m_vPosImport.size(); ++i)
				m_vPosData.push_back(m_vPosImport[i]->get_data());

		//	get connected DataExport
			for(size_t i = 0; i < m_vIDataImport.size(); ++i)
			{
				IDataExport<TAlgebra>* exp =
						dynamic_cast<IDataExport<TAlgebra>*>(m_vIDataImport[i]->get_data());

				if(exp == NULL)
					throw(UGFatalError("Something wrong in extract_exports"));

				m_vIDataExport.push_back(exp);

			//	remember func map
				FunctionIndexMapping map;
				if(!CreateFunctionIndexMapping(map,
				                               m_vIDataExport.back()->get_function_group(),
				                               m_commonFctGroup))
				{
					UG_LOG("Cannot create Function Index Mapping for disc.\n");
					throw(UGFatalError("Function Mapping not created."));
				}

				m_vMapExp.push_back(map);
			}
		}

		bool prepare_elem_loop(int id, LocalIndices& ind, IElemDiscNeed need, number time = 0.0)
		{
		//	remove ip series for all used IPData
			for(size_t i = 0; i < m_vConstData.size(); ++i)
				m_vConstData[i]->clear_ips();
			for(size_t i = 0; i < m_vPosData.size(); ++i)
				m_vPosData[i]->clear_ips();
			for(size_t i = 0; i < m_vIDataExport.size(); ++i)
				m_vIDataExport[i]->clear_ips_virt();

		// 	set elem type in elem disc
			for(size_t i = 0; i < m_vElemDisc.size(); ++i)
				if(!m_vElemDisc[i]->set_geometric_object_type(id, need))
				{
					UG_LOG("Cannot set geometric object type for Disc " << i <<".\n");
					return false;
				}

		// 	prepare loop (elem disc set local ip series here)
			for(size_t i = 0; i < m_vElemDisc.size(); ++i)
				if(!m_vElemDisc[i]->prepare_element_loop())
				{
					UG_LOG("Cannot prepare element loop.\n");
					return false;
				}

		//	evaluate constant data
			for(size_t i = 0; i < m_vConstData.size(); ++i)
				m_vConstData[i]->compute();

			for(size_t i = 0; i < m_vIDataImport.size(); ++i)
			{
			//	set id for imports
				m_vIDataImport[i]->set_geometric_object_type(id);

			//	adjust derivative array
				ind.access_by_map(m_vMapImp[i]);
				m_vIDataImport[i]->resize(ind);
			}

			for(size_t i = 0; i < m_vIDataExport.size(); ++i)
			{
			//	set id for imports
				m_vIDataExport[i]->set_geometric_object_type(id);

			//	adjust derivative array
				ind.access_by_map(m_vMapExp[i]);
				m_vIDataExport[i]->resize(ind);
			}

		//	we're done
			return true;
		}

		bool compute_elem_data(const local_vector_type & u,
		                       LocalIndices& ind,
		                       bool computeDeriv = false)
		{
		//	evaluate position data
			for(size_t i = 0; i < m_vPosData.size(); ++i)
				m_vPosData[i]->compute();

		//	evaluate export data
			for(size_t i = 0; i < m_vIDataExport.size(); ++i)
			{
				ind.access_by_map(m_vMapExp[i]);
				m_vIDataExport[i]->compute_export(u);
			}

		//	if no derivatives needed, we're done
			if(!computeDeriv) return true;

		//	compute linearized defect
			for(size_t i = 0; i < m_vIDataImport.size(); ++i)
			{
				ind.access_by_map(m_vMapImp[i]);
				m_vIDataImport[i]->compute_lin_defect(u);
			}

		//	we're done
			return true;
		}

		bool add_coupl_JA(local_matrix_type& J,
		                  LocalIndices& indRow)
		{
		//	copy local indices
			LocalIndices indCol = indRow;

		//	set different col indices
			J.set_col_indices_no_resize(indCol);

			for(size_t i = 0; i < m_vIDataImport.size(); ++i)
			{
			//	rows are given by import
				indRow.access_by_map(m_vMapImp[i]);

			//	cols are given by export
				// todo this not right iff numExp != numImp
				indCol.access_by_map(m_vMapExp[i]);

			//	add off diagonal coupling
				m_vIDataImport[i]->assemble_jacobian(J);
			}

		//	reset indices
			J.set_col_indices_no_resize(indRow);

		//	we're done
			return true;
		}

	protected:
	//	common functions
		FunctionGroup m_commonFctGroup;

	//	Function mapping for each disc
		std::vector<FunctionIndexMapping> m_vMap;

	//	current elem discs
		const std::vector<IElemDisc<TAlgebra>*>& m_vElemDisc;

	//	constant imports
		std::vector<IDataImport<TAlgebra>*> m_vConstImport;

	//	position imports
		std::vector<IDataImport<TAlgebra>*> m_vPosImport;

	//	export data imports
		std::vector<IDataImport<TAlgebra>*> m_vIDataImport;

	//	Function mapping for import
		std::vector<FunctionIndexMapping> m_vMapImp;

	//	constant data
		std::vector<IIPData*> m_vConstData;

	//	position dependend data
		std::vector<IIPData*> m_vPosData;

	//	exports
		std::vector<IDataExport<TAlgebra>*> m_vIDataExport;

	//	Function mapping for exports
		std::vector<FunctionIndexMapping> m_vMapExp;

};

} // end namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_IMPORT_EXPORT__ */
