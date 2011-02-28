/*
 * data_export.h
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
			: m_pIDataExport(NULL), m_id(-1), m_bCompLinDefect(compLinDefect)
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
		bool non_zero_derivative() const
		{
			if(m_pIDataExport == NULL) return false;
			else return m_bCompLinDefect;
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

	///	indicates iff lin defect should be computed
		bool m_bCompLinDefect;

	public:
	///	sets the geometric object type
		bool set_geometric_object_type(int id)
		{
			if(id < (int)m_vLinDefectFunc.size() && m_vLinDefectFunc[id] != NULL)
			{
				m_id = id;
				return true;
			}
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

	///	returns true if data given
		virtual bool data_given() const {return !(m_pIPData == NULL);}

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
	///	sets the elem discs to evaluate
		bool set_elem_discs(const std::vector<IElemDisc<TAlgebra>*>& vElemDisc)
		{
		//	remember current elem discs
			m_pvElemDisc = &vElemDisc;

		//	check that all disc have the correct number of functions
			for(size_t i = 0; i < m_pvElemDisc->size(); ++i)
			{
				size_t fct_given = (*m_pvElemDisc)[i]->get_function_group().num_fct();
				size_t fct_needed = (*m_pvElemDisc)[i]->num_fct();

				if(fct_given != fct_needed)
				{
					UG_LOG("ERROR in 'DataEvaluator::set_elem_discs': Exactly "
							<< fct_needed << " functions needed, but given " <<
							fct_given <<" functions.\n");
					return false;
				}

			}

		//	create list of all needed functions
			if(!CreateUnionOfFunctions(m_commonFctGroup, *m_pvElemDisc))
			{
				UG_LOG("ERROR in 'DataEvaluator::set_elem_discs':"
						"Cannot create union of all needed functions.\n");
				return false;
			}

		//	create FunctionIndexMapping for each Disc
			m_vMap.resize(m_pvElemDisc->size());
			for(size_t i = 0; i < m_pvElemDisc->size(); ++i)
			{
				if(!CreateFunctionIndexMapping(m_vMap[i],
											   (*m_pvElemDisc)[i]->get_function_group(),
											   m_commonFctGroup))
				{
					UG_LOG("ERROR in 'DataEvaluator::set_elem_discs':"
							"Cannot create Function Index Mapping for disc.\n");
					return false;
				}
			}

			return extract_imports_and_ipdata();
		}

	///	Mapping between local functions of ElemDisc i and common FunctionGroup
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
			UG_ASSERT(i < m_pvElemDisc->size(), "Invalid index");
			return (*m_pvElemDisc)[i]->get_function_group();
		}

	///	clears imports and ip data
		void clear()
		{
		//	remove imports scheduled for evaluation of lin defect
			m_vMapImp.clear();
			m_vIDataImport.clear();

		// 	remove data scheduled for evaluation
			m_vConstData.clear();
			m_vPosData.clear();

			m_vDependData.clear();
			m_vMapDepend.clear();

			m_vMapExp.clear();
			m_vIDataExport.clear();

			m_vMapLinker.clear();
			m_vLinkerDepend.clear();
			m_vLinkerData.clear();
		}

	//	this function tries to add the last entry of vTryingToAdd to the eval
	//	data
		bool add_data_to_eval_data(std::vector<IIPData*>& vEvalData,
		                           std::vector<IIPData*>& vTryingToAdd)
		{
		//	if empty, we're done
			if(vTryingToAdd.empty()) return true;

		//	search for element in already scheduled data
			std::vector<IIPData*>::iterator it, itEnd;
			it = find(vEvalData.begin(), vEvalData.end(), vTryingToAdd.back());

		//	if found, skip this data
			if(it != vEvalData.end())
			{
				vTryingToAdd.pop_back();
				return true;
			}

		//	search if element already contained in list. Then, the element
		//	did start the adding procedure before and a circle dependency
		//	is found
			itEnd = vTryingToAdd.end(); itEnd--;
			it = find(vTryingToAdd.begin(), itEnd, *itEnd);

		//	if found, return error of circle dependency
			if(it != itEnd)
			{
				UG_LOG("ERROR in add_data_to_eval_data:"
						" Circle dependency of data detected for IP Data.\n");
				return false;
			}

		//	add all dependent datas
			IIPData* data = vTryingToAdd.back();
			for(size_t i = 0; i < data->num_needed_data(); ++i)
			{
			//	add each data separately
				vTryingToAdd.push_back(data->needed_data(i));
				if(!add_data_to_eval_data(vEvalData, vTryingToAdd))
					return false;
			}

		//	add this data to the evaluation list
			vEvalData.push_back(data);

		//	pop last one, since now added to eval list
			vTryingToAdd.pop_back();

		//	we're done
			return true;
		}

		bool extract_imports_and_ipdata()
		{
		//	clear imports and ipdata
			clear();

		//	queue for all ip data needed
			std::vector<IIPData*> vEvalData;
			std::vector<IIPData*> vTryingToAdd;

		//	loop elem discs
			for(size_t d = 0; d < m_pvElemDisc->size(); ++d)
			{
			//	loop imports
				for(size_t i = 0; i < (*m_pvElemDisc)[d]->num_imports(); ++i)
				{
				//	get import
					IDataImport<TAlgebra>* iimp = (*m_pvElemDisc)[d]->import(i);

				//	skip zero data
					if(!iimp->data_given()) continue;

				//	push export on stack of needed data
					vTryingToAdd.push_back(iimp->get_data());

				//	add data and all dependency to evaluation list
					if(!add_data_to_eval_data(vEvalData, vTryingToAdd))
					{
						UG_LOG("ERROR in extract_imports_and_ipdata:"
								" Circle dependency of data detected for IP Data.\n");
						return false;
					}
					UG_ASSERT(vTryingToAdd.empty(), "All must have been added.");

				//	done iff zero-derivative
					if(!iimp->non_zero_derivative()) continue;

				//	schedule import for computation of lin defect:
					m_vIDataImport.push_back(iimp);

				//	Remember FuncMap for lin defect (the same map as for elem disc)
					m_vMapImp.push_back(m_vMap[d]);
				}
			}

		//	loop all needed ip data and sort it
			for(size_t i = 0; i < vEvalData.size(); ++i)
			{
				IIPData* ipData = vEvalData[i];

			//	detect constant import
				if(ipData->constant_data())
				{
				//	schedule for evaluation of constant data
					m_vConstData.push_back(ipData);
					continue;
				}

			//	detect position dependent import
				if(ipData->zero_derivative())
				{
				//	schedule for evaluation of position dependent data
					m_vPosData.push_back(ipData);
					continue;
				}

			//	cast to dependent data
				IDependentIPData* dependData =
						dynamic_cast<IDependentIPData*>(ipData);

			//	check success
				if(dependData == NULL)
					{return false;}

			//	create FuncMap
				FunctionIndexMapping map;
				if(!CreateFunctionIndexMapping(map,
				                               dependData->get_function_group(),
											   m_commonFctGroup))
				{
					UG_LOG("Cannot create Function Index Mapping for disc.\n");
					return false;
				}


			//	save as dependent data
				m_vDependData.push_back(ipData);
				m_vMapDepend.push_back(map);

			//	cast to data export
				IDataExport<TAlgebra>* exp =
						dynamic_cast<IDataExport<TAlgebra>*>(ipData);

			//	Data Export case
				if(exp != NULL)
				{
				//	schedule for evaluation of IDataExports
					m_vIDataExport.push_back(exp);

				// 	remember function map
					m_vMapExp.push_back(map);
				}
				else
			//	Linker case
				{
					//	schedule for evaluation of linker
						m_vLinkerData.push_back(ipData);
						m_vLinkerDepend.push_back(dependData);

					// 	remember function map
						m_vMapLinker.push_back(map);
				}
			}

		//	we're done
			return true;
		}

		bool prepare_elem_loop(int id, LocalIndices& ind, IElemDiscNeed need, number time = 0.0)
		{
		//	remove ip series for all used IPData
			for(size_t i = 0; i < m_vConstData.size(); ++i)
				m_vConstData[i]->clear_ips();
			for(size_t i = 0; i < m_vPosData.size(); ++i)
				m_vPosData[i]->clear_ips();
			for(size_t i = 0; i < m_vIDataExport.size(); ++i)
				m_vIDataExport[i]->clear_export_ips();
			for(size_t i = 0; i < m_vLinkerData.size(); ++i)
				m_vLinkerData[i]->clear_ips();

		// 	set elem type in elem disc
			for(size_t i = 0; i < m_pvElemDisc->size(); ++i)
				if(!(*m_pvElemDisc)[i]->set_geometric_object_type(id, need))
				{
					UG_LOG("In 'DataEvaluator::prepare_elem_loop':"
							"Cannot set geometric object type for Disc " << i <<".\n");
					return false;
				}

		// 	prepare loop (elem disc set local ip series here)
			for(size_t i = 0; i < m_pvElemDisc->size(); ++i)
				if(!(*m_pvElemDisc)[i]->prepare_element_loop())
				{
					UG_LOG("In 'DataEvaluator::prepare_elem_loop':"
							"Cannot prepare element loop.\n");
					return false;
				}

		//	prepare data imports
			for(size_t i = 0; i < m_vIDataImport.size(); ++i)
			{
			//	set id for imports
				m_vIDataImport[i]->set_geometric_object_type(id);

			//	adjust derivative array
				ind.access_by_map(m_vMapImp[i]);
				m_vIDataImport[i]->resize(ind);
			}

		//	prepare data exports
			for(size_t i = 0; i < m_vIDataExport.size(); ++i)
			{
			//	set id for imports
				m_vIDataExport[i]->set_geometric_object_type(id);

			//	adjust derivative array
				ind.access_by_map(m_vMapExp[i]);
				m_vIDataExport[i]->resize(ind);
			}

		//	prepare data linker
			for(size_t i = 0; i < m_vLinkerDepend.size(); ++i)
			{
			//	adjust derivative array
				ind.access_by_map(m_vMapLinker[i]);
				m_vLinkerDepend[i]->resize(ind);
			}

		//	evaluate constant data
			for(size_t i = 0; i < m_vConstData.size(); ++i)
				m_vConstData[i]->compute();

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

		// 	process dependent data
			for(size_t i = 0; i < m_vDependData.size(); ++i)
			{
			//	cast to data export
				IDataExport<TAlgebra>* exp =
						dynamic_cast<IDataExport<TAlgebra>*>(m_vDependData[i]);

			//	linker case
				if(exp == NULL)
				{
					m_vDependData[i]->compute(computeDeriv);
				}
			//	export case
				else
				{
					ind.access_by_map(m_vMapDepend[i]);
					exp->compute_export(u, computeDeriv);
				}
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
		const std::vector<IElemDisc<TAlgebra>*>* m_pvElemDisc;

	//	data imports which are connected to non-zero derivative ip data
		std::vector<IDataImport<TAlgebra>*> m_vIDataImport;

	//	Function mapping for import
		std::vector<FunctionIndexMapping> m_vMapImp;

	//	constant data
		std::vector<IIPData*> m_vConstData;

	//	position dependend data
		std::vector<IIPData*> m_vPosData;

	//	all dependent data
		std::vector<IIPData*> m_vDependData;
		std::vector<FunctionIndexMapping> m_vMapDepend;

	//	data linker
		std::vector<IIPData*> m_vLinkerData;
		std::vector<IDependentIPData*> m_vLinkerDepend;

	//	exports
		std::vector<IDataExport<TAlgebra>*> m_vIDataExport;

	//	Function mapping for exports
		std::vector<FunctionIndexMapping> m_vMapExp;

	//	Function mapping for exports
		std::vector<FunctionIndexMapping> m_vMapLinker;

};

} // end namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_IMPORT__ */
