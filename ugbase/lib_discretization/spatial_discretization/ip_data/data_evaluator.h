/*
 * data_evaluator.h
 *
 *  Created on: 17.12.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_EVALUATOR__
#define __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_EVALUATOR__

#include "user_data.h"
#include "../elem_disc/elem_disc_interface.h"
#include "data_export.h"
#include "data_import.h"

namespace ug{

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
			m_vMapImpConn.clear();
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

				//	skip non-given data (no need for evaluation)
					if(!iimp->data_given()) continue;

				//	push export on stack of needed data
					vTryingToAdd.push_back(iimp->get_data());

				//	add data and all dependency to evaluation list
					if(!add_data_to_eval_data(vEvalData, vTryingToAdd))
					{
						UG_LOG("ERROR in DataEvaluator::extract_imports_and_ipdata:"
								" Circle dependency of data detected for IP Data.\n");
						return false;
					}
					UG_ASSERT(vTryingToAdd.empty(), "All must have been added.");

				//	done iff zero-derivative
					if(iimp->zero_derivative()) continue;

				//	schedule import for computation of lin defect:
					m_vIDataImport.push_back(iimp);

				//	cast to dependent data
					IDependentIPData* dependData =
							dynamic_cast<IDependentIPData*>(iimp->get_data());

				//	check success
					if(dependData == NULL)
					{
						UG_LOG("ERROR in DataEvaluator::extract_imports_and_ipdata:"
								" Data seems dependent, but cast failed.\n");
						return false;
					}

				//	create FuncMap
					FunctionIndexMapping map;
					if(!CreateFunctionIndexMapping(map,
												   dependData->get_function_group(),
												   m_commonFctGroup))
					{
						UG_LOG("Cannot create Function Index Mapping for disc.\n");
						return false;
					}


				//	Remember FuncMap for lin defect (the same map as for elem disc)
					m_vMapImp.push_back(m_vMap[d]);

				//  Remember FuncMap for connected data
					m_vMapImpConn.push_back(map);
				}
			}

		//	loop all needed ip data and group it
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
				{
					UG_LOG("ERROR in DataEvaluator::extract_imports_and_ipdata:"
							" Data seems dependent, but cast failed.\n");
					return false;
				}

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
				if(!m_vIDataImport[i]->set_geometric_object_type(id))
				{
					UG_LOG("In 'DataEvaluator::prepare_elem_loop':"
							"Cannot set geometric object type for Import " << i <<".\n");
					return false;
				}

			//	adjust lin defect array
				ind.access_by_map(m_vMapImp[i]);
				m_vIDataImport[i]->resize(ind);
			}

		//	prepare data exports
			for(size_t i = 0; i < m_vIDataExport.size(); ++i)
			{
			//	set id for imports
				if(!m_vIDataExport[i]->set_geometric_object_type(id))
				{
					UG_LOG("In 'DataEvaluator::prepare_elem_loop':"
							"Cannot set geometric object type for Export " << i <<".\n");
					return false;
				}

			//	adjust derivative array
				ind.access_by_map(m_vMapExp[i]);
				m_vIDataExport[i]->resize(ind);
			}

		//	prepare data linker
			for(size_t i = 0; i < m_vLinkerDepend.size(); ++i)
			{
				if(!m_vLinkerDepend[i]->make_ready())
				{
					UG_LOG("In 'DataEvaluator::prepare_elem_loop':"
							"Linker not ready.\n");
					return false;
				}

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
					if(!exp->compute_export(u, computeDeriv))
					{
						UG_LOG("In 'DataEvaluator::compute_elem_data':"
								"Cannot compute data for Export " << i <<".\n");
						return false;
					}
				}
			}

		//	if no derivatives needed, we're done
			if(!computeDeriv) return true;

		//	compute linearized defect
			for(size_t i = 0; i < m_vIDataImport.size(); ++i)
			{
				ind.access_by_map(m_vMapImp[i]);
				if(!m_vIDataImport[i]->compute_lin_defect(u))
				{
					UG_LOG("In 'DataEvaluator::compute_elem_data':"
							"Cannot compute lin defect for Import " << i <<".\n");
					return false;
				}
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
				indCol.access_by_map(m_vMapImpConn[i]);

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

	//	Function mapping for import and for the connected data
		std::vector<FunctionIndexMapping> m_vMapImp;
		std::vector<FunctionIndexMapping> m_vMapImpConn;

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
		std::vector<FunctionIndexMapping> m_vMapLinker;

	//	exports
		std::vector<IDataExport<TAlgebra>*> m_vIDataExport;
		std::vector<FunctionIndexMapping> m_vMapExp;
};

} // end namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_EVALUATOR__ */
