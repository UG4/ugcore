/*
 * p1conform.h
 *
 *  Created on: 13.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__DOF_MANAGER__P1CONFORM__
#define __H__LIB_DISCRETIZATION__DOF_MANAGER__P1CONFORM__

#include <vector>

#include "lib_grid/lib_grid.h"
#include "lib_algebra/multi_index/multi_indices.h"

#include "../dof_distribution.h"
#include "../function_pattern.h"
#include "../../spacial_discretization/function_group.h"
#include "../../spacial_discretization/local_algebra.h"
#include "../../reference_element/reference_element.h"

namespace ug{

class P1StorageManager
{
//	for DofManager
	public:
		P1StorageManager() : m_pSH(NULL) {}

		void set_subset_handler(ISubsetHandler& sh)
		{
			if(m_pSH != NULL) clear();

			m_pSH = &sh;
			m_pSH->enable_subset_attachments(true);
		}

		void clear()
		{
			for(size_t si = 0; si <= m_vSubsetInfo.size(); ++si)
			{
				m_pSH->detach_from<VertexBase>(m_vSubsetInfo[si].aDoF, si);
				m_vSubsetInfo[si].aaDoFVRT.invalidate();
			}
			m_vSubsetInfo.clear();
		}

//	for the distributor
	public:
		void update_attachments()
		{
			size_t num_subsets =  m_pSH->num_subset_infos();

			// Create level dof distributors
			for(size_t si = m_vSubsetInfo.size(); si <= num_subsets; ++si)
			{
				m_vSubsetInfo.push_back(SubsetInfo());
				m_pSH->attach_to<VertexBase>(m_vSubsetInfo[si].aDoF, si);
				m_vSubsetInfo[si].aaDoFVRT.access(*m_pSH, m_vSubsetInfo[si].aDoF, si);
			}
		}

		ISubsetHandler* m_pSH;

		struct SubsetInfo
		{
			typedef ug::Attachment<size_t> ADoF;
			typedef ISubsetHandler::AttachmentAccessor<VertexBase, ADoF> attachment_accessor_type;

			attachment_accessor_type aaDoFVRT;
			ADoF aDoF;
		};

		// informations about Subsets
		std::vector<SubsetInfo> m_vSubsetInfo;
};

class P1ConformFunctionPattern : public FunctionPattern
{
	public:
		bool add_discrete_function(std::string name, LocalShapeFunctionSetID id, int dim)
		{
			// for a P1 dof manager only Lagrange P1 function space is permitted
			if(id != LSFS_LAGRANGEP1)
				{UG_LOG("P1ConformDoFDistributor: Only LSFS_LAGRANGEP1 functions are supported.\n"); return false;}

			return FunctionPattern::add_discrete_function(name, id, dim);
		}

		bool add_discrete_function(std::string name, LocalShapeFunctionSetID id, const SubsetIndexGroup& SubsetIndices, int dim)
		{
			// for a P1 dof manager only Lagrange P1 function space is permitted
			if(id != LSFS_LAGRANGEP1)
				{UG_LOG("P1ConformDoFDistributor: Only LSFS_LAGRANGEP1 functions are supported.\n"); return false;}

			return FunctionPattern::add_discrete_function(name, id, SubsetIndices, dim);
		}
};


class P1ConformDoFDistribution : public DoFDistribution
{
	public:
		// type of multiindex used
		typedef MultiIndex<2> index_type;

		// value container for element local indices
		typedef std::vector<index_type> multi_index_vector_type;

		// algebra index vector
		typedef std::vector<size_t> algebra_index_vector_type;

		// Storage type
		typedef P1StorageManager StorageManager;

	public:
		P1ConformDoFDistribution(GeometricObjectCollection goc, ISubsetHandler& sh, StorageManager& sm, FunctionPattern& dp)
		: DoFDistribution(dp), m_goc(goc), m_pISubsetHandler(&sh), m_pStorageManager(&sm), m_numDoFs(0)
		{m_vNumDoFs.clear();}

		///////////////////////////
		// Infos
		///////////////////////////

		/// number of subsets
		inline int num_subsets() const {return m_goc.num_levels();}

		/// return the number of dofs distributed
		inline size_t num_dofs() const {return m_numDoFs;}

		/// return the number of dofs distributed on subset si
		inline size_t num_dofs(int si) const {return m_vNumDoFs[si];}

		///////////////////////////////////////
		// Elements where dofs are distributed
		///////////////////////////////////////

		template<typename TElem>
		inline size_t num(int si) const {return m_goc.num<TElem>(si);}

		// iterator for elements where this grid function is defined
		template <typename TElem>
		inline typename geometry_traits<TElem>::iterator begin(int si)
			{return m_goc.begin<TElem>(si);}

		// iterator for elements where this grid function is defined
		template <typename TElem>
		inline typename geometry_traits<TElem>::iterator end(int si)
			{return m_goc.end<TElem>(si);}

		///////////////////////////////////////
		// Index access
		///////////////////////////////////////

		/// number of algebra indices for refType, subset and function group
		size_t num_indices(ReferenceObjectID refID, int si, const FunctionGroup& funcGroup) const
		{
			const ReferenceElement& refElem = ReferenceElementFactory::get_reference_element(refID);

			size_t numFct = 0;
			for(size_t fct = 0; fct < funcGroup.num_fct(); ++fct)
			{
				if(is_def_in_subset(funcGroup[fct], si))
					numFct++;
			}

			return numFct * refElem.num_obj(0);
		}

		/// fill local informations in LocalIndex
		bool prepare_indices(ReferenceObjectID refID, int si, LocalIndices& ind) const
		{
			const ReferenceElement& refElem = ReferenceElementFactory::get_reference_element(refID);

			ind.clear();
			size_t numInd = 0;
			for(size_t fct = 0; fct < ind.num_fct(); ++fct)
			{
				if(!is_def_in_subset(ind.fct_id(fct), si)) continue;

				for(size_t dof = 0; dof < refElem.num_obj(0); ++dof)
				{
					LocalIndices::dof_index_type dof_ind;
					dof_ind[0] = dof + fct * refElem.num_obj(0);
					dof_ind[1] = 0;
					ind.add_dof(fct, dof_ind);
				}
				numInd += refElem.num_obj(0);
			}
			ind.set_num_indices(numInd);
			return true;
		}

		/// fill the global algebra indices in LocalIndex
		template<typename TElem>
		void update_indices(TElem* elem, LocalIndices& ind) const
		{
			const ReferenceObjectID refID = reference_element_traits<TElem>::reference_element_type::REFERENCE_OBJECT_ID;
			const ReferenceElement& refElem = ReferenceElementFactory::get_reference_element(refID);

			for(size_t i = 0; i <  refElem.num_obj(0); ++i)
			{
				VertexBase* vrt = elem->vertex(i);
				int si = m_pISubsetHandler->get_subset_index(vrt);

				const size_t index = m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt];
				for(size_t fct = 0; fct < ind.num_fct(); ++fct)
				{
					ind.set_index(i + fct*refElem.num_obj(0), index + ind.fct_id(fct));
				}
			}
		}

		void update_indices(VertexBase* vrt, LocalIndices& ind) const
		{
			int si = m_pISubsetHandler->get_subset_index(vrt);
			const size_t index = m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt];

			for(size_t fct = 0; fct < ind.num_fct(); ++fct)
			{
				ind.set_index(fct, index + ind.fct_id(fct));
			}
		}

		template<typename TElem>
		size_t num_multi_indices(TElem* obj, size_t fct) const
		{
			typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
			return ref_elem_type::num_corners;
		}

		/// get indices of element for a function fct including boundary of element
		template<typename TElem>
		size_t get_multi_indices(TElem* elem, size_t fct, multi_index_vector_type& ind) const
		{
			typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

			const size_t numDofs = ref_elem_type::num_corners;
			ind.resize(numDofs);
			for(size_t i = 0; i < numDofs; ++i)
			{
				VertexBase* vrt = elem->vertex(i);
				int si = m_pISubsetHandler->get_subset_index(vrt);

				ind[i][0] = m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt] + fct;
				ind[i][1] = 0;
			}
			return numDofs;
		}

		template<typename TElem>
		size_t  num_multi_indices_of_geom_obj(TElem* obj, size_t fct) const
		{
			if((GeometricBaseObject)geometry_traits<TElem>::BASE_OBJECT_TYPE_ID == (GeometricBaseObject)VERTEX)
				return 1;
			else return 0;
		}

		/// get indices of element for a function fct excluding boundary of element
		template<typename TElem>
		size_t get_multi_indices_of_geom_obj(TElem* obj, size_t fct, multi_index_vector_type& ind) const
		{
			if((GeometricBaseObject)geometry_traits<TElem>::BASE_OBJECT_TYPE_ID == (GeometricBaseObject)VERTEX)
			{
				VertexBase* vrt = (VertexBase*)obj;
				int si = m_pISubsetHandler->get_subset_index(vrt);
				ind.resize(1);
				ind[0][0] = m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt] + fct;
				ind[0][1] = 0;
				return 1;
			}
			else return 0;
		}

		/// get algebra indices of element including boundary of element
		template<typename TElem>
		void get_algebra_indices(TElem* elem, algebra_index_vector_type& ind) const
		{
			typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

			ind.clear();

			const int elem_si = m_pISubsetHandler->get_subset_index(elem);
			for(size_t fct = 0; fct < num_fct(elem_si); ++fct)
			{
				for(size_t i = 0; i < ref_elem_type::num_corners; ++i)
				{
					VertexBase* vrt = elem->vertex(i);
					int si = m_pISubsetHandler->get_subset_index(vrt);
					const size_t index = m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt] + fct;
					ind.push_back(index);
				}
			}
		}

		/// get indices of element for a function fct excluding boundary of element
		template<typename TElem>
		void get_algebra_indices_of_geom_obj(TElem* obj, algebra_index_vector_type& ind) const
		{
			ind.clear();
			if((GeometricBaseObject)geometry_traits<TElem>::BASE_OBJECT_TYPE_ID == (GeometricBaseObject)VERTEX)
			{
				VertexBase* vrt = (VertexBase*)obj;
				int si = m_pISubsetHandler->get_subset_index(vrt);

				const size_t index =  m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt];
				for(size_t fct = 0; fct < num_fct(si); ++fct)
				{
					ind.push_back(index + fct);
				}
			}
		}

		///////////////////////////
		// Creation
		///////////////////////////

		// distribute dofs for given goc
		bool distribute_dofs()
		{
			// Attach indices
			m_pStorageManager->update_attachments();

			// iterators
			geometry_traits<VertexBase>::iterator iter, iterBegin, iterEnd;

			// loop subsets
			size_t i = 0;

			// reset number of dofs
			m_vNumDoFs.clear(); m_vNumDoFs.resize(num_subsets());

			// loop subsets
			for(int si = 0; si < num_subsets(); ++si)
			{
				// skip if no dofs to be distributed
				if(!(m_pFunctionPattern->num_fct(si)>0)) continue;

				// remember number of functions
				size_t num_fct =  m_pFunctionPattern->num_fct(si);

				iterBegin = m_goc.begin<VertexBase>(si);
				iterEnd =  m_goc.end<VertexBase>(si);

				// loop Verices
				for(iter = iterBegin; iter != iterEnd; ++iter)
				{
					// get vertex
					VertexBase* vrt = *iter;

					// write index
					m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt] = i;
					i += num_fct;
				}
				if(si==0) m_vNumDoFs[si] = i;
				else m_vNumDoFs[si] = i - m_vNumDoFs[si-1];
			}
			m_numDoFs = i;

			UG_LOG(std::setw(8) << m_numDoFs <<" |     " << 1 << "     | " );

			for(int si = 0; si < num_subsets(); ++si)
			{
				UG_LOG( " (" << si << ","<<1<<"," << std::setw(8) << m_vNumDoFs[si] << ") ");
			}
			UG_LOG("\n");

			return true;
		}

	protected:
		// geometric object collection for this Distributor
		GeometricObjectCollection m_goc;

		// subset handler for this distributor
		ISubsetHandler* m_pISubsetHandler;

		// Storage Manager for dofs
		StorageManager* m_pStorageManager;

		// number of distributed dofs on whole domain
		size_t m_numDoFs;

		// number of distributed dofs on each subset
		std::vector<size_t> m_vNumDoFs;
};


class GroupedP1ConformDoFDistribution : public DoFDistribution
{
	public:
		// type of multiindex used
		typedef MultiIndex<2> index_type;

		// value container for element local indices
		typedef std::vector<index_type> multi_index_vector_type;

		// algebra index vector
		typedef std::vector<size_t> algebra_index_vector_type;

		// Storage type
		typedef P1StorageManager StorageManager;

	public:
		GroupedP1ConformDoFDistribution(GeometricObjectCollection goc, ISubsetHandler& sh, StorageManager& sm, FunctionPattern& dp)
		: DoFDistribution(dp), m_goc(goc), m_pISubsetHandler(&sh), m_pStorageManager(&sm), m_numDoFs(0)
		{m_vNumDoFs.clear();}

		///////////////////////////
		// Infos
		///////////////////////////

		/// number of subsets
		inline int num_subsets() const {return m_goc.num_levels();}

		/// return the number of dofs distributed
		inline size_t num_dofs() const {return m_numDoFs;}

		/// return the number of dofs distributed on subset si
		inline size_t num_dofs(int si) const {return m_vNumDoFs[si];}

		///////////////////////////////////////
		// Elements where dofs are distributed
		///////////////////////////////////////

		template<typename TElem>
		inline size_t num(int si) const {return m_goc.num<TElem>(si);}

		// iterator for elements where this grid function is defined
		template <typename TElem>
		inline typename geometry_traits<TElem>::iterator begin(int si)
			{return m_goc.begin<TElem>(si);}

		// iterator for elements where this grid function is defined
		template <typename TElem>
		inline typename geometry_traits<TElem>::iterator end(int si)
			{return m_goc.end<TElem>(si);}

		///////////////////////////////////////
		// Index access
		///////////////////////////////////////

		/// number of algebra indices for refType, subset and function group
		size_t num_indices(ReferenceObjectID refID, int si, const FunctionGroup& funcGroup) const
		{
			const ReferenceElement& refElem = ReferenceElementFactory::get_reference_element(refID);

			for(size_t fct = 0; fct < funcGroup.num_fct(); ++fct)
			{
				if(is_def_in_subset(funcGroup[fct], si))
					return refElem.num_obj(0);
			}

			return 0;
		}

		/// fill local informations in LocalIndex
		bool prepare_indices(ReferenceObjectID refID, int si, LocalIndices& ind) const
		{
			const ReferenceElement& refElem = ReferenceElementFactory::get_reference_element(refID);

			ind.clear();
			size_t numInd = 0;
			for(size_t fct = 0; fct < ind.num_fct(); ++fct)
			{
				if(!is_def_in_subset(ind.fct_id(fct), si)) continue;

				for(size_t dof = 0; dof < refElem.num_obj(0); ++dof)
				{
					LocalIndices::dof_index_type dof_ind;
					dof_ind[0] = dof;
					dof_ind[1] = fct;
					ind.add_dof(fct, dof_ind);
				}
				numInd = refElem.num_obj(0);
			}
			ind.set_num_indices(numInd);
			return true;
		}

		/// fill the global algebra indices in LocalIndex
		template<typename TElem>
		void update_indices(TElem* elem, LocalIndices& ind) const
		{
			const ReferenceObjectID refID = reference_element_traits<TElem>::reference_element_type::REFERENCE_OBJECT_ID;
			const ReferenceElement& refElem = ReferenceElementFactory::get_reference_element(refID);

			for(size_t i = 0; i <  refElem.num_obj(0); ++i)
			{
				VertexBase* vrt = elem->vertex(i);
				int si = m_pISubsetHandler->get_subset_index(vrt);

				const size_t index = m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt];
				ind.set_index(i, index);
			}
		}

		void update_indices(VertexBase* vrt, LocalIndices& ind) const
		{
			int si = m_pISubsetHandler->get_subset_index(vrt);

			const size_t index = m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt];
			ind.set_index(0, index);
		}


		template<typename TElem>
		size_t num_multi_indices(TElem* obj, size_t fct) const
		{
			typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
			return ref_elem_type::num_corners;
		}

		/// get indices of element for a function fct including boundary of element
		template<typename TElem>
		size_t get_multi_indices(TElem* elem, size_t fct, multi_index_vector_type& ind) const
		{
			typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

			const size_t numDofs = ref_elem_type::num_corners;
			ind.resize(numDofs);
			for(size_t i = 0; i < numDofs; ++i)
			{
				VertexBase* vrt = elem->vertex(i);
				int si = m_pISubsetHandler->get_subset_index(vrt);

				ind[i][0] = m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt];
				ind[i][1] = fct;
			}
			return numDofs;
		}

		template<typename TElem>
		size_t  num_multi_indices_of_geom_obj(TElem* obj, size_t fct) const
		{
			if((GeometricBaseObject)geometry_traits<TElem>::BASE_OBJECT_TYPE_ID == (GeometricBaseObject)VERTEX)
				return 1;
			else return 0;
		}

		/// get indices of element for a function fct excluding boundary of element
		template<typename TElem>
		size_t get_multi_indices_of_geom_obj(TElem* obj, size_t fct, multi_index_vector_type& ind) const
		{
			if((GeometricBaseObject)geometry_traits<TElem>::BASE_OBJECT_TYPE_ID == (GeometricBaseObject)VERTEX)
			{
				VertexBase* vrt = (VertexBase*)obj;
				int si = m_pISubsetHandler->get_subset_index(vrt);
				ind.resize(1);
				ind[0][0] = m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt];
				ind[0][1] = fct;
				return 1;
			}
			else return 0;
		}

		/// get algebra indices of element including boundary of element
		template<typename TElem>
		void get_algebra_indices(TElem* elem, algebra_index_vector_type& ind) const
		{
			typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

			ind.clear();

			// if no functions, return
			const int elem_si = m_pISubsetHandler->get_subset_index(elem);
			if(num_fct(elem_si) == 0) return;

			for(size_t i = 0; i < ref_elem_type::num_corners; ++i)
			{
					VertexBase* vrt = elem->vertex(i);
					int si = m_pISubsetHandler->get_subset_index(vrt);
					const size_t index = m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt];
					ind.push_back(index);
			}
		}

		/// get indices of element for a function fct excluding boundary of element
		template<typename TElem>
		void get_algebra_indices_of_geom_obj(TElem* obj, algebra_index_vector_type& ind) const
		{
			ind.clear();
			if((GeometricBaseObject)geometry_traits<TElem>::BASE_OBJECT_TYPE_ID == (GeometricBaseObject)VERTEX)
			{
				VertexBase* vrt = (VertexBase*)obj;
				int si = m_pISubsetHandler->get_subset_index(vrt);
				if(num_fct(si) == 0) return;

				const size_t index =  m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt];
				ind.push_back(index);
			}
		}

		///////////////////////////
		// Creation
		///////////////////////////

		// distribute dofs for given goc
		bool distribute_dofs()
		{
			// Attach indices
			m_pStorageManager->update_attachments();

			// iterators
			geometry_traits<VertexBase>::iterator iter, iterBegin, iterEnd;

			// loop subsets
			size_t i = 0;

			// reset number of dofs
			m_vNumDoFs.clear(); m_vNumDoFs.resize(num_subsets());

			// loop subsets
			for(int si = 0; si < num_subsets(); ++si)
			{
				// skip if no dofs to be distributed
				if(!(m_pFunctionPattern->num_fct(si)>0)) continue;

				iterBegin = m_goc.begin<VertexBase>(si);
				iterEnd =  m_goc.end<VertexBase>(si);

				// loop Verices
				for(iter = iterBegin; iter != iterEnd; ++iter)
				{
					// get vertex
					VertexBase* vrt = *iter;

					// write index
					m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt] = i;
					i ++;
				}
				if(si==0) m_vNumDoFs[si] = i;
				else m_vNumDoFs[si] = i - m_vNumDoFs[si-1];
			}
			m_numDoFs = i;

			UG_LOG(std::setw(8) << m_numDoFs <<" | " << "variable" << "  | " );

			for(int si = 0; si < num_subsets(); ++si)
			{
				size_t num_fct =  m_pFunctionPattern->num_fct(si);
				UG_LOG( " (" << si << "," <<num_fct<<","<< std::setw(8) << m_vNumDoFs[si] << ") ");
			}
			UG_LOG("\n");

			return true;
		}

	protected:
		// geometric object collection for this Distributor
		GeometricObjectCollection m_goc;

		// subset handler for this distributor
		ISubsetHandler* m_pISubsetHandler;

		// Storage Manager for dofs
		StorageManager* m_pStorageManager;

		// number of distributed dofs on whole domain
		size_t m_numDoFs;

		// number of distributed dofs on each subset
		std::vector<size_t> m_vNumDoFs;
};


} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__DOF_MANAGER__P1CONFORM__ */
