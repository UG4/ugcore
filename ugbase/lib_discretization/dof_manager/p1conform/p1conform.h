/*
 * p1conform.h
 *
 *  Created on: 13.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__DOF_MANAGER__P1CONFORM__
#define __H__LIB_DISCRETIZATION__DOF_MANAGER__P1CONFORM__

#include <vector>

#include "lib_grid/lg_base.h"

#include "../dof_distribution.h"
#include "../function_pattern.h"
#include "../../reference_element/reference_element.h"

namespace ug{

class P1StorageManager
{
//	for DofManager
	public:
		P1StorageManager() : m_pSH(NULL) {m_vSubsetInfo.clear();}

	/// set subset handler
		void set_subset_handler(ISubsetHandler& sh);

	/// remove subset handler
		void clear_subset_handler();

	/// clear all dofs
		void clear();

	/// destructor
		~P1StorageManager() {clear();};

	public:
	/// attach dofs
		void update_attachments();

	/// subset handler
		ISubsetHandler* m_pSH;

	/// Subset Infos
		struct SubsetInfo
		{
			typedef ug::Attachment<size_t> ADoF;
			typedef ISubsetHandler::AttachmentAccessor<VertexBase, ADoF>
					attachment_accessor_type;

			attachment_accessor_type aaDoFVRT;
			ADoF aDoF;
		};

	// informations about Subsets
		std::vector<SubsetInfo> m_vSubsetInfo;
};

class P1ConformDoFDistribution
	: public IDoFDistribution<P1ConformDoFDistribution>
{
	public:
	///	Base class
		typedef IDoFDistribution<P1ConformDoFDistribution> base_type;

	/// \copydoc IDoFDistribution::index_type
		typedef base_type::index_type index_type;

	/// \copydoc IDoFDistribution::multi_index_vector_type
		typedef base_type::multi_index_vector_type multi_index_vector_type;

	/// \copydoc IDoFDistribution::algebra_index_vector_type
		typedef base_type::algebra_index_vector_type algebra_index_vector_type;

	/// Storage Manager type
		typedef P1StorageManager storage_manager_type;

	public:
		P1ConformDoFDistribution(GeometricObjectCollection goc,
		                         ISubsetHandler& sh, storage_manager_type& sm,
		                         FunctionPattern& fp)
		: base_type(goc, fp), m_pISubsetHandler(&sh),
		  m_pStorageManager(&sm), m_numDoFs(0)
		{
			m_vNumDoFs.clear();
			m_vNumDoFs.resize(this->num_subsets(), 0);
		}

		P1ConformDoFDistribution(GeometricObjectCollection goc,
		                         ISubsetHandler& sh, storage_manager_type& sm,
		                         FunctionPattern& fp,
		                         const SurfaceView& surfView)
		: base_type(goc, fp, surfView), m_pISubsetHandler(&sh),
		  m_pStorageManager(&sm), m_numDoFs(0)
		{
			m_vNumDoFs.clear();
			m_vNumDoFs.resize(this->num_subsets(), 0);
		}

		///////////////////////////
		// Support Info
		///////////////////////////

		static bool supports_trial_space(LocalShapeFunctionSetID& id)
		{
			return id == LocalShapeFunctionSetID(LocalShapeFunctionSetID::LAGRANGE, 1);
		}

		///////////////////////////
		// Infos
		///////////////////////////

	/// \copydoc ug::IDoFDistribution::num_dofs()
		size_t num_dofs() const {return m_numDoFs;}

	/// \copydoc ug::IDoFDistribution::num_dofs(int) const
		size_t num_dofs(int si) const {return m_vNumDoFs[si];}

	/// \copydoc IDoFDistribution::blocksize()
		int blocksize() const {return 1;}

	/// \copydoc IDoFDistribution::blocksize(int) const
		int blocksize(int si) const {return 1;}

		///////////////////////////////////////
		// LocalIndex update
		///////////////////////////////////////

	/// \copydoc IDoFDistribution::num_indices()
		size_t num_indices(ReferenceObjectID refID, int si,
		                   const FunctionGroup& funcGroup) const;

	/// \copydoc IDoFDistribution::num_inner_indices()
		size_t num_inner_indices(ReferenceObjectID refID, int si,
		                         const FunctionGroup& funcGroup) const;

	/// \copydoc IDoFDistribution::prepare_indices()
		bool prepare_indices(ReferenceObjectID refID, int si,
		                     LocalIndices& ind, bool withHanging = false) const;

	/// \copydoc IDoFDistribution::prepare_inner_indices()
		bool prepare_inner_indices(ReferenceObjectID refID, int si,
		                           LocalIndices& ind) const;

	/// \copydoc IDoFDistribution::update_indices()
		template<typename TElem>
		void update_indices(TElem* elem, LocalIndices& ind,
		                    bool withHanging = false) const;

	/// \copydoc IDoFDistribution::update_inner_indices()
		template<typename TElem>
		void update_inner_indices(TElem* elem, LocalIndices& ind) const;


		///////////////////////////////////////
		// Multi index access
		///////////////////////////////////////

	/// \copydoc IDoFDistribution::num_multi_indices()
		template<typename TElem>
		size_t num_multi_indices(TElem* elem, size_t fct) const;

	/// \copydoc IDoFDistribution::num_inner_multi_indices()
		template<typename TElem>
		size_t num_inner_multi_indices(TElem* elem, size_t fct) const;

	/// \copydoc IDoFDistribution::get_multi_indices()
		template<typename TElem>
		size_t get_multi_indices(TElem* elem, size_t fct,
		                         multi_index_vector_type& ind) const;

	/// \copydoc IDoFDistribution::get_inner_multi_indices()
		template<typename TElem>
		size_t get_inner_multi_indices(TElem* elem, size_t fct,
		                               multi_index_vector_type& ind) const;

		///////////////////////////////////////
		// Algebra index access
		///////////////////////////////////////

	/// \copydoc IDoFDistribution::num_algebra_indices()
		template<typename TElem>
		size_t num_algebra_indices(TElem* elem, size_t fct) const;

	/// \copydoc IDoFDistribution::num_inner_algebra_indices()
		template<typename TElem>
		size_t num_inner_algebra_indices(TElem* elem, size_t fct) const;

	/// \copydoc IDoFDistribution::get_algebra_indices()
		template<typename TElem>
		void get_algebra_indices(TElem* elem,
		                         algebra_index_vector_type& ind) const;

	/// \copydoc IDoFDistribution::get_inner_algebra_indices()
		template<typename TElem>
		void get_inner_algebra_indices(TElem* elem,
		                               algebra_index_vector_type& ind) const;

		///////////////////////////
		// Creation
		///////////////////////////

	/// \copydoc IDoFDistribution::distribute_dofs()
		bool distribute_dofs();

	/// \copydoc IDoFDistribution::permute_indices()
		bool permute_indices(const std::vector<size_t>& vIndNew);

	/// \copydoc IDoFDistribution::get_connections()
		bool get_connections(std::vector<std::vector<size_t> >& vvConnection);

	/// \copydoc IDoFDistribution::vertices_created()
		bool vertex_added(VertexBase* vrt);

	/// \copydoc IDoFDistribution::vertices_to_be_erased()
		bool vertex_to_be_removed(VertexBase* vrt);

	protected:
	///	returns first algebra index of a vertex
		size_t& first_index(VertexBase* vrt, size_t si)
		{
			UG_ASSERT(m_pStorageManager != NULL, "No Storage Manager");
			UG_ASSERT(m_pISubsetHandler != NULL, "No Subset Handler");
			return m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt];
		}

	///	const access to first algebra index of a vertex
		const size_t& first_index(VertexBase* vrt, size_t si) const
		{
			UG_ASSERT(m_pStorageManager != NULL, "No Storage Manager");
			UG_ASSERT(m_pISubsetHandler != NULL, "No Subset Handler");
			return m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt];
		}

	protected:
	/// subset handler for this distributor
		ISubsetHandler* m_pISubsetHandler;

	// 	Storage Manager for dofs
		storage_manager_type* m_pStorageManager;

	/// number of distributed dofs on whole domain
		size_t m_numDoFs;

	/// number offset map
		std::vector<std::vector<size_t> > m_vvOffsets;

	/// number of distributed dofs on each subset
		std::vector<size_t> m_vNumDoFs;
};


class GroupedP1ConformDoFDistribution
	: public IDoFDistribution<GroupedP1ConformDoFDistribution>
{
	public:
	///	Base class
		typedef IDoFDistribution<GroupedP1ConformDoFDistribution> base_type;

	/// \copydoc IDoFDistribution::index_type
		typedef base_type::index_type index_type;

	/// \copydoc IDoFDistribution::multi_index_vector_type
		typedef base_type::multi_index_vector_type multi_index_vector_type;

	/// \copydoc IDoFDistribution::algebra_index_vector_type
		typedef base_type::algebra_index_vector_type algebra_index_vector_type;

	/// Storage Manager type
		typedef P1StorageManager storage_manager_type;

	public:
		GroupedP1ConformDoFDistribution(GeometricObjectCollection goc,
		                                ISubsetHandler& sh,
		                                storage_manager_type& sm,
		                                FunctionPattern& dp)
		: base_type(goc, dp), m_pISubsetHandler(&sh),
		  m_pStorageManager(&sm), m_numDoFs(0)
		{
			m_vNumDoFs.clear();
			m_vNumDoFs.resize(this->num_subsets(), 0);
		}

		GroupedP1ConformDoFDistribution(GeometricObjectCollection goc,
		                                ISubsetHandler& sh,
		                                storage_manager_type& sm,
		                                FunctionPattern& dp,
		                                const SurfaceView& surfView)
		: base_type(goc, dp, surfView), m_pISubsetHandler(&sh),
		  m_pStorageManager(&sm), m_numDoFs(0)
		{
			m_vNumDoFs.clear();
			m_vNumDoFs.resize(this->num_subsets(), 0);
		}

		///////////////////////////
		// Support Info
		///////////////////////////

		static bool supports_trial_space(LocalShapeFunctionSetID& id)
		{
			return id == LocalShapeFunctionSetID(LocalShapeFunctionSetID::LAGRANGE, 1);
		}

		///////////////////////////
		// Infos
		///////////////////////////

	/// \copydoc IDoFDistribution::num_dofs()
		inline size_t num_dofs() const {return m_numDoFs;}

	/// \copydoc IDoFDistribution::num_dofs(int) const
		inline size_t num_dofs(int si) const {return m_vNumDoFs[si];}

	/// \copydoc IDoFDistribution::blocksize()
		inline int blocksize() const
		{
			if(num_subsets()==0) return -1;
			int blockSize = m_pFuncPattern->num_fct(0);

			for(int si = 1; si < num_subsets(); ++si)
			{
				const int tmpBlockSize = m_pFuncPattern->num_fct(si);
				if(tmpBlockSize != blockSize)
					return -1;
			}
			return blockSize;
		}

	/// \copydoc IDoFDistribution::blocksize(int) const
		inline int blocksize(int si) const
			{return m_pFuncPattern->num_fct(si);}

		///////////////////////////////////////
		// LocalIndex update
		///////////////////////////////////////

	/// \copydoc IDoFDistribution::num_indices()
		size_t num_indices(ReferenceObjectID refID, int si,
		                   const FunctionGroup& funcGroup) const;

	/// \copydoc IDoFDistribution::num_inner_indices()
		size_t num_inner_indices(ReferenceObjectID refID, int si,
		                         const FunctionGroup& funcGroup) const;

	/// \copydoc IDoFDistribution::prepare_indices()
		// TODO: withHanging == true is not yet implemented
		bool prepare_indices(ReferenceObjectID refID, int si,
		                     LocalIndices& ind, bool withHanging = false) const;

	/// \copydoc IDoFDistribution::prepare_inner_indices()
		bool prepare_inner_indices(ReferenceObjectID refID, int si,
		                           LocalIndices& ind) const;

	/// \copydoc IDoFDistribution::update_indices()
		// TODO: withHanging == true is not yet implemented
		template<typename TElem>
		void update_indices(TElem* elem, LocalIndices& ind,
		                    bool withHanging = false) const;

	/// \copydoc IDoFDistribution::update_inner_indices()
		template<typename TElem>
		void update_inner_indices(TElem* elem, LocalIndices& ind) const;

		///////////////////////////////////////
		// Multi index access
		///////////////////////////////////////

	/// \copydoc IDoFDistribution::num_multi_indices()
		template<typename TElem>
		size_t num_multi_indices(TElem* elem, size_t fct) const;

	/// \copydoc IDoFDistribution::num_inner_multi_indices()
		template<typename TElem>
		size_t num_inner_multi_indices(TElem* elem, size_t fct) const;

	/// \copydoc IDoFDistribution::get_multi_indices()
		template<typename TElem>
		size_t get_multi_indices(TElem* elem, size_t fct,
		                         multi_index_vector_type& ind) const;

	/// \copydoc IDoFDistribution::get_inner_multi_indices()
		template<typename TElem>
		size_t get_inner_multi_indices(TElem* elem, size_t fct,
		                               multi_index_vector_type& ind) const;

		///////////////////////////////////////
		// Algebra index access
		///////////////////////////////////////

	/// \copydoc IDoFDistribution::num_algebra_indices()
		template<typename TElem>
		size_t num_algebra_indices(TElem* elem, size_t fct) const;

	/// \copydoc IDoFDistribution::num_inner_algebra_indices()
		template<typename TElem>
		size_t num_inner_algebra_indices(TElem* elem, size_t fct) const;

	/// \copydoc IDoFDistribution::get_algebra_indices()
		template<typename TElem>
		void get_algebra_indices(TElem* elem,
		                         algebra_index_vector_type& ind) const;

	/// \copydoc IDoFDistribution::get_inner_algebra_indices()
		template<typename TElem>
		void get_inner_algebra_indices(TElem* elem,
		                               algebra_index_vector_type& ind) const;

		///////////////////////////
		// Creation
		///////////////////////////

	/// \copydoc IDoFDistribution::distribute_dofs()
		bool distribute_dofs();

	/// \copydoc IDoFDistribution::permute_indices()
		bool permute_indices(const std::vector<size_t>& vIndNew);

	/// \copydoc IDoFDistribution::get_connections()
		bool get_connections(std::vector<std::vector<size_t> >& vvConnection);

	/// \copydoc IDoFDistribution::vertices_created()
		bool vertex_added(VertexBase* vrt);

	/// \copydoc IDoFDistribution::vertices_to_be_erased()
		bool vertex_to_be_removed(VertexBase* vrt);

	/// \copydoc IDoFDistribution::compress()
		bool defragment();

	protected:
	///	returns algebra index attached to a vertex
		size_t& alg_index(VertexBase* vrt, size_t si)
		{
			UG_ASSERT(m_pStorageManager != NULL, "No Storage Manager");
			UG_ASSERT(m_pISubsetHandler != NULL, "No Subset Handler");
			return m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt];
		}

	///	const access to algebra index of a vertex
		const size_t& alg_index(VertexBase* vrt, size_t si) const
		{
			UG_ASSERT(m_pStorageManager != NULL, "No Storage Manager");
			UG_ASSERT(m_pISubsetHandler != NULL, "No Subset Handler");
			return m_pStorageManager->m_vSubsetInfo[si].aaDoFVRT[vrt];
		}

	///	returns the next free index
		size_t get_free_index(size_t si)
		{
		//	The idea is as follows:
		//	- 	If a free index is left, the a free index is returned. This
		//		changes the number of (used) dofs, but the index set remains
		//		the same. (one hole less)
		// 	-	If no free index is left (i.e. no holes in index set and therefore
		//		m_numDoFs == m_sizeIndexSet), the index set is increased and
		//		the newly created index is returned. This changes the size of
		//		the index set and the number of dofs.

		//	strat with default index to be returned
			size_t freeIndex = m_sizeIndexSet;

		//	check if free index available
			if(!m_vFreeIndex.empty())
			{
			//	return free index instead and pop index from free index list
				freeIndex = m_vFreeIndex.back(); m_vFreeIndex.pop_back();
			}
			else
			{
			//	if using new index, increase size of index set
				++m_sizeIndexSet;
			}

		//	adjust counters
			++ m_numDoFs;
			++ (m_vNumDoFs[si]);

		//	return new index
			return freeIndex;
		}

	///	remembers a free index
		void push_free_index(size_t freeIndex, size_t si)
		{
		//	remember index
			m_vFreeIndex.push_back(freeIndex);

		//	decrease number of distributed indices
			-- m_numDoFs;
			-- (m_vNumDoFs[si]);
		}

	protected:
	/// subset handler for this distributor
		ISubsetHandler* m_pISubsetHandler;

	/// Storage Manager for dofs
		storage_manager_type* m_pStorageManager;

	/// number of distributed dofs on whole domain
		size_t m_numDoFs;

	///	number of largest index used, i.e. (0, ..., m_sizeIndexSet-1) available,
	///	but maybe some indices are not used
		size_t m_sizeIndexSet;

	/// number offsetmap
		std::vector<std::vector<size_t> > m_vvOffsets;

	/// number of distributed dofs on each subset
		std::vector<size_t> m_vNumDoFs;

	///	vector to store free algebraic indices
		std::vector<size_t> m_vFreeIndex;
};


} // end namespace ug


#include "./p1conform_impl.h"

#endif /* __H__LIB_DISCRETIZATION__DOF_MANAGER__P1CONFORM__ */
