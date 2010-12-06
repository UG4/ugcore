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
#include "../../common/common.h"
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

class P1ConformFunctionPattern : public FunctionPattern
{
	public:
		virtual bool add_discrete_function(const char* name, LocalShapeFunctionSetID id, int dim);

		virtual bool add_discrete_function(const char* name, LocalShapeFunctionSetID id, const SubsetGroup& SubsetIndices, int dim);
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

	/// ordering
		bool order_cuthill_mckee(bool bReverse = false);

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
			int blockSize = m_pFunctionPattern->num_fct(0);

			for(int si = 1; si < num_subsets(); ++si)
			{
				const int tmpBlockSize = m_pFunctionPattern->num_fct(si);
				if(tmpBlockSize != blockSize)
					return -1;
			}
			return blockSize;
		}

	/// \copydoc IDoFDistribution::blocksize(int) const
		inline int blocksize(int si) const
			{return m_pFunctionPattern->num_fct(si);}

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

	/// ordering
		bool order_cuthill_mckee(bool bReverse = false);

	protected:
	/// subset handler for this distributor
		ISubsetHandler* m_pISubsetHandler;

	/// Storage Manager for dofs
		storage_manager_type* m_pStorageManager;

	/// number of distributed dofs on whole domain
		size_t m_numDoFs;

	/// number offsetmap
		std::vector<std::vector<size_t> > m_vvOffsets;

	/// number of distributed dofs on each subset
		std::vector<size_t> m_vNumDoFs;
};


} // end namespace ug


#include "./p1conform_impl.h"

#endif /* __H__LIB_DISCRETIZATION__DOF_MANAGER__P1CONFORM__ */
