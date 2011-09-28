/*
 * p1conform.h
 *
 *  Created on: 13.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__DOF_MANAGER__P1CONFORM__
#define __H__UG__LIB_DISC__DOF_MANAGER__P1CONFORM__

#include <vector>

#include "lib_grid/lg_base.h"

#include "../dof_distribution.h"
#include "../function_pattern.h"
#include "../../reference_element/reference_element.h"
#include "../dof_distribution_type.h"

namespace ug{

/// accesses storage for DoFs and handles grid attachment process
class P1StorageManager
{
	public:
	///	type of DoF attachment
		typedef ug::Attachment<size_t> ADoF;

	///	type of accessor
		typedef Grid::AttachmentAccessor<VertexBase, ADoF>
				vertex_attachment_accessor_type;

	public:
	///	Constructor
		P1StorageManager() : m_pSH(NULL), m_pGrid(NULL) {}

	/// set subset handler
		void set_subset_handler(ISubsetHandler& sh);

	/// clear all dofs
		void clear();

	/// destructor
		~P1StorageManager() {clear();};

	/// attach dofs
		bool update_attachments();

	///	returns the associated grid
		Grid* get_assigned_grid() {return m_pGrid;}

	///	returns the underlying subset handler
		ISubsetHandler* get_subset_handler() {return m_pSH;}

	///	returns the attachment accessor
		vertex_attachment_accessor_type& get_vertex_attachment_accessor()
			{return m_aaDoFVRT;}

	protected:
	/// subset handler
		ISubsetHandler* m_pSH;

	///	assosicated grid
		Grid* m_pGrid;

	///	Attachment Accessor
		vertex_attachment_accessor_type m_aaDoFVRT;

	///	Attachment (for vertices)
		ADoF m_aDoF;
};

/// DoF Manager for P1 Functions
class P1DoFDistribution
	: public IDoFDistribution<P1DoFDistribution>
{
	public:
	///	dof distribution type
		static const DofDistributionType type = DDT_P1CONFORM;

	///	own type
		typedef P1DoFDistribution this_type;

	///	Base class
		typedef IDoFDistribution<this_type> base_type;

	/// type of multi index
		typedef MultiIndex<2> index_type;

	/// type of value container for element local indices
		typedef std::vector<index_type> multi_index_vector_type;

	/// type of algebra index vector
		typedef std::vector<size_t> algebra_index_vector_type;

	/// Storage Manager type
		typedef P1StorageManager storage_manager_type;

	///	type of attachment for vertex dofs
		typedef storage_manager_type::vertex_attachment_accessor_type
				vertex_attachment_accessor_type;

		using base_type::num_fct;
		using base_type::num_subsets;
		using base_type::is_def_in_subset;
		using base_type::indices_swaped;
		using base_type::num_indices_changed;

	public:
		P1DoFDistribution(GeometricObjectCollection goc,
		                  ISubsetHandler& sh, storage_manager_type& sm,
		                  FunctionPattern& fp)
		: base_type(goc, fp), m_bGrouped(false), m_pISubsetHandler(&sh),
		  m_pStorageManager(&sm), m_raaVrtDoF(sm.get_vertex_attachment_accessor()),
		  m_numIndex(0), m_sizeIndexSet(0)
		{
			m_vNumIndex.clear();
			m_vNumIndex.resize(this->num_subsets(), 0);

		// 	Attach indices
			if(!m_pStorageManager->update_attachments())
				throw(UGFatalError("Attachment missing in DoF Storage Manager."));

		// 	create offsets
			create_offsets();
		}

		P1DoFDistribution(GeometricObjectCollection goc,
		                  ISubsetHandler& sh, storage_manager_type& sm,
		                  FunctionPattern& fp,
		                  const SurfaceView& surfView)
		: base_type(goc, fp, surfView), m_bGrouped(false), m_pISubsetHandler(&sh),
		  m_pStorageManager(&sm), m_raaVrtDoF(sm.get_vertex_attachment_accessor()),
		  m_numIndex(0), m_sizeIndexSet(0)
		{
			m_vNumIndex.clear();
			m_vNumIndex.resize(this->num_subsets(), 0);

		// 	Attach indices
			if(!m_pStorageManager->update_attachments())
				throw(UGFatalError("Attachment missing in DoF Storage Manager."));

		// 	create offsets
			create_offsets();
		}

		///////////////////////////
		// Infos
		///////////////////////////

	///	set grouped
		void set_grouping(bool bGrouped)
		{
			if(num_indices() > 0)
				throw(UGFatalError("P1DoFDistribution: Grouping flag can not be "
						"changed after initialization."));
			m_bGrouped = bGrouped;
		}

	///	returns if trial space is supported
		static bool supports_trial_space(LFEID& id)
		{
			return id == LFEID(LFEID::LAGRANGE, 1);
		}

	/// \copydoc ug::IDoFDistribution::has_indices_on()
		bool has_indices_on(ReferenceObjectID roid) const;
		bool has_indices_on(GeometricBaseObject gbo) const;

	/// \copydoc ug::IDoFDistribution::num_indices()
		size_t num_indices() const {return m_numIndex;}

	/// \copydoc ug::IDoFDistribution::num_indices(int) const
		size_t num_indices(int si) const {return m_vNumIndex[si];}

	/// \copydoc IDoFDistribution::blocksize()
		int blocksize() const
		{
			if(!m_bGrouped) return 1;
			else{
				if(this->num_subsets()==0) return -1;
				int blockSize = this->m_pFuncPattern->num_fct(0);

				for(int si = 1; si < this->num_subsets(); ++si)
				{
					const int tmpBlockSize = this->m_pFuncPattern->num_fct(si);
					if(tmpBlockSize != blockSize)
						return -1;
				}
				return blockSize;
			}
		}

	/// \copydoc IDoFDistribution::blocksize(int) const
		int blocksize(int si) const
		{
			if(!m_bGrouped) return 1;
			else return this->m_pFuncPattern->num_fct(si);
		}

	/// \copydoc IDoFDistribution::print_local_dof_statistic(int) const
		void print_local_dof_statistic(int verboseLev = 1) const
		{
			UG_LOG("P1ConformDoFManager: print_local_dof_statistic not implemented.\n");
		}

		///////////////////////////////////////
		// Index Access
		///////////////////////////////////////

	/// \copydoc IDoFDistribution::indices()
		template<typename TElem>
		void indices(TElem* elem, LocalIndices& ind, bool bHang = false) const;

	/// \copydoc IDoFDistribution::multi_indices()
		template<typename TElem>
		size_t multi_indices(TElem* elem, size_t fct,
		                         multi_index_vector_type& ind) const;

	/// \copydoc IDoFDistribution::inner_multi_indices()
		template<typename TElem>
		size_t inner_multi_indices(TElem* elem, size_t fct,
		                               multi_index_vector_type& ind) const;

	/// \copydoc IDoFDistribution::algebra_indices()
		template<typename TElem>
		size_t algebra_indices(TElem* elem,
		                         algebra_index_vector_type& ind) const;

	/// \copydoc IDoFDistribution::inner_algebra_indices()
		template<typename TElem>
		size_t inner_algebra_indices(TElem* elem,
		                               algebra_index_vector_type& ind) const;

		///////////////////////////
		// Creation
		///////////////////////////

	/// \copydoc IDoFDistribution::distribute_indices()
		bool distribute_indices();

	/// \copydoc IDoFDistribution::permute_indices()
		bool permute_indices(std::vector<size_t>& vIndNew);

	/// \copydoc IDoFDistribution::get_connections()
		bool get_connections(std::vector<std::vector<size_t> >& vvConnection);

	/// \copydoc IDoFDistribution::vertices_created()
		void grid_obj_added(VertexBase* vrt);
		void grid_obj_added(EdgeBase* edge) {}
		void grid_obj_added(Face* face) {}
		void grid_obj_added(Volume* vol) {}

	/// \copydoc IDoFDistribution::vertices_to_be_erased()
		void grid_obj_to_be_removed(VertexBase* vrt);
		void grid_obj_to_be_removed(EdgeBase* edge) {}
		void grid_obj_to_be_removed(Face* face) {}
		void grid_obj_to_be_removed(Volume* vol) {}

	/// \copydoc IDoFDistribution::grid_obj_replaced()
		void grid_obj_replaced(VertexBase* vrtNew, VertexBase* vrtOld);
		void grid_obj_replaced(EdgeBase* edgeNew, EdgeBase* edgeOld) {}
		void grid_obj_replaced(Face* faceNew, Face* faceOld) 	{}
		void grid_obj_replaced(Volume* volNew, Volume* volOld) 	{}

	/// \copydoc IDoFDistribution::compress()
		bool defragment();

	protected:
	///	creates the offset array
		void create_offsets();

	///	returns first algebra index of a vertex
		size_t& first_index(VertexBase* vrt) {return m_raaVrtDoF[vrt];}

	///	const access to first algebra index of a vertex
		const size_t& first_index(VertexBase* vrt) const {return m_raaVrtDoF[vrt];}

	///	returns the next free index
		size_t get_free_index(size_t si);

	///	remembers a free index
		void push_free_index(size_t freeIndex, size_t si);

	protected:
	///	flag to indicate if dofs should be grouped
		bool m_bGrouped;

	/// subset handler for this distributor
		ISubsetHandler* m_pISubsetHandler;

	// 	Storage Manager for dofs
		storage_manager_type* m_pStorageManager;

	///	attachment accessor for dof vertices
		vertex_attachment_accessor_type& m_raaVrtDoF;

	/// number of distributed dofs on whole domain
		size_t m_numIndex;

	///	number of largest index used, i.e. (0, ..., m_sizeIndexSet-1) available,
	///	but maybe some indices are not used
		size_t m_sizeIndexSet;

	/// number offsetmap
		std::vector<std::vector<size_t> > m_vvOffsets;

	/// number of distributed dofs on each subset
		std::vector<size_t> m_vNumIndex;

	///	vector to store free algebraic indices
		std::vector<size_t> m_vFreeIndex;
};

} // end namespace ug

#include "./p1conform_impl.h"

#endif /* __H__UG__LIB_DISC__DOF_MANAGER__P1CONFORM__ */
