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

	///	type of attachment for vertex dofs
		typedef storage_manager_type::vertex_attachment_accessor_type
				vertex_attachment_accessor_type;

	public:
		P1ConformDoFDistribution(GeometricObjectCollection goc,
		                         ISubsetHandler& sh, storage_manager_type& sm,
		                         FunctionPattern& fp)
		: base_type(goc, fp), m_pISubsetHandler(&sh),
		  m_pStorageManager(&sm), m_raaVrtDoF(sm.get_vertex_attachment_accessor()),
		  m_numDoFs(0), m_sizeIndexSet(0)
		{
			m_vNumDoFs.clear();
			m_vNumDoFs.resize(this->num_subsets(), 0);

		// 	Attach indices
			if(!m_pStorageManager->update_attachments())
				throw(UGFatalError("Attachment missing in DoF Storage Manager."));

		// 	create offsets
			create_offsets();
		}

		P1ConformDoFDistribution(GeometricObjectCollection goc,
		                         ISubsetHandler& sh, storage_manager_type& sm,
		                         FunctionPattern& fp,
		                         const SurfaceView& surfView)
		: base_type(goc, fp, surfView), m_pISubsetHandler(&sh),
		  m_pStorageManager(&sm), m_raaVrtDoF(sm.get_vertex_attachment_accessor()),
		  m_numDoFs(0), m_sizeIndexSet(0)
		{
			m_vNumDoFs.clear();
			m_vNumDoFs.resize(this->num_subsets(), 0);

		// 	Attach indices
			if(!m_pStorageManager->update_attachments())
				throw(UGFatalError("Attachment missing in DoF Storage Manager."));

		// 	create offsets
			create_offsets();
		}

		///////////////////////////
		// Support Info
		///////////////////////////

		static bool supports_trial_space(LSFSID& id)
		{
			return id == LSFSID(LSFSID::LAGRANGE, 1);
		}

		///////////////////////////
		// Infos
		///////////////////////////

	/// \copydoc ug::IDoFDistribution::has_dofs_on()
		template <typename TElem>
		bool has_dofs_on() const;

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

	/// \copydoc IDoFDistribution::indices()
		template<typename TElem>
		void indices(TElem* elem, LocalIndices& ind, bool bHang = false) const;

		///////////////////////////////////////
		// Algebra Index / Multi index access
		///////////////////////////////////////

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

	/// \copydoc IDoFDistribution::distribute_dofs()
		bool distribute_dofs();

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
		size_t& first_index(VertexBase* vrt, size_t si) {return m_raaVrtDoF[vrt];}

	///	const access to first algebra index of a vertex
		const size_t& first_index(VertexBase* vrt, size_t si) const {return m_raaVrtDoF[vrt];}

	///	returns the next free index
		size_t get_free_index(size_t si);

	///	remembers a free index
		void push_free_index(size_t freeIndex, size_t si);

	protected:
	/// subset handler for this distributor
		ISubsetHandler* m_pISubsetHandler;

	// 	Storage Manager for dofs
		storage_manager_type* m_pStorageManager;

	///	attachment accessor for dof vertices
		vertex_attachment_accessor_type& m_raaVrtDoF;

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

	///	type of attachment for vertex dofs
		typedef storage_manager_type::vertex_attachment_accessor_type
				vertex_attachment_accessor_type;

	public:
		GroupedP1ConformDoFDistribution(GeometricObjectCollection goc,
		                                ISubsetHandler& sh,
		                                storage_manager_type& sm,
		                                FunctionPattern& dp)
		: base_type(goc, dp), m_pISubsetHandler(&sh),
		  m_pStorageManager(&sm), m_raaVrtDoF(sm.get_vertex_attachment_accessor()),
		  m_numDoFs(0), m_sizeIndexSet(0)
		{
			m_vNumDoFs.clear();
			m_vNumDoFs.resize(this->num_subsets(), 0);

		// 	Attach indices
			if(!m_pStorageManager->update_attachments())
				throw(UGFatalError("Attachment missing in DoF Storage Manager."));

		// 	create offsets
			create_offsets();
		}

		GroupedP1ConformDoFDistribution(GeometricObjectCollection goc,
		                                ISubsetHandler& sh,
		                                storage_manager_type& sm,
		                                FunctionPattern& dp,
		                                const SurfaceView& surfView)
		: base_type(goc, dp, surfView), m_pISubsetHandler(&sh),
		  m_pStorageManager(&sm),
		  m_raaVrtDoF(sm.get_vertex_attachment_accessor()),
		  m_numDoFs(0), m_sizeIndexSet(0)
		{
			m_vNumDoFs.clear();
			m_vNumDoFs.resize(this->num_subsets(), 0);

		// 	Attach indices
			if(!m_pStorageManager->update_attachments())
				throw(UGFatalError("Attachment missing in DoF Storage Manager."));

		// 	create offsets
			create_offsets();
		}

		///////////////////////////
		// Support Info
		///////////////////////////

		static bool supports_trial_space(LSFSID& id)
		{
			return id == LSFSID(LSFSID::LAGRANGE, 1);
		}

		///////////////////////////
		// Infos
		///////////////////////////

	/// \copydoc ug::IDoFDistribution::has_dofs_on()
		template <typename TElem>
		bool has_dofs_on() const;

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

	/// \copydoc IDoFDistribution::indices()
		template<typename TElem>
		void indices(TElem* elem, LocalIndices& ind, bool bHang = false) const;

		///////////////////////////////////////
		// Algebra Index / Multi index access
		///////////////////////////////////////

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

	/// \copydoc IDoFDistribution::distribute_dofs()
		bool distribute_dofs();

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

	///	returns algebra index attached to a vertex
		size_t& alg_index(VertexBase* vrt, size_t si) {return m_raaVrtDoF[vrt];}

	///	const access to algebra index of a vertex
		const size_t& alg_index(VertexBase* vrt, size_t si) const {return m_raaVrtDoF[vrt];}

	///	returns the next free index
		size_t get_free_index(size_t si);

	///	remembers a free index
		void push_free_index(size_t freeIndex, size_t si);

	protected:
	/// subset handler for this distributor
		ISubsetHandler* m_pISubsetHandler;

	/// Storage Manager for dofs
		storage_manager_type* m_pStorageManager;

	///	attachment accessor for dof vertices
		vertex_attachment_accessor_type& m_raaVrtDoF;

	/// number of distributed dofs on whole domain
		size_t m_numDoFs;

	///	number of largest index used, i.e. (0, ..., m_sizeIndexSet-1) available,
	///	but maybe some indices are not used
		size_t m_sizeIndexSet;

	/// number of distributed dofs on each subset
		std::vector<size_t> m_vNumDoFs;

	/// number offsetmap
		std::vector<std::vector<size_t> > m_vvOffsets;

	///	vector to store free algebraic indices
		std::vector<size_t> m_vFreeIndex;
};


} // end namespace ug


#include "./p1conform_impl.h"

#endif /* __H__LIB_DISCRETIZATION__DOF_MANAGER__P1CONFORM__ */
