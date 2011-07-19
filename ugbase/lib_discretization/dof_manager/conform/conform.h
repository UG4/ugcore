/*
 * conform.h
 *
 *  Created on: 11.07.2011
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__DOF_MANAGER__CONFORM__
#define __H__LIB_DISCRETIZATION__DOF_MANAGER__CONFORM__

#include <vector>

#include "lib_grid/lg_base.h"

#include "../dof_distribution.h"
#include "../function_pattern.h"
#include "../../reference_element/reference_element.h"

namespace ug{

/// accesses storage for DoFs and handles grid attachment process
class ConformStorageManager
{
	public:
	///	type of DoF attachment
		typedef ug::Attachment<size_t> ADoF;

	///	type of accessor
		typedef Grid::AttachmentAccessor<VertexBase, ADoF>
				vertex_attachment_accessor_type;

	///	type of accessor
		typedef Grid::AttachmentAccessor<EdgeBase, ADoF>
				edge_attachment_accessor_type;

	///	type of accessor
		typedef Grid::AttachmentAccessor<Face, ADoF>
				face_attachment_accessor_type;

	///	type of accessor
		typedef Grid::AttachmentAccessor<Volume, ADoF>
				volume_attachment_accessor_type;

	public:
	///	Constructor
		ConformStorageManager() : m_pSH(NULL), m_pGrid(NULL) {}

	/// set subset handler
		void set_subset_handler(ISubsetHandler& sh);

	/// clear all dofs
		void clear();

	/// destructor
		~ConformStorageManager() {clear();};

	/// attach indices
		bool update_attachments();

	///	returns the associated grid
		Grid* get_assigned_grid() {return m_pGrid;}

	///	returns the underlying subset handler
		ISubsetHandler* get_subset_handler() {return m_pSH;}

	///	returns the attachment accessor
		vertex_attachment_accessor_type& vertex_attachment_accessor()
			{return m_aaIndexVRT;}

	///	returns the attachment accessor
		edge_attachment_accessor_type& edge_attachment_accessor()
			{return m_aaIndexEDGE;}

	///	returns the attachment accessor
		face_attachment_accessor_type& face_attachment_accessor()
			{return m_aaIndexFACE;}

	///	returns the attachment accessor
		volume_attachment_accessor_type& volume_attachment_accessor()
			{return m_aaIndexVOL;}

	protected:
	/// subset handler
		ISubsetHandler* m_pSH;

	///	assosicated grid
		Grid* m_pGrid;

	///	Attachment Accessor
		vertex_attachment_accessor_type m_aaIndexVRT;
		edge_attachment_accessor_type m_aaIndexEDGE;
		face_attachment_accessor_type m_aaIndexFACE;
		volume_attachment_accessor_type m_aaIndexVOL;

	///	Attachment (for vertices)
		ADoF m_aIndex;
};

/// DoF Manager for P1 Functions
/**
 * \tparam	bGrouped	true if dofs are grouped on elements
 */
class DoFDistribution
	: public IDoFDistribution<DoFDistribution>
{
	public:
	///	own type
		typedef DoFDistribution this_type;

	///	Base class
		typedef IDoFDistribution<this_type> base_type;

	/// type of multi index
		typedef MultiIndex<2> index_type;

	/// type of value container for element local indices
		typedef std::vector<index_type> multi_index_vector_type;

	/// type of algebra index vector
		typedef std::vector<size_t> algebra_index_vector_type;

	/// Storage Manager type
		typedef ConformStorageManager storage_manager_type;

	public:
	///	constructor for level DoFDistributions
		DoFDistribution(GeometricObjectCollection goc,
		                ISubsetHandler& sh, storage_manager_type& sm,
		                FunctionPattern& fp);

	///	constructor for surface DoFDistributions (passing surface view)
		DoFDistribution(GeometricObjectCollection goc,
		                  ISubsetHandler& sh, storage_manager_type& sm,
		                  FunctionPattern& fp,
		                  const SurfaceView& surfView);

		///////////////////////////
		// Infos
		///////////////////////////

	///	returns if trial space is supported
		static bool supports_trial_space(LFEID& id)
		{
			return id == LFEID(LFEID::LAGRANGE, 1);
		}

	/// \copydoc ug::IDoFDistribution::has_indices_on( ReferenceObjectID ) const
		bool has_indices_on(ReferenceObjectID roid) const;

	/// \copydoc ug::IDoFDistribution::has_indices_on( GeometricBaseObject ) const
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
		void print_local_dof_statistic(int verboseLev = 1) const;

		///////////////////////////////////////
		// Index Access
		///////////////////////////////////////

	/// \copydoc IDoFDistribution::indices()
		template<typename TElem>
		void indices(TElem* elem, LocalIndices& ind, bool bHang = false) const;

	/// \copydoc IDoFDistribution::multi_indices()
		template<typename TElem>
		size_t multi_indices(TElem* elem, size_t fct,
		                     multi_index_vector_type& ind, bool bClear=true) const;

	/// \copydoc IDoFDistribution::inner_multi_indices()
		template<typename TElem>
		size_t inner_multi_indices(TElem* elem, size_t fct,
		                           multi_index_vector_type& ind, bool bClear=true) const;

	/// \copydoc IDoFDistribution::algebra_indices()
		template<typename TElem>
		size_t algebra_indices(TElem* elem,
		                       algebra_index_vector_type& ind, bool bClear=true) const;

	/// \copydoc IDoFDistribution::inner_algebra_indices()
		template<typename TElem>
		size_t inner_algebra_indices(TElem* elem,
		                             algebra_index_vector_type& ind, bool bClear=true) const;

		///////////////////////////
		// Creation
		///////////////////////////

	/// \copydoc IDoFDistribution::distribute_indices()
		bool distribute_indices();

	/// \copydoc IDoFDistribution::compress()
		bool defragment();

	/// \copydoc IDoFDistribution::permute_indices()
		bool permute_indices(std::vector<size_t>& vIndNew);

	/// \copydoc IDoFDistribution::get_connections()
		bool get_connections(std::vector<std::vector<size_t> >& vvConnection);

	/// \copydoc IDoFDistribution::vertices_created()
		void grid_obj_added(VertexBase* vrt);
		void grid_obj_added(EdgeBase* edge);
		void grid_obj_added(Face* face);
		void grid_obj_added(Volume* vol);

	/// \copydoc IDoFDistribution::vertices_to_be_erased()
		void grid_obj_to_be_removed(VertexBase* vrt);
		void grid_obj_to_be_removed(EdgeBase* edge);
		void grid_obj_to_be_removed(Face* face);
		void grid_obj_to_be_removed(Volume* vol);

	/// \copydoc IDoFDistribution::grid_obj_replaced()
		void grid_obj_replaced(VertexBase* vrtNew, VertexBase* vrtOld);
		void grid_obj_replaced(EdgeBase* edgeNew, EdgeBase* edgeOld);
		void grid_obj_replaced(Face* faceNew, Face* faceOld);
		void grid_obj_replaced(Volume* volNew, Volume* volOld);

	protected:
		template<typename TElem, typename TBaseElem>
		void indices(TElem* elem, LocalIndices& ind,std::vector<TBaseElem*> vElem,
		             size_t numNatural, bool bHang) const;

		template<typename TBaseElem>
		size_t multi_indices(std::vector<TBaseElem*> vElem, size_t fct,
		                     multi_index_vector_type& ind) const;

		size_t inner_multi_indices(multi_index_vector_type& ind,
		                           const size_t firstIndex, const int si,
		                           const size_t fct, const ReferenceObjectID type) const;

		template<typename TBaseElem>
		size_t algebra_indices(std::vector<TBaseElem*> vElem,
		                       algebra_index_vector_type& ind) const;

		size_t	inner_algebra_indices(algebra_index_vector_type& ind,
		      	                      const size_t firstIndex, const int si,
		      	                      const ReferenceObjectID type) const;

	///	distributes the dofs on an element type
		template <typename TElem>
		bool distribute_indices();

	///	permutes indices on an element type
		template <typename TElem>
		bool permute_indices(std::vector<size_t>& vIndNew);

		template <typename TElem, typename TBaseElem>
		bool add_connections_for_adjacent(std::vector<std::vector<size_t> >& vvConnection,
		                                  TElem* elem, std::vector<TBaseElem*> vConnElem);

	///	fills the connections of an element type
		template <typename TBaseElem>
		bool get_connections(std::vector<std::vector<size_t> >& vvConnection);

	///	replaces the indices on an element type
		template <typename TElem>
		bool defragment(std::vector<std::pair<size_t, size_t> >& vReplaced);

	///	callback when an element is added to grid for an element type
		template <typename TBaseElem>
		void grid_obj_added(TBaseElem* elem);

	///	callback when an element is removed from grid for an element type
		template <typename TBaseElem>
		void grid_obj_to_be_removed(TBaseElem* elem);

	///	callback when an element is replaced in grid for an element type
		template <typename TBaseElem>
		void grid_obj_replaced(TBaseElem* elemNew, TBaseElem* elemOld);

	///	creates the offset array
		void create_offsets();

	///	creates the offset array for a reference element type
		void create_offsets(ReferenceObjectID roid);

	///	returns first algebra index of a vertex
		size_t& first_index(VertexBase* vrt) {return m_raaIndexVRT[vrt];}
		size_t& first_index(EdgeBase* ed) 	 {return m_raaIndexEDGE[ed];}
		size_t& first_index(Face* face)      {return m_raaIndexFACE[face];}
		size_t& first_index(Volume* vol)     {return m_raaIndexVOL[vol];}

	///	const access to first algebra index of a vertex
		const size_t& first_index(VertexBase* vrt) const {return m_raaIndexVRT[vrt];}
		const size_t& first_index(EdgeBase* ed)    const {return m_raaIndexEDGE[ed];}
		const size_t& first_index(Face* face)      const {return m_raaIndexFACE[face];}
		const size_t& first_index(Volume* vol)     const {return m_raaIndexVOL[vol];}

	///	returns the next free index
		size_t get_free_index(size_t si, ReferenceObjectID roid);

	///	remembers a free index
		void push_free_index(size_t freeIndex, size_t si, ReferenceObjectID roid);

	protected:
	/// subset handler for this distributor
		ISubsetHandler* m_pISubsetHandler;

	// 	Storage Manager for dofs
		storage_manager_type* m_pStorageManager;

	//	type of attachment for indices
		typedef storage_manager_type::vertex_attachment_accessor_type
				vertex_attachment_accessor_type;
		typedef storage_manager_type::edge_attachment_accessor_type
				edge_attachment_accessor_type;
		typedef storage_manager_type::face_attachment_accessor_type
				face_attachment_accessor_type;
		typedef storage_manager_type::volume_attachment_accessor_type
				volume_attachment_accessor_type;

	//	attachment accessor for vertices
		vertex_attachment_accessor_type& m_raaIndexVRT;
		edge_attachment_accessor_type& m_raaIndexEDGE;
		face_attachment_accessor_type& m_raaIndexFACE;
		volume_attachment_accessor_type& m_raaIndexVOL;

	/// number of distributed indices on whole domain
		size_t m_numIndex;

	/// number of distributed indices on each subset
		std::vector<size_t> m_vNumIndex;

	///	number of largest index used, i.e. (0, ..., m_sizeIndexSet-1) available,
	///	but maybe some indices are not used
		size_t m_sizeIndexSet;

	///	vector to store free algebraic indices
		std::vector<size_t> m_vFreeIndex;

	///	indication that function is not defined on a subset
		enum{NOT_DEF_ON_SUBSET = (size_t) -1};

	/// offset map
		std::vector<std::vector<size_t> > m_vvvOffsets[NUM_REFERENCE_OBJECTS];

	///	number of DoFs on a reference element type on a subset
		std::vector<size_t> m_vvNumDoFsOnROID[NUM_REFERENCE_OBJECTS];

	///	maximum number of DoFs on a reference type
		size_t m_vMaxDoFsOnROID[NUM_REFERENCE_OBJECTS];

	///	maximum number of DoFs on geometric objects in a dimension
		size_t m_vMaxDoFsInDim[4];

	///	flag if parameters are grouped
		bool m_bGrouped;
};

} // end namespace ug


#include "conform_impl.h"

#endif /* __H__LIB_DISCRETIZATION__DOF_MANAGER__CONFORM__ */
