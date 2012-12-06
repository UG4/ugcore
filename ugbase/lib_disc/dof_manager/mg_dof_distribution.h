/*
 * mg_dof_distribution.h
 *
 *  Created on: 29.11.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__DOF_MANAGER__MG_DOF_DISTRIBUTION__
#define __H__UG__LIB_DISC__DOF_MANAGER__MG_DOF_DISTRIBUTION__

#include "lib_grid/tools/surface_view.h"
#include "function_pattern.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/local_finite_element/local_finite_element_id.h"
#include "lib_disc/local_finite_element/local_dof_set.h"
#include "lib_disc/common/local_algebra.h"
#include "lib_disc/dof_manager/grid_level.h"

#ifdef UG_PARALLEL
#include "pcl/pcl_base.h"
#include "lib_algebra/parallelization/parallel_index_layout.h"
#endif

namespace ug{

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Lev Info
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

struct LevInfoBase
{
///	constructor
	LevInfoBase() : numIndex(0), sizeIndexSet(0) {}

/// number of distributed indices on whole domain
	size_t numIndex;

/// number of distributed indices on each subset
	std::vector<size_t> vNumIndexOnSubset;

///	number of largest index used, i.e. (0, ..., m_sizeIndexSet-1) available,
///	but maybe some indices are not used
	size_t sizeIndexSet;

#ifdef UG_PARALLEL
	public:
///	(horizontal) master index layout
	IndexLayout masterLayout;

///	(horizontal) slave index layout
	IndexLayout slaveLayout;

///	vertical master index layout
	IndexLayout verticalMasterLayout;

///	vertical slave index layout
	IndexLayout verticalSlaveLayout;

///	process communicator
	pcl::ProcessCommunicator processCommunicator;

///	communicator
	pcl::InterfaceCommunicator<IndexLayout> communicator;
#endif
};

template <typename TContainer = std::vector<size_t> >
struct LevInfo : public LevInfoBase
{
///	returns if free index available
	inline bool free_index_available() const;

///	returns number of free index available
	inline size_t num_free_index() const;

///	returns a free index
	inline size_t pop_free_index();

///	adds a free index, returns if index has not been contained before
	inline bool push_free_index(size_t index);
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// MGDoFDistribution
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

class MGDoFDistribution : public GridObserver
{
	public:
	//	type of multi index
		typedef MultiIndex<2> multi_index_type;

	public:
		MGDoFDistribution(SmartPtr<MultiGrid> spMG,
						  SmartPtr<MGSubsetHandler> spMGSH,
		                  FunctionPattern& fctPatt,
		                  bool bGrouped);

		~MGDoFDistribution();

		///	returns the multigrid
		SmartPtr<MultiGrid> multi_grid() {return m_spMG;}
		ConstSmartPtr<MultiGrid> multi_grid() const {return m_spMG;}

		///	returns the subset handler
		ConstSmartPtr<ISubsetHandler> subset_handler() const {return m_spMGSH;}

		///	returns function pattern
		const FunctionPattern& function_pattern() const {return m_rFctPatt;}

		/// number of discrete functions on subset si
		size_t num_fct() const {return m_rFctPatt.num_fct();}

		/// number of discrete functions on subset si
		size_t num_fct(int si) const {return m_rFctPatt.num_fct(si);}

		/// returns the name of the discrete function nr_fct
		std::string name(size_t fct) const {return m_rFctPatt.name(fct);}

		/// returns the dimension in which solution lives
		int dim(size_t fct) const {return m_rFctPatt.dim(fct);}

		///	returns dimension of subset
		int dim_subset(int si) const {return m_rFctPatt.dim_subset(si);}

		///	returns subset name
		std::string subset_name(int si) const {return m_rFctPatt.subset_name(si);}

		///	returns the local finite element id of a function
		const LFEID& local_finite_element_id(size_t fct) const {return m_vLFEID[fct];}

		///	returns subset group by name
		SubsetGroup subset_grp_by_name(const char* names) const;

		/// returns fct id by name
		size_t fct_id_by_name(const char* name) const{return m_rFctPatt.fct_id_by_name(name);}

		///	returns a function group to a string of functions
		FunctionGroup fct_grp_by_name(const char* names) const {return m_rFctPatt.fct_grp_by_name(names);}

		///	retruns if a function is defined on a subset
		//	\todo cache
		bool is_def_in_subset(size_t fct, int si) const {return m_rFctPatt.is_def_in_subset(fct, si);}

		/// returns true if the discrete function nr_fct is defined everywhere
		bool is_def_everywhere(size_t fct) const {return m_rFctPatt.is_def_everywhere(fct);}

		///	returns if dofs are grouped
		bool grouped() const {return m_bGrouped;}

		///	returns blocksize
		std::string blocksize() const {return "?";}

		/// returns the maximum number of dofs on grid objects in a dimension on a subset
		size_t max_dofs(const int dim, const int si) const {return m_vvMaxDoFsInDimPerSubset[dim][si];}

		/// return the maximum number of dofs on grid objects in a dimension
		size_t max_dofs(const int dim) const {return m_vMaxDoFsInDim[dim];}

		///	returns the maximum number of dofs on a Reference Object
		size_t num_dofs(const ReferenceObjectID roid, const int si) const {return m_vvNumDoFsOnROID[roid][si];}

		///	returns the number of dofs on a subelement of an element
		size_t num_dofs(size_t fct, const ReferenceObjectID roid, const ReferenceObjectID subRoid) const {return m_vNumDoFOnSubelem[fct](roid, subRoid);}

		///	returns number of subsets
		inline int num_subsets() const {return m_spMGSH->num_subsets();}

		/// returns number of levels
		inline int num_levels() const {return m_spMGSH->num_levels();}

		///	returns if indices are defined on a geometric object
		bool has_indices_on(GeometricBaseObject gbo) const {return m_vMaxDoFsInDim[gbo] > 0;}

		///	returns if indices are defined on a reference object
		bool has_indices_on(ReferenceObjectID roid) const {return m_vMaxDoFsOnROID[roid] > 0;}

		///	prints informations
		void print_local_dof_statistic(int verboseLev) const;

		/// writes the local finite element ids to LocalIndeces
		void local_finite_element_ids(LocalIndices& ind) const;

		/// extracts all indices of the element (sorted)
		/**
		 * All Indices of the element (including the subelements) are extracted
		 * and stored in the LocalIndices structure. The order of the indices
		 * is sorted, i.e. the dofs are provided as specified in the local
		 * dof set of the local finite element trial space.
		 * If bHang is set to true, also the DoFs on the Constained Objects
		 * belonging to the constraining Subelements are extracted and added at
		 * the end of the indices.
		 *
		 * \param[in]		elem		the element
		 * \param[out]		ind			Local indices
		 * \param[in]		bHang		flag if extracting of constrained dofs required
		 */
		template <typename TBaseElem>
		void indices(TBaseElem* elem, LocalIndices& ind, bool bHang = false) const;
		void indices(GeometricObject* elem, LocalIndices& ind, bool bHang = false) const;

		/// extracts all multiindices for a function (sorted)
		/**
		 * All Multi-Indices of a function living on the element (including the
		 * subelements) are extracted and stored in a std::vector. The order of
		 * the indices is sorted, i.e. the dofs are provided as specified in the
		 * local dof set of the local finite element trial space.
		 * If bHang is set to true, also the DoFs on the Constrained Objects
		 * belonging to the constraining Subelements are extracted and added at
		 * the end of the indices.
		 * If bClear is set to true, the vector is cleared before insertion.
		 *
		 * \param[in]		elem		the element
		 * \param[in]		fct			the function
		 * \param[out]		ind			vector of multi indices
		 * \param[in]		bHang		flag if extracting of constrained dofs required
		 * \param[in]		bClear		flag if vector has to be clear before insertion
		 */
		template<typename TBaseElem>
		size_t multi_indices(TBaseElem* elem, size_t fct,
		                     std::vector<multi_index_type>& ind,
		                     bool bHang = false, bool bClear = true) const;
		size_t multi_indices(GeometricObject* elem, size_t fct,
		                     std::vector<multi_index_type>& ind,
		                     bool bHang = false, bool bClear = true) const;

		/// extracts all multiindices of a function in the inner (sorted)
		/**
		 * All Multi-Indices of a function living on the element (including the
		 * subelements) are extracted and stored in a std::vector. The order of
		 * the indices is sorted, i.e. the dofs are provided as specified in the
		 * local dof set of the local finite element trial space.
		 * If bClear is set to true, the vector is cleared before insertion.
		 *
		 * \param[in]		elem		the element
		 * \param[in]		fct			the function
		 * \param[out]		ind			vector of multi indices
		 * \param[in]		bClear		flag if vector has to be clear before insertion
		 */
		template<typename TBaseElem>
		size_t inner_multi_indices(TBaseElem* elem, size_t fct,
		                           std::vector<multi_index_type>& ind,
		                           bool bClear = true) const;
		size_t inner_multi_indices(GeometricObject* elem, size_t fct,
		                           std::vector<multi_index_type>& ind,
		                           bool bClear = true) const;

		/// extracts all algebra indices of an element (not sorted)
		/**
		 * All Algebra-Indices of the element (including the subelements) are
		 * extracted and stored in a std::vector. The order of the indices is
		 * not sorted, no constrained DoFs are extracted.
		 * If bClear is set to true, the vector is cleared before insertion.
		 *
		 * \param[in]		elem		the element
		 * \param[out]		ind			vector of algebra indices
		 * \param[in]		bClear		flag if vector has to be clear before insertion
		 */
		template<typename TBaseElem>
		size_t algebra_indices(TBaseElem* elem,	std::vector<size_t>& ind,
		                       bool bClear = true) const;
		size_t algebra_indices(GeometricObject* elem,	std::vector<size_t>& ind,
		                       bool bClear = true) const;

		/// extracts all algebra indices in the inner of the element (not sorted)
		/**
		 * All Algebra-Indices of the element (excluding the subelements) are
		 * extracted and stored in a std::vector. The order of the indices is
		 * not sorted, no constrained DoFs are extracted.
		 * If bClear is set to true, the vector is cleared before insertion.
		 *
		 * \param[in]		elem		the element
		 * \param[out]		ind			vector of algebra indices
		 * \param[in]		bClear		flag if vector has to be clear before insertion
		 */
		template<typename TBaseElem>
		size_t inner_algebra_indices(TBaseElem* elem, std::vector<size_t>& ind,
		                             bool bClear = true) const;
		size_t inner_algebra_indices(GeometricObject* elem, std::vector<size_t>& ind,
		                             bool bClear = true) const;


	///	returns indices, that can be changed on the element
		template <typename TBaseElem>
		void changable_indices(std::vector<size_t>& vIndex,
		                       const std::vector<TBaseElem*>& vElem) const;

	///	returns parent != NULL, if available
		template <typename TBaseElem>
		GeometricObject* get_parent(TBaseElem* elem) const;

	///	returns parent != NULL, if is copy in sense of Multigrid
		template <typename TBaseElem>
		TBaseElem* parent_if_copy(TBaseElem* elem) const;

	///	returns parent != NULL, if of same base object type
		template <typename TBaseElem>
		TBaseElem* parent_if_same_type(TBaseElem* elem) const;

	///	returns child != NULL, if is copy in sense of Multigrid
		template <typename TBaseElem>
		TBaseElem* child_if_copy(TBaseElem* elem) const;

	protected:
		///	returns the offset for reference element, subset and function
		size_t offset(const ReferenceObjectID roid, const int si, const size_t fct) const {return m_vvvOffsets[roid][si][fct];}

		///	extracts the indices of the vertices
		template<typename TBaseElem>
		void indices_on_vertex(TBaseElem* elem, const ReferenceObjectID roid,
		                       LocalIndices& ind,
		                       const Grid::SecureVertexContainer& vElem) const;

		///	extract dofs on constrained objects
		template <typename TConstraining, typename TConstrained, typename TBaseElem>
		void constrained_indices(LocalIndices& ind,
		                         const typename Grid::traits<TBaseElem>::secure_container& vSubElem) const;

		/// extracts the indices of the subelement of an element
		template<typename TBaseElem, typename TSubBaseElem>
		void indices(TBaseElem* elem, const ReferenceObjectID roid,
		             LocalIndices& ind,
		             const typename Grid::traits<TSubBaseElem>::secure_container& vElem,
		             const Grid::SecureVertexContainer& vCorner) const;

		/// extracts the indices of a subelement of an element
		template<typename TBaseElem, typename TSubBaseElem>
		void multi_indices(TBaseElem* elem, const ReferenceObjectID roid,
		                   size_t fct, std::vector<multi_index_type>& ind,
		                   const typename Grid::traits<TSubBaseElem>::secure_container& vElem,
		                   const Grid::SecureVertexContainer& vCorner, bool bHang) const;


		/// adds all algebra indices of an geom object to the LocalIndices
		template <typename TBaseElem>
		size_t extract_inner_algebra_indices(TBaseElem* elem,
		                                     std::vector<size_t>& ind) const;

		///	adds all algebra indices of a set of geometric objects
		template<typename TBaseElem>
		void extract_inner_algebra_indices(const typename Grid::traits<TBaseElem>::secure_container& vElem,
		                                   std::vector<size_t>& ind) const;

	protected:
		/**
		 * adds indices to a geometric object.
		 *
		 * \param[in]		obj			Geometric Object
		 * \param[in]		roid		Reference Object id
		 * \param[in]		si			Subset of Geometric Object
		 * \param[in,out]	li			Level Information about Indices
		 */
		template <typename TBaseObject, typename T>
		void add(TBaseObject* obj, const ReferenceObjectID roid,
		         const int si, LevInfo<T>& li);
		template <typename T>
		void add(GeometricObject* obj, const ReferenceObjectID roid,
		         const int si, LevInfo<T>& li);

		/**
		 * adds indices to a geometric object. Tries to reuse free indices.
		 *
		 * \param[in]		obj			Geometric Object
		 * \param[in]		roid		Reference Object id
		 * \param[in]		si			Subset of Geometric Object
		 * \param[in,out]	li			Level Information about Indices
		 */
		template <typename TBaseObject, typename T>
		void add_from_free(TBaseObject* obj, const ReferenceObjectID roid,
		                   const int si, LevInfo<T>& li);
		template <typename T>
		void add_from_free(GeometricObject* obj, const ReferenceObjectID roid,
		                   const int si, LevInfo<T>& li);

		/**
		 * removes indices from the geometric object. The freed index (the
		 * produced hole in the index set) is stored in vFreeIndex
		 *
		 * \param[in]		obj					Geometric Object
		 * \param[in]		roid				Reference Object id
		 * \param[in]		si					Subset of Geometric Object
		 * \param[in,out]	li			Level Information about Indices
		 */
		template <typename TBaseObject, typename T>
		void erase(TBaseObject* obj, const ReferenceObjectID roid,
		           const int si, LevInfo<T>& li);
		template <typename T>
		void erase(GeometricObject* obj, const ReferenceObjectID roid,
		           const int si, LevInfo<T>& li);

		/**
		 * checks if the indices on the geometric object are in the range of the
		 * index set. If this is not the case, the index is replaced by a hole
		 * in the index set. The replacement of indices is remembered in the
		 * vReplaced vector by adding the pair (oldIndex, newIndex) at the end.
		 *
		 * \param[in]		obj			Geometric Object
		 * \param[in]		roid		Reference Object id
		 * \param[in]		si			Subset of Geometric Object
		 * \param[in,out]	li			Level Information about Indices
		 * \param[in,out]	vReplaced	pairs of replaced indices
		 *
		 * \returns true if index has been replaced, false else
		 */
		template <typename TBaseObject, typename T>
		void defragment(TBaseObject* obj, const ReferenceObjectID roid, const int si,
		                LevInfo<T>& li, std::vector<std::pair<size_t, size_t> >& vReplaced);
		template <typename T>
		void defragment(GeometricObject* obj, const ReferenceObjectID roid, const int si,
		                LevInfo<T>& li, std::vector<std::pair<size_t, size_t> >& vReplaced);

		/**
		 * copies the indices from one geometric object to another one. This
		 * method is usually called, if an object of the grid is replaced by a
		 * similar one.
		 */
		inline void copy(GeometricObject* objNew, GeometricObject* objOld);

		///	checks that subset assigment is ok
		void check_subsets();

		/// creates offset arrays
		void create_offsets(ReferenceObjectID roid);
		void create_offsets();

		/// initializes the attachments
		void init_attachments();

		/// removes the attachments
		void clear_attachments();

		///	registers this class as a grid observer
		void register_observer();

		///	unregisters this class as a grid observer
		void unregister_observer();

		///	returns first algebra index of a geometric object
		/// \{
		inline size_t& obj_index(GeometricObject* obj);
		inline size_t& obj_index(VertexBase* vrt) 	{return m_aaIndexVRT[vrt];}
		inline size_t& obj_index(EdgeBase* ed) 		{return m_aaIndexEDGE[ed];}
		inline size_t& obj_index(Face* face)     	{return m_aaIndexFACE[face];}
		inline size_t& obj_index(Volume* vol)     	{return m_aaIndexVOL[vol];}
		/// \}

		///	const access to first algebra index of a geometric object
		/// \{
		inline const size_t& obj_index(GeometricObject* obj) const;
		inline const size_t& obj_index(VertexBase* vrt) const {return m_aaIndexVRT[vrt];}
		inline const size_t& obj_index(EdgeBase* ed)    const {return m_aaIndexEDGE[ed];}
		inline const size_t& obj_index(Face* face)      const {return m_aaIndexFACE[face];}
		inline const size_t& obj_index(Volume* vol)     const {return m_aaIndexVOL[vol];}
		/// \}

	protected:
	///	grouping
		bool m_bGrouped;

	///	indication that function is not defined on a subset
		enum{NOT_DEF_ON_SUBSET = (size_t) -1};

	/// offset map
		std::vector<std::vector<size_t> > m_vvvOffsets[NUM_REFERENCE_OBJECTS];

	///	number of DoFs on a reference element type on a subset
		std::vector<size_t> m_vvNumDoFsOnROID[NUM_REFERENCE_OBJECTS];

	///	maximum number of DoFs on a reference type
		size_t m_vMaxDoFsOnROID[NUM_REFERENCE_OBJECTS];

	///	maximum number of DoFs on geometric objects in a dimension
		size_t m_vMaxDoFsInDim[NUM_GEOMETRIC_BASE_OBJECTS];

	///	maximum number of DoFs on geometric objects in a dimension per subset
		std::vector<size_t> m_vvMaxDoFsInDimPerSubset[NUM_GEOMETRIC_BASE_OBJECTS];

	///	local finite element id
		std::vector<LFEID> m_vLFEID;

	///	local dof sets
		std::vector<const ILocalDoFSet*> m_vLocalDoFSet[NUM_REFERENCE_OBJECTS];

	///	maximum dimensions where dofs must be ordered
		std::vector<int> m_vMaxDimToOrderDoFs;

	///	number Dofs for local DoF set and subelement of element
		std::vector<MathMatrix<NUM_REFERENCE_OBJECTS,NUM_REFERENCE_OBJECTS, int> > m_vNumDoFOnSubelem;

	///	definition of function on subset
		std::vector<std::vector<bool> > m_vvFctDefInSubset;

	///	Multi Grid
		SmartPtr<MultiGrid> m_spMG;

	///	Subset Handler
		SmartPtr<MGSubsetHandler> m_spMGSH;

	///	Function Pattern
		FunctionPattern& m_rFctPatt;

	///	cached multigrid
		MultiGrid& m_rMultiGrid;

	///	Attachment type
		typedef ug::Attachment<size_t> ADoF;
		ADoF m_aIndex;

	///	Attachment Accessors
	///	\{
		typedef Grid::AttachmentAccessor<VertexBase, ADoF> vertex_attachment_accessor_type;
		typedef Grid::AttachmentAccessor<EdgeBase, ADoF> edge_attachment_accessor_type;
		typedef Grid::AttachmentAccessor<Face, ADoF> face_attachment_accessor_type;
		typedef Grid::AttachmentAccessor<Volume, ADoF> volume_attachment_accessor_type;
	/// \}

	///	Attachments
	///	\{
		vertex_attachment_accessor_type m_aaIndexVRT;
		edge_attachment_accessor_type m_aaIndexEDGE;
		face_attachment_accessor_type m_aaIndexFACE;
		volume_attachment_accessor_type m_aaIndexVOL;
	///	\}
};


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Lev Info
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

template <>
struct LevInfo<std::vector<size_t> > : public LevInfoBase
{
	typedef std::vector<size_t>::iterator iterator;
	typedef std::vector<size_t>::const_iterator const_iterator;

///	returns if free index avaiable
	inline bool free_index_available() const
	{
		return !vFreeIndex.empty();
	}

///	returns number of free index avaiable
	inline size_t num_free_index() const
	{
		return vFreeIndex.size();
	}

///	returns a free index
	inline size_t pop_free_index()
	{
		const size_t index = vFreeIndex.back();
		vFreeIndex.pop_back();
		return index;
	}

///	adds a free index, returns if index has not been contained before
	inline bool push_free_index(size_t index)
	{
		vFreeIndex.push_back(index); return true;
	}

///	returns iterators
///	\{
	inline iterator begin() {return vFreeIndex.begin();}
	inline iterator end() {return vFreeIndex.end();}
	inline const_iterator begin() const {return vFreeIndex.begin();}
	inline const_iterator end() const {return vFreeIndex.end();}
///	\}

///	clear container
	void clear() {vFreeIndex.clear();}

	protected:
	std::vector<size_t> vFreeIndex;
};

template <>
struct LevInfo<std::set<size_t> > : public LevInfoBase
{
	typedef std::set<size_t>::iterator iterator;
	typedef std::set<size_t>::const_iterator const_iterator;

///	returns if free index avaiable
	inline bool free_index_available() const
	{
		return !vFreeIndex.empty();
	}

///	returns number of free index avaiable
	inline size_t num_free_index() const
	{
		return vFreeIndex.size();
	}

///	returns a free index
	inline size_t pop_free_index()
	{
		const size_t index = *vFreeIndex.begin();
		vFreeIndex.erase(vFreeIndex.begin());
		return index;
	}

///	adds a free index, returns if index has not been contained before
	inline bool push_free_index(size_t index)
	{
	//	already contained, just return flag
		if(vFreeIndex.find(index) != vFreeIndex.end()) return false;

		vFreeIndex.insert(index); return true;
	}

///	returns iterators
///	\{
	inline iterator begin() {return vFreeIndex.begin();}
	inline iterator end() {return vFreeIndex.end();}
	inline const_iterator begin() const {return vFreeIndex.begin();}
	inline const_iterator end() const {return vFreeIndex.end();}
///	\}

///	clear container
	void clear() {vFreeIndex.clear();}

	protected:
	std::set<size_t> vFreeIndex;
};

} // end namespace ug

//#include "mg_dof_distribution_impl.h"
#include "level_dof_distribution.h"
#include "surface_dof_distribution.h"

#endif /* __H__UG__LIB_DISC__DOF_MANAGER__MG_DOF_DISTRIBUTION__ */
