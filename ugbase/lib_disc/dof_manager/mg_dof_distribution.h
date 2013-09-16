/*
 * mg_dof_distribution.h
 *
 *  Created on: 29.11.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__DOF_MANAGER__MG_DOF_DISTRIBUTION__
#define __H__UG__LIB_DISC__DOF_MANAGER__MG_DOF_DISTRIBUTION__

#include "function_pattern.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/local_finite_element/local_finite_element_id.h"
#include "lib_disc/common/local_algebra.h"
#include "dof_distribution_info.h"
#include "lib_algebra/parallelization/algebra_layouts.h"

namespace ug{

///////////////////////////////////////////////////////////////////////////////
// Lev Info
///////////////////////////////////////////////////////////////////////////////

struct LevInfo
{
///	constructor
	LevInfo() : numIndex(0)
#ifdef UG_PARALLEL
	,spAlgebraLayouts(new AlgebraLayouts)
#endif
	{}

/// number of distributed indices on whole domain
	size_t numIndex;

/// number of distributed indices on each subset
	std::vector<size_t> vNumIndexOnSubset;

#ifdef UG_PARALLEL
	protected:
///	algebra layouts
	SmartPtr<AlgebraLayouts> spAlgebraLayouts;

	public:
///	returns algebra layouts
///	\{
	ConstSmartPtr<AlgebraLayouts> layouts() const 	{return spAlgebraLayouts;}
	SmartPtr<AlgebraLayouts> layouts() 				{return spAlgebraLayouts;}
///	\}
#endif

///	clears the struct
	void clear()
	{
#ifdef UG_PARALLEL
		spAlgebraLayouts->clear();
#endif
		numIndex = 0;
		vNumIndexOnSubset.clear();
	}
};

///////////////////////////////////////////////////////////////////////////////
// MGDoFDistribution
///////////////////////////////////////////////////////////////////////////////

class MGDoFDistribution : virtual public DoFDistributionInfoProvider, public GridObserver
{
	public:
	//	type of multi index
		typedef MultiIndex<2> multi_index_type;

	public:
		MGDoFDistribution(SmartPtr<MultiGrid> spMG,
						  SmartPtr<MGSubsetHandler> spMGSH,
						  ConstSmartPtr<DoFDistributionInfo> spDDInfo,
		                  bool bGrouped);

		~MGDoFDistribution();

		///	returns the multigrid
		SmartPtr<MultiGrid> multi_grid() {return m_spMG;}
		ConstSmartPtr<MultiGrid> multi_grid() const {return m_spMG;}


		///	returns if dofs are grouped
		bool grouped() const {return m_bGrouped;}

		///	returns blocksize
		std::string blocksize() const {return "?";}


		/// returns number of levels
		inline int num_levels() const {return m_spMGSH->num_levels();}


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
		///	extracts the indices of the vertices
		template<typename TBaseElem>
		void indices_on_vertex(TBaseElem* elem, const ReferenceObjectID roid,
		                       LocalIndices& ind,
		                       const Grid::SecureVertexContainer& vElem) const;

		///	extract dofs on constrained objects
		template <typename TConstraining, typename TConstrained, typename TBaseElem>
		void constrained_vertex_indices(LocalIndices& ind,
		                         const typename Grid::traits<TBaseElem>::secure_container& vSubElem) const;

		template <typename TBaseElem,typename TConstraining, typename TConstrained, typename TSubElem>
		void constrained_edge_indices(TBaseElem* elem,LocalIndices& ind,
		                         const typename Grid::traits<TSubElem>::secure_container& vSubElem) const;

		template <typename TBaseElem,typename TConstraining, typename TConstrained, typename TSubElem>
		void constrained_face_indices(TBaseElem* elem,LocalIndices& ind,
		                         const typename Grid::traits<TSubElem>::secure_container& vSubElem) const;

		// sorts indices on constrained edges
		template <typename TBaseElem,typename TConstraining, typename TConstrained>
		void sort_constrained_edges(std::vector<size_t>& sortedInd,TBaseElem* elem,TConstraining* constrainingObj,size_t objIndex) const;

		// sorts indices on constrained faces
		template <typename TBaseElem,typename TConstraining, typename TConstrained>
		void sort_constrained_faces(std::vector<size_t>& sortedInd,TBaseElem* elem,TConstraining* constrainingObj,size_t objIndex) const;

		/// extracts the indices of the subelement of an element
		template<typename TBaseElem, typename TSubBaseElem>
		void indices(TBaseElem* elem, const ReferenceObjectID roid,
		             LocalIndices& ind,
		             const typename Grid::traits<TSubBaseElem>::secure_container& vElem) const;

		/// extracts the indices of a subelement of an element
		template<typename TBaseElem, typename TSubBaseElem>
		void multi_indices(TBaseElem* elem, const ReferenceObjectID roid,
		                   size_t fct, std::vector<multi_index_type>& ind,
		                   const typename Grid::traits<TSubBaseElem>::secure_container& vElem) const;


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
		template <typename TBaseObject>
		bool add(TBaseObject* obj, const ReferenceObjectID roid,
		         const int si, LevInfo& li);
		bool add(GeometricObject* obj, const ReferenceObjectID roid,
		         const int si, LevInfo& li);

		///	checks that subset assigment is ok
		void check_subsets();

		/// initializes the attachments
		void init_attachments();

		/// removes the attachments
		void clear_attachments();

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

		MessageHub::SPCallbackId	m_callbackId_GridCreation;

	///	Multi Grid
		SmartPtr<MultiGrid> m_spMG;
		MultiGrid* m_pMG;
		SmartPtr<MGSubsetHandler> m_spMGSH;

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

} // end namespace ug


#endif /* __H__UG__LIB_DISC__DOF_MANAGER__MG_DOF_DISTRIBUTION__ */
