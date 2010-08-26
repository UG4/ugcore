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

#include "../dof_distribution.h"
#include "../function_pattern.h"
#include "../../common/common.h"
#include "../../reference_element/reference_element.h"

namespace ug{

class P1StorageManager
{
//	for DofManager
	public:
		P1StorageManager() : m_pSH(NULL) {}

		/// set subset handler
		void set_subset_handler(ISubsetHandler& sh);

		/// clear all dofs
		void clear();

	// TODO: Why all public ?
	public:
		/// attach dofs
		void update_attachments();

		/// subset handler
		ISubsetHandler* m_pSH;

		/// Subset Infos
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
		bool add_discrete_function(std::string name, LocalShapeFunctionSetID id, int dim);

		bool add_discrete_function(std::string name, LocalShapeFunctionSetID id, const SubsetGroup& SubsetIndices, int dim);
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
		inline size_t num() const {return m_goc.num<TElem>();}

		// iterator for elements where this grid function is defined
		template <typename TElem>
		inline typename geometry_traits<TElem>::iterator begin()
			{return m_goc.begin<TElem>();}

		// iterator for elements where this grid function is defined
		template <typename TElem>
		inline typename geometry_traits<TElem>::iterator end()
			{return m_goc.end<TElem>();}

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
		// LocalIndex update
		///////////////////////////////////////

		/// number of algebra indices for refID, subset and function group (Element + Closure of Element)
		size_t num_indices(ReferenceObjectID refID, int si, const FunctionGroup& funcGroup) const;

		/// number of algebra indices for refType, subset and function group (only inner part of Element)
		size_t num_inner_indices(ReferenceObjectID refID, int si, const FunctionGroup& funcGroup) const;

		/// fill local informations in LocalIndex (Element + Closure of Element)
		bool prepare_indices(ReferenceObjectID refID, int si, LocalIndices& ind, bool withHanging = false) const;

		/// fill local informations in LocalIndex (only inner part of Element)
		bool prepare_inner_indices(ReferenceObjectID refID, int si, LocalIndices& ind) const;

		/// fill the global algebra indices in LocalIndex (Element + Closure of Element)
		template<typename TElem>
		void update_indices(TElem* elem, LocalIndices& ind, bool withHanging = false) const;

		/// fill the global algebra indices in LocalIndex (only inner part of Element)
		template<typename TElem>
		void update_inner_indices(TElem* elem, LocalIndices& ind) const;


		///////////////////////////////////////
		// Multi index access
		///////////////////////////////////////

		/// number of multi indices on element for a function (Element + Closure of Element)
		template<typename TElem>
		size_t num_multi_indices(TElem* elem, size_t fct) const;

		/// number of multi indices on element for a function (only inner part of Element)
		template<typename TElem>
		size_t num_inner_multi_indices(TElem* elem, size_t fct) const;

		/// get multi indices of element for a function fct (Element + Closure of Element)
		template<typename TElem>
		size_t get_multi_indices(TElem* elem, size_t fct, multi_index_vector_type& ind) const;

		/// get multi indices of element for a function fct (only inner part of Element)
		template<typename TElem>
		size_t get_inner_multi_indices(TElem* elem, size_t fct, multi_index_vector_type& ind) const;

		///////////////////////////////////////
		// Algebra index access
		///////////////////////////////////////

		/// number of algebra indices on element for a function (Element + Closure of Element)
		template<typename TElem>
		size_t num_algebra_indices(TElem* elem, size_t fct) const;

		/// number of algebras indices on element for a function (only inner part of Element)
		template<typename TElem>
		size_t num_inner_algebra_indices(TElem* elem, size_t fct) const;

		/// get algebra indices of element for a function fct (Element + Closure of Element)
		template<typename TElem>
		void get_algebra_indices(TElem* elem, algebra_index_vector_type& ind) const;

		/// get algebra indices of element for a function fct (only inner part of Element)
		template<typename TElem>
		void get_inner_algebra_indices(TElem* elem, algebra_index_vector_type& ind) const;

		///////////////////////////
		// Creation
		///////////////////////////

		/// distribute dofs for given goc
		bool distribute_dofs();

	protected:
		VertexBase* get_vertex(VertexBase* vrt, size_t i) const;
		VertexBase* get_vertex(EdgeBase* edge, size_t i) const;
		VertexBase* get_vertex(Face* face, size_t i) const;
		VertexBase* get_vertex(Volume* vol, size_t i) const;

	protected:
		// geometric object collection for this Distributor
		GeometricObjectCollection m_goc;

		// subset handler for this distributor
		ISubsetHandler* m_pISubsetHandler;

		// Storage Manager for dofs
		StorageManager* m_pStorageManager;

		// number of distributed dofs on whole domain
		size_t m_numDoFs;

		// number offsetmap
		std::vector<std::vector<size_t> > m_vvOffsets;

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
		inline size_t num() const {return m_goc.num<TElem>();}

		// iterator for elements where this grid function is defined
		template <typename TElem>
		inline typename geometry_traits<TElem>::iterator begin()
			{return m_goc.begin<TElem>();}

		// iterator for elements where this grid function is defined
		template <typename TElem>
		inline typename geometry_traits<TElem>::iterator end()
			{return m_goc.end<TElem>();}

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
		// LocalIndex update
		///////////////////////////////////////

		/// number of algebra indices for refID, subset and function group (Element + Closure of Element)
		size_t num_indices(ReferenceObjectID refID, int si, const FunctionGroup& funcGroup) const;

		/// number of algebra indices for refType, subset and function group (only inner part of Element)
		size_t num_inner_indices(ReferenceObjectID refID, int si, const FunctionGroup& funcGroup) const;

		/// fill local informations in LocalIndex (Element + Closure of Element)
		bool prepare_indices(ReferenceObjectID refID, int si, LocalIndices& ind) const;

		/// fill local informations in LocalIndex (only inner part of Element)
		bool prepare_inner_indices(ReferenceObjectID refID, int si, LocalIndices& ind) const;

		/// fill the global algebra indices in LocalIndex (Element + Closure of Element)
		template<typename TElem>
		void update_indices(TElem* elem, LocalIndices& ind) const;

		/// fill the global algebra indices in LocalIndex (only inner part of Element)
		template<typename TElem>
		void update_inner_indices(TElem* elem, LocalIndices& ind) const;


		///////////////////////////////////////
		// Multi index access
		///////////////////////////////////////

		/// number of multi indices on element for a function (Element + Closure of Element)
		template<typename TElem>
		size_t num_multi_indices(TElem* elem, size_t fct) const;

		/// number of multi indices on element for a function (only inner part of Element)
		template<typename TElem>
		size_t num_inner_multi_indices(TElem* elem, size_t fct) const;

		/// get multi indices of element for a function fct (Element + Closure of Element)
		template<typename TElem>
		size_t get_multi_indices(TElem* elem, size_t fct, multi_index_vector_type& ind) const;

		/// get multi indices of element for a function fct (only inner part of Element)
		template<typename TElem>
		size_t get_inner_multi_indices(TElem* elem, size_t fct, multi_index_vector_type& ind) const;

		///////////////////////////////////////
		// Algebra index access
		///////////////////////////////////////

		/// number of algebra indices on element for a function (Element + Closure of Element)
		template<typename TElem>
		size_t num_algebra_indices(TElem* elem, size_t fct) const;

		/// number of algebras indices on element for a function (only inner part of Element)
		template<typename TElem>
		size_t num_inner_algebra_indices(TElem* elem, size_t fct) const;

		/// get algebra indices of element for a function fct (Element + Closure of Element)
		template<typename TElem>
		void get_algebra_indices(TElem* elem, algebra_index_vector_type& ind) const;

		/// get algebra indices of element for a function fct (only inner part of Element)
		template<typename TElem>
		void get_inner_algebra_indices(TElem* elem, algebra_index_vector_type& ind) const;

		///////////////////////////
		// Creation
		///////////////////////////

		/// distribute dofs for given goc
		bool distribute_dofs();

	protected:
		VertexBase* get_vertex(VertexBase* vrt, size_t i) const;
		VertexBase* get_vertex(EdgeBase* edge, size_t i) const;
		VertexBase* get_vertex(Face* face, size_t i) const;
		VertexBase* get_vertex(Volume* vol, size_t i) const;

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


#include "./p1conform_impl.h"

#endif /* __H__LIB_DISCRETIZATION__DOF_MANAGER__P1CONFORM__ */
