/*
 * ass_adapter.h
 *
 *  Created on: 04.02.2013
 *      Author: raphaelprohl, Andreas Vogel
 */

#ifndef ASS_ADAPTER_H_
#define ASS_ADAPTER_H_

#include "lib_grid/tools/bool_marker.h"
#include "lib_grid/tools/selector_grid.h"
#include "lib_disc/spatial_disc/local_to_global/local_to_global_mapper.h"

namespace ug{

template <typename TDomain, typename TAlgebra>
class IDomainConstraint;


template <typename TAlgebra>
class LocalToGlobalMapper : public ILocalToGlobalMapper<TAlgebra>
{
	public:
	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	///	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	public:
	///	default constructor
		LocalToGlobalMapper() {}

	///	adds a local vector to the global one
		void AddLocalVec(vector_type& vec, const LocalVector& lvec, ConstSmartPtr<DoFDistribution> dd)
			{ AddLocalVector(vec,lvec);}

	///	adds a local matrix to the global one
		void AddLocalMatToGlobal(matrix_type& mat, const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd)
			{ AddLocalMatrixToGlobal(mat,lmat);}

	///	destructor
		~LocalToGlobalMapper() {};
};

// AssAdapter combines tools to adapt the assemble routine
template <typename TAlgebra>
class AssAdapter
{
	public:
	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	///	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	///	Type of algebra value
		typedef typename vector_type::value_type value_type;

	public:
	/// constructor
		AssAdapter(): m_pBoolMarker(NULL), m_pSelector(NULL)
		{
			m_assIndex.index_set = false;
			m_pMapper = &m_pMapperCommon;
		}

	/// set local to global mapping
		void set_mapping(ILocalToGlobalMapper<TAlgebra>* pMapper = NULL)
		{
			if(pMapper){ m_pMapper = pMapper;}
			else{ m_pMapper = &m_pMapperCommon;}
		}

		void AddLocalVec(vector_type& vec, const LocalVector& lvec, ConstSmartPtr<DoFDistribution> dd)
		{ m_pMapper->AddLocalVec(vec, lvec, dd);}
		void AddLocalMatToGlobal(matrix_type& mat, const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd)
		{ m_pMapper->AddLocalMatToGlobal(mat, lmat, dd);}

		void resize(ConstSmartPtr<DoFDistribution> dd, vector_type& vec);
		void resize(ConstSmartPtr<DoFDistribution> dd, matrix_type& mat);

	///	sets a marker to exclude elements from assembling
	/**
	 * This methods sets a marker. Only elements that are marked will be
	 * assembled during assembling process. If no marker is set, this
	 * corresponds to a marker where all elements have been marked.
	 *
	 * \param[in]	mark	BoolMarker
	 */
		void set_marker(BoolMarker* mark = NULL){ m_pBoolMarker = mark; }
		BoolMarker* marker(){ return m_pBoolMarker;}

	///	sets a selector of elements for assembling
	/**
	 * This methods sets an element list. Only elements of this list will be
	 * assembled during assembling process. Especially the list defines the begin
	 * and end of the element-iterator in the element assembling-loop.
	 * If no element list is set, this corresponds to a assembling where the loop is
	 * carried out over all elements of a subset.
	 *
	 * \param[in]	sel		Selector
	 */
		void set_selector(Selector* sel = NULL){ m_pSelector = sel; }
		Selector* selector(){ return m_pSelector;}

	///	sets an index for which the assembling should be carried out
	/**
	 * This methods sets a boolean if an index-wise assemble routine should be used.
	 * This proceeding is e.g. useful for a nonlinear Gauss-Seidel or nonlinear
	 * Jacobi solver. The specific index is passed.
	 *
	 * \param[in]	ind			size_t
	 * \param[in]	index_set	bool
	 */
		void set_ass_index(){ set_ass_index(0, false);}
		void set_ass_index(size_t ind, bool index_set = true)
		{
			m_assIndex.index = ind; m_assIndex.index_set = index_set;
		}
		///	returns whether the assemble Index is set or not
		size_t is_ass_index_set(){ return m_assIndex.index_set;}
		///	gets assemble Index
		size_t ass_index(){ return m_assIndex.index;}

		template <typename TElem>
		void elemIter_fromSel(ConstSmartPtr<DoFDistribution> dd, int si,
				std::vector<TElem*>& elems);

		template <typename TDomain>
		void adaptConstraint(SmartPtr<IDomainConstraint<TDomain, TAlgebra> >& constraint);

		void adjust_matrix(matrix_type& mat, const size_t index);
		void adjust_vector(vector_type& vec, const size_t index, const value_type& val);
	//private:

	///	marker used to skip elements
		BoolMarker* m_pBoolMarker;

	///	selector used to set a list of elements for the assembling
		Selector* 	m_pSelector;

		struct AssIndex{
			///	should assemble be index-wise
			bool index_set;

			///	current index
			size_t index;
		};

	///	object for index-wise assemble routine
		AssIndex 	m_assIndex;

	///	LocalToGlobalMapper-Func
		ILocalToGlobalMapper<TAlgebra>* m_pMapper;
		LocalToGlobalMapper<TAlgebra> 	m_pMapperCommon;
};

} // end namespace ug

#include "ass_adapter_impl.h"

#endif /* ASS_ADAPTER_H_ */
