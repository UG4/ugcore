/*
 * local_algebra.h
 *
 *  Created on: 19.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__COMMON__LOCAL_ALGEBRA__
#define __H__LIB_DISCRETIZATION__COMMON__LOCAL_ALGEBRA__

#include <vector>

#include "./multi_index.h"
#include "./function_group.h"
#include "lib_algebra/small_algebra/small_algebra.h"

namespace ug{

/** Local index set
 *
 * The class LocalIndices is used to represent a local index set. It is intended to suit
 * for local index sets arising in Finite Element/Finite Volume assemblings. Thus, one
 * can load the LocalIndices for a FunctionGroup, representing the set of functions in
 * the discretization. It is also possible to select only a SubGroup of the Function Group.
 * This is done by adding a SubFunctionMap. By default, the SubFunctionMap is identity.
 * The access form discretization side is thus by a pair (fct, dof) representing the selected
 * function (0 <= fct <= # Selected Functions - 1) and local dof indices (0 <= dof <= #DoFs -1).
 * The access from algebra side is only an access to the algebra i-index. Using
 *  - local_index(size_t fct, size_t dof)
 *  - comp(size_t fct, size_t dof)
 *  will return the correspinding (i,alpha) entry for a (fct, dof) pair.
 */
class LocalIndices
{
	public:
		// index type used by algebra
		typedef size_t index_type;

		// component type used by algebra
		typedef size_t comp_type;

		// multi indices used by discretization
		typedef MultiIndex<2> multi_index_type;

	public:
		LocalIndices() : m_pFunctionGroup(NULL), m_bAccessAll(true), m_pSubFunctionMap(NULL)
		{
			m_vAlgIndices.clear();
			m_vvDofIndices.clear();
		};

		///////////////////////////
		// setup
		///////////////////////////

		/// adjust for the function group
		void set_function_group(const FunctionGroup& funcGroup)
		{
			m_pFunctionGroup = &funcGroup;

			// resize vectors
			resize_dof_indices();

			// create default SubFunctionMap containing all functions of function Group
			m_identitySubFunctionMap.clear();
			for(size_t fct = 0; fct < funcGroup.num_fct(); ++fct)
			{
				m_identitySubFunctionMap.select(fct);
			}

			// default function group
			m_pSubFunctionMap = &m_identitySubFunctionMap;
		}

		/// access only part of the function group
		bool access_sub_function_group(const SubFunctionMap& subFuncMap)
		{
			m_pSubFunctionMap = &subFuncMap;
			m_bAccessAll = false;
			return true;
		}

		/// access all functions in FunctionGroup
		bool access_all()
		{
			m_pSubFunctionMap = NULL;
			m_bAccessAll = true;
			return true;
//			return access_sub_function_group(m_identitySubFunctionMap);
		}

		/// get current sub function map
		const SubFunctionMap& get_sub_function_map() const
		{
			return *m_pSubFunctionMap;
		}

		/// clear all dofs and indices
		void clear()
		{
			m_vAlgIndices.clear();
			for(size_t i = 0; i < m_vvDofIndices.size(); ++i)
				m_vvDofIndices[i].clear();
		}

		/// add local dof index for local function
		void add_dof(size_t fct, const multi_index_type& dofIndex)
		{
			UG_ASSERT(fct < num_fct(), "Local function index not valid.");
			m_vvDofIndices[fct].push_back(dofIndex);
		}

		/// set number of algebra indices
		void set_num_indices(size_t numInd) {m_vAlgIndices.resize(numInd);}

		/// set algebra index
		void set_index(size_t i, index_type ind)
		{
			UG_ASSERT(i < num_indices(), "Algebra index not valid. (i="<<i<<", num_indices="<<num_indices()<<")");
			m_vAlgIndices[i] = ind;
		}

		///////////////////////////
		// discretization access
		///////////////////////////

		/// number of selected functions
		size_t num_fct() const {return m_pSubFunctionMap->num_fct();}

		/// return global function id for selected functions
		size_t fct_id(size_t fct) const
		{
			UG_ASSERT(fct < num_fct(), "local fct index '"<< fct <<"' not valid. (num_fct = "<<num_fct()<<")");
			const size_t accFct = (*m_pSubFunctionMap)[fct];
			return m_pFunctionGroup->fct_id(accFct);
		}

		/// number of dofs of selected functions (sum)
		size_t num_dofs() const
		{
			// TODO: This could be cached. Would it make sense ???
			size_t num = 0;
			for(size_t fct = 0; fct < num_fct(); ++fct)
			{
				num += num_dofs(fct);
			}
			return num;
		}

		/// number of dofs per selected function
		size_t num_dofs(size_t fct) const
		{
			UG_ASSERT(fct < num_fct(), "local fct index '"<< fct <<"' not valid. (num_fct = "<<num_fct()<<")");
			const size_t accFct = (*m_pSubFunctionMap)[fct];
			UG_ASSERT(accFct < m_vvDofIndices.size(), "accFct index " << accFct << " not valid. (Size = " << m_vvDofIndices.size() << ")\n");
			return m_vvDofIndices[accFct].size();
		}

		/// local index for (selected fct, dof)
		index_type local_index(size_t fct, size_t dof) const
		{
			UG_ASSERT(fct < num_fct(), "local fct index '"<< fct <<"' not valid. (num_fct = "<<num_fct()<<")");
			if(m_bAccessAll)
			{
				UG_ASSERT(fct < m_vvDofIndices.size(), "fct index " << fct << " not valid. (Size = " << m_vvDofIndices.size() << ")\n");
				UG_ASSERT(dof < m_vvDofIndices[fct].size(), "local dof index not valid");
				return  m_vvDofIndices[fct][dof][0];
			}

			const size_t accFct = (*m_pSubFunctionMap)[fct];
			UG_ASSERT(accFct < m_vvDofIndices.size(), "accFct index " << accFct << " not valid. (Size = " << m_vvDofIndices.size() << ")\n");
			UG_ASSERT(dof < m_vvDofIndices[accFct].size(), "local dof index not valid");

			return m_vvDofIndices[accFct][dof][0];
		}

		/// global algebra index for (selected fct, dof)
		index_type global_index(size_t fct, size_t dof) const
		{
			const size_t loc_i = local_index(fct, dof);
			return index(loc_i);
		}

		/// algebra comp for (selected fct, dof)
		comp_type comp(size_t fct, size_t dof) const
		{
			UG_ASSERT(fct < num_fct(), "local fct index '"<< fct <<"' not valid. (num_fct = "<<num_fct()<<")");
			if(m_bAccessAll)
			{
				UG_ASSERT(fct < m_vvDofIndices.size(), "fct index " << fct << " not valid. (Size = " << m_vvDofIndices.size() << ")\n");
				UG_ASSERT(dof < m_vvDofIndices[fct].size(), "local dof index not valid");
				return  m_vvDofIndices[fct][dof][1];
			}

			const size_t accFct = (*m_pSubFunctionMap)[fct];
			UG_ASSERT(accFct < m_vvDofIndices.size(), "accFct index " << accFct << " not valid. (Size = " << m_vvDofIndices.size() << ")\n");
			UG_ASSERT(dof < m_vvDofIndices[accFct].size(), "local dof index not valid");
			return m_vvDofIndices[accFct][dof][1];
		}

		multi_index_type local_multi_index(size_t fct, size_t dof) const
		{
			UG_ASSERT(fct < num_fct(), "local fct index '"<< fct <<"' not valid. (num_fct = "<<num_fct()<<")");
			if(m_bAccessAll)
			{
				UG_ASSERT(fct < m_vvDofIndices.size(), "fct index " << fct << " not valid. (Size = " << m_vvDofIndices.size() << ")\n");
				UG_ASSERT(dof < m_vvDofIndices[fct].size(), "local dof index not valid");
				return  m_vvDofIndices[fct][dof];
			}

			const size_t accFct = (*m_pSubFunctionMap)[fct];
			UG_ASSERT(accFct < m_vvDofIndices.size(), "accFct index " << accFct << " not valid. (Size = " << m_vvDofIndices.size() << ")\n");
			UG_ASSERT(dof < m_vvDofIndices[accFct].size(), "local dof index not valid");
			return m_vvDofIndices[accFct][dof];
		}

		///////////////////////////
		// algebra access
		///////////////////////////

		/// number of all dof indices used by matrix
		size_t num_indices() const {return m_vAlgIndices.size();}

		/// algebra index
		index_type index(size_t i) const
		{
			UG_ASSERT(i < num_indices(), "Algebra index not valid");
			return m_vAlgIndices[i];
		}

	protected:
		// resize the dof_index arrays for the current function group
		void resize_dof_indices()
		{
			m_vvDofIndices.resize(m_pFunctionGroup->num_fct());
			for(size_t fct = 0; fct < m_vvDofIndices.size(); ++fct)
			{
				m_vvDofIndices[fct].clear();
			}
		}

	protected:
		// current function group
		const FunctionGroup* m_pFunctionGroup;

		// access all function (flag)
		bool m_bAccessAll;

		// Mapping (fct, dof) -> local index
		std::vector<std::vector<multi_index_type> > m_vvDofIndices;

		// algebra indices
		std::vector<index_type> m_vAlgIndices;

		// sub function group, that is accessed
		const SubFunctionMap* m_pSubFunctionMap;

		// default sub function group, i.e. all functions
		SubFunctionMap m_identitySubFunctionMap;


};

template <typename TEntry>
class LocalVector
{
	public:
		typedef LocalVector<TEntry> this_type;

		// index type used by algebra
		typedef size_t index_type;

		// component type used by algebra
		typedef size_t comp_type;

		// entry type used by algebra
		typedef TEntry value_type;

	protected:
		typedef typename std::vector<value_type>::iterator iterator;

	public:
		///////////////////////////
		// setup
		///////////////////////////

		LocalVector() : m_pIndices(NULL) {m_entries.clear();}

		LocalVector(LocalIndices& ind)
			: m_pIndices(NULL)
		{m_entries.clear(); set_indices(ind);}

		/// set new local indices
		void set_indices(LocalIndices& ind)
		{
			if(ind.num_indices() != size()) m_entries.resize(ind.num_indices());
			m_pIndices = &ind;
		}

		/// get current local indices
		LocalIndices& get_indices() {return *m_pIndices;}

		template <typename TDiscreteFunction>
		void read_values(const TDiscreteFunction& v)
		{
			for(size_t i = 0; i < size(); ++i)
				m_entries[i] = v[index(i)];
		}

		///////////////////////////
		// vector functions
		///////////////////////////

		/// size
		size_t size() const {return m_entries.size();}

		/// access to entry
		value_type& operator[](size_t i)
		{
			UG_ASSERT(i < size(), "Index not valid");
			return m_entries[i];
		}

		/// access to entry
		const value_type& operator[](size_t i) const
		{
			UG_ASSERT(i < size(), "Index not valid");
			return m_entries[i];
		}

		/// set all components of the vector
		void set(number val)
		{
			iterator iterEnd = m_entries.end();
			for(iterator iter = m_entries.begin(); iter != iterEnd; ++iter)
				(*iter) = val;
		}

		/// multiply all components of the vector
		this_type& operator*(number val)
		{
			iterator iterEnd = m_entries.end();
			for(iterator iter = m_entries.begin(); iter != iterEnd; ++iter)
				(*iter) *= val;
			return *this;
		}

		/// add a vector
		this_type& operator+=(const this_type& rhs)
		{
			UG_ASSERT(size() == rhs.size(), "Size does not match.");
			for(size_t i = 0; i < size(); ++i)
				m_entries[i] += rhs[i];
			return *this;
		}

		/// subtract a vector
		this_type& operator-=(const this_type& rhs)
		{
			UG_ASSERT(size() == rhs.size(), "Size does not match.");
			for(size_t i = 0; i < size(); ++i)
				m_entries[i] -= rhs[i];
			return *this;
		}

		///////////////////////////
		// algebra access
		///////////////////////////

		/// index of component
		index_type index(size_t i) const {return m_pIndices->index(i);}

		///////////////////////////
		// discretization access
		///////////////////////////

		/// number of selected functions
		size_t num_fct() const {return m_pIndices->num_fct();}

		/// return global function id for selected functions
		size_t fct_id(size_t fct) const	{return m_pIndices->fct_id(fct);}

		/// number of dofs of selected functions (sum)
		size_t num_dofs() const {return m_pIndices->num_dofs();}

		/// number of dofs per selected function
		size_t num_dofs(size_t fct) const {return m_pIndices->num_dofs(fct);}

		/// access to dof of function fct
		number& operator()(size_t fct, size_t dof)
		{
			const index_type index = m_pIndices->local_index(fct, dof);
			const comp_type comp = m_pIndices->comp(fct, dof);
			UG_ASSERT(index < size(), "Index not valid");
			return BlockRef(m_entries[index], comp);
		}

		/// const access to dof of function fct
		const number& operator()(size_t fct, size_t dof) const
		{
			const index_type index = m_pIndices->local_index(fct, dof);
			const comp_type comp = m_pIndices->comp(fct, dof);
			UG_ASSERT(index < size(), "Index not valid");
			return BlockRef(m_entries[index], comp);
		}

	protected:
		// indices
		LocalIndices* m_pIndices;

		// entries (size = m_indices.num_indices())
		std::vector<value_type> m_entries;
};

class LocalMatrixBase{
public:
	virtual ~LocalMatrixBase(){};
};

template <typename TEntry>
class LocalMatrix : public LocalMatrixBase
{
	public:
		typedef LocalMatrix<TEntry> this_type;

		// index type used by algebra
		typedef size_t index_type;

		// component type used by algebra
		typedef size_t comp_type;

		// entry type used by algebra
		typedef TEntry value_type;

	protected:
		typedef typename std::vector<std::vector<value_type> >::iterator   row_iterator;
		typedef typename std::vector<value_type>::iterator col_iterator;

	public:
		///////////////////////////
		// setup
		///////////////////////////

		LocalMatrix() : m_pRowIndices(NULL), m_pColIndices(NULL) {}

		LocalMatrix(const LocalIndices& rowInd, const LocalIndices& colInd)
			: m_pRowIndices(NULL), m_pColIndices(NULL)
		{
			set_row_indices(rowInd);
			set_col_indices(colInd);
		}

		/// set new local row indices
		void set_row_indices(const LocalIndices& ind)
		{
			m_pRowIndices = &ind;
			resize();
		}

		/// set new local column indices
		void set_col_indices(const LocalIndices& ind)
		{
			m_pColIndices = &ind;
			resize();
		}

		/// set new local indices
		void set_indices(const LocalIndices& rowInd, const LocalIndices& colInd)
		{
			set_row_indices(rowInd);
			set_col_indices(colInd);
		}

		template <typename TMatrix>
		void read_values(const TMatrix& mat)
		{
			for(size_t i = 0; i < num_rows(); ++i)
				for(size_t j = 0; j < num_cols(); ++j)
					m_entries(i,j) = mat(row_index(i), col_index(j));
		}

		///////////////////////////
		// matrix functions
		///////////////////////////

		/// number of rows
		size_t num_rows() const
		{
			//return m_entries.size();
			return m_entries.num_rows();
		}

		/// number of columns
		size_t num_cols() const
		{
			return m_entries.num_cols();
		}

		/// access to entry
		value_type& operator() (size_t i, size_t j)
		{
			UG_ASSERT(i < num_rows(), "Row does not exist.");
			UG_ASSERT(j < num_cols(), "Column does not exist.");
			return m_entries(i,j);
		}

		/// const access to entry
		const value_type& operator() (size_t i, size_t j) const
		{
			UG_ASSERT(i < num_rows(), "Row does not exist.");
			UG_ASSERT(j < num_cols(), "Column does not exist.");
			return m_entries(i,j);
		}

		/// set all entries
		void set(number val)
		{
			m_entries = val;
		}

		/// multiply all entries
		this_type& operator*(number val)
		{
			m_entries *= val;
			return *this;
		}

		/// add matrix
		this_type& operator+=(const this_type& rhs)
		{
			UG_ASSERT(num_rows() == rhs.num_rows(), "Row size does not match");
			UG_ASSERT(num_cols() == rhs.num_cols(), "Column size does not match");

			m_entries += rhs.m_entries;
			return *this;
		}

		/// subtract matrix
		this_type& operator-=(const this_type& rhs)
		{
			UG_ASSERT(num_rows() == rhs.num_rows(), "Row size does not match");
			UG_ASSERT(num_cols() == rhs.num_cols(), "Column size does not match");

			m_entries -= rhs.m_entries;
			return *this;
		}

		///////////////////////////
		// algebra access
		///////////////////////////

		/// index of component
		index_type row_index(size_t i) const {return m_pRowIndices->index(i);}
		index_type col_index(size_t j) const {return m_pColIndices->index(j);}

		comp_type row_sub_index(size_t rowFct, size_t rowDoF) const {return m_pRowIndices->comp(rowFct, rowDoF);}
		comp_type col_sub_index(size_t colFct, size_t colDoF) const {return m_pColIndices->comp(colFct, colDoF);}

		///////////////////////////
		// discretization access
		///////////////////////////

		number& operator()(size_t rowFct, size_t rowDof, size_t colFct, size_t colDof)
		{
			const typename LocalIndices::multi_index_type rowInd = m_pRowIndices->local_multi_index(rowFct, rowDof);
			const typename LocalIndices::multi_index_type colInd = m_pColIndices->local_multi_index(colFct, colDof);
			UG_ASSERT(rowInd[0] < num_rows(), "Row does not exist.");
			UG_ASSERT(colInd[0] < num_cols(), "Column does not exist.");

			return BlockRef(m_entries(rowInd[0],colInd[0]), rowInd[1], colInd[1]);
		}

		/// const access to dof of function fct
		const number& operator()(size_t rowFct, size_t rowDof, size_t colFct, size_t colDof) const
		{
			const typename LocalIndices::multi_index_type rowInd = m_pRowIndices->local_multi_index(rowFct, rowDof);
			const typename LocalIndices::multi_index_type colInd = m_pColIndices->local_multi_index(colFct, colDof);
			UG_ASSERT(rowInd[0] < num_rows(), "Row does not exist.");
			UG_ASSERT(colInd[0] < num_cols(), "Column does not exist.");

			return BlockRef(m_entries(rowInd[0],colInd[0]), rowInd[1], colInd[1]);
		}

	protected:
		void resize()
		{
			m_entries.resize(0,0);

			if(m_pRowIndices == NULL || m_pColIndices == NULL) return;

			const size_t num_rows = m_pRowIndices->num_indices();
			const size_t num_cols = m_pColIndices->num_indices();

			m_entries.resize(num_rows, num_cols);
		}

	protected:
		// row indices
		const LocalIndices* m_pRowIndices;

		// column indices
		const LocalIndices* m_pColIndices;

		// entries
		DenseMatrix<VariableArray2<value_type> > m_entries;
};

template <typename Entry>
std::ostream& operator<< (std::ostream& outStream, const ug::LocalMatrix<Entry>& mat)
{
	for(size_t i = 0; i < mat.num_rows(); ++i)
		for(size_t j = 0; j < mat.num_cols(); ++j)
		{
			outStream << "["<<i<<","<<j<<"]: (" << mat.row_index(i) << " -> " << mat.col_index(j) << ") : " << mat(i,j) << "\n";
		}
	return outStream;
}

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__COMMON__LOCAL_ALGEBRA__ */
