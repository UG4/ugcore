/*
 * local_algebra.h
 *
 *  Created on: 19.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__LOCAL_ALGEBRA__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__LOCAL_ALGEBRA__

#include <vector>

namespace ug{

class LocalIndices
{
	public:
		// index type used by algebra
		typedef size_t index_type;

		// component type used by algebra
		typedef size_t comp_type;

		// multi indices used by discretization
		typedef MultiIndex<2> dof_index_type;

	public:
		LocalIndices() : m_pFunctionGroup(NULL)
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

			// create default SubFunctionGroup containing all functions of function Group
			m_defaultSubFunctionGroup.clear();
			for(size_t fct = 0; fct < funcGroup.num_fct(); ++fct)
			{
				m_defaultSubFunctionGroup.select(fct);
			}

			// access default function group
			access_sub_function_group(m_defaultSubFunctionGroup);
		}

		/// access only part of the function group
		bool access_sub_function_group(const SubFunctionGroup& subFuncGroup)
		{
			m_pSubFunctionGroup = &subFuncGroup;
			return true;
		}

		/// clear all dofs and indices
		void clear()
		{
			m_vAlgIndices.clear();
			for(size_t i = 0; i < m_vvDofIndices.size(); ++i)
				m_vvDofIndices[i].clear();
		}

		/// add local dof index for local function
		void add_dof(size_t fct, const dof_index_type& dofIndex)
		{
			UG_ASSERT(fct < num_fct(), "Local function index not valid.");
			m_vvDofIndices[fct].push_back(dofIndex);
		}

		/// set number of algebra indices
		void set_num_indices(size_t numInd) {m_vAlgIndices.resize(numInd);}

		/// set algebra index
		void set_index(size_t i, index_type ind)
		{
			UG_ASSERT(i < num_indices(), "Algebra index not valid.");
			m_vAlgIndices[i] = ind;
		}

		///////////////////////////
		// discretization access
		///////////////////////////

		/// number of functions handled
		size_t num_fct() const {return m_pSubFunctionGroup->num_fct();}

		/// return global function id
		size_t fct_id(size_t fct) const
		{
			UG_ASSERT(fct < num_fct(), "local fct index '"<< fct <<"' not valid. (num_fct = "<<num_fct()<<")");
			const size_t accFct = (*m_pSubFunctionGroup)[fct];
			return m_pFunctionGroup->fct_id(accFct);
		}

		/// number of dofs per local function
		size_t num_dofs(size_t fct) const
		{
			UG_ASSERT(fct < num_fct(), "local fct index '"<< fct <<"' not valid. (num_fct = "<<num_fct()<<")");
			const size_t accFct = (*m_pSubFunctionGroup)[fct];
			return m_vvDofIndices[accFct].size();
		}

		/// local index for (fct, dof)
		index_type local_index(size_t fct, size_t dof) const
		{
			UG_ASSERT(fct < num_fct(), "local fct index '"<< fct <<"' not valid. (num_fct = "<<num_fct()<<")");
			const size_t accFct = (*m_pSubFunctionGroup)[fct];
			UG_ASSERT(dof < m_vvDofIndices[accFct].size(), "local dof index not valid");
			return m_vvDofIndices[accFct][dof][0];
		}

		/// algebra comp for (fct, dof)
		comp_type comp(size_t fct, size_t dof) const
		{
			UG_ASSERT(fct < num_fct(), "local fct index '"<< fct <<"' not valid. (num_fct = "<<num_fct()<<")");
			const size_t accFct = (*m_pSubFunctionGroup)[fct];
			UG_ASSERT(dof < m_vvDofIndices[accFct].size(), "local dof index not valid");
			return m_vvDofIndices[accFct][dof][1];
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

		// sub function group, that is accessed
		const SubFunctionGroup* m_pSubFunctionGroup;

		// default sub function group, i.e. all functions
		SubFunctionGroup m_defaultSubFunctionGroup;

		// Mapping (fct, dof) -> local index
		std::vector<std::vector<dof_index_type> > m_vvDofIndices;

		// algebra indices
		std::vector<index_type> m_vAlgIndices;
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
		typedef TEntry entry_type;

	protected:
		typedef typename std::vector<entry_type>::iterator iterator;

	public:
		///////////////////////////
		// setup
		///////////////////////////

		LocalVector() : m_pIndices(NULL) {m_entries.clear();}

		LocalVector(const LocalIndices& ind){m_entries.clear(); set_indices(ind);}

		/// set new local indices
		void set_indices(const LocalIndices& ind)
		{
			if(ind.num_indices() != size()) m_entries.resize(ind.num_indices());
			m_pIndices = &ind;
		}

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
		entry_type& operator[](size_t i)
		{
			UG_ASSERT(i < size(), "Index not valid");
			return m_entries[i];
		}

		/// access to entry
		const entry_type& operator[](size_t i) const
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
		const LocalIndices* m_pIndices;

		// entries (size = m_indices.num_indices())
		std::vector<entry_type> m_entries;
};



template <typename TEntry>
class LocalMatrix
{
	public:
		typedef LocalMatrix<TEntry> this_type;

		// index type used by algebra
		typedef size_t index_type;

		// component type used by algebra
		typedef size_t comp_type;

		// entry type used by algebra
		typedef TEntry entry_type;

	protected:
		typedef typename std::vector<std::vector<entry_type> >::iterator   row_iterator;
		typedef typename std::vector<entry_type>::iterator col_iterator;

	public:
		///////////////////////////
		// setup
		///////////////////////////

		LocalMatrix() : m_pRowIndices(NULL), m_pColIndices(NULL) {}

		LocalMatrix(const LocalIndices& rowInd, const LocalIndices& colInd)
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
					m_entries[i][j] = mat(row_index(i), col_index(j));
		}

		///////////////////////////
		// matrix functions
		///////////////////////////

		/// number of rows
		size_t num_rows() const {return m_entries.size();}

		/// number of columns
		size_t num_cols() const
		{
			UG_ASSERT(num_rows() > 0, "No row set. Thus, now col size present");
			return m_entries[0].size();
		}

		/// access to entry
		entry_type& operator() (size_t i, size_t j)
		{
			UG_ASSERT(i < num_rows(), "Row does not exist.");
			UG_ASSERT(j < num_cols(), "Column does not exist.");
			return m_entries[i][j];
		}

		/// const access to entry
		const entry_type& operator() (size_t i, size_t j) const
		{
			UG_ASSERT(i < num_rows(), "Row does not exist.");
			UG_ASSERT(j < num_cols(), "Column does not exist.");
			return m_entries[i][j];
		}

		/// set all entries
		void set(number val)
		{
			for(size_t i = 0; i < num_rows(); ++i)
				for(size_t j = 0; j < num_cols(); ++j)
					m_entries[i][j] = val;
		}

		/// multiply all entries
		this_type& operator*(number val)
		{
			for(size_t i = 0; i < num_rows(); ++i)
				for(size_t j = 0; j < num_cols(); ++j)
					m_entries[i][j] *= val;
			return *this;
		}

		/// add matrix
		this_type& operator+=(const this_type& rhs)
		{
			UG_ASSERT(num_rows() == rhs.num_rows(), "Row size does not match");
			UG_ASSERT(num_cols() == rhs.num_cols(), "Column size does not match");

			for(size_t i = 0; i < num_rows(); ++i)
				for(size_t j = 0; j < num_cols(); ++j)
					m_entries[i][j] += rhs(i,j);
			return *this;
		}

		/// subtract matrix
		this_type& operator-=(const this_type& rhs)
		{
			UG_ASSERT(num_rows() == rhs.num_rows(), "Row size does not match");
			UG_ASSERT(num_cols() == rhs.num_cols(), "Column size does not match");

			for(size_t i = 0; i < num_rows(); ++i)
				for(size_t j = 0; j < num_cols(); ++j)
					m_entries[i][j] -= rhs(i,j);
			return *this;
		}

		///////////////////////////
		// algebra access
		///////////////////////////

		/// index of component
		index_type row_index(size_t i) const {return m_pRowIndices->index(i);}
		index_type col_index(size_t j) const {return m_pColIndices->index(j);}

		///////////////////////////
		// discretization access
		///////////////////////////

		number& operator()(size_t rowFct, size_t rowDof, size_t colFct, size_t colDof)
		{
			const index_type row_index = m_pRowIndices->local_index(rowFct, rowDof);
			const index_type col_index = m_pColIndices->local_index(colFct, colDof);
			const comp_type row_comp = m_pRowIndices->comp(rowFct, rowDof);
			const comp_type col_comp = m_pColIndices->comp(colFct, colDof);

			UG_ASSERT(row_index < num_rows(), "Row does not exist.");
			UG_ASSERT(col_index < num_cols(), "Column does not exist.");
			return BlockRef(m_entries[row_index][col_index], row_comp, col_comp);
		}

		/// const access to dof of function fct
		const number& operator()(size_t rowFct, size_t rowDof, size_t colFct, size_t colDof) const
		{
			const index_type row_index = m_pRowIndices->local_index(rowFct, rowDof);
			const index_type col_index = m_pColIndices->local_index(colFct, colDof);
			const comp_type row_comp = m_pRowIndices->comp(rowFct, rowDof);
			const comp_type col_comp = m_pColIndices->comp(colFct, colDof);

			UG_ASSERT(row_index < num_rows(), "Row does not exist.");
			UG_ASSERT(col_index < num_cols(), "Column does not exist.");
			return BlockRef(m_entries[row_index][col_index], row_comp, col_comp);
		}

	protected:
		void resize()
		{
			for(size_t i = 0; i < m_entries.size(); ++i)
			{
				m_entries[i].clear();
			}
			m_entries.clear();

			if(m_pRowIndices == NULL) return;
			m_entries.resize(m_pRowIndices->num_indices());

			if(m_pColIndices == NULL) return;
			for(size_t i = 0; i < m_pRowIndices->num_indices(); ++i)
			{
				m_entries[i].resize(m_pColIndices->num_indices());
			}
		}

	protected:
		// row indices
		const LocalIndices* m_pRowIndices;

		// column indices
		const LocalIndices* m_pColIndices;

		// entries
		std::vector<std::vector<entry_type> > m_entries;
};


} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__LOCAL_ALGEBRA__ */
