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
 * function (0 <= fct <= \#Selected Functions - 1) and local dof indices (0 <= dof <= \#DoFs -1).
 * The access from algebra side is only an access to the algebra i-index. Using
 *  - local_index(size_t fct, size_t dof)
 *  - comp(size_t fct, size_t dof)
 *  will return the correspinding (i,alpha) entry for a (fct, dof) pair.
 */
class LocalIndices
{
	public:
	//	Index type used by algebra
		typedef size_t index_type;

	// 	Component type used by algebra
		typedef size_t comp_type;

	//  Multi Indix type used by discretization
		typedef MultiIndex<2> multi_index_type;

	public:
	//	Constructor
		LocalIndices() :
			m_pFuncGroup(NULL), m_pFuncMap(NULL)
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
			m_pFuncGroup = &funcGroup;

		// resize vectors
			resize_dof_indices();

		//	update identity map
			for(size_t i = 0; i < m_pFuncGroup->num_fct(); ++i)
				m_IdentityMap.add(i,i);

		//	access all by default
			access_all();
		}

	/// access only part of the function group using mapping
		void access_by_map(const FunctionIndexMapping& funcMap){m_pFuncMap = &funcMap;}

	/// access all functions in FunctionGroup
		void access_all(){m_pFuncMap = &m_IdentityMap;}

	/// get current mapping
		const FunctionIndexMapping* get_map() const {return m_pFuncMap;}

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
			UG_ASSERT(i < num_indices(), "Algebra index not valid. (i="<<i<<", "
			          	  	  	  	  	  "num_indices="<<num_indices()<<")");
			m_vAlgIndices[i] = ind;
		}

		///////////////////////////
		// Discretization access
		///////////////////////////

	/// number of accessible functions
		size_t num_fct() const {return m_pFuncMap->num_fct();}

	/// return unique function id
		size_t unique_id(size_t fct) const
		{
			UG_ASSERT(fct < num_fct(), "Local fct index '"<<fct<<"' not valid."
			          	  	  	  	   " (num_fct = "<<num_fct()<<")");
			const size_t accFct = (*m_pFuncMap)[fct];
			return m_pFuncGroup->unique_id(accFct);
		}

	/// number of dofs of all accessible (sum)
		size_t num_dofs() const
		{
			size_t num = 0;
			for(size_t fct = 0; fct < num_fct(); ++fct)
				num += num_dofs(fct);
			return num;
		}

	/// number of dofs for accessible function
		size_t num_dofs(size_t fct) const
		{
			UG_ASSERT(fct < num_fct(), "local fct index '"<<fct<<"' not valid."
			          	  	  	  	   " (num_fct = "<<num_fct()<<")");
			const size_t accFct = (*m_pFuncMap)[fct];
			UG_ASSERT(accFct < m_vvDofIndices.size(), "accFct index "<<accFct<<
			          	  " not valid. (Size = "<<m_vvDofIndices.size()<<")\n");
			return m_vvDofIndices[accFct].size();
		}

	/// global algebra index for (fct, dof)
		index_type global_index(size_t fct, size_t dof) const{return index(local_index(fct, dof));}

	/// local algebra index for (fct, dof)
		index_type local_index(size_t fct, size_t dof) const {return local_multi_index(fct,dof)[0];}

	/// algebra comp for (fct, dof)
		comp_type comp(size_t fct, size_t dof) const {return local_multi_index(fct,dof)[1];}

	///	return local multi index for (fct, dof)
		const multi_index_type& local_multi_index(size_t fct, size_t dof) const
		{
			UG_ASSERT(fct < num_fct(), "local fct index '"<< fct <<
			          	  	  	  "' not valid. (num_fct = "<<num_fct()<<")");
			const size_t accFct = (*m_pFuncMap)[fct];
			UG_ASSERT(accFct < m_vvDofIndices.size(), "accFct index " << accFct
			          << " not valid. (Size = " << m_vvDofIndices.size() << ")\n");
			UG_ASSERT(dof < m_vvDofIndices[accFct].size(), "Local dof "<< dof <<
			          " not valid. (Size = "<<m_vvDofIndices[accFct].size()<<")");
			return m_vvDofIndices[accFct][dof];
		}

		///////////////////////////
		// Algebra access
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
			m_vvDofIndices.resize(m_pFuncGroup->num_fct());
			for(size_t fct = 0; fct < m_vvDofIndices.size(); ++fct)
			{
				m_vvDofIndices[fct].clear();
			}
		}

	protected:
	///	underlying function group
		const FunctionGroup* m_pFuncGroup;

	///	default Function map, accessing all
		FunctionIndexMapping m_IdentityMap;

	/// flag iff all functions are accessed
		bool m_bAccessAll;

	/// Access Mapping
		const FunctionIndexMapping* m_pFuncMap;

	// 	Mapping (fct, dof) -> local index
		std::vector<std::vector<multi_index_type> > m_vvDofIndices;

	// 	Algebra indices
		std::vector<index_type> m_vAlgIndices;
};

template <typename TEntry>
class LocalVector
{
	public:
	///	own type
		typedef LocalVector<TEntry> this_type;

	/// Index type used by algebra
		typedef size_t index_type;

	// 	Component type used by algebra
		typedef size_t comp_type;

	// 	Entry type used by algebra
		typedef TEntry value_type;

	protected:
		typedef typename std::vector<value_type>::iterator iterator;

	public:
	///	Constructor
		LocalVector() : m_pIndex(NULL) {m_entries.clear();}

	///	Constructor
		LocalVector(LocalIndices& ind)
			: m_pIndex(NULL)
		{m_entries.clear(); set_indices(ind);}

		///////////////////////////
		// setup
		///////////////////////////

	/// set new local indices
		void set_indices(LocalIndices& ind)
		{
			if(ind.num_indices() != size()) m_entries.resize(ind.num_indices());
			m_pIndex = &ind;
		}

	/// get current local indices
		LocalIndices& get_indices() {return *m_pIndex;}

	///	read values from a vector
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

	/// const access to entry
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

	/// multiply all components of the vector
		this_type& operator*=(number val)
		{
			iterator iterEnd = m_entries.end();
			for(iterator iter = m_entries.begin(); iter != iterEnd; ++iter)
				(*iter) *= val;
			return *this;
		}

	/// add a local vector
		this_type& operator+=(const this_type& rhs)
		{
			UG_ASSERT(size() == rhs.size(), "Size does not match.");
			for(size_t i = 0; i < size(); ++i)
				m_entries[i] += rhs[i];
			return *this;
		}

	/// subtract a local vector
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

	/// global index of component
		index_type index(size_t i) const {return m_pIndex->index(i);}

		///////////////////////////
		// discretization access
		///////////////////////////

	/// number of accessible functions
		size_t num_fct() const {return m_pIndex->num_fct();}

	/// return global function id for accessible functions
		size_t unique_id(size_t fct) const	{return m_pIndex->unique_id(fct);}

	/// number of dofs of accessible functions (sum)
		size_t num_dofs() const {return m_pIndex->num_dofs();}

	/// number of dofs per accessible function
		size_t num_dofs(size_t fct) const {return m_pIndex->num_dofs(fct);}

	/// access to dof of accessible function fct
		number& operator()(size_t fct, size_t dof)
		{
			const typename LocalIndices::multi_index_type& ind
										= m_pIndex->local_multi_index(fct, dof);
			UG_ASSERT(ind[0] < size(), "Index not valid");
			return BlockRef(m_entries[ind[0]], ind[1]);
		}

	/// const access to dof of accessible function fct
		const number& operator()(size_t fct, size_t dof) const
		{
			const typename LocalIndices::multi_index_type& ind
										= m_pIndex->local_multi_index(fct, dof);
			UG_ASSERT(ind[0] < size(), "Index not valid");
			return BlockRef(m_entries[ind[0]], ind[1]);
		}

	protected:
	/// Indices
		LocalIndices* m_pIndex;

	/// Entries
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

	/// Index type used by algebra
		typedef size_t index_type;

	//  Component type used by algebra
		typedef size_t comp_type;

	// 	Entry type used by algebra
		typedef TEntry value_type;

	protected:
		typedef typename std::vector<std::vector<value_type> >::iterator   row_iterator;
		typedef typename std::vector<value_type>::iterator col_iterator;

	public:
	///	Constructor
		LocalMatrix() : m_pRowIndex(NULL), m_pColIndex(NULL) {}

	///	Constructor
		LocalMatrix(const LocalIndices& rowInd, const LocalIndices& colInd)
			: m_pRowIndex(NULL), m_pColIndex(NULL)
		{
			set_row_indices(rowInd);
			set_col_indices(colInd);
		}

		///////////////////////////
		// Setup
		///////////////////////////

	/// set new local row indices
		void set_row_indices(const LocalIndices& ind)
		{
			m_pRowIndex = &ind;
			resize();
		}

	/// set new local column indices
		void set_col_indices(const LocalIndices& ind)
		{
			m_pColIndex = &ind;
			resize();
		}

	/// set new local column indices but do not resize
		void set_col_indices_no_resize(const LocalIndices& ind)
		{
			m_pColIndex = &ind;
		}

	/// set new local indices
		void set_indices(const LocalIndices& rowInd, const LocalIndices& colInd)
		{
			set_row_indices(rowInd);
			set_col_indices(colInd);
		}

	///	read entries from a matrix
		template <typename TMatrix>
		void read_values(const TMatrix& mat)
		{
			for(size_t i = 0; i < num_rows(); ++i)
				for(size_t j = 0; j < num_cols(); ++j)
					m_entries(i,j) = mat(row_index(i), col_index(j));
		}

		///////////////////////////
		// Matrix functions
		///////////////////////////

	/// number of rows
		size_t num_rows() const	{return m_entries.num_rows();}

	/// number of columns
		size_t num_cols() const{return m_entries.num_cols();}

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
		void set(number val){m_entries = val;}

	/// multiply all entries
		this_type& operator*(number val)
		{
			m_entries *= val;
			return *this;
		}

	/// multiply matrix
		this_type& operator*=(number val)
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

	/// global algebra index of local algebra component
		index_type row_index(size_t i) const {return m_pRowIndex->index(i);}

	/// global algebra index of local algebra component
		index_type col_index(size_t j) const {return m_pColIndex->index(j);}

	/// algebra component index for row (fct,dof)
		comp_type row_sub_index(size_t rowFct, size_t rowDoF) const
			{return m_pRowIndex->comp(rowFct, rowDoF);}

	/// algebra component index for column (fct,dof)
		comp_type col_sub_index(size_t colFct, size_t colDoF) const
			{return m_pColIndex->comp(colFct, colDoF);}

		///////////////////////////
		// Discretization access
		///////////////////////////

	/// access to coupling (rowFct, rowDoF) x (colFct, colDoF)
		number& operator()(size_t rowFct, size_t rowDoF, size_t colFct, size_t colDoF)
		{
			const typename LocalIndices::multi_index_type& rowInd
				= m_pRowIndex->local_multi_index(rowFct, rowDoF);
			const typename LocalIndices::multi_index_type& colInd
				= m_pColIndex->local_multi_index(colFct, colDoF);
			UG_ASSERT(rowInd[0] < num_rows(), "Row does not exist.");
			UG_ASSERT(colInd[0] < num_cols(), "Column does not exist.");

			return BlockRef(m_entries(rowInd[0],colInd[0]), rowInd[1], colInd[1]);
		}

	/// const access to coupling (rowFct, rowDoF) x (colFct, colDoF)
		const number& operator()(size_t rowFct, size_t rowDoF, size_t colFct, size_t colDoF) const
		{
			const typename LocalIndices::multi_index_type& rowInd
				= m_pRowIndex->local_multi_index(rowFct, rowDoF);
			const typename LocalIndices::multi_index_type& colInd
				= m_pColIndex->local_multi_index(colFct, colDoF);
			UG_ASSERT(rowInd[0] < num_rows(), "Row does not exist.");
			UG_ASSERT(colInd[0] < num_cols(), "Column does not exist.");

			return BlockRef(m_entries(rowInd[0],colInd[0]), rowInd[1], colInd[1]);
		}

	protected:
		void resize()
		{
			m_entries.resize(0,0);

			if(m_pRowIndex == NULL || m_pColIndex == NULL) return;

			const size_t num_rows = m_pRowIndex->num_indices();
			const size_t num_cols = m_pColIndex->num_indices();

			m_entries.resize(num_rows, num_cols);
		}

	protected:
	// 	Row indices
		const LocalIndices* m_pRowIndex;

	//	Column indices
		const LocalIndices* m_pColIndex;

	// 	Entries of local matrix
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

template <typename Entry>
std::ostream& operator<< (std::ostream& outStream, const ug::LocalVector<Entry>& vec)
{
	for(size_t i = 0; i < vec.size(); ++i)
	{
		outStream << "["<<i<<"]: (" << vec.index(i) << ") : " << vec[i] << "\n";
	}
	return outStream;
}

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__COMMON__LOCAL_ALGEBRA__ */
