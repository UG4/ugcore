/*
 * flex_local_matrix.h
 *
 *  Created on: 04.03.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_ALGEBRA__LOCAL_MATRIX_VECTOR__FLEX_LOCAL_MATRIX_VECTOR__
#define __H__LIB_ALGEBRA__LOCAL_MATRIX_VECTOR__FLEX_LOCAL_MATRIX_VECTOR__

#include <vector>

#include "common/common.h"

namespace ug{

/** Flex Local Matrix
 *
 * This is a Matrix of small range, that has components of TValue.
 * The value_type must provide the following operators:
 * 	- operator=
 *  - operator*
 * 	- operator-=
 *  - operator+=
 */
template <typename TValue>
class FlexLocalMatrix{

	public:
		// type of value entry
		typedef TValue entry_type;

		// own type
		typedef FlexLocalMatrix<TValue> this_type;

	private:
		typedef typename std::vector<std::vector<entry_type> >::iterator   row_iterator;
		typedef typename std::vector<entry_type>::iterator col_iterator;

	public:
		/// constructor for empty matrix
		FlexLocalMatrix(){m_values.clear();}

		/// constructor for matrix of size nrow x ncol
		FlexLocalMatrix(size_t nrow, size_t ncol)
		{
			m_values.resize(nrow);
			row_iterator row_iterEnd = m_values.end();
			for(row_iterator iter = m_values.begin(); iter != row_iterEnd; ++iter)
				(*iter).resize(ncol);
		}

		/// resize matrix
		void resize(size_t nrow, size_t ncol)
		{
			m_values.resize(nrow);
			row_iterator row_iterEnd = m_values.end();
			for(row_iterator iter = m_values.begin(); iter != row_iterEnd; ++iter)
				(*iter).resize(ncol);
		}

		/// number of rows
		size_t num_rows() const {return m_values.size();}

		/// number of columns
		size_t num_cols() const {return m_values[0].size();}

		/// set all entries
		void set(entry_type val)
		{
			for(size_t i = 0; i < num_rows(); ++i)
				for(size_t j = 0; j < num_cols(); ++j)
					m_values[i][j] = val;
		}

		/// multiply all entries
		this_type& operator*(entry_type val)
		{
			for(size_t i = 0; i < num_rows(); ++i)
				for(size_t j = 0; j < num_cols(); ++j)
					m_values[i][j] *= val;
			return *this;
		}

		/// add matrix
		this_type& operator+=(const this_type& rhs)
		{
			UG_ASSERT(num_rows() == rhs.num_rows(), "Row size does not match");
			UG_ASSERT(num_cols() == rhs.num_cols(), "Column size does not match");

			for(size_t i = 0; i < num_rows(); ++i)
				for(size_t j = 0; j < num_cols(); ++j)
					m_values[i][j] += rhs(i,j);
			return *this;
		}

		/// subtract matrix
		this_type& operator-=(const this_type& rhs)
		{
			UG_ASSERT(num_rows() == rhs.num_rows(), "Row size does not match");
			UG_ASSERT(num_cols() == rhs.num_cols(), "Column size does not match");

			for(size_t i = 0; i < num_rows(); ++i)
				for(size_t j = 0; j < num_cols(); ++j)
					m_values[i][j] -= rhs(i,j);
			return *this;
		}

		/// access to entry
		entry_type& operator() (size_t i, size_t j)
		{
			UG_ASSERT(i < num_rows(), "Row does not exist.");
			UG_ASSERT(j < num_cols(), "Column does not exist.");
			return m_values[i][j];
		}

		/// const access to entry
		const entry_type& operator() (size_t i, size_t j) const
		{
			UG_ASSERT(i < num_rows(), "Row does not exist.");
			UG_ASSERT(j < num_cols(), "Column does not exist.");
			return m_values[i][j];
		}

	private:
		std::vector<std::vector<number> > m_values;
};


/** Flex Local Vector
 *
 * This is a Vector of small range, that has components of TValue.
 * The value_type must provide the following operators:
 * 	- operator=
 *  - operator*
 * 	- operator-=
 *  - operator+=
 */
template <typename TValue>
class FlexLocalVector{

	public:
		// type of value entry
		typedef TValue entry_type;

		// own type
		typedef FlexLocalVector<TValue> this_type;

	private:
		typedef typename std::vector<entry_type>::iterator   iterator;

	public:
		/// constructor for vector of size 0
		FlexLocalVector(){m_values.clear();}

		/// constructor of vector of size nrow
		FlexLocalVector(size_t size){m_values.resize(size);}

		/// resize the vector, entries 0, ..., newSize-1 remain valid
		void resize(size_t newSize){m_values.resize(newSize);}

		/// size of vector
		size_t size() const {return m_values.size();}

		/// set all components of the vector
		void set(entry_type val)
		{
			iterator iterEnd = m_values.end();
			for(iterator iter = m_values.begin(); iter != iterEnd; ++iter)
				(*iter) = val;
		}

		/// multiply all components of the vector
		this_type& operator*(entry_type val)
		{
			iterator iterEnd = m_values.end();
			for(iterator iter = m_values.begin(); iter != iterEnd; ++iter)
				(*iter) *= val;
			return *this;
		}

		/// add a vector
		this_type& operator+=(const this_type& rhs)
		{
			UG_ASSERT(size() == rhs.size(), "Size does not match.");
			for(size_t i = 0; i < size(); ++i)
				m_values[i] += rhs[i];
			return *this;
		}

		/// subtract a vector
		this_type& operator-=(const this_type& rhs)
		{
			UG_ASSERT(size() == rhs.size(), "Size does not match.");
			for(size_t i = 0; i < m_values.size(); ++i)
					m_values[i] -= rhs[i];
			return *this;
		}

		/// access to component
		entry_type& operator[] (size_t i)
		{
			UG_ASSERT(i < size(), "Out of range.");
			return m_values[i];
		}

		/// const access to component
		const entry_type& operator[] (size_t i) const
		{
			UG_ASSERT(i < size(), "Out of range.");
			return m_values[i];
		}

		// depreciated
		void push_back(const entry_type& v) {m_values.push_back(v);}

	private:
		// values
		std::vector<entry_type> m_values;
};

template <typename TValue>
inline std::ostream& operator<< (std::ostream& outStream, const typename ug::FlexLocalMatrix<TValue>& m)
{
	for(size_t i = 0; i < m.num_rows(); ++i)
		for(size_t j = 0; j < m.num_cols(); ++j)
			outStream << "[" << i << ", " << j << "]: " <<  m(i,j) << std::endl;
 	return outStream;
}

template <typename TValue>
inline std::ostream& operator<< (std::ostream& outStream, const typename ug::FlexLocalVector<TValue>& v)
{
	for(size_t i = 0; i < v.size(); ++i)
		outStream << "[" << i << "]: " <<  v[i] << std::endl;
 	return outStream;
}


}

#endif /* __H__LIB_ALGEBRA__LOCAL_MATRIX_VECTOR__FLEX_LOCAL_MATRIX_VECTOR__ */
