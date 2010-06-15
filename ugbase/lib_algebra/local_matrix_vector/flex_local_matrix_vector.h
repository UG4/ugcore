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

class FlexLocalMatrix{

	public:
		// type of value entry
		typedef number value_type;

	private:
		typedef std::vector<std::vector<number> >::iterator   row_iterator;
		typedef std::vector<number>::iterator col_iterator;

	public:
	FlexLocalMatrix()
	{
		m_values.clear();
	}

	FlexLocalMatrix(size_t nrow, size_t ncol)
	{
		m_values.resize(nrow);
		row_iterator row_iterEnd = m_values.end();
		for(row_iterator iter = m_values.begin(); iter != row_iterEnd; ++iter)
		{
			(*iter).resize(ncol);
		}
	}

	void resize(size_t nrow, size_t ncol)
	{
		m_values.resize(nrow);
		row_iterator row_iterEnd = m_values.end();
		for(row_iterator iter = m_values.begin(); iter != row_iterEnd; ++iter)
		{
			(*iter).resize(ncol);
		}
	}

	size_t num_rows() const
	{
		return m_values.size();
	}

	size_t num_cols() const
	{
		return m_values[0].size();
	}

	void set(number val)
	{
		for(size_t i = 0; i < m_values.size(); ++i)
		{
			for(size_t j = 0; j < m_values[i].size(); ++j)
			{
				m_values[i][j] = val;
			}
		}
	}

	FlexLocalMatrix& operator*(number val)
	{
		for(size_t i = 0; i < m_values.size(); ++i)
		{
			for(size_t j = 0; j < m_values[i].size(); ++j)
			{
				m_values[i][j] *= val;
			}
		}
		return *this;
	}

	FlexLocalMatrix& operator+=(const FlexLocalMatrix& rhs)
	{
		assert(m_values.size() == rhs.m_values.size());

		for(size_t i = 0; i < m_values.size(); ++i)
		{
			for(size_t j = 0; j < m_values[i].size(); ++j)
			{
				m_values[i][j] += rhs.m_values[i][j];
			}
		}
		return *this;
	}

	FlexLocalMatrix& operator-=(const FlexLocalMatrix& rhs)
	{
		assert(m_values.size() == rhs.m_values.size());

		for(size_t i = 0; i < m_values.size(); ++i)
		{
			for(size_t j = 0; j < m_values[i].size(); ++j)
			{
				m_values[i][j] -= rhs.m_values[i][j];
			}
		}
		return *this;
	}

	number& operator() (size_t i, size_t j)
	{
		assert(i < m_values.size());
		assert(j < m_values[i].size());
		return m_values[i][j];
	}

	const number& operator() (size_t i, size_t j) const
	{
		assert(i < m_values.size());
		assert(j < m_values[i].size());
		return m_values[i][j];
	}

	private:
		std::vector<std::vector<number> > m_values;
};



class FlexLocalVector{

	public:
	// type of value entry
	typedef number value_type;

	private:
		typedef std::vector<number>::iterator   iterator;

	public:
	FlexLocalVector()
	{
		m_values.clear();
	}

	FlexLocalVector(size_t nrow)
	{
		m_values.resize(nrow);
	}

	void resize(size_t nrow)
	{
		m_values.resize(nrow);
	}

	size_t size() const
	{
		return m_values.size();
	}

	void set(number val)
	{
		iterator iterEnd = m_values.end();
		for(iterator i = m_values.begin(); i != iterEnd; ++i)
		{
				(*i) = val;
		}
	}

	void push_back(number val)
	{
		m_values.push_back(val);
	}

	FlexLocalVector& operator*(number val)
	{
		iterator iterEnd = m_values.end();
		for(iterator i = m_values.begin(); i != iterEnd; ++i)
		{
				(*i) *= val;
		}
		return *this;
	}

	FlexLocalVector& operator+=(const FlexLocalVector& rhs)
	{
		assert(m_values.size() == rhs.m_values.size());

		for(size_t i = 0; i < m_values.size(); ++i)
		{
				m_values[i] += rhs.m_values[i];
		}
		return *this;
	}

	FlexLocalVector& operator-=(const FlexLocalVector& rhs)
	{
		assert(m_values.size() == rhs.m_values.size());

		for(size_t i = 0; i < m_values.size(); ++i)
		{
				m_values[i] -= rhs.m_values[i];
		}
		return *this;
	}

	number& operator[] (size_t i)
	{
		assert(i < m_values.size());
		return m_values[i];
	}

	const number& operator[] (size_t i) const
	{
		assert(i < m_values.size());
		return m_values[i];
	}


	private:
		std::vector<number> m_values;
};

inline std::ostream& operator<< (std::ostream& outStream, const ug::FlexLocalMatrix& m)
{
	for(size_t i = 0; i < m.num_rows(); ++i)
	{
		for(size_t j = 0; j < m.num_cols(); ++j)
		{
			outStream << "[" << i << ", " << j << "]: " <<  m(i,j) << std::endl;
		}
	}

 	return outStream;
}


}

#endif /* __H__LIB_ALGEBRA__LOCAL_MATRIX_VECTOR__FLEX_LOCAL_MATRIX_VECTOR__ */
