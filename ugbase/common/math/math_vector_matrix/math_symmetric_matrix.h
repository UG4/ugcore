/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Christian Wehner
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef SYM_MATRIX_H_
#define SYM_MATRIX_H_

#include <cstddef>
#include <ostream>
#include <iomanip>

#include "../../types.h"
#include "math_vector.h"
#include "common/assert.h"

namespace ug {

/**
 * \defgroup math_matrix Matrix
 * \ingroup ugbase_math
 * \{
 */

//template <std::size_t N, std::size_t M, typename T = number> class MathSymmetricMatrix;

/**
 * \class MathSymmetricMatrix
 *
 * \brief A class for fixed size, dense matrices.
 *
 *	A static memory NxN symmetric matrix
 */
template <std::size_t N, typename T = number>
class MathSymmetricMatrix
{
	public:
		using value_type = T;
		using size_type = std::size_t;
		static constexpr std::size_t RowSize = N;
		static constexpr std::size_t ColSize = N;

	public:
		MathSymmetricMatrix() {}
		MathSymmetricMatrix(const MathSymmetricMatrix& v)	{assign(v);}
		
		/**
		 * \brief Assigns the elements of the given matrix to this one.
		 *
		 * \param v The matrix to be assigned.
		 * \return A reference to this matrix.
		 */
		MathSymmetricMatrix& operator = (const MathSymmetricMatrix& v)
		{
			assign(v);
			return *this;
		}

		/**
		 * \brief Adds a matrix to 'this' one: \f$ A_{this} \leftarrow A_{this} + B\f$.
		 *
		 * \param B The matrix to be added.
		 * \return A reference to this matrix.
		 */
		MathSymmetricMatrix& operator += (const MathSymmetricMatrix& B)
		{
			for(std::size_t i = 0; i < m_size; ++i)
			{
				m_data[i] += B[i];
			}
			return *this;
		}

		/**
		 * \brief Subtracts a matrix from 'this' one: \f$ A_{this} \leftarrow A_{this} - B\f$.
		 *
		 * \param B The matrix to be subtracted.
		 * \return A reference to this matrix.
		 */
		MathSymmetricMatrix& operator -= (const MathSymmetricMatrix& B)
		{
			for(std::size_t i = 0; i < m_size; ++i)
			{
				m_data[i] -= B[i];
			}
			return *this;
		}

		/**
		 * \brief Assigns the given value to all elements of the matrix.
		 *
		 * \param val The value to be assigned to the matrix.
		 * \return A reference to this matrix.
		 */
		MathSymmetricMatrix& operator = (const value_type& val)
		{
			for(std::size_t i = 0; i < m_size; ++i)
			{
				m_data[i] = val;
			}
			return *this;
		}

		/**
		 * \brief Adds the given value to all elements of the matrix.
		 *
		 * \param val The value to be added.
		 * \return A reference to this matrix.
		 */
		MathSymmetricMatrix& operator += (const value_type& val)
		{
			for(std::size_t i = 0; i < m_size; ++i)
			{
				m_data[i] += val;
			}
			return *this;
		}

		/**
		 * \brief Subtracts the given value from all elements of the matrix.
		 *
		 * \param val The value to be subtracted.
		 * \return A reference to this matrix.
		 */
		MathSymmetricMatrix& operator -= (const value_type& val)
		{
			for(std::size_t i = 0; i < m_size; ++i)
			{
				m_data[i] -= val;
			}
			return *this;
		}

		/**
		 * \brief Divides all elements of the matrix by the given value.
		 *
		 * \param val The divisor.
		 * \return A reference to this matrix.
		 */
		MathSymmetricMatrix& operator /= (const value_type& val)
		{
			for(std::size_t i = 0; i < m_size; ++i)
			{
				m_data[i] /= val;
			}
			return *this;
		}

		/**
		 * \brief Multiplies all elements of the matrix with the given value.
		 *
		 * \param val The factor.
		 * \return A reference to this matrix.
		 */
		MathSymmetricMatrix& operator *= (const value_type& val)
		{
			for(std::size_t i = 0; i < m_size; ++i)
			{
				m_data[i] *= val;
			}
			return *this;
		}
		
		/**
		 * \brief Multiplies the matrix element-wise with another matrix and sums up the entries.
		 *
		 * \param v The Matrix.
		 * \return A scalar value of the element-wise summed up products
		 */
		value_type operator * (const MathSymmetricMatrix& v) const
		{
			value_type res = 0.0;
			for(std::size_t i = 0; i < N; ++i)
			{
				for(std::size_t j = 0; j < N; ++j)
				{
					res += this->entry(i,j) * v(i,j);
				}
			}
			return res;
		}

		inline std::size_t num_rows() const {return N;}
		inline std::size_t num_cols() const {return N;}

		inline value_type& entry(std::size_t row, std::size_t col)				
		{
			UG_ASSERT(row < N && col < N, "Accessing "<<N<<"x"<<N<<"Matrix at entry ("<<row<<","<<col<<")"); 
			if (row<col) return m_data[row * N - (row - 1) * row / 2 + col - row];
				else return m_data[col * N - (col - 1) * col / 2 + row - col];
		}
		
		inline const value_type& entry(std::size_t row, std::size_t col) const
		{
			UG_ASSERT(row < N && col < N, "Accessing "<<N<<"x"<<N<<"Matrix at entry ("<<row<<","<<col<<")"); 
			if (row<col) return m_data[row * N - (row - 1) * row / 2 + col - row];
				else return m_data[col * N - (col - 1) * col / 2 + row - col];
		}

		inline value_type& operator [] (std::size_t index)				{UG_ASSERT(index < m_size, "Invalid index"); return m_data[index];}
		inline const value_type& operator [] (std::size_t index) const	{UG_ASSERT(index < m_size, "Invalid index"); return m_data[index];}

		inline value_type& operator () (std::size_t row, std::size_t col)
		{
			UG_ASSERT(row < N && col < N, "Accessing "<<N<<"x"<<N<<"Matrix at entry ("<<row<<","<<col<<")"); 
			if (row<col) return m_data[row * N - (row - 1) * row / 2 + col - row];
				else return m_data[col * N - (col - 1) * col / 2 + row - col];
		}
		
		inline const value_type& operator () (std::size_t row, std::size_t col) const
		{
			UG_ASSERT(row < N && col < N, "Accessing "<<N<<"x"<<N<<"Matrix at entry ("<<row<<","<<col<<")"); 
			if (row<col) return m_data[row * N - (row - 1) * row / 2 + col - row];
				else return m_data[col * N - (col - 1) * col / 2 + row - col];
		}

		// frobenius norm |A|_{F} = \sqrt{ 2 \sum_{i,j} A_{i,j} } 
		inline value_type fnorm(){
			value_type norm=m_data[0]*m_data[0];
			for (size_t i=1;i<m_size-1;i++){
				norm += 2*m_data[i]*m_data[i];
			}
			norm+=m_data[m_size-1]*m_data[m_size-1];
			if (N>2)
				for (size_t diag=1;diag<N-1;diag++)
					norm -= m_data[diag * N - (diag - 1) * diag / 2]*m_data[diag * N - (diag - 1) * diag / 2];
			return sqrt(2.0*norm);
		}
		
		// scale matrix by frobenius norm
		void scale_by_fnorm(){
			number norm=fnorm();
			for (size_t i=0;i<m_size;i++) m_data[i]*=norm;
		}
		
	protected:
		static constexpr size_t m_size = (N*N+N)/2;
		
		value_type m_data[m_size];

		inline void assign(const MathSymmetricMatrix& v)
		{
			for(std::size_t i = 0; i < m_size; ++i)
			{
				m_data[i] = v[i];
			}
		}
};
	

/// Print MathSymmetricMatrix<N> to standard output
template <std::size_t N>
std::ostream& operator << (std::ostream& outStream, const MathSymmetricMatrix<N>& m)
{
	for(std::size_t i = 0; i < N; ++i)
	{
		for(std::size_t j = 0; j < N; ++j)
		{
			outStream << "[" << i << "][" << j << "]: " << std::scientific << std::setprecision(8) << std::setw(15) << m.entry(i, j) << std::endl;
		}
	}
	return outStream;
}

std::ostream& operator << (std::ostream& outStream, const MathSymmetricMatrix<2>& m);
std::ostream& operator << (std::ostream& outStream, const MathSymmetricMatrix<3>& m);


} //end of namespace: lgmath


#endif