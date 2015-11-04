#ifndef MATRIX_H_
#define MATRIX_H_

#include <cstddef>
#include <ostream>
#include <iomanip>
#include "../../types.h"
#include "math_vector.h"
#include "common/assert.h"

namespace ug
{

/**
 * \defgroup math_matrix Matrix
 * \ingroup ugbase_math
 * \{
 */

template <std::size_t N, std::size_t M, typename T = number> class MathMatrix;

/**
 * \class MathMatrix
 *
 * \brief A class for fixed size, dense matrices.
 *
 *	A static memory NxM matrix
 */
template <std::size_t N, std::size_t M, typename T>
class MathMatrix
{
	public:
		typedef T value_type;
		typedef std::size_t size_type;
		static const std::size_t RowSize = N;
		static const std::size_t ColSize = M;

	public:
		MathMatrix() {}
		MathMatrix(const MathMatrix& v)	{assign(v);}

		/**
		 * \brief Assigns the elements of the given matrix to this one.
		 *
		 * \param v The matrix to be assigned.
		 * \return A reference to this matrix.
		 */
		MathMatrix& operator=  (const MathMatrix& v)
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
		MathMatrix& operator+= (const MathMatrix& B)
		{
			for(std::size_t i = 0; i < N; ++i){
				for(std::size_t j = 0; j < M; ++j){
					m_data[i][j] += B(i,j);
				}
			}
			return *this;
		}

		/**
		 * \brief Subtracts a matrix from 'this' one: \f$ A_{this} \leftarrow A_{this} - B\f$.
		 *
		 * \param B The matrix to be subtracted.
		 * \return A reference to this matrix.
		 */
		MathMatrix& operator-= (const MathMatrix& B)
		{
			for(std::size_t i = 0; i < N; ++i){
				for(std::size_t j = 0; j < M; ++j){
					m_data[i][j] -= B(i,j);
				}
			}
			return *this;
		}

		/**
		 * \brief Assigns the given value to all elements of the matrix.
		 *
		 * \param val The value to be assigned to the matrix.
		 * \return A reference to this matrix.
		 */
		MathMatrix& operator= (const value_type& val)
		{
			for(std::size_t i = 0; i < N; ++i){
				for(std::size_t j = 0; j < M; ++j){
					m_data[i][j] = val;
				}
			}
			return *this;
		}

		/**
		 * \brief Adds the given value to all elements of the matrix.
		 *
		 * \param val The value to be added.
		 * \return A reference to this matrix.
		 */
		MathMatrix& operator+= (const value_type& val)
		{
			for(std::size_t i = 0; i < N; ++i){
				for(std::size_t j = 0; j < M; ++j){
					m_data[i][j] += val;
				}
			}
			return *this;
		}

		/**
		 * \brief Subtracts the given value from all elements of the matrix.
		 *
		 * \param val The value to be subtracted.
		 * \return A reference to this matrix.
		 */
		MathMatrix& operator-= (const value_type& val)
		{
			for(std::size_t i = 0; i < N; ++i){
				for(std::size_t j = 0; j < M; ++j){
					m_data[i][j] -= val;
				}
			}
			return *this;
		}

		/**
		 * \brief Divides all elements of the matrix by the given value.
		 *
		 * \param val The divisor.
		 * \return A reference to this matrix.
		 */
		MathMatrix& operator/= (const value_type& val)
		{
			for(std::size_t i = 0; i < N; ++i){
				for(std::size_t j = 0; j < M; ++j){
					m_data[i][j] /= val;
				}
			}
			return *this;
		}

		/**
		 * \brief Multiplies all elements of the matrix with the given value.
		 *
		 * \param val The factor.
		 * \return A reference to this matrix.
		 */
		MathMatrix& operator*= (const value_type& val)
		{
			for(std::size_t i = 0; i < N; ++i){
				for(std::size_t j = 0; j < M; ++j){
					m_data[i][j] *= val;
				}
			}
			return *this;
		}

		/**
		 * \brief Multiplies the matrix element-wise with another matrix and sums up the entries.
		 *
		 * \param v The Matrix.
		 * \return A scalar value of the element-wise summed up products
		 */
		value_type operator* (const MathMatrix& v) const
		{
			value_type res = 0.0;
			for(std::size_t i = 0; i < N; ++i){
				for(std::size_t j = 0; j < M; ++j){
					res += m_data[i][j] * v.m_data[i][j];
				}
			}
			return res;
		}

		inline std::size_t num_rows() const {return N;}
		inline std::size_t num_cols() const {return M;}

		inline value_type* operator[](std::size_t index){
			UG_ASSERT(index < N, "Invalid index");
			return m_data[index];
		}
		inline const value_type* operator[](std::size_t index) const{
			UG_ASSERT(index < N, "Invalid index");
			return m_data[index];
		}

		inline value_type& operator() (std::size_t row, std::size_t col){
			UG_ASSERT(row < N && col < M, "Accessing "<<N<<"x"<<M<<"Matrix at entry ("<<row<<","<<col<<")");
			return m_data[row][col];
		}
		inline const value_type& operator() (std::size_t row, std::size_t col) const{
			UG_ASSERT(row < N && col < M, "Accessing "<<N<<"x"<<M<<"Matrix at entry ("<<row<<","<<col<<")");
			return m_data[row][col];
		}

		inline value_type& entry(std::size_t row, std::size_t col)				{return (*this)(row,col);}
		inline const value_type& entry(std::size_t row, std::size_t col) const	{return (*this)(row,col);}

		inline void assign(const MathVector<N, value_type>& vec, const std::size_t row) {
			UG_ASSERT(vec.Size == N, "Wrong vector size");
			for(std::size_t j = 0; j < N; j++)
				m_data[row][j] = vec[j];
		}

	protected:
		value_type m_data[N][M];

		inline void assign(const MathMatrix& v)
		{
			for(std::size_t i = 0; i < N; ++i){
				for(std::size_t j = 0; j < M; ++j){
					m_data[i][j] = v.m_data[i][j] ;
				}
			}
		}

};

// this are explicit instantiations to avoid compiler errors on XL
// {
template <typename T> class MathMatrix<0,0,T>{
	public:
		typedef T value_type;
		typedef std::size_t size_type;
};
template <std::size_t N, typename T> class MathMatrix<N,0,T>{
	public:
		typedef T value_type;
		typedef std::size_t size_type;
};
template <std::size_t N, typename T> class MathMatrix<0,N,T>{
	public:
		typedef T value_type;
		typedef std::size_t size_type;
};
// }

/// Print MathMatrix<N,M> to standard output
template <std::size_t N, std::size_t M>
std::ostream& operator<< (std::ostream& outStream, const ug::MathMatrix<N,M>& m)
{
	for(std::size_t i = 0; i < N; ++i)
	{
		for(std::size_t j = 0; j < M; ++j)
		{
			outStream << "[" << i << "][" << j << "]: " << std::scientific << std::setprecision(8) << std::setw(15) << m.entry(i, j) << std::endl;
		}
	}
	return outStream;
}

std::ostream& operator<< (std::ostream& outStream, const ug::MathMatrix<2,2>& m);
std::ostream& operator<< (std::ostream& outStream, const ug::MathMatrix<2,3>& m);
std::ostream& operator<< (std::ostream& outStream, const ug::MathMatrix<3,2>& m);
std::ostream& operator<< (std::ostream& outStream, const ug::MathMatrix<3,3>& m);

// end group math_matrix
/// \}

} //end of namespace: lgmath


#endif /* MathMatrix_H_ */
