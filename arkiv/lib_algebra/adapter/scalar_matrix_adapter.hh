/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Arne Nägel
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

#ifndef SCALAR_MATRIX_ADAPTER_HH_
#define SCALAR_MATRIX_ADAPTER_HH_

#include <cstdlib>

#include "lib_algebra/lib_algebra.h"
#include "lib_algebra/lib_algebra_impl.h"
#include "lib_algebra/cpu_algebra_types.h"

#include "lib_algebra/common/matrixio/matrix_io.h"
#include "lib_algebra/common/matrixio/matrix_io_mtx.h"

#include "common/assert.h"
using namespace ug;

// provides an interface for matrix of algebra type B for matrices originally of algebra type A
// allows to access a CPUBlockAlgebra (AT) as a scalar CPUAlgebra (ST)


template<typename MT>
inline void TruncateOffDiag(MT &A, size_t ncmp)
{ UG_ASSERT(0, "Implement!"); }

template<>
inline void TruncateOffDiag(CPUAlgebra::matrix_type &A, size_t ncmp)
{
	using matrix_type = CPUAlgebra::matrix_type;

	// for each row
	for (size_t i=0; i<A.num_rows(); ++i)
	{
		// for each column
		for (matrix_type::row_iterator matij = A.begin_row(i);
				matij != A.end_row(i); ++matij)
		{
			size_t j = matij.index();
			if ((i%ncmp) != (j%ncmp))
			{
			//	std::cerr << "Eliminating " << i << ", " << j<< std::endl;
				matij.value() = 0.0;
			}
		}
	}

}

template<typename AT, typename ST=CPUAlgebra>
class ScalarMatrixAdapter{

public:
	using encapsulated_matrix_type = typename AT::matrix_type;
	using value_type = typename ST::matrix_type::value_type;
	static constexpr int blockSize = AT::blockSize;

	// using const_row_iterator = typename ST::matrix_type::const_row_iterator ;

	ScalarMatrixAdapter(encapsulated_matrix_type& mat) : m_src(mat), m_const(mat) {};

	void resize_and_clear(size_t newRows, size_t newCols)
	{ m_src.resize_and_clear(newRows/blockSize, newCols/blockSize);}

	bool resize_and_keep_values(size_t newRows, size_t newCols)
	{return m_src.resize_and_keep_values(newRows/blockSize, newCols/blockSize);}


	value_type &operator () (size_t r, size_t c)
	{
		UG_ASSERT(r < num_rows(), "row index is too large");
		UG_ASSERT(c < num_cols(), "col index is too large");
		return BlockRef(m_src(r/blockSize, c/blockSize), r%blockSize, c%blockSize);
	};

	const value_type &operator () (size_t r, size_t c)  const
	{
		UG_ASSERT(r < num_rows(), "row index is too large");
		UG_ASSERT(c < num_cols(), "col index is too large");
		return BlockRef(m_src(r/blockSize, c/blockSize), r%blockSize, c%blockSize);
	}

	//! returns number of rows
	size_t num_rows() const
	{ return m_src.num_rows()*blockSize; }

	//! returns the number of cols
	size_t num_cols() const
	{ return m_src.num_cols()*blockSize; }

	//! returns the total number of connections
	size_t total_num_connections() const
	{ return m_src.total_num_connections()*blockSize*blockSize; }

	//! print (underlying) matrix
	void print(const char *text) const {m_src.print(text);}

	//! print (block) row of (underlying) matrix
	void printrow(size_t row) const {m_src.printrow(row/blockSize);}

	//! operator overloading for streams
	friend std::ostream& operator << (std::ostream& os, ScalarMatrixAdapter<AT,ST> const &a)
	{ a.outputToStream(os); return os; }

	/**
	 *  row_iterator
	 *  iterator over a row
	 */


	class row_iterator
	{
		typename encapsulated_matrix_type::row_iterator iter;
	public:
		inline void check() const {iter.check(); }
		row_iterator(typename encapsulated_matrix_type::row_iterator _iter)
		: iter(_iter) {}
		~row_iterator() = default;
		row_iterator *operator ->() { return iter.operator -> (); }
		bool operator != (const row_iterator &o) const { return *iter != o->iter;  }
		void operator ++ () { ++iter; }
		void operator += (int nr) { iter+=nr; }
		bool operator == (const row_iterator &other) const { return other->iter == *iter;}
		size_t index() const { return iter.index(); }
		value_type &value() { return BlockRef(iter.value(), blockSize, blockSize); }
	};

	class const_row_iterator
	{
		typename encapsulated_matrix_type::const_row_iterator iter;
	    public:
			inline void check() const {iter.check(); }
			const_row_iterator(typename encapsulated_matrix_type::const_row_iterator _iter)
			: iter(_iter) {}
			~const_row_iterator() = default;
			const_row_iterator *operator ->() { return iter.operator -> (); }
			bool operator != (const const_row_iterator &o) const { return iter!= o.iter;  }
			void operator ++ () { ++iter; }
			void operator += (int nr) { iter+=nr; }
			bool operator == (const const_row_iterator &other) const { return other.iter == iter;}
			size_t index() const { return iter.index(); }
			const value_type &value() const { return BlockRef(iter.value(), blockSize, blockSize); }
	  };


	row_iterator begin_row(size_t r)
	{ return row_iterator(m_src.begin_row(r)); }
	row_iterator end_row(size_t r)
	{ return row_iterator(m_src.end_row(r)); }

	const_row_iterator begin_row(size_t r) const
	{ return const_row_iterator(m_const.begin_row(r)); }
	const_row_iterator   end_row(size_t r) const
	{ return const_row_iterator(m_const.end_row(r)); }

protected:
	std::ostream& outputToStream(std::ostream& os) const
	{ return os << m_src; }

	encapsulated_matrix_type &m_src;
	const encapsulated_matrix_type &m_const;
};


// partielle Spezialisierung fuer CPUAlgebra
template<>
inline ScalarMatrixAdapter<CPUAlgebra, CPUAlgebra>::value_type&
ScalarMatrixAdapter<CPUAlgebra, CPUAlgebra>::operator () (size_t r, size_t c)
{ return m_src(r,c); };

#endif