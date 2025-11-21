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

#ifndef SCALAR_SUBMATRIX_ADAPTER_HH_
#define SCALAR_SUBMATRIX_ADAPTER_HH_

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


template<typename AT, typename ST, int R, int C>
class ScalarSubMatrixAdapter{


public:
	using encapsulated_matrix_type = typename AT::matrix_type;
	using value_type = typename ST::matrix_type::value_type;

	using const_row_iterator = typename ST::matrix_type::const_row_iterator;

	ScalarSubMatrixAdapter(encapsulated_matrix_type& mat)
	: m_src(mat), m_const(mat) {}; //, m_subr(subr), m_subc(subc) {};

	// forward
	bool resize_and_clear(size_t newRows, size_t newCols)
	{return m_src.resize_and_clear(newRows, newCols);}

	bool resize_and_keep_values(size_t newRows, size_t newCols)
	{return m_src.resize_and_keep_values(newRows, newCols);}


	value_type &operator () (size_t r, size_t c)
	{
		return BlockRef(m_src(r, c), m_subr, m_subc);
	};
	const value_type &operator () (size_t r, size_t c)  const
	{
		return BlockRef(m_src(r, c), m_subr, m_subc);
	}

	//! returns number of rows
	size_t num_rows() const
	{ return m_src.num_rows(); }

	//! returns the number of cols
	size_t num_cols() const
	{ return m_src.num_cols(); }

	inline bool is_isolated(size_t i) const
	{ return m_src.is_isolated(i); }

	bool scale(double d) {return m_src.scale(d);}
	//SparseMatrix<value_type> &operator *= (double d) { scale(d); return *this; }

	//! returns the total number of connections
	size_t total_num_connections() const
	{ return m_src.total_num_connections(); }

	//! print (underlying) matrix
	void print(const char *text) const {m_src.print(text);}

	//! print (block) row of (underlying) matrix
	void printrow(size_t row) const {m_src.printrow(row);}


	class row_iterator
		{
			typename encapsulated_matrix_type::row_iterator iter;
		public:
			inline void check() const {iter.check(); }
			row_iterator(typename encapsulated_matrix_type::row_iterator _iter)
			: iter(_iter) {}
			~row_iterator() {}
			row_iterator *operator -> () { return iter.operator -> (); }
			bool operator != (const row_iterator &o) const { return *iter != o->iter;  }
			void operator ++ () { ++iter; }
			void operator += (int nr) { iter+=nr; }
			bool operator == (const row_iterator &other) const { return other->iter == *iter;}
			size_t index() const { return iter.index(); }
			value_type &value() { return BlockRef(iter.value(), m_subr, m_subc); }
		};

		class const_row_iterator
		{
			typename encapsulated_matrix_type::const_row_iterator iter;
		    public:
				inline void check() const {iter.check(); }
				const_row_iterator(typename encapsulated_matrix_type::const_row_iterator _iter)
				: iter(_iter) {}
				~const_row_iterator() {}
				const_row_iterator *operator ->() { return iter.operator -> (); }
				bool operator != (const const_row_iterator &o) const { return iter!= o.iter;  }
				void operator ++ () { ++iter; }
				void operator += (int nr) { iter+=nr; }
				bool operator == (const const_row_iterator &other) const { return other.iter == iter;}
				size_t index() const { return iter.index(); }
				const value_type &value() const { return BlockRef(iter.value(), m_subr, m_subc); }
		  };


		row_iterator         begin_row(size_t r)
		{ return row_iterator(m_src.begin_row(r)); }
		row_iterator         end_row(size_t r)
		{ return row_iterator(m_src.end_row(r)); }

		const_row_iterator   begin_row(size_t r) const
		{ return const_row_iterator(m_const.begin_row(r)); }
		const_row_iterator   end_row(size_t r) const
		{return const_row_iterator(m_const.end_row(r)); }


		 bool has_connection(size_t r, size_t c) const
		 { return m_src.has_connection(r,c); }

		row_iterator get_iterator_or_next(size_t r, size_t c)
		{ return row_iterator(m_src.get_iterator_or_next(r,c));}

		const_row_iterator get_connection(size_t r, size_t c, bool &bFound) const
		{ return const_row_iterator(m_const.get_connection(r,c,bFound));}

		const_row_iterator get_connection(size_t r, size_t c) const
		{ return const_row_iterator(m_const.get_connection(r,c)); }

		row_iterator get_connection(size_t r, size_t c)
		{ return row_iterator (m_src.get_connection(r,c));}

		void defragment()
		{m_src.defragment();}



protected:
	encapsulated_matrix_type &m_src;
	const encapsulated_matrix_type &m_const;
	static constexpr int m_subr=R;
	static constexpr int m_subc=C;
};


/*template<>
class ScalarSubMatrixAdapter<CPUAlgebra, CPUAlgebra, 0,0> : public CPUAlgebra{};
*/

#endif
