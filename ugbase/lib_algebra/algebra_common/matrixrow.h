/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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

#ifndef __H__UG__CPU_ALGEBRA__MATRIXROW__
#define __H__UG__CPU_ALGEBRA__MATRIXROW__

namespace ug{
///////////////////////////////////////////////////////////////////

/*template<typename TValue, typename TIterator>
class AlgebraicConnectionIterator : public TIterator
{
	using TIterator::operator *;
public:
	AlgebraicConnectionIterator(TIterator &it) : TIterator(it) {}
	TValue &value() { check(); return (operator * ()).value();   }
	size_t index() const { check(); return (operator * ()).index(); }
};
template<typename TValue, typename TIterator>
class ConstAlgebraicConnectionIterator : public TIterator
{
	using TIterator::operator *;
public:
	AlgebraicConnectionIterator(TIterator &it) : TIterator(it) {}
	const TValue &value() { check(); return (operator * ()).value();   }
	size_t index() const { check(); return (operator * ()).index(); }
};*/


/// \addtogroup cpu_algebra
/// \{


template<typename TMatrix>
class MatrixRow
{
	TMatrix &A;
	size_t r;
public:
	using iterator = typename TMatrix::row_iterator;
	using const_iterator = typename TMatrix::const_row_iterator;
	using value_type = typename TMatrix::value_type;
	MatrixRow(TMatrix &_A, size_t _r) : A(_A), r(_r)
	{
	}

	iterator begin()
	{
		return A.begin_row(r);
	}

	iterator end()
	{
		return A.end_row(r);
	}

	iterator begin() const
	{
		return A.begin_row(r);
	}

	iterator end() const
	{
		return A.end_row(r);
	}
	
	value_type &operator () (size_t c)
	{
		return A(r, c);
	}
	
	value_type &operator () (size_t c) const
	{
		return A(r, c);
	}

	bool has_connection(size_t c) const
	{
		return A.has_connection(r, c);
	}

	size_t size() const
	{
		return A.num_cols();
	}

	size_t num_connections() const
	{
		return A.num_connections(r);
	}
};

template<typename TMatrix>
class ConstMatrixRow
{
	const TMatrix &A;
	size_t r;
public:
	using const_iterator = typename TMatrix::const_row_iterator;
	using value_type = typename TMatrix::value_type;

	ConstMatrixRow(const TMatrix &_A, size_t _r) : A(_A), r(_r)
	{
	}
	const_iterator begin() const
	{
		return A.begin_row(r);
	}
	const_iterator end() const
	{
		return A.end_row(r);
	}
	
	const value_type &operator () (size_t c) const
	{
		return A(r, c);
	}
	value_type &operator () (size_t c)
	{
		return A(r, c);
	}
	bool has_connection(size_t c) const
	{
		return A.has_connection(r, c);
	}

	size_t num_connections() const
	{
		return A.num_connections(r);
	}
	size_t size() const
	{
		return A.num_cols();
	}
};



// end group cpu_algebra
/// \}

} // namespace ug



#endif
