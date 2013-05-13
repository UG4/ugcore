/**
 * \file matrixrow.h
 *
 * \author Martin Rupp
 *
 * \date 18.01.10
 *
 * Goethe-Center for Scientific Computing 2010.
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
	TValue &value() { check(); return (operator*()).value();   }
	size_t index() const { check(); return (operator*()).index(); }
};
template<typename TValue, typename TIterator>
class ConstAlgebraicConnectionIterator : public TIterator
{
	using TIterator::operator *;
public:
	AlgebraicConnectionIterator(TIterator &it) : TIterator(it) {}
	const TValue &value() { check(); return (operator*()).value();   }
	size_t index() const { check(); return (operator*()).index(); }
};*/


/// \addtogroup cpu_algebra
/// \{


template<typename TMatrix>
class MatrixRow
{
	TMatrix &A;
	size_t r;
public:
	typedef typename TMatrix::row_iterator iterator;
	typedef typename TMatrix::const_row_iterator const_iterator;
	typedef typename TMatrix::value_type value_type;
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
	
	value_type &operator()(size_t c)
	{
		return A(r, c);
	}
	
	value_type &operator()(size_t c) const
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
	typedef typename TMatrix::const_row_iterator const_iterator;
	typedef typename TMatrix::value_type value_type;
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
	
	const value_type &operator()(size_t c) const
	{
		return A(r, c);
	}
	value_type &operator()(size_t c)
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
