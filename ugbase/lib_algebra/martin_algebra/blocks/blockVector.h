/*
 *  blockVector.h
 *  flexamg
 *
 *  Created by Martin Rupp on 06.01.10.
 *  Copyright 2010 G-CSC, University of Frankfurt. All rights reserved.
 *
 */

#ifndef __H__UG__MARTIN_ALGEBRA__BLOCK_VECTOR__
#define __H__UG__MARTIN_ALGEBRA__BLOCK_VECTOR__

#include "arrayStorage.h"
#include "blockDenseMatrix.h"

namespace ug{
///////////////////////////////////////////////////////////////////////////////////////
//!
//! blockVector
//! template parameters:
//! 1. value_type: i.e. float, double or blockDenseMatrix (recursive)
//! 2. storage_type: fixedStorage, arrayStorage
//! 3. n: with storage_type=fixedStorage, size of the fixed matrix
//! if storage_type=variableStorage, n is ignored
template<typename value_type, typename storage_type, int n_=0>
class blockVector
{
private:
//	- storage -
	typedef typename storage_traits<storage_type, double, n_, 0>::array_type array_type;
	typedef blockVector<value_type, storage_type, n_> vector_type;
	enum { fixed_n=n_};
	array_type values;

	friend class smallInverse<storage_type, n_, n_>;

public:
	inline void setSize(int n, bool bZero=true)
	{
		values.setSize(n, bZero);
	}

	inline int getSize() const
	{
		return values.size();
	}


public:
// access functions
	value_type &getAt(int i)
	{
		return values[i];
	}
	const value_type &getAt(int i) const
	{
		return values[i];
	}
	value_type &operator ()(int i)
	{
		return values[i];
	}
	const value_type &operator () (int i) const
	{
		return values[i];
	}

// creation
	blockVector() : values()
	{
	}

	blockVector(double d) : values(n_)
	{
		for(int i=0; i<n_; i++)
			values[i] = d;
	}

	blockVector(int n) : values(n)
	{
	}

	blockVector(const vector_type &other) : values(other.values)
	{
		//UG_ASSERT(0, "thou shall not use the copy constructor for it is slow");
	}


// algebra functions
	double operator = (double d)
	{
		for(int i=0; i<getSize(); i++)
			getAt(i) = d;
		return d;
	}

	void operator = (const vector_type &other)
	{
		values = other.values; // use assignment operator of array class
	}

// add and substract
	vector_type operator + (const vector_type &other ) const
	{
		if(other.getSize() == 0)
			return *this;
		else
		{
			UG_ASSERT(getSize() == other.getSize(), "");
			vector_type erg(getSize());
			for(int i=0; i<getSize(); i++)
				erg.values[i] = values[i] + other.values[i];
			return erg;
		}
	}

	vector_type operator - (const vector_type &other ) const
	{
		if(other.getSize() == 0)
			return *this;
		else
		{
			UG_ASSERT(getSize() == other.getSize(), "");
			vector_type erg(getSize());
			for(int i=0; i<getSize(); i++)
				erg.values[i] = values[i] - other.values[i];
			return erg;
		}
	}

	void operator += (const vector_type &other )
	{
		if(other.getSize() == 0) return;
		UG_ASSERT(getSize() == other.getSize(), "");
		for(int i=0; i<getSize(); i++)
			values[i] += other.values[i];
	}

	void operator -= (const vector_type &other )
	{
		if(other.getSize() == 0) return;

		UG_ASSERT(getSize() == other.getSize(), "");
		for(int i=0; i<getSize(); i++)
			values[i] -= other.values[i];
	}

	//! dot product
	double operator * (const vector_type &other ) const
	{
		if(other.getSize() == 0) return 0.0;
		UG_ASSERT(getSize() == other.getSize(), "");
		double s=0;
		for(int i=0; i<getSize(); i++)
			s += values[i] * other.values[i];
		return s;
	}

	//! scale vector by alpha
	vector_type operator * (double alpha) const
	{
		vector_type erg(getSize());
		for(int i=0; i<getSize(); i++)
			erg(i) = getAt(i) * alpha;
		return erg;
	}

	vector_type& operator *= (number alpha)
	{
		for(int i=0; i<this->getSize(); i++)
			getAt(i) *= alpha;
		return *this;
	}

	//! multiply with matrix, dont use (obviously wrong)
	//template<typename array_type>
	//vector_type operator * (const blockDenseMatrix<array_type> &mat )
	//{
	//return (mat * (*this));
	//}

	//! calc this = this/mat = mat^{-1} * this

	inline void operator /= (const blockDenseMatrix<value_type, storage_type, n_, n_> &mat);

	//! return mat^{-1} * this

	vector_type operator / (const blockDenseMatrix<value_type, storage_type, n_, n_> &mat) const
	{
		vector_type erg = *this;
		erg /= mat;
		return erg;
	}

	//! return sum_i this[i]^2
	double norm2() const
	{
		double s =0;
		for(int i=0; i<getSize(); i++) s += getAt(i)*getAt(i);
		return s;
	}


	//! this += alpha *vec . use this to prevent temporary variables
	void add_mult(double alpha, const vector_type &vec)
	{
		for(int i=0; i<getSize(); i++)
			//add_mult(getAt(i), alpha, vec(i));
			AddMult(getAt(i), vec(i), alpha);
	}

	//! this -= alpha *vec . use this to prevent temporary variables
	void sub_mult(const double alpha, const vector_type &vec)
	{
		for(int i=0; i<getSize(); i++)
			AddMult(getAt(i), vec(i), alpha);
	}

	//! this = alpha *vec . use this to prevent temporary variables
	void assign_mult(double alpha, const vector_type &vec)
	{
		for(int i=0; i<getSize(); i++)
			AddMult(getAt(i), vec(i), alpha);
	}


// print functions
	void p();
	void print() { p(); }
	friend std::ostream &operator << (std::ostream &out, const vector_type &v)
	{
		out << "( ";
		for(int i=0; i < v.getSize(); i++)
			out << v(i) << " ";
		out << ") ";
		return out;
	}
};

template<typename value_type, typename storage_type, int n>
inline blockVector<value_type, storage_type, n> operator * (double alpha, const blockVector<value_type, storage_type, n> &v)
{
	return v * alpha;
}


 template<typename value_type, typename storage_type, int rows_>
 inline void blockVector<value_type, storage_type, rows_>::operator /= (const blockDenseMatrix<value_type, storage_type, rows_, rows_> &mat )
 {
	 smallInverse<storage_type, rows_, rows_> inv;
	 inv.setAsInverseOf(mat);

	 inv.apply(&values[0], *this);
 }

// TODO: do this with template specialisation
// template<typename value_type, typename storage_type>
// inline void blockVector<value_type, storage_type, 1>::operator /= (const blockDenseMatrix<value_type, storage_type, 1, 1> &mat )

template<>
inline void blockVector<double, fixedStorage, 1>::operator /= (const blockDenseMatrix<double, fixedStorage, 1, 1> &mat )
{
	getAt(0) = getAt(0) / mat(0, 0);
}


template<>
inline void blockVector<double, fixedStorage, 2>::operator /= (const blockDenseMatrix<double, fixedStorage, 2, 2> &mat )
{
	double invD = 1.0/(mat(0, 0)*mat(1, 1) - mat(0, 1)*mat(1, 0));
	UG_ASSERT(invD != 0.0, "");
	double a = invD*(getAt(0) * mat(1,1) - getAt(1) * mat(0,1));
	double b = invD*(getAt(1) * mat(0,0) - getAt(0) * mat(1,0));
	getAt(0) = a;
	getAt(1) = b;
}


/*template<>
 template<typename array_type>
 inline void blockDenseMatrix<array_type, 2>::setAsInverseOf(const blockDenseMatrix<array_type, 2> &mat )
 {
 double invD = 1.0/(mat(0, 0)*mat(1, 1) - mat(0, 1)*mat(1, 0));
 ASSERT(invD != 0.0);
 getAt(0,0) = invD * mat(1,1);
 getAt(0,1) = -invD * mat(0,1);
 getAt(1,0) = -invD * mat(1,0);
 getAt(1,1) = invD * mat(0,0);
 }*/

/*
 a b   1 0
 c d   0 1

 c    bc/a  ,   c/a      0
 c    d         0        1

 c    bc/a  ,     c/a      0
 0    d-bc/a      -c/a     1


 1    b/a  ,     1/a         0
 0    1          -c/(ad-bc)  a/(ad-bc)

 1    0    ,     d/(ad-bc)      -b/(ad-bc)
 0    1          -c/(ad-bc)      a/(ad-bc)
 */



template<typename value_type, typename storage_type, int n>
blockVector<value_type, storage_type, n> operator * (double alpha, blockVector<value_type, storage_type, n> &vec)
{
	return vec*alpha;
}

} // namespace ug

#endif
