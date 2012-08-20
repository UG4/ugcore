/*
 * vector_interface.h
 *
 *  Created on: 22.03.2011
 *      Author: mrupp
 */

#ifndef VECTOR_INTERFACE_H_
#define VECTOR_INTERFACE_H_

template <typename TValueType>
class Vector
{
public:
	typedef TValueType value_type;
	//typedef subvector<value_type> subvector_type;
	typedef Vector<TValueType> vector_type;

public:
	//! constructor
	Vector();

	//! constructor with length
	Vector(size_t _length);

	//! destructor
	~Vector();

	Vector(const vector_type & v);

public:
		//! resize vector
	bool resize(size_t new_length, bool bCopyValues=true);

	//! access element i of the vector
	inline value_type &operator [] (size_t i);
	inline const value_type &operator [] (size_t i) const;

	/** add/set/get a local vector
	 *
	 * The local vector type must provide the following members:
	 * - size()					- length of local vector
	 * - index(size_t i)		- global index for component i
	 * - operator[](size_t i)	- access to value of component i
	 */

	template <typename V> bool add(const V& u);
	template <typename V> bool set(const V& u);
	template <typename V> bool get(V& u) const;


	//! returns v.T w, that is the dotprod of this vector and w
	double dotprod(const Vector &w); //const;


	//double energynorm2(const SparseMatrix &A) const;
	/*double energynorm(const SparseMatrix &A) const
	{
		return sqrt(energynorm2(A));
	}*/

	//! assign double d to whole Vector
	double operator = (double d);
	//! assign double d to whole Vector
	bool set(double d);
	bool set_random(double from, double to);


	//! assign other vector v
	void operator = (const Vector &v);
	void operator += (const Vector &v);
	void operator -= (const Vector &v);

	bool operator *= (const number &a);

	//! return sqrt(sum values[i]^2) (euclidian norm)
	inline double norm() const;

	size_t size();
	void defragment();

public:
	/*size_t begin_index() { return 0;}
	size_t end_index() { return size();}

	value_type *begin() { return values + begin_index(); }
	value_type *end() { return values + end_index(); }*/
};

template<typename TValueType>
bool CloneVector(Vector<TValueType> &dest, const Vector<TValueType> src)
{
	// clone stuff like
	dest.resize(src.size());
}

#endif /* VECTOR_INTERFACE_H_ */
