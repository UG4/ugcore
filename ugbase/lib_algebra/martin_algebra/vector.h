/*
 *  Vector.h
 *  flexamg
 *
 *  Created by Martin Rupp on 04.11.09.
 *  Copyright 2009 G-CSC, University of Frankfurt. All rights reserved.
 *
 */
#ifndef __H__UG__MARTIN_ALGEBRA__VECTOR__
#define __H__UG__MARTIN_ALGEBRA__VECTOR__

#include "sparsematrix.h"

#include "template_expressions.h"

namespace ug{
///////////////////////////////////////////////////////////////////
//							Vector
///////////////////////////////////////////////////////////////////

#define SPECIALIZE_EXPRESSION_TEMPLATES
//!
//! "big" Vector class for use with the big SparseMatrix
//! can = template expressions like x = 0.5*x - y + A*z
//! see TemplateExpressions.h
template <typename templ_entry_type>
class Vector : public TE_VEC<Vector<templ_entry_type> >
	//public XD< Vector<templ_entry_type> >
{
public:
	typedef MultiIndex<1> index_type;
	typedef FlexLocalVector<double> local_vector_type;
	typedef std::vector<MultiIndex<1> > local_index_type;

	// functions
public:
	typedef templ_entry_type entry_type;
	//typedef subvector<entry_type> subvector_type;
	typedef Vector<templ_entry_type> vector_type;

	//! constructor
	Vector(const char *_name = "");

	//! constructor with length
	Vector(size_t _length, const char *_name = "");

	//! destructor
	~Vector();

private: // forbidden functions
	Vector(Vector&); // disallow copy operator

public:
	//! create a vector with specific length
	bool create(size_t _length);
	//! create as a copy of other vector
	bool create(const Vector &v);

	bool destroy();


	//! access element i of the vector
	inline entry_type &operator [] (size_t i);
	inline const entry_type &operator [] (size_t i) const;




	//! returns v.T w, that is the dotprod of this vector and w
	double dotprod(const Vector &w); //const;

	// deprecated, use x.T() * y.
	//inline double operator *(const Vector &w); ///< shortcut for .dotprod(w)

	//double energynorm2(const SparseMatrix &A) const;
	/*double energynorm(const SparseMatrix &A) const
	{
		return sqrt(energynorm2(A));
	}*/

	//! assign double d to whole Vector
	double operator = (double d);
	//! assign double d to whole Vector
	bool set(double d) { operator = (d); return true; }


	//! add subvector
/*	void add(const subvector<entry_type> &subvec);
	//! set subvector
	void set(const subvector<entry_type> &subvec);
	//! get subvector
	void get(subvector<entry_type> &subvec) const;*/

	void add(const entry_type &d, size_t i);
	void set(const entry_type &d, size_t i);
	void get(entry_type &d, size_t i) const;

	/** Add a local vector
	 *
	 * The local vector type must provide the following members:
	 * - size()					- length of local vector
	 * - index(size_t i)		- global index for component i
	 * - operator[](size_t i)	- access to value of component i
	 */
	template <typename V>
	bool add(const V& u);

	bool add(const local_vector_type &u, const local_index_type &ind);
	bool set(const local_vector_type &u, const local_index_type &ind);
	bool get(local_vector_type &u, const local_index_type &ind) const;

	//template<typename T> inline void apply(Operation_type op, const T &t);

	//! assign other vector v
	inline void operator = (const Vector &v);
	inline bool operator *= (const number &a)
	{
		for(size_t i=0; i<size(); i++) values[i] *= a;
		return true;
	}

	//! assign this vector to another vector v
	inline void applyto(Vector &v) const;

	// template expressions for functions
	template<class Function> inline void operator = (Function &ex);

	//! template expression assignment
	template<typename Type> inline void operator = (const Type &t);
	//! template expression +=
	template<typename Type> inline void operator += (const Type &t);
	//! template expression -=
	template<typename Type> inline void operator -= (const Type &t);

	//! return sqrt(sum values[i]^2) (euclidian norm)
	inline double norm() { return ug::norm(*this); }
	inline double two_norm() { return norm(); }

	//! printofile: posx posy value
	void printtofile(const char *filename);

	size_t size() const { return length; }

	void addTo(entry_type &dest, size_t i) const
	{
		dest += values[i];
	}

	void substractFrom(entry_type &dest, size_t i) const
	{
		dest -= values[i];
	}

	void assign(entry_type &dest, size_t i) const
	{
		dest = values[i];
	}

	void preventForbiddenDestination(void *p, bool &bFirst) const
	{
		assert(bFirst == true || p != this);
		bFirst = false;
	}

public: // output functions
	//! print vector to console
	void print(const char * const text = NULL) const;
	void p() {print(); } ///< gdb shortcut for print

	//! ostream << operator
	friend ostream &operator<<(ostream &output, const Vector &v)
	{
		output << "Vector " <<  v.name << "[" << v.length << "]";
		return output;
	}

	void printtype() const;
	size_t finalize() { return true; }

public:
	int firstIndex() { return 0; }
	int lastIndex() { return size(); }

	// data
public:
	size_t level;				///< multigrid level of this vecotr
	const char *name;		///< name

	// TODO: for parallelization: auto-detect non-matching distributions
	/*enum vector_mode
	{
		DIST_ADDITIVE,
		DIST_CONSISTENT
	};

	vector_mode getDistMode()
	{
		return dist_mode;
	}

	void assureConsistent(); */


private:
	size_t length;				///< length of the vector (vector is from 0..length-1)
	entry_type *values;		///< array where the values are stored, size length

	//mutable vector_mode dist_mode;
};

} // namespace ug

#include "vector_impl.h"

#endif
