/*
 *  Vector.h
 *  flexamg
 *
 *  Created by Martin Rupp on 04.11.09.
 *  Copyright 2009 . All rights reserved.
 *
 */

///////////////////////////////////////////////////////////////////
//							Vector
///////////////////////////////////////////////////////////////////

#define SPECIALIZE_EXPRESSION_TEMPLATES
//! "big" Vector class for use with the big SparseMatrix
template <typename templ_entry_type>
class Vector : public TE_VEC<Vector<templ_entry_type> >
	//public XD< Vector<templ_entry_type> >
{
	// functions
public:
	typedef templ_entry_type entry_type;
	//typedef subvector<entry_type> subvector_type;
	typedef Vector<templ_entry_type> vector_type;

	//! constructor
	Vector(const char *_name = "");		
	
	//! constructor with length
	Vector(int _length, const char *_name = "");		
	
	//! destructor
	~Vector();
	
private: // forbidden functions
	Vector(Vector&); // disallow copy operator
	
public:
	//! create a vector with specific length
	void create(int _length);
	
	//! access element i of the vector
	inline entry_type &operator [] (int i);
	inline const entry_type &operator [] (int i) const;
	
	

	
	//! returns v.T w, that is the dotprod of this vector and w
	double dotprod(const Vector &w) const;	
		
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
	void set(double d) { operator = (d); }	
	
	
	//! add subvector
/*	void add(const subvector<entry_type> &subvec);
	//! set subvector
	void set(const subvector<entry_type> &subvec);
	//! get subvector
	void get(subvector<entry_type> &subvec) const;*/
	
	void add(const entry_type &d, int i);
	void set(const entry_type &d, int i);
	void get(entry_type &d, int i) const;
	
	
	// for Function Expression, sh. TemplateExpression.h // remove this
	template<class Function> inline void operator = (Function &ex);
	
	//! assign other vector v
	inline void operator = (const Vector &v);
	
	//! assign this vector to another vector v
	inline void applyto(Vector &v) const;
			
	//! template expression assignment
	template<typename Type> inline void operator = (const Type &t);		
	//! template expression +=
	template<typename Type> inline void operator += (const Type &t);	
	//! template expression -=
	template<typename Type> inline void operator -= (const Type &t);
	
	//! return sqrt(sum values[i]^2) (euclidian norm)
	inline double norm();		

	//! printofile: posx posy value
	void printtofile(const char *filename);
			
	int getLength() const { return length; }
	
	void addTo(entry_type &dest, int i) const
	{
		dest += values[i];
	}
	
	void substractFrom(entry_type &dest, int i) const
	{
		dest -= values[i];
	}
	
	void assign(entry_type &dest, int i) const
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
	
	// data
public:
	int level;				///< multigrid level of this vecotr
	const char *name;		///< name 
	
	// TODO: for parallelization: auto-detect non-matching distributions
	enum vector_mode		
	{
		DIST_ADDITIVE, 
		DIST_CONSISTENT
	};
	
	/*vector_mode getDistMode()
	{
		return dist_mode;
	}*/
	
	void assureConsistent();
	
	
private:
	int length;				///< length of the vector (vector is from 0..length-1)
	entry_type *values;		///< array where the values are stored, size length

	//mutable vector_mode dist_mode;

	
};


