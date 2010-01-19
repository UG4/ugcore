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
template <typename templ_entry_type>
class Vector : public XD< Vector<templ_entry_type> >
{
	// functions
public:
	typedef templ_entry_type entry_type;
	typedef subvector<entry_type> subvector_type;

	Vector(const char *_name = "");		
	Vector(int _length, const char *_name = "");		
	~Vector();
	
	Vector(Vector&); // disallow copy operator
	
	void create(int _length);
	
	inline entry_type &operator [] (int i);
	inline const entry_type &operator [] (int i) const;
	
	void print(const char * const text = NULL) const;
	void p(); // gdb
	
	friend ostream &operator<<(ostream &output, const Vector &v)
	{
		output << "Vector " <<  v.name << "[" << v.length << "]";
		return output;
	}
	
	void printtype() const; 
	
	
	double dotprod(const Vector &w) const;	
	inline double operator *(const Vector &w);
	
	//double energynorm2(const SparseMatrix &A) const;
	/*double energynorm(const SparseMatrix &A) const
	{
		return sqrt(energynorm2(A));
	}*/
	
	// assign double to whole Vector
	double operator = (double d);
	void set(double d) { operator = (d); }	
	
	
	void add(const subvector<entry_type> &subvec);
	void set(const subvector<entry_type> &subvec);
	void get(subvector<entry_type> &subvec) const;
	
	
	// für Function Expression, sh. TemplateExpression.h
	template<class Function> inline void operator = (Function &ex);
	
	inline void operator = (const Vector &v);
	inline void applyto(Vector &v) const;
			
	template<typename Type> inline void operator = (const Type &t);		
	template<typename Type> inline void operator += (const Type &t);	
	template<typename Type> inline void operator -= (const Type &t);
	
	inline double norm();		
	void printtofile(const char *filename);
	
	int getLength() const { return length; }
	// data
public:
	int length;
	int level;
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
	entry_type *values;

	//mutable vector_mode dist_mode;
public:
	const char *name;
	
};


