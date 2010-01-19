/*
 *  Vector.hpp
 *  flexamg
 *
 *  Created by Martin Rupp on 04.11.09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#include <fstream>

template<typename vec_type>
inline vec_type &Vector<vec_type>::operator [] (int i)
{
	ASSERT2(i >= 0 && i < length, *this << ": tried to access element " << i);	
	return values[i];
}

template<typename vec_type>
inline const vec_type &Vector<vec_type>::operator [] (int i) const
{
	ASSERT2(i >= 0 && i < length, *this << ": tried to access element " << i);
	return values[i];
}


// energynorm2 = x*(A*x)
/*inline double Vector<vec_type>::energynorm2(const SparseMatrix &A) const
{
	double sum=0;
	for(int i=0; i<length; i++)	sum += (A[i] * (*this)) * values[i];
	//FOR_UNROLL_FWD(i, 0, length, UNROLL, sum += A[i] * (*this) * values[i]);
	return sum;
}*/

// dotprod
template<typename vec_type>
inline double Vector<vec_type>::dotprod(const Vector &w) const
{
	ASSERT2(length == w.length,  *this << " has not same length as " << w);
	
	double sum=0;
	for(int i=0; i<length; i++)	sum += values[i] * w[i];
	//FOR_UNROLL_FWD(i, 0, length, UNROLL, sum += values[i] * w[i]);
	return sum;
}

template<typename vec_type>
inline double Vector<vec_type>::operator *(const Vector &w)
{
	return dotprod(w);	
}

// assign double to whole Vector
template<typename vec_type>
inline double Vector<vec_type>::operator = (double d)
{
	for(int i=0; i<length; i++)	
	{
		prefetchReadWrite(values+i+512);
		values[i] = d;
	}
	//FOR_UNROLL_FWD(i, 0, length, UNROLL, values[i] = d);
	return d;
}

template<typename vec_type> 
template<typename Type> inline void Vector<vec_type>::operator = (const Type &t)
{ 
	IF_PRINTLEVEL(5) cout << *this << " = " << t << " (unspecialized) " << endl;
	ASSERT2(t.getLength() == length, *this << " has not same length as " << t);
	for(int i=0; i < length; i++)
	{
		prefetchReadWrite(values+i+512);
		//values[i] = t[i];
		t.copyTo(values[i], i);
	}
	//FOR_UNROLL_FWD(i, 0, length, UNROLL, values[i] = t[i]);
}

// für Function Expression, sh. TemplateExpression.h
template<typename vec_type> 
template<class Function> inline void Vector<vec_type>::operator = (Function &ex)
{ 
	ex.applyto(*this);
}

template<typename vec_type>
inline void Vector<vec_type>::operator = (const Vector &v)
{ 
	v.applyto(*this);
}	
template<typename vec_type>
inline void Vector<vec_type>::applyto(Vector &v) const
{
	ASSERT2(v.length == length, *this << " has not same length as " << v);
	//memcpy(v.values, values, length*sizeof(vec_type));
	for(int i=0; i<length; i++)
		v.values[i] = values[i];
}


/*void operator = (const Expression<SparseMatrix, Multiply_Operator, Vector> ex)
 {
 ASSERT2(ex.getLength() == length, *this << " has not same length as " << ex);
 const matrix &m = ex.l;
 const Vector &r = ex.r;
 //for(int i=0; i < length; i++) values[i] = m[i]*r; 
 FOR_UNROLL_FWD(i, 0, length, UNROLL, values[i] = m[i]*r);
 }	*/

// v += exp
template<typename vec_type> 
template<typename Type> inline void Vector<vec_type>::operator += (const Type &t)
{ 
	IF_PRINTLEVEL(5) cout << *this << " += " << t << " (unspecialized) " << endl;
	ASSERT2(t.getLength() == length, *this << " has not same length as " << t);
	for(int i=0; i < length; i++) 
	{
		prefetchReadWrite(values+i+512);
		//values[i] += t[i]; 
		t.addTo(values[i], i);
	}
	//FOR_UNROLL_FWD(i, 0, length, UNROLL, values[i] += t[i]);
} 

// v -= exp
template<typename vec_type> 
template<typename Type> inline void Vector<vec_type>::operator -= (const Type &t)
{
	IF_PRINTLEVEL(5) cout << *this << " -= " << t << " (unspecialized) " << endl;
	ASSERT2(t.getLength() == length, *this << " has not same length as " << t);
	for(int i=0; i < length; i++) 
	{
		prefetchReadWrite(values+i+512);
		//values[i] -= t[i]; 
		t.substractFrom(values[i], i);
	}			
	//FOR_UNROLL_FWD(i, 0, length, UNROLL, values[i] -= t[i]);
}

/*
#ifdef SPECIALIZE_EXPRESSION_TEMPLATES
inline void Vector::operator = (const Expression<SparseMatrix, Multiply_Operator, Vector> ex)
{
	IF_PRINTLEVEL(5) cout << *this << " = " << ex << endl;
	const SparseMatrix &A = ex.l;
	const Vector &v = ex.r;
	for(int i=0; i<length; i++)
		values[i] = A[i] * v;
}
#endif
*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename vec_type>
Vector<vec_type>::Vector (const char *_name)
{
	if(never_happens) p(); // force creation of this rountines for gdb.
		

	length = 0; values = NULL; name = _name; level = 0;
}	

template<typename vec_type>
Vector<vec_type>::Vector(int _length, const char *_name)
{
	if(never_happens) print(); // force creation of this rountines for gdb.

	length = 0;
	create(_length);
	name = _name;
	level = 0;
}	

template<typename vec_type>
Vector<vec_type>::~Vector()
{
	delete [] values;
}

template<typename vec_type>
void Vector<vec_type>::create(int _length)
{
	ASSERT2(length == 0, *this << " already created");
	length = _length;
	values = new vec_type[length];
}




// printofile: posx posy value
template<typename vec_type>
void Vector<vec_type>::printtofile(const char *filename)
{
	fstream fil(filename, ios::out);	
	
	for(int i=0; i < length; i++)
	{
		postype pos = GetPosForIndex(i);
		fil << pos.x << " " << pos.y << " " << values[i] << endl;
	}
	
}

// print
template<typename vec_type>
void Vector<vec_type>::print(const char * const text) const
{
  
  //cout << endl << "================ " << name;
	if(text) cout << " == " << text;
	cout << " == level: " << level << " length: " << length << " =================" << endl;
	for(int i=0; i<length; i++)
		cout << values[i] << " ";
	cout << endl;
}

template<typename vec_type>
void Vector<vec_type>::printtype() const
{ 
	cout << *this;
}


