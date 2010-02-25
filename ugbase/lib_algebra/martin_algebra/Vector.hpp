/*
 *  Vector.hpp
 *  flexamg
 *
 *  Created by Martin Rupp on 04.11.09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#include <fstream>

template<typename entry_type>
inline entry_type &Vector<entry_type>::operator [] (int i)
{
	ASSERT2(i >= 0 && i < length, *this << ": tried to access element " << i);	
	return values[i];
}

template<typename entry_type>
inline const entry_type &Vector<entry_type>::operator [] (int i) const
{
	ASSERT2(i >= 0 && i < length, *this << ": tried to access element " << i);
	return values[i];
}


// energynorm2 = x*(A*x)
/*inline double Vector<entry_type>::energynorm2(const SparseMatrix &A) const
{
	double sum=0;
	for(int i=0; i<length; i++)	sum += (A[i] * (*this)) * values[i];
	//FOR_UNROLL_FWD(i, 0, length, UNROLL, sum += A[i] * (*this) * values[i]);
	return sum;
}*/

// dotprod
template<typename entry_type>
inline double Vector<entry_type>::dotprod(const Vector &w) const
{
	ASSERT2(length == w.length,  *this << " has not same length as " << w);
	
	double sum=0;
	for(int i=0; i<length; i++)	sum += values[i] * w[i];
	//FOR_UNROLL_FWD(i, 0, length, UNROLL, sum += values[i] * w[i]);
	return sum;
}

// assign double to whole Vector
template<typename entry_type>
inline double Vector<entry_type>::operator = (double d)
{
	for(int i=0; i<length; i++)	
	{
		prefetchReadWrite(values+i+512);
		values[i] = d;
	}
	//FOR_UNROLL_FWD(i, 0, length, UNROLL, values[i] = d);
	return d;
}



// für Function Expression, sh. TemplateExpression.h
template<typename entry_type> 
template<class Function> inline void Vector<entry_type>::operator = (Function &ex)
{ 
	ex.applyto(*this);
}

template<typename entry_type>
inline void Vector<entry_type>::operator = (const Vector &v)
{ 
	v.applyto(*this);
}	
template<typename entry_type>
inline void Vector<entry_type>::applyto(Vector &v) const
{
	ASSERT2(v.length == length, *this << " has not same length as " << v);
	//memcpy(v.values, values, length*sizeof(entry_type));
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

template<typename entry_type> 
template<typename Type> inline void Vector<entry_type>::operator = (const Type &t)
{ 
	//IF_PRINTLEVEL(5) cout << *this << " = " << t << " (unspecialized) " << endl;
	ASSERT2(t.getLength() == length, *this << " has not same length as " << t);
	t.preventForbiddenDestination(this);

	for(int i=0; i < length; i++)
	{
		prefetchReadWrite(values+i+512);
		//values[i] = t[i];
		t.assign(values[i], i);
		//t.copyTo(values[i], i);
	}
	//FOR_UNROLL_FWD(i, 0, length, UNROLL, values[i] = t[i]);
}

// v += exp
template<typename entry_type> 
template<typename Type> inline void Vector<entry_type>::operator += (const Type &t)
{ 
	//IF_PRINTLEVEL(5) cout << *this << " += " << t << " (unspecialized) " << endl;
	ASSERT2(t.getLength() == length, *this << " has not same length as " << t);
	t.preventForbiddenDestination(this);
	
	for(int i=0; i < length; i++) 
	{
		prefetchReadWrite(values+i+512);
		//values[i] += t[i]; 
		t.addTo(values[i], i);
	}
	//FOR_UNROLL_FWD(i, 0, length, UNROLL, values[i] += t[i]);
} 

// v -= exp
template<typename entry_type> 
template<typename Type> inline void Vector<entry_type>::operator -= (const Type &t)
{
	IF_PRINTLEVEL(5) cout << *this << " -= " << t << " (unspecialized) " << endl;
	ASSERT2(t.getLength() == length, *this << " has not same length as " << t);
	t.preventForbiddenDestination(this);
	
	for(int i=0; i < length; i++) 
	{
		prefetchReadWrite(values+i+512);
		//values[i] -= t[i]; 
		t.substractFrom(values[i], i);
	}			
	//FOR_UNROLL_FWD(i, 0, length, UNROLL, values[i] -= t[i]);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename entry_type>
Vector<entry_type>::Vector (const char *_name)
{
	if(never_happens) p(); // force creation of this rountines for gdb.
		
	length = 0; values = NULL; name = _name; level = 0;
}	

template<typename entry_type>
Vector<entry_type>::Vector(int _length, const char *_name)
{
	if(never_happens) p(); // force creation of this rountines for gdb.

	length = 0;
	create(_length);
	name = _name;
	level = 0;
}	

template<typename entry_type>
Vector<entry_type>::~Vector()
{
	delete [] values;
}

template<typename entry_type>
void Vector<entry_type>::create(int _length)
{
	ASSERT2(length == 0, *this << " already created");
	length = _length;
	values = new entry_type[length];
}




// printofile: posx posy value
template<typename entry_type>
void Vector<entry_type>::printtofile(const char *filename)
{
	fstream fil(filename, ios::out);	
	
	for(int i=0; i < length; i++)
	{
		postype pos = GetPosForIndex(i);
		fil << pos.x << " " << pos.y << " " << values[i] << endl;
	}
	
}

// print
template<typename entry_type>
void Vector<entry_type>::print(const char * const text) const
{
  
  if(name) cout << endl << "================ " << name;
	if(text) cout << " == " << text;
	cout << " == level: " << level << " length: " << length << " =================" << endl;
	for(int i=0; i<length; i++)
		//cout << values[i] << " ";
		cout << i << ": " << values[i] << endl;
	cout << endl;
}

template<typename entry_type>
void Vector<entry_type>::printtype() const
{ 
	cout << *this;
}

/*template<typename entry_type>
void Vector<entry_type>::add(const subvector<entry_type> &subvec)
{
	for(int i=0; i < subvec.getNr(); i++)
		values[subvec.getIndex(i)] += subvec(i);
}

template<typename entry_type>
void Vector<entry_type>::set(const subvector<entry_type> &subvec)
{
	for(int i=0; i < subvec.getNr(); i++)
		values[subvec.getIndex(i)] = subvec(i);	
}

template<typename entry_type>
void Vector<entry_type>::get(subvector<entry_type> &subvec) const
{
	for(int i=0; i < subvec.getNr(); i++)
		subvec(i) = values[subvec.getIndex(i)];
}*/


template<typename entry_type>
void Vector<entry_type>::add(const entry_type &d, int i)
{
	values[i] += d;
}
template<typename entry_type>
void Vector<entry_type>::set(const entry_type &d, int i)
{
	values[i] = d;
}
template<typename entry_type>
void Vector<entry_type>::get(entry_type &d, int i) const
{
	d = values[i];
}


template<typename entry_type>
double operator *(const TRANSPOSED<Vector<entry_type> > &x, const Vector<entry_type> &y)
{
	return x.T().dotprod(y);
}
