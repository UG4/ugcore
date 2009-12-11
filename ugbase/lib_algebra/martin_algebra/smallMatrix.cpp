#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
using namespace std;
#include <math.h>

#include "misc.h"
#include "smallMatrix.h"

// define this here so it is available in gdb via "call mat.p()"
//////////////////////////////////////////////////////////////

template<> void fixedMatrix<1>::p() { cout << *this << endl; }
template<> void fixedMatrix<2>::p() { cout << *this << endl; }
template<> void fixedMatrix<3>::p() { cout << *this << endl; }

//////////////////////////////////////////////////////////////

template<> void fixedVector<1>::p() { cout << *this << endl; }
template<> void fixedVector<2>::p() { cout << *this << endl; }
template<> void fixedVector<3>::p() { cout << *this << endl; }

//////////////////////////////////////////////////////////////

#if 0



class Vector;

class smallmatrix
{
	friend class Vector;
	
public:
	// Konstruktor
	smallmatrix()
	{
		mat = NULL;
		z = s = 0;
	}
	
	void create(const int zeilen, const int spalten=-1)
	{
		if(mat) delete[] mat;
		z = zeilen;
		if(spalten == -1)
			s = z;
		else
			s = spalten;
		mat = new double[z*s];
		memset(mat, 0, sizeof(double)*z*s);
	}
	
	smallmatrix(const int zeilen, const int spalten=-1)
	{
		mat = NULL;
		create(zeilen, spalten);
	}
	
	// Copykonstruktor
	smallmatrix(const smallmatrix &m)
	{
		z = m.z;
		s = m.s;
		if(m.mat)
		{
			mat = new double[z*s];
			memcpy(mat, m.mat, sizeof(double)*z*s);
		}
		else mat = NULL;
	}
	
	
	// Dekonstruktor
	~smallmatrix()
	{
		delete[] mat;
	}	
	
	// Operator zum Zugreifen auf Elemente und ver‰ndern
	inline double &operator [] (int zeile, int spalte)
	{	
		return mat[zeile*s+spalte];
	}
	
	// const Operator zum Zugreifen auf Elemente
	inline double operator [] (int zeile, int spalte) const
	{
		return mat[zeile*s+spalte];
	}
	
	// const zum Zugreifen auf Elemente
	inline a(int zeile, int spalte) const
	{
		return mat[zeile*s+spalte];
	}
	

	// getZeilen
	inline int getZeilen() const
	{
		return z;
	}
	
	// getSpalten
	inline int getSpalten() const
	{
		return s;
	}
	
	// print - Ausgabe der smallmatrix
	void print(const char *str=NULL) const
	{
		if(str)
			cout << str << " =" << endl;
		for(int i=0; i<z; i++)
		{
			cout << "|\t";
			for(int j=0; j<s; j++)
				cout << mat[i*s+j] << '\t';
			cout << '|' << endl;			
		}
		if(str)
			cout << endl;
	}
private:
	double *mat; // Zeiger auf Array mit smallmatrixelementen
	int s, z; // s = Spalten, z = Zeilen
};

class cPermutation
{
public:
	cPermutation()
	{
		length = 0;
		vals = NULL;
	}
	cPermutation(int _length)
	{
		vals = NULL;
		create(length);
	}
	
	~cPermutation()	
	{
		if(vals) delete[] vals;
	}
	
	void create(length)
	{
		if(vals) delete[] vals;
		vals = new [length];
		for(int i=0; i<length; i++) vals[i] = i;
	}
	
	void transpose(int i, int j)
	{
		int x = vals[i];
		vals[i] = vals[j];
		vals[j] = x;
	}
	
	int operator [] (int i)
	{
		return vals[i];
	}	
}


class LU_Zerlegung
{
public:
	LU_Zerlegung(const matrix &A)
	{
		create(A);
	}
	
	~LU_Zerlegung() {}
	
	void create(const matrix &A);
	void solve(Vector &x, const Vector &b);
private:
	smallmatrix LU;
	//cPermutation zeilenPermutation;
}

void LU_Zerlegung::create(const matrix &M)
{
	LU.create(M.length, M.length);
	
	for(int i=0; i<M.length; i++)
		for(cmatrixrow::iterator it(M[i]); !it.isEnd(); it++)
			LU(i, (*it).iIndex) = (*it).dValue;
	
	for(int s=0; s < A.getSpalten(); s++)
	{
		// Suchen des grˆﬂte Elements in den Zeilen
		/*int sgross = s;
		cmatrixrow::iterator it(A, s);
		for(j=s+1; j<A.getSpalten(); j++)
		{
			if(dabs(A(s, j)) > dabs(A(s, sgross)))
				sgross = j;
		}				
		
		// Vertausche die spalte sgross und s
		for(j=0; j<A.getSpalten(); j++)
		{
			x = A(j, s);
			A(j, s) = A(j, sgross);
			A(j, sgross) = x;
			
			zeilenPermutation.transpose(sgross);
		}*/
		
		// Elimination
		for(int z=s+1; z < A.length; z++)
		{
			double x;
			LU(z, s) = x = LU(z, s)/KY(s, s); // x = alpha
			for(int j=s+1; j<LU.getSpalten(); j++)
				LU(z, j) -= x*LU(s, j);
		}
	}
	return A;
}


void LU_Zerlegung::solve(Vector &x, const Vector &b)
{
	int n = LU.getSpalten();
	assert(n == LU.getZeilen());
	
	for(int i=0; i<n; i++)
	{
		double s = b[i];
		for(int k=0; k<i; k++)
			s -= LU(i, k)*x[k];
		x[i] = s;
	}
	
	for(int i=n-1; i>=0; i--)
	{
		double s = x[i];
		for(int k=i+1; k<n; k++)
			s -= LU(i, k)*x[k];
		x[i] = s/LU(i,i);
	}
}

#endif