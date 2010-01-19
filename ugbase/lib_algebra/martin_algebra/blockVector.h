/*
 *  blockVector.h
 *  flexamg
 *
 *  Created by Martin Rupp on 06.01.10.
 *  Copyright 2010 . All rights reserved.
 *
 */


template<typename storage_type, int n_>
class blockVector
{
	//storage
	typedef typename storage_traits<storage_type, double, n_, 0>::array_type array_type;
	typedef blockVector<storage_type, n_> vector_type;
	enum { fixed_n=n_};
public:
	inline void setSize(int n, bool bZero=true)
	{
		values.setSize(n, bZero);
	}
	
	inline int getSize() const
	{
		return values.size();
	}	
	//private:
	array_type values;
	
	///
	
	
public:
	double &getAt(int i)
	{
		return values[i];
	}
	double getAt(int i) const
	{
		return values[i];
	}
	double &operator ()(int i)
	{
		return values[i];
	}
	double operator () (int i) const
	{
		return values[i];
	}
	
	
	blockVector() : values()
	{		
	}
	
	blockVector(int n) : values(n)
	{
	}
	
	blockVector(const vector_type &other)
	{
		values = other.values;
	}
	
	friend ostream &operator << (ostream &out, const vector_type &v)
	{
		out << "( ";
		for(int i=0; i < v.getSize(); i++)
			out << v(i) << " ";			
		out << ") ";
		return out;
	}
	
	double operator = (double d)
	{
		for(int i=0; i<getSize(); i++)
			getAt(i) = d;
		return d;
	}
	
	void operator = (const vector_type &other)
	{
		values = other.values;
		//	memcpy(values.values, other.values, sizeof(double)*getSize());	
	}
	
	vector_type operator + (const vector_type &other ) const
	{
		if(other.getSize() == 0)
			return *this;
		else
		{
			ASSERT(getSize() == other.getSize());
			vector_type erg(getSize());
			for(int i=0; i<getSize(); i++)
				erg.values[i] = values[i] + other.values[i];
			return erg;
		}
	}
	
	void operator += (const vector_type &other )
	{
		if(other.getSize() == 0) return;
		ASSERT1(getSize() == other.getSize());
		for(int i=0; i<getSize(); i++)
			values[i] += other.values[i];		
	}
	
	vector_type operator - (const vector_type &other ) const
	{
		if(other.getSize() == 0)
			return *this;
		else
		{
			ASSERT1(getSize() == other.getSize());
			vector_type erg(getSize());
			for(int i=0; i<getSize(); i++)
				erg.values[i] = values[i] - other.values[i];
			return erg;
		}
	}
	
	void operator -= (const vector_type &other )
	{
		if(other.getSize() == 0) return;
		
		ASSERT1(getSize() == other.getSize());
		for(int i=0; i<getSize(); i++)
			values[i] -= other.values[i];		
	}
	
	double operator * (const vector_type &other ) const
	{
		if(other.getSize() == 0) return 0.0;			
		ASSERT(getSize() == other.getSize());
		double s=0;
		for(int i=0; i<getSize(); i++)
			s += values[i] * other.values[i];
		return s;
	}
	
	vector_type operator * (double alpha) const
	{
		vector_type erg(getSize());
		for(int i=0; i<getSize(); i++)
			erg(i) = getAt(i) * alpha;
		return erg;
	}
	
	template<typename array_type>
	vector_type operator * (const blockDenseMatrix<array_type> &mat )
	{
		return (mat * (*this));
	}
	
	
	template<typename array_type>
	inline void operator /= (const blockDenseMatrix<array_type> &mat);
	
	template<typename array_type>
	vector_type operator / (const blockDenseMatrix<array_type> &mat )
	{
		vector_type erg = *this;
		erg /= mat;
		return erg;
	}
	
	double norm2() const
	{
		double s =0;
		for(int i=0; i<getSize(); i++) s += getAt(i)*getAt(i);
		return s;
	}
	
	void p();
	void print() { p(); }
};

template<typename storage_type, int n>
inline blockVector<storage_type, n> operator * (double alpha, const blockVector<storage_type, n> &v) 
{
	return v * alpha;
}

/*
 template<typename storage_type>
 inline void blockVector<storage_type, 1>::operator /= (const blockDenseMatrix<storage_type, 1, 1> &mat )
 {
 getAt(0) = getAt(0) / mat(0, 0);
 }
 
 template<typename storage_type>
 inline void blockVector<storage_type, 2>::operator /= (const blockDenseMatrix<storage_type, 2> &mat )
 {
 double invD = 1.0/(mat(0, 0)*mat(1, 1) - mat(0, 1)*mat(1, 0));
 ASSERT(invD != 0.0);
 
 double a = invD*(getAt(0) * mat(1,1) - getAt(1) * mat(0,1));
 double b = invD*(getAt(1) * mat(0,0) - getAt(0) * mat(1,0));
 getAt(0) = a;
 getAt(1) = b;
 }
 
 
 template<typename storage_type>
 inline void blockVector<storage_type, rows_>::operator /= (const blockDenseMatrix<storage_type, rows_, rows_> &mat )
 {
 double densemat[n*n];	
 
 for(int r=0; r<n; r++)
 for(int c=0; c<n; c++)
 densemat[c + r*n] = mat(r, c);
 
 __CLPK_integer interchange[n];
 
 __CLPK_integer info = 0;
 __CLPK_integer dim = n;
 dgetrf_(&dim, &dim, densemat, &dim, interchange, &info);
 ASSERT2(info == 0, "info is " << info << ( info > 0 ? ": SparseMatrix singular in U(i,i)" : ": i-th argument had had illegal value"));
 
 char trans ='N';
 __CLPK_integer nrhs = 1;
 dgetrs_(&trans, &dim, &nrhs, densemat, &dim, interchange, values, &dim, &info);	
 }
 */
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



template<typename storage_type, int n>
blockVector<storage_type, n> operator * (double alpha, blockVector<storage_type, n> &vec)
{
	return vec*alpha;
}