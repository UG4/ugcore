/*
 *  AlgebraAccessor.h
 *  flexamg
 *
 *  Created by Martin Rupp on 23.02.10.
 *  Copyright 2010 . All rights reserved.
 *
 */

#include "interface.h"

typedef blockDenseMatrix<double, variableStorage> DenseMatrix;
typedef blockVector<double, variableStorage> DenseVector;

template<typename T>
class MatrixAccessorBase
{
public:
	typedef T multi_index_type;
	
	MatrixAccessorBase() {}
	virtual ~MatrixAccessorBase() {}
	virtual bool Create(const IndexInfo &i) = 0;
	
	virtual bool add(DenseMatrix &M, vector<T> &I, vector<T> &J) = 0;
	//virtual bool set(DenseMatrix &M, vector<T> &I, vector<T> &J) { return true;};
	//virtual bool get(DenseMatrix &M, vector<T> &I, vector<T> &J) {return true;};	
	virtual void print() {return;};
};


template<typename T>
class VectorAccessorBase
{
public:
	typedef T multi_index_type;
	
	VectorAccessorBase() {}
	virtual ~VectorAccessorBase() {}
	virtual bool Create(const IndexInfo &i) = 0;
	
	virtual bool add(DenseVector &v, vector<T> &I) = 0;
	//virtual bool set(DenseVector &v, vector<T> &I) {return true;};
	//virtual bool get(DenseVector &v, vector<T> &I) {return true;};
};

template<typename T>
class AlgebraAccessorBase
{
public:
	typedef T multi_index_type;
	
	AlgebraAccessorBase() {}
	virtual ~AlgebraAccessorBase() {}
	
	virtual MatrixAccessorBase<T> *newMatrix() = 0;
	virtual VectorAccessorBase<T> *newVector() = 0;	
};


/////////////////////////////
template<typename T, typename M, typename V>
class templateAlgebraAccessorBase : public AlgebraAccessorBase<T>
{
public:
	typedef T multi_index_type;
	
	templateAlgebraAccessorBase() {}
	virtual ~templateAlgebraAccessorBase() {}
	
	virtual MatrixAccessorBase<T> *newMatrix()
	{
		return new M;
	}
	virtual VectorAccessorBase<T> *newVector()
	{
		return new V;
	}
};

///////////////////////////////////////////////////////////////////////////////

template<typename T>
class MatrixAccessor_UnknownWise : public MatrixAccessorBase<T>
{
public:
	typedef T multi_index_type;
	MatrixAccessor_UnknownWise()
	{
		M = NULL;
	}
	
	virtual ~MatrixAccessor_UnknownWise()
	{
		if(M) delete[] M;
	}
	
	virtual bool Create(const IndexInfo &subdomainInfo)
	{
		ASSERT1(subdomainInfo.num_index() == 1); // subdomains
		const IndexInfo &dofgroupInfo = subdomainInfo.get_index_info(0);
		ASSERT1(dofgroupInfo.num_index() == 1); // dofgroups
		const IndexInfo &indexInfo = dofgroupInfo.get_index_info(0);
		
		m_rows = indexInfo.num_comp();
		
		M = new SparseMatrix<double> [m_rows*m_rows];
		for(int i=0;i<m_rows;i++)
			for(int j=0; j<m_rows;j++)
				M[i*m_rows+j].create(indexInfo.num_index(), indexInfo.num_index());		

		return true;
	}
	
	// s d i a
	virtual bool add(DenseMatrix &A, vector<T> &I, vector<T> &J)
	{
		for(int i=0; i<I.size(); i++)
		{
			for(int j=0; j<J.size(); j++)
			{
				multi_index_type &ind_i = I[i];
				multi_index_type &ind_j = J[j];
				double &d = A(i,j);
				M[ ind_i[3] *m_rows + ind_j[3] ].add(d, ind_i[2], ind_j[2]);
			}
		}
		
		return true;
	}
	
	virtual void print()
	{
		for(int i=0; i<m_rows; i++)
			for(int j=0; j<m_rows; j++)
			{
				cout << "(" << i << ", " << j << ")" << endl;
				M[i*m_rows + j].print();
			}
	}
	
	// calcs  b -= Ax
/*	virtual bool matmul_minus(VectorAccessorBase *x_, VectorAccessorBase *b_)
	{
		VectorAccessor_UnknownWise *x = (VectorAccessor_UnknownWise *) x_;
		VectorAccessor_UnknownWise *b = (VectorAccessor_UnknownWise *) b_;
		for(int i=0; i<m_rows; i++)
			for(int j=0; j<m_rows; j++)
				b.V[i].sub_mul(M[i*m_rows+j], x.V[j]);
	}
*/	
private:
	SparseMatrix<double> *M;
	int m_rows;
};


template<typename T>
class VectorAccessor_UnknownWise : public VectorAccessorBase<T>
{
public:
	typedef T multi_index_type;
	VectorAccessor_UnknownWise()
	{
	}
	
	virtual ~VectorAccessor_UnknownWise()
	{
	}
	
	virtual bool Create(const IndexInfo &subdomainInfo)
	{
		ASSERT1(subdomainInfo.num_index() == 1); // subdomains
		const IndexInfo &dofgroupInfo = subdomainInfo.get_index_info(0);
		ASSERT1(dofgroupInfo.num_index() == 1); // dofgroups
		const IndexInfo &indexInfo = dofgroupInfo.get_index_info(0);
		
		m_rows = indexInfo.num_comp();
		
		V = new Vector<double> [m_rows];
		for(int i=0;i<m_rows;i++)
			for(int j=0; j<m_rows;j++)
				V[i].create(indexInfo.num_index());		
	
		return true;
	}
	
	
	// s d i a
	virtual bool add(DenseVector &v, vector<multi_index_type> &I)
	{
		for(int i=0; i<I.size(); i++)
		{
			multi_index_type &ind_i = I[i];
			V[ ind_i[3] ].add(v(i), ind_i[2]);
		}		
		return true;
	}	
	
private:
	Vector<double> *V;
	int m_rows;
};

template<typename T>
class AlgebraAccessor_UnknownWise : public templateAlgebraAccessorBase<T, MatrixAccessor_UnknownWise<T>, VectorAccessor_UnknownWise<T> >
{
};


///////////////////////////////////////////////////////////////////////////////

template<typename T, int blocksize>
class MatrixAccessor_PointBlock : public MatrixAccessorBase<T>
{
public:
	typedef T multi_index_type;
	MatrixAccessor_PointBlock() : M()
	{
	}

	virtual ~MatrixAccessor_PointBlock()
	{
	}

	virtual bool Create(const IndexInfo &subdomainInfo)
	{
		ASSERT1(subdomainInfo.num_index() == 1); // subdomains
		const IndexInfo &dofgroupInfo = subdomainInfo.get_index_info(0);
		ASSERT1(dofgroupInfo.num_index() == 1); // dofgroups
		const IndexInfo &indexInfo = dofgroupInfo.get_index_info(0);
		
		ASSERT1(indexInfo.num_comp() == blocksize);
		
		M.create(indexInfo.num_index(), indexInfo.num_index());
		return true;
	}

	// s d i a
	virtual bool add(DenseMatrix &A, vector<multi_index_type> &I, vector<multi_index_type> &J)
	{
		blockDenseMatrix<double, fixedStorage, blocksize, blocksize> d;
		d = 0;
		for(int i=0; i<I.size(); i++)
		{
			for(int j=0; j<J.size(); j++)
			{
				multi_index_type &ind_i = I[i];
				multi_index_type &ind_j = J[j];
				d(ind_i[3], ind_j[3]) = A(i,j);
				M.add(d, ind_i[2], ind_j[2]);
				d(ind_i[3], ind_j[3]) = 0.0;
			}
		}
		
		return true;
	}
	
	virtual void print()
	{
		M.print();
	}

	/*// calcs  b -= Ax
	virtual bool matmul_minus(VectorAccessorBase *x_, VectorAccessorBase *b_)
	{
		VectorAccessor_PointBlock *x = (VectorAccessor_PointBlock *) x_;
		VectorAccessor_PointBlock *b = (VectorAccessor_PointBlock *) b_;
		b.V.sub_mul(M, x.V);
	}*/

private:
	SparseMatrix<blockDenseMatrix<double, fixedStorage, blocksize, blocksize> > M;
};


template<typename T, int blocksize>
class VectorAccessor_PointBlock : public VectorAccessorBase<T>
{
public:
	typedef T multi_index_type;
	VectorAccessor_PointBlock()
	{
	}
	
	virtual ~VectorAccessor_PointBlock()
	{
	}
	
	virtual bool Create(const IndexInfo &subdomainInfo)
	{
		ASSERT1(subdomainInfo.num_index() == 1); // subdomains
		const IndexInfo &dofgroupInfo = subdomainInfo.get_index_info(0);
		ASSERT1(dofgroupInfo.num_index() == 1); // dofgroups
		const IndexInfo &indexInfo = dofgroupInfo.get_index_info(0);
		
		ASSERT1(indexInfo.num_comp() == 3);
		
		V.create(indexInfo.num_index());
		return true;
	}
	
	// s d i a
	virtual bool add(DenseVector &v, vector<multi_index_type> &I)
	{
		blockVector<double, fixedStorage, blocksize> d;

		for(int i=0; i<I.size(); i++)
		{
			multi_index_type &ind_i = I[i];
			d(ind_i[3]) = v(i);
			V.add(d, ind_i[2]);
			d(ind_i[3]) = 0.0;
		}		
		return true;
	}	
	
private:
	Vector<blockVector<double, fixedStorage, blocksize> > V;
	int m_rows;
};

template<typename T, int blocksize>
class AlgebraAccessor_PointBlock : public templateAlgebraAccessorBase<T, MatrixAccessor_PointBlock<T,blocksize>, VectorAccessor_PointBlock<T,blocksize> >
{
};