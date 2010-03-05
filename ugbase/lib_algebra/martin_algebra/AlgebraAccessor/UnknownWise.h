/*
 *  UnknownWise.h
 *  flexamg
 *
 *  Created by Martin Rupp on 03.03.10.
 *  Copyright 2010 G-CSC, University of Frankfurt. All rights reserved.
 *
 */

#pragma once
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
				mat(i,j).create(indexInfo.num_index(), indexInfo.num_index());		
		
		return true;
	}
	

	// s d i a
	virtual bool add(DenseMatrix &A, vector<T> &I, vector<T> &J)
	{
		for(int i=0; i<I.size(); i++)
		{
			for(int j=0; j<J.size(); j++)
			{
				if(A(i,j) == 0.0) continue;
				multi_index_type &ind_i = I[i];
				multi_index_type &ind_j = J[j];
				double &d = A(i,j);
				mat(ind_i[3], ind_j[3]).add(d, ind_i[2], ind_j[2]);
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
				mat(i,j).print();
			}
	}
	
	const SparseMatrix<double> &mat(int i, int j) const
	{
		assert(i >=0 && j >= 0 && i < m_rows && j < m_rows);
		return M[i*m_rows+j];
	}
	
	SparseMatrix<double> &mat(int i, int j)
	{
		assert(i >=0 && j >= 0 && i < m_rows && j < m_rows);
		return M[i*m_rows+j];
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
	//friend class VectorAccessor_UnknownWise<T>;
	//private:
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
		if(V) delete[] V;
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
			V[i].create(indexInfo.num_index());		
		
		return true;
	}
	
	virtual VectorAccessorBase<T> *newClone() const
	{
		VectorAccessor_UnknownWise<T> *v = new VectorAccessor_UnknownWise<T>;
		v->V = new Vector<double> [m_rows];
		v->m_rows = m_rows;
		for(int i=0;i<m_rows;i++)
			v->V[i].create(V[i].getLength());
		return v;
	}
	
	
	// s d i a
	virtual bool add(DenseVector &v, vector<multi_index_type> &I)
	{
		for(int i=0; i<I.size(); i++)
		{
			if(v(i) == 0.0) continue;
			multi_index_type &ind_i = I[i];
			V[ ind_i[3] ].add(v(i), ind_i[2]);
		}		
		return true;
	}	
	
	virtual void print()
	{
		for(int i=0; i<m_rows; i++)
		{
			cout << "(" << i <<  ")" << endl;
			V[i].print();
		}
	}
	
	virtual bool apply(Operation_type operation, const AA_AlphaVec_Expression<T> &ex )
	{
		const VectorAccessor_UnknownWise<T> *vec = dynamic_cast<const VectorAccessor_UnknownWise<T>*>(ex.r);
		assert(vec != NULL); assert(vec->m_rows == m_rows);
		
		for(int i=0; i<m_rows; i++)
			V[i].apply(operation, ex.l*vec->V[i]);
		return true;
	}
	virtual bool apply(Operation_type operation, const AA_AlphaMatVec_Add_Expression<AA_AlphaVec_Expression<T>, AA_AlphaVec_Expression<T> > &ex)
	
	{
		const VectorAccessor_UnknownWise<T> *vec1 = dynamic_cast<const VectorAccessor_UnknownWise<T>*>(ex.l.r);
		const VectorAccessor_UnknownWise<T> *vec2 = dynamic_cast<const VectorAccessor_UnknownWise<T>*>(ex.r.r);
		assert(vec1 != NULL); assert(vec2 != NULL); assert(vec1->m_rows == m_rows); assert(vec2->m_rows == m_rows);
		
		for(int i=0; i<m_rows; i++)
			V[i].apply(operation, ex.l.l*vec1->V[i] + ex.r.l*vec2->V[i]);
		return true;
	}
	virtual bool apply(Operation_type operation, const AA_AlphaMatVec_Add_Expression<AA_AlphaVec_Expression<T>, AA_MatVec_Expression<T> > &ex)
	{
		const VectorAccessor_UnknownWise<T> *vec1 = dynamic_cast<const VectorAccessor_UnknownWise<T>*>(ex.l.r);
		const MatrixAccessor_UnknownWise<T> *mat2 = dynamic_cast<const MatrixAccessor_UnknownWise<T>*>(ex.r.l);
		const VectorAccessor_UnknownWise<T> *vec2 = dynamic_cast<const VectorAccessor_UnknownWise<T>*>(ex.r.r);
		assert(vec1 != NULL); assert(mat2 != NULL); assert(vec2 != NULL); 
		assert(vec1->m_rows == m_rows); assert(mat2->m_rows == m_rows); assert(vec2->m_rows == m_rows);
		
		if(m_rows == 1)
			V[0].apply(operation, ex.l.l*vec1->V[0] + mat2->mat(0,0)*vec2->V[0]);
		else
		{
			apply(operation, ex.l);
			apply(getAdd(operation), ex.r);			
		}
		return true;
	}
	
	virtual bool apply(Operation_type operation, const AA_MatVec_Expression<T> &ex)
	{
		const MatrixAccessor_UnknownWise<T> *mat = dynamic_cast<const MatrixAccessor_UnknownWise<T>*>(ex.l);
		const VectorAccessor_UnknownWise<T> *vec = dynamic_cast<const VectorAccessor_UnknownWise<T>*>(ex.r);
		assert(mat != NULL); assert(vec != NULL);
		assert(vec->m_rows == m_rows); assert(mat->m_rows == m_rows);
		if(m_rows == 1)
			V[0].apply(operation, mat->mat(0,0)*vec->V[0]);
		else if(m_rows == 2)
		{
			V[0].apply(operation, mat->mat(0,0)*vec->V[0] + mat->mat(0,1)*vec->V[1]);
			V[1].apply(operation, mat->mat(1,0)*vec->V[0] + mat->mat(1,1)*vec->V[1]);
		}			
		else if(m_rows == 3)
		{
			V[0].apply(operation, mat->mat(0,0)*vec->V[0] + mat->mat(0,1)*vec->V[1] + mat->mat(0,2)*vec->V[2]);
			V[1].apply(operation, mat->mat(1,0)*vec->V[0] + mat->mat(1,1)*vec->V[1] + mat->mat(1,2)*vec->V[2]);
			V[2].apply(operation, mat->mat(2,0)*vec->V[0] + mat->mat(2,1)*vec->V[1] + mat->mat(2,2)*vec->V[2]);
		}
		else if(m_rows == 4)
		{
			V[0].apply(operation, mat->mat(0,0)*vec->V[0] + mat->mat(0,1)*vec->V[1] + mat->mat(0,2)*vec->V[2] + mat->mat(0,3)*vec->V[3]);
			V[1].apply(operation, mat->mat(1,0)*vec->V[0] + mat->mat(1,1)*vec->V[1] + mat->mat(1,2)*vec->V[2] + mat->mat(1,3)*vec->V[3]);
			V[2].apply(operation, mat->mat(2,0)*vec->V[0] + mat->mat(2,1)*vec->V[1] + mat->mat(2,2)*vec->V[2] + mat->mat(2,3)*vec->V[3]);
			V[3].apply(operation, mat->mat(3,0)*vec->V[0] + mat->mat(3,1)*vec->V[1] + mat->mat(3,2)*vec->V[2] + mat->mat(3,3)*vec->V[3]);
		}
		else
		{
			for(int i=0; i<m_rows; i++)
				for(int j=0; j<m_rows; j++)
					V[i].apply(operation, mat->mat(i,j)*vec->V[j]);
		}
		return true;
	}
	
	virtual double scal_prod(const VectorAccessorBase<T> *v) const
	{
		const VectorAccessor_UnknownWise<T> *vec = dynamic_cast<const VectorAccessor_UnknownWise<T>*>(v);
		assert(vec != 0); assert(vec->m_rows == m_rows);
		double d=0;
		for(int i=0; i<m_rows; i++)
			d += V[i].T() * vec->V[i];
		return d;			
	}
	virtual double norm2() const
	{
		double d=0;
		for(int i=0; i<m_rows; i++)
			d += ::norm2(V[i]);
		return d;
	}
	
	
private:
	Vector<double> *V;
	int m_rows;
};

template<typename T>
class AlgebraAccessor_UnknownWise : public templateAlgebraAccessorBase<T, MatrixAccessor_UnknownWise<T>, VectorAccessor_UnknownWise<T> >
{
};