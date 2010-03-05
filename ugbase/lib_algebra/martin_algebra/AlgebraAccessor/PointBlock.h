/*
 *  PointBlock.h
 *  flexamg
 *
 *  Created by Martin Rupp on 03.03.10.
 *  Copyright 2010 G-CSC, University of Frankfurt. All rights reserved.
 *
 */

#pragma once
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
		d = 0.0;
		for(int i=0; i<I.size(); i++)
		{
			for(int j=0; j<J.size(); j++)
			{
				if(A(i,j) == 0.0) continue;
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
	
//private:
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
		
		ASSERT1(indexInfo.num_comp() == blocksize);
		
		V.create(indexInfo.num_index());
		return true;
	}

	virtual VectorAccessorBase<T> *newClone() const
	{
		VectorAccessor_PointBlock<T, blocksize> *v = new VectorAccessor_PointBlock<T, blocksize>;
		v->V.create(V.getLength());
		return v;
	}
	
	// s d i a
	virtual bool add(DenseVector &v, vector<multi_index_type> &I)
	{
		blockVector<double, fixedStorage, blocksize> d;
		d = 0.0;

		for(int i=0; i<I.size(); i++)
		{
			if(v(i) == 0.0) continue;
			multi_index_type &ind_i = I[i];
			d(ind_i[3]) = v(i);
			V.add(d, ind_i[2]);
			d(ind_i[3]) = 0.0;
		}
		return true;
	}	
	
	virtual void print()
	{
		V.print();
	}
	
	virtual bool apply(Operation_type operation, const AA_AlphaVec_Expression<T> &ex )
	{
		const VectorAccessor_PointBlock<T, blocksize> *vec = dynamic_cast<const VectorAccessor_PointBlock<T, blocksize>*>(ex.r);
		assert(vec != NULL);
		
		V.apply(operation, ex.l*vec->V);
		return true;
	}
	virtual bool apply(Operation_type operation, const AA_AlphaMatVec_Add_Expression<AA_AlphaVec_Expression<T>, AA_AlphaVec_Expression<T> > &ex)
	
	{
		const VectorAccessor_PointBlock<T, blocksize> *vec1 = dynamic_cast<const VectorAccessor_PointBlock<T, blocksize>*>(ex.l.r);
		const VectorAccessor_PointBlock<T, blocksize> *vec2 = dynamic_cast<const VectorAccessor_PointBlock<T, blocksize>*>(ex.r.r);
		assert(vec1 != NULL); assert(vec2 != NULL);
		
		V.apply(operation, ex.l.l*vec1->V + ex.r.l*vec2->V);
		return true;
	}
	virtual bool apply(Operation_type operation, const AA_AlphaMatVec_Add_Expression<AA_AlphaVec_Expression<T>, AA_MatVec_Expression<T> > &ex)
	{
		const VectorAccessor_PointBlock<T, blocksize> *vec1 = dynamic_cast<const VectorAccessor_PointBlock<T, blocksize>*>(ex.l.r);
		const MatrixAccessor_PointBlock<T, blocksize> *mat2 = dynamic_cast<const MatrixAccessor_PointBlock<T, blocksize>*>(ex.r.l);
		const VectorAccessor_PointBlock<T, blocksize> *vec2 = dynamic_cast<const VectorAccessor_PointBlock<T, blocksize>*>(ex.r.r);
		assert(vec1 != NULL); assert(mat2 != NULL); assert(vec2 != NULL); 
		
		V.apply(operation, ex.l.l*vec1->V + mat2->M*vec2->V);
		return true;
	}
	
	virtual bool apply(Operation_type operation, const AA_MatVec_Expression<T> &ex)
	{
		const MatrixAccessor_PointBlock<T, blocksize> *mat = dynamic_cast<const MatrixAccessor_PointBlock<T, blocksize>*>(ex.l);
		const VectorAccessor_PointBlock<T, blocksize> *vec = dynamic_cast<const VectorAccessor_PointBlock<T, blocksize>*>(ex.r);
		assert(mat != NULL); assert(vec != NULL);
		
		V.apply(operation, mat->M*vec->V);
		return true;
	}
	
	virtual double scal_prod(const VectorAccessorBase<T> *v) const
	{
		const VectorAccessor_PointBlock<T, blocksize> *vec = dynamic_cast<const VectorAccessor_PointBlock<T, blocksize>*>(v);
		assert(vec != 0);
		double d=0;
		d += V.T() * vec->V;
		return d;			
	}
	virtual double norm2() const
	{
		return ::norm2(V);
	}
	
private:
	Vector<blockVector<double, fixedStorage, blocksize> > V;
};

template<typename T, int blocksize>
class AlgebraAccessor_PointBlock : public templateAlgebraAccessorBase<T, MatrixAccessor_PointBlock<T,blocksize>, VectorAccessor_PointBlock<T,blocksize> >
{ };