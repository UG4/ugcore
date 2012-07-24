/*
 * vanka.h
 *
 *  Created on: 24.07.2012
 *      Author: Christian Wehner
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__VANKA__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__VANKA__

#include "lib_algebra/operator/interface/operator_iterator.h"

#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

namespace ug{

static const int MAXBLOCKSIZE = 19;

template<typename Matrix_type, typename Vector_type>
bool Vanka_step(const Matrix_type &A, Vector_type &x, const Vector_type &b)
{
	DenseVector< VariableArray1<number> > s;
	DenseVector< VariableArray1<number> > localx;
	DenseMatrix< VariableArray2<number> > mat;
	
	size_t blockind[MAXBLOCKSIZE];

	for(size_t i=0; i < x.size(); i++)
	{
		if (A(i,i)!=0) continue;

		size_t blocksize=0;

		for(typename Matrix_type::const_row_iterator it = A.begin_row(i); it != A.end_row(i) ; ++it){
			blockind[blocksize] = it.index();
			x[it.index()] = 0;
			blocksize++;
		};
		mat.resize(blocksize,blocksize);
		s.resize(blocksize);
		localx.resize(blocksize);
		for (size_t j=0;j<blocksize;j++){
			// fill local block matrix
			for (size_t k=0;k<blocksize;k++){
				mat.subassign(j,k,A(blockind[j],blockind[k]));
			};
			// compute rhs
			typename Vector_type::value_type sj = b[blockind[j]];
			for(typename Matrix_type::const_row_iterator it = A.begin_row(blockind[j]); it != A.end_row(blockind[j]) ; ++it){
				MatMultAdd(sj, 1.0, sj, -1.0, it.value(), x[it.index()]);
			};
			s.subassign(j,sj);
		};
		//UG_LOG("mat = " << mat << "\n\n\n");
		// solve block
		InverseMatMult(localx,1.0,mat,s);
		for (size_t j=0;j<blocksize;j++){
			x[blockind[j]] = localx[j];
		};
	}

	return true;
}

///	Vanka Preconditioner
template <typename TAlgebra>
class Vanka : public IPreconditioner<TAlgebra>
{
	public:
	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	///	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Matrix Operator type
		typedef typename IPreconditioner<TAlgebra>::matrix_operator_type matrix_operator_type;

	///	Base type
		typedef IPreconditioner<TAlgebra> base_type;

	protected:
		using base_type::set_debug;
		using base_type::debug_writer;
		using base_type::write_debug;

	public:
	///	default constructor
		Vanka() {};

	///	Clone
		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			SmartPtr<Vanka<algebra_type> > newInst(new Vanka<algebra_type>());
			newInst->set_debug(debug_writer());
			return newInst;
		}

	///	Destructor
		virtual ~Vanka()
		{};

	protected:
	///	Name of preconditioner
		virtual const char* name() const {return "Vanka";}

	///	Preprocess routine
		virtual bool preprocess(matrix_operator_type& mat)
		{
#ifdef UG_PARALLEL
			if(pcl::GetNumProcesses() > 1)
			{
				//	copy original matrix
				MakeConsistent(mat, m_A);
				//	set zero on slaves
				std::vector<IndexLayout::Element> vIndex;
				CollectUniqueElements(vIndex,  m_A.slave_layout());
				SetDirichletRow(m_A, vIndex);
			}
#endif
			return true;
		}

		virtual bool step(matrix_operator_type& mat, vector_type& c, const vector_type& d)
		{
#ifdef UG_PARALLEL
			if(pcl::GetNumProcesses() > 1)
			{
				//	make defect unique
				// todo: change that copying
				vector_type dhelp;
				dhelp.resize(d.size()); dhelp = d;
				dhelp.change_storage_type(PST_UNIQUE);

				if(!Vanka_step(m_A, c, dhelp)) return false;
				c.set_storage_type(PST_UNIQUE);
				return true;
			}
			else
#endif
			{
				if(!Vanka_step(mat, c, d)) return false;
#ifdef UG_PARALLEL
				c.set_storage_type(PST_UNIQUE);
#endif
				return true;
			}
		}

	///	Postprocess routine
		virtual bool postprocess() {return true;}

	protected:
#ifdef UG_PARALLEL
		matrix_type m_A;
#endif

};


} // end namespace ug

#endif
