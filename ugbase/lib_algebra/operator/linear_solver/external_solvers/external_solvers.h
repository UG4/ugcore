/*
 * \file external_solver.h
 *
 *
 * \date 16.01.2014
 * \author Martin Rupp
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__EXTERNAL_SOLVER_
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__EXTERNAL_SOLVER_

#include "common/common.h"
#include "lib_algebra/operator/interface/matrix_operator_inverse.h"


#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif
#include "common/progress.h"
#include "common/util/ostream_util.h"

#include "lib_algebra/operator/linear_solver/linear_solver.h"
#include "lib_algebra/cpu_algebra_types.h"


namespace ug{

class IExternalSolverImplementation
{
public:
	virtual bool init(const CPUAlgebra::matrix_type &A) = 0;
	virtual bool apply(CPUAlgebra::vector_type &c, const CPUAlgebra::vector_type &d) = 0;
	virtual const char* name() const = 0;
	virtual ~IExternalSolverImplementation() {}
};

template <typename TAlgebra>
class IExternalSolver
		: public IMatrixOperatorInverse<typename TAlgebra::matrix_type,
			  	  	  	  	  	  	  	    typename TAlgebra::vector_type>,
			  	  	  	  	  	  public VectorDebugWritingObject<typename TAlgebra::vector_type>
{
	public:
		virtual const char *double_name() const = 0;

		virtual const char* name() const
		{
			return double_name();
		}

	//	Algebra type
		typedef TAlgebra algebra_type;

	//	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	//	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	public:
	//	Constructor
		IExternalSolver()
		{
			m_size = 0;
			m_blockSize = 0;
		};

	// 	Clone

		SmartPtr<ILinearIterator<vector_type> > clone()
		{
			UG_THROW("");
			return NULL;
		}


	///	returns if parallel solving is supported
		virtual bool supports_parallel() const {return false;}


	public:

		virtual void double_init(const CPUAlgebra::matrix_type &mat) = 0;

		void mat_preprocess(const matrix_type &A)
		{
			STATIC_ASSERT(matrix_type::rows_sorted, Matrix_has_to_have_sorted_rows);

			CPUAlgebra::matrix_type mat;

			#ifdef UG_PARALLEL
				mat.set_storage_type(PST_ADDITIVE);
				mat.set_layouts(CreateLocalAlgebraLayouts());
			#endif

			m_size = GetDoubleSparseFromBlockSparse(mat, A);
			m_blockSize = mat.num_rows()/A.num_rows();

			double_init(mat);
		}

		SmartPtr<MatrixOperator<matrix_type, vector_type> > m_spOperator;
		virtual bool init(SmartPtr<MatrixOperator<matrix_type, vector_type> > Op)
		{
		// 	remember operator
			m_spOperator = Op;

		//	init LU operator
			mat_preprocess(m_spOperator->get_matrix());
		//	we're done
			return true;
		}


	protected:

	//	Preprocess routine
		virtual bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp)
		{
			matrix_type &A = *pOp;
			mat_preprocess(A);

			return true;
		}

		void get_vector(CPUAlgebra::vector_type &v_scalar, const vector_type& v)
		{
			for(size_t i=0, k=0; i<v.size(); i++)
			{
				for(size_t j=0; j<GetSize(v[i]); j++)
					v_scalar[k++] = BlockRef(v[i],j);
			}
		}

		void set_vector(CPUAlgebra::vector_type &v_scalar, vector_type& v)
		{
			for(size_t i=0, k=0; i<v.size(); i++)
			{
				for(size_t j=0; j<GetSize(v[i]); j++)
					BlockRef(v[i],j) = v_scalar[k++];
			}
		}

		virtual bool double_apply(CPUAlgebra::vector_type &c, const CPUAlgebra::vector_type &d) = 0;


	//	Stepping routine
		virtual bool apply(vector_type& c, const vector_type& d)
		{
			m_c.resize(m_size);
			m_d.resize(m_size);

#ifdef UG_PARALLEL
			m_d.set_storage_type(PST_ADDITIVE);
			m_c.set_storage_type(PST_CONSISTENT);
#endif
			get_vector(m_d, d);
			m_c.set(0.0);


			double_apply(m_c, m_d);

			set_vector(m_c, c);

#ifdef 	UG_PARALLEL
		//	Correction is always consistent
		//	todo: We set here correction to consistent, but it is not. Think about how to use ilu in parallel.
			c.set_storage_type(PST_CONSISTENT);
#endif
			return true;
		}

		virtual bool apply_return_defect(vector_type& u, vector_type& f)
		{
		//	solve u
			if(!apply(u, f)) return false;

		//	update defect
			if(!m_spOperator->matmul_minus(f, u))
			{
				UG_LOG("ERROR in 'LU::apply_return_defect': "
						"Cannot apply matmul_minus.\n");
				return false;
			}

		//	we're done
			return true;
		}

public:
		using VectorDebugWritingObject<typename TAlgebra::vector_type>::vector_debug_writer;

		int get_dim()
		{
			if(vector_debug_writer().valid() == false) return -1;
			vector_debug_writer()->update_positions();
			return vector_debug_writer()->get_dim();
		}

		template<int Tdim>
		bool get_positions(std::vector<MathVector<Tdim> > &coord)
		{
			UG_COND_THROW(vector_debug_writer().valid() == false, "no debug writer set.");
			int dim = get_dim();
			UG_COND_THROW(dim != Tdim, "wrong dimension");

			return copy_pos<Tdim, Tdim>(coord, get_positions<Tdim>());
		}

		bool get_positions3(std::vector<MathVector<3> > &coord)
		{
			UG_COND_THROW(vector_debug_writer().valid() == false, "no debug writer set.");
			int dim = get_dim();
			switch(dim)
			{
			case 1:
				return copy_pos(coord, get_positions<1>());
			case 2:
				return copy_pos(coord, get_positions<2>());
			case 3:
				return copy_pos(coord, get_positions<3>());
			case -1:
				return false;
			}
		}


		template<int dim>
		const std::vector<MathVector<dim> > &get_positions()
		{
			if(vector_debug_writer().valid())
			{
				vector_debug_writer()->update_positions();
				return vector_debug_writer()->template get_positions<dim>();
			}
			else UG_THROW("no debug_writer!");
		}
		template<int dim1, int dim2>
		bool copy_pos(std::vector<MathVector<dim1> > &dest, const std::vector<MathVector<dim2> > &src)
		{
			UG_COND_THROW(m_size == 0 || m_blockSize == 0, "not initialized");
			UG_COND_THROW(dim1 < dim2, "loss of data");

			dest.resize(m_size);
			for(size_t i=0; i<src.size(); i++)
			{
				for(size_t k=0; k<m_blockSize; k++)
				{
					dest[i*m_blockSize + k]=0;
					for(size_t j=0; j<dim2; j++)
						dest[i*m_blockSize + k][j] = src[i][j];
				}
			}
		}


	protected:
	//	Postprocess routine
		virtual bool postprocess() {return true;}

	CPUAlgebra::vector_type m_c, m_d;
	size_t m_size;
	size_t m_blockSize;
};

} // namespace ug

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__EXTERNAL_SOLVER_ */
