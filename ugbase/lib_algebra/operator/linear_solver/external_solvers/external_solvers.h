/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
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
			return SPNULL;
		}


	///	returns if parallel solving is supported
		virtual bool supports_parallel() const {return false;}


	public:

		virtual void double_init(const CPUAlgebra::matrix_type &mat) = 0;

		void mat_preprocess(const matrix_type &A)
		{
			if( A.num_rows() == 0 || A.num_cols() == 0) { m_size = 0; return; }
			STATIC_ASSERT(matrix_type::rows_sorted, Matrix_has_to_have_sorted_rows);

			CPUAlgebra::matrix_type mat;

			#ifdef UG_PARALLEL
				// add slave rows to master rows (in indices which this is possible for)
				matrix_type A_tmp; A_tmp = A;
				MatAddSlaveRowsToMasterRowOverlap0(A_tmp);

				// set zero on slaves
				std::vector<IndexLayout::Element> vIndex;
				CollectUniqueElements(vIndex, A.layouts()->slave());
				SetDirichletRow(A_tmp, vIndex);

				mat.set_storage_type(PST_ADDITIVE);
				mat.set_layouts(CreateLocalAlgebraLayouts());
			#else
				const matrix_type& A_tmp = A;
			#endif

			m_size = GetDoubleSparseFromBlockSparse(mat, A_tmp);
			m_blockSize = mat.num_rows()/A_tmp.num_rows();

			if(m_size != 0)
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
			if(m_size == 0) return true;
			m_c.resize(m_size);
			m_d.resize(m_size);

#ifdef UG_PARALLEL
			m_d.set_storage_type(PST_ADDITIVE);
			m_c.set_storage_type(PST_CONSISTENT);

			//	make defect unique
			SmartPtr<vector_type> spDtmp = d.clone();
			spDtmp->change_storage_type(PST_UNIQUE);

			get_vector(m_d, *spDtmp);
#else
			get_vector(m_d, d);
#endif

			m_c.set(0.0);


			double_apply(m_c, m_d);

			set_vector(m_c, c);

#ifdef 	UG_PARALLEL
			// correction must always be consistent (but is unique by construction)
			c.set_storage_type(PST_UNIQUE);
			c.change_storage_type(PST_CONSISTENT);
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
			return true;
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
