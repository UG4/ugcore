/*
 * parallel_lu.h
 *
 *  Created on: 16.06.2010
 *      Author: mrupp
 */

#ifndef __H__LIB_ALGEBRA__LAPACK_AGGLOMERATING_LU_OPERATOR__
#define __H__LIB_ALGEBRA__LAPACK_AGGLOMERATING_LU_OPERATOR__
#include <iostream>
#include <sstream>

#include "lib_algebra/operator/interface/operator_inverse.h"

#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/collect_matrix.h"
	#include "lib_algebra/parallelization/parallelization.h"
	#include "lib_algebra/parallelization/parallelization_util.h"
#endif

namespace ug{


#ifdef UG_PARALLEL
template<typename T>
void GatherVectorOnOne(HorizontalAlgebraLayouts &agglomerationLayout,
		ParallelVector<T> &collectedVec, const ParallelVector<T> &vec,
		ParallelStorageType type)
{
	bool bRoot = pcl::GetProcRank() == vec.layouts()->proc_comm().get_proc_id(0);
	GatherVectorOnOne(agglomerationLayout.master(), agglomerationLayout.slave(), agglomerationLayout.comm(),
			collectedVec, vec, type, bRoot);
}



template<typename T>
void BroadcastVectorFromOne(HorizontalAlgebraLayouts &agglomerationLayout,
		ParallelVector<T> &vec, const ParallelVector<T> &collectedVec,
		ParallelStorageType type)
{
	bool bRoot = pcl::GetProcRank() == vec.layouts()->proc_comm().get_proc_id(0);
	BroadcastVectorFromOne(agglomerationLayout.master(), agglomerationLayout.slave(),
			agglomerationLayout.comm(), vec, collectedVec, type, bRoot);
}

template <typename TAlgebra>
class AgglomeratingBase
{
	public:
	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	// 	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	public:
	// 	Destructor
		virtual ~AgglomeratingBase() {};

		bool i_am_root()
		{
			return m_bRoot;
			//m_bRoot = pcl::GetProcRank() == agglomerationLayout.proc_comm().get_proc_id(0);;
		}

		bool empty()
		{
			return m_bEmpty;
		}

		bool init_mat(const matrix_type &A)
		{
			m_bEmpty = A.layouts()->proc_comm().empty();
			if(m_bEmpty) return true;
			m_bRoot = pcl::GetProcRank() == A.layouts()->proc_comm().get_proc_id(0);
			PROFILE_FUNC();

			m_spCollectedOp = new MatrixOperator<matrix_type, vector_type>();
			matrix_type &collectedA = m_spCollectedOp->get_matrix();

			CollectMatrixOnOneProc(A, collectedA, agglomerationLayout.master(), agglomerationLayout.slave());
			agglomerationLayout.comm() = A.layouts()->comm();
			agglomerationLayout.proc_comm() = A.layouts()->proc_comm();

			m_spLocalAlgebraLayouts = CreateLocalAlgebraLayouts();
			collectedA.set_layouts(m_spLocalAlgebraLayouts);

			return true;
		}

		SmartPtr<MatrixOperator<matrix_type, vector_type> > get_collected_linear_op()
		{
			return m_spCollectedOp;
		}

		void init_collected_vec(vector_type &collectedX)
		{
			if(i_am_root())
			{
				matrix_type &collectedA = m_spCollectedOp->get_matrix();
				collectedX.resize(collectedA.num_rows());
				collectedX.set_layouts(m_spLocalAlgebraLayouts);
			}
		}


	//	set operator L, that will be inverted
		bool init(SmartPtr<MatrixOperator<matrix_type, vector_type> > Op)
		{
			PROFILE_FUNC();
		//	get matrix of Operator
			m_pMatrix = &Op->get_matrix();

		//	check that matrix exist
			if(m_pMatrix == NULL)
				{UG_LOG("ERROR in LUOperator::init: No Matrix given,\n"); return false;}

		//	init LU operator
			if(!init_mat(*m_pMatrix))
				{UG_LOG("ERROR in LUOperator::init: Cannot init LU Decomposition.\n"); return false;}

		//	we're done
			return true;
		}

		void gather_vector_on_one(vector_type &collectedB, const vector_type &b, ParallelStorageType type)
		{
			GatherVectorOnOne(agglomerationLayout, collectedB, b, PST_ADDITIVE);
		}

		void broadcast_vector_from_one(vector_type &x, const vector_type &collectedX, ParallelStorageType type)
		{
			BroadcastVectorFromOne(agglomerationLayout, x, collectedX, PST_CONSISTENT);
		}


	protected:
		// Operator to invert

		// matrix to invert
		matrix_type* m_pMatrix;

		HorizontalAlgebraLayouts agglomerationLayout;
		bool m_bRoot;


		SmartPtr<MatrixOperator<matrix_type, vector_type> > m_spCollectedOp;
		SmartPtr<AlgebraLayouts> m_spLocalAlgebraLayouts;
		bool m_bEmpty;

};



template <typename TAlgebra>
class AgglomeratingSolver : public IMatrixOperatorInverse<	typename TAlgebra::matrix_type,
														typename TAlgebra::vector_type>
{
	public:
	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	// 	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Base type
		typedef IMatrixOperatorInverse<matrix_type,vector_type> base_type;

	protected:
		using base_type::convergence_check;


	public:
		AgglomeratingSolver(SmartPtr<ILinearOperatorInverse<vector_type, vector_type> > linOpInverse)
		{
			UG_COND_THROW(linOpInverse.valid()==false, "linOpInverse has to be != NULL");
			m_pLinOpInverse = linOpInverse;
		};

	// 	Destructor
		virtual ~AgglomeratingSolver() {};


		virtual const char* name() const
		{
			return "AgglomeratingSolver";
		}


	//	set operator L, that will be inverted
		virtual bool init(SmartPtr<MatrixOperator<matrix_type, vector_type> > Op)
		{
			m_pMatrix = &Op->get_matrix();
			if(agglomeration.init(Op)==false) return false;

			bool bSuccess=true;
			if(agglomeration.i_am_root())
			{
				agglomeration.init_collected_vec(collectedX);
				agglomeration.init_collected_vec(collectedB);
				bSuccess = m_pLinOpInverse->init(agglomeration.get_collected_linear_op());
				UG_DLOG(LIB_ALG_LINEAR_SOLVER, 1, "Agglomerated on proc 0. Size is " << collectedX.size() << "\n");
			}
			if(pcl::AllProcsTrue(bSuccess, Op->get_matrix().layouts()->proc_comm()) == false) return false;
			return true;
		}

		virtual bool supports_parallel() const { return true; }
	// 	Compute u = L^{-1} * f
		virtual bool apply(vector_type& x, const vector_type& b)
		{
			if(agglomeration.empty()) return true;

			agglomeration.gather_vector_on_one(collectedB, b, PST_ADDITIVE);

			collectedX.set(0.0);
			if(agglomeration.i_am_root())
				m_pLinOpInverse->apply(collectedX, collectedB);

			agglomeration.broadcast_vector_from_one(x, collectedX, PST_CONSISTENT);
		//	we're done
			return true;
		}

	// 	Compute u = L^{-1} * f AND return defect f := f - L*u
		virtual bool apply_return_defect(vector_type& u, vector_type& f)
		{
			PROFILE_FUNC();
		//	solve u
			if(!apply(u, f)) return false;
		//	update defect
			if(!m_pMatrix->matmul_minus(f, u))
			{
				UG_LOG("ERROR in 'LUSolver::apply_return_defect': Cannot apply matmul_minus.\n");
				return false;
			}

		//	we're done
			return true;
		}

		virtual std::string config_string() const
		{
			std::stringstream ss;
			ss << "AgglomeratingSolver of "<< ConfigShift(m_pLinOpInverse->config_string());
			return ss.str();
		}


	protected:
		// Operator to invert
		// matrix to invert
		matrix_type* m_pMatrix;
		vector_type collectedB, collectedX;

		SmartPtr<ILinearOperatorInverse<vector_type, vector_type> > m_pLinOpInverse;
		AgglomeratingBase<TAlgebra> agglomeration;
};

///	Jacobi Preconditioner
template <typename TAlgebra>
class AgglomeratingPreconditioner : public IPreconditioner<TAlgebra>
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

	public:
	///	default constructor
		AgglomeratingPreconditioner(SmartPtr<IPreconditioner<TAlgebra> > prec)
		{
			m_prec = prec;
		}

	///	returns if parallel solving is supported
		virtual bool supports_parallel() const {return true;}

	///	Clone
		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			return new AgglomeratingPreconditioner<algebra_type>(m_prec->clone());
		}

	///	Destructor
		~AgglomeratingPreconditioner()
		{};

	protected:
	///	Name of preconditioner
		virtual const char* name() const {return "AgglomeratingPreconditioner";}

	///	Preprocess routine
		virtual bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > Op)
		{
			if(agglomeration.init(Op)==false) return false;

			bool bSuccess=true;
			if(agglomeration.i_am_root())
			{
				agglomeration.init_collected_vec(collectedC);
				agglomeration.init_collected_vec(collectedD);
				bSuccess = m_prec->init(agglomeration.get_collected_linear_op());
				UG_DLOG(LIB_ALG_LINEAR_SOLVER, 1, "Agglomerated on proc 0. Size is " << collectedC.size() << "\n");
			}
			if(pcl::AllProcsTrue(bSuccess, Op->get_matrix().layouts()->proc_comm()) == false) return false;
			return true;
		}

		virtual bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d)
		{
			if(agglomeration.empty()) return true;

			agglomeration.gather_vector_on_one(collectedD, d, PST_ADDITIVE);

			collectedC.set(0.0);
			if(agglomeration.i_am_root())
				m_prec->step(collectedC, collectedD);

			agglomeration.broadcast_vector_from_one(c, collectedC, PST_CONSISTENT);
		//	we're done
			return true;
		}

	///	Postprocess routine
		virtual bool postprocess() {return true;}


		vector_type collectedC, collectedD;

		SmartPtr<IPreconditioner<TAlgebra> > m_prec;
		AgglomeratingBase<TAlgebra> agglomeration;
};

#else


template <typename TAlgebra>
class AgglomeratingSolver : public IMatrixOperatorInverse<	typename TAlgebra::matrix_type,
														typename TAlgebra::vector_type>
{
	public:
	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	// 	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Base type
		typedef IMatrixOperatorInverse<matrix_type,vector_type> base_type;

	protected:
		using base_type::convergence_check;


	public:
		AgglomeratingSolver(SmartPtr<ILinearOperatorInverse<vector_type, vector_type> > linOpInverse)
		{
			UG_COND_THROW(linOpInverse.valid()==false, "linOpInverse has to be != NULL");
			m_pLinOpInverse = linOpInverse;
		};

	// 	Destructor
		virtual ~AgglomeratingSolver() {};


		virtual const char* name() const
		{
			return "AgglomeratingSolver";
		}


	//	set operator L, that will be inverted
		virtual bool init(SmartPtr<MatrixOperator<matrix_type, vector_type> > Op)
		{
			return m_pLinOpInverse->init(Op);
			return true;
		}

		virtual bool supports_parallel() const { return true; }

	// 	Compute u = L^{-1} * f
		virtual bool apply(vector_type& x, const vector_type& b)
		{
			m_pLinOpInverse->apply(x, b);
			return true;
		}

	// 	Compute u = L^{-1} * f AND return defect f := f - L*u
		virtual bool apply_return_defect(vector_type& u, vector_type& f)
		{
			PROFILE_FUNC();
		//	solve u
			if(!apply(u, f)) return false;
		//	update defect
			if(!m_pMatrix->matmul_minus(f, u))
			{
				UG_LOG("ERROR in 'LUSolver::apply_return_defect': Cannot apply matmul_minus.\n");
				return false;
			}

		//	we're done
			return true;
		}

		virtual std::string config_string() const
		{
			std::stringstream ss;
			ss << "AgglomeratingSolver of "<< ConfigShift(m_pLinOpInverse->config_string());
			return ss.str();
		}


	protected:
		// Operator to invert
		// matrix to invert
		matrix_type* m_pMatrix;
		SmartPtr<ILinearOperatorInverse<vector_type, vector_type> > m_pLinOpInverse;
};
#endif

} // end namespace ug

#endif /* __H__LIB_ALGEBRA__LAPACK_LU_OPERATOR__ */
