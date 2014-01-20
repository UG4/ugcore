/*
 * agglomerating_solver.h
 *
 *  Created on: 16.06.2010
 *      Author: mrupp
 */

#ifndef __H__LIB_ALGEBRA__LAPACK_AGGLOMERATING_LU_OPERATOR__
#define __H__LIB_ALGEBRA__LAPACK_AGGLOMERATING_LU_OPERATOR__
#include <iostream>
#include <sstream>

#include "lib_algebra/operator/interface/operator_inverse.h"
#include "lib_algebra/operator/interface/matrix_operator_inverse.h"
#include "lib_algebra/operator/interface/preconditioner.h"

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
#endif

template <typename TBase, typename TAlgebra>
class AgglomeratingBase : public TBase
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
			try{
#ifdef UG_PARALLEL
			m_bEmpty = A.layouts()->proc_comm().empty();
			if(m_bEmpty) return true;
			m_bRoot = pcl::GetProcRank() == A.layouts()->proc_comm().get_proc_id(0);
			PROFILE_FUNC();

			m_spCollectedOp = make_sp(new MatrixOperator<matrix_type, vector_type>());
			matrix_type &collectedA = m_spCollectedOp->get_matrix();

			CollectMatrixOnOneProc(A, collectedA, agglomerationLayout.master(), agglomerationLayout.slave());
			agglomerationLayout.comm() = A.layouts()->comm();
			agglomerationLayout.proc_comm() = A.layouts()->proc_comm();

			m_spLocalAlgebraLayouts = CreateLocalAlgebraLayouts();
			collectedA.set_layouts(m_spLocalAlgebraLayouts);
#else
			m_bEmpty = false;
			m_bRoot = true;
#endif
			}UG_CATCH_THROW("AgglomeratingBase::" << __FUNCTION__ << " failed")
			return true;
		}

#ifdef UG_PARALLEL
		void init_collected_vec(vector_type &collectedX)
		{
			if(i_am_root())
			{
				matrix_type &collectedA = m_spCollectedOp->get_matrix();
				collectedX.resize(collectedA.num_rows());
				collectedX.set_layouts(m_spLocalAlgebraLayouts);
			}
		}

		void gather_vector_on_one(vector_type &collectedB, const vector_type &b, ParallelStorageType type)
		{
			GatherVectorOnOne(agglomerationLayout, collectedB, b, PST_ADDITIVE);
		}

		void broadcast_vector_from_one(vector_type &x, const vector_type &collectedX, ParallelStorageType type)
		{
			BroadcastVectorFromOne(agglomerationLayout, x, collectedX, PST_CONSISTENT);
		}
#endif

		virtual bool init_agglomerated(SmartPtr<MatrixOperator<matrix_type, vector_type> > Op) = 0;
		virtual bool apply_agglomerated(vector_type& x, const vector_type& b) = 0;

///	Preprocess routine
		bool base_init(SmartPtr<MatrixOperator<matrix_type, vector_type> > Op)
		{
			try{
			PROFILE_FUNC();
		//	get matrix of Operator
			m_pMatrix = &Op->get_matrix();

#ifdef UG_PARALLEL
			PROFILE_BEGIN(AGGLOMERATING_Solver_BARRIER);
				m_pMatrix->layouts()->proc_comm().barrier();
			PROFILE_END();
		//	check that matrix exist
			if(m_pMatrix == NULL)
				{UG_LOG("ERROR in LUOperator::init: No Matrix given,\n"); return false;}

		//	init LU operator
			if(!init_mat(*m_pMatrix))
				{UG_LOG("ERROR in LUOperator::init: Cannot init LU Decomposition.\n"); return false;}

			bool bSuccess=true;
			if(i_am_root())
			{
				init_collected_vec(collectedX);
				init_collected_vec(collectedB);
				bSuccess = init_agglomerated(m_spCollectedOp);
				//UG_DLOG(LIB_ALG_LINEAR_SOLVER, 1,
						UG_LOG("Agglomerated on proc 0. Size is " << collectedX.size() << "(was on this proc: "
						<< m_pMatrix->num_rows() << ")\n");
			}
			if(pcl::AllProcsTrue(bSuccess, m_pMatrix->layouts()->proc_comm()) == false) return false;
			return true;
#else
			return init_agglomerated(Op);
#endif
			}UG_CATCH_THROW("AgglomeratingBase::" << __FUNCTION__ << " failed")
		}


		virtual bool base_init(SmartPtr<ILinearOperator<vector_type> > A)
		{
		//	cast operator
			SmartPtr<MatrixOperator<matrix_type,vector_type> > op =
									A.template cast_dynamic<MatrixOperator<matrix_type,vector_type> >();

		//	check if correct types are present
			if(op.invalid())
				UG_THROW("IMatrixOperatorInverse::init:"
						" Passed operator is not matrix-based.");

		//	forward request
			return base_init(op);
		}

		virtual bool apply(vector_type& x, const vector_type& b)
		{
			try{
#if UG_PARALLEL
			if(empty()) return true;

			gather_vector_on_one(collectedB, b, PST_ADDITIVE);

			collectedX.set(0.0);
			if(i_am_root())
				apply_agglomerated(collectedX, collectedB);

			broadcast_vector_from_one(x, collectedX, PST_CONSISTENT);
#else
			apply_agglomerated(x, b);
#endif
			}UG_CATCH_THROW("AgglomeratingBase::" << __FUNCTION__ << " failed")
		//	we're done
			return true;
		}

	// 	Compute u = L^{-1} * f AND return defect f := f - L*u

		virtual bool apply_return_defect(vector_type& u, vector_type& f)
		{
			PROFILE_FUNC();
		//	solve u
			if(!apply(u, f)) return false;

		//	calculate defect
			f.set(0.0);
			if(!m_pMatrix->matmul_minus(f, u))
			{
				UG_LOG("ERROR in 'LUSolver::apply_return_defect': Cannot apply matmul_minus.\n");
				return false;
			}

		//	we're done
			return true;
		}

		virtual bool apply_update_defect(vector_type& u, vector_type& f)
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

		virtual bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d)
		{
			return apply(c, d);
		}



	protected:
		// Operator to invert

		// matrix to invert
		matrix_type* m_pMatrix;
#ifdef UG_PARALLEL
		vector_type collectedB, collectedX;
		HorizontalAlgebraLayouts agglomerationLayout;
		SmartPtr<MatrixOperator<matrix_type, vector_type> > m_spCollectedOp;
		SmartPtr<AlgebraLayouts> m_spLocalAlgebraLayouts;
#endif

		bool m_bRoot;
		bool m_bEmpty;

};



template <typename TAlgebra>
class AgglomeratingSolver : public
	AgglomeratingBase<IMatrixOperatorInverse<	typename TAlgebra::matrix_type, typename TAlgebra::vector_type>, TAlgebra >
{
	public:
	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	// 	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Base type
		typedef AgglomeratingBase<IMatrixOperatorInverse<matrix_type,vector_type>, algebra_type > base_type;


	public:
		AgglomeratingSolver(SmartPtr<ILinearOperatorInverse<vector_type, vector_type> > linOpInverse)
		{
			UG_COND_THROW(linOpInverse.valid()==false, "linOpInverse has to be != NULL");
			m_pLinOpInverse = linOpInverse;
			m_name = std::string("AgglomeratingSolver(") + linOpInverse->name() + ")";
		};

		virtual bool init_agglomerated(SmartPtr<MatrixOperator<matrix_type, vector_type> > Op)
		{
			try{
				return m_pLinOpInverse->init(Op);
			}UG_CATCH_THROW("AgglomeratingSolver::" << __FUNCTION__ << " failed")
		}
		virtual bool apply_agglomerated(vector_type& x, const vector_type& b)
		{
			try{
				return m_pLinOpInverse->apply(x, b);
			}UG_CATCH_THROW("AgglomeratingSolver::" << __FUNCTION__ << " failed")
		}

	// 	Destructor
		virtual ~AgglomeratingSolver() {};


		virtual const char* name() const
		{
			return m_name.c_str();
		}

		virtual bool supports_parallel() const { return true; }

		virtual std::string config_string() const
		{
			return std::string("AgglomeratingSolver: ") + m_pLinOpInverse->config_string();
		}

		virtual bool init(SmartPtr<ILinearOperator<vector_type> > A)
		{
			return base_type::base_init(A);
		}
		virtual bool init(SmartPtr<ILinearOperator<vector_type> > A, const vector_type& u)
		{
			return base_type::base_init(A);
		}
		virtual bool init(SmartPtr<MatrixOperator<matrix_type, vector_type> > Op)
		{
			return base_type::base_init(Op);
		}


	protected:
		SmartPtr<ILinearOperatorInverse<vector_type, vector_type> > m_pLinOpInverse;
		std::string m_name;
};



template <typename TAlgebra>
class AgglomeratingIterator : public
	AgglomeratingBase<ILinearIterator<typename TAlgebra::vector_type>, TAlgebra >
{
	public:
	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	// 	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Base type
		typedef AgglomeratingBase<ILinearIterator<vector_type>, algebra_type > base_type;

	public:
		AgglomeratingIterator(SmartPtr<ILinearIterator<vector_type> > splinIt)
		{
			UG_COND_THROW(splinIt.valid()==false, "linOpInverse has to be != NULL");
			m_splinIt = splinIt;
			m_name = std::string("AgglomeratingIterator(") + splinIt->name() + std::string(")");
		}

		virtual bool init_agglomerated(SmartPtr<MatrixOperator<matrix_type, vector_type> > Op)
		{
			try{
				return m_splinIt->init(Op);
			}UG_CATCH_THROW("AgglomeratingIterator::" << __FUNCTION__ << " failed")
		}
		virtual bool apply_agglomerated(vector_type& x, const vector_type& b)
		{
			try{
				return m_splinIt->apply(x, b);
			}UG_CATCH_THROW("AgglomeratingIterator::" << __FUNCTION__ << " failed")
		}

	// 	Destructor
		virtual ~AgglomeratingIterator() {};


		virtual const char* name() const
		{
			return m_name.c_str();
		}

		virtual bool supports_parallel() const { return true; }

		///	Clone
		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			SmartPtr<ILinearIterator<vector_type> > linIt = m_splinIt->clone();
			return make_sp(new AgglomeratingIterator<algebra_type>(linIt));
		}

		virtual std::string config_string() const
		{
			return std::string("AgglomeratingIterator: ") + m_splinIt->config_string();
		}


		virtual bool init(SmartPtr<ILinearOperator<vector_type> > A)
		{
			return base_type::base_init(A);
		}
		virtual bool init(SmartPtr<ILinearOperator<vector_type> > A, const vector_type& u)
		{
			return base_type::base_init(A);
		}


	protected:
		SmartPtr<ILinearIterator<vector_type> > m_splinIt;
		std::string m_name;
};

template <typename TAlgebra>
class AgglomeratingPreconditioner: public
	AgglomeratingBase<IPreconditioner<TAlgebra>, TAlgebra >
{
	public:
	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	// 	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Base type
		typedef AgglomeratingBase<IPreconditioner<TAlgebra>, TAlgebra > base_type;

	public:
		AgglomeratingPreconditioner(SmartPtr<ILinearIterator<vector_type> > splinIt)
		{
			UG_COND_THROW(splinIt.valid()==false, "linOpInverse has to be != NULL");
			m_splinIt = splinIt;
			m_name = std::string("AgglomeratingPreconditioner(") + splinIt->name() + std::string(")");
		}

		virtual bool init_agglomerated(SmartPtr<MatrixOperator<matrix_type, vector_type> > Op)
		{
			try{
				return m_splinIt->init(Op);
			}UG_CATCH_THROW("AgglomeratingIterator::" << __FUNCTION__ << " failed")
		}
		virtual bool apply_agglomerated(vector_type& x, const vector_type& b)
		{
			try{
				return m_splinIt->apply(x, b);
			}UG_CATCH_THROW("AgglomeratingIterator::" << __FUNCTION__ << " failed")
		}

	// 	Destructor
		virtual ~AgglomeratingPreconditioner() {};


		virtual const char* name() const
		{
			return m_name.c_str();
		}

		virtual bool supports_parallel() const { return true; }

		///	Clone
		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			SmartPtr<ILinearIterator<vector_type> > linIt = m_splinIt->clone();
			return make_sp(new AgglomeratingIterator<algebra_type>(linIt));
		}

		virtual std::string config_string() const
		{
			return std::string("AgglomeratingPreconditioner: ") + m_splinIt->config_string();
		}

		virtual bool postprocess() { return true; }
		virtual bool preprocess(SmartPtr<MatrixOperator<typename TAlgebra::matrix_type, typename TAlgebra::vector_type> > pOp)
		{
			return base_type::base_init(pOp);
		}


	protected:
		SmartPtr<ILinearIterator<vector_type> > m_splinIt;
		std::string m_name;
};



} // end namespace ug

#endif /* __H__LIB_ALGEBRA__LAPACK_LU_OPERATOR__ */
