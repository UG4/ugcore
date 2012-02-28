/*
 * collected_solver.h
 *
 *  Created on: 04.05.2011
 *      Author: mrupp
 */

#ifndef COLLECTED_SOLVER_H_
#define COLLECTED_SOLVER_H_

#ifndef UG_PARALLEL
#error "This only works with a UG_PARALLEL define."
#endif


template <typename TAlgebra>
class CollectSolver: public IMatrixOperatorInverse<	typename TAlgebra::vector_type,
														typename TAlgebra::vector_type,
														typename TAlgebra::matrix_type>
{
	public:
	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	// 	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

		typedef IMatrixOperatorInverse<	typename TAlgebra::vector_type, typename TAlgebra::vector_type, typename TAlgebra::matrix_type> moi_type;

	public:
		LUSolver() :
			m_pConvCheck(NULL)
		{
			sName = "CollectedSolver(uninitialized)"
		};

		virtual const char* name() const {return name.c_str();}

		void set_convergence_check(IConvergenceCheck& convCheck)
		{
			UG_ASSERT(m_pInverse != NULL, "needs Inverse first");
			m_pInverse->set_convergence_check(convCheck);
		}
		IConvergenceCheck* get_convergence_check()
		{
			UG_ASSERT(m_pInverse != NULL, "needs Inverse first");
			return m_pInverse->get_convergence_check(convCheck);
		}

		void set_inverse(moi_type *inverse)
		{
			m_pInverse = inverse;
			sName = "CollectedSolver(";
			sName.append(inverse->name());
			sName.append(")");
		}

		moi_type *get_inverse()
		{
			return m_pInverse;
		}


		//	set operator L, that will be inverted
		virtual bool init(MatrixOperator<vector_type, vector_type, matrix_type>& Op)
		{
		// 	remember operator
			m_pOperator = &Op;

		//	get matrix of Operator
			m_pMatrix = &m_pOperator->get_matrix();

		//	check that matrix exist
			if(m_pMatrix == NULL)
				{UG_LOG("ERROR in CollectSolver::init: No Matrix given,\n"); return false;}


			// create basesolver
			collect_matrix(*m_pMatrix, collectedMatrix, masterColl, slaveColl);

			m_SMO.setmatrix(&collectedMatrix);
			m_pInverse->init(m_SMO);
		}


		virtual bool apply_return_defect(vector_type& u, vector_type& f)
		{
			vector_type collC;
			vector_type collD;
			if(pcl::GetProcRank() == 0)
			{
				size_t N = collectedBaseA.num_rows();
				collC.resize(N);
				collC.set(0.0);

				collD = d;
				collD.resize(N);
				for(size_t i=Ah.num_rows(); i<N; i++)
					collD[i] = 0.0;

			}


			// send d -> collD
			ComPol_VecAdd<vector_type > compolAdd(&collD, &d);
			pcl::ParallelCommunicator<IndexLayout> &com = m_A[level]->communicator();
			com.send_data(slaveColl, compolAdd);
			com.receive_data(masterColl, compolAdd);
			com.communicate();

			if(pcl::GetProcRank() == 0)
				m_pInverse->apply_return_defect(collC, collD);

			// send c -> collC
			ComPol_VecCopy<vector_type> compolCopy(&c, &collC);
			com.send_data(masterColl, compolCopy);
			com.receive_data(slaveColl, compolCopy);
			com.communicate();

			if(pcl::GetProcRank() == 0)
				for(size_t i=0; i<Ah.num_rows(); i++)
					c[i] = collC[i];
			d.set(0.0);

			return true;
		}

	// 	Destructor
		virtual ~LUSolver() {};

	private:
		std::string sName;
		moi_type *m_pInverse;
};

#endif /* COLLECTED_SOLVER_H_ */
