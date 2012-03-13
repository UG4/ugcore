/*
 * lu_operator.h
 *
 *  Created on: 16.06.2010
 *      Author: mrupp
 */

#ifndef __H__LIB_ALGEBRA__LAPACK_LU_OPERATOR__
#define __H__LIB_ALGEBRA__LAPACK_LU_OPERATOR__
#include <iostream>
#include <sstream>

#include "lib_algebra/operator/operator_inverse_interface.h"

#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

#define PCL_DT_BYTE 			MPI_BYTE
#define PCL_DT_PACKED 			MPI_PACKED
#define PCL_DT_CHAR 			MPI_CHAR
#define PCL_DT_SHORT 			MPI_SHORT
#define PCL_DT_INT 				MPI_INT
#define PCL_DT_LONG 			MPI_LONG
#define PCL_DT_FLOAT 			MPI_FLOAT
#define PCL_DT_DOUBLE 			MPI_DOUBLE
#define PCL_DT_LONG_DOUBLE 		MPI_LONG_DOUBLE
#define PCL_DT_UNSIGNED_CHAR 	MPI_UNSIGNED_CHAR

void PCLAllReduce(ProcessCommunicator &comm, size_t &src, size_t &dest, ReduceOperation op);
{
	comm.allreduce(&src, &dest, 1, PCL_DT_UNSIGNED_LONG, op);
}

template<typename TLayout>
void GenerateConsecutiveGlobalAlgebraIDs(std::vector<AlgebraID>& globalIDs,
		std::vector<size_t>& consecutiveIDs,
							  size_t numIDs,
							  TLayout& masterLayout,
							  TLayout& slaveLayout)
{
	size_t local = 0;
	for(size_t i=0; i<globalIDs.size(); i++)
		if(globalIDs[i].first == pcl::GetProcRank())
			local++;

	size_t global = 0;
	PCLAllReduce(local, global, PCL_RO_SUM);



}

template<typename TMatrix>
bool GetCompleteMatrixOnProcessor0(TMatrix &matrix, TMatrix &completeMatrix)
{
	std::vector<AlgebraID> globalIDs;
	GenerateGlobalAlgebraIDs(globalIDs, std::max(matrix.num_rows(), matrix.num_cols()),
								matrix.get_master_layout(), matrix.get_slave_layout());

	ProcessCommunicator &comm = matrix.get_communicator();
	size_t localIndices = 0;
		for(size_t i=0; i<globalIDs; i++)
			if(globalIDs[i].second ==


	comm.allreduce(&tNormLocal, &tNormGlobal, 1, PCL_DT_INT, PCL_RO_SUM);
	if(pcl::GetProcRank() == 0)
		completeMatrix
}

namespace ug{

template <typename TAlgebra>
class ParallelLUSolver : public IMatrixOperatorInverse<	typename TAlgebra::vector_type,
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

	///	Base type
		typedef IMatrixOperatorInverse<vector_type,vector_type,matrix_type> base_type;

	protected:
		using base_type::convergence_check;

	public:
		LUSolver() :
			m_pOperator(NULL), m_mat()
		{};

		virtual const char* name() const {return "LUSolver";}

		bool init_lu(const matrix_type &A)
		{
			matrix_type allA;
			GetCompleteMatrixOnProcessor0(allA, masterLayout, slaveLayout);

			if(pcl::GetProcId() == 0)
				lu.init_lu(allA);

		}

		bool apply_lu(vector_type &x, const vector_type &b)
		{
			GetCompleteVectorOnProcessor0(allX, x, masterLayout, slaveLayout);
			GetCompleteVectorOnProcessor0(allB, b, masterLayout, slaveLayout);

			if(pcl::GetProcId() == 0)
				lu.apply_lu(u, h);
			return true;
		}

	//	set operator L, that will be inverted
		virtual bool init(IMatrixOperator<vector_type, vector_type, matrix_type>& Op)
		{
		// 	remember operator
			m_pOperator = &Op;

		//	get matrix of Operator
			m_pMatrix = &m_pOperator->get_matrix();

		//	check that matrix exist
			if(m_pMatrix == NULL)
				{UG_LOG("ERROR in LUOperator::init: No Matrix given,\n"); return false;}

		//	init LU operator
			if(!init_lu(*m_pMatrix))
				{UG_LOG("ERROR in LUOperator::init: Cannot init LU Decomposition.\n"); return false;}

		//	we're done
			return true;
		}

	// 	Compute u = L^{-1} * f
		virtual bool apply(vector_type& u, const vector_type& f)
		{
			convergence_check()->set_symbol('%');
			convergence_check()->set_name("LU Solver");

			UG_ASSERT(f.size() == m_pMatrix->num_rows(), "Vector and Row sizes have to match!");
			UG_ASSERT(u.size() == m_pMatrix->num_cols(), "Vector and Column sizes have to match!");
			UG_ASSERT(f.size() == u.size(), "Vector sizes have to match!");

			vector_type allU, allF;
			GetCompleteVectorOnProcessor0(allU, u, masterLayout, slaveLayout);
			GetCompleteVectorOnProcessor0(allF, f, masterLayout, slaveLayout);

			if(pcl::GetProcRank() == 0)
				lu.apply_lu(allU, allF);

			DistributeCompleteVectorFromProcessor0(allU, u, masterLayout, slaveLayout);
			DistributeCompleteVectorFromProcessor0(allF, f, masterLayout, slaveLayout);

		//	we're done
			return true;
		}

	// 	Compute u = L^{-1} * f AND return defect f := f - L*u
		virtual bool apply_return_defect(vector_type& u, vector_type& f)
		{
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

	// 	Destructor
		virtual ~LUSolver() {};

	protected:
		LUSolver<TAlgebra> m_lu;
		// Operator to invert
		IMatrixOperator<vector_type, vector_type, matrix_type>* m_pOperator;

		// matrix to invert
		matrix_type* m_pMatrix;
};

} // end namespace ug

#endif /* __H__LIB_ALGEBRA__LAPACK_LU_OPERATOR__ */
