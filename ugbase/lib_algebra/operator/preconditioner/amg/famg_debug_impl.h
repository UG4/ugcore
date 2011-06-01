/**
 * \file famg_debug_impl.h
 *
 * \author Martin Rupp
 *
 * \date 14.02.2011
 *
 * implementation file for famg debug output
 *
 * Goethe-Center for Scientific Computing 2011.
 *
 */


#ifndef __H__LIB_ALGEBRA__AMG__FAMG_DEBUG_IMPL_H__
#define __H__LIB_ALGEBRA__AMG__FAMG_DEBUG_IMPL_H__


namespace ug
{


template<typename matrix_type, typename prolongation_matrix_type, typename vector_type>
void FAMGLevelCalculator<matrix_type, prolongation_matrix_type, vector_type>::write_debug_matrices()
{
	//------------------------------------
	if(m_famg.m_amghelper.positions[level].size() > 0)
	{
		m_famg.m_amghelper.positions.resize(level+2);
		m_famg.m_amghelper.positions[level+1].resize(AH.num_rows());
		for(size_t i=0; i < AH.num_rows(); i++)
			m_famg.m_amghelper.positions[level+1][i] = m_famg.m_amghelper.positions[level][m_famg.m_parentIndex[level+1][i]];
	}

	if(m_famg.m_writeMatrices && A.num_rows() < AMG_WRITE_MATRICES_MAX)
	{
		AMG_PROFILE_FUNC();
		UG_DLOG(LIB_ALG_AMG, 1, "write matrices");

		write_debug_matrix(A, level, level, "AMG_A_L");			UG_DLOG(LIB_ALG_AMG, 1, ".");
		write_debug_matrix(PoldIndices, level, level, "AMG_P_L");			UG_DLOG(LIB_ALG_AMG, 1, ".");
		write_debug_matrix(R, level, level+1, "AMG_R_L");			UG_DLOG(LIB_ALG_AMG, 1, ".");

		AMGWriteToFile(AH, level+1, level+1, GetProcFilename(m_famg.m_writeMatrixPath,
				ToString("AMG_A_L") + ToString(level+1),".mat").c_str(), m_famg.m_amghelper);

		UG_DLOG(LIB_ALG_AMG, 1, ". done.\n");
	}
}

template<typename matrix_type, typename prolongation_matrix_type, typename vector_type>
void FAMGLevelCalculator<matrix_type, prolongation_matrix_type, vector_type>::write_debug_matrix_markers()
{
	if(m_famg.m_writeMatrices)
	{
		AMG_PROFILE_FUNC();
		std::fstream ffine(GetProcFilename(m_famg.m_writeMatrixPath, std::string("AMG_fine_L") + ToString(level), ".marks").c_str(), std::ios::out);
		ffine << "1 0 0 1 0\n";
		std::fstream ffine2(GetProcFilename(m_famg.m_writeMatrixPath, std::string("AMG_aggfine_L") + ToString(level), ".marks").c_str(), std::ios::out);
		ffine2 << "1 0.2 1 1 0\n";
		std::fstream fcoarse(GetProcFilename(m_famg.m_writeMatrixPath, std::string("AMG_coarse_L") + ToString(level), ".marks").c_str(), std::ios::out);
		fcoarse << "0 0 1 1 2\n";
		std::fstream fother(GetProcFilename(m_famg.m_writeMatrixPath, std::string("AMG_other_L") + ToString(level), ".marks").c_str(), std::ios::out);
		fother << "1 1 0 1 0\n";
		std::fstream fdirichlet(GetProcFilename(m_famg.m_writeMatrixPath, std::string("AMG_dirichlet_L") + ToString(level), ".marks").c_str(), std::ios::out);
		fdirichlet << "0 1 1 1 0\n";
		for(size_t i=0; i < rating.size(); i++)
		{
			//int o = m_famg.m_amghelper.GetOriginalIndex(level, i);
			int o=i;
			if(rating[i].is_aggressive_fine()) ffine2 << o << "\n";
			else if(rating[i].is_fine()) ffine << o << "\n";
			else if(rating[i].is_coarse()) fcoarse << o << "\n";
			else if(rating[i].is_dirichlet()) fdirichlet << o << "\n";
			else fother << o << "\n";
		}
	}

}

template<typename matrix_type, typename prolongation_matrix_type, typename vector_type>
template<typename TMatrix>
void FAMGLevelCalculator<matrix_type, prolongation_matrix_type, vector_type>::write_debug_matrix(TMatrix &mat, size_t fromlevel, size_t tolevel, const char *name)
{
	if(m_famg.m_writeMatrices)
	{
		AMG_PROFILE_FUNC();
		std::string filename = GetProcFilename(m_famg.m_writeMatrixPath, ToString(name) + ToString(fromlevel),".mat");
		AMGWriteToFile(mat, fromlevel, tolevel, filename.c_str(), m_famg.m_amghelper);
		std::fstream f2(filename.c_str(), std::ios::out | std::ios::app);
		f2 << "c " << GetProcFilename(m_famg.m_writeMatrixPath, std::string("AMG_fine_L") + ToString(level), ".marks") << "\n";
		f2 << "c " << GetProcFilename(m_famg.m_writeMatrixPath, std::string("AMG_aggfine_L") + ToString(level), ".marks") << "\n";
		f2 << "c " << GetProcFilename(m_famg.m_writeMatrixPath, std::string("AMG_coarse_L") + ToString(level), ".marks") << "\n";
		f2 << "c " << GetProcFilename(m_famg.m_writeMatrixPath, std::string("AMG_other_L") + ToString(level), ".marks") << "\n";
		f2 << "c " << GetProcFilename(m_famg.m_writeMatrixPath, std::string("AMG_dirichlet_L") + ToString(level), ".marks") << "\n";
	}
}

}
#endif
