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
	if(m_famg.m_writeMatrices && A.num_rows() < AMG_WRITE_MATRICES_MAX)
	{
		AMG_PROFILE_FUNC();
		UG_DLOG(LIB_ALG_AMG, 1, "write matrices");

		write_debug_matrix(A, level, level, "AMG_A");			UG_DLOG(LIB_ALG_AMG, 1, ".");
		write_debug_matrix(P, level+1, level, "AMG_P");			UG_DLOG(LIB_ALG_AMG, 1, ".");
		write_debug_matrix(R, level, level+1, "AMG_R");			UG_DLOG(LIB_ALG_AMG, 1, ".");

		AMGWriteToFile(AH, level+1, level+1, GetProcFilename(m_famg.m_writeMatrixPath,
				ToString("AMG_A") + ToString(level+1),".mat").c_str(), m_famg.m_amghelper);

		UG_DLOG(LIB_ALG_AMG, 1, ". done.\n");
	}
}

template<typename matrix_type, typename prolongation_matrix_type, typename vector_type>
void FAMGLevelCalculator<matrix_type, prolongation_matrix_type, vector_type>::write_debug_matrix_markers()
{
	if(m_famg.m_writeMatrices)
	{
		AMG_PROFILE_FUNC();
		std::fstream ffine(GetProcFilename(m_famg.m_writeMatrixPath, std::string("AMG_fine") + ToString(level), ".marks").c_str(), std::ios::out);
		std::fstream fcoarse(GetProcFilename(m_famg.m_writeMatrixPath, std::string("AMG_coarse") + ToString(level), ".marks").c_str(), std::ios::out);
		std::fstream fother(GetProcFilename(m_famg.m_writeMatrixPath, std::string("AMG_other") + ToString(level), ".marks").c_str(), std::ios::out);
		std::fstream fdirichlet(GetProcFilename(m_famg.m_writeMatrixPath, std::string("AMG_dirichlet") + ToString(level), ".marks").c_str(), std::ios::out);
		for(size_t i=0; i < rating.size(); i++)
		{
			//int o = m_famg.m_amghelper.GetOriginalIndex(level, i);
			int o=i;
			if(rating[i].is_fine()) ffine << o << "\n";
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
		f2 << "c " << GetProcFilename(m_famg.m_writeMatrixPath, std::string("AMG_fine") + ToString(level), ".marks") << "\n";
		f2 << "c " << GetProcFilename(m_famg.m_writeMatrixPath, std::string("AMG_coarse") + ToString(level), ".marks") << "\n";
		f2 << "c " << GetProcFilename(m_famg.m_writeMatrixPath, std::string("AMG_other") + ToString(level), ".marks") << "\n";
		f2 << "c " << GetProcFilename(m_famg.m_writeMatrixPath, std::string("AMG_dirichlet") + ToString(level), ".marks") << "\n";
	}
}

}
#endif
