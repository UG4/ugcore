/*
 * super_lu.cpp
 *
 *  Created on: 16.01.2014
 *      Author: mrupp
 */

#include "super_lu.h"
#ifdef UG_SUPERLU
// fix warning
#undef TRUE
#undef FALSE
#include "slu_ddefs.h"

namespace ug{


class SuperLUImplementation : public IExternalSolverImplementation
{
	bool m_bInited;
	std::vector<double> rhs, nzval;
	/* row permutations from partial pivoting */
	/* column permutation vector */
	std::vector<int> perm_r, perm_c, colind, rowptr;

	SuperMatrix SuperLU_A, SuperLU_L, SuperLU_U, SuperLU_B;

	SuperLUConfiguration &config;

public:
	SuperLUImplementation(SuperLUConfiguration &_config) : config(_config)
	{
		m_bInited=false;
	}

	~SuperLUImplementation()
	{
		destroy();
	}
	void destroy()
	{
		if(m_bInited)
		{
			Destroy_CompCol_Matrix(&SuperLU_A);
			memset(&SuperLU_A, 0, sizeof(SuperMatrix));
			Destroy_SuperMatrix_Store(&SuperLU_B);
			memset(&SuperLU_B, 0, sizeof(SuperMatrix));
			Destroy_SuperNode_Matrix(&SuperLU_L);
			memset(&SuperLU_L, 0, sizeof(SuperMatrix));
			Destroy_CompCol_Matrix(&SuperLU_U);
			memset(&SuperLU_U, 0, sizeof(SuperMatrix));
			m_bInited = false;
		}
	}

	void get_options(superlu_options_t opt)
	{
		//opt.PrintStat = config.bPrintStat ? YES : NO;
		opt.Equil = config.equil ? YES : NO;
		switch(config.colPerm)
		{
		case SuperLUConfiguration::CPT_NATURAL: opt.ColPerm = NATURAL; break;
		case SuperLUConfiguration::CPT_MMD_ATA: opt.ColPerm = MMD_ATA; break;
		case SuperLUConfiguration::CPT_MMD_AT_PLUS_A: opt.ColPerm = MMD_AT_PLUS_A; break;
		case SuperLUConfiguration::CPT_COLAMD: opt.ColPerm = COLAMD; break;
		}
	}

	virtual void init(const CPUAlgebra::matrix_type &A)
	{

//		destroy();
		PROFILE_BEGIN_GROUP(SuperLU_Preprocess, "algebra SuperLU");
		typedef CPUAlgebra::matrix_type::const_row_iterator row_it;
		typedef CPUAlgebra::matrix_type::value_type value_type;
		size_t nnz =0;
		size_t N = A.num_rows();


		for(size_t r=0; r<N; r++)
		{
			for(row_it it = A.begin_row(r); it != A.end_row(r); ++it)
				nnz++;
		}
		if(N > 40000 && nnz > 400000) { UG_LOG("SuperLU preprocess, N = " << N << ", nnz = " << nnz << "... "); }

		nzval.resize(nnz);
		colind.resize(nnz);
		rowptr.resize(N+1);

		size_t colIndPos=0;
		for(size_t r=0; r<N; r++)
		{
			rowptr[r] = colIndPos;
			for(row_it it = A.begin_row(r); it != A.end_row(r); ++it)
			{
				nzval[colIndPos] = it.value();
				colind[colIndPos] = it.index();
				colIndPos++;
			}
			rowptr[r+1] = colIndPos;
		}

//			PrintVector(nzval, nnz, "nzval");
//			PrintVector(colind, nnz, "colind");
//			PrintVector(rowptr, N+1, "rowptr");

		superlu_options_t options;
		SuperLUStat_t stat;

		dCreate_CompRow_Matrix(&SuperLU_A, N, N, nnz, &nzval[0], &colind[0], &rowptr[0], SLU_NR, SLU_D, SLU_GE);

		//rhs = doubleMalloc(N);
		rhs.resize(N);
		for (size_t i = 0; i < N; ++i)
			rhs[i] = 0;
		dCreate_Dense_Matrix(&SuperLU_B, N, 1, &rhs[0], N, SLU_DN, SLU_D, SLU_GE);

		perm_r.resize(N+1);
		perm_c.resize(N+1);
		/* Set the default input options. */
		set_default_options(&options);

		get_options(options);

		StatInit(&stat);
		int info;


//			A.print("A");

		dgssv(&options, &SuperLU_A, &perm_c[0], &perm_r[0], &SuperLU_L, &SuperLU_U, &SuperLU_B, &stat, &info);

		dgssv_check_info(info, N);
		if(config.bPrintStat)
			StatPrint(&stat);
		StatFree(&stat);
		if(N > 40000 && nnz > 400000) { UG_LOG("done.\n"); }
		m_bInited = true;

	}

	void dgssv_check_info(int info, size_t N)
	{
		if(info > 0)
		{
			if(info < (int)N)
			{
				UG_THROW("ERROR in SuperLU: U(i,i) with i=" << info << "is exactly zero. The factorization has\
							been completed, but the factor U is exactly singular,\
							so the solution could not be computed.");
			}
			else
			{ UG_THROW("ERROR in SuperLU: memory allocation failure");}
		}
		else if(info < 0)
		{
			UG_THROW("ERROR in SuperLU: info < 0 ???");
		}
	}

	virtual void apply(CPUAlgebra::vector_type &c, const CPUAlgebra::vector_type &d)
	{
		PROFILE_BEGIN_GROUP(SuperLU_Apply, "algebra SuperLU");
		size_t N = c.size();
		double *b = (double*) ((DNformat*) SuperLU_B.Store)->nzval;

		for (size_t i = 0; i < N; ++i)
			b[i] = d[i];

		superlu_options_t options;
		set_default_options(&options);
		//options.Fact = FACTORED;
		int info;
		SuperLUStat_t stat;
		StatInit(&stat);
		dgssv(&options, &SuperLU_A, &perm_c[0], &perm_r[0], &SuperLU_L, &SuperLU_U, &SuperLU_B, &stat, &info);
		StatFree(&stat);
		dgssv_check_info(info, N);

		for (size_t i = 0; i < N; ++i)
			c[i] = b[i];
	}

	virtual const char* name() const { return "SuperLU"; }
};


IExternalSolverImplementation *CreateSuperLUImplementation(SuperLUConfiguration &config)
{
	return new SuperLUImplementation(config);
}

}
#endif
