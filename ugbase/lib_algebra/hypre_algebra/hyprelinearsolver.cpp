
#include "hyprelinearsolver.h"

namespace ug{

HYPREboomerAMG::HYPREboomerAMG()
{
	// Create BoomerAMG
	HYPRE_BoomerAMGCreate(&m_boomer);
	m_hypre_type = HYPRE_INVALID;
}

bool HYPREboomerAMG::setParameters( double tol,
									int MaxIter,
									int MaxLevels,
									double StrongThreshold,
									int CoarsenType,
									int CycleType,
									int NumSweeps,
									int RelaxType,
									int InterpType,
									int SmoothType,
									int PrintLevel	)
{
	m_tol = tol;
	m_MaxIter = MaxIter;
	m_MaxLevels = MaxLevels;
	m_StrongThreshold = StrongThreshold;
	m_CoarsenType = CoarsenType;
	m_CycleType = CycleType;
	m_NumSweeps = NumSweeps;
	m_RelaxType = RelaxType;
	m_InterpType = InterpType;
	m_SmoothType = SmoothType;
	m_PrintLevel = PrintLevel;

	return true;
}

bool HYPREboomerAMG::setType( HYPRE_TYPE type)
{
	m_hypre_type = type;
	return true;
}


bool HYPREboomerAMG::solve(HypreMatrix& A, HypreVector& x, HypreVector& b)
{
	HYPRE_ParCSRMatrix csrpar_A;
	HYPRE_ParVector par_b;
	HYPRE_ParVector par_x;
	int err;

	// Finalize A,x,b
	err = 0;
	err += HYPRE_IJVectorAssemble(b.getStorage());
	err += HYPRE_IJVectorGetObject(b.getStorage(), (void**) &par_b);
	if(err){
		std::cout << "HYPREboomerAMG.solve(): Can not get Vector b" << std::endl;
		return false;
	}
	err = 0;
	err += HYPRE_IJVectorAssemble(x.getStorage());
	err += HYPRE_IJVectorGetObject(x.getStorage(), (void**) &par_x);
	if(err){
		std::cout << "HYPREboomerAMG.solve(): Can not get Vector x" << std::endl;
		return false;
	}
	err = 0;
	err += HYPRE_IJMatrixAssemble(A.getStorage());
	err += HYPRE_IJMatrixGetObject(A.getStorage(), (void**)&csrpar_A);
	if(err){
		std::cout << "HYPREboomerAMG.solve(): Can not get Matrix A" << std::endl;
		return false;
	}


	// Setup boomer
	err = 0;
	err += HYPRE_BoomerAMGSetMaxLevels(m_boomer, m_MaxLevels);
	err += HYPRE_BoomerAMGSetStrongThreshold (m_boomer, m_StrongThreshold);
	err += HYPRE_BoomerAMGSetCoarsenType(m_boomer, m_CoarsenType);
	err += HYPRE_BoomerAMGSetCycleType(m_boomer, m_CycleType);
	err += HYPRE_BoomerAMGSetNumSweeps(m_boomer, m_NumSweeps);
	err += HYPRE_BoomerAMGSetRelaxType(m_boomer, m_RelaxType);
	err += HYPRE_BoomerAMGSetInterpType(m_boomer, m_InterpType);
	err += HYPRE_BoomerAMGSetSmoothType(m_boomer, m_SmoothType);
	err += HYPRE_BoomerAMGSetPrintLevel(m_boomer, m_PrintLevel);
	if(err){
		std::cout << "HYPREboomerAMG.solve(): Can not set up boomer" << std::endl;
		return false;
	}


	if(m_hypre_type == HYPRE_BOOMER)
	{
	err = 0;

	err += HYPRE_BoomerAMGSetTol(m_boomer, m_tol);
	err += HYPRE_BoomerAMGSetMaxIter(m_boomer, m_MaxIter);

	err += HYPRE_BoomerAMGSetup(m_boomer, csrpar_A, par_b, par_x);
	err += HYPRE_BoomerAMGSolve(m_boomer, csrpar_A, par_b, par_x);
	}
	else if(m_hypre_type == HYPRE_PCG)
	// boomer AMG as Proconditioner in pcg
	{
	HYPRE_Solver pcg;

	err += HYPRE_BoomerAMGSetPrintLevel(m_boomer, 1);
	err += HYPRE_BoomerAMGSetTol(m_boomer, 0.0);
	err += HYPRE_BoomerAMGSetMaxIter(m_boomer, 1);

	err += 	HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &pcg);

	err +=	HYPRE_PCGSetMaxIter(pcg, m_MaxIter);
	err += 	HYPRE_PCGSetTol(pcg, m_tol);

	err += 	HYPRE_PCGSetTwoNorm(pcg, 1);
	err += 	HYPRE_PCGSetPrintLevel(pcg, 5);
	err += 	HYPRE_PCGSetLogging(pcg, 1);

		/* Set the PCG preconditioner */
	err += 	HYPRE_PCGSetPrecond(pcg, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
				    (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, m_boomer);

		/* Now setup and solve! */
	err += 	HYPRE_ParCSRPCGSetup(pcg, csrpar_A, par_b, par_x);
	err += 	HYPRE_ParCSRPCGSolve(pcg, csrpar_A, par_b, par_x);

	err += 	HYPRE_ParCSRPCGDestroy(pcg);
	}
	else if(m_hypre_type == HYPRE_BICGSTAB)
	// boomer AMG as Preconditioner in BiCGStab
	{
	HYPRE_Solver bicgstab;

	err += 	HYPRE_BoomerAMGSetTol(m_boomer, 0.0);
	err += 	HYPRE_BoomerAMGSetMaxIter(m_boomer, 1);
	err +=  HYPRE_BoomerAMGSetPrintLevel(m_boomer, 1);

		/* Create solver */
	err += 	HYPRE_ParCSRBiCGSTABCreate(MPI_COMM_WORLD, &bicgstab);

	err += 	HYPRE_BiCGSTABSetMaxIter(bicgstab, m_MaxIter);
	err += 	HYPRE_BiCGSTABSetTol(bicgstab, m_tol);

	err += 	HYPRE_BiCGSTABSetPrintLevel(bicgstab, 5);
	err += 	HYPRE_BiCGSTABSetLogging(bicgstab, 1);

		/* Set the BiCGSTAB preconditioner */
	err += 	HYPRE_BiCGSTABSetPrecond(bicgstab, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
							(HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, m_boomer);


		/* Now setup and solve! */
	err += 	HYPRE_ParCSRBiCGSTABSetup(bicgstab, csrpar_A, par_b, par_x);
	err += 	HYPRE_ParCSRBiCGSTABSolve(bicgstab, csrpar_A, par_b, par_x);

	err += 	HYPRE_ParCSRBiCGSTABDestroy(bicgstab);
	}
	else if(m_hypre_type == HYPRE_GMRES)
	{
		HYPRE_Solver gmres;
		err += 	HYPRE_BoomerAMGSetPrintLevel(m_boomer, 1);
		err += 	HYPRE_BoomerAMGSetTol(m_boomer, 0.0);
		err += 	HYPRE_BoomerAMGSetMaxIter(m_boomer, 1);

		/* Create solver */
		err += 	HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, &gmres);

		err += 	HYPRE_GMRESSetMaxIter(gmres, m_MaxIter);
		err += 	HYPRE_GMRESSetTol(gmres, m_tol);

		err += 	HYPRE_GMRESSetPrintLevel(gmres, 2);
		err += 	HYPRE_GMRESSetLogging(gmres, 1);

		/* Set the GMRES preconditioner */
		err += 	HYPRE_GMRESSetPrecond(gmres, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
							(HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, m_boomer);


        /* Now setup and solve! */
		err += 	HYPRE_ParCSRGMRESSetup(gmres, csrpar_A, par_b, par_x);
		err += 	HYPRE_ParCSRGMRESSolve(gmres, csrpar_A, par_b, par_x);

		err += 	HYPRE_ParCSRGMRESDestroy(gmres);
	}
	else
	{
		std::cout << "HYPREboomerAMG.solve(): Type not specified" << std::endl;
		return false;

	}

	if(err){
		std::cout << "HYPREboomerAMG.solve(): Error while solving" << std::endl;
		return false;
	}


	return !(bool)err;
}

HYPREboomerAMG::~HYPREboomerAMG()
{
	// Delete BoomerAMG
	HYPRE_BoomerAMGDestroy(m_boomer);
}

}
