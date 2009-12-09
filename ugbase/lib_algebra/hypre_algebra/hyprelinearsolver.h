/*
 * linearsolver.h
 *
 *  Created on: 04.07.2009
 *      Author: andreasvogel
 */

#ifndef LINEARSOLVER_H_
#define LINEARSOLVER_H_

#include <HYPRE.h>
#include <_hypre_utilities.h>
#include <HYPRE_krylov.h>
#include <HYPRE_parcsr_ls.h>

#include "hyprevector.h"
#include "hyprematrix.h"

namespace ug{

enum HYPRE_TYPE
{
	HYPRE_INVALID = 0,
	HYPRE_BOOMER,
	HYPRE_PCG,
	HYPRE_BICGSTAB,
	HYPRE_GMRES
};

class HYPREboomerAMG{

	public:
	HYPREboomerAMG();

	bool setParameters( double tol = 1.e-8,
						int MaxIter = 100,
						int MaxLevels = 5,
						double StrongThreshold = 0.25,
						int CoarsenType = 6,
						int CycleType = 1,
						int NumSweeps = 1,
						int RelaxType = 3,
						int InterpType = 0,
						int SmoothType = 6,
						int PrintLevel = 3	);

	bool setType(HYPRE_TYPE type);

	bool solve(HypreMatrix& A, HypreVector& x, HypreVector& b);
	~HYPREboomerAMG();

	private:
		HYPRE_Solver m_boomer;

		double m_tol;
		int m_MaxIter;
		int m_MaxLevels;
		double m_StrongThreshold;
		int m_CoarsenType;
		int m_CycleType;
		int m_NumSweeps;
		int m_RelaxType;
		int m_InterpType;
		int m_SmoothType;
		int m_PrintLevel;

		HYPRE_TYPE m_hypre_type;

};


}

#endif /* LINEARSOLVER_H_ */
