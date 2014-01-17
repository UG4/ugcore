/*
 * \file super_lu.h
 *
 * http://crd-legacy.lbl.gov/~xiaoye/SuperLU
 *
 * \date 20.07.2013
 * \author Martin Rupp
 */
#ifdef UG_SUPERLU
#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__SUPER_LU_SOLVER_
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__SUPER_LU_SOLVER_

#include "external_solvers.h"
#include <vector>


namespace ug{

struct SuperLUConfiguration
{
	bool bPrintStat;
	bool equil;
	enum
	{
	CPT_NATURAL, CPT_MMD_ATA, CPT_MMD_AT_PLUS_A, CPT_COLAMD
	} colPerm;
};


IExternalSolverImplementation *CreateSuperLUImplementation(SuperLUConfiguration &config);


template<typename TAlgebra>
class SuperLUSolver : public IExternalSolver<TAlgebra>
{
	SuperLUConfiguration config;
	using IExternalSolver<TAlgebra>::impl;
public:
	SuperLUSolver()
	{
		config.bPrintStat = false;
		config.equil = true;
		config.colPerm = SuperLUConfiguration::CPT_COLAMD;

		impl = CreateSuperLUImplementation(config);
	}

	void print_stat(bool b)
	{
		config.bPrintStat = b;
	}

	void equil(bool b)
	{
		config.equil = b;
	}

	void col_perm_natural()
	{
		config.colPerm = SuperLUConfiguration::CPT_NATURAL;
	}

	void col_perm_mdo_ATA()
	{
		config.colPerm = SuperLUConfiguration::CPT_MMD_ATA;
	}

	void col_perm_mdo_AT_plus_A()
	{
		config.colPerm = SuperLUConfiguration::CPT_MMD_AT_PLUS_A;
	}

	void col_perm_approximate()
	{
		config.colPerm = SuperLUConfiguration::CPT_COLAMD;
	}
};

} // end namespace ug

#endif
#endif
