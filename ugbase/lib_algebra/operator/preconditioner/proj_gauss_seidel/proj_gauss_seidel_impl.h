/*
 * proj_gauss_seidel_impl.h
 *
 *  Created on: 10.10.2013
 *      Author: raphaelprohl
 */

#ifndef __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJ_GAUSS_SEIDEL__PROJ_GAUSS_SEIDEL_IMPL__
#define __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJ_GAUSS_SEIDEL__PROJ_GAUSS_SEIDEL_IMPL__

#include "proj_gauss_seidel.h"

#define PROFILE_PROJ_GS
#ifdef PROFILE_PROJ_GS
	#define PROJ_GS_PROFILE_FUNC()		PROFILE_FUNC()
	#define PROJ_GS_PROFILE_BEGIN(name)	PROFILE_BEGIN_GROUP(name, "projGS")
	#define PROJ_GS_PROFILE_END()		PROFILE_END()
#else
	#define PROJ_GS_PROFILE_FUNC()
	#define PROJ_GS_PROFILE_BEGIN(name)
	#define PROJ_GS_PROFILE_END()
#endif

namespace ug{

template <typename TAlgebra>
bool
ProjGaussSeidel<TAlgebra>::
preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp)
{
#ifdef UG_PARALLEL
	if(pcl::GetNumProcesses() > 1)
	{
		//	copy original matrix
		MakeConsistent(*pOp, *m_spMat);
		//	set zero on slaves
		std::vector<IndexLayout::Element> vIndex;
		CollectUniqueElements(vIndex, (*m_spMat).layouts()->slave());
		SetDirichletRow(*m_spMat, vIndex);
	}
	matrix_type &A = *m_spMat;
#else
	matrix_type &A = *pOp;
#endif

	CheckDiagonalInvertible(A);
	return true;
}

template <typename TAlgebra>
bool
ProjGaussSeidel<TAlgebra>::
init(SmartPtr<ILinearOperator<vector_type> > J, const vector_type& u)
{
	PROFILE_FUNC_GROUP("projGS");

	//try{

// 	cast operator and remember it
	m_spMat = J.template cast_dynamic<matrix_type>();

//	check that operator type is correct
	if(m_spMat.invalid())
		UG_THROW("ProjGaussSeidel:init: Can not cast Operator to Matrix.");

//	preprocess
	if(!preprocess(m_spMat))
	{
		UG_LOG("ERROR in '"<<name()<<"::init': preprocess failed.\n");
		return false;
	}

//	remember solution //TODO: resize, ect. necessary?
	m_lastSol = u;

//	remember, that operator has been initialized
	m_bInit = true;

	//} UG_CATCH_THROW("ProjGaussSeidel: Init failure for init(u)");
	return true;
}

template <typename TAlgebra>
bool
ProjGaussSeidel<TAlgebra>::
init(SmartPtr<ILinearOperator<vector_type> > L)
{
	PROFILE_FUNC_GROUP("projGS");
	UG_THROW("Solution u is not set in ProjGaussSeidel:init(L)! "
			"Use ProjGaussSeidel:init(J,u) instead! \n");

	return false;
}

template <typename TAlgebra>
void
ProjGaussSeidel<TAlgebra>::
gs_step_with_projection(vector_type& c, const matrix_type& A, const vector_type& d)
{
	//	TODO: matrix_type correct or MatrixOperator<matrix_type, vector_type>?
	typename vector_type::value_type s;
	typename vector_type::value_type tmpSol;

	for(size_t i = 0; i < c.size(); i++)
	{
		s = d[i];

		for(typename matrix_type::const_row_iterator it = A.begin_row(i);
				it != A.end_row(i) && it.index() < i; ++it)
			// s -= it.value() * x[it.index()];
			MatMultAdd(s, 1.0, s, -1.0, it.value(), c[it.index()]);

		//	c[i] = s / A(i,i)
		InverseMatMult(c[i], 1.0, A(i,i), s);

		//	compute temporary solution tmpSol := u_{s-1/2} = u_{s-1} + c
		tmpSol = m_lastSol[i] + c[i];

		//	perform projection: check whether the temporary solution u_{s-1/2}
		//	fulfills the underlying constraint or not
		if (tmpSol < 0.0)
		{
			//	a constraint is not fulfilled

			//	adjust correction c := u_s - u_{s-1} = 0 - u_{s-1}
			c[i] = - m_lastSol[i];

			//	set new solution u_s to zero
			//	and store the current index in a vector for further treatment
			m_lastSol[i] = 0.0;
			m_vActiveIndices.push_back(i);
		}
		else
		{
			//	the 'tmpSol' is valid with respect to all constraints
			m_lastSol[i] = tmpSol;
			m_vInactiveIndices.push_back(i);
		}
	}
}

template <typename TAlgebra>
bool
ProjGaussSeidel<TAlgebra>::
apply(vector_type &c, const vector_type& d)
{
	PROFILE_FUNC_GROUP("projGS");
//	Check that operator is initialized
	if(!m_bInit)
	{
		UG_LOG("ERROR in '"<<name()<<"::apply': Iterator not initialized.\n");
		return false;
	}

	//	Check parallel status
#ifdef UG_PARALLEL
	if(!d.has_storage_type(PST_ADDITIVE))
		UG_THROW(name() << "::apply: Wrong parallel "
					   "storage format. Defect must be additive.");
#endif

	//	check sizes
	if(d.size() != m_spMat->num_rows())
		UG_THROW("Vector [size= "<<d.size()<<"] and Row [size= "
					   <<m_spMat->num_rows()<<"] sizes have to match!");
	if(c.size() != m_spMat->num_cols())
		UG_THROW("Vector [size= "<<c.size()<<"] and Column [size= "
					   <<m_spMat->num_cols()<<"] sizes have to match!");
	if(d.size() != c.size())
		UG_THROW("Vector [d size= "<<d.size()<<", c size = "
					   <<c.size()<< "] sizes have to match!");
	if(d.size() != m_lastSol.size())
			UG_THROW("Vector [d size= "<<d.size()<<", m_lastSol size = "
					   <<m_lastSol.size()<< "] sizes have to match!");
	//	reset vectors
	m_vActiveIndices.resize(0); m_vInactiveIndices.resize(0);

	//	perform a forward GaussSeidel-step: c = (D-L)^-1 * d, project on the underlying constraint
	//	and store the indices, which satisfy the constraint with equality, in the vector vActiveIndices

#ifdef UG_PARALLEL
	if(pcl::GetNumProcesses() > 1)
	{
	//	make defect unique
		SmartPtr<vector_type> spDtmp = d.clone();
		spDtmp->change_storage_type(PST_UNIQUE);

		gs_step_with_projection(c, *m_spMat, *spDtmp);
		c.set_storage_type(PST_UNIQUE);
	}
	else
#endif
	{
		gs_step_with_projection(c, *m_spMat, d);
#ifdef UG_PARALLEL
		c.set_storage_type(PST_UNIQUE);
#endif
	}

//	apply scaling
	const number kappa = this->damping()->damping(c, d, m_spMat.template cast_dynamic<ILinearOperator<vector_type> >());
	if(kappa != 1.0){
		//	damp only those indices which are inactive
		for (std::vector<size_t>::iterator itInactiveInd = m_vInactiveIndices.begin();
				itInactiveInd <  m_vInactiveIndices.end(); ++itInactiveInd)
			c[*itInactiveInd] *= kappa;
	}

//	Correction is always consistent
#ifdef 	UG_PARALLEL
	if(!c.change_storage_type(PST_CONSISTENT))
		UG_THROW(name() << "::apply': Cannot change "
				"parallel storage type of correction to consistent.");
#endif

	return true;
}

template <typename TAlgebra>
bool
ProjGaussSeidel<TAlgebra>::
apply_update_defect(vector_type &c, vector_type& d)
{
	PROFILE_FUNC_GROUP("projGS");

	//	set correction to zero
	c.set(0.0);

	//	compute new correction and perform the projection
	//	(adjusting the solution m_lastSol to the underlying constraint)
	if(!apply(c, d)) return false;

	// 	update defect d := d - A*c = b - A*(x+c) (= b - A*x_new)
	if(!m_spMat->matmul_minus(d, c))
	{
		UG_LOG("ERROR in '"<<name()<<"::apply_update_defect': "
				"Cannot execute matmul_minus to compute d:=d-A*c.\n");
		return false;
	}

	//	adjust defect due to constraint/active indices
	for (std::vector<size_t>::iterator itActiveInd = m_vActiveIndices.begin();
			itActiveInd < m_vActiveIndices.end(); ++itActiveInd)
	{
		if (d[*itActiveInd] > 0.0)
			d[*itActiveInd] = 0.0;
	}

	return true;
}

} // end namespace ug

#endif /* __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJ_GAUSS_SEIDEL__PROJ_GAUSS_SEIDEL_IMPL__ */
