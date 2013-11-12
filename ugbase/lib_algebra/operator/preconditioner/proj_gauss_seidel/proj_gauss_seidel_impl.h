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

// 	init vector of obstacle values and init values with zero
	if (!m_bObs)
		m_spVecOfObsValues = u.clone_without_values();

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
	value_type s, tmpSol, obsVal;
	const number relaxFactor = m_relax;

	for(size_t i = 0; i < c.size(); i++)
	{
		s = d[i];

		for(typename matrix_type::const_row_iterator it = A.begin_row(i);
				it != A.end_row(i) && it.index() < i; ++it)
			// s -= it.value() * x[it.index()];
			MatMultAdd(s, 1.0, s, -1.0, it.value(), c[it.index()]);

		//	c[i] = relaxFactor * s / A(i,i)
		InverseMatMult(c[i], relaxFactor, A(i,i), s);

		//	compute temporary solution (solution of a common (forward) GaussSeidel-step)
		//	tmpSol := u_{s-1/2} = u_{s-1} + c
		tmpSol = m_lastSol[i] + c[i];

		//	get i-th obstacle value
		obsVal = (*m_spVecOfObsValues)[i];

		//	perform projection: check whether the temporary solution u_{s-1/2}
		//	fulfills the underlying constraint or not
		if ((tmpSol - obsVal) < 0.0)
		{
			//	a constraint is not fulfilled

			//	adjust correction c := u_s - u_{s-1} = m_obsVal - u_{s-1}
			// TODO: need to be replaced by c[i] = ConsValue - m_lastSol[i];
			c[i] = (*m_spVecOfObsValues)[i] - m_lastSol[i];

			//	set new solution u_s to the obstacle value
			//	and store the current index in a vector for further treatment
			m_lastSol[i] = (*m_spVecOfObsValues)[i];
			m_vActiveIndices.push_back(i);
		}
		else
		{
			//	the 'tmpSol' is valid with respect to all constraints
			m_lastSol[i] = tmpSol;
			m_vInactiveIndices.push_back(i);
		}

		//UG_LOG("m_lastSol[" << i << "]: " << m_lastSol[i] << "\n");
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

	//	perform a forward GaussSeidel-step: c = (D-L)^-1 * d, then project on the underlying constraint
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

	/*for(size_t i = 0; i < m_spMat->num_rows(); i++)
	{
		for(size_t j = 0; j < m_spMat->num_cols(); j++)
		{
			UG_LOG("A(" << i << "," << j << "]: " << (*m_spMat)(i,j) << "\n");
		}
	}*/

	//	set correction to zero
	c.set(0.0);

	/*for(size_t i = 0; i < m_lastSol.size(); i++)
	{
		UG_LOG("m_lastSol[" << i << "]: " << m_lastSol[i] << "\n");
	}*/

	//	compute new correction and perform the projection
	//	(and update of the solution m_lastSol under consideration of the constraint)
	if(!apply(c, d)) return false;

	/*for(size_t i = 0; i < c.size(); i++)
	{
		UG_LOG("c[" << i << "]: " << c[i] << "\n");
	}

	for(size_t i = 0; i < d.size(); i++)
	{
		UG_LOG("d[" << i << "]: " << d[i] << "\n");
	}*/

	// 	update defect d := d - A*c (= b - A*(x+c) = b - A*x_new)
	if(!m_spMat->matmul_minus(d, c))
	{
		UG_LOG("ERROR in '"<<name()<<"::apply_update_defect': "
				"Cannot execute matmul_minus to compute d:=d-A*c.\n");
		return false;
	}

	/*for(size_t i = 0; i < d.size(); i++)
	{
		UG_LOG("d_nachUpdate[" << i << "]: " << d[i] << "\n");
	}*/

	//	adjust defect of the active indices due to the constraint.
	for (std::vector<size_t>::iterator itActiveInd = m_vActiveIndices.begin();
			itActiveInd < m_vActiveIndices.end(); ++itActiveInd)
	{
		UG_LOG("activeIndex: " << *itActiveInd << "\n");

		//	check, if Ax <= b. For that case the new defect is set to zero,
		//	since all equations/constraints are fulfilled
		if (d[*itActiveInd] < 0.0)
			d[*itActiveInd] = 0.0;
	}

	/*for(size_t i = 0; i < d.size(); i++)
	{
		UG_LOG("d_nachActiveIndUpdate[" << i << "]: " << d[i] << "\n");
	}*/

	return true;
}

} // end namespace ug

#endif /* __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJ_GAUSS_SEIDEL__PROJ_GAUSS_SEIDEL_IMPL__ */
