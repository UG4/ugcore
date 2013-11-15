/*
 * proj_preconditioners_impl.h
 *
 *  Created on: 13.11.2013
 *      Author: raphaelprohl
 */

#ifndef __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_PRECONDS__PROJ_PRECOND_INTERFACE_IMPL__
#define __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_PRECONDS__PROJ_PRECOND_INTERFACE_IMPL__

#include "proj_precond_interface.h"

#define PROFILE_PROJ_PRECOND
#ifdef PROFILE_PROJ_PRECOND
	#define PROJ_PRECOND_PROFILE_FUNC()		PROFILE_FUNC()
	#define PROJ_PRECOND_PROFILE_BEGIN(name)	PROFILE_BEGIN_GROUP(name, "IProjPrecond")
	#define PROJ_PRECOND_PROFILE_END()		PROFILE_END()
#else
	#define PROJ_PRECOND_PROFILE_FUNC()
	#define PROJ_PRECOND_PROFILE_BEGIN(name)
	#define PROJ_PRECOND_PROFILE_END()
#endif

namespace ug{

template <typename TAlgebra>
bool
IProjPreconditioner<TAlgebra>::
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
IProjPreconditioner<TAlgebra>::
init(SmartPtr<ILinearOperator<vector_type> > J, const vector_type& u)
{
	PROFILE_FUNC_GROUP("IProjPrecond");

	//try{

// 	cast operator and remember it
	m_spMat = J.template cast_dynamic<matrix_type>();

//	check that operator type is correct
	if(m_spMat.invalid())
		UG_THROW("ProjGaussSeidel:init: Can not cast Operator to Matrix.");

//	preprocess
	if(!preprocess(m_spMat))
	{
		UG_LOG("ERROR in 'IProjPreconditioner::init': preprocess failed.\n");
		return false;
	}

//	remember solution
	m_lastSol = u.clone();

//	remember, that operator has been initialized
	m_bInit = true;

// 	init vector of obstacle values and init values with zero
	if (!m_bLowerObs)
		m_spVecOfLowObsValues = u.clone_without_values();
	if (!m_bUpperObs)
		m_spVecOfUpObsValues = u.clone_without_values();

//	check, that lower obstacle is <= upper obstacle (for all indices)
	if (m_bLowerObs && m_bUpperObs)
	{
		if ((*m_spVecOfLowObsValues).size() != (*m_spVecOfUpObsValues).size())
			UG_THROW("In IProjPreconditioner::init(J,u) :Vector of lower obstacle values [size= "
					<<(*m_spVecOfLowObsValues).size()<<"] and "
					" Vector of upper obstacle values [size= "
					<<(*m_spVecOfUpObsValues).size()<<"] sizes have to match!");

		for(size_t i = 0; i < (*m_spVecOfLowObsValues).size(); i++)
		{
			value_type lowerObsVal = (*m_spVecOfLowObsValues)[i];
			value_type upperObsVal = (*m_spVecOfUpObsValues)[i];
			for(size_t j = 0; j < GetSize(lowerObsVal); j++)
			{
				if (BlockRef(lowerObsVal, j) - BlockRef(upperObsVal, j) > 0.0)
					UG_THROW("In IProjPreconditioner::init(J,u) " <<i<<"-th index and "<<j<<"-th"
						" component of vector of lower obstacle [value= "<<lowerObsVal<<"] needs "
						"to be lower equal the "<<i<<"-th value of vector of upper obstacle "
						"[value= "<<upperObsVal<<"]!");
			}
		}
	}

//	(ugly) hint, that usual damping (x += damp * c) does not make sense for the projected
//	GaussSeidel-method.
	const number kappa = this->damping()->damping(u, u, m_spMat.template cast_dynamic<ILinearOperator<vector_type> >());
	if(kappa != 1.0){
		UG_THROW("IProjPreconditioner::set_damp': Ususal damping is not possible "
				"for IProjPreconditioner! Use 'set_sor_relax' instead!");
	}

	//} UG_CATCH_THROW("ProjGaussSeidel: Init failure for init(u)");
	return true;
}

template <typename TAlgebra>
bool
IProjPreconditioner<TAlgebra>::
init(SmartPtr<ILinearOperator<vector_type> > L)
{
	PROFILE_FUNC_GROUP("IProjPrecond");
	UG_THROW("Solution u is not set in ProjGaussSeidel::init(L)! "
			"Use ProjGaussSeidel::init(J,u) instead! \n");

	return false;
}

template <typename TAlgebra>
void
IProjPreconditioner<TAlgebra>::
correction_for_lower_obs(vector_type& c, const size_t index, const value_type& tmpSol)
{
	//	get index-th lower obstacle value
	value_type lowerObsVal = (*m_spVecOfLowObsValues)[index];

	UG_ASSERT(GetSize(tmpSol) == GetSize(lowerObsVal), "size of tmpSol and size "
			"of lowerObsVal need to be the same");

	for(size_t j = 0; j < GetSize(tmpSol); j++)
	{
		if ( (BlockRef(tmpSol, j) - BlockRef(lowerObsVal, j)) < 0.0)
		{
			//	u_{s-1/2} < lowerObsValue (:the lower constraint is not fulfilled)

			//	adjust correction c := u_s - u_{s-1} = m_obsVal - u_{s-1}
			BlockRef(c[index], j) = BlockRef(lowerObsVal, j) - BlockRef((*m_lastSol)[index], j);

			//	set new solution u_s to the obstacle value
			//	and store the current index in a vector for further treatment
			BlockRef((*m_lastSol)[index], j) = BlockRef(lowerObsVal, j);
			m_vActiveIndicesLow.push_back(MultiIndex<2>(index, j) );
		}
		else
		{
			//	the 'tmpSol' is valid with respect to the lower constraints
			BlockRef((*m_lastSol)[index], j) = BlockRef(tmpSol, j);
			m_vInactiveIndices.push_back(MultiIndex<2>(index, j) );
		}
	}
}

template <typename TAlgebra>
void
IProjPreconditioner<TAlgebra>::
correction_for_upper_obs(vector_type& c, const size_t index, const value_type& tmpSol)
{
	//	get index-th upper obstacle value
	value_type upperObsVal = (*m_spVecOfUpObsValues)[index];

	UG_ASSERT(GetSize(tmpSol) == GetSize(upperObsVal), "size of tmpSol and size "
			"of upperObsVal need to be the same");

	for(size_t j = 0; j < GetSize(tmpSol); j++)
	{
		if ( (BlockRef(tmpSol, j) - BlockRef(upperObsVal, j)) > 0.0)
		{
			//	u_{s-1/2} > upperObsValue (:the upper constraint is not fulfilled)

			//	adjust correction c := u_s - u_{s-1} = m_obsVal - u_{s-1}
			BlockRef(c[index], j) = BlockRef(upperObsVal, j) - BlockRef((*m_lastSol)[index], j);

			//	set new solution u_s to the obstacle value
			//	and store the current index in a vector for further treatment
			BlockRef((*m_lastSol)[index], j) = BlockRef(upperObsVal, j);
			m_vActiveIndicesUp.push_back(MultiIndex<2>(index, j) );
		}
		else
		{
			//	the 'tmpSol' is valid with respect to the upper constraints
			BlockRef((*m_lastSol)[index], j) = BlockRef(tmpSol, j);
			m_vInactiveIndices.push_back(MultiIndex<2>(index, j) );
		}
	}
}

template <typename TAlgebra>
void
IProjPreconditioner<TAlgebra>::
correction_for_lower_and_upper_obs(vector_type& c, const size_t index, const value_type& tmpSol)
{
	//	get index-th lower obstacle value
	value_type upperObsVal = (*m_spVecOfUpObsValues)[index];
	value_type lowerObsVal = (*m_spVecOfLowObsValues)[index];

	UG_ASSERT(GetSize(tmpSol) == GetSize(upperObsVal), "size of tmpSol and size "
			"of upperObsVal need to be the same");
	UG_ASSERT(GetSize(tmpSol) == GetSize(lowerObsVal), "size of tmpSol and size "
			"of lowerObsVal need to be the same");

	for(size_t j = 0; j < GetSize(tmpSol); j++)
	{
		if ( (BlockRef(tmpSol, j) - BlockRef(upperObsVal, j)) > 0.0)
		{
			//	u_{s-1/2} > upperObsValue (:the upper constraint is not fulfilled)

			//	adjust correction c := u_s - u_{s-1} = m_obsVal - u_{s-1}
			BlockRef(c[index], j)  = BlockRef(upperObsVal, j) - BlockRef((*m_lastSol)[index], j);

			//	set new solution u_s to the obstacle value
			//	and store the current index in a vector for further treatment
			BlockRef((*m_lastSol)[index], j) = BlockRef(upperObsVal, j);
			m_vActiveIndicesUp.push_back(MultiIndex<2>(index, j) );
		}
		else
		{
			if ( (BlockRef(tmpSol, j) - BlockRef(lowerObsVal, j)) < 0.0)
			{
				//	u_{s-1/2} < lowerObsValue (:the lower constraint is not fulfilled)

				//	adjust correction c := u_s - u_{s-1} = m_obsVal - u_{s-1}
				BlockRef(c[index], j) = BlockRef(lowerObsVal, j) - BlockRef((*m_lastSol)[index], j);

				//	set new solution u_s to the obstacle value
				//	and store the current index in a vector for further treatment
				BlockRef((*m_lastSol)[index], j) = BlockRef(lowerObsVal, j);
				m_vActiveIndicesLow.push_back(MultiIndex<2>(index, j) );
			}
			else
			{
				//	the 'tmpSol' is valid with respect to all constraints
				BlockRef((*m_lastSol)[index], j) = BlockRef(tmpSol, j);
				m_vInactiveIndices.push_back(MultiIndex<2>(index, j) );
			}
		}
	}
}

template <typename TAlgebra>
void
IProjPreconditioner<TAlgebra>::
adjust_defect(vector_type& d)
{
	for (std::vector<MultiIndex<2> >::iterator itActiveInd = m_vActiveIndicesLow.begin();
					itActiveInd < m_vActiveIndicesLow.end(); ++itActiveInd)
	{
		//	check, if Ax <= b. For that case the new defect is set to zero,
		//	since all equations/constraints are fulfilled
		if (BlockRef(d[(*itActiveInd)[0]], (*itActiveInd)[1]) < 0.0)
			BlockRef(d[(*itActiveInd)[0]], (*itActiveInd)[1]) = 0.0;
	}

	for (std::vector<MultiIndex<2> >::iterator itActiveInd = m_vActiveIndicesUp.begin();
				itActiveInd < m_vActiveIndicesUp.end(); ++itActiveInd)
	{
		//	check, if Ax >= b. For that case the new defect is set to zero,
		//	since all equations/constraints are fulfilled
		if (BlockRef(d[(*itActiveInd)[0]], (*itActiveInd)[1]) > 0.0)
			BlockRef(d[(*itActiveInd)[0]], (*itActiveInd)[1]) = 0.0;
	}
}

template <typename TAlgebra>
bool
IProjPreconditioner<TAlgebra>::
apply(vector_type &c, const vector_type& d)
{
	PROFILE_FUNC_GROUP("IProjPrecond");
//	Check that operator is initialized
	if(!m_bInit)
	{
		UG_LOG("ERROR in 'IProjPreconditioner::apply': Iterator not initialized.\n");
		return false;
	}

	//	Check parallel status
#ifdef UG_PARALLEL
	if(!d.has_storage_type(PST_ADDITIVE))
		UG_THROW("IProjPreconditioner::apply: Wrong parallel "
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
	if(d.size() != (*m_lastSol).size())
			UG_THROW("Vector [d size= "<<d.size()<<", m_lastSol size = "
					   <<(*m_lastSol).size()<< "] sizes have to match!");
	//	reset vectors
	m_vActiveIndicesLow.resize(0); m_vActiveIndicesUp.resize(0); m_vInactiveIndices.resize(0);

	//	perform a forward GaussSeidel-step: c = (D-L)^-1 * d, then project on the underlying constraint
	//	and store the indices, which satisfy the constraint with equality, in the vector vActiveIndices

#ifdef UG_PARALLEL
	if(pcl::GetNumProcesses() > 1)
	{
	//	make defect unique
		SmartPtr<vector_type> spDtmp = d.clone();
		spDtmp->change_storage_type(PST_UNIQUE);

		projected_precond_step(c, *m_spMat, *spDtmp);
		c.set_storage_type(PST_UNIQUE);
	}
	else
#endif
	{
		projected_precond_step(c, *m_spMat, d);
#ifdef UG_PARALLEL
		c.set_storage_type(PST_UNIQUE);
#endif
	}

//	Correction is always consistent
#ifdef 	UG_PARALLEL
	if(!c.change_storage_type(PST_CONSISTENT))
		UG_THROW("IProjPreconditioner::apply': Cannot change "
				"parallel storage type of correction to consistent.");
#endif

	return true;
}

template <typename TAlgebra>
bool
IProjPreconditioner<TAlgebra>::
apply_update_defect(vector_type &c, vector_type& d)
{
	PROFILE_FUNC_GROUP("IProjPrecond");

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
		UG_LOG("ERROR in 'IProjPreconditioner::apply_update_defect': "
				"Cannot execute matmul_minus to compute d:=d-A*c.\n");
		return false;
	}

	/*for(size_t i = 0; i < d.size(); i++)
	{
		UG_LOG("d_nachUpdate[" << i << "]: " << d[i] << "\n");
	}*/

	//	adjust defect of the active indices for the case that a constraint is set
	if(m_bLowerObs || m_bUpperObs)
		adjust_defect(d);

	/*for(size_t i = 0; i < d.size(); i++)
	{
		UG_LOG("d_nachActiveIndUpdate[" << i << "]: " << d[i] << "\n");
	}*/

	return true;
}

} // end namespace ug

#endif /* __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_PRECONDS__PROJ_PRECOND_INTERFACE_IMPL__ */
