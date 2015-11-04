/*
 * proj_gauss_seidel_interface_impl.h
 *
 *  Created on: 13.11.2013
 *      Author: raphaelprohl
 */

#ifndef __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__PROJ_GAUSS_SEIDEL_INTERFACE_IMPL__
#define __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__PROJ_GAUSS_SEIDEL_INTERFACE_IMPL__

#include "proj_gauss_seidel_interface.h"

#define PROFILE_PROJ_GAUSS_SEIDEL
#ifdef PROFILE_PROJ_GAUSS_SEIDEL
	#define PROJ_GAUSS_SEIDEL_PROFILE_FUNC()		PROFILE_FUNC()
	#define PROJ_GAUSS_SEIDEL_PROFILE_BEGIN(name)	PROFILE_BEGIN_GROUP(name, "IProjGaussSeidel")
	#define PROJ_GAUSS_SEIDEL_PROFILE_END()		PROFILE_END()
#else
	#define PROJ_GAUSS_SEIDEL_PROFILE_FUNC()
	#define PROJ_GAUSS_SEIDEL_PROFILE_BEGIN(name)
	#define PROJ_GAUSS_SEIDEL_PROFILE_END()
#endif

namespace ug{

template <typename TDomain, typename TAlgebra>
void
IProjGaussSeidel<TDomain,TAlgebra>::
truncateVec(vector_type& vec, vector<DoFIndex>& vInd)
{
	typedef typename vector<DoFIndex>::iterator iter_type;
	iter_type dofIter = vInd.begin();
	iter_type dofIterEnd = vInd.end();
	for( ; dofIter != dofIterEnd; dofIter++)
	{
		if (vec.size() <= (*dofIter)[0])
			UG_THROW("vec size is to small in IProjGaussSeidel::truncateVec \n");

		UG_LOG("truncateVec: " <<*dofIter<<"\n");
		//vec[(*dofIter)[0]] = 0.0;
		DoFRef(vec, *dofIter) = 0.0;
	}
}

template <typename TDomain, typename TAlgebra>
void
IProjGaussSeidel<TDomain,TAlgebra>::
truncateMat(matrix_type& mat, vector<DoFIndex>& vInd)
{
	typedef typename vector<DoFIndex>::iterator iter_type;
	iter_type dofIter = vInd.begin();
	iter_type dofIterEnd = vInd.end();
	for( ; dofIter != dofIterEnd; dofIter++)
	{
		UG_LOG("activeDof : " <<*dofIter<< "\n");

		//	set row to zero (for dof '(*dofIter)[0]' and its first function)
		SetRow(mat, (*dofIter)[0], 0.0); //(*dofIter)[1]);
		//	set col to zero
		//BlockRef(mat, , (*dofIter)[0]) = 0.0;
		UG_LOG("truncateMat: mat(" <<(*dofIter)[1]<<","<<(*dofIter)[0]<<") \n");
		mat((*dofIter)[1], (*dofIter)[0]) = 0.0;
	}
}

template <typename TDomain, typename TAlgebra>
bool
IProjGaussSeidel<TDomain,TAlgebra>::
init(SmartPtr<ILinearOperator<vector_type> > J, const vector_type& u)
{
	PROFILE_FUNC_GROUP("IProjGaussSeidel");

//	call GaussSeidelBase-init
	base_type::init(J,u);

//	remember solution
	m_spSol = u.clone();

//	remember, that operator has been initialized
	m_bInit = true;

	UG_LOG("In IProjGaussSeidel::init u hat "<<(*m_spSol).size()<<"Eintraege \n");
	UG_LOG("\n");

	typedef typename vector<SmartPtr<IObstacleConstraint<TDomain,TAlgebra> > >::iterator iter_type;
	iter_type iter = m_spvObsConstraint.begin();
	iter_type iterEnd = m_spvObsConstraint.end();
	for( ; iter != iterEnd; iter++)
		(*iter)->preprocess();

//	(ugly) hint, that usual damping (x += damp * c) does not make sense for the projected
//	GaussSeidel-method.
/*	const number kappa = this->damping()->damping(u, u, m_spOperator);
	if(kappa != 1.0){
		UG_THROW("IProjGaussSeidel::set_damp': Ususal damping is not possible "
				"for IProjGaussSeidel! Use 'set_sor_relax' instead!");
	}*/

	return true;
}

template <typename TDomain, typename TAlgebra>
void
IProjGaussSeidel<TDomain,TAlgebra>::
project_correction(value_type& c_i, const size_t i)
{
	if(!m_bObsCons)
		return;

	typedef typename vector<SmartPtr<IObstacleConstraint<TDomain,TAlgebra> > >::iterator iter_type;
	iter_type iterEnd = m_spvObsConstraint.end();

	for(size_t comp = 0; comp < GetSize(c_i); comp++)
	{
		DoFIndex dof = DoFIndex(i, comp);

		//	loop all obstacle constraint, which are set
		//	& perform a projection: check whether the temporary solution u_{s-1/2}
		//	fulfills the underlying constraint(s) or not
		bool dofIsActive = false;
		bool dofIsObsDoF = false;
		//UG_LOG("dof "<<dof<<"\n");
		//	set iterator to the first obstacle constraint
		iter_type iter = m_spvObsConstraint.begin();
		for( ; iter != iterEnd; iter++)
		{
			//	check, if the dof lies in an obstacle subset: if not -> continue!
			if (!((*iter)->is_obs_dof(dof)))
				continue;

			//UG_LOG("IS IN OBS SUBSET \n");
			dofIsObsDoF = true;

			(*iter)->adjust_sol_and_cor((*m_spSol)[i], c_i, dofIsActive, dof);
		}

		if (dofIsObsDoF && (!dofIsActive))
		{
			// 	dof is admissible -> do regular solution update intern of the
			//	Projected GaussSeidel class
			BlockRef((*m_spSol)[i], comp) += BlockRef(c_i, comp);
		}
	}
}

template <typename TDomain, typename TAlgebra>
bool
IProjGaussSeidel<TDomain,TAlgebra>::
apply(vector_type &c, const vector_type& d)
{
	PROFILE_FUNC_GROUP("IProjGaussSeidel");
//	Check that operator is initialized
	if(!m_bInit)
	{
		UG_LOG("ERROR in 'IProjGaussSeidel::apply': Iterator not initialized.\n");
		return false;
	}

	//	loop all obstacle constraints, which are set & reset its active dofs
	if(m_bObsCons)
	{
		typedef typename vector<SmartPtr<IObstacleConstraint<TDomain,TAlgebra> > >::iterator iter_type;
		iter_type iter = m_spvObsConstraint.begin();
		iter_type iterEnd = m_spvObsConstraint.end();

		for( ; iter != iterEnd; iter++)
			(*iter)->reset_active_dofs();
	}

	base_type::apply(c, d);

	// TODO: in case of using the projected GaussSeidel as smoother in a GMG: TRUNCATION
	const GF& def = dynamic_cast<const GF&>(d);
	ConstSmartPtr<GF> spD = def.clone_without_values();
	if(spD.valid())
	{
		//	check if the DofDistribution is a MGDofDistribution
		//if (spD->dof_distribution()->multi_grid().valid())
		int surfaceLev = spD->dof_distribution()->grid_level().level();

		UG_LOG("NumIndices :" <<spD->dof_distribution()->num_indices() << "\n");
		UG_LOG("numLevels: " << spD->approx_space()->num_levels() << "\n");

		if (spD->dof_distribution()->multi_grid().valid())
		{
			size_t topLev = spD->dof_distribution()->multi_grid()->top_level();
			UG_LOG("surfaceLev: " << surfaceLev << "!\n");
			UG_LOG("topLev: " << topLev << "!\n");

			if((size_t)surfaceLev == topLev)
			{
				UG_LOG("topLev gleich surfaceLev!\n");
				//TODO: for all obstacle constraints:
				if(m_bObsCons)
				{
					typedef typename vector<SmartPtr<IObstacleConstraint<TDomain,TAlgebra> > >::iterator iter_type;
					iter_type iter = m_spvObsConstraint.begin();
					iter_type iterEnd = m_spvObsConstraint.end();

					for( ; iter != iterEnd; iter++)
					{
						//	1. get all active indices
						vector<DoFIndex> vActiveDoFs;
						(*iter)->active_dofs(vActiveDoFs);
						UG_LOG("vActiveDoFs.size() : " <<vActiveDoFs.size()<< "\n");

						typedef typename vector<DoFIndex>::iterator dof_iter_type;
						dof_iter_type dofIter = vActiveDoFs.begin();
						dof_iter_type dofIterEnd = vActiveDoFs.end();
						for( ; dofIter != dofIterEnd; dofIter++)
							UG_LOG("activeDof : " <<*dofIter<< "\n");

						/*UG_LOG("\n");
						//	2. truncation call!
						vector<DoFIndex> testVec;
						DoFIndex dofIndex1(0,1); DoFIndex dofIndex2(1,0); DoFIndex dofIndex3(2,1);
						testVec.push_back(dofIndex1);
						testVec.push_back(dofIndex2);
						testVec.push_back(dofIndex3);

						vector_type& d_top = const_cast<vector_type&>(d);
						d_top[0] = 1.0; d_top[1] = 3.0; d_top[2] = 5.0;
						UG_LOG("d_top(0): "<<d_top[0]<< "\n");
						UG_LOG("d_top(1): "<<d_top[1]<< "\n");
						UG_LOG("d_top(2): "<<d_top[2]<< "\n");
						UG_LOG("d_top(3): "<<d_top[3]<< "\n");
						truncateVec(d_top, testVec);
						UG_LOG("d_top(0): "<<d_top[0]<< "\n");
						UG_LOG("d_top(1): "<<d_top[1]<< "\n");
						UG_LOG("d_top(2): "<<d_top[2]<< "\n");
						UG_LOG("d_top(3): "<<d_top[3]<< "\n");

						matrix_type mat_top;
						mat_top.resize_and_clear(4, 4);
						UG_LOG("#rows(mat): "<<mat_top.num_rows()<<"\n");
						for(size_t i=0; i < mat_top.num_rows(); i++)
						{
							size_t num_connections = mat_top.num_connections(i);
							UG_LOG("#connections in "<<i<<"-th row: "<<num_connections<<"\n");
						}
						mat_top(0,0) = 1.0; mat_top(0,1) = 1.5; mat_top(1,1) = 3.0; mat_top(2,2) = 5.0;
						for(size_t i=0; i < mat_top.num_rows(); i++)
							for(typename matrix_type::row_iterator conn = mat_top.begin_row(i); conn != mat_top.end_row(i); ++conn)
							{
								size_t num_connections = mat_top.num_connections(i);
								UG_LOG("#connections in "<<i<<"-th row: "<<num_connections<<"\n");
								UG_LOG("mat: "<<i<<"-th row: "<< conn.value() << "\n");
							}

						truncateMat(mat_top, testVec);

						for(size_t i=0; i < mat_top.num_rows(); i++)
							for(typename matrix_type::row_iterator conn = mat_top.begin_row(i); conn != mat_top.end_row(i); ++conn)
							{
								UG_LOG("mat: "<<i<<"-th row: "<< conn.value() << "\n");
							}
						UG_LOG("\n");*/

					}
				} //end(if(m_bObsCons))


			}
		}
	}

	return true;
}

template <typename TDomain, typename TAlgebra>
bool
IProjGaussSeidel<TDomain,TAlgebra>::
apply_update_defect(vector_type &c, vector_type& d)
{
	PROFILE_FUNC_GROUP("IProjGaussSeidel");

	//	by calling 'apply_update_defect' the projected Gauss Seidel preconditioner cannot be
	//	a smoother within a multigrid method
	//bIsASmoother = false;

	base_type::apply_update_defect(c, d);

	//	adjust defect for the active dofs
	if(m_bObsCons)
	{
		typedef typename vector<SmartPtr<IObstacleConstraint<TDomain,TAlgebra> > >::iterator iter_type;
		iter_type iter = m_spvObsConstraint.begin();
		iter_type iterEnd = m_spvObsConstraint.end();

		for( ; iter != iterEnd; iter++)
			(*iter)->adjust_defect_to_constraint(d);
	}

	return true;
}


} // end namespace ug

#endif /* __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__PROJ_GAUSS_SEIDEL_INTERFACE_IMPL__ */
