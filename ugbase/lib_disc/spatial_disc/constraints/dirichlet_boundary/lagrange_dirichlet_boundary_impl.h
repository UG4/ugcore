/*
 * lagrange_dirichlet_boundary_impl.h
 *
 *  Created on: 08.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__LAGRANGE_DIRICHLET_BOUNDARY_IMPL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__LAGRANGE_DIRICHLET_BOUNDARY_IMPL__

#include "lagrange_dirichlet_boundary.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/function_spaces/dof_position_util.h"

#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{


////////////////////////////////////////////////////////////////////////////////
//	setup
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
set_approximation_space(SmartPtr<ApproximationSpace<TDomain> > approxSpace)
{
	base_type::set_approximation_space(approxSpace);
	m_spApproxSpace = approxSpace;
	m_spDomain = approxSpace->domain();
	m_aaPos = m_spDomain->position_accessor();
}

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
clear()
{
	m_vBNDNumberData.clear();
	m_vNumberData.clear();
	m_vConstNumberData.clear();
	m_vVectorData.clear();
}

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
add(SmartPtr<UserData<number, dim, bool> > func, const char* function, const char* subsets)
{
	m_vBNDNumberData.push_back(CondNumberData(func, function, subsets));
}

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
add(SmartPtr<UserData<number, dim> > func, const char* function, const char* subsets)
{
	m_vNumberData.push_back(NumberData(func, function, subsets));
}

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
add(number value, const char* function, const char* subsets)
{
	m_vConstNumberData.push_back(ConstNumberData(value, function, subsets));
}

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
add(SmartPtr<UserData<MathVector<dim>, dim> > func, const char* functions, const char* subsets)
{
	m_vVectorData.push_back(VectorData(func, functions, subsets));
}

#ifdef UG_FOR_LUA
template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
add(const char* name, const char* function, const char* subsets)
{
	if(LuaUserData<number, dim>::check_callback_returns(name)){
		SmartPtr<UserData<number, dim> > sp =
							LuaUserDataFactory<number, dim>::create(name);
		add(sp, function, subsets);
		return;
	}
	if(LuaUserData<number, dim, bool>::check_callback_returns(name)){
		SmartPtr<UserData<number, dim, bool> > sp =
						LuaUserDataFactory<number, dim, bool>::create(name);
		add(sp, function, subsets);
		return;
	}
	if(LuaUserData<MathVector<dim>, dim>::check_callback_returns(name)){
		SmartPtr<UserData<MathVector<dim>, dim> > sp =
				LuaUserDataFactory<MathVector<dim>, dim>::create(name);
		add(sp, function, subsets);
		return;
	}

//	no match found
	if(!CheckLuaCallbackName(name))
		UG_THROW("LagrangeDirichlet::add: Lua-Callback with name '"<<name<<
		               "' does not exist.");

//	name exists but wrong signature
	UG_THROW("LagrangeDirichlet::add: Cannot find matching callback "
					"signature. Use one of:\n"
					"a) Number - Callback\n"
					<< (LuaUserData<number, dim>::signature()) << "\n" <<
					"b) Conditional Number - Callback\n"
					<< (LuaUserData<number, dim, bool>::signature()) << "\n" <<
					"c) "<<dim<<"d Vector - Callback\n"
					<< (LuaUserData<MathVector<dim>, dim>::signature()));
}
#endif

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
check_functions_and_subsets(FunctionGroup& functionGroup, SubsetGroup& subsetGroup, size_t numFct) const
{
//	only number of functions allowed
	if(functionGroup.size() != numFct)
		UG_THROW("DirichletBoundary:extract_data:"
					" Only "<<numFct<<" function(s) allowed in specification of a"
					" Dirichlet Value, but the following functions given:"
					<<functionGroup);

//	get subsethandler
	ConstSmartPtr<ISubsetHandler> pSH = m_spApproxSpace->subset_handler();

// 	loop subsets
	for(size_t si = 0; si < subsetGroup.size(); ++si)
	{
	//	get subset index
		const int subsetIndex = subsetGroup[si];

	//	check that subsetIndex is valid
		if(subsetIndex < 0 || subsetIndex >= pSH->num_subsets())
			UG_THROW("DirichletBoundary:extract_data:"
							" Invalid Subset Index " << subsetIndex << ". (Valid is"
							" 0, .. , " << pSH->num_subsets() <<").");

	//	check all functions
		for(size_t i=0; i < functionGroup.size(); ++i)
		{
			const size_t fct = functionGroup[i];

		// 	check if function exist
			if(fct >= m_spApproxSpace->num_fct())
				UG_THROW("DirichletBoundary:extract_data:"
							" Function "<< fct << " does not exist in pattern.");

		// 	check that function is defined for segment
			if(!m_spApproxSpace->is_def_in_subset(fct, subsetIndex))
				UG_THROW("DirichletBoundary:extract_data:"
								" Function "<<fct<<" not defined on subset "<<subsetIndex);
		}
	}

//	everything ok
}

template <typename TDomain, typename TAlgebra>
template <typename TUserData, typename TScheduledUserData>
void DirichletBoundary<TDomain, TAlgebra>::
extract_data(std::map<int, std::vector<TUserData*> >& mvUserDataBndSegment,
             std::vector<TScheduledUserData>& vUserData)
{
//	clear the extracted data
	mvUserDataBndSegment.clear();

	for(size_t i = 0; i < vUserData.size(); ++i)
	{
	//	create Function Group and Subset Group
		try{
			vUserData[i].ssGrp = m_spApproxSpace->subset_grp_by_name(vUserData[i].ssName.c_str());
		}UG_CATCH_THROW(" Subsets '"<<vUserData[i].ssName<<"' not"
		                " all contained in ApproximationSpace.");

		FunctionGroup fctGrp;
		try{
			fctGrp = m_spApproxSpace->fct_grp_by_name(vUserData[i].fctName.c_str());
		}UG_CATCH_THROW(" Functions '"<<vUserData[i].fctName<<"' not"
		                " all contained in ApproximationSpace.");

	//	check functions and subsets
		check_functions_and_subsets(fctGrp, vUserData[i].ssGrp, TUserData::numFct);

	//	set functions
		UG_ASSERT(fctGrp.size() == TUserData::numFct, "wrong number of fct");
		for(size_t fct = 0; fct < TUserData::numFct; ++fct)
		{
			vUserData[i].fct[fct] = fctGrp[fct];
		}

	// 	loop subsets
		for(size_t si = 0; si < vUserData[i].ssGrp.size(); ++si)
		{
		//	get subset index and function
			const int subsetIndex = vUserData[i].ssGrp[si];

		//	remember functor and function
			mvUserDataBndSegment[subsetIndex].push_back(&vUserData[i]);
		}
	}
}


template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
extract_data()
{
//	check that function pattern exists
	if(!m_spApproxSpace.valid())
		UG_THROW("DirichletBoundary:extract_data: "
				" Approximation Space not set.");

	extract_data(m_mNumberBndSegment, m_vNumberData);
	extract_data(m_mBNDNumberBndSegment, m_vBNDNumberData);
	extract_data(m_mConstNumberBndSegment, m_vConstNumberData);
	extract_data(m_mVectorBndSegment, m_vVectorData);
}

////////////////////////////////////////////////////////////////////////////////
//	assemble_dirichlet_rows
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
assemble_dirichlet_rows(matrix_type& mat, ConstSmartPtr<DoFDistribution> dd, number time)
{
	extract_data();

//	loop boundary subsets
	typename std::map<int, std::vector<CondNumberData*> >::const_iterator iter;
	for(iter = m_mBNDNumberBndSegment.begin(); iter != m_mBNDNumberBndSegment.end(); ++iter)
	{
		int si = (*iter).first;
		const std::vector<CondNumberData*>& userData = (*iter).second;

		DoFDistribution::traits<VertexBase>::const_iterator iterBegin 	= dd->begin<VertexBase>(si);
		DoFDistribution::traits<VertexBase>::const_iterator iterEnd 	= dd->end<VertexBase>(si);

	//	create Multiindex
		std::vector<DoFIndex> multInd;

	//	for readin
		MathVector<1> val;
		position_type corner;

	//	loop vertices
		for(DoFDistribution::traits<VertexBase>::const_iterator iter = iterBegin; iter != iterEnd; iter++)
		{
		//	get vertex
			VertexBase* vertex = *iter;

		//	get corner position
			corner = m_aaPos[vertex];

		//	loop dirichlet functions on this segment
			for(size_t i = 0; i < userData.size(); ++i)
			{
			// 	check if function is dirichlet
				if(!(*userData[i])(val, corner, time, si)) continue;

			//	get function index
				const size_t fct = userData[i]->fct[0];

			//	get multi indices
				if(dd->inner_dof_indices(vertex, fct, multInd) != 1)
					return;

				this->m_spAssTuner->set_dirichlet_row(mat, multInd[0]);
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
//	adjust TRANSFER
////////////////////////////////////////////////////////////////////////////////


template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_prolongation(matrix_type& P,
                    ConstSmartPtr<DoFDistribution> ddFine,
                    ConstSmartPtr<DoFDistribution> ddCoarse,
                    number time)
{
	extract_data();

	adjust_prolongation<CondNumberData>(m_mBNDNumberBndSegment, P, ddFine, ddCoarse, time);
	adjust_prolongation<NumberData>(m_mNumberBndSegment, P, ddFine, ddCoarse, time);
	adjust_prolongation<ConstNumberData>(m_mConstNumberBndSegment, P, ddFine, ddCoarse, time);

	adjust_prolongation<VectorData>(m_mVectorBndSegment, P, ddFine, ddCoarse, time);
}

template <typename TDomain, typename TAlgebra>
template <typename TUserData>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_prolongation(const std::map<int, std::vector<TUserData*> >& mvUserData,
                    matrix_type& P,
                    ConstSmartPtr<DoFDistribution> ddFine,
                    ConstSmartPtr<DoFDistribution> ddCoarse,
                    number time)
{
//	loop boundary subsets
	typename std::map<int, std::vector<TUserData*> >::const_iterator iter;
	for(iter = mvUserData.begin(); iter != mvUserData.end(); ++iter)
	{
	//	get subset index
		const int si = (*iter).first;

	//	get vector of scheduled dirichlet data on this subset
		const std::vector<TUserData*>& vUserData = (*iter).second;

	//	adapt jacobian for dofs in each base element type
		try
		{
		if(ddFine->max_dofs(VERTEX)) adjust_prolongation<Vertex, TUserData>(vUserData, si, P, ddFine, ddCoarse, time);
		if(ddFine->max_dofs(EDGE))   adjust_prolongation<EdgeBase, TUserData>(vUserData, si, P, ddFine, ddCoarse, time);
		if(ddFine->max_dofs(FACE))   adjust_prolongation<Face, TUserData>(vUserData, si, P, ddFine, ddCoarse, time);
		if(ddFine->max_dofs(VOLUME)) adjust_prolongation<Volume, TUserData>(vUserData, si, P, ddFine, ddCoarse, time);
		}
		UG_CATCH_THROW("DirichletBoundary::adjust_prolongation:"
						" While calling 'adapt_prolongation' for TUserData, aborting.");
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TBaseElem, typename TUserData>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_prolongation(const std::vector<TUserData*>& vUserData, int si,
                    matrix_type& P,
                    ConstSmartPtr<DoFDistribution> ddFine,
                    ConstSmartPtr<DoFDistribution> ddCoarse,
                    number time)
{
//	create Multiindex
	std::vector<DoFIndex> vFineDoF, vCoarseDoF;

//	dummy for readin
	typename TUserData::value_type val;

//	position of dofs
	std::vector<position_type> vPos;

//	iterators
	typename DoFDistribution::traits<TBaseElem>::const_iterator iter, iterEnd;
	iter = ddFine->begin<TBaseElem>(si);
	iterEnd = ddFine->end<TBaseElem>(si);

//	loop elements
	for( ; iter != iterEnd; iter++)
	{
	//	get vertex
		TBaseElem* elem = *iter;
		GeometricObject* parent = m_spDomain->grid()->get_parent(elem);
		if(!parent) continue;

	//	loop dirichlet functions on this segment
		for(size_t i = 0; i < vUserData.size(); ++i)
		{
			for(size_t f = 0; f < TUserData::numFct; ++f)
			{
			//	get function index
				const size_t fct = vUserData[i]->fct[f];

			//	get local finite element id
				const LFEID& lfeID = ddFine->local_finite_element_id(fct);

			//	get multi indices
				ddFine->inner_dof_indices(elem, fct, vFineDoF);
				ddCoarse->inner_dof_indices(parent, fct, vCoarseDoF);

			//	get dof position
				if(TUserData::isConditional){
					InnerDoFPosition<TDomain>(vPos, elem, *m_spDomain, lfeID);
					UG_ASSERT(vFineDoF.size() == vPos.size(), "Size mismatch");
				}

			//	loop dofs on element
				for(size_t j = 0; j < vFineDoF.size(); ++j)
				{
				// 	check if function is dirichlet
					if(TUserData::isConditional){
						if(!(*vUserData[i])(val, vPos[j], time, si)) continue;
					}

					SetRow(P, vFineDoF[j], 0.0);
				}

				if(vFineDoF.size() > 0){
					for(size_t k = 0; k < vCoarseDoF.size(); ++k){
						DoFRef(P, vFineDoF[0], vCoarseDoF[k]) = 1.0;
					}
				}
			}
		}
	}
}

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_restriction(matrix_type& R,
					ConstSmartPtr<DoFDistribution> ddCoarse,
					ConstSmartPtr<DoFDistribution> ddFine,
					number time)
{
	extract_data();

	adjust_restriction<CondNumberData>(m_mBNDNumberBndSegment, R, ddCoarse, ddFine, time);
	adjust_restriction<NumberData>(m_mNumberBndSegment, R, ddCoarse, ddFine, time);
	adjust_restriction<ConstNumberData>(m_mConstNumberBndSegment, R, ddCoarse, ddFine, time);

	adjust_restriction<VectorData>(m_mVectorBndSegment, R, ddCoarse, ddFine, time);
}

template <typename TDomain, typename TAlgebra>
template <typename TUserData>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_restriction(const std::map<int, std::vector<TUserData*> >& mvUserData,
                   matrix_type& R,
                   ConstSmartPtr<DoFDistribution> ddCoarse,
                   ConstSmartPtr<DoFDistribution> ddFine,
                   number time)
{
//	loop boundary subsets
	typename std::map<int, std::vector<TUserData*> >::const_iterator iter;
	for(iter = mvUserData.begin(); iter != mvUserData.end(); ++iter)
	{
	//	get subset index
		const int si = (*iter).first;

	//	get vector of scheduled dirichlet data on this subset
		const std::vector<TUserData*>& vUserData = (*iter).second;

	//	adapt jacobian for dofs in each base element type
		try
		{
		if(ddFine->max_dofs(VERTEX)) adjust_restriction<Vertex, TUserData>(vUserData, si, R, ddCoarse, ddFine, time);
		if(ddFine->max_dofs(EDGE))   adjust_restriction<EdgeBase, TUserData>(vUserData, si, R, ddCoarse, ddFine, time);
		if(ddFine->max_dofs(FACE))   adjust_restriction<Face, TUserData>(vUserData, si, R, ddCoarse, ddFine, time);
		if(ddFine->max_dofs(VOLUME)) adjust_restriction<Volume, TUserData>(vUserData, si, R, ddCoarse, ddFine, time);
		}
		UG_CATCH_THROW("DirichletBoundary::adjust_restriction:"
						" While calling 'adjust_restriction' for TUserData, aborting.");
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TBaseElem, typename TUserData>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_restriction(const std::vector<TUserData*>& vUserData, int si,
                   matrix_type& R,
                   ConstSmartPtr<DoFDistribution> ddCoarse,
                   ConstSmartPtr<DoFDistribution> ddFine,
                   number time)
{
//	create Multiindex
	std::vector<DoFIndex> vFineDoF, vCoarseDoF;

//	dummy for readin
	typename TUserData::value_type val;

//	position of dofs
	std::vector<position_type> vPos;

//	iterators
	typename DoFDistribution::traits<TBaseElem>::const_iterator iter, iterEnd;
	iter = ddFine->begin<TBaseElem>(si);
	iterEnd = ddFine->end<TBaseElem>(si);

//	loop elements
	for( ; iter != iterEnd; iter++)
	{
	//	get vertex
		TBaseElem* elem = *iter;
		GeometricObject* parent = m_spDomain->grid()->get_parent(elem);
		if(!parent) continue;

	//	loop dirichlet functions on this segment
		for(size_t i = 0; i < vUserData.size(); ++i)
		{
			for(size_t f = 0; f < TUserData::numFct; ++f)
			{
			//	get function index
				const size_t fct = vUserData[i]->fct[f];

			//	get local finite element id
				const LFEID& lfeID = ddFine->local_finite_element_id(fct);

			//	get multi indices
				ddFine->inner_dof_indices(elem, fct, vFineDoF);
				ddCoarse->inner_dof_indices(parent, fct, vCoarseDoF);

			//	get dof position
				if(TUserData::isConditional){
					InnerDoFPosition<TDomain>(vPos, parent, *m_spDomain, lfeID);
					UG_ASSERT(vCoarseDoF.size() == vPos.size(), "Size mismatch");
				}

			//	loop dofs on element
				for(size_t j = 0; j < vCoarseDoF.size(); ++j)
				{
				// 	check if function is dirichlet
					if(TUserData::isConditional){
						if(!(*vUserData[i])(val, vPos[j], time, si)) continue;
					}

					SetRow(R, vCoarseDoF[j], 0.0);
				}

				if(vFineDoF.size() > 0){
					for(size_t k = 0; k < vCoarseDoF.size(); ++k){
						DoFRef(R, vCoarseDoF[k], vFineDoF[0]) = 1.0;
					}
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
//	adjust JACOBIAN
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_jacobian(matrix_type& J, const vector_type& u,
		ConstSmartPtr<DoFDistribution> dd, number time,
		ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		const number s_a0)
{
	extract_data();

	adjust_jacobian<CondNumberData>(m_mBNDNumberBndSegment, J, u, dd, time);
	adjust_jacobian<NumberData>(m_mNumberBndSegment, J, u, dd, time);
	adjust_jacobian<ConstNumberData>(m_mConstNumberBndSegment, J, u, dd, time);

	adjust_jacobian<VectorData>(m_mVectorBndSegment, J, u, dd, time);
}


template <typename TDomain, typename TAlgebra>
template <typename TUserData>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_jacobian(const std::map<int, std::vector<TUserData*> >& mvUserData,
                matrix_type& J, const vector_type& u,
           	    ConstSmartPtr<DoFDistribution> dd, number time)
{
//	loop boundary subsets
	typename std::map<int, std::vector<TUserData*> >::const_iterator iter;
	for(iter = mvUserData.begin(); iter != mvUserData.end(); ++iter)
	{
	//	get subset index
		const int si = (*iter).first;

	//	get vector of scheduled dirichlet data on this subset
		const std::vector<TUserData*>& vUserData = (*iter).second;

	//	adapt jacobian for dofs in each base element type
		try
		{
		if(dd->max_dofs(VERTEX))
			adjust_jacobian<Vertex, TUserData>(vUserData, si, J, u, dd, time);
		if(dd->max_dofs(EDGE))
			adjust_jacobian<EdgeBase, TUserData>(vUserData, si, J, u, dd, time);
		if(dd->max_dofs(FACE))
			adjust_jacobian<Face, TUserData>(vUserData, si, J, u, dd, time);
		if(dd->max_dofs(VOLUME))
			adjust_jacobian<Volume, TUserData>(vUserData, si, J, u, dd, time);
		}
		UG_CATCH_THROW("DirichletBoundary::adjust_jacobian:"
						" While calling 'adapt_jacobian' for TUserData, aborting.");
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TBaseElem, typename TUserData>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_jacobian(const std::vector<TUserData*>& vUserData, int si,
                matrix_type& J, const vector_type& u,
           	    ConstSmartPtr<DoFDistribution> dd, number time)
{
//	create Multiindex
	std::vector<DoFIndex> multInd;

//	dummy for readin
	typename TUserData::value_type val;

//	position of dofs
	std::vector<position_type> vPos;

//	iterators
	typename DoFDistribution::traits<TBaseElem>::const_iterator iter, iterEnd;
	iter = dd->begin<TBaseElem>(si);
	iterEnd = dd->end<TBaseElem>(si);

//	loop elements
	for( ; iter != iterEnd; iter++)
	{
	//	get vertex
		TBaseElem* elem = *iter;

	//	loop dirichlet functions on this segment
		for(size_t i = 0; i < vUserData.size(); ++i)
		{
			for(size_t f = 0; f < TUserData::numFct; ++f)
			{
			//	get function index
				const size_t fct = vUserData[i]->fct[f];

			//	get local finite element id
				const LFEID& lfeID = dd->local_finite_element_id(fct);

			//	get multi indices
				dd->inner_dof_indices(elem, fct, multInd);

			//	get dof position
				if(TUserData::isConditional){
					InnerDoFPosition<TDomain>(vPos, elem, *m_spDomain, lfeID);
					UG_ASSERT(multInd.size() == vPos.size(), "Size mismatch");
				}

			//	loop dofs on element
				for(size_t j = 0; j < multInd.size(); ++j)
				{
				// 	check if function is dirichlet
					if(TUserData::isConditional){
						if(!(*vUserData[i])(val, vPos[j], time, si)) continue;
					}

					this->m_spAssTuner->set_dirichlet_row(J, multInd[j]);
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
//	adjust DEFECT
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_defect(vector_type& d, const vector_type& u,
              ConstSmartPtr<DoFDistribution> dd, number time,
              ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
			  const std::vector<number>* vScaleMass,
			  const std::vector<number>* vScaleStiff)
{
	extract_data();

	adjust_defect<CondNumberData>(m_mBNDNumberBndSegment, d, u, dd, time);
	adjust_defect<NumberData>(m_mNumberBndSegment, d, u, dd, time);
	adjust_defect<ConstNumberData>(m_mConstNumberBndSegment, d, u, dd, time);

	adjust_defect<VectorData>(m_mVectorBndSegment, d, u, dd, time);
}


template <typename TDomain, typename TAlgebra>
template <typename TUserData>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_defect(const std::map<int, std::vector<TUserData*> >& mvUserData,
               vector_type& d, const vector_type& u,
               ConstSmartPtr<DoFDistribution> dd, number time)
{
//	loop boundary subsets
	typename std::map<int, std::vector<TUserData*> >::const_iterator iter;
	for(iter = mvUserData.begin(); iter != mvUserData.end(); ++iter)
	{
	//	get subset index
		const int si = (*iter).first;

	//	get vector of scheduled dirichlet data on this subset
		const std::vector<TUserData*>& vUserData = (*iter).second;

	//	adapt jacobian for dofs in each base element type
		try
		{
		if(dd->max_dofs(VERTEX))
			adjust_defect<Vertex, TUserData>(vUserData, si, d, u, dd, time);
		if(dd->max_dofs(EDGE))
			adjust_defect<EdgeBase, TUserData>(vUserData, si, d, u, dd, time);
		if(dd->max_dofs(FACE))
			adjust_defect<Face, TUserData>(vUserData, si, d, u, dd, time);
		if(dd->max_dofs(VOLUME))
			adjust_defect<Volume, TUserData>(vUserData, si, d, u, dd, time);
		}
		UG_CATCH_THROW("DirichletBoundary::adjust_defect:"
						" While calling 'adjust_defect' for TUserData, aborting.");
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TBaseElem, typename TUserData>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_defect(const std::vector<TUserData*>& vUserData, int si,
              vector_type& d, const vector_type& u,
              ConstSmartPtr<DoFDistribution> dd, number time)
{
//	create Multiindex
	std::vector<DoFIndex> multInd;

//	dummy for readin
	typename TUserData::value_type val;

//	position of dofs
	std::vector<position_type> vPos;

//	iterators
	typename DoFDistribution::traits<TBaseElem>::const_iterator iter, iterEnd;
	iter = dd->begin<TBaseElem>(si);
	iterEnd = dd->end<TBaseElem>(si);

//	loop elements
	for( ; iter != iterEnd; iter++)
	{
	//	get vertex
		TBaseElem* elem = *iter;

	//	loop dirichlet functions on this segment
		for(size_t i = 0; i < vUserData.size(); ++i)
		{
			for(size_t f = 0; f < TUserData::numFct; ++f)
			{
			//	get function index
				const size_t fct = vUserData[i]->fct[f];

			//	get local finite element id
				const LFEID& lfeID = dd->local_finite_element_id(fct);

			//	get multi indices
				dd->inner_dof_indices(elem, fct, multInd);

			//	get dof position
				if(TUserData::isConditional){
					InnerDoFPosition<TDomain>(vPos, elem, *m_spDomain, lfeID);
					UG_ASSERT(multInd.size() == vPos.size(), "Size mismatch. (multInd.size()="<<
					          multInd.size()<<", vPos.size()="<<vPos.size()<<")");
				}

			//	loop dofs on element
				for(size_t j = 0; j < multInd.size(); ++j)
				{
				// 	check if function is dirichlet
					if(TUserData::isConditional){
						if(!(*vUserData[i])(val, vPos[j], time, si)) continue;
					}

					//	set zero for dirichlet values
					this->m_spAssTuner->set_dirichlet_val(d, multInd[j], 0.0);
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
//	adjust SOLUTION
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_solution(vector_type& u, ConstSmartPtr<DoFDistribution> dd, number time)
{
	extract_data();

	adjust_solution<CondNumberData>(m_mBNDNumberBndSegment, u, dd, time);
	adjust_solution<NumberData>(m_mNumberBndSegment, u, dd, time);
	adjust_solution<ConstNumberData>(m_mConstNumberBndSegment, u, dd, time);

	adjust_solution<VectorData>(m_mVectorBndSegment, u, dd, time);
}


template <typename TDomain, typename TAlgebra>
template <typename TUserData>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_solution(const std::map<int, std::vector<TUserData*> >& mvUserData,
                vector_type& u, ConstSmartPtr<DoFDistribution> dd, number time)
{
//	loop boundary subsets
	typename std::map<int, std::vector<TUserData*> >::const_iterator iter;
	for(iter = mvUserData.begin(); iter != mvUserData.end(); ++iter)
	{
	//	get subset index
		const int si = (*iter).first;

	//	get vector of scheduled dirichlet data on this subset
		const std::vector<TUserData*>& vUserData = (*iter).second;

	//	adapt jacobian for dofs in each base element type
		try
		{
		if(dd->max_dofs(VERTEX))
			adjust_solution<Vertex, TUserData>(vUserData, si, u, dd, time);
		if(dd->max_dofs(EDGE))
			adjust_solution<EdgeBase, TUserData>(vUserData, si, u, dd, time);
		if(dd->max_dofs(FACE))
			adjust_solution<Face, TUserData>(vUserData, si, u, dd, time);
		if(dd->max_dofs(VOLUME))
			adjust_solution<Volume, TUserData>(vUserData, si, u, dd, time);
		}
		UG_CATCH_THROW("DirichletBoundary::adjust_solution:"
						" While calling 'adjust_solution' for TUserData, aborting.");
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TBaseElem, typename TUserData>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_solution(const std::vector<TUserData*>& vUserData, int si,
                vector_type& u, ConstSmartPtr<DoFDistribution> dd, number time)
{
//	create Multiindex
	std::vector<DoFIndex> multInd;

//	value readin
	typename TUserData::value_type val;

//	position of dofs
	std::vector<position_type> vPos;

//	iterators
	typename DoFDistribution::traits<TBaseElem>::const_iterator iter, iterEnd;
	iter = dd->begin<TBaseElem>(si);
	iterEnd = dd->end<TBaseElem>(si);

//	loop elements
	for( ; iter != iterEnd; iter++)
	{
	//	get vertex
		TBaseElem* elem = *iter;

	//	loop dirichlet functions on this segment
		for(size_t i = 0; i < vUserData.size(); ++i)
		{
			for(size_t f = 0; f < TUserData::numFct; ++f)
			{
			//	get function index
				const size_t fct = vUserData[i]->fct[f];

			//	get local finite element id
				const LFEID& lfeID = dd->local_finite_element_id(fct);

			//	get dof position
				InnerDoFPosition<TDomain>(vPos, elem, *m_spDomain, lfeID);

			//	get multi indices
				dd->inner_dof_indices(elem, fct, multInd);

				UG_ASSERT(multInd.size() == vPos.size(), "Size mismatch");

			//	loop dofs on element
				for(size_t j = 0; j < multInd.size(); ++j)
				{
				//  get dirichlet value
					if(!(*vUserData[i])(val, vPos[j], time, si)) continue;

					this->m_spAssTuner->set_dirichlet_val(u, multInd[j], val[f]);
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
//	adjust LINEAR
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_linear(matrix_type& A, vector_type& b,
              ConstSmartPtr<DoFDistribution> dd, number time)
{
	extract_data();

	adjust_linear<CondNumberData>(m_mBNDNumberBndSegment, A, b, dd, time);
	adjust_linear<NumberData>(m_mNumberBndSegment, A, b, dd, time);
	adjust_linear<ConstNumberData>(m_mConstNumberBndSegment, A, b, dd, time);

	adjust_linear<VectorData>(m_mVectorBndSegment, A, b, dd, time);
}


template <typename TDomain, typename TAlgebra>
template <typename TUserData>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_linear(const std::map<int, std::vector<TUserData*> >& mvUserData,
              matrix_type& A, vector_type& b,
           	  ConstSmartPtr<DoFDistribution> dd, number time)
{
//	loop boundary subsets
	typename std::map<int, std::vector<TUserData*> >::const_iterator iter;
	for(iter = mvUserData.begin(); iter != mvUserData.end(); ++iter)
	{
	//	get subset index
		const int si = (*iter).first;

	//	get vector of scheduled dirichlet data on this subset
		const std::vector<TUserData*>& vUserData = (*iter).second;

	//	adapt jacobian for dofs in each base element type
		try
		{
		if(dd->max_dofs(VERTEX))
			adjust_linear<Vertex, TUserData>(vUserData, si, A, b, dd, time);
		if(dd->max_dofs(EDGE))
			adjust_linear<EdgeBase, TUserData>(vUserData, si, A, b, dd, time);
		if(dd->max_dofs(FACE))
			adjust_linear<Face, TUserData>(vUserData, si, A, b, dd, time);
		if(dd->max_dofs(VOLUME))
			adjust_linear<Volume, TUserData>(vUserData, si, A, b, dd, time);
		}
		UG_CATCH_THROW("DirichletBoundary::adjust_linear:"
						" While calling 'adjust_linear' for TUserData, aborting.");
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TBaseElem, typename TUserData>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_linear(const std::vector<TUserData*>& vUserData, int si,
              matrix_type& A, vector_type& b,
              ConstSmartPtr<DoFDistribution> dd, number time)
{
//	create Multiindex
	std::vector<DoFIndex> multInd;

//	readin value
	typename TUserData::value_type val;

//	position of dofs
	std::vector<position_type> vPos;

//	iterators
	typename DoFDistribution::traits<TBaseElem>::const_iterator iter, iterEnd;
	iter = dd->begin<TBaseElem>(si);
	iterEnd = dd->end<TBaseElem>(si);

//	loop elements
	for( ; iter != iterEnd; iter++)
	{
	//	get vertex
		TBaseElem* elem = *iter;

	//	loop dirichlet functions on this segment
		for(size_t i = 0; i < vUserData.size(); ++i)
		{
			for(size_t f = 0; f < TUserData::numFct; ++f)
			{
			//	get function index
				const size_t fct = vUserData[i]->fct[f];

			//	get local finite element id
				const LFEID& lfeID = dd->local_finite_element_id(fct);

			//	get dof position
				InnerDoFPosition<TDomain>(vPos, elem, *m_spDomain, lfeID);

			//	get multi indices
				dd->inner_dof_indices(elem, fct, multInd);

				UG_ASSERT(multInd.size() == vPos.size(),
						  "Mismatch: numInd="<<multInd.size()<<", numPos="
						  <<vPos.size()<<" on "<<elem->reference_object_id());

			//	loop dofs on element
				for(size_t j = 0; j < multInd.size(); ++j)
				{
				// 	check if function is dirichlet and read value
					if(!(*vUserData[i])(val, vPos[j], time, si)) continue;

					this->m_spAssTuner->set_dirichlet_row(A, multInd[j]);
					this->m_spAssTuner->set_dirichlet_val(b, multInd[j], val[f]);
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
//	adjust RHS
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_rhs(vector_type& b, const vector_type& u,
           ConstSmartPtr<DoFDistribution> dd, number time)
{
	extract_data();

	adjust_rhs<CondNumberData>(m_mBNDNumberBndSegment, b, u, dd, time);
	adjust_rhs<NumberData>(m_mNumberBndSegment, b, u, dd, time);
	adjust_rhs<ConstNumberData>(m_mConstNumberBndSegment, b, u, dd, time);

	adjust_rhs<VectorData>(m_mVectorBndSegment, b, u, dd, time);
}


template <typename TDomain, typename TAlgebra>
template <typename TUserData>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_rhs(const std::map<int, std::vector<TUserData*> >& mvUserData,
           vector_type& b, const vector_type& u,
           ConstSmartPtr<DoFDistribution> dd, number time)
{
//	loop boundary subsets
	typename std::map<int, std::vector<TUserData*> >::const_iterator iter;
	for(iter = mvUserData.begin(); iter != mvUserData.end(); ++iter)
	{
	//	get subset index
		const int si = (*iter).first;

	//	get vector of scheduled dirichlet data on this subset
		const std::vector<TUserData*>& vUserData = (*iter).second;

	//	adapt jacobian for dofs in each base element type
		try
		{
		if(dd->max_dofs(VERTEX))
			adjust_rhs<Vertex, TUserData>(vUserData, si, b, u, dd, time);
		if(dd->max_dofs(EDGE))
			adjust_rhs<EdgeBase, TUserData>(vUserData, si, b, u, dd, time);
		if(dd->max_dofs(FACE))
			adjust_rhs<Face, TUserData>(vUserData, si, b, u, dd, time);
		if(dd->max_dofs(VOLUME))
			adjust_rhs<Volume, TUserData>(vUserData, si, b, u, dd, time);
		}
		UG_CATCH_THROW("DirichletBoundary::adjust_rhs:"
						" While calling 'adjust_rhs' for TUserData, aborting.");
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TBaseElem, typename TUserData>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_rhs(const std::vector<TUserData*>& vUserData, int si,
           vector_type& b, const vector_type& u,
           ConstSmartPtr<DoFDistribution> dd, number time)
{
//	create Multiindex
	std::vector<DoFIndex> multInd;

//	readin value
	typename TUserData::value_type val;

//	position of dofs
	std::vector<position_type> vPos;

//	iterators
	typename DoFDistribution::traits<TBaseElem>::const_iterator iter, iterEnd;
	iter = dd->begin<TBaseElem>(si);
	iterEnd = dd->end<TBaseElem>(si);

//	loop elements
	for( ; iter != iterEnd; iter++)
	{
	//	get vertex
		TBaseElem* elem = *iter;

	//	loop dirichlet functions on this segment
		for(size_t i = 0; i < vUserData.size(); ++i)
		{
			for(size_t f = 0; f < TUserData::numFct; ++f)
			{
			//	get function index
				const size_t fct = vUserData[i]->fct[f];

			//	get local finite element id
				const LFEID& lfeID = dd->local_finite_element_id(fct);

			//	get dof position
				InnerDoFPosition<TDomain>(vPos, elem, *m_spDomain, lfeID);

			//	get multi indices
				dd->inner_dof_indices(elem, fct, multInd);

				UG_ASSERT(multInd.size() == vPos.size(), "Size mismatch");

			//	loop dofs on element
				for(size_t j = 0; j < multInd.size(); ++j)
				{
				// 	check if function is dirichlet and read value
					if(!(*vUserData[i])(val, vPos[j], time, si)) continue;

					this->m_spAssTuner->set_dirichlet_val(b, multInd[j], val[f]);

				}
			}
		}
	}
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__LAGRANGE_DIRICHLET_BOUNDARY_IMPL__ */
