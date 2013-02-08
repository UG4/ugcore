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
			if(fct >= m_spApproxSpace->function_pattern()->num_fct())
				UG_THROW("DirichletBoundary:extract_data:"
							" Function "<< fct << " does not exist in pattern.");

		// 	check that function is defined for segment
			if(!m_spApproxSpace->function_pattern()->is_def_in_subset(fct, subsetIndex))
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
template <typename TDD>
void DirichletBoundary<TDomain, TAlgebra>::
assemble_dirichlet_rows(matrix_type& mat, ConstSmartPtr<TDD> dd, number time)
{
	extract_data();

//	loop boundary subsets
	typename std::map<int, std::vector<CondNumberData*> >::const_iterator iter;
	for(iter = m_mBNDNumberBndSegment.begin(); iter != m_mBNDNumberBndSegment.end(); ++iter)
	{
		int si = (*iter).first;
		const std::vector<CondNumberData*>& userData = (*iter).second;

		typename TDD::template traits<VertexBase>::const_iterator iterBegin 	= dd->template begin<VertexBase>(si);
		typename TDD::template traits<VertexBase>::const_iterator iterEnd 	= dd->template end<VertexBase>(si);

	//	create Multiindex
		std::vector<MultiIndex<2> >  multInd;

	//	for readin
		MathVector<1> val;
		position_type corner;

	//	loop vertices
		for(typename TDD::template traits<VertexBase>::const_iterator iter = iterBegin; iter != iterEnd; iter++)
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
				if(dd->inner_multi_indices(vertex, fct, multInd) != 1)
					return;

				const size_t index = multInd[0][0];
				const size_t alpha = multInd[0][1];

			//	set dirichlet row
				SetDirichletRow(mat, index, alpha);
			}
		}
	}
}


////////////////////////////////////////////////////////////////////////////////
//	adjust JACOBIAN
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
template <typename TDD>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_jacobian(matrix_type& J, const vector_type& u, ConstSmartPtr<TDD> dd, number time)
{
	extract_data();

	adjust_jacobian<CondNumberData, TDD>(m_mBNDNumberBndSegment, J, u, dd, time);
	adjust_jacobian<NumberData, TDD>(m_mNumberBndSegment, J, u, dd, time);
	adjust_jacobian<ConstNumberData, TDD>(m_mConstNumberBndSegment, J, u, dd, time);

	adjust_jacobian<VectorData, TDD>(m_mVectorBndSegment, J, u, dd, time);
}


template <typename TDomain, typename TAlgebra>
template <typename TUserData, typename TDD>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_jacobian(const std::map<int, std::vector<TUserData*> >& mvUserData,
                matrix_type& J, const vector_type& u,
           	    ConstSmartPtr<TDD> dd, number time)
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
		if(dd->has_indices_on(VERTEX))
			adjust_jacobian<Vertex, TUserData, TDD>(vUserData, si, J, u, dd, time);
		if(dd->has_indices_on(EDGE))
			adjust_jacobian<EdgeBase, TUserData, TDD>(vUserData, si, J, u, dd, time);
		if(dd->has_indices_on(FACE))
			adjust_jacobian<Face, TUserData, TDD>(vUserData, si, J, u, dd, time);
		if(dd->has_indices_on(VOLUME))
			adjust_jacobian<Volume, TUserData, TDD>(vUserData, si, J, u, dd, time);
		}
		UG_CATCH_THROW("DirichletBoundary::adjust_jacobian:"
						" While calling 'adapt_jacobian' for TUserData, aborting.");
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TBaseElem, typename TUserData, typename TDD>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_jacobian(const std::vector<TUserData*>& vUserData, int si,
                matrix_type& J, const vector_type& u,
           	    ConstSmartPtr<TDD> dd, number time)
{
//	create Multiindex
	std::vector<MultiIndex<2> >  multInd;

//	dummy for readin
	typename TUserData::value_type val;

//	position of dofs
	std::vector<position_type> vPos;

//	iterators
	typename TDD::template traits<TBaseElem>::const_iterator iter, iterEnd;
	iter = dd->template begin<TBaseElem>(si);
	iterEnd = dd->template end<TBaseElem>(si);

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

			//	dimension of fct
				const int dim = dd->dim(fct);

			//	get multi indices
				dd->inner_multi_indices(elem, fct, multInd);

			//	get dof position
				if(TUserData::isConditional){
					InnerDoFPosition<TDomain>(vPos, elem, *m_spDomain, lfeID, dim);
					UG_ASSERT(multInd.size() == vPos.size(), "Size mismatch");
				}

			//	loop dofs on element
				for(size_t j = 0; j < multInd.size(); ++j)
				{
				// 	check if function is dirichlet
					if(TUserData::isConditional){
						if(!(*vUserData[i])(val, vPos[j], time, si)) continue;
					}

				//	set dirichlet row
					SetDirichletRow(J, multInd[j][0], multInd[j][1]);
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
//	adjust DEFECT
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
template <typename TDD>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_defect(vector_type& d, const vector_type& u,
              ConstSmartPtr<TDD> dd, number time)
{
	extract_data();

	adjust_defect<CondNumberData, TDD>(m_mBNDNumberBndSegment, d, u, dd, time);
	adjust_defect<NumberData, TDD>(m_mNumberBndSegment, d, u, dd, time);
	adjust_defect<ConstNumberData, TDD>(m_mConstNumberBndSegment, d, u, dd, time);

	adjust_defect<VectorData, TDD>(m_mVectorBndSegment, d, u, dd, time);
}


template <typename TDomain, typename TAlgebra>
template <typename TUserData, typename TDD>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_defect(const std::map<int, std::vector<TUserData*> >& mvUserData,
               vector_type& d, const vector_type& u,
               ConstSmartPtr<TDD> dd, number time)
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
		if(dd->has_indices_on(VERTEX))
			adjust_defect<Vertex, TUserData, TDD>(vUserData, si, d, u, dd, time);
		if(dd->has_indices_on(EDGE))
			adjust_defect<EdgeBase, TUserData, TDD>(vUserData, si, d, u, dd, time);
		if(dd->has_indices_on(FACE))
			adjust_defect<Face, TUserData, TDD>(vUserData, si, d, u, dd, time);
		if(dd->has_indices_on(VOLUME))
			adjust_defect<Volume, TUserData, TDD>(vUserData, si, d, u, dd, time);
		}
		UG_CATCH_THROW("DirichletBoundary::adjust_defect:"
						" While calling 'adjust_defect' for TUserData, aborting.");
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TBaseElem, typename TUserData, typename TDD>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_defect(const std::vector<TUserData*>& vUserData, int si,
              vector_type& d, const vector_type& u,
              ConstSmartPtr<TDD> dd, number time)
{
//	create Multiindex
	std::vector<MultiIndex<2> >  multInd;

//	dummy for readin
	typename TUserData::value_type val;

//	position of dofs
	std::vector<position_type> vPos;

//	iterators
	typename TDD::template traits<TBaseElem>::const_iterator iter, iterEnd;
	iter = dd->template begin<TBaseElem>(si);
	iterEnd = dd->template end<TBaseElem>(si);

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
				dd->inner_multi_indices(elem, fct, multInd);

			//	dimension of fct
				const int dim = dd->dim(fct);

			//	get dof position
				if(TUserData::isConditional){
					InnerDoFPosition<TDomain>(vPos, elem, *m_spDomain, lfeID, dim);
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
					BlockRef(d[multInd[j][0]], multInd[j][1]) = 0.0;
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
//	adjust SOLUTION
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
template <typename TDD>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_solution(vector_type& u, ConstSmartPtr<TDD> dd, number time)
{
	extract_data();

	adjust_solution<CondNumberData, TDD>(m_mBNDNumberBndSegment, u, dd, time);
	adjust_solution<NumberData, TDD>(m_mNumberBndSegment, u, dd, time);
	adjust_solution<ConstNumberData, TDD>(m_mConstNumberBndSegment, u, dd, time);

	adjust_solution<VectorData, TDD>(m_mVectorBndSegment, u, dd, time);
}


template <typename TDomain, typename TAlgebra>
template <typename TUserData, typename TDD>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_solution(const std::map<int, std::vector<TUserData*> >& mvUserData,
                vector_type& u, ConstSmartPtr<TDD> dd, number time)
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
		if(dd->has_indices_on(VERTEX))
			adjust_solution<Vertex, TUserData, TDD>(vUserData, si, u, dd, time);
		if(dd->has_indices_on(EDGE))
			adjust_solution<EdgeBase, TUserData, TDD>(vUserData, si, u, dd, time);
		if(dd->has_indices_on(FACE))
			adjust_solution<Face, TUserData, TDD>(vUserData, si, u, dd, time);
		if(dd->has_indices_on(VOLUME))
			adjust_solution<Volume, TUserData, TDD>(vUserData, si, u, dd, time);
		}
		UG_CATCH_THROW("DirichletBoundary::adjust_solution:"
						" While calling 'adjust_solution' for TUserData, aborting.");
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TBaseElem, typename TUserData, typename TDD>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_solution(const std::vector<TUserData*>& vUserData, int si,
                vector_type& u, ConstSmartPtr<TDD> dd, number time)
{
//	create Multiindex
	std::vector<MultiIndex<2> >  multInd;

//	value readin
	typename TUserData::value_type val;

//	position of dofs
	std::vector<position_type> vPos;

//	iterators
	typename TDD::template traits<TBaseElem>::const_iterator iter, iterEnd;
	iter = dd->template begin<TBaseElem>(si);
	iterEnd = dd->template end<TBaseElem>(si);

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

			//	dimension of fct
				const int dim = dd->dim(fct);

			//	get dof position
				InnerDoFPosition<TDomain>(vPos, elem, *m_spDomain, lfeID, dim);

			//	get multi indices
				dd->inner_multi_indices(elem, fct, multInd);

				UG_ASSERT(multInd.size() == vPos.size(), "Size mismatch");

			//	loop dofs on element
				for(size_t j = 0; j < multInd.size(); ++j)
				{
				//  get dirichlet value
					if(!(*vUserData[i])(val, vPos[j], time, si)) continue;

				//	set zero for dirichlet values
					BlockRef(u[multInd[j][0]], multInd[j][1]) = val[f];
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
//	adjust LINEAR
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
template <typename TDD>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_linear(matrix_type& A, vector_type& b,
              ConstSmartPtr<TDD> dd, number time)
{
	extract_data();

	adjust_linear<CondNumberData, TDD>(m_mBNDNumberBndSegment, A, b, dd, time);
	adjust_linear<NumberData, TDD>(m_mNumberBndSegment, A, b, dd, time);
	adjust_linear<ConstNumberData, TDD>(m_mConstNumberBndSegment, A, b, dd, time);

	adjust_linear<VectorData, TDD>(m_mVectorBndSegment, A, b, dd, time);
}


template <typename TDomain, typename TAlgebra>
template <typename TUserData, typename TDD>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_linear(const std::map<int, std::vector<TUserData*> >& mvUserData,
              matrix_type& A, vector_type& b,
           	  ConstSmartPtr<TDD> dd, number time)
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
		if(dd->has_indices_on(VERTEX))
			adjust_linear<Vertex, TUserData, TDD>(vUserData, si, A, b, dd, time);
		if(dd->has_indices_on(EDGE))
			adjust_linear<EdgeBase, TUserData, TDD>(vUserData, si, A, b, dd, time);
		if(dd->has_indices_on(FACE))
			adjust_linear<Face, TUserData, TDD>(vUserData, si, A, b, dd, time);
		if(dd->has_indices_on(VOLUME))
			adjust_linear<Volume, TUserData, TDD>(vUserData, si, A, b, dd, time);
		}
		UG_CATCH_THROW("DirichletBoundary::adjust_linear:"
						" While calling 'adjust_linear' for TUserData, aborting.");
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TBaseElem, typename TUserData, typename TDD>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_linear(const std::vector<TUserData*>& vUserData, int si,
              matrix_type& A, vector_type& b,
              ConstSmartPtr<TDD> dd, number time)
{
//	create Multiindex
	std::vector<MultiIndex<2> >  multInd;

//	readin value
	typename TUserData::value_type val;

//	position of dofs
	std::vector<position_type> vPos;

//	iterators
	typename TDD::template traits<TBaseElem>::const_iterator iter, iterEnd;
	iter = dd->template begin<TBaseElem>(si);
	iterEnd = dd->template end<TBaseElem>(si);

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

			//	dimension of fct
				const int dim = dd->dim(fct);

			//	get dof position
				InnerDoFPosition<TDomain>(vPos, elem, *m_spDomain, lfeID, dim);

			//	get multi indices
				dd->inner_multi_indices(elem, fct, multInd);

				UG_ASSERT(multInd.size() == vPos.size(),
						  "Mismatch: numInd="<<multInd.size()<<", numPos="<<vPos.size());

			//	loop dofs on element
				for(size_t j = 0; j < multInd.size(); ++j)
				{
				// 	check if function is dirichlet and read value
					if(!(*vUserData[i])(val, vPos[j], time, si)) continue;

					const size_t index = multInd[j][0];
					const size_t alpha = multInd[j][1];

				//	set dirichlet row
					SetDirichletRow(A, index, alpha);

				//	set dirichlet value in rhs
					BlockRef(b[index], alpha) = val[f];
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
//	adjust RHS
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
template <typename TDD>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_rhs(vector_type& b, const vector_type& u,
           ConstSmartPtr<TDD> dd, number time)
{
	extract_data();

	adjust_rhs<CondNumberData, TDD>(m_mBNDNumberBndSegment, b, u, dd, time);
	adjust_rhs<NumberData, TDD>(m_mNumberBndSegment, b, u, dd, time);
	adjust_rhs<ConstNumberData, TDD>(m_mConstNumberBndSegment, b, u, dd, time);

	adjust_rhs<VectorData, TDD>(m_mVectorBndSegment, b, u, dd, time);
}


template <typename TDomain, typename TAlgebra>
template <typename TUserData, typename TDD>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_rhs(const std::map<int, std::vector<TUserData*> >& mvUserData,
           vector_type& b, const vector_type& u,
           ConstSmartPtr<TDD> dd, number time)
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
		if(dd->has_indices_on(VERTEX))
			adjust_rhs<Vertex, TUserData, TDD>(vUserData, si, b, u, dd, time);
		if(dd->has_indices_on(EDGE))
			adjust_rhs<EdgeBase, TUserData, TDD>(vUserData, si, b, u, dd, time);
		if(dd->has_indices_on(FACE))
			adjust_rhs<Face, TUserData, TDD>(vUserData, si, b, u, dd, time);
		if(dd->has_indices_on(VOLUME))
			adjust_rhs<Volume, TUserData, TDD>(vUserData, si, b, u, dd, time);
		}
		UG_CATCH_THROW("DirichletBoundary::adjust_rhs:"
						" While calling 'adjust_rhs' for TUserData, aborting.");
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TBaseElem, typename TUserData, typename TDD>
void DirichletBoundary<TDomain, TAlgebra>::
adjust_rhs(const std::vector<TUserData*>& vUserData, int si,
           vector_type& b, const vector_type& u,
           ConstSmartPtr<TDD> dd, number time)
{
//	create Multiindex
	std::vector<MultiIndex<2> >  multInd;

//	readin value
	typename TUserData::value_type val;

//	position of dofs
	std::vector<position_type> vPos;

//	iterators
	typename TDD::template traits<TBaseElem>::const_iterator iter, iterEnd;
	iter = dd->template begin<TBaseElem>(si);
	iterEnd = dd->template end<TBaseElem>(si);

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

			//	dimension of fct
				const int dim = dd->dim(fct);

			//	get dof position
				InnerDoFPosition<TDomain>(vPos, elem, *m_spDomain, lfeID, dim);

			//	get multi indices
				dd->inner_multi_indices(elem, fct, multInd);

				UG_ASSERT(multInd.size() == vPos.size(), "Size mismatch");

			//	loop dofs on element
				for(size_t j = 0; j < multInd.size(); ++j)
				{
				// 	check if function is dirichlet and read value
					if(!(*vUserData[i])(val, vPos[j], time, si)) continue;

					const size_t index = multInd[j][0];
					const size_t alpha = multInd[j][1];

				//	set dirichlet value in rhs
					BlockRef(b[index], alpha) = val[f];
				}
			}
		}
	}
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__LAGRANGE_DIRICHLET_BOUNDARY_IMPL__ */
