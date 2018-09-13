/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG__LIB_DISC__FUNCTION_SPACES__INTEGRATE__
#define __H__UG__LIB_DISC__FUNCTION_SPACES__INTEGRATE__

#include <cmath>

#include <boost/function.hpp>

#include "common/common.h"

#include "lib_grid/tools/subset_group.h"

#include "lib_disc/common/function_group.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/domain_util.h"  // for CollectCornerCoordinates
#include "lib_disc/quadrature/quadrature_provider.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/user_data/user_data.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"

#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{


////////////////////////////////////////////////////////////////////////////////
// Object encapsulating (scalar) GridFunction related data
////////////////////////////////////////////////////////////////////////////////
template <typename TGridFunction>
class ScalarGridFunctionData
{
public:
	typedef typename TGridFunction::domain_type domain_type;

	ScalarGridFunctionData(TGridFunction& gridFct, size_t cmp)
	: m_gridFct(gridFct), m_fct(cmp),
	 m_id(m_gridFct.local_finite_element_id(m_fct))
	{}

	TGridFunction& grid_function()
	{return  m_gridFct;}

	const TGridFunction& grid_function() const
	{return  m_gridFct;}

	size_t fct()
	{return m_fct;}

	const LFEID &id() const
	{return m_id;}

	//! returns true, iff scalar function is defined in subset si
	bool is_def_in_subset(int si) const
	{ return m_gridFct.is_def_in_subset(m_fct, si); }

	///	returns domain (forward)
	SmartPtr<domain_type> domain() {return m_gridFct.domain();}

	///	returns const domain (forward)
	ConstSmartPtr<domain_type> domain() const {return m_gridFct.domain();}

	template <typename TElem>
	size_t dof_indices(TElem* elem, std::vector<DoFIndex>& ind, bool bHang = false, bool bClear = true) const
	{return m_gridFct.dof_indices(elem, m_fct, ind, bHang, bClear);}

private:
	/// grid function
		TGridFunction& m_gridFct;

	///	component of function
		size_t m_fct;

	///	local finite element id
		LFEID m_id;
};

////////////////////////////////////////////////////////////////////////////////
// Integrand Interface
////////////////////////////////////////////////////////////////////////////////


/// Abstract integrand interface
/*! An integrand is a short-living (temporary) object that is generated/used in various integration functions.*/
template <typename TData, int TWorldDim>
class IIntegrand
{
	public:
	///	world dimension
		static const int worldDim = TWorldDim;

	///	data type
		typedef TData data_type;

	///	returns the values of the integrand for a bunch of ips
	/**
	 *
	 * @param vValue[out]		the value of the integrand at the ips
	 * @param vGlobIP[in]		global integration point positions
	 * @param pElem[in] 		the element to integrate
	 * @param vCornerCoords[in]	corner coordinates of the element
	 * @param vLocIP[in]		local integration point positions
	 * @param vJT[in]			jacobian transposed at integration point
	 * @param numIP[in]			number of integration points
	 */
	/// \{
		virtual void values(TData vValue[],
		                    const MathVector<worldDim> vGlobIP[],
		                    GridObject* pElem,
	                        const MathVector<worldDim> vCornerCoords[],
		                    const MathVector<1> vLocIP[],
		                    const MathMatrix<1, worldDim> vJT[],
		                    const size_t numIP) = 0;
		virtual void values(TData vValue[],
		                    const MathVector<worldDim> vGlobIP[],
		                    GridObject* pElem,
	                        const MathVector<worldDim> vCornerCoords[],
		                    const MathVector<2> vLocIP[],
		                    const MathMatrix<2, worldDim> vJT[],
		                    const size_t numIP) = 0;
		virtual void values(TData vValue[],
		                    const MathVector<worldDim> vGlobIP[],
		                    GridObject* pElem,
	                        const MathVector<worldDim> vCornerCoords[],
		                    const MathVector<3> vLocIP[],
		                    const MathMatrix<3, worldDim> vJT[],
		                    const size_t numIP) = 0;
	/// \}

		virtual ~IIntegrand() {}


	///	sets the subset
		virtual void set_subset(int si) {m_si = si;}

	///	returns the subset
		inline int subset() const {return m_si;}

	protected:
	///	subset
		int m_si;
};

/// Abstract integrand interface (using CRTP)
template <typename TData, int TWorldDim, typename TImpl>
class StdIntegrand : public IIntegrand<TData, TWorldDim>
{
	public:
	///	world dimension
		static const int worldDim = TWorldDim;

	///	data type
		typedef TData data_type;

	/// \copydoc IIntegrand::values
	/// \{
		virtual void values(TData vValue[],
		                    const MathVector<worldDim> vGlobIP[],
		                    GridObject* pElem,
	                        const MathVector<worldDim> vCornerCoords[],
		                    const MathVector<1> vLocIP[],
		                    const MathMatrix<1, worldDim> vJT[],
		                    const size_t numIP)
		{
			getImpl().template evaluate<1>(vValue,vGlobIP,pElem,vCornerCoords,vLocIP,vJT,numIP);
		}
		virtual void values(TData vValue[],
		                    const MathVector<worldDim> vGlobIP[],
		                    GridObject* pElem,
	                        const MathVector<worldDim> vCornerCoords[],
		                    const MathVector<2> vLocIP[],
		                    const MathMatrix<2, worldDim> vJT[],
		                    const size_t numIP)
		{
			getImpl().template evaluate<2>(vValue,vGlobIP,pElem,vCornerCoords,vLocIP,vJT,numIP);
		}
		virtual void values(TData vValue[],
		                    const MathVector<worldDim> vGlobIP[],
		                    GridObject* pElem,
	                        const MathVector<worldDim> vCornerCoords[],
		                    const MathVector<3> vLocIP[],
		                    const MathMatrix<3, worldDim> vJT[],
		                    const size_t numIP)
		{
			getImpl().template evaluate<3>(vValue,vGlobIP,pElem,vCornerCoords,vLocIP,vJT,numIP);
		}
	/// \}

	protected:
	///	access to implementation
		TImpl& getImpl() {return static_cast<TImpl&>(*this);}

	///	const access to implementation
		const TImpl& getImpl() const {return static_cast<const TImpl&>(*this);}
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Generic Volume Integration Routine
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/// integrates on the whole domain
/**
 * This function integrates an arbitrary integrand over the whole domain.
 * Note:
 *  - only grid elements of the same dimension as the world dimension of the
 *    domain are integrated. Thus, no manifolds.
 *  - The implementation is using virtual functions. Thus, there is a small
 *    performance drawback compared to hard coding everything, but we gain
 *    flexibility. In addition all virtual calls compute for the whole set of
 *    integration points to avoid many virtual calls, i.e. only one virtual
 *    call for all integration points is needed.
 *
 * \param[in]		iterBegin	iterator to first geometric object to integrate
 * \param[in]		iterBegin	iterator to last geometric object to integrate
 * \param[in]		integrand	Integrand
 * \param[in]		quadOrder	order of quadrature rule
 * \param[in]		quadType
 * \param[in]		paaElemContribs	(optional). If != NULL, the method will store
 *									the contribution of each element in the
 *									associated attachment entry.
 * \returns			value of the integral
 */
template <int WorldDim, int dim, typename TConstIterator>
number Integrate(TConstIterator iterBegin,
                 TConstIterator iterEnd,
                 typename domain_traits<WorldDim>::position_accessor_type& aaPos,
                 IIntegrand<number, WorldDim>& integrand,
                 int quadOrder, std::string quadType,
                 Grid::AttachmentAccessor<
                 	typename domain_traits<dim>::grid_base_object, ANumber>
                 	*paaElemContribs = NULL
                 )
{
	PROFILE_FUNC();

//	reset the result
	number integral = 0;

//	note: this iterator is for the base elements, e.g. Face and not
//			for the special type, e.g. Triangle, Quadrilateral
	TConstIterator iter = iterBegin;

//	this is the base element type (e.g. Face). This is the type when the
//	iterators above are dereferenciated.
	typedef typename domain_traits<dim>::grid_base_object grid_base_object;

//	get quad type
	if(quadType.empty()) quadType = "best";
	QuadType type = GetQuadratureType(quadType);

//	accessing without dereferencing a pointer first is simpler...
	Grid::AttachmentAccessor<grid_base_object, ANumber> aaElemContribs;
	if(paaElemContribs)
		aaElemContribs = *paaElemContribs;

//	We'll reuse containers to avoid reallocations
	std::vector<MathVector<WorldDim> > vCorner;
	std::vector<MathVector<WorldDim> > vGlobIP;
	std::vector<MathMatrix<dim, WorldDim> > vJT;
	std::vector<number> vValue;

// 	iterate over all elements
	for(; iter != iterEnd; ++iter)
	{
	//	get element
		grid_base_object* pElem = *iter;

	//	get reference object id (i.e. Triangle, Quadrilateral, Tetrahedron, ...)
		ReferenceObjectID roid = (ReferenceObjectID) pElem->reference_object_id();

	//	get quadrature Rule for reference object id and order
		try{
		const QuadratureRule<dim>& rQuadRule
					= QuadratureRuleProvider<dim>::get(roid, quadOrder, type);

	//	get reference element mapping by reference object id
		DimReferenceMapping<dim, WorldDim>& mapping
							= ReferenceMappingProvider::get<dim, WorldDim>(roid);

	//	number of integration points
		const size_t numIP = rQuadRule.size();

	//	get all corner coordinates
		CollectCornerCoordinates(vCorner, *pElem, aaPos, true);

	//	update the reference mapping for the corners
		mapping.update(&vCorner[0]);

	//	compute global integration points
		vGlobIP.resize(numIP);
		mapping.local_to_global(&(vGlobIP[0]), rQuadRule.points(), numIP);

	//	compute transformation matrices
		vJT.resize(numIP);
		mapping.jacobian_transposed(&(vJT[0]), rQuadRule.points(), numIP);

	//	compute integrand values at integration points
		vValue.resize(numIP);
		try
		{
			integrand.values(&(vValue[0]), &(vGlobIP[0]),
			                 pElem, &vCorner[0], rQuadRule.points(),
							 &(vJT[0]),
							 numIP);
		}
		UG_CATCH_THROW("Unable to compute values of integrand at integration point.");

	//	reset contribution of this element
		number intValElem = 0;

	//	loop integration points
		for(size_t ip = 0; ip < numIP; ++ip)
		{
		//	get quadrature weight
			const number weightIP = rQuadRule.weight(ip);

		//	get determinate of mapping
			const number det = SqrtGramDeterminant(vJT[ip]);

		//	add contribution of integration point
			intValElem += vValue[ip] * weightIP * det;
		}

	//	add to global sum
		integral += intValElem;
		if(aaElemContribs.valid())
			aaElemContribs[pElem] = intValElem;

		}UG_CATCH_THROW("SumValuesOnElems failed.");
	} // end elem

//	return the summed integral contributions of all elements
	return integral;
}

template <typename TGridFunction, int dim>
number IntegrateSubset(IIntegrand<number, TGridFunction::dim> &spIntegrand,
                       TGridFunction& spGridFct,
                       int si, int quadOrder, std::string quadType)
{
//	integrate elements of subset
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object grid_base_object;
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator const_iterator;

	spIntegrand.set_subset(si);
	
	return Integrate<TGridFunction::dim,dim,const_iterator>
					(spGridFct.template begin<grid_base_object>(si),
	                 spGridFct.template end<grid_base_object>(si),
					 spGridFct.domain()->position_accessor(),
	                 spIntegrand,
	                 quadOrder, quadType);
}


template <typename TGridFunction>
number IntegrateSubsets(IIntegrand<number, TGridFunction::dim> &spIntegrand,
                        TGridFunction& spGridFct,
                        const char* subsets, int quadOrder,
                        std::string quadType = std::string())
{
//	world dimensions
	static const int dim = TGridFunction::dim;

//	read subsets
	SubsetGroup ssGrp(spGridFct.domain()->subset_handler());
	if(subsets != NULL)
	{
		ssGrp.add(TokenizeString(subsets));
		if(!SameDimensionsInAllSubsets(ssGrp))
			UG_THROW("IntegrateSubsets: Subsets '"<<subsets<<"' do not have same dimension."
			         "Can not integrate on subsets of different dimensions.");
	}
	else
	{
	//	add all subsets and remove lower dim subsets afterwards
		ssGrp.add_all();
		RemoveLowerDimSubsets(ssGrp);
	}

//	reset value
	number value = 0;

//	loop subsets
	for(size_t i = 0; i < ssGrp.size(); ++i)
	{
	//	get subset index
		const int si = ssGrp[i];

	//	check dimension
		if(ssGrp.dim(i) > dim)
			UG_THROW("IntegrateSubsets: Dimension of subset is "<<ssGrp.dim(i)<<", but "
			         " World Dimension is "<<dim<<". Cannot integrate this.");

	//	integrate elements of subset
		try{
		switch(ssGrp.dim(i))
		{
			case DIM_SUBSET_EMPTY_GRID: break;
			case 1: value += IntegrateSubset<TGridFunction, 1>(spIntegrand, spGridFct, si, quadOrder, quadType); break;
			case 2: value += IntegrateSubset<TGridFunction, 2>(spIntegrand, spGridFct, si, quadOrder, quadType); break;
			case 3: value += IntegrateSubset<TGridFunction, 3>(spIntegrand, spGridFct, si, quadOrder, quadType); break;
			default: UG_THROW("IntegrateSubsets: Dimension "<<ssGrp.dim(i)<<" not supported. "
			                  " World dimension is "<<dim<<".");
		}
		}
		UG_CATCH_THROW("IntegrateSubsets: Integration failed on subset "<<si);
	}

#ifdef UG_PARALLEL
	// sum over processes
	if(pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		number local = value;
		com.allreduce(&local, &value, 1, PCL_DT_DOUBLE, PCL_RO_SUM);
	}
#endif

//	return the result
	return value;
}


////////////////////////////////////////////////////////////////////////////////
// UserData Integrand
////////////////////////////////////////////////////////////////////////////////

//! For arbitrary UserData \f$\rho\f$, this class defines the integrand \f$\rho(u)\f$.
template <typename TData, typename TGridFunction>
class UserDataIntegrand
	: public StdIntegrand<TData, TGridFunction::dim, UserDataIntegrand<TData, TGridFunction> >
{
	public:
	//	world dimension of grid function
		static const int worldDim = TGridFunction::dim;

	//	data type
		typedef TData data_type;

	private:
	//  data to integrate
		SmartPtr<UserData<TData, worldDim> > m_spData;

	// 	grid function
		TGridFunction* m_spGridFct;

	//	time
		number m_time;

	public:
	/// constructor
		UserDataIntegrand(SmartPtr<UserData<TData, worldDim> > spData,
		                      TGridFunction* spGridFct,
		                      number time)
		: m_spData(spData), m_spGridFct(spGridFct), m_time(time)
		{
			m_spData->set_function_pattern(spGridFct->function_pattern());
		};

	/// constructor
		UserDataIntegrand(SmartPtr<UserData<TData, worldDim> > spData,
							  number time)
		: m_spData(spData), m_spGridFct(NULL), m_time(time)
		{
			UG_COND_THROW(m_spData->requires_grid_fct(),
						"UserDataIntegrand: Missing GridFunction, but data requires grid function.");
		};

	/// \copydoc IIntegrand::values
		template <int elemDim>
		void evaluate(TData vValue[],
		              const MathVector<worldDim> vGlobIP[],
		              GridObject* pElem,
		              const MathVector<worldDim> vCornerCoords[],
		              const MathVector<elemDim> vLocIP[],
		              const MathMatrix<elemDim, worldDim> vJT[],
		              const size_t numIP)
		{
		//	get local solution if needed
			if(m_spData->requires_grid_fct())
			{
			//	create storage
				LocalIndices ind;
				LocalVector u;

			// 	get global indices
				m_spGridFct->indices(pElem, ind);

			// 	adapt local algebra
				u.resize(ind);

			// 	read local values of u
				GetLocalVector(u, *m_spGridFct);

			//	compute data
				try{
					(*m_spData)(vValue, vGlobIP, m_time, this->m_si, pElem,
								vCornerCoords, vLocIP, numIP, &u, &vJT[0]);
				}
				UG_CATCH_THROW("UserDataIntegrand: Cannot evaluate data.");
			}
			else
			{
			//	compute data
				try{
					(*m_spData)(vValue, vGlobIP, m_time, this->m_si, numIP);
				}
				UG_CATCH_THROW("UserDataIntegrand: Cannot evaluate data.");
			}

		};
};




//! For arbitrary UserData \f$f\f$ (of type TData), this class defines the integrand \f$f^2(u)\f$.
template <typename TData, typename TGridFunction>
class UserDataIntegrandSq
	: public StdIntegrand<number, TGridFunction::dim, UserDataIntegrandSq<TData, TGridFunction> >
{
	public:
	//	world dimension of grid function
		static const int worldDim = TGridFunction::dim;

	//	data type
		typedef TData data_type;

	private:
	//  data to integrate
		SmartPtr<UserData<TData, worldDim> > m_spData;

	// 	grid function
		const TGridFunction* m_pGridFct;

	//	time
		number m_time;

	public:
	/// constructor
		UserDataIntegrandSq(SmartPtr<UserData<TData, worldDim> > spData,
		                      const TGridFunction* pGridFct,
		                      number time)
		: m_spData(spData), m_pGridFct(pGridFct), m_time(time)
		{
			m_spData->set_function_pattern(pGridFct->function_pattern());
		};

	/// constructor
		UserDataIntegrandSq(SmartPtr<UserData<TData, worldDim> > spData,
							  number time)
		: m_spData(spData), m_pGridFct(NULL), m_time(time)
		{
			if(m_spData->requires_grid_fct())
				UG_THROW("UserDataIntegrand: Missing GridFunction, but "
						" data requires grid function.")
		};


	/// \copydoc IIntegrand::values
		template <int elemDim>
		void evaluate(number vValue[],
		              const MathVector<worldDim> vGlobIP[],
		              GridObject* pElem,
		              const MathVector<worldDim> vCornerCoords[],
		              const MathVector<elemDim> vLocIP[],
		              const MathMatrix<elemDim, worldDim> vJT[],
		              const size_t numIP)
		{

			std::vector<TData> tmpValues(numIP);

		//	get local solution if needed
			if(m_spData->requires_grid_fct())
			{
				UG_ASSERT(m_pGridFct!=NULL, "Error: Requires valid pointer!")
			//	create storage
				LocalIndices ind;
				LocalVector u;

			// 	get global indices
				m_pGridFct->indices(pElem, ind);

			// 	adapt local algebra
				u.resize(ind);

			// 	read local values of u
				GetLocalVector(u, *m_pGridFct);

			//	compute data
				try{
					(*m_spData)(&tmpValues.front(), vGlobIP, m_time, this->m_si, pElem,
								vCornerCoords, vLocIP, numIP, &u, &vJT[0]);
				}
				UG_CATCH_THROW("UserDataIntegrand: Cannot evaluate data.");
			}
			else
			{
			//	compute data
				try{
					(*m_spData)(&tmpValues.front(), vGlobIP, m_time, this->m_si, numIP);
				}
				UG_CATCH_THROW("UserDataIntegrand: Cannot evaluate data.");
			}

			for (size_t i=0; i<numIP; ++i)
			{
				vValue[i]=inner_prod(tmpValues[i], tmpValues[i]);
			}

		};
	protected:

		number inner_prod(const number &d1, const number &d2)
		{return d1*d2;}

		number inner_prod(const MathVector<worldDim> &d1, const MathVector<worldDim> &d2)
		{return VecDot(d1, d2);}

		template<typename T>
		number inner_prod(const T &d1, const T &d2)
		{ UG_ASSERT(0, "NOT IMPLEMENTED"); return 0.0;}
};




//! For arbitrary UserData $f$ and grid functions u_1 and u_2, this class (should) define the integrand $ (f(u_1)- f(u_2))^2$.
template <typename TData, typename TGridFunction>
class UserDataDistIntegrandSq
		: public StdIntegrand<number, TGridFunction::dim, UserDataDistIntegrandSq<TData, TGridFunction> >
{
	public:
	///	world dimension of grid function
		static const int worldDim = TGridFunction::dim;
		// typedef typename L2Integrand<TGridFunction>::weight_type weight_type;

	protected:
		ScalarGridFunctionData<TGridFunction> m_fineData;
		const int m_fineTopLevel;

		ScalarGridFunctionData<TGridFunction> m_coarseData;
		const int m_coarseTopLevel;

	///	multigrid
		SmartPtr<MultiGrid> m_spMG;



		SmartPtr<UserData<TData, worldDim> > m_spData;
		// ConstSmartPtr<weight_type> m_spWeight;

		double m_time;
	public:

		/// constructor (1st is fine grid function)
		UserDataDistIntegrandSq(SmartPtr<UserData<TData, worldDim> > spData, TGridFunction& fineGridFct, size_t fineCmp,
					TGridFunction& coarseGridFct, size_t coarseCmp)
		: m_fineData(fineGridFct, fineCmp), m_fineTopLevel(fineGridFct.dof_distribution()->grid_level().level()),
		  m_coarseData(coarseGridFct, coarseCmp), m_coarseTopLevel(coarseGridFct.dof_distribution()->grid_level().level()),
		  m_spMG(m_fineData.domain()->grid()), m_spData(spData), m_time(0.0) /*, m_spWeight(make_sp(new ConstUserNumber<TGridFunction::dim>(1.0)))*/
		{
			if(m_fineTopLevel < m_coarseTopLevel)
				UG_THROW("UserDataDistIntegrandSq: fine and top level inverted.");

			if(m_fineData.domain().get() != m_coarseData.domain().get())
				UG_THROW("UserDataDistIntegrandSq: grid functions defined on different domains.");
		};

		virtual ~UserDataDistIntegrandSq() {}

		///	sets subset
		virtual void set_subset(int si)
		{
			if(!m_fineData.is_def_in_subset(si))
				UG_THROW("UserDataDistIntegrandSq: Grid function component"
						<<m_fineData.fct()<<" not defined on subset "<<si);
			if(!m_coarseData.is_def_in_subset(si))
				UG_THROW("UserDataDistIntegrandSq: Grid function component"
						<<m_coarseData.fct()<<" not defined on subset "<<si);
			IIntegrand<number, worldDim>::set_subset(si);
		}

	/// \copydoc IIntegrand::values
		template <int elemDim>
		void evaluate(number vValue[],
		              const MathVector<worldDim> vGlobIP[],
		              GridObject* pFineElem,
		              const MathVector<worldDim> vCornerCoords[],
		              const MathVector<elemDim> vFineLocIP[],
		              const MathMatrix<elemDim, worldDim> vJT[],
		              const size_t numIP)
		{


		// must return 0.0, if m_spData is independent of grid functions
			if(!m_spData->requires_grid_fct()) {
				for (size_t i=0; i<numIP; ++i) { vValue[i]=0.0; }
				return;
			}

		//	get coarse element
			GridObject* pCoarseElem = pFineElem;
			if(m_coarseTopLevel < m_fineTopLevel){
				int parentLevel = m_spMG->get_level(pCoarseElem);
				while(parentLevel > m_coarseTopLevel){
					pCoarseElem = m_spMG->get_parent(pCoarseElem);
					parentLevel = m_spMG->get_level(pCoarseElem);
				}
			}

		//	get corner coordinates
			typedef typename TGridFunction::template dim_traits<elemDim>::grid_base_object TElem;
			std::vector<MathVector<worldDim> > vCornerCoarse;
			CollectCornerCoordinates(vCornerCoarse, *static_cast<TElem*>(pCoarseElem), *m_coarseData.domain());

		//	get Reference Mapping
			const ReferenceObjectID coarseROID = pCoarseElem->reference_object_id();
			DimReferenceMapping<elemDim, worldDim>& map
				= ReferenceMappingProvider::get<elemDim, worldDim>(coarseROID, vCornerCoarse);

			std::vector<MathVector<elemDim> > vCoarseLocIP;
			vCoarseLocIP.resize(numIP);
			for(size_t ip = 0; ip < vCoarseLocIP.size(); ++ip) VecSet(vCoarseLocIP[ip], 0.0);
			map.global_to_local(&vCoarseLocIP[0], vGlobIP, numIP);

		// element weights
		/*	typedef typename weight_type::data_type ipdata_type;
			std::vector<ipdata_type> fineElemWeights(numIP, 1.0);
			UG_ASSERT(m_spWeight.valid(), "L2DistIntegrand::evaluate requires valid weights! ");
			(*m_spWeight)(&fineElemWeights[0], vFineGlobIP, 0.0, IIntegrand<number, worldDim>::subset(), numIP);*/


		//	get trial space
		/*	const LocalShapeFunctionSet<elemDim>& rFineLSFS =
					LocalFiniteElementProvider::get<elemDim>(fineROID, m_fineData.id());
			const LocalShapeFunctionSet<elemDim>& rCoarseLSFS =
					LocalFiniteElementProvider::get<elemDim>(coarseROID, m_coarseData.id());*/

		//	get multiindices of element
			/*std::vector<DoFIndex> vFineMI, vCoarseMI;
			m_fineData.dof_indices(pFineElem, vFineMI);
			m_coarseData.dof_indices(pCoarseElem, vCoarseMI);*/

			TData fineValues[numIP];
			TData coarseValues[numIP];



		//	compute coarse data
			try{
				LocalVector uCoarse;
				LocalIndices indCoarse;

				m_coarseData.grid_function().indices(pCoarseElem, indCoarse);
				uCoarse.resize(indCoarse);
				GetLocalVector(uCoarse, m_coarseData.grid_function());

				(*m_spData)(coarseValues, vGlobIP, m_time, this->m_si, pCoarseElem,
						&vCornerCoarse[0], &vCoarseLocIP[0], numIP, &uCoarse, &vJT[0]);
			} UG_CATCH_THROW("UserDataDistIntegrandSq: Cannot evaluate coarse data.");


		//	compute fine data
			try{
				LocalVector uFine;
				LocalIndices indFine;

				m_fineData.grid_function().indices(pFineElem, indFine);
				uFine.resize(indFine);
				GetLocalVector(uFine, m_fineData.grid_function());

				(*m_spData)(fineValues, vGlobIP, m_time, this->m_si, pFineElem,
						vCornerCoords, vFineLocIP, numIP, &uFine, &vJT[0]);
			} UG_CATCH_THROW("UserDataDistIntegrandSq: Cannot evaluate fine data.");

		//	loop all integration points
			for(size_t ip = 0; ip < numIP; ++ip)
			{
				vValue[ip] = inner_dist2(fineValues[ip], coarseValues[ip]);
			}

		};


	protected:

			number inner_prod(const number &d1, const number &d2)
			{return d1*d2;}

			number inner_prod(const MathVector<worldDim> &d1, const MathVector<worldDim> &d2)
			{return VecDot(d1, d2);}

			template<typename T>
			number inner_prod(const T &d1, const T &d2)
			{ UG_ASSERT(0, "NOT IMPLEMENTED"); return 0.0;}



			number inner_dist2(const number &v1, const number &v2)
			{ return (v1-v2)*(v1-v2); }

			number inner_dist2(const MathVector<worldDim> &v1, const MathVector<worldDim> &v2)
			{ return VecDistanceSq(v1, v2); }

			template<typename T>
			number inner_dist2(const T &d1, const T &d2)
			{ UG_ASSERT(0, "NOT IMPLEMENTED"); return 0.0;}

};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Volume Integrand implementations
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <typename TGridFunction>
number Integral(SmartPtr<UserData<number, TGridFunction::dim> > spData,
                TGridFunction& spGridFct,
                const char* subsets, number time,
                int quadOrder, std::string quadType)
{
	UserDataIntegrand<number, TGridFunction> spIntegrand(spData, &spGridFct, time);
	return IntegrateSubsets(spIntegrand, spGridFct, subsets, quadOrder, quadType);
}

template <typename TGridFunction>
number Integral(SmartPtr<UserData<number, TGridFunction::dim> > spData, SmartPtr<TGridFunction> spGridFct,
                const char* subsets, number time, int quadOrder, std::string quadType)
{ return Integral(spData, *spGridFct, subsets, time, quadOrder, quadType); }

template <typename TGridFunction>
number Integral(SmartPtr<UserData<number, TGridFunction::dim> > spData, SmartPtr<TGridFunction> spGridFct,const char* subsets,number time, int order)
{return Integral(spData, spGridFct, subsets, time, order, "best");}

template <typename TGridFunction>
number Integral(SmartPtr<UserData<number, TGridFunction::dim> > spData, SmartPtr<TGridFunction> spGridFct,const char* subsets,number time)
{return Integral(spData, spGridFct, subsets, time, 1, "best");}

template <typename TGridFunction>
number Integral(SmartPtr<UserData<number, TGridFunction::dim> > spData, SmartPtr<TGridFunction> spGridFct,number time)
{return Integral(spData, spGridFct, NULL, time, 1, "best");}

template <typename TGridFunction>
number Integral(SmartPtr<UserData<number, TGridFunction::dim> > spData, SmartPtr<TGridFunction> spGridFct,const char* subsets)
{return Integral(spData, spGridFct, subsets, 0.0, 1, "best");}

template <typename TGridFunction>
number Integral(SmartPtr<UserData<number, TGridFunction::dim> > spData, SmartPtr<TGridFunction> spGridFct)
{return Integral(spData, spGridFct, NULL, 0.0, 1, "best");}

///////////////
// const data
///////////////

template <typename TGridFunction>
number Integral(number val, SmartPtr<TGridFunction> spGridFct,
                const char* subsets,
                number time, int quadOrder)
{
	static const int dim = TGridFunction::dim;
	SmartPtr<UserData<number, dim> > sp =
			make_sp(new ConstUserNumber<dim>(val));
	return Integral(sp, spGridFct, subsets, time, quadOrder);
}

template <typename TGridFunction>
number Integral(number val, SmartPtr<TGridFunction> spGridFct,const char* subsets,number time)
{return Integral(val, spGridFct, subsets, time, 1);}

template <typename TGridFunction>
number Integral(number val, SmartPtr<TGridFunction> spGridFct,number time)
{return Integral(val, spGridFct, NULL, time, 1);}

template <typename TGridFunction>
number Integral(number val, SmartPtr<TGridFunction> spGridFct,const char* subsets)
{return Integral(val, spGridFct, subsets, 0.0, 1);}

template <typename TGridFunction>
number Integral(number val, SmartPtr<TGridFunction> spGridFct)
{return Integral(val, spGridFct, NULL, 0.0, 1);}

///////////////
// lua data
///////////////

#ifdef UG_FOR_LUA
template <typename TGridFunction>
number Integral(const char* luaFct,
                SmartPtr<TGridFunction> spGridFct,
                const char* subsets, number time,
                int quadOrder, std::string quadType)
{
	static const int dim = TGridFunction::dim;
	SmartPtr<UserData<number, dim> > sp =
			LuaUserDataFactory<number, dim>::create(luaFct);
	return Integral(sp, spGridFct, subsets, time, quadOrder, quadType);
}

template <typename TGridFunction>
number Integral(const char* luaFct, SmartPtr<TGridFunction> spGridFct,const char* subsets,number time, int quadOrder)
{return Integral(luaFct, spGridFct, subsets, time, quadOrder, "best");}

template <typename TGridFunction>
number Integral(const char* luaFct, SmartPtr<TGridFunction> spGridFct,const char* subsets,number time)
{return Integral(luaFct, spGridFct, subsets, time, 1, "best");}

template <typename TGridFunction>
number Integral(const char* luaFct, SmartPtr<TGridFunction> spGridFct,number time)
{return Integral(luaFct, spGridFct, NULL, time, 1, "best");}

template <typename TGridFunction>
number Integral(const char* luaFct, SmartPtr<TGridFunction> spGridFct,const char* subsets)
{return Integral(luaFct, spGridFct, subsets, 0.0, 1, "best");}

template <typename TGridFunction>
number Integral(const char* luaFct, SmartPtr<TGridFunction> spGridFct)
{return Integral(luaFct, spGridFct, NULL, 0.0, 1, "best");}

#endif


////////////////////////////////////////////////////////////////////////////////
// L2 Error Integrand
////////////////////////////////////////////////////////////////////////////////

template <typename TGridFunction>
class L2ErrorIntegrand
	: public StdIntegrand<number, TGridFunction::dim, L2ErrorIntegrand<TGridFunction> >
{
	public:
	///	world dimension of grid function
		static const int worldDim = TGridFunction::dim;

	private:
		ScalarGridFunctionData<TGridFunction> m_scalarData;

	///  exact solution
		SmartPtr<UserData<number, worldDim> > m_spExactSolution;

	///	time
		number m_time;

	public:
	/// constructor
		L2ErrorIntegrand(SmartPtr<UserData<number, worldDim> > spExactSol,
		                 TGridFunction& gridFct, size_t cmp,
		                 number time)
		: m_scalarData(gridFct, cmp),
		  m_spExactSolution(spExactSol), m_time(time)
		{};

		virtual ~L2ErrorIntegrand() {};

	///	sets subset
		virtual void set_subset(int si)
		{
			if(!m_scalarData.is_def_in_subset(si))
				UG_THROW("L2ErrorIntegrand: Grid function component"
						<<m_scalarData.fct()<<" not defined on subset "<<si);
			IIntegrand<number, worldDim>::set_subset(si);
		}

	/// \copydoc IIntegrand::values
		template <int elemDim>
		void evaluate(number vValue[],
		              const MathVector<worldDim> vGlobIP[],
		              GridObject* pElem,
		              const MathVector<worldDim> vCornerCoords[],
		              const MathVector<elemDim> vLocIP[],
		              const MathMatrix<elemDim, worldDim> vJT[],
		              const size_t numIP)
		{
		//	get reference object id (i.e. Triangle, Quadrilateral, Tetrahedron, ...)
			const ReferenceObjectID roid = pElem->reference_object_id();

			try{
		//	get trial space
			const LocalShapeFunctionSet<elemDim>& rTrialSpace =
							LocalFiniteElementProvider::get<elemDim>(roid, m_scalarData.id());

		//	number of dofs on element
			const size_t num_sh = rTrialSpace.num_sh();

		//	get multiindices of element
			std::vector<DoFIndex> ind;  // 	aux. index array
			m_scalarData.dof_indices(pElem, ind);

		//	check multi indices
			if(ind.size() != num_sh)
				UG_THROW("L2ErrorIntegrand::evaluate: Wrong number of"
						" multi indices.");

		//	loop all integration points
			for(size_t ip = 0; ip < numIP; ++ip)
			{
			//	compute exact solution at integration point
				number exactSolIP;
				(*m_spExactSolution)(exactSolIP, vGlobIP[ip], m_time, this->subset());

			// 	compute approximated solution at integration point
				number approxSolIP = 0.0;
				for(size_t sh = 0; sh < num_sh; ++sh)
				{
				//	get value at shape point (e.g. corner for P1 fct)
					const number valSH = DoFRef(m_scalarData.grid_function(), ind[sh]);

				//	add shape fct at ip * value at shape
					approxSolIP += valSH * rTrialSpace.shape(sh, vLocIP[ip]);
				}

			//	get squared of difference
				vValue[ip] = (exactSolIP - approxSolIP);
				vValue[ip] *= vValue[ip];
			}

			}
			UG_CATCH_THROW("L2ErrorIntegrand::evaluate: trial space missing.");
		};
};

/// computes the l2 error function on the whole domain or on some subsets
/**
 * This function computes the L2-difference between a given analytic function
 * and a grid function.
 *
 * \param[in]		ExactSol	analytic function
 * \param[in]		spGridFct	grid function
 * \param[in]		cmp			symbolic name of component function
 * \param[in]		time		time point
 * \param[in]		quadOrder	order of quadrature rule
 * \param[in]		subsets		subsets, where to compute
 * \returns			number 		l2-norm of difference
 */
template <typename TGridFunction>
number L2Error(SmartPtr<UserData<number, TGridFunction::dim> > spExactSol,
               TGridFunction& gridFct, const char* cmp,
               number time, int quadOrder, const char* subsets)
{
//	get function id of name
	const size_t fct = gridFct.fct_id_by_name(cmp);

//	check that function exists
	UG_COND_THROW(fct >= gridFct.num_fct(),
				"L2Error: Function space does not contain a function with name " << cmp << ".");

	L2ErrorIntegrand<TGridFunction> spIntegrand(spExactSol, gridFct, fct, time);
	return sqrt(IntegrateSubsets(spIntegrand, gridFct, subsets, quadOrder));
}

template <typename TGridFunction>
number L2Error(SmartPtr<UserData<number, TGridFunction::dim> > spExactSol,
               SmartPtr<TGridFunction> spGridFct, const char* cmp,
               number time, int quadOrder, const char* subsets)
{ return L2Error(spExactSol, *spGridFct, cmp, time, quadOrder, subsets); }

template <typename TGridFunction>
number L2Error(SmartPtr<UserData<number, TGridFunction::dim> > spExactSol,
               SmartPtr<TGridFunction> spGridFct, const char* cmp,
               number time, int quadOrder)
{ return L2Error(spExactSol, *spGridFct, cmp, time, quadOrder, NULL); }




#ifdef UG_FOR_LUA
template <typename TGridFunction>
number L2Error(const char* ExactSol,
               SmartPtr<TGridFunction> spGridFct, const char* cmp,
               number time, int quadOrder, const char* subsets)
{
	SmartPtr<UserData<number, TGridFunction::dim> > spExactSol
	 = make_sp(new LuaUserData<number, TGridFunction::domain_type::dim>(ExactSol));
	return L2Error(spExactSol, spGridFct, cmp, time, quadOrder, subsets);
}

template <typename TGridFunction>
number L2Error(const char* ExactSol,
               SmartPtr<TGridFunction> spGridFct, const char* cmp,
               number time, int quadOrder)
{
	return L2Error(ExactSol, spGridFct, cmp, time, quadOrder, NULL);
}
#endif




////////////////////////////////////////////////////////////////////////////////
// H1 Error Integrand
////////////////////////////////////////////////////////////////////////////////

template <typename TGridFunction>
class H1ErrorIntegrand
	: public StdIntegrand<number, TGridFunction::dim, H1ErrorIntegrand<TGridFunction> >
{
	public:
	///	world dimension of grid function
		static const int worldDim = TGridFunction::dim;

	private:
	/// grid function
		ScalarGridFunctionData<TGridFunction> m_scalarData;

	///	exact solution
		SmartPtr<UserData<number, worldDim> > m_spExactSolution;

	///	exact gradient
		SmartPtr<UserData<MathVector<worldDim>, worldDim> > m_spExactGrad;

	///	time
		number m_time;

	public:
	/// constructor
		H1ErrorIntegrand(SmartPtr<UserData<number, worldDim> > spExactSol,
		                 SmartPtr<UserData<MathVector<worldDim>, worldDim> > spExactGrad,
		                 TGridFunction& gridFct, size_t cmp,
		                 number time)
		: m_scalarData (gridFct, cmp),
		  m_spExactSolution(spExactSol),
		  m_spExactGrad(spExactGrad),
		  m_time(time)
		{}

	///	sets subset
		virtual void set_subset(int si)
		{
			if(!m_scalarData.is_def_in_subset(si))
				UG_THROW("H1Error: Grid function component"
						<<m_scalarData.fct()<<" not defined on subset "<<si);
			IIntegrand<number, worldDim>::set_subset(si);
		}

	/// \copydoc IIntegrand::values
		template <int elemDim>
		void evaluate(number vValue[],
		              const MathVector<worldDim> vGlobIP[],
		              GridObject* pElem,
		              const MathVector<worldDim> vCornerCoords[],
		              const MathVector<elemDim> vLocIP[],
		              const MathMatrix<elemDim, worldDim> vJT[],
		              const size_t numIP)
		{
		//	get reference object id (i.e. Triangle, Quadrilateral, Tetrahedron, ...)
			const ReferenceObjectID roid = pElem->reference_object_id();

			try{
		//	get trial space
			const LocalShapeFunctionSet<elemDim>& rTrialSpace =
							LocalFiniteElementProvider::get<elemDim>(roid, m_scalarData.id());

		//	number of dofs on element
			const size_t num_sh = rTrialSpace.num_sh();

		//	get multiindices of element
			std::vector<DoFIndex> ind;  // 	aux. index array
			m_scalarData.dof_indices(pElem, ind);

		//	check multi indices
			if(ind.size() != num_sh)
				UG_THROW("H1ErrorIntegrand::evaluate: Wrong number of"
						" multi indices.");

		//	loop all integration points
			std::vector<MathVector<elemDim> > vLocGradient(num_sh);
			for(size_t ip = 0; ip < numIP; ++ip)
			{
			//	compute exact solution at integration point
				number exactSolIP;
				(*m_spExactSolution)(exactSolIP, vGlobIP[ip], m_time, this->subset());

			//	compute exact gradient at integration point
				MathVector<worldDim> exactGradIP;
				(*m_spExactGrad)(exactGradIP, vGlobIP[ip], m_time, this->subset());

			//	compute shape gradients at ip
				rTrialSpace.grads(&vLocGradient[0], vLocIP[ip]);

			// 	compute approximated solution at integration point
				number approxSolIP = 0.0;
				MathVector<elemDim> locTmp; VecSet(locTmp, 0.0);
				for(size_t sh = 0; sh < num_sh; ++sh)
				{
				//	get value at shape point (e.g. corner for P1 fct)
					const number valSH = DoFRef(m_scalarData.grid_function(), ind[sh]);

				//	add shape fct at ip * value at shape
					approxSolIP += valSH * rTrialSpace.shape(sh, vLocIP[ip]);

				//	add gradient at ip
					VecScaleAppend(locTmp, valSH, vLocGradient[sh]);
				}

			//	compute global gradient
				MathVector<worldDim> approxGradIP;
				MathMatrix<worldDim, elemDim> JTInv;
				Inverse(JTInv, vJT[ip]);
				MatVecMult(approxGradIP, JTInv, locTmp);

			//	get squared of difference
				vValue[ip] = (exactSolIP - approxSolIP) * (exactSolIP - approxSolIP);
				vValue[ip] += VecDistanceSq(approxGradIP, exactGradIP);
			}

			}
			UG_CATCH_THROW("H1ErrorIntegrand::evaluate: trial space missing.");
		};
};

/// compute H1 error of a function on the whole domain or on some subsets
/**
 * This function computes the H1-difference between a given analytic function
 * and a grid function.
 *
 * \param[in]		spExactSol	analytic function
 * \param[in]		spExactGrad	analytic gradient
 * \param[in]		spGridFct	grid function
 * \param[in]		cmp			symbolic name of component function
 * \param[in]		time		time point
 * \param[in]		quadOrder	order of quadrature rule
 * \param[in]		subsets		subsets, where to compute
 * \returns			number 		l2-norm of difference
 */
template <typename TGridFunction>
number H1Error(SmartPtr<UserData<number, TGridFunction::dim> > spExactSol,
               SmartPtr<UserData<MathVector<TGridFunction::dim>, TGridFunction::dim> > spExactGrad,
			   SmartPtr<TGridFunction> spGridFct, const char* cmp,
			   number time, int quadOrder, const char* subsets)
{
//	get function id of name
	const size_t fct = spGridFct->fct_id_by_name(cmp);

//	check that function exists
	if(fct >= spGridFct->num_fct())
		UG_THROW("H1Error: Function space does not contain"
				" a function with name " << cmp << ".");

	H1ErrorIntegrand<TGridFunction> spIntegrand(spExactSol, spExactGrad, *spGridFct, fct, time);
	return sqrt(IntegrateSubsets(spIntegrand, *spGridFct, subsets, quadOrder));
}

template <typename TGridFunction>
number H1Error(SmartPtr<UserData<number, TGridFunction::dim> > spExactSol,
               SmartPtr<UserData<MathVector<TGridFunction::dim>, TGridFunction::dim> > spExactGrad,
			   SmartPtr<TGridFunction> spGridFct, const char* cmp,
			   number time, int quadOrder)
{
	return H1Error(spExactSol, spExactGrad, spGridFct, cmp, time, quadOrder, NULL);
}

#ifdef UG_FOR_LUA
template <typename TGridFunction>
number H1Error(const char* ExactSol, const char* ExactGrad,
			   SmartPtr<TGridFunction> spGridFct, const char* cmp,
			   number time, int quadOrder, const char* subsets)
{
	static const int dim = TGridFunction::domain_type::dim;
	SmartPtr<UserData<number, dim> > spExactSol
	 = make_sp(new LuaUserData<number, dim>(ExactSol));
	SmartPtr<UserData<MathVector<dim>, dim> > spExactGrad
	 = make_sp(new LuaUserData<MathVector<dim>, dim>(ExactGrad));
	return H1Error(spExactSol, spExactGrad, spGridFct, cmp, time, quadOrder, subsets);
}

template <typename TGridFunction>
number H1Error(const char* ExactSol, const char* ExactGrad,
			   SmartPtr<TGridFunction> spGridFct, const char* cmp,
			   number time, int quadOrder)
{
	return H1Error(ExactSol, ExactGrad, spGridFct, cmp, time, quadOrder, NULL);
}
#endif



/// computes an (abstract) distance between two functions
/**
 * This function computes the (abstract) TDiffError-difference between two grid functions that
 * may be defined on different grids. The element loop is performed over the
 * finer level.
 *
 * \param[in]		spGridFct1	grid function 1
 * \param[in]		cmp1		symbolic name of component function
 * \param[in]		spGridFct2	grid function 2
 * \param[in]		cmp2		symbolic name of component function
 * \param[in]		quadOrder	order of quadrature rule
 * \param[in]		subsets		subsets, where to compute
 * \returns			number 		H1-norm of difference
 */
template <typename TDistIntegrand, typename TGridFunction>
number GridFunctionDistance2(TGridFunction& spGridFct1, const char* cmp1,
							TGridFunction& spGridFct2, const char* cmp2,
							int quadOrder, const char* subsets)
{
//	get function id of name
	const size_t fct1 = spGridFct1.fct_id_by_name(cmp1);
	const size_t fct2 = spGridFct2.fct_id_by_name(cmp2);

//	check that function exists
	if(fct1 >= spGridFct1.num_fct())
		UG_THROW("GridFunctionDistance: Function space does not contain"
				" a function with name " << cmp1 << ".");
	if(fct2 >= spGridFct2.num_fct())
		UG_THROW("GridFunctionDistance: Function space does not contain"
				" a function with name " << cmp2 << ".");

//	get top level of grid functions
	const int level1 = spGridFct1.dof_distribution()->grid_level().level();
	const int level2 = spGridFct2.dof_distribution()->grid_level().level();

	if(level1 > level2){
		// level check
		TDistIntegrand spIntegrand(spGridFct1, fct1, spGridFct2, fct2);
		return IntegrateSubsets(spIntegrand, spGridFct1, subsets, quadOrder);
	}else{
		TDistIntegrand spIntegrand(spGridFct2, fct2, spGridFct1, fct1);
		return IntegrateSubsets(spIntegrand, spGridFct2, subsets, quadOrder);
	}

}

//! Computes (weighted) distance
template <typename TDistIntegrand, typename TGridFunction>
number GridFunctionDistance2(TGridFunction& spGridFct1, const char* cmp1,
               TGridFunction& spGridFct2, const char* cmp2,
               int quadOrder, const char* subsets, ConstSmartPtr<typename TDistIntegrand::weight_type> spWeights)
{
//	get function id of name
	const size_t fct1 = spGridFct1.fct_id_by_name(cmp1);
	const size_t fct2 = spGridFct2.fct_id_by_name(cmp2);

//	check that function exists
	if(fct1 >= spGridFct1.num_fct())
		UG_THROW("GridFunctionDistance: Function space does not contain"
				" a function with name " << cmp1 << ".");
	if(fct2 >= spGridFct2.num_fct())
		UG_THROW("GridFunctionDistance: Function space does not contain"
				" a function with name " << cmp2 << ".");

//	get top level of gridfunctions
	const int level1 = spGridFct1.dof_distribution()->grid_level().level();
	const int level2 = spGridFct2.dof_distribution()->grid_level().level();

	// w/ weights
	if(level1 > level2){
		TDistIntegrand spIntegrand(spGridFct1, fct1, spGridFct2, fct2, spWeights);
		return IntegrateSubsets(spIntegrand, spGridFct1, subsets, quadOrder);
	}else{
		TDistIntegrand spIntegrand(spGridFct2, fct2, spGridFct1, fct1, spWeights);
		return IntegrateSubsets(spIntegrand, spGridFct2, subsets, quadOrder);
	}
}


//! Computes (weighted) distance with shift for averages
template <typename TDistIntegrand, typename TGridFunction>
number GridFunctionDistance2(TGridFunction& spGridFct1, const char* cmp1,
               TGridFunction& spGridFct2, const char* cmp2,
               int quadOrder, const char* subsets,
			   ConstSmartPtr<typename TDistIntegrand::weight_type> spWeights,
			   number distAvg12)
{
//	get function id of name
	const size_t fct1 = spGridFct1.fct_id_by_name(cmp1);
	const size_t fct2 = spGridFct2.fct_id_by_name(cmp2);

//	check that function exists
	if(fct1 >= spGridFct1.num_fct())
		UG_THROW("GridFunctionDistance: Function space does not contain"
				" a function with name " << cmp1 << ".");
	if(fct2 >= spGridFct2.num_fct())
		UG_THROW("GridFunctionDistance: Function space does not contain"
				" a function with name " << cmp2 << ".");

//	get top level of gridfunctions
	const int level1 = spGridFct1.dof_distribution()->grid_level().level();
	const int level2 = spGridFct2.dof_distribution()->grid_level().level();

	// w/ weights
	if(level1 > level2){
		TDistIntegrand spIntegrand(spGridFct1, fct1, spGridFct2, fct2, spWeights, distAvg12);
		return IntegrateSubsets(spIntegrand, spGridFct1, subsets, quadOrder);
	}else{
		TDistIntegrand spIntegrand(spGridFct2, fct2, spGridFct1, fct1, spWeights, -distAvg12);
		return IntegrateSubsets(spIntegrand, spGridFct2, subsets, quadOrder);
	}
}



////////////////////////////////////////////////////////////////////////////////
// L2 Integrand
////////////////////////////////////////////////////////////////////////////////

/// Grid function as L2 integrand
template <typename TGridFunction>
class L2Integrand
	: public StdIntegrand<number, TGridFunction::dim, L2Integrand<TGridFunction> >
{
	public:
	//	world dimension of grid function
		static const int worldDim = TGridFunction::dim;
		typedef UserData<number, worldDim> weight_type;

	protected:
	// grid function data
		ScalarGridFunctionData<TGridFunction> m_scalarData;

	/// scalar weight (optional, default is 1.0)
		ConstSmartPtr<weight_type> m_spWeight;

	public:
	/// CTOR
		L2Integrand(TGridFunction& spGridFct, size_t cmp)
		: m_scalarData(spGridFct, cmp), m_spWeight(make_sp(new ConstUserNumber<worldDim>(1.0)))
		{};

		L2Integrand(TGridFunction& spGridFct, size_t cmp, ConstSmartPtr<weight_type> spWeight)
		: m_scalarData(spGridFct, cmp), m_spWeight(spWeight)
		{};

	/// DTOR
		virtual ~L2Integrand() {};

	///	sets subset
		virtual void set_subset(int si)
		{
			if(!m_scalarData.is_def_in_subset(si))
				UG_THROW("L2ErrorIntegrand: Grid function component" <<m_scalarData.fct() <<" not defined on subset "<<si);
			IIntegrand<number, worldDim>::set_subset(si);
		}

	/// \copydoc IIntegrand::values
		template <int elemDim>
		void evaluate(number vValue[],
		              const MathVector<worldDim> vGlobIP[],
		              GridObject* pElem,
		              const MathVector<worldDim> vCornerCoords[],
		              const MathVector<elemDim> vLocIP[],
		              const MathMatrix<elemDim, worldDim> vJT[],
		              const size_t numIP)
		{
		//	get reference object id (i.e. Triangle, Quadrilateral, Tetrahedron, ...)
			ReferenceObjectID roid = (ReferenceObjectID) pElem->reference_object_id();

		// element weights
			typedef typename weight_type::data_type ipdata_type;
			std::vector<ipdata_type> locElemWeights(numIP, 1.0);
			UG_ASSERT(m_spWeight.valid(), "L2Integrand::evaluate requires valid weights!");
			(*m_spWeight)(&locElemWeights[0], vGlobIP, 0.0, IIntegrand<number, worldDim>::subset(), numIP);

			try{
		//	get trial space
			const LocalShapeFunctionSet<elemDim>& rTrialSpace =
							LocalFiniteElementProvider::get<elemDim>(roid, m_scalarData.id());

		//	number of dofs on element
			const size_t num_sh = rTrialSpace.num_sh();

		//	get multiindices of element
			std::vector<DoFIndex> ind;  // 	aux. index array
			m_scalarData.dof_indices(pElem, ind);

		//	check multi indices
			if(ind.size() != num_sh)
				UG_THROW("L2Integrand::evaluate: Wrong number of multi indices.");

		//	loop all integration points
			for(size_t ip = 0; ip < numIP; ++ip)
			{

			// 	compute approximated solution at integration point
				number approxSolIP = 0.0;
				for(size_t sh = 0; sh < num_sh; ++sh)
				{
					//	get value at shape point (e.g. corner for P1 fct)
					//	and add shape fct at ip * value at shape
					const number valSH = DoFRef(m_scalarData.grid_function(), ind[sh]);
					approxSolIP += valSH * rTrialSpace.shape(sh, vLocIP[ip]);
				}

				//	get square
				vValue[ip] = locElemWeights[ip]*approxSolIP*approxSolIP;

			}

			}
			UG_CATCH_THROW("L2FuncIntegrand::values: trial space missing.");
		};
};

/**
 * This function computes the square of the L2-norm of a grid function.
 *
 * \param[in]		spGridFct	grid function
 * \param[in]		cmp			symbolic name of function
 * \param[in]		quadOrder	order of quadrature rule
 * \param[in]		subsets		subsets, where to interpolate
 * 								(NULL indicates that all full-dimensional subsets
 * 								shall be considered)
 * \returns			number 		l2-norm
 */

template <typename TGridFunction>
number L2Norm2(TGridFunction& u, const char* cmp,
              int quadOrder, const char* subsets,
			  ConstSmartPtr<typename L2Integrand<TGridFunction>::weight_type> spWeight)
{
//	get function id of name
	const size_t fct = u.fct_id_by_name(cmp);

//	check that function exists
	if(fct >= u.num_fct())
		UG_THROW("L2Norm: Function space does not contain"
				" a function with name " << cmp << ".");

	L2Integrand<TGridFunction> integrandL2(u, fct, spWeight);
	return IntegrateSubsets(integrandL2, u, subsets, quadOrder);
}

template <typename TGridFunction>
number L2Norm2(TGridFunction& u, const char* cmp,
              int quadOrder, const char* subsets)
{
//	get function id of name
	const size_t fct = u.fct_id_by_name(cmp);

//	check that function exists
	if(fct >= u.num_fct())
		UG_THROW("L2Norm: Function space does not contain"
				" a function with name " << cmp << ".");

	L2Integrand<TGridFunction> integrandL2(u, fct);
	return IntegrateSubsets(integrandL2, u, subsets, quadOrder);
}

template <typename TGridFunction>
number L2Norm(TGridFunction& u, const char* cmp,
              int quadOrder, const char* subsets)
{
	return sqrt(L2Norm2(u, cmp, quadOrder, subsets));
}
/**
 * This function computes the L2-norm of a grid function on all full-dim subsets.
 *
 * \param[in]		spGridFct	grid function
 * \param[in]		cmp			symbolic name of function
 * \param[in]		quadOrder	order of quadrature rule
 * \returns			number 		l2-norm
 */
template <typename TGridFunction>
number L2Norm(TGridFunction& gridFct, const char* cmp, int quadOrder)
{ return L2Norm(gridFct, cmp, quadOrder, NULL); }

template <typename TGridFunction>
number L2Norm(SmartPtr<TGridFunction> spGridFct, const char* cmp, int quadOrder, const char* subsets)
{ return L2Norm(*spGridFct, cmp, quadOrder, subsets); }

template <typename TGridFunction>
number L2Norm(SmartPtr<TGridFunction> spGridFct, const char* cmp, int quadOrder)
{ return L2Norm(spGridFct, cmp, quadOrder, NULL); }


/// Integrand for the distance of two grid functions - evaluated in the (weighted) H1-semi norm
template <typename TGridFunction>
class L2DistIntegrand
		: public StdIntegrand<number, TGridFunction::dim, L2DistIntegrand<TGridFunction> >
{
	public:
	///	world dimension of grid function
		static const int worldDim = TGridFunction::dim;
		typedef typename L2Integrand<TGridFunction>::weight_type weight_type;

	protected:
		ScalarGridFunctionData<TGridFunction> m_fineData;
		const int m_fineTopLevel;

		ScalarGridFunctionData<TGridFunction> m_coarseData;
		const int m_coarseTopLevel;

	///	multigrid
		SmartPtr<MultiGrid> m_spMG;

		ConstSmartPtr<weight_type> m_spWeight;

	/// shift
		double m_deltaFineCoarse;

	public:

		/// constructor (1st is fine grid function)
		L2DistIntegrand(TGridFunction& fineGridFct, size_t fineCmp,
					TGridFunction& coarseGridFct, size_t coarseCmp)
		: m_fineData(fineGridFct, fineCmp), m_fineTopLevel(fineGridFct.dof_distribution()->grid_level().level()),
		  m_coarseData(coarseGridFct, coarseCmp), m_coarseTopLevel(coarseGridFct.dof_distribution()->grid_level().level()),
		  m_spMG(m_fineData.domain()->grid()), m_spWeight(make_sp(new ConstUserNumber<TGridFunction::dim>(1.0))),
		  m_deltaFineCoarse(0.0)
		{
			UG_COND_THROW(m_fineTopLevel < m_coarseTopLevel,
				"L2DiffIntegrand: fine and top level inverted.");
			UG_COND_THROW(m_fineData.domain().get() != m_coarseData.domain().get(),
				"L2DiffIntegrand: grid functions defined on different domains.");
		};

		/// constructor (1st is fine grid function)
		L2DistIntegrand(TGridFunction& fineGridFct, size_t fineCmp,
					TGridFunction& coarseGridFct, size_t coarseCmp, ConstSmartPtr<weight_type> spWeight)
		: m_fineData(fineGridFct, fineCmp), m_fineTopLevel(fineGridFct.dof_distribution()->grid_level().level()),
		  m_coarseData(coarseGridFct, coarseCmp), m_coarseTopLevel(coarseGridFct.dof_distribution()->grid_level().level()),
		  m_spMG(m_fineData.domain()->grid()), m_spWeight(spWeight),
		  m_deltaFineCoarse(0.0)
		{
			UG_COND_THROW(m_fineTopLevel < m_coarseTopLevel,
					"L2DiffIntegrand: fine and top level inverted.");
			UG_COND_THROW(m_fineData.domain().get() != m_coarseData.domain().get(),
					"L2DiffIntegrand: grid functions defined on different domains.");
		};

		/// constructor (1st is fine grid function)
		L2DistIntegrand(TGridFunction& fineGridFct, size_t fineCmp,
				TGridFunction& coarseGridFct, size_t coarseCmp, ConstSmartPtr<weight_type> spWeight, number dist12)
		: m_fineData(fineGridFct, fineCmp), m_fineTopLevel(fineGridFct.dof_distribution()->grid_level().level()),
		  m_coarseData(coarseGridFct, coarseCmp), m_coarseTopLevel(coarseGridFct.dof_distribution()->grid_level().level()),
		  m_spMG(m_fineData.domain()->grid()), m_spWeight(spWeight),
		  m_deltaFineCoarse(dist12)
		{
			UG_COND_THROW(m_fineTopLevel < m_coarseTopLevel,
					"L2DiffIntegrand: fine and top level inverted.");
			UG_COND_THROW(m_fineData.domain().get() != m_coarseData.domain().get(),
					"L2DiffIntegrand: grid functions defined on different domains.");
		};


		virtual ~L2DistIntegrand() {}

		///	sets subset
		virtual void set_subset(int si)
		{
			UG_COND_THROW(!m_fineData.is_def_in_subset(si),
				"L2DiffIntegrand: Grid function component" <<m_fineData.fct()<<" not defined on subset "<<si);
			UG_COND_THROW(!m_coarseData.is_def_in_subset(si),
				"L2DiffIntegrand: Grid function component" <<m_coarseData.fct()<<" not defined on subset "<<si);
			IIntegrand<number, worldDim>::set_subset(si);
		}

	/// \copydoc IIntegrand::values
		template <int elemDim>
		void evaluate(number vValue[],
		              const MathVector<worldDim> vFineGlobIP[],
		              GridObject* pFineElem,
		              const MathVector<worldDim> vCornerCoords[],
		              const MathVector<elemDim> vFineLocIP[],
		              const MathMatrix<elemDim, worldDim> vJT[],
		              const size_t numIP)
		{
			typedef typename TGridFunction::template dim_traits<elemDim>::grid_base_object Element;

			//	get coarse element
			GridObject* pCoarseElem = pFineElem;
			if(m_coarseTopLevel < m_fineTopLevel){
				int parentLevel = m_spMG->get_level(pCoarseElem);
				while(parentLevel > m_coarseTopLevel){
					pCoarseElem = m_spMG->get_parent(pCoarseElem);
					parentLevel = m_spMG->get_level(pCoarseElem);
				}
			}

		//	get reference object id (i.e. Triangle, Quadrilateral, Tetrahedron, ...)
			const ReferenceObjectID fineROID = pFineElem->reference_object_id();
			const ReferenceObjectID coarseROID = pCoarseElem->reference_object_id();

		//	get corner coordinates
			std::vector<MathVector<worldDim> > vCornerCoarse;
			CollectCornerCoordinates(vCornerCoarse, *static_cast<Element*>(pCoarseElem), *m_coarseData.domain());

		//	get Reference Mapping
			DimReferenceMapping<elemDim, worldDim>& map
				= ReferenceMappingProvider::get<elemDim, worldDim>(coarseROID, vCornerCoarse);

			std::vector<MathVector<elemDim> > vCoarseLocIP;
			vCoarseLocIP.resize(numIP);
			for(size_t ip = 0; ip < vCoarseLocIP.size(); ++ip) VecSet(vCoarseLocIP[ip], 0.0);
			map.global_to_local(&vCoarseLocIP[0], vFineGlobIP, numIP);

		// element weights
			typedef typename weight_type::data_type ipdata_type;
			std::vector<ipdata_type> fineElemWeights(numIP, 1.0);
			UG_ASSERT(m_spWeight.valid(), "L2DistIntegrand::evaluate requires valid weights! ");
			(*m_spWeight)(&fineElemWeights[0], vFineGlobIP, 0.0, IIntegrand<number, worldDim>::subset(), numIP);

		try{
		//	get trial space
			const LocalShapeFunctionSet<elemDim>& rFineLSFS =
					LocalFiniteElementProvider::get<elemDim>(fineROID, m_fineData.id());
			const LocalShapeFunctionSet<elemDim>& rCoarseLSFS =
					LocalFiniteElementProvider::get<elemDim>(coarseROID, m_coarseData.id());

		//	get multiindices of element
			std::vector<DoFIndex> vFineMI, vCoarseMI;
			m_fineData.dof_indices(pFineElem, vFineMI);
			m_coarseData.dof_indices(pCoarseElem, vCoarseMI);

		//	loop all integration points
			for(size_t ip = 0; ip < numIP; ++ip)
			{
			// 	compute approximated solution at integration point
				number fineSolIP = 0.0;
				for(size_t sh = 0; sh < vFineMI.size(); ++sh)
				{
				//	get value at shape point (e.g. corner for P1 fct)
					const number val = DoFRef(m_fineData.grid_function(), vFineMI[sh]);

				//	add shape fct at ip * value at shape
					fineSolIP += val * rFineLSFS.shape(sh, vFineLocIP[ip]);
				}
				number coarseSolIP = 0.0;
				for(size_t sh = 0; sh < vCoarseMI.size(); ++sh)
				{
				//	get value at shape point (e.g. corner for P1 fct)
					const number val = DoFRef(m_coarseData.grid_function(), vCoarseMI[sh]);

				//	add shape fct at ip * value at shape
					coarseSolIP += val * rCoarseLSFS.shape(sh, vCoarseLocIP[ip]);
				}

			//	get squared of difference
				vValue[ip] = fineElemWeights[ip]*(fineSolIP - coarseSolIP -m_deltaFineCoarse)*(fineSolIP-coarseSolIP-m_deltaFineCoarse);
			}

			}
			UG_CATCH_THROW("L2DistIntegrand::evaluate: trial space missing.");
		};
};


/// computes the squared l2 distance between two functions
template <typename TGridFunction>
number L2Distance2(TGridFunction& spGridFct1, const char* cmp1,
               TGridFunction& spGridFct2, const char* cmp2,
               int quadOrder, const char* subsets,
			   ConstSmartPtr<typename L2Integrand<TGridFunction>::weight_type> spWeight, number avgDist12=0.0)
{
	return GridFunctionDistance2<L2DistIntegrand<TGridFunction>, TGridFunction>
		(spGridFct1, cmp1, spGridFct2, cmp2, quadOrder, subsets, spWeight, avgDist12);
}


/// computes the squared l2 distance between two functions
template <typename TGridFunction>
number L2Distance2(TGridFunction& spGridFct1, const char* cmp1,
               TGridFunction& spGridFct2, const char* cmp2,
               int quadOrder, const char* subsets)
{
	return GridFunctionDistance2<L2DistIntegrand<TGridFunction>, TGridFunction>
		(spGridFct1, cmp1, spGridFct2, cmp2, quadOrder, subsets);
}

/// computes the l2 distance between two functions
template <typename TGridFunction>
number L2Distance(TGridFunction& spGridFct1, const char* cmp1,
               TGridFunction& spGridFct2, const char* cmp2,
               int quadOrder, const char* subsets)
{
	return sqrt(L2Distance2(spGridFct1, cmp1, spGridFct2, cmp2, quadOrder, subsets));
}


template <typename TGridFunction>
number L2Error(SmartPtr<TGridFunction> spGridFct1, const char* cmp1,
               SmartPtr<TGridFunction> spGridFct2, const char* cmp2,
               int quadOrder)
{
	return L2Distance(*spGridFct1, cmp1, *spGridFct2, cmp2, quadOrder, NULL);
}

template <typename TGridFunction>
number L2Error(SmartPtr<TGridFunction> spGridFct1, const char* cmp1,
               SmartPtr<TGridFunction> spGridFct2, const char* cmp2,
               int quadOrder, const char* subsets)
{
	return L2Distance(*spGridFct1, cmp1, *spGridFct2, cmp2, quadOrder, subsets);
}

////////////////////////////////////////////////////////////////////////////////
// H1 semi-norm Integrand
////////////////////////////////////////////////////////////////////////////////



//! Norm of a grid function, evaluated in (weighted) H1-semi norm
template <typename TGridFunction>
class H1SemiIntegrand
	: public StdIntegrand<number, TGridFunction::dim, H1SemiIntegrand<TGridFunction> >
{
	public:
	///	world dimension of grid function
		static const int worldDim = TGridFunction::dim;
		typedef UserData<MathMatrix<worldDim, worldDim>, worldDim> weight_type;

	protected:
	/// grid function data
		ScalarGridFunctionData<TGridFunction> m_scalarData;

	/// scalar weight (optional)
		ConstSmartPtr<weight_type> m_spWeight;

	public:
	/// constructor
		H1SemiIntegrand(TGridFunction& gridFct, size_t cmp)
		: m_scalarData(gridFct, cmp), m_spWeight(make_sp(new ConstUserMatrix<worldDim>(1.0))) {}

	/// constructor
		H1SemiIntegrand(TGridFunction& gridFct, size_t cmp, ConstSmartPtr<weight_type> spWeight)
		: m_scalarData(gridFct, cmp), m_spWeight(spWeight) {}

	/// DTOR
		virtual ~H1SemiIntegrand() {};

	///	sets subset
		virtual void set_subset(int si)
		{
			if(!m_scalarData.is_def_in_subset(si))
				UG_THROW("H1Error: Grid function component"
						<<m_scalarData.fct()<<" not defined on subset "<<si);
			IIntegrand<number, worldDim>::set_subset(si);
		}

	/// \copydoc IIntegrand::values
		template <int elemDim>
		void evaluate(number vValue[],
		              const MathVector<worldDim> vGlobIP[],
		              GridObject* pElem,
		              const MathVector<worldDim> vCornerCoords[],
		              const MathVector<elemDim> vLocIP[],
		              const MathMatrix<elemDim, worldDim> vJT[],
		              const size_t numIP)
		{
		//	get reference object id (i.e. Triangle, Quadrilateral, Tetrahedron, ...)
			const ReferenceObjectID roid = pElem->reference_object_id();
			const TGridFunction &gridFct= m_scalarData.grid_function();

			typedef typename weight_type::data_type ipdata_type;

			std::vector<ipdata_type> elemWeights(numIP, MathMatrix<worldDim, worldDim>());
			UG_ASSERT(m_spWeight.valid(), "H1SemiIntegrand::evaluate requires valid weights!");


			if(m_spWeight->requires_grid_fct())
			{
				//	get local solution if needed
				LocalIndices ind;
				LocalVector uloc;
				gridFct.indices(pElem, ind);	// 	get global indices
				uloc.resize(ind);				// 	adapt local algebra
				GetLocalVector(uloc, gridFct);	// 	read local values of u

				//	compute data
				try{
					(*m_spWeight)(&elemWeights[0], vGlobIP, 0.0, IIntegrand<number, worldDim>::subset(), pElem,
							vCornerCoords, vLocIP, numIP, &uloc, NULL);
				} UG_CATCH_THROW("H1SemiIntegrand: Cannot evaluate weight data.");
			}
			else
			{
				//	compute data
				try{
					(*m_spWeight)(&elemWeights[0], vGlobIP, 0.0, IIntegrand<number, worldDim>::subset(), numIP);
				} UG_CATCH_THROW("H1SemiIntegrand: Cannot evaluate weight data.");
			}


			try{




		//	get trial space
			const LocalShapeFunctionSet<elemDim>& rTrialSpace =
							LocalFiniteElementProvider::get<elemDim>(roid, m_scalarData.id());

		//	number of dofs on element
			const size_t num_sh = rTrialSpace.num_sh();

		//	get multiindices of element
			std::vector<DoFIndex> ind;  // 	aux. index array
			gridFct.dof_indices(pElem, m_scalarData.fct(), ind);

		//	check multi indices
			if(ind.size() != num_sh)
				UG_THROW("H1SemiNormFuncIntegrand::evaluate: Wrong number of multi-)indices.");

		//	loop all integration points
			std::vector<MathVector<elemDim> > vLocGradient(num_sh);
			for(size_t ip = 0; ip < numIP; ++ip)
			{
			//	compute shape gradients at ip
				rTrialSpace.grads(&vLocGradient[0], vLocIP[ip]);

			// 	compute approximated solution at integration point
				number approxSolIP = 0.0;
				MathVector<elemDim> tmpVec(0.0);
				for(size_t sh = 0; sh < num_sh; ++sh)
				{
				//	get value at shape point (e.g. corner for P1 fct)
					const number valSH = DoFRef(gridFct, ind[sh]);

				//	add shape fct at ip * value at shape
					approxSolIP += valSH * rTrialSpace.shape(sh, vLocIP[ip]);

				//	add gradient at ip
					VecScaleAppend(tmpVec, valSH, vLocGradient[sh]);
				}

			//	compute gradient
				MathVector<worldDim> approxGradIP;
				MathMatrix<worldDim, elemDim> JTInv;
				Inverse(JTInv, vJT[ip]);
				MatVecMult(approxGradIP, JTInv, tmpVec);

			//	get norm squared
				MathVector<worldDim> approxDGradIP;
				MatVecMult(approxDGradIP, elemWeights[ip], approxGradIP);
				vValue[ip] = VecDot(approxDGradIP, approxGradIP);
			}

			}
			UG_CATCH_THROW("H1SemiIntegrand::evaluate: trial space missing.");
		};
};


/// compute H1 semi-norm of a function on the whole domain (or on selected subsets)
/**
 * This function computes the H1 semi-norm of a grid function
 *
 * \param[in]		spGridFct	grid function
 * \param[in]		cmp			symbolic name of component function
 * \param[in]		time		time point
 * \param[in]		quadOrder	order of quadrature rule
 * \param[in]		subsets		subsets, where to compute (OPTIONAL)
 * \param[in]		weights		element-wise weights (OPTIONAL)
 * \returns			number 		l2-norm
 */
template <typename TGridFunction>
number H1SemiNorm2(TGridFunction& gridFct, const char* cmp, int quadOrder, const char* subsets=NULL,
				ConstSmartPtr<typename H1SemiIntegrand<TGridFunction>::weight_type> weights = SPNULL)
{
//	get function id of name
	const size_t fct = gridFct.fct_id_by_name(cmp);

//	check that function exists
	if(fct >= gridFct.num_fct())
		UG_THROW("H1SemiNorm: Function space does not contain"
				" a function with name " << cmp << ".");
	if  (weights.invalid()) {
		H1SemiIntegrand<TGridFunction> integrand(gridFct, fct);
		return IntegrateSubsets(integrand, gridFct, subsets, quadOrder);
	} else {
		H1SemiIntegrand<TGridFunction> integrand(gridFct, fct, weights);
		return IntegrateSubsets(integrand, gridFct, subsets, quadOrder);
	}
}

/// Computes the H1SemiNorm
template <typename TGridFunction>
number H1SemiNorm(TGridFunction& gridFct, const char* cmp, int quadOrder, const char* subsets=NULL,
				ConstSmartPtr<typename H1SemiIntegrand<TGridFunction>::weight_type> weights = SPNULL)
{
	return (sqrt(H1SemiNorm2(gridFct, cmp, quadOrder, subsets, weights)));
}

// Delegating to H1SemiNorm
template <typename TGridFunction>
number H1SemiNorm(SmartPtr<TGridFunction> spGridFct, const char* cmp, int quadOrder, const char* subsets)
{ return H1SemiNorm(*spGridFct, cmp, quadOrder, subsets); }

template <typename TGridFunction>
number H1SemiNorm(SmartPtr<TGridFunction> spGridFct, const char* cmp, int quadOrder,
		const char* subsets, ConstSmartPtr<typename H1SemiIntegrand<TGridFunction>::weight_type> weights = SPNULL)
{ return H1SemiNorm(*spGridFct, cmp, quadOrder, subsets, weights); }

template <typename TGridFunction>
number H1SemiNorm( SmartPtr<TGridFunction> spGridFct, const char* cmp, int quadOrder)
{ return H1SemiNorm(*spGridFct, cmp, quadOrder, NULL); }

template <typename TGridFunction>
number H1SemiNorm( SmartPtr<TGridFunction> spGridFct, const char* cmp, int quadOrder,
		ConstSmartPtr<typename H1SemiIntegrand<TGridFunction>::weight_type> weights)
{ return H1SemiNorm(*spGridFct, cmp, quadOrder, NULL, weights); }

/// Integrand for the distance of two grid functions - evaluated in the (weighted) H1-semi norm
template <typename TGridFunction>
class H1SemiDistIntegrand : public StdIntegrand<number, TGridFunction::dim, H1SemiDistIntegrand<TGridFunction> >
{
	public:
	///	world dimension of grid function
		static const int worldDim = TGridFunction::dim;
		typedef typename H1SemiIntegrand<TGridFunction>::weight_type weight_type;

	private:
		ScalarGridFunctionData<TGridFunction> m_fineData;
		const int m_fineTopLevel;

		ScalarGridFunctionData<TGridFunction> m_coarseData;
		const int m_coarseTopLevel;

	///	multigrid
		SmartPtr<MultiGrid> m_spMG;

	/// scalar weight (optional)
		ConstSmartPtr<weight_type> m_spWeight;

	public:
		/// constructor
		H1SemiDistIntegrand(TGridFunction& fineGridFct, size_t fineCmp,
							TGridFunction& coarseGridFct, size_t coarseCmp)
		: m_fineData(fineGridFct, fineCmp),
		  m_fineTopLevel(fineGridFct.dof_distribution()->grid_level().level()),
		  m_coarseData(coarseGridFct, coarseCmp),
		  m_coarseTopLevel(coarseGridFct.dof_distribution()->grid_level().level()),
		  m_spMG(fineGridFct.domain()->grid()),
		  m_spWeight(new ConstUserNumber<TGridFunction::dim>(1.0))
		{
			if(m_fineTopLevel < m_coarseTopLevel)
				UG_THROW("H1SemiDiffIntegrand: fine and top level inverted.");

			if(m_fineData.domain().get() != m_coarseData.domain().get())
				UG_THROW("H1SemiDiffIntegrand: grid functions defined on different domains.");
		};

		/// constructor
		H1SemiDistIntegrand(TGridFunction& fineGridFct, size_t fineCmp,
							TGridFunction& coarseGridFct, size_t coarseCmp,
							ConstSmartPtr<weight_type> spWeight)
		: m_fineData(fineGridFct, fineCmp),
		  m_fineTopLevel(fineGridFct.dof_distribution()->grid_level().level()),
		  m_coarseData(coarseGridFct, coarseCmp),
		  m_coarseTopLevel(coarseGridFct.dof_distribution()->grid_level().level()),
		  m_spMG(fineGridFct.domain()->grid()),
		  m_spWeight(spWeight)
		{
			if(m_fineTopLevel < m_coarseTopLevel)
				UG_THROW("H1SemiDiffIntegrand: fine and top level inverted.");

			if(m_fineData.domain().get() != m_coarseData.domain().get())
				UG_THROW("H1SemiDiffIntegrand: grid functions defined on different domains.");
		}
		virtual ~H1SemiDistIntegrand(){}

	///	sets subset
		virtual void set_subset(int si)
		{
			if(!m_fineData.is_def_in_subset(si))
				UG_THROW("H1SemiDiffIntegrand: Grid function component"
						<<m_fineData.fct()<<" not defined on subset "<<si);
			if(!m_coarseData.is_def_in_subset(si))
				UG_THROW("H1SemiDiffIntegrand: Grid function component"
						<<m_coarseData.fct()<<" not defined on subset "<<si);
			IIntegrand<number, worldDim>::set_subset(si);
		}

	/// \copydoc IIntegrand::values
		template <int elemDim>
		void evaluate(number vValue[],
		              const MathVector<worldDim> vFineGlobIP[],
		              GridObject* pFineElem,
		              const MathVector<worldDim> vCornerCoords[],
		              const MathVector<elemDim> vFineLocIP[],
		              const MathMatrix<elemDim, worldDim> vJT[],
		              const size_t numIP)
		{
			typedef typename TGridFunction::template dim_traits<elemDim>::grid_base_object Element;

			const TGridFunction &fineGridFct  = m_fineData.grid_function();
			const TGridFunction &coarseGridFct  = m_coarseData.grid_function();

			//	get coarse element
			GridObject* pCoarseElem = pFineElem;
			if(m_coarseTopLevel < m_fineTopLevel){
				int parentLevel = m_spMG->get_level(pCoarseElem);
				while(parentLevel > m_coarseTopLevel){
					pCoarseElem = m_spMG->get_parent(pCoarseElem);
					parentLevel = m_spMG->get_level(pCoarseElem);
				}
			}

		//	get reference object id (i.e. Triangle, Quadrilateral, Tetrahedron, ...)
			const ReferenceObjectID fineROID = pFineElem->reference_object_id();
			const ReferenceObjectID coarseROID = pCoarseElem->reference_object_id();

		//	get corner coordinates
			std::vector<MathVector<worldDim> > vCornerCoarse;
			CollectCornerCoordinates(vCornerCoarse, *static_cast<Element*>(pCoarseElem), *m_coarseData.domain());

		//	get reference Mapping
			DimReferenceMapping<elemDim, worldDim>& map
				= ReferenceMappingProvider::get<elemDim, worldDim>(coarseROID, vCornerCoarse);

			std::vector<MathVector<elemDim> > vCoarseLocIP;
			vCoarseLocIP.resize(numIP);
			for(size_t ip = 0; ip < vCoarseLocIP.size(); ++ip) VecSet(vCoarseLocIP[ip], 0.0);
			map.global_to_local(&vCoarseLocIP[0], vFineGlobIP, numIP);


			// determine weights
			std::vector<typename weight_type::data_type> elemWeights(numIP, MathMatrix<worldDim, worldDim>());
			UG_ASSERT(m_spWeight.valid(), "H1SemiDistIntegrand::evaluate requires valid weights!");

			//	get local solution (if required)
			if(m_spWeight->requires_grid_fct())
			{

				LocalIndices ind;
				LocalVector u;
				fineGridFct.indices(pFineElem, ind);	// 	get global indices
				u.resize(ind);							// 	adapt local algebra
				GetLocalVector(u, fineGridFct);			// 	read local values of u

				//	compute data
				try{
					(*m_spWeight)(&elemWeights[0], vFineGlobIP, 0.0, IIntegrand<number, worldDim>::subset(), pFineElem,
							vCornerCoords, vFineLocIP, numIP, &u, NULL);
				} UG_CATCH_THROW("H1SemiDistIntegrand: Cannot evaluate data.");
			}
			else
			{
				//	compute data
				try{
					(*m_spWeight)(&elemWeights[0], vFineGlobIP, 0.0, IIntegrand<number, worldDim>::subset(), numIP);
				} UG_CATCH_THROW("H1SemiDistIntegrand: Cannot evaluate data.");
			}

			try{


		//	get trial space
			const LocalShapeFunctionSet<elemDim>& rFineLSFS =
					LocalFiniteElementProvider::get<elemDim>(fineROID, m_fineData.id());
			const LocalShapeFunctionSet<elemDim>& rCoarseLSFS =
					LocalFiniteElementProvider::get<elemDim>(coarseROID, m_coarseData.id());

		//	get multiindices of element
			std::vector<DoFIndex> vFineMI, vCoarseMI;
			m_fineData.dof_indices(pFineElem, vFineMI);
			m_coarseData.dof_indices(pCoarseElem, vCoarseMI);

			std::vector<MathVector<elemDim> > vFineLocGradient(vFineMI.size());
			std::vector<MathVector<elemDim> > vCoarseLocGradient(vCoarseMI.size());

		//	loop all integration points
			for(size_t ip = 0; ip < numIP; ++ip)
			{
			//	compute shape gradients at ip
				rFineLSFS.grads(&vFineLocGradient[0], vFineLocIP[ip]);
				rCoarseLSFS.grads(&vCoarseLocGradient[0], vCoarseLocIP[ip]);

			// 	compute approximated solutions at integration point
				number fineSolIP = 0.0;
				MathVector<elemDim> fineLocTmp(0.0);
				for(size_t sh = 0; sh < vFineMI.size(); ++sh)
				{
					const number val = DoFRef(fineGridFct, vFineMI[sh]);
					fineSolIP += val * rFineLSFS.shape(sh, vFineLocIP[ip]);
					VecScaleAppend(fineLocTmp, val, vFineLocGradient[sh]);
				}

				number coarseSolIP = 0.0;
				MathVector<elemDim> coarseLocTmp(0.0);
				for(size_t sh = 0; sh < vCoarseMI.size(); ++sh)
				{
					const number val = DoFRef(coarseGridFct, vCoarseMI[sh]);
					coarseSolIP += val * rCoarseLSFS.shape(sh, vCoarseLocIP[ip]);
					VecScaleAppend(coarseLocTmp, val, vCoarseLocGradient[sh]);
				}

			//	compute global gradient
				MathVector<worldDim> fineGradIP;
				MathMatrix<worldDim, elemDim> fineJTInv;
				Inverse(fineJTInv, vJT[ip]);
				MatVecMult(fineGradIP, fineJTInv, fineLocTmp);

			//	compute global gradient
				MathVector<worldDim> coarseGradIP;
				MathMatrix<worldDim, elemDim> coarseJTInv;
				map.jacobian_transposed_inverse(coarseJTInv, vCoarseLocIP[ip]);
				MatVecMult(coarseGradIP, coarseJTInv, coarseLocTmp);

			//	get squared of difference
				/*vValue[ip] = (coarseSolIP - fineSolIP);
				vValue[ip] *= vValue[ip]; */
				vValue[ip] = VecDistanceSq(coarseGradIP, fineGradIP, elemWeights[ip]);
			}

			}
			UG_CATCH_THROW("H1SemiDiffIntegrand::evaluate: trial space missing.");
		};
};


//! Distance in H1 semi norm (with subset selection)
template <typename TGridFunction>
number H1SemiError2(SmartPtr<TGridFunction> spGridFct1, const char* cmp1,
               SmartPtr<TGridFunction> spGridFct2, const char* cmp2,
               int quadOrder, const char* subsets)
{
	return GridFunctionDistance2<H1SemiDistIntegrand<TGridFunction>, TGridFunction>
		(*spGridFct1, cmp1, *spGridFct2, cmp2, quadOrder, subsets);
}

//! Distance in H1 semi norm (with subset selection)
template <typename TGridFunction>
number H1SemiError(SmartPtr<TGridFunction> spGridFct1, const char* cmp1,
               SmartPtr<TGridFunction> spGridFct2, const char* cmp2,
               int quadOrder, const char* subsets)
{
	return sqrt(H1SemiError2(*spGridFct1, cmp1, *spGridFct2, cmp2, quadOrder, subsets));
}

//! Distance in H1 semi norm (all subsets)
template <typename TGridFunction>
number H1SemiError(SmartPtr<TGridFunction> spGridFct1, const char* cmp1,
               SmartPtr<TGridFunction> spGridFct2, const char* cmp2,
               int quadOrder)
{
	return H1SemiError(*spGridFct1, cmp1, *spGridFct2, cmp2, quadOrder, NULL);
}

//! Squared distance in H1 semi norm (with select subsets & weights)
template <typename TGridFunction>
number H1SemiDistance2(TGridFunction& spGridFct1, const char* cmp1,
		TGridFunction& spGridFct2, const char* cmp2,
        int quadOrder, const char* subsets,
		ConstSmartPtr<typename H1SemiDistIntegrand<TGridFunction>::weight_type> weights)
{
	return GridFunctionDistance2<H1SemiDistIntegrand<TGridFunction>, TGridFunction>
		(spGridFct1, cmp1, spGridFct2, cmp2, quadOrder, subsets, weights);
}

//! Distance in H1 semi norm (with select subsets & weights)
template <typename TGridFunction>
number H1SemiDistance(TGridFunction& spGridFct1, const char* cmp1,
		TGridFunction& spGridFct2, const char* cmp2,
        int quadOrder, const char* subsets,
		ConstSmartPtr<typename H1SemiDistIntegrand<TGridFunction>::weight_type> weights)
{ return sqrt(H1SemiDistance2(spGridFct1, cmp1, spGridFct2, cmp2, quadOrder, subsets, weights)); }

//! Squared distance in H1 semi norm (all subsets, with weights)
template <typename TGridFunction>
number H1SemiDistance2(TGridFunction& spGridFct1, const char* cmp1,
					  TGridFunction& spGridFct2, const char* cmp2,
               int quadOrder, ConstSmartPtr<typename H1SemiDistIntegrand<TGridFunction>::weight_type > weights)
{ return H1SemiDistance2(spGridFct1, cmp1, spGridFct2, cmp2, quadOrder, NULL, weights); }

//! Distance in H1 semi norm (all subsets, with weights)
template <typename TGridFunction>
number H1SemiDistance(TGridFunction& spGridFct1, const char* cmp1,
					  TGridFunction& spGridFct2, const char* cmp2,
               int quadOrder, ConstSmartPtr<typename H1SemiDistIntegrand<TGridFunction>::weight_type > weights)
{ return sqrt(H1SemiDistance2(spGridFct1, cmp1, spGridFct2, cmp2, quadOrder, NULL, weights)); }




////////////////////////////////////////////////////////////////////////////////
// H1 semi-norm Integrand
////////////////////////////////////////////////////////////////////////////////

//! Norm of a grid function, evaluated in (weighted) H1-semi norm
template <typename TGridFunction>
class H1EnergyIntegrand
	: public StdIntegrand<number, TGridFunction::dim, H1EnergyIntegrand<TGridFunction> >
{
	public:
	///	world dimension of grid function
		static const int worldDim = TGridFunction::dim;
		typedef UserData<MathMatrix<worldDim, worldDim>, worldDim> weight_type;

	protected:
	/// grid function data
		ScalarGridFunctionData<TGridFunction> m_scalarData;

	/// scalar weight (optional)
		ConstSmartPtr<weight_type> m_spWeight;

	public:
	/// constructor
		H1EnergyIntegrand(TGridFunction& gridFct, size_t cmp)
		: m_scalarData(gridFct, cmp), m_spWeight(make_sp(new ConstUserMatrix<worldDim>(1.0))) {}

	/// constructor
		H1EnergyIntegrand(TGridFunction& gridFct, size_t cmp, ConstSmartPtr<weight_type> spWeight)
		: m_scalarData(gridFct, cmp), m_spWeight(spWeight) {}

	/// DTOR
		virtual ~H1EnergyIntegrand() {};

	///	sets subset
		virtual void set_subset(int si)
		{
			if(!m_scalarData.is_def_in_subset(si))
				UG_THROW("H1EnergyIntegrand: Grid function component"
						<<m_scalarData.fct()<<" not defined on subset "<<si);
			IIntegrand<number, worldDim>::set_subset(si);
		}

	/// \copydoc IIntegrand::values
		template <int elemDim>
		void evaluate(number vValue[],
		              const MathVector<worldDim> vGlobIP[],
		              GridObject* pElem,
		              const MathVector<worldDim> vCornerCoords[],
		              const MathVector<elemDim> vLocIP[],
		              const MathMatrix<elemDim, worldDim> vJT[],
		              const size_t numIP)
		{
		//	get reference object id (i.e. Triangle, Quadrilateral, Tetrahedron, ...)
			const ReferenceObjectID roid = pElem->reference_object_id();
			const TGridFunction &gridFct= m_scalarData.grid_function();

			typedef typename weight_type::data_type ipdata_type;

			std::vector<ipdata_type> elemWeights(numIP, MathMatrix<worldDim, worldDim>());
			UG_ASSERT(m_spWeight.valid(), "H1EnergyIntegrand::evaluate requires valid weights!");


			if(m_spWeight->requires_grid_fct())
			{
				//	get local solution if needed
				LocalIndices ind;
				LocalVector uloc;
				gridFct.indices(pElem, ind);	// 	get global indices
				uloc.resize(ind);				// 	adapt local algebra
				GetLocalVector(uloc, gridFct);	// 	read local values of u

				//	compute data
				try{
					(*m_spWeight)(&elemWeights[0], vGlobIP, 0.0, IIntegrand<number, worldDim>::subset(), pElem,
							vCornerCoords, vLocIP, numIP, &uloc, NULL);
				} UG_CATCH_THROW("H1EnergyIntegrand: Cannot evaluate weight data.");
			}
			else
			{
				//	compute data
				try{
					(*m_spWeight)(&elemWeights[0], vGlobIP, 0.0, IIntegrand<number, worldDim>::subset(), numIP);
				} UG_CATCH_THROW("H1EnergyIntegrand: Cannot evaluate weight data.");
			}


			try{




		//	get trial space
			const LocalShapeFunctionSet<elemDim>& rTrialSpace =
							LocalFiniteElementProvider::get<elemDim>(roid, m_scalarData.id());

		//	number of dofs on element
			const size_t num_sh = rTrialSpace.num_sh();

		//	get multiindices of element
			std::vector<DoFIndex> ind;  // 	aux. index array
			gridFct.dof_indices(pElem, m_scalarData.fct(), ind);

		//	check multi indices
			if(ind.size() != num_sh)
				UG_THROW("H1EnergyIntegrand::evaluate: Wrong number of multi-)indices.");

		//	loop all integration points
			std::vector<MathVector<elemDim> > vLocGradient(num_sh);
			for(size_t ip = 0; ip < numIP; ++ip)
			{
			//	compute shape gradients at ip
				rTrialSpace.grads(&vLocGradient[0], vLocIP[ip]);

			// 	compute approximated solution at integration point
				number approxSolIP = 0.0;
				MathVector<elemDim> tmpVec(0.0);
				for(size_t sh = 0; sh < num_sh; ++sh)
				{
				//	get value at shape point (e.g. corner for P1 fct)
					const number valSH = DoFRef(gridFct, ind[sh]);

				//	add shape fct at ip * value at shape
					approxSolIP += valSH * rTrialSpace.shape(sh, vLocIP[ip]);

				//	add gradient at ip
					VecScaleAppend(tmpVec, valSH, vLocGradient[sh]);
				}

			//	compute gradient
				MathVector<worldDim> approxGradIP;
				MathMatrix<worldDim, elemDim> JTInv;
				Inverse(JTInv, vJT[ip]);
				MatVecMult(approxGradIP, JTInv, tmpVec);

			//	get norm squared
				MathVector<worldDim> approxDGradIP;
				MatVecMult(approxDGradIP, elemWeights[ip], approxGradIP);
				vValue[ip] = VecTwoNormSq(approxDGradIP);
			}

			}
			UG_CATCH_THROW("H1EnergyIntegrand::evaluate: trial space missing.");
		};
};


/// compute energy -norm of a function on the whole domain (or on selected subsets)
/**
 * This function computes the integral over  \f$ \| q  \|^2 \f where the velocity is given by \f$ q:= \kappa \nabla u\f$ of a grid function u.
 *
 * \param[in]		spGridFct	grid function
 * \param[in]		cmp			symbolic name of component function
 * \param[in]		time		time point
 * \param[in]		quadOrder	order of quadrature rule
 * \param[in]		subsets		subsets, where to compute (OPTIONAL)
 * \param[in]		weights		element-wise matrix-valued weights kappa (OPTIONAL)
 * \returns			number 		l2-norm
 */
template <typename TGridFunction>
number H1EnergyNorm2(TGridFunction& gridFct, const char* cmp, int quadOrder, const char* subsets=NULL,
				ConstSmartPtr<typename H1SemiIntegrand<TGridFunction>::weight_type> weights = SPNULL)
{
//	get function id of name
	const size_t fct = gridFct.fct_id_by_name(cmp);

//	check that function exists
	if(fct >= gridFct.num_fct())
		UG_THROW("H1SemiNorm: Function space does not contain"
				" a function with name " << cmp << ".");
	if  (weights.invalid()) {
		H1EnergyIntegrand<TGridFunction> integrand(gridFct, fct);
		return IntegrateSubsets(integrand, gridFct, subsets, quadOrder);
	} else {
		H1EnergyIntegrand<TGridFunction> integrand(gridFct, fct, weights);
		return IntegrateSubsets(integrand, gridFct, subsets, quadOrder);
	}
}
template <typename TGridFunction>
number H1EnergyNorm(TGridFunction& gridFct, const char* cmp, int quadOrder, const char* subsets=NULL,
				ConstSmartPtr<typename H1SemiIntegrand<TGridFunction>::weight_type> weights = SPNULL)
{
	return (sqrt(H1EnergyNorm2(gridFct, cmp, quadOrder, subsets, weights)));
}

template <typename TGridFunction>
number H1EnergyNorm(SmartPtr<TGridFunction> spGridFct, const char* cmp, int quadOrder,
		const char* subsets, ConstSmartPtr<typename H1SemiIntegrand<TGridFunction>::weight_type> weights = SPNULL)
{ return H1EnergyNorm(*spGridFct, cmp, quadOrder, subsets, weights); }

template <typename TGridFunction>
number H1EnergyNorm( SmartPtr<TGridFunction> spGridFct, const char* cmp, int quadOrder)
{ return H1EnergyNorm(spGridFct, cmp, quadOrder, NULL); }

template <typename TGridFunction>
number H1EnergyNorm( SmartPtr<TGridFunction> spGridFct, const char* cmp, int quadOrder,
		ConstSmartPtr<typename H1SemiIntegrand<TGridFunction>::weight_type> weights)
{ return H1EnergyNorm(spGridFct, cmp, quadOrder, NULL, weights); }




/// Integrand for the distance of two grid functions - evaluated in the norm |D \nabla u|^2
template <typename TGridFunction>
class H1EnergyDistIntegrand
		: public StdIntegrand<number, TGridFunction::dim, H1EnergyDistIntegrand<TGridFunction> >
{
	public:
	///	world dimension of grid function
		static const int worldDim = TGridFunction::dim;
		typedef typename H1SemiIntegrand<TGridFunction>::weight_type weight_type;

	private:
		ScalarGridFunctionData<TGridFunction> m_fineData;
		const int m_fineTopLevel;

		ScalarGridFunctionData<TGridFunction> m_coarseData;
		const int m_coarseTopLevel;

	///	multigrid
		SmartPtr<MultiGrid> m_spMG;

	/// scalar weight (optional)
		ConstSmartPtr<weight_type> m_spWeight;

	public:
		/// constructor
		H1EnergyDistIntegrand(TGridFunction& fineGridFct, size_t fineCmp,
							TGridFunction& coarseGridFct, size_t coarseCmp)
		: m_fineData(fineGridFct, fineCmp),
		  m_fineTopLevel(fineGridFct.dof_distribution()->grid_level().level()),
		  m_coarseData(coarseGridFct, coarseCmp),
		  m_coarseTopLevel(coarseGridFct.dof_distribution()->grid_level().level()),
		  m_spMG(fineGridFct.domain()->grid()),
		  m_spWeight(new ConstUserNumber<TGridFunction::dim>(1.0))
		{
			if(m_fineTopLevel < m_coarseTopLevel)
				UG_THROW("H1EnergyDistIntegrand: fine and top level inverted.");

			if(m_fineData.domain().get() != m_coarseData.domain().get())
				UG_THROW("H1EnergyDistIntegrand: grid functions defined on different domains.");
		};

		/// constructor
		H1EnergyDistIntegrand(TGridFunction& fineGridFct, size_t fineCmp,
							TGridFunction& coarseGridFct, size_t coarseCmp,
							ConstSmartPtr<weight_type> spWeight)
		: m_fineData(fineGridFct, fineCmp),
		  m_fineTopLevel(fineGridFct.dof_distribution()->grid_level().level()),
		  m_coarseData(coarseGridFct, coarseCmp),
		  m_coarseTopLevel(coarseGridFct.dof_distribution()->grid_level().level()),
		  m_spMG(fineGridFct.domain()->grid()),
		  m_spWeight(spWeight)
		{
			if(m_fineTopLevel < m_coarseTopLevel)
				UG_THROW("H1EnergyDistIntegrand: fine and top level inverted.");

			if(m_fineData.domain().get() != m_coarseData.domain().get())
				UG_THROW("H1EnergyDistIntegrand: grid functions defined on different domains.");
		}
		virtual ~H1EnergyDistIntegrand(){}

	///	sets subset
		virtual void set_subset(int si)
		{
			if(!m_fineData.is_def_in_subset(si))
				UG_THROW("H1EnergyDistIntegrand: Grid function component"
						<<m_fineData.fct()<<" not defined on subset "<<si);
			if(!m_coarseData.is_def_in_subset(si))
				UG_THROW("H1EnergyDistIntegrand: Grid function component"
						<<m_coarseData.fct()<<" not defined on subset "<<si);
			IIntegrand<number, worldDim>::set_subset(si);
		}

	/// \copydoc IIntegrand::values
		template <int elemDim>
		void evaluate(number vValue[],
		              const MathVector<worldDim> vFineGlobIP[],
		              GridObject* pFineElem,
		              const MathVector<worldDim> vCornerCoords[],
		              const MathVector<elemDim> vFineLocIP[],
		              const MathMatrix<elemDim, worldDim> vJT[],
		              const size_t numIP)
		{
			typedef typename TGridFunction::template dim_traits<elemDim>::grid_base_object Element;

			const TGridFunction &fineGridFct  = m_fineData.grid_function();
			const TGridFunction &coarseGridFct  = m_coarseData.grid_function();

			//	get coarse element
			GridObject* pCoarseElem = pFineElem;
			if(m_coarseTopLevel < m_fineTopLevel){
				int parentLevel = m_spMG->get_level(pCoarseElem);
				while(parentLevel > m_coarseTopLevel){
					pCoarseElem = m_spMG->get_parent(pCoarseElem);
					parentLevel = m_spMG->get_level(pCoarseElem);
				}
			}

		//	get reference object id (i.e. Triangle, Quadrilateral, Tetrahedron, ...)
			const ReferenceObjectID fineROID = pFineElem->reference_object_id();
			const ReferenceObjectID coarseROID = pCoarseElem->reference_object_id();

		//	get corner coordinates
			std::vector<MathVector<worldDim> > vCornerCoarse;
			CollectCornerCoordinates(vCornerCoarse, *static_cast<Element*>(pCoarseElem), *m_coarseData.domain());

		//	get reference Mapping
			DimReferenceMapping<elemDim, worldDim>& map
				= ReferenceMappingProvider::get<elemDim, worldDim>(coarseROID, vCornerCoarse);

			std::vector<MathVector<elemDim> > vCoarseLocIP;
			vCoarseLocIP.resize(numIP);
			for(size_t ip = 0; ip < vCoarseLocIP.size(); ++ip) VecSet(vCoarseLocIP[ip], 0.0);
			map.global_to_local(&vCoarseLocIP[0], vFineGlobIP, numIP);


			// determine weights
			std::vector<typename weight_type::data_type> elemWeights(numIP, MathMatrix<worldDim, worldDim>());
			UG_ASSERT(m_spWeight.valid(), "H1SemiDistIntegrand::evaluate requires valid weights!");

			//	get local solution if needed
			if(m_spWeight->requires_grid_fct())
			{

				LocalIndices ind;
				LocalVector u;
				fineGridFct.indices(pFineElem, ind);	// 	get global indices
				u.resize(ind);							// 	adapt local algebra
				GetLocalVector(u, fineGridFct);			// 	read local values of u

				//	compute data
				try{
					(*m_spWeight)(&elemWeights[0], vFineGlobIP, 0.0, IIntegrand<number, worldDim>::subset(), pFineElem,
							vCornerCoords, vFineLocIP, numIP, &u, NULL);
				} UG_CATCH_THROW("H1SemiDistIntegrand: Cannot evaluate data.");
			}
			else
			{
				//	compute data
				try{
					(*m_spWeight)(&elemWeights[0], vFineGlobIP, 0.0, IIntegrand<number, worldDim>::subset(), numIP);
				} UG_CATCH_THROW("H1SemiDistIntegrand: Cannot evaluate data.");
			}

			try{


		//	get trial space
			const LocalShapeFunctionSet<elemDim>& rFineLSFS =
					LocalFiniteElementProvider::get<elemDim>(fineROID, m_fineData.id());
			const LocalShapeFunctionSet<elemDim>& rCoarseLSFS =
					LocalFiniteElementProvider::get<elemDim>(coarseROID, m_coarseData.id());

		//	get multiindices of element
			std::vector<DoFIndex> vFineMI, vCoarseMI;
			m_fineData.dof_indices(pFineElem, vFineMI);
			m_coarseData.dof_indices(pCoarseElem, vCoarseMI);

			std::vector<MathVector<elemDim> > vFineLocGradient(vFineMI.size());
			std::vector<MathVector<elemDim> > vCoarseLocGradient(vCoarseMI.size());

		//	loop all integration points
			for(size_t ip = 0; ip < numIP; ++ip)
			{
			//	compute shape gradients at ip
				rFineLSFS.grads(&vFineLocGradient[0], vFineLocIP[ip]);
				rCoarseLSFS.grads(&vCoarseLocGradient[0], vCoarseLocIP[ip]);

			// 	compute approximated solutions at integration point
				number fineSolIP = 0.0;
				MathVector<elemDim> fineLocTmp(0.0);
				for(size_t sh = 0; sh < vFineMI.size(); ++sh)
				{
					const number val = DoFRef(fineGridFct, vFineMI[sh]);
					fineSolIP += val * rFineLSFS.shape(sh, vFineLocIP[ip]);
					VecScaleAppend(fineLocTmp, val, vFineLocGradient[sh]);
				}

				number coarseSolIP = 0.0;
				MathVector<elemDim> coarseLocTmp(0.0);
				for(size_t sh = 0; sh < vCoarseMI.size(); ++sh)
				{
					const number val = DoFRef(coarseGridFct, vCoarseMI[sh]);
					coarseSolIP += val * rCoarseLSFS.shape(sh, vCoarseLocIP[ip]);
					VecScaleAppend(coarseLocTmp, val, vCoarseLocGradient[sh]);
				}

			//	compute global D*gradient
				MathVector<worldDim> fineGradIP;
				MathMatrix<worldDim, elemDim> fineJTInv;
				Inverse(fineJTInv, vJT[ip]);
				MatVecMult(fineGradIP, fineJTInv, fineLocTmp);
				MatVecMult(fineLocTmp, elemWeights[ip], fineGradIP);

			//	compute global D*gradient
				MathVector<worldDim> coarseGradIP;
				MathMatrix<worldDim, elemDim> coarseJTInv;
				map.jacobian_transposed_inverse(coarseJTInv, vCoarseLocIP[ip]);
				MatVecMult(coarseGradIP, coarseJTInv, coarseLocTmp);
				MatVecMult(coarseLocTmp, elemWeights[ip], coarseGradIP);

			//	get squared difference
				vValue[ip] = VecDistanceSq(fineLocTmp, coarseLocTmp);

			}

			}
			UG_CATCH_THROW("H1EnergyDiffIntegrand::evaluate: trial space missing.");
		};
};


//! Squared distance in H1 semi norm (with select subsets & weights)
template <typename TGridFunction>
number H1EnergyDistance2(TGridFunction& spGridFct1, const char* cmp1,
		TGridFunction& spGridFct2, const char* cmp2,
        int quadOrder, const char* subsets,
		ConstSmartPtr<typename H1SemiDistIntegrand<TGridFunction>::weight_type> weights)
{
	return GridFunctionDistance2<H1EnergyDistIntegrand<TGridFunction>, TGridFunction>
		(spGridFct1, cmp1, spGridFct2, cmp2, quadOrder, subsets, weights);
}

//! Distance in H1 semi norm (with select subsets & weights)
template <typename TGridFunction>
number H1EnergyDistance(TGridFunction& spGridFct1, const char* cmp1,
		TGridFunction& spGridFct2, const char* cmp2,
        int quadOrder, const char* subsets,
		ConstSmartPtr<typename H1SemiDistIntegrand<TGridFunction>::weight_type> weights)
{ return sqrt(H1EnergyDistance2(spGridFct1, cmp1, spGridFct2, cmp2, quadOrder, subsets, weights)); }

//! Squared distance in H1 semi norm (all subsets, with weights)
template <typename TGridFunction>
number H1EnergyDistance2(TGridFunction& spGridFct1, const char* cmp1,
					  TGridFunction& spGridFct2, const char* cmp2,
               int quadOrder, ConstSmartPtr<typename H1SemiDistIntegrand<TGridFunction>::weight_type > weights)
{ return H1EnergyDistance2(spGridFct1, cmp1, spGridFct2, cmp2, quadOrder, NULL, weights); }

//! Distance in H1 semi norm (all subsets, with weights)
template <typename TGridFunction>
number H1EnergyDistance(TGridFunction& spGridFct1, const char* cmp1,
					  TGridFunction& spGridFct2, const char* cmp2,
               int quadOrder, ConstSmartPtr<typename H1SemiDistIntegrand<TGridFunction>::weight_type > weights)
{ return sqrt(H1EnergyDistance2(spGridFct1, cmp1, spGridFct2, cmp2, quadOrder, NULL, weights)); }




////////////////////////////////////////////////////////////////////////////////
// H1 norm integrand
////////////////////////////////////////////////////////////////////////////////
template <typename TGridFunction>
class H1NormIntegrand
	: public StdIntegrand<number, TGridFunction::dim, H1NormIntegrand<TGridFunction> >
{
	public:
	///	world dimension of grid function
		static const int worldDim = TGridFunction::dim;

	private:
		ScalarGridFunctionData<TGridFunction> m_scalarData;

	public:
	/// CTOR
		H1NormIntegrand(TGridFunction& gridFct, size_t cmp)
		: m_scalarData(gridFct, cmp) {}

	/// DTOR
		virtual ~H1NormIntegrand() {}

	///	sets subset
		virtual void set_subset(int si)
		{
			if(!m_scalarData.is_def_in_subset(si))
				UG_THROW("H1Norm: Grid function component"
						<<m_scalarData.fct()<<" not defined on subset "<<si);
			IIntegrand<number, worldDim>::set_subset(si);
		}

	/// \copydoc IIntegrand::values
		template <int elemDim>
		void evaluate(number vValue[],
		              const MathVector<worldDim> vGlobIP[],
		              GridObject* pElem,
		              const MathVector<worldDim> vCornerCoords[],
		              const MathVector<elemDim> vLocIP[],
		              const MathMatrix<elemDim, worldDim> vJT[],
		              const size_t numIP)
		{
		//	get reference object id (i.e. Triangle, Quadrilateral, Tetrahedron, ...)
			const ReferenceObjectID roid = pElem->reference_object_id();

			try{
		//	get trial space
			const LocalShapeFunctionSet<elemDim>& rTrialSpace =
							LocalFiniteElementProvider::get<elemDim>(roid, m_scalarData.id());

		//	number of dofs on element
			const size_t num_sh = rTrialSpace.num_sh();

		//	get multiindices of element
			std::vector<DoFIndex> ind;  // 	aux. index array
			m_scalarData.dof_indices(pElem, ind);

		//	check multi indices
			if(ind.size() != num_sh)
				UG_THROW("H1ErrorIntegrand::evaluate: Wrong number of"
						" multi indices.");

		//	loop all integration points
			std::vector<MathVector<elemDim> > vLocGradient(num_sh);
			for(size_t ip = 0; ip < numIP; ++ip)
			{

			//	compute shape gradients at ip
				rTrialSpace.grads(&vLocGradient[0], vLocIP[ip]);

			// 	compute approximated solution at integration point
				number approxSolIP = 0.0;
				MathVector<elemDim> locTmp; VecSet(locTmp, 0.0);
				for(size_t sh = 0; sh < num_sh; ++sh)
				{
				//	get value at shape point (e.g. corner for P1 fct)
					const number valSH = DoFRef(m_scalarData.grid_function(), ind[sh]);

				//	add shape fct at ip * value at shape
					approxSolIP += valSH * rTrialSpace.shape(sh, vLocIP[ip]);

				//	add gradient at ip
					VecScaleAppend(locTmp, valSH, vLocGradient[sh]);
				}

			//	compute global gradient
				MathVector<worldDim> approxGradIP;
				MathMatrix<worldDim, elemDim> JTInv;
				Inverse(JTInv, vJT[ip]);
				MatVecMult(approxGradIP, JTInv, locTmp);

			//	get squared of difference
				vValue[ip] = approxSolIP * approxSolIP;
				vValue[ip] += VecDot(approxGradIP, approxGradIP);
			}

			}
			UG_CATCH_THROW("H1SemiNormFuncIntegrand::evaluate: trial space missing.");
		};
};




template <typename TGridFunction>
number H1Norm2(TGridFunction& u, const char* cmp,
			   int quadOrder, const char* subsets=NULL)
{
//	get function id of name
	const size_t fct = u.fct_id_by_name(cmp);

//	check that function exists
	if(fct >= u.num_fct())
		UG_THROW("H1Norm: Function space does not contain"
				" a function with name " << cmp << ".");

	H1NormIntegrand<TGridFunction> spIntegrand(u, fct);
	return IntegrateSubsets(spIntegrand, u, subsets, quadOrder);
}


template <typename TGridFunction>
number H1Norm(TGridFunction& u, const char* cmp,
			   int quadOrder, const char* subsets=NULL)
{
	return sqrt(H1Norm2(u, cmp, quadOrder,subsets));
}

template <typename TGridFunction>
number H1Norm(SmartPtr<TGridFunction> spGridFct, const char* cmp,
			   int quadOrder, const char* subsets)
{
	return H1Norm(*spGridFct, cmp, quadOrder, subsets);
}

template <typename TGridFunction>
number H1Norm(SmartPtr<TGridFunction> spGridFct, const char* cmp, int quadOrder)
{
	return H1Norm(*spGridFct, cmp, quadOrder, NULL);
}



/// Integrand for the distance of two grid functions - evaluated in the H1 norm
template <typename TGridFunction>
class H1DistIntegrand
	: public StdIntegrand<number, TGridFunction::dim, H1DistIntegrand<TGridFunction> >
{
	public:
	///	world dimension of grid function
		static const int worldDim = TGridFunction::dim;

	private:
		ScalarGridFunctionData<TGridFunction> m_fineData;
		const int m_fineTopLevel;

		ScalarGridFunctionData<TGridFunction> m_coarseData;
		const int m_coarseTopLevel;

	///	multigrid
		SmartPtr<MultiGrid> m_spMG;

	public:
	/// constructor (1 is fine grid function)
		H1DistIntegrand(TGridFunction& fineGridFct, size_t fineCmp,
		                TGridFunction& coarseGridFct, size_t coarseCmp)
		: m_fineData(fineGridFct, fineCmp),
		  m_fineTopLevel(fineGridFct.dof_distribution()->grid_level().level()),
		  m_coarseData(coarseGridFct, coarseCmp),
		  m_coarseTopLevel(coarseGridFct.dof_distribution()->grid_level().level()),
		  m_spMG(fineGridFct.domain()->grid())
		{
			if(m_fineTopLevel < m_coarseTopLevel)
				UG_THROW("H1DiffIntegrand: fine and top level inverted.");

			if(fineGridFct.domain().get() !=
					coarseGridFct.domain().get())
				UG_THROW("H1DiffIntegrand: grid functions defined on different domains.");
		};

	/// DTOR
		virtual ~H1DistIntegrand() {};

	///	sets subset
		virtual void set_subset(int si)
		{
			if(!m_fineData.is_def_in_subset(si))
				UG_THROW("H1DiffIntegrand: Grid function component"
						<<m_fineData.fct()<<" not defined on subset "<<si);
			if(!m_coarseData.is_def_in_subset(si))
				UG_THROW("H1DiffIntegrand: Grid function component"
						<<m_coarseData.fct()<<" not defined on subset "<<si);
			IIntegrand<number, worldDim>::set_subset(si);
		}

	/// \copydoc IIntegrand::values
		template <int elemDim>
		void evaluate(number vValue[],
		              const MathVector<worldDim> vFineGlobIP[],
		              GridObject* pFineElem,
		              const MathVector<worldDim> vCornerCoords[],
		              const MathVector<elemDim> vFineLocIP[],
		              const MathMatrix<elemDim, worldDim> vJT[],
		              const size_t numIP)
		{
			typedef typename TGridFunction::template dim_traits<elemDim>::grid_base_object Element;

			//	get coarse element
			GridObject* pCoarseElem = pFineElem;
			if(m_coarseTopLevel < m_fineTopLevel){
				int parentLevel = m_spMG->get_level(pCoarseElem);
				while(parentLevel > m_coarseTopLevel){
					pCoarseElem = m_spMG->get_parent(pCoarseElem);
					parentLevel = m_spMG->get_level(pCoarseElem);
				}
			}

		//	get reference object id (i.e. Triangle, Quadrilateral, Tetrahedron, ...)
			const ReferenceObjectID fineROID = pFineElem->reference_object_id();
			const ReferenceObjectID coarseROID = pCoarseElem->reference_object_id();

		//	get corner coordinates
			std::vector<MathVector<worldDim> > vCornerCoarse;
			CollectCornerCoordinates(vCornerCoarse, *static_cast<Element*>(pCoarseElem), *m_coarseData.domain());

		//	get Reference Mapping
			DimReferenceMapping<elemDim, worldDim>& map
				= ReferenceMappingProvider::get<elemDim, worldDim>(coarseROID, vCornerCoarse);

			std::vector<MathVector<elemDim> > vCoarseLocIP;
			vCoarseLocIP.resize(numIP);
			for(size_t ip = 0; ip < vCoarseLocIP.size(); ++ip) VecSet(vCoarseLocIP[ip], 0.0);
			map.global_to_local(&vCoarseLocIP[0], vFineGlobIP, numIP);

			try{
		//	get trial space
			const LocalShapeFunctionSet<elemDim>& rFineLSFS =
					LocalFiniteElementProvider::get<elemDim>(fineROID, m_fineData.id());
			const LocalShapeFunctionSet<elemDim>& rCoarseLSFS =
					LocalFiniteElementProvider::get<elemDim>(coarseROID, m_coarseData.id());

		//	get multiindices of element
			std::vector<DoFIndex> vFineMI, vCoarseMI;
			m_fineData.dof_indices(pFineElem, vFineMI);
			m_coarseData.dof_indices(pCoarseElem, vCoarseMI);

			std::vector<MathVector<elemDim> > vFineLocGradient(vFineMI.size());
			std::vector<MathVector<elemDim> > vCoarseLocGradient(vCoarseMI.size());

		//	loop all integration points
			for(size_t ip = 0; ip < numIP; ++ip)
			{
			//	compute shape gradients at ip
				rFineLSFS.grads(&vFineLocGradient[0], vFineLocIP[ip]);
				rCoarseLSFS.grads(&vCoarseLocGradient[0], vCoarseLocIP[ip]);

			// 	compute approximated solution at integration point
				number fineSolIP = 0.0;
				MathVector<elemDim> fineLocTmp; VecSet(fineLocTmp, 0.0);
				for(size_t sh = 0; sh < vFineMI.size(); ++sh)
				{
					const number val = DoFRef(m_fineData.grid_function(), vFineMI[sh]);
					fineSolIP += val * rFineLSFS.shape(sh, vFineLocIP[ip]);
					VecScaleAppend(fineLocTmp, val, vFineLocGradient[sh]);
				}
				number coarseSolIP = 0.0;
				MathVector<elemDim> coarseLocTmp; VecSet(coarseLocTmp, 0.0);
				for(size_t sh = 0; sh < vCoarseMI.size(); ++sh)
				{
					const number val = DoFRef(m_coarseData.grid_function(), vCoarseMI[sh]);
					coarseSolIP += val * rCoarseLSFS.shape(sh, vCoarseLocIP[ip]);
					VecScaleAppend(coarseLocTmp, val, vCoarseLocGradient[sh]);
				}

			//	compute global gradient
				MathVector<worldDim> fineGradIP;
				MathMatrix<worldDim, elemDim> fineJTInv;
				Inverse(fineJTInv, vJT[ip]);
				MatVecMult(fineGradIP, fineJTInv, fineLocTmp);

			//	compute global gradient
				MathVector<worldDim> coarseGradIP;
				MathMatrix<worldDim, elemDim> coarseJTInv;
				map.jacobian_transposed_inverse(coarseJTInv, vCoarseLocIP[ip]);
				MatVecMult(coarseGradIP, coarseJTInv, coarseLocTmp);

			//	get squared of difference
				vValue[ip] = (coarseSolIP - fineSolIP);
				vValue[ip] *= vValue[ip];
				vValue[ip] += VecDistanceSq(coarseGradIP, fineGradIP);
			}

			}
			UG_CATCH_THROW("H1DiffIntegrand::evaluate: trial space missing.");
		};
};



template <typename TGridFunction>
number H1Distance2(TGridFunction& spGridFct1, const char* cmp1,
               TGridFunction& spGridFct2, const char* cmp2,
               int quadOrder, const char* subsets=NULL)
{
	return GridFunctionDistance2<H1DistIntegrand<TGridFunction>, TGridFunction>
		(spGridFct1, cmp1, spGridFct2, cmp2, quadOrder, subsets);
}


template <typename TGridFunction>
number H1Distance(TGridFunction& spGridFct1, const char* cmp1,
               TGridFunction& spGridFct2, const char* cmp2,
               int quadOrder, const char* subsets=NULL)
{
	return sqrt(H1Distance2(spGridFct1, cmp1, spGridFct2, cmp2, quadOrder, subsets));
}
/// for lua shell
template <typename TGridFunction>
number H1Error(SmartPtr<TGridFunction> spGridFct1, const char* cmp1,
               SmartPtr<TGridFunction> spGridFct2, const char* cmp2,
               int quadOrder, const char* subsets)
{
	return H1Distance(*spGridFct1, cmp1, *spGridFct2, cmp2, quadOrder, subsets);
}

/// for lua shell
template <typename TGridFunction>
number H1Error(SmartPtr<TGridFunction> spGridFct1, const char* cmp1,
               SmartPtr<TGridFunction> spGridFct2, const char* cmp2,
               int quadOrder)
{
	return H1Distance(*spGridFct1, cmp1, *spGridFct2, cmp2, quadOrder, NULL);
}


////////////////////////////////////////////////////////////////////////////////
// Standard Integrand
////////////////////////////////////////////////////////////////////////////////

template <typename TGridFunction>
class StdFuncIntegrand
	: public StdIntegrand<number, TGridFunction::dim, StdFuncIntegrand<TGridFunction> >
{
	public:
	//	world dimension of grid function
		static const int worldDim = TGridFunction::dim;

	private:
	// grid function
		TGridFunction* m_pGridFct;

	//	component of function
		const size_t m_fct;

	public:
	/// constructor
		StdFuncIntegrand(TGridFunction* pGridFct, size_t cmp)
		: m_pGridFct(pGridFct), m_fct(cmp)
		{};

		virtual ~StdFuncIntegrand(){}

	///	sets subset
		virtual void set_subset(int si)
		{
			if(!m_pGridFct->is_def_in_subset(m_fct, si))
				UG_THROW("L2ErrorIntegrand: Grid function component"
						<<m_fct<<" not defined on subset "<<si);
			IIntegrand<number, worldDim>::set_subset(si);
		}

	/// \copydoc IIntegrand::values
		template <int elemDim>
		void evaluate(number vValue[],
		              const MathVector<worldDim> vGlobIP[],
		              GridObject* pElem,
		              const MathVector<worldDim> vCornerCoords[],
		              const MathVector<elemDim> vLocIP[],
		              const MathMatrix<elemDim, worldDim> vJT[],
		              const size_t numIP)
		{
		//	get reference object id (i.e. Triangle, Quadrilateral, Tetrahedron, ...)
			ReferenceObjectID roid = (ReferenceObjectID) pElem->reference_object_id();

			const LFEID m_id = m_pGridFct->local_finite_element_id(m_fct);

			try{
		//	get trial space
			const LocalShapeFunctionSet<elemDim>& rTrialSpace =
							LocalFiniteElementProvider::get<elemDim>(roid, m_id);

		//	number of dofs on element
			const size_t num_sh = rTrialSpace.num_sh();

		//	get multiindices of element

			std::vector<DoFIndex> ind;  // 	aux. index array
			m_pGridFct->dof_indices(pElem, m_fct, ind);

		//	check multi indices
			if(ind.size() != num_sh)
				UG_THROW("StdFuncIntegrand::evaluate: Wrong number of"
						" multi indices.");

		//	loop all integration points
			for(size_t ip = 0; ip < numIP; ++ip)
			{

			// 	compute approximated solution at integration point
				number approxSolIP = 0.0;
				for(size_t sh = 0; sh < num_sh; ++sh)
				{
					//	get value at shape point (e.g. corner for P1 fct)
					//	and add shape fct at ip * value at shape
					const number valSH = DoFRef((*m_pGridFct), ind[sh]);
					approxSolIP += valSH * rTrialSpace.shape(sh, vLocIP[ip]);
				}

				//	get function value at ip
				vValue[ip] = approxSolIP;

			}

			}
			UG_CATCH_THROW("StdFuncIntegrand::evaluate: trial space missing.");
		};
};


template <typename TGridFunction>
number StdFuncIntegralOnVertex(TGridFunction& gridFct,
							   size_t fct,
							   int si)
{
//	integrate elements of subset
	typedef typename TGridFunction::template dim_traits<0>::grid_base_object grid_base_object;
	typedef typename TGridFunction::template dim_traits<0>::const_iterator const_iterator;

	//	reset the result
	number integral = 0;

//	note: this iterator is for the base elements, e.g. Face and not
//			for the special type, e.g. Triangle, Quadrilateral
	const_iterator iter = gridFct.template begin<grid_base_object>(si);
	const_iterator iterEnd = gridFct.template end<grid_base_object>(si);

// 	iterate over all elements
	for(; iter != iterEnd; ++iter)
	{
	//	get element
		grid_base_object* pElem = *iter;

		std::vector<DoFIndex> ind;  // 	aux. index array
		gridFct.dof_indices(pElem, fct, ind);

	// 	compute approximated solution at integration point
		number value = 0.0;
		for(size_t sh = 0; sh < ind.size(); ++sh)
		{
			value += DoFRef(gridFct, ind[sh]);
		}

	//	add to global sum
		integral += value;

	} // end elem

//	return the summed integral contributions of all elements
	return integral;
}

template <typename TGridFunction>
number StdFuncIntegralOnVertex(SmartPtr<TGridFunction> spGridFct, size_t fct, int si)
{ return StdFuncIntegralOnVertex(*spGridFct, fct, si); }


template <typename TGridFunction>
number Integral(TGridFunction& gridFct, const char* cmp,
                const char* subsets, int quadOrder)
{
//	get function id of name
	const size_t fct = gridFct.fct_id_by_name(cmp);

//	check that function exists
	if(fct >= gridFct.num_fct())
		UG_THROW("L2Norm: Function space does not contain"
				" a function with name " << cmp << ".");

//	read subsets
	SubsetGroup ssGrp(gridFct.domain()->subset_handler());
	if(subsets != NULL)
	{
		ssGrp.add(TokenizeString(subsets));
		if(!SameDimensionsInAllSubsets(ssGrp))
			UG_THROW("IntegrateSubsets: Subsets '"<<subsets<<"' do not have same dimension."
					 "Can not integrate on subsets of different dimensions.");
	}
	else
	{
	//	add all subsets and remove lower dim subsets afterwards
		ssGrp.add_all();
		RemoveLowerDimSubsets(ssGrp);
	}

	// \TODO: This should be generalite in IntegrateSubset instead of doing it directly here
	bool bOnlyVertex = true;
	for(size_t s = 0; s < ssGrp.size(); ++s)
		if(ssGrp.dim(s) != 0) bOnlyVertex = false;

	if(bOnlyVertex)
	{
		number value = 0;
		for(size_t s = 0; s < ssGrp.size(); ++s)
			value += StdFuncIntegralOnVertex(gridFct, fct, ssGrp[s]);

#ifdef UG_PARALLEL
	// sum over processes
	if(pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		number local = value;
		com.allreduce(&local, &value, 1, PCL_DT_DOUBLE, PCL_RO_SUM);
	}
#endif
		return value;
	}

	StdFuncIntegrand<TGridFunction> integrand(&gridFct, fct);
	return IntegrateSubsets(integrand, gridFct, subsets, quadOrder);
}

template <typename TGridFunction>
number Integral(SmartPtr<TGridFunction> spGridFct, const char* cmp,
                const char* subsets, int quadOrder)
{ return Integral(*spGridFct, cmp, subsets, quadOrder); }


template <typename TGridFunction>
number Integral(SmartPtr<TGridFunction> spGridFct, const char* cmp,
                const char* subsets)
{
	return Integral(spGridFct, cmp, subsets, 1);
}
template <typename TGridFunction>
number Integral(SmartPtr<TGridFunction> spGridFct, const char* cmp)
{
	return Integral(spGridFct, cmp, NULL, 1);
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Generic Boundary Integration Routine
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <int WorldDim, int dim, typename TConstIterator>
number IntegralNormalComponentOnManifoldUsingFV1Geom(TConstIterator iterBegin,
                                       TConstIterator iterEnd,
                                       typename domain_traits<WorldDim>::position_accessor_type& aaPos,
                                       const ISubsetHandler* ish,
                                       IIntegrand<MathVector<WorldDim>, WorldDim>& integrand,
                                       const SubsetGroup& bndSSGrp)
{
//	reset the result
	number integral = 0;

//	note: this iterator is for the base elements, e.g. Face and not
//			for the special type, e.g. Triangle, Quadrilateral
	TConstIterator iter = iterBegin;

//	this is the base element type (e.g. Face). This is the type when the
//	iterators above are dereferenciated.
	typedef typename domain_traits<dim>::grid_base_object grid_base_object;

//	create a FV1 Geometry
	DimFV1Geometry<dim> geo;

//	specify, which subsets are boundary
	for(size_t s = 0; s < bndSSGrp.size(); ++s)
	{
	//	get subset index
		const int bndSubset = bndSSGrp[s];

	//	request this subset index as boundary subset. This will force the
	//	creation of boundary subsets when calling geo.update
		geo.add_boundary_subset(bndSubset);
	}

//	vector of corner coordinates of element corners (to be filled for each elem)
	std::vector<MathVector<WorldDim> > vCorner;

// 	iterate over all elements
	for(; iter != iterEnd; ++iter)
	{
	//	get element
		grid_base_object* pElem = *iter;

	//	get all corner coordinates
		CollectCornerCoordinates(vCorner, *pElem, aaPos, true);

	//	compute bf and grads at bip for element
		try{
			geo.update(pElem, &vCorner[0], ish);
		}
		UG_CATCH_THROW("IntegralNormalComponentOnManifold: "
					"Cannot update Finite Volume Geometry.");

	//	specify, which subsets are boundary
		for(size_t s = 0; s < bndSSGrp.size(); ++s)
		{
		//	get subset index
			const int bndSubset = bndSSGrp[s];

		//	get all bf of this subset
			typedef typename DimFV1Geometry<dim>::BF BF;
			const std::vector<BF>& vBF = geo.bf(bndSubset);

		//	loop boundary faces
			for(size_t b = 0; b < vBF.size(); ++b)
			{
			//	get bf
				const BF& bf = vBF[b];

			//	value
				MathVector<WorldDim> value;

			//	jacobian
				MathMatrix<dim, WorldDim> JT;
				Inverse(JT, bf.JTInv());

			//	compute integrand values at integration points
				try
				{
					integrand.values(&value, &(bf.global_ip()),
									 pElem, &vCorner[0], &(bf.local_ip()),
									 &(JT),
									 1);
				}
				UG_CATCH_THROW("IntegralNormalComponentOnManifold: Unable to compute values of "
								"integrand at integration point.");

			//	sum up contribution (normal includes area)
				integral += VecDot(value, bf.normal());

			} // end bf
		} // end bnd subsets
	} // end elem

//	return the summed integral contributions of all elements
	return integral;
}


template <int WorldDim, int dim, typename TConstIterator>
number IntegralNormalComponentOnManifoldGeneral(
		TConstIterator iterBegin,
		TConstIterator iterEnd,
		typename domain_traits<WorldDim>::position_accessor_type& aaPos,
		const ISubsetHandler* ish,
		IIntegrand<MathVector<WorldDim>, WorldDim>& integrand,
		const SubsetGroup& bndSSGrp,
		int quadOrder,
		Grid& grid)
{
//	reset the result
	number integral = 0;

//	note: this iterator is for the base elements, e.g. Face and not
//			for the special type, e.g. Triangle, Quadrilateral
	TConstIterator iter = iterBegin;

//	this is the base element type (e.g. Face). This is the type when the
//	iterators above are dereferenciated.
	typedef typename domain_traits<dim>::element_type Element;
	typedef typename domain_traits<dim>::side_type Side;

//	vector of corner coordinates of element corners (to be filled for each elem)
	std::vector<MathVector<WorldDim> > vCorner;
	std::vector<int> vSubsetIndex;

// 	iterate over all elements
	for(; iter != iterEnd; ++iter)
	{
	//	get element
		Element* pElem = *iter;

	//	get all corner coordinates
		CollectCornerCoordinates(vCorner, *pElem, aaPos, true);

	//	get reference object id
		const ReferenceObjectID elemRoid = pElem->reference_object_id();

	//	get sides
		typename Grid::traits<Side>::secure_container vSide;
		grid.associated_elements_sorted(vSide, pElem);
		vSubsetIndex.resize(vSide.size());
		for(size_t i = 0; i < vSide.size(); ++i)
			vSubsetIndex[i] = ish->get_subset_index(vSide[i]);

		DimReferenceMapping<dim, WorldDim>& rMapping
			= ReferenceMappingProvider::get<dim, WorldDim>(elemRoid, vCorner);

		const DimReferenceElement<dim>& rRefElem
			= ReferenceElementProvider::get<dim>(elemRoid);

	//	loop sub elements
		for(size_t side = 0; side < vSide.size(); ++side)
		{
		//	check if side used
			if(!bndSSGrp.contains(vSubsetIndex[side])) continue;

		//	get side
			Side* pSide = vSide[side];

			std::vector<MathVector<WorldDim> > vSideCorner(rRefElem.num(dim-1, side, 0));
			std::vector<MathVector<dim> > vLocalSideCorner(rRefElem.num(dim-1, side, 0));
			for(size_t co = 0; co < vSideCorner.size(); ++co){
				vSideCorner[co] = vCorner[rRefElem.id(dim-1, side, 0, co)];
				vLocalSideCorner[co] = rRefElem.corner(rRefElem.id(dim-1, side, 0, co));
			}

		//	side quad rule
			const ReferenceObjectID sideRoid = pSide->reference_object_id();
			const QuadratureRule<dim-1>& rSideQuadRule
					= QuadratureRuleProvider<dim-1>::get(sideRoid, quadOrder);

		// 	normal
			MathVector<WorldDim> Normal;
			ElementNormal<WorldDim>(sideRoid, Normal, &vSideCorner[0]);

		//	quadrature points
			const number* vWeight = rSideQuadRule.weights();
			const size_t nip = rSideQuadRule.size();
			std::vector<MathVector<dim> > vLocalIP(nip);
			std::vector<MathVector<dim> > vGlobalIP(nip);

			DimReferenceMapping<dim-1, dim>& map
				= ReferenceMappingProvider::get<dim-1, dim>(sideRoid, vLocalSideCorner);

			for(size_t ip = 0; ip < nip; ++ip)
				map.local_to_global(vLocalIP[ip], rSideQuadRule.point(ip));

			for(size_t ip = 0; ip < nip; ++ip)
				rMapping.local_to_global(vGlobalIP[ip], vLocalIP[ip]);

		//	compute transformation matrices
			std::vector<MathMatrix<dim-1, WorldDim> > vJT(nip);
			map.jacobian_transposed(&(vJT[0]), rSideQuadRule.points(), nip);

			std::vector<MathMatrix<dim, WorldDim> > vElemJT(nip);
			rMapping.jacobian_transposed(&(vElemJT[0]), &vLocalIP[0], nip);

			std::vector<MathVector<WorldDim> > vValue(nip);

		//	compute integrand values at integration points
			try
			{
				integrand.values(&vValue[0], &vGlobalIP[0],
								 pElem, &vCorner[0], &vLocalIP[0],
								 &(vElemJT[0]),
								 nip);
			}
			UG_CATCH_THROW("IntegralNormalComponentOnManifold: Unable to compute values of "
							"integrand at integration point.");

		//	loop integration points
			for(size_t ip = 0; ip < nip; ++ip)
			{
			//	get quadrature weight
				const number weightIP = vWeight[ip];

			//	get determinate of mapping
				const number det = SqrtGramDeterminant(vJT[ip]);

			//	add contribution of integration point
				integral +=  VecDot(vValue[ip], Normal) * weightIP * det;
			}
		} // end bf
	} // end elem

//	return the summed integral contributions of all elements
	return integral;
}

template <typename TGridFunction, int dim>
number IntegralNormalComponentOnManifoldSubset(SmartPtr<IIntegrand<MathVector<TGridFunction::dim>, TGridFunction::dim> > spIntegrand,
                                 SmartPtr<TGridFunction> spGridFct,
                                 int si, const SubsetGroup& bndSSGrp, int quadOrder)
{
//	integrate elements of subset
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object grid_base_object;
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator const_iterator;
	static const int WorldDim = TGridFunction::dim;

	spIntegrand->set_subset(si);

	if(quadOrder == 1)
		return IntegralNormalComponentOnManifoldUsingFV1Geom<WorldDim,dim,const_iterator>
					(spGridFct->template begin<grid_base_object>(si),
	                 spGridFct->template end<grid_base_object>(si),
	                 spGridFct->domain()->position_accessor(),
	                 spGridFct->domain()->subset_handler().get(),
	                 *spIntegrand, bndSSGrp);
	else{
		UG_LOG(" #### IntegralNormalComponentOnManifoldSubset ####:\n")
		return IntegralNormalComponentOnManifoldGeneral<WorldDim,dim,const_iterator>
					(spGridFct->template begin<grid_base_object>(si),
	                 spGridFct->template end<grid_base_object>(si),
	                 spGridFct->domain()->position_accessor(),
	                 spGridFct->domain()->subset_handler().get(),
	                 *spIntegrand, bndSSGrp, quadOrder, *spGridFct->domain()->grid());
	}
}

template <typename TGridFunction>
number IntegralNormalComponentOnManifoldSubsets(
		SmartPtr<IIntegrand<MathVector<TGridFunction::dim>, TGridFunction::dim> > spIntegrand,
		SmartPtr<TGridFunction> spGridFct,
		const char* BndSubsets, const char* InnerSubsets,
		int quadOrder)
{
//	world dimensions
	static const int dim = TGridFunction::dim;

//	read subsets
	SubsetGroup innerSSGrp(spGridFct->domain()->subset_handler());
	if(InnerSubsets != NULL)
	{
		innerSSGrp.add(TokenizeString(InnerSubsets));
		if(!SameDimensionsInAllSubsets(innerSSGrp))
			UG_THROW("IntegralNormalComponentOnManifold: Subsets '"<<InnerSubsets<<"' do not have same dimension."
			         "Can not integrate on subsets of different dimensions.");
	}
	else
	{
	//	add all subsets and remove lower dim subsets afterwards
		innerSSGrp.add_all();
		RemoveLowerDimSubsets(innerSSGrp);
	}

//	read subsets
	SubsetGroup bndSSGrp(spGridFct->domain()->subset_handler());
	if(BndSubsets != NULL)
		bndSSGrp.add(TokenizeString(BndSubsets));
	else
		UG_THROW("IntegralNormalComponentOnManifold: No boundary subsets passed.");

//	reset value
	number value = 0;

//	loop subsets
	for(size_t i = 0; i < innerSSGrp.size(); ++i)
	{
	//	get subset index
		const int si = innerSSGrp[i];

	//	check dimension
		if(innerSSGrp.dim(i) != dim && innerSSGrp.dim(i) != DIM_SUBSET_EMPTY_GRID)
			UG_THROW("IntegralNormalComponentOnManifold: Dimension of inner subset is "<<
			         innerSSGrp.dim(i)<<", but only World Dimension "<<dim<<
			         " subsets can be used for inner subsets.");

	//	integrate elements of subset
		switch(innerSSGrp.dim(i))
		{
			case DIM_SUBSET_EMPTY_GRID: break;
			case dim: value += IntegralNormalComponentOnManifoldSubset<TGridFunction, dim>(spIntegrand, spGridFct, si, bndSSGrp, quadOrder); break;
			default: UG_THROW("IntegralNormalComponentOnManifold: Dimension "<<innerSSGrp.dim(i)<<" not supported. "
			                  " World dimension is "<<dim<<".");
		}
	}

#ifdef UG_PARALLEL
	// sum over processes
	if(pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		number local = value;
		com.allreduce(&local, &value, 1, PCL_DT_DOUBLE, PCL_RO_SUM);
	}
#endif

//	return the result
	return value;
}

template <typename TGridFunction>
number IntegralNormalComponentOnManifold(
		SmartPtr<UserData<MathVector<TGridFunction::dim>, TGridFunction::dim> > spData,
		SmartPtr<TGridFunction> spGridFct,
		const char* BndSubset, const char* InnerSubset,
		number time,
		int quadOrder)
{
	SmartPtr<IIntegrand<MathVector<TGridFunction::dim>, TGridFunction::dim> > spIntegrand
		= make_sp(new UserDataIntegrand<MathVector<TGridFunction::dim>, TGridFunction>(spData, &(*spGridFct), time));

	return IntegralNormalComponentOnManifoldSubsets(spIntegrand, spGridFct, BndSubset, InnerSubset, quadOrder);
}

template <typename TGridFunction>
number IntegralNormalComponentOnManifold(
		SmartPtr<UserData<MathVector<TGridFunction::dim>, TGridFunction::dim> > spData,
		SmartPtr<TGridFunction> spGridFct,
		const char* BndSubset, const char* InnerSubset,
		number time)
{return IntegralNormalComponentOnManifold(spData, spGridFct, BndSubset, InnerSubset, time, 1);}

template <typename TGridFunction>
number IntegralNormalComponentOnManifold(
		SmartPtr<UserData<MathVector<TGridFunction::dim>, TGridFunction::dim> > spData,
		SmartPtr<TGridFunction> spGridFct,
		const char* BndSubset,
		number time)
{return IntegralNormalComponentOnManifold(spData, spGridFct, BndSubset, NULL, time, 1);}

template <typename TGridFunction>
number IntegralNormalComponentOnManifold(
		SmartPtr<UserData<MathVector<TGridFunction::dim>, TGridFunction::dim> > spData,
		SmartPtr<TGridFunction> spGridFct,
		const char* BndSubset, const char* InnerSubset)
{return IntegralNormalComponentOnManifold(spData, spGridFct, BndSubset, InnerSubset, 0.0, 1);}

template <typename TGridFunction>
number IntegralNormalComponentOnManifold(
		SmartPtr<UserData<MathVector<TGridFunction::dim>, TGridFunction::dim> > spData,
		SmartPtr<TGridFunction> spGridFct,
		const char* BndSubset)
{return IntegralNormalComponentOnManifold(spData, spGridFct, BndSubset, NULL, 0.0, 1);}

///////////////
// lua data
///////////////

#ifdef UG_FOR_LUA
template <typename TGridFunction>
number IntegralNormalComponentOnManifold(
		const char* luaFct, SmartPtr<TGridFunction> spGridFct,
		const char* BndSubset, const char* InnerSubset,
		number time, int quadOrder)
{
	static const int dim = TGridFunction::dim;
	SmartPtr<UserData<MathVector<dim>, dim> > sp =
			LuaUserDataFactory<MathVector<dim>, dim>::create(luaFct);
	return IntegralNormalComponentOnManifold(sp, spGridFct, BndSubset, InnerSubset, time, quadOrder);
}

template <typename TGridFunction>
number IntegralNormalComponentOnManifold(
		const char* luaFct, SmartPtr<TGridFunction> spGridFct,
		const char* BndSubset, const char* InnerSubset,
		number time)
{return IntegralNormalComponentOnManifold(luaFct, spGridFct, BndSubset, InnerSubset, time, 1);}

template <typename TGridFunction>
number IntegralNormalComponentOnManifold(
		const char* luaFct, SmartPtr<TGridFunction> spGridFct,
		const char* BndSubset,
		number time)
{return IntegralNormalComponentOnManifold(luaFct, spGridFct, BndSubset, NULL, time, 1);}

template <typename TGridFunction>
number IntegralNormalComponentOnManifold(
		const char* luaFct, SmartPtr<TGridFunction> spGridFct,
		const char* BndSubset, const char* InnerSubset)
{return IntegralNormalComponentOnManifold(luaFct, spGridFct, BndSubset, InnerSubset, 0.0, 1);}

template <typename TGridFunction>
number IntegralNormalComponentOnManifold(
		const char* luaFct, SmartPtr<TGridFunction> spGridFct,
		const char* BndSubset)
{return IntegralNormalComponentOnManifold(luaFct, spGridFct, BndSubset, NULL, 0.0, 1);}
#endif

////////////////////////////////////////////////////////////////////////////////
// Boundary Integral
////////////////////////////////////////////////////////////////////////////////

/// Integrates the Flux of a component over some boundary subsets
/**
 * This function integrates \f$ - \nabla c \cdot \vec{n} \f$ over some selected
 * boundary subsets. In order to compute the gradient of a function a world-
 * dimension element must be given and, thus, some "inner" subsets have to be
 * specified to indicate the full-dimensional subsets. The integral sum is then
 * taken over all boundary subset manifold geometric objects, that are part of
 * an element of the scheduled "inner" elements
 *
 * @param u				a grid function
 * @param cmp			the component, whose gradient should be integrated
 * @param BndSubset		a comma-separated string of symbolic boundary subset names
 * @param InnerSubset	a comma-separated string of symbolic inner subset names
 * @return	the integral
 */
template <typename TGridFunction>
number IntegrateNormalGradientOnManifold(TGridFunction& u, const char* cmp,
                                   const char* BndSubset, const char* InnerSubset = NULL)
{
//	get function id of name
	const size_t fct = u.fct_id_by_name(cmp);

//	check that function exists
	if(fct >= u.num_fct())
		UG_THROW("IntegrateNormalGradientOnManifold: Function space does not contain"
				" a function with name " << cmp << ".");

	if(u.local_finite_element_id(fct) != LFEID(LFEID::LAGRANGE, TGridFunction::dim, 1))
		UG_THROW("IntegrateNormalGradientOnManifold:"
				 "Only implemented for Lagrange P1 functions.");

//	read subsets
	SubsetGroup innerSSGrp(u.domain()->subset_handler());
	if(InnerSubset != NULL)
		innerSSGrp.add(TokenizeString(InnerSubset));
	else // add all if no subset specified
		innerSSGrp.add_all();

//	read bnd subsets
	SubsetGroup bndSSGrp(u.domain()->subset_handler());
	if(InnerSubset != NULL){
		bndSSGrp.add(TokenizeString(BndSubset));
	}
	else{
		UG_THROW("IntegrateNormalGradientOnManifold: No boundary subsets specified. Aborting.");
	}

//	reset value
	number value = 0;

//	loop subsets
	for(size_t i = 0; i < innerSSGrp.size(); ++i)
	{
	//	get subset index
		const int si = innerSSGrp[i];

	//	skip if function is not defined in subset
		if(!u.is_def_in_subset(fct, si)) continue;


		if (innerSSGrp.dim(i) != TGridFunction::dim)
			UG_THROW("IntegrateNormalGradientOnManifold: Element dimension does not match world dimension!");

	//	create integration kernel
		static const int dim = TGridFunction::dim;

	//	integrate elements of subset
		typedef typename domain_traits<dim>::grid_base_object grid_base_object;

	//	get iterators for all elems on subset
		typedef typename TGridFunction::template dim_traits<dim>::const_iterator const_iterator;
		const_iterator iter = u.template begin<grid_base_object>();
		const_iterator iterEnd = u.template end<grid_base_object>();

	//	create a FV1 Geometry
		DimFV1Geometry<dim> geo;

	//	specify, which subsets are boundary
		for(size_t s = 0; s < bndSSGrp.size(); ++s)
		{
		//	get subset index
			const int bndSubset = bndSSGrp[s];

		//	request this subset index as boundary subset. This will force the
		//	creation of boundary subsets when calling geo.update
			geo.add_boundary_subset(bndSubset);
		}

	//	vector of corner coordinates of element corners (to be filled for each elem)
		std::vector<MathVector<dim> > vCorner;

	//	loop elements of subset
		for( ; iter != iterEnd; ++iter)
		{
		//	get element
			grid_base_object* elem = *iter;

		//	get all corner coordinates
			CollectCornerCoordinates(vCorner, *elem, u.domain()->position_accessor(), true);

		//	compute bf and grads at bip for element
			try{
				geo.update(elem, &vCorner[0], u.domain()->subset_handler().get());
			}
			UG_CATCH_THROW("IntegrateNormalGradientOnManifold: "
						"Cannot update Finite Volume Geometry.");

		//	get fct multi-indeces of element
			std::vector<DoFIndex> ind;
			u.dof_indices(elem, fct, ind);

		//	specify, which subsets are boundary
			for(size_t s = 0; s < bndSSGrp.size(); ++s)
			{
			//	get subset index
				const int bndSubset = bndSSGrp[s];

			//	get all bf of this subset
				typedef typename DimFV1Geometry<dim>::BF BF;
				const std::vector<BF>& vBF = geo.bf(bndSubset);

			//	loop boundary faces
				for(size_t b = 0; b < vBF.size(); ++b)
				{
				//	get bf
					const BF& bf = vBF[b];

				//	get normal on bf
					const MathVector<dim>& normal = bf.normal();

				//	check multi indices
					UG_ASSERT(ind.size() == bf.num_sh(),
					          "IntegrateNormalGradientOnManifold::values: Wrong number of"
							  " multi indices, ind: "<<ind.size() << ", bf.num_sh: "
							  << bf.num_sh());

				// 	compute gradient of solution at bip
					MathVector<dim> grad; VecSet(grad, 0.0);
					for(size_t sh = 0; sh < bf.num_sh(); ++sh)
					{
						const number fctVal = DoFRef(u, ind[sh]);

						VecScaleAdd(grad, 1.0, grad, fctVal, bf.global_grad(sh));
					}

				//	sum up contributions
					value -= VecDot(grad, normal);
				}
			}
		}
	}

#ifdef UG_PARALLEL
	// sum over processes
	if(pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		number local = value;
		com.allreduce(&local, &value, 1, PCL_DT_DOUBLE, PCL_RO_SUM);
	}
#endif

//	return the result
	return value;
}

} // namespace ug

#endif /*__H__UG__LIB_DISC__FUNCTION_SPACES__INTEGRATE__*/
