/*
 * integrate.h
 *
 *  Created on: 04.04.2011
 *      Author: avogel, anaegel (some parts)
 */

#ifndef __H__UG__LIB_DISC__FUNCTION_SPACES__INTEGRATE__
#define __H__UG__LIB_DISC__FUNCTION_SPACES__INTEGRATE__

#include <cmath>

#include "common/common.h"

#include "lib_disc/common/subset_group.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/quadrature/quadrature.h"
#include "lib_disc/local_finite_element/local_shape_function_set.h"
#include "lib_disc/spatial_disc/disc_util/finite_volume_geometry.h"
#include "lib_disc/spatial_disc/user_data/user_data.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"
#include <boost/function.hpp>

#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// Integrand Interface
////////////////////////////////////////////////////////////////////////////////

/// Abstract integrand interface
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
		                    GeometricObject* pElem,
	                        const MathVector<worldDim> vCornerCoords[],
		                    const MathVector<1> vLocIP[],
		                    const MathMatrix<1, worldDim> vJT[],
		                    const size_t numIP) = 0;
		virtual void values(TData vValue[],
		                    const MathVector<worldDim> vGlobIP[],
		                    GeometricObject* pElem,
	                        const MathVector<worldDim> vCornerCoords[],
		                    const MathVector<2> vLocIP[],
		                    const MathMatrix<2, worldDim> vJT[],
		                    const size_t numIP) = 0;
		virtual void values(TData vValue[],
		                    const MathVector<worldDim> vGlobIP[],
		                    GeometricObject* pElem,
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

/// Abstract integrand interface
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
		                    GeometricObject* pElem,
	                        const MathVector<worldDim> vCornerCoords[],
		                    const MathVector<1> vLocIP[],
		                    const MathMatrix<1, worldDim> vJT[],
		                    const size_t numIP)
		{
			getImpl().template evaluate<1>(vValue,vGlobIP,pElem,vCornerCoords,vLocIP,vJT,numIP);
		}
		virtual void values(TData vValue[],
		                    const MathVector<worldDim> vGlobIP[],
		                    GeometricObject* pElem,
	                        const MathVector<worldDim> vCornerCoords[],
		                    const MathVector<2> vLocIP[],
		                    const MathMatrix<2, worldDim> vJT[],
		                    const size_t numIP)
		{
			getImpl().template evaluate<2>(vValue,vGlobIP,pElem,vCornerCoords,vLocIP,vJT,numIP);
		}
		virtual void values(TData vValue[],
		                    const MathVector<worldDim> vGlobIP[],
		                    GeometricObject* pElem,
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
 * \returns			value of the integral
 */
template <int WorldDim, int dim, typename TConstIterator>
number Integrate(TConstIterator iterBegin,
                 TConstIterator iterEnd,
                 typename domain_traits<WorldDim>::position_accessor_type& aaPos,
                 IIntegrand<number, WorldDim>& integrand,
                 int quadOrder)
{
//	reset the result
	number integral = 0;

//	note: this iterator is for the base elements, e.g. Face and not
//			for the special type, e.g. Triangle, Quadrilateral
	TConstIterator iter = iterBegin;

//	this is the base element type (e.g. Face). This is the type when the
//	iterators above are dereferenciated.
	typedef typename domain_traits<dim>::geometric_base_object geometric_base_object;

// 	iterate over all elements
	for(; iter != iterEnd; ++iter)
	{
	//	get element
		geometric_base_object* pElem = *iter;

	//	get reference object id (i.e. Triangle, Quadrilateral, Tetrahedron, ...)
		ReferenceObjectID roid = (ReferenceObjectID) pElem->reference_object_id();

	//	get quadrature Rule for reference object id and order
		try{
		const QuadratureRule<dim>& rQuadRule
					= QuadratureRuleProvider<dim>::get_rule(roid, quadOrder);

	//	number of integration points
		const size_t numIP = rQuadRule.size();

	//	get reference element mapping by reference object id
		try{
		DimReferenceMapping<dim, WorldDim>& mapping
							= ReferenceMappingProvider::get<dim, WorldDim>(roid);

	//	get all corner coordinates
		std::vector<MathVector<WorldDim> > vCorner;
		CollectCornerCoordinates(vCorner, *pElem, aaPos, true);

	//	update the reference mapping for the corners
		mapping.update(&vCorner[0]);

	//	compute global integration points
		std::vector<MathVector<WorldDim> > vGlobIP(numIP);
		mapping.local_to_global(&(vGlobIP[0]), rQuadRule.points(), numIP);

	//	compute transformation matrices
		std::vector<MathMatrix<dim, WorldDim> > vJT(numIP);
		mapping.jacobian_transposed(&(vJT[0]), rQuadRule.points(), numIP);

	//	compute integrand values at integration points
		std::vector<number> vValue(numIP);
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

		}catch(UG_ERROR_QuadratureRuleNotRegistered& ex){
			UG_THROW("SumValuesOnElems: " << ex.get_msg() << ".");
		}
		}catch(UGError_ReferenceMappingMissing& ex){
			UG_THROW("SumValuesOnElems: " << ex.get_msg() << ".");
		}
	} // end elem

//	return the summed integral contributions of all elements
	return integral;
}

template <typename TGridFunction, int dim>
number IntegrateSubset(SmartPtr<IIntegrand<number, TGridFunction::dim> > spIntegrand,
                       SmartPtr<TGridFunction> spGridFct,
                       int si, int quadOrder)
{
//	integrate elements of subset
	typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object geometric_base_object;
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator const_iterator;

	spIntegrand->set_subset(si);

	return Integrate<TGridFunction::dim,dim,const_iterator>
					(spGridFct->template begin<geometric_base_object>(si),
	                 spGridFct->template end<geometric_base_object>(si),
	                 spGridFct->domain()->position_accessor(),
	                 *spIntegrand,
	                 quadOrder);
}


template <typename TGridFunction>
number IntegrateSubsets(SmartPtr<IIntegrand<number, TGridFunction::dim> > spIntegrand,
                        SmartPtr<TGridFunction> spGridFct,
                        const char* subsets, int quadOrder)
{
//	world dimensions
	static const int dim = TGridFunction::dim;

//	read subsets
	SubsetGroup ssGrp(spGridFct->domain()->subset_handler());
	if(subsets != NULL)
	{
		ConvertStringToSubsetGroup(ssGrp, subsets);
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
			case 1: value += IntegrateSubset<TGridFunction, 1>(spIntegrand, spGridFct, si, quadOrder); break;
			case 2: value += IntegrateSubset<TGridFunction, 2>(spIntegrand, spGridFct, si, quadOrder); break;
			case 3: value += IntegrateSubset<TGridFunction, 3>(spIntegrand, spGridFct, si, quadOrder); break;
			default: UG_THROW("IntegrateSubsets: Dimension "<<ssGrp.dim(i)<<" not supported. "
			                  " World dimension is "<<dim<<".");
		}
		}
		UG_CATCH_THROW("IntegrateSubsets: Integration failed on subset "<<si);
	}

#ifdef UG_PARALLEL
	// sum over processes
	if(pcl::GetNumProcesses() > 1)
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
		SmartPtr<TGridFunction> m_spGridFct;

	//	time
		number m_time;

	public:
	/// constructor
		UserDataIntegrand(SmartPtr<UserData<TData, worldDim> > spData,
		                      SmartPtr<TGridFunction> spGridFct,
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
			if(m_spData->requires_grid_fct())
				UG_THROW("DirectUserDataIntegrand: Missing GridFunction, but "
						" data requires grid function.")
		};

	/// \copydoc IIntegrand::values
		template <int elemDim>
		void evaluate(TData vValue[],
		              const MathVector<worldDim> vGlobIP[],
		              GeometricObject* pElem,
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
					(*m_spData)(vValue, vGlobIP, m_time, this->m_si, u, pElem,
								vCornerCoords, vLocIP, numIP, &vJT[0]);
				}
				UG_CATCH_THROW("DirectUserDataIntegrand: Cannot evaluate data.");
			}
			else
			{
			//	compute data
				try{
					(*m_spData)(vValue, vGlobIP, m_time, this->m_si, numIP);
				}
				UG_CATCH_THROW("DirectUserDataIntegrand: Cannot evaluate data.");
			}

		};
};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Volume Integrand implementations
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

///////////////
// DirectUserData
///////////////

template <typename TGridFunction>
number Integral(SmartPtr<UserData<number, TGridFunction::dim> > spData,
                SmartPtr<TGridFunction> spGridFct,
                const char* subsets, number time,
                int quadOrder)
{
	SmartPtr<IIntegrand<number, TGridFunction::dim> > spIntegrand
		= CreateSmartPtr(new UserDataIntegrand<number, TGridFunction>(spData, spGridFct, time));

	return IntegrateSubsets(spIntegrand, spGridFct, subsets, quadOrder);
}

template <typename TGridFunction>
number Integral(SmartPtr<UserData<number, TGridFunction::dim> > spData, SmartPtr<TGridFunction> spGridFct,const char* subsets,number time)
{return Integral(spData, spGridFct, subsets, time, 1);}

template <typename TGridFunction>
number Integral(SmartPtr<UserData<number, TGridFunction::dim> > spData, SmartPtr<TGridFunction> spGridFct,number time)
{return Integral(spData, spGridFct, NULL, time, 1);}

template <typename TGridFunction>
number Integral(SmartPtr<UserData<number, TGridFunction::dim> > spData, SmartPtr<TGridFunction> spGridFct,const char* subsets)
{return Integral(spData, spGridFct, subsets, 0.0, 1);}

template <typename TGridFunction>
number Integral(SmartPtr<UserData<number, TGridFunction::dim> > spData, SmartPtr<TGridFunction> spGridFct)
{return Integral(spData, spGridFct, NULL, 0.0, 1);}

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
			CreateSmartPtr(new ConstUserNumber<dim>(val));
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
                int quadOrder)
{
	static const int dim = TGridFunction::dim;
	SmartPtr<UserData<number, dim> > sp =
			LuaUserDataFactory<number, dim>::create(luaFct);
	return Integral(sp, spGridFct, subsets, time, quadOrder);
}

template <typename TGridFunction>
number Integral(const char* luaFct, SmartPtr<TGridFunction> spGridFct,const char* subsets,number time)
{return Integral(luaFct, spGridFct, subsets, time, 1);}

template <typename TGridFunction>
number Integral(const char* luaFct, SmartPtr<TGridFunction> spGridFct,number time)
{return Integral(luaFct, spGridFct, NULL, time, 1);}

template <typename TGridFunction>
number Integral(const char* luaFct, SmartPtr<TGridFunction> spGridFct,const char* subsets)
{return Integral(luaFct, spGridFct, subsets, 0.0, 1);}

template <typename TGridFunction>
number Integral(const char* luaFct, SmartPtr<TGridFunction> spGridFct)
{return Integral(luaFct, spGridFct, NULL, 0.0, 1);}

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
	/// grid function
		SmartPtr<TGridFunction> m_spGridFct;

	///	component of function
		const size_t m_fct;

	///  exact solution
		SmartPtr<UserData<number, worldDim> > m_spExactSolution;

	///	time
		number m_time;

	public:
	/// constructor
		L2ErrorIntegrand(SmartPtr<UserData<number, worldDim> > spExactSol,
		                 SmartPtr<TGridFunction> gridFct, size_t cmp,
		                 number time)
		: m_spGridFct(gridFct), m_fct(cmp),
		  m_spExactSolution(spExactSol), m_time(time)
		{};

	///	sets subset
		virtual void set_subset(int si)
		{
			if(!m_spGridFct->is_def_in_subset(m_fct, si))
				UG_THROW("L2ErrorIntegrand: Grid function component"
						<<m_fct<<" not defined on subset "<<si);
			IIntegrand<number, worldDim>::set_subset(si);
		}

	/// \copydoc IIntegrand::values
		template <int elemDim>
		void evaluate(number vValue[],
		              const MathVector<worldDim> vGlobIP[],
		              GeometricObject* pElem,
		              const MathVector<worldDim> vCornerCoords[],
		              const MathVector<elemDim> vLocIP[],
		              const MathMatrix<elemDim, worldDim> vJT[],
		              const size_t numIP)
		{
		//	get reference object id (i.e. Triangle, Quadrilateral, Tetrahedron, ...)
			ReferenceObjectID roid = (ReferenceObjectID) pElem->reference_object_id();

			//	local finite element id
			const LFEID m_id = m_spGridFct->local_finite_element_id(m_fct);

			try{
		//	get trial space
			const LocalShapeFunctionSet<elemDim>& rTrialSpace =
							LocalShapeFunctionSetProvider::get<elemDim>(roid, m_id);

		//	number of dofs on element
			const size_t num_sh = rTrialSpace.num_sh();

		//	get multiindices of element

			std::vector<MultiIndex<2> > ind;  // 	aux. index array
			m_spGridFct->multi_indices(pElem, m_fct, ind);

		//	check multi indices
			if(ind.size() != num_sh)
				UG_THROW("L2ErrorIntegrand::evaluate: Wrong number of"
						" multi indices.");

		//	loop all integration points
			for(size_t ip = 0; ip < numIP; ++ip)
			{
				//value[ip] = ipvalueFct(vLocIP[ip], vGlobIP[ip], vJT[ip], ind)

			//	compute exact solution at integration point
				number exactSolIP;
				(*m_spExactSolution)(exactSolIP, vGlobIP[ip], m_time, this->subset());

			// 	compute approximated solution at integration point
				number approxSolIP = 0.0;
				for(size_t sh = 0; sh < num_sh; ++sh)
				{
				//	get value at shape point (e.g. corner for P1 fct)
					const number valSH = DoFRef(*m_spGridFct, ind[sh]);

				//	add shape fct at ip * value at shape
					approxSolIP += valSH * rTrialSpace.shape(sh, vLocIP[ip]);
				}

			//	get squared of difference
				vValue[ip] = (exactSolIP - approxSolIP);
				vValue[ip] *= vValue[ip];
			}

			}catch(UGError_LocalShapeFunctionSetNotRegistered& ex)
			{
				UG_THROW("L2ErrorIntegrand::evaluate: "<<ex.get_msg());
			}
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
               SmartPtr<TGridFunction> spGridFct, const char* cmp,
               number time, int quadOrder, const char* subsets)
{
//	get function id of name
	const size_t fct = spGridFct->fct_id_by_name(cmp);

//	check that function exists
	if(fct >= spGridFct->num_fct())
		UG_THROW("L2Error: Function space does not contain"
				" a function with name " << cmp << ".");

	SmartPtr<IIntegrand<number, TGridFunction::dim> > spIntegrand
		= CreateSmartPtr(new L2ErrorIntegrand<TGridFunction>(spExactSol, spGridFct, fct, time));

	return sqrt(IntegrateSubsets(spIntegrand, spGridFct, subsets, quadOrder));
}

#ifdef UG_FOR_LUA
template <typename TGridFunction>
number L2Error(const char* ExactSol,
               SmartPtr<TGridFunction> spGridFct, const char* cmp,
               number time, int quadOrder, const char* subsets)
{
	SmartPtr<UserData<number, TGridFunction::dim> > spExactSol
	 = CreateSmartPtr(new LuaUserData<number, TGridFunction::domain_type::dim>(ExactSol));
	return L2Error(spExactSol, spGridFct, cmp, time, quadOrder, subsets);
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
		SmartPtr<TGridFunction> m_spGridFct;

	///	component of function
		size_t m_fct;

	///	local finite element id
		LFEID m_id;

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
		                 SmartPtr<TGridFunction> gridFct, size_t cmp,
		                 number time)
		: m_spGridFct(gridFct), m_fct(cmp),
		  m_id(m_spGridFct->local_finite_element_id(m_fct)),
		  m_spExactSolution(spExactSol),
		  m_spExactGrad(spExactGrad),
		  m_time(time)
		{}

	///	sets subset
		virtual void set_subset(int si)
		{
			if(!m_spGridFct->is_def_in_subset(m_fct, si))
				UG_THROW("H1Error: Grid function component"
						<<m_fct<<" not defined on subset "<<si);
			IIntegrand<number, worldDim>::set_subset(si);
		}

	/// \copydoc IIntegrand::values
		template <int elemDim>
		void evaluate(number vValue[],
		              const MathVector<worldDim> vGlobIP[],
		              GeometricObject* pElem,
		              const MathVector<worldDim> vCornerCoords[],
		              const MathVector<elemDim> vLocIP[],
		              const MathMatrix<elemDim, worldDim> vJT[],
		              const size_t numIP)
		{
		//	get reference object id (i.e. Triangle, Quadrilateral, Tetrahedron, ...)
			ReferenceObjectID roid = (ReferenceObjectID) pElem->reference_object_id();

			try{
		//	get trial space
			const LocalShapeFunctionSet<elemDim>& rTrialSpace =
							LocalShapeFunctionSetProvider::get<elemDim>(roid, m_id);

		//	number of dofs on element
			const size_t num_sh = rTrialSpace.num_sh();

		//	get multiindices of element
			std::vector<MultiIndex<2> > ind;  // 	aux. index array
			m_spGridFct->multi_indices(pElem, m_fct, ind);

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
					const number valSH = DoFRef(*m_spGridFct, ind[sh]);

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

			}catch(UGError_LocalShapeFunctionSetNotRegistered& ex)
			{
				UG_THROW("H1ErrorIntegrand::evaluate: "<<ex.get_msg());
			}
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

	SmartPtr<IIntegrand<number, TGridFunction::dim> > spIntegrand
		= CreateSmartPtr(new H1ErrorIntegrand<TGridFunction>(spExactSol, spExactGrad, spGridFct, fct, time));

	return sqrt(IntegrateSubsets(spIntegrand, spGridFct, subsets, quadOrder));
}

#ifdef UG_FOR_LUA
template <typename TGridFunction>
number H1Error(const char* ExactSol, const char* ExactGrad,
			   SmartPtr<TGridFunction> spGridFct, const char* cmp,
			   number time, int quadOrder, const char* subsets)
{
	static const int dim = TGridFunction::domain_type::dim;
	SmartPtr<UserData<number, dim> > spExactSol
	 = CreateSmartPtr(new LuaUserData<number, dim>(ExactSol));
	SmartPtr<UserData<MathVector<dim>, dim> > spExactGrad
	 = CreateSmartPtr(new LuaUserData<MathVector<dim>, dim>(ExactGrad));
	return H1Error(spExactSol, spExactGrad, spGridFct, cmp, time, quadOrder, subsets);
}
#endif

////////////////////////////////////////////////////////////////////////////////
// L2 Integrand
////////////////////////////////////////////////////////////////////////////////

template <typename TGridFunction>
class L2FuncIntegrand
	: public StdIntegrand<number, TGridFunction::dim, L2FuncIntegrand<TGridFunction> >
{
	public:
	//	world dimension of grid function
		static const int worldDim = TGridFunction::dim;

	private:
	// grid function
		SmartPtr<TGridFunction> m_spGridFct;

	//	component of function
		const size_t m_fct;

	public:
	/// constructor
		L2FuncIntegrand(SmartPtr<TGridFunction> spGridFct, size_t cmp)
		: m_spGridFct(spGridFct), m_fct(cmp)
		{};

	///	sets subset
		virtual void set_subset(int si)
		{
			if(!m_spGridFct->is_def_in_subset(m_fct, si))
				UG_THROW("L2ErrorIntegrand: Grid function component"
						<<m_fct<<" not defined on subset "<<si);
			IIntegrand<number, worldDim>::set_subset(si);
		}

	/// \copydoc IIntegrand::values
		template <int elemDim>
		void evaluate(number vValue[],
		              const MathVector<worldDim> vGlobIP[],
		              GeometricObject* pElem,
		              const MathVector<worldDim> vCornerCoords[],
		              const MathVector<elemDim> vLocIP[],
		              const MathMatrix<elemDim, worldDim> vJT[],
		              const size_t numIP)
		{
		//	get reference object id (i.e. Triangle, Quadrilateral, Tetrahedron, ...)
			ReferenceObjectID roid = (ReferenceObjectID) pElem->reference_object_id();

			const LFEID m_id = m_spGridFct->local_finite_element_id(m_fct);

			try{
		//	get trial space
			const LocalShapeFunctionSet<elemDim>& rTrialSpace =
							LocalShapeFunctionSetProvider::get<elemDim>(roid, m_id);

		//	number of dofs on element
			const size_t num_sh = rTrialSpace.num_sh();

		//	get multiindices of element

			std::vector<MultiIndex<2> > ind;  // 	aux. index array
			m_spGridFct->multi_indices(pElem, m_fct, ind);

		//	check multi indices
			if(ind.size() != num_sh)
				UG_THROW("L2FuncIntegrand::values: Wrong number of"
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
					const number valSH = BlockRef((*m_spGridFct)[ind[sh][0]], ind[sh][1]);
					approxSolIP += valSH * rTrialSpace.shape(sh, vLocIP[ip]);
				}

				//	get square
				vValue[ip] = approxSolIP*approxSolIP;

			}

			}catch(UGError_LocalShapeFunctionSetNotRegistered& ex)
			{
				UG_THROW("L2FuncIntegrand::values: "<<ex.get_msg());
			}
		};
};

/// interpolates a function on the whole domain or on some subsets
/**
 * This function computes the L2-norm of a grid function.
 *
 * \param[in]		u1			grid function
 * \param[in]		u2			grid function
 * \param[in]		name		symbolic name of function
 * \param[in]		time		time point
 * \param[in]		quadOrder	order of quadrature rule
 * \param[in]		subsets		subsets, where to interpolate
 * \returns			number 		l2-norm of difference
 */
template <typename TGridFunction>
number L2Norm(SmartPtr<TGridFunction> spGridFct, const char* cmp,
              int quadOrder, const char* subsets)
{
//	get function id of name
	const size_t fct = spGridFct->fct_id_by_name(cmp);

//	check that function exists
	if(fct >= spGridFct->num_fct())
		UG_THROW("L2Norm: Function space does not contain"
				" a function with name " << cmp << ".");

	SmartPtr<IIntegrand<number, TGridFunction::dim> > spIntegrand
		= CreateSmartPtr(new L2FuncIntegrand<TGridFunction>(spGridFct, fct));

	return sqrt(IntegrateSubsets(spIntegrand, spGridFct, subsets, quadOrder));
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
		SmartPtr<TGridFunction> m_spGridFct;

	//	component of function
		const size_t m_fct;

	public:
	/// constructor
		StdFuncIntegrand(SmartPtr<TGridFunction> spGridFct, size_t cmp)
		: m_spGridFct(spGridFct), m_fct(cmp)
		{};

	///	sets subset
		virtual void set_subset(int si)
		{
			if(!m_spGridFct->is_def_in_subset(m_fct, si))
				UG_THROW("L2ErrorIntegrand: Grid function component"
						<<m_fct<<" not defined on subset "<<si);
			IIntegrand<number, worldDim>::set_subset(si);
		}

	/// \copydoc IIntegrand::values
		template <int elemDim>
		void evaluate(number vValue[],
		              const MathVector<worldDim> vGlobIP[],
		              GeometricObject* pElem,
		              const MathVector<worldDim> vCornerCoords[],
		              const MathVector<elemDim> vLocIP[],
		              const MathMatrix<elemDim, worldDim> vJT[],
		              const size_t numIP)
		{
		//	get reference object id (i.e. Triangle, Quadrilateral, Tetrahedron, ...)
			ReferenceObjectID roid = (ReferenceObjectID) pElem->reference_object_id();

			const LFEID m_id = m_spGridFct->local_finite_element_id(m_fct);

			try{
		//	get trial space
			const LocalShapeFunctionSet<elemDim>& rTrialSpace =
							LocalShapeFunctionSetProvider::get<elemDim>(roid, m_id);

		//	number of dofs on element
			const size_t num_sh = rTrialSpace.num_sh();

		//	get multiindices of element

			std::vector<MultiIndex<2> > ind;  // 	aux. index array
			m_spGridFct->multi_indices(pElem, m_fct, ind);

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
					const number valSH = BlockRef((*m_spGridFct)[ind[sh][0]], ind[sh][1]);
					approxSolIP += valSH * rTrialSpace.shape(sh, vLocIP[ip]);
				}

				//	get squared of difference
				vValue[ip] = approxSolIP;

			}

			}catch(UGError_LocalShapeFunctionSetNotRegistered& ex)
			{
				UG_THROW("StdFuncIntegrand::evaluate: "<<ex.get_msg());
			}
		};
};


template <typename TGridFunction>
number StdFuncIntegralOnVertex(SmartPtr<TGridFunction> spGridFct,
							   size_t fct,
							   int si)
{
//	integrate elements of subset
	typedef typename TGridFunction::template dim_traits<0>::geometric_base_object geometric_base_object;
	typedef typename TGridFunction::template dim_traits<0>::const_iterator const_iterator;

	//	reset the result
	number integral = 0;

//	note: this iterator is for the base elements, e.g. Face and not
//			for the special type, e.g. Triangle, Quadrilateral
	const_iterator iter = spGridFct->template begin<geometric_base_object>(si);
	const_iterator iterEnd = spGridFct->template end<geometric_base_object>(si);

// 	iterate over all elements
	for(; iter != iterEnd; ++iter)
	{
	//	get element
		geometric_base_object* pElem = *iter;

		std::vector<MultiIndex<2> > ind;  // 	aux. index array
		spGridFct->multi_indices(pElem, fct, ind);

	// 	compute approximated solution at integration point
		number value = 0.0;
		for(size_t sh = 0; sh < ind.size(); ++sh)
		{
			value += BlockRef((*spGridFct)[ind[sh][0]], ind[sh][1]);
		}

	//	add to global sum
		integral += value;

	} // end elem

//	return the summed integral contributions of all elements
	return integral;
}


template <typename TGridFunction>
number Integral(SmartPtr<TGridFunction> spGridFct, const char* cmp,
                const char* subsets, int quadOrder)
{
//	get function id of name
	const size_t fct = spGridFct->fct_id_by_name(cmp);

//	check that function exists
	if(fct >= spGridFct->num_fct())
		UG_THROW("L2Norm: Function space does not contain"
				" a function with name " << cmp << ".");

//	read subsets
	SubsetGroup ssGrp(spGridFct->domain()->subset_handler());
	if(subsets != NULL)
	{
		ConvertStringToSubsetGroup(ssGrp, subsets);
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
			value += StdFuncIntegralOnVertex(spGridFct, fct, ssGrp[s]);

#ifdef UG_PARALLEL
	// sum over processes
	if(pcl::GetNumProcesses() > 1)
	{
		pcl::ProcessCommunicator com;
		number local = value;
		com.allreduce(&local, &value, 1, PCL_DT_DOUBLE, PCL_RO_SUM);
	}
#endif
		return value;
	}

	SmartPtr<IIntegrand<number, TGridFunction::dim> > spIntegrand
		= CreateSmartPtr(new StdFuncIntegrand<TGridFunction>(spGridFct, fct));

	return IntegrateSubsets(spIntegrand, spGridFct, subsets, quadOrder);
}

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
	typedef typename domain_traits<dim>::geometric_base_object geometric_base_object;

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
		geometric_base_object* pElem = *iter;

	//	get all corner coordinates
		CollectCornerCoordinates(vCorner, *pElem, aaPos, true);

	//	compute bf and grads at bip for element
		if(!geo.update(pElem, &vCorner[0], ish))
			UG_THROW("IntegralNormalComponentOnManifold: "
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

template <typename TGridFunction, int dim>
number IntegralNormalComponentOnManifoldSubset(SmartPtr<IIntegrand<MathVector<TGridFunction::dim>, TGridFunction::dim> > spIntegrand,
                                 SmartPtr<TGridFunction> spGridFct,
                                 int si, const SubsetGroup& bndSSGrp, int quadOrder)
{
//	integrate elements of subset
	typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object geometric_base_object;
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator const_iterator;

	spIntegrand->set_subset(si);

	if(quadOrder > 2)
		UG_THROW("IntegrateOverManifold: Currently only middle point rule implemented.");

	return IntegralNormalComponentOnManifoldUsingFV1Geom<TGridFunction::dim,dim,const_iterator>
					(spGridFct->template begin<geometric_base_object>(si),
	                 spGridFct->template end<geometric_base_object>(si),
	                 spGridFct->domain()->position_accessor(),
	                 spGridFct->domain()->subset_handler().get(),
	                 *spIntegrand, bndSSGrp);
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
		ConvertStringToSubsetGroup(innerSSGrp, InnerSubsets);
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
		ConvertStringToSubsetGroup(bndSSGrp, BndSubsets);
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
	if(pcl::GetNumProcesses() > 1)
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
		= CreateSmartPtr(new UserDataIntegrand<MathVector<TGridFunction::dim>, TGridFunction>(spData, spGridFct, time));

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

	if(u.local_finite_element_id(fct) != LFEID(LFEID::LAGRANGE, 1))
		UG_THROW("IntegrateNormalGradientOnManifold:"
				 "Only implemented for Lagrange P1 functions.");

//	read subsets
	SubsetGroup innerSSGrp; innerSSGrp.set_subset_handler(u.domain()->subset_handler());
	if(InnerSubset != NULL)
		ConvertStringToSubsetGroup(innerSSGrp, u.domain()->subset_handler(), InnerSubset);
	else // add all if no subset specified
		innerSSGrp.add_all();

//	read bnd subsets
	SubsetGroup bndSSGrp; bndSSGrp.set_subset_handler(u.domain()->subset_handler());
	if(InnerSubset != NULL){
		ConvertStringToSubsetGroup(bndSSGrp, u.domain()->subset_handler(), BndSubset);
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
		typedef typename domain_traits<dim>::geometric_base_object geometric_base_object;

	//	get iterators for all elems on subset
		typedef typename TGridFunction::template dim_traits<dim>::const_iterator const_iterator;
		const_iterator iter = u.template begin<geometric_base_object>();
		const_iterator iterEnd = u.template end<geometric_base_object>();

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
			geometric_base_object* elem = *iter;

		//	get all corner coordinates
			CollectCornerCoordinates(vCorner, *elem, u.domain()->position_accessor(), true);

		//	compute bf and grads at bip for element
			if(!geo.update(elem, &vCorner[0], u.domain()->subset_handler().get()))
				UG_THROW("IntegrateNormalGradientOnManifold: "
						"Cannot update Finite Volume Geometry.");

		//	get fct multi-indeces of element
			std::vector<MultiIndex<2> > ind;
			u.multi_indices(elem, fct, ind);

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
						const number fctVal = BlockRef(u[ind[sh][0]], ind[sh][1]);

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
	if(pcl::GetNumProcesses() > 1)
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
