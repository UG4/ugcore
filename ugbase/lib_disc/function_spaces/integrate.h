/*
 * integrate.h
 *
 *  Created on: 04.04.2011
 *      Author: avogel, anaegel (some parts)
 */

#ifndef __H__UG__LIB_DISC__FUNCTION_SPACES__INTEGRATEDRAFT__
#define __H__UG__LIB_DISC__FUNCTION_SPACES__INTEGRATEDRAFT__

#include <cmath>

#include "common/common.h"

#include "lib_disc/common/subset_group.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/quadrature/quadrature.h"
#include "lib_disc/local_finite_element/local_shape_function_set.h"
#include "lib_disc/spatial_disc/disc_util/finite_volume_geometry.h"
#include "lib_disc/spatial_disc/ip_data/ip_data.h"
#include "lib_disc/spatial_disc/ip_data/const_user_data.h"
#include <boost/function.hpp>

#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// Integrand Interface
////////////////////////////////////////////////////////////////////////////////

/// Abstract integrand interface
template <int TWorldDim, int TElemDim>
class IIntegrand
{
	public:
	///	reference dimension of elements
		static const int elemDim = TElemDim;

	///	world dimension
		static const int worldDim = TWorldDim;

	public:
	///	returns the values of the integrand for a bunch of ips
	/**
	 * \param[in]	pElem	the element to integrate
	 * \param[out]	value	the value of the integrand at the ips
	 * \param[in]	vLocIP	local integration point positions
	 * \param[in]	vGlobIP	global integration point positions
	 * \param[in]	numIP	number of integration points
	 */
		virtual void values(number vValue[],
		                    const MathVector<worldDim> vGlobIP[],
		                    GeometricObject* pElem,
	                        const MathVector<worldDim> vCornerCoords[],
		                    const MathVector<elemDim> vLocIP[],
		                    const MathMatrix<elemDim, worldDim> vJT[],
		                    const size_t numIP) = 0;

		virtual ~IIntegrand() {}
};

////////////////////////////////////////////////////////////////////////////////
// Generic Integration Routine
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
                 IIntegrand<WorldDim,dim>& integrand,
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
	typedef typename domain_traits<dim>::position_type position_type;

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
		std::vector<position_type> vCorner;
		CollectCornerCoordinates(vCorner, *pElem, aaPos, true);

	//	update the reference mapping for the corners
		mapping.update(&vCorner[0]);

	//	compute global integration points
		std::vector<position_type> vGlobIP(numIP);
		mapping.local_to_global(&(vGlobIP[0]), rQuadRule.points(), numIP);

	//	compute transformation matrices
		std::vector<MathMatrix<WorldDim, dim> > vJT(numIP);
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
			const number det = Determinant(vJT[ip]);

		//	add contribution of integration point
			intValElem += vValue[ip] * weightIP * det;
		}

	//	add to global sum
		integral += intValElem;

		}catch(UG_ERROR_QuadratureRuleNotRegistered& ex){
			UG_THROW("SumValuesOnElems: " << ex.get_msg() << ".");
		}
		}catch(UG_ERROR_ReferenceMappingMissing& ex){
			UG_THROW("SumValuesOnElems: " << ex.get_msg() << ".");
		}
	} // end elem

//	return the summed integral contributions of all elements
	return integral;
}

////////////////////////////////////////////////////////////////////////////////
// IPData Integrand
////////////////////////////////////////////////////////////////////////////////

template <typename TGridFunction, int TDim = TGridFunction::dim>
class DirectIPDataIntegrand : public IIntegrand<TGridFunction::dim, TDim>
{
	public:
	//	world dimension of grid function
		static const int worldDim = TGridFunction::dim;

	//	reference element dimension
		static const int elemDim = TDim;

	private:
	//  data to integrate
		SmartPtr<IDirectIPData<number, worldDim> > m_spData;

	// 	grid function
		SmartPtr<TGridFunction> m_spGridFct;

	//	time
		number m_time;

	//	subset
		int m_si;

	public:
	/// constructor
		DirectIPDataIntegrand(SmartPtr<IDirectIPData<number, worldDim> > spData,
		                      SmartPtr<TGridFunction> spGridFct,
		                      number time, int si)
		: m_spData(spData), m_spGridFct(spGridFct), m_time(time), m_si(si)
		{
			m_spData->set_function_pattern(spGridFct->function_pattern());
		};

	/// constructor
		DirectIPDataIntegrand(SmartPtr<IDirectIPData<number, worldDim> > spData,
							  number time, int si)
		: m_spData(spData), m_spGridFct(NULL), m_time(time), m_si(si)
		{
			if(m_spData->requires_grid_fct())
				UG_THROW("DirectIPDataIntegrand: Missing GridFunction, but "
						" data requires grid function.")
		};

	/// \copydoc IIntegrand::getValues
		virtual void values(number vValue[],
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
					(*m_spData)(vValue, vGlobIP, m_time, m_si, u, pElem,
								vCornerCoords, vLocIP, numIP, &vJT[0]);
				}
				UG_CATCH_THROW("DirectIPDataIntegrand: Cannot evaluate data.");
			}
			else
			{
			//	compute data
				try{
					(*m_spData)(vValue, vGlobIP, m_time, m_si, numIP);
				}
				UG_CATCH_THROW("DirectIPDataIntegrand: Cannot evaluate data.");
			}

		};
};

template <typename TGridFunction, int dim>
number IntegrateDirectIPDataInDim(SmartPtr<IDirectIPData<number, TGridFunction::dim> > spData,
                                  SmartPtr<TGridFunction> spGridFct,
                                  number time,
                                  int quadOrder, int si)
{
	DirectIPDataIntegrand<TGridFunction, dim> integrand(spData, spGridFct, time, si);

//	integrate elements of subset
	typedef typename domain_traits<dim>::geometric_base_object geometric_base_object;

	return Integrate(spGridFct->template begin<geometric_base_object>(si),
	                   spGridFct->template end<geometric_base_object>(si),
	                   spGridFct->domain()->position_accessor(),
	                   integrand,
	                   quadOrder);
}

template <typename TGridFunction>
number Integral(SmartPtr<IDirectIPData<number, TGridFunction::dim> > spData,
                SmartPtr<TGridFunction> spGridFct,
                const char* subsets, number time,
                int quadOrder)
{
//	world dimensions
	static const int dim = TGridFunction::dim;

//	read subsets
	SubsetGroup ssGrp(spGridFct->domain()->subset_handler());
	if(subsets != NULL)
	{
		ConvertStringToSubsetGroup(ssGrp, subsets);
		if(!SameDimensionsInAllSubsets(ssGrp))
			UG_THROW("Integral: Subsets '"<<subsets<<"' do not have same dimension."
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
	for(size_t i = 0; i < ssGrp.num_subsets(); ++i)
	{
	//	get subset index
		const int si = ssGrp[i];

	//	check dimension
		if(ssGrp.dim(si) > dim)
			UG_THROW("Integral: Dimension of subset is "<<ssGrp.dim(si)<<", but "
			         " World Dimension is "<<dim<<". Cannot integrate this.");

	//	integrate elements of subset
		switch(ssGrp.dim(si))
		{
			case 1: value += IntegrateDirectIPDataInDim<TGridFunction, dim>(spData, spGridFct, time, quadOrder, si); break;
			case 2: value += IntegrateDirectIPDataInDim<TGridFunction, dim>(spData, spGridFct, time, quadOrder, si); break;
			case 3: value += IntegrateDirectIPDataInDim<TGridFunction, dim>(spData, spGridFct, time, quadOrder, si); break;
			default: UG_THROW("Integral: Dimension "<<ssGrp.dim(si)<<" not supported. "
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


///////////////
// const data
///////////////

template <typename TGridFunction>
number Integral(number val, SmartPtr<TGridFunction> spGridFct,
                const char* subsets,
                number time, int quadOrder)
{
	static const int dim = TGridFunction::dim;
	SmartPtr<IDirectIPData<number, dim> > sp =
			CreateSmartPtr(new ConstUserNumber<dim>(val));
	return Integral(sp, spGridFct, subsets, time, quadOrder);
}

template <typename TGridFunction>
number Integral(number val, SmartPtr<TGridFunction> spGridFct,const char* subsets,number time)
{return Integral(val, spGridFct, subsets, time, 1);}

template <typename TGridFunction>
number Integral(number val, SmartPtr<TGridFunction> spGridFct,const char* subsets)
{return Integral(val, spGridFct, subsets, 0.0, 1);}

template <typename TGridFunction>
number Integral(number val, SmartPtr<TGridFunction> spGridFct,number time)
{return Integral(val, spGridFct, NULL, time, 1);}

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
	SmartPtr<IDirectIPData<number, dim> > sp =
			LuaUserDataFactory<number, dim>::create(luaFct);
	return Integral(sp, spGridFct, subsets, time, quadOrder);
}
#endif

////////////////////////////////////////////////////////////////////////////////
// L2 Error Integrand
////////////////////////////////////////////////////////////////////////////////

template <typename TGridFunction, int TDim = TGridFunction::dim>
class L2ErrorIntegrand : public IIntegrand<TGridFunction::dim, TDim>
{
	public:
	//	world dimension of grid function
		static const int worldDim = TGridFunction::dim;

	//	reference element dimension
		static const int elemDim = TDim;

	private:
	// grid function
		TGridFunction& m_rGridFct;

	//	component of function
		const size_t m_fct;

	//  exact solution
		IPData<number, worldDim>& m_ExactSolution;

	//	time
		number m_time;

	//	subset
		int m_si;

	public:
	/// constructor
		L2ErrorIntegrand(IPData<number, worldDim>& exactSol,
		                 TGridFunction& gridFct, size_t cmp,
		                 number time, int si)
		: m_rGridFct(gridFct), m_fct(cmp),
		  m_ExactSolution(exactSol), m_time(time), m_si(si)
		{};

	/// \copydoc IIntegrand::getValues
		virtual void values(number vValue[],
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
			const LFEID m_id = m_rGridFct.local_finite_element_id(m_fct);

			try{
		//	get trial space
			const DimLocalShapeFunctionSet<elemDim>& rTrialSpace =
							LocalShapeFunctionSetProvider::get<elemDim>(roid, m_id);

		//	number of dofs on element
			const size_t num_sh = rTrialSpace.num_sh();

		//	get multiindices of element

			std::vector<MultiIndex<2> > ind;  // 	aux. index array
			m_rGridFct.multi_indices(pElem, m_fct, ind);

		//	check multi indices
			if(ind.size() != num_sh)
				UG_THROW("L2ErrorIntegrand::values: Wrong number of"
						" multi indices.");

		//	loop all integration points
			for(size_t ip = 0; ip < numIP; ++ip)
			{
				//value[ip] = ipvalueFct(vLocIP[ip], vGlobIP[ip], vJT[ip], ind)

			//	compute exact solution at integration point
				number exactSolIP;
				m_ExactSolution(exactSolIP, vGlobIP[ip], m_time, m_si);

			// 	compute approximated solution at integration point
				number approxSolIP = 0.0;
				for(size_t sh = 0; sh < num_sh; ++sh)
				{
				//	get value at shape point (e.g. corner for P1 fct)
					const number valSH = BlockRef(m_rGridFct[ind[sh][0]], ind[sh][1]);

				//	add shape fct at ip * value at shape
					approxSolIP += valSH * rTrialSpace.shape(sh, vLocIP[ip]);
				}

			//	get squared of difference
				vValue[ip] = (exactSolIP - approxSolIP);
				vValue[ip] *= vValue[ip];
			}

			}catch(UG_ERROR_LocalShapeFunctionSetNotRegistered& ex)
			{
				UG_THROW("L2ErrorIntegrand::getValues: "<<ex.get_msg());
			}
		};
};

/// interpolates a function on the whole domain or on some subsets
/**
 * This function computes the L2-difference between a given analytic function
 * and a grid function.
 *
 * \param[in]		ExactSol	analytic function
 * \param[in]		u			grid function
 * \param[in]		name		symbolic name of function
 * \param[in]		time		time point
 * \param[in]		quadOrder	order of quadrature rule
 * \param[in]		subsets		subsets, where to interpolate
 * \returns			number 		l2-norm of difference
 */
template <typename TGridFunction>
number L2Error(IPData<number, TGridFunction::dim>& ExactSol,
                    TGridFunction& u, const char* name, number time, int quadOrder, const char* subsets)
{
//	get function id of name
	const size_t fct = u.fct_id_by_name(name);

//	check that function exists
	if(fct >= u.num_fct())
		UG_THROW("L2Error: Function space does not contain"
				" a function with name " << name << ".");

//	create subset group
	SubsetGroup ssGrp; ssGrp.set_subset_handler(u.domain()->subset_handler());

//	read subsets
	if(subsets != NULL)
		ConvertStringToSubsetGroup(ssGrp, u.domain()->subset_handler(), subsets);
	else // add all if no subset specified
		ssGrp.add_all();

//	reset l2norm
	number l2norm2 = 0;

//	loop subsets
	for(size_t i = 0; i < ssGrp.num_subsets(); ++i)
	{
	//	get subset index
		const int si = ssGrp[i];

	//	skip if function is not defined in subset
		if(!u.is_def_in_subset(fct, si)) continue;


		if (ssGrp.dim(i) != TGridFunction::dim)
			UG_THROW("L2Error: Element dimension does not match world dimension!");


	//	create integration kernel
		static const int dim = TGridFunction::dim;
		L2ErrorIntegrand<TGridFunction, dim> integrandKernel(ExactSol, u, fct, time, si);

	//	integrate elements of subset
		typedef typename domain_traits<dim>::geometric_base_object geometric_base_object;
		l2norm2 += Integrate(u.template begin<geometric_base_object>(si),
							 u.template end<geometric_base_object>(si),
							 u.domain()->position_accessor(),
							 integrandKernel,
							 quadOrder);
	}

//	return the sqrt of the result
	return sqrt(l2norm2);
}

#ifdef UG_FOR_LUA
template <typename TGridFunction>
number L2Error(const char* ExactSol,
                         TGridFunction& u, const char* name, number time, int quadOrder, const char* subsets)
{
	LuaUserData<number, TGridFunction::domain_type::dim> p(ExactSol);
	return L2Error(p, u, name, time, quadOrder, subsets);
}
#endif

////////////////////////////////////////////////////////////////////////////////
// H1 Error Integrand
////////////////////////////////////////////////////////////////////////////////

template <typename TGridFunction, int TDim = TGridFunction::dim>
class H1ErrorIntegrand : IIntegrand<TGridFunction::dim, TDim>
{
	public:
	//	world dimension of grid function
		static const int worldDim = TGridFunction::dim;

	//	reference element dimension
		static const int elemDim = TDim;

	private:
	// grid function
		TGridFunction& m_rGridFct;

	//	component of function
		size_t m_fct;

	//	local finite element id
		LFEID m_id;

	// 	exact solution
		IPData<number, worldDim>& m_ExactSolution;

	// 	exact gradient
		IPData<MathVector<worldDim>, worldDim>& m_ExactGrad;

	//	time
		number m_time;

	public:
	/// constructor
		H1ErrorIntegrand(IPData<number, worldDim>& exactSol,
		                 IPData<MathVector<worldDim>, worldDim>& exactGrad,
		                 TGridFunction& gridFct, size_t cmp,
		                 number time)
		: m_rGridFct(gridFct), m_fct(cmp),
		  m_id(m_rGridFct.local_finite_element_id(m_fct)),
		  m_ExactSolution(exactSol),
		  m_ExactGrad(exactGrad),
		  m_time(time)
		{}

	/// \copydoc IIntegrand::getValues
		virtual void values(number vValue[],
		                    const MathVector<worldDim> vGlobIP[],
		                    GeometricObject* pElem,
	                        const MathVector<worldDim> vCornerCoords[],
		                    const MathVector<elemDim> vLocIP[],
		                    const MathMatrix<elemDim, worldDim> vJT[],
		                    const size_t numIP)
		{
		//	get reference object id (i.e. Triangle, Quadrilateral, Tetrahedron, ...)
			ReferenceObjectID roid = (ReferenceObjectID) pElem->reference_object_id();

			// 	aux. index array
					typename TGridFunction::multi_index_vector_type ind;

			try{
		//	get trial space
			const DimLocalShapeFunctionSet<elemDim>& rTrialSpace =
							LocalShapeFunctionSetProvider::get<elemDim>(roid, m_id);

		//	number of dofs on element
			const size_t num_sh = rTrialSpace.num_sh();

		//	get multiindices of element
			m_rGridFct.multi_indices(pElem, m_fct, ind);

		//	check multi indices
			if(ind.size() != num_sh)
				UG_THROW("L2ErrorIntegrand::values: Wrong number of"
						" multi indices.");

		//	loop all integration points
			std::vector<MathVector<elemDim> > vLocGradient(num_sh);
			for(size_t ip = 0; ip < numIP; ++ip)
			{
			//	compute exact solution at integration point
				number exactSolIP;
				m_ExactSolution(exactSolIP, vGlobIP[ip], m_time);

			//	compute exact gradient at integration point
				MathVector<worldDim> exactGradIP;
				m_ExactGrad(exactGradIP, vGlobIP[ip], m_time);

			//	compute shape gradients at ip
				rTrialSpace.grads(&vLocGradient[0], vGlobIP[ip]);

			// 	compute approximated solution at integration point
				number approxSolIP = 0.0;
				MathVector<elemDim> locTmp = 0.0;
				for(size_t sh = 0; sh < num_sh; ++sh)
				{
				//	get value at shape point (e.g. corner for P1 fct)
					const number valSH = BlockRef(m_rGridFct[ind[sh][0]], ind[sh][1]);

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

			}catch(UG_ERROR_LocalShapeFunctionSetNotRegistered& ex)
			{
				UG_THROW("L2ErrorIntegrand::getValues: "<<ex.get_msg());
			}
		};
};

////////////////////////////////////////////////////////////////////////////////
// L2 Integrand
////////////////////////////////////////////////////////////////////////////////

template <typename TGridFunction, int TDim = TGridFunction::dim>
class L2FuncIntegrand : public IIntegrand<TGridFunction::dim, TDim>
{
	public:
	//	world dimension of grid function
		static const int worldDim = TGridFunction::dim;

	//	reference element dimension
		static const int elemDim = TDim;

	private:
	// grid function
		TGridFunction& m_rGridFct;

	//	component of function
		const size_t m_fct;

	public:
	/// constructor
		L2FuncIntegrand(TGridFunction& gridFct, size_t cmp)
		: m_rGridFct(gridFct), m_fct(cmp)
		{};

	/// \copydoc IIntegrand::getValues
		virtual void values(number vValue[],
		                    const MathVector<worldDim> vGlobIP[],
		                    GeometricObject* pElem,
	                        const MathVector<worldDim> vCornerCoords[],
		                    const MathVector<elemDim> vLocIP[],
		                    const MathMatrix<elemDim, worldDim> vJT[],
		                    const size_t numIP)
		{
		//	get reference object id (i.e. Triangle, Quadrilateral, Tetrahedron, ...)
			ReferenceObjectID roid = (ReferenceObjectID) pElem->reference_object_id();

			const LFEID m_id = m_rGridFct.local_finite_element_id(m_fct);

			try{
		//	get trial space
			const DimLocalShapeFunctionSet<elemDim>& rTrialSpace =
							LocalShapeFunctionSetProvider::get<elemDim>(roid, m_id);

		//	number of dofs on element
			const size_t num_sh = rTrialSpace.num_sh();

		//	get multiindices of element

			std::vector<MultiIndex<2> > ind;  // 	aux. index array
			m_rGridFct.multi_indices(pElem, m_fct, ind);

		//	check multi indices
			if(ind.size() != num_sh)
				UG_THROW("L2ErrorIntegrand::values: Wrong number of"
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
					const number valSH = BlockRef(m_rGridFct[ind[sh][0]], ind[sh][1]);
					approxSolIP += valSH * rTrialSpace.shape(sh, vLocIP[ip]);
				}

				//	get square
				vValue[ip] = approxSolIP*approxSolIP;

			}

			}catch(UG_ERROR_LocalShapeFunctionSetNotRegistered& ex)
			{
				UG_THROW("L2ErrorIntegrand::values: "<<ex.get_msg());
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
number L2Norm(TGridFunction& u, const char* name, int quadOrder, const char* subsets)
{
//	get function id of name
	const size_t fct = u.fct_id_by_name(name);

//	check that function exists
	if(fct >= u.num_fct())
		UG_THROW("L2ErrorFunc: Function space does not contain"
				" a function with name " << name << ".");

//	create subset group
	SubsetGroup ssGrp; ssGrp.set_subset_handler(u.domain()->subset_handler());

//	read subsets
	if(subsets != NULL)
		ConvertStringToSubsetGroup(ssGrp, u.domain()->subset_handler(), subsets);
	else // add all if no subset specified
		ssGrp.add_all();

//	reset l2norm
	number l2norm2 = 0;

//	loop subsets
	for(size_t i = 0; i < ssGrp.num_subsets(); ++i)
	{
	//	get subset index
		const int si = ssGrp[i];

	//	skip if function is not defined in subset
		if(!u.is_def_in_subset(fct, si)) continue;


		if (ssGrp.dim(i) != TGridFunction::dim)
			UG_THROW("L2Norm: Element dimension does not match world dimension!");

	//	create integration kernel
		static const int dim = TGridFunction::dim;
		L2FuncIntegrand<TGridFunction, dim> integrand(u, fct);

	//	integrate elements of subset
		typedef typename domain_traits<dim>::geometric_base_object geometric_base_object;
		l2norm2 += Integrate(u.template begin<geometric_base_object>(si),
							 u.template end<geometric_base_object>(si),
							 u.domain()->position_accessor(),
							 integrand,
							 quadOrder);
	}

#ifdef UG_PARALLEL
	// sum over all processes
	if(pcl::GetNumProcesses() > 1)
	{
		pcl::ProcessCommunicator com;
		number local = l2norm2;
		com.allreduce(&local, &l2norm2, 1, PCL_DT_DOUBLE, PCL_RO_SUM);
	}
#endif
//	return the sqrt of the result
	return sqrt(l2norm2);
}

////////////////////////////////////////////////////////////////////////////////
// Standard Integrand
////////////////////////////////////////////////////////////////////////////////

template <typename TGridFunction, int TDim = TGridFunction::dim>
class StdFuncIntegrand : public IIntegrand<TGridFunction::dim, TDim>
{
	public:
	//	world dimension of grid function
		static const int worldDim = TGridFunction::dim;

	//	reference element dimension
		static const int elemDim = TDim;

	private:
	// grid function
		TGridFunction& m_rGridFct;

	//	component of function
		const size_t m_fct;

	public:
	/// constructor
		StdFuncIntegrand(TGridFunction& gridFct, size_t cmp)
		: m_rGridFct(gridFct), m_fct(cmp)
		{};

	/// \copydoc IIntegrand::getValues
		virtual void values(number vValue[],
		                    const MathVector<worldDim> vGlobIP[],
		                    GeometricObject* pElem,
	                        const MathVector<worldDim> vCornerCoords[],
		                    const MathVector<elemDim> vLocIP[],
		                    const MathMatrix<elemDim, worldDim> vJT[],
		                    const size_t numIP)
		{
		//	get reference object id (i.e. Triangle, Quadrilateral, Tetrahedron, ...)
			ReferenceObjectID roid = (ReferenceObjectID) pElem->reference_object_id();

			const LFEID m_id = m_rGridFct.local_finite_element_id(m_fct);

			try{
		//	get trial space
			const DimLocalShapeFunctionSet<elemDim>& rTrialSpace =
							LocalShapeFunctionSetProvider::get<elemDim>(roid, m_id);

		//	number of dofs on element
			const size_t num_sh = rTrialSpace.num_sh();

		//	get multiindices of element

			std::vector<MultiIndex<2> > ind;  // 	aux. index array
			m_rGridFct.multi_indices(pElem, m_fct, ind);

		//	check multi indices
			if(ind.size() != num_sh)
				UG_THROW("StdFuncIntegrand::values: Wrong number of"
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
					const number valSH = BlockRef(m_rGridFct[ind[sh][0]], ind[sh][1]);
					approxSolIP += valSH * rTrialSpace.shape(sh, vLocIP[ip]);
				}

				//	get squared of difference
				vValue[ip] = approxSolIP;

			}

			}catch(UG_ERROR_LocalShapeFunctionSetNotRegistered& ex)
			{
				UG_THROW("StdFuncIntegrand::values: "<<ex.get_msg());
			}
		};
};


template <typename TGridFunction>
number StdFuncIntegral(TGridFunction& u, const char* name, int quadOrder, const char* subsets)
{
//	get function id of name
	const size_t fct = u.fct_id_by_name(name);

//	check that function exists
	if(fct >= u.num_fct())
		UG_THROW("StdIntegral: Function space does not contain"
				" a function with name " << name << ".");

//	create subset group
	SubsetGroup ssGrp; ssGrp.set_subset_handler(u.domain()->subset_handler());

//	read subsets
	if(subsets != NULL)
		ConvertStringToSubsetGroup(ssGrp, u.domain()->subset_handler(), subsets);
	else // add all if no subset specified
		ssGrp.add_all();

//	reset value
	number value = 0;

//	loop subsets
	for(size_t i = 0; i < ssGrp.num_subsets(); ++i)
	{
	//	get subset index
		const int si = ssGrp[i];

	//	skip if function is not defined in subset
		if(!u.is_def_in_subset(fct, si)) continue;


		if (ssGrp.dim(i) != TGridFunction::dim)
			UG_THROW("StdIntegral: Element dimension does not match world dimension!");

	//	create integration kernel
		static const int dim = TGridFunction::dim;
		StdFuncIntegrand<TGridFunction, dim> integrand(u, fct);

	//	integrate elements of subset
		typedef typename domain_traits<dim>::geometric_base_object geometric_base_object;
		value += Integrate(u.template begin<geometric_base_object>(si),
							 u.template end<geometric_base_object>(si),
							 u.domain()->position_accessor(),
							 integrand,
							 quadOrder);
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
number IntegrateFluxOnBoundary(TGridFunction& u, const char* cmp,
                               const char* BndSubset, const char* InnerSubset = NULL)
{
//	get function id of name
	const size_t fct = u.fct_id_by_name(cmp);

//	check that function exists
	if(fct >= u.num_fct())
		UG_THROW("IntegrateFluxOnBoundary: Function space does not contain"
				" a function with name " << cmp << ".");

	if(u.local_finite_element_id(fct) != LFEID(LFEID::LAGRANGE, 1))
		UG_THROW("IntegrateFluxOnBoundary:"
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
		UG_THROW("IntegrateFluxOnBoundary: No boundary subsets specified. Aborting.");
	}

//	reset value
	number value = 0;

//	loop subsets
	for(size_t i = 0; i < innerSSGrp.num_subsets(); ++i)
	{
	//	get subset index
		const int si = innerSSGrp[i];

	//	skip if function is not defined in subset
		if(!u.is_def_in_subset(fct, si)) continue;


		if (innerSSGrp.dim(i) != TGridFunction::dim)
			UG_THROW("IntegrateFluxOnBoundary: Element dimension does not match world dimension!");

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
		for(size_t s = 0; s < bndSSGrp.num_subsets(); ++s)
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
				UG_THROW("IntegrateFluxOnBoundary: "
						"Cannot update Finite Volume Geometry.");

		//	get fct multi-indeces of element
			std::vector<MultiIndex<2> > ind;
			u.multi_indices(elem, fct, ind);

		//	specify, which subsets are boundary
			for(size_t s = 0; s < bndSSGrp.num_subsets(); ++s)
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
					          "IntegrateFluxOnBoundary::values: Wrong number of"
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


/// Integrates the AceticAcid Flux over some boundary subsets
/**
 * This function integrates \f$ \rho \omega \vec{q} \cdot \vec{n} \f$m, where
 * \f$ \vec{q} := -\frac{K}{\mu}(\nabla p - \rho \vec{g}) \f$ is the
 * Darcy velocity over some selected
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
number IntegrateAceticAcidFluxOnBoundary(TGridFunction& u, const char* PressCmp,
                                         const char* OmegaCmp,
                                         number Permeability,
                                         number Viscosity,
                                         number Density,
                                         number GravityNorm,
                                         const char* BndSubset, const char* InnerSubset = NULL)
{
//	get function id of name
	const size_t pressFct = u.fct_id_by_name(PressCmp);
	const size_t omegaFct = u.fct_id_by_name(OmegaCmp);

//	check that function exists
	if(pressFct >= u.num_fct())
		UG_THROW("IntegrateFluxOnBoundary: Function space does not contain"
				" a function with name " << PressCmp << ".");
	if(omegaFct >= u.num_fct())
		UG_THROW("IntegrateFluxOnBoundary: Function space does not contain"
				" a function with name " << omegaFct << ".");

	if(u.local_finite_element_id(pressFct) != LFEID(LFEID::LAGRANGE, 1))
		UG_THROW("IntegrateFluxOnBoundary:"
				"Only implemented for Lagrange P1 functions.");
	if(u.local_finite_element_id(omegaFct) != LFEID(LFEID::LAGRANGE, 1))
		UG_THROW("IntegrateFluxOnBoundary:"
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
		UG_THROW("IntegrateFluxOnBoundary: No boundary subsets specified. Aborting.");
	}

//	reset value
	number value = 0;

//	loop subsets
	for(size_t i = 0; i < innerSSGrp.num_subsets(); ++i)
	{
	//	get subset index
		const int si = innerSSGrp[i];

	//	skip if function is not defined in subset
		if(!u.is_def_in_subset(pressFct, si)) continue;


		if (innerSSGrp.dim(i) != TGridFunction::dim)
			UG_THROW("IntegrateFluxOnBoundary: Element dimension does not match world dimension!");

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
		for(size_t s = 0; s < bndSSGrp.num_subsets(); ++s)
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
				UG_THROW("IntegrateFluxOnBoundary: "
								"Cannot update Finite Volume Geometry.");

		//	get fct multi-indeces of element
			std::vector<MultiIndex<2> > indPressure;
			u.multi_indices(elem, pressFct, indPressure);
			std::vector<MultiIndex<2> > indOmega;
			u.multi_indices(elem, omegaFct, indOmega);

		//	specify, which subsets are boundary
			for(size_t s = 0; s < bndSSGrp.num_subsets(); ++s)
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
					UG_ASSERT(indPressure.size() != bf.num_sh(),
					          "IntegrateFluxOnBoundary::values: Wrong number of"
										" multi indices.");

				// 	compute gradient of solution at bip
					MathVector<dim> gradPressure; VecSet(gradPressure, 0.0);
					for(size_t sh = 0; sh < bf.num_sh(); ++sh)
					{
						const number fctVal = BlockRef(u[indPressure[sh][0]], indPressure[sh][1]);

						VecScaleAdd(gradPressure, 1.0, gradPressure, fctVal, bf.global_grad(sh));
					}

				// 	compute omega at bip
					number omega = 0.0;
					for(size_t sh = 0; sh < bf.num_sh(); ++sh)
					{
						const number fctVal = BlockRef(u[indOmega[sh][0]], indOmega[sh][1]);

						omega += fctVal * bf.shape(sh);
					}

					MathVector<dim> flux; VecSet(flux, 0.0);
					flux[dim-1] = GravityNorm * Density;

					VecScaleAppend(flux, -1.0, gradPressure);

					VecScale(flux, flux, Density* omega* Permeability/Viscosity);

				//	sum up contributions
					value += VecDot(flux, normal);
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

#endif /*__H__UG__LIB_DISC__FUNCTION_SPACES__INTEGRATEDRAFT__*/
