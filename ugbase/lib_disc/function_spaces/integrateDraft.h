/*
 * integrate.h
 *
 *  Created on: 04.04.2011
 *      Author: kxylouris, avogel
 */

#ifndef __H__UG__LIB_DISC__FUNCTION_SPACES__INTEGRATEDRAFT__
#define __H__UG__LIB_DISC__FUNCTION_SPACES__INTEGRATEDRAFT__

#include <cmath>

#include "common/common.h"

#include "lib_disc/common/subset_group.h"
#include "lib_disc/quadrature/quadrature.h"
#include "lib_disc/local_finite_element/local_shape_function_set.h"
#include <boost/function.hpp>

namespace ug{


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
		virtual void values(GeometricObject* pElem,
		                    const MathVector<elemDim> vLocIP[],
		                    const MathVector<worldDim> vGlobIP[],
		                    const MathMatrix<worldDim, elemDim> vJT[],
		                    const size_t numIP,
		                    number value[]) = 0;

		virtual ~IIntegrand() {}
};


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
			integrand.values(pElem, rQuadRule.points(),
							&(vGlobIP[0]), &(vJT[0]),
							numIP, &(vValue[0]));
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
			UG_THROW_FATAL("SumValuesOnElems: " << ex.get_msg() << ".");
		}
		}catch(UG_ERROR_ReferenceMappingMissing& ex){
			UG_THROW_FATAL("SumValuesOnElems: " << ex.get_msg() << ".");
		}
	} // end elem

//	return the summed integral contributions of all elements
	return integral;
}



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

	public:
	/// constructor
		L2ErrorIntegrand(IPData<number, worldDim>& exactSol,
		                 TGridFunction& gridFct, size_t cmp,
		                 number time)
		: m_rGridFct(gridFct), m_fct(cmp),
		  m_ExactSolution(exactSol), m_time(time)
		{};

	/// \copydoc IIntegrand::getValues
		virtual void values(GeometricObject* pElem,
		                    const MathVector<elemDim> vLocIP[],
		                    const MathVector<worldDim> vGlobIP[],
		                    const MathMatrix<worldDim, elemDim> vJT[],
		                    const size_t numIP,
		                    number value[])
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
				UG_THROW_FATAL("L2ErrorIntegrand::values: Wrong number of"
						" multi indices.");

		//	loop all integration points
			for(size_t ip = 0; ip < numIP; ++ip)
			{
				//value[ip] = ipvalueFct(vLocIP[ip], vGlobIP[ip], vJT[ip], ind)

			//	compute exact solution at integration point
				number exactSolIP;
				m_ExactSolution(exactSolIP, vGlobIP[ip], m_time);

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
				value[ip] = (exactSolIP - approxSolIP);
				value[ip] *= value[ip];
			}

			}catch(UG_ERROR_LocalShapeFunctionSetNotRegistered& ex)
			{
				UG_THROW_FATAL("L2ErrorIntegrand::getValues: "<<ex.get_msg());
			}
		};
};


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
		virtual void values(GeometricObject* pElem,
		                    const MathVector<elemDim> vLocIP[],
		                    const MathVector<worldDim> vGlobIP[],
		                    const MathMatrix<worldDim, elemDim> vJT[],
		                    const size_t numIP,
		                    number value[])
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
				UG_THROW_FATAL("L2ErrorIntegrand::values: Wrong number of"
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
				MathMatrix<elemDim,worldDim> JTInv;
				Inverse(JTInv, vJT[ip]);
				MatVecMult(approxGradIP, JTInv, locTmp);

			//	get squared of difference
				value[ip] = (exactSolIP - approxSolIP) * (exactSolIP - approxSolIP);
				value[ip] += VecDistanceSq(approxGradIP, exactGradIP);
			}

			}catch(UG_ERROR_LocalShapeFunctionSetNotRegistered& ex)
			{
				UG_THROW_FATAL("L2ErrorIntegrand::getValues: "<<ex.get_msg());
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
number L2ErrorDraft(IPData<number, TGridFunction::dim>& ExactSol,
                    TGridFunction& u, const char* name, number time, int quadOrder, const char* subsets)
{
//	get function id of name
	const size_t fct = u.fct_id_by_name(name);

//	check that function exists
	if(fct >= u.num_fct())
		UG_THROW_FATAL("L2ErrorDraft: Function space does not contain"
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
			UG_THROW_FATAL("L2ErrorDraft: Element dimension does not match world dimension!");


	//	create integration kernel
		static const int dim = TGridFunction::dim;
		L2ErrorIntegrand<TGridFunction, dim> integrandKernel(ExactSol, u, fct, time);

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
		virtual void values(GeometricObject* pElem,
		                    const MathVector<elemDim> vLocIP[],
		                    const MathVector<worldDim> vGlobIP[],
		                    const MathMatrix<worldDim, elemDim> vJT[],
		                    const size_t numIP,
		                    number value[])
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
				UG_THROW_FATAL("L2ErrorIntegrand::values: Wrong number of"
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
				value[ip] = approxSolIP*approxSolIP;

			}

			}catch(UG_ERROR_LocalShapeFunctionSetNotRegistered& ex)
			{
				UG_THROW_FATAL("L2ErrorIntegrand::values: "<<ex.get_msg());
			}
		};
};


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
		virtual void values(GeometricObject* pElem,
		                    const MathVector<elemDim> vLocIP[],
		                    const MathVector<worldDim> vGlobIP[],
		                    const MathMatrix<worldDim, elemDim> vJT[],
		                    const size_t numIP,
		                    number value[])
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
				UG_THROW_FATAL("StdFuncIntegrand::values: Wrong number of"
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
				value[ip] = approxSolIP;

			}

			}catch(UG_ERROR_LocalShapeFunctionSetNotRegistered& ex)
			{
				UG_THROW_FATAL("StdFuncIntegrand::values: "<<ex.get_msg());
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
		UG_THROW_FATAL("L2ErrorFunc: Function space does not contain"
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
			UG_THROW_FATAL("L2Norm: Element dimension does not match world dimension!");

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

template <typename TGridFunction>
number StdFuncIntegral(TGridFunction& u, const char* name, int quadOrder, const char* subsets)
{
//	get function id of name
	const size_t fct = u.fct_id_by_name(name);

//	check that function exists
	if(fct >= u.num_fct())
		UG_THROW_FATAL("StdIntegral: Function space does not contain"
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
			UG_THROW_FATAL("StdIntegral: Element dimension does not match world dimension!");

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

} // namespace ug

#endif /*__H__UG__LIB_DISC__FUNCTION_SPACES__INTEGRATEDRAFT__*/
