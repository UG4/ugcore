/*
 * integrate.h
 *
 *  Created on: 04.04.2011
 *      Author: kxylouris, avogel
 */

#ifndef __H__LIBDISCRETIZATION__FUNCTION_SPACES__INTEGRATEDRAFT__
#define __H__LIBDISCRETIZATION__FUNCTION_SPACES__INTEGRATEDRAFT__

#include <cmath>

#include "common/common.h"

#include "lib_discretization/common/subset_group.h"
#include "lib_discretization/domain_util.h"
#include "lib_discretization/quadrature/quadrature.h"
#include "lib_discretization/local_finite_element/local_shape_function_set.h"
#include <boost/function.hpp>

namespace ug{


/// Abstract integrand interface
template <int TWorldDim, int TDim>
class IIntegrand
{
	public:
	///	reference dimension of elements
		static const int dim = TDim;

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
		virtual bool values(GeometricObject* pElem,
		                    number value[],
		                    const MathVector<dim> vLocIP[],
		                    const MathVector<worldDim> vGlobIP[],
		                    const MathMatrix<worldDim, dim> vJT[],
		                    size_t numIP) = 0;
};


template <typename TGridFunction, int TDim = TGridFunction::dim>
class L2ErrorIntegrand : public IIntegrand<TGridFunction::dim, TDim>
{
	public:
	//	world dimension of grid function
		static const int worldDim = TGridFunction::dim;

	//	reference element dimension
		static const int dim = TDim;

	private:
	// grid function
		TGridFunction& m_rGridFct;

	//	component of function
		size_t m_fct;

	//	local finite element id
		LFEID m_id;

	//  exact solution
		typedef boost::function<void (number&, const MathVector<dim>&, number)> ExactSolFunctor;
		ExactSolFunctor m_ExactSolution;

	//	time
		number m_time;

	// 	aux. index array
		typename TGridFunction::multi_index_vector_type ind;

	public:

	/// constructor
		L2ErrorIntegrand(const ExactSolFunctor& exactSol,
		                 TGridFunction& gridFct, size_t cmp,
		                 number time)
		: m_rGridFct(gridFct), m_fct(cmp),
		  m_id(m_rGridFct.local_finite_element_id(m_fct)),
		  m_ExactSolution(exactSol), m_time(time)
		{};

	/// \copydoc IIntegrand::getValues
		virtual bool values(GeometricObject* pElem,
		                    number value[],
		                    const MathVector<dim> vLocIP[],
		                    const MathVector<worldDim> vGlobIP[],
		                    const MathMatrix<worldDim, dim> vJT[],
		                    size_t numIP)
		{
		//	get reference object id (i.e. Triangle, Quadrilateral, Tetrahedron, ...)
			ReferenceObjectID roid = (ReferenceObjectID) pElem->reference_object_id();

			try{
		//	get trial space
			const DimLocalShapeFunctionSet<dim>& rTrialSpace =
							LocalShapeFunctionSetProvider::get<dim>(roid, m_id);

		//	number of dofs on element
			const size_t num_sh = rTrialSpace.num_sh();

		//	get multiindices of element
			m_rGridFct.multi_indices(pElem, m_fct, ind);

		//	check multi indices
			if(ind.size() != num_sh)
			{
				UG_LOG("ERROR in 'L2ErrorIntegrand::values': Wrong number of"
						" multi indices.\n");
				return false;
			}

		//	loop all integration points
			for(size_t ip = 0; ip < numIP; ++ip)
			{
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
				UG_LOG("ERROR in 'L2ErrorIntegrand::getValues': "<<ex.get_msg()<<"\n");
				return false;
			}

			return true;
		};
};


template <typename TGridFunction, int TDim = TGridFunction::dim>
class H1ErrorIntegrand : IIntegrand<TGridFunction::dim, TDim>
{
	public:
	//	world dimension of grid function
		static const int worldDim = TGridFunction::dim;

	//	reference element dimension
		static const int dim = TDim;

	private:
	// grid function
		TGridFunction& m_rGridFct;

	//	component of function
		size_t m_fct;

	//	local finite element id
		LFEID m_id;

	// 	exact solution
		typedef boost::function<void (number&, const MathVector<dim>&, number)> ExactSolFunctor;
		ExactSolFunctor m_ExactSolution;

	// 	exact gradient
		typedef boost::function<void (MathVector<dim>&, const MathVector<dim>&, number)> ExactGradFunctor;
		ExactGradFunctor m_ExactGrad;

	//	time
		number m_time;

	// 	aux. index array
		typename TGridFunction::multi_index_vector_type ind;

	public:

	/// constructor
		H1ErrorIntegrand(ExactSolFunctor& exactSol,
		                 ExactGradFunctor& exactGrad,
		                 TGridFunction& gridFct, size_t cmp,
		                 number time)
		: m_rGridFct(gridFct), m_fct(cmp),
		  m_id(m_rGridFct.local_finite_element_id(m_fct)),
		  m_ExactSolution(exactSol),
		  m_ExactGrad(exactGrad),
		  m_time(time)
		{}

	/// \copydoc IIntegrand::getValues
		virtual bool values(GeometricObject* pElem,
		                    number value[],
		                    const MathVector<dim> vLocIP[],
		                    const MathVector<worldDim> vGlobIP[],
		                    const MathMatrix<worldDim, dim> vJT[],
		                    size_t numIP)
		{
		//	get reference object id (i.e. Triangle, Quadrilateral, Tetrahedron, ...)
			ReferenceObjectID roid = (ReferenceObjectID) pElem->reference_object_id();

			try{
		//	get trial space
			const DimLocalShapeFunctionSet<dim>& rTrialSpace =
							LocalShapeFunctionSetProvider::get<dim>(roid, m_id);

		//	number of dofs on element
			const size_t num_sh = rTrialSpace.num_sh();

		//	get multiindices of element
			m_rGridFct.multi_indices(pElem, m_fct, ind);

		//	check multi indices
			if(ind.size() != num_sh)
			{
				UG_LOG("ERROR in 'L2ErrorIntegrand::values': Wrong number of"
						" multi indices.\n");
				return false;
			}

		//	loop all integration points
			std::vector<MathVector<dim> > vLocGradient(num_sh);
			for(size_t ip = 0; ip < numIP; ++ip)
			{
			//	compute exact solution at integration point
				number exactSolIP;
				m_ExactSolution(exactSolIP, vGlobIP[ip], m_time);

			//	compute exact gradient at integration point
				MathVector<dim> exactGradIP;
				m_ExactGrad(exactGradIP, vGlobIP[ip], m_time);

			//	compute shape gradients at ip
				rTrialSpace.grads(&vLocGradient[0], vGlobIP[ip]);

			// 	compute approximated solution at integration point
				number approxSolIP = 0.0;
				MathVector<dim> locTmp = 0.0;
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
				MathMatrix<dim,worldDim> JTInv;
				Inverse(JTInv, vJT[ip]);
				MatVecMult(approxGradIP, JTInv, locTmp);

			//	get squared of difference
				value[ip] = (exactSolIP - approxSolIP) * (exactSolIP - approxSolIP);
				value[ip] += VecDistanceSq(approxGradIP, exactGradIP);
			}

			}catch(UG_ERROR_LocalShapeFunctionSetNotRegistered& ex)
			{
				UG_LOG("ERROR in 'L2ErrorIntegrand::getValues': "<<ex.get_msg()<<"\n");
				return false;
			}

			return true;
		};
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
template <int WorldDim, int dim>
number Integrate(typename domain_traits<dim>::const_iterator iterBegin,
                 typename domain_traits<dim>::const_iterator iterEnd,
                 typename domain_traits<WorldDim>::position_accessor_type aaPos,
                 IIntegrand<WorldDim,dim>& integrand, int quadOrder)
{
//	reset the result
	number integral = 0;

//	note: this iterator is for the base elements, e.g. Face and not
//			for the special type, e.g. Triangle, Quadrilateral
	typename domain_traits<dim>::const_iterator iter = iterBegin;

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
		CollectCornerCoordinates(vCorner, *pElem, aaPos);

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
		if(!integrand.values(pElem, &(vValue[0]), rQuadRule.points(), &(vGlobIP[0]), &(vJT[0]), numIP))
		{
			UG_LOG("Unable to compute values of integrand at integration point.\n");
			throw(UGFatalError("Cannot compute integration values."));
		}

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
			UG_LOG("ERROR in 'SumValuesOnElems': " << ex.get_msg() << ".\n");
			return false;
		}
		}catch(UG_ERROR_ReferenceMappingMissing& ex){
			UG_LOG("ERROR in 'SumValuesOnElems': " << ex.get_msg() << ".\n");
			return false;
		}
	} // end elem

//	return the summed integral contributions of all elements
	return integral;
}



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
number L2ErrorDraft(
	const boost::function<void (number& res, const MathVector<TGridFunction::dim>& x, number m_time)>& ExactSol,
	TGridFunction& u, const char* name, number time, int quadOrder, const char* subsets)
{
//	get Function Pattern
	const typename TGridFunction::approximation_space_type& approxSpace
				= u.get_approximation_space();

//	get function id of name
	const size_t fct = approxSpace.fct_id_by_name(name);

//	check that function exists
	if(fct >= u.num_fct())
	{
		UG_LOG("ERROR in L2Error: Function space does not contain"
				" a function with name " << name << ".\n");
		return false;
	}

//	create subset group
	SubsetGroup ssGrp; ssGrp.set_subset_handler(*approxSpace.get_subset_handler());

//	read subsets
	if(subsets != NULL)
		ConvertStringToSubsetGroup(ssGrp, *approxSpace.get_subset_handler(), subsets);
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

	//	switch dimensions
		switch(ssGrp.dim(i))
		{
/*			case 1:
			{
			//	create integration Kernel
				L2ErrorIntegrand<TGridFunction, 1> integrandKernel(ExactSol, u, fct, time);

			//	integrate elements of subset
				l2norm2 += Integrate(u.template begin<EdgeBase>(si),
									 u.template end<EdgeBase>(si),
									 u.get_domain().get_position_accessor(),
									 integrandKernel,
									 quadOrder);
			}
*/			case 2:
			{
			//	create integration Kernel
				L2ErrorIntegrand<TGridFunction, 2> integrandKernel(ExactSol, u, fct, time);

			//	integrate elements of subset
				l2norm2 += Integrate(u.template begin<Face>(si),
									 u.template end<Face>(si),
									 u.get_domain().get_position_accessor(),
									 integrandKernel,
									 quadOrder);
			}
/*			case 3:
			{
			//	create integration Kernel
				L2ErrorIntegrand<TGridFunction, 3> integrandKernel(ExactSol, u, fct, time);

			//	integrate elements of subset
				l2norm2 += Integrate(u.template begin<Volume>(si),
									 u.template end<Volume>(si),
									 u.get_domain().get_position_accessor(),
									 integrandKernel,
									 quadOrder);
			}
*/			default:
				UG_LOG("ERROR in L2Error: Dimension "<<ssGrp.dim(i) <<
				       " not supported in world dimension "<<TGridFunction::dim<<".");
				throw(UGFatalError("Dimension not supported."));
		}
	}

//	return the sqrt of the result
	return sqrt(l2norm2);
}

} // namespace ug

#endif /*__H__LIBDISCRETIZATION__FUNCTION_SPACES__INTEGRATEDRAFT__*/
