/*
 * lagrange.h
 *
 *  Created on: 17.11.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET__LAGRANGE__LAGRANGE__
#define __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET__LAGRANGE__LAGRANGE__

#include "../common/lagrange1d.h"
#include "../local_shape_function_set.h"
#include "../local_dof_pattern.h"
#include "lib_discretization/common/multi_index.h"
#include "common/util/metaprogramming_util.h"

namespace ug{

template <typename TRefElem, size_t TOrder>
struct LagrangeLSFS{};

template <>
template <size_t TOrder>
class LagrangeLSFS<ReferenceEdge, TOrder>
{
	public:
	///	Reference Element type
		typedef ReferenceEdge reference_element_type;

	///	Order of Shape functions
		static const size_t p = TOrder;

	///	Dimension, where shape functions are defined
		static const int dim = reference_element_type::dim;

	///	Domain position type
		typedef MathVector<dim> position_type;

	///	Shape type
		typedef number shape_type;

	///	Gradient type
		typedef MathVector<dim> grad_type;

	/// Number of shape functions
		static const size_t nsh = p+1;

	///	Multi Index type
		typedef MultiIndex<dim> multi_index_type;

	public:
	///	Constructor
		LagrangeLSFS()
		{
			m_vPolynom.resize(nsh);
			m_vDPolynom.resize(nsh);

			for(size_t i = 0; i < nsh; ++i)
			{
				m_vPolynom[i] = EquidistantLagrange1D(i, p);
				m_vDPolynom[i] = m_vPolynom[i].derivative();
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		bool position(size_t i, position_type& pos) const
		{
			pos = EquidistantLagrange1D::position(i, p);
			return true;
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		shape_type shape(size_t i, const position_type& x) const
		{
			return m_vPolynom[i].value(x[0]);
		}

	///	\copydoc ug::LocalShapeFunctionSet::shapes()
		void shapes(shape_type* sOut, const position_type& x) const
		{
			for(size_t sh = 0; sh < num_sh(); ++sh)
				sOut[sh] = shape(sh, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		grad_type grad(size_t i, const position_type& x) const
		{
			grad_type tmpGrad;
			eval_grad(tmpGrad, i, x);
			return tmpGrad;
		}

	///	\copydoc ug::LocalShapeFunctionSet::grads()
		void grads(grad_type* gOut, const position_type& x) const
		{
			for(size_t sh = 0; sh < num_sh(); ++sh)
				eval_grad(gOut[sh], sh, x);
		}

	///	given an index i of the element, we return the multiindex
		inline static MultiIndex<1> multi_index(size_t i)
		{
			UG_ASSERT(i < nsh, "i must be smaller than Number of DoFs.");
			return MultiIndex<1>(i);
		}

	protected:
		void eval_grad(grad_type& grad, size_t i, const position_type& x) const
		{
			grad[0] = m_vDPolynom[i].value(x[0]);
		}

	protected:
		std::vector<Polynomial1D> m_vPolynom;	///< Shape Polynomials
		std::vector<Polynomial1D> m_vDPolynom;	///< Derivative of Shape Polynomial
};

template <>
template <size_t TOrder>
class LagrangeLSFS<ReferenceTriangle, TOrder>
{
	public:
	///	dimension of reference element
		static const int dim = ReferenceTriangle::dim;

	///	order of lagrange function
		static const size_t p = TOrder;

	///	number of shapes
		static const size_t nsh = BinomialCoefficient<dim + p, p>::value;

	///	MultiIndex type
		typedef MultiIndex<dim> multi_index_type;

	///	return the index for a multi_index
		inline static size_t index(const MultiIndex<dim>& ind)
		{
			UG_ASSERT(ind[0] <= p, "Wrong Multiindex.");
			UG_ASSERT(ind[1] <= p, "Wrong Multiindex.");
			UG_ASSERT(ind[0] + ind[1] <= p, "Wrong Multiindex.");

		//	todo: Replace loop by explicit formula
			size_t res = ind[0];
			for(size_t i = 0; i < ind[1]; ++i)
				res += (p+1-i);

			UG_ASSERT(res < nsh, "Wrong index computation. Internal error.");

			return res;
		}

	///	return the multi_index for an index
		inline static MultiIndex<dim> multi_index(size_t i)
		{
			UG_ASSERT(i < nsh, "i must be smaller than Number of DoFs.");

		//	todo: replace by explicit formula (iff possible)
			int i0 = i;
			int i1;
			for(i1 = 0; i1 < p; ++i1)
			{
				if(i0 - (p+1-i1) < 0) break;
				i0 -= (p+1-i1);
			}

			return MultiIndex<dim>( i0, i1 );
		}

	///	return shape at point x for index i
		inline number shape(const size_t i, const MathVector<dim>& x)
		{
			UG_ASSERT(i < nsh, "i must be smaller than Number of DoFs.");

		//	get corresponding multi index
			MultiIndex<dim> ind = multi_index(i);

		//	forward
			return shape(ind, x);
		}

	///	return shape at point x for a MultiIndex
		inline number shape(const MultiIndex<dim>& ind, const MathVector<dim>& x)
		{
			UG_ASSERT(ind[0] <= p, "Wrong Multiindex.");
			UG_ASSERT(ind[1] <= p, "Wrong Multiindex.");
			UG_ASSERT(ind[0] + ind[1] <= p, "Wrong Multiindex.");

		//	get adjoint barycentric index
			const size_t i0 = p + 1 + ind[0] + ind[1];

			return m_vPolynom[ ind[0] ].value(x[0])
					* m_vPolynom[ ind[1] ].value(x[1])
					* m_vPolynom[ i0 ].value(1-x[0]-x[1]);
		}

	private:
		std::vector<Polynomial1D> m_vPolynom;
		std::vector<Polynomial1D> m_vDPolynom;

};

template <>
template <size_t TOrder>
struct LagrangeLSFS<ReferenceTetrahedron, TOrder>
{
	static const int dim = ReferenceTetrahedron::dim;
	static const size_t nsh = BinomialCoefficient<dim + TOrder, TOrder>::value;
};

template <>
template <size_t TOrder>
struct LagrangeLSFS<ReferenceQuadrilateral, TOrder>
{
	static const int dim = ReferenceQuadrilateral::dim;
	static const size_t p = TOrder;
	static const size_t nsh = (p+1)*(p+1);

	typedef MultiIndex<dim> multi_index_type;

//	return the index for a multi_index
	inline static size_t index(const MultiIndex<dim>& ind)
	{
		UG_ASSERT(ind[0] < p, "Wrong Multiindex.");
		UG_ASSERT(ind[1] < p, "Wrong Multiindex.");

		return ind[1] * p + ind[0];
	}

//	return the multi_index for an index
	inline static MultiIndex<dim> multi_index(size_t i)
	{
		UG_ASSERT(i < nsh, "i must be smaller than Number of DoFs.");

		return MultiIndex<dim>( i%p, i/p );
	}
};

template <>
template <size_t TOrder>
struct LagrangeLSFS<ReferenceHexahedron, TOrder>
{
	static const int dim = ReferenceHexahedron::dim;
	static const size_t nsh = (TOrder+1)*(TOrder+1);
};


template <typename TImpl>
class LocalShapeFunctionSetWrapper
	: public ug::LocalShapeFunctionSet<typename TImpl::reference_element_type>,
	  public TImpl
{
	/// Implementation
		typedef TImpl ImplType;

	public:
	///	Reference Element type
		typedef typename ImplType::reference_element_type reference_element_type;

	///	Order of Shape functions
		static const size_t order = ImplType::p;

	///	Dimension, where shape functions are defined
		static const int dim = ImplType::dim;

	///	Domain position type
		typedef typename ImplType::position_type position_type;

	///	Shape type
		typedef typename ImplType::shape_type shape_type;

	///	Gradient type
		typedef typename ImplType::grad_type grad_type;

	/// Number of shape functions
		static const size_t nsh = ImplType::nsh;

	public:
	///	constructor
		LocalShapeFunctionSetWrapper(){}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		virtual size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		virtual bool position(size_t i, position_type& pos) const
		{
			return ImplType::position(i, pos);
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		virtual shape_type shape(size_t i, const position_type& x) const
		{
			return ImplType::shape(i, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::shapes()
		virtual void shapes(shape_type* sOut, const position_type& x) const
		{
			for(size_t sh = 0; sh < num_sh(); ++sh)
				sOut[sh] = ImplType::shape(sh, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		virtual grad_type grad(size_t i, const position_type& x) const
		{
			return ImplType::grad(i, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grads()
		virtual void grads(grad_type* gOut, const position_type& x) const
		{
			for(size_t sh = 0; sh < num_sh(); ++sh)
				gOut[sh] = ImplType::grad(sh, x);
		}

		const LocalDoFPattern<reference_element_type>& local_dof_pattern() const
		{
			return m_ElementDoFPattern;
		}

	private:
		LocalDoFPattern<reference_element_type> m_ElementDoFPattern;
};

} //namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET__LAGRANGE__LAGRANGE__ */
