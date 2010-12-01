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

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

/// Lagrange Shape Function Set without virtual functions
template <typename TRefElem, size_t TOrder>
struct LagrangeLSFS{};

/// specialization for Edges
/**
 * Lagrange shape function of any order for the Reference Edge
 * \tparam 	TOrder		requested order
 */
template <>
template <size_t TOrder>
class LagrangeLSFS<ReferenceEdge, TOrder>
{
	private:
	///	abbreviation for order
		static const size_t p = TOrder;

	public:
	///	Reference Element type
		typedef ReferenceEdge reference_element_type;

	///	Order of Shape functions
		static const size_t order = TOrder;

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
			//	create equidistant polynomials and its derivatives
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

	///	return Multi index for index i
		inline MultiIndex<dim> multi_index(size_t i) const
		{
			UG_ASSERT(i < nsh, "i must be smaller than Number of DoFs.");
			return MultiIndex<1>(i);
		}

	///	return the index for a multi_index
		inline size_t index(const MultiIndex<dim>& ind) const
		{
			UG_ASSERT(ind[0] < nsh, "ind[0] must be smaller than Number of DoFs.");
			return ind[0];
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
	private:
	///	abbreviation for order
		static const size_t p = TOrder;

	public:
	///	Reference Element type
		typedef ReferenceTriangle reference_element_type;

	///	Order of Shape functions
		static const size_t order = TOrder;

	///	Dimension, where shape functions are defined
		static const int dim = reference_element_type::dim;

	///	Domain position type
		typedef MathVector<dim> position_type;

	///	Shape type
		typedef number shape_type;

	///	Gradient type
		typedef MathVector<dim> grad_type;

	/// Number of shape functions
		static const size_t nsh = BinomialCoefficient<dim + p, p>::value;

	///	Multi Index type
		typedef MultiIndex<dim> multi_index_type;

	public:
	///	Constructor
		LagrangeLSFS()
		{
			m_vPolynom.resize(p+1);
			m_vDPolynom.resize(p+1);

			for(size_t i = 0; i <= p; ++i)
			{
			//	create trancated equidistant polynomials and its derivatives
				m_vPolynom[i] = TruncatedEquidistantLagrange1D(i, p);
				m_vDPolynom[i] = m_vPolynom[i].derivative();
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		bool position(size_t i, position_type& pos) const
		{
		//	get Multi Index
			MultiIndex<dim> ind = multi_index(i);

		//	set position
			for(int d = 0; d < dim; ++d)
				pos[d] = TruncatedEquidistantLagrange1D::position(ind[d], p);

			return true;
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline number shape(const size_t i, const MathVector<dim>& x) const
		{
		//	forward
			return shape(multi_index(i), x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::shapes()
		void shapes(shape_type* sOut, const position_type& x) const
		{
			for(size_t sh = 0; sh < num_sh(); ++sh)
				sOut[sh] = shape(sh, x);
		}

	///	shape value for a Multi Index
		inline number shape(const MultiIndex<dim>& ind, const MathVector<dim>& x) const
		{
			check_multi_index(ind);

		//	get adjoint barycentric index
			const size_t i0 = p - ind[0] - ind[1];
			const number x0 = 1.0 - x[0] - x[1];

			return    m_vPolynom[ ind[0] ].value(x[0])
					* m_vPolynom[ ind[1] ].value(x[1])
					* m_vPolynom[ i0     ].value(x0);
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


	///	return the index for a multi_index
		inline size_t index(const MultiIndex<dim>& ind) const
		{
			check_multi_index(ind);

		//	todo: Replace loop by explicit formula
			size_t res = ind[0];
			for(size_t i = 0; i < ind[1]; ++i)
				res += (p+1-i);

			UG_ASSERT(res < nsh, "Wrong index computation. Internal error.");

			return res;
		}

	///	return the multi_index for an index
		inline MultiIndex<dim> multi_index(size_t i) const
		{
			UG_ASSERT(i < nsh, "i must be smaller than Number of DoFs.");

		//	todo: replace by explicit formula (iff possible)
			int i0 = i;
			int i1;
			for(i1 = 0; i1 < (int)p; ++i1)
			{
				const int diff = i0 - (p+1-i1);
				if(diff < 0) break;
				i0 = diff;
			}

			UG_ASSERT(i0 >= 0, "i0 is negative ("<<i0<<")");
			UG_ASSERT(i1 >= 0, "i1 is negative ("<<i1<<")");

			return MultiIndex<dim>( i0, i1 );
		}

	protected:
		inline void check_multi_index(const MultiIndex<dim>& ind) const
		{
			UG_ASSERT(ind[0] <= p, "Wrong Multiindex.");
			UG_ASSERT(ind[1] <= p, "Wrong Multiindex.");
			UG_ASSERT(ind[0] + ind[1] <= p, "Wrong Multiindex.");
		}

		void eval_grad(grad_type& grad, const size_t i,
		               	   	   	   	   	   const position_type& x) const
		{
			return eval_grad(grad, multi_index(i), x);
		}

		void eval_grad(grad_type& grad, const MultiIndex<dim> ind,
		               	   	   	   	   	   const position_type& x) const
		{
			check_multi_index(ind);

		//	get adjoint barycentric index and position
			const size_t i0 = p - ind[0] - ind[1];
			const number x0 = 1.0 - x[0] - x[1];

			UG_ASSERT(i0 <= p && i0 >= 0, "Wrong Multiindex.");
			UG_ASSERT(x0 <= 1.0 && x0 >= 0.0, "Wrong Position.");

		//	loop dimensions
			for(int d = 0; d < dim; ++d)
			{
				grad[d] = m_vDPolynom[ind[d]].value(x[d])
						* m_vPolynom[i0].value(x0);
				grad[d] += (-1) * m_vDPolynom[i0].value(x0)
						   * m_vPolynom[ind[d]].value(x[d]);

			//	multiply by all functions not depending on x[d]
				for(int d2 = 0; d2 < dim; ++d2)
				{
				// 	skip own value
					if(d2 == d) continue;

					grad[d] *= m_vPolynom[ind[d2]].value(x[d2]);
				}
			}
		}

	private:
		std::vector<Polynomial1D> m_vPolynom;
		std::vector<Polynomial1D> m_vDPolynom;

};


template <>
template <size_t TOrder>
class LagrangeLSFS<ReferenceQuadrilateral, TOrder>
{
	private:
	///	abbreviation for order
		static const size_t p = TOrder;

	public:
	///	Reference Element type
		typedef ReferenceQuadrilateral reference_element_type;

	///	Order of Shape functions
		static const size_t order = TOrder;

	///	Dimension, where shape functions are defined
		static const int dim = reference_element_type::dim;

	///	Domain position type
		typedef MathVector<dim> position_type;

	///	Shape type
		typedef number shape_type;

	///	Gradient type
		typedef MathVector<dim> grad_type;

	/// Number of shape functions
		static const size_t nsh = (p+1)*(p+1);

	///	Multi Index type
		typedef MultiIndex<dim> multi_index_type;

	public:
	///	Constructor
		LagrangeLSFS()
		{
			m_vPolynom.resize(p+1);
			m_vDPolynom.resize(p+1);

			for(size_t i = 0; i <= p; ++i)
			{
			//	create trancated equidistant polynomials and its derivatives
				m_vPolynom[i] = EquidistantLagrange1D(i, p);
				m_vDPolynom[i] = m_vPolynom[i].derivative();
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		bool position(size_t i, position_type& pos) const
		{
		//	get Multi Index
			MultiIndex<dim> ind = multi_index(i);

		//	set position
			for(int d = 0; d < dim; ++d)
				pos[d] = EquidistantLagrange1D::position(ind[d], p);

			return true;
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline number shape(const size_t i, const MathVector<dim>& x) const
		{
		//	forward
			return shape(multi_index(i), x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::shapes()
		void shapes(shape_type* sOut, const position_type& x) const
		{
			for(size_t sh = 0; sh < num_sh(); ++sh)
				sOut[sh] = shape(sh, x);
		}

	///	shape value for a Multi Index
		inline number shape(const MultiIndex<dim>& ind, const MathVector<dim>& x) const
		{
			check_multi_index(ind);

			return    m_vPolynom[ ind[0] ].value(x[0])
					* m_vPolynom[ ind[1] ].value(x[1]);
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


	///	return the index for a multi_index
		inline size_t index(const MultiIndex<dim>& ind) const
		{
			check_multi_index(ind);

			return ind[1] * (p+1) + ind[0];
		}

	///	return the multi_index for an index
		inline MultiIndex<dim> multi_index(size_t i) const
		{
			UG_ASSERT(i < nsh, "i must be smaller than Number of DoFs.");

			return MultiIndex<dim>( i%(p+1), i/(p+1) );
		}

	protected:
		inline void check_multi_index(const MultiIndex<dim>& ind) const
		{
			UG_ASSERT(ind[0] <= p && ind[0] >= 0, "Wrong Multiindex.");
			UG_ASSERT(ind[1] <= p && ind[1] >= 0, "Wrong Multiindex.");
		}

		void eval_grad(grad_type& grad, const size_t i,
		               	   	   	   	    const position_type& x) const
		{
			return eval_grad(grad, multi_index(i), x);
		}

		void eval_grad(grad_type& grad, const MultiIndex<dim> ind,
		               	   	   	   	   	   const position_type& x) const
		{
			check_multi_index(ind);

		//	loop dimensions
			for(int d = 0; d < dim; ++d)
			{
				grad[d] = m_vDPolynom[ind[d]].value(x[d]);

			//	multiply by all functions not depending on x[d]
				for(int d2 = 0; d2 < dim; ++d2)
				{
				// 	skip own value
					if(d2 == d) continue;

					grad[d] *= m_vPolynom[ind[d2]].value(x[d2]);
				}
			}
		}

	private:
		std::vector<Polynomial1D> m_vPolynom;
		std::vector<Polynomial1D> m_vDPolynom;
};


template <>
template <size_t TOrder>
class LagrangeLSFS<ReferenceTetrahedron, TOrder>
{
	private:
	///	abbreviation for order
		static const size_t p = TOrder;

	public:
	///	Reference Element type
		typedef ReferenceTetrahedron reference_element_type;

	///	Order of Shape functions
		static const size_t order = TOrder;

	///	Dimension, where shape functions are defined
		static const int dim = reference_element_type::dim;

	///	Domain position type
		typedef MathVector<dim> position_type;

	///	Shape type
		typedef number shape_type;

	///	Gradient type
		typedef MathVector<dim> grad_type;

	/// Number of shape functions
		static const size_t nsh = BinomialCoefficient<dim + p, p>::value;

	///	Multi Index type
		typedef MultiIndex<dim> multi_index_type;

	public:
	///	Constructor
		LagrangeLSFS()
		{
			m_vPolynom.resize(p+1);
			m_vDPolynom.resize(p+1);

			for(size_t i = 0; i <= p; ++i)
			{
			//	create trancated equidistant polynomials and its derivatives
				m_vPolynom[i] = TruncatedEquidistantLagrange1D(i, p);
				m_vDPolynom[i] = m_vPolynom[i].derivative();
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		bool position(size_t i, position_type& pos) const
		{
		//	get Multi Index
			MultiIndex<dim> ind = multi_index(i);

		//	set position
			for(int d = 0; d < dim; ++d)
				pos[d] = TruncatedEquidistantLagrange1D::position(ind[d], p);

			return true;
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline number shape(const size_t i, const MathVector<dim>& x) const
		{
		//	forward
			return shape(multi_index(i), x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::shapes()
		void shapes(shape_type* sOut, const position_type& x) const
		{
			for(size_t sh = 0; sh < num_sh(); ++sh)
				sOut[sh] = shape(sh, x);
		}

	///	shape value for a Multi Index
		inline number shape(const MultiIndex<dim>& ind, const MathVector<dim>& x) const
		{
			check_multi_index(ind);

		//	get adjoint barycentric index
			const size_t i0 = p - ind[0] - ind[1] - ind[2];
			const number x0 = 1.0 - x[0] - x[1] - x[2];

			return    m_vPolynom[ ind[0] ].value(x[0])
					* m_vPolynom[ ind[1] ].value(x[1])
					* m_vPolynom[ ind[2] ].value(x[2])
					* m_vPolynom[ i0 ].value(x0);
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


	///	return the index for a multi_index
		inline size_t index(const MultiIndex<dim>& ind) const
		{
			check_multi_index(ind);

		//	todo: Replace loop by explicit formula
		//	add x range
			size_t res = ind[0];

		//	add y range
			for(size_t i = 0; i < ind[1]; ++i)
				res += (p+1-ind[2]-i);

		//	add z range
			for(size_t i2 = 0; i2 < ind[2]; ++i2)
				res += BinomCoeff(p+2-i2, p-i2);

			UG_ASSERT(res < nsh, "Wrong index computation. Internal error.");

			return res;
		}

	///	return the multi_index for an index
		inline MultiIndex<dim> multi_index(size_t i) const
		{
			UG_ASSERT(i < nsh, "i must be smaller than Number of DoFs.");

		//	todo: replace by explicit formula (iff possible)
			int i0 = i;
			int i1, i2;
			for(i2 = 0; i2 <= (int)p; ++i2)
			{
				const int binom = BinomCoeff(p+2-i2, p-i2);

				// if i2 is correct
				const int diff = i0 - binom;
				if(diff < 0)
				{
					for(i1 = 0; i1 <= (int)p; ++i1)
					{
						// if i1 is correct return values
						const int diff =  i0 - (p+1-i2-i1);
						if(diff < 0)
							return MultiIndex<dim>( i0, i1, i2);

						// else decrease i1
						i0 = diff;
					}
				}
				// else go one level lower
				else
					i0 = diff;
			}

			UG_ASSERT(i0 >= 0, "i0 is negative ("<<i0<<")");
			UG_ASSERT(i1 >= 0, "i1 is negative ("<<i1<<")");
			UG_ASSERT(i2 >= 0, "i1 is negative ("<<i2<<")");

			UG_ASSERT(0, "Should not reach this line.");
			return MultiIndex<dim>( i0, i1, i2);
		}

	protected:
		inline void check_multi_index(const MultiIndex<dim>& ind) const
		{
			UG_ASSERT(ind[0] <= p, "Wrong Multiindex.");
			UG_ASSERT(ind[1] <= p, "Wrong Multiindex.");
			UG_ASSERT(ind[2] <= p, "Wrong Multiindex.");
			UG_ASSERT(ind[0] + ind[1] + ind[2] <= p, "Wrong Multiindex.");
		}

		void eval_grad(grad_type& grad, const size_t i,
		               	   	   	   	    const position_type& x) const
		{
			return eval_grad(grad, multi_index(i), x);
		}

		void eval_grad(grad_type& grad, const MultiIndex<dim> ind,
		               	   	   	   	   	   const position_type& x) const
		{
			check_multi_index(ind);

		//	get adjoint barycentric index and position
			const size_t i0 = p - ind[0] - ind[1] - ind[2];
			const number x0 = 1 - x[0] - x[1] - x[2];

			UG_ASSERT(i0 <= p && i0 >= 0, "Wrong Multiindex.");
			UG_ASSERT(x0 <= 1.0 && x0 >= 0.0, "Wrong Position.");

		//	loop dimensions
			for(int d = 0; d < dim; ++d)
			{
				grad[d] = m_vDPolynom[ind[d]].value(x[d])
						* m_vPolynom[i0].value(x0);
				grad[d] += (-1) * m_vDPolynom[i0].value(x0)
						   * m_vPolynom[ind[d]].value(x[d]);

			//	multiply by all functions not depending on x[d]
				for(int d2 = 0; d2 < dim; ++d2)
				{
				// 	skip own value
					if(d2 == d) continue;

					grad[d] *= m_vPolynom[ind[d2]].value(x[d2]);
				}
			}
		}

	private:
		std::vector<Polynomial1D> m_vPolynom;
		std::vector<Polynomial1D> m_vDPolynom;
};


template <>
template <size_t TOrder>
class LagrangeLSFS<ReferencePrism, TOrder>
{
	public:
	///	Reference Element type
		typedef ReferencePrism reference_element_type;

	///	Order of Shape functions
		static const size_t order = TOrder;

	///	Dimension, where shape functions are defined
		static const int dim = reference_element_type::dim;

	private:
	///	abbreviation for order
		static const size_t p = TOrder;

	/// dofs per layer
		static const size_t dofPerLayer = BinomialCoefficient<2 + p, p>::value;

	public:
	///	Domain position type
		typedef MathVector<dim> position_type;

	///	Shape type
		typedef number shape_type;

	///	Gradient type
		typedef MathVector<dim> grad_type;

	/// Number of shape functions
		static const size_t nsh = dofPerLayer*(p+1);

	///	Multi Index type
		typedef MultiIndex<dim> multi_index_type;

	public:
	///	Constructor
		LagrangeLSFS()
		{
			m_vPolynom.resize(p+1);
			m_vDPolynom.resize(p+1);
			m_vTruncPolynom.resize(p+1);
			m_vDTruncPolynom.resize(p+1);

			for(size_t i = 0; i <= p; ++i)
			{
			//	create truncated equidistant polynomials and its derivatives
				m_vTruncPolynom[i] = TruncatedEquidistantLagrange1D(i, p);
				m_vDTruncPolynom[i] = m_vTruncPolynom[i].derivative();

			//	create equidistant polynomials and its derivatives
				m_vPolynom[i] = EquidistantLagrange1D(i, p);
				m_vDPolynom[i] = m_vPolynom[i].derivative();
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		bool position(size_t i, position_type& pos) const
		{
		//	get Multi Index
			MultiIndex<dim> ind = multi_index(i);

		//	set position
			for(int d = 0; d < 2; ++d)
				pos[d] = TruncatedEquidistantLagrange1D::position(ind[d], p);

			pos[2] = EquidistantLagrange1D::position(ind[2], p);

			return true;
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline number shape(const size_t i, const MathVector<dim>& x) const
		{
		//	forward
			return shape(multi_index(i), x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::shapes()
		void shapes(shape_type* sOut, const position_type& x) const
		{
			for(size_t sh = 0; sh < num_sh(); ++sh)
				sOut[sh] = shape(sh, x);
		}

	///	shape value for a Multi Index
		inline number shape(const MultiIndex<dim>& ind, const MathVector<dim>& x) const
		{
			check_multi_index(ind);

		//	get adjoint barycentric index
			const size_t i0 = p - ind[0] - ind[1];
			const number x0 = 1.0 - x[0] - x[1];

				//	x-y direction
			return    m_vTruncPolynom[ ind[0] ].value(x[0])
					* m_vTruncPolynom[ ind[1] ].value(x[1])
					* m_vTruncPolynom[   i0   ].value( x0 )
				//	z direction
					* m_vPolynom[ ind[2] ].value(x[2]);
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


	///	return the index for a multi_index
		inline size_t index(const MultiIndex<dim>& ind) const
		{
			check_multi_index(ind);

			size_t res = ind[0];
			for(size_t i = 0; i < ind[1]; ++i)
				res += (p+1-i);

			// add height
			res += ind[2] * dofPerLayer;

			return res;
		}

	///	return the multi_index for an index
		inline MultiIndex<dim> multi_index(size_t i) const
		{
			UG_ASSERT(i < nsh, "i must be smaller than Number of DoFs.");

			const size_t i2 = i / dofPerLayer;

			int i0 = i - i2*dofPerLayer;
			int i1;
			for(i1 = 0; i1 < (int)p; ++i1)
			{
				const int diff = i0 - (p+1-i1);
				if(diff < 0)
					break;
				i0 = diff;
			}

			return MultiIndex<dim>( i0, i1, i2);
		}

	protected:
		inline void check_multi_index(const MultiIndex<dim>& ind) const
		{
			UG_ASSERT(ind[0] <= p, "Wrong Multiindex.");
			UG_ASSERT(ind[1] <= p, "Wrong Multiindex.");
			UG_ASSERT(ind[0] + ind[1] <= p, "Wrong Multiindex.");
			UG_ASSERT(ind[2] <= p && ind[2] >= 0, "Wrong Multiindex.");
		}

		void eval_grad(grad_type& grad, const size_t i,
		               	   	   	   	    const position_type& x) const
		{
			return eval_grad(grad, multi_index(i), x);
		}

		void eval_grad(grad_type& grad, const MultiIndex<dim> ind,
		               	   	   	   	   	   const position_type& x) const
		{
			check_multi_index(ind);

		//	get adjoint barycentric index and position
			const size_t i0 = p - ind[0] - ind[1];
			const number x0 = 1 - x[0] - x[1];

			UG_ASSERT(i0 <= p && i0 >= 0, "Wrong Multiindex.");
			UG_ASSERT(x0 <= 1.0 && x0 >= 0.0, "Wrong Position.");

		//	x-y gradient
			for(size_t d = 0; d < 2; ++d)
			{
				grad[d] = m_vDTruncPolynom[ind[d]].value(x[d])
						* m_vTruncPolynom[i0].value(x0);
				grad[d] += (-1) * m_vDTruncPolynom[i0].value(x0)
						   * m_vTruncPolynom[ind[d]].value(x[d]);

			//	multiply by all functions not depending on x[d]
				for(size_t d2 = 0; d2 < 2; ++d2)
				{
				// 	skip own value
					if(d2 == d) continue;

					grad[d] *= m_vTruncPolynom[ind[d2]].value(x[d2]);
				}

			//	multiply by z coordinate
				grad[d] *= m_vPolynom[ind[2]].value(x[2]);
			}

		//	z gradient
			grad[2] = m_vDPolynom[ind[2]].value(x[2]);
			grad[2] *=   m_vTruncPolynom[ ind[0] ].value(x[0])
					   * m_vTruncPolynom[ ind[1] ].value(x[1])
					   * m_vTruncPolynom[   i0   ].value( x0 );
		}

	private:
		std::vector<Polynomial1D> m_vPolynom;
		std::vector<Polynomial1D> m_vDPolynom;
		std::vector<Polynomial1D> m_vTruncPolynom;
		std::vector<Polynomial1D> m_vDTruncPolynom;
};

namespace {

template <int p>
struct NumberOfDoFsOfPyramid
{
	enum
	{
		value = NumberOfDoFsOfPyramid<p-1>::value
				+ (p+1)*(p+1)
	};
};

// specialization to end recursion
template <>
struct NumberOfDoFsOfPyramid<1>
{
	enum
	{
		value = 5
	};
};

} // end empty namespace

// todo: Implement
template <>
template <size_t TOrder>
class LagrangeLSFS<ReferencePyramid, TOrder>
{
	private:
	///	abbreviation for order
		static const size_t p = TOrder;

	public:
	///	Reference Element type
		typedef ReferencePyramid reference_element_type;

	///	Order of Shape functions
		static const size_t order = TOrder;

	///	Dimension, where shape functions are defined
		static const int dim = reference_element_type::dim;

	///	Domain position type
		typedef MathVector<dim> position_type;

	///	Shape type
		typedef number shape_type;

	///	Gradient type
		typedef MathVector<dim> grad_type;

	/// Number of shape functions
		static const size_t nsh = NumberOfDoFsOfPyramid<p>::value;

	///	Multi Index type
		typedef MultiIndex<dim> multi_index_type;

	public:
	///	Constructor
		LagrangeLSFS()
		{
			m_vvPolynom.resize(p+1);
			m_vvDPolynom.resize(p+1);

			for(size_t i2 = 0; i2 <= p; ++i2)
			{
				m_vvPolynom[i2].resize(p+1);
				m_vvDPolynom[i2].resize(p+1);

				for(size_t i = 0; i <= p-i2; ++i)
				{
					m_vvPolynom[i2][i] = BoundedEquidistantLagrange1D(i, p, p-i2);
					m_vvDPolynom[i2][i] = m_vvPolynom[i2][i].derivative();
				}
			}
			UG_LOG("Pol (0,0) = " << m_vvPolynom[0][0] << "\n");
			UG_LOG("Pol (0,1) = " << m_vvPolynom[0][1] << "\n");
			UG_LOG("Pol (1,0) = " << m_vvPolynom[1][0] << "\n");
		}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		bool position(size_t i, position_type& pos) const
		{
		//	get Multi Index
			MultiIndex<dim> ind = multi_index(i);

		//	set position
			for(int d = 0; d < dim; ++d)
				pos[d] = EquidistantLagrange1D::position(ind[d], p);

			return true;
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline number shape(const size_t i, const MathVector<dim>& x) const
		{
		//	forward
			return shape(multi_index(i), x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::shapes()
		void shapes(shape_type* sOut, const position_type& x) const
		{
			for(size_t sh = 0; sh < num_sh(); ++sh)
				sOut[sh] = shape(sh, x);
		}

	///	shape value for a Multi Index
		inline number shape(const MultiIndex<dim>& ind, const MathVector<dim>& x) const
		{
			check_multi_index(ind);

			if(ind[2] == 0 && ind[0]==0 && ind[1] == 0)
				return 	  m_vvPolynom[ 0 ][ ind[0] ].value(x[0])
						* m_vvPolynom[ 0 ][ ind[1] ].value(x[1])
						- m_vvPolynom[ 0 ][ 1 ].value(x[2]);
			else if(ind[2] == 0)
				return    m_vvPolynom[ 0 ][ ind[0] ].value(x[0])
						* m_vvPolynom[ 0 ][ ind[1] ].value(x[1]);
			else
				return m_vvPolynom[ 0 ][ 1 ].value(x[2]);
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


	///	return the index for a multi_index
		inline size_t index(const MultiIndex<dim>& ind) const
		{
			check_multi_index(ind);

			size_t res = 0;

		//	add layers that are completely filled
			for(size_t i2 = 0; i2 < ind[2]; ++i2)
				res += (p+1-i2)*(p+1-i2);

		//	add dofs of top z-layer
			res += ind[1] * (p+1-ind[2]) + ind[0];

			UG_ASSERT(res < nsh, "Wrong index computation. Internal error.");

			return res;
		}

	///	return the multi_index for an index
		inline MultiIndex<dim> multi_index(size_t i) const
		{
			UG_ASSERT(i < nsh, "i must be smaller than Number of DoFs.");

		//	get z layer
			int iTmp = i;
			int i2;
			for(i2 = 0; i2 < (int)p; ++i2)
			{
				const int diff = iTmp - (p+1-i2)*(p+1-i2);
				if(diff < 0) break;

				iTmp = diff;
			}

			return MultiIndex<dim>( iTmp%(p+1-i2), iTmp/(p+1-i2), i2);
		}

	protected:
		inline void check_multi_index(const MultiIndex<dim>& ind) const
		{
			UG_ASSERT(ind[0] <= p-ind[2] && ind[0] >= 0, "Wrong Multiindex.");
			UG_ASSERT(ind[1] <= p-ind[2] && ind[1] >= 0, "Wrong Multiindex.");
			UG_ASSERT(ind[2] <= p && ind[2] >= 0, "Wrong Multiindex.");
		}

		void eval_grad(grad_type& grad, const size_t i,
		               	   	   	   	    const position_type& x) const
		{
			return eval_grad(grad, multi_index(i), x);
		}

		void eval_grad(grad_type& grad, const MultiIndex<dim> ind,
		               	   	   	   	   	   const position_type& x) const
		{
			check_multi_index(ind);

		//	loop x,y
			if(ind[2] == 0)
			{
				for(int d = 0; d < 2; ++d)
				{
					grad[d] = m_vvDPolynom[ 0 ][ ind[d] ].value(x[d]);

				//	multiply by all functions not depending on x[d]
					for(int d2 = 0; d2 < 2; ++d2)
					{
					// 	skip own value
						if(d2 == d) continue;

						grad[d] *= m_vvPolynom[ 0 ][ ind[d2] ].value(x[d2]);
					}

				//	multiply by z factor
					//grad[d] *= m_vvPolynom[ 0 ][ ind[2] ].value(x[2]);
				}

				if(ind[2] == 0 && ind[0]==0 && ind[1] == 0)
					grad[2] = - m_vvDPolynom[ 0 ][ 1 ].value(x[2]);
				else
					grad[2] = 0.0;
			}
			else
			{
				grad[0] = 0.0; grad[1] = 0.0; grad[2] = 1.0; return;
			}
		//	do z
/*			for(int d2 = 0; d2 < 2; ++d2)
			{
				grad[2] *= m_vvPolynom[ ind[2] ][ ind[d2] ].value(x[d2]);
			}
*/		}

	private:
		std::vector<std::vector<Polynomial1D> > m_vvPolynom;
		std::vector<std::vector<Polynomial1D> > m_vvDPolynom;
};


template <>
template <size_t TOrder>
class LagrangeLSFS<ReferenceHexahedron, TOrder>
{
	private:
	///	abbreviation for order
		static const size_t p = TOrder;

	public:
	///	Reference Element type
		typedef ReferenceHexahedron reference_element_type;

	///	Order of Shape functions
		static const size_t order = TOrder;

	///	Dimension, where shape functions are defined
		static const int dim = reference_element_type::dim;

	///	Domain position type
		typedef MathVector<dim> position_type;

	///	Shape type
		typedef number shape_type;

	///	Gradient type
		typedef MathVector<dim> grad_type;

	/// Number of shape functions
		static const size_t nsh = (p+1)*(p+1)*(p+1);

	///	Multi Index type
		typedef MultiIndex<dim> multi_index_type;

	public:
	///	Constructor
		LagrangeLSFS()
		{
			m_vPolynom.resize(p+1);
			m_vDPolynom.resize(p+1);

			for(size_t i = 0; i <= p; ++i)
			{
			//	create trancated equidistant polynomials and its derivatives
				m_vPolynom[i] = EquidistantLagrange1D(i, p);
				m_vDPolynom[i] = m_vPolynom[i].derivative();
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		bool position(size_t i, position_type& pos) const
		{
		//	get Multi Index
			MultiIndex<dim> ind = multi_index(i);

		//	set position
			for(int d = 0; d < dim; ++d)
				pos[d] = EquidistantLagrange1D::position(ind[d], p);

			return true;
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline number shape(const size_t i, const MathVector<dim>& x) const
		{
		//	forward
			return shape(multi_index(i), x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::shapes()
		void shapes(shape_type* sOut, const position_type& x) const
		{
			for(size_t sh = 0; sh < num_sh(); ++sh)
				sOut[sh] = shape(sh, x);
		}

	///	shape value for a Multi Index
		inline number shape(const MultiIndex<dim>& ind, const MathVector<dim>& x) const
		{
			check_multi_index(ind);

			return 	  m_vPolynom[ ind[0] ].value(x[0])
					* m_vPolynom[ ind[1] ].value(x[1])
					* m_vPolynom[ ind[2] ].value(x[2]);
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


	///	return the index for a multi_index
		inline size_t index(const MultiIndex<dim>& ind) const
		{
			check_multi_index(ind);

			return ind[2] * (p+1)*(p+1) + ind[1] * (p+1) + ind[0];
		}

	///	return the multi_index for an index
		inline MultiIndex<dim> multi_index(size_t i) const
		{
			UG_ASSERT(i < nsh, "i must be smaller than Number of DoFs.");

			return MultiIndex<dim>( i%(p+1), i/(p+1)%(p+1), i/((p+1)*(p+1)));
		}

	protected:
		inline void check_multi_index(const MultiIndex<dim>& ind) const
		{
			UG_ASSERT(ind[0] <= p && ind[0] >= 0, "Wrong Multiindex.");
			UG_ASSERT(ind[1] <= p && ind[1] >= 0, "Wrong Multiindex.");
			UG_ASSERT(ind[2] <= p && ind[2] >= 0, "Wrong Multiindex.");
		}

		void eval_grad(grad_type& grad, const size_t i,
		               	   	   	   	    const position_type& x) const
		{
			return eval_grad(grad, multi_index(i), x);
		}

		void eval_grad(grad_type& grad, const MultiIndex<dim> ind,
		               	   	   	   	   	   const position_type& x) const
		{
			check_multi_index(ind);

		//	loop dimensions
			for(int d = 0; d < dim; ++d)
			{
				grad[d] = m_vDPolynom[ind[d]].value(x[d]);

			//	multiply by all functions not depending on x[d]
				for(int d2 = 0; d2 < dim; ++d2)
				{
				// 	skip own value
					if(d2 == d) continue;

					grad[d] *= m_vPolynom[ind[d2]].value(x[d2]);
				}
			}
		}

	private:
		std::vector<Polynomial1D> m_vPolynom;
		std::vector<Polynomial1D> m_vDPolynom;
};


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

/// wrapper class implementing the LocalShapeFunctionSet interface
/**
 * This class wrappes a class passed by the template arguement into the
 * virtual ILocalShapeFunctionSet interface and makes it thus usable in that
 * context on the price of virtual functions.
 * \tparam 	TImpl		Implementation of a Local Shape Function Set
 */
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
		static const size_t order = ImplType::order;

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
			ImplType::shapes(sOut, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		virtual grad_type grad(size_t i, const position_type& x) const
		{
			return ImplType::grad(i, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grads()
		virtual void grads(grad_type* gOut, const position_type& x) const
		{
			ImplType::grads(gOut, x);
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
