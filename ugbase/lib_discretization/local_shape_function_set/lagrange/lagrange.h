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
#include "../local_dof.h"
#include "lib_discretization/common/multi_index.h"
#include "common/util/metaprogramming_util.h"
#include "lib_grid/grid/geometric_base_objects.h"

namespace ug{

/// Lagrange Shape Function Set without virtual functions
template <typename TRefElem, int TOrder>
struct LagrangeLSFS{};

///////////////////////////////////////////////////////////////////////////////
// Edge
///////////////////////////////////////////////////////////////////////////////

/// specialization for Edges
/**
 * Lagrange shape function of any order for the Reference Edge
 * \tparam 	TOrder		requested order
 */
template <>
template <int TOrder>
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
			for(size_t i = 0; i < nsh; ++i)
			{
			//	create equidistant polynomials and its derivatives
				m_vPolynom[i] = EquidistantLagrange1D(i, p);
				m_vDPolynom[i] = m_vPolynom[i].derivative();
			}
		}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline static size_t num_sh() {return nsh;}

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
			grad_type tmpGrad; grad(tmpGrad, i, x); return tmpGrad;
		}

	///	\copydoc ug::LocalShapeFunctionSet::grads()
		void grads(grad_type* gOut, const position_type& x) const
		{
			for(size_t sh = 0; sh < num_sh(); ++sh)
				grad(gOut[sh], sh, x);
		}

	///	evaluates the gradient
		void grad(grad_type& g, size_t i, const position_type& x) const
		{
			g[0] = m_vDPolynom[i].value(x[0]);
		}

	///	return Multi index for index i
		inline static MultiIndex<dim> multi_index(size_t i)
		{
			check_index(i);
			return MultiIndex<1>(i);
		}

	///	return the index for a multi_index
		inline static size_t index(const MultiIndex<dim>& ind)
		{
			check_multi_index(ind);
			return ind[0];
		}

	///	returns if DoFs are assigned to geometric objects of the dimension
		bool has_sh_on(int d) const
		{
			return (d==0) || (d==1 && p>1);
		}

	///	returns the number of DoFs on a given sub-geometric object
		size_t num_sh(int d, size_t id) const
		{
			if(d==0) return 1;
			if(d==1) return (p-1);
			else return 0;
		}

	///	returns the number of DoFs on a sub-geometric object type
		size_t num_sh(ReferenceObjectID type) const
		{
			if(type == ROID_VERTEX) return 1;
			if(type == ROID_EDGE) return (p-1);
			else return 0;
		}

	///	checks in debug mode that index is valid
		inline static void check_index(size_t i)
		{
			UG_ASSERT(i < nsh, "i must be smaller than Number of DoFs.");
		}

	///	checks in debug mode that multi-index is valid
		inline static void check_multi_index(const MultiIndex<dim>& ind)
		{
			UG_ASSERT(ind[0] < nsh, "ind[0] must be smaller than Number of DoFs.");
		}

	protected:
		Polynomial1D m_vPolynom[p+1];	///< Shape Polynomials
		Polynomial1D m_vDPolynom[p+1];	///< Derivative of Shape Polynomial
};

///////////////////////////////////////////////////////////////////////////////
// Triangle
///////////////////////////////////////////////////////////////////////////////

template <>
template <int TOrder>
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
			ReferenceTriangle::check_position(x);

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
			grad_type tmpGrad; grad(tmpGrad, i, x); return tmpGrad;
		}

	///	\copydoc ug::LocalShapeFunctionSet::grads()
		void grads(grad_type* gOut, const position_type& x) const
		{
			for(size_t sh = 0; sh < num_sh(); ++sh)
				grad(gOut[sh], sh, x);
		}

	///	evaluates the gradient
		void grad(grad_type& g, const size_t i,	const position_type& x) const
		{
			grad(g, multi_index(i), x);
		}

	///	evaluates the gradient
		void grad(grad_type& g, const MultiIndex<dim> ind,
		               	   	   	   	   	   const position_type& x) const
		{
			check_multi_index(ind);
			ReferenceTriangle::check_position(x);

		//	get adjoint barycentric index and position
			const size_t i0 = p - ind[0] - ind[1];
			const number x0 = 1.0 - x[0] - x[1];

			UG_ASSERT(i0 <= p && i0 >= 0, "Wrong Multiindex.");
			UG_ASSERT(x0 <= 1.0 && x0 >= 0.0, "Wrong Position.");

		//	loop dimensions
			for(int d = 0; d < dim; ++d)
			{
				g[d] = m_vDPolynom[ind[d]].value(x[d])
						* m_vPolynom[i0].value(x0);
				g[d] += (-1) * m_vDPolynom[i0].value(x0)
						   * m_vPolynom[ind[d]].value(x[d]);

			//	multiply by all functions not depending on x[d]
				for(int d2 = 0; d2 < dim; ++d2)
				{
				// 	skip own value
					if(d2 == d) continue;

					g[d] *= m_vPolynom[ind[d2]].value(x[d2]);
				}
			}
		}

	///	return the index for a multi_index
		inline size_t index(const MultiIndex<dim>& ind) const
		{
			check_multi_index(ind);

			size_t res = ind[0];
			for(size_t i = 0; i < ind[1]; ++i)
				res += (p+1-i);

			check_index(res);
			return res;
		}

	///	return the multi_index for an index
		inline MultiIndex<dim> multi_index(size_t i) const
		{
			check_index(i);

			int i0 = i, i1;
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

	///	returns if DoFs are assigned to geometric objects of the dimension
		bool has_sh_on(int d) const
		{
			return (d==0) || (d==1 && p>1) || (d==2 && p>2);
		}

	///	returns the number of DoFs on a given sub-geometric object
		size_t num_sh(int d, size_t id) const
		{
			if(d==0) return 1;
			if(d==1) return (p-1);
			if(d==2) return ((p>2) ? (BinomialCoefficient<dim + p-3, p-3>::value) : 0);
			else return 0;
		}

	///	returns the number of DoFs on a sub-geometric object type
		size_t num_sh(ReferenceObjectID type) const
		{
			if(type == ROID_VERTEX) return 1;
			if(type == ROID_EDGE) return (p-1);
			if(type == ROID_TRIANGLE)
				return ((p>2) ? (BinomialCoefficient<dim + p-3, p-3>::value) : 0);
			else return 0;
		}

	///	checks in debug mode that index is valid
		inline static void check_index(size_t i)
		{
			UG_ASSERT(i < nsh, "i must be smaller than Number of DoFs.");
		}

	///	checks in debug mode that multi-index is valid
		inline static void check_multi_index(const MultiIndex<dim>& ind)
		{
			UG_ASSERT(ind[0] <= p, "Wrong Multiindex.");
			UG_ASSERT(ind[1] <= p, "Wrong Multiindex.");
			UG_ASSERT(ind[0] + ind[1] <= p, "Wrong Multiindex.");
		}

	private:
		Polynomial1D m_vPolynom[p+1];
		Polynomial1D m_vDPolynom[p+1];
};

///////////////////////////////////////////////////////////////////////////////
// Quadrilateral
///////////////////////////////////////////////////////////////////////////////

template <>
template <int TOrder>
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
			ReferenceQuadrilateral::check_position(x);

			return    m_vPolynom[ ind[0] ].value(x[0])
					* m_vPolynom[ ind[1] ].value(x[1]);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		grad_type grad(size_t i, const position_type& x) const
		{
			grad_type tmpGrad; grad(tmpGrad, i, x); return tmpGrad;
		}

	///	\copydoc ug::LocalShapeFunctionSet::grads()
		void grads(grad_type* gOut, const position_type& x) const
		{
			for(size_t sh = 0; sh < num_sh(); ++sh)
				grad(gOut[sh], sh, x);
		}

	/// evaluates the gradient
		void grad(grad_type& g, const size_t i, const position_type& x) const
		{
			grad(g, multi_index(i), x);
		}

	/// evaluates the gradient
		void grad(grad_type& g, const MultiIndex<dim> ind,
		               	   	   	const position_type& x) const
		{
			check_multi_index(ind);
			ReferenceQuadrilateral::check_position(x);

		//	loop dimensions
			for(int d = 0; d < dim; ++d)
			{
				g[d] = m_vDPolynom[ind[d]].value(x[d]);

			//	multiply by all functions not depending on x[d]
				for(int d2 = 0; d2 < dim; ++d2)
				{
				// 	skip own value
					if(d2 == d) continue;

					g[d] *= m_vPolynom[ind[d2]].value(x[d2]);
				}
			}
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
			check_index(i);

			return MultiIndex<dim>( i%(p+1), i/(p+1) );
		}

	///	returns if DoFs are assigned to geometric objects of the dimension
		bool has_sh_on(int d) const
		{
			return (d==0) || (d==1 && p>1) || (d==2 && p>1);
		}

	///	returns the number of DoFs on a given sub-geometric object
		size_t num_sh(int d, size_t id) const
		{
			if(d==0) return 1;
			if(d==1) return (p-1);
			if(d==2) return (p-1)*(p-1);
			else return 0;
		}

	///	returns the number of DoFs on a sub-geometric object type
		size_t num_sh(ReferenceObjectID type) const
		{
			if(type == ROID_VERTEX) return 1;
			if(type == ROID_EDGE) return (p-1);
			if(type == ROID_QUADRILATERAL) return (p-1)*(p-1);
			else return 0;
		}

	///	checks in debug mode that index is valid
		inline static void check_index(size_t i)
		{
			UG_ASSERT(i < nsh, "i must be smaller than Number of DoFs.");
		}

	///	checks in debug mode that multi-index is valid
		inline static void check_multi_index(const MultiIndex<dim>& ind)
		{
			UG_ASSERT(ind[0] <= p && ind[0] >= 0, "Wrong Multiindex.");
			UG_ASSERT(ind[1] <= p && ind[1] >= 0, "Wrong Multiindex.");
		}

	private:
		Polynomial1D m_vPolynom[p+1];
		Polynomial1D m_vDPolynom[p+1];
};

///////////////////////////////////////////////////////////////////////////////
// Tetrahedron
///////////////////////////////////////////////////////////////////////////////

template <>
template <int TOrder>
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
			ReferenceTetrahedron::check_position(x);

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
			grad_type tmpGrad; grad(tmpGrad, i, x); return tmpGrad;
		}

	///	\copydoc ug::LocalShapeFunctionSet::grads()
		void grads(grad_type* gOut, const position_type& x) const
		{
			for(size_t sh = 0; sh < num_sh(); ++sh)
				grad(gOut[sh], sh, x);
		}

	///	evaluates the gradient
		void grad(grad_type& g, const size_t i, const position_type& x) const
		{
			grad(g, multi_index(i), x);
		}

	///	evaluates the gradient
		void grad(grad_type& g, const MultiIndex<dim> ind,
		               	   	   	   	    const position_type& x) const
		{
			check_multi_index(ind);
			ReferenceTetrahedron::check_position(x);

		//	get adjoint barycentric index and position
			const size_t i0 = p - ind[0] - ind[1] - ind[2];
			const number x0 = 1 - x[0] - x[1] - x[2];

			UG_ASSERT(i0 <= p && i0 >= 0, "Wrong Multiindex.");
			UG_ASSERT(x0 <= 1.0 && x0 >= 0.0, "Wrong Position.");

		//	loop dimensions
			for(int d = 0; d < dim; ++d)
			{
				g[d] = m_vDPolynom[ind[d]].value(x[d])
						* m_vPolynom[i0].value(x0);
				g[d] += (-1) * m_vDPolynom[i0].value(x0)
						   * m_vPolynom[ind[d]].value(x[d]);

			//	multiply by all functions not depending on x[d]
				for(int d2 = 0; d2 < dim; ++d2)
				{
				// 	skip own value
					if(d2 == d) continue;

					g[d] *= m_vPolynom[ind[d2]].value(x[d2]);
				}
			}
		}

	///	return the index for a multi_index
		inline size_t index(const MultiIndex<dim>& ind) const
		{
			check_multi_index(ind);

		//	add x range
			size_t res = ind[0];

		//	add y range
			for(size_t i = 0; i < ind[1]; ++i)
				res += (p+1-ind[2]-i);

		//	add z range
			for(size_t i2 = 0; i2 < ind[2]; ++i2)
				res += BinomCoeff(p+2-i2, p-i2);

			check_index(res);
			return res;
		}

	///	return the multi_index for an index
		inline MultiIndex<dim> multi_index(size_t i) const
		{
			check_index(i);

			int i0 = i, i1 = 0, i2 = 0;
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

		///	returns if DoFs are assigned to geometric objects of the dimension
			bool has_sh_on(int d) const
			{
				return (d==0) || (d==1 && p>1) || (d==2 && p>2) || (d==3 && p>3);
			}

		///	returns the number of DoFs on a given sub-geometric object
			size_t num_sh(int d, size_t id) const
			{
				if(d==0) return 1;
			//	number of shapes on edge is same as for edge with p-1
				if(d==1) return (p-1);
			//	number of shapes on faces is same as for triangles in 2d
				if(d==2) return ((p>2) ? (BinomialCoefficient<dim-1 + p-3, p-3>::value) : 0);
			//	number of shapes on interior is same as for tetrahedra with p-4
				if(d==3) return ((p>3) ? (BinomialCoefficient<dim + p-4, p-4>::value) : 0);
				else return 0;
			}

		///	returns the number of DoFs on a sub-geometric object type
			size_t num_sh(ReferenceObjectID type) const
			{
				if(type == ROID_VERTEX) return 1;
				if(type == ROID_EDGE) return (p-1);
				if(type == ROID_TRIANGLE)
					return ((p>2) ? (BinomialCoefficient<dim-1 + p-3, p-3>::value) : 0);
				if(type == ROID_TETRAHEDRON)
					return ((p>3) ? (BinomialCoefficient<dim + p-4, p-4>::value) : 0);
				else return 0;
			}

	///	checks in debug mode that index is valid
		inline static void check_index(size_t i)
		{
			UG_ASSERT(i < nsh, "i must be smaller than Number of DoFs.");
		}

	///	checks in debug mode that multi-index is valid
		inline static void check_multi_index(const MultiIndex<dim>& ind)
		{
			UG_ASSERT(ind[0] <= p, "Wrong Multiindex.");
			UG_ASSERT(ind[1] <= p, "Wrong Multiindex.");
			UG_ASSERT(ind[2] <= p, "Wrong Multiindex.");
			UG_ASSERT(ind[0] + ind[1] + ind[2] <= p, "Wrong Multiindex.");
		}

	private:
		Polynomial1D m_vPolynom[p+1];
		Polynomial1D m_vDPolynom[p+1];
};

///////////////////////////////////////////////////////////////////////////////
// Prism
///////////////////////////////////////////////////////////////////////////////

template <>
template <int TOrder>
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
			ReferencePrism::check_position(x);

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
			grad_type tmpGrad; grad(tmpGrad, i, x); return tmpGrad;
		}

	///	\copydoc ug::LocalShapeFunctionSet::grads()
		void grads(grad_type* gOut, const position_type& x) const
		{
			for(size_t sh = 0; sh < num_sh(); ++sh)
				grad(gOut[sh], sh, x);
		}

	///	evaluates the gradient
		void grad(grad_type& g, const size_t i, const position_type& x) const
		{
			grad(g, multi_index(i), x);
		}

	///	evaluates the gradient
		void grad(grad_type& g, const MultiIndex<dim> ind,
		               	   	   	   	   	const position_type& x) const
		{
			check_multi_index(ind);
			ReferencePrism::check_position(x);

		//	get adjoint barycentric index and position
			const size_t i0 = p - ind[0] - ind[1];
			const number x0 = 1 - x[0] - x[1];

			UG_ASSERT(i0 <= p && i0 >= 0, "Wrong Multiindex.");
			UG_ASSERT(x0 <= 1.0 && x0 >= 0.0, "Wrong Position.");

		//	x-y gradient
			for(size_t d = 0; d < 2; ++d)
			{
				g[d] = m_vDTruncPolynom[ind[d]].value(x[d])
						* m_vTruncPolynom[i0].value(x0);
				g[d] += (-1) * m_vDTruncPolynom[i0].value(x0)
						   * m_vTruncPolynom[ind[d]].value(x[d]);

			//	multiply by all functions not depending on x[d]
				for(size_t d2 = 0; d2 < 2; ++d2)
				{
				// 	skip own value
					if(d2 == d) continue;

					g[d] *= m_vTruncPolynom[ind[d2]].value(x[d2]);
				}

			//	multiply by z coordinate
				g[d] *= m_vPolynom[ind[2]].value(x[2]);
			}

		//	z gradient
			g[2] = m_vDPolynom[ind[2]].value(x[2]);
			g[2] *=   m_vTruncPolynom[ ind[0] ].value(x[0])
					   * m_vTruncPolynom[ ind[1] ].value(x[1])
					   * m_vTruncPolynom[   i0   ].value( x0 );
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
			check_index(i);

			const size_t i2 = i / dofPerLayer;

			int i0 = i - i2*dofPerLayer, i1;
			for(i1 = 0; i1 < (int)p; ++i1)
			{
				const int diff = i0 - (p+1-i1);
				if(diff < 0)
					break;
				i0 = diff;
			}

			return MultiIndex<dim>( i0, i1, i2);
		}

	///	returns if DoFs are assigned to geometric objects of the dimension
		bool has_sh_on(int d) const
		{
			return (d==0) || (d==1 && p>1) || (d==2 && p>1) || (d==3 && p >2);
		}

	///	returns the number of DoFs on a given sub-geometric object
		size_t num_sh(int d, size_t id) const
		{
		//	get reference element
			static const ReferencePrism& rRef
						= ReferenceElementProvider::get<ReferencePrism>();

		//	get type of subelement
			const ReferenceObjectID type = rRef.ref_elem_type(d, id);

		//	forward request
			return num_sh(type);
		}

	///	returns the number of DoFs on a sub-geometric object type
		size_t num_sh(ReferenceObjectID type) const
		{
			if(type == ROID_VERTEX) return 1;
			if(type == ROID_EDGE) return (p-1);
		//	same as for a 2d triangle of order p-3
			if(type == ROID_TRIANGLE)
				return ((p>2) ? (BinomialCoefficient<dim + p-3, p-3>::value) : 0);
		//	same as for a 2d quadrilateral of order p-2
			if(type == ROID_QUADRILATERAL)
				return ((p>1) ? ((p-1)*(p-1)) : 0);
		//	same as for a 3d prism of order p-2
			if(type == ROID_PRISM)
				return ((p>2) ? (LagrangeLSFS<ReferencePrism, p-2>::nsh) : 0);
			else return 0;
		}

	///	checks in debug mode that index is valid
		inline static void check_index(size_t i)
		{
			UG_ASSERT(i < nsh, "i must be smaller than Number of DoFs.");
		}

	///	checks in debug mode that multi-index is valid
		inline static void check_multi_index(const MultiIndex<dim>& ind)
		{
			UG_ASSERT(ind[0] <= p, "Wrong Multiindex.");
			UG_ASSERT(ind[1] <= p, "Wrong Multiindex.");
			UG_ASSERT(ind[0] + ind[1] <= p, "Wrong Multiindex.");
			UG_ASSERT(ind[2] <= p && ind[2] >= 0, "Wrong Multiindex.");
		}

	private:
		Polynomial1D m_vPolynom[p+1];
		Polynomial1D m_vDPolynom[p+1];
		Polynomial1D m_vTruncPolynom[p+1];
		Polynomial1D m_vDTruncPolynom[p+1];
};

///////////////////////////////////////////////////////////////////////////////
// Pyramid
///////////////////////////////////////////////////////////////////////////////

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
	enum {value = 5};
};
template <>
struct NumberOfDoFsOfPyramid<0>
{
	enum {value = 0};
};
template <>
struct NumberOfDoFsOfPyramid<-1>
{
	enum {value = 0};
};

} // end empty namespace

// todo: Implement higher order (impossible ?)
//	NOTE:	Currently only 1st order is implemented. There is no shape function
//			set for pyramids, that is continuous and allows a continuous
//			derivative in the inner of the pyramid. This is basically, since
//			one reguards the pyramid as two tetrahedrons, glued together.
template <>
template <int TOrder>
class LagrangeLSFS<ReferencePyramid, TOrder>
{
	private:
	///	abbreviation for order
		static const size_t p = 1;

	public:
	///	Reference Element type
		typedef ReferencePyramid reference_element_type;

	///	Order of Shape functions
		static const size_t order = 1;

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
			ReferencePyramid::check_position(x);

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
			grad_type tmpGrad; grad(tmpGrad, i, x); return tmpGrad;
		}

	///	\copydoc ug::LocalShapeFunctionSet::grads()
		void grads(grad_type* gOut, const position_type& x) const
		{
			for(size_t sh = 0; sh < num_sh(); ++sh)
				grad(gOut[sh], sh, x);
		}

	///	evaluates the gradient
		void grad(grad_type& g, const size_t i, const position_type& x) const
		{
			grad(g, multi_index(i), x);
		}

	///	evaluates the gradient
		void grad(grad_type& g, const MultiIndex<dim> ind,
		               	   	   	   	   	   const position_type& x) const
		{
			check_multi_index(ind);
			ReferencePyramid::check_position(x);

		//	loop x,y
			if(ind[2] == 0)
			{
				for(int d = 0; d < 2; ++d)
				{
					g[d] = m_vvDPolynom[ 0 ][ ind[d] ].value(x[d]);

				//	multiply by all functions not depending on x[d]
					for(int d2 = 0; d2 < 2; ++d2)
					{
					// 	skip own value
						if(d2 == d) continue;

						g[d] *= m_vvPolynom[ 0 ][ ind[d2] ].value(x[d2]);
					}

				//	multiply by z factor
					//grad[d] *= m_vvPolynom[ 0 ][ ind[2] ].value(x[2]);
				}

				if(ind[2] == 0 && ind[0]==0 && ind[1] == 0)
					g[2] = - m_vvDPolynom[ 0 ][ 1 ].value(x[2]);
				else
					g[2] = 0.0;
			}
			else
			{
				g[0] = 0.0; g[1] = 0.0; g[2] = 1.0; return;
			}
		//	do z
/*			for(int d2 = 0; d2 < 2; ++d2)
			{
				grad[2] *= m_vvPolynom[ ind[2] ][ ind[d2] ].value(x[d2]);
			}
*/		}

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

			check_index(res);
			return res;
		}

	///	return the multi_index for an index
		inline MultiIndex<dim> multi_index(size_t i) const
		{
			check_index(i);

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

	///	returns if DoFs are assigned to geometric objects of the dimension
		bool has_sh_on(int d) const
		{
			return (d==0) || (d==1 && p>1) || (d==2 && p>1) || (d==3 && p >2);
		}

	///	returns the number of DoFs on a given sub-geometric object
		size_t num_sh(int d, size_t id) const
		{
		//	get reference element
			static const ReferencePrism& rRef
						= ReferenceElementProvider::get<ReferencePrism>();

		//	get type of subelement
			const ReferenceObjectID type = rRef.ref_elem_type(d, id);

		//	forward request
			return num_sh(type);
		}

	///	returns the number of DoFs on a sub-geometric object type
		size_t num_sh(ReferenceObjectID type) const
		{
			if(type == ROID_VERTEX) return 1;
			if(type == ROID_EDGE) return (p-1);
		//	same as for a 2d triangle of order p-3
			if(type == ROID_TRIANGLE)
				return ((p>2) ? (BinomialCoefficient<dim + p-3, p-3>::value) : 0);
		//	same as for a 2d quadrilateral of order p-2
			if(type == ROID_QUADRILATERAL)
				return ((p>1) ? ((p-1)*(p-1)) : 0);
		//	same as for a 3d pyramid of order p-2
			if(type == ROID_PYRAMID)
				return ((p>2) ? (NumberOfDoFsOfPyramid<p-2>::value) : 0);
			else return 0;
		}

	///	checks in debug mode that index is valid
		inline static void check_index(size_t i)
		{
			UG_ASSERT(i < nsh, "i must be smaller than Number of DoFs.");
		}

	///	checks in debug mode that multi-index is valid
		inline static void check_multi_index(const MultiIndex<dim>& ind)
		{
			UG_ASSERT(ind[0] <= p-ind[2] && ind[0] >= 0, "Wrong Multiindex.");
			UG_ASSERT(ind[1] <= p-ind[2] && ind[1] >= 0, "Wrong Multiindex.");
			UG_ASSERT(ind[2] <= p && ind[2] >= 0, "Wrong Multiindex.");
		}

	private:
		std::vector<std::vector<Polynomial1D> > m_vvPolynom;
		std::vector<std::vector<Polynomial1D> > m_vvDPolynom;
};

///////////////////////////////////////////////////////////////////////////////
// Hexahedron
///////////////////////////////////////////////////////////////////////////////

template <>
template <int TOrder>
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
			ReferenceHexahedron::check_position(x);

			return 	  m_vPolynom[ ind[0] ].value(x[0])
					* m_vPolynom[ ind[1] ].value(x[1])
					* m_vPolynom[ ind[2] ].value(x[2]);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		grad_type grad(size_t i, const position_type& x) const
		{
			grad_type tmpGrad; grad(tmpGrad, i, x); return tmpGrad;
		}

	///	\copydoc ug::LocalShapeFunctionSet::grads()
		void grads(grad_type* gOut, const position_type& x) const
		{
			for(size_t sh = 0; sh < num_sh(); ++sh)
				grad(gOut[sh], sh, x);
		}

	///	evaluates the gradient
		void grad(grad_type& g, const size_t i, const position_type& x) const
		{
			grad(g, multi_index(i), x);
		}

	///	evaluates the gradient
		void grad(grad_type& g, const MultiIndex<dim> ind,
		          	  	  	  	const position_type& x) const
		{
			check_multi_index(ind);
			ReferenceHexahedron::check_position(x);

		//	loop dimensions
			for(int d = 0; d < dim; ++d)
			{
				g[d] = m_vDPolynom[ind[d]].value(x[d]);

			//	multiply by all functions not depending on x[d]
				for(int d2 = 0; d2 < dim; ++d2)
				{
				// 	skip own value
					if(d2 == d) continue;

					g[d] *= m_vPolynom[ind[d2]].value(x[d2]);
				}
			}
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
			check_index(i);

			return MultiIndex<dim>( i%(p+1), i/(p+1)%(p+1), i/((p+1)*(p+1)));
		}

	///	returns if DoFs are assigned to geometric objects of the dimension
		bool has_sh_on(int d) const
		{
			return (d==0) || (d==1 && p>1) || (d==2 && p>1) || (d==3 && p>1);
		}

	///	returns the number of DoFs on a given sub-geometric object
		size_t num_sh(int d, size_t id) const
		{
			if(d==0) return 1;
			if(d==1) return (p-1);
			if(d==2) return (p-1)*(p-1);
			if(d==3) return (p-1)*(p-1)*(p-1);
			else return 0;
		}

	///	returns the number of DoFs on a sub-geometric object type
		size_t num_sh(ReferenceObjectID type) const
		{
			if(type == ROID_VERTEX) return 1;
			if(type == ROID_EDGE) return (p-1);
			if(type == ROID_QUADRILATERAL) return (p-1)*(p-1);
			if(type == ROID_HEXAHEDRON) return (p-1)*(p-1)*(p-1);
			else return 0;
		}

	///	checks in debug mode that index is valid
		inline static void check_index(size_t i)
		{
			UG_ASSERT(i < nsh, "i must be smaller than Number of DoFs.");
		}

	///	checks in debug mode that multi-index is valid
		inline static void check_multi_index(const MultiIndex<dim>& ind)
		{
			UG_ASSERT(ind[0] <= p && ind[0] >= 0, "Wrong Multiindex.");
			UG_ASSERT(ind[1] <= p && ind[1] >= 0, "Wrong Multiindex.");
			UG_ASSERT(ind[2] <= p && ind[2] >= 0, "Wrong Multiindex.");
		}

	private:
		Polynomial1D m_vPolynom[p+1];
		Polynomial1D m_vDPolynom[p+1];
};

} //namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET__LAGRANGE__LAGRANGE__ */
