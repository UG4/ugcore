/*
 * lagrange.h
 *
 *  Created on: 17.11.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__LAGRANGE__LAGRANGE__
#define __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__LAGRANGE__LAGRANGE__

#include "../common/lagrange1d.h"
#include "../local_finite_element_provider.h"
#include "../local_dof_set.h"
#include "lagrange_local_dof.h"
#include "lib_disc/common/multi_index.h"
#include "common/util/metaprogramming_util.h"
#include "lib_grid/grid/grid_base_objects.h"

namespace ug{

/// Lagrange Shape Function Set without virtual functions and fixed order
template <typename TRefElem, int TOrder>
class LagrangeLSFS;

/// Lagrange Shape Function Set without virtual functions and flexible order
template <typename TRefElem>
class FlexLagrangeLSFS;


///////////////////////////////////////////////////////////////////////////////
// Vertex
///////////////////////////////////////////////////////////////////////////////

/// specialization for Vertex
/**
 * Lagrange shape function of any order for the Reference Vertex
 * \tparam 	TOrder		requested order
 */
template <int TOrder>
class LagrangeLSFS<ReferenceVertex, TOrder>
	: public LagrangeLDS<ReferenceVertex>,
	  public BaseLSFS<LagrangeLSFS<ReferenceVertex, TOrder>, 0>
{
	private:
	///	abbreviation for order
		static const size_t p = TOrder;

	///	base class
		typedef BaseLSFS<LagrangeLSFS<ReferenceVertex, TOrder>, 0> base_type;

	public:
	///	Shape type
		typedef typename base_type::shape_type shape_type;

	///	Gradient type
		typedef typename base_type::grad_type grad_type;

	public:
	///	Order of Shape functions
		static const size_t order = TOrder;

	///	Dimension, where shape functions are defined
		static const int dim = ReferenceVertex::dim;

	/// Number of shape functions
		static const size_t nsh = 1;

	public:
	///	Constructor
		LagrangeLSFS() {}

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline bool continuous() const {return true;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
		{
			return true;
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline shape_type shape(size_t i, const MathVector<dim>& x) const
		{
			return 1.0;
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		inline void grad(grad_type& g, size_t i, const MathVector<dim>& x) const
		{
		}
};

/// specialization for Edges
/**
 * Lagrange shape function of any order for the Reference Edge
 */
template <>
class FlexLagrangeLSFS<ReferenceVertex>
	: public LagrangeLDS<ReferenceVertex>,
	  public BaseLSFS<FlexLagrangeLSFS<ReferenceVertex>, 0>
{
	public:
	///	Dimension, where shape functions are defined
		static const int dim = ReferenceVertex::dim;

	public:
	///	default Constructor
		FlexLagrangeLSFS() {}

	///	Constructor
		FlexLagrangeLSFS(size_t order) {p = order;}

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline bool continuous() const {return true;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
		{
			return true;
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline shape_type shape(size_t i, const MathVector<dim>& x) const
		{
			return 1.0;
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		inline void grad(grad_type& g, size_t i, const MathVector<dim>& x) const
		{
		}

	protected:
	///	order
		size_t p;

	/// Number of shape functions
		static const size_t nsh = 1;
};

///////////////////////////////////////////////////////////////////////////////
// Edge
///////////////////////////////////////////////////////////////////////////////

/// specialization for Edges
/**
 * Lagrange shape function of any order for the Reference Edge
 * \tparam 	TOrder		requested order
 */
template <int TOrder>
class LagrangeLSFS<ReferenceEdge, TOrder>
	: public LagrangeLDS<ReferenceEdge>,
	  public BaseLSFS<LagrangeLSFS<ReferenceEdge, TOrder>, 1>
{
	private:
	///	abbreviation for order
		static const size_t p = TOrder;

	///	base class
		typedef BaseLSFS<LagrangeLSFS<ReferenceEdge, TOrder>, 1> base_type;

	public:
	///	Shape type
		typedef typename base_type::shape_type shape_type;

	///	Gradient type
		typedef typename base_type::grad_type grad_type;

	public:
	///	Order of Shape functions
		static const size_t order = TOrder;

	///	Dimension, where shape functions are defined
		static const int dim = ReferenceEdge::dim;

	/// Number of shape functions
		static const size_t nsh = p+1;

	public:
	///	Constructor
		LagrangeLSFS();

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline bool continuous() const {return true;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
		{
			pos = EquidistantLagrange1D::position(multi_index(i)[0], p);
			return true;
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline shape_type shape(size_t i, const MathVector<dim>& x) const
		{
			return m_vPolynom[multi_index(i)[0]].value(x[0]);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		inline void grad(grad_type& g, size_t i, const MathVector<dim>& x) const
		{
			g[0] = m_vDPolynom[multi_index(i)[0]].value(x[0]);
		}

	///	return Multi index for index i
		inline const MathVector<dim,int>& multi_index(size_t i) const
		{
			check_index(i);
			return m_vMultiIndex[i];
		}

	///	return the index for a multi_index
		inline size_t index(const MathVector<dim,int>& ind) const
		{
			check_multi_index(ind);
			for(size_t i=0; i<nsh; ++i)
				if(multi_index(i) == ind) return i;
			UG_THROW("Index not found in LagrangeLSFS");
		}

	///	return Multi index for index i
		inline MathVector<dim,int> mapped_multi_index(size_t i) const
		{
			check_index(i);
			return MathVector<1,int>(i);
		}

	///	return the index for a multi_index
		inline size_t mapped_index(const MathVector<dim,int>& ind) const
		{
			check_multi_index(ind);
			return ind[0];
		}

	///	checks in debug mode that index is valid
		inline void check_index(size_t i) const
		{
			UG_ASSERT(i < nsh, "Wrong index.");
		}

	///	checks in debug mode that multi-index is valid
		inline void check_multi_index(const MathVector<dim,int>& ind) const
		{
			UG_ASSERT(ind[0] < (int)nsh && ind[0] >= 0, "Wrong MultiIndex");
		}

	protected:
		Polynomial1D m_vPolynom[p+1];	///< Shape Polynomials
		Polynomial1D m_vDPolynom[p+1];	///< Derivative of Shape Polynomial

		MathVector<dim,int> m_vMultiIndex[nsh];
};

/// specialization for Edges
/**
 * Lagrange shape function of any order for the Reference Edge
 */
template <>
class FlexLagrangeLSFS<ReferenceEdge>
	: public LagrangeLDS<ReferenceEdge>,
	  public BaseLSFS<FlexLagrangeLSFS<ReferenceEdge>, 1>
{
	public:
	///	Dimension, where shape functions are defined
		static const int dim = ReferenceEdge::dim;

	public:
	///	default Constructor
		FlexLagrangeLSFS() {set_order(1);}

	///	Constructor
		FlexLagrangeLSFS(size_t order) {set_order(order);}

	///	sets the order
		void set_order(size_t order);

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline bool continuous() const {return true;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
		{
			pos = EquidistantLagrange1D::position(multi_index(i)[0], p);
			return true;
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline shape_type shape(size_t i, const MathVector<dim>& x) const
		{
			return m_vPolynom[multi_index(i)[0]].value(x[0]);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		inline void grad(grad_type& g, size_t i, const MathVector<dim>& x) const
		{
			g[0] = m_vDPolynom[multi_index(i)[0]].value(x[0]);
		}

	///	return Multi index for index i
		inline const MathVector<dim,int>& multi_index(size_t i) const
		{
			check_index(i);
			return m_vMultiIndex[i];
		}

	///	return the index for a multi_index
		inline size_t index(const MathVector<dim,int>& ind) const
		{
			check_multi_index(ind);
			for(size_t i=0; i<nsh; ++i)
				if(multi_index(i) == ind) return i;
			UG_THROW("Index not found in LagrangeLSFS");
		}

	///	return Multi index for index i
		inline MathVector<dim,int> mapped_multi_index(size_t i) const
		{
			check_index(i);
			return MathVector<1,int>(i);
		}

	///	return the index for a multi_index
		inline size_t mapped_index(const MathVector<dim,int>& ind) const
		{
			check_multi_index(ind);
			return ind[0];
		}

	///	checks in debug mode that index is valid
		inline void check_index(size_t i) const
		{
			UG_ASSERT(i < nsh, "Wrong index.");
		}

	///	checks in debug mode that multi-index is valid
		inline void check_multi_index(const MathVector<dim,int>& ind) const
		{
			UG_ASSERT(ind[0] < (int)nsh && ind[0] >= 0, "Wrong MultiIndex");
		}

	protected:
	///	order
		size_t p;

	/// Number of shape functions
		size_t nsh;

		std::vector<Polynomial1D> m_vPolynom;	///< Shape Polynomials
		std::vector<Polynomial1D> m_vDPolynom;	///< Derivative of Shape Polynomial

		std::vector<MathVector<dim,int> > m_vMultiIndex;
};

///////////////////////////////////////////////////////////////////////////////
// Triangle
///////////////////////////////////////////////////////////////////////////////

template <int TOrder>
class LagrangeLSFS<ReferenceTriangle, TOrder>
	: public LagrangeLDS<ReferenceTriangle>,
	  public BaseLSFS<LagrangeLSFS<ReferenceTriangle, TOrder>, 2>
{
	private:
	///	abbreviation for order
		static const size_t p = TOrder;

	///	base class
		typedef BaseLSFS<LagrangeLSFS<ReferenceTriangle, TOrder>, 2> base_type;

	public:
	///	Shape type
		typedef typename base_type::shape_type shape_type;

	///	Gradient type
		typedef typename base_type::grad_type grad_type;

	public:
	///	Order of Shape functions
		static const size_t order = TOrder;

	///	Dimension, where shape functions are defined
		static const int dim = ReferenceTriangle::dim;

	/// Number of shape functions
		static const size_t nsh = BinomialCoefficient<dim + p, p>::value;

	public:
	///	Constructor
		LagrangeLSFS();

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline bool continuous() const {return true;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
		{
		//	get Multi Index
			MathVector<dim,int> ind = multi_index(i);

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

	///	shape value for a Multi Index
		inline number shape(const MathVector<dim,int>& ind, const MathVector<dim>& x) const
		{
			check_multi_index(ind);
			//ReferenceTriangle::check_position(x);

		//	get adjoint barycentric index
			const size_t i0 = p - ind[0] - ind[1];
			const number x0 = 1.0 - x[0] - x[1];

			return    m_vPolynom[ ind[0] ].value(x[0])
					* m_vPolynom[ ind[1] ].value(x[1])
					* m_vPolynom[ i0     ].value(x0);
		}


	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline void grad(grad_type& g, const size_t i,	const MathVector<dim>& x) const
		{
			grad(g, multi_index(i), x);
		}

	///	evaluates the gradient
		inline void grad(grad_type& g, const MathVector<dim,int> ind,
		               	   	   	   	   	   const MathVector<dim>& x) const
		{
			check_multi_index(ind);
			//ReferenceTriangle::check_position(x);

		//	get adjoint barycentric index and position
			const size_t i0 = p - ind[0] - ind[1];
			const number x0 = 1.0 - x[0] - x[1];

			UG_ASSERT(i0 <= (int)p, "Wrong Multiindex.");
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

	///	return Multi index for index i
		inline const MathVector<dim,int>& multi_index(size_t i) const
		{
			check_index(i);
			return m_vMultiIndex[i];
		}

	///	return the index for a multi_index
		inline size_t index(const MathVector<dim,int>& ind) const
		{
			check_multi_index(ind);
			for(size_t i=0; i<nsh; ++i)
				if(multi_index(i) == ind) return i;
			UG_THROW("Index not found in LagrangeLSFS");
		}

	///	return the index for a multi_index
		inline size_t mapped_index(const MathVector<dim,int>& ind) const
		{
			check_multi_index(ind);

			size_t res = ind[0];
			for(int i = 0; i < ind[1]; ++i)
				res += (p+1-i);

			check_index(res);
			return res;
		}

	///	return the multi_index for an index
		inline MathVector<dim,int> mapped_multi_index(size_t i) const
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
			return MathVector<dim,int>( i0, i1 );
		}

	///	checks in debug mode that index is valid
		inline void check_index(size_t i) const
		{
			UG_ASSERT(i < nsh, "Wrong index.");
		}

	///	checks in debug mode that multi-index is valid
		inline void check_multi_index(const MathVector<dim,int>& ind) const
		{
			UG_ASSERT(ind[0] <= (int)p && ind[0]>=0, "Wrong Multiindex.");
			UG_ASSERT(ind[1] <= (int)p && ind[0]>=0, "Wrong Multiindex.");
			UG_ASSERT(ind[0] + ind[1] <= (int)p && ind[0]>=0, "Wrong Multiindex.");
		}

	private:
		Polynomial1D m_vPolynom[p+1];
		Polynomial1D m_vDPolynom[p+1];

		MathVector<dim,int> m_vMultiIndex[nsh];
};


template <>
class FlexLagrangeLSFS<ReferenceTriangle>
	: public LagrangeLDS<ReferenceTriangle>,
	  public BaseLSFS<FlexLagrangeLSFS<ReferenceTriangle>, 2>
{
	public:
	///	Dimension, where shape functions are defined
		static const int dim = ReferenceTriangle::dim;

	public:
	///	default Constructor
		FlexLagrangeLSFS() {set_order(1);}

	///	Constructor
		FlexLagrangeLSFS(size_t order) {set_order(order);}

	///	sets the order
		void set_order(size_t order);

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline bool continuous() const {return true;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
		{
		//	get Multi Index
			MathVector<dim,int> ind = multi_index(i);

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

	///	shape value for a Multi Index
		inline number shape(const MathVector<dim,int>& ind, const MathVector<dim>& x) const
		{
			check_multi_index(ind);
			//ReferenceTriangle::check_position(x);

		//	get adjoint barycentric index
			const size_t i0 = p - ind[0] - ind[1];
			const number x0 = 1.0 - x[0] - x[1];

			return    m_vPolynom[ ind[0] ].value(x[0])
					* m_vPolynom[ ind[1] ].value(x[1])
					* m_vPolynom[ i0     ].value(x0);
		}


	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline void grad(grad_type& g, const size_t i,	const MathVector<dim>& x) const
		{
			grad(g, multi_index(i), x);
		}

	///	evaluates the gradient
		inline void grad(grad_type& g, const MathVector<dim,int> ind,
		               	   	   	   	   	   const MathVector<dim>& x) const
		{
			check_multi_index(ind);
			//ReferenceTriangle::check_position(x);

		//	get adjoint barycentric index and position
			const int i0 = p - ind[0] - ind[1];
			const number x0 = 1.0 - x[0] - x[1];

			UG_ASSERT(i0 <= (int)p && i0 >= 0, "Wrong Multiindex.");
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

	///	return Multi index for index i
		inline const MathVector<dim,int>& multi_index(size_t i) const
		{
			check_index(i);
			return m_vMultiIndex[i];
		}

	///	return the index for a multi_index
		inline size_t index(const MathVector<dim,int>& ind) const
		{
			check_multi_index(ind);
			for(size_t i=0; i<nsh; ++i)
				if(multi_index(i) == ind) return i;
			UG_THROW("Index not found in LagrangeLSFS");
		}

	///	return the index for a multi_index
		inline size_t mapped_index(const MathVector<dim,int>& ind) const
		{
			check_multi_index(ind);

			size_t res = ind[0];
			for(int i = 0; i < ind[1]; ++i)
				res += (p+1-i);

			check_index(res);
			return res;
		}

	///	return the multi_index for an index
		inline MathVector<dim,int> mapped_multi_index(size_t i) const
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
			return MathVector<dim,int>( i0, i1 );
		}

	///	checks in debug mode that index is valid
		inline void check_index(size_t i) const
		{
			UG_ASSERT(i < nsh, "Wrong index.");
		}

	///	checks in debug mode that multi-index is valid
		inline void check_multi_index(const MathVector<dim,int>& ind) const
		{
			UG_ASSERT(ind[0] <= (int)p && ind[0]>=0, "Wrong Multiindex.");
			UG_ASSERT(ind[1] <= (int)p && ind[0]>=0, "Wrong Multiindex.");
			UG_ASSERT(ind[0] + ind[1] <= (int)p && ind[0]>=0, "Wrong Multiindex.");
		}

	private:
	///	order
		size_t p;

	/// Number of shape functions
		size_t nsh;

		std::vector<Polynomial1D> m_vPolynom;	///< Shape Polynomials
		std::vector<Polynomial1D> m_vDPolynom;	///< Derivative of Shape Polynomial

		std::vector<MathVector<dim,int> > m_vMultiIndex;
};

///////////////////////////////////////////////////////////////////////////////
// Quadrilateral
///////////////////////////////////////////////////////////////////////////////

template <int TOrder>
class LagrangeLSFS<ReferenceQuadrilateral, TOrder>
	: public LagrangeLDS<ReferenceQuadrilateral>,
	  public BaseLSFS<LagrangeLSFS<ReferenceQuadrilateral, TOrder>, 2>
{
	private:
	///	abbreviation for order
		static const size_t p = TOrder;

	///	base class
		typedef BaseLSFS<LagrangeLSFS<ReferenceQuadrilateral, TOrder>, 2> base_type;

	public:
	///	Shape type
		typedef typename base_type::shape_type shape_type;

	///	Gradient type
		typedef typename base_type::grad_type grad_type;

	public:
	///	Order of Shape functions
		static const size_t order = TOrder;

	///	Dimension, where shape functions are defined
		static const int dim = ReferenceQuadrilateral::dim;

	/// Number of shape functions
		static const size_t nsh = (p+1)*(p+1);

	public:
	///	Constructor
		LagrangeLSFS();

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline bool continuous() const {return true;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
		{
		//	get Multi Index
			MathVector<dim,int> ind = multi_index(i);

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

	///	shape value for a Multi Index
		inline number shape(const MathVector<dim,int>& ind, const MathVector<dim>& x) const
		{
			check_multi_index(ind);
			//ReferenceQuadrilateral::check_position(x);

			return    m_vPolynom[ ind[0] ].value(x[0])
					* m_vPolynom[ ind[1] ].value(x[1]);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		inline void grad(grad_type& g, const size_t i, const MathVector<dim>& x) const
		{
			grad(g, multi_index(i), x);
		}

	/// evaluates the gradient
		inline void grad(grad_type& g, const MathVector<dim,int> ind,
		               	   	   	const MathVector<dim>& x) const
		{
			check_multi_index(ind);
			//ReferenceQuadrilateral::check_position(x);

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

	///	return Multi index for index i
		inline const MathVector<dim,int>& multi_index(size_t i) const
		{
			check_index(i);
			return m_vMultiIndex[i];
		}

	///	return the index for a multi_index
		inline size_t index(const MathVector<dim,int>& ind) const
		{
			check_multi_index(ind);
			for(size_t i=0; i<nsh; ++i)
				if(multi_index(i) == ind) return i;
			UG_THROW("Index not found in LagrangeLSFS");
		}

	///	return the index for a multi_index
		inline size_t mapped_index(const MathVector<dim,int>& ind) const
		{
			check_multi_index(ind);

			return ind[1] * (p+1) + ind[0];
		}

	///	return the multi_index for an index
		inline MathVector<dim,int> mapped_multi_index(size_t i) const
		{
			check_index(i);

			return MathVector<dim,int>( i%(p+1), i/(p+1) );
		}

	///	checks in debug mode that index is valid
		inline void check_index(size_t i) const
		{
			UG_ASSERT(i < nsh, "Wrong index.");
		}

	///	checks in debug mode that multi-index is valid
		inline void check_multi_index(const MathVector<dim,int>& ind) const
		{
			UG_ASSERT(ind[0] <= (int)p && ind[0] >= 0, "Wrong Multiindex.");
			UG_ASSERT(ind[1] <= (int)p && ind[1] >= 0, "Wrong Multiindex.");
		}

	private:
		Polynomial1D m_vPolynom[p+1];
		Polynomial1D m_vDPolynom[p+1];

		MathVector<dim,int> m_vMultiIndex[nsh];
};


template <>
class FlexLagrangeLSFS<ReferenceQuadrilateral>
	: public LagrangeLDS<ReferenceQuadrilateral>,
	  public BaseLSFS<FlexLagrangeLSFS<ReferenceQuadrilateral>, 2>
{
	public:
	///	Dimension, where shape functions are defined
		static const int dim = ReferenceQuadrilateral::dim;

	public:
	///	default Constructor
		FlexLagrangeLSFS() {set_order(1);}

	///	Constructor
		FlexLagrangeLSFS(size_t order) {set_order(order);}

	///	sets the order
		void set_order(size_t order);

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline bool continuous() const {return true;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
		{
		//	get Multi Index
			MathVector<dim,int> ind = multi_index(i);

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

	///	shape value for a Multi Index
		inline number shape(const MathVector<dim,int>& ind, const MathVector<dim>& x) const
		{
			check_multi_index(ind);
			//ReferenceQuadrilateral::check_position(x);

			return    m_vPolynom[ ind[0] ].value(x[0])
					* m_vPolynom[ ind[1] ].value(x[1]);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		inline void grad(grad_type& g, const size_t i, const MathVector<dim>& x) const
		{
			grad(g, multi_index(i), x);
		}

	/// evaluates the gradient
		inline void grad(grad_type& g, const MathVector<dim,int> ind,
		               	   	   	const MathVector<dim>& x) const
		{
			check_multi_index(ind);
			//ReferenceQuadrilateral::check_position(x);

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

	///	return Multi index for index i
		inline const MathVector<dim,int>& multi_index(size_t i) const
		{
			check_index(i);
			return m_vMultiIndex[i];
		}

	///	return the index for a multi_index
		inline size_t index(const MathVector<dim,int>& ind) const
		{
			check_multi_index(ind);
			for(size_t i=0; i<nsh; ++i)
				if(multi_index(i) == ind) return i;
			UG_THROW("Index not found in LagrangeLSFS");
		}

	///	return the index for a multi_index
		inline size_t mapped_index(const MathVector<dim,int>& ind) const
		{
			check_multi_index(ind);

			return ind[1] * (p+1) + ind[0];
		}

	///	return the multi_index for an index
		inline MathVector<dim,int> mapped_multi_index(size_t i) const
		{
			check_index(i);

			return MathVector<dim,int>( i%(p+1), i/(p+1) );
		}

	///	checks in debug mode that index is valid
		inline void check_index(size_t i) const
		{
			UG_ASSERT(i < nsh, "Wrong index.");
		}

	///	checks in debug mode that multi-index is valid
		inline void check_multi_index(const MathVector<dim,int>& ind) const
		{
			UG_ASSERT(ind[0] <= (int)p && ind[0] >= 0, "Wrong Multiindex.");
			UG_ASSERT(ind[1] <= (int)p && ind[1] >= 0, "Wrong Multiindex.");
		}

	private:
	///	order
		size_t p;

	/// Number of shape functions
		size_t nsh;

		std::vector<Polynomial1D> m_vPolynom;	///< Shape Polynomials
		std::vector<Polynomial1D> m_vDPolynom;	///< Derivative of Shape Polynomial

		std::vector<MathVector<dim,int> > m_vMultiIndex;
};

///////////////////////////////////////////////////////////////////////////////
// Tetrahedron
///////////////////////////////////////////////////////////////////////////////

template <int TOrder>
class LagrangeLSFS<ReferenceTetrahedron, TOrder>
	: public LagrangeLDS<ReferenceTetrahedron>,
	  public BaseLSFS<LagrangeLSFS<ReferenceTetrahedron, TOrder>, 3>
{
	private:
	///	abbreviation for order
		static const size_t p = TOrder;

	///	base class
		typedef BaseLSFS<LagrangeLSFS<ReferenceTetrahedron, TOrder>, 3> base_type;

	public:
	///	Shape type
		typedef typename base_type::shape_type shape_type;

	///	Gradient type
		typedef typename base_type::grad_type grad_type;

	public:
	///	Order of Shape functions
		static const size_t order = TOrder;

	///	Dimension, where shape functions are defined
		static const int dim = ReferenceTetrahedron::dim;

	/// Number of shape functions
		static const size_t nsh = BinomialCoefficient<dim + p, p>::value;

	public:
	///	Constructor
		LagrangeLSFS();

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline bool continuous() const {return true;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
		{
		//	get Multi Index
			MathVector<dim,int> ind = multi_index(i);

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

	///	shape value for a Multi Index
		inline number shape(const MathVector<dim,int>& ind, const MathVector<dim>& x) const
		{
			check_multi_index(ind);
			//ReferenceTetrahedron::check_position(x);

		//	get adjoint barycentric index
			const size_t i0 = p - ind[0] - ind[1] - ind[2];
			const number x0 = 1.0 - x[0] - x[1] - x[2];

			return    m_vPolynom[ ind[0] ].value(x[0])
					* m_vPolynom[ ind[1] ].value(x[1])
					* m_vPolynom[ ind[2] ].value(x[2])
					* m_vPolynom[ i0 ].value(x0);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		inline void grad(grad_type& g, const size_t i, const MathVector<dim>& x) const
		{
			grad(g, multi_index(i), x);
		}

	///	evaluates the gradient
		inline void grad(grad_type& g, const MathVector<dim,int> ind,
		               	   	   	   	    const MathVector<dim>& x) const
		{
			check_multi_index(ind);
			//ReferenceTetrahedron::check_position(x);

		//	get adjoint barycentric index and position
			const size_t i0 = p - ind[0] - ind[1] - ind[2];
			const number x0 = 1 - x[0] - x[1] - x[2];

			UG_ASSERT(i0 <= p, "Wrong Multiindex.");
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

	///	return Multi index for index i
		inline const MathVector<dim,int>& multi_index(size_t i) const
		{
			check_index(i);
			return m_vMultiIndex[i];
		}

	///	return the index for a multi_index
		inline size_t index(const MathVector<dim,int>& ind) const
		{
			check_multi_index(ind);
			for(size_t i=0; i<nsh; ++i)
				if(multi_index(i) == ind) return i;
			UG_THROW("Index not found in LagrangeLSFS");
		}

	///	return the index for a multi_index
		inline size_t mapped_index(const MathVector<dim,int>& ind) const
		{
			check_multi_index(ind);

		//	add x range
			size_t res = ind[0];

		//	add y range
			for(int i = 0; i < ind[1]; ++i)
				res += (p+1-ind[2]-i);

		//	add z range
			for(int i2 = 0; i2 < ind[2]; ++i2)
				res += BinomCoeff(p+2-i2, p-i2);

			check_index(res);
			return res;
		}

	///	return the multi_index for an index
		inline MathVector<dim,int> mapped_multi_index(size_t i) const
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
							return MathVector<dim,int>( i0, i1, i2);

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
			return MathVector<dim,int>( i0, i1, i2);
		}

	///	checks in debug mode that index is valid
		inline void check_index(size_t i) const
		{
			UG_ASSERT(i < nsh, "Wrong index.");
		}

	///	checks in debug mode that multi-index is valid
		inline void check_multi_index(const MathVector<dim,int>& ind) const
		{
			UG_ASSERT(ind[0] <= (int)p && ind[0] >= 0, "Wrong Multiindex.");
			UG_ASSERT(ind[1] <= (int)p && ind[1] >= 0, "Wrong Multiindex.");
			UG_ASSERT(ind[2] <= (int)p && ind[2] >= 0, "Wrong Multiindex.");
			UG_ASSERT(ind[0] + ind[1] + ind[2] <= (int)p, "Wrong Multiindex.");
		}

	private:
		Polynomial1D m_vPolynom[p+1];
		Polynomial1D m_vDPolynom[p+1];

		MathVector<dim,int> m_vMultiIndex[nsh];
};


template <>
class FlexLagrangeLSFS<ReferenceTetrahedron>
	: public LagrangeLDS<ReferenceTetrahedron>,
	  public BaseLSFS<FlexLagrangeLSFS<ReferenceTetrahedron>, 3>
{
	public:
	///	Dimension, where shape functions are defined
		static const int dim = ReferenceTetrahedron::dim;

	public:
	///	default Constructor
		FlexLagrangeLSFS() {set_order(1);}

	///	Constructor
		FlexLagrangeLSFS(size_t order) {set_order(order);}

	///	sets the order
		void set_order(size_t order);

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline bool continuous() const {return true;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
		{
		//	get Multi Index
			MathVector<dim,int> ind = multi_index(i);

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

	///	shape value for a Multi Index
		inline number shape(const MathVector<dim,int>& ind, const MathVector<dim>& x) const
		{
			check_multi_index(ind);
			//ReferenceTetrahedron::check_position(x);

		//	get adjoint barycentric index
			const size_t i0 = p - ind[0] - ind[1] - ind[2];
			const number x0 = 1.0 - x[0] - x[1] - x[2];

			return    m_vPolynom[ ind[0] ].value(x[0])
					* m_vPolynom[ ind[1] ].value(x[1])
					* m_vPolynom[ ind[2] ].value(x[2])
					* m_vPolynom[ i0 ].value(x0);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		inline void grad(grad_type& g, const size_t i, const MathVector<dim>& x) const
		{
			grad(g, multi_index(i), x);
		}

	///	evaluates the gradient
		inline void grad(grad_type& g, const MathVector<dim,int> ind,
		               	   	   	   	    const MathVector<dim>& x) const
		{
			check_multi_index(ind);
			//ReferenceTetrahedron::check_position(x);

		//	get adjoint barycentric index and position
			const size_t i0 = p - ind[0] - ind[1] - ind[2];
			const number x0 = 1 - x[0] - x[1] - x[2];

			UG_ASSERT(i0 <= p, "Wrong Multiindex.");
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

	///	return Multi index for index i
		inline const MathVector<dim,int>& multi_index(size_t i) const
		{
			check_index(i);
			return m_vMultiIndex[i];
		}

	///	return the index for a multi_index
		inline size_t index(const MathVector<dim,int>& ind) const
		{
			check_multi_index(ind);
			for(size_t i=0; i<nsh; ++i)
				if(multi_index(i) == ind) return i;
			UG_THROW("Index not found in LagrangeLSFS");
		}

	///	return the index for a multi_index
		inline size_t mapped_index(const MathVector<dim,int>& ind) const
		{
			check_multi_index(ind);

		//	add x range
			size_t res = ind[0];

		//	add y range
			for(int i = 0; i < ind[1]; ++i)
				res += (p+1-ind[2]-i);

		//	add z range
			for(int i2 = 0; i2 < ind[2]; ++i2)
				res += BinomCoeff(p+2-i2, p-i2);

			check_index(res);
			return res;
		}

	///	return the multi_index for an index
		inline MathVector<dim,int> mapped_multi_index(size_t i) const
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
							return MathVector<dim,int>( i0, i1, i2);

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
			return MathVector<dim,int>( i0, i1, i2);
		}

	///	checks in debug mode that index is valid
		inline void check_index(size_t i) const
		{
			UG_ASSERT(i < nsh, "Wrong index.");
		}

	///	checks in debug mode that multi-index is valid
		inline void check_multi_index(const MathVector<dim,int>& ind) const
		{
			UG_ASSERT(ind[0] <= (int)p && ind[0] >= 0, "Wrong Multiindex.");
			UG_ASSERT(ind[1] <= (int)p && ind[1] >= 0, "Wrong Multiindex.");
			UG_ASSERT(ind[2] <= (int)p && ind[2] >= 0, "Wrong Multiindex.");
			UG_ASSERT(ind[0] + ind[1] + ind[2] <= (int)p, "Wrong Multiindex.");
		}

	private:
	///	order
		size_t p;

	/// Number of shape functions
		size_t nsh;

		std::vector<Polynomial1D> m_vPolynom;	///< Shape Polynomials
		std::vector<Polynomial1D> m_vDPolynom;	///< Derivative of Shape Polynomial

		std::vector<MathVector<dim,int> > m_vMultiIndex;
};

///////////////////////////////////////////////////////////////////////////////
// Prism
///////////////////////////////////////////////////////////////////////////////

template <int TOrder>
class LagrangeLSFS<ReferencePrism, TOrder>
	: public LagrangeLDS<ReferencePrism>,
	  public BaseLSFS<LagrangeLSFS<ReferencePrism, TOrder>, 3>
{
	private:
	///	abbreviation for order
		static const size_t p = TOrder;

	/// dofs per layer
		static const size_t dofPerLayer = BinomialCoefficient<2 + p, p>::value;

	///	base class
		typedef BaseLSFS<LagrangeLSFS<ReferencePrism, TOrder>, 3> base_type;

	public:
	///	Shape type
		typedef typename base_type::shape_type shape_type;

	///	Gradient type
		typedef typename base_type::grad_type grad_type;

	public:
	///	Order of Shape functions
		static const size_t order = TOrder;

	///	Dimension, where shape functions are defined
		static const int dim = ReferencePrism::dim;

	/// Number of shape functions
		static const size_t nsh = dofPerLayer*(p+1);

	public:
	///	Constructor
		LagrangeLSFS();

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline bool continuous() const {return true;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
		{
		//	get Multi Index
			MathVector<dim,int> ind = multi_index(i);

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

	///	shape value for a Multi Index
		inline number shape(const MathVector<dim,int>& ind, const MathVector<dim>& x) const
		{
			check_multi_index(ind);
			//ReferencePrism::check_position(x);

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
		void grad(grad_type& g, const size_t i, const MathVector<dim>& x) const
		{
			grad(g, multi_index(i), x);
		}

	///	evaluates the gradient
		void grad(grad_type& g, const MathVector<dim,int> ind,
		               	   	   	   	   	const MathVector<dim>& x) const
		{
			check_multi_index(ind);
			//ReferencePrism::check_position(x);

		//	get adjoint barycentric index and position
			const size_t i0 = p - ind[0] - ind[1];
			const number x0 = 1 - x[0] - x[1];

			UG_ASSERT(i0 <= p, "Wrong Multiindex.");
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

	///	return Multi index for index i
		inline const MathVector<dim,int>& multi_index(size_t i) const
		{
			check_index(i);
			return m_vMultiIndex[i];
		}

	///	return the index for a multi_index
		inline size_t index(const MathVector<dim,int>& ind) const
		{
			check_multi_index(ind);
			for(size_t i=0; i<nsh; ++i)
				if(multi_index(i) == ind) return i;
			UG_THROW("Index not found in LagrangeLSFS");
		}

	///	return the index for a multi_index
		inline size_t mapped_index(const MathVector<dim,int>& ind) const
		{
			check_multi_index(ind);

			size_t res = ind[0];
			for(int i = 0; i < ind[1]; ++i)
				res += (p+1-i);

			// add height
			res += ind[2] * dofPerLayer;

			return res;
		}

	///	return the multi_index for an index
		inline MathVector<dim,int> mapped_multi_index(size_t i) const
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

			return MathVector<dim,int>( i0, i1, i2);
		}

	///	checks in debug mode that index is valid
		inline void check_index(size_t i) const
		{
			UG_ASSERT(i < nsh, "Wrong index.");
		}

	///	checks in debug mode that multi-index is valid
		inline void check_multi_index(const MathVector<dim,int>& ind) const
		{
			UG_ASSERT(ind[0] <= (int)p && ind[0] >= 0, "Wrong Multiindex.");
			UG_ASSERT(ind[1] <= (int)p && ind[1] >= 0, "Wrong Multiindex.");
			UG_ASSERT(ind[0] + ind[1] <= (int)p, "Wrong Multiindex.");
			UG_ASSERT(ind[2] <= (int)p && ind[2] >= 0, "Wrong Multiindex.");
		}

	private:
		Polynomial1D m_vPolynom[p+1];
		Polynomial1D m_vDPolynom[p+1];
		Polynomial1D m_vTruncPolynom[p+1];
		Polynomial1D m_vDTruncPolynom[p+1];

		MathVector<dim,int> m_vMultiIndex[nsh];
};


template <>
class FlexLagrangeLSFS<ReferencePrism>
	: public LagrangeLDS<ReferencePrism>,
	  public BaseLSFS<FlexLagrangeLSFS<ReferencePrism>, 3>
{
	public:
	///	Dimension, where shape functions are defined
		static const int dim = ReferencePrism::dim;

	public:
	///	default Constructor
		FlexLagrangeLSFS() {set_order(1);}

	///	Constructor
		FlexLagrangeLSFS(size_t order) {set_order(order);}

	///	sets the order
		void set_order(size_t order);

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline bool continuous() const {return true;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
		{
		//	get Multi Index
			MathVector<dim,int> ind = multi_index(i);

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

	///	shape value for a Multi Index
		inline number shape(const MathVector<dim,int>& ind, const MathVector<dim>& x) const
		{
			check_multi_index(ind);
			//ReferencePrism::check_position(x);

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
		void grad(grad_type& g, const size_t i, const MathVector<dim>& x) const
		{
			grad(g, multi_index(i), x);
		}

	///	evaluates the gradient
		void grad(grad_type& g, const MathVector<dim,int> ind,
		               	   	   	   	   	const MathVector<dim>& x) const
		{
			check_multi_index(ind);
			//ReferencePrism::check_position(x);

		//	get adjoint barycentric index and position
			const size_t i0 = p - ind[0] - ind[1];
			const number x0 = 1 - x[0] - x[1];

			UG_ASSERT(i0 <= p, "Wrong Multiindex.");
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

	///	return Multi index for index i
		inline const MathVector<dim,int>& multi_index(size_t i) const
		{
			check_index(i);
			return m_vMultiIndex[i];
		}

	///	return the index for a multi_index
		inline size_t index(const MathVector<dim,int>& ind) const
		{
			check_multi_index(ind);
			for(size_t i=0; i<nsh; ++i)
				if(multi_index(i) == ind) return i;
			UG_THROW("Index not found in LagrangeLSFS");
		}

	///	return the index for a multi_index
		inline size_t mapped_index(const MathVector<dim,int>& ind) const
		{
			check_multi_index(ind);

			size_t res = ind[0];
			for(int i = 0; i < ind[1]; ++i)
				res += (p+1-i);

			// add height
			res += ind[2] * dofPerLayer;

			return res;
		}

	///	return the multi_index for an index
		inline MathVector<dim,int> mapped_multi_index(size_t i) const
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

			return MathVector<dim,int>( i0, i1, i2);
		}

	///	checks in debug mode that index is valid
		inline void check_index(size_t i) const
		{
			UG_ASSERT(i < nsh, "Wrong index.");
		}

	///	checks in debug mode that multi-index is valid
		inline void check_multi_index(const MathVector<dim,int>& ind) const
		{
			UG_ASSERT(ind[0] <= (int)p && ind[0] >= 0, "Wrong Multiindex.");
			UG_ASSERT(ind[1] <= (int)p && ind[1] >= 0, "Wrong Multiindex.");
			UG_ASSERT(ind[0] + ind[1] <= (int)p, "Wrong Multiindex.");
			UG_ASSERT(ind[2] <= (int)p && ind[2] >= 0, "Wrong Multiindex.");
		}

	private:
	///	order
		size_t p;

	/// Number of shape functions
		size_t nsh;

	///	number of dofs per layer
		size_t dofPerLayer;

		std::vector<Polynomial1D> m_vPolynom;
		std::vector<Polynomial1D> m_vDPolynom;
		std::vector<Polynomial1D> m_vTruncPolynom;
		std::vector<Polynomial1D> m_vDTruncPolynom;

		std::vector<MathVector<dim,int> > m_vMultiIndex;
};

///////////////////////////////////////////////////////////////////////////////
// Pyramid
///////////////////////////////////////////////////////////////////////////////

namespace {
template <int p>
struct NumberOfDoFsOfPyramid{
	enum{value = NumberOfDoFsOfPyramid<p-1>::value + (p+1)*(p+1)};
};

// specialization to end recursion
template <> struct NumberOfDoFsOfPyramid<1>{enum {value = 5};};
template <> struct NumberOfDoFsOfPyramid<0>{enum {value = 0};};
template <> struct NumberOfDoFsOfPyramid<-1>{enum {value = 0};};
} // end empty namespace

// todo: Implement higher order (impossible ?)
//	NOTE:	Currently only 1st order is implemented. There is no shape function
//			set for pyramids, that is continuous and allows a continuous
//			derivative in the inner of the pyramid. This is basically, since
//			one regards the pyramid as two tetrahedrons, glued together.
template <int TOrder>
class LagrangeLSFS<ReferencePyramid, TOrder>
	: public LagrangeLDS<ReferencePyramid>,
	  public BaseLSFS<LagrangeLSFS<ReferencePyramid, TOrder>, 3>
{
	private:
	///	abbreviation for order
		static const size_t p = TOrder;

	///	base class
		typedef BaseLSFS<LagrangeLSFS<ReferencePyramid, TOrder>, 3> base_type;

	public:
	///	Shape type
		typedef typename base_type::shape_type shape_type;

	///	Gradient type
		typedef typename base_type::grad_type grad_type;

	public:
	///	Order of Shape functions
		static const size_t order = TOrder;

	///	Dimension, where shape functions are defined
		static const int dim = 3;	//reference_element_type::dim; (compile error on OSX 10.5)

	/// Number of shape functions
		static const size_t nsh = NumberOfDoFsOfPyramid<p>::value;

	public:
	///	Constructor
		LagrangeLSFS();

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline bool continuous() const {return true;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
		{
		//	get Multi Index
			MathVector<dim,int> ind = multi_index(i);

		//	set position
			for(int d = 0; d < dim; ++d)
				pos[d] = EquidistantLagrange1D::position(ind[d], p);

			return true;
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline number shape(const size_t i, const MathVector<dim>& x) const
		{
		//	only first order
			if(p != 1) UG_THROW("Only 1. order Lagrange Pyramid implemented.");

		//	smaller value of x and y
			number m = x[0];
			if (x[0] > x[1]) m = x[1];

			switch(i)
			{
			  case 0 : return((1.0-x[0])*(1.0-x[1]) + x[2]*(m-1.0));
			  case 1 : return(x[0]*(1.0-x[1])       - x[2]*m);
			  case 2 : return(x[0]*x[1]             + x[2]*m);
			  case 3 : return((1.0-x[0])*x[1]       - x[2]*m);
			  case 4 : return(x[2]);
			  default: UG_THROW("Wrong index "<< i<<" in Pyramid");
			}
		}

	///	shape value for a Multi Index
		inline number shape(const MathVector<dim,int>& ind, const MathVector<dim>& x) const
		{
			check_multi_index(ind);
			//ReferencePyramid::check_position(x);

		//	forward
			return shape(index(ind), x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		inline void grad(grad_type& g, const size_t i, const MathVector<dim>& x) const
		{
		//	only first order
			if(p != 1) UG_THROW("Only 1. order Lagrange Pyramid implemented.");

			int m = 0;
			if (x[0] > x[1]) m = 1;

			switch(i)
			{
			  case 0:
				g[0] = -(1.0-x[1]);
				g[1] = -(1.0-x[0]) + x[2];
				g[2] = -(1.0-x[1]);
				g[m] += x[2];
				break;
			  case 1:
				g[0] = (1.0-x[1]);
				g[1] = -x[0];
				g[2] = -x[1];
				g[m] -= x[2];
				break;
			  case 2:
				g[0] = x[1];
				g[1] = x[0];
				g[2] = x[1];
				g[m] += x[2];
				break;
			  case 3:
				g[0] = -x[1];
				g[1] = 1.0-x[0];
				g[2] = -x[1];
				g[m] -= x[2];
				break;
		      case 4:
				g[0] = 0.0;
				g[1] = 0.0;
				g[2] = 1.0;
				break;
		      default: UG_THROW("Wrong index "<< i<<" in Pyramid");
			}

		}

	///	evaluates the gradient
		inline void grad(grad_type& g, const MathVector<dim,int> ind,
		               	   	   	   	   	   const MathVector<dim>& x) const
		{
			grad(g, index(ind), x);
		}

	///	return Multi index for index i
		inline const MathVector<dim,int>& multi_index(size_t i) const
		{
			check_index(i);
			return m_vMultiIndex[i];
		}

	///	return the index for a multi_index
		inline size_t index(const MathVector<dim,int>& ind) const
		{
			check_multi_index(ind);
			for(size_t i=0; i<nsh; ++i)
				if(multi_index(i) == ind) return i;
			UG_THROW("Index not found in LagrangeLSFS");
		}

	///	return the index for a multi_index
		inline size_t mapped_index(const MathVector<dim,int>& ind) const
		{
			check_multi_index(ind);

			size_t res = 0;

		//	add layers that are completely filled
			for(int i2 = 0; i2 < ind[2]; ++i2)
				res += (p+1-i2)*(p+1-i2);

		//	add dofs of top z-layer
			res += ind[1] * (p+1-ind[2]) + ind[0];

			check_index(res);
			return res;
		}

	///	return the multi_index for an index
		inline MathVector<dim,int> mapped_multi_index(size_t i) const
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

			return MathVector<dim,int>( iTmp%(p+1-i2), iTmp/(p+1-i2), i2);
		}

	///	checks in debug mode that index is valid
		inline void check_index(size_t i) const
		{
			UG_ASSERT(i < nsh, "Wrong index.");
		}

	///	checks in debug mode that multi-index is valid
		inline void check_multi_index(const MathVector<dim,int>& ind) const
		{
			UG_ASSERT(ind[0] <= (int)p-ind[2] && ind[0] >= 0, "Wrong Multiindex.");
			UG_ASSERT(ind[1] <= (int)p-ind[2] && ind[1] >= 0, "Wrong Multiindex.");
			UG_ASSERT(ind[2] <= (int)p && ind[2] >= 0, "Wrong Multiindex.");
		}

	private:
		std::vector<std::vector<Polynomial1D> > m_vvPolynom;
		std::vector<std::vector<Polynomial1D> > m_vvDPolynom;

		MathVector<dim,int> m_vMultiIndex[nsh];
};

///////////////////////////////////////////////////////////////////////////////
// Hexahedron
///////////////////////////////////////////////////////////////////////////////

template <int TOrder>
class LagrangeLSFS<ReferenceHexahedron, TOrder>
	: public LagrangeLDS<ReferenceHexahedron>,
	  public BaseLSFS<LagrangeLSFS<ReferenceHexahedron, TOrder>, 3>
{
	private:
	///	abbreviation for order
		static const size_t p = TOrder;

	///	base class
		typedef BaseLSFS<LagrangeLSFS<ReferenceHexahedron, TOrder>, 3> base_type;

	public:
	///	Shape type
		typedef typename base_type::shape_type shape_type;

	///	Gradient type
		typedef typename base_type::grad_type grad_type;

	public:
	///	Order of Shape functions
		static const size_t order = TOrder;

	///	Dimension, where shape functions are defined
		static const int dim = ReferenceHexahedron::dim;

	/// Number of shape functions
		static const size_t nsh = (p+1)*(p+1)*(p+1);

	public:
	///	Constructor
		LagrangeLSFS();

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline bool continuous() const {return true;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
		{
		//	get Multi Index
			MathVector<dim,int> ind = multi_index(i);

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

	///	shape value for a Multi Index
		inline number shape(const MathVector<dim,int>& ind, const MathVector<dim>& x) const
		{
			check_multi_index(ind);
			//ReferenceHexahedron::check_position(x);

			return 	  m_vPolynom[ ind[0] ].value(x[0])
					* m_vPolynom[ ind[1] ].value(x[1])
					* m_vPolynom[ ind[2] ].value(x[2]);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		inline void grad(grad_type& g, const size_t i, const MathVector<dim>& x) const
		{
			grad(g, multi_index(i), x);
		}

	///	evaluates the gradient
		inline void grad(grad_type& g, const MathVector<dim,int> ind,
		          	  	  	  	const MathVector<dim>& x) const
		{
			check_multi_index(ind);
			//ReferenceHexahedron::check_position(x);

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

	///	return Multi index for index i
		inline const MathVector<dim,int>& multi_index(size_t i) const
		{
			check_index(i);
			return m_vMultiIndex[i];
		}

	///	return the index for a multi_index
		inline size_t index(const MathVector<dim,int>& ind) const
		{
			check_multi_index(ind);
			for(size_t i=0; i<nsh; ++i)
				if(multi_index(i) == ind) return i;
			UG_THROW("Index not found in LagrangeLSFS");
		}

	///	return the index for a multi_index
		inline size_t mapped_index(const MathVector<dim,int>& ind) const
		{
			check_multi_index(ind);

			return ind[2] * (p+1)*(p+1) + ind[1] * (p+1) + ind[0];
		}

	///	return the multi_index for an index
		inline MathVector<dim,int> mapped_multi_index(size_t i) const
		{
			check_index(i);

			return MathVector<dim,int>( i%(p+1), i/(p+1)%(p+1), i/((p+1)*(p+1)));
		}

	///	checks in debug mode that index is valid
		inline void check_index(size_t i) const
		{
			UG_ASSERT(i < nsh, "Wrong index.");
		}

	///	checks in debug mode that multi-index is valid
		inline void check_multi_index(const MathVector<dim,int>& ind) const
		{
			UG_ASSERT(ind[0] <= (int)p && ind[0] >= 0, "Wrong Multiindex.");
			UG_ASSERT(ind[1] <= (int)p && ind[1] >= 0, "Wrong Multiindex.");
			UG_ASSERT(ind[2] <= (int)p && ind[2] >= 0, "Wrong Multiindex.");
		}

	private:
		Polynomial1D m_vPolynom[p+1];
		Polynomial1D m_vDPolynom[p+1];

		MathVector<dim,int> m_vMultiIndex[nsh];
};

template <>
class FlexLagrangeLSFS<ReferenceHexahedron>
	: public LagrangeLDS<ReferenceHexahedron>,
	  public BaseLSFS<FlexLagrangeLSFS<ReferenceHexahedron>, 3>
{
	public:
	///	Dimension, where shape functions are defined
		static const int dim = ReferenceHexahedron::dim;

	public:
	///	default Constructor
		FlexLagrangeLSFS() {set_order(1);}

	///	Constructor
		FlexLagrangeLSFS(size_t order) {set_order(order);}

	///	sets the order
		void set_order(size_t order);

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline bool continuous() const {return true;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
		{
		//	get Multi Index
			MathVector<dim,int> ind = multi_index(i);

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

	///	shape value for a Multi Index
		inline number shape(const MathVector<dim,int>& ind, const MathVector<dim>& x) const
		{
			check_multi_index(ind);
			//ReferenceHexahedron::check_position(x);

			return 	  m_vPolynom[ ind[0] ].value(x[0])
					* m_vPolynom[ ind[1] ].value(x[1])
					* m_vPolynom[ ind[2] ].value(x[2]);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		inline void grad(grad_type& g, const size_t i, const MathVector<dim>& x) const
		{
			grad(g, multi_index(i), x);
		}

	///	evaluates the gradient
		inline void grad(grad_type& g, const MathVector<dim,int> ind,
		          	  	  	  	const MathVector<dim>& x) const
		{
			check_multi_index(ind);
			//ReferenceHexahedron::check_position(x);

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

	///	return Multi index for index i
		inline const MathVector<dim,int>& multi_index(size_t i) const
		{
			check_index(i);
			return m_vMultiIndex[i];
		}

	///	return the index for a multi_index
		inline size_t index(const MathVector<dim,int>& ind) const
		{
			check_multi_index(ind);
			for(size_t i=0; i<nsh; ++i)
				if(multi_index(i) == ind) return i;
			UG_THROW("Index not found in LagrangeLSFS");
		}

	///	return the index for a multi_index
		inline size_t mapped_index(const MathVector<dim,int>& ind) const
		{
			check_multi_index(ind);

			return ind[2] * (p+1)*(p+1) + ind[1] * (p+1) + ind[0];
		}

	///	return the multi_index for an index
		inline MathVector<dim,int> mapped_multi_index(size_t i) const
		{
			check_index(i);

			return MathVector<dim,int>( i%(p+1), i/(p+1)%(p+1), i/((p+1)*(p+1)));
		}

	///	checks in debug mode that index is valid
		inline void check_index(size_t i) const
		{
			UG_ASSERT(i < nsh, "Wrong index.");
		}

	///	checks in debug mode that multi-index is valid
		inline void check_multi_index(const MathVector<dim,int>& ind) const
		{
			UG_ASSERT(ind[0] <= (int)p && ind[0] >= 0, "Wrong Multiindex.");
			UG_ASSERT(ind[1] <= (int)p && ind[1] >= 0, "Wrong Multiindex.");
			UG_ASSERT(ind[2] <= (int)p && ind[2] >= 0, "Wrong Multiindex.");
		}

	private:
	///	order
		size_t p;

	/// Number of shape functions
		size_t nsh;

		std::vector<Polynomial1D> m_vPolynom;	///< Shape Polynomials
		std::vector<Polynomial1D> m_vDPolynom;	///< Derivative of Shape Polynomial

		std::vector<MathVector<dim,int> > m_vMultiIndex;
};

///////////////////////////////////////////////////////////////////////////////
// Octahedron
///////////////////////////////////////////////////////////////////////////////

//	NOTE:	Currently only 1st order is implemented. There is no shape function
//			set for octahedrons, that is continuous and allows a continuous
//			derivative in the inner of the pyramid. This is basically, since
//			one regards the octahedron as 4 tetrahedrons, glued together.
template <int TOrder>
class LagrangeLSFS<ReferenceOctahedron, TOrder>
	: public LagrangeLDS<ReferenceOctahedron>,
	  public BaseLSFS<LagrangeLSFS<ReferenceOctahedron, TOrder>, 3>
{
	private:
	///	abbreviation for order
		static const size_t p = TOrder;

	///	base class
		typedef BaseLSFS<LagrangeLSFS<ReferenceOctahedron, TOrder>, 3> base_type;

	public:
	///	Shape type
		typedef typename base_type::shape_type shape_type;

	///	Gradient type
		typedef typename base_type::grad_type grad_type;

	public:
	///	Order of Shape functions
		static const size_t order = TOrder;

	///	Dimension, where shape functions are defined
		static const int dim = 3;	//reference_element_type::dim; (compile error on OSX 10.5)

	/// Number of shape functions
		static const size_t nsh = 6;

	public:
	///	Constructor
		LagrangeLSFS();

	///	\copydoc ug::LocalShapeFunctionSet::continuous()
		inline bool continuous() const {return true;}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
		{
		//	get Multi Index
			MathVector<dim,int> ind = multi_index(i);

		//	set position
			for(int d = 0; d < dim; ++d)
				pos[d] = EquidistantLagrange1D::position(ind[d], p);

			return true;
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline number shape(const size_t i, const MathVector<dim>& x) const
		{
		//	only first order
			if(p != 1) UG_THROW("Only 1. order Lagrange Octahedron implemented.");

		//	shape analogously to pyramidal case introducing additional distinction of cases
		//	z >= 0 and z < 0
			const number z_sgn 	= (x[2] < 0) ? -1.0 : 1.0;

			switch(i)
			{
			  case 0 :
				if (x[2] < 0)
				  return(-1.0*x[2]);
				else
				  return(0.0);
			  case 1 :
				if (x[0] > x[1])
				  return((1.0-x[0])*(1.0-x[1]) + z_sgn*x[2]*(x[1]-1.0));
				else
				  return((1.0-x[0])*(1.0-x[1]) + z_sgn*x[2]*(x[0]-1.0));
			  case 2 :
				if (x[0] > x[1])
				  return(x[0]*(1.0-x[1])       - z_sgn*x[2]*x[1]);
				else
				  return(x[0]*(1.0-x[1])       - z_sgn*x[2]*x[0]);
			  case 3 :
				if (x[0] > x[1])
				  return(x[0]*x[1]             + z_sgn*x[2]*x[1]);
				else
				  return(x[0]*x[1]             + z_sgn*x[2]*x[0]);
			  case 4 :
				if (x[0] > x[1])
				  return((1.0-x[0])*x[1]       - z_sgn*x[2]*x[1]);
				else
				  return((1.0-x[0])*x[1]       - z_sgn*x[2]*x[0]);
			  case 5 :
				if (x[2] < 0)
				  return(0.0);
				else
				  return(x[2]);
			  default: UG_THROW("Wrong index "<< i<<" in Octahedron.");
			}
		}

	///	shape value for a Multi Index
		inline number shape(const MathVector<dim,int>& ind, const MathVector<dim>& x) const
		{
			check_multi_index(ind);

		//	forward
			return shape(index(ind), x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		inline void grad(grad_type& g, const size_t i, const MathVector<dim>& x) const
		{
		//	only first order
			if(p != 1) UG_THROW("Only 1. order Lagrange Octahedron implemented.");

			//	shape analogously to pyramidal case introducing additional distinction of cases
			//	z >= 0 and z < 0
				const number z_sgn 	= (x[2] < 0) ? -1.0 : 1.0;

				switch(i)
				{
				  case 0:
					if (x[2] < 0.0)
					{
						g[0] = 0.0;
						g[1] = 0.0;
						g[2] = -1.0;
						break;
					}
					else
					{
						g[0] = 0.0;
						g[1] = 0.0;
						g[2] = 0.0;
						break;
					}
				  case 1:
					if (x[0] > x[1])
					  {
						g[0] = -(1.0-x[1]);
						g[1] = -(1.0-x[0]) + z_sgn*x[2];
						g[2] = -z_sgn*(1.0-x[1]);
						break;
					  }
					else
					  {
						g[0] = -(1.0-x[1]) + z_sgn*x[2];
						g[1] = -(1.0-x[0]);
						g[2] = -z_sgn*(1.0-x[0]);
						break;
					  }
				  case 2:
					if (x[0] > x[1])
					  {
						g[0] = (1.0-x[1]);
						g[1] = -x[0] - z_sgn*x[2];
						g[2] = -z_sgn*x[1];
						break;
					  }
					else
					  {
						g[0] = (1.0-x[1]) - z_sgn*x[2];
						g[1] = -x[0];
						g[2] = -z_sgn*x[0];
						break;
					  }
				  case 3:
					if (x[0] > x[1])
					  {
						g[0] = x[1];
						g[1] = x[0] + z_sgn*x[2];
						g[2] = z_sgn*x[1];
						break;
					  }
					else
					  {
						g[0] = x[1] + z_sgn*x[2];
						g[1] = x[0];
						g[2] = z_sgn*x[0];
						break;
					  }
				  case 4:
					if (x[0] > x[1])
					  {
						g[0] = -x[1];
						g[1] = 1.0-x[0] - z_sgn*x[2];
						g[2] = -z_sgn*x[1];
						break;
					  }
					else
					  {
						g[0] = -x[1] - z_sgn*x[2];
						g[1] = 1.0-x[0];
						g[2] = -z_sgn*x[0];
						break;
					  }
			      case 5:
			        if (x[2] < 0.0)
			        {
			        	g[0] = 0.0;
			        	g[1] = 0.0;
			        	g[2] = 0.0;
			        	break;
			        }
			        else
			        {
			        	g[0] = 0.0;
			        	g[1] = 0.0;
			        	g[2] = 1.0;
			        	break;
			        }
			      default: UG_THROW("Wrong index "<< i<<" in Octahedron.");
				}
		}

	///	evaluates the gradient
		inline void grad(grad_type& g, const MathVector<dim,int> ind,
		               	   	   	   	   	   const MathVector<dim>& x) const
		{
			grad(g, index(ind), x);
		}

	///	return Multi index for index i
		inline const MathVector<dim,int>& multi_index(size_t i) const
		{
			check_index(i);
			return m_vMultiIndex[i];
		}

	///	return the index for a multi_index
		inline size_t index(const MathVector<dim,int>& ind) const
		{
			check_multi_index(ind);
			for(size_t i=0; i<nsh; ++i)
				if(multi_index(i) == ind) return i;
			UG_THROW("Index not found in LagrangeLSFS");
		}

	///	checks in debug mode that index is valid
		inline void check_index(size_t i) const
		{
			UG_ASSERT(i < nsh, "Wrong index.");
		}

	///	checks in debug mode that multi-index is valid
		inline void check_multi_index(const MathVector<dim,int>& ind) const
		{
			UG_ASSERT(ind[0] <= (int)p-ind[2] && ind[0] >= 0, "Wrong Multiindex.");
			UG_ASSERT(ind[1] <= (int)p-ind[2] && ind[1] >= 0, "Wrong Multiindex.");
			UG_ASSERT(ind[2] <= (int)p && ind[2] >= 0, "Wrong Multiindex.");
		}

	private:

		MathVector<dim,int> m_vMultiIndex[nsh];
};


} //namespace ug

#endif /* __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__LAGRANGE__LAGRANGE__ */

