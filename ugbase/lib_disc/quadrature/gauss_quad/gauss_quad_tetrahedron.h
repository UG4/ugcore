//  This file is parsed from UG 3.9.
//  It provides the Gauss Quadratures for a reference tetrahedron.


#ifndef __H__UG__LIB_DISCRETIZATION__QUADRATURE__GAUSS_QUADRATURE__ReferenceTetrahedron__
#define __H__UG__LIB_DISCRETIZATION__QUADRATURE__GAUSS_QUADRATURE__ReferenceTetrahedron__

#include "gauss_quad.h"

namespace ug{

template <>
class GaussQuadrature<ReferenceTetrahedron, 0>
{
	public:
	/// Dimension of integration domain
		static const size_t dim = 3;

	/// Order of quadrature rule
		static const size_t p = 0;

	/// Number of integration points
		static const size_t nip = 7;

	/// Constructor
		GaussQuadrature();

	/// number of integration points
		inline size_t size() const {return nip;}

	/// returns i'th integration point
		inline const MathVector<dim>& point(size_t i) const
			{UG_ASSERT(i < size(), "Wrong index"); return m_vPoint[i];}

	/// returns all positions in an array of size()
		inline const MathVector<dim>* points() const {return m_vPoint;}

	/// return the i'th weight
		inline number weight(size_t i) const
			{UG_ASSERT(i < size(), "Wrong index"); return m_vWeight[i];}

	/// returns all weights in an array of size()
		inline const number* weights() const {return m_vWeight;}

	/// returns the order
		inline size_t order() const {return p;}

	protected:
	/// integration points
		MathVector<dim> m_vPoint[nip];

	/// weights
		number m_vWeight[nip];
};

template <>
class GaussQuadrature<ReferenceTetrahedron, 1>
{
	public:
	/// Dimension of integration domain
		static const size_t dim = 3;

	/// Order of quadrature rule
		static const size_t p = 1;

	/// Number of integration points
		static const size_t nip = 1;

	/// Constructor
		GaussQuadrature();

	/// number of integration points
		inline size_t size() const {return nip;}

	/// returns i'th integration point
		inline const MathVector<dim>& point(size_t i) const
			{UG_ASSERT(i < size(), "Wrong index"); return m_vPoint[i];}

	/// returns all positions in an array of size()
		inline const MathVector<dim>* points() const {return m_vPoint;}

	/// return the i'th weight
		inline number weight(size_t i) const
			{UG_ASSERT(i < size(), "Wrong index"); return m_vWeight[i];}

	/// returns all weights in an array of size()
		inline const number* weights() const {return m_vWeight;}

	/// returns the order
		inline size_t order() const {return p;}

	protected:
	/// integration points
		MathVector<dim> m_vPoint[nip];

	/// weights
		number m_vWeight[nip];
};

template <>
class GaussQuadrature<ReferenceTetrahedron, 2>
{
	public:
	/// Dimension of integration domain
		static const size_t dim = 3;

	/// Order of quadrature rule
		static const size_t p = 2;

	/// Number of integration points
		static const size_t nip = 4;

	/// Constructor
		GaussQuadrature();

	/// number of integration points
		inline size_t size() const {return nip;}

	/// returns i'th integration point
		inline const MathVector<dim>& point(size_t i) const
			{UG_ASSERT(i < size(), "Wrong index"); return m_vPoint[i];}

	/// returns all positions in an array of size()
		inline const MathVector<dim>* points() const {return m_vPoint;}

	/// return the i'th weight
		inline number weight(size_t i) const
			{UG_ASSERT(i < size(), "Wrong index"); return m_vWeight[i];}

	/// returns all weights in an array of size()
		inline const number* weights() const {return m_vWeight;}

	/// returns the order
		inline size_t order() const {return p;}

	protected:
	/// integration points
		MathVector<dim> m_vPoint[nip];

	/// weights
		number m_vWeight[nip];
};

template <>
class GaussQuadrature<ReferenceTetrahedron, 3>
{
	public:
	/// Dimension of integration domain
		static const size_t dim = 3;

	/// Order of quadrature rule
		static const size_t p = 3;

	/// Number of integration points
		static const size_t nip = 8;

	/// Constructor
		GaussQuadrature();

	/// number of integration points
		inline size_t size() const {return nip;}

	/// returns i'th integration point
		inline const MathVector<dim>& point(size_t i) const
			{UG_ASSERT(i < size(), "Wrong index"); return m_vPoint[i];}

	/// returns all positions in an array of size()
		inline const MathVector<dim>* points() const {return m_vPoint;}

	/// return the i'th weight
		inline number weight(size_t i) const
			{UG_ASSERT(i < size(), "Wrong index"); return m_vWeight[i];}

	/// returns all weights in an array of size()
		inline const number* weights() const {return m_vWeight;}

	/// returns the order
		inline size_t order() const {return p;}

	protected:
	/// integration points
		MathVector<dim> m_vPoint[nip];

	/// weights
		number m_vWeight[nip];
};

template <>
class GaussQuadrature<ReferenceTetrahedron, 5>
{
	public:
	/// Dimension of integration domain
		static const size_t dim = 3;

	/// Order of quadrature rule
		static const size_t p = 5;

	/// Number of integration points
		static const size_t nip = 15;

	/// Constructor
		GaussQuadrature();

	/// number of integration points
		inline size_t size() const {return nip;}

	/// returns i'th integration point
		inline const MathVector<dim>& point(size_t i) const
			{UG_ASSERT(i < size(), "Wrong index"); return m_vPoint[i];}

	/// returns all positions in an array of size()
		inline const MathVector<dim>* points() const {return m_vPoint;}

	/// return the i'th weight
		inline number weight(size_t i) const
			{UG_ASSERT(i < size(), "Wrong index"); return m_vWeight[i];}

	/// returns all weights in an array of size()
		inline const number* weights() const {return m_vWeight;}

	/// returns the order
		inline size_t order() const {return p;}

	protected:
	/// integration points
		MathVector<dim> m_vPoint[nip];

	/// weights
		number m_vWeight[nip];
};

template <>
class GaussQuadrature<ReferenceTetrahedron, 6>
{
	public:
	/// Dimension of integration domain
		static const size_t dim = 3;

	/// Order of quadrature rule
		static const size_t p = 6;

	/// Number of integration points
		static const size_t nip = 24;

	/// Constructor
		GaussQuadrature();

	/// number of integration points
		inline size_t size() const {return nip;}

	/// returns i'th integration point
		inline const MathVector<dim>& point(size_t i) const
			{UG_ASSERT(i < size(), "Wrong index"); return m_vPoint[i];}

	/// returns all positions in an array of size()
		inline const MathVector<dim>* points() const {return m_vPoint;}

	/// return the i'th weight
		inline number weight(size_t i) const
			{UG_ASSERT(i < size(), "Wrong index"); return m_vWeight[i];}

	/// returns all weights in an array of size()
		inline const number* weights() const {return m_vWeight;}

	/// returns the order
		inline size_t order() const {return p;}

	protected:
	/// integration points
		MathVector<dim> m_vPoint[nip];

	/// weights
		number m_vWeight[nip];
};

template <>
class GaussQuadrature<ReferenceTetrahedron, 7>
{
	public:
	/// Dimension of integration domain
		static const size_t dim = 3;

	/// Order of quadrature rule
		static const size_t p = 7;

	/// Number of integration points
		static const size_t nip = 31;

	/// Constructor
		GaussQuadrature();

	/// number of integration points
		inline size_t size() const {return nip;}

	/// returns i'th integration point
		inline const MathVector<dim>& point(size_t i) const
			{UG_ASSERT(i < size(), "Wrong index"); return m_vPoint[i];}

	/// returns all positions in an array of size()
		inline const MathVector<dim>* points() const {return m_vPoint;}

	/// return the i'th weight
		inline number weight(size_t i) const
			{UG_ASSERT(i < size(), "Wrong index"); return m_vWeight[i];}

	/// returns all weights in an array of size()
		inline const number* weights() const {return m_vWeight;}

	/// returns the order
		inline size_t order() const {return p;}

	protected:
	/// integration points
		MathVector<dim> m_vPoint[nip];

	/// weights
		number m_vWeight[nip];
};

template <>
class GaussQuadrature<ReferenceTetrahedron, 8>
{
	public:
	/// Dimension of integration domain
		static const size_t dim = 3;

	/// Order of quadrature rule
		static const size_t p = 8;

	/// Number of integration points
		static const size_t nip = 43;

	/// Constructor
		GaussQuadrature();

	/// number of integration points
		inline size_t size() const {return nip;}

	/// returns i'th integration point
		inline const MathVector<dim>& point(size_t i) const
			{UG_ASSERT(i < size(), "Wrong index"); return m_vPoint[i];}

	/// returns all positions in an array of size()
		inline const MathVector<dim>* points() const {return m_vPoint;}

	/// return the i'th weight
		inline number weight(size_t i) const
			{UG_ASSERT(i < size(), "Wrong index"); return m_vWeight[i];}

	/// returns all weights in an array of size()
		inline const number* weights() const {return m_vWeight;}

	/// returns the order
		inline size_t order() const {return p;}

	protected:
	/// integration points
		MathVector<dim> m_vPoint[nip];

	/// weights
		number m_vWeight[nip];
};

template <>
class GaussQuadrature<ReferenceTetrahedron, 4> : public GaussQuadrature<ReferenceTetrahedron, 5>{};

}; // namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__QUADRATURE__GAUSS_QUADRATURE__ReferenceTetrahedron__ */

