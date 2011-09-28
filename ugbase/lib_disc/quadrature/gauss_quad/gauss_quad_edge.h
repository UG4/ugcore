//  This file is parsed from UG 3.9.
//  It provides the Gauss Quadratures for a reference edge.


#ifndef __H__UG__LIB_DISC__QUADRATURE__GAUSS_QUADRATURE__ReferenceEdge__
#define __H__UG__LIB_DISC__QUADRATURE__GAUSS_QUADRATURE__ReferenceEdge__

#include "gauss_quad.h"

namespace ug{

template <>
class GaussQuadrature<ReferenceEdge, 1>
{
	public:
	/// Dimension of integration domain
		static const size_t dim = 1;

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
class GaussQuadrature<ReferenceEdge, 3>
{
	public:
	/// Dimension of integration domain
		static const size_t dim = 1;

	/// Order of quadrature rule
		static const size_t p = 3;

	/// Number of integration points
		static const size_t nip = 2;

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
class GaussQuadrature<ReferenceEdge, 5>
{
	public:
	/// Dimension of integration domain
		static const size_t dim = 1;

	/// Order of quadrature rule
		static const size_t p = 5;

	/// Number of integration points
		static const size_t nip = 3;

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
class GaussQuadrature<ReferenceEdge, 7>
{
	public:
	/// Dimension of integration domain
		static const size_t dim = 1;

	/// Order of quadrature rule
		static const size_t p = 7;

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
class GaussQuadrature<ReferenceEdge, 9>
{
	public:
	/// Dimension of integration domain
		static const size_t dim = 1;

	/// Order of quadrature rule
		static const size_t p = 9;

	/// Number of integration points
		static const size_t nip = 5;

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
class GaussQuadrature<ReferenceEdge, 11>
{
	public:
	/// Dimension of integration domain
		static const size_t dim = 1;

	/// Order of quadrature rule
		static const size_t p = 11;

	/// Number of integration points
		static const size_t nip = 6;

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
class GaussQuadrature<ReferenceEdge, 13>
{
	public:
	/// Dimension of integration domain
		static const size_t dim = 1;

	/// Order of quadrature rule
		static const size_t p = 13;

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
class GaussQuadrature<ReferenceEdge, 15>
{
	public:
	/// Dimension of integration domain
		static const size_t dim = 1;

	/// Order of quadrature rule
		static const size_t p = 15;

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
class GaussQuadrature<ReferenceEdge, 17>
{
	public:
	/// Dimension of integration domain
		static const size_t dim = 1;

	/// Order of quadrature rule
		static const size_t p = 17;

	/// Number of integration points
		static const size_t nip = 9;

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
class GaussQuadrature<ReferenceEdge, 19>
{
	public:
	/// Dimension of integration domain
		static const size_t dim = 1;

	/// Order of quadrature rule
		static const size_t p = 19;

	/// Number of integration points
		static const size_t nip = 10;

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
class GaussQuadrature<ReferenceEdge, 18> : public GaussQuadrature<ReferenceEdge, 19>{};

template <>
class GaussQuadrature<ReferenceEdge, 16> : public GaussQuadrature<ReferenceEdge, 17>{};

template <>
class GaussQuadrature<ReferenceEdge, 14> : public GaussQuadrature<ReferenceEdge, 15>{};

template <>
class GaussQuadrature<ReferenceEdge, 12> : public GaussQuadrature<ReferenceEdge, 13>{};

template <>
class GaussQuadrature<ReferenceEdge, 10> : public GaussQuadrature<ReferenceEdge, 11>{};

template <>
class GaussQuadrature<ReferenceEdge, 8> : public GaussQuadrature<ReferenceEdge, 9>{};

template <>
class GaussQuadrature<ReferenceEdge, 6> : public GaussQuadrature<ReferenceEdge, 7>{};

template <>
class GaussQuadrature<ReferenceEdge, 4> : public GaussQuadrature<ReferenceEdge, 5>{};

template <>
class GaussQuadrature<ReferenceEdge, 2> : public GaussQuadrature<ReferenceEdge, 3>{};

template <>
class GaussQuadrature<ReferenceEdge, 0> : public GaussQuadrature<ReferenceEdge, 1>{};

}; // namespace ug

#endif /* __H__UG__LIB_DISC__QUADRATURE__GAUSS_QUADRATURE__ReferenceEdge__ */

