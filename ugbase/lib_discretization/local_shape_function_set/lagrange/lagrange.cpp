/*
 * lagrange.cpp
 *
 *  Created on: 01.07.2011
 *      Author: andreasvogel
 */

#include "lagrange.h"

namespace ug{

///////////////////////////////////////////////////////////////////////////////
// Edge
///////////////////////////////////////////////////////////////////////////////

template <>
template <int TOrder>
LagrangeLSFS<ReferenceEdge, TOrder>::LagrangeLSFS()
{
//	init polynomials
	for(size_t i = 0; i < nsh; ++i)
	{
	//	create equidistant polynomials and its derivatives
		m_vPolynom[i] = EquidistantLagrange1D(i, p);
		m_vDPolynom[i] = m_vPolynom[i].derivative();
	}

//	reference element
	const ReferenceEdge& rRef =
			ReferenceElementProvider::get<ReferenceEdge>();

//	init shape -> multi-index mapping
	size_t index = 0;

//	vertices
	for(size_t i = 0; i< rRef.num(0); ++i)
		m_vLocalDoF[index++] = LocalDoF(0, i, 0);
	m_vMultiIndex[0] = MultiIndex<dim>(0);
	m_vMultiIndex[1] = MultiIndex<dim>(p);

//	edge
	for(size_t i = 1; i < p; ++i)
	{
		m_vLocalDoF[index] = LocalDoF(1, 0, i-1);
		m_vMultiIndex[index++] = MultiIndex<dim>(i);
	}

	UG_ASSERT(index == nsh, "Not all indices distributed.");
}

template class LagrangeLSFS<ReferenceEdge, 1>;
template class LagrangeLSFS<ReferenceEdge, 2>;
template class LagrangeLSFS<ReferenceEdge, 3>;
template class LagrangeLSFS<ReferenceEdge, 4>;
template class LagrangeLSFS<ReferenceEdge, 5>;

///////////////////////////////////////////////////////////////////////////////
// Triangle
///////////////////////////////////////////////////////////////////////////////

template <>
template <int TOrder>
LagrangeLSFS<ReferenceTriangle, TOrder>::LagrangeLSFS()
{
//	init polynomials
	for(size_t i = 0; i <= p; ++i)
	{
	//	create trancated equidistant polynomials and its derivatives
		m_vPolynom[i] = TruncatedEquidistantLagrange1D(i, p);
		m_vDPolynom[i] = m_vPolynom[i].derivative();
	}

//	reference element
	const ReferenceTriangle& rRef =
			ReferenceElementProvider::get<ReferenceTriangle>();

//	init shape -> multi-index mapping
	size_t index = 0;

	int vCo[3][dim] = {{0,0}, {1,0}, {0,1}};

//	vertices
	for(size_t i = 0; i< rRef.num(0); ++i)
	{
		m_vLocalDoF[index] = LocalDoF(0, i, 0);
		m_vMultiIndex[index++] = MultiIndex<dim>(p*vCo[i][0],p*vCo[i][1]);
	}

//	edge
	for(size_t e = 0; e< rRef.num(1); ++e)
	{
	//	get ids of corners of edge
		const size_t co0 = rRef.id(1,e, 0,0);
		const size_t co1 = rRef.id(1,e, 0,1);

	//	add dofs on the edge
		for(size_t i = 1; i < p; ++i)
		{
		//	set: dim=1, id=e, offset=i-1
			m_vLocalDoF[index] = LocalDoF(1, e, i-1);

		//	compute multi index
			const size_t m0 = i*vCo[co0][0] + (p-i)*vCo[co1][0];
			const size_t m1 = i*vCo[co0][1] + (p-i)*vCo[co1][1];

			m_vMultiIndex[index++] = MultiIndex<dim>(m0, m1);
		}
	}

//	add dof on triangle
	const size_t used = index;
	for(size_t m1 = 1; m1 < p; ++m1)
	{
		for(size_t m0 = 1; m0 < p-m1; ++m0)
		{
		//	set: dim=2, id=0, offset=i
			m_vLocalDoF[index] = LocalDoF(2, 0, index-used);

		//	use regular mapping for inner DoFs
			m_vMultiIndex[index++] = MultiIndex<dim>(m0,m1);
		}
	}

	UG_ASSERT(index == nsh, "Not all indices distributed.");
}

template class LagrangeLSFS<ReferenceTriangle, 1>;
template class LagrangeLSFS<ReferenceTriangle, 2>;
template class LagrangeLSFS<ReferenceTriangle, 3>;
template class LagrangeLSFS<ReferenceTriangle, 4>;
template class LagrangeLSFS<ReferenceTriangle, 5>;

///////////////////////////////////////////////////////////////////////////////
// Quadrilateral
///////////////////////////////////////////////////////////////////////////////

template <>
template <int TOrder>
LagrangeLSFS<ReferenceQuadrilateral, TOrder>::LagrangeLSFS()
{
//	init polynomials
	for(size_t i = 0; i <= p; ++i)
	{
	//	create trancated equidistant polynomials and its derivatives
		m_vPolynom[i] = EquidistantLagrange1D(i, p);
		m_vDPolynom[i] = m_vPolynom[i].derivative();
	}

//	reference element
	const ReferenceQuadrilateral& rRef =
			ReferenceElementProvider::get<ReferenceQuadrilateral>();

//	init shape -> multi-index mapping
	size_t index = 0;

	int vCo[4][dim] = {{0,0}, {1,0}, {1,1}, {0,1}};

//	vertices
	for(size_t i = 0; i< rRef.num(0); ++i)
	{
		m_vLocalDoF[index] = LocalDoF(0, i, 0);
		m_vMultiIndex[index++] = MultiIndex<dim>(p*vCo[i][0],p*vCo[i][1]);
	}

//	edge
	for(size_t e = 0; e< rRef.num(1); ++e)
	{
	//	get ids of corners of edge
		const size_t co0 = rRef.id(1,e, 0,0);
		const size_t co1 = rRef.id(1,e, 0,1);

	//	add dofs on the edge
		for(size_t i = 1; i < p; ++i)
		{
		//	set: dim=1, id=e, offset=i-1
			m_vLocalDoF[index] = LocalDoF(1, e, i-1);

		//	compute multi index
			const size_t m0 = i*vCo[co0][0] + (p-i)*vCo[co1][0];
			const size_t m1 = i*vCo[co0][1] + (p-i)*vCo[co1][1];

			m_vMultiIndex[index++] = MultiIndex<dim>(m0, m1);
		}
	}

//	add dof on triangle
	const size_t used = index;
	for(size_t m1 = 1; m1 < p; ++m1)
	{
		for(size_t m0 = 1; m0 < p; ++m0)
		{
		//	set: dim=2, id=0, offset=i
			m_vLocalDoF[index] = LocalDoF(2, 0, index-used);

		//	use regular mapping for inner DoFs
			m_vMultiIndex[index++] = MultiIndex<dim>(m0,m1);
		}
	}

	UG_ASSERT(index == nsh, "Not all indices distributed.");
}

template class LagrangeLSFS<ReferenceQuadrilateral, 1>;
template class LagrangeLSFS<ReferenceQuadrilateral, 2>;
template class LagrangeLSFS<ReferenceQuadrilateral, 3>;
template class LagrangeLSFS<ReferenceQuadrilateral, 4>;
template class LagrangeLSFS<ReferenceQuadrilateral, 5>;

///////////////////////////////////////////////////////////////////////////////
// Tetrahedron
///////////////////////////////////////////////////////////////////////////////

template <>
template <int TOrder>
LagrangeLSFS<ReferenceTetrahedron, TOrder>::LagrangeLSFS()
{
	for(size_t i = 0; i <= p; ++i)
	{
	//	create trancated equidistant polynomials and its derivatives
		m_vPolynom[i] = TruncatedEquidistantLagrange1D(i, p);
		m_vDPolynom[i] = m_vPolynom[i].derivative();
	}
}

template class LagrangeLSFS<ReferenceTetrahedron, 1>;
template class LagrangeLSFS<ReferenceTetrahedron, 2>;
template class LagrangeLSFS<ReferenceTetrahedron, 3>;
template class LagrangeLSFS<ReferenceTetrahedron, 4>;
template class LagrangeLSFS<ReferenceTetrahedron, 5>;

///////////////////////////////////////////////////////////////////////////////
// Prism
///////////////////////////////////////////////////////////////////////////////

template <>
template <int TOrder>
LagrangeLSFS<ReferencePrism, TOrder>::LagrangeLSFS()
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

template class LagrangeLSFS<ReferencePrism, 1>;
template class LagrangeLSFS<ReferencePrism, 2>;
template class LagrangeLSFS<ReferencePrism, 3>;
template class LagrangeLSFS<ReferencePrism, 4>;
template class LagrangeLSFS<ReferencePrism, 5>;

///////////////////////////////////////////////////////////////////////////////
// Pyramid
///////////////////////////////////////////////////////////////////////////////

template <>
template <int TOrder>
LagrangeLSFS<ReferencePyramid, TOrder>::LagrangeLSFS()
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

template class LagrangeLSFS<ReferencePyramid, 1>;
template class LagrangeLSFS<ReferencePyramid, 2>;

///////////////////////////////////////////////////////////////////////////////
// Hexahedron
///////////////////////////////////////////////////////////////////////////////

template <>
template <int TOrder>
LagrangeLSFS<ReferenceHexahedron, TOrder>::LagrangeLSFS()
{
//	init polynomials
	for(size_t i = 0; i <= p; ++i)
	{
	//	create trancated equidistant polynomials and its derivatives
		m_vPolynom[i] = EquidistantLagrange1D(i, p);
		m_vDPolynom[i] = m_vPolynom[i].derivative();
	}

//	reference element
	const ReferenceHexahedron& rRef =
			ReferenceElementProvider::get<ReferenceHexahedron>();

//	init shape -> multi-index mapping
	size_t index = 0;

	int vCo[8][dim] = {{0,0,0}, {1,0,0}, {1,1,0}, {0,1,0},
	                   {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}};

//	vertices
	for(size_t i = 0; i< rRef.num(0); ++i)
	{
		m_vLocalDoF[index] = LocalDoF(0, i, 0);
		m_vMultiIndex[index++] = MultiIndex<dim>(p*vCo[i][0],p*vCo[i][1],p*vCo[i][2]);
	}

//	edge
	for(size_t e = 0; e< rRef.num(1); ++e)
	{
	//	get ids of corners of edge
		const size_t co0 = rRef.id(1,e, 0,0);
		const size_t co1 = rRef.id(1,e, 0,1);

	//	add dofs on the edge
		for(size_t i = 1; i < p; ++i)
		{
		//	set: dim=1, id=e, offset=i-1
			m_vLocalDoF[index] = LocalDoF(1, e, i-1);

		//	compute multi index
			const size_t m0 = i*vCo[co0][0] + (p-i)*vCo[co1][0];
			const size_t m1 = i*vCo[co0][1] + (p-i)*vCo[co1][1];
			const size_t m2 = i*vCo[co0][2] + (p-i)*vCo[co1][2];

			m_vMultiIndex[index++] = MultiIndex<dim>(m0, m1, m2);
		}
	}

//	add dof on quadrilateral
	for(size_t f = 0; f< rRef.num(2); ++f)
	{
	//	get ids of corners of edge
		const size_t co0 = rRef.id(2,f, 0,0);
		const size_t co1 = rRef.id(2,f, 0,1);
		const size_t co2 = rRef.id(2,f, 0,2);
		//const size_t co3 = rRef.id(2,f, 0,3);

		int p0[dim], d1[dim], d2[dim];

	//	origin of counting
		p0[0] = vCo[co0][0];
		p0[1] = vCo[co0][1];
		p0[2] = vCo[co0][2];

	//	directions of counting
		d1[0] = vCo[co1][0] - vCo[co0][0];
		d1[1] = vCo[co1][1] - vCo[co0][1];
		d1[2] = vCo[co1][2] - vCo[co0][2];

		d2[0] = vCo[co2][0] - vCo[co0][0];
		d2[1] = vCo[co2][1] - vCo[co0][1];
		d2[2] = vCo[co2][2] - vCo[co0][2];

	//	add dofs on the edge
		size_t cnt = 0;
		for(size_t j = 1; j < p; ++j)
		{
			for(size_t i = 1; i < p; ++i)
			{
			//	set: dim=2, id=f, offset=cnt
				m_vLocalDoF[index] = LocalDoF(2, f, cnt++);

			//	compute multi index
				const size_t m0 = p0[0] + i*d1[0] + j*d2[0];
				const size_t m1 = p0[1] + i*d1[1] + j*d2[1];
				const size_t m2 = p0[2] + i*d1[2] + j*d2[2];

				m_vMultiIndex[index++] = MultiIndex<dim>(m0, m1, m2);
			}
		}
	}

//	add dof on hexahedron
	const size_t used = index;
	for(size_t m2 = 1; m2 < p; ++m2)
	{
		for(size_t m1 = 1; m1 < p; ++m1)
		{
			for(size_t m0 = 1; m0 < p; ++m0)
			{
			//	set: dim=2, id=0, offset=i
				m_vLocalDoF[index] = LocalDoF(3, 0, index-used);

			//	use regular mapping for inner DoFs
				m_vMultiIndex[index++] = MultiIndex<dim>(m0,m1,m2);
			}
		}
	}

	UG_ASSERT(index == nsh, "Not all indices distributed.");
}

template class LagrangeLSFS<ReferenceHexahedron, 1>;
template class LagrangeLSFS<ReferenceHexahedron, 2>;
template class LagrangeLSFS<ReferenceHexahedron, 3>;
template class LagrangeLSFS<ReferenceHexahedron, 4>;
template class LagrangeLSFS<ReferenceHexahedron, 5>;

} // end namespace ug
