/*
 * lagrange.cpp
 *
 *  Created on: 01.07.2011
 *      Author: andreasvogel
 */

#include "lagrange.h"
#include <sstream>
#include "common/util/provider.h"

namespace ug{

template <typename TRefElem>
void SetLagrangeVertexMultiIndex(	MathVector<TRefElem::dim,int>* vMultiIndex,
                        	const TRefElem& rRef,
                        	size_t p,
                        	size_t& index)
{
//	dimension of Reference element
	static const int dim = TRefElem::dim;

//	get corner position integer
	const MathVector<dim,int>* vCo = rRef.corner();

//	loop all vertices
	for(size_t i = 0; i< rRef.num(0); ++i)
	{
	//	set multi index
		for(int d = 0; d<dim; ++d)
		{
			UG_ASSERT(((long)p)*vCo[i][d] >= 0,
			          "Wrong multi index m["<<d<<"]="<<p*vCo[i][d]);
			vMultiIndex[index][d] = p*vCo[i][d];
		}

	//	next index
		++index;
	}
}

template <typename TRefElem>
void SetLagrangeEdgeMultiIndex(	MathVector<TRefElem::dim,int>* vMultiIndex,
                      	const TRefElem& rRef,
                      	size_t p,
                      	size_t& index)
{
//	dimension of Reference element
	static const int dim = TRefElem::dim;

//	only for 2d,3d elems we do something
	if(dim < 1) return;

//	get corner position integer
	const MathVector<dim,int>* vCo = rRef.corner();

//	loop all edges
	for(size_t e = 0; e< rRef.num(1); ++e)
	{
	//	get ids of corners of edge
		const size_t co0 = rRef.id(1,e, 0,0);
		const size_t co1 = rRef.id(1,e, 0,1);

	//	add dofs on the edge
		for(size_t i = 1; i < p; ++i)
		{
		//	compute multi index
			MathVector<dim,int> m;
			VecScaleAdd(m, i, vCo[co1], (p-i), vCo[co0]);

		//	set multi index
			for(int d = 0; d<dim; ++d)
			{
				UG_ASSERT(m[d] >= 0, "Wrong multi index m["<<d<<"]="<<m[d]);
				vMultiIndex[index][d] = m[d];
			}

		//	next index
			index++;
		}
	}
}

template <typename TRefElem>
void SetLagrangeFaceMultiIndex(	MathVector<TRefElem::dim,int>* vMultiIndex,
                      	const TRefElem& rRef,
                      	size_t p,
                      	size_t& index)
{
//	dimension of Reference element
	static const int dim = TRefElem::dim;

//	only for 2d,3d elems we do something
	if(dim < 2) return;

//	get corner position integer
	const MathVector<dim,int>* vCo = rRef.corner();

//	add dof on quadrilateral
	for(size_t f = 0; f< rRef.num(2); ++f)
	{
	//	get ids of corners of edge
		const size_t co0 = rRef.id(2,f, 0,0);
		const size_t co1 = rRef.id(2,f, 0,1);
		size_t co2 = rRef.id(2,f, 0,2);
		if(rRef.num(2,f,0)==4) co2 = rRef.id(2,f, 0,3);

	//	directions of counting
		MathVector<dim,int> d1, d2;
		VecSubtract(d1, vCo[co1], vCo[co0]);
		VecSubtract(d2, vCo[co2], vCo[co0]);

	//	loop 'y'-direction
		for(size_t j = 1; j < p; ++j)
		{
		//	for a quadrilateral we have a quadratic loop, but for a
		//	triangle we need to stop at the diagonal
			const size_t off = ((rRef.num(2,f,0)==3) ? j : 0);

		//	loop 'x'-direction
			for(size_t i = 1; i < p-off; ++i)
			{
			//	compute multi index
				MathVector<dim,int> m = vCo[co0]; m *= p;
				VecScaleAppend(m, i, d1, j, d2);

			//	set multi index
				for(int d = 0; d<dim; ++d)
				{
					UG_ASSERT(m[d] >= 0, "Wrong multi index m["<<d<<"]="<<m[d]);
					vMultiIndex[index][d] = m[d];
				}
				
			//	next index
				++index;
			}
		}
	}
}

template <typename TRefElem>
void SetLagrangeVolumeMultiIndex(	MathVector<TRefElem::dim,int>* vMultiIndex,
                                 	const TRefElem& rRef,
                                 	size_t p,
                                 	size_t& index)
{
//	dimension of Reference element
	static const int dim = TRefElem::dim;

//	only for 3d elems we do something
	if(dim < 3) return;

//	get corner position integer
//	const MathVector<dim,int>* vCo = rRef.corner();

//	get type of reference element
	ReferenceObjectID type = rRef.roid(dim, 0);

//	handle elems
	switch(type)
	{
	case ROID_TETRAHEDRON:
		for(size_t m2 = 1; m2 < p; ++m2)
			for(size_t m1 = 1; m1 < p-m2; ++m1)
				for(size_t m0 = 1; m0 < p-m2-m1; ++m0)
				{
				//	use regular mapping for inner DoFs
					vMultiIndex[index][0] = m0;
					vMultiIndex[index][1] = m1;
					vMultiIndex[index++][2] = m2;
				}
		break;

	case ROID_PYRAMID:
		if(p>1) UG_THROW("LagrangeLSFS: Higher order Pyramid not implemented.");
		break;

	case ROID_PRISM:
		for(size_t m2 = 1; m2 < p; ++m2)
			for(size_t m1 = 1; m1 < p; ++m1)
				for(size_t m0 = 1; m0 < p-m1; ++m0)
				{
				//	use regular mapping for inner DoFs
					vMultiIndex[index][0] = m0;
					vMultiIndex[index][1] = m1;
					vMultiIndex[index++][2] = m2;
				}
		break;

	case ROID_HEXAHEDRON:
		for(size_t m2 = 1; m2 < p; ++m2)
			for(size_t m1 = 1; m1 < p; ++m1)
				for(size_t m0 = 1; m0 < p; ++m0)
				{
				//	use regular mapping for inner DoFs
					vMultiIndex[index][0] = m0;
					vMultiIndex[index][1] = m1;
					vMultiIndex[index++][2] = m2;
				}
		break;
	default: UG_THROW("LagrangeLSFS: Missing 3d mapping for type '"<<type<<"'."
	                        " roid="<<rRef.reference_object_id());
	}
}

template <typename TRefElem>
void SetLagrangeMultiIndex(	MathVector<TRefElem::dim,int>* vMultiIndex,
                           	const TRefElem& rRef,
                           	size_t p)
{
//	init shape -> multi-index mapping
	size_t index = 0;

//	vertices
	SetLagrangeVertexMultiIndex(vMultiIndex, rRef, p, index);

//	edge
	SetLagrangeEdgeMultiIndex(vMultiIndex, rRef, p, index);

//	faces
	SetLagrangeFaceMultiIndex(vMultiIndex, rRef, p, index);

//	volumes
	SetLagrangeVolumeMultiIndex(vMultiIndex, rRef, p, index);
}

///////////////////////////////////////////////////////////////////////////////
// Edge
///////////////////////////////////////////////////////////////////////////////

template <int TOrder>
LagrangeLSFS<ReferenceEdge, TOrder>::LagrangeLSFS()
	: LagrangeLDS<ReferenceEdge>(p)
{
//	init polynomials
	for(size_t i = 0; i < nsh; ++i)
	{
	//	create equidistant polynomials and its derivatives
		m_vPolynom[i] = EquidistantLagrange1D(i, p);
		m_vDPolynom[i] = m_vPolynom[i].derivative();
	}

//	reference element
	const ReferenceEdge& rRef =	Provider<ReferenceEdge>::get();

//	init shape -> multi-index mapping
	SetLagrangeMultiIndex(m_vMultiIndex, rRef, p);
}

void FlexLagrangeLSFS<ReferenceEdge>::set_order(size_t order)
{
	LagrangeLDS<ReferenceEdge>::set_order(order);

//	resize
	p = order;
	nsh = p+1;
	m_vPolynom.resize(nsh);
	m_vDPolynom.resize(nsh);
	m_vMultiIndex.resize(nsh);

//	init polynomials
	for(size_t i = 0; i < nsh; ++i)
	{
	//	create equidistant polynomials and its derivatives
		m_vPolynom[i] = EquidistantLagrange1D(i, p);
		m_vDPolynom[i] = m_vPolynom[i].derivative();
	}

//	reference element
	const ReferenceEdge& rRef = Provider<ReferenceEdge>::get();

//	init shape -> multi-index mapping
	SetLagrangeMultiIndex(&m_vMultiIndex[0], rRef, p);
}

template class LagrangeLSFS<ReferenceEdge, 1>;
template class LagrangeLSFS<ReferenceEdge, 2>;
template class LagrangeLSFS<ReferenceEdge, 3>;
template class LagrangeLSFS<ReferenceEdge, 4>;
template class LagrangeLSFS<ReferenceEdge, 5>;

//template class FlexLagrangeLSFS<ReferenceEdge>;

///////////////////////////////////////////////////////////////////////////////
// Triangle
///////////////////////////////////////////////////////////////////////////////

template <int TOrder>
LagrangeLSFS<ReferenceTriangle, TOrder>::LagrangeLSFS()
	: LagrangeLDS<ReferenceTriangle>(p)
{
//	init polynomials
	for(size_t i = 0; i <= p; ++i)
	{
	//	create trancated equidistant polynomials and its derivatives
		m_vPolynom[i] = TruncatedEquidistantLagrange1D(i, p);
		m_vDPolynom[i] = m_vPolynom[i].derivative();
	}

//	reference element
	const ReferenceTriangle& rRef = Provider<ReferenceTriangle>::get();

//	init shape -> multi-index mapping
	SetLagrangeMultiIndex(m_vMultiIndex, rRef, p);
}

void FlexLagrangeLSFS<ReferenceTriangle>::set_order(size_t order)
{
	LagrangeLDS<ReferenceTriangle>::set_order(order);

//	resize
	p = order;
	nsh = BinomCoeff(dim+p, p);

	const size_t polSize = p+1;
	m_vPolynom.resize(polSize);
	m_vDPolynom.resize(polSize);
	m_vMultiIndex.resize(nsh);

//	init polynomials
	for(size_t i = 0; i <= p; ++i)
	{
	//	create trancated equidistant polynomials and its derivatives
		m_vPolynom[i] = TruncatedEquidistantLagrange1D(i, p);
		m_vDPolynom[i] = m_vPolynom[i].derivative();
	}

//	reference element
	const ReferenceTriangle& rRef = Provider<ReferenceTriangle>::get();

//	init shape -> multi-index mapping
	SetLagrangeMultiIndex(&m_vMultiIndex[0], rRef, p);
}

template class LagrangeLSFS<ReferenceTriangle, 1>;
template class LagrangeLSFS<ReferenceTriangle, 2>;
template class LagrangeLSFS<ReferenceTriangle, 3>;
template class LagrangeLSFS<ReferenceTriangle, 4>;
template class LagrangeLSFS<ReferenceTriangle, 5>;

//template class FlexLagrangeLSFS<ReferenceTriangle>;

///////////////////////////////////////////////////////////////////////////////
// Quadrilateral
///////////////////////////////////////////////////////////////////////////////

template <int TOrder>
LagrangeLSFS<ReferenceQuadrilateral, TOrder>::LagrangeLSFS()
	: LagrangeLDS<ReferenceQuadrilateral>(p)
{
//	init polynomials
	for(size_t i = 0; i <= p; ++i)
	{
	//	create truncated equidistant polynomials and its derivatives
		m_vPolynom[i] = EquidistantLagrange1D(i, p);
		m_vDPolynom[i] = m_vPolynom[i].derivative();
	}

//	reference element
	const ReferenceQuadrilateral& rRef = Provider<ReferenceQuadrilateral>::get();

//	init shape -> multi-index mapping
	SetLagrangeMultiIndex(m_vMultiIndex, rRef, p);
}

void FlexLagrangeLSFS<ReferenceQuadrilateral>::set_order(size_t order)
{
	LagrangeLDS<ReferenceQuadrilateral>::set_order(order);

//	resize
	p = order;
	nsh = (p+1)*(p+1);

	const size_t polSize = p+1;
	m_vPolynom.resize(polSize);
	m_vDPolynom.resize(polSize);
	m_vMultiIndex.resize(nsh);

//	init polynomials
	for(size_t i = 0; i <= p; ++i)
	{
	//	create truncated equidistant polynomials and its derivatives
		m_vPolynom[i] = EquidistantLagrange1D(i, p);
		m_vDPolynom[i] = m_vPolynom[i].derivative();
	}

//	reference element
	const ReferenceQuadrilateral& rRef = Provider<ReferenceQuadrilateral>::get();

//	init shape -> multi-index mapping
	SetLagrangeMultiIndex(&m_vMultiIndex[0], rRef, p);
}

template class LagrangeLSFS<ReferenceQuadrilateral, 1>;
template class LagrangeLSFS<ReferenceQuadrilateral, 2>;
template class LagrangeLSFS<ReferenceQuadrilateral, 3>;
template class LagrangeLSFS<ReferenceQuadrilateral, 4>;
template class LagrangeLSFS<ReferenceQuadrilateral, 5>;

//template class FlexLagrangeLSFS<ReferenceQuadrilateral>;

///////////////////////////////////////////////////////////////////////////////
// Tetrahedron
///////////////////////////////////////////////////////////////////////////////

template <int TOrder>
LagrangeLSFS<ReferenceTetrahedron, TOrder>::LagrangeLSFS()
	: LagrangeLDS<ReferenceTetrahedron>(p)
{
	for(size_t i = 0; i <= p; ++i)
	{
	//	create trancated equidistant polynomials and its derivatives
		m_vPolynom[i] = TruncatedEquidistantLagrange1D(i, p);
		m_vDPolynom[i] = m_vPolynom[i].derivative();
	}

//	reference element
	const ReferenceTetrahedron& rRef = Provider<ReferenceTetrahedron>::get();

//	init shape -> multi-index mapping
	SetLagrangeMultiIndex(m_vMultiIndex, rRef, p);
}

void FlexLagrangeLSFS<ReferenceTetrahedron>::set_order(size_t order)
{
	LagrangeLDS<ReferenceTetrahedron>::set_order(order);

//	resize
	p = order;
	nsh = BinomCoeff(dim + p, p);

	const size_t polSize = p+1;
	m_vPolynom.resize(polSize);
	m_vDPolynom.resize(polSize);
	m_vMultiIndex.resize(nsh);

//	init polynomials
	for(size_t i = 0; i <= p; ++i)
	{
	//	create trancated equidistant polynomials and its derivatives
		m_vPolynom[i] = TruncatedEquidistantLagrange1D(i, p);
		m_vDPolynom[i] = m_vPolynom[i].derivative();
	}

//	reference element
	const ReferenceTetrahedron& rRef = Provider<ReferenceTetrahedron>::get();

//	init shape -> multi-index mapping
	SetLagrangeMultiIndex(&m_vMultiIndex[0], rRef, p);
}

template class LagrangeLSFS<ReferenceTetrahedron, 1>;
template class LagrangeLSFS<ReferenceTetrahedron, 2>;
template class LagrangeLSFS<ReferenceTetrahedron, 3>;
template class LagrangeLSFS<ReferenceTetrahedron, 4>;
template class LagrangeLSFS<ReferenceTetrahedron, 5>;

//template class FlexLagrangeLSFS<ReferenceTetrahedron>;

///////////////////////////////////////////////////////////////////////////////
// Prism
///////////////////////////////////////////////////////////////////////////////

template <int TOrder>
LagrangeLSFS<ReferencePrism, TOrder>::LagrangeLSFS()
	: LagrangeLDS<ReferencePrism>(p)
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

//	reference element
	const ReferencePrism& rRef = Provider<ReferencePrism>::get();

//	init shape -> multi-index mapping
	SetLagrangeMultiIndex(m_vMultiIndex, rRef, p);
}

void FlexLagrangeLSFS<ReferencePrism>::set_order(size_t order)
{
	LagrangeLDS<ReferencePrism>::set_order(order);

//	resize
	p = order;
	dofPerLayer = BinomCoeff(2 + p, p);
	nsh = dofPerLayer * (p+1);

	const size_t polSize = p+1;
	m_vPolynom.resize(polSize);
	m_vDPolynom.resize(polSize);
	m_vMultiIndex.resize(nsh);

//	init polynomials
	for(size_t i = 0; i <= p; ++i)
	{
	//	create truncated equidistant polynomials and its derivatives
		m_vTruncPolynom[i] = TruncatedEquidistantLagrange1D(i, p);
		m_vDTruncPolynom[i] = m_vTruncPolynom[i].derivative();

	//	create equidistant polynomials and its derivatives
		m_vPolynom[i] = EquidistantLagrange1D(i, p);
		m_vDPolynom[i] = m_vPolynom[i].derivative();
	}

//	reference element
	const ReferencePrism& rRef = Provider<ReferencePrism>::get();

//	init shape -> multi-index mapping
	SetLagrangeMultiIndex(&m_vMultiIndex[0], rRef, p);
}

template class LagrangeLSFS<ReferencePrism, 1>;
template class LagrangeLSFS<ReferencePrism, 2>;
template class LagrangeLSFS<ReferencePrism, 3>;
template class LagrangeLSFS<ReferencePrism, 4>;
template class LagrangeLSFS<ReferencePrism, 5>;

//template class FlexLagrangeLSFS<ReferencePrism>;

///////////////////////////////////////////////////////////////////////////////
// Pyramid
///////////////////////////////////////////////////////////////////////////////

template <int TOrder>
LagrangeLSFS<ReferencePyramid, TOrder>::LagrangeLSFS()
	: LagrangeLDS<ReferencePyramid>(p)
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

	//	reference element
		const ReferencePyramid& rRef =
				Provider<ReferencePyramid>::get();

	//	init shape -> multi-index mapping
		SetLagrangeMultiIndex(m_vMultiIndex, rRef, p);
}

template class LagrangeLSFS<ReferencePyramid, 1>;
template class LagrangeLSFS<ReferencePyramid, 2>;

///////////////////////////////////////////////////////////////////////////////
// Hexahedron
///////////////////////////////////////////////////////////////////////////////

template <int TOrder>
LagrangeLSFS<ReferenceHexahedron, TOrder>::LagrangeLSFS()
	: LagrangeLDS<ReferenceHexahedron>(p)
{
//	init polynomials
	for(size_t i = 0; i <= p; ++i)
	{
	//	create trancated equidistant polynomials and its derivatives
		m_vPolynom[i] = EquidistantLagrange1D(i, p);
		m_vDPolynom[i] = m_vPolynom[i].derivative();
	}

//	reference element
	const ReferenceHexahedron& rRef = Provider<ReferenceHexahedron>::get();

//	init shape -> multi-index mapping
	SetLagrangeMultiIndex(m_vMultiIndex, rRef, p);
}

void FlexLagrangeLSFS<ReferenceHexahedron>::set_order(size_t order)
{
	LagrangeLDS<ReferenceHexahedron>::set_order(order);

//	resize
	p = order;
	nsh = (p+1) * (p+1) * (p+1);

	const size_t polSize = p+1;
	m_vPolynom.resize(polSize);
	m_vDPolynom.resize(polSize);
	m_vMultiIndex.resize(nsh);

//	init polynomials
	for(size_t i = 0; i <= p; ++i)
	{
	//	create trancated equidistant polynomials and its derivatives
		m_vPolynom[i] = EquidistantLagrange1D(i, p);
		m_vDPolynom[i] = m_vPolynom[i].derivative();
	}

//	reference element
	const ReferenceHexahedron& rRef = Provider<ReferenceHexahedron>::get();

//	init shape -> multi-index mapping
	SetLagrangeMultiIndex(&m_vMultiIndex[0], rRef, p);
}

template class LagrangeLSFS<ReferenceHexahedron, 1>;
template class LagrangeLSFS<ReferenceHexahedron, 2>;
template class LagrangeLSFS<ReferenceHexahedron, 3>;
template class LagrangeLSFS<ReferenceHexahedron, 4>;
template class LagrangeLSFS<ReferenceHexahedron, 5>;

} // end namespace ug
