/*
 * lagrange.cpp
 *
 *  Created on: 01.07.2011
 *      Author: andreasvogel
 */

#include "lagrange.h"
#include <sstream>

namespace ug{

template <typename TRefElem>
void SetVertexLocalDoFs(	MultiIndex<TRefElem::dim>* vMultiIndex,
                        	LocalDoF* vLocalDoF,
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
		vLocalDoF[index] = LocalDoF(0, i, 0);

	//	set multi index
		for(int d = 0; d<dim; ++d)
			vMultiIndex[index][d] = p*vCo[i][d];

	//	next index
		++index;
	}
}

template <typename TRefElem>
void SetEdgeLocalDoFs(	MultiIndex<TRefElem::dim>* vMultiIndex,
                      	LocalDoF* vLocalDoF,
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
		//	set: dim=1, id=e, offset=i-1
			vLocalDoF[index] = LocalDoF(1, e, i-1);

		//	compute multi index
			MathVector<dim,int> m;
			VecScaleAdd(m, i, vCo[co1], (p-i), vCo[co0]);

		//	set multi index
			for(int d = 0; d<dim; ++d)
				vMultiIndex[index][d] = m[d];

		//	next index
			index++;
		}
	}
}

template <typename TRefElem>
void SetFaceLocalDoFs(	MultiIndex<TRefElem::dim>* vMultiIndex,
                      	LocalDoF* vLocalDoF,
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

	//	reset counter
		size_t cnt = 0;

	//	loop 'y'-direction
		for(size_t j = 1; j < p; ++j)
		{
		//	for a quadrilateral we have a quadratic loop, but for a
		//	triangle we need to stop at the diagonal
			const size_t off = ((rRef.num(2,f,0)==3) ? j : 0);

		//	loop 'x'-direction
			for(size_t i = 1; i < p-off; ++i)
			{
			//	set: dim=2, id=f, offset=cnt
				vLocalDoF[index] = LocalDoF(2, f, cnt++);

			//	compute multi index
				MathVector<dim,int> m = vCo[co0];
				VecScaleAppend(m, i, d1, j, d2);

			//	set multi index
				for(int d = 0; d<dim; ++d)
					vMultiIndex[index][d] = m[d];

			//	next index
				++index;
			}
		}
	}
}

template <typename TRefElem>
void SetVolumeLocalDoFs(	MultiIndex<TRefElem::dim>* vMultiIndex,
                        	LocalDoF* vLocalDoF,
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
	ReferenceObjectID type = rRef.ref_elem_type(dim, 0);

//	handle elems
	size_t cnt = 0;
	switch(type)
	{
	case ROID_TETRAHEDRON:
		for(size_t m2 = 1; m2 < p; ++m2)
			for(size_t m1 = 1; m1 < p-m2; ++m1)
				for(size_t m0 = 1; m0 < p-m2-m1; ++m0)
				{
				//	set: dim=2, id=0, offset=i
					vLocalDoF[index] = LocalDoF(3, 0, cnt++);

				//	use regular mapping for inner DoFs
					vMultiIndex[index][0] = m0;
					vMultiIndex[index][1] = m1;
					vMultiIndex[index++][2] = m2;
				}
		break;

	case ROID_PYRAMID:
		if(p>1) throw(UGFatalError("LagrangeLSFS: Higher order Pyramid not implemented."));
		break;

	case ROID_PRISM:
		for(size_t m2 = 1; m2 < p; ++m2)
			for(size_t m1 = 1; m1 < p; ++m1)
				for(size_t m0 = 1; m0 < p-m1; ++m0)
				{
				//	set: dim=2, id=0, offset=i
					vLocalDoF[index] = LocalDoF(3, 0, cnt++);

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
				//	set: dim=2, id=0, offset=i
					vLocalDoF[index] = LocalDoF(3, 0, cnt++);

				//	use regular mapping for inner DoFs
					vMultiIndex[index][0] = m0;
					vMultiIndex[index][1] = m1;
					vMultiIndex[index++][2] = m2;
				}
		break;
	default: std::stringstream ss;
		ss << "LagrangeLSFS: Missing 3d mapping for type '"<<type<<"'.";
		throw(UGFatalError(ss.str().c_str()));
	}
}

template <typename TRefElem>
void SetLocalDoFs(	MultiIndex<TRefElem::dim>* vMultiIndex,
                  	LocalDoF* vLocalDoF,
                  	const TRefElem& rRef,
                  	size_t p)
{
//	init shape -> multi-index mapping
	size_t index = 0;

//	vertices
	SetVertexLocalDoFs(vMultiIndex, vLocalDoF, rRef, p, index);

//	edge
	SetEdgeLocalDoFs(vMultiIndex, vLocalDoF, rRef, p, index);

//	faces
	SetFaceLocalDoFs(vMultiIndex, vLocalDoF, rRef, p, index);

//	volumes
	SetVolumeLocalDoFs(vMultiIndex, vLocalDoF, rRef, p, index);
}

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
	SetLocalDoFs(m_vMultiIndex, m_vLocalDoF, rRef, p);
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
	SetLocalDoFs(m_vMultiIndex, m_vLocalDoF, rRef, p);
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
	SetLocalDoFs(m_vMultiIndex, m_vLocalDoF, rRef, p);
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

//	reference element
	const ReferenceTetrahedron& rRef =
			ReferenceElementProvider::get<ReferenceTetrahedron>();

//	init shape -> multi-index mapping
	SetLocalDoFs(m_vMultiIndex, m_vLocalDoF, rRef, p);
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

//	reference element
	const ReferencePrism& rRef =
			ReferenceElementProvider::get<ReferencePrism>();

//	init shape -> multi-index mapping
	SetLocalDoFs(m_vMultiIndex, m_vLocalDoF, rRef, p);
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

	//	reference element
		const ReferencePyramid& rRef =
				ReferenceElementProvider::get<ReferencePyramid>();

	//	init shape -> multi-index mapping
		SetLocalDoFs(m_vMultiIndex, m_vLocalDoF, rRef, p);
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
	SetLocalDoFs(m_vMultiIndex, m_vLocalDoF, rRef, p);
}

template class LagrangeLSFS<ReferenceHexahedron, 1>;
template class LagrangeLSFS<ReferenceHexahedron, 2>;
template class LagrangeLSFS<ReferenceHexahedron, 3>;
template class LagrangeLSFS<ReferenceHexahedron, 4>;
template class LagrangeLSFS<ReferenceHexahedron, 5>;

} // end namespace ug
