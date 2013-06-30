/*
 * lagrange_local_dof.cpp
 *
 *  Created on: 30.06.2013
 *      Author: andreasvogel
 */

#include "lagrange_local_dof.h"
#include "common/util/provider.h"
#include "lib_disc/common/multi_index.h"

namespace ug{

///////////////////////////////////////////////////////////////////////////////
// Help Functions to create LocalDoFs
///////////////////////////////////////////////////////////////////////////////

void SetLagrangeVertexLocalDoFs(std::vector<LocalDoF>& vLocalDoF,
                                const ReferenceElement& rRef,
                                const size_t p)
{
//	loop all vertices
	for(size_t i = 0; i< rRef.num(0); ++i)
		vLocalDoF.push_back(LocalDoF(0, i, 0));
}

void SetLagrangeEdgeLocalDoFs(std::vector<LocalDoF>& vLocalDoF,
                              const ReferenceElement& rRef,
                              const size_t p)
{
//	only for 2d,3d elems we do something
	if(rRef.dimension() < 1) return;

//	loop all edges
	for(size_t e = 0; e< rRef.num(1); ++e)
	{
	//	add dofs on the edge
		for(size_t i = 1; i < p; ++i)
		{
		//	set: dim=1, id=e, offset=i-1
			vLocalDoF.push_back(LocalDoF(1, e, i-1));
		}
	}
}

void SetLagrangeFaceLocalDoFs(std::vector<LocalDoF>& vLocalDoF,
                              const ReferenceElement& rRef,
                              const size_t p)
{
//	only for 2d,3d elems we do something
	if(rRef.dimension() < 2) return;

//	add dof on quadrilateral
	for(size_t f = 0; f< rRef.num(2); ++f)
	{
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
				vLocalDoF.push_back(LocalDoF(2, f, cnt++));
			}
		}
	}
}

void SetLagrangeVolumeLocalDoFs(std::vector<LocalDoF>& vLocalDoF,
                                const ReferenceElement& rRef,
                                const size_t p)
{
//	only for 3d elems we do something
	if(rRef.dimension() < 3) return;

//	get type of reference element
	const ReferenceObjectID roid = rRef.roid();

//	handle elems
	size_t cnt = 0;
	switch(roid)
	{
	case ROID_TETRAHEDRON:
		for(size_t m2 = 1; m2 < p; ++m2)
			for(size_t m1 = 1; m1 < p-m2; ++m1)
				for(size_t m0 = 1; m0 < p-m2-m1; ++m0)
				{
				//	set: dim=2, id=0, offset=i
					vLocalDoF.push_back(LocalDoF(3, 0, cnt++));
				}
		break;

	case ROID_PYRAMID:
		//\todo:order dofs
		{
			size_t numInnerDoF = 0;
			for(int i=1; i <= (int)p -2; ++i) numInnerDoF += i*i;

			for(size_t i = 0; i < numInnerDoF; ++i)
				vLocalDoF.push_back(LocalDoF(3, 0, i));
		}
		break;

	case ROID_PRISM:
		for(size_t m2 = 1; m2 < p; ++m2)
			for(size_t m1 = 1; m1 < p; ++m1)
				for(size_t m0 = 1; m0 < p-m1; ++m0)
				{
				//	set: dim=2, id=0, offset=i
					vLocalDoF.push_back(LocalDoF(3, 0, cnt++));
				}
		break;

	case ROID_HEXAHEDRON:
		for(size_t m2 = 1; m2 < p; ++m2)
			for(size_t m1 = 1; m1 < p; ++m1)
				for(size_t m0 = 1; m0 < p; ++m0)
				{
				//	set: dim=2, id=0, offset=i
					vLocalDoF.push_back(LocalDoF(3, 0, cnt++));
				}
		break;
	default: UG_THROW("SetLagrangeVolumeLocalDoFs: Missing 3d mapping "
							"for type '"<<roid<<"'.");
	}
}

void SetLagrangeLocalDoFs(	std::vector<LocalDoF>& vLocalDoF,
                          	const ReferenceElement& rRef,
                          	const size_t p)
{
	SetLagrangeVertexLocalDoFs(vLocalDoF, rRef, p);
	SetLagrangeEdgeLocalDoFs(vLocalDoF, rRef, p);
	SetLagrangeFaceLocalDoFs(vLocalDoF, rRef, p);
	SetLagrangeVolumeLocalDoFs(vLocalDoF, rRef, p);

	if(vLocalDoF.size() != LagrangeNumDoFs(rRef.roid(), p))
		UG_THROW("Wrong number of LocalDoFs ("<<vLocalDoF.size()<<") distributed, "
		         "correct is "<<LagrangeNumDoFs(rRef.roid(), p));
}

size_t GetNumberOfDoFsOfPyramid(int p)
{
	if(p <= 0) return 0;
	if(p == 0) return 1;
	if(p == 1) return 5;
	else return GetNumberOfDoFsOfPyramid(p-1) + (p+1)*(p+1);
}

size_t LagrangeNumDoFOnSub(const ReferenceObjectID elem,
                           const ReferenceObjectID sub, const size_t p)
{
	switch(elem){
		case ROID_VERTEX:
			if(sub == ROID_VERTEX) return 1;
			else return 0;
		case ROID_EDGE:
				 if(sub == ROID_VERTEX) 	return 1;
			else if(sub == ROID_EDGE) 		return p-1;
			else return 0;
		case ROID_TRIANGLE:
			if(sub == ROID_VERTEX)   return 1;
			if(sub == ROID_EDGE) 	  return (p-1);
			if(sub == ROID_TRIANGLE) return ((p>2) ? BinomCoeff(p-1, p-3) : 0);
			else return 0;
		case ROID_QUADRILATERAL:
			if(sub == ROID_VERTEX)		   return 1;
			if(sub == ROID_EDGE) 		   return (p-1);
			if(sub == ROID_QUADRILATERAL) return (p-1)*(p-1);
			else return 0;
		case ROID_TETRAHEDRON:
			if(sub == ROID_VERTEX)      return 1;
			if(sub == ROID_EDGE) 	     return (p-1);
			if(sub == ROID_TRIANGLE)    return ((p>2) ? BinomCoeff(p-1, p-3) : 0);
			if(sub == ROID_TETRAHEDRON) return ((p>3) ? BinomCoeff(p-1, p-4) : 0);
			else return 0;
		case ROID_PRISM:
			if(sub == ROID_VERTEX)        return 1;
			if(sub == ROID_EDGE)          return (p-1);
		//	same as for a 2d triangle of order p-3
			if(sub == ROID_TRIANGLE)      return ((p>2) ? BinomCoeff(p-1, p-3) : 0);
		//	same as for a 2d quadrilateral of order p-2
			if(sub == ROID_QUADRILATERAL) return (p-1)*(p-1);
		//	same as for a 3d prism of order p-2
			if(sub == ROID_PRISM)		   return ((p>2) ? BinomCoeff(p-1, p-3)*(p-1) : 0);
			else return 0;
		case ROID_PYRAMID:
			if(sub == ROID_VERTEX)			return 1;
			if(sub == ROID_EDGE) 			return (p-1);
		//	same as for a 2d triangle of order p-3
			if(sub == ROID_TRIANGLE)   	return ((p>2) ? BinomCoeff(p-1, p-3) : 0);
		//	same as for a 2d quadrilateral of order p-2
			if(sub == ROID_QUADRILATERAL)	return (p-1)*(p-1);
		//	same as for a 3d pyramid of order p-2
			if(sub == ROID_PYRAMID)		return ((p>2) ? GetNumberOfDoFsOfPyramid(p-3) : 0);
			else return 0;
		case ROID_HEXAHEDRON:
			if(sub == ROID_VERTEX)		   return 1;
			if(sub == ROID_EDGE) 	 	   return (p-1);
			if(sub == ROID_QUADRILATERAL) return (p-1)*(p-1);
			if(sub == ROID_HEXAHEDRON)    return (p-1)*(p-1)*(p-1);
			else return 0;
		default: UG_THROW("LagrangeLDS: Invalid ReferenceObjectID: "<<elem);
	}
}

size_t LagrangeNumDoFs(const ReferenceObjectID elem, const size_t p)
{
	switch(elem){
		case ROID_VERTEX: 			return 1;
		case ROID_EDGE: 			return p+1;
		case ROID_TRIANGLE: 		return BinomCoeff(2 + p, p);
		case ROID_QUADRILATERAL: 	return (p+1)*(p+1);
		case ROID_TETRAHEDRON: 		return BinomCoeff(3 + p, p);
		case ROID_PRISM: 			return BinomCoeff(2+p,p) * (p+1);
		case ROID_PYRAMID: 			return GetNumberOfDoFsOfPyramid(p);
		case ROID_HEXAHEDRON: 		return (p+1)*(p+1)*(p+1);
		default: UG_THROW("LagrangeLDS: Invalid ReferenceObjectID: "<<elem);
	}
}


///////////////////////////////////////////////////////////////////////////////
// LagrangeLDS
///////////////////////////////////////////////////////////////////////////////

template <typename TRefElem>
LagrangeLDS<TRefElem>::LagrangeLDS(size_t order)
{
	set_order(order);
}

template <typename TRefElem>
void LagrangeLDS<TRefElem>::set_order(size_t order)
{
	p = order;
	m_vLocalDoF.clear();
	SetLagrangeLocalDoFs(m_vLocalDoF, Provider<TRefElem>::get(), p);
}

template <typename TRefElem>
size_t LagrangeLDS<TRefElem>::num_dof(ReferenceObjectID type) const
{
	return LagrangeNumDoFOnSub(roid(), type, p);
}

template class LagrangeLDS<ReferenceVertex>;
template class LagrangeLDS<ReferenceEdge>;
template class LagrangeLDS<ReferenceTriangle>;
template class LagrangeLDS<ReferenceQuadrilateral>;
template class LagrangeLDS<ReferenceTetrahedron>;
template class LagrangeLDS<ReferencePrism>;
template class LagrangeLDS<ReferencePyramid>;
template class LagrangeLDS<ReferenceHexahedron>;

} // end namespace ug
