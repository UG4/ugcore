/*
 * lagrange_local_dof.h
 *
 *  Created on: 27.06.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__LAGRANGE__LAGRANGE_LOCAL_DOF__
#define __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__LAGRANGE__LAGRANGE_LOCAL_DOF__

#include "common/util/provider.h"
#include "../local_dof_set.h"
#include "lib_disc/common/multi_index.h"

namespace ug{

///////////////////////////////////////////////////////////////////////////////
// Help Functions to create LocalDoFs
///////////////////////////////////////////////////////////////////////////////

template <typename TRefElem>
void SetLagrangeVertexLocalDoFs(std::vector<LocalDoF>& vLocalDoF,
                                const TRefElem& rRef,
                                size_t p,
                                size_t& index)
{
//	loop all vertices
	for(size_t i = 0; i< rRef.num(0); ++i)
		vLocalDoF[index++] = LocalDoF(0, i, 0);
}

template <typename TRefElem>
void SetLagrangeEdgeLocalDoFs(std::vector<LocalDoF>& vLocalDoF,
                              const TRefElem& rRef,
                              size_t p,
                              size_t& index)
{
//	only for 2d,3d elems we do something
	if(TRefElem::dim < 1) return;

//	loop all edges
	for(size_t e = 0; e< rRef.num(1); ++e)
	{
	//	add dofs on the edge
		for(size_t i = 1; i < p; ++i)
		{
		//	set: dim=1, id=e, offset=i-1
			vLocalDoF[index++] = LocalDoF(1, e, i-1);
		}
	}
}

template <typename TRefElem>
void SetLagrangeFaceLocalDoFs(std::vector<LocalDoF>& vLocalDoF,
                              const TRefElem& rRef,
                              size_t p,
                              size_t& index)
{
//	only for 2d,3d elems we do something
	if(TRefElem::dim < 2) return;

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
				vLocalDoF[index++] = LocalDoF(2, f, cnt++);
			}
		}
	}
}

template <typename TRefElem>
void SetLagrangeVolumeLocalDoFs(std::vector<LocalDoF>& vLocalDoF,
                                const TRefElem& rRef,
                                size_t p,
                                size_t& index)
{
//	only for 3d elems we do something
	if(TRefElem::dim < 3) return;

//	get type of reference element
	ReferenceObjectID type = TRefElem::REFERENCE_OBJECT_ID;

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
					vLocalDoF[index++] = LocalDoF(3, 0, cnt++);
				}
		break;

	case ROID_PYRAMID:
		//\todo:order dofs
		{
			size_t numInnerDoF = 0;
			for(int i=1; i <= (int)p -2; ++i) numInnerDoF += i*i;

			for(size_t i = 0; i < numInnerDoF; ++i)
				vLocalDoF[index++] = LocalDoF(3, 0, i);
		}
		break;

	case ROID_PRISM:
		for(size_t m2 = 1; m2 < p; ++m2)
			for(size_t m1 = 1; m1 < p; ++m1)
				for(size_t m0 = 1; m0 < p-m1; ++m0)
				{
				//	set: dim=2, id=0, offset=i
					vLocalDoF[index++] = LocalDoF(3, 0, cnt++);
				}
		break;

	case ROID_HEXAHEDRON:
		for(size_t m2 = 1; m2 < p; ++m2)
			for(size_t m1 = 1; m1 < p; ++m1)
				for(size_t m0 = 1; m0 < p; ++m0)
				{
				//	set: dim=2, id=0, offset=i
					vLocalDoF[index++] = LocalDoF(3, 0, cnt++);
				}
		break;
	default: UG_THROW("SetLagrangeVolumeLocalDoFs: Missing 3d mapping "
							"for type '"<<type<<"'.");
	}
}

template <typename TRefElem>
void SetLagrangeLocalDoFs(	std::vector<LocalDoF>& vLocalDoF,
                          	const TRefElem& rRef,
                          	size_t p)
{
//	init shape -> multi-index mapping
	size_t index = 0;

//	vertices
	SetLagrangeVertexLocalDoFs(vLocalDoF, rRef, p, index);

//	edge
	SetLagrangeEdgeLocalDoFs(vLocalDoF, rRef, p, index);

//	faces
	SetLagrangeFaceLocalDoFs(vLocalDoF, rRef, p, index);

//	volumes
	SetLagrangeVolumeLocalDoFs(vLocalDoF, rRef, p, index);

	UG_ASSERT(index == vLocalDoF.size(), "Wrong number of LocalDoFs ("<<index<<
	         					") distributed, correct is "<<vLocalDoF.size());
}

///////////////////////////////////////////////////////////////////////////////
// Lagrange LocalDoFSets
///////////////////////////////////////////////////////////////////////////////

/// Lagrange DoF Set
template <typename TRefElem>
class LagrangeLDS{};


///////////////////////////////////////////////////////////////////////////////
// Edge
///////////////////////////////////////////////////////////////////////////////

/// specialization for Edges
template <>
class LagrangeLDS<ReferenceVertex>
	: public LocalDoFSet
{
	public:
	///	constructor
		LagrangeLDS() {set_order(1);}

	///	constructor
		LagrangeLDS(size_t order) {set_order(order);}

	///	sets the order
		void set_order(size_t order)
		{
			p = order;
			nsh = 1;
			m_vLocalDoF.resize(nsh);
			SetLagrangeLocalDoFs(m_vLocalDoF, Provider<ReferenceVertex>::get(), p);
		}

	///	returns the reference dimension
		int dim() const {return ReferenceVertex::dim;}

	///	returns the type of reference element
		ReferenceObjectID roid() const {return ReferenceVertex::REFERENCE_OBJECT_ID;}

	///	returns the total number of DoFs on the finite element
		size_t num_dof() const {return nsh;};

	///	returns the number of DoFs on a sub-geometric object type
		int num_dof(ReferenceObjectID type) const
		{
			if(type == ROID_VERTEX) return 1;
			else return -1;
		}

	///	returns the number of DoFs on sub-geometric object in dimension and id
		size_t num_dof(int d, size_t id) const
		{
			return num_dof(Provider<ReferenceVertex>::get().roid(d, id));
		}

	///	returns if the storage needs objects of a given dimension
		size_t max_num_dof(int d) const
		{
			if(d == 0) return 1;
			else return 0;
		}

	///	returns the dof storage
		const LocalDoF& local_dof(size_t dof) const {return m_vLocalDoF[dof];}

	protected:
	///	number of shapes
		size_t nsh;

	///	order
		size_t p;

	///	association to elements
		std::vector<LocalDoF> m_vLocalDoF;
};

///////////////////////////////////////////////////////////////////////////////
// Edge
///////////////////////////////////////////////////////////////////////////////

/// specialization for Edges
template <>
class LagrangeLDS<ReferenceEdge>
: public LocalDoFSet
{
	public:
	///	constructor
		LagrangeLDS() {set_order(1);}

	///	constructor
		LagrangeLDS(size_t order) {set_order(order);}

	///	sets the order
		void set_order(size_t order)
		{
			p = order;
			nsh = p+1;
			m_vLocalDoF.resize(nsh);
			SetLagrangeLocalDoFs(m_vLocalDoF, Provider<ReferenceEdge>::get(), p);
		}

	///	returns the reference dimension
		int dim() const {return ReferenceEdge::dim;}

	///	returns the type of reference element
		ReferenceObjectID roid() const {return ReferenceEdge::REFERENCE_OBJECT_ID;}

	///	returns the total number of DoFs on the finite element
		size_t num_dof() const {return nsh;};

	///	returns the number of DoFs on a sub-geometric object type
		int num_dof(ReferenceObjectID type) const
		{
				 if(type == ROID_VERTEX) return 1;
			else if(type == ROID_EDGE) return p-1;
			else return -1;
		}

	///	returns the number of DoFs on sub-geometric object in dimension and id
		size_t num_dof(int d, size_t id) const
		{
			return num_dof(Provider<ReferenceEdge>::get().roid(d, id));
		}

	///	returns if the storage needs objects of a given dimension
		size_t max_num_dof(int d) const
		{
				 if(d == 0) return 1;
			else if(d == 1) return p-1;
			else return 0;
		}

	///	returns the dof storage
		const LocalDoF& local_dof(size_t dof) const {return m_vLocalDoF[dof];}

	protected:
	///	number of shapes
		size_t nsh;

	///	order
		size_t p;

	///	association to elements
		std::vector<LocalDoF> m_vLocalDoF;
};


///////////////////////////////////////////////////////////////////////////////
// Triangle
///////////////////////////////////////////////////////////////////////////////

/// specialization for Triangles
template <>
class LagrangeLDS<ReferenceTriangle>
: public LocalDoFSet
{
	public:
	///	constructor
		LagrangeLDS() {set_order(1);}

	///	constructor
		LagrangeLDS(size_t order) {set_order(order);}

	///	sets the order
		void set_order(size_t order)
		{
			p = order;
			nsh = BinomCoeff(2 + p, p);
			m_vLocalDoF.resize(nsh);
			SetLagrangeLocalDoFs(m_vLocalDoF, Provider<ReferenceTriangle>::get(), p);
		}

	///	returns the reference dimension
		int dim() const {return ReferenceTriangle::dim;}

	///	returns the type of reference element
		ReferenceObjectID roid() const {return ReferenceTriangle::REFERENCE_OBJECT_ID;}

	///	returns the total number of DoFs on the finite element
		size_t num_dof() const {return nsh;};

	///	returns the number of DoFs on a sub-geometric object type
		int num_dof(ReferenceObjectID type) const
		{
			if(type == ROID_VERTEX)   return 1;
			if(type == ROID_EDGE) 	  return (p-1);
			if(type == ROID_TRIANGLE) return ((p>2) ? BinomCoeff(p-1, p-3) : 0);
			else return -1;
		}

	///	returns the number of DoFs on sub-geometric object in dimension and id
		size_t num_dof(int d, size_t id) const
		{
			return num_dof(Provider<ReferenceTriangle>::get().roid(d, id));
		}

	///	returns if the storage needs objects of a given dimension
		size_t max_num_dof(int d) const
		{
			if(d==0) return 1;
			if(d==1) return (p-1);
			if(d==2) return ((p>2) ? BinomCoeff(p-1, p-3) : 0);
			else return 0;
		}

	///	returns the dof storage
		const LocalDoF& local_dof(size_t dof) const {return m_vLocalDoF[dof];}

	protected:
	///	number of shapes
		size_t nsh;

	///	order
		size_t p;

	///	association to elements
		std::vector<LocalDoF> m_vLocalDoF;
};


///////////////////////////////////////////////////////////////////////////////
// Quadrilateral
///////////////////////////////////////////////////////////////////////////////

/// specialization for Quadrilateral
template <>
class LagrangeLDS<ReferenceQuadrilateral>
: public LocalDoFSet
{
	public:
	///	constructor
		LagrangeLDS() {set_order(1);}

	///	constructor
		LagrangeLDS(size_t order) {set_order(order);}

	///	sets the order
		void set_order(size_t order)
		{
			p = order;
			nsh = (p+1)*(p+1);
			m_vLocalDoF.resize(nsh);
			SetLagrangeLocalDoFs(m_vLocalDoF, Provider<ReferenceQuadrilateral>::get(), p);
		}

	///	returns the reference dimension
		int dim() const {return ReferenceQuadrilateral::dim;}

	///	returns the type of reference element
		ReferenceObjectID roid() const {return ReferenceQuadrilateral::REFERENCE_OBJECT_ID;}

	///	returns the total number of DoFs on the finite element
		size_t num_dof() const {return nsh;};

	///	returns the number of DoFs on a sub-geometric object type
		int num_dof(ReferenceObjectID type) const
		{
			if(type == ROID_VERTEX)		   return 1;
			if(type == ROID_EDGE) 		   return (p-1);
			if(type == ROID_QUADRILATERAL) return (p-1)*(p-1);
			else return -1;
		}

	///	returns if the storage needs objects of a given dimension
		size_t max_num_dof(int d) const
		{
			if(d==0) return 1;
			if(d==1) return (p-1);
			if(d==2) return (p-1)*(p-1);
			else return 0;
		}

	///	returns the number of DoFs on sub-geometric object in dimension and id
		size_t num_dof(int d, size_t id) const
		{
			return num_dof(Provider<ReferenceQuadrilateral>::get().roid(d, id));
		}

	///	returns the dof storage
		const LocalDoF& local_dof(size_t dof) const {return m_vLocalDoF[dof];}

	protected:
	///	number of shapes
		size_t nsh;

	///	order
		size_t p;

	///	association to elements
		std::vector<LocalDoF> m_vLocalDoF;
};

///////////////////////////////////////////////////////////////////////////////
// Tetrahedron
///////////////////////////////////////////////////////////////////////////////

/// specialization for Tetrahedron
template <>
class LagrangeLDS<ReferenceTetrahedron>
: public LocalDoFSet
{
	public:
	///	constructor
		LagrangeLDS() {set_order(1);}

	///	constructor
		LagrangeLDS(size_t order) {set_order(order);}

	///	sets the order
		void set_order(size_t order)
		{
			p = order;
			nsh = BinomCoeff(3 + p, p);
			m_vLocalDoF.resize(nsh);
			SetLagrangeLocalDoFs(m_vLocalDoF, Provider<ReferenceTetrahedron>::get(), p);
		}

	///	returns the reference dimension
		int dim() const {return ReferenceTetrahedron::dim;}

	///	returns the type of reference element
		ReferenceObjectID roid() const {return ReferenceTetrahedron::REFERENCE_OBJECT_ID;}

	///	returns the total number of DoFs on the finite element
		size_t num_dof() const {return nsh;};

	///	returns the number of DoFs on a sub-geometric object type
		int num_dof(ReferenceObjectID type) const
		{
			if(type == ROID_VERTEX)      return 1;
			if(type == ROID_EDGE) 	     return (p-1);
			if(type == ROID_TRIANGLE)    return ((p>2) ? BinomCoeff(p-1, p-3) : 0);
			if(type == ROID_TETRAHEDRON) return ((p>3) ? BinomCoeff(p-1, p-4) : 0);
			else return -1;
		}

	///	returns the number of DoFs on sub-geometric object in dimension and id
		size_t num_dof(int d, size_t id) const
		{
			return num_dof(Provider<ReferenceTetrahedron>::get().roid(d, id));
		}

	///	returns if the storage needs objects of a given dimension
		size_t max_num_dof(int d) const
		{
			if(d==0) return 1;
		//	number of shapes on edge is same as for edge with p-1
			if(d==1) return (p-1);
		//	number of shapes on faces is same as for triangles in 2d
			if(d==2) return ((p>2) ? BinomCoeff(p-1, p-3) : 0);
		//	number of shapes on interior is same as for tetrahedra with p-4
			if(d==3) return ((p>3) ? BinomCoeff(p-1, p-4) : 0);
			else return 0;
		}

	///	returns the dof storage
		const LocalDoF& local_dof(size_t dof) const {return m_vLocalDoF[dof];}

	protected:
	///	number of shapes
		size_t nsh;

	///	order
		size_t p;

	///	association to elements
		std::vector<LocalDoF> m_vLocalDoF;
};

///////////////////////////////////////////////////////////////////////////////
// Prism
///////////////////////////////////////////////////////////////////////////////

/// specialization for Prism
template <>
class LagrangeLDS<ReferencePrism>
: public LocalDoFSet
{
	public:
	///	constructor
		LagrangeLDS() {set_order(1);}

	///	constructor
		LagrangeLDS(size_t order) {set_order(order);}

	///	sets the order
		void set_order(size_t order)
		{
			p = order;
			nsh = BinomCoeff(2+p,p) * (p+1);
			m_vLocalDoF.resize(nsh);
			SetLagrangeLocalDoFs(m_vLocalDoF, Provider<ReferencePrism>::get(), p);
		}

	///	returns the reference dimension
		int dim() const {return ReferencePrism::dim;}

	///	returns the type of reference element
		ReferenceObjectID roid() const {return ReferencePrism::REFERENCE_OBJECT_ID;}

	///	returns the total number of DoFs on the finite element
		size_t num_dof() const {return nsh;};

	///	returns the number of DoFs on a sub-geometric object type
		int num_dof(ReferenceObjectID type) const
		{
			if(type == ROID_VERTEX)        return 1;
			if(type == ROID_EDGE)          return (p-1);
		//	same as for a 2d triangle of order p-3
			if(type == ROID_TRIANGLE)      return ((p>2) ? BinomCoeff(p-1, p-3) : 0);
		//	same as for a 2d quadrilateral of order p-2
			if(type == ROID_QUADRILATERAL) return (p-1)*(p-1);
		//	same as for a 3d prism of order p-2
			if(type == ROID_PRISM)		   return ((p>2) ? BinomCoeff(p-1, p-3)*(p-1) : 0);
			else return -1;
		}

	///	returns the number of DoFs on sub-geometric object in dimension and id
		size_t num_dof(int d, size_t id) const
		{
			return num_dof(Provider<ReferencePrism>::get().roid(d, id));
		}

	///	returns if the storage needs objects of a given dimension
		size_t max_num_dof(int d) const
		{
			if(d==0) return 1;
			if(d==1) return (p-1);
			if(d==2) return num_dof(ROID_QUADRILATERAL);
			if(d==3) return num_dof(ROID_PRISM);
			else return 0;
		}

	///	returns the dof storage
		const LocalDoF& local_dof(size_t dof) const {return m_vLocalDoF[dof];}

	protected:
	///	number of shapes
		size_t nsh;

	///	order
		size_t p;

	///	association to elements
		std::vector<LocalDoF> m_vLocalDoF;
};

///////////////////////////////////////////////////////////////////////////////
// Pyramid
///////////////////////////////////////////////////////////////////////////////

size_t GetNumberOfDoFsOfPyramid(int p)
{
	if(p <= 0) return 0;
	if(p == 0) return 1;
	if(p == 1) return 5;
	else return GetNumberOfDoFsOfPyramid(p-1) + (p+1)*(p+1);
}


/// specialization for Pyramid
template <>
class LagrangeLDS<ReferencePyramid>
: public LocalDoFSet
{
	public:
	///	constructor
		LagrangeLDS() {set_order(1);}

	///	constructor
		LagrangeLDS(size_t order) {set_order(order);}

	///	sets the order
		void set_order(size_t order)
		{
			p = order;
			nsh = GetNumberOfDoFsOfPyramid(p);
			m_vLocalDoF.resize(nsh);
			SetLagrangeLocalDoFs(m_vLocalDoF, Provider<ReferencePyramid>::get(), p);
		}

	///	returns the reference dimension
		int dim() const {return ReferencePyramid::dim;}

	///	returns the type of reference element
		ReferenceObjectID roid() const {return ReferencePyramid::REFERENCE_OBJECT_ID;}

	///	returns the total number of DoFs on the finite element
		size_t num_dof() const {return nsh;};

	///	returns the number of DoFs on a sub-geometric object type
		int num_dof(ReferenceObjectID type) const
		{
			if(type == ROID_VERTEX)			return 1;
			if(type == ROID_EDGE) 			return (p-1);
		//	same as for a 2d triangle of order p-3
			if(type == ROID_TRIANGLE)   	return ((p>2) ? BinomCoeff(p-1, p-3) : 0);
		//	same as for a 2d quadrilateral of order p-2
			if(type == ROID_QUADRILATERAL)	return (p-1)*(p-1);
		//	same as for a 3d pyramid of order p-2
			if(type == ROID_PYRAMID)		return ((p>2) ? GetNumberOfDoFsOfPyramid(p-3) : 0);
			else return -1;
		}

	///	returns the number of DoFs on sub-geometric object in dimension and id
		size_t num_dof(int d, size_t id) const
		{
			return num_dof(Provider<ReferencePyramid>::get().roid(d, id));
		}

	///	returns if the storage needs objects of a given dimension
		size_t max_num_dof(int d) const
		{
			if(d==0) return 1;
			if(d==1) return (p-1);
			if(d==2) return num_dof(ROID_QUADRILATERAL);
			if(d==3) return num_dof(ROID_PYRAMID);
			else return 0;
		}

	///	returns the dof storage
		const LocalDoF& local_dof(size_t dof) const {return m_vLocalDoF[dof];}

	protected:
	///	number of shapes
		size_t nsh;

	///	order
		size_t p;

	///	association to elements
		std::vector<LocalDoF> m_vLocalDoF;
};


///////////////////////////////////////////////////////////////////////////////
// Hexahedron
///////////////////////////////////////////////////////////////////////////////

/// specialization for Hexahedron
template <>
class LagrangeLDS<ReferenceHexahedron>
: public LocalDoFSet
{
	public:
	///	constructor
		LagrangeLDS() {set_order(1);}

	///	constructor
		LagrangeLDS(size_t order) {set_order(order);}

	///	sets the order
		void set_order(size_t order)
		{
			p = order;
			nsh = (p+1)*(p+1)*(p+1);
			m_vLocalDoF.resize(nsh);
			SetLagrangeLocalDoFs(m_vLocalDoF, Provider<ReferenceHexahedron>::get(), p);
		}

	///	returns the reference dimension
		int dim() const {return ReferenceHexahedron::dim;}

	///	returns the type of reference element
		ReferenceObjectID roid() const {return ReferenceHexahedron::REFERENCE_OBJECT_ID;}

	///	returns the total number of DoFs on the finite element
		size_t num_dof() const {return nsh;};

	///	returns the number of DoFs on a sub-geometric object type
		int num_dof(ReferenceObjectID type) const
		{
			if(type == ROID_VERTEX)		   return 1;
			if(type == ROID_EDGE) 	 	   return (p-1);
			if(type == ROID_QUADRILATERAL) return (p-1)*(p-1);
			if(type == ROID_HEXAHEDRON)    return (p-1)*(p-1)*(p-1);
			else return -1;
		}

	///	returns the number of DoFs on sub-geometric object in dimension and id
		size_t num_dof(int d, size_t id) const
		{
			return num_dof(Provider<ReferenceHexahedron>::get().roid(d, id));
		}

	///	returns if the storage needs objects of a given dimension
		size_t max_num_dof(int d) const
		{
			if(d==0) return 1;
			if(d==1) return (p-1);
			if(d==2) return (p-1)*(p-1);
			if(d==3) return (p-1)*(p-1)*(p-1);
			else return 0;
		}

	///	returns the dof storage
		const LocalDoF& local_dof(size_t dof) const {return m_vLocalDoF[dof];}

	protected:
	///	number of shapes
		size_t nsh;

	///	order
		size_t p;

	///	association to elements
		std::vector<LocalDoF> m_vLocalDoF;
};

} //namespace ug

#endif /* __H__UG__LIB_DISC__LOCAL_SHAPE_FUNCTION_SET__LAGRANGE__LAGRANGE_LOCAL_DOF__ */
