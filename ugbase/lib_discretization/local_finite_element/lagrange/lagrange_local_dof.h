/*
 * lagrange_local_dof.h
 *
 *  Created on: 27.06.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET__LAGRANGE__LAGRANGE_LOCAL_DOF__
#define __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET__LAGRANGE__LAGRANGE_LOCAL_DOF__

#include "common/util/provider.h"
#include "lagrange.h"
#include "../local_dof_set.h"
#include "lib_discretization/common/multi_index.h"

namespace ug{

///////////////////////////////////////////////////////////////////////////////
// Help Functions to create LocalDoFs
///////////////////////////////////////////////////////////////////////////////

template <typename TRefElem>
void SetLagrangeVertexLocalDoFs(LocalDoF* vLocalDoF,
                                const TRefElem& rRef,
                                size_t p,
                                size_t& index)
{
//	loop all vertices
	for(size_t i = 0; i< rRef.num(0); ++i)
		vLocalDoF[index++] = LocalDoF(0, i, 0);
}

template <typename TRefElem>
void SetLagrangeEdgeLocalDoFs(LocalDoF* vLocalDoF,
                              const TRefElem& rRef,
                              size_t p,
                              size_t& index)
{
//	dimension of Reference element
	static const int dim = TRefElem::dim;

//	only for 2d,3d elems we do something
	if(dim < 1) return;

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
void SetLagrangeFaceLocalDoFs(LocalDoF* vLocalDoF,
                              const TRefElem& rRef,
                              size_t p,
                              size_t& index)
{
//	dimension of Reference element
	static const int dim = TRefElem::dim;

//	only for 2d,3d elems we do something
	if(dim < 2) return;

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
void SetLagrangeVolumeLocalDoFs(LocalDoF* vLocalDoF,
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
	default: std::stringstream ss;
		ss << "SetLagrangeVolumeLocalDoFs: Missing 3d mapping for type '"<<type<<"'.";
		throw(UGFatalError(ss.str().c_str()));
	}
}

template <typename TRefElem>
void SetLagrangeLocalDoFs(	LocalDoF* vLocalDoF,
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
}

///////////////////////////////////////////////////////////////////////////////
// Lagrange LocalDoFSets
///////////////////////////////////////////////////////////////////////////////

/// Lagrange DoF Set
template <typename TRefElem, int TOrder>
struct LagrangeLDS{};


///////////////////////////////////////////////////////////////////////////////
// Edge
///////////////////////////////////////////////////////////////////////////////

/// specialization for Edges
template <>
template <int TOrder>
class LagrangeLDS<ReferenceVertex, TOrder>
{
	public:
	///	number of shapes
		static const size_t nsh = 1;

	///	dimension of reference element
		const static int Dim = 0;

	///	order
		static const size_t order = TOrder;

	public:
	///	constructor
		LagrangeLDS() : m_rRef(Provider::get<ReferenceVertex>())
		{
		//	create LocalDoF vector
			SetLagrangeLocalDoFs(m_vLocalDoF, m_rRef, order);
		}

	///	returns the reference dimension
		static int dim() {return Dim;}

	///	returns the type of reference element
		static ReferenceObjectID roid() {return ReferenceVertex::REFERENCE_OBJECT_ID;}

	///	returns the total number of DoFs on the finite element
		static size_t num_dof() {return nsh;};

	///	returns the number of DoFs on a sub-geometric object type
		static int num_dof(ReferenceObjectID type)
		{
			if(type == ROID_VERTEX) return 1;
			else return -1;
		}

	///	returns the number of DoFs on sub-geometric object in dimension and id
		size_t num_dof(int d, size_t id) const
		{
			return num_dof(m_rRef.ref_elem_type(d, id));
		}

	///	returns if the storage needs objects of a given dimension
		static size_t max_num_dof(int d)
		{
			if(d == 1) return 1;
			else return 0;
		}

	///	returns the dof storage
		const LocalDoF& local_dof(size_t dof) const {return m_vLocalDoF[dof];}

	protected:
	///	association to elements
		LocalDoF m_vLocalDoF[nsh];

	//	reference element
		const ReferenceVertex& m_rRef;
};

///////////////////////////////////////////////////////////////////////////////
// Edge
///////////////////////////////////////////////////////////////////////////////

/// specialization for Edges
template <>
template <int TOrder>
class LagrangeLDS<ReferenceEdge, TOrder>
{
	protected:
	///	corresponding local shape function set
		typedef LagrangeLSFS<ReferenceEdge, TOrder> LSFS;

	///	abbreviation for order
		static const size_t p = TOrder;

	public:
	///	number of shapes
		static const size_t nsh = LSFS::nsh;

	///	dimension of reference element
		const static int Dim = LSFS::dim;

	///	order
		static const size_t order = TOrder;

	public:
	///	constructor
		LagrangeLDS() : m_rRef(Provider::get<ReferenceEdge>())
		{
		//	create LocalDoF vector
			SetLagrangeLocalDoFs(m_vLocalDoF, m_rRef, p);
		}

	///	returns the reference dimension
		static int dim() {return Dim;}

	///	returns the type of reference element
		static ReferenceObjectID roid() {return ReferenceEdge::REFERENCE_OBJECT_ID;}

	///	returns the total number of DoFs on the finite element
		static size_t num_dof() {return nsh;};

	///	returns the number of DoFs on a sub-geometric object type
		static int num_dof(ReferenceObjectID type)
		{
				 if(type == ROID_VERTEX) return 1;
			else if(type == ROID_EDGE) return p-1;
			else return -1;
		}

	///	returns the number of DoFs on sub-geometric object in dimension and id
		size_t num_dof(int d, size_t id) const
		{
			return num_dof(m_rRef.ref_elem_type(d, id));
		}

	///	returns if the storage needs objects of a given dimension
		static size_t max_num_dof(int d)
		{
				 if(d == 1) return 1;
			else if(d == 2) return p-1;
			else return 0;
		}

	///	returns the dof storage
		const LocalDoF& local_dof(size_t dof) const {return m_vLocalDoF[dof];}

	protected:
	///	association to elements
		LocalDoF m_vLocalDoF[nsh];

	//	reference element
		const ReferenceEdge& m_rRef;
};


///////////////////////////////////////////////////////////////////////////////
// Triangle
///////////////////////////////////////////////////////////////////////////////

/// specialization for Triangles
template <>
template <int TOrder>
class LagrangeLDS<ReferenceTriangle, TOrder>
{
	protected:
	///	corresponding local shape function set
		typedef LagrangeLSFS<ReferenceTriangle, TOrder> LSFS;

	///	abbreviation for order
		static const size_t p = TOrder;

	public:
	///	number of shapes
		static const size_t nsh = LSFS::nsh;

	///	dimension of reference element
		const static int Dim = LSFS::dim;

	///	order
		static const size_t order = TOrder;

	public:
	///	constructor
		LagrangeLDS() : m_rRef(Provider::get<ReferenceTriangle>())
		{
		//	create LocalDoF vector
			SetLagrangeLocalDoFs(m_vLocalDoF, m_rRef, p);
		}

	///	returns the reference dimension
		static int dim() {return Dim;}

	///	returns the type of reference element
		static ReferenceObjectID roid() {return ReferenceTriangle::REFERENCE_OBJECT_ID;}

	///	returns the total number of DoFs on the finite element
		static size_t num_dof() {return nsh;};

	///	returns the number of DoFs on a sub-geometric object type
		static int num_dof(ReferenceObjectID type)
		{
			if(type == ROID_VERTEX) return 1;
			if(type == ROID_EDGE) return (p-1);
			if(type == ROID_TRIANGLE)
				return ((p>2) ? (BinomialCoefficient<Dim + p-3, p-3>::value) : 0);
			else return -1;
		}

	///	returns the number of DoFs on sub-geometric object in dimension and id
		size_t num_dof(int d, size_t id) const
		{
			return num_dof(m_rRef.ref_elem_type(d, id));
		}

	///	returns if the storage needs objects of a given dimension
		static size_t max_num_dof(int d)
		{
			if(d==0) return 1;
			if(d==1) return (p-1);
			if(d==2) return ((p>2) ? (BinomialCoefficient<Dim + p-3, p-3>::value) : 0);
			else return 0;
		}

	///	returns the dof storage
		const LocalDoF& local_dof(size_t dof) const {return m_vLocalDoF[dof];}

	protected:
	///	association to elements
		LocalDoF m_vLocalDoF[nsh];

	///	reference element
		const ReferenceTriangle& m_rRef;
};


///////////////////////////////////////////////////////////////////////////////
// Quadrilateral
///////////////////////////////////////////////////////////////////////////////

/// specialization for Quadrilateral
template <>
template <int TOrder>
class LagrangeLDS<ReferenceQuadrilateral, TOrder>
{
	protected:
	///	corresponding local shape function set
		typedef LagrangeLSFS<ReferenceQuadrilateral, TOrder> LSFS;

	///	abbreviation for order
		static const size_t p = TOrder;

	public:
	///	number of shapes
		static const size_t nsh = LSFS::nsh;

	///	dimension of reference element
		const static int Dim = LSFS::dim;

	///	order
		static const size_t order = TOrder;

	public:
	///	constructor
		LagrangeLDS() : m_rRef(Provider::get<ReferenceQuadrilateral>())
		{
		//	create LocalDoF vector
			SetLagrangeLocalDoFs(m_vLocalDoF, m_rRef, p);
		}

	///	returns the reference dimension
		static int dim() {return Dim;}

	///	returns the type of reference element
		static ReferenceObjectID roid() {return ReferenceQuadrilateral::REFERENCE_OBJECT_ID;}

	///	returns the total number of DoFs on the finite element
		static size_t num_dof() {return nsh;};

	///	returns the number of DoFs on a sub-geometric object type
		static int num_dof(ReferenceObjectID type)
		{
			if(type == ROID_VERTEX) return 1;
			if(type == ROID_EDGE) return (p-1);
			if(type == ROID_QUADRILATERAL) return (p-1)*(p-1);
			else return -1;
		}

	///	returns if the storage needs objects of a given dimension
		static size_t max_num_dof(int d)
		{
			if(d==0) return 1;
			if(d==1) return (p-1);
			if(d==2) return (p-1)*(p-1);
			else return 0;
		}

	///	returns the number of DoFs on sub-geometric object in dimension and id
		size_t num_dof(int d, size_t id) const
		{
			return num_dof(m_rRef.ref_elem_type(d, id));
		}

	///	returns the dof storage
		const LocalDoF& local_dof(size_t dof) const {return m_vLocalDoF[dof];}

	protected:
	///	association to elements
		LocalDoF m_vLocalDoF[nsh];

	//	reference element
		const ReferenceQuadrilateral& m_rRef;
};

///////////////////////////////////////////////////////////////////////////////
// Tetrahedron
///////////////////////////////////////////////////////////////////////////////

/// specialization for Tetrahedron
template <>
template <int TOrder>
class LagrangeLDS<ReferenceTetrahedron, TOrder>
{
	protected:
	///	corresponding local shape function set
		typedef LagrangeLSFS<ReferenceTetrahedron, TOrder> LSFS;

	///	abbreviation for order
		static const size_t p = TOrder;

	public:
	///	number of shapes
		static const size_t nsh = LSFS::nsh;

	///	dimension of reference element
		const static int Dim = LSFS::dim;

	///	order
		static const size_t order = TOrder;

	public:
	///	constructor
		LagrangeLDS() : m_rRef(Provider::get<ReferenceTetrahedron>())
		{
		//	create LocalDoF vector
			SetLagrangeLocalDoFs(m_vLocalDoF, m_rRef, p);
		}

	///	returns the reference dimension
		static int dim() {return Dim;}

	///	returns the type of reference element
		static ReferenceObjectID roid() {return ReferenceTetrahedron::REFERENCE_OBJECT_ID;}

	///	returns the total number of DoFs on the finite element
		static size_t num_dof() {return nsh;};

	///	returns the number of DoFs on a sub-geometric object type
		static int num_dof(ReferenceObjectID type)
		{
			if(type == ROID_VERTEX) return 1;
			if(type == ROID_EDGE) return (p-1);
			if(type == ROID_TRIANGLE)
				return ((p>2) ? (BinomialCoefficient<Dim-1 + p-3, p-3>::value) : 0);
			if(type == ROID_TETRAHEDRON)
				return ((p>3) ? (BinomialCoefficient<Dim + p-4, p-4>::value) : 0);
			else return -1;
		}

	///	returns the number of DoFs on sub-geometric object in dimension and id
		size_t num_dof(int d, size_t id) const
		{
			return num_dof(m_rRef.ref_elem_type(d, id));
		}

	///	returns if the storage needs objects of a given dimension
		static size_t max_num_dof(int d)
		{
			if(d==0) return 1;
		//	number of shapes on edge is same as for edge with p-1
			if(d==1) return (p-1);
		//	number of shapes on faces is same as for triangles in 2d
			if(d==2) return ((p>2) ? (BinomialCoefficient<Dim-1 + p-3, p-3>::value) : 0);
		//	number of shapes on interior is same as for tetrahedra with p-4
			if(d==3) return ((p>3) ? (BinomialCoefficient<Dim + p-4, p-4>::value) : 0);
			else return 0;
		}

	///	returns the dof storage
		const LocalDoF& local_dof(size_t dof) const {return m_vLocalDoF[dof];}

	protected:
	///	association to elements
		LocalDoF m_vLocalDoF[nsh];

	///	reference element
		const ReferenceTetrahedron& m_rRef;
};

///////////////////////////////////////////////////////////////////////////////
// Prism
///////////////////////////////////////////////////////////////////////////////

/// specialization for Prism
template <>
template <int TOrder>
class LagrangeLDS<ReferencePrism, TOrder>
{
	protected:
	///	corresponding local shape function set
		typedef LagrangeLSFS<ReferencePrism, TOrder> LSFS;

	///	abbreviation for order
		static const size_t p = TOrder;

	public:
	///	number of shapes
		static const size_t nsh = LSFS::nsh;

	///	dimension of reference element
		const static int Dim = LSFS::dim;

	///	order
		static const size_t order = TOrder;

	public:
	///	constructor
		LagrangeLDS() : m_rRef(Provider::get<ReferencePrism>())
		{
		//	create LocalDoF vector
			SetLagrangeLocalDoFs(m_vLocalDoF, m_rRef, p);
		}

	///	returns the reference dimension
		static int dim() {return Dim;}

	///	returns the type of reference element
		static ReferenceObjectID roid() {return ReferencePrism::REFERENCE_OBJECT_ID;}

	///	returns the total number of DoFs on the finite element
		static size_t num_dof() {return nsh;};

	///	returns the number of DoFs on a sub-geometric object type
		static int num_dof(ReferenceObjectID type)
		{
			if(type == ROID_VERTEX) return 1;
			if(type == ROID_EDGE) return (p-1);
		//	same as for a 2d triangle of order p-3
			if(type == ROID_TRIANGLE)
				return ((p>2) ? (BinomialCoefficient<Dim + p-3, p-3>::value) : 0);
		//	same as for a 2d quadrilateral of order p-2
			if(type == ROID_QUADRILATERAL)
				return ((p>1) ? ((p-1)*(p-1)) : 0);
		//	same as for a 3d prism of order p-2
			if(type == ROID_PRISM)
				return ((p>2) ? (BinomialCoefficient<2 + p-2, p-2>::value)*(p-1) : 0);
			else return -1;
		}

	///	returns the number of DoFs on sub-geometric object in dimension and id
		size_t num_dof(int d, size_t id) const
		{
			return num_dof(m_rRef.ref_elem_type(d, id));
		}

	///	returns if the storage needs objects of a given dimension
		static size_t max_num_dof(int d)
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
	///	association to elements
		LocalDoF m_vLocalDoF[nsh];

	///	reference element
		const ReferencePrism& m_rRef;
};

///////////////////////////////////////////////////////////////////////////////
// Pyramid
///////////////////////////////////////////////////////////////////////////////

/// specialization for Pyramid
template <>
template <int TOrder>
class LagrangeLDS<ReferencePyramid, TOrder>
{
	protected:
	///	corresponding local shape function set
		typedef LagrangeLSFS<ReferencePyramid, TOrder> LSFS;

	///	abbreviation for order
		static const size_t p = TOrder;

	public:
	///	number of shapes
		static const size_t nsh = LSFS::nsh;

	///	dimension of reference element
		const static int Dim = LSFS::dim;

	///	order
		static const size_t order = TOrder;

	public:
	///	constructor
		LagrangeLDS() : m_rRef(Provider::get<ReferencePyramid>())
		{
		//	create LocalDoF vector
			SetLagrangeLocalDoFs(m_vLocalDoF, m_rRef, p);
		}

	///	returns the reference dimension
		static int dim() {return Dim;}

	///	returns the type of reference element
		static ReferenceObjectID roid() {return ReferencePyramid::REFERENCE_OBJECT_ID;}

	///	returns the total number of DoFs on the finite element
		static size_t num_dof() {return nsh;};

	///	returns the number of DoFs on a sub-geometric object type
		static int num_dof(ReferenceObjectID type)
		{
			if(type == ROID_VERTEX) return 1;
			if(type == ROID_EDGE) return (p-1);
		//	same as for a 2d triangle of order p-3
			if(type == ROID_TRIANGLE)
				return ((p>2) ? (BinomialCoefficient<Dim + p-3, p-3>::value) : 0);
		//	same as for a 2d quadrilateral of order p-2
			if(type == ROID_QUADRILATERAL)
				return ((p>1) ? ((p-1)*(p-1)) : 0);
		//	same as for a 3d pyramid of order p-2
			if(type == ROID_PYRAMID)
				return ((p>2) ? (NumberOfDoFsOfPyramid<p-2>::value) : 0);
			else return -1;
		}

	///	returns the number of DoFs on sub-geometric object in dimension and id
		size_t num_dof(int d, size_t id) const
		{
			return num_dof(m_rRef.ref_elem_type(d, id));
		}

	///	returns if the storage needs objects of a given dimension
		static size_t max_num_dof(int d)
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
	///	association to elements
		LocalDoF m_vLocalDoF[nsh];

	///	reference element
		const ReferencePyramid& m_rRef;
};


///////////////////////////////////////////////////////////////////////////////
// Hexahedron
///////////////////////////////////////////////////////////////////////////////

/// specialization for Hexahedron
template <>
template <int TOrder>
class LagrangeLDS<ReferenceHexahedron, TOrder>
{
	protected:
	///	corresponding local shape function set
		typedef LagrangeLSFS<ReferenceHexahedron, TOrder> LSFS;

	///	abbreviation for order
		static const size_t p = TOrder;

	public:
	///	number of shapes
		static const size_t nsh = LSFS::nsh;

	///	dimension of reference element
		const static int Dim = LSFS::dim;

	///	order
		static const size_t order = TOrder;

	public:
	///	constructor
		LagrangeLDS() : m_rRef(Provider::get<ReferenceHexahedron>())
		{
		//	create LocalDoF vector
			SetLagrangeLocalDoFs(m_vLocalDoF, m_rRef, p);
		}

	///	returns the reference dimension
		static int dim() {return Dim;}

	///	returns the type of reference element
		static ReferenceObjectID roid() {return ReferenceHexahedron::REFERENCE_OBJECT_ID;}

	///	returns the total number of DoFs on the finite element
		static size_t num_dof() {return nsh;};

	///	returns the number of DoFs on a sub-geometric object type
		static int num_dof(ReferenceObjectID type)
		{
			if(type == ROID_VERTEX) return 1;
			if(type == ROID_EDGE) return (p-1);
			if(type == ROID_QUADRILATERAL) return (p-1)*(p-1);
			if(type == ROID_HEXAHEDRON) return (p-1)*(p-1)*(p-1);
			else return -1;
		}

	///	returns the number of DoFs on sub-geometric object in dimension and id
		size_t num_dof(int d, size_t id) const
		{
			return num_dof(m_rRef.ref_elem_type(d, id));
		}

	///	returns if the storage needs objects of a given dimension
		static size_t max_num_dof(int d)
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
	///	association to elements
		LocalDoF m_vLocalDoF[nsh];

	///	reference element
		const ReferenceHexahedron& m_rRef;
};

} //namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__LOCAL_SHAPE_FUNCTION_SET__LAGRANGE__LAGRANGE_LOCAL_DOF__ */
