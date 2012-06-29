/*
 * reference_mapping_provider.cpp
 *
 *  Created on: 21.07.2011
 *      Author: andreasvogel
 */

#include "common/util/provider.h"
#include "reference_mapping_provider.h"
#include "reference_mapping.h"

namespace ug{


/// wrapper of a ReferenceElementMapping into the virtual base class
template <typename TRefMapping>
class DimReferenceMappingWrapper
	: public DimReferenceMapping<TRefMapping::dim, TRefMapping::worldDim>,
	  	  	  TRefMapping
{
	public:
	///	world dimension (range space dimension)
		static const int worldDim = TRefMapping::worldDim;

	///	reference dimension (domain space dimension)
		static const int dim = TRefMapping::dim;

	public:
	///	returns if mapping is affine
		virtual bool is_linear() const {return TRefMapping::isLinear;}

	///	refresh mapping for new set of corners
		virtual void update(const MathVector<worldDim>* vCorner)
		{
			TRefMapping::update(vCorner);
		}

	///	map local coordinate to global coordinate
		virtual void local_to_global(MathVector<worldDim>& globPos,
		                             const MathVector<dim> locPos) const
		{
			TRefMapping::local_to_global(globPos, locPos);
		}

	///	map n local coordinate to global coordinate
		virtual void local_to_global(MathVector<worldDim>* globPos,
									 const MathVector<dim>* locPos, size_t n) const
		{
			for(size_t ip = 0; ip < n; ++ip)
				TRefMapping::local_to_global(globPos[ip], locPos[ip]);
		}

	///	returns transposed of jacobian
		virtual void jacobian_transposed(MathMatrix<dim, worldDim>& JT,
		                                 const MathVector<dim> locPos) const
		{
			TRefMapping::jacobian_transposed(JT, locPos);
		}

	///	returns transposed of jacobian for n local positions
		virtual void jacobian_transposed(MathMatrix<dim, worldDim>* JT,
										 const MathVector<dim>* locPos, size_t n) const
		{
			if(TRefMapping::isLinear)
			{
				UG_ASSERT(n > 0, "No point specified in mapping.");
				TRefMapping::jacobian_transposed(JT[0], locPos[0]);
				for(size_t ip = 1; ip < n; ++ip)
					JT[ip] = JT[0];
			}
			else
			{
				for(size_t ip = 0; ip < n; ++ip)
					TRefMapping::jacobian_transposed(JT[ip], locPos[ip]);
			}
		}

	///	returns transposed of the inverse of the jacobian
		virtual void jacobian_transposed_inverse(MathMatrix<worldDim, dim>& JTInv,
		                                         const MathVector<dim> locPos) const
		{
			TRefMapping::jacobian_transposed_inverse(JTInv, locPos);
		}

	///	returns transposed of the inverse of the jacobian for n local positions
		virtual void jacobian_transposed_inverse(MathMatrix<worldDim, dim>* JTInv,
												 const MathVector<dim>* locPos, size_t n) const
		{
			if(TRefMapping::isLinear)
			{
				UG_ASSERT(n > 0, "No point specified in mapping.");
				TRefMapping::jacobian_transposed_inverse(JTInv[0], locPos[0]);
				for(size_t ip = 1; ip < n; ++ip)
					JTInv[ip] = JTInv[0];
			}
			else
			{
				for(size_t ip = 0; ip < n; ++ip)
					TRefMapping::jacobian_transposed_inverse(JTInv[ip], locPos[ip]);
			}
		}

	///	returns the determinate of the jacobian
		virtual number jacobian_det(const MathVector<dim> locPos) const
		{
			return TRefMapping::jacobian_det(locPos);
		}

	///	returns the determinate of the jacobian for n local positions
		virtual void jacobian_det(number* det, const MathVector<dim>* locPos, size_t n) const
		{
			if(TRefMapping::isLinear)
			{
				UG_ASSERT(n > 0, "No point specified in mapping.");
				det[0] = TRefMapping::jacobian_det(locPos[0]);
				for(size_t ip = 1; ip < n; ++ip)
					det[ip] = det[0];
			}
			else
			{
				for(size_t ip = 0; ip < n; ++ip)
					det[ip] = TRefMapping::jacobian_det(locPos[ip]);
			}
		}

	///	virtual destructor
		virtual ~DimReferenceMappingWrapper() {}
};


ReferenceMappingProvider::
ReferenceMappingProvider()
{
//	clear mappings
	for(int d = 0; d < 4; ++d)
		for(int rd = 0; rd < 4; ++rd)
			for(int roid = 0; roid < NUM_REFERENCE_OBJECTS; ++roid)
				m_vvvMapping[d][rd][roid] = NULL;

//	set mappings

//	edge
	set_mapping<1,1>(ROID_EDGE, Provider<DimReferenceMappingWrapper<ReferenceMapping<ReferenceEdge, 1> > >::get());
	set_mapping<1,2>(ROID_EDGE, Provider<DimReferenceMappingWrapper<ReferenceMapping<ReferenceEdge, 2> > >::get());
	set_mapping<1,3>(ROID_EDGE, Provider<DimReferenceMappingWrapper<ReferenceMapping<ReferenceEdge, 3> > >::get());

//	triangle
	set_mapping<2,2>(ROID_TRIANGLE, Provider<DimReferenceMappingWrapper<ReferenceMapping<ReferenceTriangle, 2> > >::get());
	set_mapping<2,3>(ROID_TRIANGLE, Provider<DimReferenceMappingWrapper<ReferenceMapping<ReferenceTriangle, 3> > >::get());

//	quadrilateral
	set_mapping<2,2>(ROID_QUADRILATERAL, Provider<DimReferenceMappingWrapper<ReferenceMapping<ReferenceQuadrilateral, 2> > >::get());
	set_mapping<2,3>(ROID_QUADRILATERAL, Provider<DimReferenceMappingWrapper<ReferenceMapping<ReferenceQuadrilateral, 3> > >::get());

//	3d elements
	set_mapping<3,3>(ROID_TETRAHEDRON, Provider<DimReferenceMappingWrapper<ReferenceMapping<ReferenceTetrahedron, 3> > >::get());
	set_mapping<3,3>(ROID_PRISM, Provider<DimReferenceMappingWrapper<ReferenceMapping<ReferencePrism, 3> > >::get());
	set_mapping<3,3>(ROID_PYRAMID, Provider<DimReferenceMappingWrapper<ReferenceMapping<ReferencePyramid, 3> > >::get());
	set_mapping<3,3>(ROID_HEXAHEDRON, Provider<DimReferenceMappingWrapper<ReferenceMapping<ReferenceHexahedron, 3> > >::get());
}


} // end namespace ug
