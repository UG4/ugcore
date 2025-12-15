/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include "reference_mapping_provider.h"

#include "common/util/provider.h"
#include "reference_mapping.h"

namespace ug {


/// wrapper of a ReferenceElementMapping into the virtual base class
template <typename TRefMapping>
class DimReferenceMappingWrapper
	: public DimReferenceMapping<TRefMapping::dim, TRefMapping::worldDim>,
	  public TRefMapping
{
	public:
	///	world dimension (range space dimension)
		static constexpr int worldDim = TRefMapping::worldDim;

	///	reference dimension (domain space dimension)
		static constexpr int dim = TRefMapping::dim;

	public:
	///	returns if mapping is affine
		bool is_linear() const override {return TRefMapping::isLinear;}

	///	refresh mapping for new set of corners
		void update(const MathVector<worldDim>* vCorner) override {
			TRefMapping::update(vCorner);
		}

	///	refresh mapping for new set of corners
		void update(const std::vector<MathVector<worldDim> >& vCorner) override {
			TRefMapping::update(vCorner);
		}

	///	map local coordinate to global coordinate
		void local_to_global(MathVector<worldDim>& globPos,
	                     const MathVector<dim>& locPos) const override {
			TRefMapping::local_to_global(globPos, locPos);
		}

	///	map n local coordinate to global coordinate
		void local_to_global(MathVector<worldDim>* vGlobPos,
	                     const MathVector<dim>* vLocPos, size_t n) const override {
			TRefMapping::local_to_global(vGlobPos, vLocPos, n);
		}

	///	map local coordinate to global coordinate for a vector of local positions
		void local_to_global(std::vector<MathVector<worldDim> >& vGlobPos,
	                     const std::vector<MathVector<dim> >& vLocPos) const override {
			TRefMapping::local_to_global(vGlobPos, vLocPos);
		}

	///	map global coordinate to local coordinate
		void global_to_local(MathVector<dim>& locPos,
							 const MathVector<worldDim>& globPos,
							 const size_t maxIter = 1000,
							 const number tol = 1e-10) const override {
			TRefMapping::global_to_local(locPos, globPos, maxIter, tol);
		}

	///	map global coordinate to local coordinate for n local positions
		void global_to_local(MathVector<dim>* vLocPos,
							 const MathVector<worldDim>* vGlobPos, size_t n,
							 const size_t maxIter = 1000,
							 const number tol = 1e-10) const override {
			TRefMapping::global_to_local(vLocPos, vGlobPos, n, maxIter, tol);
		}

	///	map global coordinate to local coordinate for a vector of local positions
		void global_to_local(std::vector<MathVector<dim> >& vLocPos,
							 const std::vector<MathVector<worldDim> >& vGlobPos,
							 const size_t maxIter = 1000,
							 const number tol = 1e-10) const override {
			TRefMapping::global_to_local(vLocPos, vGlobPos, maxIter, tol);
		}

	///	returns jacobian
		void jacobian(MathMatrix<worldDim, dim>& J,
	              const MathVector<dim>& locPos) const override {
			TRefMapping::jacobian(J, locPos);
		}

	///	returns jacobian for n local positions
		void jacobian(MathMatrix<worldDim, dim>* vJ,
	              const MathVector<dim>* vLocPos, size_t n) const override {
			TRefMapping::jacobian(vJ, vLocPos, n);
		}

	///	returns jacobian for a vector of local positions
		void jacobian(std::vector<MathMatrix<worldDim, dim> >& vJ,
	              const std::vector<MathVector<dim> >& vLocPos) const override {
			TRefMapping::jacobian(vJ, vLocPos);
		}

	///	returns transposed of jacobian
		void jacobian_transposed(MathMatrix<dim, worldDim>& JT,
	                         const MathVector<dim>& locPos) const override {
			TRefMapping::jacobian_transposed(JT, locPos);
		}

	///	returns transposed of jacobian for n local positions
		void jacobian_transposed(MathMatrix<dim, worldDim>* vJT,
	                         const MathVector<dim>* vLocPos, size_t n) const override {
			TRefMapping::jacobian_transposed(vJT, vLocPos, n);
		}

	///	returns transposed of jacobian for a vector of positions
		void jacobian_transposed(std::vector<MathMatrix<dim, worldDim> >& vJT,
	                         const std::vector<MathVector<dim> >& vLocPos) const override {
			TRefMapping::jacobian_transposed(vJT, vLocPos);
		}

	///	returns transposed of the inverse of the jacobian
		number jacobian_transposed_inverse(MathMatrix<worldDim, dim>& JTInv,
	                                   const MathVector<dim>& locPos) const override {
			return TRefMapping::jacobian_transposed_inverse(JTInv, locPos);
		}

	///	returns transposed of the inverse of the jacobian for n local positions
		void jacobian_transposed_inverse(MathMatrix<worldDim, dim>* vJTInv,
	                                 const MathVector<dim>* vLocPos, size_t n) const override {
			TRefMapping::jacobian_transposed_inverse(vJTInv, vLocPos, n);
		}

	///	returns transposed of the inverse of the jacobian for n local positions
		void jacobian_transposed_inverse(MathMatrix<worldDim, dim>* vJTInv,
	                                 number* vDet,
	                                 const MathVector<dim>* vLocPos, size_t n) const override {
			TRefMapping::jacobian_transposed_inverse(vJTInv, vDet, vLocPos, n);
		}

	///	returns transposed of the inverse of the jacobian for a vector of positions
		void jacobian_transposed_inverse(std::vector<MathMatrix<worldDim, dim> >& vJTInv,
	                                 const std::vector<MathVector<dim> >& vLocPos) const override {
			TRefMapping::jacobian_transposed_inverse(vJTInv, vLocPos);
		}

	///	returns transposed of the inverse of the jacobian for a vector of positions
		void jacobian_transposed_inverse(std::vector<MathMatrix<worldDim, dim> >& vJTInv,
	                                 std::vector<number>& vDet,
	                                 const std::vector<MathVector<dim> >& vLocPos) const override {
			TRefMapping::jacobian_transposed_inverse(vJTInv, vDet, vLocPos);
		}

	///	returns the determinate of the jacobian
		[[nodiscard]] number sqrt_gram_det(const MathVector<dim>& locPos) const override {
			return TRefMapping::sqrt_gram_det(locPos);
		}

	///	returns the determinate of the jacobian for n local positions
		void sqrt_gram_det(number* vDet,
	                   const MathVector<dim>* vLocPos, size_t n) const override {
			TRefMapping::sqrt_gram_det(vDet, vLocPos, n);
		}

	///	returns the determinate of the jacobian for a vector of local positions
		void sqrt_gram_det(std::vector<number>& vDet,
	                   const std::vector<MathVector<dim> >& vLocPos) const override {
			TRefMapping::sqrt_gram_det(vDet, vLocPos);
		}

	///	virtual destructor
		~DimReferenceMappingWrapper() override = default;
};


ReferenceMappingProvider::
ReferenceMappingProvider()
{
//	clear mappings
	for(int d = 0; d < 4; ++d)
		for(int rd = 0; rd < 4; ++rd)
			for(int roid = 0; roid < NUM_REFERENCE_OBJECTS; ++roid)
				m_vvvMapping[d][rd][roid] = nullptr;

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
	set_mapping<3,3>(ROID_OCTAHEDRON, Provider<DimReferenceMappingWrapper<ReferenceMapping<ReferenceOctahedron, 3> > >::get());
}


} // end namespace ug
