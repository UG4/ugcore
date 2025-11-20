/*
 * Copyright (c) 2016:  G-CSC, Goethe University Frankfurt
 * Author: Markus Breit
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
 *
 *  Created on: 2016-12-19
 */

#ifndef UG__LIB_GRID__REFINEMENT__PROJECTORS__NEURITE_PROJECTOR_H
#define UG__LIB_GRID__REFINEMENT__PROJECTORS__NEURITE_PROJECTOR_H

#include "common/types.h"
#include "lib_grid/refinement/projectors/refinement_projector.h"

#include <boost/serialization/split_member.hpp> // for separate load/save methods

namespace ug {

class NeuriteProjector
: public RefinementProjector
{
	public:
		NeuriteProjector();
		NeuriteProjector(SPIGeometry3d geometry);

		virtual ~NeuriteProjector();

		virtual void set_geometry(SPIGeometry3d geometry);

		/// called when a new vertex was created from an old vertex
		virtual number new_vertex(Vertex* vrt, Vertex* parent);

		/// called when a new vertex was created from an old edge
		virtual number new_vertex(Vertex* vrt, Edge* parent);

		/// called when a new vertex was created from an old face
		virtual number new_vertex(Vertex* vrt, Face* parent);

		/// called when a new vertex was created from an old volume
		virtual number new_vertex(Vertex* vrt, Volume* parent);

		/// project a vertex to its model position
		void project(Vertex* vrt);

		/// spline direction at some grid object
		void direction_at_grid_object(vector3& dirOut, GridObject* o) const;


	protected:
		void attach_surf_params();

	public:
		struct Section
		{
			Section() : endParam(0) {}                // constructor for serialization
			Section(number _endParam)   // constructor for search with CompareSections
			: endParam(_endParam) {}

			number endParam;

			number splineParamsX[4];
			number splineParamsY[4];
			number splineParamsZ[4];
			number splineParamsR[4];

			template <typename Archive>
			void serialize(Archive& ar, const unsigned int version)
			{
				ar & endParam;
				for (size_t i = 0; i < 4; ++i)
				{
					ar & splineParamsX[i];
					ar & splineParamsY[i];
					ar & splineParamsZ[i];
					ar & splineParamsR[i];
				}
			}
		};

		struct BranchingPoint;
		struct BranchingRegion
		{
			BranchingRegion() : t(0.0), bp(nullptr) {} // constructor for serialization
			BranchingRegion(number _t)       // constructor for search with CompareBranchingPointEnds
			: t(_t) {}

			/// the axial parameter where other neurite(s) branch off
			number t;

			/// pointer to whole branching point
			SmartPtr<BranchingPoint> bp;

			template <typename Archive>
			void save(Archive & ar, const unsigned int version) const
			{
				ar << t;

				bool owns_bp = bp->vRegions[0] == this;
				ar << owns_bp;
				if (owns_bp)
					ar << *bp;
			}

			template <typename Archive>
			void load(Archive & ar, const unsigned int version)
			{
				// invoke serialization of the base class
				ar >> t;

				bool owns_bp;
				ar >> owns_bp;
				if (owns_bp)
				{
					bp = make_sp(new BranchingPoint());
					ar >> *bp;
				}
			}

			BOOST_SERIALIZATION_SPLIT_MEMBER()
		};


		struct BranchingPoint
		{
			std::vector<uint32_t> vNid;
			std::vector<BranchingRegion*> vRegions;

			template <typename Archive>
			void serialize(Archive& ar, const unsigned int version)
			{
				// please note: We could simply use
				// ar & vNid
				// but this somehow takes up a huge number of template recursions
				size_t sz = vNid.size(); // this does not hurt in the load-case
				ar & sz;
				vNid.resize(sz);    // this does not hurt in the save-case
				for (size_t i = 0; i < sz; ++i)
					ar & vNid[i];
			}
		};

		struct SomaPoint {
			vector3 soma;
			number radius;
			SomaPoint(const vector3& soma, number radius) : soma(soma), radius(radius) {}
		};

		struct SomaBranchingRegion {
			SmartPtr<SomaPoint> somaPt;
			number radius;
			vector3 center;
			number t;
			vector3 bp;

			template <typename Archive>
			void serialize(Archive& ar, const unsigned int version) {
				ar & radius;
				ar & t;
			}

			template <typename Archive>
			void load(Archive & ar, const unsigned int version)
			{
				ar >> radius;
				ar >> t;
			}

			SomaBranchingRegion(const vector3& center, number radius, number t)
	  		 	  : radius(radius),
				  t(t),
				  bp(center)
			{}

			SomaBranchingRegion(number t)
				 : radius(-1),
				  t(t)
			{}

			SomaBranchingRegion()
				 : radius(-1),
				  t(-1)
			{}
		};

		/*!
		 * \brief Mapping of model to surface vertices
		 */
		struct Mapping {
			/// TODO: Mapping should store unique indices of 1d vertices from SWC file
			vector3 v1; /// start vertex of edge
			vector3 v2; /// end vertex of edge
			number lambda; /// projection parameter lambda
		};

		struct Neurite
		{
			vector3 refDir;
			std::vector<Section> vSec;
			std::vector<BranchingRegion> vBR;
			std::vector<Section> vSomaSec;
			std::vector<SomaBranchingRegion> vSBR;

			template <typename Archive>
			void serialize(Archive& ar, const unsigned int version)
			{
				/// neurite information
				ar & refDir;

				size_t sz = vSec.size();
				ar & sz;
				vSec.resize(sz);
				for (size_t i = 0; i < sz; ++i)
					ar & vSec[i];

				sz = vBR.size();
				ar & sz;
				vBR.resize(sz);
				for (size_t i = 0; i < sz; ++i)
					ar & vBR[i];

				/// soma information
				sz = vSomaSec.size();
				ar & sz;
				vSomaSec.resize(sz);
				for (size_t i = 0; i < sz; ++i)
					ar & vSomaSec[i];

				sz = vSBR.size();
				ar & sz;
				vSBR.resize(sz);
				for (size_t i = 0; i < sz; ++i) {
					ar & vSBR[i];
				}
			}
		};

		struct SurfaceParams
		{
			/**
			 * @brief Neurite ID a vertex belongs to
			 *
			 * Only the low 20 bits are used for ordinary vertices, i.e.,
			 * vertices that do not belong to two neurites at a branching point.
			 * In branching points, the high 4 bits encode which of the child
			 * neurites share the vertex: bit 28 for child 0, bit 29 for child 1,
			 * bit 30 for child 2 and bit 31 for child 3.
			 * Bits 20-27 are used to encode the branching region on the parent
			 * neurite.
			 * There cannot be more than 256 branching regions per neurite.
			 * There cannot be more than 4 children per branching point.
			 */
			uint32 neuriteID;
			float axial;
			float angular;
			float radial;
		};


	public:
		struct CompareSections
		{
			bool operator () (const Section& a, const Section& b)
			{return a.endParam < b.endParam;}
		};
		struct CompareBranchingRegionEnds
		{
			bool operator () (const BranchingRegion& a, const BranchingRegion& b)
			{return a.t < b.t;}
		};

		// helper struct for branching point calculations
		struct BPProjectionHelper
		{
			number start;
			number end;
			std::vector<Section>::const_iterator sec_start;
			const SomaBranchingRegion* sr;
		};

		void debug_neurites() const;

		struct CompareSomaBranchingRegionsEnd {
			bool operator () (const SomaBranchingRegion& a, const SomaBranchingRegion& b) {
				return a.t < b.t;
			}
		};

	public:
		std::vector<Neurite>& neurites();
		std::vector<SomaBranchingRegion>& somata();
		const Neurite& neurite(uint32_t nid) const;
		Grid::VertexAttachmentAccessor<Attachment<SurfaceParams> >& surface_params_accessor();
		const Grid::VertexAttachmentAccessor<Attachment<SurfaceParams> >& surface_params_accessor() const;
		const std::vector<std::pair<number, number> >& quadrature_points() const;
		void average_pos_from_parent(vector3& posOut, const IVertexGroup* const parent) const;
		vector3 position(Vertex* vrt) const;

		number axial_range_around_branching_region
		(
			uint32_t nid,
			size_t brInd,
			number numberOfRadii = 5.0
		) const;

		number axial_range_around_soma_region
		(
			const SomaBranchingRegion& sr,
			size_t numRadii,
			size_t nid,
			Vertex* vrt
		) const;

		bool is_in_axial_range_around_soma_region
		(
			const SomaBranchingRegion& sr,
			size_t nid,
			Vertex* vrt
		) const;

		void print_surface_params(const Vertex* const v) const;

	protected:
		std::vector<Section>::const_iterator get_section_iterator(uint32_t nid, float t) const;
		std::vector<NeuriteProjector::Section>::const_iterator get_soma_section_iterator(uint32_t nid, float t) const;

		void prepare_quadrature();

		void average_params
		(
			uint32_t& neuriteID,
			float& t,
			float& angle,
			float& radius,
			const IVertexGroup* parent
		) const;

		number push_into_place(Vertex* vrt, const IVertexGroup* parent);

	private:
		Attachment<SurfaceParams> m_aSurfParams;
		Grid::VertexAttachmentAccessor<Attachment<SurfaceParams> > m_aaSurfParams;

		/**
		 * @brief storage for cubic splines:
		 * Vector of sorted vectors of sections; each inner vector represents a single neurite.
		 * Any complete neurite is parameterized by t in [0,1]. Each section in the neurite vector
		 * consists of the parameter t at which the section ends and the sixteen coefficients
		 * describing the spline in each dimension (monomial basis {(t_(i+1) - t)^i}_i).
		**/
		std::vector<Neurite> m_vNeurites; //!< spline information for neurites

		std::vector<std::pair<number, number> > m_qPoints; //!< for quadrature when projecting within branching points

		/*!
		 * \brief check if an element is inside a sphere
		 * \param[in] elem
		 * \param[in] center
		 * \param[in] radius
		 * \return \c true if element is inside or on sphere else false
		 */
		template <typename TElem>
		bool IsElementInsideSphere
		(
			const TElem* elem,
			const vector3& center,
			const number radius
		) const
		{
			Grid::VertexAttachmentAccessor aaPos(this->geometry()->grid(), aPosition);
			vector3 c = CalculateCenter(elem, aaPos);
			vector3 diff;
			VecSubtract(diff, c, center);
			VecPow(diff, diff, 2.0);
			number s = diff.x() + diff.y() + diff.z();
			return (s <= (sq(radius) + SMALL));
		}

		friend class boost::serialization::access;
		template<typename Archive>
		void save(Archive & ar, const unsigned int version) const
		{
			UG_EMPTY_BASE_CLASS_SERIALIZATION(NeuriteProjector, RefinementProjector);

			// only write if data is to be written
			if(ArchiveInfo<Archive>::TYPE == AT_DATA)
			{
				// neurites
				size_t sz = m_vNeurites.size();
				ar << sz;
				for (size_t i = 0; i < sz; ++i) {
					ar << m_vNeurites[i];
				}
			}

			// do not do anything otherwise
			else if (ArchiveInfo<Archive>::TYPE == AT_GUI)
			{}
		}

		template<typename Archive>
		void load(Archive & ar, const unsigned int version)
		{
			UG_EMPTY_BASE_CLASS_SERIALIZATION(NeuriteProjector, RefinementProjector);

			// neurites
			size_t sz;
			ar >> sz;
			m_vNeurites.resize(sz);
			for (size_t i = 0; i < sz; ++i) {
				ar >> m_vNeurites[i];
			}

			// reconstruct uninitialized pointers in branching points/ranges
			size_t nNeurites = m_vNeurites.size();
			for (size_t n = 0; n < nNeurites; ++n)
			{
				Neurite& neurite = m_vNeurites[n];
				size_t nBR = neurite.vBR.size();
				for (size_t b = 0; b < nBR; ++b)
				{
					BranchingRegion& br = neurite.vBR[b];
					SmartPtr<BranchingPoint> bp = br.bp;
					if (bp.valid() && !bp->vRegions.size())
					{
						size_t nReg = bp->vNid.size();
						bp->vRegions.resize(nReg);
						bp->vRegions[0] = &br;
						for (size_t r = 1; r < nReg; ++r)
						{
							UG_ASSERT(m_vNeurites[bp->vNid[r]].vBR.size(),
								"No branching region in neurite where there should be one.")
							bp->vRegions[r] = &m_vNeurites[bp->vNid[r]].vBR[0];
							bp->vRegions[r]->bp = bp;
						}
					}
				}
			}
			// debug_neurites();
		}
		BOOST_SERIALIZATION_SPLIT_MEMBER()
};

// DO NOT CHANGE LINES BELOW! Needed for serialization! //
std::ostream& operator << (std::ostream& os, const NeuriteProjector::SurfaceParams& surfParams);
std::istream& operator >> (std::istream& in, NeuriteProjector::SurfaceParams& surfParams);
DECLARE_ATTACHMENT_INFO_TRAITS(Attachment<NeuriteProjector::SurfaceParams>, "NeuriteProjectorSurfaceParams");
std::ostream& operator << (std::ostream& os, const NeuriteProjector::Mapping&  mapping);
std::istream& operator >> (std::istream& in, NeuriteProjector::Mapping& mapping);
DECLARE_ATTACHMENT_INFO_TRAITS(Attachment<NeuriteProjector::Mapping>, "Mapping");
} // namespace ug

#endif