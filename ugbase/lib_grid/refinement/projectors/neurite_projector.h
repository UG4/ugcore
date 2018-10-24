/*
 * neurite_projector.h
 *
 *  Created on: 19.12.2016
 *      Author: mbreit
 */

#ifndef UG__LIB_GRID__REFINEMENT__PROJECTORS__NEURITE_PROJECTOR_H
#define UG__LIB_GRID__REFINEMENT__PROJECTORS__NEURITE_PROJECTOR_H

#include "common/types.h"
#include "lib_grid/refinement/projectors/refinement_projector.h"
#include "lib_grid/tools/copy_attachment_handler.h"
//#include "lib_disc/quadrature/gauss_legendre/gauss_legendre.h"

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

		///	called when a new vertex was created from an old edge.
		virtual number new_vertex(Vertex* vrt, Edge* parent);

		///	called when a new vertex was created from an old face.
		virtual number new_vertex(Vertex* vrt, Face* parent);

		///	called when a new vertex was created from an old face.
		virtual number new_vertex(Vertex* vrt, Volume* parent);

		/// spline direction at some grid object
		void direction_at_grid_object(vector3& dirOut, GridObject* o) const;

	protected:
		void attach_surf_params();

	public:
		struct Section
		{
		    Section() {}                // constructor for serialization
		    Section(number _endParam)   // constructor for search with CompareSections
		    : endParam(_endParam) {}

		    number endParam;

            number splineParamsX[4];
            number splineParamsY[4];
            number splineParamsZ[4];
            number splineParamsR[4];

            template <class Archive>
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
            BranchingRegion() : tstart(0.0), tend(0.0), bp(SPNULL) {} // constructor for serialization
            BranchingRegion(number _tend)       // constructor for search with CompareBranchingPointEnds
		    : tstart(_tend), tend(_tend) {}

		    /// the parameter range for this branching point
            number tstart;
            number tend;

            /// pointer to whole branching point
            SmartPtr<BranchingPoint> bp;

            template<class Archive>
            void save(Archive & ar, const unsigned int version) const
            {
                ar << tstart;
                ar << tend;

                // the following won't work:
                //bool owns_bp = bp->vRegions[0] == this;
                // pointers are not identical if constructed outside of projector
            	// but this is unsafe:
                bool owns_bp = tstart == bp->vRegions[0]->tstart && tend == bp->vRegions[0]->tend;
                ar << owns_bp;
                if (owns_bp)
                	ar << *bp;
            }

            template<class Archive>
            void load(Archive & ar, const unsigned int version)
            {
                // invoke serialization of the base class
                ar >> tstart;
                ar >> tend;

                bool owns_bp;
                ar >> owns_bp;
                if (owns_bp)
                {
                    bp = make_sp(new BranchingPoint());
                    ar >> *bp;
                }
            }

            BOOST_SERIALIZATION_SPLIT_MEMBER();
		};


        struct BranchingPoint
        {
	        std::vector<size_t> vNid;
            std::vector<BranchingRegion*> vRegions;

            template <class Archive>
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


		struct Neurite
        {
		    vector3 refDir;
            std::vector<Section> vSec;
            std::vector<BranchingRegion> vBR;

            template <class Archive>
            void serialize(Archive& ar, const unsigned int version)
            {
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
            }
        };

		struct SurfaceParams
		{
		    uint32 neuriteID;
		    float axial;
		    float angular;
		    float radial;
		};


    public:
        struct CompareSections
        {
            bool operator()(const Section& a, const Section& b)
            {return a.endParam < b.endParam;}
        };
        struct CompareBranchingRegionEnds
        {
            bool operator()(const BranchingRegion& a, const BranchingRegion& b)
            {return a.tend < b.tend;}
        };

        // helper struct for branching point calculations
        struct BPProjectionHelper
        {
            number start;
            number end;
            std::vector<Section>::const_iterator sec_start;
        };

        void debug_neurites() const;

	public:
        void add_neurite(const Neurite& n);
        const Neurite& neurite(size_t nid) const;

        Grid::VertexAttachmentAccessor<Attachment<SurfaceParams> >& surface_params_accessor();
        const Grid::VertexAttachmentAccessor<Attachment<SurfaceParams> >& surface_params_accessor() const;

        const std::vector<std::pair<number, number> >& quadrature_points() const;

        void average_pos_from_parent(vector3& posOut, const IVertexGroup* parent) const;

        vector3 position(Vertex* vrt) const;

    protected:
        const Section& get_section(uint32 nid, float t) const;

        void prepare_quadrature();

		void average_params
		(
			size_t& neuriteID,
			float& t,
			float& angle,
			float& radius,
			const IVertexGroup* parent
		) const;

		number push_into_place(Vertex* vrt, const IVertexGroup* parent);

	private:
		Attachment<SurfaceParams> m_aSurfParams;
        Grid::VertexAttachmentAccessor<Attachment<SurfaceParams> > m_aaSurfParams;

        /// handles propagation of surface param attachment to children on higher levels
        CopyAttachmentHandler<Vertex, Attachment<SurfaceParams> > m_cah;

        /**
		 * @brief storage for cubic splines:
		 * Vector of sorted vectors of sections; each inner vector represents a single neurite.
		 * Any complete neurite is parameterized by t in [0,1]. Each section in the neurite vector
		 * consists of the parameter t at which the section ends and the sixteen coefficients
		 * describing the spline in each dimension (monomial basis {(t_(i+1) - t)^i}_i).
		**/
		std::vector<Neurite> m_vNeurites;

		/// for quadrature when projecting within branching points
		//size_t m_quadOrder;

        std::vector<std::pair<number, number> > m_qPoints;



		friend class boost::serialization::access;

        template<class Archive>
        void save(Archive & ar, const unsigned int version) const
        {
            UG_EMPTY_BASE_CLASS_SERIALIZATION(NeuriteProjector, RefinementProjector);

            // only write if data is to be written
            if(ArchiveInfo<Archive>::TYPE == AT_DATA)
            {
				size_t sz = m_vNeurites.size();
				ar << sz;
				for (size_t i = 0; i < sz; ++i)
					ar << m_vNeurites[i];
            }

            // do not do anything otherwise
            else if (ArchiveInfo<Archive>::TYPE == AT_GUI)
            {}

            //ar << m_quadOrder;
        }

        template<class Archive>
        void load(Archive & ar, const unsigned int version)
        {
            UG_EMPTY_BASE_CLASS_SERIALIZATION(NeuriteProjector, RefinementProjector);

            size_t sz;
            ar >> sz;
            m_vNeurites.resize(sz);
            for (size_t i = 0; i < sz; ++i)
                ar >> m_vNeurites[i];

            //ar >> m_quadOrder;

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

            //debug_neurites();
        }

        BOOST_SERIALIZATION_SPLIT_MEMBER();
};

// DO NOT CHANGE! Needed for serialization! //
std::ostream& operator<<(std::ostream &os, const NeuriteProjector::SurfaceParams& surfParams);
std::istream& operator>>(std::istream& in, NeuriteProjector::SurfaceParams& surfParams);
DECLARE_ATTACHMENT_INFO_TRAITS(Attachment<NeuriteProjector::SurfaceParams>, "NeuriteProjectorSurfaceParams");

} // namespace ug


#endif // UG__LIB_GRID__REFINEMENT__PROJECTORS__NEURITE_PROJECTOR_H
