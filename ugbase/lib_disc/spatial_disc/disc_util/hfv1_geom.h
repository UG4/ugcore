/*
 * hfv1_geom.h
 *
 *  Created on: 08.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__HANGING_FINITE_VOLUME_GEOMETRY__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__HANGING_FINITE_VOLUME_GEOMETRY__

// extern libraries
#include <cmath>
#include <vector>

// other ug4 modules
#include "common/common.h"

// library intern includes
#include "lib_disc/reference_element/reference_element.h"
#include "lib_disc/local_finite_element/local_shape_function_set.h"
#include "fv_util.h"
#include "fv1_geom.h"

namespace ug{

template <	typename TElem, int TWorldDim>
class HFV1Geometry : public FVGeometryBase{
	private:
	/// type of reference element
		typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

	/// number of SubControlVolumes
		static const size_t m_numNaturalSCV = ref_elem_type::numCorners;

	/// number of SubControlVolumeFaces
		static const size_t m_numNaturalSCVF = ref_elem_type::numEdges;

	public:
	/// dimension of reference element
		static const int dim = ref_elem_type::dim;

	/// dimension of world
		static const int worldDim = TWorldDim;

	/// Hanging node flag: this Geometry does support hanging nodes
		static const bool usesHangingNodes = true;

	/// flag indicating if local data may change
		static const bool staticLocalData = true;

	protected:
		struct MidID
		{
				MidID() : dim(0), id(0) {};
				MidID(size_t dim_, size_t id_) : dim(dim_), id(id_) {};
				size_t dim;
				size_t id;
		};

	public:
		class SCVF
		{
			public:
			/// dimension of reference element
				static const int dim = ref_elem_type::dim;

			/// dimension of world
				static const int worldDim = TWorldDim;

			private:
			/// let outer class access private members
				friend class HFV1Geometry<TElem, TWorldDim>;

			/// number of integration points
				static const size_t m_numIP = 1;

			/// Number of corners of scvf
				static const size_t m_numCorners = hfv1_traits<ref_elem_type, TWorldDim>::NumCornersOfSCVF;

			public:
				SCVF() {};

			/// index of SubControlVolume on one side of the scvf
				inline size_t from() const {return m_from;}

			/// index of SubControlVolume on one side of the scvf
				inline size_t to() const {return m_to;}

			/// number of integration points on scvf
				inline size_t num_ip() const {return m_numIP;}

			/// local integration point of scvf
				inline const MathVector<dim>& local_ip() const {return localIP;}

			/// global integration point of scvf
				inline const MathVector<worldDim>& global_ip() const {return globalIP;}

			/// normal on scvf (points direction "from"->"to"). Norm is equal to area
				inline const MathVector<worldDim>& normal() const {return Normal;} // includes area

			/// Transposed Inverse of Jacobian in integration point
				inline const MathMatrix<worldDim,dim>& JTInv() const {return JtInv;}

			/// Determinante of Jacobian in integration point
				inline number detJ() const {return detj;}

			/// number of shape functions
				inline size_t num_sh() const {return vShape.size();}

			/// value of shape function i in integration point
				inline number shape(size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid ip index"); return vShape[sh];}

			/// value of local gradient of shape function i in integration point
				inline const MathVector<dim>& local_grad(size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid ip index"); return localGrad[sh];}

			/// vector of local gradients in ip point
				inline const MathVector<dim>* local_grad_vector() const {return &localGrad[0];}

			/// value of global gradient of shape function i in integration point
				inline const MathVector<worldDim>& global_grad(size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid ip index"); return globalGrad[sh];}

			/// vector of global gradients in ip point
				inline const MathVector<worldDim>* global_grad_vector() const {return &globalGrad[0];}

			/// number of corners, that bound the scvf
				inline size_t num_corners() const {return m_numCorners;}

			/// return local corner number i
				inline const MathVector<dim>& local_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid corner index."); return m_vLocPos[co];}

			/// return glbal corner number i
				inline const MathVector<worldDim>& global_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid corner index."); return m_vGloPos[co];}

			private:
				// This scvf separates the scv with the ids given in "from" and "to"
				// The computed normal points in direction from->to
				size_t m_from, m_to;

				// ordering is:
				// 2D: edgeMidPoint, CenterOfElement
				// 3D: edgeMidPoint, Side one, CenterOfElement, Side two
				MathVector<dim> m_vLocPos[m_numCorners]; // local corners of scvf
				MathVector<worldDim> m_vGloPos[m_numCorners]; // global corners of scvf
				MidID m_midId[m_numCorners]; // dimension and id of object, that's midpoint bounds the scvf

				// scvf part
				MathVector<dim> localIP; // local integration point
				MathVector<worldDim> globalIP; // global intergration point
				MathVector<worldDim> Normal; // normal (incl. area)

				// shapes and derivatives
				std::vector<number> vShape; // shapes at ip
				std::vector<MathVector<dim> > localGrad; // local grad at ip
				std::vector<MathVector<worldDim> > globalGrad; // global grad at ip
				MathMatrix<worldDim,dim> JtInv; // Jacobian transposed at ip
				number detj; // Jacobian det at ip
		};

		class SCV
		{
			private:
			//  let outer class access private members
				friend class HFV1Geometry<TElem, TWorldDim>;

			/// number of integration points
				static const size_t m_numIP = 1;

			/// Number of corners of scvf
				static const size_t m_maxNumCorners = hfv1_traits<ref_elem_type, TWorldDim>::MaxNumCornersOfSCV;

			/// type of element the subcontrol volume represents
				typedef typename hfv1_traits<ref_elem_type, TWorldDim>::scv_type scv_type;

			public:
				SCV() : m_numCorners(m_maxNumCorners) {};

			/// node id that this scv is associated to
				inline size_t node_id() const {return nodeId;}

			/// number of integration points
				inline size_t num_ip() const {return m_numIP;}

			/// local integration point of scv
				inline const MathVector<dim>& local_ip() const {return m_vLocPos[0];}

			/// global integration point
				inline const MathVector<worldDim>& global_ip() const {return m_vGloPos[0];}

			/// volume of scv
				inline number volume() const {return vol;}

			/// number of corners, that bound the scvf
				inline size_t num_corners() const {return m_numCorners;}

			/// return local corner number i
				inline const MathVector<dim>& local_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid corner index."); return m_vLocPos[co];}

			/// return global corner number i
				inline const MathVector<worldDim>& global_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid corner index."); return m_vGloPos[co];}

			private:
				size_t nodeId; // node id of associated node
				size_t m_numCorners;

				// The ordering is: Corner, ...
				MathVector<dim> m_vLocPos[m_maxNumCorners]; // local position of node
				MathVector<worldDim> m_vGloPos[m_maxNumCorners]; // global position of node
				MidID m_midId[m_maxNumCorners]; // dimension and id of object, that's midpoint bounds the scv
				number vol;
		};


	public:
	///	constructor
		HFV1Geometry();

	///	update values for an element
		void update(TElem* elem, const MathVector<worldDim>* vCornerCoords,
		            			 const ISubsetHandler* ish = NULL);

	// 	debug output
		void print();

	public:
	/// number of SubControlVolumeFaces
		inline size_t num_scvf() const {return m_vSCVF.size();}

	/// const access to SubControlVolumeFace number i
		inline const SCVF& scvf(size_t i) const
			{UG_ASSERT(i < num_scvf(), "Invalid Index."); return m_vSCVF[i];}

	/// number of SubControlVolumes
		size_t num_scv() const {return m_vSCV.size();}

	/// const access to SubControlVolume number i
		inline const SCV& scv(size_t i) const
			{UG_ASSERT(i < num_scv(), "Invalid Index."); return m_vSCV[i];}

	/// number of shape functions
		//\TODO:	This is not yet like it should be and only conserving behavior as it has been.
		// 			In fact, this should of course return the number of shapes, which is
		//			generally NOT equal to the number of SCVs.
		inline size_t num_sh() const
		{
			//if (m_pElem->reference_object_id() == ROID_PYRAMID)
			//	UG_THROW("Function num_sh() for hanging nodes not yet implemented for pyramids.");
			return num_scv();
		};

	public:
	/// returns number of all scv ips
		size_t num_scvf_ips() const {return m_vGlobSCVFIP.size();}

	/// returns all ips of scv as they appear in scv loop
		const MathVector<worldDim>* scvf_global_ips() const {return &m_vGlobSCVFIP[0];}

	/// returns all ips of scv as they appear in scv loop
		const MathVector<dim>* scvf_local_ips() const {return &m_vLocSCVFIP[0];}

	/// returns number of all scv ips
		size_t num_scv_ips() const {return m_vGlobSCVIP.size();}

	/// returns all ips of scv as they appear in scv loop
		const MathVector<worldDim>* scv_global_ips() const {return &m_vGlobSCVIP[0];}

	/// returns all ips of scv as they appear in scv loop
		const MathVector<dim>* scv_local_ips() const {return &m_vLocSCVIP[0];}

	protected:
		std::vector<MathVector<worldDim> > m_vGlobSCVFIP;
		std::vector<MathVector<dim> > m_vLocSCVFIP;
		std::vector<MathVector<worldDim> > m_vGlobSCVIP;
		std::vector<MathVector<dim> > m_vLocSCVIP;

	protected:
		void copy_local_corners(SCVF& scvf)
		{
			for(size_t i = 0; i < scvf.num_corners(); ++i)
			{
				const size_t dim = scvf.m_midId[i].dim;
				const size_t id = scvf.m_midId[i].id;
				scvf.m_vLocPos[i] = m_locMid[dim][id];
			}
		}

		void copy_global_corners(SCVF& scvf)
		{
			for(size_t i = 0; i < scvf.num_corners(); ++i)
			{
				const size_t dim = scvf.m_midId[i].dim;
				const size_t id = scvf.m_midId[i].id;
				scvf.m_vGloPos[i] = m_gloMid[dim][id];
			}
		}

		void copy_local_corners(SCV& scv)
		{
			for(size_t i = 0; i < scv.num_corners(); ++i)
			{
				const size_t dim_ = scv.m_midId[i].dim;
				const size_t id = scv.m_midId[i].id;
				UG_ASSERT(dim_ <= (size_t)dim, "Dimension wrong");
				UG_ASSERT(id < m_locMid[dim_].size(), "id " << id << " in dim="
				          <<dim_<<" wrong. (size is "<< m_locMid[dim_].size()<<")\n");
				scv.m_vLocPos[i] = m_locMid[dim_][id];
			}
		}

		void copy_global_corners(SCV& scv)
		{
			for(size_t i = 0; i < scv.num_corners(); ++i)
			{
				const size_t dim_ = scv.m_midId[i].dim;
				const size_t id = scv.m_midId[i].id;
				UG_ASSERT(dim_ <= (size_t)dim, "Dimension wrong");
				UG_ASSERT(id < m_gloMid[dim_].size(), "id " << id << " in dim="
				          <<dim_<<" wrong. (size is "<< m_gloMid[dim_].size()<<")\n");
				scv.m_vGloPos[i] = m_gloMid[dim_][id];
			}
		}

		// i = number of side
		void compute_side_midpoints(size_t i,
		                            MathVector<dim>& locSideMid,
		                            MathVector<worldDim>& gloSideMid)
		{
			if(dim>1)
			{
				const size_t coID0 = m_rRefElem.id(2, i, 0, 0);
				locSideMid = m_locMid[0][coID0];
				gloSideMid = m_gloMid[0][coID0];

				// add corner coordinates of the corners of the geometric object
				for(size_t j = 1; j < m_rRefElem.num(2, i, 0); ++j)
				{
					const size_t coID = m_rRefElem.id(2, i, 0, j);
					locSideMid += m_locMid[0][coID];
					gloSideMid += m_gloMid[0][coID];
				}

				// scale for correct averaging
				locSideMid *= 1./(m_rRefElem.num(2, i, 0));
				gloSideMid *= 1./(m_rRefElem.num(2, i, 0));
			}
		}

		// i, j, k, l = number nodes
		void compute_side_midpoints(size_t i, size_t j, size_t k, size_t l,
									MathVector<dim>& locSideMid,
									MathVector<worldDim>& gloSideMid)
		{
			VecAdd(locSideMid, m_locMid[0][i], m_locMid[0][j], m_locMid[0][k], m_locMid[0][l]);
			VecAdd(gloSideMid, m_gloMid[0][i], m_gloMid[0][j], m_gloMid[0][k], m_gloMid[0][l]);

			// scale for correct averaging
			locSideMid *= 0.25;
			gloSideMid *= 0.25;
		}

		// i, j, k = number nodes
		void compute_side_midpoints(size_t i, size_t j, size_t k,
									MathVector<dim>& locSideMid, MathVector<worldDim>& gloSideMid)
		{
			VecAdd(locSideMid, m_locMid[0][i], m_locMid[0][j], m_locMid[0][k]);
			VecAdd(gloSideMid, m_gloMid[0][i], m_gloMid[0][j], m_gloMid[0][k]);

			// scale for correct averaging
			locSideMid *= 1./3.;
			gloSideMid *= 1./3.;
		}

		// returns edgeID of child edge, that has corner co. i is the natural edge index
		size_t get_child_edge_of_corner(size_t i, size_t co)
        {
         	for(size_t e = 0; e < m_vNatEdgeInfo[i].num_child_edges(); ++e)
         	{
         		const size_t childId = m_vNatEdgeInfo[i].child_edge(e);
         		if(m_vNewEdgeInfo[childId].from() == co)
         			return childId;
         		if(m_vNewEdgeInfo[childId].to() == co)
         			return childId;
         	}
         	return -1;
        }

	private:
		struct NewEdgeInfo
		{
			friend class HFV1Geometry<TElem, TWorldDim>;
		public:
			NewEdgeInfo() : m_from(-1), m_to(-1){}
            size_t from() const {return m_from;}
            size_t to() const {return m_to;}
		private:
			size_t m_from;
			size_t m_to;
		};

        struct NatEdgeInfo
        {
                friend class HFV1Geometry<TElem, TWorldDim>;
            public:
            	NatEdgeInfo() : nodeId(-1), numChildEdges(0) {for(size_t i=0; i<2; ++i) childEdge[i] = -1;}
                bool is_hanging() const {return numChildEdges == 2;}
                size_t node_id() const {UG_ASSERT(is_hanging(), "Should only be called, if edge is hanging."); return nodeId;}
                size_t num_child_edges() const {return numChildEdges;}
                size_t child_edge(size_t i) const {UG_ASSERT(i < num_child_edges(), "Wrong id."); return childEdge[i];}

            private:
                size_t nodeId;
                size_t numChildEdges;
                size_t childEdge[2];
        };

    //  help array: maps NaturalEdge -> NodeOnEdge (-1 if no node on edge)
        std::vector<NatEdgeInfo> m_vNatEdgeInfo;
        std::vector<NewEdgeInfo> m_vNewEdgeInfo;

	private:
	// 	pointer to current element
		TElem* m_pElem;

		std::vector<MathVector<dim> > m_locMid[dim+1];
		std::vector<MathVector<worldDim> > m_gloMid[dim+1];

	// 	SubControlVolumeFaces
		std::vector<SCVF> m_vSCVF;

	// 	SubControlVolumes
		std::vector<SCV> m_vSCV;

	// 	Reference Mapping
		ReferenceMapping<ref_elem_type, worldDim> m_rMapping;

	// 	Reference Element
		const ref_elem_type& m_rRefElem;
};

}

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__HANGING_FINITE_VOLUME_GEOMETRY__ */
