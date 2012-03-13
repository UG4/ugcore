/*
 * local_transfer.h
 *
 *  Created on: 07.03.2012
 *      Author: andreasvogel
 */

#ifndef LOCAL_TRANSFER_H_
#define LOCAL_TRANSFER_H_

#include "lib_disc/function_spaces/local_transfer_interface.h"

namespace ug{

template <typename TAlgebra>
class P1LocalTransfer : public ILocalTransferAlgebra<TAlgebra>
{
	public:
		typedef typename TAlgebra::vector_type vector_type;
		using ILocalTransferAlgebra<TAlgebra>::m_pVec;

	public:
	///	Constructor, indicating to transfer only one component
		P1LocalTransfer(size_t fct) : m_fct(fct) {}

	///	returns if prolongation is performed on type
		virtual bool prolongation_needed(GeometricBaseObject gbo) const
		{
			return (gbo == VERTEX);
		}

	///	returns if restriction is performed on type
		virtual bool restriction_needed(GeometricBaseObject gbo) const
		{
			return (gbo == VERTEX);
		}

		void prolongate_values(VertexBase* vrt, GeometricObject* parent, const MGDoFDistribution& mgDD) const
		{
			std::vector<MultiIndex<2> > vFineMI;
			std::vector<MultiIndex<2> > vCoarseMI;
			switch(parent->reference_object_id())
			{
				case ROID_VERTEX:
				{
					VertexBase* pParent = static_cast<VertexBase*>(parent);
					mgDD.inner_multi_indices(vrt, m_fct, vFineMI);
					mgDD.inner_multi_indices(pParent, m_fct, vCoarseMI);

					for(size_t i = 0; i < vFineMI.size(); ++i)
						BlockRef((*m_pVec)[ vFineMI[i][0] ], vFineMI[i][1]) =
								BlockRef((*m_pVec)[ vCoarseMI[i][0] ], vCoarseMI[i][1]);
				}
				break;
				case ROID_EDGE:
				{
					mgDD.inner_multi_indices(vrt, m_fct, vFineMI);
					for(size_t i = 0; i < vFineMI.size(); ++i)
						BlockRef((*m_pVec)[ vFineMI[i][0] ], vFineMI[i][1]) = 0.0;

					EdgeBase* pParent = static_cast<EdgeBase*>(parent);
					for(size_t i = 0; i < pParent->num_vertices(); ++i)
					{
						VertexBase* edgeVrt = pParent->vertex(i);
						mgDD.inner_multi_indices(edgeVrt, m_fct, vCoarseMI);

						for(size_t i = 0; i < vFineMI.size(); ++i)
							BlockRef((*m_pVec)[ vFineMI[i][0] ], vFineMI[i][1])
							+= 0.5 * BlockRef((*m_pVec)[ vCoarseMI[i][0] ], vCoarseMI[i][1]);
					}
				}
				break;
				case ROID_QUADRILATERAL:
				{
					mgDD.inner_multi_indices(vrt, m_fct, vFineMI);
					for(size_t i = 0; i < vFineMI.size(); ++i)
						BlockRef((*m_pVec)[ vFineMI[i][0] ], vFineMI[i][1]) = 0.0;

					Face* pParent = static_cast<Face*>(parent);
					for(size_t i = 0; i < pParent->num_vertices(); ++i)
					{
						VertexBase* faceVrt = pParent->vertex(i);
						mgDD.inner_multi_indices(faceVrt, m_fct, vCoarseMI);

						for(size_t i = 0; i < vFineMI.size(); ++i)
							BlockRef((*m_pVec)[ vFineMI[i][0] ], vFineMI[i][1])
							+= 0.25 * BlockRef((*m_pVec)[ vCoarseMI[i][0] ], vCoarseMI[i][1]);
					}
				}
				break;
				case ROID_HEXAHEDRON:
				{
					mgDD.inner_multi_indices(vrt, m_fct, vFineMI);
					for(size_t i = 0; i < vFineMI.size(); ++i)
						BlockRef((*m_pVec)[ vFineMI[i][0] ], vFineMI[i][1]) = 0.0;

					Volume* pParent = static_cast<Volume*>(parent);
					for(size_t i = 0; i < pParent->num_vertices(); ++i)
					{
						VertexBase* hexVrt = pParent->vertex(i);
						mgDD.inner_multi_indices(hexVrt, m_fct, vCoarseMI);

						for(size_t i = 0; i < vFineMI.size(); ++i)
							BlockRef((*m_pVec)[ vFineMI[i][0] ], vFineMI[i][1])
							+= 0.125 * BlockRef((*m_pVec)[ vCoarseMI[i][0] ], vCoarseMI[i][1]);
					}
				}
				break;
				case ROID_TRIANGLE:
				case ROID_TETRAHEDRON:
				case ROID_PRISM:
				case ROID_PYRAMID: /*nothing to do in those cases */ break;
				default: UG_THROW_FATAL("Reference Object type not found.");
			}
		}
		void prolongate_values(EdgeBase* elem, GeometricObject* parent, const MGDoFDistribution& mgDD) const {}
		void prolongate_values(Face* elem, GeometricObject* parent, const MGDoFDistribution& mgDD) const {}
		void prolongate_values(Volume* elem, GeometricObject* parent, const MGDoFDistribution& mgDD) const {}

		void restrict_values(VertexBase* vrt, GeometricObject* parent, const MGDoFDistribution& mgDD) const
		{
			std::vector<MultiIndex<2> > vFineMI;
			std::vector<MultiIndex<2> > vCoarseMI;
			switch(parent->reference_object_id())
			{
				case ROID_VERTEX:
				{
					VertexBase* pParent = static_cast<VertexBase*>(parent);
					mgDD.inner_multi_indices(vrt, m_fct, vFineMI);
					mgDD.inner_multi_indices(pParent, m_fct, vCoarseMI);

					for(size_t i = 0; i < vFineMI.size(); ++i)
						BlockRef((*m_pVec)[ vCoarseMI[i][0] ], vCoarseMI[i][1]) =
								BlockRef((*m_pVec)[ vFineMI[i][0] ], vFineMI[i][1]);
				}
				break;
				case ROID_EDGE:
				case ROID_TRIANGLE:
				case ROID_QUADRILATERAL:
				case ROID_TETRAHEDRON:
				case ROID_PRISM:
				case ROID_PYRAMID:
				case ROID_HEXAHEDRON: /*nothing to do in those cases */ break;
				default: UG_THROW_FATAL("Reference Object type not found.");
			}
		}
		void restrict_values(EdgeBase* elem, GeometricObject* parent, const MGDoFDistribution& mgDD) const {}
		void restrict_values(Face* elem, GeometricObject* parent, const MGDoFDistribution& mgDD) const {}
		void restrict_values(Volume* elem, GeometricObject* parent, const MGDoFDistribution& mgDD) const {}

	protected:
		size_t m_fct;
};


} // end namespace ug

#endif /* LOCAL_TRANSFER_H_ */
