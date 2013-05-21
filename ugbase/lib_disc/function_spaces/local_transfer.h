/*
 * local_transfer.h
 *
 *  Created on: 07.03.2012
 *      Author: Andreas Vogel, Christian Wehner
 */

#ifndef LOCAL_TRANSFER_H_
#define LOCAL_TRANSFER_H_

#include "lib_disc/function_spaces/local_transfer_interface.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"
#include "lib_disc/local_finite_element/local_shape_function_set.h"
#include "lib_disc/function_spaces/grid_function_util.h"

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
				default: UG_THROW("Reference Object type not found.");
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
				default: UG_THROW("Reference Object type not found.");
			}
		}
		void restrict_values(EdgeBase* elem, GeometricObject* parent, const MGDoFDistribution& mgDD) const {}
		void restrict_values(Face* elem, GeometricObject* parent, const MGDoFDistribution& mgDD) const {}
		void restrict_values(Volume* elem, GeometricObject* parent, const MGDoFDistribution& mgDD) const {}

	protected:
		size_t m_fct;
};

template <typename TAlgebra,typename TGridFunction>
class ElementLocalTransfer : public ILocalTransferAlgebra<TAlgebra>
{
		///	domain type
		typedef typename TGridFunction::domain_type domain_type;
	
		///	world dimension
		static const int dim = domain_type::dim;
	
		/// element type
		typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object elem_type;
	
		///	grid type
		typedef typename domain_type::grid_type grid_type;
	public:
		typedef typename TAlgebra::vector_type vector_type;
		using ILocalTransferAlgebra<TAlgebra>::m_pVec;
	private:
		elem_type* m_coarseElem;
		ReferenceObjectID m_coarseRoid;
		LocalVector m_coarseLocU;
		typedef DimReferenceMapping<dim, dim> DimReferenceMapping_type;
		DimReferenceMapping<dim, dim>* m_coarseMap;
		ConstSmartPtr<domain_type> m_spDomain;
		std::vector<LFEID> m_vLFEID;
		typedef LocalShapeFunctionSet<dim> lsfs_type;
		grid_type* m_grid;
	public:
		///	Constructor
		ElementLocalTransfer(SmartPtr<TGridFunction> spGridFct)
		{
			m_spDomain = spGridFct->domain();
			m_vLFEID.resize(spGridFct->num_fct());
			for(size_t fct = 0; fct < spGridFct->num_fct(); ++fct)
					m_vLFEID[fct] = spGridFct->local_finite_element_id(fct);
			domain_type& domain = *spGridFct->domain().get();
			grid_type& grid = *domain.grid();
			m_grid = &grid;
		}
		
		///	returns if prolongation is performed on type
		virtual bool prolongation_needed(GeometricBaseObject gbo) const
		{
			return true;
		}
		
		///	returns if restriction is performed on type
		virtual bool restriction_needed(GeometricBaseObject gbo) const
		{
			return true;
		}

		void prolongate_values_general(GeometricObject* child, GeometricObject* parent, const MGDoFDistribution& mgDD) const
		{		
			std::vector<MultiIndex<2> > vFineMI;
			std::vector<MultiIndex<2> > vCoarseMI;
			typename grid_type::template traits<elem_type>::secure_container assoElements;
			elem_type* coarseElem;
			ReferenceObjectID coarseRoid;
			LocalVector coarseLocU;
			typedef DimReferenceMapping<dim, dim> DimReferenceMapping_type;
			DimReferenceMapping<dim, dim>* coarseMap;
			bool isElemType = true;
			if (dim==2){
					switch(parent->reference_object_id()){
						case ROID_VERTEX: 
							isElemType = false;
							break;
						case ROID_EDGE: 
							isElemType = false;
							break;
						case ROID_UNKNOWN:break;
						case ROID_TRIANGLE:break;
						case ROID_PRISM:break;
						case ROID_PYRAMID:break;
						case NUM_REFERENCE_OBJECTS:break;
						case ROID_QUADRILATERAL:break;
						case ROID_TETRAHEDRON:break;
						case ROID_HEXAHEDRON:break;
						default: UG_THROW("Reference Object type not found.");
					} 
			} else 
				if (dim==3){
					switch(parent->reference_object_id()){
						case ROID_VERTEX: 
							isElemType = false;
							break;
						case ROID_EDGE: 
							isElemType = false;
							break;
						case ROID_QUADRILATERAL: 
							isElemType = false;
							break;
						case ROID_TRIANGLE: 
							isElemType = false;
							break;
						case ROID_UNKNOWN:break;
						case ROID_TETRAHEDRON:break;
						case ROID_HEXAHEDRON:break;
						case ROID_PRISM:break;
						case ROID_PYRAMID:break;
						case NUM_REFERENCE_OBJECTS:break;
						default: UG_THROW("Reference Object type not found.");
					} 
				}
			if (isElemType==false){
				m_grid->associated_elements(assoElements,parent);
				coarseElem =  assoElements[0];
			} else {
				coarseElem = dynamic_cast<elem_type*>(parent);
			}
			coarseRoid = (ReferenceObjectID) coarseElem->reference_object_id();
			// 	get local values for parent (coarse element)
			LocalIndices ind;
			mgDD.indices(coarseElem, ind, false);
			coarseLocU.resize(ind);
			// 	read local values of u in coarse element
			GetLocalVector(coarseLocU, *m_pVec);
			//	get corner coordinates
			std::vector<MathVector<dim> > vCornerCoarse;
			CollectCornerCoordinates(vCornerCoarse, *coarseElem, *m_spDomain);
			//	get Reference Mapping
			coarseMap = &ReferenceMappingProvider::get<dim, dim>(coarseRoid, vCornerCoarse);
			for(size_t fct = 0; fct < mgDD.num_fct(); fct++){
				const LocalShapeFunctionSet<dim>& lsfs = LocalShapeFunctionSetProvider::get<dim>(coarseRoid, m_vLFEID[fct]);
				std::vector<MultiIndex<2> > vCoarseMultInd, vFineMultInd;
				mgDD.inner_multi_indices(child, fct, vFineMultInd);
				if (vFineMultInd.size()==0) continue;
				std::vector<MathVector<dim> > vDoFPos, vLocPos;
				DoFPosition(vDoFPos, child, *m_spDomain, m_vLFEID[fct], dim);
				//	get local position of DoF
				vLocPos.resize(vDoFPos.size());
				for(size_t ip = 0; ip < vLocPos.size(); ++ip) VecSet(vLocPos[ip], 0.0);
				coarseMap->global_to_local(vLocPos, vDoFPos);
				std::vector<std::vector<number> > vvShape;
				//	evaluate coarse shape fct at fine local point
				lsfs.shapes(vvShape, vLocPos);
				for(size_t ip = 0; ip < vFineMultInd.size(); ++ip){
					number interValue = 0;
					for(size_t sh = 0; sh < vvShape[ip].size(); ++sh){
						interValue += vvShape[ip][sh] * coarseLocU(fct,sh);
					}
					BlockRef((*m_pVec)[ vFineMultInd[ip][0] ], vFineMultInd[ip][1]) = interValue;
				}
			}
		}
	
		void prolongate_values(VertexBase* vrt, GeometricObject* parent, const MGDoFDistribution& mgDD) const{ prolongate_values_general(vrt,parent,mgDD); }
		void prolongate_values(EdgeBase* elem, GeometricObject* parent, const MGDoFDistribution& mgDD) const {prolongate_values_general(elem,parent,mgDD); }
		void prolongate_values(Face* elem, GeometricObject* parent, const MGDoFDistribution& mgDD) const { prolongate_values_general(elem,parent,mgDD);}
		void prolongate_values(Volume* elem, GeometricObject* parent, const MGDoFDistribution& mgDD) const {prolongate_values_general(elem,parent,mgDD);}
		void restrict_values(VertexBase* vrt, GeometricObject* parent, const MGDoFDistribution& mgDD) const{}
		void restrict_values(EdgeBase* elem, GeometricObject* parent, const MGDoFDistribution& mgDD) const {}
		void restrict_values(Face* elem, GeometricObject* parent, const MGDoFDistribution& mgDD) const {}
		void restrict_values(Volume* elem, GeometricObject* parent, const MGDoFDistribution& mgDD) const {}
	};

} // end namespace ug

#endif /* LOCAL_TRANSFER_H_ */
