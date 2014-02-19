/*
 * adaption_surface_grid_function.h
 *
 *  Created on: 14.06.2013
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__FUNCTION_SPACE__ADAPTION_SURFACE_GRID_FUNCTION__
#define __H__UG__LIB_DISC__FUNCTION_SPACE__ADAPTION_SURFACE_GRID_FUNCTION__

#include "grid_function.h"
#include "local_transfer_interface.h"

namespace ug{

template <typename TDomain>
class AdaptionSurfaceGridFunction : public GridObserver
{
	public:
		typedef std::vector<std::vector<number> > Values;
		typedef Attachment<Values>	AValues;

	public:
		AdaptionSurfaceGridFunction(SmartPtr<TDomain> spDomain,
				 bool bObserveStorage = true);

	public:
		template <typename TAlgebra>
		void copy_from_surface(const GridFunction<TDomain,TAlgebra>& rSurfaceFct);

		template <typename TAlgebra>
		void copy_to_surface(GridFunction<TDomain,TAlgebra>& rSurfaceFct);

		AValues value_attachment()	{return m_aValue;}

	protected:
		template <typename TElem, typename TAlgebra>
		void copy_from_surface(const GridFunction<TDomain,TAlgebra>& rSurfaceFct, TElem* elem);
		template <typename TElem, typename TAlgebra>
		void copy_from_surface(const GridFunction<TDomain,TAlgebra>& rSurfaceFct);

		template <typename TElem, typename TAlgebra>
		void copy_to_surface(GridFunction<TDomain,TAlgebra>& rSurfaceFct, TElem* elem);
		template <typename TElem, typename TAlgebra>
		void copy_to_surface(GridFunction<TDomain,TAlgebra>& rSurfaceFct);

	public:
		void prolongate(const GridMessage_Adaption& msg);
		void do_restrict(const GridMessage_Adaption& msg);

	protected:
		template <typename TBaseElem>
		void prolongate(const GridMessage_Adaption& msg, const size_t lvl);

		template <typename TBaseElem>
		void do_restrict(const MGSelector& sel, const GridMessage_Adaption& msg);

		template <typename TBaseElem>
		void select_parents(MGSelector& sel, const GridMessage_Adaption& msg);

	protected:
		class ValueAccessor : public TransferValueAccessor
		{
			public:
				ValueAccessor(AdaptionSurfaceGridFunction<TDomain>& rASGF,
				              size_t fct);

				void access_inner(GridObject* elem);

				void access_closure(GridObject* elem);

			protected:
				template <typename TBaseElem, typename TSubBaseElem>
				void access_closure(TBaseElem* elem);

				template <typename TBaseElem>
				void access_closure(TBaseElem* elem);

			protected:
				AdaptionSurfaceGridFunction<TDomain>& m_rASGF;
				const size_t m_fct;
				bool m_HasDoFs[NUM_GEOMETRIC_BASE_OBJECTS];
		};
		friend class ValueAccessor;

	public:
		///	creates storage when object created
		template <typename TBaseElem>
		inline void obj_created(TBaseElem* elem);

		/// grid observer callbacks
		/// \{
		virtual void vertex_created(Grid* grid, Vertex* vrt, GridObject* pParent = NULL, bool replacesParent = false){obj_created(vrt);}
		virtual void edge_created(Grid* grid, Edge* e, GridObject* pParent = NULL, bool replacesParent = false){obj_created(e);}
		virtual void face_created(Grid* grid, Face* f, GridObject* pParent = NULL, bool replacesParent = false){obj_created(f);}
		virtual void volume_created(Grid* grid, Volume* vol, GridObject* pParent = NULL, bool replacesParent = false){obj_created(vol);}
		/// \}

	protected:
		SmartPtr<TDomain> m_spDomain;
		SmartPtr<MultiGrid> m_spGrid;
		ConstSmartPtr<DoFDistributionInfo> m_spDDInfo;
		int m_ParallelStorageType;
		bool m_bObserveStorage;

		void attach_entries(ConstSmartPtr<DoFDistributionInfo> spDDInfo);

		template <typename TElem> void detach_entries();
		void detach_entries();

		AValues									m_aValue;
		MultiElementAttachmentAccessor<AValues>	m_aaValue;

		std::vector<SmartPtr<IElemProlongation<TDomain> > > m_vpProlong;
		std::vector<SmartPtr<IElemRestriction<TDomain> > > m_vpRestrict;

};

} // end namespace ug

#include "adaption_surface_grid_function_impl.h"

#endif /* __H__UG__LIB_DISC__FUNCTION_SPACE__ADAPTION_SURFACE_GRID_FUNCTION__ */
