/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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
