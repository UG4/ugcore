/*
 * local_transfer_interface.h
 *
 *  Created on: 07.03.2012
 *      Author: andreasvogel
 */

#ifndef LOCAL_TRANSFER_INTERFACE_H_
#define LOCAL_TRANSFER_INTERFACE_H_

#include "lib_disc/dof_manager/mg_dof_distribution.h"

namespace ug{

class ILocalTransfer : public GridObserver
{
	public:
	///	constructor
		ILocalTransfer() : m_spMGDD(NULL) {}

	///	virtual destructor
		virtual ~ILocalTransfer() {set_dof_distribution(NULL);}

	///	sets the dof distribution
		void set_dof_distribution(SmartPtr<MGDoFDistribution> spMGDD)
		{
		//	NOTE: 	it is very important, that the GridObservation of the
		//			DoFDistribution is registered BEFORE this Observation,
		//			since the DoFDistribution attaches new indices to the new
		//			element, that are needed by the prolongation / restriction.
		//			Therefore we register at the grid when the DoFDistr is passed.
		//			This ensures, that the observer of the DoFDistr has already
		//			been registered in the constructor of the DoFDistr. The
		//			same way this observer is unregistered before the DD may
		//			be destroyed, since a smart pointer is stored.
			if(m_spMGDD != spMGDD)
			{
				if(m_spMGDD.is_valid()) unregister_observer();
				m_spMGDD = spMGDD;
				if(m_spMGDD.is_valid()) register_observer();
			}
		}

	///	returns if prolongation is performed on type
		virtual bool prolongation_needed(GeometricBaseObject gbo) const = 0;

	///	returns if restriction is performed on type
		virtual bool restriction_needed(GeometricBaseObject gbo) const = 0;

	protected:
		void register_observer()
		{
			int type = OT_GRID_OBSERVER;

			if(prolongation_needed(VERTEX)) type |= OT_VERTEX_OBSERVER;
			if(prolongation_needed(EDGE)) type |= OT_EDGE_OBSERVER;
			if(prolongation_needed(FACE)) type |= OT_FACE_OBSERVER;
			if(prolongation_needed(VOLUME)) type |= OT_VOLUME_OBSERVER;

			if(restriction_needed(VERTEX)) type |= OT_VERTEX_OBSERVER;
			if(restriction_needed(EDGE)) type |= OT_EDGE_OBSERVER;
			if(restriction_needed(FACE)) type |= OT_FACE_OBSERVER;
			if(restriction_needed(VOLUME)) type |= OT_VOLUME_OBSERVER;

			m_spMGDD->multi_grid().register_observer(this, type);
		}

		void unregister_observer()
		{
			m_spMGDD->multi_grid().unregister_observer(this);
		}

	///	underlying DoF Distribution
		SmartPtr<MGDoFDistribution> m_spMGDD;
};

template <typename TImpl>
class ILocalTransferImpl : 	public ILocalTransfer
{
	public:
	///	prolongate
		void prolongate_values(VertexBase* elem, GeometricObject* parent) const {}
		void prolongate_values(EdgeBase* elem, GeometricObject* parent) const {}
		void prolongate_values(Face* elem, GeometricObject* parent) const {}
		void prolongate_values(Volume* elem, GeometricObject* parent) const {}

	///	restrict
		void restrict_values(VertexBase* elem, GeometricObject* parent) const {}
		void restrict_values(EdgeBase* elem, GeometricObject* parent) const {}
		void restrict_values(Face* elem, GeometricObject* parent) const {}
		void restrict_values(Volume* elem, GeometricObject* parent) const {}

	protected:
	protected:
	///	access to implementation
		TImpl& getImpl() {return static_cast<TImpl&>(*this);}

	///	const access to implementation
		const TImpl& getImpl() const {return static_cast<const TImpl&>(*this);}

		template <typename TBaseElem>
		void obj_created(TBaseElem* elem, GeometricObject* pParent, bool replacesParent)
		{
		//	only replacement nothing to do
			if(replacesParent) return;

		//	forward to implementation
			getImpl().prolongate_values(elem, pParent);
		}

		template <typename TBaseElem>
		void obj_to_be_erased(TBaseElem* elem, TBaseElem* replacedBy)
		{
		//	only replacement, nothing to do
			if(replacedBy) return;

		//	identical on coarser grid, indices are copied. nothing to do
			if(m_spMGDD->parent_if_copy(elem)) return;

		//	get parent
			GeometricObject* pParent = m_spMGDD->get_parent(elem);

		//	forward to implementation
			getImpl().restrict_values(elem, pParent);
		}

	///	grid callback invoked on object creation
	///	\{
		virtual void vertex_created(Grid* grid, VertexBase* vrt, GeometricObject* pParent = NULL, bool replacesParent = false) {obj_created<VertexBase>(vrt, pParent, replacesParent);}
		virtual void edge_created(Grid* grid, EdgeBase* e, GeometricObject* pParent = NULL, bool replacesParent = false) {obj_created<EdgeBase>(e, pParent, replacesParent);}
		virtual void face_created(Grid* grid, Face* f, GeometricObject* pParent = NULL, bool replacesParent = false) {obj_created<Face>(f, pParent, replacesParent);}
		virtual void volume_created(Grid* grid, Volume* vol, GeometricObject* pParent = NULL, bool replacesParent = false) {obj_created<Volume>(vol, pParent, replacesParent);}
	///	\}

	///	grid callbacks invoked on erase
	///	\{
		virtual void vertex_to_be_erased(Grid* grid, VertexBase* vrt, VertexBase* replacedBy = NULL) {obj_to_be_erased<VertexBase>(vrt, replacedBy);}
		virtual void edge_to_be_erased(Grid* grid, EdgeBase* e, EdgeBase* replacedBy = NULL) {obj_to_be_erased<EdgeBase>(e, replacedBy);}
		virtual void face_to_be_erased(Grid* grid, Face* f, Face* replacedBy = NULL) {obj_to_be_erased<Face>(f, replacedBy);}
		virtual void volume_to_be_erased(Grid* grid, Volume* vol, Volume* replacedBy = NULL) {obj_to_be_erased<Volume>(vol, replacedBy);}
	///	\}
};

} // end namespace ug

#endif /* LOCAL_TRANSFER_INTERFACE_H_ */
