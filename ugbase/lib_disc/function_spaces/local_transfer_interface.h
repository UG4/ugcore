/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__FUNCTION_SPACS__LOCAL_TRANSFER_INTERFACE__
#define __H__UG__LIB_DISC__FUNCTION_SPACS__LOCAL_TRANSFER_INTERFACE__

#include "lib_disc/domain.h"
#include "lib_disc/local_finite_element/local_finite_element_id.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// 	Value Access Base class
////////////////////////////////////////////////////////////////////////////////

class TransferValueAccessor
{
	public:
		virtual void access_inner(GridObject* elem) = 0;

		virtual void access_closure(GridObject* elem) = 0;

		virtual ~TransferValueAccessor() = default;

		const number& operator [] (size_t i) const {
			UG_ASSERT(i < m_Val.size(), "Wrong index "<<i<<" (size: "<<m_Val.size()<<")");
			return *m_Val[i];
		}

		number& operator [] (size_t i) {
			UG_ASSERT(i < m_Val.size(), "Wrong index "<<i<<" (size: "<<m_Val.size()<<")");
			return *m_Val[i];
		}

		size_t size() const {return m_Val.size();}

	protected:
		std::vector<number*> m_Val;
};

////////////////////////////////////////////////////////////////////////////////
// 	Element Prolongation
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain>
class IElemProlongation
{
	public:
		virtual void init(ConstSmartPtr<TDomain> spDomain,
		                  SmartPtr<TransferValueAccessor> vValueChild,
		          		  SmartPtr<TransferValueAccessor> vValueParent)
		{
			m_spDomain = spDomain;
			m_spGrid = m_spDomain->grid();
			m_vValueChild = vValueChild;
			m_vValueParent = vValueParent;
		}

		virtual bool perform_prolongation_on(GridBaseObjectId gbo) = 0;

		virtual void prolongate(Vertex* parent) = 0;
		virtual void prolongate(Edge* parent) = 0;
		virtual void prolongate(Face* parent) = 0;
		virtual void prolongate(Volume* parent) = 0;

		IElemProlongation() = default;
		virtual ~IElemProlongation() = default;

	protected:
		ConstSmartPtr<TDomain> m_spDomain;
		ConstSmartPtr<MultiGrid> m_spGrid;

        SmartPtr<TransferValueAccessor> m_vValueChild;
        SmartPtr<TransferValueAccessor> m_vValueParent;
};


template <typename TDomain, typename TImpl>
class ElemProlongationBase : public IElemProlongation<TDomain>
{
	public:
		ElemProlongationBase() = default;
		~ElemProlongationBase() override = default;

		void prolongate(Vertex* parent) override {
			this->getImpl().prolongate(parent, *this->m_vValueChild, *this->m_vValueParent);
		}

		void prolongate(Edge* parent) override {
			this->getImpl().prolongate(parent, *this->m_vValueChild, *this->m_vValueParent);
		}

		void prolongate(Face* parent) override {
			this->getImpl().prolongate(parent, *this->m_vValueChild, *this->m_vValueParent);
		}

		void prolongate(Volume* parent) override {
			this->getImpl().prolongate(parent, *this->m_vValueChild, *this->m_vValueParent);
		}

	protected:
	///	access to implementation
		TImpl& getImpl() {return static_cast<TImpl&>(*this);}

	///	const access to implementation
		const TImpl& getImpl() const {return static_cast<const TImpl&>(*this);}
};


template <typename TDomain>
SmartPtr<IElemProlongation<TDomain> >
GetStandardElementProlongation(const LFEID& lfeid);


////////////////////////////////////////////////////////////////////////////////
// 	Element Restriction
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain>
class IElemRestriction
{
	public:
		virtual void init(ConstSmartPtr<TDomain> spDomain,
		                  SmartPtr<TransferValueAccessor> vValueChild,
		          		  SmartPtr<TransferValueAccessor> vValueParent)
		{
			m_spDomain = spDomain;
			m_spGrid = m_spDomain->grid();
			m_vValueChild = vValueChild;
			m_vValueParent = vValueParent;
		}

		virtual bool perform_restriction_on(GridBaseObjectId gbo) = 0;

		virtual void do_restrict(Vertex* parent) = 0;
		virtual void do_restrict(Edge* parent) = 0;
		virtual void do_restrict(Face* parent) = 0;
		virtual void do_restrict(Volume* parent) = 0;

		virtual ~IElemRestriction() {}

	protected:
		ConstSmartPtr<TDomain> m_spDomain;
		ConstSmartPtr<MultiGrid> m_spGrid;

        SmartPtr<TransferValueAccessor> m_vValueChild;
        SmartPtr<TransferValueAccessor> m_vValueParent;
};


template <typename TDomain, typename TImpl>
class ElemRestrictionBase : public IElemRestriction<TDomain>
{
	public:
		void do_restrict(Vertex* parent) override {
			this->getImpl().do_restrict(parent, *this->m_vValueChild, *this->m_vValueParent);
		}

		void do_restrict(Edge* parent) override {
			this->getImpl().do_restrict(parent, *this->m_vValueChild, *this->m_vValueParent);
		}

		void do_restrict(Face* parent) override {
			this->getImpl().do_restrict(parent, *this->m_vValueChild, *this->m_vValueParent);
		}

		void do_restrict(Volume* parent) override {
			this->getImpl().do_restrict(parent, *this->m_vValueChild, *this->m_vValueParent);
		}

	protected:
	///	access to implementation
		TImpl& getImpl() {return static_cast<TImpl&>(*this);}

	///	const access to implementation
		const TImpl& getImpl() const {return static_cast<const TImpl&>(*this);}
};


template <typename TDomain>
SmartPtr<IElemRestriction<TDomain> >
GetStandardElementRestriction(const LFEID& lfeid);


} // end namespace ug

#endif