/*
 * local_transfer_interface.h
 *
 *  Created on: 07.03.2012
 *      Author: andreasvogel
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
		virtual void access_inner(GeometricObject* elem) = 0;

		virtual void access_closure(GeometricObject* elem) = 0;

		virtual ~TransferValueAccessor() {}

		const number& operator[](size_t i) const {
			UG_ASSERT(i < m_Val.size(), "Wront index "<<i<<" (size: "<<m_Val.size()<<")");
			return *m_Val[i];
		}

		number& operator[](size_t i) {
			UG_ASSERT(i < m_Val.size(), "Wront index "<<i<<" (size: "<<m_Val.size()<<")");
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

		virtual bool perform_prolongation_on(GeometricBaseObject gbo) = 0;

		virtual void prolongate(VertexBase* parent) = 0;
		virtual void prolongate(EdgeBase* parent) = 0;
		virtual void prolongate(Face* parent) = 0;
		virtual void prolongate(Volume* parent) = 0;

		virtual ~IElemProlongation() {}

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
		virtual void prolongate(VertexBase* parent){
			this->getImpl().prolongate(parent, *this->m_vValueChild, *this->m_vValueParent);
		}
		virtual void prolongate(EdgeBase* parent){
			this->getImpl().prolongate(parent, *this->m_vValueChild, *this->m_vValueParent);
		}
		virtual void prolongate(Face* parent){
			this->getImpl().prolongate(parent, *this->m_vValueChild, *this->m_vValueParent);
		}
		virtual void prolongate(Volume* parent){
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

		virtual bool perform_restriction_on(GeometricBaseObject gbo) = 0;

		virtual void do_restrict(VertexBase* parent) = 0;
		virtual void do_restrict(EdgeBase* parent) = 0;
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
		virtual void do_restrict(VertexBase* parent){
			this->getImpl().do_restrict(parent, *this->m_vValueChild, *this->m_vValueParent);
		}
		virtual void do_restrict(EdgeBase* parent){
			this->getImpl().do_restrict(parent, *this->m_vValueChild, *this->m_vValueParent);
		}
		virtual void do_restrict(Face* parent){
			this->getImpl().do_restrict(parent, *this->m_vValueChild, *this->m_vValueParent);
		}
		virtual void do_restrict(Volume* parent){
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

#endif /* __H__UG__LIB_DISC__FUNCTION_SPACS__LOCAL_TRANSFER_INTERFACE__ */
