/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Arne Nägel
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

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PRODUCT__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PRODUCT__

#include "lib_algebra/operator/interface/linear_iterator.h"
#include "common/util/smart_pointer.h"
#include <vector>

#ifdef UG_PARALLEL
#include "lib_algebra/parallelization/parallelization.h"
#endif

namespace ug{

/** Base class for ILinearIterators build from other ILinearIterators */
template <typename X, typename Y>
class CombinedLinearIterator : public ILinearIterator<X,Y>
{
	public:
		CombinedLinearIterator() {};

		CombinedLinearIterator(const std::vector<SmartPtr<ILinearIterator<X,Y> > >& vIterator)
			: m_vIterator(vIterator)
		{};

		//	Name of Iterator
		virtual const char* name() const = 0;

		///	returns if parallel solving is supported
		virtual bool supports_parallel() const
		{
			for(size_t i = 0; i < m_vIterator.size(); ++i)
				if(!m_vIterator[i]->supports_parallel())
					return false;
			return true;
		}

		// 	Prepare for Operator J(u) and linearization point u (current solution)
		virtual bool init(SmartPtr<ILinearOperator<Y,X> > J, const Y& u) = 0;

		// Prepare for Linear Operator L
		virtual bool init(SmartPtr<ILinearOperator<Y,X> > L) = 0;

		// Compute correction
		virtual bool apply(Y& c, const X& d) = 0;

		// Compute correction and update defect
		virtual bool apply_update_defect(Y& c, X& d) = 0;

		void add_iterator(SmartPtr<ILinearIterator<X,Y> > I) {m_vIterator.push_back(I);}

		void add_iterator(SmartPtr<ILinearIterator<X,Y> > I,size_t nr) { for (size_t i=0;i<nr;i++) m_vIterator.push_back(I);}

		//	Clone
		virtual SmartPtr<ILinearIterator<X,Y> > clone() = 0;

	protected:
		std::vector<SmartPtr<ILinearIterator<X,Y> > > m_vIterator;
};


/** This operator is a product of ILinearIterator (multiplicative composition).*/
template <typename X, typename Y>
class LinearIteratorProduct : public CombinedLinearIterator<X,Y>
{
	protected:
		typedef CombinedLinearIterator<X,Y> base_type;
		using base_type::m_vIterator;

	public:
		LinearIteratorProduct() {};

		LinearIteratorProduct(const std::vector<SmartPtr<ILinearIterator<X,Y> > >& vIterator)
			: CombinedLinearIterator<X,Y>(vIterator)
		{};

		//	Name of Iterator
		virtual const char* name() const {return "IteratorProduct";}

		// 	Prepare for Operator J(u) and linearization point u (current solution)
		virtual bool init(SmartPtr<ILinearOperator<Y,X> > J, const Y& u) {
			bool bRes = true;
			for(size_t i = 0; i < m_vIterator.size(); i++){
				if ((i>0) && (m_vIterator[i]==m_vIterator[i-1])) continue;
				bRes &= m_vIterator[i]->init(J,u);
			}
			return bRes;
		}

		// Prepare for Linear Operator L
		virtual bool init(SmartPtr<ILinearOperator<Y,X> > L)
		{
			bool bRes = true;
			for(size_t i = 0; i < m_vIterator.size(); i++) {
				if ((i>0) && (m_vIterator[i]==m_vIterator[i-1])) continue;
				bRes &= m_vIterator[i]->init(L);
			}
			return bRes;
		}

		// Compute correction
		virtual bool apply(Y& c, const X& d)
		{
			// create temporary defect and forward request
			SmartPtr<X> spDTmp = d.clone();
			return apply_update_defect(c, *spDTmp);
		}

		// Compute correction and update defect
		virtual bool apply_update_defect(Y& c, X& d)
		{
			SmartPtr<Y> spCTmp = c.clone_without_values();

			bool bRes = true;
			c.set(0.0);
			for(size_t i = 0; i < m_vIterator.size(); i++)
			{
				bRes &= m_vIterator[i]->apply_update_defect(*spCTmp,d);
				c += (*spCTmp);
			}
			return bRes;
		}

		//	Clone
		virtual SmartPtr<ILinearIterator<X,Y> > clone() {
			return make_sp(new LinearIteratorProduct(*this));
		}
};

/** This operator is a sum of ILinearIterator (additive composition). */
template <typename X, typename Y>
class LinearIteratorSum : public CombinedLinearIterator<X,Y>
{
	protected:
		typedef CombinedLinearIterator<X,Y> base_type;
		using base_type::m_vIterator;

	public:
		LinearIteratorSum() {};

		LinearIteratorSum(const std::vector<SmartPtr<ILinearIterator<X,Y> > >& vIterator)
			: CombinedLinearIterator<X,Y>(vIterator)
		{};

		//	Name of Iterator
		virtual const char* name() const {return "IteratorProduct";}

		// 	Prepare for Operator J(u) and linearization point u (current solution)
		virtual bool init(SmartPtr<ILinearOperator<Y,X> > J, const Y& u) {
			m_spOp = J;
			bool bRes = true;
			for(size_t i = 0; i < m_vIterator.size(); i++){
				if ((i>0) && (m_vIterator[i]==m_vIterator[i-1])) continue;
				bRes &= m_vIterator[i]->init(J,u);
			}
			return bRes;
		}

		// Prepare for Linear Operator L
		virtual bool init(SmartPtr<ILinearOperator<Y,X> > L)
		{
			m_spOp = L;
			bool bRes = true;
			for(size_t i = 0; i < m_vIterator.size(); i++) {
				if ((i>0) && (m_vIterator[i]==m_vIterator[i-1])) continue;
				bRes &= m_vIterator[i]->init(L);
			}
			return bRes;
		}

		// Compute correction
		virtual bool apply(Y& c, const X& d)
		{
			// create temporary correction
			SmartPtr<Y> spCTmp = c.clone_without_values();

			// this part computes all corrections independently
			c.set(0.0);

			bool bRes = true;
			for(size_t i = 0; i < m_vIterator.size(); i++)
			{
				bRes &= m_vIterator[i]->apply(*spCTmp,d);
				c += (*spCTmp);
			}

			return bRes;
		}

		//	Compute correction and update defect
		virtual bool apply_update_defect(Y& c, X& d)
		{
			bool bRet = apply(c, d);
			m_spOp->apply_sub(d, c);
			return bRet;
		}

		//	Clone
		virtual SmartPtr<ILinearIterator<X,Y> > clone() {
			return make_sp(new LinearIteratorSum(*this));
		}

	protected:
		SmartPtr<ILinearOperator<Y,X> > m_spOp;
};


} // end namespace ug

#endif
