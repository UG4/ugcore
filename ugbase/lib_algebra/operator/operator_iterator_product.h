/*
 * mi.h
 *
 *  Created on: 04.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__PRODUCT__
#define __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__PRODUCT__

#include "lib_algebra/operator/operator.h"
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

namespace ug{

/** This operator is a product of ILinearIterator (multiplicative composition).*/
template <typename X, typename Y>
class LinearIteratorProduct : public ILinearIterator<X,Y>
{
	public:
	// 	Domain space
		typedef X domain_function_type;

	// 	Range space
		typedef Y codomain_function_type;

		typedef ILinearIterator<X,Y> iterator_type;

	public:

		LinearIteratorProduct() {};

	//	Name of Iterator
		virtual const char* name() const {return "IteratorProduct";}

	// 	Prepare for Operator J(u) and linearization point u (current solution)
		virtual bool init(ILinearOperator<Y,X>& J, const Y& u) {

			const size_t max = m_vIterator.size();
			for(size_t i=0; i<max; i++){
					m_vIterator[i]->init(J,u);
			}
			return true;

		}

	//	Prepare for Linear Operartor L
		virtual bool init(ILinearOperator<Y,X>& L)
		{

			const size_t max = m_vIterator.size();
			for(size_t i=0; i<max; i++) {
				m_vIterator[i]->init(L);
			}
			return true;

		}

	//	Compute new correction c = B*d
		virtual bool apply(Y& c, const X& d)
		{
			// create temporary defect and forward request
			X dTmp;
			dTmp.resize(d.size()); dTmp= d; // maybe resize needed ?

		   return apply_update_defect(c, dTmp);
		}

	//	Compute new correction c = B*d and return new defect d := d - A*c
		virtual bool apply_update_defect(Y& c, X& d)
		{
			Y cTmp;
			cTmp.resize(c.size());
			cTmp = c; // TODO: obsolete, replace by clone

			const size_t max = m_vIterator.size();
			for(size_t i=0; i<max; i++)
			{
				m_vIterator[i]->apply_update_defect(cTmp,d);
				c+=cTmp;
			}
			return true;
		}

		void add_iterator(iterator_type &I) {m_vIterator.push_back(&I);}

	//	Clone
		virtual ILinearIterator<X,Y>* clone() {UG_ASSERT(0,"Implement!"); return 0;}

	// 	Destructor
		virtual ~LinearIteratorProduct() {};

	private:

		std::vector<iterator_type*> m_vIterator;
};

/** This operator is a sum of ILinearIterator (additive composition). */
template <typename X, typename Y>
class LinearIteratorSum: public LinearIteratorProduct<X,Y>
{
	public:
	// 	Domain space
	//	typedef LinearIteratorProduct::X X;

	// 	Range space
///		typedef Y LinearIteratorProduct::codomain_function_type;

		//typedef ILinearIterator<X,Y> iterator_type;

	public:

		LinearIteratorSum() {};
		virtual const char* name() const {return "IteratorSum";}

		//	Compute new correction c = B*d
		virtual bool apply(Y& c, const X& d)
		{
			// create temporary defect and forward request
			Y cTmp;
			cTmp.resize(c.size());
			cTmp = c; // TODO: obsolete, replace by clone

			const size_t max = LinearIteratorProduct<X,Y>::m_vIterator.size();
			for(size_t i=0; i<max; i++)
			{
				LinearIteratorProduct<X,Y>::m_vIterator[i]->apply_update_defect(cTmp,d);
					c+=cTmp;
			}

			return true;
		}

		//	Compute new correction c = B*d and return new defect d := d - A*c
		virtual bool apply_update_defect(Y& c, X& d)
		{
			UG_ASSERT(0, name() << "Please implement me!");
			return false;
		}
};


} // end namespace ug

#endif
