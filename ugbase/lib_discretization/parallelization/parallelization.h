/*
 * parallelization.h
 *
 *  Created on: 26.5.2010
 *      Author: A. Vogel, S.Reiter
 */

#ifndef __H__LIB_DISCRETIZATION__PARALLELIZATION__PARALLELIZATION__
#define __H__LIB_DISCRETIZATION__PARALLELIZATION__PARALLELIZATION__

#include "parallel_dof_manager.h"
#include "parallel_index_layout.h"
#include "parallelization_util.h"
#include "pcl/pcl.h"

namespace ug
{
//TODO: Move the following methods to an appropriate place
template <class TVector>
class ComPol_VecCopy : public pcl::ICommunicationPolicy<IndexLayout>
{
	public:
		ComPol_VecCopy() : m_pVec(NULL)	{}
		ComPol_VecCopy(TVector* pVec) : m_pVec(pVec)	{}
		
		void set_vector(TVector* pVec)	{m_pVec = pVec;}
		
		virtual bool
		collect(std::ostream& buff, Interface& interface)
		{
			typename TVector::local_vector_type u(1);
			typename TVector::local_index_type i(1);
			TVector& v = *m_pVec;
			
			for(typename Interface::iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				i[0][0] = interface.get_element(iter);
				v.get(u, i);
				buff.write((char*)&u[0], sizeof(typename TVector::
												local_vector_type::
												value_type));
			}
			return true;
		}
		
		virtual bool
		extract(std::istream& buff, Interface& interface)
		{
			typename TVector::local_vector_type u(1);
			typename TVector::local_index_type i(1);
			TVector& v = *m_pVec;
			
			for(typename Interface::iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				i[0][0] = interface.get_element(iter);
				buff.read((char*)&u[0], sizeof(typename TVector::
												local_vector_type::
												value_type));
				v.set(u, i);
			}
			return true;
		}
		
	private:
		TVector* m_pVec;
};

//TODO: Move the following methods to an appropriate place
template <class TVector>
class ComPol_VecAdd : public pcl::ICommunicationPolicy<IndexLayout>
{
	public:
		ComPol_VecAdd() : m_pVec(NULL)	{}
		ComPol_VecAdd(TVector* pVec) : m_pVec(pVec)	{}
		
		void set_vector(TVector* pVec)	{m_pVec = pVec;}
		
		virtual bool
		collect(std::ostream& buff, Interface& interface)
		{
			typename TVector::local_vector_type u(1);
			typename TVector::local_index_type i(1);
			TVector& v = *m_pVec;
			
			for(typename Interface::iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				i[0][0] = interface.get_element(iter);
				v.get(u, i);
				buff.write((char*)&u[0], sizeof(typename TVector::
												local_vector_type::
												value_type));
			}
			return true;
		}
		
		virtual bool
		extract(std::istream& buff, Interface& interface)
		{
			typename TVector::local_vector_type u(1);
			typename TVector::local_index_type i(1);
			TVector& v = *m_pVec;
			
			for(typename Interface::iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				i[0][0] = interface.get_element(iter);
				buff.read((char*)&u[0], sizeof(typename TVector::
												local_vector_type::
												value_type));
				v.add(u, i);
			}
			return true;
		}
		
	private:
		TVector* m_pVec;
};

}

#endif
