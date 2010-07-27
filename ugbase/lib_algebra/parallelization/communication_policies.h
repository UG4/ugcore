/*
 * communication_policies.h
 *
 *  Created on: 14.6.2010
 *      Author: A. Vogel, S.Reiter
 */

#ifndef __H__LIB_ALGEBRA__PARALLELIZATION__COMMUNICATION_POLICIES__
#define __H__LIB_ALGEBRA__PARALLELIZATION__COMMUNICATION_POLICIES__

#include "parallel_index_layout.h"
#include "pcl/pcl.h"


namespace ug{

template <class TVector>
class ComPol_VecCopy : public pcl::ICommunicationPolicy<IndexLayout>
{
	public:
		ComPol_VecCopy() : m_pVec(NULL)	{}
		ComPol_VecCopy(TVector* pVec) : m_pVec(pVec)	{}

		void set_vector(TVector* pVec)	{m_pVec = pVec;}

		virtual int
		get_required_buffer_size(Interface& interface)
		{
			return interface.size() * sizeof(typename TVector::
												local_vector_type::
												entry_type);
		}

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
												entry_type));
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
												entry_type));
				v.set(u, i);
			}
			return true;
		}

	private:
		TVector* m_pVec;
};

template <class TVector>
class ComPol_VecAdd : public pcl::ICommunicationPolicy<IndexLayout>
{
	public:
		ComPol_VecAdd() : m_pVec(NULL)	{}
		ComPol_VecAdd(TVector* pVec) : m_pVec(pVec)	{}

		void set_vector(TVector* pVec)	{m_pVec = pVec;}

		virtual int
		get_required_buffer_size(Interface& interface)
		{
			return interface.size() * sizeof(typename TVector::
												local_vector_type::
												entry_type);
		}

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
												entry_type));
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
												entry_type));
				v.add(u, i);
			}
			return true;
		}

	private:
		TVector* m_pVec;
};

template <class TVector>
class ComPol_VecAddSetZero : public pcl::ICommunicationPolicy<IndexLayout>
{
	public:
		ComPol_VecAddSetZero() : m_pVec(NULL)	{}
		ComPol_VecAddSetZero(TVector* pVec) : m_pVec(pVec)	{}

		void set_vector(TVector* pVec)	{m_pVec = pVec;}

		virtual int
		get_required_buffer_size(Interface& interface)
		{
			return interface.size() * sizeof(typename TVector::
												local_vector_type::
												entry_type);
		}

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
												entry_type));
				u[0] = 0.0;
				v.set(u, i);
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
												entry_type));
				v.add(u, i);
			}
			return true;
		}

	private:
		TVector* m_pVec;
};


}//	end of namespace


#endif /* __H__LIB_ALGEBRA__PARALLELIZATION__COMMUNICATION_POLICIES__ */
