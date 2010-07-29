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
			return interface.size() * sizeof(typename TVector::entry_type);
		}

		virtual bool
		collect(std::ostream& buff, Interface& interface)
		{
			TVector& v = *m_pVec;

			for(typename Interface::iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				const size_t index = interface.get_element(iter);
				typename TVector::entry_type entry = v[index];
				buff.write((char*)&entry, sizeof(typename TVector::entry_type));
			}
			return true;
		}

		virtual bool
		extract(std::istream& buff, Interface& interface)
		{
			typename TVector::entry_type entry;
			TVector& v = *m_pVec;

			for(typename Interface::iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				const size_t index = interface.get_element(iter);
				buff.read((char*)&entry, sizeof(typename TVector::entry_type));
				v[index] = entry;
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
			return interface.size() * sizeof(typename TVector::entry_type);
		}

		virtual bool
		collect(std::ostream& buff, Interface& interface)
		{
			TVector& v = *m_pVec;

			for(typename Interface::iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				const size_t index = interface.get_element(iter);
				typename TVector::entry_type entry = v[index];
				buff.write((char*)&entry, sizeof(typename TVector::entry_type));
			}
			return true;
		}

		virtual bool
		extract(std::istream& buff, Interface& interface)
		{
			typename TVector::entry_type entry;
			TVector& v = *m_pVec;

			for(typename Interface::iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				const size_t index = interface.get_element(iter);
				buff.read((char*)&entry, sizeof(typename TVector::entry_type));
				v[index] += entry;
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
			return interface.size() * sizeof(typename TVector::entry_type);
		}

		virtual bool
		collect(std::ostream& buff, Interface& interface)
		{
			TVector& v = *m_pVec;

			for(typename Interface::iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				const size_t index = interface.get_element(iter);
				typename TVector::entry_type entry = v[index];
				buff.write((char*)&entry, sizeof(typename TVector::entry_type));
				v[index] = 0.0;
			}
			return true;
		}

		virtual bool
		extract(std::istream& buff, Interface& interface)
		{
			typename TVector::entry_type entry;
			TVector& v = *m_pVec;

			for(typename Interface::iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				const size_t index = interface.get_element(iter);
				buff.read((char*)&entry, sizeof(typename TVector::entry_type));
				v[index] += entry;
			}
			return true;
		}

	private:
		TVector* m_pVec;
};


}//	end of namespace


#endif /* __H__LIB_ALGEBRA__PARALLELIZATION__COMMUNICATION_POLICIES__ */
