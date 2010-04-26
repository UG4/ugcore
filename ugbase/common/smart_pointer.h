//	smart_pointer.h
//	created by Sebastian Reiter

#ifndef __SMART_POINTER__
#define __SMART_POINTER__

#include <functional>

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	Policies

////////////////////////////////////////////////////////////////////////
//	FreeDelete
template <typename T>
class FreeDelete
{
	public:
		void free(T* data)	{delete data;}
};
////////////////////////////////////////////////////////////////////////
//	FreeRelease
template <typename T>
class FreeRelease
{
	public:
		void free(T* data)	{data->Release;}
};


////////////////////////////////////////////////////////////////////////
//	SmartPtr
/**
*	The FreePolicy has to feature the method free().
*/
template <typename T, class FreePolicy = FreeDelete<T> >
class SmartPtr : public FreePolicy
{
	public:
		SmartPtr() : m_ptr(0), m_refCount(0)	{}
		SmartPtr(T* ptr) : m_ptr(ptr), m_refCount(0)	{if(ptr) m_refCount = new int(1);}
		SmartPtr(const SmartPtr<T>& sp) : m_ptr(sp.m_ptr), m_refCount(sp.m_refCount)
		{
			if(m_refCount) (*m_refCount)++;
		}

	/**	this template method allows to assign smart-pointers that encapsulate
	 *	derivates of T. Make sure that the pointer-type of TSmartPtr is castable
	 *	to T*.*/
		template <class TPtr>
		SmartPtr(const SmartPtr<TPtr>& sp) :
			m_ptr(static_cast<T*>(sp.get_impl())),
			m_refCount(sp.get_refcount_ptr())
		{
			if(m_refCount) (*m_refCount)++;
		}

		~SmartPtr() {release();}

		T* operator->() const {return m_ptr;}
		T& operator*() const {return *m_ptr;}
		SmartPtr<T>& operator=(const SmartPtr& sp)	{
			if(m_ptr)
				release();
			m_ptr = sp.m_ptr;
			m_refCount = sp.m_refCount;
			if(m_refCount)
				(*m_refCount)++;
			return *this;
		}

		T* get_impl() const	{return m_ptr;}
		int get_refcount() const {if(m_refCount) return *m_refCount; return 0;}
		
	///	WARNING: this method is dangerous!
	/**	This method should never be used since it may be removed in future
	 *	versions of the SmartPtr class.
	 *	It is featured in order to allow to implement a template-constructor
	 *	that casts element-pointers of a smart pointer.*/
		int* get_refcount_ptr() const {return m_refCount;}

	private:
		void release() {
			if(m_refCount)
			{
				(*m_refCount)--;
				if((*m_refCount) < 1)
				{
					delete m_refCount;
					//delete m_ptr;
					free(m_ptr);
				}
			}
		}

		T*		m_ptr;
		int*	m_refCount;
};

namespace std
{
	template <class T>
	struct less<SmartPtr<T> > : public binary_function<SmartPtr<T>, SmartPtr<T>, bool>
	{
		bool operator()(const SmartPtr<T>& lhs, const SmartPtr<T>& rhs) const
		{
			return less<T*>()(lhs.get_impl(), rhs.get_impl());
		}
	};
}

#endif