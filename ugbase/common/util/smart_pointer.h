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
		static void free(T* data)	{delete data;}
};

////////////////////////////////////////////////////////////////////////
//	FreeRelease
template <typename T>
class FreeRelease
{
	public:
		static void free(T* data)	{data->Release;}
};


////////////////////////////////////////////////////////////////////////
//	SmartPtr
/**
 * A Smart Pointer stores a pointer encapsulates a pointer. With each copy
 * of the SmartPointer a shared reference count is increased and with
 * each destructor call to a SmartPointer this shared reference count is
 * decreased. If the shared reference count reaches 0 free as specified in
 * the FreePolicy is executed on the pointer.
 *
 *	The FreePolicy has to feature the method free().
 */
template <typename T, class FreePolicy = FreeDelete<T> >
class SmartPtr
{
	friend class SmartPtr<void>;

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

	///	copies a SmartPtr and applies a reinterpret_cast on its impl.
		template <class TIn>
		static SmartPtr<T> from_reinterpret_cast(SmartPtr<TIn>& in)	{
			return SmartPtr<T>(reinterpret_cast<T*>(in.m_ptr),
								in.get_refcount_ptr());
		}

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

	///	returns true if the pointer is valid, false if not.
		inline bool is_valid()		{return m_ptr != NULL;}
		
	///	WARNING: this method is dangerous!
	/**	This method should never be used since it may be removed in future
	 *	versions of the SmartPtr class.
	 *	It is featured in order to allow to implement a template-constructor
	 *	that casts element-pointers of a smart pointer.*/
		int* get_refcount_ptr() const {return m_refCount;}

	private:
	///	decrements the refCount and frees the encapsulated pointer if required.
		void release() {
			if(m_refCount)
			{
				(*m_refCount)--;
				if((*m_refCount) < 1)
				{
					delete m_refCount;
					//delete m_ptr;
					FreePolicy::free(m_ptr);
				}
			}
		}

	////////////////////////////////
	//	The following methods are required for SmartPtr<void>
	///	you should only use this constructor if you really know what you're doing!
		SmartPtr(T* ptr, int* refCount) : m_ptr(ptr), m_refCount(refCount)
		{
			if(m_refCount)
				(*m_refCount)++;
		}

	//	this release method is required by SmartPtr<void>
		static void free_void_ptr(void* ptr){
			FreePolicy::free(reinterpret_cast<T*>(ptr));
		}

	private:
		T*		m_ptr;
		int*	m_refCount;
};

/**	The SmartPtr<void> is a specialization of the SmartPtr class.
 * It can only be constructed from an existing SmartPtr or from an
 * existing SmartPtr<void>. This is crucial to guarantee save
 * release methods.
 *
 * In contrary to the original SmartPtr class, the void specialization
 * does not feature the -> and * operators. One also can't directly access
 * the encapsulated pointer.
 *
 * \todo	add to_smart_ptr_dynamic
 */
template <>
class SmartPtr<void>
{
	public:
		SmartPtr() : m_ptr(0), m_refCountPtr(0), m_freeFunc(0) {}

		SmartPtr(const SmartPtr<void>& sp) :
			m_ptr(sp.m_ptr),
			m_refCountPtr(sp.m_refCountPtr),
			m_freeFunc(sp.m_freeFunc)
		{
			if(m_refCountPtr) (*m_refCountPtr)++;
		}

		template <class T>
		SmartPtr(SmartPtr<T>& sp) :
			m_ptr((void*)sp.get_impl()),
			m_refCountPtr(sp.get_refcount_ptr()),
			m_freeFunc(&SmartPtr<T>::free_void_ptr)
		{
			if(m_refCountPtr) (*m_refCountPtr)++;
		}

		SmartPtr<void>& operator=(const SmartPtr<void>& sp)
		{
			if(m_ptr)
				release();
			m_ptr = sp.m_ptr;
			m_refCountPtr = sp.m_refCountPtr;
			if(m_refCountPtr)
				(*m_refCountPtr)++;
			m_freeFunc = sp.m_freeFunc;
			return *this;
		}

		template <class T>
		SmartPtr<void>& operator=(const SmartPtr<T>& sp)
		{
			if(m_ptr)
				release();
			m_ptr = sp.m_ptr;
			m_refCountPtr = sp.m_refCountPtr;
			if(m_refCountPtr)
				(*m_refCountPtr)++;
			m_freeFunc = &SmartPtr<T>::free_void_ptr;
			return *this;
		}

	///	Returns a SmartPtr with the specified type and shared reference counting.
	/**	USE WITH CARE! ONLY COMPATIBLE TYPES SHOULD BE USED*/
		template <class T>
		SmartPtr<T> to_smart_pointer_reinterpret(){
			return SmartPtr<T>(reinterpret_cast<T*>(m_ptr), m_refCountPtr);
		}

		inline bool is_valid()		{return m_ptr != NULL;}

	private:
		void release() {
			if(m_refCountPtr)
			{
				(*m_refCountPtr)--;
				if((*m_refCountPtr) < 1)
				{
					delete m_refCountPtr;
					m_freeFunc(m_ptr);
				}
			}
		}

		void* m_ptr;
		int* m_refCountPtr;
		void (*m_freeFunc)(void*);
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
