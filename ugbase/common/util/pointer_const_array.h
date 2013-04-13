// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 06.08.2012 (m,d,y)

#ifndef __H__UG__pointer_const_array__
#define __H__UG__pointer_const_array__

#include <vector>
#include <common/assert.h>
#include <common/types.h>

namespace ug
{

/// \addtogroup ugbase_common_util
/// \{

///	Container which holds an array of pointers
/**	Gives access to an array of elements of type TPtr. TPtr is assumed to be a
 * pointer type and treated as thus. This means it has to support assignment
 * of NULL. Furthermore no destructors are called when an array is cleared or
 * shrinks. Operators +(int) and ++ also have to be supported.
 *
 * The array provided by PointerArray is either managed by the PointerArray class
 * itself or can be set from outside.
 *
 * If an external array is set from outside, the contents will not be copied, unless
 * indicated otherwise. Instead only the pointer is stored and used.
 * The class will not take responsibility over the lifetime of that array.
 * That means that the user is responsible to make sure that the external array,
 * to which the class points to, is valid when accessed.
 *
 * If the user however uses the methods reserve or push_back, then
 * dedicated memory is allocated and managed by the class itself. This memory
 * is freed when no longer needed.
 *
 * If an array was set from outside and push_back is called, then the contents
 * of the external array are copied to the local memory of the class.
 *
 * Note that no memory is freed when clear is called. To make sure that
 * all internal memory is freed, you have to assign an empty instance to the
 * class.
 *
 * Note that entries of this container may not be changed once they have been added!
 *
 * If in doubt whether this class should be used, or if a std::vector<TPtr> would
 * be more appropriate, you should most likely choose the std::vector.
 * The main purpose of this class is to allow optional access to an array existing
 * outside of the class itself.
 *
 * \todo:	Add iterators
 * \todo:	Implement analogous PointerArray with additional resize method, which
 * 			allows to change existing entries.
 */
template <class TPtr>
class PointerConstArray{
	typedef TPtr const*	ConstPtrArray;
	typedef TPtr*		PtrArray;

	public:
		PointerConstArray();
		PointerConstArray(const PointerConstArray& pa);
		~PointerConstArray();

		PointerConstArray& operator=(const PointerConstArray& pa);

	///	returns the size of the associated array.
		inline size_t size() const;

	///	returns true if the associated array is empty
		inline bool empty() const;

	///	returns the i-th entry of the array. Make sure that i < size().
	/**	Note that this is a read only operation and that entries in the associated
	 * array may not be changed once they were added.*/
		inline TPtr const operator[](size_t i) const;

	///	set the array on which the container operates.
	/**	The container will not take ownership over the array until bCopy is set
	 * to true. The caller is thus responsible to free the array if necessary.
	 * \note	entries in the external entry will not be changed at any time by
	 * 			this class or through use of this class.
	 * \note	Setting an external array does not affect the capacity of this
	 * 			container, unless bCopy is set to true.*/
		void set_external_array(TPtr const* array, size_t size, bool bCopy = false);

	///	reserves memory but does not alter the size
	/**	Note that the method does not free memory if the new capacity is lower
	 * than the current capacity.*/
		void reserve(size_t capacity);

	///	clears the container
	/**	The methods sets the internal size to 0, but does not change the capacity.
	 * If you want to free the memory of a container, either delete it or assign
	 * an empty one.*/
		void clear();

	///	appends the given pointer at the end of the array.
	/**	If the capacity is too small, it is doubled before the element is appended.*/
		void push_back(TPtr p);

	private:
	///	copies contents from the given PointerConstArray
	/**	The method assumes, that all internal memory has been freed prior to
	 * calling this method!*/
		void assign_pointer_const_array(const PointerConstArray& pa);

	///	reserves the specified memory and optionally copies old data
		void reserve(size_t capacity, bool copyOldData);

	private:
		ConstPtrArray		m_array;
		PtrArray			m_data;
		size_t				m_size;
		size_t				m_capacity;
};

// end group ugbase_common_util
/// \}

}//	end of namespace

////////////////////////////////////////
//	include implementation
#include "pointer_const_array_impl.hpp"

#endif
