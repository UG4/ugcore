// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y2010

#include <cstring>
#include "class_name_provider.h"
#include "function_traits.h"
#include "common/common.h"
#include "common/util/smart_pointer.h"
#include "stdvectorwrap.h"

#ifndef __H__UG_BRIDGE__PARAMETER_STACK__
#define __H__UG_BRIDGE__PARAMETER_STACK__


#define PUSH_PARAM_TO_STACK(paramVar, val, paramType, clName)	{m_entries[m_numEntries].param.paramVar = (val);\
																m_entries[m_numEntries].type = (paramType);\
																m_entries[m_numEntries].pClassNameNode = (clName);\
																++m_numEntries;}

#define PUSH_STD_STRING_TO_STACK(val, clName)		{m_entries[m_numEntries].param.m_stdString = new std::string(val);\
													 m_entries[m_numEntries].type = PT_STD_STRING;\
													 m_entries[m_numEntries].pClassNameNode = (clName);\
													 ++m_numEntries;}

//	call the constructor and assign the smart-ptr afterwards.
#define PUSH_SP_TO_STACK(val, clName)				{m_entries[m_numEntries].param.m_smartPtrWrapper = new SmartPtr<void>(val);\
													 m_entries[m_numEntries].type = PT_SMART_POINTER;\
													 m_entries[m_numEntries].pClassNameNode = (clName);\
													 ++m_numEntries;}

//	call the constructor and assign the smart-ptr afterwards.
#define PUSH_CSP_TO_STACK(val, clName)				{m_entries[m_numEntries].param.m_constSmartPtrWrapper = new ConstSmartPtr<void>(val);\
													 m_entries[m_numEntries].type = PT_CONST_SMART_POINTER;\
													 m_entries[m_numEntries].pClassNameNode = (clName);\
													 ++m_numEntries;}

namespace ug
{
namespace bridge
{


struct ERROR_BadIndex{
	ERROR_BadIndex(int index) : m_index(index)	{}
	int m_index;
};

struct ERROR_BadConversion{
	ERROR_BadConversion(int index, int from, int to) :
		m_index(index), m_from(from), m_to(to)	{}
	
	int	m_index;
	int m_from;
	int m_to;
};

struct ERROR_IncompatibleClasses{
	ERROR_IncompatibleClasses(int index, const std::string& from, const std::string& to) :
		m_index(index), m_from(from), m_to(to)	{}

	int	m_index;
	std::string m_from;
	std::string m_to;
};



static inline int ARRAY_INDEX_TO_STACK_INDEX(int index, int stackSize)
{
	int nIndex = index;
	if(nIndex < 0)
		nIndex = stackSize + nIndex;
	
	if(nIndex < 0 || nIndex >= stackSize)
		throw(ERROR_BadIndex(index));
	
	return nIndex;
}


const int PARAMETER_STACK_MAX_SIZE = UG_REGISTRY_MAX_NUM_ARGS;

// CAUTION:
// Type values must not be changed! Bindings rely on the exact values.
// Append new types at the end and update bindings.
// If in doubt contact binding developers!
enum ParameterTypes
{
	PT_UNKNOWN = 0,
	PT_BOOL = 1,
	PT_INTEGER = 2,
	PT_NUMBER = 3,
	PT_CSTRING = 4,
	PT_STD_STRING = 5,
	PT_POINTER = 6,
	PT_CONST_POINTER = 7,
	PT_SMART_POINTER = 8,
	PT_CONST_SMART_POINTER = 9
};

////////////////////////////////////////////////////////////////////////
///	A stack that can hold values together with their type-id.
/**
 * This class is mainly used as an intermediate parameter storage during
 * calls to ugbridge methods. Its focus is on being leightweight and fast,
 * which makes it a little unflexible at times.
 * Note that the maximal number of parameters is specified by the constant
 * PARAMETER_STACK_MAX_SIZE. This is set to 10 by default. Please note, that
 * this value should not be unnecessarily high. This wouldn't make sense anyway,
 * since the template-method-wrappers can't take any more parameters.
 *
 * Supported types are integer, number, string, reference, pointer and smart-pointer.
 * References and pointers are stored in a void*. The user is responsible to
 * associate the correct types.
 *
 * Use push_* to add new parameters to the stack - * stands for one of
 * the parameters above.
 *
 * Use to_* to retrieve a value.
 *
 * Use set_* to set a value in an existing entry.
 *
 * Indices start with 0. Negative indices can be used to start indexing from
 * the top of the stack.
 *
 * If a bad index is passed, an instance of ERROR_BadIndex is thrown.
 * 
 * If set_* or get_* are performed on wrong types, an instance of
 * ERROR_BadConversion is thrown.
 */
class ParameterStack
{
	public:
	///	default constructor
		ParameterStack() :
			m_numEntries(0), m_bHasSmartPtrs(false), m_bHasStringCopies(false)
		{}
		
	///	destructor
		~ParameterStack()
		{
		//	we have to release all string copies and smart pointers.
			if(m_bHasStringCopies){
				for(int i = 0; i < m_numEntries; ++i){
					if(m_entries[i].type == PT_STD_STRING)
						delete m_entries[i].param.m_stdString;
				}
			}

			if(m_bHasSmartPtrs){
				for(int i = 0; i < m_numEntries; ++i){
					if(m_entries[i].type == PT_SMART_POINTER)
						delete (SmartPtr<void>*)m_entries[i].param.m_smartPtrWrapper;
					else if(m_entries[i].type == PT_CONST_SMART_POINTER)
						delete (ConstSmartPtr<void>*)m_entries[i].param.m_constSmartPtrWrapper;
				}
			}
		}

	////////////////////////////////
	//	info
	////////////////////////////////
		inline int size() const		{return m_numEntries;}
		
	////////////////////////////////
	//	push
	////////////////////////////////

	///	push bool
		inline void push_bool(bool val = true)
			{PUSH_PARAM_TO_STACK(m_bool, val, PT_BOOL,
			  &ClassNameProvider<bool>::class_name_node());}

	///	push integer
		inline void push_integer(int val = 0)
			{PUSH_PARAM_TO_STACK(m_int, val, PT_INTEGER,
			  &ClassNameProvider<int>::class_name_node());}

	///	push a number (double/float)
		inline void push_number(number val = 0)
			{PUSH_PARAM_TO_STACK(m_number, val, PT_NUMBER,
			  &ClassNameProvider<number>::class_name_node());}
		
	///	strings are not bufferd.
		inline void push_cstring(const char* str = "")
		{
			PUSH_PARAM_TO_STACK(m_cstring, str, PT_CSTRING,
								&ClassNameProvider<char>::class_name_node());
		}

		inline void push_std_string(const char* str = "")
		{
		//	push a std string to the stack. this performs a copy of the string
			PUSH_STD_STRING_TO_STACK(str, &ClassNameProvider<char>::class_name_node());
			m_bHasStringCopies = true;
		}

		inline void push_std_string(const std::string& str)
		{
		//	push a std string to the stack. this performs a copy of the string
			PUSH_STD_STRING_TO_STACK(str, &ClassNameProvider<char>::class_name_node());
			m_bHasStringCopies = true;
		}

	/// push user defined classes
		template<class T>
		inline void push_pointer(T* ptr = NULL)
			{PUSH_PARAM_TO_STACK(	m_ptr, (void*)ptr, PT_POINTER,
		                     		&ClassNameProvider<T>::class_name_node());}

		inline void push_pointer(void* ptr, const ClassNameNode* classNameNode)
			{PUSH_PARAM_TO_STACK(	m_ptr, ptr, PT_POINTER, classNameNode);}

		
		template<class T>
		inline void push_smart_pointer_std_vector(const SmartPtr<std::vector<T> >& ptr = SmartPtr<std::vector<T> >(NULL))
			{PUSH_SP_TO_STACK(ptr, &ClassNameProvider<std_vector_wrap<T> >::class_name_node());
			m_bHasSmartPtrs = true;}
			
		
	/// user defined classes
		template<class T>
		inline void push_const_pointer(const T* ptr = NULL)
			{PUSH_PARAM_TO_STACK(	m_constPtr, (const void*)ptr, PT_CONST_POINTER,
			                     	&ClassNameProvider<T>::class_name_node());}

		inline void push_const_pointer(const void* ptr, const ClassNameNode* classNameNode)
			{PUSH_PARAM_TO_STACK(m_constPtr, ptr, PT_CONST_POINTER, classNameNode);}

	/// SmartPtrs to user defined classes
		template<class T>
		inline void push_smart_pointer(const SmartPtr<T>& ptr = SmartPtr<T>(NULL))
			{PUSH_SP_TO_STACK(ptr, &ClassNameProvider<T>::class_name_node());
				m_bHasSmartPtrs = true;}

		inline void push_smart_pointer(const SmartPtr<void>& ptr,
		                               const ClassNameNode* classNameNode)
			{PUSH_SP_TO_STACK(ptr, classNameNode);
				m_bHasSmartPtrs = true;}

	/// ConstSmartPtrs to user defined classes
		template<class T>
		inline void push_const_smart_pointer(const ConstSmartPtr<T>& ptr = ConstSmartPtr<T>(NULL))
			{PUSH_CSP_TO_STACK(ptr, &ClassNameProvider<T>::class_name_node());
				m_bHasSmartPtrs = true;}

		inline void push_const_smart_pointer(const ConstSmartPtr<void>& ptr,
		                                     const ClassNameNode* classNameNode)
			{PUSH_CSP_TO_STACK(ptr, classNameNode);
				m_bHasSmartPtrs = true;}

	////////////////////////////////
	//	get
	////////////////////////////////

	///	returns ParameterType enum of data type for a stack entry
		uint get_type(int index) const
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);
			return m_entries[index].type;
		}
		
	//\todo: return const std::string&
	///	returns the class name for an element in the param stack
		const char* class_name(int index) const
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);
			if(m_entries[index].pClassNameNode == NULL)
				throw(UGError("ClassNameNode missing in Parameter stack."));
			return (*m_entries[index].pClassNameNode).name().c_str();
		}

	///	returns the class name node for an element in the param stack
		const ClassNameNode* class_name_node(int index) const
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);
			if(m_entries[index].pClassNameNode == NULL)
				throw(UGError("ClassNameNode missing in Parameter stack."));
			return m_entries[index].pClassNameNode;
		}

	///	return element in param stack casted to bool
		bool to_bool(int index) const
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);

			const Entry& e = m_entries[index];
			if(e.type == PT_BOOL)
				return e.param.m_bool;
			else if(e.type == PT_INTEGER)
				return (bool) e.param.m_int;

			throw(ERROR_BadConversion(index, e.type, PT_BOOL));
		}

	///	return element in param stack casted to int
		int to_integer(int index) const
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);
			
			const Entry& e = m_entries[index];
			if(e.type == PT_INTEGER)
				return e.param.m_int;
			else if(e.type == PT_NUMBER)
				return (int) e.param.m_number;
			
			throw(ERROR_BadConversion(index, e.type, PT_INTEGER));
		}

	///	return element in param stack casted to number (double/float)
		number to_number(int index) const
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);
			
			const Entry& e = m_entries[index];
			if(e.type == PT_INTEGER)
				return (number)e.param.m_int;
			else if(e.type == PT_NUMBER)
				return e.param.m_number;
			
			throw(ERROR_BadConversion(index, e.type, PT_NUMBER));
		}

	///	return element in param stack casted to const char*
		const char* to_cstring(int index) const
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);
			
			const Entry& e = m_entries[index];

			if(e.type == PT_CSTRING)
				return e.param.m_cstring;
			else if(e.type == PT_STD_STRING)
				return e.param.m_stdString->c_str();
			
			throw(ERROR_BadConversion(index, e.type, PT_CSTRING));
		}

	///	return const reference to std::string
		const std::string& to_std_string(int index) const
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);

			const Entry& e = m_entries[index];

			if(e.type == PT_STD_STRING)
				return *e.param.m_stdString;

			throw(ERROR_BadConversion(index, e.type, PT_STD_STRING));
		}

	///	return element in param stack casted to user defined type
		template <class T>
		T* to_pointer(int index) const
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);
			
			const Entry& e = m_entries[index];
			if(e.type == PT_POINTER)
			{
			//	copy pointer; only copy will be changed
				const ClassNameNode* pClassNameNode = e.pClassNameNode;
				void* ptr = ClassCastProvider::cast_to_base_class(
						e.param.m_ptr, pClassNameNode, ClassNameProvider<T>::name());

				if(ptr != NULL)
					return reinterpret_cast<T*>(ptr);
				else
					throw(ERROR_IncompatibleClasses(index, class_name(index),
					                                ClassNameProvider<T>::name()));
			}
			
			throw(ERROR_BadConversion(index, e.type, PT_POINTER));
		}

	///	return element in param stack casted to void*
		void* to_pointer(int index) const
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);

			const Entry& e = m_entries[index];
			if(e.type == PT_POINTER)
				return e.param.m_ptr;

			throw(ERROR_BadConversion(index, e.type, PT_POINTER));
		}

	///	return element in param stack casted to user defined type
		template <class T>
		const T* to_const_pointer(int index) const
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);
			
			const Entry& e = m_entries[index];
			if(e.type == PT_CONST_POINTER)
			{
			//	copy pointer; only copy will be changed
				const ClassNameNode* pClassNameNode = e.pClassNameNode;
				void* ptr = ClassCastProvider::cast_to_base_class(
						e.param.m_ptr, pClassNameNode, ClassNameProvider<T>::name());

				if(ptr != NULL)
					return reinterpret_cast<const T*>(ptr);
				else
					throw(ERROR_IncompatibleClasses(index, class_name(index),
													ClassNameProvider<T>::name()));
			}
			
			throw(ERROR_BadConversion(index, e.type, PT_CONST_POINTER));
		}

	///	return element in param stack casted to void*
		const void* to_const_pointer(int index) const
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);

			const Entry& e = m_entries[index];
			if(e.type == PT_CONST_POINTER)
				return e.param.m_constPtr;

			throw(ERROR_BadConversion(index, e.type, PT_CONST_POINTER));
		}

	///	return element in param stack casted to user defined type in SmartPtr
	/**	\todo	currently only works for FreePolicy FreeDelete...*/
		template <class T>
		SmartPtr<T> to_smart_pointer(int index) const
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);

			const Entry& e = m_entries[index];
			if(e.type == PT_SMART_POINTER)
			{
			//	copy pointer; only copy will be changed
				const ClassNameNode* pClassNameNode = e.pClassNameNode;

			//	copy smart pointer
				SmartPtr<void> smartPtr =  *(e.param.m_smartPtrWrapper);

			//	cast raw pointer
				void* rawPtr = smartPtr.get();

			// 	cast raw pointer to desired class
				rawPtr = ClassCastProvider::cast_to_base_class(
						rawPtr, pClassNameNode, ClassNameProvider<T>::name());

			//	set raw ptr into smart pointer
				smartPtr.set_impl<T, FreeDelete>(rawPtr);

				if(rawPtr != NULL)
					return smartPtr.cast_reinterpret<T, FreeDelete>();
				else
					throw(ERROR_IncompatibleClasses(index, class_name(index),
					                                ClassNameProvider<T>::name()));
			}

			throw(ERROR_BadConversion(index, e.type, PT_SMART_POINTER));
		}

	///	return element in param stack casted to user defined type in SmartPtr
		SmartPtr<void> to_smart_pointer(int index) const
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);

			const Entry& e = m_entries[index];
			if(e.type == PT_SMART_POINTER)
				return *e.param.m_smartPtrWrapper;

			throw(ERROR_BadConversion(index, e.type, PT_SMART_POINTER));
		}

	///	return element in param stack casted to user defined type in SmartPtr
	/**	\todo	currently only works for FreePolicy FreeDelete...*/
		template <class T>
		ConstSmartPtr<T> to_const_smart_pointer(int index) const
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);

			const Entry& e = m_entries[index];
			if(e.type == PT_CONST_SMART_POINTER)
			{
			//	copy pointer; only copy will be changed
				const ClassNameNode* pClassNameNode = e.pClassNameNode;

			//	copy smart pointer
				ConstSmartPtr<void> smartPtr =  *(e.param.m_constSmartPtrWrapper);

			//	cast raw pointer
				const void* rawPtrConst = smartPtr.get();
				void* rawPtr = const_cast<void*>(rawPtrConst);

			// 	cast raw pointer to desired class
				rawPtr = ClassCastProvider::cast_to_base_class(
						rawPtr, pClassNameNode, ClassNameProvider<T>::name());

			//	set raw ptr into smart pointer
				smartPtr.set_impl<T, FreeDelete>(const_cast<const void*>(rawPtr));

				if(rawPtr != NULL)
					return smartPtr.cast_reinterpret<T, FreeDelete>();
				else
					throw(ERROR_IncompatibleClasses(index, class_name(index),
					                                ClassNameProvider<T>::name()));
			}

			throw(ERROR_BadConversion(index, e.type, PT_CONST_SMART_POINTER));
		}

	///	return element in param stack casted to user defined type in SmartPtr
		ConstSmartPtr<void> to_const_smart_pointer(int index) const
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);

			const Entry& e = m_entries[index];
			if(e.type == PT_CONST_SMART_POINTER)
				return *e.param.m_constSmartPtrWrapper;

			throw(ERROR_BadConversion(index, e.type, PT_CONST_SMART_POINTER));
		}

	////////////////////////////////
	//	set		
	////////////////////////////////

		void set_bool(int index, bool val)
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);

			Entry& e = m_entries[index];
			if(e.type == PT_BOOL)
				e.param.m_bool = val;
			else if(e.type == PT_INTEGER)
				e.param.m_int = (bool)val;
			else
				throw(ERROR_BadConversion(index, e.type, PT_BOOL));
		}

		void set_integer(int index, int val)
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);
			
			Entry& e = m_entries[index];
			if(e.type == PT_INTEGER)
				e.param.m_int = val;
			else if(e.type == PT_NUMBER)
				e.param.m_number = (number)val;
			else
				throw(ERROR_BadConversion(index, e.type, PT_INTEGER));
		}

		void set_number(int index, number val)
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);
			
			Entry& e = m_entries[index];
			if(e.type == PT_INTEGER)
				e.param.m_int = (int)val;
			else if(e.type == PT_NUMBER)
				e.param.m_number = val;
			else
				throw(ERROR_BadConversion(index, e.type, PT_NUMBER));
		}

		void set_cstring(int index, const char* str)
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);

			Entry& e = m_entries[index];

			if(e.type == PT_CSTRING)
				e.param.m_cstring = str;
			else
				throw(ERROR_BadConversion(index, e.type, PT_CSTRING));
		}

		void set_std_string(int index, const char* str)
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);

			Entry& e = m_entries[index];
			if(e.type == PT_STD_STRING){
				*e.param.m_stdString = str;
			}
			else
				throw(ERROR_BadConversion(index, e.type, PT_STD_STRING));
		}

		void set_std_string(int index, const std::string& str)
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);
			
			Entry& e = m_entries[index];
			if(e.type == PT_STD_STRING){
				*e.param.m_stdString = str;
			}
			else
				throw(ERROR_BadConversion(index, e.type, PT_STD_STRING));
		}

		template <class T>
		void set_pointer(int index, T* ptr)
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);
			
			Entry& e = m_entries[index];
			if(e.type == PT_POINTER)
				e.param.m_ptr = (void*)ptr;
			else
				throw(ERROR_BadConversion(index, e.type, PT_POINTER));
		}
		
		template <class T>
		void set_const_pointer(int index, const T* ptr)
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);
			
			Entry& e = m_entries[index];
			if(e.type == PT_CONST_POINTER)
				e.param.m_constPtr = (void*)ptr;
			else
				throw(ERROR_BadConversion(index, e.type, PT_CONST_POINTER));
		}
		
		template <class T>
		void set_smart_pointer(int index, const SmartPtr<T>& ptr)
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);

			Entry& e = m_entries[index];
			if(e.type == PT_SMART_POINTER)
				*e.param.m_smartPtrWrapper = ptr;
			else
				throw(ERROR_BadConversion(index, e.type, PT_SMART_POINTER));
		}

		template <class T>
		void set_const_smart_pointer(int index, const ConstSmartPtr<T>& ptr)
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);

			Entry& e = m_entries[index];
			if(e.type == PT_CONST_SMART_POINTER)
				*e.param.m_constSmartPtrWrapper = ptr;
			else
				throw(ERROR_BadConversion(index, e.type, PT_CONST_SMART_POINTER));
		}

	///	returns true if a parameter of the stack has been named by user
		bool parameter_named(int index) const
		{
			return class_name_node(index)->named();
		}

	private:
		union Parameter{
			bool m_bool;
			int	m_int;
			number m_number;
			const char* m_cstring;
			std::string* m_stdString;
			void* m_ptr;
			const void* m_constPtr;
			SmartPtr<void>* m_smartPtrWrapper;
			ConstSmartPtr<void>* m_constSmartPtrWrapper;
		};
		
	///	structure to store a data entry with additional information
		struct Entry{			
			Parameter param;  	//< Storage of data (native data or pointer)
			uint type;			//<	enum ParameterTypes indicating stored type
			const ClassNameNode* pClassNameNode; //< class name for user defined data
		};

	//	This array is of fixed size, since we want to introduce a minimal
	//	overhead during argument assignment.
		Entry m_entries[PARAMETER_STACK_MAX_SIZE];

	///	number of currently stored entries
		int m_numEntries;

	//	variables that tell whether m_entries contains string copies or smart pointers
		bool m_bHasSmartPtrs;
		bool m_bHasStringCopies;
};

} // end namespace bridge
} // end namespace ug

#endif
