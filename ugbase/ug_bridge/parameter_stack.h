// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y2010

#include <cstring>
#include "class_name_provider.h"
#include "common/common.h"
#include "common/util/smart_pointer.h"

#ifndef __H__UG_BRIDGE__PARAMETER_STACK__
#define __H__UG_BRIDGE__PARAMETER_STACK__


#define PUSH_PARAM_TO_STACK(paramVar, val, paramType, clName)	{m_entries[m_numEntries].param.paramVar = (val);\
																m_entries[m_numEntries].type = (paramType);\
																m_entries[m_numEntries].pClassNames = (clName);\
																++m_numEntries;}

//	call the constructor and assign the smart-ptr afterwards.
#define PUSH_SP_TO_STACK(val, clName)				{new((SmartPtr<void>*)m_entries[m_numEntries].param.m_smartPtrWrapper) SmartPtr<void>;\
													 *(SmartPtr<void>*)m_entries[m_numEntries].param.m_smartPtrWrapper = (val);\
													 m_entries[m_numEntries].type = PT_SMART_POINTER;\
													 m_entries[m_numEntries].pClassNames = (clName);\
													 ++m_numEntries;}

//	call the constructor and assign the smart-ptr afterwards.
#define PUSH_CSP_TO_STACK(val, clName)				{new((ConstSmartPtr<void>*)m_entries[m_numEntries].param.m_constSmartPtrWrapper) ConstSmartPtr<void>;\
													 *(ConstSmartPtr<void>*)m_entries[m_numEntries].param.m_constSmartPtrWrapper = (val);\
													 m_entries[m_numEntries].type = PT_CONST_SMART_POINTER;\
													 m_entries[m_numEntries].pClassNames = (clName);\
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
	ERROR_IncompatibleClasses(int index, const char* from, const char* to) :
		m_index(index), m_from(from), m_to(to)	{}

	int	m_index;
	const char* m_from;
	const char* m_to;
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




const int PARAMETER_STACK_MAX_SIZE = 10;

enum ParameterTypes
{
	PT_UNKNOWN = 0,
	PT_BOOL,
	PT_INTEGER,
	PT_NUMBER,
	PT_STRING,
	PT_POINTER,
	PT_CONST_POINTER,
	PT_SMART_POINTER,
	PT_CONST_SMART_POINTER,
	PT_RANGE = 0xFFFF
};

enum ParameterFlags
{
	PF_STRING_COPY = 0x10000,
	PF_RANGE = 0xFFFF0000
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
		ParameterStack() : m_numEntries(0), m_bHasSmartPtrs(false), m_bHasStringCopies(false)	{}
		
		~ParameterStack()
		{
		//	we have to release all string copies and smart pointers.
			if(m_bHasStringCopies){
				for(int i = 0; i < m_numEntries; ++i){
					if((m_entries[i].type & PF_RANGE) == PF_STRING_COPY)
						delete[] m_entries[i].param.m_string;
				}
			}

			if(m_bHasSmartPtrs){
				for(int i = 0; i < m_numEntries; ++i){
					if(m_entries[i].type == PT_SMART_POINTER)
						((SmartPtr<void>*)m_entries[i].param.m_smartPtrWrapper)->invalidate();
					else if(m_entries[i].type == PT_CONST_SMART_POINTER)
						((ConstSmartPtr<void>*)m_entries[i].param.m_constSmartPtrWrapper)->invalidate();
				}
			}
		}

	////////////////////////////////
	//	info
		inline int size() const		{return m_numEntries;}
		
	////////////////////////////////
	//	push
		inline void push_bool(bool val = true)			{PUSH_PARAM_TO_STACK(m_bool, val, PT_BOOL, NULL);}
		inline void push_integer(int val = 0)			{PUSH_PARAM_TO_STACK(m_int, val, PT_INTEGER, NULL);}
		inline void push_number(number val = 0)			{PUSH_PARAM_TO_STACK(m_number, val, PT_NUMBER, NULL);}
		
	///	strings are not bufferd.
		inline void push_string(const char* str = "", bool bCopy = false)
		{
			if(bCopy){
			//	copy the string and store it
				int strSize = strlen(str) + 1;	// don't forget terminating 0
				char* tstr = new char[strSize];
				memcpy(tstr, str, strSize);
				PUSH_PARAM_TO_STACK(m_string, tstr, PT_STRING | PF_STRING_COPY, NULL);
				m_bHasStringCopies = true;
			}
			else{
				PUSH_PARAM_TO_STACK(m_string, str, PT_STRING, NULL);
			}
		}

	/// user defined classes
		template<class T>
		inline void push_pointer(T* ptr = NULL)			{PUSH_PARAM_TO_STACK(	m_ptr, (void*)ptr, PT_POINTER,
																				&ClassNameProvider<T>::names());}

		inline void push_pointer(void* ptr, const std::vector<const char*>* classNames)
														{PUSH_PARAM_TO_STACK(	m_ptr, ptr, PT_POINTER, classNames);}

	/// user defined classes
		template<class T>
		inline void push_const_pointer(const T* ptr = NULL)	{PUSH_PARAM_TO_STACK(	m_constPtr, (const void*)ptr, PT_CONST_POINTER,
																				&ClassNameProvider<T>::names());}

		inline void push_const_pointer(const void* ptr, const std::vector<const char*>* classNames)
														{PUSH_PARAM_TO_STACK(m_constPtr, ptr, PT_CONST_POINTER, classNames);}

	/// SmartPtrs to user defined classes
		template<class T>
		inline void push_smart_pointer(const SmartPtr<T>& ptr = SmartPtr<T>(NULL))
			{
				PUSH_SP_TO_STACK(ptr, &ClassNameProvider<T>::names());
				m_bHasSmartPtrs = true;
			}

		inline void push_smart_pointer(const SmartPtr<void>& ptr, const std::vector<const char*>* classNames)
			{
				PUSH_SP_TO_STACK(ptr, classNames);
				m_bHasSmartPtrs = true;
			}

	/// ConstSmartPtrs to user defined classes
		template<class T>
		inline void push_const_smart_pointer(const ConstSmartPtr<T>& ptr = ConstSmartPtr<T>(NULL))
			{
				PUSH_CSP_TO_STACK(ptr, &ClassNameProvider<T>::names());
				m_bHasSmartPtrs = true;
			}

		inline void push_const_smart_pointer(const ConstSmartPtr<void>& ptr, const std::vector<const char*>* classNames)
			{
				PUSH_CSP_TO_STACK(ptr, classNames);
				m_bHasSmartPtrs = true;
			}

	////////////////////////////////
	//	get
		uint get_type(int index) const
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);
		//	skip the flags by binary and operation
			return m_entries[index].type & PT_RANGE;
		}
		
		const char* class_name(int index) const
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);
			if(m_entries[index].pClassNames == NULL) return "";
			else if (m_entries[index].pClassNames->empty()) return "";
			return (*m_entries[index].pClassNames)[0];
		}

		const std::vector<const char*>* class_names(int index) const
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);
			return m_entries[index].pClassNames;
		}

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

		const char* to_string(int index) const
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);
			
			const Entry& e = m_entries[index];

			if((e.type & PT_RANGE) == PT_STRING)
				return e.param.m_string;
			
			throw(ERROR_BadConversion(index, e.type, PT_STRING));
		}

		template <class T>
		T* to_pointer(int index) const
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);
			
			const Entry& e = m_entries[index];
			if(e.type == PT_POINTER)
			{
				if(ClassNameVecContains(*e.pClassNames, ClassNameProvider<T>::name()))
					return reinterpret_cast<T*>(e.param.m_ptr);
				else
					throw(ERROR_IncompatibleClasses(index, class_name(index), ClassNameProvider<T>::name()));
			}
			
			throw(ERROR_BadConversion(index, e.type, PT_POINTER));
		}

		void* to_pointer(int index) const
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);

			const Entry& e = m_entries[index];
			if(e.type == PT_POINTER)
				return e.param.m_ptr;

			throw(ERROR_BadConversion(index, e.type, PT_POINTER));
		}

		template <class T>
		const T* to_const_pointer(int index) const
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);
			
			const Entry& e = m_entries[index];
			if(e.type == PT_CONST_POINTER)
			{
				if(ClassNameVecContains(*e.pClassNames, ClassNameProvider<T>::name()))
					return reinterpret_cast<const T*>(e.param.m_constPtr);
				else
					throw(ERROR_IncompatibleClasses(index, class_name(index), ClassNameProvider<T>::name()));
			}
			
			throw(ERROR_BadConversion(index, e.type, PT_CONST_POINTER));
		}

		const void* to_const_pointer(int index) const
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);

			const Entry& e = m_entries[index];
			if(e.type == PT_CONST_POINTER)
				return e.param.m_constPtr;

			throw(ERROR_BadConversion(index, e.type, PT_CONST_POINTER));
		}

		template <class T>
		SmartPtr<T> to_smart_pointer(int index) const
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);

			const Entry& e = m_entries[index];
			if(e.type == PT_SMART_POINTER)
			{
				if(ClassNameVecContains(*e.pClassNames, ClassNameProvider<T>::name()))
					return ((SmartPtr<void>*)e.param.m_smartPtrWrapper)->
													to_smart_pointer_reinterpret<T>();
				else
					throw(ERROR_IncompatibleClasses(index, class_name(index), ClassNameProvider<T>::name()));
			}

			throw(ERROR_BadConversion(index, e.type, PT_SMART_POINTER));
		}

		SmartPtr<void> to_smart_pointer(int index) const
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);

			const Entry& e = m_entries[index];
			if(e.type == PT_SMART_POINTER)
				return *(SmartPtr<void>*)e.param.m_smartPtrWrapper;

			throw(ERROR_BadConversion(index, e.type, PT_SMART_POINTER));
		}

		template <class T>
		ConstSmartPtr<T> to_const_smart_pointer(int index) const
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);

			const Entry& e = m_entries[index];
			if(e.type == PT_CONST_SMART_POINTER)
			{
				if(ClassNameVecContains(*e.pClassNames, ClassNameProvider<T>::name()))
					return ((ConstSmartPtr<void>*)e.param.m_constSmartPtrWrapper)->
													to_smart_pointer_reinterpret<T>();
				else
					throw(ERROR_IncompatibleClasses(index, class_name(index), ClassNameProvider<T>::name()));
			}

			throw(ERROR_BadConversion(index, e.type, PT_CONST_SMART_POINTER));
		}

		ConstSmartPtr<void> to_const_smart_pointer(int index) const
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);

			const Entry& e = m_entries[index];
			if(e.type == PT_CONST_SMART_POINTER)
				return *(ConstSmartPtr<void>*)e.param.m_constSmartPtrWrapper;

			throw(ERROR_BadConversion(index, e.type, PT_CONST_SMART_POINTER));
		}

	////////////////////////////////
	//	set		
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

		void set_string(int index, const char* str, bool bCopy = false)
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);
			
			Entry& e = m_entries[index];
		//	first check whether the old string has to be cleared
			if(e.type == (PT_STRING | PF_STRING_COPY)){
				delete[] e.param.m_string;
				e.type = PT_STRING;
			}

			if(e.type == PT_STRING){
				if(bCopy){
					int strSize = strlen(str) + 1;	// don't forget terminating 0
					char* tstr = new char[strSize];
					memcpy(tstr, str, strSize);
					e.param.m_string = tstr;
					e.type = PT_STRING | PF_STRING_COPY;
					m_bHasStringCopies = true;
				}
				else{
					e.param.m_string = str;
					e.type = PT_STRING;				
				}
			}
			else
				throw(ERROR_BadConversion(index, e.type, PT_STRING));
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
				*(SmartPtr<void>*)e.param.m_smartPtrWrapper = ptr;
			else
				throw(ERROR_BadConversion(index, e.type, PT_SMART_POINTER));
		}

		template <class T>
		void set_const_smart_pointer(int index, const ConstSmartPtr<T>& ptr)
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);

			Entry& e = m_entries[index];
			if(e.type == PT_CONST_SMART_POINTER)
				*(ConstSmartPtr<void>*)e.param.m_constSmartPtrWrapper = ptr;
			else
				throw(ERROR_BadConversion(index, e.type, PT_CONST_SMART_POINTER));
		}

		bool is_parameter_undeclared(int index) const
		{
			if(class_name(index) == NULL || strlen(class_name(index)) == 0)
				return false;
			return class_name(index)[0] == '[';
		}

	private:
		union Parameter{
			bool m_bool;
			int	m_int;
			number m_number;
			const char* m_string;
			void* m_ptr;
			const void* m_constPtr;
			byte m_smartPtrWrapper[sizeof(SmartPtr<void>)];
			byte m_constSmartPtrWrapper[sizeof(ConstSmartPtr<void>)];
		};
		
		struct Entry{			
			Parameter param;
			uint type;
			const std::vector<const char*>*	pClassNames;
		};

	//	This array is of fixed size, since we want to introduce a minimal
	//	overhead during argument assignment.
		Entry m_entries[PARAMETER_STACK_MAX_SIZE];
		int m_numEntries;

	//	variables that tell whether m_entries contains string copies or smart pointers
		bool m_bHasSmartPtrs;
		bool m_bHasStringCopies;
};

} // end namespace bridge
} // end namespace ug

#endif
