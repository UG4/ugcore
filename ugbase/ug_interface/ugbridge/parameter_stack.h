#ifndef __H__UG_BRIDGE__PARAMETER_STACK__
#define __H__UG_BRIDGE__PARAMETER_STACK__

typedef double number;

#define PUSH_PARAM_TO_STACK(paramVar, val, paramType)	{m_entries[m_numEntries].param.paramVar = (val);\
													 	m_entries[m_numEntries].type = (paramType);\
													 	++m_numEntries;}

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

static inline int ARRAY_INDEX_TO_STACK_INDEX(int index, int stackSize)
{
	int nIndex = index;
	if(nIndex < 0)
		nIndex = stackSize + nIndex;
	
	if(nIndex < 0 || nIndex >= stackSize)
		throw(ERROR_BadIndex(index));
}




const int PARAMETER_STACK_MAX_SIZE = 10;

enum ParameterTypes
{
	PT_UNKNOWN = 0,
	PT_INTEGER = 1,
	PT_NUMBER = 2,
	PT_STRING = 3,
	PT_REFERENCE = 4,
	PT_POINTER = 5
};

////////////////////////////////////////////////////////////////////////
///	A stack that can hold values together with their type-id.
/**
 * This class is mainly used as a intermediate parameter storage during
 * calls to ugbridge methods. Its focus is on being leightweight and fast,
 * which makes it a little unflexible at times.
 * Note that the maximal number of parameters is specified by the constant
 * PARAMETER_STACK_MAX_SIZE. This is set to 10 by default. Please note, that
 * this value should not be unnecessarily high. This wouldn't make sense anyway,
 * since the template-method-wrappers can't take any more parameters.
 *
 * Supported types are integer, number, string, reference and pointer.
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
		ParameterStack() : m_numEntries(0)			{}
		
	////////////////////////////////
	//	info
		inline int size() const		{return m_numEntries;}
		
	////////////////////////////////
	//	push
		inline void push_integer(int val)			{PUSH_PARAM_TO_STACK(m_int, val, PT_INTEGER);}
		inline void push_number(number val)			{PUSH_PARAM_TO_STACK(m_number, val, PT_NUMBER);}
		
	///	strings are not bufferd.
		inline void push_string(const char* str)	{PUSH_PARAM_TO_STACK(m_string, str, PT_STRING);}
		
		template<class T>
		inline void push_reference(T& ref)			{PUSH_PARAM_TO_STACK(m_ptr, (void*)&ref, PT_REFERENCE);}
		
		template<class T>
		inline void push_pointer(T* ptr)			{PUSH_PARAM_TO_STACK(m_ptr, (void*)ptr, PT_POINTER);}
		
		
	////////////////////////////////
	//	get
		int get_type(int index) const
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);
			return m_entries[index].type;
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
			if(e.type == PT_STRING)
				return e.param.m_string;
			
			throw(ERROR_BadConversion(index, e.type, PT_STRING));
		}
		
		template <class T>
		T& to_reference(int index) const
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);
			
			const Entry& e = m_entries[index];
			if(e.type == PT_REFERENCE)
				return *reinterpret_cast<T*>(e.param.m_ptr);
			
			throw(ERROR_BadConversion(index, e.type, PT_REFERENCE));
		}

		template <class T>
		T* to_pointer(int index) const
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);
			
			const Entry& e = m_entries[index];
			if(e.type == PT_POINTER)
				return reinterpret_cast<T*>(e.param.m_ptr);
			
			throw(ERROR_BadConversion(index, e.type, PT_POINTER));
		}

	////////////////////////////////
	//	set		
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

		void set_string(int index, const char* str)
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);
			
			Entry& e = m_entries[index];
			if(e.type == PT_STRING)
				e.param.m_string = str;
			else
				throw(ERROR_BadConversion(index, e.type, PT_STRING));
		}
		
		template <class T>
		void set_reference(int index, T& ref)
		{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);
			
			Entry& e = m_entries[index];
			if(e.type == PT_REFERENCE)
				e.param.m_ptr = (void*)&ref;
			else
				throw(ERROR_BadConversion(index, e.type, PT_REFERENCE));
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
		
	private:
		union Parameter{
			int	m_int;
			number m_number;
			const char* m_string;
			void* m_ptr;
		};
		
		struct Entry{
			Parameter param;
			int type;
		};

	//	This array is of fixed size, since we want to introduce a minimal
	//	overhead during argument assignment.
		Entry m_entries[PARAMETER_STACK_MAX_SIZE];
		int m_numEntries;
};

#endif
