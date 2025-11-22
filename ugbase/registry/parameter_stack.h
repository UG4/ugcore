/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter, Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */
#ifndef IG_UGBASE_REGISTRY_PARAMETER_STACK_H
#define IG_UGBASE_REGISTRY_PARAMETER_STACK_H

#include <cstring>
#include "class_name_provider.h"
#include "function_traits.h"
#include "common/common.h"
#include "common/util/smart_pointer.h"
#include "common/util/variant.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_function_handle.h"
#endif


#include <iostream>
#define untested() ( std::cerr <<  "@@#\n@@@:"<< __FILE__ << ":"<< __LINE__ \
          <<":" << __func__ << "\n" )

namespace ug
{
namespace bridge
{

/// \addtogroup registry
/// \{

/// a stack holding parameter infos about a parameter stack
/**
 * This class is used to store type information about the entries in a parameter
 * list.
 *
 * Note that the maximal number of parameters is specified by the constant
 * PARAMETER_STACK_MAX_SIZE. Please note, that this value should not be
 * unnecessarily high. The appropriate choice is UG_REGISTRY_MAX_NUM_ARGS,
 * since the template-method-wrappers can't take any more parameters.
 *
 * Supported types are bool, integer, number, const char*, std::string,
 * reference, pointer and smart-pointer. References and pointers are stored in
 * a void*. The user is responsible to associate the correct types.
 */
class ParameterInfo
{
	protected:
	///	maximal number of parameter in a parameter list
		static constexpr int PARAMETER_STACK_MAX_SIZE = UG_REGISTRY_MAX_NUM_ARGS;

	///	help function to compute correct parameter index
		static inline int ARRAY_INDEX_TO_STACK_INDEX(int index, int stackSize)
		{
			int nIndex = index;
			if(nIndex < 0)
				nIndex = stackSize + nIndex;

			if(nIndex < 0 || nIndex >= stackSize)
				UG_THROW("Invalid index "<<nIndex<<" used in Parameter Stack.");

			return nIndex;
		}

	public:
	///	default constructor
		ParameterInfo() : m_numEntries(0) {}

	////////////////////////////////
	//	info
	////////////////////////////////

	///	returns number of parameters in the param stack
		inline int size() const		{return m_numEntries;}

	///	returns if index is a std::vector
		inline bool is_vector(int index) const{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);
			return m_vEntryType[index].bVector;
		}

	///	returns ParameterType enum of data type for a stack entry
		Variant::Type type(int index) const{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);
			return m_vEntryType[index].type;
		}

	///	returns the class name node for an element in the param stack
		const ClassNameNode* class_name_node(int index) const{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);
			if(m_vEntryType[index].pClassNameNode == nullptr)
				UG_THROW("ClassNameNode missing in Parameter stack.");
			return m_vEntryType[index].pClassNameNode;
		}

	///	returns the class name for an element in the param stack
		const char* class_name(int index) const{
			return class_name_node(index)->name().c_str();
		}

	///	returns true if a parameter of the stack has been named by user
		bool parameter_named(int index) const{
			return class_name_node(index)->named();
		}

	///////////////////////////////////////////////////////////////////////
	//	push_type
	///////////////////////////////////////////////////////////////////////
	protected:
		template <typename TType, typename TNode>
		inline void _push_type(){
			m_vEntryType[m_numEntries].type = Variant::type<TType>();
			m_vEntryType[m_numEntries].pClassNameNode = &ClassNameProvider<TNode>::class_name_node();
			m_vEntryType[m_numEntries].bVector = false;
			++m_numEntries;
		}

		template <typename TNative>
		inline void _push_type(){_push_type<TNative,TNative>();}

		template <typename TType, typename TNode>
		inline void _push_vector_type(){
			m_vEntryType[m_numEntries].type = Variant::type<TType>();
			m_vEntryType[m_numEntries].pClassNameNode = &ClassNameProvider<TNode>::class_name_node();
			m_vEntryType[m_numEntries].bVector = true;
			++m_numEntries;
		}

		template <typename TNative>
		inline void _push_vector_type(){_push_vector_type<TNative,TNative>();}

		template <typename T>
		struct PushType{
			static void push(ParameterInfo* This){
				T::___UG_REGISTRY_ERROR___FUNCTION_OR_METHOD_PARAMETERS_RESTRICTED_to__NATIVE_TYPES__or__POINTER_resp_SMARTPOINTER_to_registered_types____();
			}
		};

	public:
	///	pushes a type to the parameter stack
		template <typename T>
		inline void push_type(){PushType<T>::push(this);}

	protected:
	///	structure to store a data entry with additional information
		struct EntryType{
			EntryType() :
				type(Variant::VT_INVALID), pClassNameNode(nullptr), bVector(false)
			{}
			Variant::Type type;						//<	enum ParameterTypes indicating stored type
			const ClassNameNode* pClassNameNode; 	//< class name for user defined data
			bool bVector;							//< boolean if vector data
		};

	//	This array is of fixed size, since we want to introduce a minimal
	//	overhead during argument assignment.
		EntryType m_vEntryType[PARAMETER_STACK_MAX_SIZE];

	///	number of currently stored entries
		int m_numEntries;
};

// implementation native types
template <> struct ParameterInfo::PushType<bool>				{static void push(ParameterInfo* This){This->_push_type<bool>();}};
template <> struct ParameterInfo::PushType<int>					{static void push(ParameterInfo* This){This->_push_type<int>();}};
template <> struct ParameterInfo::PushType<size_t>				{static void push(ParameterInfo* This){This->_push_type<size_t>();}};
template <> struct ParameterInfo::PushType<float>				{static void push(ParameterInfo* This){This->_push_type<float>();}};
template <> struct ParameterInfo::PushType<double>				{static void push(ParameterInfo* This){This->_push_type<double>();}};
template <> struct ParameterInfo::PushType<const char*>			{static void push(ParameterInfo* This){This->_push_type<const char*>();}};
template <> struct ParameterInfo::PushType<std::string>			{static void push(ParameterInfo* This){This->_push_type<std::string>();}};
template <> struct ParameterInfo::PushType<const std::string&>	{static void push(ParameterInfo* This){This->_push_type<std::string>();}};
#ifdef UG_FOR_LUA
template <> struct ParameterInfo::PushType<LuaFunctionHandle>	{static void push(ParameterInfo* This){This->_push_type<LuaFunctionHandle>();}};
template <> struct ParameterInfo::PushType<LuaTableHandle>	{static void push(ParameterInfo* This){ This->_push_type<LuaTableHandle>();}};
#endif

// implementation pointers and references to registered types
template <typename TClass> struct ParameterInfo::PushType<TClass*>				{static void push(ParameterInfo* This){This->_push_type<void*, TClass>();}};
template <typename TClass> struct ParameterInfo::PushType<const TClass*>		{static void push(ParameterInfo* This){This->_push_type<const void*, TClass>();}};
template <typename TClass> struct ParameterInfo::PushType<TClass&>				{static void push(ParameterInfo* This){This->_push_type<void*, TClass>();}};
template <typename TClass> struct ParameterInfo::PushType<const TClass&>			{static void push(ParameterInfo* This){This->_push_type<const void*, TClass>();}};
template <typename TClass> struct ParameterInfo::PushType<SmartPtr<TClass> >		{static void push(ParameterInfo* This){This->_push_type<SmartPtr<void>, TClass>();}};
template <typename TClass> struct ParameterInfo::PushType<ConstSmartPtr<TClass> >	{static void push(ParameterInfo* This){This->_push_type<ConstSmartPtr<void>, TClass>();}};

// implementation for std::vector, std::vector& and const std::vector&  (native types)
template <> struct ParameterInfo::PushType<std::vector<bool> >			{static void push(ParameterInfo* This){This->_push_vector_type<bool>();}};
template <> struct ParameterInfo::PushType<std::vector<int> >			{static void push(ParameterInfo* This){This->_push_vector_type<int>();}};
template <> struct ParameterInfo::PushType<std::vector<size_t> >		{static void push(ParameterInfo* This){This->_push_vector_type<size_t>();}};
template <> struct ParameterInfo::PushType<std::vector<float> >			{static void push(ParameterInfo* This){This->_push_vector_type<float>();}};
template <> struct ParameterInfo::PushType<std::vector<double> >		{static void push(ParameterInfo* This){This->_push_vector_type<double>();}};
template <> struct ParameterInfo::PushType<std::vector<const char*> >	{static void push(ParameterInfo* This){This->_push_vector_type<const char*>();}};
template <> struct ParameterInfo::PushType<std::vector<std::string> >	{static void push(ParameterInfo* This){This->_push_vector_type<std::string>();}};

/* Note: we do not support non-const references, since the bindings do not support the back-copy into the the bound language
template <> struct ParameterInfo::PushType<std::vector<bool>&>			{static void push(ParameterInfo* This){This->_push_vector_type<bool>();}};
template <> struct ParameterInfo::PushType<std::vector<int>&>			{static void push(ParameterInfo* This){This->_push_vector_type<int>();}};
template <> struct ParameterInfo::PushType<std::vector<size_t>&>		{static void push(ParameterInfo* This){This->_push_vector_type<size_t>();}};
template <> struct ParameterInfo::PushType<std::vector<float>&>			{static void push(ParameterInfo* This){This->_push_vector_type<float>();}};
template <> struct ParameterInfo::PushType<std::vector<double>&>		{static void push(ParameterInfo* This){This->_push_vector_type<double>();}};
template <> struct ParameterInfo::PushType<std::vector<const char*>&>	{static void push(ParameterInfo* This){This->_push_vector_type<const char*>();}};
template <> struct ParameterInfo::PushType<std::vector<std::string>&>	{static void push(ParameterInfo* This){This->_push_vector_type<std::string>();}};
*/

template <> struct ParameterInfo::PushType<const std::vector<bool>&>		{static void push(ParameterInfo* This){This->_push_vector_type<bool>();}};
template <> struct ParameterInfo::PushType<const std::vector<int>&>			{static void push(ParameterInfo* This){This->_push_vector_type<int>();}};
template <> struct ParameterInfo::PushType<const std::vector<size_t>&>		{static void push(ParameterInfo* This){This->_push_vector_type<size_t>();}};
template <> struct ParameterInfo::PushType<const std::vector<float>&>		{static void push(ParameterInfo* This){This->_push_vector_type<float>();}};
template <> struct ParameterInfo::PushType<const std::vector<double>&>		{static void push(ParameterInfo* This){This->_push_vector_type<double>();}};
template <> struct ParameterInfo::PushType<const std::vector<const char*>&>	{static void push(ParameterInfo* This){This->_push_vector_type<const char*>();}};
template <> struct ParameterInfo::PushType<const std::vector<std::string>&>	{static void push(ParameterInfo* This){This->_push_vector_type<std::string>();}};

// implementation for std::vector, std::vector& and const std::vector&  (registered types)
template <typename TClass> struct ParameterInfo::PushType<std::vector<TClass*> >	{static void push(ParameterInfo* This){This->_push_vector_type<void*, TClass>();}};
//template <typename TClass> struct ParameterInfo::PushType<std::vector<TClass*>& >	{static void push(ParameterInfo* This){This->_push_vector_type<void*, TClass>();}};
template <typename TClass> struct ParameterInfo::PushType<const std::vector<TClass*>&>	{static void push(ParameterInfo* This){This->_push_vector_type<void*, TClass>();}};

template <typename TClass> struct ParameterInfo::PushType<std::vector<const TClass*> >	{static void push(ParameterInfo* This){This->_push_vector_type<const void*, TClass>();}};
//template <typename TClass> struct ParameterInfo::PushType<std::vector<const TClass*>& >	{static void push(ParameterInfo* This){This->_push_vector_type<const void*, TClass>();}};
template <typename TClass> struct ParameterInfo::PushType<const std::vector<const TClass*>&>	{static void push(ParameterInfo* This){This->_push_vector_type<const void*, TClass>();}};

template <typename TClass> struct ParameterInfo::PushType<std::vector<SmartPtr<TClass> > >	{static void push(ParameterInfo* This){This->_push_vector_type<SmartPtr<void>, TClass>();}};
//template <typename TClass> struct ParameterInfo::PushType<std::vector<SmartPtr<TClass> >& >	{static void push(ParameterInfo* This){This->_push_vector_type<SmartPtr<void>, TClass>();}};
template <typename TClass> struct ParameterInfo::PushType<const std::vector<SmartPtr<TClass> >&>	{static void push(ParameterInfo* This){This->_push_vector_type<SmartPtr<void>, TClass>();}};

template <typename TClass> struct ParameterInfo::PushType<std::vector<ConstSmartPtr<TClass> > >	{static void push(ParameterInfo* This){This->_push_vector_type<ConstSmartPtr<void>, TClass>();}};
//template <typename TClass> struct ParameterInfo::PushType<std::vector<ConstSmartPtr<TClass> >& >	{static void push(ParameterInfo* This){This->_push_vector_type<ConstSmartPtr<void>, TClass>();}};
template <typename TClass> struct ParameterInfo::PushType<const std::vector<ConstSmartPtr<TClass> >&>	{static void push(ParameterInfo* This){This->_push_vector_type<ConstSmartPtr<void>, TClass>();}};


////////////////////////////////////////////////////////////////////////
///	A stack that can hold values together with their type-id.
/**
 * This class is mainly used as an intermediate parameter storage during
 * calls to ugbridge methods. Its focus is on being lightweight and fast.
 *
 * Use push() to add new parameters to the stack.
 * Use to() to retrieve a value.
 * Use set() to set a value in an existing entry.
 *
 * Indices start with 0. Negative indices can be used to start indexing from
 * the top of the stack.
 */
class ParameterStack : public ParameterInfo
{

	///////////////////////////////////////////////////////////////////////
	//	push
	///////////////////////////////////////////////////////////////////////
	protected:
		/**
		 * pushes a native type or a ptr/smartptr to a registered type to
		 * the stack. The value is stored in a Variant (thus, ptr are stored
		 * as void*, const void*, SmartPtr<void>, ConstSmartPtr<void>)
		 *
		 * @param val	the parameter to push
		 */
		template <typename T>
		inline void _push_native(const T& val){
			this->push_type<T>();
			m_vEntry[m_numEntries-1] = Variant(val);
		}

		template <typename TPtr, typename TType>
		inline void _push_pointer(TPtr val){
			this->push_type<TType>();
			m_vEntry[m_numEntries-1] = Variant(val);
		}

		/**
		 * pushes a native type to a ptr/smartptr. In order to keep track of the
		 * concrete type of the object, the classNameNode is stored as well
		 *
		 * @param val		the value to push
		 * @param classNameNode		the values classNameNode
		 * \tparam T void-ptr-type (one of void*, const void*, SmartPtr<void>, ConstSmartPtr<void>)
		 */
		template <typename T>
		inline void _push_void_pointer(T val, const ClassNameNode* classNameNode){
			m_vEntry[m_numEntries] = Variant(val);
			m_vEntryType[m_numEntries].type = Variant::type<T>();
			m_vEntryType[m_numEntries].pClassNameNode = classNameNode;
			m_vEntryType[m_numEntries].bVector = false;
			++m_numEntries;
		}

		/**
		 * pushes an std::vector to the stack, holding native type entries.
		 *
		 * @param spVec a smart-ptr to the std::vector
		 */
		template<typename T>
		inline void _push_vector(SmartPtr<std::vector<T> > spVec)
		{
				this->push_type<std::vector<T> >();
				SmartPtr<void> sp = spVec;
				m_vEntry[m_numEntries-1] = Variant(sp);
		}

		/**
		 * pushes an std::vector to the stack, holding ptr/smartptr to
		 * user-defined (registered) type. Ptrs are cast to void, and in order
		 * to get the concrete type (for casting back), the ClassNameNode is
		 * stored.
		 *
		 * @param spVec a smart-ptr to the std::vector
		 * \tparam TVoid void-ptr-type (one of void*, const void*, SmartPtr<void>, ConstSmartPtr<void>)
		 */
		template <typename TVoid>
		inline void _push_void_pointer_vector(SmartPtr<std::vector<std::pair<TVoid, const ClassNameNode*> > > spVec,
		                                      const ClassNameNode* baseNameNode = nullptr){
			SmartPtr<void> sp = spVec;
			m_vEntry[m_numEntries] = Variant(sp);
			m_vEntryType[m_numEntries].type = Variant::type<TVoid>();
			m_vEntryType[m_numEntries].pClassNameNode = baseNameNode;
			m_vEntryType[m_numEntries].bVector = true;
			++m_numEntries;
		}

		/**
		 * pushes an std::vector to the stack, holding ptr/smartptr to
		 * user-defined (registered) type. Ptrs are cast to void, and in order
		 * to get the concrete type (for casting back), the ClassNameNode is
		 * stored.
		 *
		 * @param vec a smart-ptr to the std::vector
		 */
		template <typename TVoid, typename TPtr, typename TNode>
		inline void _push_pointer_vector(const std::vector<TPtr>& vec){
			SmartPtr<std::vector<std::pair<TVoid, const ClassNameNode*> > > spVec
				= SmartPtr<std::vector<std::pair<TVoid, const ClassNameNode*> > >(new std::vector<std::pair<TVoid, const ClassNameNode*> >());

			for(size_t i = 0; i < vec.size(); ++i){
				spVec->push_back(std::pair<TVoid, const ClassNameNode*>(vec[i], &ClassNameProvider<TNode>::class_name_node()));
			}
			_push_void_pointer_vector(spVec, &ClassNameProvider<TNode>::class_name_node());
		}

	public:
	/// push user defined classes casted to void
	///	\{
		inline void push(void* ptr, const ClassNameNode* classNameNode){_push_void_pointer<void*>(ptr, classNameNode);}
		inline void push(const void* ptr, const ClassNameNode* classNameNode){_push_void_pointer<const void*>(ptr, classNameNode);}
		inline void push(SmartPtr<void> ptr, const ClassNameNode* classNameNode){_push_void_pointer<SmartPtr<void> >(ptr, classNameNode);}
		inline void push(ConstSmartPtr<void> ptr, const ClassNameNode* classNameNode){_push_void_pointer<ConstSmartPtr<void> >(ptr, classNameNode);}
	///	\}

	///	push array type
	///	\{
		inline void push(SmartPtr<std::vector<std::pair<void*, const ClassNameNode*> > > spVec){_push_void_pointer_vector<void*>(spVec);}
		inline void push(SmartPtr<std::vector<std::pair<const void*, const ClassNameNode*> > > spVec){_push_void_pointer_vector<const void*>(spVec);}
		inline void push(SmartPtr<std::vector<std::pair<SmartPtr<void>, const ClassNameNode*> > > spVec){_push_void_pointer_vector<SmartPtr<void> >(spVec);}
		inline void push(SmartPtr<std::vector<std::pair<ConstSmartPtr<void>, const ClassNameNode*> > > spVec){_push_void_pointer_vector<ConstSmartPtr<void> >(spVec);}
	/// \}

	///	push native array type
	///	\{
		inline void push(SmartPtr<std::vector<bool> > spVec){_push_vector<bool>(spVec);}
		inline void push(SmartPtr<std::vector<size_t> > spVec){_push_vector<size_t>(spVec);}
		inline void push(SmartPtr<std::vector<int> > spVec){_push_vector<int>(spVec);}
		inline void push(SmartPtr<std::vector<float> > spVec){_push_vector<float>(spVec);}
		inline void push(SmartPtr<std::vector<double> > spVec){_push_vector<double>(spVec);}
		inline void push(SmartPtr<std::vector<const char*> > spVec){_push_vector<const char*>(spVec);}
		inline void push(SmartPtr<std::vector<std::string> > spVec){_push_vector<std::string>(spVec);}
	/// \}

	protected:
		template <typename T>
		struct PushType{
			static void push(ParameterStack* This, T data){
				T::___UG_REGISTRY_ERROR___FUNCTION_OR_METHOD_PARAMETERS_RESTRICTED_to__NATIVE_TYPES__or__POINTER_resp_SMARTPOINTER_to_registered_types____();
		}};

	public:
	///	return element in param stack cast to type
		template <typename T>
		inline void push(T data) {PushType<T>::push(this, data);}

	///////////////////////////////////////////////////////////////////////
	//	to
	///////////////////////////////////////////////////////////////////////
	protected:
	///	return element in param stack cast to native type
		template <typename T>
		inline T _to_native(int index) const{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);
			return m_vEntry[index].to<T>();
		}

	///	returns element in param stack cast to pointer type
	/**
	 * returns the element at index in the stack cast to a pointer type
	 *
	 * @param index		the stack index
	 * @return	the cast pointer
	 * \tparam	TPtr concrete pointer type
	 * \tparam	TVoid ptr-type (one of void*, const void*, SmartPtr<void>, ConstSmartPtr<void>)
	 */
		template <typename T, typename TPtr, typename TVoid>
		inline TPtr _to_pointer(int index) const{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);
			const ClassNameNode* pClassNameNode = m_vEntryType[index].pClassNameNode;
			TPtr ptr = ClassCastProvider::cast_to<T>(m_vEntry[index].to<TVoid>(), pClassNameNode);
			return ptr;
		}

	///	return element in param stack casted to native type vector
		template <typename T>
		inline std::vector<T>& _to_native_vector(int index) const{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);
			SmartPtr<void> smartPtr = m_vEntry[index].to<SmartPtr<void> >();
			SmartPtr<std::vector<T> > spVec = smartPtr.cast_reinterpret<std::vector<T>, FreeDelete>();
			if(spVec.invalid()) UG_THROW("Cannot cast back to std::vector<T> for native type.");
			return *spVec;
		}

	///	return element in param stack cast to native type vector
		template <typename T, typename TPtr, typename TVoid>
		inline std::vector<TPtr>& _to_pointer_vector(int index) const{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);
			SmartPtr<void> smartPtr = m_vEntry[index].to<SmartPtr<void> >();
			SmartPtr<std::vector<std::pair<TVoid, const ClassNameNode*> > > spVec = smartPtr.cast_reinterpret<std::vector<std::pair<TVoid, const ClassNameNode*> > , FreeDelete>();
			if(spVec.invalid()) UG_THROW("Cannot cast back to std::vector<T> for native type.");

			SmartPtr<std::vector<TPtr> > sp = SmartPtr<std::vector<TPtr> >(new std::vector<TPtr>());

			for(size_t i = 0; i < spVec->size(); ++i){
				sp->push_back(ClassCastProvider::cast_to<T>((*spVec)[i].first, (*spVec)[i].second));
			}

			 const_cast<std::vector<SmartPtr<void> >*>(&m_vStoredSmartPtr)->push_back(sp);
			return *sp;
		}
		std::vector<SmartPtr<void> > m_vStoredSmartPtr;

	///	return element in param stack cast to native type vector
		template <typename TPtr>
		inline SmartPtr<std::vector<std::pair<TPtr, const ClassNameNode*> > > _to_void_pointer_vector(int index) const{
			index = ARRAY_INDEX_TO_STACK_INDEX(index, m_numEntries);
			SmartPtr<void> smartPtr = m_vEntry[index].to<SmartPtr<void> >();
			SmartPtr<std::vector<std::pair<TPtr, const ClassNameNode*> > > spVec = smartPtr.cast_reinterpret<std::vector<std::pair<TPtr, const ClassNameNode*> > , FreeDelete>();
			if(spVec.invalid()) UG_THROW("Cannot cast back to std::vector<T> for native type.");

			return spVec;
		}

		template <typename T>
		struct ToType{
			static T to(const ParameterStack* This, int index){
				return T::___UG_REGISTRY_ERROR___FUNCTION_OR_METHOD_PARAMETERS_RESTRICTED_to__NATIVE_TYPES__or__POINTER_resp_SMARTPOINTER_to_registered_types____();
		}};

	public:
	///	return element in param stack cast to type
		template <typename T>
		inline T to(int index) const {return ToType<T>::to(this, index);}

	///	return element in param stack as plain variant
		const Variant& get(int index) const	{return m_vEntry[index];}
		
	private:
	///	fixed size array storing the data for a stack entry
		Variant m_vEntry[PARAMETER_STACK_MAX_SIZE];
};

////////////////////////////////////////////////////////////////////////////////
//	PushType
////////////////////////////////////////////////////////////////////////////////

// convert to native types
template <> struct ParameterStack::PushType<bool>			{static void push(ParameterStack* This, bool data)					{This->_push_native<bool>(data);}};
template <> struct ParameterStack::PushType<int>			{static void push(ParameterStack* This, int data)					{This->_push_native<int>(data);}};
template <> struct ParameterStack::PushType<size_t>			{static void push(ParameterStack* This, size_t data)				{This->_push_native<size_t>(data);}};
template <> struct ParameterStack::PushType<float>			{static void push(ParameterStack* This, float data)					{This->_push_native<float>(data);}};
template <> struct ParameterStack::PushType<double>			{static void push(ParameterStack* This, double data)				{This->_push_native<double>(data);}};
template <> struct ParameterStack::PushType<const char*>		{static void push(ParameterStack* This, const char* data)		{This->_push_native<const char*>(data);}};
template <> struct ParameterStack::PushType<std::string>		{static void push(ParameterStack* This, std::string data){This->_push_native<std::string>(data);}};
template <> struct ParameterStack::PushType<const std::string&>	{static void push(ParameterStack* This, const std::string& data){This->_push_native<std::string>(data);}};
#ifdef UG_FOR_LUA
template <> struct ParameterStack::PushType<LuaFunctionHandle>	{static void push(ParameterStack* This, LuaFunctionHandle data)	{This->_push_native<LuaFunctionHandle>(data);}};
template <> struct ParameterStack::PushType<LuaTableHandle>	{static void push(ParameterStack* This, LuaTableHandle data)	{This->_push_native<LuaTableHandle>(data);}};
#endif

// convert pointers to native types
template <> struct ParameterStack::PushType<std::vector<bool> >			{static void push(ParameterStack* This, const std::vector<bool>& spVec)			{This->push(SmartPtr<std::vector<bool> >(new std::vector<bool>(spVec)));}};
template <> struct ParameterStack::PushType<std::vector<int> >			{static void push(ParameterStack* This, const std::vector<int>& spVec)			{This->push(SmartPtr<std::vector<int> >(new std::vector<int>(spVec)));}};
template <> struct ParameterStack::PushType<std::vector<size_t> >		{static void push(ParameterStack* This, const std::vector<size_t>& spVec)		{This->push(SmartPtr<std::vector<size_t> >(new std::vector<size_t>(spVec)));}};
template <> struct ParameterStack::PushType<std::vector<float> >		{static void push(ParameterStack* This, const std::vector<float>& spVec)		{This->push(SmartPtr<std::vector<float> >(new std::vector<float>(spVec)));}};
template <> struct ParameterStack::PushType<std::vector<double> >		{static void push(ParameterStack* This, const std::vector<double>& spVec)		{This->push(SmartPtr<std::vector<double> >(new std::vector<double>(spVec)));}};
template <> struct ParameterStack::PushType<std::vector<const char*> >	{static void push(ParameterStack* This, const std::vector<const char*>& spVec)	{This->push(SmartPtr<std::vector<const char*> >(new std::vector<const char*>(spVec)));}};
template <> struct ParameterStack::PushType<std::vector<std::string> >	{static void push(ParameterStack* This, const std::vector<std::string>& spVec)	{This->push(SmartPtr<std::vector<std::string> >(new std::vector<std::string>(spVec)));}};

/* Note: we do not support non-const references, since the bindings do not support the back-copy into the the bound language
template <> struct ParameterStack::PushType<std::vector<bool>& >		{static void push(ParameterStack* This, const std::vector<bool>& spVec)			{This->push(SmartPtr<std::vector<bool> >(new std::vector<bool>(spVec)));}};
template <> struct ParameterStack::PushType<std::vector<int>& >			{static void push(ParameterStack* This, const std::vector<int>& spVec)			{This->push(SmartPtr<std::vector<int> >(new std::vector<int>(spVec)));}};
template <> struct ParameterStack::PushType<std::vector<size_t>& >		{static void push(ParameterStack* This, const std::vector<size_t>& spVec)		{This->push(SmartPtr<std::vector<size_t> >(new std::vector<size_t>(spVec)));}};
template <> struct ParameterStack::PushType<std::vector<float>& >		{static void push(ParameterStack* This, const std::vector<float>& spVec)		{This->push(SmartPtr<std::vector<float> >(new std::vector<float>(spVec)));}};
template <> struct ParameterStack::PushType<std::vector<double>& >		{static void push(ParameterStack* This, const std::vector<double>& spVec)		{This->push(SmartPtr<std::vector<double> >(new std::vector<double>(spVec)));}};
template <> struct ParameterStack::PushType<std::vector<const char*>& >	{static void push(ParameterStack* This, const std::vector<const char*>& spVec)	{This->push(SmartPtr<std::vector<const char*> >(new std::vector<const char*>(spVec)));}};
template <> struct ParameterStack::PushType<std::vector<std::string>& >	{static void push(ParameterStack* This, const std::vector<std::string>& spVec)	{This->push(SmartPtr<std::vector<std::string> >(new std::vector<std::string>(spVec)));}};
*/

template <> struct ParameterStack::PushType<const std::vector<bool>& >			{static void push(ParameterStack* This, const std::vector<bool>& spVec)			{This->push(SmartPtr<std::vector<bool> >(new std::vector<bool>(spVec)));}};
template <> struct ParameterStack::PushType<const std::vector<int>& >			{static void push(ParameterStack* This, const std::vector<int>& spVec)			{This->push(SmartPtr<std::vector<int> >(new std::vector<int>(spVec)));}};
template <> struct ParameterStack::PushType<const std::vector<size_t>& >		{static void push(ParameterStack* This, const std::vector<size_t>& spVec)		{This->push(SmartPtr<std::vector<size_t> >(new std::vector<size_t>(spVec)));}};
template <> struct ParameterStack::PushType<const std::vector<float>& >			{static void push(ParameterStack* This, const std::vector<float>& spVec)		{This->push(SmartPtr<std::vector<float> >(new std::vector<float>(spVec)));}};
template <> struct ParameterStack::PushType<const std::vector<double>& >		{static void push(ParameterStack* This, const std::vector<double>& spVec)		{This->push(SmartPtr<std::vector<double> >(new std::vector<double>(spVec)));}};
template <> struct ParameterStack::PushType<const std::vector<const char*>& >	{static void push(ParameterStack* This, const std::vector<const char*>& spVec)	{This->push(SmartPtr<std::vector<const char*> >(new std::vector<const char*>(spVec)));}};
template <> struct ParameterStack::PushType<const std::vector<std::string>& >	{static void push(ParameterStack* This, const std::vector<std::string>& spVec)	{This->push(SmartPtr<std::vector<std::string> >(new std::vector<std::string>(spVec)));}};

// convert push concrete pointer types
template <typename T> struct ParameterStack::PushType<T*>					{static void push(ParameterStack* This, T* data)				{This->_push_pointer<void*, T*>(data);}};
template <typename T> struct ParameterStack::PushType<T&>					{static void push(ParameterStack* This, T& data)				{PushType<T*>::push(This, &data);}};
template <typename T> struct ParameterStack::PushType<const T*>			{static void push(ParameterStack* This, const T* data)			{This->_push_pointer<const void*, const T*>(data);}};
template <typename T> struct ParameterStack::PushType<const T&>			{static void push(ParameterStack* This, const T& data)			{PushType<const T*>::push(This, &data);}};
template <typename T> struct ParameterStack::PushType<SmartPtr<T> >		{static void push(ParameterStack* This, SmartPtr<T> data)		{This->_push_pointer<SmartPtr<void> , SmartPtr<T> >(data);}};
template <typename T> struct ParameterStack::PushType<ConstSmartPtr<T> >	{static void push(ParameterStack* This, ConstSmartPtr<T> data)	{This->_push_pointer<ConstSmartPtr<void>, ConstSmartPtr<T> >(data);}};

// convert to std::vector, std::vector& and const std::vector& (registered types)
template<typename T> struct ParameterStack::PushType<std::vector<T*> >		 	{static void push(ParameterStack* This, const std::vector<T*>& data)	{This->_push_pointer_vector<void*, T*, T>(data);}};
//template<typename T> struct ParameterStack::PushType<std::vector<T*>& >		{static void push(ParameterStack* This, const std::vector<T*>& data)	{This->_push_pointer_vector<void*, T*, T>(data);}};
template<typename T> struct ParameterStack::PushType<const std::vector<T*>&>	{static void push(ParameterStack* This, const std::vector<T*>& data)	{This->_push_pointer_vector<void*, T*, T>(data);}};

template<typename T> struct ParameterStack::PushType<std::vector<const T*> >		{static void push(ParameterStack* This, const std::vector<const T*>& data)	{This->_push_pointer_vector<const void*, const T*, T>(data);}};
//template<typename T> struct ParameterStack::PushType<std::vector<const T*>& >		{static void push(ParameterStack* This, const std::vector<const T*>& data)	{This->_push_pointer_vector<const void*, const T*, T>(data);}};
template<typename T> struct ParameterStack::PushType<const std::vector<const T*>&>	{static void push(ParameterStack* This, const std::vector<const T*>& data)	{This->_push_pointer_vector<const void*, const T*, T>(data);}};

template<typename T> struct ParameterStack::PushType<std::vector<SmartPtr<T> > >		{static void push(ParameterStack* This, const std::vector<SmartPtr<T> >& data)	{This->_push_pointer_vector<SmartPtr<void>, SmartPtr<T>, T>(data);}};
//template<typename T> struct ParameterStack::PushType<std::vector<SmartPtr<T> >& >		{static void push(ParameterStack* This, const std::vector<SmartPtr<T> >& data)	{This->_push_pointer_vector<SmartPtr<void>, SmartPtr<T>, T>(data);}};
template<typename T> struct ParameterStack::PushType<const std::vector<SmartPtr<T> >&>	{static void push(ParameterStack* This, const std::vector<SmartPtr<T> >& data)	{This->_push_pointer_vector<SmartPtr<void>, SmartPtr<T>, T>(data);}};

template<typename T> struct ParameterStack::PushType<std::vector<ConstSmartPtr<T> > >		{static void push(ParameterStack* This, const std::vector<ConstSmartPtr<T> >& data)	{This->_push_pointer_vector<ConstSmartPtr<void>, ConstSmartPtr<T>, T>(data);}};
//template<typename T> struct ParameterStack::PushType<std::vector<ConstSmartPtr<T> >& >		{static void push(ParameterStack* This, const std::vector<ConstSmartPtr<T> >& data)	{This->_push_pointer_vector<ConstSmartPtr<void>, ConstSmartPtr<T>, T>(data);}};
template<typename T> struct ParameterStack::PushType<const std::vector<ConstSmartPtr<T> >&>{static void push(ParameterStack* This, const std::vector<ConstSmartPtr<T> >& data)	{This->_push_pointer_vector<ConstSmartPtr<void>, ConstSmartPtr<T>, T>(data);}};


////////////////////////////////////////////////////////////////////////////////
//	ToType
////////////////////////////////////////////////////////////////////////////////

// convert to native types
template <> struct ParameterStack::ToType<bool>				{static bool to(const ParameterStack* This, int index)					{return This->_to_native<bool>(index);}};
template <> struct ParameterStack::ToType<int>				{static int to(const ParameterStack* This, int index)					{return This->_to_native<int>(index);}};
template <> struct ParameterStack::ToType<size_t>			{static size_t to(const ParameterStack* This, int index)				{return This->_to_native<size_t>(index);}};
template <> struct ParameterStack::ToType<float>			{static float to(const ParameterStack* This, int index)					{return This->_to_native<float>(index);}};
template <> struct ParameterStack::ToType<double>			{static double to(const ParameterStack* This, int index)				{return This->_to_native<double>(index);}};
template <> struct ParameterStack::ToType<const char*>		{static const char* to(const ParameterStack* This, int index)			{return This->_to_native<const char*>(index);}};
template <> struct ParameterStack::ToType<std::string>		{static std::string to(const ParameterStack* This, int index)			{return This->_to_native<std::string>(index);}};
template <> struct ParameterStack::ToType<const std::string&>{static const std::string& to(const ParameterStack* This, int index)	{return This->_to_native<const std::string&>(index);}};
#ifdef UG_FOR_LUA
template <> struct ParameterStack::ToType<LuaFunctionHandle>	{static LuaFunctionHandle to(const ParameterStack* This, int index)		{return This->_to_native<LuaFunctionHandle>(index);}};
template <> struct ParameterStack::ToType<LuaTableHandle>	{static LuaTableHandle to(const ParameterStack* This, int index)		{return This->_to_native<LuaTableHandle>(index);}};
#endif

// convert to void types
template <> struct ParameterStack::ToType<void*>			{static void* to(const ParameterStack* This, int index)						{return This->_to_native<void*>(index);}};
template <> struct ParameterStack::ToType<const void*>		{static const void* to(const ParameterStack* This, int index)				{return This->_to_native<const void*>(index);}};
template <> struct ParameterStack::ToType<SmartPtr<void> >	{static SmartPtr<void> to(const ParameterStack* This, int index)			{return This->_to_native<SmartPtr<void> >(index);}};
template <> struct ParameterStack::ToType<ConstSmartPtr<void> >	{static ConstSmartPtr<void> to(const ParameterStack* This, int index)	{return This->_to_native<ConstSmartPtr<void> >(index);}};

// convert to concrete pointer types
template <typename T> struct ParameterStack::ToType<T*>				{static T* to(const ParameterStack* This, int index){return This->_to_pointer<T, T*, void*>(index);}};
template <typename T> struct ParameterStack::ToType<T&>				{static T& to(const ParameterStack* This, int index){return *ToType<T*>::to(This, index);}};
template <typename T> struct ParameterStack::ToType<const T*>			{static const T* to(const ParameterStack* This, int index){return This->_to_pointer<T, const T*, const void*>(index);}};
template <typename T> struct ParameterStack::ToType<const T&>			{static const T& to(const ParameterStack* This, int index){return *ToType<const T*>::to(This, index);}};
template <typename T> struct ParameterStack::ToType<SmartPtr<T> >		{static SmartPtr<T> to(const ParameterStack* This, int index){return This->_to_pointer<T, SmartPtr<T>, SmartPtr<void> >(index);}};
template <typename T> struct ParameterStack::ToType<ConstSmartPtr<T> >	{static ConstSmartPtr<T> to(const ParameterStack* This, int index){return This->_to_pointer<T, ConstSmartPtr<T>, ConstSmartPtr<void> >(index);}};

// convert to std::vector, std::vector& and const std::vector& (native types)
template<> struct ParameterStack::ToType<std::vector<bool> >		{static std::vector<bool> to(const ParameterStack* This, int index)			{return This->_to_native_vector<bool>(index);}};
template<> struct ParameterStack::ToType<std::vector<int> >			{static std::vector<int> to(const ParameterStack* This, int index)			{return This->_to_native_vector<int>(index);}};
template<> struct ParameterStack::ToType<std::vector<size_t> >		{static std::vector<size_t> to(const ParameterStack* This, int index)		{return This->_to_native_vector<size_t>(index);}};
template<> struct ParameterStack::ToType<std::vector<float> >		{static std::vector<float> to(const ParameterStack* This, int index)		{return This->_to_native_vector<float>(index);}};
template<> struct ParameterStack::ToType<std::vector<double> >		{static std::vector<double> to(const ParameterStack* This, int index)		{return This->_to_native_vector<double>(index);}};
template<> struct ParameterStack::ToType<std::vector<const char*> >	{static std::vector<const char*> to(const ParameterStack* This, int index)	{return This->_to_native_vector<const char*>(index);}};
template<> struct ParameterStack::ToType<std::vector<std::string> >	{static std::vector<std::string> to(const ParameterStack* This, int index)	{return This->_to_native_vector<std::string>(index);}};

/* Note: we do not support non-const references, since the bindings do not support the back-copy into the the binded language
template<> struct ParameterStack::ToType<std::vector<bool>&>		{static std::vector<bool>& to(const ParameterStack* This, int index)			{return This->_to_native_vector<bool>(index);}};
template<> struct ParameterStack::ToType<std::vector<int>&>			{static std::vector<int>& to(const ParameterStack* This, int index)			{return This->_to_native_vector<int>(index);}};
template<> struct ParameterStack::ToType<std::vector<size_t>&>		{static std::vector<size_t>& to(const ParameterStack* This, int index)		{return This->_to_native_vector<size_t>(index);}};
template<> struct ParameterStack::ToType<std::vector<float>&>		{static std::vector<float>& to(const ParameterStack* This, int index)			{return This->_to_native_vector<float>(index);}};
template<> struct ParameterStack::ToType<std::vector<double>&>		{static std::vector<double>& to(const ParameterStack* This, int index)		{return This->_to_native_vector<double>(index);}};
template<> struct ParameterStack::ToType<std::vector<const char*>&>	{static std::vector<const char*>& to(const ParameterStack* This, int index)	{return This->_to_native_vector<const char*>(index);}};
template<> struct ParameterStack::ToType<std::vector<std::string>&>	{static std::vector<std::string>& to(const ParameterStack* This, int index)	{return This->_to_native_vector<std::string>(index);}};
*/

template<> struct ParameterStack::ToType<const std::vector<bool>&>			{static const std::vector<bool>& to(const ParameterStack* This, int index)			{return This->_to_native_vector<bool>(index);}};
template<> struct ParameterStack::ToType<const std::vector<int>&>			{static const std::vector<int>& to(const ParameterStack* This, int index)			{return This->_to_native_vector<int>(index);}};
template<> struct ParameterStack::ToType<const std::vector<size_t>&>		{static const std::vector<size_t>& to(const ParameterStack* This, int index)		{return This->_to_native_vector<size_t>(index);}};
template<> struct ParameterStack::ToType<const std::vector<float>&>			{static const std::vector<float>& to(const ParameterStack* This, int index)			{return This->_to_native_vector<float>(index);}};
template<> struct ParameterStack::ToType<const std::vector<double>&>		{static const std::vector<double>& to(const ParameterStack* This, int index)		{return This->_to_native_vector<double>(index);}};
template<> struct ParameterStack::ToType<const std::vector<const char*>&>	{static const std::vector<const char*>& to(const ParameterStack* This, int index)	{return This->_to_native_vector<const char*>(index);}};
template<> struct ParameterStack::ToType<const std::vector<std::string>&>	{static const std::vector<std::string>& to(const ParameterStack* This, int index)	{return This->_to_native_vector<std::string>(index);}};

// convert to std::vector, std::vector& and const std::vector& (registered types)
template<typename T> struct ParameterStack::ToType<std::vector<T*> >{static std::vector<T*> to(const ParameterStack* This, int index){return This->_to_pointer_vector<T, T*, void*>(index);}};
//template<typename T> struct ParameterStack::ToType<std::vector<T*>& >{static std::vector<T*>& to(const ParameterStack* This, int index){return This->_to_pointer_vector<T, T*, void*>(index);}};
template<typename T> struct ParameterStack::ToType<const std::vector<T*>&>{static const std::vector<T*>& to(const ParameterStack* This, int index){return This->_to_pointer_vector<T, T*, void*>(index);}};

template<typename T> struct ParameterStack::ToType<std::vector<const T*> >{static std::vector<const T*> to(const ParameterStack* This, int index){return This->_to_pointer_vector<T, const T*, const void*>(index);}};
//template<typename T> struct ParameterStack::ToType<std::vector<const T*>& >{static std::vector<const T*>& to(const ParameterStack* This, int index){return This->_to_pointer_vector<T, const T*, const void*>(index);}};
template<typename T> struct ParameterStack::ToType<const std::vector<const T*>&>{static const std::vector<const T*>& to(const ParameterStack* This, int index){return This->_to_pointer_vector<T, const T*, const void*>(index);}};

template<typename T> struct ParameterStack::ToType<std::vector<SmartPtr<T> > >{static std::vector<SmartPtr<T> > to(const ParameterStack* This, int index){return This->_to_pointer_vector<T, SmartPtr<T>, SmartPtr<void> >(index);}};
//template<typename T> struct ParameterStack::ToType<std::vector<SmartPtr<T> >& >{static std::vector<SmartPtr<T> >& to(const ParameterStack* This, int index){return This->_to_pointer_vector<T, SmartPtr<T>, SmartPtr<void> >(index);}};
template<typename T> struct ParameterStack::ToType<const std::vector<SmartPtr<T> >&>{static const std::vector<SmartPtr<T> >& to(const ParameterStack* This, int index){return This->_to_pointer_vector<T, SmartPtr<T>, SmartPtr<void> >(index);}};

template<typename T> struct ParameterStack::ToType<std::vector<ConstSmartPtr<T> > >{static std::vector<ConstSmartPtr<T> > to(const ParameterStack* This, int index){return This->_to_pointer_vector<T, ConstSmartPtr<T>, ConstSmartPtr<void> >(index);}};
//template<typename T> struct ParameterStack::ToType<std::vector<ConstSmartPtr<T> >& >{static std::vector<ConstSmartPtr<T> >& to(const ParameterStack* This, int index){return This->_to_pointer_vector<T, ConstSmartPtr<T>, ConstSmartPtr<void> >(index);}};
template<typename T> struct ParameterStack::ToType<const std::vector<ConstSmartPtr<T> >&>{static const std::vector<ConstSmartPtr<T> >& to(const ParameterStack* This, int index){return This->_to_pointer_vector<T, ConstSmartPtr<T>, ConstSmartPtr<void> >(index);}};

// convert to std::vector for void pointer (registered types)
template<> struct ParameterStack::ToType<SmartPtr<std::vector<std::pair<void*, const ClassNameNode*> > > >{static SmartPtr<std::vector<std::pair<void*, const ClassNameNode*> > >  to(const ParameterStack* This, int index){return This->_to_void_pointer_vector<void*>(index);}};
template<> struct ParameterStack::ToType<SmartPtr<std::vector<std::pair<const void*, const ClassNameNode*> > > >{static SmartPtr<std::vector<std::pair<const void*, const ClassNameNode*> > >  to(const ParameterStack* This, int index){return This->_to_void_pointer_vector<const void*>(index);}};
template<> struct ParameterStack::ToType<SmartPtr<std::vector<std::pair<SmartPtr<void>, const ClassNameNode*> > > >{static SmartPtr<std::vector<std::pair<SmartPtr<void>, const ClassNameNode*> > >  to(const ParameterStack* This, int index){return This->_to_void_pointer_vector<SmartPtr<void> >(index);}};
template<> struct ParameterStack::ToType<SmartPtr<std::vector<std::pair<ConstSmartPtr<void>, const ClassNameNode*> > > >{static SmartPtr<std::vector<std::pair<ConstSmartPtr<void>, const ClassNameNode*> > >  to(const ParameterStack* This, int index){return This->_to_void_pointer_vector<ConstSmartPtr<void> >(index);}};

// end group registry
/// \}

} // end namespace bridge
} // end namespace ug

#endif