//	authors: Sebastian Reiter, Andreas Vogel


#ifndef __H__UG_BRIDGE__CLASS_NAME_PROVIDER__
#define __H__UG_BRIDGE__CLASS_NAME_PROVIDER__

#include "common/assert.h"
#include <vector>
#include <cstring>
#include <string>
#include <algorithm>
#include <iostream>
#include <typeinfo>
#include <map>
#include "common/common.h"
#include "common/util/smart_pointer.h"
#include "common/ug_config.h"
#include "error.h"

namespace ug
{
namespace bridge
{

/// \addtogroup registry
/// \{

/// node for class names
/**
 * A ClassNameNode stores the name of a registered class and pointers to
 * the ClassNameNodes of the direct base classes of this class. By traversing
 * the tree of ClassNameNodes all parent classes of a class can be found.
 */
class UG_API ClassNameNode
{
	public:
	///	constructor
		ClassNameNode();

	///	set name
		void set_name(const std::string& name);

	///	add a base class
		void add_base_class(const ClassNameNode& node);

	///	returns own name
		const std::string& name() const {return m_name;}

	///	returns if a name has been set
		bool empty() const {return m_name.empty();}

	///	returns if a name has been set by the user
		bool named() const;

	///	returns number of parents
		size_t num_base_classes() const {return m_vBaseClass.size();}

	///	return a base class
		const ClassNameNode& base_class(size_t i) const {return *m_vBaseClass.at(i);}

	protected:
	///	own name
		std::string m_name;

	///	base classes
		std::vector<const ClassNameNode*> m_vBaseClass;
};

///	provides the name for a class
template <typename TClass>
class UG_API ClassNameProvider
{
	public:
	/// set name of class and copy parent names
		static void set_name(const std::string& nameIn, const std::string& group,
		                     bool newName = false);

	/// set name of class and copy parent names
		template <typename TParent1>
		static void set_name(const std::string& nameIn, const std::string& group,
		                     bool newName = false);

	/// set name of class and copy parent names
		template <typename TParent1, typename TParent2>
		static void set_name(const std::string& nameIn, const std::string& group,
		                     bool newName = false);

	/** check if a given class name is equal
	 * This function checks if a given class name is equal to this class name.
	 * If the flag 'strict' is set to true, the class name must match exactly.
	 * If it is set to false, also parent names (of the class hierarchy) are checked
	 * and if this class inherits from the parent class true is returned.
	 */
		static bool is_a(const std::string& parent, bool strict = false);

	/// name of this class
		static const std::string& name();

	/// groups
		static const std::string& group(){return m_group;}

	/// returns vector of all names including parent class names
		static const std::vector<const char*>& names();

	///	return the class name node in the class hierarchy
		static const ClassNameNode& class_name_node() {return m_ClassNameNode;}

	///	returns if class name is forward declared
		static bool forward_declared() {return m_bForwardDeclared;}

	///	returns if the class has been named by user
		static bool named() {return !m_bForwardDeclared && !m_ClassNameNode.empty();}

	protected:
	///	sets a temporary name to the class
		static void set_foreward_declared();

	private:
	///	vector of parent class names (depreciated)
		static std::vector<const char*> m_names;

	///	Name of group, we're this class is sorted
		static std::string m_group;

	///	set to true if class has not been named by user, but a default name given
		static bool m_bForwardDeclared;

	///	class name node holding own name and pointers to base class name nodes
		static ClassNameNode m_ClassNameNode;
};

template <typename TClass>
const char* GetClassName(){
	return ClassNameProvider<TClass>::name().c_str();
}

///	static cast function for two classes
template <typename TBase, typename TDerived>
void* StaticVoidCast(void* DerivVoidPtr);


struct UGError_ClassCastFailed : public UGError{
	UGError_ClassCastFailed(const std::string& from, const std::string& to) :
		UGError("Class cast failed"), m_from(from), m_to(to)	{}

	std::string m_from;
	std::string m_to;
};

///	provides castings from derived classes to base classes
class ClassCastProvider
{
	public:
	///	add a cast function to the registry: Casts: Derived -> Base
		template <typename TBase, typename TDerived>
		static void add_cast_func();

	///	cast a pointer to the desired base class
	/**
	 * This method casts a void pointer to a given derived class to the void
	 * pointer of a base class. If conversion fails, an exception of type
	 * UGError_ClassCastFailed is thrown.
	 *
	 * \param[in]	pDerivVoid		void pointer to Derived object
	 * \param[in,out]	node		on entry: class name node corresponding to pDerivVoid
	 * 								on exit:  class name node corresponding to baseName
	 * \param[in]		baseName	name of base class the pointer should be casted to
	 * \returns		void* to base class
	 */
	/// \{
		static void* cast_to_base_class(void* pDerivVoid,
		                                const ClassNameNode*& node,
		                                const std::string& baseName);
		static const void* cast_to_base_class(const void* pDerivVoid,
		                                      const ClassNameNode*& node,
		                                      const std::string& baseName);
	//// \}

	///	casts a void pointer to a concrete class
	/**
	 * This method casts a void pointer to a given derived classed and returns it
	 * as a reinterpreted cast to the type specified by the template argument.
	 * If conversion fails, an exception of type UGError_ClassCastFailed is thrown.
	 *
	 * \param[in]	pDerivVoid		void pointer to Derived object
	 * \param[in,out]	node		on entry: class name node corresponding to pDerivVoid
	 * 								on exit:  class name node corresponding to baseName
	 * \returns		void* to base class
	 */
	/// \{
		template <typename T>
		static T* cast_to(void* pDerivVoid, const ClassNameNode*& node);
		template <typename T>
		static const T* cast_to(const void* pDerivVoid, const ClassNameNode*& node);
		template <typename T>
		static SmartPtr<T> cast_to(SmartPtr<void> pDerivVoid, const ClassNameNode*& node);
		template <typename T>
		static ConstSmartPtr<T> cast_to(ConstSmartPtr<void> pDerivVoid, const ClassNameNode*& node);
	///	\}

	protected:
	//	type of cast pointer
		typedef void* (*CastFunc)(void*);

	//	cast map
		static std::map<std::pair<const ClassNameNode*, const ClassNameNode*>, CastFunc> m_mmCast;
};


/// returns the vector containing all names in the name tree for node and its base classes
void ExtractClassNameVec(std::vector<const char*>& names,
                         const ClassNameNode& node,
                         bool clearVec = true);

///	returns if a name is contained in the name vector
bool ClassNameVecContains(const std::vector<const char*>& names, const std::string& name);

///	returns if a name is contained in the name tree at node or in base classes
bool ClassNameTreeContains(const ClassNameNode& node, const std::string& name);

/// returns an std::vector that contains in reverse order the base class that
///	must be used in the Class Hierarchy to get to the base class
bool ClassNameTreeWay(std::vector<size_t>& vWay, const ClassNameNode& node, const std::string& name);

// end group registry
/// \}

} // end namespace
} // end namespace

// include implementation
#include "class_name_provider_impl.h"

#endif /* __H__UG_BRIDGE__CLASS_NAME_PROVIDER__ */
