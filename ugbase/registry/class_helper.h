/**
 * \file class_helper.h
 *
 * \author Martin Rupp
 *
 * \date 20.10.2010
 *
 * ClassHierarchy implementation, GetClassHierarchy, FindClass
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */

#ifndef __H__UG_BRIDGE__CLASS_HELPER__
#define __H__UG_BRIDGE__CLASS_HELPER__

#include <string>
#include "common/ug_config.h"

#include "registry.h"
#include "class.h"
#include "global_function.h"

namespace ug
{
namespace bridge
{

// ClassHierarchy
//--------------
/**
 * \brief Class Hierarchy Helper Class for UG Registry
 * This class stores class names and their subclasses
 * \sa GetClassHierarchy
 */
class UG_API ClassHierarchy
{
	public:
		ClassHierarchy() : name(), bGroup(false), subclasses() {}

		/**
		 * adds the class c to the class hierarchy by attaching it to its base
		 * hierarchy (base hierarchy taken from c->class_names()). automatically
		 * creates nonexisting base hierarchy.
		 * \param	c		Class to be inserted
		 */
		void insert_class(const IExportedClass &c);

		/**
		 * searches the hierarchy for the classname name.
		 * \param	name	class name to be searched
		 * \return NULL if class name not found,
		 * 				otherwise ClassHierarchy with the class as base
		 * 				(find_class(name)->name() == name)
		 */
		ClassHierarchy *find_class(const char *name);

		bool operator < (const ClassHierarchy &other) const
		{
			return name < other.name;
		}

		void sort();

		std::string name;
		bool bGroup;
		std::vector<ClassHierarchy> subclasses;
};
/**
 * inits hierarchy with all classes of UGBridge
 */
UG_API void GetClassHierarchy(ClassHierarchy &hierarchy, const Registry &reg);

/**
 * Finds the class classname in the default ug registry and returns
 * IExportedClass pointer if found, otherwise NULL
 */
UG_API bool PrintFunctionInfo(Registry &reg, const char *functionname);

UG_API void PrintFunctionInfo(const ExportedFunctionBase &thefunc, bool isConst=false,
                       const char *classname=NULL, const char *highlightclassname=NULL);
UG_API void PrintConstructorInfo(const ExportedConstructor &constr, const char *classname,
		const char *highlightclassname=NULL);

UG_API const IExportedClass *FindClass(bridge::Registry &reg, const char* classname);
UG_API bool PrintClassHierarchy(Registry &reg, const char *classname);
UG_API void PrintClassInfo(const IExportedClass &c);
UG_API bool PrintClassInfo(Registry &reg, const char *classname);
UG_API bool ClassUsageExact(Registry &reg, const char *classname, bool OutParameters);
UG_API const IExportedClass *FindClass(Registry &reg, const char* classname);
UG_API const ExportedFunction *FindFunction(Registry &reg, const char *functionname);
UG_API std::string ParameterToString(const ParameterStack &par, int i);

} // end namespace
} // end namespace

#endif
