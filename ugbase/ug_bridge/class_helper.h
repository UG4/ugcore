/**
 * \file class_helper.cpp
 *
 * \author Martin Rupp
 *
 * \date 20.10.2010
 *
 * ClassHierarchy implementation, GetClassHierarchy, FindClass
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */

#include <string>

#ifndef __H__UG_BRIDGE__CLASS_HELPER__
#define __H__UG_BRIDGE__CLASS_HELPER__

namespace ug
{
namespace bridge
{

class IExportedClass;
class ExportedFunctionBase;
class ExportedFunction;
class Registry;

// ClassHierarchy
//--------------
/**
 * \brief Class Hierarchy Helper Class for UG Registry
 * This class stores class names and their subclasses
 * \sa GetClassHierarchy
 */
class ClassHierarchy
{
public:
	ClassHierarchy() : name(), bGroup(false), subclasses() {}
	/**
	 * adds the class c to the class hierarchy by attaching it to its base hierarchy
	 * (base hierarchy taken from c->class_names()). automatically creates nonexisting base hierarchy.
	 * \param	c		Class to be inserted
	 */
	void insert_class(const IExportedClass &c);
	/**
	 * searches the hierarchy for the classname name.
	 * \param	name	class name to be searched
	 * \return NULL if class name not found, otherwise ClassHierarchy with the class as base (find_class(name)->name() == name)
	 */
	ClassHierarchy *find_class(const char *name);

	bool operator < (const ClassHierarchy &other) const
	{
		return name < other.name;
	}

	void sort()
	{
		std::sort(subclasses.begin(), subclasses.end());
		for(size_t i=0; i<subclasses.size(); i++)
			subclasses[i].sort();
	}

	std::string name;
	bool bGroup;
	std::vector<ClassHierarchy> subclasses;
};
/**
 * inits hierarchy with all classes of UGBridge
 */
void GetClassHierarchy(ClassHierarchy &hierarchy, const Registry &reg);

/**
 * Finds the class classname in the default ug registry and returns IExportedClass pointer if found, otherwise NULL
 */
bool PrintFunctionInfo(Registry &reg, const char *functionname);
void PrintFunctionInfo(const ExportedFunctionBase &thefunc, bool isConst=false, const char *classname=NULL, const char *highlightclassname=NULL);

const IExportedClass *FindClass(bridge::Registry &reg, const char* classname);
bool PrintClassHierarchy(Registry &reg, const char *classname);
void PrintClassInfo(const IExportedClass &c);
bool PrintClassInfo(Registry &reg, const char *classname);
bool ClassUsageExact(Registry &reg, const char *classname, bool OutParameters);
const IExportedClass *FindClass(Registry &reg, const char* classname);
const ExportedFunction *FindFunction(Registry &reg, const char *functionname);
}
}

#endif
