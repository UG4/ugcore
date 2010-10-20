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
#include "registry.h"

using namespace std;

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
class ClassHierarchy
{
public:
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

	string name;
	vector<ClassHierarchy> subclasses;
};
/**
 * inits hierarchy with all classes of UGBridge
 */
void GetClassHierarchy(ClassHierarchy &hierarchy, const bridge::Registry &reg);

/**
 * Finds the class classname in the default ug registry and returns IExportedClass pointer if found, otherwise NULL
 */
const IExportedClass *FindClass(const char* classname);

}
}
