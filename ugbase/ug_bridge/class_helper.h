#include <string>
#include "registry.h"

using namespace std;

namespace ug
{
namespace bridge
{

class ClassHierarchy
{
public:
	void insert_class(const IExportedClass &c);
	ClassHierarchy *find_class(const char *name);

	string name;
	vector<ClassHierarchy> subclasses;
};

void GetClassHierarchy(ClassHierarchy &hierarchy, const bridge::Registry &reg);
const IExportedClass *FindClass(const char* classname);

}
}
