// author: Martin Rupp

#include "class_helper.h"
#include "ug_script/ug_script.h"

namespace ug
{
namespace bridge
{
extern Registry& GetUGRegistry();

void ClassHierarchy::insert_class(const IExportedClass &c)
{
	const vector<const char *> *pNames = c.class_names();

	if(pNames == NULL) return;

	ClassHierarchy *base = this;

	for(vector<const char*>::const_reverse_iterator rit = pNames->rbegin(); rit < pNames->rend(); ++rit)
	{
		const char *thename = (*rit);
		size_t j;
		for(j=0; j<base->subclasses.size(); j++)
		{
			if(base->subclasses.at(j).name.compare(thename) == 0)
			{
				base = &base->subclasses.at(j);
				break;
			}
		}

		if(j == base->subclasses.size())
		{
			ClassHierarchy newclass;
			newclass.name = thename;
			base->subclasses.push_back(newclass);
		}
	}
}


ClassHierarchy *ClassHierarchy::find_class(const char *classname)
{
	if(name.compare(classname) == 0)
		return this;
	for(size_t i=0; i < subclasses.size(); i++)
	{
		ClassHierarchy *c = subclasses[i].find_class(classname);
		if(c) return c;
	}
	return NULL;
}

void GetClassHierarchy(ClassHierarchy &hierarchy, const bridge::Registry &reg)
{
	hierarchy.subclasses.clear();
	hierarchy.name = "UGBase";
	for(size_t i=0; i<reg.num_classes(); ++i)
		hierarchy.insert_class(reg.get_class(i));
}

const IExportedClass *FindClass(const char* classname)
{
	bridge::Registry &reg = ug::GetUGRegistry();
	for(size_t j=0; j<reg.num_classes(); ++j)
		if(strcmp(classname, reg.get_class(j).name()) == 0)
		{
			return &reg.get_class(j);
			break;
		}
	return NULL;
}

} // namespace bridge
} // namespace ug

