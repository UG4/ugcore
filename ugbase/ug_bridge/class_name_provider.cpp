#include "class_name_provider.h"

namespace ug
{
namespace bridge
{

bool ClassNameVecContains(const std::vector<const char*>& names, const char* name)
{
	//  return true if pointers are equal
		for(size_t i = 0; i < names.size(); ++i)
			if(name == names[i]) return true;

	//  security fallback: if pointers not equal, compare also strings
		for(size_t i = 0; i < names.size(); ++i)
			if(strcmp(name, names[i]) == 0) return true;

	//  both comparisons fail. Name is not this class, nor on of its parents
		return false;
}

}//	end of namespace
}//	end of namespace
