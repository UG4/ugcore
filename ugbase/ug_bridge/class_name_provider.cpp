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

void ExtractClassNameVec(std::vector<const char*>& names, const ClassNameNode& node, bool clearVec)
{
//	clear vector
	if(clearVec)
		names.clear();

//	add node name
	names.push_back(node.name().c_str());

//	add all base classes
	for(size_t i = 0; i < node.num_base_classes(); ++i)
		ExtractClassNameVec(names, node.base_class(i), false);
}

bool ClassNameTreeContains(const ClassNameNode& node, const char* name)
{
//	if the node is the name, return true
	if(strcmp(node.name().c_str(), name) == 0) return true;

//	else search in parents
	bool bContains = false;

	for(size_t i = 0; i < node.num_base_classes(); ++i)
		bContains |= ClassNameTreeContains(node.base_class(i), name);

//	return if found in parents
	return bContains;
}

bool ClassNameTreeWay(std::vector<size_t>& vWay, const ClassNameNode& node, const char* name)
{
//	if the node is the name, return true
	if(strcmp(node.name().c_str(), name) == 0) return true;

//	look in parents
	for(size_t i = 0; i < node.num_base_classes(); ++i)
	{
		if(ClassNameTreeWay(vWay, node.base_class(i), name))
		{
			vWay.push_back(i);
			return true;
		}
	}
//	return if found in parents
	return false;
}

std::map<std::pair<const ClassNameNode*, const ClassNameNode*>, void* (*)(void*)>
	ClassCastProvider::m_mmCast = std::map<std::pair<const ClassNameNode*, const ClassNameNode*>,  void* (*)(void*)> ();

}//	end of namespace
}//	end of namespace
