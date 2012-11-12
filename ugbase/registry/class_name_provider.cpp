#include "class_name_provider.h"

namespace ug
{
namespace bridge
{

ClassNameNode::ClassNameNode()
{
	m_name.clear();
}

void ClassNameNode::set_name(const std::string& name)
{
//	set name
	m_name = std::string(name);

//	check size
	if(m_name.empty())
		UG_THROW_REGISTRY_MSG("Name must be longer than 0 characters.");
}

void ClassNameNode::add_base_class(const ClassNameNode& node)
{
	std::vector<const ClassNameNode*>::iterator it
		= std::find(m_vBaseClass.begin(), m_vBaseClass.end(), &node);
	if(it == m_vBaseClass.end())
		m_vBaseClass.push_back(&node);
}

bool ClassNameNode::named() const
{
//	if string is empty, then has not been named
	if(empty()) return false;

//	if string starts with '[' then has only been predeclared
	if(m_name.c_str()[0] == '[') return false;

//	has been named
	return true;
}


bool ClassNameVecContains(const std::vector<const char*>& names, const std::string& name)
{
	//  return true if pointers are equal
		for(size_t i = 0; i < names.size(); ++i)
			if(name.c_str() == names[i]) return true;

	//  security fallback: if pointers not equal, compare also strings
		for(size_t i = 0; i < names.size(); ++i)
			if(name == names[i]) return true;

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

bool ClassNameTreeContains(const ClassNameNode& node, const std::string& name)
{
//	if the node is the name, return true
	if(node.name() == name) return true;

//	else search in parents
	bool bContains = false;

	for(size_t i = 0; i < node.num_base_classes(); ++i)
		bContains |= ClassNameTreeContains(node.base_class(i), name);

//	return if found in parents
	return bContains;
}

bool ClassNameTreeWay(std::vector<size_t>& vWay, const ClassNameNode& node,
                      const std::string& name)
{
//	if the node is the name, return true
	if(node.name() == name) return true;

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

void* ClassCastProvider::
cast_to_base_class(void* pDerivVoid, const ClassNameNode*& node, const std::string& baseName)
{
//	find way to base class
	std::vector<size_t> vWay;
	if(!ClassNameTreeWay(vWay, *node, baseName))
	{
		UG_ERR_LOG("ERROR in ClassCastProvider::cast_to_base_class: Request"
				" to cast from derived class '"<< node->name()<<"' to "
				" base class '"<<baseName<<"', but no such base class in"
				" registered class hierarchy.");
		throw new UGError_ClassCastFailed(node->name(), baseName);
	}

	void* currPtr = pDerivVoid;
	const ClassNameNode* pCurrNode = node;

//	cast all the way down
	while(!vWay.empty())
	{
	//	get base class to cast to
		const ClassNameNode* pBaseClassNode = &pCurrNode->base_class(vWay.back());

	//	get name pair
		std::pair<const ClassNameNode*, const ClassNameNode*> namePair(pBaseClassNode, pCurrNode);

	//	find in map
		std::map<std::pair<const ClassNameNode*, const ClassNameNode*>, CastFunc>::iterator it;
		it = m_mmCast.find(namePair);

		if(it == m_mmCast.end())
		{
			UG_ERR_LOG("ERROR in ClassCastProvider::cast_to_base_class:"
					" Request intermediate cast from derived class '" <<
					pCurrNode->name() <<"' to direct base class '"
					<<pBaseClassNode->name()<<"', but no such cast "
					" function registered.");
			throw new UGError_ClassCastFailed(node->name(), baseName);
		}

	//	get cast function
		CastFunc pCastFunc = it->second;

	//	cast
		currPtr = (*pCastFunc)(currPtr);

	//	set node to base class
		pCurrNode = pBaseClassNode;

	//	pop way
		vWay.pop_back();
	}

//	write current node on exit
	node = pCurrNode;

//	return current pointer
	return currPtr;
}


const void* ClassCastProvider::
cast_to_base_class(const void* pDerivVoid, const ClassNameNode*& node, const std::string& baseName)
{
	return const_cast<const void*>(cast_to_base_class(const_cast<void*>(pDerivVoid), node, baseName));
}




std::map<std::pair<const ClassNameNode*, const ClassNameNode*>, void* (*)(void*)>
	ClassCastProvider::m_mmCast = std::map<std::pair<const ClassNameNode*, const ClassNameNode*>,  void* (*)(void*)> ();

}//	end of namespace
}//	end of namespace
