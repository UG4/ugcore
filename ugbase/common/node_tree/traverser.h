#ifndef __H__UG__NODE_TREE__TRAVERSER__
#define __H__UG__NODE_TREE__TRAVERSER__

#include <vector>
#include "object.h"
#include "node.h"

namespace ug{
namespace node_tree
{

class Traverser;
class GroupNode;
class BoxedGroupNode;

////////////////////////////////////////////////////////////////////////
//	Traverser
///	Derivates of a Traverser can be used to traverse a scenegraph.
/**
*/
class Traverser
{
	public:
		Traverser();
		virtual ~Traverser();

	protected:
		void apply(SPNode& node);

		template<typename HandlerType>
		void register_handler_function(unsigned int oc, HandlerType func);

		void traverse_object(Object* obj);

		virtual void handle_group(GroupNode* group);
		virtual void handle_boxed_group(BoxedGroupNode* boxedGroup);

	private:
		bool handler_function_registered(unsigned int oc);

	private:
		typedef void (Traverser::*HandlerFunc)(Object* obj);
		std::vector<HandlerFunc>	m_vHandlerFuncs;
};


template<typename HandlerType>
void Traverser::register_handler_function(unsigned int oc, HandlerType func)
{
//	make sure that there is enough space
	if(oc >= m_vHandlerFuncs.size())
		m_vHandlerFuncs.resize(oc+1, 0);

	m_vHandlerFuncs[oc] = (HandlerFunc)func;
}


}//	end of namespace node_tree
}//	end of namespace ug

#endif
