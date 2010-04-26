//	object.h
//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com

#ifndef __H__UG__NODE_TREE__OBJECT__
#define __H__UG__NODE_TREE__OBJECT__

namespace ug{
namespace node_tree
{

////////////////////////////////////////////////////////////////////////
//	Object-Codes
///	ids associated with objects. Values shouldn't be unnessesarily high.
enum ObjectCode
{
	OC_INVALID = 0xFFFFFFFF,
	OC_OBJECT = 0,
//	nodes
	OC_NODE = 1,
	OC_GROUP_NODE,
	OC_BOXED_GROUP_NODE,
	OC_COLLISION_TREE_ROOT_NODE,
	OC_COLLISION_EDGES_NODE,
	OC_COLLISION_TRIANGLES_NODE,

	OC_NODES_END,

//	custom objects
	OC_CUSTOM_OBJECT
};


////////////////////////////////////////////////////////////////////////
//	Object
///	An Object serves as the base-class for most of the polymorphic node-tree objects
/**
 * Each derivative should have its own ObjectCode, with which it can
 * be uniquely identified.
 */
class Object
{
	public:
		virtual ~Object()	{}

		inline unsigned int getObjectCode()	{return m_objectCode;}

	protected:
		Object()	{};
		Object(const Object& obj)	{};

	protected:
		unsigned int m_objectCode;

};

}//	end of namespace node_tree
}//	end of namespace ug

#endif
