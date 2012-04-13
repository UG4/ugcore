//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m04 d26

#ifndef __H__UG__NODE_TREE__COLLSISION_ELEMENT_INFO__
#define __H__UG__NODE_TREE__COLLSISION_ELEMENT_INFO__

namespace ug{
namespace node_tree
{

////////////////////////////////////////////////////////////////////////
//	CollisionElementID
/**
 * This struct allows to specify a pointer or an integer-id that shall be
 * stored with each element (edge or triangle) in a collision-tree.
 *
 * Note that only a pointer or an integer may be stored not both at the
 * same time.
 *
 * If you construct an empty CollisionElementID, it will be set to
 * invalid by default.
 * Note that an integer value of -1 marks the CollisionElementID as invalid.
 */
struct CollisionElementID{
	CollisionElementID() : m_intID(-1)			{}
	CollisionElementID(void* p)	: m_ptrID(p)	{}
	CollisionElementID(int i)	: m_intID(i)	{}

	inline bool is_valid()	{return m_intID != -1;}
	
	bool operator == (const CollisionElementID& id) const {return m_ptrID == id.m_ptrID;}

	union{
		void* 	m_ptrID;
		int		m_intID;
	};
};

}//	end of namespace node_tree
}//	end of namespace ug

#endif

