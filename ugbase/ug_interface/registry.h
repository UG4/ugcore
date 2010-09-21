//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m09 d16

#ifndef __H__UG__BINDINGS__REGISTRY__
#define __H__UG__BINDINGS__REGISTRY__

#include <vector>
#include "common/common.h"
#include "interface_base.h"

namespace ug
{
namespace interface
{

/**
 * The Registry singleton allows registration of global interface functions and
 * of interface objects.
 *
 * The Registry can be used to create bindings to different languages and programs, such
 * as lua-scripting, vrl-visual-programming or ProMesh
 */
class Registry
{
	public:
		static Registry& inst()
		{
			static Registry reg;
			return reg;
		}
		
	///	registers an object-type. Make sure that TObject derives from IObject.
		template <class TObject>
		bool register_object(){
		//todo:	check whether an object in the objects group
		//		has already been registered with the same name.
		//		If so - don't register the object and return false.
			m_objects.push_back(new TObject);
			return true;
		}

	///	registers a global function. Make sure that TFunc derives from IGlobalFunction.
		template <class TFunc>
		bool register_global_function(){
		//todo:	check whether a function in the functions group
		//		has already been registered with the same name.
		//		If so - don't register the function and return false.
			m_functions.push_back(new TFunc);
			return true;
		}
		
	///	Pass a prototype object to the register method.
	/**	Make sure that the prototype has been created with new.
	 *	This method is only useful if your object can not be created with
	 *	the default constructor or if you want to initialize it with special values.
	 *	Please note that you pass control over obj to the ObjectRegistry
	 *	when calling this method. This means that you neither should use nor
	 *	delete the object after you passed it to this method.*/
		bool register_object(IObject* obj){
		//todo:	check whether an object in the objects group
		//		has already been registered with the same name.
		//		If so - don't register the object and return false.
			m_objects.push_back(obj);
			return true;
		}
		
		size_t num_objects()				{return m_objects.size();}
		IObject* get_object(size_t index)	{return m_objects.at(index);}
		
		size_t num_functions()							{return m_functions.size();}
		IGlobalFunction* get_function(size_t index)	{return m_functions.at(index);}
		
	private:
		Registry();
		Registry(const Registry& reg)	{}
		~Registry()
		{
		//	delete all objects
			for(size_t i = 0; i < m_objects.size(); ++i)
				delete m_objects[i];
				
		//	delete all functions
			for(size_t i = 0; i < m_functions.size(); ++i)
				delete m_functions[i];
		}
		
	protected:
		std::vector<IObject*>			m_objects;
		std::vector<IGlobalFunction*>	m_functions;
};

}//	end of namespace
}//	end of namespace

#endif
