/* 
 * File:   stdvectorwrap_register.h
 * Author: mrupp
 *
 * Created on 30. Oktober 2012, 14:58
 */


#ifndef STDVECTORWRAP_REGISTER_H
#define	STDVECTORWRAP_REGISTER_H

#include "stdvectorwrap.h"
#include "registry.h"
#include <string>

namespace ug{ namespace bridge{
	

template<typename TValueType>
void RegStdVectorWrap(Registry &registry, std::string grp, const char *pname=NULL)
{		
	std::string name;
	if(pname == NULL)
	{
		name = "stdvector_";
		name.append(ClassNameProvider<TValueType>::name());
	}
	else
		name = pname;
	typedef std_vector_wrap<TValueType> T;
	registry.add_class_< T > (name, grp)
		.add_constructor()
		.template add_constructor<void (*)(size_t s)>("size")
		.template add_constructor<void (*)(size_t s, TValueType)>("size#initVal")
		.add_method("get", &T::index, "element at position 'index'", "index")
		.add_method("__tostring", &T::tostring, "")
		.add_method("size", &T::size, "size")
	
		.add_method("resize", static_cast<void (T::*)(size_t)>(&T::resize), "", "size")
		.add_method("resize", static_cast<void (T::*)(size_t, TValueType)>(&T::resize), "","size#initVal")
		.add_method("set", &T::set, "", "index#value")
		.add_method("push_back", &T::push_back, "", "value")
		.add_method("reserve", &T::reserve, "", "size")
		.add_method("clear", &T::clear, "", "")
		;
}

} }
#endif	/* STDVECTORWRAP_REGISTER_H */

