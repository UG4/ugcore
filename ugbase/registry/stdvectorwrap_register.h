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
	

template<typename T>
void RegStdVectorWrap(Registry &registry, std::string grp, const char *pname=NULL)
{		
	std::string name;
	if(pname == NULL)
	{
		name = "stdvector_";
		name.append(ClassNameProvider<T>::name());
	}
	else
		name = pname;
	registry.add_class_< std_vector_wrap<T> > (name, grp)
		.add_method("get", &std_vector_wrap<T>::index, "")
		.add_method("__tostring", &std_vector_wrap<T>::tostring, "")
		.add_method("size", &std_vector_wrap<T>::size, "");
}

} }
#endif	/* STDVECTORWRAP_REGISTER_H */

