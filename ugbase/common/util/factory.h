/*
 * Copyright (c) 2016:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG_factory
#define __H__UG_factory

#include <boost/static_assert.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include "../boost_serialization.h"
#include "detail/register_type_pair_functor.h"

namespace ug{

namespace detail{
namespace factory{
	///	used internally to construct derived classes of a given type
	template <class TBase, class TDerived>
	TBase* DerivedClassFactory ()
	{
		return new TDerived;
	}
}
}


///	A factory class which creates instances given a class-name
/** Given a base-class 'TBase' and an optional boost::mpl sequence of pairs of
 * derived classes and strings 'TPairSeq', the ug::Factory class allows for the
 * creation of instances ofderived classes given a name. Additionally to the
 * classes contained in 'TPairSeq', one may register derived classes at runtime
 * through the 'register_class' method.
 *
 * Besides creation of registered classes through the 'create' method,
 * one can iterate over all registered classes or obtain a registered class-name
 * given a pointer to the base class.
 *
 * \note	This factory class is somewhat limited compared to e.g. to boost::factory and
 * 			is not meant as a replacement. It's main purpose is to allow iteration
 *			over registered classes and to perform automatic registration through
 *			a type list.
 */
template <class TBase, class TPairSeq = boost::mpl::vector<> >
class Factory
{
public:
	Factory()
	{
		detail::RegisterTypePairFunctor<Factory> func(this);
		boost::mpl::for_each<TPairSeq>(func);
	}

	SmartPtr<TBase> create(const std::string& className)
	{
		return create(className.c_str());
	}

	SmartPtr<TBase> create(const char* className)
	{
		return SmartPtr<TBase>(create_raw(className));
	}

	TBase* create_raw(const std::string& className)
	{
		return create_raw(className.c_str());
	}

	TBase* create_raw(const char* className)
	{
		typename class_map_t::iterator iclass = m_classMap.find(std::string(className));
		UG_COND_THROW(iclass == m_classMap.end(),
					  "Unregistered class-name used in 'Factory::create_raw': " << className);
		return iclass->second.factory();
	}

	template <class TDerived>
	void register_class(const char* name)
	{
		// UG_LOG("registering derived class: " << className << std::endl);
		BOOST_STATIC_ASSERT((boost::is_base_of<TBase, TDerived>::value));

		std::string className(name);

		m_classMap[className] = ClassInfo(
				name,
				&detail::factory::DerivedClassFactory<TBase, TDerived>);
		
		std::string typeName(typeid(TDerived).name());
		m_classNameMap[typeName] = className;

		m_classNames.push_back(className);
	}

	size_t num_classes() const								{return m_classNames.size();}
	const std::string& class_name(const size_t i) const		{return m_classNames.at(i);}
	const std::string& class_name(const TBase& cls) const
	{
		std::string typeName = typeid(cls).name();
		class_name_map_t::const_iterator i = m_classNameMap.find(typeName);
		if(i != m_classNameMap.end())
			return i->second;
		static const std::string defaultName ("");
		return defaultName;
	}

private:
	typedef TBase* (*factory_sig)();

	struct ClassInfo {
		ClassInfo ()	{}
		ClassInfo (const char* _name,
				   factory_sig _factory) :
			name(_name),
			factory(_factory)
		{}

		const char*			name;
		factory_sig			factory;
	};

	typedef std::map<std::string, ClassInfo>	class_map_t;
	typedef std::map<std::string, std::string>	class_name_map_t;
	class_map_t					m_classMap;///< contains ClassInfo objects accessible by class-name.
	class_name_map_t			m_classNameMap;///< key: type-name, value: class-names
	std::vector<std::string>	m_classNames;
};

}//	end of namespace

#endif	//__H__UG_factory
