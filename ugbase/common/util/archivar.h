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

#ifndef __H__UG_archivar
#define __H__UG_archivar

#include <boost/static_assert.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include "../boost_serialization.h"
#include "detail/register_type_pair_functor.h"

namespace ug{


namespace detail{
namespace archivar{
	template <typename TArchive, typename TBase, typename TDerived>
	void CallArchiveOnDerivedClass (TArchive& ar, TBase& base, const char* name)
	{
		BOOST_STATIC_ASSERT((boost::is_base_of<TBase, TDerived>::value));
		TDerived& derived = dynamic_cast<TDerived&>(base);
		ar & make_nvp(name, derived);
	}
}
}

template <class TArchive, class TBase, class TPairSeq>
class Archivar
{
public:
	Archivar()
	{
		detail::RegisterTypePairFunctor<Archivar> func(this);
		boost::mpl::for_each<TPairSeq>(func);
	}

	template <class T>
	void archive(TArchive& ar, T& t)
	{
		archive(ar, t, "");
	}

	template <class T>
	void archive(TArchive& ar, T& t, const char* name)
	{
		std::string typeName(typeid(t).name());

		typename callback_map_t::iterator icallback = m_callbackMap.find(typeName);
		UG_COND_THROW(icallback == m_callbackMap.end(),
					  "Unregistered class used in 'Archivar::archive': " << typeName);

		icallback->second(ar, t, name);
	}

	template <class TDerived>
	void register_class(const char* name)
	{
		BOOST_STATIC_ASSERT((boost::is_base_of<TBase, TDerived>::value));
		std::string typeName(typeid(TDerived).name());
		m_callbackMap[typeName] =
				&detail::archivar::CallArchiveOnDerivedClass<TArchive, TBase, TDerived>;
	}

private:
	using archive_sig = void(*)(TArchive&, TBase&, const char* name);

	using callback_map_t = std::map<std::string, archive_sig>;
	callback_map_t	m_callbackMap;
};

}//	end of namespace

#endif	//__H__UG_archivar
