/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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

#ifndef __H__UG__END_BOOST_LIST_H_
#define __H__UG__END_BOOST_LIST_H_


#include <boost/mpl/list.hpp>

/**
 * This structure is a helper structure to end an boost::mpl::list
 * \code
 * boost::mpl::list<int, double, char, end_boost_list>
 * \endcode
 * is the same as
 * \code
 * boost::mpl::list<int, double, char>
 * \endcode
 *
 * This comes in handy when you are using defines to enable/disable parts of the list:
 * \code
 * boost::mpl::list<
 * 	#ifdef USE_INT
 * 	  int,
 * 	#endif
 * 	#ifdef USE_DOUBLE
 * 	  double,
 * 	#endif
 * 	#ifdef USE_CHAR
 * 	  char,
 * 	#endif
 * 	end_boost_list
 * 	>
 * \endcode
 *  Without end_boost_list, we would need to add some other #ifdefs because
 *  of the extra Comma , at the end (imagine only USE_INT is defined).
 *
 *  you can use that list now as normal:
 *  static const bool isEmpty = boost::mpl::empty<List>::value;
	typename boost::mpl::if_c<isEmpty, ListProcessEnd, ListProcessNext>::type (reg,grp);
 */
namespace ug{
	struct end_boost_list{};
}

namespace boost
{
namespace mpl
{
/// specializing boost::mpl::empty so that an list<end_boost_list> appears to be empty
template<>
struct empty<list1<ug::end_boost_list> >
{
	enum { value = true };
};

/// lists that contain nothing but ug::end_boost_list from the start need also be regarded as empty
template<>
struct empty<list<ug::end_boost_list> >
{
	enum { value = true };
};

}
}

#endif
