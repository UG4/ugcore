
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

}
}

#endif /* __H__UG__END_BOOST_LIST_H_ */
