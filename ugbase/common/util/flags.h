// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Jul 31, 2013 (d,m,y)

#ifndef __H__UG__flags__
#define __H__UG__flags__

namespace ug{

///	Helps maintaining, activating and deactivating a set of flags from an enum.
/**	Given an enum
 * \code
 * enum SomeEnum{
 * 	E1 = 0,
 * 	E2 = 1,
 * 	E3 = 1 << 1,
 * 	E4 = 1 << 2
 * };
 *
 * \endcode
 * One can use the Flag class as follows:
 *
 * \code
 * typedef Flag<SomeEnum> SomeFlag;
 * ...
 * SomeFlag f1(E2), f2(E4);
 * SomeFlag f3(f1 | f2);
 * if(f3.contains(E2))
 * 	f3.remove(E2);
 * ...
 * \endcode
 *
 */
template <class TEnum, class TStorageType = unsigned int, TStorageType defaultValue = 0>
class Flag{
	public:
		Flag()					: m_value(defaultValue)	{}
		Flag(TStorageType flag)	: m_value(flag)			{}
		Flag(const Flag& flag)	: m_value(flag.m_value)	{}

		bool contains(TStorageType flag) const	{return (m_value & flag) == flag;}
		bool contains(const Flag& flag) const	{return (m_value & flag.m_value) == flag.m_value;}

		bool partially_contains(TStorageType flag) const	{return (m_value & flag) != 0;}
		bool partially_contains(const Flag& flag) const		{return (m_value & flag.m_value) != 0;}

		Flag& set(TStorageType flag)			{m_value = flag; return *this;}
		Flag& add(TStorageType flag)			{m_value |= flag; return *this;}
		Flag& remove(TStorageType flag)			{m_value &= (~flag); return *this;}

		Flag operator& (const Flag& flag) const	{return Flag(m_value & flag.m_value);}
		Flag operator&= (const Flag& flag)		{m_value &= flag.m_value; return *this;}
		Flag operator| (const Flag& flag) const	{return Flag(m_value | flag.m_value);}
		Flag operator|= (const Flag& flag)		{m_value |= flag.m_value; return *this;}
		Flag operator= (const Flag& flag)		{m_value = flag.m_value; return *this;}
		Flag operator= (TStorageType val)		{m_value = val; return *this;}

		TStorageType operator()() const			{return m_value;}
		TStorageType get() const				{return m_value;}

		bool operator== (const Flag& flag) const	{return m_value == flag.m_value;}
		bool operator== (TStorageType val) const	{return m_value == val;}

		bool operator!= (const Flag& flag) const	{return m_value != flag.m_value;}
		bool operator!= (TStorageType val) const	{return m_value != val;}

	private:
		TStorageType	m_value;
};

}// end of namespace

#endif
