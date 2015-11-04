#ifndef __H__UG__generic_grid_object_iterator__
#define __H__UG__generic_grid_object_iterator__

namespace ug
{

////////////////////////////////////////////////////////////////////////////////////////////////
//	GenericGridObjectIterator
///	Use this class as a tool to create iterators to your own geometric objects.
template <class TValue, class TBaseIterator>
class GenericGridObjectIterator : public TBaseIterator
{
	friend class Grid;
	template <class TIterDest, class TIterSrc> friend TIterDest iterator_cast(const TIterSrc& iter);

	public:
		typedef TValue	value_type;

	public:
		GenericGridObjectIterator()	{}

		GenericGridObjectIterator(const GenericGridObjectIterator& iter) :
			TBaseIterator(iter)	{}

	///	note that the * operator is read only.
		inline TValue operator* () const	{return static_cast<TValue>(TBaseIterator::operator*());}

	protected:
		GenericGridObjectIterator(const TBaseIterator& iter) :
			TBaseIterator(iter)	{}
};

////////////////////////////////////////////////////////////////////////////////////////////////
//	ConstGenericGridObjectIterator
///	Use this class as a tool to create const_iterators to your own geometric objects.
template <class TValue, class TBaseIterator, class TConstBaseIterator>
class ConstGenericGridObjectIterator : public TConstBaseIterator
{
	friend class Grid;
	template <class TIterDest, class TIterSrc> friend TIterDest iterator_cast(const TIterSrc& iter);

	public:
		typedef TValue	value_type;

	public:
		ConstGenericGridObjectIterator()	{}

		ConstGenericGridObjectIterator(const ConstGenericGridObjectIterator& iter) :
			TConstBaseIterator(iter)	{}

		ConstGenericGridObjectIterator(const GenericGridObjectIterator<TValue, TBaseIterator>& iter) :
			TConstBaseIterator(iter)	{}

	///	note that the * operator is read only.
		inline TValue operator* () const	{return static_cast<TValue>(TConstBaseIterator::operator*());}

	protected:
		ConstGenericGridObjectIterator(const TBaseIterator& iter) :
			TConstBaseIterator(iter)	{}

		ConstGenericGridObjectIterator(const TConstBaseIterator& iter) :
			TConstBaseIterator(iter)	{}
};

////////////////////////////////////////////////////////////////////////
//	iterator_cast
///	You should avoid casting whenever possible!
template <class TIterDest, class TIterSrc>
inline TIterDest
iterator_cast(const TIterSrc& iter)
{
	return TIterDest(iter);
}

}//	end of namespace

#endif
