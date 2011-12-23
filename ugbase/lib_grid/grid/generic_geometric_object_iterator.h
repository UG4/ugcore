// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 22.12.2011 (m,d,y)

#ifndef __H__UG__generic_geometric_object_iterator__
#define __H__UG__generic_geometric_object_iterator__

namespace ug
{

////////////////////////////////////////////////////////////////////////////////////////////////
//	GenericGeometricObjectIterator
///	Use this class as a tool to create iterators to your own geometric objects.
template <class TValue, class TBaseIterator>
class GenericGeometricObjectIterator : public TBaseIterator
{
	friend class Grid;
	template <class TIterDest, class TIterSrc> friend TIterDest iterator_cast(const TIterSrc& iter);

	public:
		typedef TValue	value_type;

	public:
		GenericGeometricObjectIterator()	{}

		GenericGeometricObjectIterator(const GenericGeometricObjectIterator& iter) :
			TBaseIterator(iter)	{}

	///	note that the * operator is read only.
		inline TValue operator* () const	{return static_cast<TValue>(TBaseIterator::operator*());}

	protected:
		GenericGeometricObjectIterator(const TBaseIterator& iter) :
			TBaseIterator(iter)	{}
};

////////////////////////////////////////////////////////////////////////////////////////////////
//	ConstGenericGeometricObjectIterator
///	Use this class as a tool to create const_iterators to your own geometric objects.
template <class TValue, class TBaseIterator, class TConstBaseIterator>
class ConstGenericGeometricObjectIterator : public TConstBaseIterator
{
	friend class Grid;
	template <class TIterDest, class TIterSrc> friend TIterDest iterator_cast(const TIterSrc& iter);

	public:
		typedef TValue	value_type;

	public:
		ConstGenericGeometricObjectIterator()	{}

		ConstGenericGeometricObjectIterator(const ConstGenericGeometricObjectIterator& iter) :
			TConstBaseIterator(iter)	{}

	///	note that the * operator is read only.
		inline TValue operator* () const	{return static_cast<TValue>(TConstBaseIterator::operator*());}

	protected:
		ConstGenericGeometricObjectIterator(const TBaseIterator& iter) :
			TConstBaseIterator(iter)	{}

		ConstGenericGeometricObjectIterator(const TConstBaseIterator& iter) :
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
