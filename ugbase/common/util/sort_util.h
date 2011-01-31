// created by Martin Rupp 31.01.2011

#ifndef __H__SORT_UTIL__
#define __H__SORT_UTIL__

template<typename TIndex, typename TValue>
struct SortStruct
{
	TIndex index;
	TValue value;

	void set_value(TValue &val)
	{
		value = val;
	}

	bool operator < (const SortStruct<TIndex, TValue> &other) const
	{
		return value < other.value;
	}
};


template<typename TIndex, typename TValue>
struct SortStructIndirect
{
	TIndex index;
	TValue *value;

	void set_value(TValue &val)
	{
		value = &val;
	}

	bool operator < (const SortStructIndirect<TIndex, TValue> &other) const
	{
		return *value < *(other.value);
	}
};


template<typename TSortType, typename TIterator, typename TPivot>
void IndicesSortHelper(const TIterator &begin, const TIterator &end, std::vector<TPivot> &vElements)
{
	std::vector<TSortType> vSort;
	size_t N = end-begin;
	vSort.resize(N);

	size_t i=0;
	for(TIterator it = begin; it != end; ++it, ++i)
	{
		size_t index = *it;
		assert(index >= 0 && index < vElements.size());
		vSort[i].set_value(vElements[index]);
		vSort[i].index = index;
	}
	sort(vSort.begin(), vSort.end());
	i=0;
	for(TIterator it = begin; it != end; ++it, ++i)
		*it = vSort[i].index;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//	IndicesSort
/**
 * Allows sorting a range between two iterators begin and end which contain indices in a vector
 * whose elements are compareable. To do this, we create a list with copies of those elements
 */
template<typename TIterator, typename TPivot>
void IndicesSort(const TIterator &begin, const TIterator &end, std::vector<TPivot> &vElements)
{
	IndicesSortHelper<SortStruct<typename TIterator::value_type, TPivot> , TIterator, TPivot> (begin, end, vElements);
}

////////////////////////////////////////////////////////////////////////////////////////////////
//	IndicesSort
/**
 * Allows sorting a range between two iterators begin and end which contain indices in a vector
 * whose elements are compareable. To do this, we create a list pointers to those elements
 * use this if elements of vElements are large.
 */
template<typename TIterator, typename TPivot>
void IndicesSortIndirect(const TIterator &begin, const TIterator &end, std::vector<TPivot> &vElements)
{
	IndicesSortHelper<SortStructIndirect<typename TIterator::value_type, TPivot> , TIterator, TPivot>(begin, end, vElements);
}
#endif
