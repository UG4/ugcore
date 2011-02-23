// created by Martin Rupp 31.01.2011

#include <utility>

#ifndef __H__SORT_UTIL__
#define __H__SORT_UTIL__

template<typename TIndex, typename TValue>
struct SortStruct
{
	SortStruct() { }
	SortStruct(TIndex &i, TValue &v) : index(i), value(v) { }

	void set_value(TValue &val)
	{
		value = val;
	}

	TValue &get_value()
	{
		return value;
	}

	bool operator < (const SortStruct<TIndex, TValue> &other) const
	{
		return value < other.value;
	}

	TIndex index;
	TValue value;
};



template<typename TCompareValues>
class CompareIndicesByClass
{
public:
	CompareIndicesByClass(const TCompareValues &compareValues) : m_compareValues(compareValues) { }

	template<typename TIndex>
	bool operator () (const TIndex &i, const TIndex &j) const
	{
		return m_compareValues[i] < m_compareValues[j];
	}
private:
	const TCompareValues &m_compareValues;
};

template<typename TCompareValues, typename TCompare>
class CompareIndicesByClass2
{
public:
	CompareIndicesByClass2(const TCompareValues &compareValues, TCompare comp) : m_compareValues(compareValues), m_comp(comp) { }

	template<typename TIndex>
	bool operator () (const TIndex &i, const TIndex &j) const
	{
		return m_comp(m_compareValues[i], m_compareValues[j]);
	}
private:
	const TCompareValues &m_compareValues;
	TCompare m_comp;
};

////////////////////////////////////////////////////////////////////////////////////////////////
//	CompareIndicesBy
/**
 * Allows sorting of indices by their values in another array.
 * example:
 *  int v[] = {1, 3, 2};
	int pivot[] = {5, 2, 3, 1};
	sort(v, v+3, CompareIndicesBy(pivot));
	for(size_t i=0; i<3; i++)
		cout << v[i] << " = " << pivot[v[i]] << endl;
 */
template<typename TCompareValues>
inline CompareIndicesByClass<TCompareValues> CompareIndicesBy(const TCompareValues &values)
{
	return CompareIndicesByClass<TCompareValues>(values);
}

template<typename TCompareValues, typename TCompare>
inline CompareIndicesByClass2<TCompareValues, TCompare> CompareIndicesBy(const TCompareValues &values, TCompare comp)
{
	return CompareIndicesByClass2<TCompareValues, TCompare>(values, comp);
}


inline bool boolstrcmp(const char *a, const char *b)
{
	return strcmp(a, b) < 0;
}


#endif


