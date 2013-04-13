// created by Martin Rupp 31.01.2011

#include <utility>

#ifndef __H__SORT_UTIL__
#define __H__SORT_UTIL__

/// \addtogroup ugbase_common_util
/// \{

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

/**
 * Allows sorting of indices by their values in another array with a comp function
 */
template<typename TCompareValues, typename TCompare>
inline CompareIndicesByClass2<TCompareValues, TCompare> CompareIndicesBy(const TCompareValues &values, TCompare comp)
{
	return CompareIndicesByClass2<TCompareValues, TCompare>(values, comp);
}

/**
 * example:
 * ind = GetSortedIndices(v);
 * then v[ ind[0] ] <= v[ ind[1] ] <= ... <= v [ ind[v.size() -1] ]; ,
 * @param values array where indices should be sorted after
 * @return the indices list sorted with respect to the ordering in values
 */
template<typename TCompareValues>
std::vector<size_t> GetSortedIndices(const TCompareValues &values)
{
	std::vector<size_t> indices(values.size());
	for(size_t i=0; i<indices.size(); i++) indices[i] = i;
	std::sort(indices.begin(), indices.end(), CompareIndicesByClass<TCompareValues>(values));
	return indices;
}

/**
 * example:
 * ind = GetSortedIndices(v);
 * then v[ ind[0] ] <= v[ ind[1] ] <= ... <= v [ ind[v.size() -1] ]; ,
 * @param values array where indices should be sorted after
 * @param comp compare function
 * @return the indices list sorted with respect to the ordering in values (with compare function comp)
 */
template<typename TCompareValues, typename TCompare>
std::vector<size_t> GetSortedIndices(const TCompareValues &values, TCompare comp)
{
	std::vector<size_t> indices(values.size());
	for(size_t i=0; i<indices.size(); i++) indices[i] = i;
	std::sort(indices.begin(), indices.end(), CompareIndicesByClass2<TCompareValues, TCompare>(values, comp));
	return indices;
}

inline bool boolstrcmp(const char *a, const char *b)
{
	return strcmp(a, b) < 0;
}

// end group ugbase_common_util
/// \}

#endif


