/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

	bool operator < (const SortStruct &other) const
	{
		return value < other.value;
	}

	TIndex index;
	TValue value;
};



template<typename TCompareValues, bool bReverse>
class CompareIndicesByClass
{
public:
	CompareIndicesByClass(const TCompareValues &compareValues) : m_compareValues(compareValues) { }

	template<typename TIndex>
	bool operator () (const TIndex &i, const TIndex &j) const
	{
		if(bReverse)
			return m_compareValues[i] > m_compareValues[j];
		else
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
inline CompareIndicesByClass<TCompareValues, false> CompareIndicesBy(const TCompareValues &values)
{
	return CompareIndicesByClass<TCompareValues, false>(values);
}

template<typename TCompareValues>
inline CompareIndicesByClass<TCompareValues, true> CompareIndicesReversedBy(const TCompareValues &values)
{
	return CompareIndicesByClass<TCompareValues, true>(values);
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
 * std::vector<size_t> v; // = { 45, 3, 6, 4, 9};
 * ind = GetSortedIndices(ind, v);
 * then v[ ind[0] ] <= v[ ind[1] ] <= ... <= v [ ind[v.size() -1] ]; ,
 * and here ind = 1, 3, 2, 4, 0, because e.g. v[ind[0]] = v[1] = 3.
 * @param values array where indices should be sorted after
 * @return the indices list sorted with respect to the ordering in values
 */
template<typename TCompareValues>
void GetSortedIndices(std::vector<size_t> &indices, const TCompareValues &values, bool bReverse=false)
{
	indices.resize(values.size());
	for(size_t i=0; i<indices.size(); i++) indices[i] = i;
	if(bReverse)
		std::sort(indices.begin(), indices.end(), CompareIndicesByClass<TCompareValues, true>(values));
	else
		std::sort(indices.begin(), indices.end(), CompareIndicesByClass<TCompareValues, false>(values));
}

/**
 * example:
 * std::vector<const char*> v; // some strings
 * GetSortedIndices(ind, v, boolstrcmp);
 * then v[ ind[0] ] <= v[ ind[1] ] <= ... <= v [ ind[v.size() -1] ]; ,
 * @param indices array with indices
 * @param values array where indices should be sorted after
 * @param comp compare function
 */
template<typename TCompareValues, typename TCompare>
void GetSortedIndices(std::vector<size_t> &indices, const TCompareValues &values, TCompare comp)
{
	indices.resize(values.size());
	for(size_t i=0; i<indices.size(); i++) indices[i] = i;
	std::sort(indices.begin(), indices.end(), CompareIndicesByClass2<TCompareValues, TCompare>(values, comp));
}

template<typename TCompareValues, typename TCompare>
std::vector<size_t> GetSortedIndices(const TCompareValues &values, TCompare comp)
{
	std::vector<size_t> indices(values.size());
	GetSortedIndices(indices, values, comp);
	return indices;
}

template<typename TCompareValues>
std::vector<size_t> GetSortedIndices(const TCompareValues &values)
{
	std::vector<size_t> indices(values.size());
	GetSortedIndices(indices, values);
	return indices;
}



inline bool boolstrcmp(const char *a, const char *b)
{
	return strcmp(a, b) < 0;
}

inline bool boolstrcmpReversed(const char *a, const char *b)
{
	return strcmp(a, b) > 0;
}

// end group ugbase_common_util
/// \}

#endif


