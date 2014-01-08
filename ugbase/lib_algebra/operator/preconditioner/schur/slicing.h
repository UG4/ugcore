/*
 * slicing.h
 *
 *  Created on: 08.01.2014
 *      Author: anaegel
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__SCHUR_SLICING_H_
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__SCHUR_SLICING_H_


#ifdef UG_PARALLEL

#include <iostream>
#include <sstream>
#include <string>
#include <set>

#include "lib_algebra/parallelization/parallelization.h"
#include "pcl/pcl.h"

#include "common/log.h"


namespace ug{

extern DebugID SchurDebug;

class SlicingOp
{

};

class SlicingData{

public:
	enum slice_desc_type {SD_INNER=0, SD_SKELETON, SLICE_DESC_SIZE};
	typedef std::vector<slice_desc_type> slice_desc_type_vector;
	//typedef std::set<int> slice_desc_set;
	typedef std::vector<int> slice_desc_set;

protected:
	slice_desc_type_vector m_slice_types;
	slice_desc_set m_slice_set[(int) SLICE_DESC_SIZE]; //!< mapping islice -> iglobal


public:
	/// constructor
	SlicingData(const slice_desc_type_vector &types)
	: m_slice_types(types)
	{
		auto_fill_sets();
	}

	/// copy types
	void set_types(const slice_desc_type_vector &types)
	{ m_slice_types = types; }

	/// auto fill for sets
	/// assigns every i=0.. m_slice_types.size()-1 to exactly one set
	void auto_fill_sets()
	{
		const size_t ntypes = m_slice_types.size();

		for (size_t i=0; i< ntypes; ++i)
		{
			slice_desc_type tt = get_type(i);
			slice_desc_set &set =slice(tt);

			// if sd is ordered, then this is constant
			//set.insert(set.end(), i);
			set.push_back(i);
		}

		UG_DLOG(SchurDebug, 5,"SlicingData::auto_fill_sets:" << ntypes << " "<< slice(SD_INNER).size() << " "<< slice(SD_SKELETON).size() << std::endl);

		slice_desc_set::const_iterator it;

		UG_DEBUG_BEGIN(SchurDebug, 5)

		{
			UG_LOG("Skeleton:");
			const slice_desc_set &myset=slice(SD_SKELETON);
			for (it=myset.begin(); it!=myset.end(); ++it)  UG_LOG(*it << " ");
		}

		{
			UG_LOG("\nInner:");
	    	const slice_desc_set &myset=slice(SD_INNER);
	    	for (it=myset.begin(); it!=myset.end(); ++it) UG_LOG(*it << " ");
		}
		UG_DEBUG_END(SchurDebug, 5)

	}


	template<class VT>
	SmartPtr<VT> slice_clone_without_values(const VT &full_src, slice_desc_type type) const
	{
		const slice_desc_set &slice_desc = slice(type);
		// SmartPtr<VT> slice_clone = full_src.clone_without_values();
		//slice_clone->resize(slice_desc.size());

		SmartPtr<VT> slice_clone = new VT(slice_desc.size());

		// convert layouts (vector->slice)
		SmartPtr<AlgebraLayouts> slice_layouts = new AlgebraLayouts(*full_src.layouts());
		replace_indices_in_layout(type, slice_layouts->master());
		replace_indices_in_layout(type, slice_layouts->slave());
		slice_clone->set_layouts(slice_layouts);
		//UG_LOG(*slice_layouts);

		return slice_clone;

	}

    /// copy: slice of vector -> small vector
	template<class VT>
	void get_vector_slice(const VT &full_src, slice_desc_type desc, VT &small_dst) const
	{
		const slice_desc_set &slice_desc = slice(desc);
		small_dst.resize(slice_desc.size());
		slice_desc_set::const_iterator elem = slice_desc.begin();
		for (size_t i=0; i<slice_desc.size(); ++i, ++elem)
			small_dst[i] = full_src[*elem];
	}



	/// copy: small vector -> slice of a vector
	template<class VT>
	void set_vector_slice(const VT &small_src, VT &full_dst, slice_desc_type desc) const
	{
		const slice_desc_set &slice_desc = slice(desc);
		slice_desc_set::const_iterator elem = slice_desc.begin();
		for (size_t i=0; i<slice_desc.size(); ++i, ++elem)
				full_dst[*elem] = small_src[i];
	}

	 /// copy: slice of vector -> small vector
		template<class VT>
		void subtract_vector_slice(const VT &full_src, slice_desc_type desc, VT &small_dst) const
		{
				const slice_desc_set &slice_desc = slice(desc);
				small_dst.resize(slice_desc.size());
				slice_desc_set::const_iterator elem = slice_desc.begin();
				for (size_t i=0; i<slice_desc.size(); ++i, ++elem)
					small_dst[i] -= full_src[*elem];
		}

		/// copy: small vector -> slice of a vector
			template<class VT>
			void subtract_vector_slice(const VT &small_src, VT &full_dst, slice_desc_type desc) const
			{
				const slice_desc_set &slice_desc = slice(desc);
				slice_desc_set::const_iterator elem = slice_desc.begin();
				for (size_t i=0; i<slice_desc.size(); ++i, ++elem)
						full_dst[*elem] -= small_src[i];
			}


	// Extracts a slice from a (full) matrix
	template<class MT>
	void get_matrix(const MT &A, slice_desc_type row_type, slice_desc_type col_type, MT &Aslice) const
	{
		const slice_desc_set &row_slice = slice(row_type);
		const slice_desc_set &col_slice = slice(col_type);
		UG_DLOG(SchurDebug, 5, "SlicingData::get_matrix:" << row_slice.size() << "x" << col_slice.size()<< std::endl)
		Aslice.resize_and_clear(row_slice.size(), col_slice.size());

		int ii=0;
		for (slice_desc_set::const_iterator elem = row_slice.begin();
			elem!=row_slice.end(); ++elem, ++ii)
		{
			const int i = *elem; // global index
			for(typename MT::const_row_iterator it = A.begin_row(i);
				it != A.end_row(i); ++it)

			{
				const int j=it.index();
				int jj;
				// if (get_type(j)!=col_type) continue;
				if (find_index(col_type, j, jj))
				{
					Aslice(ii, jj) = it.value();
				}
			 }
		}

	}

	size_t get_num_elems(slice_desc_type type) const
	{return slice(type).size();}

protected:

	/// returns type for a global index
	slice_desc_type get_type(int index)
	{return m_slice_types[index];}


	/// returns the set of global indices for a given type
	const slice_desc_set &slice(slice_desc_type type) const
	{return m_slice_set[type];}

	slice_desc_set &slice(slice_desc_type type)
	{return m_slice_set[type];}


	/// returns local index for a global index
	bool find_index(slice_desc_type type, int gindex, int &index) const
	{
		bool found=false;

		const slice_desc_set &myset=slice(type);
		index = myset.size();

		//slice_desc_set::const_iterator it = myset.find(gindex);
		slice_desc_set::const_iterator it = lower_bound(myset.begin(), myset.end(), gindex);
		if (it != myset.end() && *it == gindex) {
			//index =// *it;
			index = std::distance(myset.begin(), it);
			found = true;
		}
	//	if (found && index >=myset.size()) {
		UG_ASSERT( (!found || index<myset.size()) , "Invalid index found!");
	//	}
		return found;
	}

	void replace_indices_in_layout(slice_desc_type type, IndexLayout &il) const
	{
		IndexLayout::iterator iter;
		for (iter = il.begin(); iter!=il.end(); ++iter)
		{
			// iterate over interfaces
			// i.e., pcl::OrderedInterface<size_t, std::vector>
			IndexLayout::Interface &interf=il.interface(iter);

			IndexLayout::Interface::iterator eiter;
			for (eiter = interf.begin(); eiter!=interf.end(); ++eiter)
			{
				// replace elements in interface
				size_t &elem = interf.get_element(eiter);
				int newind;

				bool found=find_index(type, elem, newind);
				UG_ASSERT(found, "SlicingData:: Did not find index???");
				interf.get_element(eiter) = newind;

			}
		}
	}

};

}

#endif /* UG_PARALLEL */
#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__SCHUR_SLICING_H_ */