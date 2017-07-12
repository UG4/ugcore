/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
 * Author: Arne Nägel
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

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__SCHUR_SLICING_H_
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__SCHUR_SLICING_H_




#include <iostream>
#include <sstream>
#include <string>
#include <set>

#ifdef UG_PARALLEL
#include "lib_algebra/parallelization/parallelization.h"
#include "pcl/pcl.h"
#endif

#include "common/log.h"


namespace ug{


//! todo: replace DebugID
extern DebugID SchurDebug;


/*
 *
 * Allows splitting index sets for vectors/operator into different slices,
 * e.g., for schur complement, component-wise splitting, etc.
 *
 *
 * **/

template <class TVec, size_t N>
class SlicingData{

public:
	typedef std::vector<int> slice_desc_set;
	typedef TVec slice_desc_type_vector;
	typedef typename TVec::value_type slice_desc_type;

protected:
	bool m_valid;
	slice_desc_type_vector m_slice_types;
	slice_desc_set m_slice_set[N]; //!< N mappings: islice -> iglobal


public:
	//! Constructor

	SlicingData() : m_valid(false)
	{}

	/*! Builds index mappings based on types */
	SlicingData(const slice_desc_type_vector &types)
	: m_valid(true), m_slice_types(types)
	{
		reset_set_mappings();
	}

	/// copy types
	void set_types(const slice_desc_type_vector &types, bool bClear=false)
	{
		m_slice_types = types;
		reset_set_mappings();
	}

	bool is_valid()
	{ return m_valid;}


protected:

	void clear_set_mappings()
	{
		const size_t ntypes = m_slice_types.size();
		for (size_t i=0; i< ntypes; ++i)
		{
			slice_desc_type tt = get_type(i);
				slice_desc_set &set =slice(tt);
				set.clear();
		}
	}

	/// auto fill for sets
	/// assigns every index i=0.. m_slice_types.size()-1 to exactly one set
	void fill_set_mappings()
	{
		const size_t ntypes = m_slice_types.size();
		for (size_t i=0; i< ntypes; ++i)
		{
			slice_desc_type tt = get_type(i);
			slice_desc_set &set =slice(tt);

			set.push_back(i);
		}

		UG_DLOG(SchurDebug, 5,"SlicingData::fill_set_mappings:" << ntypes);
		//UG_LOG("SlicingData::fill_set_mappings:" << ntypes);
		for (size_t i=0; i<N; ++i)
			UG_DLOG(SchurDebug, 0, "  " << m_slice_set[i].size());
			//UG_LOG("  " << slice(get_type(i)).size());
		UG_DLOG(SchurDebug, 0, 	std::endl);
		//UG_LOG(std::endl);


		m_valid = true;
	}


	void reset_set_mappings()
	{
		clear_set_mappings();
		fill_set_mappings();
	}
public:




    /// copy: slice of vector -> small vector
	template<class VT>
	void get_vector_slice(const VT &full_src, slice_desc_type desc, VT &small_dst) const
	{

		const slice_desc_set &slice_desc = slice(desc);

		// UG_DLOG("get_vector_slice:" << slice_desc.size());
		small_dst.resize(slice_desc.size());
		slice_desc_set::const_iterator elem = slice_desc.begin();
		for (size_t i=0; i<slice_desc.size(); ++i, ++elem)
			{ small_dst[i] = full_src[*elem]; }
	}



	/// copy: small vector -> slice of a vector
	template<class VT>
	void set_vector_slice(const VT &small_src, VT &full_dst, slice_desc_type desc) const
	{
		const slice_desc_set &slice_desc = slice(desc);
		//UG_LOG("get_vector_slice:" << slice_desc.size());
		slice_desc_set::const_iterator elem = slice_desc.begin();
		for (size_t i=0; i<slice_desc.size(); ++i, ++elem)
			{ full_dst[*elem] = small_src[i]; }
	}

	 /// add: slice of vector -> small vector
	template<class VT>
	void add_vector_slice(const VT &full_src, slice_desc_type desc, VT &small_dst, double sigma=1.0) const
	{
		const slice_desc_set &slice_desc = slice(desc);
		small_dst.resize(slice_desc.size());
		slice_desc_set::const_iterator elem = slice_desc.begin();
		for (size_t i=0; i<slice_desc.size(); ++i, ++elem)
			{ small_dst[i] += sigma*full_src[*elem]; }
	}

	/// add: small vector -> slice of a vector
		template<class VT>
		void add_vector_slice(const VT &small_src, VT &full_dst, slice_desc_type desc, double sigma=1.0) const
		{
			const slice_desc_set &slice_desc = slice(desc);
			slice_desc_set::const_iterator elem = slice_desc.begin();
			for (size_t i=0; i<slice_desc.size(); ++i, ++elem)
				{ full_dst[*elem] += sigma*small_src[i]; }
		}


	 /// substract: slice of vector -> small vector
	template<class VT>
	void subtract_vector_slice(const VT &full_src, slice_desc_type desc, VT &small_dst) const
	{ add_vector_slice(full_src, desc, small_dst, -1.0); }

	/// substract: small vector -> slice of a vector
	template<class VT>
	void subtract_vector_slice(const VT &small_src, VT &full_dst, slice_desc_type desc) const
	{ add_vector_slice(small_src, full_dst, desc, -1.0); }


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

	//! Number of elements for each type
	size_t get_num_elems(slice_desc_type type) const
	{return slice(type).size();}


	//! Create a (partial) clone, without copying values
	template<class VT>
	SmartPtr<VT> slice_clone_without_values(const VT &full_src, slice_desc_type type) const
	{
		const slice_desc_set &slice_desc = slice(type);

		SmartPtr<VT> clone(new VT(slice_desc.size()));
#ifdef UG_PARALLEL
		// set layout
		clone->set_layouts(get_slice_layouts(full_src.layouts(), type));

		// set mask
		uint mask = full_src.get_storage_mask();
		clone->set_storage_type(mask);
		//std::cout << "Storage Type:" << full_src.get_storage_type() <<", mask=" << mask << std::endl;
#endif
		return clone;

	}

	//! Create a (partial) clone
	template<class VT>
	SmartPtr<VT> slice_clone(const VT &full_src, slice_desc_type type) const
	{
		SmartPtr<VT> clone = slice_clone_without_values(full_src, type);
		get_vector_slice(full_src, type, *clone);
		return clone;
	}


protected:

	//! returns type for a global index
	slice_desc_type get_type(size_t index)
	{return m_slice_types[index];}


	//! returns the set of global indices for a given type
	const slice_desc_set &slice(slice_desc_type type) const
	{return m_slice_set[type];}

	slice_desc_set &slice(slice_desc_type type)
	{return m_slice_set[type];}


	//! returns local index for a global index
	bool find_index(slice_desc_type type, int gindex, int &index) const
	{
		// WARNING int index < size_t myset.size() WARNING
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
		UG_ASSERT( (!found || index<(int)myset.size()) , "Invalid index found!");
	//	}
		return found;
	}
#ifdef UG_PARALLEL
public:
	SmartPtr<AlgebraLayouts> get_slice_layouts(ConstSmartPtr<AlgebraLayouts> layouts, slice_desc_type type) const
	{
		// convert layouts (vector->slice)
		SmartPtr<AlgebraLayouts> slice_layouts(new AlgebraLayouts(*layouts));
		replace_indices_in_layout(type, slice_layouts->master());
		replace_indices_in_layout(type, slice_layouts->slave());


		//UG_LOG(*slice_layouts);
		return slice_layouts;
	}

protected:
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
				UG_COND_THROW(!found, "SlicingData:: Did not find index???");
				interf.get_element(eiter) = newind;

			}
		}
	}

#endif /* UG_PARALLEL */

};

}


#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__SCHUR_SLICING_H_ */
