/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#ifndef __H__UG_PCL_COMPOL_ATTACHMENT_REDUCE__
#define __H__UG_PCL_COMPOL_ATTACHMENT_REDUCE__

#include <algorithm>
#include "pcl/pcl_methods.h"
#include "common/serialization.h"

namespace ug{

///	methods defined in those traits are used by ComPol_AttachmentReduce
/**	A default implementation is provided which works for integer types*/
template <typename TValue>
struct attachment_reduce_traits
{
	using value_t = TValue;
	static inline value_t min(value_t v1, value_t v2)	{return std::min(v1, v2);}
	static inline value_t max(value_t v1, value_t v2)	{return std::max(v1, v2);}
	static inline value_t sum(value_t v1, value_t v2)	{return v1 + v2;}
	static inline value_t prod(value_t v1, value_t v2)	{return v1 * v2;}
	static inline value_t land(value_t v1, value_t v2)	{return v1 && v2;}
	static inline value_t band(value_t v1, value_t v2)	{return v1 & v2;}
	static inline value_t lor(value_t v1, value_t v2)	{return v1 || v2;}
	static inline value_t bor(value_t v1, value_t v2)	{return v1 | v2;}
};

/**	Specialization for float. No band and bor operations are supported.*/
template <>
struct attachment_reduce_traits<float>
{
	using value_t = float;
	static inline value_t min(value_t v1, value_t v2)	{return std::min(v1, v2);}
	static inline value_t max(value_t v1, value_t v2)	{return std::max(v1, v2);}
	static inline value_t sum(value_t v1, value_t v2)	{return v1 + v2;}
	static inline value_t prod(value_t v1, value_t v2)	{return v1 * v2;}
	static inline value_t land(value_t v1, value_t v2)	{return v1 && v2;}
	static inline value_t band(value_t v1, value_t v2)	{UG_THROW("floats do not support a binary and operation.");}
	static inline value_t lor(value_t v1, value_t v2)	{return v1 || v2;}
	static inline value_t bor(value_t v1, value_t v2)	{UG_THROW("floats do not support a binary or operation.");}
};

/**	Specialization for double. No band and bor operations are supported.*/
template <>
struct attachment_reduce_traits<double>
{
	using value_t = double;
	static inline value_t min(value_t v1, value_t v2)	{return std::min(v1, v2);}
	static inline value_t max(value_t v1, value_t v2)	{return std::max(v1, v2);}
	static inline value_t sum(value_t v1, value_t v2)	{return v1 + v2;}
	static inline value_t prod(value_t v1, value_t v2)	{return v1 * v2;}
	static inline value_t land(value_t v1, value_t v2)	{return v1 && v2;}
	static inline value_t band(value_t v1, value_t v2)	{UG_THROW("doubles do not support a binary and operation.");}
	static inline value_t lor(value_t v1, value_t v2)	{return v1 || v2;}
	static inline value_t bor(value_t v1, value_t v2)	{UG_THROW("doubles do not support a binary or operation.");}
};

template <int dim>
struct vector_attachment_reduce_traits
{
	using value_t = MathVector<dim>;
	static inline value_t min(value_t v1, value_t v2)
	{
		value_t v;
		for(int i = 0; i < dim; ++i)
			v[i] = std::min(v1[i], v2[i]);
		return v;
	}
	static inline value_t max(value_t v1, value_t v2)
	{
		value_t v;
		for(int i = 0; i < dim; ++i)
			v[i] = std::max(v1[i], v2[i]);
		return v;
	}
	static inline value_t sum(value_t v1, value_t v2)
	{
		value_t v = v1;
		v += v2;
		return v;
	}
	static inline value_t prod(value_t v1, value_t v2)
	{
		value_t v;
		for(int i = 0; i < dim; ++i)
			v[i] = v1[i] * v2[i];
		return v;
	}
	static inline value_t land(value_t v1, value_t v2)
	{
		value_t v;
		for(int i = 0; i < dim; ++i)
			v[i] = v1[i] && v2[i];
		return v;
	}
	static inline value_t band(value_t v1, value_t v2)	{UG_THROW("vectors do not support a binary and operation.");}
	static inline value_t lor(value_t v1, value_t v2)
	{
		value_t v;
		for(int i = 0; i < dim; ++i)
			v[i] = v1[i] || v2[i];
		return v;
	}
	static inline value_t bor(value_t v1, value_t v2)	{UG_THROW("vectors do not support a binary or operation.");}
};

template <>
struct attachment_reduce_traits<MathVector<1> > :
		vector_attachment_reduce_traits<1>	{};

template <>
struct attachment_reduce_traits<MathVector<2> > :
		vector_attachment_reduce_traits<2>	{};

template <>
struct attachment_reduce_traits<MathVector<3> > :
		vector_attachment_reduce_traits<3>	{};

template <>
struct attachment_reduce_traits<MathVector<4> > :
		vector_attachment_reduce_traits<4>	{};


// implementation for a std::vector<number> of arbitrary size
struct std_number_vector_attachment_reduce_traits
{
	using value_t = std::vector<number>;

	// we have to make sure the two vectors are of the same length
	// if they are not, the smaller one will be resized to fit the larger one's size
	// fill-in with zeros
	static inline void check_length(value_t& v1, value_t& v2)
	{
		if (v1.size() == v2.size()) return;

		if (v1.size() > v2.size())
		{
			v2.resize(v1.size(), 0.0);
			return;
		}

		v1.resize(v2.size(), 0.0);
	}

	static inline value_t min(value_t v1, value_t v2)
	{
		check_length(v1,v2);
		size_t sz = v1.size();
		value_t v(sz);
		for (size_t i = 0; i < sz; i++)
			v[i] = std::min(v1[i], v2[i]);
		return v;
	}

	static inline value_t max(value_t v1, value_t v2)
	{
		check_length(v1,v2);
		size_t sz = v1.size();
		value_t v(sz);
		for (size_t i = 0; i < sz; i++)
			v[i] = std::max(v1[i], v2[i]);
		return v;
	}

	static inline value_t sum(value_t v1, value_t v2)
	{
		check_length(v1,v2);
		size_t sz = v1.size();
		value_t v(sz);
		for (size_t i = 0; i < sz; i++)
			v[i] = v1[i] + v2[i];
		return v;
	}

	static inline value_t prod(value_t v1, value_t v2)
	{
		check_length(v1,v2);
		size_t sz = v1.size();
		value_t v(sz);
		for (size_t i = 0; i < sz; i++)
			v[i] = v1[i] * v2[i];
		return v;
	}

	static inline value_t land(value_t v1, value_t v2)
	{
		check_length(v1,v2);
		size_t sz = v1.size();
		value_t v(sz);
		for (size_t i = 0; i < sz; i++)
			v[i] = v1[i] && v2[i];
		return v;
	}

	static inline value_t band(value_t v1, value_t v2)
	{
		UG_THROW("Vectors of number do not support a binary and operation.");
	}

	static inline value_t lor(value_t v1, value_t v2)
	{
		check_length(v1,v2);
		size_t sz = v1.size();
		value_t v(sz);
		for (size_t i = 0; i < sz; i++)
			v[i] = v1[i] || v2[i];
		return v;
	}

	static inline value_t bor(value_t v1, value_t v2)
	{
		UG_THROW("Vectors of number do not support a binary or operation.");
	}
};

template <>
struct attachment_reduce_traits<std::vector<number> > :
	public std_number_vector_attachment_reduce_traits {};


///	Performs reduce operations on the specified attachment
/**	Currently the following reduce operations are supported:
 *		- PCL_RO_MAX
 *		- PCL_RO_MIN
 *		- PCL_RO_SUM
 *		- PCL_RO_PROD
 *		- PCL_RO_LAND
 *		- PCL_RO_BAND
 *		- PCL_RO_LOR
 *		- PCL_RO_BOR
 *
 * The reduce operation is performed on the attachment of each grid-object
 * separately during extraction.
 *
 * Reduce operations are locally performed using the operations defined in
 * attachment_reduce_traits for the given type.
 *
 * If you want to perform an 'allreduce' on the
 * attached values you have to perform a reduce first, followed by a copy through
 * ComPol_CopyAttachment
 */
template <typename TLayout, typename TAttachment>
class ComPol_AttachmentReduce : public pcl::ICommunicationPolicy<TLayout>
{
	public:
		using Layout = TLayout;
		using GeomObj = typename Layout::Type;
		using Element = typename Layout::Element;
		using Interface = typename Layout::Interface;
		using Value = typename TAttachment::ValueType;
		using iiter_t = typename Interface::const_iterator;
		using art = attachment_reduce_traits<Value>;

		enum ReduceOperation{
			NONE,
			MAX,
			MIN,
			SUM,
			PROD,
			LAND,
			BAND,
			LOR,
			BOR
		};

		ComPol_AttachmentReduce(Grid& grid, TAttachment& attachment,
								pcl::ReduceOperation op) :
			m_grid(grid),
			m_aa(grid, attachment),
			m_op(NONE)
		{
			if(op == PCL_RO_MAX)		m_op = MAX;
			else if(op == PCL_RO_MIN)	m_op = MIN;
			else if(op == PCL_RO_SUM)	m_op = SUM;
			else if(op == PCL_RO_PROD)	m_op = PROD;
			else if(op == PCL_RO_LAND)	m_op = LAND;
			else if(op == PCL_RO_BAND)	m_op = BAND;
			else if(op == PCL_RO_LOR)	m_op = LOR;
			else if(op == PCL_RO_BOR)	m_op = BOR;
		}

		ComPol_AttachmentReduce(Grid& grid, TAttachment& attachment,
								ReduceOperation op) :
			m_grid(grid),
			m_aa(grid, attachment),
			m_op(op)
		{}

	///	writes the data for the given interface to the buffer.
		bool collect(BinaryBuffer& buff, const Interface& interface);

	///	reads the data from the buffer to the given interface .
		bool extract(BinaryBuffer& buff, const Interface& interface);

	protected:
		bool extract_max(BinaryBuffer& buff, const Interface& interface);
		bool extract_min(BinaryBuffer& buff, const Interface& interface);
		bool extract_sum(BinaryBuffer& buff, const Interface& interface);
		bool extract_prod(BinaryBuffer& buff, const Interface& interface);
		bool extract_land(BinaryBuffer& buff, const Interface& interface);
		bool extract_band(BinaryBuffer& buff, const Interface& interface);
		bool extract_lor(BinaryBuffer& buff, const Interface& interface);
		bool extract_bor(BinaryBuffer& buff, const Interface& interface);

	private:
		Grid& m_grid;
		Grid::AttachmentAccessor<GeomObj, TAttachment>	m_aa;
		ReduceOperation m_op;
};


template <typename TLayout, typename TAttachment>
bool ComPol_AttachmentReduce<TLayout, TAttachment>::
collect(BinaryBuffer& buff, const Interface& interface)
{
	for(iiter_t iter = interface.begin(); iter != interface.end(); ++iter)
		Serialize(buff, m_aa[interface.get_element(iter)]);
	return true;
}

template <typename TLayout, typename TAttachment>
bool ComPol_AttachmentReduce<TLayout, TAttachment>::
extract(BinaryBuffer& buff, const Interface& interface)
{
	switch(m_op){
		case MAX:	return extract_max(buff, interface);
		case MIN:	return extract_min(buff, interface);
		case SUM:	return extract_sum(buff, interface);
		case PROD:	return extract_prod(buff, interface);
		case LAND:	return extract_land(buff, interface);
		case BAND:	return extract_band(buff, interface);
		case LOR:	return extract_lor(buff, interface);
		case BOR:	return extract_bor(buff, interface);
		default:
			UG_THROW("Unsupported reduce operation in ComPol_AttachmentReduce::extract:"
					 << m_op);
			break;
	}
	return false;
}

template <typename TLayout, typename TAttachment>
bool ComPol_AttachmentReduce<TLayout, TAttachment>::
extract_max(BinaryBuffer& buff, const Interface& interface)
{
	using std::max;
	Value v;
	for(iiter_t iter = interface.begin(); iter != interface.end(); ++iter){
		Deserialize(buff, v);
		Element e = interface.get_element(iter);
		m_aa[e] = art::max(v, m_aa[e]);
	}
	return true;
}

template <typename TLayout, typename TAttachment>
bool ComPol_AttachmentReduce<TLayout, TAttachment>::
extract_min(BinaryBuffer& buff, const Interface& interface)
{
	using std::min;
	Value v;
	for(iiter_t iter = interface.begin(); iter != interface.end(); ++iter){
		Deserialize(buff, v);
		Element e = interface.get_element(iter);
		m_aa[e] = art::min(v, m_aa[e]);
	}
	return true;
}

template <typename TLayout, typename TAttachment>
bool ComPol_AttachmentReduce<TLayout, TAttachment>::
extract_sum(BinaryBuffer& buff, const Interface& interface)
{
	Value v;
	for(iiter_t iter = interface.begin(); iter != interface.end(); ++iter){
		Deserialize(buff, v);
		Element e = interface.get_element(iter);
		m_aa[e] = art::sum(v, m_aa[e]);
	}
	return true;
}

template <typename TLayout, typename TAttachment>
bool ComPol_AttachmentReduce<TLayout, TAttachment>::
extract_prod(BinaryBuffer& buff, const Interface& interface)
{
	Value v;
	for(iiter_t iter = interface.begin(); iter != interface.end(); ++iter){
		Deserialize(buff, v);
		Element e = interface.get_element(iter);
		m_aa[e] = art::prod(v, m_aa[e]);
	}
	return true;
}

template <typename TLayout, typename TAttachment>
bool ComPol_AttachmentReduce<TLayout, TAttachment>::
extract_land(BinaryBuffer& buff, const Interface& interface)
{
	Value v;
	for(iiter_t iter = interface.begin(); iter != interface.end(); ++iter){
		Deserialize(buff, v);
		Element e = interface.get_element(iter);
		m_aa[e] = art::land(v, m_aa[e]);
	}
	return true;
}

template <typename TLayout, typename TAttachment>
bool ComPol_AttachmentReduce<TLayout, TAttachment>::
extract_band(BinaryBuffer& buff, const Interface& interface)
{
	Value v;
	for(iiter_t iter = interface.begin(); iter != interface.end(); ++iter){
		Deserialize(buff, v);
		Element e = interface.get_element(iter);
		m_aa[e] = art::band(v, m_aa[e]);
	}
	return true;
}

template <typename TLayout, typename TAttachment>
bool ComPol_AttachmentReduce<TLayout, TAttachment>::
extract_lor(BinaryBuffer& buff, const Interface& interface)
{
	Value v;
	for(iiter_t iter = interface.begin(); iter != interface.end(); ++iter){
		Deserialize(buff, v);
		Element e = interface.get_element(iter);
		m_aa[e] = art::lor(v, m_aa[e]);
	}
	return true;
}

template <typename TLayout, typename TAttachment>
bool ComPol_AttachmentReduce<TLayout, TAttachment>::
extract_bor(BinaryBuffer& buff, const Interface& interface)
{
	Value v;
	for(iiter_t iter = interface.begin(); iter != interface.end(); ++iter){
		Deserialize(buff, v);
		Element e = interface.get_element(iter);
		m_aa[e] = art::bor(v, m_aa[e]);
	}
	return true;
}

}//	end of namespace

#endif
