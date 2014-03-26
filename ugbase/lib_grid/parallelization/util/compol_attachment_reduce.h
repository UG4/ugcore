// created by Sebastian Reiter
// s.b.reiter@gmail.com
// march 2014

#ifndef __H__UG_PCL_COMPOL_ATTACHMENT_REDUCE__
#define __H__UG_PCL_COMPOL_ATTACHMENT_REDUCE__

#include <algorithm>
#include "pcl/pcl_methods.h"
#include "common/serialization.h"

namespace ug{

///	methods defined in those traits are used by ComPol_AttachmentReduce
/**	A default implementation is provided which works for integer types*/
template <class TValue>
struct attachment_reduce_traits
{
	typedef TValue value_t;
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
	typedef float value_t;
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
	typedef double value_t;
	static inline value_t min(value_t v1, value_t v2)	{return std::min(v1, v2);}
	static inline value_t max(value_t v1, value_t v2)	{return std::max(v1, v2);}
	static inline value_t sum(value_t v1, value_t v2)	{return v1 + v2;}
	static inline value_t prod(value_t v1, value_t v2)	{return v1 * v2;}
	static inline value_t land(value_t v1, value_t v2)	{return v1 && v2;}
	static inline value_t band(value_t v1, value_t v2)	{UG_THROW("doubles do not support a binary and operation.");}
	static inline value_t lor(value_t v1, value_t v2)	{return v1 || v2;}
	static inline value_t bor(value_t v1, value_t v2)	{UG_THROW("doubles do not support a binary or operation.");}
};


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
template <class TLayout, class TAttachment>
class ComPol_AttachmentReduce : public pcl::ICommunicationPolicy<TLayout>
{
	public:
		typedef TLayout							Layout;
		typedef typename Layout::Type			GeomObj;
		typedef typename Layout::Element		Element;
		typedef typename Layout::Interface		Interface;
		typedef typename TAttachment::ValueType Value;
		typedef typename Interface::const_iterator	iiter_t;
		typedef attachment_reduce_traits<Value>	art;

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


template <class TLayout, class TAttachment>
bool ComPol_AttachmentReduce<TLayout, TAttachment>::
collect(BinaryBuffer& buff, const Interface& interface)
{
	for(iiter_t iter = interface.begin(); iter != interface.end(); ++iter)
		Serialize(buff, m_aa[interface.get_element(iter)]);
	return true;
}

template <class TLayout, class TAttachment>
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

template <class TLayout, class TAttachment>
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

template <class TLayout, class TAttachment>
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

template <class TLayout, class TAttachment>
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

template <class TLayout, class TAttachment>
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

template <class TLayout, class TAttachment>
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

template <class TLayout, class TAttachment>
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

template <class TLayout, class TAttachment>
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

template <class TLayout, class TAttachment>
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
