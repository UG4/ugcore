#ifndef __H__UG_attachment_io_traits
#define __H__UG_attachment_io_traits

namespace ug{


template <class TAttachment>
struct attachment_io_traits{
	typedef typename TAttachment::ValueType									value_type;
	typedef typename attachment_value_traits<value_type>::reference			reference_type;
	typedef typename attachment_value_traits<value_type>::const_reference 	const_reference_type;

	static void write_value (std::ostream& out, const_reference_type v)	{out << v;}
	static void read_value(std::istream& in, reference_type v)			{in >> v;}
};

template <>
struct attachment_io_traits<Attachment<bool> > {
	typedef bool													value_type;
	typedef attachment_value_traits<value_type>::reference			reference_type;
	typedef attachment_value_traits<value_type>::const_reference 	const_reference_type;

	static void write_value (std::ostream& out, const_reference_type v)	{out << v;}
	static void read_value(std::istream& in, reference_type v)			{value_type tmp; in >> tmp; v = tmp;}
};

template <>
struct attachment_io_traits<Attachment<vector1> > {
	typedef vector1													value_type;
	typedef attachment_value_traits<value_type>::reference			reference_type;
	typedef attachment_value_traits<value_type>::const_reference 	const_reference_type;

	static void write_value (std::ostream& out, const_reference_type v)
	{
		out << v[0];
	}
	static void read_value(std::istream& in, reference_type v)	
	{
		in >> v[0];
	}
};


template <>
struct attachment_io_traits<Attachment<vector2> > {
	typedef vector2													value_type;
	typedef attachment_value_traits<value_type>::reference			reference_type;
	typedef attachment_value_traits<value_type>::const_reference 	const_reference_type;

	static void write_value (std::ostream& out, const_reference_type v)
	{
		out << v[0] << " " << v[1];
	}
	static void read_value(std::istream& in, reference_type v)	
	{
		in >> v[0] >> v[1];
	}
};


template <>
struct attachment_io_traits<Attachment<vector3> > {
	typedef vector3													value_type;
	typedef attachment_value_traits<value_type>::reference			reference_type;
	typedef attachment_value_traits<value_type>::const_reference 	const_reference_type;

	static void write_value (std::ostream& out, const_reference_type v)
	{
		out << v[0] << " " << v[1] << " " << v[2];
	}
	static void read_value(std::istream& in, reference_type v)	
	{
		in >> v[0] >> v[1] >> v[2];
	}
};


template <>
struct attachment_io_traits<Attachment<vector4> > {
	typedef vector4													value_type;
	typedef attachment_value_traits<value_type>::reference			reference_type;
	typedef attachment_value_traits<value_type>::const_reference 	const_reference_type;

	static void write_value (std::ostream& out, const_reference_type v)
	{
		out << v[0] << " " << v[1] << " " << v[2] << " " << v[3];
	}
	static void read_value(std::istream& in, reference_type v)	
	{
		in >> v[0] >> v[1] >> v[2] >> v[3];
	}
};


////////////////////////////////////////////////////////////////////////////////
/// serialization for std::vector<T> with type T - e. g. std::vector<bool>
////////////////////////////////////////////////////////////////////////////////
template <typename T>
struct attachment_io_traits<Attachment<std::vector<T> > >
/*!
 * \brief serializes/deserializes an Attachment of std::vector<T>
 * \note  this is possible whenever the corresponding Attachment<T>
 *        can be serialized/deserialized (cf. above)
 */
{
	typedef typename std::vector<T>											value_type;
	typedef typename attachment_value_traits<value_type>::reference			reference_type;
	typedef typename attachment_value_traits<value_type>::const_reference 	const_reference_type;

	static void write_value (std::ostream& out, const_reference_type v)
	{
		out << v.size() << " ";
		for (size_t i = 0; i < v.size(); ++i)
			attachment_io_traits<Attachment<T> >::write_value(out, v[i]);
	}
	static void read_value(std::istream& in, reference_type v)
	{
		size_t sz;
		in >> sz;
		v.resize(sz);
		for (size_t i = 0; i < sz; ++i)
			attachment_io_traits<Attachment<T> >::read_value(in, v[i]);
	}
};


}//	end of namespace

#endif	//__H__UG_attachment_io_traits
