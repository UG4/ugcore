// created by Sebastian Reiter
// s.b.reiter@gmail.com

#ifndef __H__UG_attachment_io_traits
#define __H__UG_attachment_io_traits

namespace ug{

template <class TAttachment>
struct attachment_io_traits {
	typedef typename TAttachment::ValueType									value_type;
	typedef typename attachment_value_traits<value_type>::reference			reference_type;
	typedef typename attachment_value_traits<value_type>::const_reference 	const_reference_type;

	static void write_value (std::ostream& out, const_reference_type v)	{out << v;}
	static void read_value(std::istream& in, reference_type v)			{in >> v;}
};


template <>
struct attachment_io_traits<Attachment<bool> > {
	typedef bool															value_type;
	typedef typename attachment_value_traits<value_type>::reference			reference_type;
	typedef typename attachment_value_traits<value_type>::const_reference 	const_reference_type;

	static void write_value (std::ostream& out, const_reference_type v)	{out << v;}
	static void read_value(std::istream& in, reference_type v)			{value_type tmp; in >> tmp; v = tmp;}
};


template <>
struct attachment_io_traits<Attachment<vector1> > {
	typedef vector1															value_type;
	typedef typename attachment_value_traits<value_type>::reference			reference_type;
	typedef typename attachment_value_traits<value_type>::const_reference 	const_reference_type;

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
	typedef vector2															value_type;
	typedef typename attachment_value_traits<value_type>::reference			reference_type;
	typedef typename attachment_value_traits<value_type>::const_reference 	const_reference_type;

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
	typedef vector3															value_type;
	typedef typename attachment_value_traits<value_type>::reference			reference_type;
	typedef typename attachment_value_traits<value_type>::const_reference 	const_reference_type;

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
	typedef vector4															value_type;
	typedef typename attachment_value_traits<value_type>::reference			reference_type;
	typedef typename attachment_value_traits<value_type>::const_reference 	const_reference_type;

	static void write_value (std::ostream& out, const_reference_type v)
	{
		out << v[0] << " " << v[1] << " " << v[2] << " " << v[3];
	}
	static void read_value(std::istream& in, reference_type v)	
	{
		in >> v[0] >> v[1] >> v[2] >> v[3];
	}
};

}//	end of namespace

#endif	//__H__UG_attachment_io_traits
