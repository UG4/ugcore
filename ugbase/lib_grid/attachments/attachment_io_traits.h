/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG_attachment_io_traits
#define __H__UG_attachment_io_traits

namespace ug{


template <typename TAttachment>
struct attachment_io_traits{
	using value_type = typename TAttachment::ValueType;
	using reference_type = typename attachment_value_traits<value_type>::reference;
	using const_reference_type = typename attachment_value_traits<value_type>::const_reference;

	static void write_value (std::ostream& out, const_reference_type v)	{out << v;}
	static void read_value(std::istream& in, reference_type v) {in >> v;}
};

template <>
struct attachment_io_traits<Attachment<bool> > {
	using value_type = bool;
	using reference_type = attachment_value_traits<value_type>::reference;
	using const_reference_type = attachment_value_traits<value_type>::const_reference;

	static void write_value (std::ostream& out, const_reference_type v)	{out << v;}
	static void read_value(std::istream& in, reference_type v) {value_type tmp; in >> tmp; v = tmp;}
};

template <>
struct attachment_io_traits<Attachment<vector1> > {
	using value_type = vector1;
	using reference_type = attachment_value_traits<value_type>::reference;
	using const_reference_type = attachment_value_traits<value_type>::const_reference;

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
	using value_type = vector2;
	using reference_type = attachment_value_traits<value_type>::reference;
	using const_reference_type = attachment_value_traits<value_type>::const_reference;

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
	using value_type = vector3;
	using reference_type = attachment_value_traits<value_type>::reference;
	using const_reference_type = attachment_value_traits<value_type>::const_reference;

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
	using value_type = vector4;
	using reference_type = attachment_value_traits<value_type>::reference;
	using const_reference_type = attachment_value_traits<value_type>::const_reference;

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
	using value_type = std::vector<T>;
	using reference_type = typename attachment_value_traits<value_type>::reference;
	using const_reference_type = typename attachment_value_traits<value_type>::const_reference;

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

#endif