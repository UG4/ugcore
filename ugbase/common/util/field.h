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

#ifndef __H__UG_field__
#define __H__UG_field__

namespace ug{

template <class T>
class Field{
	public:
		Field();
		Field(size_t width, size_t height);
		Field(size_t width, size_t height, const T& value);
		Field(const Field& f);
		~Field();

		Field& operator=(const Field& field);

		void		resize_no_copy(size_t width, size_t height);

		inline T&		at(size_t x, size_t y);
		inline const T&	at(size_t x, size_t y) const;

		size_t 		width() const		{return m_width;}
		size_t 		height() const		{return m_height;}
		size_t 		size() const		{return m_width * m_height;}
		size_t		capacity() const	{return m_capacity;}
		T*			data()				{return m_data;}
		const T*	data() const		{return m_data;}

		void fill(size_t x, size_t y, size_t w, size_t h, const T& value);
		void fill_all(const T& value);
		void copy(size_t x, size_t y, const Field& f);
		void swap(Field& f);

	private:
		inline size_t array_index(size_t x, size_t y) const;

		size_t	m_width;
		size_t	m_height;
		size_t	m_capacity;
		T*		m_data;
};

}//	end of namespace


////////////////////////////////////////
//	include implementation
#include "field_impl.hpp"

#endif	//__H__UG_field__
