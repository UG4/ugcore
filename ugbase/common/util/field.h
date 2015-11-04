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

		T& operator=(const Field& field);

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
