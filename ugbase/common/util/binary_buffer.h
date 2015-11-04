#ifndef __H__UG__binary_buffer__
#define __H__UG__binary_buffer__

#include <vector>
#include "common/types.h"

namespace ug
{

/// \addtogroup ugbase_common_io
/// \{

///	A Buffer for binary data.
/**	The BinaryBuffer allows read and write access, which mimics the
 * behavior of std::iostream. However, in contrary to ug::BinaryStream,
 * this class is not a specialization of std::iostream. This is crucial
 * to achieve maximal performance.
 *
 * Furthermore BinaryBuffer gives access to its internal buffer,
 * which can be handy in some situations. However, this feature
 * has to be used with extreme care!
 */
class BinaryBuffer
{
	public:
		BinaryBuffer();

	///	creates a binary buffer and reserves bufSize bytes.
		BinaryBuffer(size_t bufSize);

	///	clears the buffer
	/**	This method does not free associated memory. It only
	 * resets the read and write positions. To free the memory
	 * you may create a new instance of ug::BinaryBuffer and
	 * assign it by value to your current instance.*/
		void clear();

	///	resizes the associated buffer to the given size.
	/**	You can check the number of bytes reserved in the buffer through
	 * ug::BinaryBuffer::capacity*/
		void reserve(size_t newSize);

	///	returns the capacity (reserved memory) of the buffer
		inline size_t capacity() const;

	///	returns the current read-pos (in bytes)
		inline size_t read_pos() const;

	///	returns the current write-pos (in bytes)
		inline size_t write_pos() const;

	///	reads data of the given size (in bytes)
	/**	This automatically advances the read position.*/
		inline void read(char* buf, size_t size);

	/// writes data of the given size (in bytes)
	/**	This automatically advances the write position.*/
		inline void write(const char* buf, size_t size);

	///	returns the raw buffer pointer or NULL if the buffer is empty (capacity() == 0)
		inline char* buffer();

	///	returns true if the read-position reached the write-position
		inline bool eof();

	///	sets the read position (in bytes).
		void set_read_pos(size_t pos);

	///	sets the write position.
		void set_write_pos(size_t pos);

	private:
		std::vector<char>	m_data;
		size_t				m_readPos;
		size_t				m_writePos;
};

// end group ugbase_common_io
/// \}

}//	end of namespace


////////////////////////////////
//	include implementation
#include "binary_buffer_impl.h"

#endif
