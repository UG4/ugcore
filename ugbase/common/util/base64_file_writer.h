/**
 * \file base64_file_writer.h
 * \author Torbjoern Klatt
 * \author Martin Scherer
 * \date 2013-02-25
 */

#ifndef __H__UG__COMMON__UTIL__BASE64_FILE_WRITER__
#define __H__UG__COMMON__UTIL__BASE64_FILE_WRITER__

#include <sstream>
#include <fstream>
#include <vector>

namespace ug {

/// \addtogroup ugbase_common_io
/// \{

/**
 * \brief File writer allowing selective base64 encoding of arbitrary data
 * \details This provides a convenient and well tested way of writing selectivly
 *   base64 encoded data to an arbitrary file.
 *
 *   After initialization with a file name
 *   \code{.cpp}
 *     Base64FileWriter writer( "output_file.txt" );
 *   \endcode
 *   arbitrary data can be written to it solely using the <tt>&lt;&lt;</tt>-operator:
 *   \code{.cpp}
 *     int myInt = 3;
 *     writer << "My int is: " << myInt;
 *   \endcode
 *   The default behaviour is to write base64 binary or ascii encoded data from
 *   the very beginning.
 *   This can be overridden by passing <tt>Base64FileWriter::normal</tt> as data
 *   to it:
 *   \code{.cpp}
 *     writer << Base64FileWriter::normal << "Base64 is not always useful.";
 *   \endcode
 *   Switching back to base64 encoding is simple as
 *   \code{.cpp}
 *     writer << Base64FileWriter::base64_binary << "and back to base64 encoding!";
 *   \endcode
 *   To finish off writing, just close the writer:
 *   \code{.cpp}
 *   writer.close();
 *   \endcode
 */
class Base64FileWriter {
public:

	/**
	 * \brief Format flags to enable deactivation of base64 encoding selectivly
	 */
	enum fmtflag {
		//! enables base64 encoding of all following data until Base64FileWriter::normal is set
		base64_ascii,
		//! encode given values in binary base64
		base64_binary,
		//! behaves as usual std::ofstream
		normal
	};

	//////////////////////////////////////////////////////////////////////////
	Base64FileWriter();

	/**
	 * \brief Constructor with name of file to write to
	 * \param[in] filename name of the output file
	 * \param[in] mode     openmode for the file as defined in `std::ios_base`.
	 *                     Defaults to `out | trunc`, i.e. create new
	 *                     file or truncate existing.
	 * \throws UGError if \c filename can not be opened or is not writeable
	 */
	Base64FileWriter(const char* filename,
				const std::ios_base::openmode mode = (std::ios_base::out |
														std::ios_base::trunc));

	/**
	 * \brief Destructor, which properly flushs encoder's internal buffer and closes file stream
	 */
	~Base64FileWriter();

public:

	void open(const char *filename,
			const std::ios_base::openmode mode = std::ios_base::out );

	/**
	 * \brief gets the current set format
	 */
	fmtflag format() const;

	/**
	 * \brief Closes the file writer properly and encodes any remaining buffer content
	 */
	void close();

	/**
	 * \brief Switch between normal and base64 encoded output
	 * \param format one of the values defined in Base64FileWriter::fmtflag
	 */
	Base64FileWriter& operator<<(const fmtflag format);

	// insert plain standard types to this filewriter
	Base64FileWriter& operator<<(int i);
	Base64FileWriter& operator<<(char c);
	Base64FileWriter& operator<<(const char* cstr);
	Base64FileWriter& operator<<(const std::string& str);

	Base64FileWriter& operator<<(float f);
	Base64FileWriter& operator<<(double d);
	Base64FileWriter& operator<<(long l);
	Base64FileWriter& operator<<(size_t s);

private:
	/**
	 * \brief Writes given data to the output file and encodes it if Base64FileWriter::base64 is set
	 * \param value arbitrary data to write to the output file
	 * \tparam T type of the data to write
	 * \throws UGError if Base64FileWriter::fmtflag is not set to Base64FileWriter::normal
	 *                 or Base64FileWriter::base64
	 */
	template <typename T>
	void dispatch(const T& value);

	/**
	 * \brief File stream to write everything to
	 */
	std::fstream m_fStream;
	/**
	 * \brief Current write format (\c base64 or \c normal)
	 */
	fmtflag m_currFormat;
	/**
	 * \brief Internal input buffer for the encoder
	 */
	std::stringstream m_inBuffer;

	/*
	 * input byte size during a format block
	 */
	size_t m_lastInputByteSize;

	/**
	 * buffer used while flushing
	 */
	std::vector<char> m_tmpBuff;

	/**
	 * number of bytes written between flushes
	 */
	size_t m_numBytesWritten;

	/**
	 * \brief Flushes input buffer
	 * \param force whether to forcefully flush the buffer
	 */
	void flushInputBuffer(bool force = false);

	/**
	 * \brief Check on readiness of the file stream
	 * \throws UGError if either the filestream's \c badbit is set or the file
	 *                 stream is not open.
	 */
	inline void assertFileOpen();
};

// end group ugbase_common_io

} // namespace: ug

#endif // __H__UG__COMMON__UTIL__BASE64_FILE_WRITER__
