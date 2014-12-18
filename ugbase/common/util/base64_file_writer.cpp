#include "common/util/base64_file_writer.h"

// for base64 encoding with boost
#include <boost/archive/iterators/transform_width.hpp>
#include <boost/archive/iterators/base64_from_binary.hpp>
#include <boost/archive/iterators/ostream_iterator.hpp>
//#include <boost/filesystem.hpp>

// debug includes!!
#include "common/profiler/profiler.h"
#include "common/error.h"
#include "common/assert.h"

using namespace std;

/**
 * \brief This is the actual encoder using Boost iterators
 * note that final padding to triplet boundary has to be performed manually!
 * see: http://www.boost.org/doc/libs/1_48_0/libs/serialization/doc/dataflow.html
 */
typedef boost::archive::iterators::base64_from_binary<
		// convert binary values to base64 characters
			boost::archive::iterators::transform_width<
		// retrieve 6 bit integers from a sequence of 8 bit bytes
			const char *, 6, 8>
		// compose all the above operations in to a new iterator
		> base64_text;

namespace ug {

////////////////////////////////////////////////////////////////////////////////
// PUBLIC FUNCTIONS

Base64FileWriter::fmtflag Base64FileWriter::format() const {
	return m_currFormat;
}

Base64FileWriter& Base64FileWriter::operator<<(const fmtflag format)
{
	PROFILE_FUNC();

	// forceful flushing of encoder's internal input buffer is necessary
	// if we are switching formats.
	if (format != m_currFormat && m_numBytesWritten > 0) {
		flushInputBuffer(true);
	}
	m_currFormat = format;
	return *this;
}

Base64FileWriter& Base64FileWriter::operator<<(int value)
{
	dispatch(value);
	return *this;
}

Base64FileWriter& Base64FileWriter::operator<<(char value) {
	dispatch(value);
	return *this;
}

Base64FileWriter& Base64FileWriter::operator<<(const char* value)
{
	dispatch(value);
	return *this;
}

Base64FileWriter& Base64FileWriter::operator<<(const string& value)
{
	dispatch(value);
	return *this;
}

Base64FileWriter& Base64FileWriter::operator<<(float value)
{
	dispatch(value);
	return *this;
}

Base64FileWriter& Base64FileWriter::operator<<(double value)
{
	dispatch(value);
	return *this;
}

Base64FileWriter& Base64FileWriter::operator<<(long value)
{
	dispatch(value);
	return *this;
}

Base64FileWriter& Base64FileWriter::operator<<(size_t value)
{
	dispatch(value);
	return *this;
}

Base64FileWriter::Base64FileWriter() :
	m_currFormat(base64_ascii),
	m_inBuffer(ios_base::binary | ios_base::out | ios_base::in),
	m_lastInputByteSize(0),
	m_numBytesWritten(0)
{}

Base64FileWriter::Base64FileWriter(const char* filename,
		const ios_base::openmode mode) :
	m_currFormat(base64_ascii),
	m_inBuffer(ios_base::binary | ios_base::out | ios_base::in),
	m_lastInputByteSize(0),
	m_numBytesWritten(0)
{
	PROFILE_FUNC();

	open(filename, mode);
}

Base64FileWriter::~Base64FileWriter()
{
	flushInputBuffer(true);
	m_fStream.close();
}

void Base64FileWriter::open(const char *filename,
								const ios_base::openmode mode)
{
	// TODO: Create non-existing subdirectories in a platform-independent way
	// like in the following code. (Problem: no header-only boost implementation!)
	/*
	if (mode & std::ios_base::out)
	{
		boost::filesystem::path p(filename);
		boost::filesystem::create_directories(p.parent_path());
	}
	*/

	m_fStream.open(filename, mode);
	if (!m_fStream.is_open()) {
		UG_THROW( "Could not open output file: " << filename);
	} else if (!m_fStream.good()) {
		UG_THROW( "Can not write to output file: " << filename);
	}
}

////////////////////////////////////////////////////////////////////////////////
// PRIVATE FUNCTIONS

// this function performs conversion to plain const char* and stores it in
// m_inBuffer, either raw for binary mode or as string conversion in base64_ascii.
template <typename T>
void Base64FileWriter::dispatch(const T& value)
{
//	PROFILE_FUNC(); // this profile node is too small
	assertFileOpen();

	switch ( m_currFormat ) {
		case base64_ascii:
			// create string representation of value and store it in buffer
			m_inBuffer << value;
			// written bytes is equivalent to current write pointer pos
			m_numBytesWritten = m_inBuffer.tellp();
			m_lastInputByteSize = sizeof(char);
			// check if a buffer flush is needed
			flushInputBuffer();
			break;
		case base64_binary: {
			// write the value in binary mode to the input buffer
			UG_ASSERT(m_inBuffer.good(), "can not write to buffer")
			m_inBuffer.write(reinterpret_cast<const char*>(&value), sizeof(T));
			UG_ASSERT(m_inBuffer.good(), "write failed")

			m_numBytesWritten += sizeof(T);
			m_lastInputByteSize = sizeof(T);
			// check if a buffer flush is needed
			flushInputBuffer();
			break;
		}
		case normal:
			// nothing to do here, almost
			m_fStream << value;
			break;
	}
}

inline void Base64FileWriter::assertFileOpen()
{
	if (m_fStream.bad() || !m_fStream.is_open()) {
		UG_THROW( "File stream is not open." );
	}
}

void Base64FileWriter::flushInputBuffer(bool force)
{
//	PROFILE_FUNC(); // this profile node is too small

	size_t buff_len = 0;
	// amount of elements to flush at once
	const uint elements_to_flush = 12;
	// in case of normal format, no input size is known, so this evals to zero.
	const uint bytes_to_flush = 3 * m_lastInputByteSize * elements_to_flush;
	// in case of forced flush we have to add padding
	uint paddChars = 0;

	// if force, flush all bytes in input stream
	if (force) {
		paddChars = (3 - m_numBytesWritten % 3) % 3;
		buff_len = m_numBytesWritten;
	} else if (bytes_to_flush <= m_numBytesWritten) {
		buff_len = bytes_to_flush;
	} else {
		// buff_len == 0
		return;
	}

	m_tmpBuff.clear();
	// enlarge the input stream by # paddChars because boost reads beyond buffer
	m_tmpBuff.resize(buff_len + paddChars, 0x0);

	if(!m_tmpBuff.empty()){
		char* buff = &m_tmpBuff[0];

		if (!m_inBuffer.good()) {
			UG_THROW("input buffer not good before read");
		}

		// read the buffer
		m_inBuffer.read(buff, buff_len);

		if (!m_inBuffer.good()) {
			UG_THROW("failed to read from input buffer");
		}

		// encode buff in base64
		copy(base64_text(buff), base64_text(buff + buff_len),
				boost::archive::iterators::ostream_iterator<char>(m_fStream));
	}

	size_t rest_len = m_numBytesWritten - buff_len;

	if (rest_len > 0) {
		m_tmpBuff.clear();
		m_tmpBuff.resize(rest_len);
		char* rest = &m_tmpBuff[0];
		// read the rest
		m_inBuffer.read(rest, rest_len);
		if (!m_inBuffer.good()) {
			UG_THROW("failed to read from input buffer");
		}

		// reset the internal string
		m_inBuffer.str("");
		// set write pos to beginning
		m_inBuffer.seekp(0, ios_base::beg);

		// and write the rest
		m_inBuffer.write(rest, rest_len);
		if (!m_inBuffer.good()) {
			UG_THROW("failed to write from input buffer");
		}
	} else {
		// reset buffer
		m_inBuffer.str("");
		m_inBuffer.seekp(0, ios_base::beg);
	}

	// set read and write position to beginning
	m_inBuffer.seekg(0, ios_base::beg);

	if (force) {
		for(uint i = 0; i < paddChars; ++i)
			m_fStream << '=';

		// resetting num bytes written and bytes in block
		m_numBytesWritten = 0;
	} else {
		// update amount of bytes in input buffer
		m_numBytesWritten -= buff_len;
	}
}

void Base64FileWriter::close()
{
	PROFILE_FUNC();

	// make sure all remaining content of the input buffer is encoded and flushed
	flushInputBuffer(true);

	// only when this is done, close the file stream
	m_fStream.close();
	UG_ASSERT(m_fStream.good(), "could not close output file.");
}

}	// namespace ug
