/**
 * \file base64_file_writer.h
 * \author Torbjoern Klatt
 * \author Martin Scherer
 * \date 2013-02-25
 */

#ifndef __H__UG__LIB_DISC__IO_BASE64_FILE_WRITER__
#define __H__UG__LIB_DISC__IO_BASE64_FILE_WRITER__

#include <fstream>
#include <cstring>

#include "common/assert.h"
#include "common/error.h"
#include "common/log.h"
#include "common/profiler/profiler.h"

namespace ug {
	
	/**
	 */
	class Base64FileWriter {
		private:
			/**
			 */
			const static size_t INPUT_BUFFER_SIZE = 64;
			/**
			 */
			const static size_t ENCODE_TRIPLET_SIZE = 24;
		
		public:
			/**
			 */
			enum fmtflag {
				base64,
				normal
			};
			
			/**
			 */
			explicit Base64FileWriter( const char* filename ) :
				m_fStream( NULL ),
				m_currFormat( normal ),
				m_inBuffer( new char[INPUT_BUFFER_SIZE] ),
				m_inBufferSize(0)
			{
				PROFILE_FUNC()
				m_fStream.open( filename, std::ios_base::out | std::ios_base::trunc );
				if (!m_fStream.is_open() ) {
					UG_THROW( "Could not open output file: " << filename );
				} else if(! m_fStream.good()) {
					UG_THROW( "Can not write to output file: " << filename );
				}
			}
			
			/**
			 */
			~Base64FileWriter()
			{
				m_fStream.close();
				delete( m_inBuffer );
			}
			
			/**
			 */
			template <typename T>
			Base64FileWriter& operator<<( const T& value ) {
				PROFILE_FUNC()
				UG_ASSERT( ( m_fStream.good() && m_fStream.is_open() ), 
				           "File stream is not open." );
				
				if ( m_currFormat == base64 ) {
					// we should know how big the data is
					size_t valueSize = sizeof( value );
					
					// copy the data to the end of our input buffer
					std::memcpy( m_inBuffer + m_inBufferSize, &value, valueSize );
					// and update the end of it
					m_inBufferSize += valueSize;
					
					// now encode as much as possible from the front of our input buffer
					flushInputBuffer();
					
				} else if ( m_currFormat == normal ) {
					// nothing to do here, almost
					m_fStream << value;
					
				} else {
					UG_THROW( "Output format mode not set." );
				}
				return *this;
			}
			
			/**
			 */
			Base64FileWriter& operator<<( const fmtflag format ) {
				PROFILE_FUNC()
				if ( format != m_currFormat && format == normal ) {
					// base64 encoded stream should end here
					// thus, fill the rest of the buffer with zero-bytes until
					// a triplet is full and encode it
					std::memset( m_inBuffer + m_inBufferSize, 0, ENCODE_TRIPLET_SIZE - m_inBufferSize );
					m_inBufferSize = ENCODE_TRIPLET_SIZE;
					flushInputBuffer();
				}
				return *this;
			}

			void close()
			{
				m_fStream.close();
				UG_ASSERT( m_fStream.good(), "Base64FileWriter: Can not close output file." );
			}
			
			/**
			 */
			static void encodeTriplet( char *in, char *out, int n )
			{
				static char digits[] =
					"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
				unsigned char *_in = (unsigned char *)in;
				
				out[0] = digits[_in[0] >> 2];
				out[1] = digits[((_in[0] & 0x03) << 4) | (_in[1] >> 4)];
				out[2] = n > 1? digits[((_in[1] & 0x0F) << 2) | (_in[2] >> 6)] : '=';
				out[3] = n > 2? digits[_in[2] & 0x3F] : '=';
			}
			
		private:
			/**
			 */
			std::fstream m_fStream;
			/**
			 */
			fmtflag m_currFormat;
			/**
			 */
			char *m_inBuffer;
			/**
			 */
			size_t m_inBufferSize;
			
			/**
			 */
			void flushInputBuffer()
			{
				PROFILE_FUNC()
				// encode as many triplets as there are in our input buffer
				while ( m_inBufferSize >= ENCODE_TRIPLET_SIZE ) {
					char out[ENCODE_TRIPLET_SIZE];
					encodeTriplet( m_inBuffer, out, ENCODE_TRIPLET_SIZE );
					
					// save the encoded triplet to the output stream
					m_fStream << out;
					
					// move the remaining part of the input buffer to the start of it
					std::memmove( m_inBuffer, 
					              m_inBuffer + ENCODE_TRIPLET_SIZE, 
					              m_inBufferSize - ENCODE_TRIPLET_SIZE );
					
					// and update the size of the input buffer
					m_inBufferSize -= ENCODE_TRIPLET_SIZE;
				}
			}
	};
	
} // namespace: ug

#endif // __H__UG__LIB_DISC__IO_BASE64_FILE_WRITER__

// EOF
