/**
 * \file matrix_io.cpp
 * \author Torbjoern Klatt
 * \date 2012-05-06
 */

#include "matrix_io.h"

namespace ug
{

// /////////////////////////////////////////////////////////////////////////////
// Constructors / Destructor
MatrixIO::MatrixIO() :
  m_pMatFileName( NULL ), m_matFileStream(), m_matFileType( 0 )
//  ,m_rows(0), m_cols(0), m_lines(0)
{
}

MatrixIO::MatrixIO( std::string mFile, int openMode ) :
  m_matFileStream(), m_matFileType( 0 )//, m_rows(0), m_cols(0), m_lines(0)
{
  set_mat_file_name( mFile, openMode );
}

MatrixIO::~MatrixIO()
{
  close_file();
  delete m_pMatFileName;
}

// /////////////////////////////////////////////////////////////////////////////
// Public Member Functions
void MatrixIO::set_mat_file_name( std::string mFile, int openMode )
{
  if ( !mFile.empty() ) {
    if ( openMode == EXISTING ) {
      UG_ASSERT( FileExists( mFile.c_str() ),
                "File " << mFile.c_str() << " could not be found." );
    } else if( openMode == NEW ) {
      std::ofstream createFile;
      createFile.open( mFile.c_str(), std::ios_base::out );
      UG_ASSERT( createFile.is_open(), "File could not be created." );
      createFile.close();
    } else {
      UG_THROW( "Invalid open mode specified: " << openMode );
    }
    m_pMatFileName = new std::string( mFile );
  }
}

std::string MatrixIO::get_mat_file_name() const
{
  return std::string( *m_pMatFileName );
}

// /////////////////////////////////////////////////////////////////////////////
// Private Member Functions
void MatrixIO::open_file( std::ios_base::openmode mode )
{
  UG_ASSERT( !m_pMatFileName->empty(), "Matrix File not set." );

  m_matFileStream.open( m_pMatFileName->c_str(), mode );
  UG_ASSERT( m_matFileStream.fail() || m_matFileStream.bad(),
             "Matrix File could not be opend for reading/writing.\n"
             << "iostate: " << m_matFileStream.rdstate() );
}

void MatrixIO::close_file()
{
  if ( m_matFileStream.is_open() ) {
//     m_pMatFileStream->flush();
    m_matFileStream.close();
  }
}



} // namespace ug

// EOF
