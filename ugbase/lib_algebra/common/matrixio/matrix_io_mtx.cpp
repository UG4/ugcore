
#include "matrix_io_mtx.h"

namespace ug
{

// /////////////////////////////////////////////////////////////////////////////
// Constructor / Destructor
MatrixIOMtx::MatrixIOMtx() :
  m_pMatFileName( NULL ), m_matFileStream(), m_firstDataLine(0),
  m_matFileType( MatrixFileType::MATRIX_MARKET ),
  m_mmTypeCode( MMTypeCode() ),
  m_rows( 0 ), m_cols( 0 ), m_lines( 0 )
{}

MatrixIOMtx::MatrixIOMtx( std::string mFile, int openMode ) :
  m_matFileStream(), m_firstDataLine(0),
  m_matFileType( MatrixFileType::MATRIX_MARKET ),
  m_mmTypeCode( MMTypeCode() ), m_rows( 0 ), m_cols( 0 ), m_lines( 0 )
{
  PROFILE_FUNC();
  // set the file name
  set_mat_file_name( mFile, openMode );

  if ( openMode == MatrixIO::EXISTING ) {
    // read matrix type information
    query_matrix_type();

    // read matrix size information
    query_matrix_characteristics();
  }
}

MatrixIOMtx::~MatrixIOMtx()
{
  close_file();
  delete m_pMatFileName;
}

// /////////////////////////////////////////////////////////////////////////////
// Public Functions
void MatrixIOMtx::set_mat_file_name( std::string mFile, int openMode )
{
  PROFILE_FUNC();
  if( !mFile.empty() ) {
    if ( openMode == MatrixIO::EXISTING ) {
      UG_ASSERT( FileExists( mFile.c_str() ),
                "File " << mFile.c_str() << " could not be found." );
    } else if( openMode == MatrixIO::NEW ) {
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

std::string MatrixIOMtx::get_mat_file_name() const
{
  return std::string( *m_pMatFileName );
}

void MatrixIOMtx::set_mat_dims( size_t rows, size_t cols, size_t lines )
{
  PROFILE_FUNC();
  UG_ASSERT( rows > 0, "Number rows must be positive." );
  UG_ASSERT( cols > 0, "Number columns must be positive." );
  UG_ASSERT( lines > 0, "Number data lines must be positive." );

  m_rows = rows;
  m_cols = cols;
  m_lines = lines;
}

size_t MatrixIOMtx::get_num_rows() const
{
  return (size_t)m_rows;
}

size_t MatrixIOMtx::get_num_cols() const
{
  return (size_t)m_cols;
}

size_t MatrixIOMtx::get_num_lines() const
{
  return (size_t)m_lines;
}

bool MatrixIOMtx::is_sparse() const
{
  return m_mmTypeCode.is_sparse();
}


// /////////////////////////////////////////////////////////////////////////////
// Private Functions
void MatrixIOMtx::open_file( std::ios_base::openmode mode )
{
  UG_ASSERT( !m_pMatFileName->empty(), "Matrix File not set." );

  UG_ASSERT( !m_matFileStream.is_open(), "File alread open." );

  m_matFileStream.open( m_pMatFileName->c_str(), mode );
  UG_ASSERT( !m_matFileStream.fail() || !m_matFileStream.bad(),
             "Matrix File could not be opend for reading/writing.\n"
             << "iostate: " << m_matFileStream.rdstate() );
}

void MatrixIOMtx::close_file()
{
  if( m_matFileStream.is_open() ) {
//     m_pMatFileStream->flush();
    m_matFileStream.close();
  }
}

void MatrixIOMtx::query_matrix_type()
{
  PROFILE_FUNC();
  // open the file
  open_file();

  std::stringstream first_line;
  // make sure we are at the beginning of the file
  m_matFileStream.seekg( 0 );

  // get the first line
  std::string s;
  std::getline( m_matFileStream, s );
  first_line.str( s );

  // split the line into its specified parts
  std::string buffer_str = "";
  std::vector<std::string> banner_items;
  while( first_line >> buffer_str ) {
    banner_items.push_back( buffer_str );
  };

  // and make sure the file is a valid MatrixMarket file
  // a) Banner
  UG_ASSERT( banner_items.at( 0 ).compare( MM_BANNER_STR ) == 0,
             "Given file is not a valid Matrix Market file.\n\t"
             << "First line must start with '" << MM_BANNER_STR
             << "' and not with '" << banner_items.at( 0 ) << "'." );
  // b) matrix
  UG_ASSERT( banner_items.at( 1 ).compare( MM_MTX_STR ) == 0,
             "Given file is not a valid Matrix Market file.\n\t"
             << "First line must contain '" << MM_MTX_STR
             << "' as second element and not '" << banner_items.at( 1 ) << "'." );

  // c) coordinate / array
  m_mmTypeCode.set_class_type( banner_items.at( 2 ) );

  // d) data type
  m_mmTypeCode.set_numeric_type( banner_items.at( 3 ) );

  // e) algebraic type
  m_mmTypeCode.set_algebraic_type( banner_items.at( 4 ) );

  // close it
  close_file();
}

void MatrixIOMtx::query_matrix_characteristics()
{
  PROFILE_FUNC();
  // open the file
  open_file();

  // reach end of comments
  std::string str;
  do {
    UG_ASSERT( !m_matFileStream.eof(), "Unexpected end of file." );
    std::getline( m_matFileStream, str );
    m_firstDataLine++;
  } while( str.at( 0 ) == '%' );

  // get next non-empty line
  while( str.empty() ) {
    UG_ASSERT( !m_matFileStream.eof(), "Unexpected end of file." );
    std::getline( m_matFileStream, str );
    m_firstDataLine++;
  }
  
  // split the line
  std::vector<std::string> entriesVec;
  if( m_mmTypeCode.is_sparse() ) {
    boost::algorithm::split( entriesVec, str, boost::is_any_of( " " ),
                             boost::algorithm::token_compress_on );
  } else {
    UG_THROW( "Other than sparse MatrixMarket matrices are not yet implemented." );
  }

  set_mat_dims( boost::lexical_cast<int>( entriesVec.at( 0 ) ),
                boost::lexical_cast<int>( entriesVec.at( 1 ) ),
                boost::lexical_cast<int>( entriesVec.at( 2 ) ) );

  // close the file
  close_file();
}


//template<typename matrix_type>
//  std::vector< std::vector<size_t> > determine_matrix_characteristics(const matrix_type &matrix );

void MatrixIOMtx::read_entry( size_t *m, size_t *n,
                              CPUAlgebra::matrix_type::value_type *val )
{
  PROFILE_FUNC();
  // make sure, the file is open
  UG_ASSERT( m_matFileStream.is_open(), "File is not open." );

  // get the data line and split it
  std::string line;
  std::getline( m_matFileStream, line );
  std::vector<std::string> elements;
  boost::algorithm::split( elements, line, boost::is_any_of(" "), boost::algorithm::token_compress_on );

  // parse it according to the matrix type
  if( is_sparse() ) {
    UG_ASSERT( elements.size() == 3,
               "Sparse matrix requires three values per line. Found: "
               << elements.size() );
    *m = boost::lexical_cast<size_t>( elements.at( 0 ) );
    *n = boost::lexical_cast<size_t>( elements.at( 1 ) );
    *val = boost::lexical_cast<CPUAlgebra::matrix_type::value_type>( elements.at( 2 ) );
  } else {
    UG_THROW( "Other than sparse MatrixMarket matrices are not yet implemented." );
  }
}

void MatrixIOMtx::write_banner()
{
  // Write the first lines (banner, comment, characteristics)
  open_file( std::ios_base::out | std::ios_base::trunc );

  // print out the banner
  m_matFileStream << MM_BANNER_STR << " " << MM_MTX_STR << " "
                  << MM_COORDINATE_STR << " " << MM_REAL_STR << " ";
  if ( m_mmTypeCode.is_general() ) {
    m_matFileStream << MM_GENERAL_STR;
  } else if ( m_mmTypeCode.is_symmetric() ) {
    m_matFileStream << MM_SYMMETRIC_STR;
  } else if ( m_mmTypeCode.is_skew_symmetric() ) {
    m_matFileStream << MM_SKEW_STR;
  }
  m_matFileStream << "\n";

  close_file();
}

void MatrixIOMtx::write_entry( size_t m, size_t n ,
                               CPUAlgebra::matrix_type::value_type val )
{
  PROFILE_FUNC();
  UG_ASSERT( m_matFileStream.is_open(), "File is not open." );
  UG_ASSERT( m > 0, "Row index not positive." );
  UG_ASSERT( n > 0, "Column index not positive." );
  m_matFileStream.unsetf( std::ios_base::scientific );
  m_matFileStream << m << " " << n;
  m_matFileStream << ( (val < 0) ? " " : "  " );
  m_matFileStream.setf( std::ios_base::scientific);
  m_matFileStream << std::setprecision(13);
  m_matFileStream << val << "\n";
  m_matFileStream.flush();
}

} // namespace ug

// EOF
