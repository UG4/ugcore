/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Torbjörn Klatt
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

#include "matrix_io.h"

namespace ug
{

// /////////////////////////////////////////////////////////////////////////////
// Constructors / Destructor
MatrixIO::MatrixIO() :
  m_pMatFileName( nullptr ), m_matFileStream(), m_matFileType( 0 )
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
