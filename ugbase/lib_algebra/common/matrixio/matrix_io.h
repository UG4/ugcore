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

#ifndef __H__UG__LIB_ALGEBRA__MATRIX_IO_H
#define __H__UG__LIB_ALGEBRA__MATRIX_IO_H

#include <string>
#include <fstream>

#include "common/util/file_util.h"
#include "common/assert.h"
#include "common/error.h"
#include "common/log.h"
#include "common/profiler/profiler.h"
#include "lib_algebra/lib_algebra.h"
#include "lib_algebra/cpu_algebra_types.h"

namespace ug
{

/**
 * \defgroup matrixio MatrixIO
 * \ingroup lib_algebra_common
 * \{
 */

/**
 * \brief Representation of a matrix exchange file format
 *
 * Supported matrix exchange file formats are:
 * - MatrixMarket (http://math.nist.gov/MatrixMarket/formats.html#MMformat)
 *
 * \note The Harwell-Boing (http://math.nist.gov/MatrixMarket/formats.html#hb)
 * format is not yet implemented.
 */
class MatrixFileType
{
  public:
    /// Supported matrix exchange file formats
    enum Type {
      /// undefined format
      UNDEFINED = 0,
      /// MatrixMarket format (by NIS)
      MATRIX_MARKET = 1,
      /// Harwell-Boing format (by NIS)
      HARWELL_BOING = 2   // TODO: HarwellBoing format not yet implemented
    };

    /**
     * \brief Default Constructor
     *
     * Initializes file format type to \c 0 (thus undefined).
     */
    MatrixFileType() : m_type( UNDEFINED ) {}
    /**
     * \brief Constructor specifying file type
     *
     * Initializes file format type to the given type.
     *
     * \param[in] type  a valid value of MatrixFileType::Type
     */
    MatrixFileType( int type ) : m_type() {
      UG_ASSERT( type >= 0 && type <= 2, "Unknown MatrixFileType code: " << type );
      m_type = type;
    }
    /**
     * \brief Tells whether this is a MatrixMarket type
     * \return \c true if it is a MatrixMarket file type, \c false otherwise.
     */
    bool is_mm() const {
      return m_type == MATRIX_MARKET;
    }
    /**
     * \brief Tells whether this is a Harwell-Boing type
     * \return \c true if it is a Harwell-Boing file type, \c false otherwise.
     */
    bool is_hb() const {
      return m_type == HARWELL_BOING;
    }

  private:
    /// Holds the file type value (which is one of MatrixFileType::Type)
    int m_type;
};

/**
 * \brief Generic matrix I/O functionality
 *
 * This is the generic class for matrix I/O functionality.
 * For actual matrix I/O operations use one of it's subclasses specifying the
 * specific behaviour for one of the supported matrix exchange file formats.
 *
 * \note Please include this header file, when you want to use the I/O functions
 * for a specific exchange file format.
 */
class MatrixIO
{
  public:
    /// Specifies how file names pointing to non-existing files should be handeld
    enum OpenMode {
      /// Only existing files are valid files (i.e. non-existing file result in error)
      EXISTING = 1,
      /// Non-existing files are created with the specified file name
      NEW = 2
    };
    
  private:
    /// Full path name of the matrix exchange file
    std::string *m_pMatFileName;
    /// Internal file stream for reading from and writing into the matrix exchange file
    std::fstream m_matFileStream;
    /// Matrix exchange file type
    MatrixFileType m_matFileType;
    /// Number of rows as specified in the matrix exchange file
    //int m_rows;
    /// Number of columns as specified in the matrix exchange file
    //int m_cols;
    /// Number of data lines as specified in the matrix exchange file
    //int m_lines;

  public:
    /**
     * \brief Default constructor
     */
    MatrixIO();
    /**
     * \brief Constructor specifying matrix exchange file path
     *
     * \param[in] mFile     path and name of matrix exchange file.
     * \param[in] openMode  how to deal with non-existing files (a value of
     *                      MatrixIO::OpenMode)
     */
    MatrixIO( std::string mFile, int openMode=EXISTING );
    /**
     * \brief Destructor
     *
     * The destructor calls MatrixIO::close_file to make sure, that write
     * operations to the associated exchange file are flushed and terminated
     * correctly.
     */
    ~MatrixIO();

    /**
     * \brief Sets associated file to another one
     * 
     * \param[in] mFile Full path and name of matrix exchange file.
     * \param[in] openMode  how to deal with non-existing files
     *
     * \throws std::runtime_error if exchange file is not found and
     *                            <tt>openMode=OpenMode::Existing</tt>.
     * \throws UGError            if \c openMode is non of OpenMode
     */
    void set_mat_file_name( std::string mFile, int openMode=EXISTING );
    /**
     * \brief Retreive associated exchange file path and name
     *
     * \return String representation of exchange file with full path as set by
     * constructor or MatrixIO::set_mat_file_name.
     */
    std::string get_mat_file_name() const;

  private:
    /**
     * \brief Opens the associated file in \c m_matFileStream
     *
     * \param[in] mode  Modus of the resulting file stream (a value of
     *                  std::ios_base::openmode). Defaults to std::ios_base::in.
     *
     * \throws std::runtime_error if \c m_pMatFileName is empty
     * \throws std::runtime_error if file stream can not be opened (i.e.
     *                            \c failbit or \c badbit are set)
     */
    void open_file( std::ios_base::openmode mode=std::ios_base::in );
    /**
     * \brief Closes file stream
     */
    void close_file();
};

// end group matrixio
/// \}

} // namespace ug

// include specializations
#include "matrix_io_mtx.h"

#endif // __H__UG__LIB_ALGEBRA__MATRIXIO_H

// EOF
