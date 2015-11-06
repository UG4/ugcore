/*
 * Copyright (c) 2012-2013:  G-CSC, Goethe University Frankfurt
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

/**
 * \file matrix_io_mtx.h
 * \author Torbjoern Klatt
 * \date 2012-05-06
 *
 * \details For some parts the implementation of a pure C library for the
 * MatrixMarket exchange format was used as a source of inspiration.
 * These can be found at http://math.nist.gov/MatrixMarket/mmio-c.html
 */

#ifndef __H__UG__LIB_ALGEBRA__MATRIX_IO_MTX_H
#define __H__UG__LIB_ALGEBRA__MATRIX_IO_MTX_H

#include <sstream>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "lib_algebra/common/matrixio/matrix_io.h"
#include "lib_algebra/common/matrixio/mm_type_code.h"

namespace ug
{

/// \addtogroup matrixio
/// \{

/**
 * \brief Provides I/O functionality for MatrixMarket exchange file format
 *
 * Given the path to a \c .mtx file it queries the therein stored matrix type
 * and characteristics.
 * With <tt>readInto(matrixType &matrix)</tt> an easy way of reading the
 * represented matrix into a <tt>CPUAlgebra::matrix_type</tt> object is
 * provided.
 *
 * \note So far, only real sparse matrices either symmetric or skew-symmetric or
 * none of both are supported.
 * Support for dense matrices (either really dense or just in array format) will
 * come soon.
 * Non-real (i.e. complex, integer or pattern) matrices and hermitian matrices
 * might be supported in the future.
 *
 * The specification of the MatrixMarket matrix exchange format can be found
 * here: http://math.nist.gov/MatrixMarket/formats.html#MMformat
 *
 * \note Please include the generic MatrixIO header file
 * (lib_algebra/common/matrixio/matrix_io.h) for usage and not the header file
 * of this specific exchange file format.
 */
class MatrixIOMtx : private MatrixIO
{
  private:
    // Full path name of the matrix exchange file
    // (docu in matrix_io.h)
    std::string *m_pMatFileName;
    // Internal file stream for reading from and writing into the matrix
    // exchange file (docu in matrix_io.h)
    std::fstream m_matFileStream;
    /// Line number of the first data line (0 based)
    size_t m_firstDataLine;
    // Matrix exchange file type (set to constant MatrixFileType::MatrixMarket)
    // (docu in matrix_io.h)
    MatrixFileType m_matFileType;
    /// Characteristics of MatrixMarket matrix
    MMTypeCode m_mmTypeCode;
    // Number of rows as specified in the matrix exchange file
    // (docu in matrix_io.h)
    size_t m_rows;
    // Number of columns as specified in the matrix exchange file
    // (docu in matrix_io.h)
    size_t m_cols;
    // Number of data lines as specified in the matrix exchange file
    // (docu in matrix_io.h)
    size_t m_lines;

  public:
    /**
     * \brief Default constructor
     */
    MatrixIOMtx();
    /**
     * \brief Constructor specifying matrix exchange file path
     *
     * After successful call to MatrixIOMtx::set_mat_file_name it calls
     * MatrixIOMtx::query_matrix_type and
     * MatrixIOMtx::query_matrix_characteristics as long as openMode is
     * MatrixIO::OpenMode::EXISTING.
     *
     * \param[in] mFile     Full path and name of matrix exchange file.
     * \param[in] openMode  how to deal with non-existing files (a value of
     *                      MatrixIO::OpenMode)
     */
    MatrixIOMtx( std::string mFile, int openMode=MatrixIO::EXISTING );
    
    /**
     * \brief Destructor
     *
     * The destructor calls MatrixIOMtx::close_file to make sure, that write
     * operations to the associated exchange file are flushed and terminated
     * correctly.
     */
    ~MatrixIOMtx();

    /**
     * \brief Sets associated file to another one
     *
     * \note It does not check, whether the given file is a MatrixMarket
     *       exchange file.
     *       Only it's existance is verified.
     *
     * \param[in] mFile Full path and name of matrix exchange file.
     * \param[in] openMode  how to deal with non-existing files
     *
     * \throws std::runtime_error if exchange file is not found and
     *                            <tt>openMode=MatrixIO::OpenMode::Existing</tt>.
     * \throws UGError            if \c openMode is non of MatrixIO::OpenMode
     */
    void set_mat_file_name( std::string mFile, int openMode=MatrixIO::EXISTING );
    
    /**
     * \brief Retreive associated exchange file path and name
     *
     * \return String representation of exchange file with full path as set by
     * constructor or MatrixIOMtx::set_mat_file_name.
     */
    std::string get_mat_file_name() const;

    /**
     * \brief Sets the size of the matrix and number of data lines
     *
     * In case of reading from a matrix exchange file, this function is called
     * by MatrixIOMtx::query_matrix_characteristics.
     *
     * \param[in] rows  Number of rows of the matrix
     * \param[in] cols  Number of columns of the matrix
     * \param[in] lines Number of data lines in the exchange file
     * 
     * \throws std::runtime_error if \c rows is not positive
     * \throws std::runtime_error if \c cols is not positive
     * \throws std::runtime_error if \c lines is not positive
     */
    void set_mat_dims( size_t rows, size_t cols, size_t lines );

    /**
     * \brief Retrieves number of rows as specified in the exchange file
     *
     * \return Number of rows as \c size_t
     */
    size_t get_num_rows() const;

    /**
     * \brief Retrieves number of columns as specified in the exchange file
     *
     * \return Number of columns as \c size_t
     */
    size_t get_num_cols() const;

    /**
     * \brief Retrieves number of data lines as specified in the exchange file
     *
     * \return Number of data lines as \c size_t
     */
    size_t get_num_lines() const;

    /**
     * \brief Tells whether the associated exchange file holds a sparse matrix
     *
     * \return \c true if it is in coordinate format, \c false otherwise
     */
    bool is_sparse() const;

    /**
     * \brief Reads data from associated exchange file into ug matrix
     *
     * The given ug-matrix is first resized.
     * Each data line of the underlying exchange file is read separately and
     * independently.
     * If the matrix is (skew-)symmetric as defined in the exchange file, this
     * is considered during writing the values to the correct positions.
     *
     * \note This function does not take care about memory alignement of the
     *       given ug-matrix (e.g. SparseMatrix::defragment ).
     *
     * \param[in,out] matrix  Matrix of CPUAlgebra::matrix_type to read into
     *
     * \throws std::runtime_error if one of \c m_rows, \c m_cols or \c m_lines
     *                            is not positive.
     * \throws UGError            if exchange file's matrix is not in coordinate
     *                            format (as only sparse matrices are supported yet)
     */
    template<typename matrix_type>
    void read_into( matrix_type &matrix )
    {
      PROFILE_FUNC();
      UG_ASSERT( m_rows > 0 || m_cols > 0 || m_lines > 0,
                 "MatrixMarket matrix dimensions seem not be reasonable: (m, n, nnz)=("
                 << m_rows << "," << m_cols << "," << m_lines << ")" );

      // open the file and go one before first data line
      open_file();

      std::string dummy;
      for( size_t i = 0; i < m_firstDataLine; i++ ) {
        std::getline( m_matFileStream, dummy );
      }

      matrix.resize_and_clear( m_rows, m_cols );
      size_t x, y;
      double val;
      if ( m_mmTypeCode.is_sparse() ) {
        for( size_t i = 0; i < m_lines; i++ ) {
          read_entry( &x, &y, &val );
          // we need to substract 1 from row and column indices as MM is 1-based
          matrix( x - 1, y - 1 ) = val;
          if ( m_mmTypeCode.is_symmetric() && x != y ) {
            matrix( y - 1, x - 1 ) = val;
          } else if ( m_mmTypeCode.is_skew_symmetric() && x != y ) {
            matrix( y - 1, x - 1 ) = -val;
          }
        }
      } else {
        UG_THROW( "Other than sparse MatrixMarket matrices are not yet implemented." );
      }

      // close the file
      close_file();
    }


    /**
     * \brief Writes a ug matrix into the associated MatrixMarket exchange file
     *
     * From this analysis the MatrixMarket exchange file banner is constructed
     * and the characteristics and data lines written to that file consecutively.
     *
     * A comment is added right after the banner.
     * The default comment is <tt>Generated with ug4.</tt>.
     * To suppress the comment line completely, pass an empty string as the
     * \c comment parameter to this function.
     * 
     * \note Please make sure, the lines of the given comment are not longer
     *       than 1025 characters and all liness start with a \c %.
     *
     * \todo Implement validity check for the comment.
     *
     * The resulting MatrixMarket exchange file is in line with the
     * specifications.
     *
     * \param[in] matrix  Matrix of CPUAlgebra::matrix_type to read from
     * \param[in] comment Optional comment to be put atop of the MatrixMarket
     *                    matrix exchange file.
     */
    template<typename matrix_type>
    void write_from( matrix_type &matrix,
                     std::string comment="%Generated with ug4." )
    {
	  PROFILE_FUNC();

	  // analyze the given matrix
	  std::vector< std::vector<size_t> > rowIndexPerCol = determine_matrix_characteristics( matrix );

	  // write the MatrixMarket banner
	  write_banner();

	  // open the file for appending
	  open_file( std::ios_base::out | std::ios_base::app );

	  // add a comment if it's not empty
	  if ( !comment.empty() ) {
		if ( comment.find_first_of( '%' ) != 0 ) {
		  UG_WARNING( "Given comment did not start with '%'. Prepending it to make it valid." );
		  comment.insert( 0, "%" );
		}
		m_matFileStream << comment << "\n";
	  }

	  // write characteristics
	  m_matFileStream << m_rows << " " << m_cols << " " << m_lines << "\n";


	  // write entries to the file
	  for ( size_t col = 0; col < m_cols; col++ ) {
		for (size_t row = 0; row < rowIndexPerCol.at(col).size(); row++ ) {
		  // we need to add 1 to row and column index, as MM indices are 1-based
		  write_entry( rowIndexPerCol.at(col).at(row) + 1, col + 1,
					   matrix(rowIndexPerCol.at(col).at(row), col) );
		}
	  }
	  close_file();
    }


  private:
    // Opens file as fstream. (docu in matrix_io.h)
    void open_file( std::ios_base::openmode mode=std::ios_base::in );
    
    // Closes fstream. (docu in matrix_io.h)
    void close_file();

    /**
     * \brief Parses banner of MatrixMarket file
     *
     * The banner is the first line of a MatrixMarket file with the syntax:
     *
     * <tt>%%MatrixMarket matrix [format] [numeric] [type]</tt>
     *
     * <dl>
     * <dt><tt>format</tt></dt>
     * <dd>Either \c coordinate or \c array, representing a sparse or dense
     * matrix respectively</dd>
     * <dt><tt>numeric</tt></dt>
     * <dd>Either \c real, \c complex, \c integer or \c pattern,
     * representing the numerical type of the values</dd>
     * <dt><tt>type</tt></dt>
     * <dd>Either \c general, \c symmetric, \c skew-symmetric or
     * \c hermitian, representing the algebraic type of the matrix</dd>
     * </dl>
     * 
     * \throws std::runtime_error if the banner starts not with
     *                            <tt>%%MatrixMarket</tt>
     * \throws std::runtime_error if the second word of the banner is not
     *                            \c matrix
     * \throws UGError            if the \c format identifier is invalid
     * \throws UGError            if the \c numeric identifier is invalid
     * \throws UGError            if the \c type identifier is invalid
     */
    void query_matrix_type();

    /**
     * \brief Retreives matrix size and number non-zeros from exchange file
     *
     * The first non-empty line after the comments (lines starting with \c %)
     * consists of three values for the coordinate format:
     *
     * <tt>nRows nCols nLines</tt>
     *
     * <dl>
     * <dt><tt>nRows</tt></dt><dd>Number of rows of the described matrix</dd>
     * <dt><tt>nCols</tt></dt><dd>Number of columns of the described matrix</dd>
     * <dt><tt>nLines</tt></dt><dd>Number of data lines. For types \c symmetric,
     * \c skew-symmetric and \c hermitian, this is not equal to the number of
     * non-zero elements.</dd>
     * </dl>
     *
     * \throws std::runtime_error if the file stream ended unexpectetly (e.g.
     *                            the exchange file does not contain data lines)
     * \throws UGError            if the exchange file does not describe a
     *                            coordinate/sparse matrix (as only sparse
     *                            matrices are supported yet)
     */
    void query_matrix_characteristics();

    /**
     * \brief Determines characteristics of given matrix
     *
     * Symmetries of the given matrix are determined and saved.
     *
     * \param[in] matrix  Matrix of CPUAlgebra::matrix_type to analyse
     * \return 2D vector with coordinates of non-zero entries of the given matrix
     *         with column indices as first and row indices as second dimension.
     *
     * \throws std::runtime_error if the matrix was determined to be symmetric
     *                            and skew-symmetric at the same time
     */


    template<typename matrix_type>
    std::vector< std::vector<size_t> > determine_matrix_characteristics(const matrix_type &matrix )
    {
      PROFILE_FUNC();

      // read matrix dimensions
      size_t rows = matrix.num_rows();
      size_t cols = matrix.num_cols();
      size_t offDiagEntries = 0;
      size_t diagEntries = 0;

      // As MatrixMarket specifies a column-first-ordering, we need to get to know
      // in which rows of each column non-zero entries are.
      // During this, we also get to know how many non-zero entries there are in
      // total and can detect symmetries of the matrix.
      std::vector< std::vector<size_t> > rowIndexPerCol;
      rowIndexPerCol.resize(cols);
      bool isSymmetric = true;
      bool isSkew = true;
      bool changed = false;
      for ( size_t r = 0; r < rows; r++ ) { // iterate rows
        for (typename matrix_type::const_row_iterator conn = matrix.begin_row(r);
              conn != matrix.end_row(r); ++conn ) {
          if ( conn.value() != 0.0 ) {
            // first add index to list
            if ( isSymmetric || isSkew ) {
              if ( conn.index() <= r ) {
                rowIndexPerCol.at( conn.index() ).push_back(r);
              }
            } else {
              rowIndexPerCol.at( conn.index() ).push_back(r);
            }

            // increment counters
            (conn.index() == r) ? diagEntries++ : offDiagEntries++ ;

            // check, whether it's still a symmetric or skew-symmetric matrix
            if ( r != conn.index() ) {
              if ( matrix(r, conn.index() ) != matrix( conn.index(), r ) ) {
                if ( isSymmetric ) {
                  changed = true;
                }
                isSymmetric = false;
              }
              if ( matrix(r, conn.index() ) != -1.0 * matrix( conn.index(), r ) ) {
                if ( isSkew ) {
                  changed = true;
                }
                isSkew = false;
              }
            }

            // We assumed the matrix to be symmetric or skew-symmetric, but it's
            // not. Thus we need to redo everything done before.
            if ( changed ) {
              offDiagEntries = 0;
              diagEntries = 0;
              rowIndexPerCol.clear();
              rowIndexPerCol.resize( cols );
              r = -1; // at the end of the current row-loop this is incremented again
              changed = false;
              break;
            }
          }
        }
      }

      // make sure the matrix is not both, symmetric and skew-symmetric
      UG_ASSERT( ( ( isSymmetric && !isSkew ) ||
                   ( isSkew && !isSymmetric ) ||
                   ( !isSkew && !isSymmetric ) ),
                 "Error on detecting symmetry of matrix. (skew:" << isSkew
                 << " sym:" << isSymmetric << ")" );

      // set MMTypeCode
      size_t entries = 0;
      m_mmTypeCode.set_class_type( MMTypeCode::COORDINATE );
      m_mmTypeCode.set_numeric_type( MMTypeCode::REAL );
      if ( isSymmetric ) {
        entries = ( offDiagEntries / 2 )+ diagEntries;
        m_mmTypeCode.set_algebraic_type( MMTypeCode::SYMMETRIC );
      } else if ( isSkew ) {
        entries = ( offDiagEntries / 2 ) + diagEntries;
        m_mmTypeCode.set_algebraic_type( MMTypeCode::SKEW );
      } else {
        entries = offDiagEntries + diagEntries;
        m_mmTypeCode.set_algebraic_type( MMTypeCode::GENERAL );
      }

      // now we can set the dimensions
      set_mat_dims( rows, cols, entries );

      return rowIndexPerCol;
    }

    /**
     * \brief Reads and parses the next data line
     *
     * \note This function does not open the filestream by itself, but requires
     *       it to be already open.
     * 
     * \param[out]  m   pointer to store the row index to (1 based)
     * \param[out]  n   pointer to store the column index to (1 based)
     * \param[out]  val pointer to store the value to
     *
     * \throws std::runtime_error if the filestream is not open
     * \throws std::runtime_error if the matrix is sparse but there were less or
     *                            more than 3 elements in the line
     * \throws ug::UGError        if the exchange file does not describe a
     *                            coordinate/sparse matrix (as only sparse
     *                            matrices are supported yet)
     */
    void read_entry( size_t *m, size_t *n, CPUAlgebra::matrix_type::value_type *val );

    /**
     * \brief Writes the banner to the associated MatrixMarket file
     *
     * \note Any existing content of the associated file will be delted.
     */
    void write_banner();
    
    /**
     * \brief Appends a data line to the associated file
     *
     * The value is formatted to scientific notation with 13 digits of precision.
     *
     * \param[out]  m   row index (1 based)
     * \param[out]  n   column index (1 based)
     * \param[out]  val value at [row,col]
     *
     * \throws std::runtime_error if the filestream is not open
     * \throws std::runtime_error if either the specified row or column index
     *                            is not positive
     */
    void write_entry( size_t m, size_t n, CPUAlgebra::matrix_type::value_type val );
};

// end group matrixio
/// \}

} // namespace ug

#endif // __H__UG__LIB_ALGEBRA__MATRIX_IO_MTX_H

// EOF
