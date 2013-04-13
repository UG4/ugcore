/**
 * \file mm_type_code.h
 * \author Torbjoern Klatt
 * \date 2012-05-06
 *
 * \details For some parts the implementation of a pure C library for the
 * MatrixMarket exchange format was used as a source of inspiration.
 * These can be found at http://math.nist.gov/MatrixMarket/mmio-c.html
 */

#ifndef __H__UG__LIB_ALGEBRA__MM_TYPE_CODE_H
#define __H__UG__LIB_ALGEBRA__MM_TYPE_CODE_H

#include <string>

#include "common/error.h"

namespace ug
{

/// \addtogroup matrixio
/// \{

/**
 * \brief Type representation for MatrixMarket matrix exchange files
 */
class MMTypeCode
{
  public:
    /**
     * \brief Class type of the described matrix
     *
     * \note Only matrices in coordinate format are supported yet.
     */
    enum ClassType {
      COORDINATE        = 1,
      ARRAY             = 2     // TODO array/dense matrix not yet implemented
    };

    /**
     * \brief Numeric type of the described matrix
     *
     * \note Only real matrices are supported yet.
     */
    enum NumericType {
      REAL              = 1,
      COMPLEX           = 2,    // TODO complex matrices not yet implemented
      INTEGER           = 3,    // TODO integer matrices not yet implemented
      PATTERN           = 4     // TODO pattern matrices not yet implemented
    };

    /**
     * \brief Algebraic type of the described matrix
     *
     * \note Hermitian matrices are not yet supported.
     */
    enum AlgebraicType {
      GENERAL           = 1,
      SYMMETRIC         = 2,
      SKEW              = 3,
      HERMITIAN         = 4     // TODO hermitian matrices not yet implemented
    };

    /// Maximum line length in characters as defined by the MatrixMarket specifications
    static const int MM_LINE_LENGTH = 1025;
#define MM_BANNER_STR      "%%MatrixMarket"
#define MM_MTX_STR         "matrix"
#define MM_COORDINATE_STR  "coordinate"
#define MM_SPARSE_STR      "coordinate"
#define MM_ARRAY_STR       "array"
#define MM_DENSE_STR       "array"
#define MM_REAL_STR        "real"
#define MM_COMPLEX_STR     "complex"
#define MM_INTEGER_STR     "integer"
#define MM_PATTERN_STR     "pattern"
#define MM_GENERAL_STR     "general"
#define MM_SYMMETRIC_STR   "symmetric"
#define MM_SKEW_STR        "skew-symmetric"
#define MM_HERMITIAN_STR   "hermitian"

  private:
    /// Holds a value of MMTypeCode::ClassType (initialized with 0)
    int m_class_type;
    /// Holds a value of MMTypeCode::NumericType (initialized with 0)
    int m_numeric_type;
    /// Holds a value of MMTypeCode::AlgebraicType (initialized with 0)
    int m_algebraic_type;

  public:
    /**
     * \brief Default constructor
     *
     * Initializes all three types to \c 0, thus undefined.
     */
    MMTypeCode() :
      m_class_type( 0 ), m_numeric_type( 0 ), m_algebraic_type( 0 )
    {}

    /**
     * \brief Pretty prints the three classes and their current value.
     *
     * <tt>MatrixMarket Type Codes:
     *   ClassType:     {0,1,2}
     *   NumericType:   {0,1,2,3,4}
     *   AlgebraicType: {0,1,2,3,4}
     * </tt>
     * 
     * \return String
     *
     * \todo Display corresponding strings instead of internal numeric codes.
     */
    std::string to_string() {
      std::stringstream out;
      out << "MatrixMarket Type Codes:\n";
      out << "  ClassType:     " << m_class_type << "\n";
      out << "  NumericType:   " << m_numeric_type << "\n";
      out << "  AlgebraicType: " << m_algebraic_type;
      return out.str();
    }

    /**
     * \brief Tells whether MMTypeCode is sparse
     * 
     * \return \c true if \c ClassType is sparse (i.e. 1), \c false otherwise.
     */
    bool is_sparse() const {
      return ( m_class_type == COORDINATE );
    }

    /**
     * \brief Tells whether MMTypeCode is real
     *
     * \return \c true if \c NumericType is real (i.e. 1), \c false otherwise.
     */
    bool is_real() const {
      return ( m_numeric_type == REAL );
    }

    /**
     * \brief Tells whether MMTypeCode is general
     *
     * \return \c true if \c AlgebraicType is general (i.e. 1), \c false otherwise.
     */
    bool is_general() const {
      return ( m_algebraic_type == GENERAL );
    }
    /**
     * \brief Tells whether MMTypeCode is symmetric
     *
     * \return \c true if \c AlgebraicType is symmetric (i.e. 2), \c false otherwise.
     */
    bool is_symmetric() const {
      return ( m_algebraic_type == SYMMETRIC );
    }
    /**
     * \brief Tells whether MMTypeCode is skew-symmetric
     *
     * \return \c true if \c AlgebraicType is skew-symmetric (i.e. 3), \c false otherwise.
     */
    bool is_skew_symmetric() const {
      return ( m_algebraic_type == SKEW );
    }

    /**
     * \brief Sets a new class type from an enum value
     *
     * \param[in] type One value of MMTypeCode::ClassType
     * \throws UGError if \c type is not a valid value
     */
    void set_class_type( int type ) {
      switch( type ) {
        case 1:
          m_class_type = COORDINATE;
          break;
        case 2:
          m_class_type = ARRAY;
          break;
        default:
          UG_THROW( "MatrixMarket Type code is invalid: " << type );
      }
    }
    /**
     * \brief Sets a new class type from a string
     *
     * \param[in] type a string either \c coordinate or \c array
     * \throws UGError if \c type is another string
     */
    void set_class_type( std::string type ) {
      if( type.compare( MM_COORDINATE_STR ) == 0 || type.compare( MM_SPARSE_STR ) == 0 ) {
        set_class_type( COORDINATE );
      } else if( type.compare( MM_ARRAY_STR ) == 0 || type.compare( MM_DENSE_STR ) == 0 ) {
        set_class_type( ARRAY );
      } else {
        UG_THROW( "MatrixMarket Type is invalid: " << type );
      }
    }

    /**
     * \brief Sets a new numeric type from an enum value
     *
     * \param[in] type One value of MMTypeCode::NumericType
     * \throws UGError if \c type is not a valid value
     */
    void set_numeric_type( int type ) {
      switch( type ) {
        case 1:
          m_numeric_type = REAL;
          break;
        case 2:
          m_numeric_type = COMPLEX;
          break;
        case 3:
          m_numeric_type = INTEGER;
          break;
        case 4:
          m_numeric_type = PATTERN;
          break;
        default:
          UG_THROW( "MatrixMarket numeric type code is invalid: " << type );
      }
    }
    /**
     * \brief Sets a new numeric type from a string
     *
     * \param[in] type a string either \c real, \c complex, \c integer or
     *                 \c pattern
     * \throws UGError if \c type is another string
     */
    void set_numeric_type( std::string type ) {
      if( type.compare( MM_REAL_STR ) == 0 ) {
        set_numeric_type( REAL );
      } else if( type.compare( MM_COMPLEX_STR ) == 0 ) {
        set_numeric_type( COMPLEX );
      } else if( type.compare( MM_INTEGER_STR ) == 0 ) {
        set_numeric_type( INTEGER );
      } else if( type.compare( MM_PATTERN_STR ) == 0 ) {
        set_numeric_type( PATTERN );
      } else {
        UG_THROW( "MatrixMarket numeric type is invalid: " << type );
      }
    }

    /**
     * \brief Sets a new algebraic type from an enum value
     *
     * \param[in] type One value of MMTypeCode::AlgebraicType
     * \throws UGError if \c type is not a valid value
     */
    void set_algebraic_type( int type ) {
      switch( type ) {
        case 1:
          m_algebraic_type = GENERAL;
          break;
        case 2:
          m_algebraic_type = SYMMETRIC;
          break;
        case 3:
          m_algebraic_type = SKEW;
          break;
        case 4:
          m_algebraic_type = HERMITIAN;
          break;
        default:
          UG_THROW( "MatrixMarket algebraic type code is invalid: " << type );
      }
    }
    /**
     * \brief Sets a new algebraic type from a string
     *
     * \param[in] type a string either \c general, \c symmetric,
     *                 \c skew-symmetric or \c hermitian
     * \throws UGError if \c type is another string
     */
    void set_algebraic_type( std::string type ) {
      if( type.compare( MM_GENERAL_STR ) == 0 ) {
        set_algebraic_type( GENERAL );
      } else if( type.compare( MM_SYMMETRIC_STR ) == 0 ) {
        set_algebraic_type( SYMMETRIC );
      } else if( type.compare( MM_SKEW_STR ) == 0 ) {
        set_algebraic_type( SKEW );
      } else if( type.compare( MM_HERMITIAN_STR ) == 0 ) {
        set_algebraic_type( HERMITIAN );
      } else {
        UG_THROW( "MatrixMarket algebraic type is invalid: " << type );
      }
    }
};

// end group matrixio
/// \}

} // namespace ug

#endif // __H__UG__LIB_ALGEBRA__MM_TYPE_CODE_H

// EOF
