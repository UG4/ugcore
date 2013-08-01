/*
 * algebra_types.h
 *
 *  Created on: 01.04.2011
 *      Author: mrupp
 */

#ifndef __H__UG__LIB_ALGEBRA__ALGEBRA_TYPE__
#define __H__UG__LIB_ALGEBRA__ALGEBRA_TYPE__

#include <ostream>

namespace ug{

/// \addtogroup lib_algebra lib_algebra
/// \{

////////////////////////////////////////////////////////////////////////////////
//   Algebra Types
////////////////////////////////////////////////////////////////////////////////

/// class describing the type of an algebra
class AlgebraType
{
	public:
	///	types of algebra
		enum Type
		{
			CPU = 0,
			GPU = 1
		};

	///	indicating variable block size
		enum {VariableBlockSize = -1};

	public:
	///	empty constructor
		AlgebraType();

	///	constructor for fix blocksize
		AlgebraType(Type type, int blockSize);

	///	constructor for fix blocksize
		AlgebraType(const char* type, int blockSize);

	///	constructor for variable block size
		AlgebraType(const char* type);

	///	returns the type
		int type() const {return m_type;}

	///	returns the blocksize
		int blocksize() const {return m_blockSize;}

	protected:
		int m_type;
		int m_blockSize;
};

/// writes the Identifier to the output stream
std::ostream& operator<<(std::ostream& out,	const AlgebraType& v);


/// Singleton, providing the current default algebra.
class DefaultAlgebra
{
	public:
		static AlgebraType get() {return m_default;}
		static void set(const AlgebraType& defaultType) {m_default = defaultType;}

	protected:
		static AlgebraType m_default;
};

// end group lib_algebra
/// \}

} // end namespace ug

#endif /* __H__UG__LIB_ALGEBRA__ALGEBRA_TYPE__ */
