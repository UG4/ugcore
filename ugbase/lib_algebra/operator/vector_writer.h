/*
 * vector_writer.h
 *
 *  Created on: 10.03.2011
 */

#ifndef __H__LIB_ALGEBRA__OPERATOR__VECTOR_WRITER__
#define __H__LIB_ALGEBRA__OPERATOR__VECTOR_WRITER__

namespace ug{

/// Interface for modifying a vector (e.g, setting Dirichlet values, ...)
/** sa GridFunctionVectorWriterDirichlet0*/
template <typename vector_type>
class IVectorWriter
{
	public:
	///	write vector
		virtual bool update(vector_type &vec) = 0;

		/// virtual destructor
		virtual ~IVectorWriter(){}

	/*private:
		std::vector<vector_type> m_vVec;
		std::vector<number> m_vWeights;*/
};





} // end namespace ug

#endif /* __H__LIB_ALGEBRA__OPERATOR__VECTOR_WRITER__ */
