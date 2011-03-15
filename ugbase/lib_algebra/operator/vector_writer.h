/*
 * vector_writer.h
 *
 *  Created on: 10.03.2011
 */

#ifndef __H__LIB_ALGEBRA__OPERATOR__VECTOR_WRITER__
#define __H__LIB_ALGEBRA__OPERATOR__VECTOR_WRITER__

namespace ug{

template <typename vector_type>
class IVectorWriter
{
	public:
	///	write vector
		virtual bool update(vector_type &vec) = 0;

		/*virtual bool update() = 0;

		void clear();

		vector_type& get_vector(size_t i) {return m_vVec[i];}

		size_t num_vector() {return m_vVec.size();}

		number weight(size_t i) {..}*/

		/// virtual destructor
		virtual ~IVectorWriter(){}

	/*private:
		std::vector<vector_type> m_vVec;
		std::vector<number> m_vWeights;*/
};


template <size_t dim>
class IPositionProvider
{
public:
	virtual bool get_positions(std::vector<MathVector<dim> >&vec) = 0;
	virtual ~IPositionProvider() {}
};



} // end namespace ug

#endif /* __H__LIB_ALGEBRA__OPERATOR__VECTOR_WRITER__ */
