/*
 * constraints_post_process_interface.h
 *
 *  Created on: 01.03.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__POST_PROCESS__CONSTRAINTS__CONSTRAINTS_POST_PROCESS_INTERFACE__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__POST_PROCESS__CONSTRAINTS__CONSTRAINTS_POST_PROCESS_INTERFACE__

#include "lib_discretization/assemble.h"

namespace ug {

template <	typename TDiscreteFunction,
			typename TAlgebra = typename TDiscreteFunction::algebra_type >
class IConstraintsPostProcess{
	public:
		// discrete function type
		typedef TDiscreteFunction discrete_function_type;

		// algebra type
		typedef TAlgebra algebra_type;

		// type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

		// type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	public:
		virtual IAssembleReturn post_process_jacobian(matrix_type& J, const discrete_function_type& u)
		{return IAssemble_NOT_IMPLEMENTED;}

		virtual IAssembleReturn post_process_defect(vector_type& d, const discrete_function_type& u)
		{return IAssemble_NOT_IMPLEMENTED;}

		virtual IAssembleReturn post_process_linear(matrix_type& mat, vector_type& rhs, const discrete_function_type& u)
		{return IAssemble_NOT_IMPLEMENTED;}

		virtual ~IConstraintsPostProcess() {};
};

template <	typename TDiscreteFunction,
			typename TAlgebra = typename TDiscreteFunction::algebra_type >
class SymP1ConstraintsPostProcess : public IConstraintsPostProcess<TDiscreteFunction, TAlgebra> {
	public:
		// discrete function type
		typedef TDiscreteFunction discrete_function_type;

		// algebra type
		typedef TAlgebra algebra_type;

		// type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

		// type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	public:
		virtual IAssembleReturn post_process_linear(matrix_type& mat, vector_type& rhs, const discrete_function_type& u)
		{
			typename geometry_traits<ConstrainingEdge>::iterator iter, iterBegin, iterEnd;

			iterBegin = u.template begin<ConstrainingEdge>();
			iterEnd = u.template end<ConstrainingEdge>();

			for(iter = iterBegin; iter != iterEnd; ++iter)
			{
				ConstrainingEdge* bigEdge = *iter;

				VertexBase* vrt1 = bigEdge->vertex(0);
				VertexBase* vrt2 = bigEdge->vertex(1);
				VertexBase* vrt3 = NULL;

				for(EdgeBaseIterator edIter = bigEdge->constrained_edges_begin(); edIter !=  bigEdge->constrained_edges_end(); ++edIter)
				{
					ConstrainedEdge* smallEdge = dynamic_cast<ConstrainedEdge*>(*edIter);
					if(smallEdge == NULL)
					{
						UG_LOG("Cannot find Constrained Edge. Aborting.\n"); return IAssemble_ERROR;
					}

					vrt3 = smallEdge->vertex(0);
					if(vrt3 != vrt2 && vrt3 != vrt1) break;

					vrt3 = smallEdge->vertex(1);
					if(vrt3 != vrt2 && vrt3 != vrt1) break;
				}

				typename discrete_function_type::algebra_index_vector_type ind1, ind2, ind3;
				u.get_inner_algebra_indices(vrt1, ind1);
				u.get_inner_algebra_indices(vrt2, ind2);
				u.get_inner_algebra_indices(vrt3, ind3);

				if(false)
				{
					UG_LOG("ind1 = " << ind1[0] << "\n");
					UG_LOG("ind2 = " << ind2[0] << "\n");
					UG_LOG("ind3 = " << ind3[0] << "\n");
				}

				if(!SplitAddRow(mat, ind1, ind2, ind3))
					{UG_LOG("ERROR while splitting rows. Aborting.\n"); return IAssemble_ERROR;}

				if(!SetInterpolation(mat, ind1, ind2, ind3))
					{UG_LOG("ERROR while setting interpolation. Aborting.\n"); return IAssemble_ERROR;}

				if(!HandleRhs(rhs, ind1, ind2, ind3))
					{UG_LOG("ERROR while setting interpolation. Aborting.\n"); return IAssemble_ERROR;}

			}

			return IAssemble_OK;
		}

	protected:
		template <typename T>
		bool SplitAddRow(SparseMatrix<T>& A,  std::vector<size_t>& ind1, std::vector<size_t>& ind2, std::vector<size_t>& ind3)
		{
			if(ind1.size() != ind2.size() || ind1.size() != ind3.size())
				{UG_LOG("Wring number of indices. Cannot split row.\n"); return false;}

			for(size_t i = 0; i < ind1.size(); ++i)
			{
				for(typename SparseMatrix<T>::rowIterator conn = A.beginRow(ind3[i]); !conn.isEnd(); ++conn)
				{
					typename SparseMatrix<T>::entry_type block = (*conn).dValue;
					block *= 1./2.;
					const size_t j = (*conn).iIndex;

					A(ind1[i], j) += block;
					A(ind2[i], j) += block;
					A(ind3[i], j) = 0.0;
				}
			}
			return true;
		}

		template <typename T>
		bool SetInterpolation(SparseMatrix<T>& A,  std::vector<size_t>& ind1, std::vector<size_t>& ind2, std::vector<size_t>& ind3)
		{
			if(ind1.size() != ind2.size() || ind1.size() != ind3.size())
				{UG_LOG("Wring number of indices. Cannot split row.\n"); return false;}

			for(size_t i = 0; i < ind1.size(); ++i)
			{
					A(ind3[i], ind3[i]) = 1.0;
					A(ind3[i], ind1[i]) = -1./2.;
					A(ind3[i], ind2[i]) = -1./2.;
			}
			return true;
		}

		template <typename T>
		bool HandleRhs(Vector<T>& rhs,  std::vector<size_t>& ind1, std::vector<size_t>& ind2, std::vector<size_t>& ind3)
		{
			if(ind1.size() != ind2.size() || ind1.size() != ind3.size())
				{UG_LOG("Wring number of indices. Cannot split row.\n"); return false;}

			for(size_t i = 0; i < ind1.size(); ++i)
			{
				typename Vector<T>::entry_type& val = rhs[ind3[i]];

				rhs[ind1[i]] += val * 1./2.;
				rhs[ind2[i]] += val * 1./2.;
				val = 0.0;
			}
			return true;
		}


};


}; // namespace ug



#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__POST_PROCESS__CONSTRAINTS__CONSTRAINTS_POST_PROCESS_INTERFACE__ */
