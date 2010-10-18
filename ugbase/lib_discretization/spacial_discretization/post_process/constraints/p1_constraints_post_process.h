/*
 * p1_constraints_post_process.h
 *
 *  Created on: 01.03.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__POST_PROCESS__CONSTRAINTS__P1_CONSTRAINTS_POST_PROCESS__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__POST_PROCESS__CONSTRAINTS__P1_CONSTRAINTS_POST_PROCESS__

#include "lib_discretization/assemble.h"

namespace ug {

template <	typename TDoFDistribution,
			typename TAlgebra>
class SymP1ConstraintsPostProcess : public IPostProcess<TDoFDistribution, TAlgebra> {
	public:
	// 	DoF Distribution Type
		typedef TDoFDistribution dof_distribution_type;

	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	// 	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	public:
		virtual int type() {return PPT_CONSTRAINTS;}
		virtual IAssembleReturn post_process_linear(matrix_type& mat, vector_type& rhs, const vector_type& u, const dof_distribution_type& dofDistr, number time = 0.0)
		{
			typename geometry_traits<ConstrainingEdge>::iterator iter, iterBegin, iterEnd;

			iterBegin = dofDistr.template begin<ConstrainingEdge>();
			iterEnd = dofDistr.template end<ConstrainingEdge>();

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

				typename dof_distribution_type::algebra_index_vector_type ind1, ind2, ind3;
				dofDistr.get_inner_algebra_indices(vrt1, ind1);
				dofDistr.get_inner_algebra_indices(vrt2, ind2);
				dofDistr.get_inner_algebra_indices(vrt3, ind3);

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



template <	typename TDoFDistribution,
			typename TAlgebra>
class OneSideP1ConstraintsPostProcess : public IPostProcess<TDoFDistribution, TAlgebra> {
	public:
	// 	DoF Distribution Type
		typedef TDoFDistribution dof_distribution_type;

	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	// 	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	// 	Type of algebra index vector
		typedef typename dof_distribution_type::algebra_index_vector_type algebra_index_vector_type;

	public:
		virtual int type() {return PPT_CONSTRAINTS;}

		virtual IAssembleReturn post_process_linear(matrix_type& mat, vector_type& rhs, const vector_type& u, const dof_distribution_type& dofDistr, number time = 0.0)
		{
			std::vector<algebra_index_vector_type> vConstrainingIndices;
			algebra_index_vector_type constrainedIndex;
			std::vector<VertexBase*> vConstrainingVertices;
			VertexBase* constrainedVertex = NULL;

			//////////////////////////
			// Constrained Edges
			//////////////////////////
			typename geometry_traits<ConstrainingEdge>::const_iterator iter, iterBegin, iterEnd;

			iterBegin = dofDistr.template begin<ConstrainingEdge>();
			iterEnd = dofDistr.template end<ConstrainingEdge>();

			vConstrainingIndices.resize(2);
			vConstrainingVertices.resize(2);
			for(iter = iterBegin; iter != iterEnd; ++iter)
			{
				ConstrainingEdge* bigEdge = *iter;

				for(size_t i=0; i < 2; ++i)
					vConstrainingVertices[i] = bigEdge->vertex(i);

				for(EdgeBaseIterator edIter = bigEdge->constrained_edges_begin(); edIter !=  bigEdge->constrained_edges_end(); ++edIter)
				{
					ConstrainedEdge* smallEdge = dynamic_cast<ConstrainedEdge*>(*edIter);
					if(smallEdge == NULL)
					{
						UG_LOG("Cannot find Constrained Edge. Aborting.\n"); return IAssemble_ERROR;
					}

					constrainedVertex = smallEdge->vertex(0);
					if(constrainedVertex != vConstrainingVertices[1] && constrainedVertex != vConstrainingVertices[0]) break;

					constrainedVertex = smallEdge->vertex(1);
					if(constrainedVertex != vConstrainingVertices[1] && constrainedVertex != vConstrainingVertices[0]) break;
				}

				// get algebra indices
				for(size_t i=0; i < 2; ++i)
					dofDistr.get_inner_algebra_indices(vConstrainingVertices[i], vConstrainingIndices[i]);

				dofDistr.get_inner_algebra_indices(constrainedVertex, constrainedIndex);

				// Split using indices
				if(!SplitAddRow(mat, constrainedIndex, vConstrainingIndices))
					{UG_LOG("ERROR while splitting rows. Aborting.\n"); return IAssemble_ERROR;}

				if(!SetInterpolation(mat, constrainedIndex, vConstrainingIndices))
					{UG_LOG("ERROR while setting interpolation. Aborting.\n"); return IAssemble_ERROR;}

				if(!HandleRhs(rhs, constrainedIndex, vConstrainingIndices))
					{UG_LOG("ERROR while setting interpolation. Aborting.\n"); return IAssemble_ERROR;}

			}

			//////////////////////////
			// Constrained Quads
			//////////////////////////
			typename geometry_traits<ConstrainingQuadrilateral>::const_iterator iterQuad, iterQuadBegin, iterQuadEnd;

			iterQuadBegin = dofDistr.template begin<ConstrainingQuadrilateral>();
			iterQuadEnd = dofDistr.template end<ConstrainingQuadrilateral>();

			vConstrainingIndices.resize(4);
			vConstrainingVertices.resize(4);
			for(iterQuad = iterQuadBegin; iterQuad != iterQuadEnd; ++iterQuad)
			{
				ConstrainingQuadrilateral* bigQuad = *iterQuad;

				for(size_t i=0; i < 4; ++i)
					vConstrainingVertices[i] = bigQuad->vertex(i);

				UG_ASSERT(bigQuad->num_constrained_vertices() == 1, "Should only be one hanging vertex in current implementation");

				constrainedVertex = *(bigQuad->constrained_vertices_begin());

				// get algebra indices
				for(size_t i=0; i < 4; ++i)
					dofDistr.get_inner_algebra_indices(vConstrainingVertices[i], vConstrainingIndices[i]);

				dofDistr.get_inner_algebra_indices(constrainedVertex, constrainedIndex);

				// Split using indices
				if(!SplitAddRow(mat, constrainedIndex, vConstrainingIndices))
					{UG_LOG("ERROR while splitting rows. Aborting.\n"); return IAssemble_ERROR;}

				if(!SetInterpolation(mat, constrainedIndex, vConstrainingIndices))
					{UG_LOG("ERROR while setting interpolation. Aborting.\n"); return IAssemble_ERROR;}

				if(!HandleRhs(rhs, constrainedIndex, vConstrainingIndices))
					{UG_LOG("ERROR while setting interpolation. Aborting.\n"); return IAssemble_ERROR;}

			}

			return IAssemble_OK;
		}

	protected:
		template <typename T>
		bool SplitAddRow(SparseMatrix<T>& A,  algebra_index_vector_type& constrainedIndex, std::vector<algebra_index_vector_type>& vConstrainingIndices)
		{
			for(size_t i = 0; i < vConstrainingIndices.size(); ++i)
			{
				if(vConstrainingIndices[i].size() != constrainedIndex.size())
					{UG_LOG("Wring number of indices. Cannot split row.\n"); return false;}
			}

			for(size_t i = 0; i < constrainedIndex.size(); ++i)
			{
				for(typename SparseMatrix<T>::rowIterator conn = A.beginRow(constrainedIndex[i]); !conn.isEnd(); ++conn)
				{
					typename SparseMatrix<T>::entry_type block = (*conn).dValue;
					const size_t j = (*conn).iIndex;

					// choose randomly the first dof to add whole row
					A(vConstrainingIndices[0][i], j) += block;
					A(constrainedIndex[i], j) = 0.0;
				}
			}
			return true;
		}

		template <typename T>
		bool SetInterpolation(SparseMatrix<T>& A,  algebra_index_vector_type& constrainedIndex, std::vector<algebra_index_vector_type>& vConstrainingIndices)
		{
			for(size_t i = 0; i < vConstrainingIndices.size(); ++i)
			{
				if(vConstrainingIndices[i].size() != constrainedIndex.size())
					{UG_LOG("Wring number of indices. Cannot split row.\n"); return false;}
			}

			const number scale = -1./(vConstrainingIndices.size());
			for(size_t i = 0; i < constrainedIndex.size(); ++i)
			{
					A(constrainedIndex[i], constrainedIndex[i]) = 1.0;
					for(size_t j = 0; j < vConstrainingIndices.size(); ++j)
					{
						A(constrainedIndex[i], vConstrainingIndices[j][i]) = scale;
					}
			}
			return true;
		}

		template <typename T>
		bool HandleRhs(Vector<T>& rhs,  algebra_index_vector_type& constrainedIndex, std::vector<algebra_index_vector_type>& vConstrainingIndices)
		{
			for(size_t i = 0; i < vConstrainingIndices.size(); ++i)
			{
				if(vConstrainingIndices[i].size() != constrainedIndex.size())
					{UG_LOG("Wrong number of indices. Cannot split row.\n"); return false;}
			}

			for(size_t i = 0; i < constrainedIndex.size(); ++i)
			{
				typename Vector<T>::entry_type& val = rhs[constrainedIndex[i]];

				// choose randomly the first dof to add whole rhs (must be the same as for row)
				rhs[vConstrainingIndices[0][i]] += val;
				val = 0.0;
			}
			return true;
		}


};


}; // namespace ug



#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__POST_PROCESS__CONSTRAINTS__P1_CONSTRAINTS_POST_PROCESS__ */
