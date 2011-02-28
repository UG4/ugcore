/*
 * p1_constraints_post_process.h
 *
 *  Created on: 01.03.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__POST_PROCESS__CONSTRAINTS__P1_CONSTRAINTS_POST_PROCESS__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__POST_PROCESS__CONSTRAINTS__P1_CONSTRAINTS_POST_PROCESS__

#include "lib_discretization/assemble.h"

namespace ug {

template <	typename TDoFDistribution,
			typename TAlgebra>
class SymP1ConstraintsPostProcess : public IPostProcess<TDoFDistribution, TAlgebra> {
	public:
	// 	DoF Distribution Type
		typedef IDoFDistribution<TDoFDistribution> dof_distribution_type;

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

		virtual IAssembleReturn post_process_jacobian(matrix_type& J,
		                                              const vector_type& u,
		                                              const dof_distribution_type& dofDistr,
		                                              number time = 0.0)
		{
		//  \todo: Implement correctly
		//	dummy for rhs
			vector_type rhsDummy; rhsDummy.resize(u.size());


			return post_process_linear(J, rhsDummy, u, dofDistr, time);
		}

		virtual IAssembleReturn post_process_linear(matrix_type& mat,
		                                            vector_type& rhs,
		                                            const vector_type& u,
		                                            const dof_distribution_type& dofDistr,
		                                            number time = 0.0)
		{
		//	algebra indices of constraining vertex
			std::vector<algebra_index_vector_type> vConstrainingInd;

		//	algebra indices of constrained vertex
			algebra_index_vector_type constrainedInd;

		//	vector of constraining vertices
			std::vector<VertexBase*> vConstrainingVrt;

		//	iterators for hanging vertices
			typename geometry_traits<HangingVertex>::const_iterator iter, iterBegin, iterEnd;

		//	get begin end of hanging vertices
			iterBegin = dofDistr.template begin<HangingVertex>();
			iterEnd = dofDistr.template end<HangingVertex>();

		//	loop constraining edges
			for(iter = iterBegin; iter != iterEnd; ++iter)
			{
			//	resize tmp arrays
				vConstrainingInd.clear();
				vConstrainingVrt.clear();

			//	get hanging vert
				HangingVertex* hgVrt = *iter;

			//	switch constraining parent
				switch(hgVrt->get_parent_base_object_type_id())
				{
				case EDGE:
				{
				//	cast to constraining edge
					ConstrainingEdge* constrainingEdge =
							dynamic_cast<ConstrainingEdge*>(hgVrt->get_parent());

				//	check that edge is correct
					if(constrainingEdge == NULL)
					{
						UG_LOG("ERROR in 'OneSideP1ConstraintsPostProcess::post_process_linear:'"
								" Parent element should be constraining edge, but is not.\n");
						return IAssemble_ERROR;
					}

				//	get constraining vertices
					EdgeBaseIterator edgeIter = constrainingEdge->constrained_edges_begin();
					for(; edgeIter != constrainingEdge->constrained_edges_end(); ++edgeIter)
					{
					//	get constrained edge
						ConstrainedEdge* constrainedEdge = dynamic_cast<ConstrainedEdge*>(*edgeIter);

					//	check
						if(constrainedEdge == NULL)
						{
							UG_LOG("ERROR in 'OneSideP1ConstraintsPostProcess::post_process_linear:'"
									" Child element should be constrained edge, but is not.\n");
							return IAssemble_ERROR;
						}

					//	get non-hanging vertex
						VertexBase* vrt = GetConnectedVertex(constrainedEdge, hgVrt);

					//	push back in list of interpolation vertices
						vConstrainingVrt.push_back(vrt);
					}
				}
					break;
				case FACE:
				{
				//	cast to constraining quadrilateral
					ConstrainingQuadrilateral* bigQuad =
							dynamic_cast<ConstrainingQuadrilateral*>(hgVrt->get_parent());

				//	check that quad is correct
					if(bigQuad == NULL)
					{
						UG_LOG("ERROR in 'OneSideP1ConstraintsPostProcess::post_process_linear:'"
								" Parent element should be constraining quad, but is not.\n");
						return IAssemble_ERROR;
					}

				//	get constraining vertices
				//	\todo: This is only valid for a surface grid!!!
					for(size_t i=0; i < bigQuad->num_vertices(); ++i)
						vConstrainingVrt.push_back(bigQuad->vertex(i));
				}
					break;
				default: UG_LOG("ERROR in 'OneSideP1ConstraintsPostProcess::post_process_linear:'"
								" Parent element of hang. vertex wrong.\n");
						return IAssemble_ERROR;
				}

			//	resize constraining indices
				vConstrainingInd.resize(vConstrainingVrt.size());

			// 	get algebra indices for constraining vertices
				for(size_t i=0; i < vConstrainingVrt.size(); ++i)
					dofDistr.get_inner_algebra_indices(vConstrainingVrt[i], vConstrainingInd[i]);

			// 	get algebra indices constrained vertices
				dofDistr.get_inner_algebra_indices(hgVrt, constrainedInd);

			// 	Split using indices
				if(!SplitAddRow(mat, constrainedInd, vConstrainingInd))
				{
					UG_LOG("ERROR while splitting rows. Aborting.\n");
					return IAssemble_ERROR;
				}

			//	Set interpolation
				if(!SetInterpolation(mat, constrainedInd, vConstrainingInd))
				{
					UG_LOG("ERROR while setting interpolation. Aborting.\n");
					return IAssemble_ERROR;
				}

			//	adapt rhs
				if(!HandleRhs(rhs, constrainedInd, vConstrainingInd))
				{
					UG_LOG("ERROR while setting interpolation. Aborting.\n");
					return IAssemble_ERROR;
				}
			}

		//  we're done
			return IAssemble_OK;
		}

	protected:
		template <typename T>
		bool SplitAddRow(SparseMatrix<T>& A	,
		                 algebra_index_vector_type& constrainedIndex,
		                 std::vector<algebra_index_vector_type>& vConstrainingIndices)
		{
			for(size_t i = 0; i < vConstrainingIndices.size(); ++i)
			{
				if(vConstrainingIndices[i].size() != constrainedIndex.size())
				{
					UG_LOG("Wrong number of indices. Cannot split row.\n");
					return false;
				}
			}

			for(size_t i = 0; i < constrainedIndex.size(); ++i)
			{
				for(typename SparseMatrix<T>::rowIterator conn = A.beginRow(constrainedIndex[i]);
						!conn.isEnd(); ++conn)
				{
					typename SparseMatrix<T>::value_type block = (*conn).dValue;
					const size_t j = (*conn).iIndex;

					A(constrainedIndex[i], j) = 0.0;

					block *= 1./(vConstrainingIndices.size());

					for(size_t k = 0; k < vConstrainingIndices.size(); ++k)
						A(vConstrainingIndices[k][i], j) += block;
				}

				A(constrainedIndex[i], constrainedIndex[i]) = 1.0;

				number frac = -1.0/(vConstrainingIndices.size());
				for(size_t j=0; j < vConstrainingIndices.size();++j)
					A(constrainedIndex[i], vConstrainingIndices[j][i]) = frac;
			}

			return true;
		}

		template <typename T>
		bool SetInterpolation(SparseMatrix<T>& A	,
		                      algebra_index_vector_type& constrainedIndex,
		                      std::vector<algebra_index_vector_type>& vConstrainingIndices)
		{
			for(size_t i = 0; i < vConstrainingIndices.size(); ++i)
			{
				if(vConstrainingIndices[i].size() != constrainedIndex.size())
				{
					UG_LOG("Wrong number of indices. Cannot split row.\n");
					return false;
				}
			}

			for(size_t i = 0; i < constrainedIndex.size(); ++i)
			{
				A(constrainedIndex[i], constrainedIndex[i]) = 1.0;

				number frac = -1.0/(vConstrainingIndices.size());
				for(size_t j=0; j < vConstrainingIndices.size();++j)
					A(constrainedIndex[i], vConstrainingIndices[j][i]) = frac;
			}
			return true;
		}

		template <typename T>
		bool HandleRhs(Vector<T>& rhs,
		               algebra_index_vector_type& constrainedIndex,
		               std::vector<algebra_index_vector_type>& vConstrainingIndices)
		{
			for(size_t i = 0; i < vConstrainingIndices.size(); ++i)
			{
				if(vConstrainingIndices[i].size() != constrainedIndex.size())
				{
					UG_LOG("Wrong number of indices. Cannot split row.\n");
					return false;
				}
			}

			for(size_t i = 0; i < constrainedIndex.size(); ++i)
			{
				typename Vector<T>::value_type& val = rhs[constrainedIndex[i]];
				val *= 1./(vConstrainingIndices.size());

				// split equally on all constraining indices
				for(size_t j=0; j < vConstrainingIndices.size(); ++j)
					rhs[vConstrainingIndices[j][i]] += val;
			}
			return true;
		}


};



template <	typename TDoFDistribution,
			typename TAlgebra>
class OneSideP1ConstraintsPostProcess : public IPostProcess<TDoFDistribution, TAlgebra> {
	public:
	// 	DoF Distribution Type
		typedef IDoFDistribution<TDoFDistribution> dof_distribution_type;

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

		virtual IAssembleReturn post_process_jacobian(matrix_type& J,
		                                              const vector_type& u,
		                                              const dof_distribution_type& dofDistr,
		                                              number time = 0.0)
		{
		//  \todo: Implement correctly
		//	dummy for rhs
			vector_type rhsDummy; rhsDummy.resize(u.size());


			return post_process_linear(J, rhsDummy, u, dofDistr, time);
		}

		virtual IAssembleReturn post_process_linear(matrix_type& mat,
		                                            vector_type& rhs,
		                                            const vector_type& u,
		                                            const dof_distribution_type& dofDistr,
		                                            number time = 0.0)
		{
		//	algebra indices of constraining vertex
			std::vector<algebra_index_vector_type> vConstrainingInd;

		//	algebra indices of constrained vertex
			algebra_index_vector_type constrainedInd;

		//	vector of constraining vertices
			std::vector<VertexBase*> vConstrainingVrt;

		//	iterators for hanging vertices
			typename geometry_traits<HangingVertex>::const_iterator iter, iterBegin, iterEnd;

		//	get begin end of hanging vertices
			iterBegin = dofDistr.template begin<HangingVertex>();
			iterEnd = dofDistr.template end<HangingVertex>();

		//	loop constraining edges
			for(iter = iterBegin; iter != iterEnd; ++iter)
			{
			//	resize tmp arrays
				vConstrainingInd.clear();
				vConstrainingVrt.clear();

			//	get hanging vert
				HangingVertex* hgVrt = *iter;

			//	switch constraining parent
				switch(hgVrt->get_parent_base_object_type_id())
				{
				case EDGE:
				{
				//	cast to constraining edge
					ConstrainingEdge* constrainingEdge =
							dynamic_cast<ConstrainingEdge*>(hgVrt->get_parent());

				//	check that edge is correct
					if(constrainingEdge == NULL)
					{
						UG_LOG("ERROR in 'OneSideP1ConstraintsPostProcess::post_process_linear:'"
								" Parent element should be constraining edge, but is not.\n");
						return IAssemble_ERROR;
					}

				//	get constraining vertices
					EdgeBaseIterator edgeIter = constrainingEdge->constrained_edges_begin();
					for(; edgeIter != constrainingEdge->constrained_edges_end(); ++edgeIter)
					{
					//	get constrained edge
						ConstrainedEdge* constrainedEdge = dynamic_cast<ConstrainedEdge*>(*edgeIter);

					//	check
						if(constrainedEdge == NULL)
						{
							UG_LOG("ERROR in 'OneSideP1ConstraintsPostProcess::post_process_linear:'"
									" Child element should be constrained edge, but is not.\n");
							return IAssemble_ERROR;
						}

					//	get non-hanging vertex
						VertexBase* vrt = GetConnectedVertex(constrainedEdge, hgVrt);

					//	push back in list of interpolation vertices
						vConstrainingVrt.push_back(vrt);
					}
				}
					break;
				case FACE:
				{
				//	cast to constraining quadrilateral
					ConstrainingQuadrilateral* bigQuad =
							dynamic_cast<ConstrainingQuadrilateral*>(hgVrt->get_parent());

				//	check that quad is correct
					if(bigQuad == NULL)
					{
						UG_LOG("ERROR in 'OneSideP1ConstraintsPostProcess::post_process_linear:'"
								" Parent element should be constraining quad, but is not.\n");
						return IAssemble_ERROR;
					}

				//	get constraining vertices
				//	\todo: This is only valid for a surface grid!!!
					for(size_t i=0; i < bigQuad->num_vertices(); ++i)
						vConstrainingVrt.push_back(bigQuad->vertex(i));
				}
					break;
				default: UG_LOG("ERROR in 'OneSideP1ConstraintsPostProcess::post_process_linear:'"
								" Parent element of hang. vertex wrong.\n");
						return IAssemble_ERROR;
				}

			//	resize constraining indices
				vConstrainingInd.resize(vConstrainingVrt.size());

			// 	get algebra indices for constraining vertices
				for(size_t i=0; i < vConstrainingVrt.size(); ++i)
					dofDistr.get_inner_algebra_indices(vConstrainingVrt[i], vConstrainingInd[i]);

			// 	get algebra indices constrained vertices
				dofDistr.get_inner_algebra_indices(hgVrt, constrainedInd);

			// 	Split using indices
				if(!SplitAddRow(mat, constrainedInd, vConstrainingInd))
				{
					UG_LOG("ERROR while splitting rows. Aborting.\n");
					return IAssemble_ERROR;
				}

			//	Set interpolation
				if(!SetInterpolation(mat, constrainedInd, vConstrainingInd))
				{
					UG_LOG("ERROR while setting interpolation. Aborting.\n");
					return IAssemble_ERROR;
				}

			//	adapt rhs
				if(!HandleRhs(rhs, constrainedInd, vConstrainingInd))
				{
					UG_LOG("ERROR while setting interpolation. Aborting.\n");
					return IAssemble_ERROR;
				}
			}

		//  we're done
			return IAssemble_OK;
		}

	protected:
		template <typename T>
		bool SplitAddRow(SparseMatrix<T>& A,
		                 algebra_index_vector_type& constrainedIndex,
		                 std::vector<algebra_index_vector_type>& vConstrainingIndices)
		{
			for(size_t i = 0; i < vConstrainingIndices.size(); ++i)
			{
				if(vConstrainingIndices[i].size() != constrainedIndex.size())
				{
					UG_LOG("Wring number of indices. Cannot split row.\n");
					return false;
				}
			}

			for(size_t i = 0; i < constrainedIndex.size(); ++i)
			{
				for(typename SparseMatrix<T>::rowIterator conn = A.beginRow(constrainedIndex[i]);
						!conn.isEnd(); ++conn)
				{
					typename SparseMatrix<T>::value_type block = (*conn).dValue;
					const size_t j = (*conn).iIndex;

					// choose randomly the first dof to add whole row
					A(vConstrainingIndices[0][i], j) += block;
					A(constrainedIndex[i], j) = 0.0;
				}
			}
			return true;
		}

		template <typename T>
		bool SetInterpolation(SparseMatrix<T>& A,
		                      algebra_index_vector_type& constrainedIndex,
		                      std::vector<algebra_index_vector_type>& vConstrainingIndices)
		{
			for(size_t i = 0; i < vConstrainingIndices.size(); ++i)
			{
				if(vConstrainingIndices[i].size() != constrainedIndex.size())
				{
					UG_LOG("Wring number of indices. Cannot split row.\n");
					return false;
				}
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
		bool HandleRhs(Vector<T>& rhs,
		               algebra_index_vector_type& constrainedIndex,
		               std::vector<algebra_index_vector_type>& vConstrainingIndices)
		{
			for(size_t i = 0; i < vConstrainingIndices.size(); ++i)
			{
				if(vConstrainingIndices[i].size() != constrainedIndex.size())
				{
					UG_LOG("Wrong number of indices. Cannot split row.\n");
					return false;
				}
			}

			for(size_t i = 0; i < constrainedIndex.size(); ++i)
			{
				typename Vector<T>::value_type& val = rhs[constrainedIndex[i]];

				// choose randomly the first dof to add whole rhs (must be the same as for row)
				rhs[vConstrainingIndices[0][i]] += val;
				val = 0.0;
			}
			return true;
		}


};


}; // namespace ug



#endif /* __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__POST_PROCESS__CONSTRAINTS__P1_CONSTRAINTS_POST_PROCESS__ */
