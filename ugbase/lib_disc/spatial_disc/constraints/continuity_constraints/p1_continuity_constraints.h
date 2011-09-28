/*
 * p1_continuity_constraints.h
 *
 *  Created on: 01.03.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__CONTINUITY_CONSTRAINTS__P1_CONTINUITY_CONSTRAINTS__
#define __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__CONTINUITY_CONSTRAINTS__P1_CONTINUITY_CONSTRAINTS__

#include "lib_disc/assemble_interface.h"
#include "lib_grid/algorithms/geom_obj_util/vertex_util.h"

namespace ug {

template <	typename TDoFDistribution,
			typename TAlgebra>
class SymP1ConstraintsPostProcess : public IConstraint<TDoFDistribution, TAlgebra> {
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
		virtual int type() {return CT_CONSTRAINTS;}

		virtual bool adjust_defect(vector_type& d, const vector_type& u,
		                           const dof_distribution_type& dofDistr,
		                           number time = 0.0)
		{
		//	not implemented
			return false;
		}

		virtual bool adjust_rhs(vector_type& rhs, const vector_type& u,
		                        const dof_distribution_type& dofDistr,
		                        number time = 0.0)
		{
		//	not implemented
			return false;
		}

		virtual bool adjust_jacobian(matrix_type& J,
		                                              const vector_type& u,
		                                              const dof_distribution_type& dofDistr,
		                                              number time = 0.0)
		{
		//  \todo: Implement correctly
		//	dummy for rhs
			vector_type rhsDummy; rhsDummy.resize(u.size());


			return adjust_linear(J, rhsDummy, u, dofDistr, time);
		}

		virtual bool adjust_linear(matrix_type& mat,
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

		//	loop subsets
			for(int si = 0; si < dofDistr.num_subsets(); ++si)
			{
			//	get begin end of hanging vertices
				iterBegin = dofDistr.template begin<HangingVertex>(si);
				iterEnd = dofDistr.template end<HangingVertex>(si);

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
							UG_LOG("ERROR in 'OneSideP1ConstraintsPostProcess::adjust_linear:'"
									" Parent element should be constraining edge, but is not.\n");
							return false;
						}

					//	get constraining vertices
						for(size_t i_cde = 0; i_cde < constrainingEdge->num_constrained_edges(); ++i_cde)
						{
						//	get constrained edge
							ConstrainedEdge* constrainedEdge = dynamic_cast<ConstrainedEdge*>(
																constrainingEdge->constrained_edge(i_cde));

						//	check
							if(constrainedEdge == NULL)
							{
								UG_LOG("ERROR in 'OneSideP1ConstraintsPostProcess::adjust_linear:'"
										" Child element should be constrained edge, but is not.\n");
								return false;
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
							UG_LOG("ERROR in 'OneSideP1ConstraintsPostProcess::adjust_linear:'"
									" Parent element should be constraining quad, but is not.\n");
							return false;
						}

					//	get constraining vertices
					//	\todo: This is only valid for a surface grid!!!
						for(size_t i_cf=0; i_cf < bigQuad->num_constrained_faces(); ++i_cf)
						{
							Face* face = bigQuad->constrained_face(i_cf);

							VertexBase* vrt = NULL;
							size_t i_vrt = 0;
							for(i_vrt = 0; i_vrt < face->num_vertices(); ++i_vrt)
							{
								vrt = face->vertex(i_vrt);
								if(hgVrt != vrt && dynamic_cast<HangingVertex*>(vrt) == NULL)
									break;
							}
							if(i_vrt == face->num_vertices()) {UG_LOG("ERROR: Vertex not detected.\n"); return false;}

							vConstrainingVrt.push_back(vrt);
						}
					}
						break;
					default: UG_LOG("ERROR in 'OneSideP1ConstraintsPostProcess::adjust_linear:'"
									" Parent element of hang. vertex wrong.\n");
							return false;
					}

				//	resize constraining indices
					vConstrainingInd.resize(vConstrainingVrt.size());

				// 	get algebra indices for constraining vertices
					for(size_t i=0; i < vConstrainingVrt.size(); ++i)
						dofDistr.inner_algebra_indices(vConstrainingVrt[i], vConstrainingInd[i]);

				// 	get algebra indices constrained vertices
					dofDistr.inner_algebra_indices(hgVrt, constrainedInd);

				// 	Split using indices
					if(!SplitAddRow(mat, constrainedInd, vConstrainingInd))
					{
						UG_LOG("ERROR while splitting rows. Aborting.\n");
						return false;
					}

				//	adapt rhs
					if(!HandleRhs(rhs, constrainedInd, vConstrainingInd))
					{
						UG_LOG("ERROR while setting interpolation. Aborting.\n");
						return false;
					}
				}

			//	second loop to set the constraints
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
							UG_LOG("ERROR in 'OneSideP1ConstraintsPostProcess::adjust_linear:'"
									" Parent element should be constraining edge, but is not.\n");
							return false;
						}

					//	get constraining vertices
						for(size_t i_cde = 0; i_cde < constrainingEdge->num_constrained_edges(); ++i_cde)
						{
						//	get constrained edge
							ConstrainedEdge* constrainedEdge = dynamic_cast<ConstrainedEdge*>(
																constrainingEdge->constrained_edge(i_cde));

						//	check
							if(constrainedEdge == NULL)
							{
								UG_LOG("ERROR in 'OneSideP1ConstraintsPostProcess::adjust_linear:'"
										" Child element should be constrained edge, but is not.\n");
								return false;
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
							UG_LOG("ERROR in 'OneSideP1ConstraintsPostProcess::adjust_linear:'"
									" Parent element should be constraining quad, but is not.\n");
							return false;
						}

					//	get constraining vertices
					//	\todo: This is only valid for a surface grid!!!
					//	since then the indices in shadowing vertex and shadow are
					//	the same. In general, on a level, we will get the wrong
					//	vertex from the coarser level
						for(size_t i_cf=0; i_cf < bigQuad->num_constrained_faces(); ++i_cf)
						{
							Face* face = bigQuad->constrained_face(i_cf);

							VertexBase* vrt = NULL;
							size_t i_vrt = 0;
							for(i_vrt = 0; i_vrt < face->num_vertices(); ++i_vrt)
							{
								vrt = face->vertex(i_vrt);
								if(hgVrt != vrt && dynamic_cast<HangingVertex*>(vrt) == NULL)
									break;
							}
							if(i_vrt == face->num_vertices()) {UG_LOG("ERROR: Vertex not detected.\n"); return false;}

							vConstrainingVrt.push_back(vrt);
						}
					}
						break;
					default: UG_LOG("ERROR in 'OneSideP1ConstraintsPostProcess::adjust_linear:'"
									" Parent element of hang. vertex wrong.\n");
							return false;
					}

				//	resize constraining indices
					vConstrainingInd.resize(vConstrainingVrt.size());

				// 	get algebra indices for constraining vertices
					for(size_t i=0; i < vConstrainingVrt.size(); ++i)
						dofDistr.inner_algebra_indices(vConstrainingVrt[i], vConstrainingInd[i]);

				// 	get algebra indices constrained vertices
					dofDistr.inner_algebra_indices(hgVrt, constrainedInd);

				//	Set interpolation
					if(!SetInterpolation(mat, constrainedInd, vConstrainingInd))
					{
						UG_LOG("ERROR while setting interpolation. Aborting.\n");
						return false;
					}
				}
			}

		//  we're done
			return true;
		}

		virtual bool adjust_solution(vector_type& u,
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

		//	loop subsets
			for(int si = 0; si < dofDistr.num_subsets(); ++si)
			{
			//	get begin end of hanging vertices
				iterBegin = dofDistr.template begin<HangingVertex>(si);
				iterEnd = dofDistr.template end<HangingVertex>(si);

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
							UG_LOG("ERROR in 'OneSideP1ConstraintsPostProcess::adjust_linear:'"
									" Parent element should be constraining edge, but is not.\n");
							return false;
						}

					//	get constraining vertices
						for(size_t i_cde = 0; i_cde < constrainingEdge->num_constrained_edges(); ++i_cde)
						{
						//	get constrained edge
							ConstrainedEdge* constrainedEdge = dynamic_cast<ConstrainedEdge*>(
																constrainingEdge->constrained_edge(i_cde));

						//	check
							if(constrainedEdge == NULL)
							{
								UG_LOG("ERROR in 'OneSideP1ConstraintsPostProcess::adjust_linear:'"
										" Child element should be constrained edge, but is not.\n");
								return false;
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
							UG_LOG("ERROR in 'OneSideP1ConstraintsPostProcess::adjust_linear:'"
									" Parent element should be constraining quad, but is not.\n");
							return false;
						}

					//	get constraining vertices
					//	\todo: This is only valid for a surface grid!!!
						for(size_t i=0; i < bigQuad->num_vertices(); ++i)
							vConstrainingVrt.push_back(bigQuad->vertex(i));
					}
						break;
					default: UG_LOG("ERROR in 'OneSideP1ConstraintsPostProcess::adjust_linear:'"
									" Parent element of hang. vertex wrong.\n");
							return false;
					}

				//	resize constraining indices
					vConstrainingInd.resize(vConstrainingVrt.size());

				// 	get algebra indices for constraining vertices
					for(size_t i=0; i < vConstrainingVrt.size(); ++i)
						dofDistr.inner_algebra_indices(vConstrainingVrt[i], vConstrainingInd[i]);

				// 	get algebra indices constrained vertices
					dofDistr.inner_algebra_indices(hgVrt, constrainedInd);

				// 	Interpolate values
					if(!InterpolateValues(u, constrainedInd, vConstrainingInd))
					{
						UG_LOG("ERROR while interpolating values rows. Aborting.\n");
						return false;
					}
				}
			}

		//  we're done
			return true;
		}

	protected:
		bool SplitAddRow(matrix_type& A	,
		                 algebra_index_vector_type& constrainedIndex,
		                 std::vector<algebra_index_vector_type>& vConstrainingIndices)
		{
		//	check number of indices passed
			for(size_t i = 0; i < vConstrainingIndices.size(); ++i)
			{
				if(vConstrainingIndices[i].size() != constrainedIndex.size())
				{
					UG_LOG("Wrong number of indices. Cannot split row.\n");
					return false;
				}
			}

		//	handle each contrained index
			for(size_t i = 0; i < constrainedIndex.size(); ++i)
			{
			//	add coupling constrained dof -> constrained dof
			//	get entry
				typename matrix_type::value_type& block
						= A(constrainedIndex[i], constrainedIndex[i]);

			//	scale by weight
				block *= (1./(vConstrainingIndices.size()))*
							(1./(vConstrainingIndices.size()));

			//	add coupling
				for(size_t k = 0; k < vConstrainingIndices.size(); ++k)
					for(size_t m = 0; m < vConstrainingIndices.size(); ++m)
				{
					A(vConstrainingIndices[k][i],
					  vConstrainingIndices[m][i]) += block;
				}

			//	reset block
				 A(constrainedIndex[i], constrainedIndex[i]) = 0.0;

			//	loop coupling between constrained dof -> constraining dof
				for(typename matrix_type::row_iterator conn = A.begin_row(constrainedIndex[i]);
						conn != A.end_row(constrainedIndex[i]); ++conn)
				{
				//	skip self-coupling (already handled)
					const size_t j = conn.index();
					if(j == constrainedIndex[i]) continue;

				//	get coupling entry
					typename matrix_type::value_type block = conn.value();

				//	get transposed coupling entry
					typename matrix_type::value_type blockT = A(j, constrainedIndex[i]);

				//	multiply the cpl value by the inverse number of constraining
				//	indices
					block *= 1./(vConstrainingIndices.size());
					blockT *= 1./(vConstrainingIndices.size());

				//	add the coupling to the constraining indices rows
					for(size_t k = 0; k < vConstrainingIndices.size(); ++k)
					{
						A(vConstrainingIndices[k][i], j) += block;

						A(j, vConstrainingIndices[k][i]) += blockT;
					}

				//	set the splitted coupling to zero
					conn.value() = 0.0;
					A(j, constrainedIndex[i]) = 0.0;
				}
			}

		//	we're done
			return true;
		}

		bool SetInterpolation(matrix_type& A	,
		                      algebra_index_vector_type& constrainedIndex,
		                      std::vector<algebra_index_vector_type>& vConstrainingIndices)
		{
		//	check number of indices passed
			for(size_t i = 0; i < vConstrainingIndices.size(); ++i)
			{
				if(vConstrainingIndices[i].size() != constrainedIndex.size())
				{
					UG_LOG("Wrong number of indices. Cannot split row.\n");
					return false;
				}
			}

		//	loop all constrained dofs
			for(size_t i = 0; i < constrainedIndex.size(); ++i)
			{
			//	remove all couplings
				for(typename matrix_type::row_iterator conn = A.begin_row(constrainedIndex[i]);
						conn != A.end_row(constrainedIndex[i]); ++conn)
				{
					conn.value() = 0.0;
				}

			//	set diag of row to identity
				A(constrainedIndex[i], constrainedIndex[i]) = 1.0;

			//	set coupling to all contraining dofs the inverse of the
			//	number of contraining dofs
				number frac = -1.0/(vConstrainingIndices.size());
				for(size_t j=0; j < vConstrainingIndices.size();++j)
					A(constrainedIndex[i], vConstrainingIndices[j][i]) = frac;
			}

		//	we're done
			return true;
		}

		bool HandleRhs(vector_type& rhs,
		               algebra_index_vector_type& constrainedIndex,
		               std::vector<algebra_index_vector_type>& vConstrainingIndices)
		{
		//	check number of indices passed
			for(size_t i = 0; i < vConstrainingIndices.size(); ++i)
			{
				if(vConstrainingIndices[i].size() != constrainedIndex.size())
				{
					UG_LOG("Wrong number of indices. Cannot split row.\n");
					return false;
				}
			}

		//	loop constrained indices
			for(size_t i = 0; i < constrainedIndex.size(); ++i)
			{
			//	get constrained rhs
				typename vector_type::value_type& val = rhs[constrainedIndex[i]];
				val *= 1./(vConstrainingIndices.size());

			// 	split equally on all constraining indices
				for(size_t j=0; j < vConstrainingIndices.size(); ++j)
					rhs[vConstrainingIndices[j][i]] += val;

			//	set rhs to zero for contrained index
				val = 0.0;
			}
			return true;
		}

		bool InterpolateValues(vector_type& u,
		                       algebra_index_vector_type& constrainedIndex,
		                       std::vector<algebra_index_vector_type>& vConstrainingIndices)
		{
		//	check number of indices passed
			for(size_t i = 0; i < vConstrainingIndices.size(); ++i)
			{
				if(vConstrainingIndices[i].size() != constrainedIndex.size())
				{
					UG_LOG("Wrong number of indices. Cannot split row.\n");
					return false;
				}
			}

		//	loop constrained indices
			for(size_t i = 0; i < constrainedIndex.size(); ++i)
			{
			//	get constrained rhs
				typename vector_type::value_type& val = u[constrainedIndex[i]];
				const number scale = 1./(vConstrainingIndices.size());

				val = 0.0;

			// 	split equally on all constraining indices
				for(size_t j=0; j < vConstrainingIndices.size(); ++j)
				{
					typename vector_type::value_type entry = u[vConstrainingIndices[j][i]];
					entry *= scale;
					val += entry;
				}
			}
			return true;
		}

};



template <	typename TDoFDistribution,
			typename TAlgebra>
class OneSideP1ConstraintsPostProcess : public IConstraint<TDoFDistribution, TAlgebra> {
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
		virtual int type() {return CT_CONSTRAINTS;}

		virtual bool adjust_jacobian(matrix_type& J,
		                                              const vector_type& u,
		                                              const dof_distribution_type& dofDistr,
		                                              number time = 0.0)
		{
		//  \todo: Implement correctly
		//	dummy for rhs
			vector_type rhsDummy; rhsDummy.resize(u.size());


			return adjust_linear(J, rhsDummy, u, dofDistr, time);
		}

		virtual bool adjust_defect(vector_type& d, const vector_type& u,
		                           const dof_distribution_type& dofDistr,
		                           number time = 0.0)
		{
		//	not implemented
			return false;
		}

		virtual bool adjust_rhs(vector_type& rhs, const vector_type& u,
		                        const dof_distribution_type& dofDistr,
		                        number time = 0.0)
		{
		//	not implemented
			return false;
		}

		virtual bool adjust_solution(vector_type& u,
		                                   const dof_distribution_type& dofDistr,
		                                   number time = 0.0)
		{
		//	not implemented
			return false;
		}

		virtual bool adjust_linear(matrix_type& mat,
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

		//	loop subsets
			for(int si = 0; si < dofDistr.num_subsets(); ++si)
			{

			//	get begin end of hanging vertices
				iterBegin = dofDistr.template begin<HangingVertex>(si);
				iterEnd = dofDistr.template end<HangingVertex>(si);

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
							UG_LOG("ERROR in 'OneSideP1ConstraintsPostProcess::adjust_linear:'"
									" Parent element should be constraining edge, but is not.\n");
							return false;
						}

					//	get constraining vertices
						for(size_t i_cde = 0; i_cde != constrainingEdge->num_constrained_edges(); ++i_cde)
						{
						//	get constrained edge
							ConstrainedEdge* constrainedEdge = dynamic_cast<ConstrainedEdge*>(
																	constrainingEdge->constrained_edge(i_cde));

						//	check
							if(constrainedEdge == NULL)
							{
								UG_LOG("ERROR in 'OneSideP1ConstraintsPostProcess::adjust_linear:'"
										" Child element should be constrained edge, but is not.\n");
								return false;
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
							UG_LOG("ERROR in 'OneSideP1ConstraintsPostProcess::adjust_linear:'"
									" Parent element should be constraining quad, but is not.\n");
							return false;
						}

					//	get constraining vertices
					//	\todo: This is only valid for a surface grid!!!
						for(size_t i=0; i < bigQuad->num_vertices(); ++i)
							vConstrainingVrt.push_back(bigQuad->vertex(i));
					}
						break;
					default: UG_LOG("ERROR in 'OneSideP1ConstraintsPostProcess::adjust_linear:'"
									" Parent element of hang. vertex wrong.\n");
							return false;
					}

				//	resize constraining indices
					vConstrainingInd.resize(vConstrainingVrt.size());

				// 	get algebra indices for constraining vertices
					for(size_t i=0; i < vConstrainingVrt.size(); ++i)
						dofDistr.inner_algebra_indices(vConstrainingVrt[i], vConstrainingInd[i]);

				// 	get algebra indices constrained vertices
					dofDistr.inner_algebra_indices(hgVrt, constrainedInd);

				// 	Split using indices
					if(!SplitAddRow(mat, constrainedInd, vConstrainingInd))
					{
						UG_LOG("ERROR while splitting rows. Aborting.\n");
						return false;
					}

				//	Set interpolation
					if(!SetInterpolation(mat, constrainedInd, vConstrainingInd))
					{
						UG_LOG("ERROR while setting interpolation. Aborting.\n");
						return false;
					}

				//	adapt rhs
					if(!HandleRhs(rhs, constrainedInd, vConstrainingInd))
					{
						UG_LOG("ERROR while setting interpolation. Aborting.\n");
						return false;
					}
				}
			}
		//  we're done
			return true;
		}

	protected:
		bool SplitAddRow(matrix_type& A,
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
				for(typename matrix_type::row_iterator conn = A.begin_row(constrainedIndex[i]);
						conn != A.end_row(constrainedIndex[i]); ++conn)
				{
					typename matrix_type::value_type block = conn.value();
					const size_t j = conn.index();

					// choose randomly the first dof to add whole row
					A(vConstrainingIndices[0][i], j) += block;
					A(constrainedIndex[i], j) = 0.0;
				}
			}
			return true;
		}

		bool SetInterpolation(matrix_type& A,
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

		bool HandleRhs(vector_type& rhs,
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
				typename vector_type::value_type& val = rhs[constrainedIndex[i]];

				// choose randomly the first dof to add whole rhs (must be the same as for row)
				rhs[vConstrainingIndices[0][i]] += val;
				val = 0.0;
			}
			return true;
		}


};


}; // namespace ug



#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__CONTINUITY_CONSTRAINTS__P1_CONTINUITY_CONSTRAINTS__ */
