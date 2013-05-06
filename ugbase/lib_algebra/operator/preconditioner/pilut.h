/*
 * \file ilut.h
 *
 * \date 30.04.2013
 * \author Martin Rupp
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PILUT__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PILUT__

#include "common/util/smart_pointer.h"
#include "lib_algebra/operator/interface/preconditioner.h"
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif
#include "progress.h"
#include "common/util/ostream_util.h"
#include "lib_algebra/parallelization/collect_matrix.h"
namespace ug{

template<typename matrix_type>
void SerializeGlobalRow(BinaryBuffer &stream, const matrix_type &mat, size_t rowIndex)
{
	PROFILE_FUNC_GROUP("algebra parallelization");

	size_t numConnections = mat.num_connections(rowIndex);
	Serialize(stream, numConnections);
		for(typename matrix_type::const_row_iterator conn = mat.begin_row(rowIndex);
						conn != mat.end_row(rowIndex); ++conn)
	{
		size_t i =conn.index();
		typename matrix_type::value_type v = conn.value();
		Serialize(stream, i);
		Serialize(stream, v);
	}
}

template<typename matrix_type>
void DeserializeGlobalRow(BinaryBuffer &stream, matrix_type &mat, size_t rowIndex)
{
	PROFILE_FUNC_GROUP("algebra parallelization");

	size_t numConnections = Deserialize<size_t>(stream);
	for(size_t i=0; i<numConnections; i++)
	{
		size_t j = Deserialize<size_t>(stream);
		typename matrix_type::value_type v;
		Deserialize(stream, v);
		mat(rowIndex, j) = v;
	}
}

template <typename TAlgebra>
class PILUTPreconditioner : public IPreconditioner<TAlgebra>
{
	public:
	//	Algebra type
		typedef TAlgebra algebra_type;

	//	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	//	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Matrix Operator type
		typedef typename IPreconditioner<TAlgebra>::matrix_operator_type matrix_operator_type;

	private:
		typedef typename matrix_type::value_type block_type;
		using IPreconditioner<TAlgebra>::debug_writer;
		using IPreconditioner<TAlgebra>::set_debug;

		typedef typename matrix_type::connection connection;
		typedef typename matrix_type::row_iterator row_iterator;
		typedef typename matrix_type::value_type value_type;

	public:
	//	Constructor
		PILUTPreconditioner(double eps=1e-6) :
			m_eps(eps), m_info(false), m_groupSize(20)
		{};

	// 	Clone

		SmartPtr<ILinearIterator<vector_type> > clone()
		{
			SmartPtr<PILUTPreconditioner<algebra_type> > newInst(new PILUTPreconditioner<algebra_type>(m_eps));
			newInst->set_debug(debug_writer());
			newInst->set_damp(this->damping());
			newInst->set_info(m_info);
			newInst->set_group_size(m_groupSize);
			return newInst;
		}

	// 	Destructor
		virtual ~PILUTPreconditioner()
		{
		};

		void set_threshold(number thresh)
		{
			m_eps = thresh;
		}

		void set_group_size(size_t i)
		{
			m_groupSize = i;
		}

	///	sets storage information output to true or false
		void set_info(bool info)
		{
			m_info = info;
		}

	protected:
	//	Name of preconditioner
		virtual const char* name() const {return "PILUTPreconditioner";}

		void get_distribution(std::vector<size_t> &dist, size_t size)
		{
			PROFILE_FUNC_GROUP("pilut algebra");
			srand(0);
			dist.resize(size);
			size_t groupsize = 5;
			size_t j = 0;
			for(size_t i=0; i<size; i++)
			{
				if(i % groupsize == groupsize-1) j = rand()%pcl::GetNumProcesses();
				dist[i] = j;
				//UG_LOG("index " << i << " is on processor " << j << "\n");
			}

		}


		void get_distribution_and_parallel_matrix(std::vector<size_t> &dist, matrix_type &mat, matrix_type &cMat)
		{
			PROFILE_FUNC_GROUP("pilut algebra");
			const pcl::ProcessCommunicator &pc = mat.layouts()->proc_comm();
			pcl::InterfaceCommunicator<IndexLayout> &ic = mat.layouts()->comm();

			int size;
			if(pcl::GetNumProcesses()>1)
			{
				if(pcl::GetProcRank() == 0)
				{

					collect_matrix(mat, cMat, masterLayout, slaveLayout);
					size = cMat.num_rows();

					pc.broadcast(size);
					get_distribution(dist, size);

					std::vector<ug::BinaryBuffer> bb;
					bb.resize(pcl::GetNumProcesses());
					for(size_t i=0; i<size; i++)
						if(dist[i]!=0)
							SerializeGlobalRow(bb[dist[i]], cMat, i);
					for(int i=1; i<pcl::GetNumProcesses(); i++)
						ic.send_raw(i, bb[i].buffer(), bb[i].write_pos(), false);
					ic.communicate();
				}
				else
				{
					collect_matrix(mat, cMat, masterLayout, slaveLayout);

					pc.broadcast(size);
					cMat.resize_and_clear(size, size);
					get_distribution(dist, size);

					ug::BinaryBuffer buf;
					ic.receive_raw(0, buf);
					ic.communicate();

					for(size_t i=0; i<size; i++)
						if(dist[i]==pcl::GetProcRank())
							DeserializeGlobalRow(buf, cMat, i);
				}
				// todo: those parallel dirichlet-point are a problem... ask andreas
				for(size_t i=0; i<cMat.num_rows(); i++)
				{
					if(cMat.is_isolated(i) == false) continue;
					bool bFound;
					row_iterator it = cMat.get_connection(i, i, bFound);
					if(bFound) it.value() = 1.0;
				}
			}
			else
			{
				cMat = mat;
				dist.resize(mat.num_rows(), 0);
			}
			//cMat.print();
		}

		void ClearRow(matrix_type &mat, matrix_type &matT, size_t i)
		{
			PROFILE_FUNC_GROUP("pilut algebra");
			for(row_iterator i_it = mat.begin_row(i); i_it != mat.end_row(i); ++i_it)
			{
				i_it.value() = 0.0;
				matT(i_it.index(), i) = 0.0;
			}
		}

		void send_to_0(matrix_type &mat, size_t k, const std::vector<size_t> &dist)
		{
			PROFILE_FUNC_GROUP("pilut algebra");
			BinaryBuffer buf;
			if(pcl::GetNumProcesses() <= 1) return;

			pcl::InterfaceCommunicator<IndexLayout> &ic = mat.layouts()->comm();
			if(pcl::GetProcRank() == dist[k])
			{
				if(pcl::GetProcRank() == 0) return;
				SerializeGlobalRow(buf, mat, k);
				ic.send_raw(0, buf.buffer(), buf.write_pos(), false);
				ic.communicate();

			}
			else if(pcl::GetProcRank() == 0)
			{
				ic.receive_raw(dist[k], buf);
				ic.communicate();
				DeserializeGlobalRow(buf, mat, k);
			}
		}

		void get_consistent_row(matrix_type &mat, size_t k, const std::vector<size_t> &dist)
		{
			PROFILE_FUNC_GROUP("pilut algebra");
			if(pcl::GetNumProcesses() <= 1) return;

			const pcl::ProcessCommunicator &pc = mat.layouts()->proc_comm();
			BinaryBuffer buf;
			if(pcl::GetProcRank() == dist[k])
				SerializeGlobalRow(buf, mat, k);
			pc.broadcast(buf, dist[k]);
			if(pcl::GetProcRank() != dist[k])
				DeserializeGlobalRow(buf, mat, k);
		}

		vector_type collC;
		vector_type collD;
		void get_vectors(matrix_type &cMat)
		{
			if(pcl::GetProcRank() == 0)
			{
				size_t N = cMat.num_rows();
				collC.resize(N);
				collD.resize(N);
			}
		}

		std::vector<connection> con;

		void eliminate_row_in_window(size_t i, size_t K1, size_t K2, matrix_type &L, matrix_type &U)
		{
			// eliminate all entries A(i, k) with k<i with rows A(k, .) and k<i

			con.clear();

			for(row_iterator it = L.begin_row(i); it != L.end_row(i); ++it)
			{
				if(it.index() >= K1)
					con.push_back(connection(it.index(), it.value()));
			}
			for(row_iterator it =U.begin_row(i); it != U.end_row(i); ++it)
				con.push_back(connection(it.index(), it.value()));

			size_t u_part;
			// eliminate all entries A(i, k) with k<i with rows U(k, .) and k<i
			for(size_t i_it = 0; i_it < con.size(); ++i_it)
			{
				size_t k = con[i_it].iIndex;
				if(k >= i || k >= K2)
				{
					// safe where U begins / L ends in con
					u_part = i_it;
					break;
				}
				if(con[i_it].dValue == 0.0) continue;
				value_type &ukk = U.begin_row(k).value();

				// add row k to row i by A(i, .) -= U(k,.)  A(i,k) / U(k,k)
				// so that A(i,k) is zero.
				// safe A(i,k)/U(k,k) in con, (later L(i,k) )
				value_type &d = con[i_it].dValue = con[i_it].dValue / ukk;

				typename matrix_type::row_iterator k_it = U.begin_row(k); // upper row iterator
				++k_it; // skip diag
				size_t j = i_it+1;
				while(k_it != U.end_row(k) && j < con.size())
				{
					// (since con and U[k] is sorted, we can do sth like a merge on the two lists)
					if(k_it.index() == con[j].iIndex)
					{
						// match
						con[j].dValue -= k_it.value() * d;
						++k_it;	++j;
					}
					else if(k_it.index() < con[j].iIndex)
					{
						// we have a value in U(k, (*k_it).iIndex), but not in A.
						// check tolerance criteria

						typename matrix_type::connection
							c(k_it.index(), k_it.value() * d * -1.0);

						// insert sorted
						con.insert(con.begin()+j, c);
						++j;
						++k_it;
					}
					else
					{
						// we have a value in A(k, con[j].iIndex), but not in U.
						++j;
					}
				}
				// insert new connections after last connection of row i
				if (k_it!=U.end_row(k)){
					for (;k_it!=U.end_row(k);++k_it)
					{
						typename matrix_type::connection c(k_it.index(),-k_it.value()*d);
							con.push_back(c);
					}
				}
			}
			L.set_matrix_row(i, &con[0], u_part);
			U.set_matrix_row(i, &con[u_part], con.size()-u_part);
		}
		void distributed_ILU(matrix_type &mat, const std::vector<size_t> &dist, matrix_type &L, matrix_type &U)
		{
			PROFILE_FUNC_GROUP("pilut algebra");
			//	Prepare Inverse Matrix

			size_t N = mat.num_rows();
			if(pcl::GetProcRank()==0)
			{
				matrix_type A;
							A= mat;

				size_t n = A.num_rows();
				for(size_t k=0; k<n; k++)
				{
					for(size_t i=k+1; i<n; i++)
					{
						if(A.has_connection(i,k) == false) continue;
						A(i,k) = A(i,k)/A(k,k);
						for(size_t j=k+1; j<n; j++)
							A(i,j) = A(i,j) - A(i,k)*A(k,j);
					}
				}
				UG_LOG("normal LU:\n");
				A.print();
				UG_LOG("----------\n");
			}


			size_t totalentries=0;
			size_t maxentries=0;
			bool bUseProgress = true; N > 20000;
			progress p;
			if(bUseProgress)
			{
				UG_LOG("Using PILUT(" << m_eps << ") on " << mat.num_rows() << " x " << mat.num_rows() << " matrix...");
				p.start(N);
			}

			assert(dist[0] == 0);

			matrix_type LT;
			U.resize_and_clear(N, N);
			L.resize_and_clear(N, N);
			LT.resize_and_clear(N, N);

			for(size_t k=0; k<N; k++)
				for(row_iterator it = mat.begin_row(k); it != mat.end_row(k); ++it)
				{
					if(it.index() >= k)
						U(k, it.index()) = it.value();
					else
					{
						if(dist[k] != pcl::GetProcRank()) continue;
						LT(it.index(), k) = it.value(); // (i,j) mit i < j
						L(k, it.index()) = it.value(); // (i,j) mit i < j
					}

				}

			L.print();
			U.print();

#if 0
			PROFILE_BEGIN(dILU1);

			for(size_t k=0; k<N; k++)
			{
				if(bUseProgress) p.set(k);

				get_consistent_row(U, k, dist);
				send_to_0(L, k, dist);
				//if(dabs(A(k,k)) < 1e-10)
					//	return false;

				for(row_iterator i_itT = LT.begin_row(k); i_itT != LT.end_row(k); ++i_itT)
				{
					//A(i,k) = A(i,k)/A(k,k);
					size_t i = i_itT.index();
					if(dist[i] != pcl::GetProcRank()) continue;

					if(L.has_connection(i,k) == false) continue;
					value_type &Aik = L(i,k);

					//for(size_t j=k+1; j<n; j++)
						//A(i,j) = A(i,j) - A(i,k)*A(k,j);
					row_iterator j_it = U.begin_row(k);
					if(j_it == U.end_row(k)) continue;
					Aik /= j_it.value();

					for(++j_it; // j = k+1
							j_it != U.end_row(k); ++j_it)
					{
						size_t j = j_it.index();
						if(i<=j)
						{
							value_type &Uij = U(i,j);
							Uij = Uij - Aik * j_it.value();
						}
						else
						{
							value_type &Lij = L(i,j);
							Lij = Lij - Aik * j_it.value();
							LT(j,i) = 1.0;
						}
					}
				}

			}
#else
			PROFILE_BEGIN(dILU1);




			size_t globalK=0;
			size_t K1, K2;
			std::vector<bool> done;
			done.resize(mat.num_rows(), false);
			std::vector<size_t> vdone;
			while(true)
			{
				K1=globalK;
				for(K2=K1; dist[K2] == dist[K1] && K2<N; K2++)  { }


				if(dist[globalK] != pcl::GetProcRank())
				{
					for(int i=K1; i<K2; i++)
						get_consistent_row(U, i, dist);
				}
				else
				{
					get_consistent_row(U, K1, dist);
					for(size_t i=K1+1; i<K2; i++)
					{
						eliminate_row_in_window(i, K1, K2, L, U);
						get_consistent_row(U, i, dist);
					}
				}


				if(K2 == N) break;

				vdone.clear();
				for(size_t ii=K1; ii<K2; ii++)
				{
					for(row_iterator it_i = LT.get_iterator_or_next(ii, K2); it_i != LT.end_row(ii); ++it_i)
					{
						size_t i=it_i.index();
						if(done[i]) continue;
						done[i] = true;
						vdone.push_back(i);
						eliminate_row_in_window(i, K1, K2, L, U);
					}
				}

				for(size_t i=0; i<vdone.size(); i++)
					done[i] = false;
				globalK = K2;
			}
			//send_to_0(L, k, dist);
#endif
			if(bUseProgress) p.stop();

			p.start(N);
			if(pcl::GetProcRank() == 0)
			{


				/*int totalentries = 0;
				for(size_t i=0; i<N; i++)
				{
					p.set(i);
					for(row_iterator it = mat.begin_row(i); it != mat.end_row(i); ++it)
					{
						size_t j = it.index();
						if(j >= i)
							U(i, j) = it.value();
						else
							L(i, j) = it.value();
						totalentries++;
					}
				}*/
				U.defragment();
				L.defragment();
				UG_LOG(reset_floats << "Total entries: " << totalentries << " (" << ((double)totalentries) / (mat.num_rows()*mat.num_rows()) << "% of dense)");
			}


			p.stop();

			UG_LOG("mat");
			L.print();
			U.print();
		}


	//	Preprocess routine
		virtual bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp)
		{
			PROFILE_FUNC_GROUP("pilut algebra");

			std::vector<size_t> dist;
			matrix_type &mat = *pOp;
			matrix_type cMat;
			get_distribution_and_parallel_matrix(dist, mat, cMat);

			distributed_ILU(cMat, dist, m_L, m_U);
			get_vectors(cMat);

			return true;
		}


		void solve_ILU(vector_type &c, vector_type &d)
		{
			PROFILE_FUNC_GROUP("pilut algebra");
			// L
			for(size_t i=0; i < m_L.num_rows(); i++)
			{
				// c[i] = d[i] - m_L[i]*c;
				c[i] = d[i];
				for(typename matrix_type::row_iterator it = m_L.begin_row(i); it != m_L.end_row(i); ++it)
					MatMultAdd(c[i], 1.0, c[i], -1.0, it.value(), c[it.index()] );
				// lii = 1.0.
			}
			// U
			//
			// last row diagonal U entry might be close to zero with corresponding zero rhs
			// when solving Navier Stokes system, therefore handle separately
			{
				size_t i=m_U.num_rows()-1;
				typename matrix_type::row_iterator it = m_U.begin_row(i);
				UG_ASSERT(it.index() == i, "");
				block_type &uii = it.value();
				typename vector_type::value_type s = c[i];
				// check if close to zero
				if (BlockNorm(uii)<m_small){
					// set correction to zero
					c[i] = 0;
					if (BlockNorm(s)>m_small){
						UG_LOG("Warning: near-zero diagonal entry in last row of U with corresponding non-near-zero rhs entry (" << BlockNorm(s) << ")\n");
					}
				} else {
					// c[i] = s/uii;
					InverseMatMult(c[i], 1.0, uii, s);
				}
			}
			// handle all other rows
			for(size_t i=m_U.num_rows()-2; ; i--)
			{
				typename matrix_type::row_iterator it = m_U.begin_row(i);
				UG_ASSERT(it.index() == i, "");
				block_type &uii = it.value();

				typename vector_type::value_type s = c[i];
				++it; // skip diag
				for(; it != m_U.end_row(i); ++it)
					// s -= it.value() * c[it.index()];
					MatMultAdd(s, 1.0, s, -1.0, it.value(), c[it.index()] );

					// c[i] = s/uii;
					InverseMatMult(c[i], 1.0, uii, s);

					if(i==0) break;
			}

		}

		void pvec(const vector_type &v, const char *name)
		{
			PROFILE_FUNC_GROUP("pilut algebra");
			IF_DEBUG(LIB_ALG_MATRIX, 3)
			{
				for(size_t i=0; i<v.size(); i++){	UG_LOG(name << "[" << i << "] = " << v[i] << "\n"); } UG_LOG("\n");
			}
		}

	//	Stepping routine
		virtual bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d)
		{
			PROFILE_FUNC_GROUP("pilut algebra");
			matrix_type &mat = *pOp;

			pvec(d, "d");

			if(pcl::GetProcRank() == 0)
			{
				collD.set(0.0);
				for(size_t i=0; i<d.size(); i++)
					collD[i] = d[i];

				pvec(collD, "collD pre");
			}

			// send d -> collD
			ComPol_VecAdd<vector_type > compolAdd(&collD, &d);
			pcl::InterfaceCommunicator<IndexLayout> &ic = mat.layouts()->comm();
			ic.receive_data(masterLayout, compolAdd);
			ic.send_data(slaveLayout, compolAdd);
			ic.communicate();

			if(pcl::GetProcRank() == 0)
			{
				pvec(collD, "collD p1");
				solve_ILU(collC, collD);
				pvec(collD, "collD p2");
				pvec(collC, "collC");
			}
			// send collC -> c
			ComPol_VecCopy<vector_type> compolCopy(&c, &collC);
			ic.receive_data(slaveLayout, compolCopy);
			ic.send_data(masterLayout, compolCopy);
			ic.communicate();

			pvec(c, "c");
			if(pcl::GetProcRank() == 0)
			{
				for(size_t i=0; i<c.size(); i++)
					c[i] = collC[i];
				pvec(c, "c");
			}
			c.set_storage_type(PST_CONSISTENT);
			return true;
		}

	//	Postprocess routine
		virtual bool postprocess() {return true;}

	protected:
		matrix_type m_L;
		matrix_type m_U;
		double m_eps;
		bool m_info;
		static const number m_small;
		size_t m_groupSize;
		IndexLayout masterLayout, slaveLayout;
};

// define constant
template <typename TAlgebra>
const number PILUTPreconditioner<TAlgebra>::m_small = 1e-8;

} // end namespace ug

#endif
