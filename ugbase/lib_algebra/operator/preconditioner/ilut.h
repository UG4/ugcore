
#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__ILUT__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__ILUT__

#include "common/util/smart_pointer.h"
#include "lib_algebra/operator/interface/preconditioner.h"
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif
#include "common/progress.h"
#include "common/util/ostream_util.h"
#include "common/profiler/profiler.h"

#include "lib_algebra/algebra_common/vector_util.h"

#include "lib_algebra/ordering_strategies/algorithms/IOrderingAlgorithm.h"
#include "lib_algebra/ordering_strategies/algorithms/native_cuthill_mckee.h" // for backward compatibility

#include "lib_algebra/algebra_common/permutation_util.h"

namespace ug{

template <typename TAlgebra>
class ILUTPreconditioner : public IPreconditioner<TAlgebra>
{
	public:
	//	Algebra type
		typedef TAlgebra algebra_type;

	//	Vector type
		typedef typename TAlgebra::vector_type vector_type;
		typedef typename vector_type::value_type vector_value;

	//	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

		typedef typename matrix_type::row_iterator matrix_row_iterator;
		typedef typename matrix_type::const_row_iterator const_matrix_row_iterator;
		typedef typename matrix_type::connection matrix_connection;

	///	Ordering type
		typedef std::vector<size_t> ordering_container_type;
		typedef IOrderingAlgorithm<TAlgebra, ordering_container_type> ordering_algo_type;

	///	Matrix Operator type
		typedef typename IPreconditioner<TAlgebra>::matrix_operator_type matrix_operator_type;
		using IPreconditioner<TAlgebra>::set_debug;

	protected:
		typedef typename matrix_type::value_type block_type;

		using IPreconditioner<TAlgebra>::debug_writer;
		using IPreconditioner<TAlgebra>::write_debug;

	private:
		typedef IPreconditioner<TAlgebra> base_type;

	public:
	///	Constructor
		ILUTPreconditioner(double eps=1e-6)
			: m_eps(eps), m_info(false), m_show_progress(true), m_bSortIsIdentity(false)
		{
			//default was set true
			m_spOrderingAlgo = make_sp(new NativeCuthillMcKeeOrdering<TAlgebra, ordering_container_type>());
		};

	/// clone constructor
		ILUTPreconditioner(const ILUTPreconditioner<TAlgebra> &parent)
			: base_type(parent), m_spOrderingAlgo(parent.m_spOrderingAlgo)
		{
			m_eps = parent.m_eps;
			set_info(parent.m_info);
			m_bSortIsIdentity = parent.m_bSortIsIdentity;
		}

	///	Clone
		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			return make_sp(new ILUTPreconditioner<algebra_type>(*this));
		}


	///	Destructor
		virtual ~ILUTPreconditioner()
		{};

	///	returns if parallel solving is supported
		virtual bool supports_parallel() const {return true;}

	///	sets threshold for incomplete LU factorisation (added 01122010ih)
		void set_threshold(number thresh)
		{
			m_eps = thresh;
		}
		
	///	sets storage information output to true or false
		void set_info(bool info)
		{
			m_info = info;
		}
		
	///	set whether the progress meter should be shown
		void set_show_progress(bool s)
		{
			m_show_progress = s;
		}

		virtual std::string config_string() const
		{
			std::stringstream ss;
			ss << "ILUT(threshold = " << m_eps << ", sort = " << (m_spOrderingAlgo.valid()?"true":"false") << ")";
			if(m_eps == 0.0) ss << " = Sparse LU";
			return ss.str();
		}

	/// 	sets an ordering algorithm
		void set_ordering_algorithm(SmartPtr<ordering_algo_type> ordering_algo){
			m_spOrderingAlgo = ordering_algo;
		}

	/// set cuthill-mckee sort on/off
		void set_sort(bool b)
		{
			if(b){
				m_spOrderingAlgo = make_sp(new NativeCuthillMcKeeOrdering<TAlgebra, ordering_container_type>());
			}
			else{
				m_spOrderingAlgo = SPNULL;
			}

			UG_LOG("\nILUT: please use 'set_ordering_algorithm(..)' in the future\n");
		}


	protected:
	//	Name of preconditioner
		virtual const char* name() const {return "ILUT";}

	protected:
		virtual bool init(SmartPtr<ILinearOperator<vector_type> > J,
		                  const vector_type& u)
		{
		//	cast to matrix based operator
			SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp =
					J.template cast_dynamic<MatrixOperator<matrix_type, vector_type> >();

		//	Check that matrix if of correct type
			if(pOp.invalid())
				UG_THROW(name() << "::init': Passed Operator is "
						"not based on matrix. This Preconditioner can only "
						"handle matrix-based operators.");

			m_u = &u;

		//	forward request to matrix based implementation
			return base_type::init(pOp);
		}

		bool init(SmartPtr<ILinearOperator<vector_type> > L)
		{
		//	cast to matrix based operator
			SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp =
					L.template cast_dynamic<MatrixOperator<matrix_type, vector_type> >();

		//	Check that matrix if of correct type
			if(pOp.invalid())
				UG_THROW(name() << "::init': Passed Operator is "
						"not based on matrix. This Preconditioner can only "
						"handle matrix-based operators.");

			m_u = NULL;

		//	forward request to matrix based implementation
			return base_type::init(pOp);
		}

		bool init(SmartPtr<MatrixOperator<matrix_type, vector_type> > Op)
		{
			m_u = NULL;

			return base_type::init(Op);
		}

	//	Preprocess routine
		virtual bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp)
		{
			return preprocess_mat(*pOp);
		}

	public:
		virtual bool preprocess_mat(matrix_type &mat)
		{
#ifdef 	UG_PARALLEL
			matrix_type m2;
			m2 = mat;

			MatAddSlaveRowsToMasterRowOverlap0(m2);

		//	set zero on slaves
			std::vector<IndexLayout::Element> vIndex;
			CollectUniqueElements(vIndex,  mat.layouts()->slave());
			SetDirichletRow(m2, vIndex);

			// Even after this setting of Dirichlet rows, it is possible that there are
			// zero rows on a proc because of the distribution:
			// For example, if one has a horizontal grid interface between two SHADOW_RIM_COPY
			// vertices, but the shadowing element for the hSlave side is vMaster (without being in
			// any horizontal interface). Then, as the horizontal interface (on the shadowed level)
			// is not part of the algebraic layouts, the hSlave is not converted into a Dirichlet row
			// by the previous commands.
			// As a first aid, we will simply convert any zero row on the current proc into a
			// Dirichlet row.
			// TODO: The corresponding rhs vector entry could be non-zero!
			// It is definitely not set to zero by change_storage_type(PST_UNIQUE) as the index is not contained
			// in the vector layouts either. Still, the defect assembling process might contain a vertex
			// loop and assemble something that is not solution-dependent! What do we do then?
			size_t nInd = m2.num_rows();
			size_t cnt = 0;
			for (size_t i = 0; i < nInd; ++i)
			{
				if (!m2.num_connections(i))
				{
					m2(i,i) = 1.0;
					++cnt;
				}
			}
#ifndef NDEBUG
			if (cnt) UG_LOG_ALL_PROCS("Converted "<<cnt<<" zero rows into Dirichlet rows.\n");
#endif

			return preprocess_mat2(m2);
#else
			return preprocess_mat2(mat);
#endif
		}

		virtual bool preprocess_mat2(matrix_type &mat)
		{
			PROFILE_BEGIN_GROUP(ILUT_preprocess, "ilut algebra");
			//matrix_type &mat = *pOp;
			STATIC_ASSERT(matrix_type::rows_sorted, Matrix_has_to_have_sorted_rows);
			write_debug(mat, "ILUT_PreprocessIn");

			matrix_type* A;
			matrix_type permA;

			if(m_spOrderingAlgo.valid() && !m_useOverlap){
				if(m_u){
					m_spOrderingAlgo->init(&m_ILU, *m_u);
				}
				else{
					m_spOrderingAlgo->init(&m_ILU);
				}
			}

			if(m_spOrderingAlgo.valid())
			{
				m_spOrderingAlgo->compute();
				m_ordering = m_spOrderingAlgo->ordering();

				m_bSortIsIdentity = GetInversePermutation(m_ordering, m_old_ordering);

				if(!m_bSortIsIdentity){
					SetMatrixAsPermutation(permA, mat, m_ordering);
					A = &permA;
				}
				else{
					A = &mat;
				}
			}
			else
			{
				A = &mat;
			}

			m_L.resize_and_clear(A->num_rows(), A->num_cols());
			m_U.resize_and_clear(A->num_rows(), A->num_cols());

			size_t totalentries=0;
			size_t maxentries=0;

			if((A->num_rows() > 0) && (A->num_cols() > 0)){
				// con is the current line of L/U
				// i also tried using std::list here or a special custom vector-based linked list
				// but vector is fastest, even with the insert operation.
				std::vector<matrix_connection> con;
				con.reserve(300);
				con.resize(0);

				// init row 0 of U
				for(matrix_row_iterator i_it = A->begin_row(0); i_it != A->end_row(0); ++i_it)
					con.push_back(matrix_connection(i_it.index(), i_it.value()));
				m_U.set_matrix_row(0, &con[0], con.size());

				Progress prog;
				if(m_show_progress)
					PROGRESS_START_WITH(prog, A->num_rows(),
						"Using ILUT(" << m_eps << ") on " << A->num_rows() << " x " << A->num_rows() << " matrix...");
				
				for(size_t i=1; i<A->num_rows(); i++)
				{
					if(m_show_progress) {PROGRESS_UPDATE(prog, i);}
					con.resize(0);
					size_t u_part=0;

					UG_COND_THROW(A->num_connections(i) == 0, "row " << i << " has no connections");

					// get the row A(i, .) into con
					double dmax=0;
					for(matrix_row_iterator i_it = A->begin_row(i); i_it != A->end_row(i); ++i_it)
					{
						con.push_back(matrix_connection(i_it.index(), i_it.value()));
						if(dmax < BlockNorm(i_it.value()))
							dmax = BlockNorm(i_it.value());
					}

					// eliminate all entries A(i, k) with k<i with rows U(k, .) and k<i
					for(size_t i_it = 0; i_it < con.size(); ++i_it)
					{
						size_t k = con[i_it].iIndex;
						if(k >= i)
						{
							// safe where U begins / L ends in con
							u_part = i_it;
							break;
						}
						if(con[i_it].dValue == 0.0) continue;
						UG_COND_THROW(!(m_U.num_connections(k) != 0 && m_U.begin_row(k).index() == k), "");
						block_type &ukk = m_U.begin_row(k).value();

						// add row k to row i by A(i, .) -= U(k,.)  A(i,k) / U(k,k)
						// so that A(i,k) is zero.
						// safe A(i,k)/U(k,k) in con, (later L(i,k) )
						con[i_it].dValue = con[i_it].dValue / ukk;
						block_type d = con[i_it].dValue;
						UG_COND_THROW(!BlockMatrixFiniteAndNotTooBig(d, 1e40), "i = " << i << " " << d);

						matrix_row_iterator k_it = m_U.begin_row(k); // upper row iterator
						++k_it; // skip diag
						size_t j = i_it+1;
						while(k_it != m_U.end_row(k) && j < con.size())
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

								matrix_connection c(k_it.index(), k_it.value() * d * -1.0);

								UG_COND_THROW(!BlockMatrixFiniteAndNotTooBig(c.dValue, 1e40), "i = " << i << " " << c.dValue);
								if(BlockNorm(c.dValue) > dmax * m_eps)
								{
									// insert sorted
									con.insert(con.begin()+j, c);
									++j; // don't do this when using iterators
								}
								// else do some lumping
								++k_it;
							}
							else
							{
								// we have a value in A(k, con[j].iIndex), but not in U.
								++j;
							}
						}
						// insert new connections after last connection of row i
						if (k_it!=m_U.end_row(k)){
							for (;k_it!=m_U.end_row(k);++k_it){
								matrix_connection c(k_it.index(),-k_it.value()*d);
						        if(BlockNorm(c.dValue) > dmax * m_eps)
						        {
						        	con.push_back(c);
						        }
						    }
						};
					}


					totalentries+=con.size();
					if(maxentries < con.size()) maxentries = con.size();

					// safe L and U
					m_L.set_matrix_row(i, &con[0], u_part);
					m_U.set_matrix_row(i, &con[u_part], con.size()-u_part);

				}
				if(m_show_progress) {PROGRESS_FINISH(prog);}

				m_L.defragment();
				m_U.defragment();
			}

			if (m_info==true)
			{
				m_L.print("L");
				m_U.print("U");
				UG_LOG("\n	ILUT storage information:\n");
				UG_LOG("	A is " << A->num_rows() << " x " << A->num_cols() << " matrix.\n");
				UG_LOG("	A nr of connections: " << A->total_num_connections()  << "\n");
				UG_LOG("	L+U nr of connections: " << m_L.total_num_connections()+m_U.total_num_connections() << "\n");
				UG_LOG("	Increase factor: " << (float)(m_L.total_num_connections() + m_U.total_num_connections() )/A->total_num_connections() << "\n");
				UG_LOG(reset_floats << "	Total entries: " << totalentries << " (" << ((double)totalentries) / (A->num_rows()*A->num_rows()) << "% of dense)\n");
				if(m_spOrderingAlgo.valid())
				{
					UG_LOG("	Using " <<  m_spOrderingAlgo->name() << "\n");
				}
				else
				{
					UG_LOG("Not using sort.");
				}
			}

			return true;
		}

	//	Stepping routine
		virtual bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d)
		{
#ifdef UG_PARALLEL
			SmartPtr<vector_type> spDtmp = d.clone();
			spDtmp->change_storage_type(PST_UNIQUE);
			bool b = solve(c, *spDtmp);

			c.set_storage_type(PST_ADDITIVE);
			c.change_storage_type(PST_CONSISTENT);
			return b;
#else
			return solve(c, d);
#endif
			return true;
		}

		virtual bool solve(vector_type& c, const vector_type& d)
		{
			if(m_spOrderingAlgo.invalid()|| m_bSortIsIdentity)
			{
				if(applyLU(c, d) == false) return false;
			}
			else
			{
				PROFILE_BEGIN_GROUP(ILUT_StepWithReorder, "ilut algebra");

				SetVectorAsPermutation(c, d, m_ordering);
				c2.resize(c.size());
				if(applyLU(c2, c) == false) return false;
				SetVectorAsPermutation(c, c2, m_old_ordering);
			}
			return true;
		}


		virtual bool applyLU(vector_type& c, const vector_type& d)
		{
			PROFILE_BEGIN_GROUP(ILUT_step, "ilut algebra");
			// apply iterator: c = LU^{-1}*d (damp is not used)
			// L
			for(size_t i=0; i < m_L.num_rows(); i++)
			{
				// c[i] = d[i] - m_L[i]*c;
				c[i] = d[i];
				for(matrix_row_iterator it = m_L.begin_row(i); it != m_L.end_row(i); ++it)
					MatMultAdd(c[i], 1.0, c[i], -1.0, it.value(), c[it.index()] );
				// lii = 1.0.
			}

			// U
			//
			// last row diagonal U entry might be close to zero with corresponding zero rhs 
			// when solving Navier Stokes system, therefore handle separately
			if(m_U.num_rows() > 0)
			{
				size_t i=m_U.num_rows()-1;
				matrix_row_iterator it = m_U.begin_row(i);
				UG_ASSERT(it != m_U.end_row(i), i);
				UG_ASSERT(it.index() == i, i);
				block_type &uii = it.value();
				vector_value s = c[i];
				// check if diag part is significantly smaller than rhs
				// This may happen when matrix is indefinite with one eigenvalue
				// zero. In that case, the factorization on the last row is
				// nearly zero due to round-off errors. In order to allow ill-
				// scaled matrices (i.e. small matrix entries row-wise) this
				// is compared to the rhs, that is small in this case as well.
				if (false && BlockNorm(uii) <= m_small * m_eps * BlockNorm(s)){
					UG_LOG("ILUT("<<m_eps<<") Warning: Near-zero diagonal entry "
							"with norm "<<BlockNorm(uii)<<" in last row of U "
							" with corresponding non-near-zero rhs with norm "
							<< BlockNorm(s) << ". Setting rhs to zero.\n");
					// set correction to zero
					c[i] = 0;
				} else {
					// c[i] = s/uii;
					InverseMatMult(c[i], 1.0, uii, s);
				}
			}

			// handle all other rows
			if(m_U.num_rows() > 1){
				for(size_t i=m_U.num_rows()-2; ; i--)
				{
					matrix_row_iterator it = m_U.begin_row(i);
					UG_ASSERT(it != m_U.end_row(i), i);
					UG_ASSERT(it.index() == i, i);
					block_type &uii = it.value();

					vector_value s = c[i];
					++it; // skip diag
					for(; it != m_U.end_row(i); ++it){
						// s -= it.value() * c[it.index()];
						MatMultAdd(s, 1.0, s, -1.0, it.value(), c[it.index()] );
					}
					// c[i] = s/uii;
					InverseMatMult(c[i], 1.0, uii, s);

					if(i==0) break;
				}
			}
			return true;
		}

		virtual bool multi_apply(std::vector<vector_type> &vc, const std::vector<vector_type> &vd)
		{
			PROFILE_BEGIN_GROUP(ILUT_step, "ilut algebra");
			// apply iterator: c = LU^{-1}*d (damp is not used)
			// L
			for(size_t i=0; i < m_L.num_rows(); i++)
			{

				for(size_t e=0; e<vc.size(); e++)
				{
					vector_type &c = vc[e];
					const vector_type &d = vd[e];
					// c[i] = d[i] - m_L[i]*c;
					c[i] = d[i];
					for(matrix_row_iterator it = m_L.begin_row(i); it != m_L.end_row(i); ++it)
						MatMultAdd(c[i], 1.0, c[i], -1.0, it.value(), c[it.index()] );
					// lii = 1.0.
				}
			}

			// U
			//
			// last row diagonal U entry might be close to zero with corresponding zero rhs
			// when solving Navier Stokes system, therefore handle separately
			if(m_U.num_rows() > 0)
			{
				for(size_t e=0; e<vc.size(); e++)
				{
					vector_type &c = vc[e];
					size_t i=m_U.num_rows()-1;
					matrix_row_iterator it = m_U.begin_row(i);
					UG_ASSERT(it != m_U.end_row(i), i);
					UG_ASSERT(it.index() == i, i);
					block_type &uii = it.value();
					vector_value s = c[i];
					// check if diag part is significantly smaller than rhs
					// This may happen when matrix is indefinite with one eigenvalue
					// zero. In that case, the factorization on the last row is
					// nearly zero due to round-off errors. In order to allow ill-
					// scaled matrices (i.e. small matrix entries row-wise) this
					// is compared to the rhs, that is small in this case as well.
					if (false && BlockNorm(uii) <= m_small * m_eps * BlockNorm(s)){
						UG_LOG("ILUT("<<m_eps<<") Warning: Near-zero diagonal entry "
								"with norm "<<BlockNorm(uii)<<" in last row of U "
								" with corresponding non-near-zero rhs with norm "
								<< BlockNorm(s) << ". Setting rhs to zero.\n");
						// set correction to zero
						c[i] = 0;
					} else {
						// c[i] = s/uii;
						InverseMatMult(c[i], 1.0, uii, s);
					}
				}
			}

			// handle all other rows
			if(m_U.num_rows() > 1){
				for(size_t i=m_U.num_rows()-2; ; i--)
				{
					for(size_t e=0; e<vc.size(); e++)
					{
						vector_type &c = vc[e];
						matrix_row_iterator it = m_U.begin_row(i);
						UG_ASSERT(it != m_U.end_row(i), i);
						UG_ASSERT(it.index() == i, i);
						block_type &uii = it.value();

						vector_value s = c[i];
						++it; // skip diag
						for(; it != m_U.end_row(i); ++it){
							// s -= it.value() * c[it.index()];
							MatMultAdd(s, 1.0, s, -1.0, it.value(), c[it.index()] );
						}
						// c[i] = s/uii;
						InverseMatMult(c[i], 1.0, uii, s);

						if(i==0) break;
					}
				}
			}
			return true;
		}

	//	Postprocess routine
		virtual bool postprocess() {return true;}

	protected:
		vector_type c2;
		matrix_type m_L;
		matrix_type m_U;
		double m_eps;
		bool m_info;
		bool m_show_progress;
		static const number m_small;
		std::vector<size_t> newIndex, oldIndex;

	/// for ordering algorithms
		SmartPtr<ordering_algo_type> m_spOrderingAlgo;
		ordering_container_type m_ordering, m_old_ordering;

		bool m_bSortIsIdentity;

		const vector_type* m_u;
};

// define constant
template <typename TAlgebra>
const number ILUTPreconditioner<TAlgebra>::m_small = 1e-8;

} // end namespace ug

#endif
