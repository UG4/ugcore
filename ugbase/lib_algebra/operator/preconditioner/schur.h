/*
 * schubt
 r.h
 *
 *  Created on: 18.12.2013
 *      Author: anaegel
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__SCHUR__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__SCHUR__


#ifdef UG_PARALLEL

#include <iostream>
#include <sstream>
#include <string>
#include <set>

#include "lib_algebra/operator/interface/preconditioner.h"
#include "lib_algebra/operator/interface/linear_operator.h"
#include "lib_algebra/operator/interface/linear_operator_inverse.h"
#include "lib_algebra/operator/interface/matrix_operator.h"
#include "lib_algebra/operator/interface/matrix_operator_inverse.h"
#include "lib_algebra/parallelization/parallelization.h"
#include "lib_algebra/operator/debug_writer.h"
#include "pcl/pcl.h"



namespace ug{


class SlicingOp
{

};

class SlicingData{

public:
	enum slice_desc_type {SD_INNER=0, SD_SKELETON, SLICE_DESC_SIZE};
	typedef std::vector<slice_desc_type> slice_desc_type_vector;
	//typedef std::set<int> slice_desc_set;
	typedef std::vector<int> slice_desc_set;

protected:
	slice_desc_type_vector m_slice_types;
	slice_desc_set m_slice_set[(int) SLICE_DESC_SIZE]; //!< mapping islice -> iglobal


public:
	/// constructor
	SlicingData(const slice_desc_type_vector &types)
	: m_slice_types(types)
	{
		auto_fill_sets();
	}

	/// copy types
	void set_types(const slice_desc_type_vector &types)
	{ m_slice_types = types; }

	/// auto fill for sets
	/// assigns every i=0.. m_slice_types.size()-1 to exactly one set
	void auto_fill_sets()
	{

		const size_t ntypes = m_slice_types.size();


		for (size_t i=0; i< ntypes; ++i)
		{
			slice_desc_type tt = get_type(i);
			slice_desc_set &set =slice(tt);

			// if sd is ordered, then this is constant
			//set.insert(set.end(), i);
			set.push_back(i);
		}



		UG_LOG("SlicingData::auto_fill_sets:" << ntypes << " "<< slice(SD_INNER).size() << " "<< slice(SD_SKELETON).size() << std::endl);

		slice_desc_set::const_iterator it;

		{
			UG_LOG("Skeleton:")
			const slice_desc_set &myset=slice(SD_SKELETON);

				for (it=myset.begin(); it!=myset.end(); ++it)   UG_LOG(*it << " ");
		}

		{

		UG_LOG("\nInner:")
	    const slice_desc_set &myset=slice(SD_INNER);
		for (it=myset.begin(); it!=myset.end(); ++it)
		    UG_LOG(*it << " ");
		}

	}


	template<class VT>
	SmartPtr<VT> slice_clone_without_values(const VT &full_src, slice_desc_type type) const
	{
		const slice_desc_set &slice_desc = slice(type);
		// SmartPtr<VT> slice_clone = full_src.clone_without_values();
		//slice_clone->resize(slice_desc.size());

		SmartPtr<VT> slice_clone = new VT(slice_desc.size());

		SmartPtr<AlgebraLayouts> slice_layouts = new AlgebraLayouts(*full_src.layouts());
		replace_indices_in_layout(type, slice_layouts->master());
		replace_indices_in_layout(type, slice_layouts->slave());
		slice_clone->set_layouts(slice_layouts);

		//UG_LOG(*slice_layouts);

		return slice_clone;

	}

    /// copy: slice of vector -> small vector
	template<class VT>
	void get_vector_slice(const VT &full_src, slice_desc_type desc, VT &small_dst) const
	{
		const slice_desc_set &slice_desc = slice(desc);
		small_dst.resize(slice_desc.size());
		slice_desc_set::const_iterator elem = slice_desc.begin();
		for (size_t i=0; i<slice_desc.size(); ++i, ++elem)
			small_dst[i] = full_src[*elem];
	}



	/// copy: small vector -> slice of a vector
	template<class VT>
	void set_vector_slice(const VT &small_src, VT &full_dst, slice_desc_type desc) const
	{
		const slice_desc_set &slice_desc = slice(desc);
		slice_desc_set::const_iterator elem = slice_desc.begin();
		for (size_t i=0; i<slice_desc.size(); ++i, ++elem)
				full_dst[*elem] = small_src[i];
	}

	 /// copy: slice of vector -> small vector
		template<class VT>
		void subtract_vector_slice(const VT &full_src, slice_desc_type desc, VT &small_dst) const
		{
				const slice_desc_set &slice_desc = slice(desc);
				small_dst.resize(slice_desc.size());
				slice_desc_set::const_iterator elem = slice_desc.begin();
				for (size_t i=0; i<slice_desc.size(); ++i, ++elem)
					small_dst[i] -= full_src[*elem];
		}

		/// copy: small vector -> slice of a vector
			template<class VT>
			void subtract_vector_slice(const VT &small_src, VT &full_dst, slice_desc_type desc) const
			{
				const slice_desc_set &slice_desc = slice(desc);
				slice_desc_set::const_iterator elem = slice_desc.begin();
				for (size_t i=0; i<slice_desc.size(); ++i, ++elem)
						full_dst[*elem] -= small_src[i];
			}


	// Extracts a slice from a (full) matrix
	template<class MT>
	void get_matrix(const MT &A, slice_desc_type row_type, slice_desc_type col_type, MT &Aslice) const
	{
		const slice_desc_set &row_slice = slice(row_type);
		const slice_desc_set &col_slice = slice(col_type);
		UG_LOG("SlicingData::get_matrix:" << row_slice.size() << "x" << col_slice.size()<< std::endl)
		Aslice.resize_and_clear(row_slice.size(), col_slice.size());

		int ii=0;
		for (slice_desc_set::const_iterator elem = row_slice.begin();
			elem!=row_slice.end(); ++elem, ++ii)
		{
			const int i = *elem; // global index
			for(typename MT::const_row_iterator it = A.begin_row(i);
				it != A.end_row(i); ++it)

			{
				const int j=it.index();
				int jj;
				// if (get_type(j)!=col_type) continue;
				if (find_index(col_type, j, jj))
					Aslice(ii, jj) = it.value();
			 }
		}

	}

	size_t get_num_elems(slice_desc_type type) const
	{return slice(type).size();}

protected:

	/// returns type for a global index
	slice_desc_type get_type(int index)
	{return m_slice_types[index];}


	/// returns the set of global indices for a given type
	const slice_desc_set &slice(slice_desc_type type) const
	{return m_slice_set[type];}

	slice_desc_set &slice(slice_desc_type type)
	{return m_slice_set[type];}


	/// returns local index for a global index
	bool find_index(slice_desc_type type, int gindex, int &index) const
	{
		bool found=false;

		const slice_desc_set &myset=slice(type);
		index = myset.size();

		//slice_desc_set::const_iterator it = myset.find(gindex);
		slice_desc_set::const_iterator it = lower_bound(myset.begin(), myset.end(), gindex);
		if (it != myset.end()) {
			//index =// *it;
			index = std::distance(myset.begin(), it);
			found = true;
		}
		if (found && index >=myset.size()) {
		UG_ASSERT( (!found || index<myset.size()) , "Invalid index found!");
		}
		return found;
	}

	void replace_indices_in_layout(slice_desc_type type, IndexLayout &il) const
	{
		// UG_LOG("ilpre:\n"<< il);
		IndexLayout::iterator iter;
		for (iter = il.begin(); iter!=il.end(); ++iter)
		{
			// iterate over interfaces
			// i.e., pcl::OrderedInterface<size_t, std::vector>
			IndexLayout::Interface &interf=il.interface(iter);

			IndexLayout::Interface::iterator eiter;
			for (eiter = interf.begin(); eiter!=interf.end(); ++eiter)
			{
				// replace elements in interface
				size_t &elem = interf.get_element(eiter);
				int newind;

				bool found=find_index(type, elem, newind);
				UG_ASSERT(found, "SlicingData:: Did not find index???");
				//UG_LOG(eiter->elem <<  ", " << eiter->localID <<  "=>");
				interf.get_element(eiter) = newind;
				//UG_LOG(eiter->elem <<  ", " << eiter->localID <<  std::endl);
				//UG_LOG("->" << interf.get_element(eiter) << std::endl);
			}
		}
		// UG_LOG("ilpost:\n"<<il);

	}

};





template <typename TAlgebra>
class SchurComplementOperator
	: public ILinearOperator<	typename TAlgebra::vector_type,
	  	  	  	  	  	  	  	typename TAlgebra::vector_type>,
	  public DebugWritingObject<TAlgebra>
{
	public:

	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	// 	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	protected:
		using DebugWritingObject<TAlgebra>::write_debug;
		using DebugWritingObject<TAlgebra>::debug_writer;

	public:
	///	constructor
	SchurComplementOperator(SmartPtr<MatrixOperator<matrix_type, vector_type> > Alocal,
							SlicingData::slice_desc_type_vector &sdv)
	: m_spOperator(Alocal),
	  m_slicing(sdv)
	{
		m_op[0][0] = new MatrixOperator<matrix_type, vector_type>();
		m_op[0][1] = new MatrixOperator<matrix_type, vector_type>();
		m_op[1][0] = new MatrixOperator<matrix_type, vector_type>();
		m_op[1][1] = new MatrixOperator<matrix_type, vector_type>();
	}

	// destructor
	virtual ~SchurComplementOperator() {};

	///	name of solver
	virtual const char* name() const {return "My local Schur complement Solver";}


	/// implementation of the operator for the solution dependent initialization.
	void init(const vector_type& u) {init();}

	///	initializes the solver for operator A
	virtual void init();

	///	applies the Schur complement built from matrix operator set via 'set_matrix()'
	/// to 'u' and returns the result 'f := S times u'
	virtual void apply(vector_type& f, const vector_type& u);

	///	solves the system
	virtual void apply_sub(vector_type& f, const vector_type& u);



	//	save current operator
	void set_matrix(SmartPtr<MatrixOperator<matrix_type, vector_type> > A)
	{ m_spOperator = A; }

	///	sets a Dirichlet solver
	void set_dirichlet_solver(SmartPtr<ILinearOperatorInverse<vector_type> > dirichletSolver)
	{ m_spDirichletSolver = dirichletSolver; }


	matrix_type &sub_matrix(int r, int c)
	{return sub_operator(r,c)->get_matrix();}

	SmartPtr<MatrixOperator<matrix_type, vector_type> > sub_operator(int r, int c)
	{return m_op[r][c];}

	size_t sub_size(SlicingData::slice_desc_type type)
	{return m_slicing.get_num_elems(type);}

	const SlicingData &slicing() const
	{return m_slicing;}

protected:
	// 	Operator that is inverted by this Inverse Operator
	SmartPtr<MatrixOperator<matrix_type,vector_type> > m_spOperator;

	// slices from matrix
	const SlicingData m_slicing;

	// 	Linear Solver to invert the local Dirichlet problems
	SmartPtr<ILinearOperatorInverse<vector_type> > m_spDirichletSolver;



	// sub matrices/operator (will be set by init)
	SmartPtr<MatrixOperator<matrix_type,vector_type> > m_op[2][2];

	int m_applyCnt;

}; /* end class 'LocalSchurComplement' */



/// operator implementation of the FETI-DP solver
/**
 * This operator implements a Schur complement solver */
//template <typename TAlgebra>
//class SchurSolver : public IMatrixOperatorInverse<	typename TAlgebra::matrix_type,
//													typename TAlgebra::vector_type>,
//	public DebugWritingObject<TAlgebra>

template <typename TAlgebra>
class SchurPrecond: public IPreconditioner<TAlgebra>
{
	public:
	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	// 	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Base type
		typedef IPreconditioner<TAlgebra> base_type;

	protected:
		//using base_type::set_debug;
		using base_type::write_debug;
		using base_type::debug_writer;

	public:
	///	constructor
		SchurPrecond();

		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			SmartPtr<SchurPrecond<algebra_type> > newInst(new SchurPrecond<algebra_type>());
			UG_ASSERT(0, "Implement SchurPrecond::clone()!")
			return newInst;
		}

	protected:
	///	name of solver
		virtual const char* name() const {return "Schur complement";}

		//	Preprocess routine
		virtual bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp);

		//	Stepping routine
		virtual bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d);

		//	Postprocess routine
		virtual bool postprocess();//  {return true;}

	///	returns if parallel solving is supported
		virtual bool supports_parallel() const
		{
			if(m_spDirichletSolver.valid()
				&& (!m_spDirichletSolver->supports_parallel()))
					return false;

			if(m_spSkeletonSolver.valid()
					&& (!m_spSkeletonSolver->supports_parallel()))
					return false;

			return true;
		}

public:
		void set_schur_complement_operator(SmartPtr<SchurComplementOperator<algebra_type> > scop)
		{ m_spSchurComplementOp = scop; }

		//	define an approximation for schur complement
		void set_schur_complement_approx(SmartPtr<MatrixOperator<matrix_type, vector_type> > S)
		{ m_spSkeletonMatrix = S; }

	///	sets the Dirichlet solver (forward to Schur complement)
		void set_dirichlet_solver(SmartPtr<ILinearOperatorInverse<vector_type> > dirichletSolver)
		{ m_spDirichletSolver = dirichletSolver; }

	///	sets the coarse problem solver
		void set_skeleton_solver(SmartPtr<ILinearOperatorInverse<vector_type> > skeletonSolver)
		{ m_spSkeletonSolver = skeletonSolver; }

	//	set debug output
		void set_debug(SmartPtr<IDebugWriter<algebra_type> > spDebugWriter)
		{
			base_type::set_debug(spDebugWriter);
			//m_spSchurComplementOp->set_debug(spDebugWriter);
		}

	/*///	initializes the solver for operator A
		virtual bool init(SmartPtr<MatrixOperator<matrix_type, vector_type> > A);

	///	solves the reduced system \f$F \lambda = d\f$ with preconditioned cg method
		virtual bool apply_return_defect(vector_type& x, vector_type& d);

	///	solves the system
		virtual bool apply(vector_type& x, const vector_type& b)
		{
		//	copy defect
			vector_type d; d.resize(b.size());
			d = b;

		//	solve on copy of defect
			return apply_return_defect(x, d);
		}*/

		// destructor
//		virtual ~SchurPrecond() {};


			/*void set_domain_decomp_info(pcl::IDomainDecompositionInfo& ddInfo)
			{
				m_pDDInfo = &ddInfo;
			}*/




	protected:
	// 	Reference to operator that is inverted by this Inverse Operator
	//	SmartPtr<MatrixOperator<matrix_type,vector_type> > m_spOperator;

	//	Local Schur complement for each subdomain
		SmartPtr<SchurComplementOperator<algebra_type> > m_spSchurComplementOp;

		// approximation of Schur complement
		SmartPtr<MatrixOperator<matrix_type,vector_type> > m_spSkeletonMatrix;


	//	Dirichlet solver for inverse of \f$A_{II}\f$ in local Schur complement
		SmartPtr<ILinearOperatorInverse<vector_type> > m_spDirichletSolver;

	//	Solver used in solving coarse problem on root.
	// 	It solves \f$S_{\Pi \Pi} u_{\Pi} = \tilde{f}_{\Pi}\f$ 
		SmartPtr<ILinearOperatorInverse<vector_type> > m_spSkeletonSolver;

		// temporary vectors for correction/defect
		SmartPtr<vector_type> m_aux_rhs[2];
		SmartPtr<vector_type> m_aux_sol[2];

	//	pointer to Domain decomposition info object
	//	pcl::IDomainDecompositionInfo* m_pDDInfo;


		int m_iterCnt;

}; /* end class 'SchurPrecond' */

} // end namespace ug

#endif /* UG_PARALLEL */

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__FETI__ */
