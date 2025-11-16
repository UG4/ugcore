/*
 * Copyright (c) 2018:  G-CSC, Goethe University Frankfurt
 * Authors: Arne Naegel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

/*
 *  Remark: Based on implementation of element/component Seidel
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__SEQUENTIAL_SSC__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__SEQUENTIAL_SSC__

#include "lib_algebra/operator/interface/preconditioner.h"

#include <vector>
#include <algorithm>

#ifdef UG_PARALLEL
	#include "pcl/pcl_util.h"
	#include "lib_algebra/parallelization/parallelization_util.h"
	#include "lib_algebra/parallelization/parallelization.h"
	#include "lib_algebra/parallelization/matrix_overlap.h"
	#include "lib_algebra/parallelization/parallel_matrix_overlap_impl.h"
#endif

namespace ug{


//! Abstract definition for subspace V_k
template <typename TDomain, typename TAlgebra, typename TObject>
class ILocalSubspace
{
public:
	///	Algebra type
		using algebra_type = TAlgebra;

	///	Vector type
		using vector_type = typename TAlgebra::vector_type;

	///	Matrix type
		using matrix_type = typename TAlgebra::matrix_type;

	/// Grid function type
		using grid_function_type = GridFunction<TDomain, TAlgebra>;
		using vector_type_local = DenseVector< VariableArray1<number> >;
		using matrix_type_local = DenseMatrix< VariableArray2<number> >;

		/// virtual DTOR
		virtual ~ILocalSubspace() = default;

		/// Called once.
		virtual bool preprocess(const vector_type &c)
		{ return true; }


		/// Extract local data (based on group obj).
		virtual void init(TObject*, const vector_type &) = 0;

		/// Extract matrix (on local index set).
		virtual void extract_matrix(const matrix_type &A) = 0;

		/// Extract rhs (on local index set) for parallel subspace correction.
		virtual void extract_rhs(const vector_type &d) = 0;

		/// Extract rhs (on local index set) for sequential subspace correction.
		virtual void extract_rhs(const vector_type &d, const matrix_type &A, const vector_type &c) = 0;

		/// u = u + omega*ck.
		virtual void update_solution(vector_type &u, double omega=1.0) = 0;

		virtual size_t size() { return 0; }
};



//! Concrete definition of subspace V_k (based on size_t)
template <typename TDomain, typename TAlgebra, typename TObject>
class LocalIndexSubspace : public ILocalSubspace<TDomain, TAlgebra, TObject>
{
public:
	///	Algebra type
	using algebra_type = TAlgebra;

	///	Vector type
	using vector_type = typename TAlgebra::vector_type;

	///	Matrix type
	using matrix_type = typename TAlgebra::matrix_type;

	/// Local data
	using index_vector = std::vector<size_t>;
	using vector_type_local = DenseVector< VariableArray1<number> >;
	using matrix_type_local = DenseMatrix< VariableArray2<number> >;

	/// virtual DTOR
	virtual ~LocalIndexSubspace() = default;

	/// Called once.
	/*virtual bool preprocess(const vector_type &c, obj_iterator_type &objBegin, obj_iterator_type &objEnd)
	{
		objBegin = c.template begin<TGroupObj>();
		objEnd = c.template end<TGroupObj>();
		return true;
	}
*/
	/// Extract local data (based on group obj)
	virtual void init(TObject*, const vector_type &) = 0;

	/// Extract matrix (on local index set)
	virtual void extract_matrix(const matrix_type &A)
	{
		using const_row_iterator = typename TAlgebra::matrix_type::const_row_iterator;
		static constexpr int blockSize = TAlgebra::blockSize;
		const size_t numIndex  = this->size();

		// fill local block matrix
		bool bFound;
		m_Aloc.resize(numIndex, numIndex);
		m_Aloc = 0.0;
		for (size_t j = 0; j<numIndex; j++)
		{
			for (size_t k=0;k<numIndex;k++)
			{
				const_row_iterator it = A.get_connection(m_vInd[j],m_vInd[k], bFound);
				if(bFound){
					m_Aloc.subassign(j*blockSize,k*blockSize,it.value());
				}
			}
		}
	}


	/// Extract rhs (on local index set) for parallel subspace correction
	virtual void extract_rhs(const vector_type &d)
	{
		const size_t numIndex  = this->size();
		static constexpr int blockSize = TAlgebra::blockSize;

		// compute s[j] := d[j] - sum_k A(j,k)*c[k]
		// note: the loop over k is the whole matrix row (not only selected indices)
		m_dloc.resize(numIndex);
		for (size_t j = 0; j<numIndex; j++)
		{
			typename TAlgebra::vector_type::value_type sumj = d[m_vInd[j]];
			m_dloc.subassign(j*blockSize,sumj);
		}
	}

	/// Extract rhs (on local index set) for sequential subspace correction
	virtual void extract_rhs(const vector_type &d, const matrix_type &A, const vector_type &c)
	{
		using const_row_iterator = typename TAlgebra::matrix_type::const_row_iterator;
		static constexpr int blockSize = TAlgebra::blockSize;
		const size_t numIndex  = this->size();

		// compute s[j] := d[j] - sum_k A(j,k)*c[k]
		// note: the loop over k is the whole matrix row (not only selected indices)
		m_dloc.resize(numIndex);
		for (size_t j = 0; j<numIndex; j++)
		{
			typename TAlgebra::vector_type::value_type sumj = d[m_vInd[j]];
			for(const_row_iterator it = A.begin_row(m_vInd[j]); it != A.end_row(m_vInd[j]) ; ++it)
			{
				MatMultAdd(sumj, 1.0, sumj, -1.0, it.value(), c[it.index()]);
			};

			m_dloc.subassign(j*blockSize,sumj);

		};
	};

	/// u = u + ck
	virtual void update_solution(vector_type &u, double omega=1.0)
	{
		const size_t numIndex  = this->size();

		// solve block
		m_uloc.resize(numIndex);
		InverseMatMult(m_uloc, 1.0, m_Aloc, m_dloc);
		for (size_t j=0;j<numIndex;j++)
		{
			u[m_vInd[j]] += omega*m_uloc[j];
		}

	}


	virtual size_t size() { return m_vInd.size(); }

protected:

	/// Algebraic indices.
	index_vector m_vInd;

	/// Memory for local algebra
	vector_type_local m_uloc;
	vector_type_local m_dloc;
	matrix_type_local m_Aloc;

};


//! Abstract definition for subspace V_k (based on DoFIndex)
template <typename TDomain, typename TAlgebra, typename TObject>
class LocalDoFSubspace : public ILocalSubspace<TDomain, TAlgebra, TObject>
{
public:
	///	Algebra type
	using algebra_type = TAlgebra;

	///	Vector type
	using vector_type = typename TAlgebra::vector_type;

	///	Matrix type
	using matrix_type = typename TAlgebra::matrix_type;

	/// Local data
	using index_vector = std::vector<DoFIndex>;
	using vector_type_local = DenseVector< VariableArray1<number> >;
	using matrix_type_local = DenseMatrix< VariableArray2<number> >;

	/// virtual DTOR
	virtual ~LocalDoFSubspace() = default;

	bool check(void *obj) const
	{ return (dynamic_cast<TObject*> (obj) != nullptr); }

	/// Extract local data (based on group obj)
	virtual void init(TObject*, const vector_type &) = 0;

	/// Extract matrix (on local index set)
	virtual void extract_matrix(const matrix_type &A)
	{
		// using const_row_iterator = typename TAlgebra::matrix_type::const_row_iterator;
		// static constexpr int blockSize = TAlgebra::blockSize;
		const size_t numIndex  = this->size();

		// fill local block matrix
		m_Aloc.resize(numIndex, numIndex);
		m_Aloc = 0.0;
		for (size_t j = 0; j<numIndex; ++j)
		{
			for (size_t k=0; k<numIndex; ++k)
			{
				m_Aloc.entry(j,k) = DoFRef(A,m_vInd[j],m_vInd[k]);
			}
		}

	}


	/// Extract rhs (on local index set) for parallel subspace correction
	virtual void extract_rhs(const vector_type &d)
	{
		const size_t numIndex  = this->size();
		m_dloc.resize(numIndex);

		for (size_t j = 0; j<numIndex; j++)
		{
			m_dloc[j] = DoFRef(d, m_vInd[j]);
		}
	}

	/// Extract rhs (on local index set) for sequential subspace correction
	virtual void extract_rhs(const vector_type &d, const matrix_type &A, const vector_type &c)
	{
		using const_row_iterator = typename TAlgebra::matrix_type::const_row_iterator;
		static constexpr int blockSize = TAlgebra::blockSize;
		const size_t numIndex  = this->size();

		// compute s[j] := d[j] - sum_k A(j,k)*c[k]
		// note: the loop over k is the whole matrix row (not only selected indices)
		// (code taken from component Gauss-Seidel)
		m_dloc.resize(numIndex);
		for (size_t j = 0; j<numIndex; j++)
		{
			typename TAlgebra::vector_type::value_type sumdj = d[m_vInd[j][0]];
			for(const_row_iterator it = A.begin_row(m_vInd[j][0]); it != A.end_row(m_vInd[j][0]) ; ++it)
			{ MatMultAdd(sumdj, 1.0, sumdj, -1.0, it.value(), c[it.index()]); }

			// TODO: CHECK!!!
			m_dloc.subassign(j*blockSize, sumdj);
		};
	};



	/// u = u + ck
	virtual void update_solution(vector_type &u, double omega=1.0)
	{
		const size_t numIndex  = this->size();

		// solve block
		m_uloc.resize(numIndex);
		InverseMatMult(m_uloc, 1.0, m_Aloc, m_dloc);
		for (size_t j=0;j<numIndex;j++)
		{
			DoFRef(u, m_vInd[j]) += omega*m_uloc[j];
		}

	}


	virtual size_t size() { return m_vInd.size(); }

protected:
	/// Algebraic indices.
	std::vector<DoFIndex> m_vInd;

	// Memory for local algebra
	vector_type_local m_uloc;
	vector_type_local m_dloc;
	matrix_type_local m_Aloc;

};



//! Collects indices on all elements with v \in Vtx(elem)
template <typename TDomain, typename TAlgebra>
class VertexBasedSubspace : public LocalIndexSubspace<TDomain, TAlgebra,Vertex>
{
protected:
	using TObject = Vertex;

public:
	///	Algebra type
	using algebra_type = TAlgebra;

	///	Vector type
	using vector_type = typename TAlgebra::vector_type;

	///	Matrix type
	using matrix_type = typename TAlgebra::matrix_type;

	/// Base type
	using base_type = LocalIndexSubspace<TDomain, TAlgebra, Vertex>;

	/// CTOR
	VertexBasedSubspace() = default;

	/// virtual DTOR
	virtual ~VertexBasedSubspace() = default;

	bool check(void *obj) const
	{ return true; /*return (dynamic_cast<TObject*> (obj) != nullptr);*/ }

	/// Extract indices for local DoFs.
	void init(void *obj, const vector_type &cvec)
	{
		UG_ASSERT(check(obj), "HUHH: Expecting vertex!");
		Vertex *groupObj = reinterpret_cast<Vertex*> (obj);

		// We will modify index list of base class.
		std::vector<size_t> &vInd = base_type::m_vInd;
		vInd.clear();

		// Union of elements (associated with grouping object).
		using TGridFunction = GridFunction<TDomain, TAlgebra>;
		using TElement = typename GridFunction<TDomain, TAlgebra>::element_type;
		const TGridFunction *c = dynamic_cast<const TGridFunction *> (&cvec);  // Need a grid function here!
		std::vector<TElement*> vElem;
		c->collect_associated(vElem, groupObj);

		// Union of algebraic indices.
		for(size_t i = 0; i < vElem.size(); ++i)
		{
			c->algebra_indices(vElem[i], vInd, false);
		}

		// Remove dublicates.
		if(vElem.size() > 1)
		{
			std::sort(vInd.begin(), vInd.end());
			vInd.erase(std::unique(vInd.begin(), vInd.end()), vInd.end());
		}
	}

};


//! Collects indices on all elements with v \in Vtx(elem)
template <typename TDomain, typename TAlgebra>
class VertexCenteredVankaSubspace : public LocalDoFSubspace<TDomain, TAlgebra, Vertex>
{
public:
	///	Algebra type
	using algebra_type = TAlgebra;

	///	Vector type
	using vector_type = typename TAlgebra::vector_type;

	///	Matrix type
	using matrix_type = typename TAlgebra::matrix_type;

	using grid_function_type = GridFunction<TDomain, TAlgebra>;

	/// Base type
	using base_type = LocalDoFSubspace<TDomain, TAlgebra, Vertex>;

	/// CTOR
	VertexCenteredVankaSubspace(const std::vector<std::string> &vVtxCmp, const std::vector<std::string> &vElemCmp)
		: m_vVtxCmp(vVtxCmp), m_vElemCmp(vElemCmp) {}

	/// virtual DTOR
	virtual ~VertexCenteredVankaSubspace() = default;

	/// Extracts function IDs.
protected:
	void extract_ids(const grid_function_type *c)
	{
		ConstSmartPtr<DoFDistributionInfo> ddinfo =
						c->approx_space()->dof_distribution_info();
		UG_COND_THROW(ddinfo.invalid(), "Requiring valid ddinfo!");
		
		// Vertex functions.
		m_vVtxFct.reserve(m_vVtxCmp.size());
		m_vVtxFct.clear();
		for(size_t i = 0; i < m_vVtxCmp.size(); ++i)
			m_vVtxFct.push_back(ddinfo->fct_id_by_name(m_vVtxCmp[i].c_str()));
		
		// Element functions.
		m_vElemFct.reserve(m_vElemCmp.size());
		m_vElemFct.clear();
		for(size_t i = 0; i < m_vElemCmp.size(); ++i)
			m_vElemFct.push_back(ddinfo->fct_id_by_name(m_vElemCmp[i].c_str()));
	}
public:
	/// OVERRIDE
	bool preprocess(const vector_type &cvec)
	{
		using TGridFunction = GridFunction<TDomain, TAlgebra>;
		const TGridFunction *c = dynamic_cast<const TGridFunction *> (&cvec);  // Need a grid function here!

		UG_COND_THROW(c==nullptr, "Requiring a grid function here!");
		extract_ids(c);

		return true;
	}

	/// Extract indices for single vertex.
	void init(Vertex *groupObj, const vector_type &cvec)
	{
		// We will modify index list of base class.
		typename base_type::index_vector &vInd = base_type::m_vInd;
		vInd.clear();

		// Union of elements (associated with grouping object).
		using TGridFunction = GridFunction<TDomain, TAlgebra>;
		using TElement = typename GridFunction<TDomain, TAlgebra>::element_type;
		const TGridFunction *c = dynamic_cast<const TGridFunction *> (&cvec);  // Need a grid function here!
		std::vector<TElement*> vElem;
		c->collect_associated(vElem, groupObj);
		
		// Collect associated indices.
		for(size_t i = 0; i < vElem.size(); ++i)
		{
			for(size_t f = 0; f < m_vElemFct.size(); ++f)
			c->dof_indices(vElem[i], m_vElemFct[f], vInd, false, false);
		}

		// Collect vertex indices.
		for(size_t f = 0; f < m_vVtxFct.size(); ++f)
			c->dof_indices(groupObj, m_vVtxFct[f], vInd, false, false);
		
		// Remove dublicates.
		if(vElem.size() > 1)
		{
			std::sort(vInd.begin(), vInd.end());
			vInd.erase(std::unique(vInd.begin(), vInd.end()), vInd.end());
		}
		
		/*
		UG_DLOG("VertexBasedVankaSubspace [for " << groupObj << "]:")
		for (typename base_type::index_vector::iterator it = vInd.begin() ; it != vInd.end(); ++it)
		{ UG_DLOG(*it <<","); }
		
		UG_DLOG(std::endl);
		*/
	}



protected:
	std::vector<std::string> m_vVtxCmp;  // primary (vertex) components
	std::vector<std::string> m_vElemCmp; // secondary (element) components

	std::vector<size_t> m_vVtxFct; 		 // mapping to fct IDs
	std::vector<size_t> m_vElemFct;		 // mapping to fct IDs

};


/*
//! Correction for single subspace
template <typename TDomain, typename TAlgebra>
class LocalSubspaceCorrection :
		public ILinearIterator<typename TAlgebra::vector_type>,
		public DebugWritingObject<TAlgebra>
{
public:
	///	Algebra type
	using algebra_type = TAlgebra;

	///	Vector type
	using vector_type = typename TAlgebra::vector_type;

	///	Matrix type
	using matrix_type = typename TAlgebra::matrix_type;


	using TLinearOperator = MatrixOperator<matrix_type, vector_type>;

	virtual ~LocalSubspaceCorrection<TDomain, TAlgebra>() {}

	///	returns the name of iterator
	virtual const char* name() const = 0;

	///	returns if parallel solving is supported
	bool supports_parallel() {return false};

	///	initialize for operator J(u) and linearization point u
	bool init(SmartPtr<TLinearOperator> J, const vector_type& u) {

	};

		///	initialize for linear operator L
		virtual bool init(SmartPtr<TLinearOperator> L) = 0;

		///	Compute new correction c = B*d (w/o updating the defect).
		virtual bool apply(vector_type& c, const vector_type& d)
		{

		}

		///	compute new correction c = B*d and update defect d := d - A*c
		virtual bool apply_update_defect(vector_type& c, vector_type& d)
		{

		}

};


*/


//! Abstract loop
template<typename TGroupObj, typename TDomain, typename TAlgebra>
void ParallelSubspaceCorrectionLoop(const typename TAlgebra::matrix_type& A,
                            GridFunction<TDomain, TAlgebra>& c,
                            const typename TAlgebra::vector_type& d,
                            number omega_relax,
							ILocalSubspace<TDomain, TAlgebra, TGroupObj> &subspace,
							typename GridFunction<TDomain, TAlgebra>::template traits<TGroupObj>::const_iterator objIterBegin,
							typename GridFunction<TDomain, TAlgebra>::template traits<TGroupObj>::const_iterator objIterEnd
)
{
	// Loop over all grouping objects.
	using GroupObjIter = typename GridFunction<TDomain, TAlgebra>::template traits<TGroupObj>::const_iterator;
	for(GroupObjIter iter = objIterBegin; iter != objIterEnd; ++iter)
	{
		// Apply subspace correction (w.r.t. obj.)
		TGroupObj* groupObj = *iter;
		subspace.init(groupObj, c);
		subspace.extract_matrix(A);
		subspace.extract_rhs(d);  // w/o updates in c => parallel
		subspace.update_solution(c, omega_relax);
	}
}

//! Abstract loop over TGroupObj (forward)
template<typename TGroupObj, typename TGroupObjIter, typename TDomain, typename TAlgebra>
void SequentialSubspaceCorrectionLoop(const typename TAlgebra::matrix_type& A,
                            GridFunction<TDomain, TAlgebra>& c,
                            const typename TAlgebra::vector_type& d,
                            number omega_relax,
							ILocalSubspace<TDomain, TAlgebra,TGroupObj> &subspace,
							TGroupObjIter objIterBegin, TGroupObjIter objIterEnd)
{
	// Loop over all grouping objects.
	size_t count=0;
	for(TGroupObjIter iter = objIterBegin; iter != objIterEnd; ++iter)
	{
		// Apply subspace correction (w.r.t. obj.)
		TGroupObj* groupObj = *iter;
		if (groupObj == nullptr)
		{
			UG_LOG("WARNING: Found empty object!: "<< count);
		}
		else
		{
			subspace.init(groupObj, c);
			subspace.extract_matrix(A);
			subspace.extract_rhs(d,A,c); // w/ updates for c => sequential
			subspace.update_solution(c, omega_relax);
		}
		count++;
	}
}

//! Abstract loop (backward)
template<typename TGroupObj, typename TDomain, typename TAlgebra>
void SequentialSubspaceCorrectionLoopBackward(const typename TAlgebra::matrix_type& A,
                            GridFunction<TDomain, TAlgebra>& c,
                            const typename TAlgebra::vector_type& d,
                            number omega_relax,
							ILocalSubspace<TDomain, TAlgebra,TGroupObj> &subspace,
							typename GridFunction<TDomain, TAlgebra>::template traits<TGroupObj>::const_iterator objIterBegin,
							typename GridFunction<TDomain, TAlgebra>::template traits<TGroupObj>::const_iterator objIterEnd
)
{
	// Loop over all grouping objects.
	using GroupObjIter = typename GridFunction<TDomain, TAlgebra>::template traits<TGroupObj>::const_iterator;

	int ncount = 0;
	for(GroupObjIter iter = objIterBegin; iter != objIterEnd; ++iter)
	{ ncount++; }

	// std::cerr << "Found " << ncount << "objects!" << std::endl;

	std::vector<TGroupObj*> objects(ncount);
	for(GroupObjIter iter = objIterBegin; iter != objIterEnd; ++iter)
	{
		// std::cerr << "Found1 " <<  *iter << std::endl;
		objects.push_back(*iter);
	}

	// std::cerr << "Created vector !" << std::endl;

	for(typename std::vector<TGroupObj*>::reverse_iterator riter = objects.rbegin(); riter != objects.rend(); ++riter)
	{
		// Apply subspace correction (w.r.t. obj.)
		TGroupObj* groupObj = *riter;
		if (groupObj == nullptr) continue;

		//std::cerr << "Found2 " <<  groupObj << std::endl;

		subspace.init(groupObj, c);
		subspace.extract_matrix(A);
		subspace.extract_rhs(d,A,c); // w/ updates for c => sequential
		subspace.update_solution(c, omega_relax);
	}
}

template <class TDomain>
class customLexLess{
		public:
			customLexLess(const typename TDomain::position_accessor_type& aaPos): _aaPos(aaPos){}
			bool operator()(Vertex* a, Vertex *b) const
			{ return _aaPos[a] < _aaPos[b]; }
		protected:
			const typename TDomain::position_accessor_type& _aaPos;
};


/// Abstract
template <typename TDomain, typename TAlgebra>
class ISpaceDecomposition {

public:
	///	Algebra type
		using algebra_type = TAlgebra;

	///	Vector type
		using vector_type = typename TAlgebra::vector_type;

	///	Matrix type
		using matrix_type = typename TAlgebra::matrix_type;

	/// Grid function type
		using grid_function_type = GridFunction<TDomain, TAlgebra>;


		virtual ~ISpaceDecomposition()= default;

	/// Initialize decomposition
		virtual void init(const vector_type &x) = 0;

		///	World Dimension
		static constexpr int dim = TDomain::dim;

};


template <typename TDomain, typename TAlgebra>
class FullVertexCover: public ISpaceDecomposition<TDomain, TAlgebra>
{
public:
	using base_type = ISpaceDecomposition<TDomain, TAlgebra>;

	using TGridFunction = GridFunction<TDomain, TAlgebra>;

	static constexpr int dim = TDomain::dim;

	using TElement = typename TGridFunction::element_type;

	FullVertexCover() : m_pSol(nullptr) {}

		/// virtual DTOR
	virtual ~FullVertexCover()= default;

	virtual void init(const typename base_type::vector_type &x)
	{
		auto* m_pSol = dynamic_cast<TGridFunction*>(&x);
		UG_ASSERT(m_pSol !=nullptr, "Error: Cast failed!");
	}

	template <typename TElem>
	typename TGridFunction::template traits<TElem>::const_iterator
	begin() const
	{ return m_pSol->template begin<Vertex>();}

	template <typename TElem>
	typename TGridFunction::template traits<TElem>::const_iterator
	end() const
	{ return m_pSol->template end<Vertex>(); }
protected:
	TGridFunction* m_pSol;

};


///	Sequential subspace correction preconditioner
template <typename TDomain, typename TAlgebra>
class SequentialSubspaceCorrection : public IPreconditioner<TAlgebra>
{
public:
	///	Algebra type
	using algebra_type = TAlgebra;

	///	Vector type
	using vector_type = typename TAlgebra::vector_type;

	///	Matrix type
	using matrix_type = typename TAlgebra::matrix_type;

	///	Matrix Operator type
	using matrix_operator_type = typename IPreconditioner<TAlgebra>::matrix_operator_type;

	/// Subspace for vertices
	using TVertexSubspace = ILocalSubspace <TDomain, TAlgebra, Vertex>;

	///	Base type
	using base_type = IPreconditioner<TAlgebra>;

protected:
	using base_type::set_debug;
	using base_type::debug_writer;
	using base_type::write_debug;

	using grid_function_type = GridFunction<TDomain, TAlgebra>;

public:
	///	default constructor
	SequentialSubspaceCorrection() : m_relax(1.0), m_type("vertex") {};

	///	constructor setting relaxation
	SequentialSubspaceCorrection(number relax) : m_relax(relax), m_type("vertex") {};


	///	Clone
	virtual SmartPtr<ILinearIterator<vector_type> > clone()
	{
		SmartPtr<SequentialSubspaceCorrection<TDomain, TAlgebra> >
		newInst(new SequentialSubspaceCorrection<TDomain, TAlgebra>());
		newInst->set_debug(debug_writer());
		newInst->set_damp(this->damping());
		newInst->set_relax(m_relax);
		newInst->set_type(m_type);
		newInst->set_vertex_subspace(m_spVertexSubspace);
		return newInst;
	}

	///	Destructor
	virtual ~SequentialSubspaceCorrection() = default;

	///	returns if parallel solving is supported
	virtual bool supports_parallel() const {return true;}

	/// set relaxation parameter
	void set_relax(number omega){ m_relax=omega; }

	/// set type
	void set_type(const std::string& type){ m_type=type; }

	/// set subspace
	void set_vertex_subspace(SmartPtr<TVertexSubspace> spVertexSubspace)
	{ m_spVertexSubspace = spVertexSubspace; }


protected:
	///	Name of preconditioner
	virtual const char* name() const
	{return "SequentialSubspaceCorrection";}

	///	Preprocess routine
	virtual bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp)
	{
		PROFILE_BEGIN_GROUP(SSC_preprocess, "algebra ssc");

		// Creating overlap 1 matrix and vectors.
		matrix_type *pA=nullptr;
#ifdef UG_PARALLEL
		if(pcl::NumProcs() > 1)
		{
			UG_ASSERT(0, "SequentialSubspaceCorrection not implemented in parallel. Need to think about this");
			// We should work with element overlap 1 here!
			m_A = *pOp;
			CreateOverlap(m_A);
			m_oD.set_layouts(m_A.layouts());
			m_oC.set_layouts(m_A.layouts());
			m_oD.resize(m_A.num_rows(), false);
			m_oC.resize(m_A.num_rows(), false);
			pA = &m_A;
		}
		else
#endif
		pA = &(*pOp);



		// Checking.
		THROW_IF_NOT_EQUAL(pA->num_rows(), pA->num_cols());
		// UG_COND_THROW(CheckDiagonalInvertible(*pA) == false, name() << ": A has noninvertible diagonal");

		return true;
	}

	virtual bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d)
	{
		PROFILE_BEGIN_GROUP(SSC_step, "algebra ssc");

		using TGridFunction = GridFunction<TDomain, TAlgebra>;
		auto* pC = dynamic_cast<TGridFunction*>(&c);

		UG_COND_THROW(pC == nullptr, "SequentialSubspaceCorrection expects correction to be a GridFunction.");

		//using Element = typename GridFunction<TDomain, TAlgebra>::element_type;
		//using Side = typename GridFunction<TDomain, TAlgebra>::side_type;


		// Set all vector entries to zero.
		pC->set(0.0);

#ifdef UG_PARALLEL
		if(pcl::NumProcs() > 1){
			UG_ASSERT(0, "SequentialSubspaceCorrection not implemented in parallel. Need to think about this");

			// The key problem is that we are iterating over grid functions, not vectors !!!
		}
		else
#endif
		{
			matrix_type &A=*pOp;
			if	(m_type == "vertex")
			{
				// Prepare
				using TVertexIter = typename GridFunction<TDomain, TAlgebra>::template traits<Vertex>::const_iterator;
				TVertexIter objIterBegin = pC->template begin<Vertex>();
				TVertexIter objIterEnd = pC->template end<Vertex>();
			
				m_spVertexSubspace->preprocess(*pC);

				if (true)
				{
					int ncount = 0;
					for(TVertexIter iter = objIterBegin; iter != objIterEnd; ++iter)
					{ ncount++; }

					// std::cerr << "Found " << ncount << "objects!" << std::endl;

					// Create bidirectional container.
					using TObjVector = std::vector<Vertex*>;
					TObjVector objects(ncount);
					objects.clear();
					for(TVertexIter iter = objIterBegin; iter != objIterEnd; ++iter)
					{
							// std::cerr << "Found1 " <<  *iter << std::endl;
							objects.push_back(*iter);
					}

				// Forward loop.
				 SequentialSubspaceCorrectionLoop<Vertex,TObjVector::const_iterator, TDomain,TAlgebra>
						(A, *pC, d, m_relax, *m_spVertexSubspace, objects.begin(), objects.end());

				// Backward loop.
				 SequentialSubspaceCorrectionLoop<Vertex,TObjVector::const_reverse_iterator, TDomain,TAlgebra>
						(A, *pC, d, m_relax, *m_spVertexSubspace, objects.rbegin(), objects.rend());

				} // true


			} else if(m_type == "element") {
/*				using TElemIter = typename GridFunction<TDomain, TAlgebra>::template traits<Element>::const_iterator;
				SequentialSubspaceCorrectionLoop<Element,TElemIter,TDomain,TAlgebra>
						(A, *pC, d, m_relax, *m_spElemSubspace, pC->template begin<Element>(), pC->template end<Element>());
*/
			}
			else UG_THROW("SequentialSubspaceCorrectionStep: wrong patch type '"<<m_type<<"'."
					" Options: vertex.")
#ifdef UG_PARALLEL
			pC->set_storage_type(PST_CONSISTENT);
#endif

			return true;
		}
		return false;
	}

	///	Postprocess routine
	virtual bool postprocess() {return true;}


public:
	using Element = typename GridFunction<TDomain, TAlgebra>::element_type;
protected:
	number m_relax;
	std::string m_type;


	SmartPtr<TVertexSubspace> m_spVertexSubspace;
	
#ifdef UG_PARALLEL

	/// matrix with overlap
	matrix_type m_A;

	///	for overlaps only
	vector_type m_oD;
	vector_type m_oC;
#endif


};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__ELEMENT_GAUSS_SEIDEL__ */
