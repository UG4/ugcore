/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Authors: Christian Wehner, Andreas Vogel
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
 *  Remark: Base on Vanka idea and implementation
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__ELEMENT_GAUSS_SEIDEL__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__ELEMENT_GAUSS_SEIDEL__

#include "lib_algebra/operator/interface/preconditioner.h"

#include <vector>
#include <algorithm>

#ifdef UG_PARALLEL
	#include "pcl/pcl_util.h"
	#include "lib_algebra/parallelization/parallelization_util.h"
	#include "lib_algebra/parallelization/parallel_matrix_overlap_impl.h"
#endif

namespace ug{

template<typename TGroupObj, typename TDomain, typename TAlgebra>
void ElementGaussSeidelStep(const typename TAlgebra::matrix_type& A,
                            GridFunction<TDomain, TAlgebra>& c,
                            const typename TAlgebra::vector_type& d,
                            number relax)
{
	typedef typename TAlgebra::matrix_type::const_row_iterator const_row_iterator;
	const static int blockSize = TAlgebra::blockSize;

	// memory for local algebra
	DenseVector< VariableArray1<number> > s;
	DenseVector< VariableArray1<number> > x;
	DenseMatrix< VariableArray2<number> > mat;
	std::vector<size_t> vInd;
	
	// set all vector entries to zero
	c.set(0.0);
#ifdef UG_PARALLEL
	c.set_storage_type(PST_ADDITIVE);
#endif
	typedef typename GridFunction<TDomain, TAlgebra>::element_type Element;
	std::vector<Element*> vElem;

	// loop all grouping objects
	typedef typename GridFunction<TDomain, TAlgebra>::template traits<TGroupObj>::const_iterator GroupObjIter;
	for(GroupObjIter iter = c.template begin<TGroupObj>();
					 iter != c.template end<TGroupObj>(); ++iter){

		// get grouping obj
		TGroupObj* groupObj = *iter;

		// collect elems associated to grouping object
		c.collect_associated(vElem, groupObj);

		// get all algebraic indices on element
		vInd.clear();
		for(size_t i = 0; i < vElem.size(); ++i)
			c.algebra_indices(vElem[i], vInd, false);

		// check for doublicates
		if(vElem.size() > 1){
		    std::sort(vInd.begin(), vInd.end());
		    vInd.erase(std::unique(vInd.begin(), vInd.end()), vInd.end());
		}

		// get number of indices on patch
		const size_t numIndex = vInd.size();

		// fill local block matrix
		bool bFound;
		mat.resize(numIndex, numIndex);
		mat = 0.0;
		for (size_t j = 0; j<numIndex; j++){
			for (size_t k=0;k<numIndex;k++){
				const_row_iterator it = A.get_connection(vInd[j],vInd[k], bFound);
				if(bFound){
					mat.subassign(j*blockSize,k*blockSize,it.value());
				}
			};
		}

		// compute s[j] := d[j] - sum_k A(j,k)*c[k]
		// note: the loop over k is the whole matrix row (not only selected indices)
		s.resize(numIndex);
		for (size_t j = 0; j<numIndex; j++){
			typename TAlgebra::vector_type::value_type sj = d[vInd[j]];
			for(const_row_iterator it = A.begin_row(vInd[j]);
									it != A.end_row(vInd[j]) ; ++it){
				MatMultAdd(sj, 1.0, sj, -1.0, it.value(), c[it.index()]);
			};
			s.subassign(j*blockSize,sj);
		};

		// solve block
		x.resize(numIndex);
		InverseMatMult(x,1,mat,s);
		for (size_t j=0;j<numIndex;j++){
			c[vInd[j]] += relax*x[j];
		};
	}
}


///	ElementGaussSeidel Preconditioner
template <typename TDomain, typename TAlgebra>
class ElementGaussSeidel : public IPreconditioner<TAlgebra>
{
	public:
	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	///	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Matrix Operator type
		typedef typename IPreconditioner<TAlgebra>::matrix_operator_type matrix_operator_type;

	///	Base type
		typedef IPreconditioner<TAlgebra> base_type;

	protected:
		using base_type::set_debug;
		using base_type::debug_writer;
		using base_type::write_debug;

	public:
	///	default constructor
		ElementGaussSeidel() : m_relax(1), m_type("element"){};

	///	constructor setting relaxation
		ElementGaussSeidel(number relax) : m_relax(relax), m_type("element") {};

	///	constructor setting type
		ElementGaussSeidel(const std::string& type) : m_relax(1), m_type(type) {};

	///	constructor setting relaxation and type
		ElementGaussSeidel(number relax, const std::string& type) : m_relax(relax), m_type(type) {};

	///	Clone
		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			SmartPtr<ElementGaussSeidel<TDomain, TAlgebra> >
							newInst(new ElementGaussSeidel<TDomain, TAlgebra>());
			newInst->set_debug(debug_writer());
			newInst->set_damp(this->damping());
			newInst->set_relax(m_relax);
			newInst->set_type(m_type);
			return newInst;
		}

	///	Destructor
		virtual ~ElementGaussSeidel()
		{};

	///	returns if parallel solving is supported
		virtual bool supports_parallel() const {return true;}

	/// set relaxation parameter
		void set_relax(number omega){m_relax=omega;};

	/// set type
		void set_type(const std::string& type){m_type=type;};

	protected:
	///	Name of preconditioner
		virtual const char* name() const {return "ElementGaussSeidel";}

	///	Preprocess routine
		virtual bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp)
		{
#ifdef UG_PARALLEL
			if(pcl::NumProcs() > 1)
			{
				//	copy original matrix
				MakeConsistent(*pOp, m_A);
				//	set zero on slaves
				//	std::vector<IndexLayout::Element> vIndex;
				//CollectUniqueElements(vIndex,  m_A.layouts()->slave());
				//SetDirichletRow(m_A, vIndex);
			}
#endif
			return true;
		}

		virtual bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d)
		{
			GridFunction<TDomain, TAlgebra>* pC
							= dynamic_cast<GridFunction<TDomain, TAlgebra>*>(&c);
			if(pC == NULL)
				UG_THROW("ElementGaussSeidel expects correction to be a GridFunction.");


			typedef typename GridFunction<TDomain, TAlgebra>::element_type Element;
			typedef typename GridFunction<TDomain, TAlgebra>::side_type Side;

			
#ifdef UG_PARALLEL
			if(pcl::NumProcs() > 1){
			         // make defect unique
				SmartPtr<vector_type> spDtmp = d.clone();
				spDtmp->change_storage_type(PST_UNIQUE);
				
				// execute step
				if (m_type == "element") ElementGaussSeidelStep<Element,TDomain,TAlgebra>(m_A, *pC, *spDtmp, m_relax);
				else UG_THROW("ElementGaussSeidel: wrong patch type '"<<m_type<<"'."
					      " Options: element, side, face, edge, vertex.");
				
				// make correction consistent
				// pC->set_storage_type(PST_CONSISTENT);
				pC->change_storage_type(PST_CONSISTENT);
				return true;
			}
			else
#endif
			  {
			    matrix_type &A=*pOp; 
			    if (m_type == "element") ElementGaussSeidelStep<Element,TDomain,TAlgebra>(A, *pC, d, m_relax);
			    else if	(m_type == "side") ElementGaussSeidelStep<Side,TDomain,TAlgebra>(A, *pC, d, m_relax);
			    else if	(m_type == "face") ElementGaussSeidelStep<Face,TDomain,TAlgebra>(A, *pC, d, m_relax);
			    else if	(m_type == "edge") ElementGaussSeidelStep<Edge,TDomain,TAlgebra>(A, *pC, d, m_relax);
			    else if	(m_type == "vertex") ElementGaussSeidelStep<Vertex,TDomain,TAlgebra>(A, *pC, d, m_relax);
			    else UG_THROW("ElementGaussSeidel: wrong patch type '"<<m_type<<"'."
					  " Options: element, side, face, edge, vertex.")
#ifdef UG_PARALLEL
				   pC->set_storage_type(PST_CONSISTENT);
#endif

			    return true;
			  }
		}

	///	Postprocess routine
		virtual bool postprocess() {return true;}

	protected:
#ifdef UG_PARALLEL
		matrix_type m_A;
#endif

		number m_relax;

		std::string m_type;
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__ELEMENT_GAUSS_SEIDEL__ */
