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


#define SCHUR_MOD
namespace ug{

template<typename TGroupObj, typename TDomain, typename TAlgebra>
void ElementGaussSeidelStep(const typename TAlgebra::matrix_type& A,
                            GridFunction<TDomain, TAlgebra>& c,
                            const typename TAlgebra::vector_type& d,
                            number relax
#ifdef SCHUR_MOD
							, const std::vector<std::string>& cmp, number alpha, bool elim_off_diag=false
#endif
							)
{
	using const_row_iterator = typename TAlgebra::matrix_type::const_row_iterator;
	static constexpr int blockSize = TAlgebra::blockSize;

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
	using Element = typename GridFunction<TDomain, TAlgebra>::element_type;
	std::vector<Element*> vElem;

	// loop all grouping objects
	using GroupObjIter = typename GridFunction<TDomain, TAlgebra>::template traits<TGroupObj>::const_iterator;
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
		for (size_t j = 0; j<numIndex; j++)
		{
			typename TAlgebra::vector_type::value_type sj = d[vInd[j]];
			for(const_row_iterator it = A.begin_row(vInd[j]); it != A.end_row(vInd[j]) ; ++it)
			{
				MatMultAdd(sj, 1.0, sj, -1.0, it.value(), c[it.index()]);
			};

			s.subassign(j*blockSize,sj);
		};


#ifdef SCHUR_MOD

		// tokenize passed functions
		ConstSmartPtr<DoFDistributionInfo> ddinfo =
			c.approx_space()->dof_distribution_info();

		// get fct IDs of selected functions
		std::vector<size_t> vFct;
		for(size_t i = 0; i < cmp.size(); ++i)
			vFct.push_back(ddinfo->fct_id_by_name(cmp[i].c_str()));


		DenseVector< VariableArray1<number> > weights;
		weights.resize(numIndex);

		if (vFct.size()==1)
		{
			std::vector<DoFIndex> vSchurDoFIndex;
			c.dof_indices(groupObj, vFct[0], vSchurDoFIndex, true);

			// identify schur rows
			UG_ASSERT(blockSize==1, "Element GS: Elimination does only work for blockSize==1!")

			std::vector<size_t> vIndElim; 		// global indices for elimination
			// std::vector<size_t> vIndKeep;	// global indices kept

			const size_t numElim = vSchurDoFIndex.size();
			for (size_t j = 0; j<numElim; j++)
			{
				vIndElim.push_back(vSchurDoFIndex[j][0]);
			}

			std::vector<size_t> mapElim; 		// local indices for elimination (j_elim -> j_local)
			std::vector<size_t> mapKeep;	    // local indices for elimination (j_nonelim -> j_local)

			// compute weights & fill map
			for (size_t j = 0; j<numIndex; j++)
			{
				std::vector<size_t>::iterator it = find (vIndElim.begin(), vIndElim.end(), vInd[j]);
				if (it != vIndElim.end()) {
					// eliminate this row
					mapElim.push_back(j);
				} else {
					// keep this row
					mapKeep.push_back(j); // vIndKeep.push_back(vInd[j]);
				}
			}

			const size_t numKeep = mapKeep.size();// vIndKeep.size();
			UG_ASSERT((numElim+numKeep == numIndex), "Map elim does not match!");
			UG_ASSERT(mapElim.size()==numElim, "Map elim does not match!");
			UG_ASSERT(mapKeep.size()==numKeep, "Map mon elim does not match!");


			// extract mat_keep (needed for inversion)
			DenseMatrix< VariableArray2<number> > matKeep;
			matKeep.resize(numKeep, numKeep);
			matKeep = 0.0;
			for (size_t j = 0; j<numKeep; j++)
			{
				for (size_t k=0;k<numKeep;k++)
				{
					matKeep(j,k) = mat(mapKeep[j], mapKeep[k]);
				}
			}

			if (false) {

				// compute mat elim
				// compute contribution Bi^T A^-1 Cj to schur complement
				// std::cout << "S (" << numElim << "):" << std::endl;

				// compute
				DenseVector< VariableArray1<number> > schur_cj; schur_cj.resize(numKeep);
				DenseVector< VariableArray1<number> > schur_yj; schur_yj.resize(numKeep);

				DenseMatrix< VariableArray2<number> > matElim; matElim.resize(numElim, numElim);
				matElim = 0.0;


				for (size_t j = 0; j<numElim; j++)
				{
					for (size_t k=0; k<numKeep; k++)
					{
						schur_cj[k]  = mat(mapKeep[k], mapElim[j]);
					}

					InverseMatMult(schur_yj, 1.0, matKeep, schur_cj);

					for (size_t i = 0; i<numElim; i++)
					{
						// compute elim_ij
						matElim(i,j) = 0.0;
						for (size_t k=0; k<numKeep; k++)
						{
							number schur_bik = mat(mapElim[i], mapKeep[k]);
							matElim(i,j) += schur_bik*schur_yj[k];

						}

						// replace/update matrix
						// std::cout << mat(mapElim[i], mapElim[j]) << "+" << matElim(i,j) << "=";
						mat(mapElim[i], mapElim[j]) -= alpha*matElim(i,j); //* alpha
						// std::cout << mat(mapElim[i], mapElim[j]) << std::endl;

						// std::cout << matElim;
					}



				}
			}


			if (false) // eliminate off-diag rows B, C?
			{
				for (size_t j = 0; j<numElim; j++)
				{
					for (size_t k=0; k<numKeep; k++)
					{
						mat(mapElim[j], mapKeep[k]) = 0.0;
						mat(mapKeep[k], mapElim[j]) = 0.0;
					}

					for (size_t k=0; k<numElim; k++)
					{

						if (j==k) continue;
							mat(mapElim[j], mapElim[k]) = 0.0;
							mat(mapElim[k], mapElim[j]) = 0.0;
					}
				}
			}



			if (false) // rescale elim part of matrix?
			{
				for (size_t i = 0; i<numElim; i++)
				{
					for (size_t j=0; j<numElim; j++)
					{
						mat(mapElim[i], mapElim[j]) *= alpha;

					}
				}
			}

			// std::cout << mat;

			// std::cout << std::endl;
		} // (vFct.size()==1)

#endif
		// solve block
		x.resize(numIndex);
		InverseMatMult(x, 1.0, mat, s);
		for (size_t j=0;j<numIndex;j++)
		{
			c[vInd[j]] += relax*x[j];
		}


	}
}


///	ElementGaussSeidel Preconditioner
template <typename TDomain, typename TAlgebra>
class ElementGaussSeidel : public IPreconditioner<TAlgebra>
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

	///	Base type
		using base_type = IPreconditioner<TAlgebra>;

	protected:
		using base_type::set_debug;
		using base_type::debug_writer;
		using base_type::write_debug;

		using grid_function_type = GridFunction<TDomain, TAlgebra>;

	public:
	///	default constructor
		ElementGaussSeidel() : m_relax(1.0), m_type("element"), m_schur_alpha(1.0), m_elim_off_diag(false) {};

	///	constructor setting relaxation
		explicit ElementGaussSeidel(number relax) : m_relax(relax), m_type("element"), m_schur_alpha(1.0), m_elim_off_diag(false) {};

	///	constructor setting type
		explicit ElementGaussSeidel(const std::string& type) : m_relax(1.0), m_type(type), m_schur_alpha(1.0), m_elim_off_diag(false) {};

	///	constructor setting relaxation and type
		ElementGaussSeidel(number relax, const std::string& type) : m_relax(relax), m_type(type), m_schur_alpha(1.0), m_elim_off_diag(false) {};

	///	Clone
		SmartPtr<ILinearIterator<vector_type> > clone() override {
			SmartPtr<ElementGaussSeidel >
							newInst(new ElementGaussSeidel());
			newInst->set_debug(debug_writer());
			newInst->set_damp(this->damping());
			newInst->set_relax(m_relax);
			newInst->set_type(m_type);

			newInst->m_schur_cmp =m_schur_cmp;
			newInst->m_schur_alpha = m_schur_alpha;
			return newInst;
		}

	///	Destructor
		~ElementGaussSeidel() override = default;

	///	returns if parallel solving is supported
		bool supports_parallel() const override {return true;}

	/// set relaxation parameter
		void set_relax(number omega){m_relax=omega;};

	/// set type
		void set_type(const std::string& type){m_type=type;};

		void select_schur_cmp(const std::vector<std::string>& cmp, number alpha)
		{
			m_schur_cmp =cmp;
			m_schur_alpha = alpha;
		};

		void set_elim_offdiag(bool elim) { m_elim_off_diag=elim;};

	protected:
	///	Name of preconditioner
		const char* name() const override {return "ElementGaussSeidel";}

	///	Preprocess routine
		bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp) override {
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

		bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d) override {
			GridFunction<TDomain, TAlgebra>* pC
							= dynamic_cast<GridFunction<TDomain, TAlgebra>*>(&c);
			if(pC == nullptr)
				UG_THROW("ElementGaussSeidel expects correction to be a GridFunction.");


			using Element = typename GridFunction<TDomain, TAlgebra>::element_type;
			using Side = typename GridFunction<TDomain, TAlgebra>::side_type;

			
#ifdef UG_PARALLEL
			if(pcl::NumProcs() > 1){
			         // make defect unique
				SmartPtr<vector_type> spDtmp = d.clone();
				spDtmp->change_storage_type(PST_UNIQUE);
				
				// execute step
				if (m_type == "element") ElementGaussSeidelStep<Element,TDomain,TAlgebra>(m_A, *pC, *spDtmp, m_relax, m_schur_cmp, m_schur_alpha);
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
			    if (m_type == "element") ElementGaussSeidelStep<Element,TDomain,TAlgebra>(A, *pC, d, m_relax, m_schur_cmp, m_schur_alpha);
			    else if	(m_type == "side") ElementGaussSeidelStep<Side,TDomain,TAlgebra>(A, *pC, d, m_relax, m_schur_cmp, m_schur_alpha);
			    else if	(m_type == "face") ElementGaussSeidelStep<Face,TDomain,TAlgebra>(A, *pC, d, m_relax, m_schur_cmp, m_schur_alpha);
			    else if	(m_type == "edge") ElementGaussSeidelStep<Edge,TDomain,TAlgebra>(A, *pC, d, m_relax, m_schur_cmp, m_schur_alpha);
			    else if	(m_type == "vertex") ElementGaussSeidelStep<Vertex,TDomain,TAlgebra>(A, *pC, d, m_relax, m_schur_cmp, m_schur_alpha);
			    else UG_THROW("ElementGaussSeidel: wrong patch type '"<<m_type<<"'."
					  " Options: element, side, face, edge, vertex.")
#ifdef UG_PARALLEL
				   pC->set_storage_type(PST_CONSISTENT);
#endif

			    return true;
			  }
		}

	///	Postprocess routine
		bool postprocess() override {return true;}

	protected:
#ifdef UG_PARALLEL
		matrix_type m_A;
#endif

		number m_relax;
		std::string m_type;

#ifdef SCHUR_MOD
		std::vector<std::string> m_schur_cmp;
		number m_schur_alpha;
		bool m_elim_off_diag;
#endif

};

} // end namespace ug

#endif