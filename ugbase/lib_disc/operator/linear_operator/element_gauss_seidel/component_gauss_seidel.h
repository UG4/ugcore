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

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__COMPONENT_GAUSS_SEIDEL__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__COMPONENT_GAUSS_SEIDEL__

#include "lib_algebra/operator/interface/preconditioner.h"

#include <vector>
#include <algorithm>

#ifdef UG_PARALLEL
	#include "pcl/pcl_util.h"
	#include "lib_algebra/parallelization/parallelization_util.h"
	#include "lib_algebra/parallelization/parallel_matrix_overlap_impl.h"
#endif

namespace ug{

///	ComponentGaussSeidel Preconditioner
template <typename TDomain, typename TAlgebra>
class ComponentGaussSeidel : public IPreconditioner<TAlgebra>
{
	public:
	/// World dimension
		typedef GridFunction<TDomain, TAlgebra> GF;
		typedef typename GF::element_type Element;
		typedef typename GF::side_type Side;
		static const int dim = TDomain::dim;

	///	Algebra types
	/// \{
		typedef TAlgebra algebra_type;
		typedef typename TAlgebra::vector_type vector_type;
		typedef typename TAlgebra::matrix_type matrix_type;
	/// \}

	protected:
		typedef IPreconditioner<TAlgebra> base_type;
		using base_type::set_debug;
		using base_type::debug_writer;
		using base_type::write_debug;

	public:
	///	default constructor
		ComponentGaussSeidel(const std::vector<std::string>& vFullRowCmp)
			: m_bInit(false), m_relax(1), m_alpha(1.0), m_beta(1.0), m_bWeighted(false), m_vFullRowCmp(vFullRowCmp) {}

	///	constructor setting relaxation and type
		ComponentGaussSeidel(number relax, const std::vector<std::string>& vFullRowCmp)
			: m_bInit(false), m_relax(relax),  m_alpha(1.0), m_beta(1.0), m_bWeighted(false), m_vFullRowCmp(vFullRowCmp){}

	///	constructor setting relaxation and type
		ComponentGaussSeidel(number relax, const std::vector<std::string>& vFullRowCmp,
		                     const std::vector<int>& vSmooth, const std::vector<number>& vDamp)
		: m_bInit(false), m_relax(relax),  m_alpha(1.0), m_beta(1.0), m_bWeighted(false), m_vFullRowCmp(vFullRowCmp),
		  m_vGroupObj(vSmooth), m_vDamp(vDamp) {};

	///	Clone
		virtual SmartPtr<ILinearIterator<vector_type> > clone();

	///	returns if parallel solving is supported
		virtual bool supports_parallel() const {return true;}

		void set_alpha(number alpha)
		{m_alpha = alpha;}

		void set_beta(number beta)
		{m_beta = beta;}

		void set_weights(bool b)
		{m_bWeighted = b;}

	protected:
	///	Name of preconditioner
		virtual const char* name() const {return "ComponentGaussSeidel";}

	///	Preprocess routine
		virtual bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp);

		struct DimCache{
			std::vector<std::vector<DoFIndex> > vvDoFIndex;
			std::vector<DenseMatrix< VariableArray2<number> > > vBlockInverse;
		};

		void apply_blocks(const matrix_type& A, GF& c,
		                  const vector_type& d, number relax,
		                  const DimCache& dimCache,
		                  bool bReverse);

		void apply_blocks_weighted(const matrix_type& A, GF& c,
				                  const vector_type& d, number relax,
				                  const DimCache& dimCache,
				                  bool bReverse);

		template<typename TGroupObj>
		void extract_by_grouping(std::vector<std::vector<DoFIndex> >& vvDoFIndex,
		                         const GF& c,
		                         const std::vector<size_t>& vFullRowCmp,
		                         const std::vector<size_t>& vRemainCmp);

		void extract_blocks(const matrix_type& A, const GF& c);

	///	step method
		virtual bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d);

	///	Postprocess routine
		virtual bool postprocess() {return true;}

	protected:
#ifdef UG_PARALLEL
		matrix_type m_A;
#endif

	///	init flag
		bool m_bInit;

	///	relaxing parameter
		number m_relax;
		number m_alpha;
		number m_beta;

		bool m_bWeighted;
		SmartPtr<GF> m_weight;

	///	components, whose matrix row must be fulfilled completely on gs-block
		std::vector<std::string> m_vFullRowCmp;

	/// smooth order
		std::vector<int> m_vGroupObj;
		std::vector<number> m_vDamp;

	///	caching storage
		DimCache m_vDimCache[VOLUME+1];
};

template <typename TDomain, typename TAlgebra>
bool ComponentGaussSeidel<TDomain, TAlgebra>::
preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp)
{
#ifdef UG_PARALLEL
	if(pcl::NumProcs() > 1)
	{
		//	copy original matrix
		MakeConsistent(*pOp, m_A);
		//	set zero on slaves
		std::vector<IndexLayout::Element> vIndex;
		CollectUniqueElements(vIndex,  m_A.layouts()->slave());
		SetDirichletRow(m_A, vIndex);
	}
#endif

	m_bInit = false;

	return true;
}

template <typename TDomain, typename TAlgebra>
void ComponentGaussSeidel<TDomain, TAlgebra>::
apply_blocks(const matrix_type& A, GF& c,
             const vector_type& d, number relax,
             const DimCache& dimCache,
             bool bReverse)
{
// 	memory for local algebra
	DenseVector< VariableArray1<number> > s;
	DenseVector< VariableArray1<number> > x;

//	get caches
	const std::vector<std::vector<DoFIndex> >& vvDoFIndex = dimCache.vvDoFIndex;
	const std::vector<DenseMatrix< VariableArray2<number> > >& vBlockInverse = dimCache.vBlockInverse;

//	loop blocks
	for(size_t b = 0; b < vvDoFIndex.size(); ++b)
	{
		size_t block = (!bReverse) ? b : (vvDoFIndex.size()-1 - b);

	//	get storage
		const std::vector<DoFIndex>& vDoFIndex = vvDoFIndex[block];
		const DenseMatrix<VariableArray2<number> >& BlockInv = vBlockInverse[block];
		const size_t numIndex = vDoFIndex.size();

	// 	compute s[j] := d[j] - sum_k A(j,k)*c[k]
	// 	note: the loop over k is the whole matrix row (not only selected indices)
		typedef typename matrix_type::const_row_iterator const_row_iterator;
		const static int blockSize = TAlgebra::blockSize;
		s.resize(numIndex);
		for (size_t j = 0; j<numIndex; j++){
			typename vector_type::value_type sj = d[vDoFIndex[j][0]];
			for(const_row_iterator it = A.begin_row(vDoFIndex[j][0]);
									it != A.end_row(vDoFIndex[j][0]) ; ++it){
				MatMultAdd(sj, 1.0, sj, -1.0, it.value(), c[it.index()]);
			};
			s.subassign(j*blockSize,sj);
		};

	// 	solve block
		x.resize(numIndex);
		MatMult(x, relax, BlockInv, s);

	//	add to global correction
		for (size_t j=0;j<numIndex;j++)
			DoFRef(c, vDoFIndex[j]) += x[j];
	}
}


template <typename TDomain, typename TAlgebra>
void ComponentGaussSeidel<TDomain, TAlgebra>::
apply_blocks_weighted(const matrix_type& A, GF& c,
             const vector_type& d, number relax,
             const DimCache& dimCache,
             bool bReverse)
{
// 	memory for local algebra
	DenseVector< VariableArray1<number> > s;
	DenseVector< VariableArray1<number> > x;

//	get caches
	const std::vector<std::vector<DoFIndex> >& vvDoFIndex = dimCache.vvDoFIndex;
	const std::vector<DenseMatrix< VariableArray2<number> > >& vBlockInverse = dimCache.vBlockInverse;

//	loop blocks
	for(size_t b = 0; b < vvDoFIndex.size(); ++b)
	{
		size_t block = (!bReverse) ? b : (vvDoFIndex.size()-1 - b);

	//	get storage
		const std::vector<DoFIndex>& vDoFIndex = vvDoFIndex[block];
		const DenseMatrix<VariableArray2<number> >& BlockInv = vBlockInverse[block];
		const size_t numIndex = vDoFIndex.size();

	// 	compute s[i] := d[i] - sum_k A(i,k)*c[k]
	// 	note: the loop over k is the whole matrix row (not only selected indices)
		typedef typename matrix_type::const_row_iterator const_row_iterator;
		const static int blockSize = TAlgebra::blockSize;
		s.resize(numIndex);
		for (size_t i = 0; i<numIndex; i++)
		{
			typename vector_type::value_type si = d[vDoFIndex[i][0]];
			for(const_row_iterator it = A.begin_row(vDoFIndex[i][0]);
									it != A.end_row(vDoFIndex[i][0]) ; ++it)
			{
				MatMultAdd(si, 1.0, si, -1.0, it.value(), c[it.index()]);
			}
			s.subassign(i*blockSize, si);
		}

	// restrict
		for (size_t i = 0; i<numIndex; ++i)
		{
			number wi = sqrt(DoFRef(*m_weight, vDoFIndex[i]));
			UG_ASSERT(wi!=0.0, "Huhh: Zero weights found!" << i << "dofIndex" << vDoFIndex[i])
			s[i] /= wi;
		}

	// 	solve block
		x.resize(numIndex);
		MatMult(x, relax, BlockInv, s);

	//	add to global correction
		for (size_t i=0; i<numIndex;i++)
		{
			number wi = sqrt(DoFRef(*m_weight, vDoFIndex[i]));
			DoFRef(c, vDoFIndex[i]) += x[i]/wi;
		}
	}
}

template <typename TDomain, typename TAlgebra>
template<typename TGroupObj>
void ComponentGaussSeidel<TDomain, TAlgebra>::
extract_by_grouping(std::vector<std::vector<DoFIndex> >& vvDoFIndex,
                    const GF& c,
                    const std::vector<size_t>& vFullRowCmp,
                    const std::vector<size_t>& vRemainCmp)
{
// 	memory for local algebra
	std::vector<DoFIndex> vFullRowDoFIndex;
	std::vector<Element*> vElem;

//	clear indices
//	\todo: improve this by precomputing size
	vvDoFIndex.clear();

// loop all grouping objects
	typedef typename GF::template traits<TGroupObj>::const_iterator GroupObjIter;
	for(GroupObjIter iter = c.template begin<TGroupObj>();
					 iter != c.template end<TGroupObj>(); ++iter)
	{
	// 	get grouping obj
		TGroupObj* groupObj = *iter;

	// 	get all dof indices on obj associated to cmps
		vFullRowDoFIndex.clear();
		for(size_t f = 0; f < vFullRowCmp.size(); ++f)
			c.inner_dof_indices(groupObj, vFullRowCmp[f], vFullRowDoFIndex, false);

	// 	check if equation present
		if(vFullRowDoFIndex.empty())
			UG_THROW("CGS: Should not happen.");

	// 	get all algebraic indices on element
		if(TGroupObj::dim <= VERTEX){
			std::vector<Edge*> vSub;
			c.collect_associated(vSub, groupObj);
			for(size_t i = 0; i < vSub.size(); ++i)
				for(size_t f = 0; f < vFullRowCmp.size(); ++f)
					c.inner_dof_indices(vSub[i], vFullRowCmp[f], vFullRowDoFIndex, false);
		}
		if(TGroupObj::dim <= EDGE){
			std::vector<Face*> vSub;
			c.collect_associated(vSub, groupObj);
			for(size_t i = 0; i < vSub.size(); ++i)
				for(size_t f = 0; f < vFullRowCmp.size(); ++f)
					c.inner_dof_indices(vSub[i], vFullRowCmp[f], vFullRowDoFIndex, false);
		}

	// 	collect elems associated to grouping object
		c.collect_associated(vElem, groupObj);

	// 	get all algebraic indices on element
		std::vector<DoFIndex> vDoFIndex;
		for(size_t i = 0; i < vElem.size(); ++i)
			for(size_t f = 0; f < vRemainCmp.size(); ++f)
				c.dof_indices(vElem[i], vRemainCmp[f], vDoFIndex, false, false);

	// 	check for dublicates
		if(vElem.size() > 1){
		    std::sort(vDoFIndex.begin(), vDoFIndex.end());
		    vDoFIndex.erase(std::unique(vDoFIndex.begin(), vDoFIndex.end()), vDoFIndex.end());
		}

	// 	concat all indices
		vDoFIndex.insert(vDoFIndex.end(), vFullRowDoFIndex.begin(), vFullRowDoFIndex.end());

	//	add
		vvDoFIndex.push_back(vDoFIndex);
	}
}

template <typename TDomain, typename TAlgebra>
void ComponentGaussSeidel<TDomain, TAlgebra>::
extract_blocks(const matrix_type& A, const GF& c)
{

	ConstSmartPtr<DoFDistributionInfo> ddinfo =
							c.approx_space()->dof_distribution_info();

	// tokenize passed functions
	if(m_vFullRowCmp.empty())
		UG_THROW("ComponentGaussSeidel: No components set.")

	// get ids of selected functions
	std::vector<size_t> vFullRowCmp;
	for(size_t i = 0; i < m_vFullRowCmp.size(); ++i)
		vFullRowCmp.push_back(ddinfo->fct_id_by_name(m_vFullRowCmp[i].c_str()));

	// compute remaining functions
	std::vector<size_t> vRemainCmp;
	for(size_t f = 0; f < ddinfo->num_fct(); ++f)
		if(std::find(vFullRowCmp.begin(), vFullRowCmp.end(), f) == vFullRowCmp.end())
			vRemainCmp.push_back(f);

	if(m_vGroupObj.empty()){
		for(int d = VERTEX; d <= VOLUME; ++d){
			bool bCarryDoFs = false;
			for(size_t f = 0; f < vFullRowCmp.size(); ++f)
				if(ddinfo->max_fct_dofs(vFullRowCmp[f], d) > 0)
					bCarryDoFs = true;
			if(bCarryDoFs)
				m_vGroupObj.push_back(d);
		}
	}

//	extract for each dim-grouping objects
	for(int d = VERTEX; d <= VOLUME; ++d)
	{
	//	only extract if needed
		if(std::find(m_vGroupObj.begin(), m_vGroupObj.end(), d) == m_vGroupObj.end())
			continue;

	//	only extract if carrying dofs
		bool bCarryDoFs = false;
		for(size_t f = 0; f < vFullRowCmp.size(); ++f)
			if(ddinfo->max_fct_dofs(vFullRowCmp[f], d) > 0)
				bCarryDoFs = true;
		if(!bCarryDoFs)
			continue;

	//	extract
		std::vector<std::vector<DoFIndex> >& vvDoFIndex = m_vDimCache[d].vvDoFIndex;
		switch(d){
			case VERTEX: extract_by_grouping<Vertex>(vvDoFIndex, c, vFullRowCmp, vRemainCmp); break;
			case EDGE:   extract_by_grouping<Edge>(vvDoFIndex, c, vFullRowCmp, vRemainCmp); break;
			case FACE:   extract_by_grouping<Face>(vvDoFIndex, c, vFullRowCmp, vRemainCmp); break;
			case VOLUME: extract_by_grouping<Volume>(vvDoFIndex, c, vFullRowCmp, vRemainCmp); break;
			default: UG_THROW("wrong dim");
		}
	}

	// compute weights v(i) (= number of groups j that i belongs to)

	if (m_bWeighted)
	{
		m_weight = c.clone_without_values();
		m_weight->set(0.0);
		for(int d = VERTEX; d <= VOLUME; ++d)
		{
			std::vector<std::vector<DoFIndex> >& vvDoFIndex = m_vDimCache[d].vvDoFIndex;

			for(size_t j = 0; j < vvDoFIndex.size(); ++j)
			{

				std::vector<DoFIndex>& vDoFIndex = vvDoFIndex[j];
				const size_t numIndex = vDoFIndex.size();

					// count number of connections for this row
					DoFRef(*m_weight, vDoFIndex[numIndex-1]) = 1.0;
					for (size_t i =0; i<numIndex-1; ++i)
						DoFRef(*m_weight, vDoFIndex[i]) += 1.0;
			}

		}
	}
	//	extract local matrices
	for(int d = VERTEX; d <= VOLUME; ++d)
	{
		std::vector<std::vector<DoFIndex> >& vvDoFIndex = m_vDimCache[d].vvDoFIndex;
		std::vector<DenseMatrix<VariableArray2<number> > >& vBlockInverse = m_vDimCache[d].vBlockInverse;

		vBlockInverse.resize(vvDoFIndex.size());
		for(size_t b = 0; b < vvDoFIndex.size(); ++b)
		{
		//	get storage
			std::vector<DoFIndex>& vDoFIndex = vvDoFIndex[b];
			DenseMatrix<VariableArray2<number> >& BlockInv = vBlockInverse[b];

		// 	get number of indices on patch
			const size_t numIndex = vDoFIndex.size();
			//std::cerr << "locSize" << numIndex << std::endl;

		// 	fill local block matrix
			BlockInv.resize(numIndex, numIndex);
			BlockInv = 0.0;

		// 	copy matrix rows (only including cols of selected indices)
			for (size_t i = 0; i < numIndex; i++)
			{
				// Diag (A)
				BlockInv(i,i) = m_alpha*DoFRef(A, vDoFIndex[i], vDoFIndex[i]);
				BlockInv(numIndex-1,i) = DoFRef(A, vDoFIndex[numIndex-1], vDoFIndex[i]);
				BlockInv(i,numIndex-1) = DoFRef(A, vDoFIndex[i], vDoFIndex[numIndex-1]);
			}


			// schur complement correction


			number schur = 0.0;
			if (m_weight.invalid())
			{
				for (size_t i = 0; i < numIndex-1; i++)
				{
					schur += BlockInv(numIndex-1,i)/BlockInv(i,i)*BlockInv(i,numIndex-1);
				}
				// std::cerr << "unweighted:" << schur << " " << BlockInv(numIndex-1, numIndex-1)<< m_beta<<  "=>";

				BlockInv(numIndex-1, numIndex-1) = schur + (BlockInv(numIndex-1, numIndex-1)-schur)/m_beta;

				//std::cerr << BlockInv(numIndex-1, numIndex-1) << std::endl;

			} else {
				for (size_t i = 0; i < numIndex-1; i++)
				{
					number wi2 = DoFRef(*m_weight, vDoFIndex[i]);
					schur += wi2*BlockInv(numIndex-1,i)/BlockInv(i,i)*BlockInv(i,numIndex-1);
				}
			   //std::cerr <<  "weighted:"<< schur << " " << BlockInv(numIndex-1, numIndex-1)<< m_beta<<  "=>";

			   //  BlockInv(numIndex-1, numIndex-1) = - schur - (BlockInv(numIndex-1, numIndex-1)-schur)/m_beta; // worked with alpha=1, beta=-2.0
			   BlockInv(numIndex-1, numIndex-1) =  schur + (BlockInv(numIndex-1, numIndex-1) - schur)/m_beta;

			   //std::cerr << BlockInv(numIndex-1, numIndex-1) << std::endl;
			}


			//BlockInv(numIndex-1, numIndex-1) = - 2.0*schur;



				/*for (size_t k = 0; k < numIndex; k++)
					BlockInv(j,k) = DoFRef(A, vDoFIndex[j], vDoFIndex[k]);*/

		//	get inverse
			if(!Invert(BlockInv)){
				std::stringstream ss;
				ss << d <<"dim - Elem-Mat "<<b<<" singular. size: "<<numIndex<<"\n";
				for (size_t j = 0; j < numIndex; j++)
					ss << vDoFIndex[j][0] << ", ";
				ss << "\n";

				BlockInv = 0.0;
				for (size_t j = 0; j < numIndex; j++)
					for (size_t k = 0; k < numIndex; k++)
						BlockInv(j,k) = DoFRef(A, vDoFIndex[j], vDoFIndex[k]);

				ss << BlockInv;

				UG_THROW(ss.str());
			}
		}
	}

	m_bInit = true;
}

template <typename TDomain, typename TAlgebra>
bool ComponentGaussSeidel<TDomain, TAlgebra>::
step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d)
{
//	check that grid funtion passed
	GridFunction<TDomain, TAlgebra>* pC
					= dynamic_cast<GridFunction<TDomain, TAlgebra>*>(&c);
	if(pC == NULL)
		UG_THROW("ComponentGaussSeidel: expects correction to be a GridFunction.");

	const vector_type* pD = &d;
	const matrix_type* pMat = pOp.get();
#ifdef UG_PARALLEL
	SmartPtr<vector_type> spDtmp;
	if(pcl::NumProcs() > 1){
	//	make defect unique
		spDtmp = d.clone();
		spDtmp->change_storage_type(PST_UNIQUE);
		pD = spDtmp.get();
		pMat = &m_A;
	}
#endif

//	check if initialized
	if(!m_bInit)
		extract_blocks(*pMat, *pC);

// 	set correction to zero
	pC->set(0.0);

//	loop Grouping objects
	for(size_t i = 0; i < m_vGroupObj.size(); ++i)
	{
	//	get its dimension
		const int d = m_vGroupObj[i];

	//	get relax param
		number damp = m_relax;
		if(d < (int)m_vDamp.size()) damp *= m_vDamp[d];

	//	apply
		if (m_weight.invalid()) {
			apply_blocks(*pMat, *pC, *pD, damp, m_vDimCache[d], false);
			apply_blocks(*pMat, *pC, *pD, damp, m_vDimCache[d], true);
		} else {
			apply_blocks_weighted(*pMat, *pC, *pD, damp, m_vDimCache[d], false);
			apply_blocks_weighted(*pMat, *pC, *pD, damp, m_vDimCache[d], true);
		}

	}

#ifdef UG_PARALLEL
	 pC->set_storage_type(PST_UNIQUE);
#endif

	return true;
}


template <typename TDomain, typename TAlgebra>
SmartPtr<ILinearIterator<typename TAlgebra::vector_type> >
ComponentGaussSeidel<TDomain, TAlgebra>::clone()
{
	SmartPtr<ComponentGaussSeidel<TDomain, TAlgebra> >
					newInst(new ComponentGaussSeidel<TDomain, TAlgebra>(m_relax, m_vFullRowCmp, m_vGroupObj, m_vDamp));
	newInst->set_debug(debug_writer());
	newInst->set_damp(this->damping());
	newInst->set_alpha(this->m_alpha);
	newInst->set_beta(this->m_beta);
	newInst->set_weights(this->m_bWeighted);
	return newInst;
}


} // end namespace ug

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__COMPONENT_GAUSS_SEIDEL__ */
