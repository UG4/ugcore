/*
 * component_gauss_seidel.h
 *
 *  Created on: 24.07.2012
 *      Author: Christian Wehner, Andreas Vogel
 *
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

template<typename TGroupObj, typename TDomain, typename TAlgebra>
void ComponentGaussSeidelStep(const typename TAlgebra::matrix_type& A,
                            GridFunction<TDomain, TAlgebra>& c,
                            const typename TAlgebra::vector_type& d,
                            number relax,
                            const std::vector<size_t>& vCmp,
							const std::vector<size_t>& vRemainCmp)
{
	// memory for local algebra
	DenseVector< VariableArray1<number> > s;
	DenseVector< VariableArray1<number> > x;
	DenseMatrix< VariableArray2<number> > mat;
	std::vector<DoFIndex> vCmpInd, vInd;
	
	typedef typename GridFunction<TDomain, TAlgebra>::element_type Element;
	std::vector<Element*> vElem;

	// loop all grouping objects
	typedef typename GridFunction<TDomain, TAlgebra>::template traits<TGroupObj>::const_iterator GroupObjIter;
	for(GroupObjIter iter = c.template begin<TGroupObj>();
					 iter != c.template end<TGroupObj>(); ++iter){

		// get grouping obj
		TGroupObj* groupObj = *iter;

		// get all dof indices on obj associated to cmps
		vCmpInd.clear();
		for(size_t f = 0; f < vCmp.size(); ++f)
			c.inner_dof_indices(groupObj, vCmp[f], vCmpInd, false);

		// check if equation present
		if(vCmpInd.empty()) continue;

		// collect elems associated to grouping object
		c.collect_associated(vElem, groupObj);

		// get all algebraic indices on element
		if(TGroupObj::dim <= VERTEX){
			std::vector<EdgeBase*> vSub;
			c.collect_associated(vSub, groupObj);
			for(size_t i = 0; i < vSub.size(); ++i)
				for(size_t f = 0; f < vCmp.size(); ++f)
					c.inner_dof_indices(vSub[i], vCmp[f], vCmpInd, false);
		}
		if(TGroupObj::dim <= EDGE){
			std::vector<Face*> vSub;
			c.collect_associated(vSub, groupObj);
			for(size_t i = 0; i < vSub.size(); ++i)
				for(size_t f = 0; f < vCmp.size(); ++f)
					c.inner_dof_indices(vSub[i], vCmp[f], vCmpInd, false);
		}

		// get all algebraic indices on element
		vInd.clear();
		for(size_t i = 0; i < vElem.size(); ++i)
			for(size_t f = 0; f < vRemainCmp.size(); ++f)
				c.dof_indices(vElem[i], vRemainCmp[f], vInd, false, false);

		// check for doublicates
		if(vElem.size() > 1){
		    std::sort(vInd.begin(), vInd.end());
		    vInd.erase(std::unique(vInd.begin(), vInd.end()), vInd.end());
		}

		// get number of indices on patch
		const size_t numIndex = vCmpInd.size() + vInd.size();

//		UG_LOG("CGS in dim: "<<TGroupObj::dim<<", #p: "<<vCmpInd.size()<<
//		       ", v: "<<vInd.size()<<" ("<<vInd.size()/2<<" each)\n")

		// concat all indices
		vInd.insert(vInd.end(), vCmpInd.begin(), vCmpInd.end());

		// fill local block matrix
		mat.resize(numIndex, numIndex);
		mat = 0.0;

		// copy matrix rows (only including cols of selected indices)
		for (size_t j = 0; j < numIndex; j++)
			for (size_t k = 0; k < numIndex; k++)
				mat(j,k) = DoFRef(A, vInd[j], vInd[k]);

		// compute s[j] := d[j] - sum_k A(j,k)*c[k]
		// note: the loop over k is the whole matrix row (not only selected indices)
		typedef typename TAlgebra::matrix_type::const_row_iterator const_row_iterator;
		const static int blockSize = TAlgebra::blockSize;
		s.resize(numIndex);
		for (size_t j = 0; j<numIndex; j++){
			typename TAlgebra::vector_type::value_type sj = d[vInd[j][0]];
			for(const_row_iterator it = A.begin_row(vInd[j][0]);
									it != A.end_row(vInd[j][0]) ; ++it){
				MatMultAdd(sj, 1.0, sj, -1.0, it.value(), c[it.index()]);
			};
			s.subassign(j*blockSize,sj);
		};

		// solve block
		x.resize(numIndex);
		InverseMatMult(x,1,mat,s);
		for (size_t j=0;j<numIndex;j++){
			DoFRef(c, vInd[j]) += relax*x[j];
		};
	}
}


///	ComponentGaussSeidel Preconditioner
template <typename TDomain, typename TAlgebra>
class ComponentGaussSeidel : public IPreconditioner<TAlgebra>
{
	public:
	/// World dimension
		typedef typename GridFunction<TDomain, TAlgebra>::element_type Element;
		typedef typename GridFunction<TDomain, TAlgebra>::side_type Side;
		static const int dim = TDomain::dim;

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
		ComponentGaussSeidel(const std::string& cmps) : m_relax(1), m_cmps(cmps){};

	///	constructor setting relaxation and type
		ComponentGaussSeidel(number relax, const std::string& cmps) : m_relax(relax), m_cmps(cmps) {};

	///	constructor setting relaxation and type
		ComponentGaussSeidel(number relax, const std::string& cmps,
		                     const std::vector<int>& vSmooth, const std::vector<number>& vDamp)
		: m_relax(relax), m_cmps(cmps), m_vSmooth(vSmooth), m_vDamp(vDamp) {};

	///	Clone
		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			SmartPtr<ComponentGaussSeidel<TDomain, TAlgebra> >
							newInst(new ComponentGaussSeidel<TDomain, TAlgebra>(m_relax, m_cmps, m_vSmooth, m_vDamp));
			newInst->set_debug(debug_writer());
			newInst->set_damp(this->damping());
			return newInst;
		}

	///	Destructor
		virtual ~ComponentGaussSeidel()
		{};

	///	returns if parallel solving is supported
		virtual bool supports_parallel() const {return true;}

	/// set relaxation parameter
		void set_relax(number omega){m_relax=omega;};

	/// set type
		void set_cmps(const std::string& cmps){m_cmps=cmps;};

	protected:
	///	Name of preconditioner
		virtual const char* name() const {return "ComponentGaussSeidel";}

	///	Preprocess routine
		virtual bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp)
		{
#ifdef UG_PARALLEL
			if(pcl::GetNumProcesses() > 1)
			{
				//	copy original matrix
				MakeConsistent(*pOp, m_A);
				//	set zero on slaves
				std::vector<IndexLayout::Element> vIndex;
				CollectUniqueElements(vIndex,  m_A.layouts()->slave());
				SetDirichletRow(m_A, vIndex);
			}
#endif
			return true;
		}

		virtual bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d)
		{
			GridFunction<TDomain, TAlgebra>* pC
							= dynamic_cast<GridFunction<TDomain, TAlgebra>*>(&c);
			if(pC == NULL)
				UG_THROW("ComponentGaussSeidel: expects correction to be a GridFunction.");

			ConstSmartPtr<DoFDistributionInfo> ddinfo =
									pC->approx_space()->dof_distribution_info();

			// tokenize passed functions
			std::vector<std::string> vCmpString = TokenizeTrimString(m_cmps);
			if(vCmpString.empty())
				UG_THROW("ComponentGaussSeidel: No components set.")

			// get ids of selected functions
			std::vector<size_t> vCmp;
			for(size_t i = 0; i < vCmpString.size(); ++i)
				vCmp.push_back(ddinfo->fct_id_by_name(vCmpString[i].c_str()));

			// compute remaining functions
			std::vector<size_t> vRemainCmp;
			for(size_t f = 0; f < ddinfo->num_fct(); ++f)
				if(std::find(vCmp.begin(), vCmp.end(), f) == vCmp.end())
					vRemainCmp.push_back(f);

			const vector_type* pD = &d;
			const matrix_type* pMat = pOp.get();

#ifdef UG_PARALLEL
			if(pcl::GetNumProcesses() > 1){
			//	make defect unique
				SmartPtr<vector_type> spDtmp = d.clone();
				spDtmp->change_storage_type(PST_UNIQUE);
				pD = spDtmp.get();
				pMat = &m_A;
			}
#endif

			// check which dimension to loop
			std::vector<bool> vLoopDim(dim+1, false);
			for(int d = 0; d <= dim; ++d)
				for(size_t f = 0; f < vCmp.size(); ++f)
					if(ddinfo->max_fct_dofs(vCmp[f], d) > 0)
						vLoopDim[d] = true;

			// set all vector entries to zero
			pC->set(0.0);

			if(m_vSmooth.empty()){

				if(vLoopDim[VERTEX]) ComponentGaussSeidelStep<VertexBase,TDomain,TAlgebra>(*pMat, *pC, *pD, m_relax, vCmp, vRemainCmp);
				if(vLoopDim[EDGE]) ComponentGaussSeidelStep<EdgeBase,TDomain,TAlgebra>(*pMat, *pC, *pD, m_relax, vCmp, vRemainCmp);
				if(vLoopDim[FACE]) ComponentGaussSeidelStep<Face,TDomain,TAlgebra>(*pMat, *pC, *pD, m_relax, vCmp, vRemainCmp);
				if(vLoopDim[VOLUME]) ComponentGaussSeidelStep<Volume,TDomain,TAlgebra>(*pMat, *pC, *pD, m_relax, vCmp, vRemainCmp);
			}
			else{
				for(size_t i = 0; i < m_vSmooth.size(); ++i){

					const int d = m_vSmooth[i];
					number damp = m_relax;

					if(d < (int)m_vDamp.size())
						damp *= m_vDamp[d];

					if(!vLoopDim[d]) continue;

					switch(d){
						case VERTEX: ComponentGaussSeidelStep<VertexBase,TDomain,TAlgebra>(*pMat, *pC, *pD, damp, vCmp, vRemainCmp); break;
						case EDGE: ComponentGaussSeidelStep<EdgeBase,TDomain,TAlgebra>(*pMat, *pC, *pD, damp, vCmp, vRemainCmp); break;
						case FACE: ComponentGaussSeidelStep<Face,TDomain,TAlgebra>(*pMat, *pC, *pD, damp, vCmp, vRemainCmp); break;
						case VOLUME: ComponentGaussSeidelStep<Volume,TDomain,TAlgebra>(*pMat, *pC, *pD, damp, vCmp, vRemainCmp); break;
						default: UG_THROW("Not implemented.")
					}
				}
			}


#ifdef UG_PARALLEL
			 pC->set_storage_type(PST_UNIQUE);
#endif
			return true;
		}

	///	Postprocess routine
		virtual bool postprocess() {return true;}

	protected:
#ifdef UG_PARALLEL
		matrix_type m_A;
#endif

	///	relaxing parameter
		number m_relax;

	///	components, whose matrix row must be fulfilled completely on gs-block
		std::string m_cmps;

	/// smooth order
		std::vector<int> m_vSmooth;
		std::vector<number> m_vDamp;
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__COMPONENT_GAUSS_SEIDEL__ */
