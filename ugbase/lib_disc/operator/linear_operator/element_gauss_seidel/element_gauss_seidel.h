/*
 * element_gauss_seidel.h
 *
 *  Created on: 24.07.2012
 *      Author: Christian Wehner, Andreas Vogel
 *
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

		// set correction values for element indices to zero
		for(size_t j = 0; j < numIndex; ++j){
			//c[vInd[j]] = 0;
		}

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
		// note: element indices in the sum are excluded by setting c to zero
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
				UG_THROW("ElementGaussSeidel expects correction to be a GridFunction.");


			typedef typename GridFunction<TDomain, TAlgebra>::element_type Element;
			typedef typename GridFunction<TDomain, TAlgebra>::side_type Side;

#ifdef UG_PARALLEL
			if(pcl::GetNumProcesses() > 1)
			{
			//	make defect unique
				SmartPtr<vector_type> spDtmp = d.clone();
				spDtmp->change_storage_type(PST_UNIQUE);

				if		(m_type == "element") ElementGaussSeidelStep<Element,TDomain,TAlgebra>(m_A, *pC, *spDtmp, m_relax);
				else if	(m_type == "side") ElementGaussSeidelStep<Side,TDomain,TAlgebra>(m_A, *pC, *spDtmp, m_relax);
				else if	(m_type == "face") ElementGaussSeidelStep<Face,TDomain,TAlgebra>(m_A, *pC, *spDtmp, m_relax);
				else if	(m_type == "edge") ElementGaussSeidelStep<EdgeBase,TDomain,TAlgebra>(m_A, *pC, *spDtmp, m_relax);
				else if	(m_type == "vertex") ElementGaussSeidelStep<VertexBase,TDomain,TAlgebra>(m_A, *pC, *spDtmp, m_relax);
				else UG_THROW("ElementGaussSeidel: wrong patch type '"<<m_type<<"'."
					         " Options: element, side, face, edge, vertex.")

				pC->set_storage_type(PST_UNIQUE);
				return true;
			}
			else
#endif
			{
				matrix_type &mat = *pOp;
				if		(m_type == "element") ElementGaussSeidelStep<Element,TDomain,TAlgebra>(mat, *pC, d, m_relax);
				else if	(m_type == "side") ElementGaussSeidelStep<Side,TDomain,TAlgebra>(mat, *pC, d, m_relax);
				else if	(m_type == "face") ElementGaussSeidelStep<Face,TDomain,TAlgebra>(mat, *pC, d, m_relax);
				else if	(m_type == "edge") ElementGaussSeidelStep<EdgeBase,TDomain,TAlgebra>(mat, *pC, d, m_relax);
				else if	(m_type == "vertex") ElementGaussSeidelStep<VertexBase,TDomain,TAlgebra>(mat, *pC, d, m_relax);
				else UG_THROW("ElementGaussSeidel: wrong patch type '"<<m_type<<"'."
					         " Options: element, side, face, edge, vertex.")

#ifdef UG_PARALLEL
				c.set_storage_type(PST_UNIQUE);
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
