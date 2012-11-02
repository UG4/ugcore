/*
 * line_smoothers.h
 *
 *  Created on: 27.08.2012
 *      Author: Christian Wehner
 *
 *  Gauss-Seidel type smoothers using alternating directions
 *  for handling of anisotropic problems 
 *	Steps forward and backward in all coordinate directions
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__LINE_SMOOTHERS__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__LINE_SMOOTHERS__

#include "lib_algebra/operator/interface/operator.h"
#include "lib_algebra/operator/interface/operator_iterator.h"
#include "lib_algebra/cpu_algebra/core_smoothers.h"
#include "lib_disc/dof_manager/lexorder.h"
#include "lib_grid/algorithms/attachment_util.h"
#include "lib_disc/common/groups_util.h"
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
	#include "lib_algebra/parallelization/parallel_matrix_overlap_impl.h"
#endif

namespace ug{

// ORDER IN Y DIRECTION

// Order for 1D
template<int dim>
bool ComparePosDimYDir(const std::pair<MathVector<dim>, size_t> &p1,
					   const std::pair<MathVector<dim>, size_t> &p2)
{return false;}

template<>
bool ComparePosDimYDir<1>(const std::pair<MathVector<1>, size_t> &p1,
						  const std::pair<MathVector<1>, size_t> &p2)
{
	return (p1.first[0]<p2.first[0]);
};

// Order for 2D
template<>
bool ComparePosDimYDir<2>(const std::pair<MathVector<2>, size_t> &p1,
						  const std::pair<MathVector<2>, size_t> &p2)
{
	if (p1.first[0]==p2.first[0]) {
		return (p1.first[1] < p2.first[1]);
	}
	else {
		return (p1.first[0] < p2.first[0]);
	}
};

// Order for 3D
template<>
bool ComparePosDimYDir<3>(const std::pair<MathVector<3>, size_t> &p1,
						  const std::pair<MathVector<3>, size_t> &p2)
{
	if (p1.first[2]==p2.first[2]){
		if (p1.first[0]==p2.first[0]) {
			return (p1.first[1] < p2.first[1]);
		}
		else {
			return (p1.first[0] < p2.first[0]);
		}
	}
	else{
		return (p1.first[2] < p2.first[2]);
	}
};

template<int dim>
void ComputeDirectionYOrder(std::vector<std::pair<MathVector<dim>, size_t> >& vPos,std::vector<size_t>& indY)
{
	indY.resize(indY.size()+vPos.size());
	//  sort indices based on their position
	std::sort(vPos.begin(), vPos.end(), ComparePosDimYDir<dim>);
	//	write mapping
	for (size_t i=0; i < vPos.size(); ++i){
			indY[i] = vPos[i].second;
	};
}

template <typename TDD, typename TDomain>
void OrderDirectionYForDofDist(SmartPtr<TDD> dd,
                        ConstSmartPtr<TDomain> domain,std::vector<size_t>& indY)
{
	//	Lex Ordering is only possible in this cases:
	//	b) Same number of DoFs on each geometric object (or no DoFs on object)
	//		--> in this case we can order all dofs
	//	a) different trial spaces, but DoFs for each trial spaces only on separate
	//	   geometric objects. (e.g. one space only vertices, one space only on edges)
	//		--> in this case we can order all geometric objects separately

	//	a) check for same number of DoFs on every geometric object
		bool bEqualNumDoFOnEachGeomObj = true;
		int numDoFOnGeomObj = -1;
		for(int si = 0; si < dd->num_subsets(); ++si){
			for(int roid = 0; roid < NUM_REFERENCE_OBJECTS; ++roid){
				const int numDoF = dd->num_dofs((ReferenceObjectID)roid, si);

				if(numDoF == 0) continue;

				if(numDoFOnGeomObj == -1){
					numDoFOnGeomObj = numDoF;
				}
				else{
					if(numDoFOnGeomObj != numDoF)
						bEqualNumDoFOnEachGeomObj = false;
				}
			}
		}

	//	b) check for non-mixed spaces
		std::vector<bool> bSingleSpaceUsage(NUM_REFERENCE_OBJECTS, true);
		std::vector<bool> vHasDoFs(NUM_REFERENCE_OBJECTS, false);
		for(size_t fct = 0; fct < dd->num_fct(); ++fct){

			LFEID lfeID = dd->local_finite_element_id(fct);
			const int lfeDim = dd->dim(fct);
			const CommonLocalDoFSet& locDoF = LocalDoFSetProvider::get(lfeDim, lfeID);

			for(int roid = 0; roid < NUM_REFERENCE_OBJECTS; ++roid){
				const int numDoF = locDoF.num_dof((ReferenceObjectID)roid);

				if(numDoF <= 0) continue;

				if(vHasDoFs[roid] == false){
					vHasDoFs[roid] = true;
				}
				else{
					bSingleSpaceUsage[roid] = false;
				}
			}
		}
		std::vector<bool> bSortableComp(dd->num_fct(), true);
		for(size_t fct = 0; fct < dd->num_fct(); ++fct){

			LFEID lfeID = dd->local_finite_element_id(fct);
			const int lfeDim = dd->dim(fct);
			const CommonLocalDoFSet& locDoF = LocalDoFSetProvider::get(lfeDim, lfeID);

			for(int roid = 0; roid < NUM_REFERENCE_OBJECTS; ++roid){
				if(locDoF.num_dof((ReferenceObjectID)roid) != 0){
					if(bSingleSpaceUsage[roid] == false)
						bSortableComp[fct] = false;
				}
			}
		}

	//	get position attachment
		typedef typename std::pair<MathVector<TDomain::dim>, size_t> pos_type;

	//	get positions of indices
		std::vector<pos_type> vPositions;

	//	a) we can order globally
		if(bEqualNumDoFOnEachGeomObj)
		{
			ExtractPositions<TDomain, TDD>(domain, dd, vPositions);

			ComputeDirectionYOrder<TDomain::dim>(vPositions, indY);
		}
	//	b) we can only order some spaces
		else
		{
//			UG_LOG("OrderLex: Cannot order globally, trying to order some components:\n");
			for(size_t fct = 0; fct < dd->num_fct(); ++fct){
				if(bSortableComp[fct] == false){
//					UG_LOG("OrderLex: "<<dd->name(fct)<<" NOT SORTED.\n");
					continue;
				}

				ExtractPositions<TDomain, TDD>(domain, dd, fct, vPositions);

				ComputeDirectionYOrder<TDomain::dim>(vPositions, indY);
			}
		}
};


// ORDER IN Z DIRECTION

// Order for 1D
template<int dim>
bool ComparePosDimZDir(const std::pair<MathVector<dim>, size_t> &p1,
					   const std::pair<MathVector<dim>, size_t> &p2)
{return false;}


template<>
bool ComparePosDimZDir<1>(const std::pair<MathVector<1>, size_t> &p1,
						  const std::pair<MathVector<1>, size_t> &p2)
{
	return (p1.first[0]<p2.first[0]);
};

// Order for 2D
template<>
bool ComparePosDimZDir<2>(const std::pair<MathVector<2>, size_t> &p1,
						  const std::pair<MathVector<2>, size_t> &p2)
{
	if (p1.first[0]==p2.first[0]) {
		return (p1.first[1] < p2.first[1]);
	}
	else {
		return (p1.first[0] < p2.first[0]);
	}
};

// Order for 3D
template<>
bool ComparePosDimZDir<3>(const std::pair<MathVector<3>, size_t> &p1,
						  const std::pair<MathVector<3>, size_t> &p2)
{
	if (p1.first[1]==p2.first[1]){
		if (p1.first[0]==p2.first[0]) {
			return (p1.first[2] < p2.first[2]);
		}
		else {
			return (p1.first[0] < p2.first[0]);
		}
	}
	else{
		return (p1.first[1] < p2.first[1]);
	}
};

template<int dim>
void ComputeDirectionZOrder(std::vector<std::pair<MathVector<dim>, size_t> >& vPos,std::vector<size_t>& indZ)
{
	indZ.resize(indZ.size()+vPos.size());
	std::sort(vPos.begin(), vPos.end(), ComparePosDimZDir<dim>);
	for (size_t i=0; i < vPos.size(); ++i){
		indZ[i] = vPos[i].second;
	};
}

template <typename TDD, typename TDomain>
void OrderDirectionZForDofDist(SmartPtr<TDD> dd,
                        ConstSmartPtr<TDomain> domain,std::vector<size_t>& indZ)
{
	//	Lex Ordering is only possible in this cases:
	//	b) Same number of DoFs on each geometric object (or no DoFs on object)
	//		--> in this case we can order all dofs
	//	a) different trial spaces, but DoFs for each trial spaces only on separate
	//	   geometric objects. (e.g. one space only vertices, one space only on edges)
	//		--> in this case we can order all geometric objects separately

	//	a) check for same number of DoFs on every geometric object
		bool bEqualNumDoFOnEachGeomObj = true;
		int numDoFOnGeomObj = -1;
		for(int si = 0; si < dd->num_subsets(); ++si){
			for(int roid = 0; roid < NUM_REFERENCE_OBJECTS; ++roid){
				const int numDoF = dd->num_dofs((ReferenceObjectID)roid, si);

				if(numDoF == 0) continue;

				if(numDoFOnGeomObj == -1){
					numDoFOnGeomObj = numDoF;
				}
				else{
					if(numDoFOnGeomObj != numDoF)
						bEqualNumDoFOnEachGeomObj = false;
				}
			}
		}

	//	b) check for non-mixed spaces
		std::vector<bool> bSingleSpaceUsage(NUM_REFERENCE_OBJECTS, true);
		std::vector<bool> vHasDoFs(NUM_REFERENCE_OBJECTS, false);
		for(size_t fct = 0; fct < dd->num_fct(); ++fct){

			LFEID lfeID = dd->local_finite_element_id(fct);
			const int lfeDim = dd->dim(fct);
			const CommonLocalDoFSet& locDoF = LocalDoFSetProvider::get(lfeDim, lfeID);

			for(int roid = 0; roid < NUM_REFERENCE_OBJECTS; ++roid){
				const int numDoF = locDoF.num_dof((ReferenceObjectID)roid);

				if(numDoF <= 0) continue;

				if(vHasDoFs[roid] == false){
					vHasDoFs[roid] = true;
				}
				else{
					bSingleSpaceUsage[roid] = false;
				}
			}
		}
		std::vector<bool> bSortableComp(dd->num_fct(), true);
		for(size_t fct = 0; fct < dd->num_fct(); ++fct){

			LFEID lfeID = dd->local_finite_element_id(fct);
			const int lfeDim = dd->dim(fct);
			const CommonLocalDoFSet& locDoF = LocalDoFSetProvider::get(lfeDim, lfeID);

			for(int roid = 0; roid < NUM_REFERENCE_OBJECTS; ++roid){
				if(locDoF.num_dof((ReferenceObjectID)roid) != 0){
					if(bSingleSpaceUsage[roid] == false)
						bSortableComp[fct] = false;
				}
			}
		}

	//	get position attachment
		typedef typename std::pair<MathVector<TDomain::dim>, size_t> pos_type;

	//	get positions of indices
		std::vector<pos_type> vPositions;

	//	a) we can order globally
		if(bEqualNumDoFOnEachGeomObj)
		{
			ExtractPositions<TDomain, TDD>(domain, dd, vPositions);

		//	get mapping: old -> new index
			ComputeDirectionZOrder<TDomain::dim>(vPositions, indZ);
		}
	//	b) we can only order some spaces
		else
		{
//			UG_LOG("OrderLex: Cannot order globally, trying to order some components:\n");
			for(size_t fct = 0; fct < dd->num_fct(); ++fct){
				if(bSortableComp[fct] == false){
					// UG_LOG("OrderLex: "<<dd->name(fct)<<" NOT SORTED.\n");
					continue;
				}

				ExtractPositions<TDomain, TDD>(domain, dd, fct, vPositions);

			//	get mapping: old -> new index
				ComputeDirectionZOrder<TDomain::dim>(vPositions, indZ);
			}
		}
};


template <typename TDomain, typename TDD, typename TBaseElem>
void collectStretchedElementIndices(ConstSmartPtr<TDomain> domain,
                          ConstSmartPtr<TDD> dd,
                          std::vector<size_t>& indarray,number alpha){
	typename TDD::template traits<TBaseElem>::const_iterator iter, iterEnd;
//	vector for positions
	std::vector<MathVector<TDomain::dim> > vElemPos;

//	algebra indices vector
	std::vector<MultiIndex<2> > ind;
	
	const typename TDomain::position_accessor_type& aaPos = domain->position_accessor();

//	loop all subsets
	for(int si = 0; si < dd->num_subsets(); ++si)
	{
	//	get iterators
		iter = dd->template begin<TBaseElem>(si);
		iterEnd = dd->template end<TBaseElem>(si);

	//	loop all elements
		for(;iter != iterEnd; ++iter)
		{
		//	get vertex
			TBaseElem* elem = *iter;

		//	loop all functions
			for(size_t fct = 0; fct < dd->num_fct(); ++fct)
			{
				//	skip non-used function
				if(!dd->is_def_in_subset(fct,si)) continue;
				std::vector<EdgeBase*> vEdge;
				CollectEdgesSorted(vEdge, domain->grid, elem);
				std::vector<number> edgeLength(vEdge.size());
				std::vector<MultiIndex<2> > ind;
				MathVector<TDomain::dim> vCoord[2];
				for (size_t i=0;i<vEdge.size();i++){
					for (size_t j=0;j<2;j++) vCoord[i] = aaPos[vEdge[i]->vertex(i)];
					edgeLength[i] = VecDistance(vCoord[0],vCoord[1]);
				};
				number minedgelength=edgeLength[0];
				number maxedgelength=edgeLength[1];
				for (size_t i=1;i<vEdge.size();i++){
					if (edgeLength[i]<minedgelength) minedgelength=edgeLength[i];
					if (edgeLength[i]>maxedgelength) maxedgelength=edgeLength[i];
				};
				if (maxedgelength/minedgelength>alpha){
					//	load indices associated with element function
					dd->inner_multi_indices(elem, fct, ind);

					//	load positions associated with element and function
					InnerDoFPosition(vElemPos, elem, *(const_cast<TDomain*>(domain.get())),
				                 dd->local_finite_element_id(fct), dd->dim(fct));

					//	check correct size
					UG_ASSERT(ind.size() == vElemPos.size(), "Num MultiIndex ("<<ind.size()
						  <<") and Num Position ("<<vElemPos.size()<<") must match."
						 "GeomObject dim="<<geometry_traits<TBaseElem>::BASE_OBJECT_ID);

					//	write position
					for(size_t sh = 0; sh < ind.size(); ++sh)
					{
						const size_t index = ind[sh][0];
						indarray.push_back(index);
					}
				};
			};
		};
	};
};



template <typename TDomain,typename TAlgebra>
class LineGaussSeidel : public IPreconditioner<TAlgebra>
{
	public:
	///	Domain
		typedef TDomain domain_type;
		
	//	Algebra type
		typedef TAlgebra algebra_type;

	//	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	//	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Matrix Operator type
		typedef typename IPreconditioner<TAlgebra>::matrix_operator_type matrix_operator_type;

	///	Base type
		typedef IPreconditioner<TAlgebra> base_type;

	protected:
		using base_type::set_debug;
		using base_type::debug_writer;
		using base_type::write_debug;
		
	protected:
		//	approximation space for level and surface grid
		SmartPtr<ApproximationSpace<TDomain> > m_spApproxSpace;

		// index sets
		std::vector<size_t> indY;
		std::vector<size_t> indZ;
		size_t m_ind_end;
		
		size_t m_nr_forwardx;
		size_t m_nr_backwardx;
		size_t m_nr_forwardy;
		size_t m_nr_backwardy;
		size_t m_nr_forwardz;
		size_t m_nr_backwardz;

		///	world dimension
		static const int dim = domain_type::dim;

		bool m_init;

	public:
	//	Constructor
		LineGaussSeidel(SmartPtr<ApproximationSpace<TDomain> > approxSpace){
			m_spApproxSpace = approxSpace;
			m_init = false;
			m_nr_forwardx=1;
			m_nr_backwardx=1;
			m_nr_forwardy=1;
			m_nr_backwardy=1;
			m_nr_forwardz=1;
			m_nr_backwardz=1;
			OrderLex<TDomain>(*m_spApproxSpace,"lr");
		};

	// 	Update ordering indices
		bool update(size_t xsize){
			indY.resize(0);
			indZ.resize(0);
			// lexicographic ordering of unknowns
			if (m_nr_forwardx+m_nr_backwardx>0){
//				OrderLex<TDomain>(*m_spApproxSpace,"lr");
			}
			if (m_nr_forwardy+m_nr_backwardy+m_nr_forwardz+m_nr_backwardz>0){
				size_t level=3289578756;
				for (size_t i=0;i<m_spApproxSpace->num_levels();i++){
					if (m_spApproxSpace->level_dof_distribution(i)->num_indices()==xsize){
						level = i;
						break;
					};
				};
				if (level==3289578756){
					return false;
				}
				if ((dim>1)&&(m_nr_forwardy+m_nr_backwardy>0)){
					OrderDirectionYForDofDist<LevelDoFDistribution, TDomain>(m_spApproxSpace->level_dof_distribution(level), m_spApproxSpace->domain(),indY);
				}
				if ((dim>2)&&(m_nr_forwardz+m_nr_backwardz>0)){
					OrderDirectionZForDofDist<LevelDoFDistribution, TDomain>(m_spApproxSpace->level_dof_distribution(level), m_spApproxSpace->domain(),indZ);
				}
			};
			m_ind_end = indY.size();
			m_init = true;
			return true;
		}
		
	//  Destructor
		~LineGaussSeidel(){
			indY.clear();
			indZ.clear();
		};
		
		void set_num_steps(size_t forwardx,size_t backwardx,size_t forwardy,size_t backwardy,size_t forwardz,size_t backwardz){
			m_nr_forwardx=forwardx;
			m_nr_backwardx=backwardx;
			m_nr_forwardy=forwardy;
			m_nr_backwardy=backwardy;
			m_nr_forwardz=forwardz;
			m_nr_backwardz=backwardz;
		};

		void set_num_steps(size_t forwardx,size_t backwardx,size_t forwardy,size_t backwardy){
			m_nr_forwardx=forwardx;
			m_nr_backwardx=backwardx;
			m_nr_forwardy=forwardy;
			m_nr_backwardy=backwardy;
			m_nr_forwardz=0;
			m_nr_backwardz=0;
		};

		void set_num_steps(size_t forwardx,size_t backwardx){
			m_nr_forwardx=forwardx;
			m_nr_backwardx=backwardx;
			m_nr_forwardy=0;
			m_nr_backwardy=0;
			m_nr_forwardz=0;
			m_nr_backwardz=0;
		};

	// 	Clone
		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			SmartPtr<LineGaussSeidel<domain_type,algebra_type> > newInst(new LineGaussSeidel<domain_type,algebra_type>(m_spApproxSpace));
			newInst->set_debug(debug_writer());
			newInst->set_damp(this->damping());
			newInst->set_num_steps(this->get_num_forwardx(),this->get_num_backwardx(),this->get_num_forwardy(),this->get_num_backwardy(),
								   this->get_num_forwardz(),this->get_num_backwardz());
			return newInst;
		}
		
	public:
		size_t get_num_forwardx(){return m_nr_forwardx;}
		size_t get_num_backwardx(){return m_nr_backwardx;}
		size_t get_num_forwardy(){return m_nr_forwardy;}
		size_t get_num_backwardy(){return m_nr_backwardy;}
		size_t get_num_forwardz(){return m_nr_forwardz;}
		size_t get_num_backwardz(){return m_nr_backwardz;}

	protected:
	//	Name of preconditioner
		virtual const char* name() const {return "Line Gauss-Seidel";}

	//	Preprocess routine
		virtual bool preprocess(matrix_operator_type& mat)
		{
#ifdef UG_PARALLEL
			if(pcl::GetNumProcesses() > 1)
			{
				//	copy original matrix
				MakeConsistent(mat, m_A);
				//	set zero on slaves
				std::vector<IndexLayout::Element> vIndex;
				CollectUniqueElements(vIndex,  m_A.slave_layout());
				SetDirichletRow(m_A, vIndex);
			}
#endif
			return true;
		}

	//	Postprocess routine
		virtual bool postprocess() {return true;}

	//	Stepping routine
		virtual bool step(matrix_operator_type& mat, vector_type& c, const vector_type& d)
		{
#ifdef UG_PARALLEL
			if(pcl::GetNumProcesses() > 1)
			{
				//	make defect unique
				// todo: change that copying
				vector_type dhelp;
				dhelp.resize(d.size()); dhelp = d;
				dhelp.change_storage_type(PST_UNIQUE);
				if (!linegs_step(m_A, c, dhelp)) return false;
//				if(!gs_step_LL(m_A, c, dhelp)) return false;
				c.set_storage_type(PST_UNIQUE);
				return true;
			}
			else
#endif
			{
				if(!linegs_step(mat, c, d)) return false;
			//	if(!gs_step_LL(mat, c, d)) return false;
#ifdef UG_PARALLEL
				c.set_storage_type(PST_UNIQUE);
#endif
				return true;
			}
		}

	protected:
#ifdef UG_PARALLEL
		matrix_type m_A;
#endif

	bool linegs_step(const matrix_type &A, vector_type &x, const vector_type &b)
	{
		// gs LL has preconditioning matrix
		// (D-L)^{-1}
		if (m_init==false) update(x.size());

		size_t i;

		typename vector_type::value_type s;
		
		// forward in x direction
		if (m_nr_forwardx>0){
			
			for(i=0; i < x.size(); i++)
			{
				s = b[i];
			
				for(typename matrix_type::const_row_iterator it = A.begin_row(i); it != A.end_row(i) && it.index() < i; ++it)
					// s -= it.value() * x[it.index()];
					MatMultAdd(s, 1.0, s, -1.0, it.value(), x[it.index()]);
			
				// x[i] = s/A(i,i)
				InverseMatMult(x[i], 1.0, A(i,i), s);
			}
			for (size_t count=1;count<m_nr_forwardx;count++){
				for(i=0; i < x.size(); i++)
				{
					s = b[i];
					for(typename matrix_type::const_row_iterator it = A.begin_row(i); it != A.end_row(i); ++it){
						if (it.index()==i) continue;
						// s -= it.value() * x[it.index()];
						MatMultAdd(s, 1.0, s, -1.0, it.value(), x[it.index()]);
					}
					// x[i] = s/A(i,i)
					InverseMatMult(x[i], 1.0, A(i,i), s);
				}
			};
		};

		// backward in x direction
		for (size_t count=0;count<m_nr_backwardx;count++){
			for	(i=x.size()-1; (int)i>= 0; i--)
			{
				s = b[i];

				for(typename matrix_type::const_row_iterator it = A.begin_row(i); it != A.end_row(i); ++it){
					if (it.index()==i) continue;
					// s -= it.value() * x[it.index()];
					MatMultAdd(s, 1.0, s, -1.0, it.value(), x[it.index()]);
				}
				// x[i] = s/A(i,i)
				InverseMatMult(x[i], 1.0, A(i,i), s);
				if (i==0) break;
			};
		};

		if (dim==1) return true;
		
		if (m_nr_forwardy+m_nr_backwardy+m_nr_forwardz+m_nr_backwardz==0) return true;
		
		// forward in y direction
		for (size_t count=0;count<m_nr_forwardy;count++){
		for (size_t j=0;j < m_ind_end; j++){
			i = indY[j];

			s = b[i];

			for(typename matrix_type::const_row_iterator it = A.begin_row(i); it != A.end_row(i); ++it){
				if (it.index()==i) continue;
				// s -= it.value() * x[it.index()];
				MatMultAdd(s, 1.0, s, -1.0, it.value(), x[it.index()]);
			}
			// x[i] = s/A(i,i)
			InverseMatMult(x[i], 1.0, A(i,i), s);
		}
		};

		// backward in y direction
		for (size_t count=0;count<m_nr_backwardy;count++){
		for (size_t j=m_ind_end-1;(int)j >= 0; j--){
			i = indY[j];

			s = b[i];

			for(typename matrix_type::const_row_iterator it = A.begin_row(i); it != A.end_row(i); ++it){
				if (it.index()==i) continue;
				// s -= it.value() * x[it.index()];
				MatMultAdd(s, 1.0, s, -1.0, it.value(), x[it.index()]);
			}
			// x[i] = s/A(i,i)
			InverseMatMult(x[i], 1.0, A(i,i), s);
			if (j==0) break;
		}
		};
		if (dim==2) return true;

		// forward in z direction
		for (size_t count=0;count<m_nr_forwardz;count++){
		for (size_t j=0;j < m_ind_end; j++){
			i = indZ[j];

			s = b[i];

			for(typename matrix_type::const_row_iterator it = A.begin_row(i); it != A.end_row(i) ; ++it){
				if (it.index()==i) continue;
				// s -= it.value() * x[it.index()];
				MatMultAdd(s, 1.0, s, -1.0, it.value(), x[it.index()]);
			};
			// x[i] = s/A(i,i)
			InverseMatMult(x[i], 1.0, A(i,i), s);
		}
		}

		// backward in z direction
		for (size_t count=0;count<m_nr_backwardz;count++){
		for (size_t j=m_ind_end-1;(int)j >= 0; j--){
			i = indZ[j];

			s = b[i];

			for(typename matrix_type::const_row_iterator it = A.begin_row(i); it != A.end_row(i) ; ++it){
				if (it.index()==i) continue;
				// s -= it.value() * x[it.index()];
				MatMultAdd(s, 1.0, s, -1.0, it.value(), x[it.index()]);
			};
			// x[i] = s/A(i,i)
			InverseMatMult(x[i], 1.0, A(i,i), s);
			if (j==0) break;
		}
		}
		return true;
	}

};

template <typename TDomain,typename TAlgebra>
class LineVanka : public IPreconditioner<TAlgebra>
{
	public:
	///	Domain
		typedef TDomain domain_type;
		
	//	Algebra type
		typedef TAlgebra algebra_type;

	//	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	//	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Matrix Operator type
		typedef typename IPreconditioner<TAlgebra>::matrix_operator_type matrix_operator_type;

	///	Base type
		typedef IPreconditioner<TAlgebra> base_type;

	protected:
		using base_type::set_debug;
		using base_type::debug_writer;
		using base_type::write_debug;
		
	protected:
		//	approximation space for level and surface grid
		SmartPtr<ApproximationSpace<TDomain> > m_spApproxSpace;

		// index sets
		std::vector<size_t> indY;
		std::vector<size_t> indZ;
		size_t m_ind_end;
		
		size_t m_nr_forwardx;
		size_t m_nr_backwardx;
		size_t m_nr_forwardy;
		size_t m_nr_backwardy;
		size_t m_nr_forwardz;
		size_t m_nr_backwardz;

		///	world dimension
		static const int dim = domain_type::dim;

		bool m_init;
	
	public:
		/// set m_relaxation parameter
		void set_relax(number omega){m_relax=omega;};
		number relax(){return m_relax;};

	protected:
		number m_relax;

	public:
	//	Constructor
		LineVanka(SmartPtr<ApproximationSpace<TDomain> > approxSpace){
			m_spApproxSpace = approxSpace;
			m_init = false;
			m_nr_forwardx=1;
			m_nr_backwardx=1;
			m_nr_forwardy=1;
			m_nr_backwardy=1;
			m_nr_forwardz=1;
			m_nr_backwardz=1;
			m_relax=1;
//			OrderLex<TDomain>(*m_spApproxSpace,"lr");
		};

	// 	Update ordering indices
		bool update(size_t xsize){
			indY.resize(0);
			indZ.resize(0);
			// lexicographic ordering of unknowns
			if (m_nr_forwardx+m_nr_backwardx>0){
//				OrderLex<TDomain>(*m_spApproxSpace,"lr");
			}
			if (m_nr_forwardy+m_nr_backwardy+m_nr_forwardz+m_nr_backwardz>0){
				size_t level=3289578756;
				for (size_t i=0;i<m_spApproxSpace->num_levels();i++){
					if (m_spApproxSpace->level_dof_distribution(i)->num_indices()==xsize){
						level = i;
						break;
					};
				};
				if (level==3289578756){
					return false;
				}
				if ((dim>1)&&(m_nr_forwardy+m_nr_backwardy>0)){
					OrderDirectionYForDofDist<LevelDoFDistribution, TDomain>(m_spApproxSpace->level_dof_distribution(level), m_spApproxSpace->domain(),indY);
				}
				if ((dim>2)&&(m_nr_forwardz+m_nr_backwardz>0)){
					OrderDirectionZForDofDist<LevelDoFDistribution, TDomain>(m_spApproxSpace->level_dof_distribution(level), m_spApproxSpace->domain(),indZ);
				}
			};
			m_ind_end = indY.size();
			m_init = true;
			return true;
		}
				
	//  Destructor
		~LineVanka(){
			indY.clear();
			indZ.clear();
		};

		void set_num_steps(size_t forwardx,size_t backwardx,size_t forwardy,size_t backwardy,size_t forwardz,size_t backwardz){
			m_nr_forwardx=forwardx;
			m_nr_backwardx=backwardx;
			m_nr_forwardy=forwardy;
			m_nr_backwardy=backwardy;
			m_nr_forwardz=forwardz;
			m_nr_backwardz=backwardz;
		};

		void set_num_steps(size_t forwardx,size_t backwardx,size_t forwardy,size_t backwardy){
			m_nr_forwardx=forwardx;
			m_nr_backwardx=backwardx;
			m_nr_forwardy=forwardy;
			m_nr_backwardy=backwardy;
			m_nr_forwardz=0;
			m_nr_backwardz=0;
		};

		void set_num_steps(size_t forwardx,size_t backwardx){
			m_nr_forwardx=forwardx;
			m_nr_backwardx=backwardx;
			m_nr_forwardy=0;
			m_nr_backwardy=0;
			m_nr_forwardz=0;
			m_nr_backwardz=0;
		};

	// 	Clone
		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			SmartPtr<LineVanka<domain_type,algebra_type> > newInst(new LineVanka<domain_type,algebra_type>(m_spApproxSpace));
			newInst->set_debug(debug_writer());
			newInst->set_damp(this->damping());
			newInst->set_num_steps(this->get_num_forwardx(),this->get_num_backwardx(),this->get_num_forwardy(),this->get_num_backwardy(),
								   this->get_num_forwardz(),this->get_num_backwardz());
			newInst->set_relax(this->relax());
			return newInst;
		}
		
	public:
		size_t get_num_forwardx(){return m_nr_forwardx;}
		size_t get_num_backwardx(){return m_nr_backwardx;}
		size_t get_num_forwardy(){return m_nr_forwardy;}
		size_t get_num_backwardy(){return m_nr_backwardy;}
		size_t get_num_forwardz(){return m_nr_forwardz;}
		size_t get_num_backwardz(){return m_nr_backwardz;}

	protected:
	//	Name of preconditioner
		virtual const char* name() const {return "Line Vanka";}

	//	Preprocess routine
		virtual bool preprocess(matrix_operator_type& mat)
		{
#ifdef UG_PARALLEL
			if(pcl::GetNumProcesses() > 1)
			{
				//	copy original matrix
				MakeConsistent(mat, m_A);
				//	set zero on slaves
				std::vector<IndexLayout::Element> vIndex;
				CollectUniqueElements(vIndex,  m_A.slave_layout());
				SetDirichletRow(m_A, vIndex);
			}
#endif
			return true;
		}

	//	Postprocess routine
		virtual bool postprocess() {return true;}

	//	Stepping routine
		virtual bool step(matrix_operator_type& mat, vector_type& c, const vector_type& d)
		{
#ifdef UG_PARALLEL
			if(pcl::GetNumProcesses() > 1)
			{
				//	make defect unique
				// todo: change that copying
				vector_type dhelp;
				dhelp.resize(d.size()); dhelp = d;
				dhelp.change_storage_type(PST_UNIQUE);
				if (!linevanka_step(m_A, c, dhelp)) return false;
//				if(!gs_step_LL(m_A, c, dhelp)) return false;
				c.set_storage_type(PST_UNIQUE);
				return true;
			}
			else
#endif
			{
				if(!linevanka_step(mat, c, d)) return false;
			//	if(!gs_step_LL(mat, c, d)) return false;
#ifdef UG_PARALLEL
				c.set_storage_type(PST_UNIQUE);
#endif
				return true;
			}
		}

	protected:
#ifdef UG_PARALLEL
		matrix_type m_A;
#endif

	bool linevanka_step(const matrix_type &A, vector_type &x, const vector_type &b)
	{
		DenseVector< VariableArray1<number> > s;
		DenseVector< VariableArray1<number> > localx;
		DenseMatrix< VariableArray2<number> > mat;
		
		static const int MAXBLOCKSIZE = 19;
		size_t blockind[MAXBLOCKSIZE];
		size_t blocksize;
		// gs LL has preconditioning matrix
		// (D-L)^{-1}
		if (m_init==false) update(x.size());

		size_t i;
		size_t nrofelements=0;
		
		for(i=0; i < x.size(); i++)
		{
			x[i]=0;
			if (A(i,i)==0) nrofelements++;
		};
		
		// forward in x direction
		for (size_t count=0;count<m_nr_forwardx;count++){
			for(i=0; i < x.size(); i++){
				if (A(i,i)!=0) continue;
				blocksize=0;
				for(typename matrix_type::const_row_iterator it = A.begin_row(i); it != A.end_row(i) ; ++it){
					blockind[blocksize] = it.index();
					x[it.index()] = 0;
					blocksize++;
				};
				mat.resize(blocksize,blocksize);
				s.resize(blocksize);
				localx.resize(blocksize);
				for (size_t j=0;j<blocksize;j++){
					// fill local block matrix
					for (size_t k=0;k<blocksize;k++){
						mat.subassign(j,k,A(blockind[j],blockind[k]));
					};
					// compute rhs
					typename vector_type::value_type sj = b[blockind[j]];
					for(typename matrix_type::const_row_iterator it = A.begin_row(blockind[j]); it != A.end_row(blockind[j]) ; ++it){
						MatMultAdd(sj, 1.0, sj, -1.0, it.value(), x[it.index()]);
					};
					s.subassign(j,sj);
				};
				// solve block
				InverseMatMult(localx,1,mat,s);
				for (size_t j=0;j<blocksize;j++){
					x[blockind[j]] = m_relax*localx[j];
				};
			};
		};
		// backward in x direction
		for (size_t count=0;count<m_nr_backwardx;count++){
			for	(i=x.size()-1;(int)i>= 0; i--)
			{
				if (A(i,i)==0){ 
				blocksize=0;
				for(typename matrix_type::const_row_iterator it = A.begin_row(i); it != A.end_row(i) ; ++it){
					blockind[blocksize] = it.index();
					x[it.index()] = 0;
					blocksize++;
				};
				mat.resize(blocksize,blocksize);
				s.resize(blocksize);
				localx.resize(blocksize);
				for (size_t j=0;j<blocksize;j++){
					// fill local block matrix
					for (size_t k=0;k<blocksize;k++){
						mat.subassign(j,k,A(blockind[j],blockind[k]));
					};
					// compute rhs
					typename vector_type::value_type sj = b[blockind[j]];
					for(typename matrix_type::const_row_iterator it = A.begin_row(blockind[j]); it != A.end_row(blockind[j]) ; ++it){
						MatMultAdd(sj, 1.0, sj, -1.0, it.value(), x[it.index()]);
					};
					s.subassign(j,sj);
				};
				// solve block
				InverseMatMult(localx,1,mat,s);
				for (size_t j=0;j<blocksize;j++){
					x[blockind[j]] = m_relax*localx[j];
				};
				};
				if (i==0) break;
			};
		};

		if (dim==1) return true;
		
		if (m_nr_forwardy+m_nr_backwardy+m_nr_forwardz+m_nr_backwardz==0) return true;
		
		// forward in y direction
		for (size_t count=0;count<m_nr_forwardy;count++){
		for (size_t sortedi=0;sortedi < m_ind_end; sortedi++){
			i = indY[sortedi];
				blocksize=0;
				for(typename matrix_type::const_row_iterator it = A.begin_row(i); it != A.end_row(i) ; ++it){
					blockind[blocksize] = it.index();
					x[it.index()] = 0;
					blocksize++;
				};
				mat.resize(blocksize,blocksize);
				s.resize(blocksize);
				localx.resize(blocksize);
				for (size_t j=0;j<blocksize;j++){
					// fill local block matrix
					for (size_t k=0;k<blocksize;k++){
						mat.subassign(j,k,A(blockind[j],blockind[k]));
					};
					// compute rhs
					typename vector_type::value_type sj = b[blockind[j]];
					for(typename matrix_type::const_row_iterator it = A.begin_row(blockind[j]); it != A.end_row(blockind[j]) ; ++it){
						MatMultAdd(sj, 1.0, sj, -1.0, it.value(), x[it.index()]);
					};
					s.subassign(j,sj);
				};
				// solve block
				InverseMatMult(localx,1,mat,s);
				for (size_t j=0;j<blocksize;j++){
					x[blockind[j]] = m_relax*localx[j];
				};
		}
		};

		// backward in y direction
		for (size_t count=0;count<m_nr_backwardy;count++){
		for (size_t sortedi=m_ind_end-1;(int)sortedi >= 0; sortedi--){
			i = indY[sortedi];
				blocksize=0;
				for(typename matrix_type::const_row_iterator it = A.begin_row(i); it != A.end_row(i) ; ++it){
					blockind[blocksize] = it.index();
					x[it.index()] = 0;
					blocksize++;
				};
				mat.resize(blocksize,blocksize);
				s.resize(blocksize);
				localx.resize(blocksize);
				for (size_t j=0;j<blocksize;j++){
					// fill local block matrix
					for (size_t k=0;k<blocksize;k++){
						mat.subassign(j,k,A(blockind[j],blockind[k]));
					};
					// compute rhs
					typename vector_type::value_type sj = b[blockind[j]];
					for(typename matrix_type::const_row_iterator it = A.begin_row(blockind[j]); it != A.end_row(blockind[j]) ; ++it){
						MatMultAdd(sj, 1.0, sj, -1.0, it.value(), x[it.index()]);
					};
					s.subassign(j,sj);
				};
				// solve block
				InverseMatMult(localx,1,mat,s);
				for (size_t j=0;j<blocksize;j++){
					x[blockind[j]] = m_relax*localx[j];
				};
			if (sortedi==0) break;
		}
		};
		if (dim==2) return true;

		// forward in z direction
		for (size_t count=0;count<m_nr_forwardz;count++){
		for (size_t sortedi=0;sortedi < m_ind_end; sortedi++){
				i = indZ[sortedi];
				blocksize=0;
				for(typename matrix_type::const_row_iterator it = A.begin_row(i); it != A.end_row(i) ; ++it){
					blockind[blocksize] = it.index();
					x[it.index()] = 0;
					blocksize++;
				};
				mat.resize(blocksize,blocksize);
				s.resize(blocksize);
				localx.resize(blocksize);
				for (size_t j=0;j<blocksize;j++){
					// fill local block matrix
					for (size_t k=0;k<blocksize;k++){
						mat.subassign(j,k,A(blockind[j],blockind[k]));
					};
					// compute rhs
					typename vector_type::value_type sj = b[blockind[j]];
					for(typename matrix_type::const_row_iterator it = A.begin_row(blockind[j]); it != A.end_row(blockind[j]) ; ++it){
						MatMultAdd(sj, 1.0, sj, -1.0, it.value(), x[it.index()]);
					};
					s.subassign(j,sj);
				};
				// solve block
				InverseMatMult(localx,1,mat,s);
				for (size_t j=0;j<blocksize;j++){
					x[blockind[j]] = m_relax*localx[j];
				};
		}
		}

		// backward in z direction
		for (size_t count=0;count<m_nr_backwardz;count++){
		for (size_t sortedi=m_ind_end-1;(int)sortedi >= 0; sortedi--){
			i = indZ[sortedi];
				blocksize=0;
				for(typename matrix_type::const_row_iterator it = A.begin_row(i); it != A.end_row(i) ; ++it){
					blockind[blocksize] = it.index();
					x[it.index()] = 0;
					blocksize++;
				};
				mat.resize(blocksize,blocksize);
				s.resize(blocksize);
				localx.resize(blocksize);
				for (size_t j=0;j<blocksize;j++){
					// fill local block matrix
					for (size_t k=0;k<blocksize;k++){
						mat.subassign(j,k,A(blockind[j],blockind[k]));
					};
					// compute rhs
					typename vector_type::value_type sj = b[blockind[j]];
					for(typename matrix_type::const_row_iterator it = A.begin_row(blockind[j]); it != A.end_row(blockind[j]) ; ++it){
						MatMultAdd(sj, 1.0, sj, -1.0, it.value(), x[it.index()]);
					};
					s.subassign(j,sj);
				};
				// solve block
				InverseMatMult(localx,1,mat,s);
				for (size_t j=0;j<blocksize;j++){
					x[blockind[j]] = m_relax*localx[j];
				};
			if (sortedi==0) break;
		}
		}
		return true;
	}

};

} // end namespace ug

#endif // __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__LINE_SMOOTHERS__
