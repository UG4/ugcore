/*
 * Copyright (c) 2017-now:  G-CSC, Goethe University Frankfurt
 * Authors: Arne Nägel
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

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__UZAWA__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__UZAWA__

#include <vector>
#include <string>


#include "lib_algebra/operator/algebra_debug_writer.h"
#include "lib_algebra/adapter/slicing.h"

#include "common/util/string_util.h"  // string tokenizer


namespace ug{

	using binary_grouping_vector = std::vector<bool>;


template<typename TGroupObj, typename TGridFunction>
void ExtractByObject(std::vector<DoFIndex>& vIndex,
	                 const TGridFunction& c,
	                 const std::vector<size_t>& vFullRowCmp,
	                 const std::vector<size_t>& vRemainCmp)
{
	// 	memory for local algebra
	using Element = typename TGridFunction::element_type;
	std::vector<Element*> vElem;


	// loop all grouping objects
	using GroupObjIter = typename TGridFunction::template traits<TGroupObj>::const_iterator;
	for(GroupObjIter iter = c.template begin<TGroupObj>();
			iter != c.template end<TGroupObj>(); ++iter)
	{
		// 	get grouping obj
		TGroupObj* groupObj = *iter;

		// 	get all dof indices on obj associated to cmps

		for(size_t f = 0; f < vFullRowCmp.size(); ++f)
			c.inner_dof_indices(groupObj, vFullRowCmp[f], vIndex, false);
	}
}

template <typename TGridFunction>
class UzawaSlicing : public SlicingData<binary_grouping_vector, 2>  {
public:
	using base_type = SlicingData;

	// build an object with
	UzawaSlicing(const std::vector<std::string>& vSchurCmps) // ø todo remove unused variable from interface
	{}

	void init(const TGridFunction &u, const std::vector<std::string>& vSchurCmps);

protected:

};



template <typename TGridFunction>
void UzawaSlicing <TGridFunction>::
init(const TGridFunction &u, const std::vector<std::string>& vSchurCmps)
{
	UG_DLOG(SchurDebug, 3, "UzawaSlicing::init" << std::endl);

	ConstSmartPtr<DoFDistributionInfo> ddinfo =
			u.approx_space()->dof_distribution_info();

	// fill this vector
	std::vector<DoFIndex> vIndex;

	// ids of selected functions
	std::vector<size_t> vFullRowCmp;

	// complementing functions
	std::vector<size_t> vRemainCmp;

	// tokenize passed functions
	UG_ASSERT(!vSchurCmps.empty(), "UzawaBase::init: No components set.");

	// get ids of selected functions
	for(size_t i = 0; i < vSchurCmps.size(); ++i)
		vFullRowCmp.push_back(ddinfo->fct_id_by_name(vSchurCmps[i].c_str()));


	// Obtain remaining (non-Schur) functions.
	for(size_t f = 0; f < ddinfo->num_fct(); ++f)
		if(std::find(vFullRowCmp.begin(), vFullRowCmp.end(), f) == vFullRowCmp.end())
			vRemainCmp.push_back(f);



	// Extract for each dim-grouping objects.
	for(int d = VERTEX; d <= VOLUME; ++d)
	{

		//	only extract if carrying dofs
		int iCarryDoFs = 0;
		for(size_t f = 0; f < vFullRowCmp.size(); ++f)
			iCarryDoFs += ddinfo->max_fct_dofs(vFullRowCmp[f], d);
		if(iCarryDoFs == 0) continue;

		//	extract
		switch(d){
		case VERTEX: ExtractByObject<Vertex, TGridFunction>(vIndex, u, vFullRowCmp, vRemainCmp); break;
		case EDGE:   ExtractByObject<Edge, TGridFunction>(vIndex, u, vFullRowCmp, vRemainCmp); break;
		case FACE:   ExtractByObject<Face, TGridFunction>(vIndex, u, vFullRowCmp, vRemainCmp); break;
		case VOLUME: ExtractByObject<Volume, TGridFunction>(vIndex, u, vFullRowCmp, vRemainCmp); break;
		default: UG_THROW("wrong dim");
		}

		UG_DLOG(SchurDebug, 4, "Found "<< vIndex.size() << " indices ( out of "<< u.size() << ") for Schur block after dimension "<< d << std::endl) ;
	}


	// todo: replace by InputIterator
	base_type::slice_desc_type_vector mapping(u.size(), false);

	for (size_t i=0; i<vIndex.size(); ++i)
	{
		UG_ASSERT(vIndex[i][1]==0, "Assuming CPUAlgebra...")
		mapping[vIndex[i][0]] = true;
	}

	UG_DLOG(SchurDebug, 3, "UzawaSlicing::init::set_types" << std::endl);
	base_type::set_types(mapping, true);

}




/*! Base class for an Uzawa iteration */

template <typename TDomain, typename TAlgebra>
class UzawaBase : public IPreconditioner<TAlgebra>
{

	static constexpr bool UZAWA_CMP_SCHUR = true;
	static constexpr bool UZAWA_CMP_DEFAULT = false;

	public:
	/// World dimension
		static constexpr int dim = TDomain::dim;
		using TGridFunction = GridFunction<TDomain, TAlgebra>;
		using TElement = typename TGridFunction::element_type;
		using TSide = typename TGridFunction::side_type;


	///	Algebra types
	/// \{
		using algebra_type = TAlgebra;
		using vector_type = typename TAlgebra::vector_type;
		using matrix_type = typename TAlgebra::matrix_type;
		using matrix_operator_type = MatrixOperator<matrix_type, vector_type>;

	/// \}

	protected:
		using base_type = IPreconditioner<TAlgebra>;

		enum BLOCK { AUX_A11, B12, B21, AUX_C22, AUX_M22, AUX_ARRAY_SIZE};


	public:
	///	default constructor
		UzawaBase(const std::vector<std::string>& vSchurCmp)
		:  m_bInit(false), m_vSchurCmp(vSchurCmp), m_slicing(vSchurCmp), m_dSchurUpdateWeight(0.0)
		{
			init_block_operators();
			for (size_t i=0; i<vSchurCmp.size(); ++i)
			{ std::cout << "Comp = " << vSchurCmp[i] << std::endl; }
		}

		UzawaBase(const char *sSchurCmp)
		:  m_bInit(false), m_vSchurCmp(TokenizeTrimString(sSchurCmp)), m_slicing(m_vSchurCmp), m_dSchurUpdateWeight(0.0)
		{
			init_block_operators();
			{ std::cout << "Comp = " << sSchurCmp << std::endl; }
		}


		~UzawaBase() override = default;

		/// Overriding base type
		virtual bool init(SmartPtr<ILinearOperator<vector_type> > J, const vector_type& u)
		{
			UG_DLOG(SchurDebug, 4, "UzawaBase::init(J,u)")

			SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp =
					J.template cast_dynamic<MatrixOperator<matrix_type, vector_type> >();

			const TGridFunction* pVecU = dynamic_cast<const TGridFunction* >(&u);

			UG_ASSERT(pOp.valid(),  "Need a matrix based operator here!");
			UG_ASSERT(pVecU!=nullptr,  "Need a GridFunction here!");
			base_type::init(J,u);   // calls preprocess

			// Extract sub-matrices afterwards.
			init_in_first_step(*pOp, *pVecU);
			m_bInit = true;

			return true;


		}

		SmartPtr<ILinearIterator<typename TAlgebra::vector_type> > clone() override;
protected:
		// Interface for IPreconditioner

		///	name
		const char* name() const override {return "IUzawaBase";}

		///	initializes the preconditioner
		bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp) override;

		///	computes a new correction c = B*d
		bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d) override;


		///	cleans the operator
		bool postprocess() override;

public:

			///	returns if parallel solving is supported
		bool supports_parallel() const override {
				UG_ASSERT(m_spForwardInverse.valid() && m_SpSchurComplementInverse.valid(), "Need valid iter");
				return (m_spForwardInverse->supports_parallel() && m_SpSchurComplementInverse->supports_parallel());

			}

			// forward approximate inverse
			void set_forward_iter(SmartPtr<ILinearIterator<vector_type> > iter)
			{ m_spForwardInverse = iter; }

			// Schur approximate inverse
			void set_schur_iter(SmartPtr<ILinearIterator<vector_type> > iter)
			{ m_SpSchurComplementInverse = iter; }

			// backward approximate inverse
			void set_backward_iter(SmartPtr<ILinearIterator<vector_type> > iter)
			{ m_spBackwardInverse = iter; }

			// assembly for update of Schur operator
			void set_schur_operator_update(SmartPtr<AssembledLinearOperator<TAlgebra> > spSchurUpdateOp, double theta=0.0)
			{ m_spSchurUpdateOp = spSchurUpdateOp; m_dSchurUpdateWeight=theta;}


			/// allocate block matrix operators
			void init_block_operators()
			{
				m_auxMat[AUX_A11] = make_sp(new matrix_operator_type());
				m_auxMat[B12]  	  = make_sp(new matrix_operator_type());
				m_auxMat[B21]  	  = make_sp(new matrix_operator_type());
				m_auxMat[AUX_C22] = make_sp(new matrix_operator_type());

				m_auxMat[AUX_M22] = make_sp(new matrix_operator_type());
			}

			/// extract block matrix operators (called once)
			void extract_sub_matrices(const matrix_type& K, const TGridFunction& c);

			/// Update C22 block by matrix.
			void extract_schur_update(const matrix_type& K, const TGridFunction& c)
			{
				if (m_spSchurUpdateOp.invalid()) return;


				const GridLevel clevel =  c.grid_level();
				my_write_debug(K, "init_KFull_ForSchurUpdate", clevel, clevel);

				//Assembling auxiliary matrix
				SmartPtr<IAssemble<TAlgebra> > m_spAss = m_spSchurUpdateOp->discretization();
				matrix_type tmpM;

				//m_spAss->ass_tuner()->set_force_regular_grid(true);
				//m_spAss->assemble_jacobian(tmpM, c, clevel);
				m_spAss->assemble_mass_matrix(tmpM, c, clevel);
				//m_spAss->ass_tuner()->set_force_regular_grid(false);

				my_write_debug(tmpM, "init_MFull_ForSchurUpdate", clevel, clevel);
				UG_DLOG(SchurDebug, 5, "extract_schur_update on level "<<  clevel << ", " << tmpM);


				/* Matrices C22 and M22 are both additive. Thus, we add both.
				 * (Note, that preconds requiring a consistent diagonal must be called afterwards!)*/
				m_slicing.get_matrix(tmpM, UZAWA_CMP_SCHUR,  UZAWA_CMP_SCHUR,
									*(m_auxMat[AUX_M22].template cast_static<matrix_type>()) );

				if (m_spSliceDebugWriter[UZAWA_CMP_SCHUR].valid()) {
					m_spSliceDebugWriter[UZAWA_CMP_SCHUR]->write_matrix(*m_auxMat[AUX_M22],
															"UZAWA_init_M22_ForSchurUpdate.mat");
				}

				UG_DLOG(SchurDebug, 4, "AUX_C22:"); CheckRowIterators(*m_auxMat[AUX_C22]);
				UG_DLOG(SchurDebug, 4, "AUX_M22:"); CheckRowIterators(*m_auxMat[AUX_M22]);

				// add matrix
				MatAddNonDirichlet<matrix_type>(*m_auxMat[AUX_C22], 1.0, *m_auxMat[AUX_C22], m_dSchurUpdateWeight, *m_auxMat[AUX_M22]);

				if (m_spSliceDebugWriter[UZAWA_CMP_SCHUR].valid()) {
					m_spSliceDebugWriter[UZAWA_CMP_SCHUR]->write_matrix(*m_auxMat[AUX_C22],
																		"UZAWA_init_C22_AfterSchurUpdate.mat");
				}
			}

			void init_block_iterations()
			{
				if (m_spForwardInverse.valid()) { m_spForwardInverse->init(m_auxMat[AUX_A11]); }
				if (m_SpSchurComplementInverse.valid()) { m_SpSchurComplementInverse->init(m_auxMat[AUX_C22]); }
				if (m_spBackwardInverse.valid()) { m_spBackwardInverse->init(m_auxMat[AUX_A11]); }
			}

/*
			void power_iteration(ConstSmartPtr<vector_type> usol, ConstSmartPtr<vector_type> dsol)
			{
				SmartPtr<vector_type> uschur(m_slicing.template slice_clone<vector_type>(*usol, UZAWA_CMP_SCHUR) );
				SmartPtr<vector_type> dschur(m_slicing.template slice_clone<vector_type>(*dsol, UZAWA_CMP_SCHUR) )
				SmartPtr<vector_type> ueigen(m_slicing.template slice_clone<vector_type>(*usol, UZAWA_CMP_DEFAULT) );
				SmartPtr<vector_type> deigen(m_slicing.template slice_clone<vector_type>(*usol, UZAWA_CMP_DEFAULT) );

				// Initialize with random guess
				ueigen->set_random(-0.5, 0.5);

				// dschur = B12*ueigen


				//
				m_spForwardInverse->apply(*uschur, *dschur);

			    // deigen = B21 * uschur

			}
*/

			void postprocess_block_iterations()
			{

			}

protected:
			void init_in_first_step(const matrix_type &pMat, const TGridFunction &pC)
			{
				UG_DLOG(SchurDebug, 4, "step-init: Size=" << m_vSchurCmp.size() << std::endl);
				m_slicing.init(pC, m_vSchurCmp);

				if (debug_writer().valid())
				{

					if (m_spGridFunctionDebugWriter.valid())
					{
						UG_LOG("Valid grid function writer for "<< m_spGridFunctionDebugWriter->grid_level() << " on level " << pC.grid_level());

						GridLevel gLevel = m_spGridFunctionDebugWriter->grid_level();
						m_spGridFunctionDebugWriter->set_grid_level(pC.grid_level());
						m_spGridFunctionDebugWriter->update_positions();
						m_spGridFunctionDebugWriter->set_grid_level(gLevel);
					}


					switch(debug_writer()->get_dim())
					{
						case 1:	reset_slice_debug_writers<1>(); break;
						case 2:	reset_slice_debug_writers<2>(); break;
						case 3:	reset_slice_debug_writers<3>(); break;
						default: UG_LOG("Invalid dimension for debug_write ???");
					}
				}

				extract_sub_matrices(pMat, pC);
				extract_schur_update(pMat, pC);

				// Initialize preconditioners.
				init_block_iterations();
			}
protected:


		void my_write_debug(const matrix_type& mat, std::string name, const TGridFunction& rTo, const TGridFunction& rFrom)
		{
			my_write_debug(mat, name, rTo.grid_level(), rFrom.grid_level());
		}

		void my_write_debug(const matrix_type& mat, std::string name, const GridLevel& glTo, const GridLevel& glFrom)
		{
			PROFILE_FUNC_GROUP("debug");

			if(m_spGridFunctionDebugWriter.invalid()) return;

		//	build name
			std::stringstream ss;
			ss << "UZAWA_" << name << GridLevelAppendix(glTo);
			if(glFrom != glTo) ss << GridLevelAppendix(glFrom);
			ss << ".mat";

		//	write
			GridLevel currGL = m_spGridFunctionDebugWriter->grid_level();
			m_spGridFunctionDebugWriter->set_grid_levels(glFrom, glTo);
			m_spGridFunctionDebugWriter->write_matrix(mat, ss.str().c_str());
			m_spGridFunctionDebugWriter->set_grid_level(currGL);
		}

		virtual void my_write_debug(const TGridFunction& rGF, std::string name)
		{
			int m_dbgIterCnt = 0;
			PROFILE_FUNC_GROUP("debug");

			if(m_spGridFunctionDebugWriter.invalid()) return;

		//	build name
			GridLevel gl = rGF.grid_level();
			std::stringstream ss;
			ss << "UZAWA_" << name << GridLevelAppendix(gl);
			ss << "_i" << std::setfill('0') << std::setw(3) << m_dbgIterCnt << ".vec";

		//	write
			GridLevel currGL = m_spGridFunctionDebugWriter->grid_level();
			m_spGridFunctionDebugWriter->set_grid_level(gl);
			m_spGridFunctionDebugWriter->write_vector(rGF, ss.str().c_str());
			m_spGridFunctionDebugWriter->set_grid_level(currGL);
		};

protected:
	using base_type::debug_writer;

	void set_debug(SmartPtr<IDebugWriter<algebra_type> > spDebugWriter);

	void create_slice_debug_writers();
	template <int d> void reset_slice_debug_writers();


	/// flag indicating, whether operator must be initialized
	bool m_bInit;

	/// vector of strings identifying components used for Schur complement
	std::vector<std::string> m_vSchurCmp;

	/// object for slicing routines
	UzawaSlicing<TGridFunction> m_slicing;

	///	iteration for forward system
	SmartPtr<ILinearIterator<vector_type> >  m_spForwardInverse;

	///	iteration for forward system
	SmartPtr<ILinearIterator<vector_type> > m_SpSchurComplementInverse;

	///	iteration for forward system
	SmartPtr<ILinearIterator<vector_type> > m_spBackwardInverse;

	///	assembly for (additive) Schur complement update
	SmartPtr<AssembledLinearOperator<TAlgebra> > m_spSchurUpdateOp;

	///	scaling factor for (additive) Schur complement update
	double m_dSchurUpdateWeight;
	SmartPtr<matrix_operator_type> m_auxMat[AUX_ARRAY_SIZE];			/// auxiliary matrices (not cloned!)


	SmartPtr<GridFunctionDebugWriter<TDomain, TAlgebra> >m_spGridFunctionDebugWriter;
	SmartPtr<IDebugWriter<algebra_type> > m_spSliceDebugWriter[2];

#ifdef UG_PARALLEL
		matrix_type m_A;
#endif
};


template <typename TDomain, typename TAlgebra>
bool UzawaBase<TDomain, TAlgebra>::
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

	// more preprocessing is based on grid functions
	// thus, it must be performed later...
	//m_bInit = false;

	return true;
}


//! Single step of an (inexact) Uzawa iteration
template <typename TDomain, typename TAlgebra>
bool UzawaBase<TDomain, TAlgebra>::
step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d)
{
	// Check that grid function was passed.
	GridFunction<TDomain, TAlgebra>* pC = dynamic_cast<GridFunction<TDomain, TAlgebra>*>(&c);
	if(pC == nullptr) UG_THROW("UzawaBase: expects correction to be a GridFunction.");

	const vector_type* pD = &d;
	const matrix_type* pMat = pOp.get();

#ifdef UG_PARALLEL
	SmartPtr<vector_type> spDtmp;
	if(pcl::NumProcs() > 1)
	{
		// Make defect unique
		spDtmp = d.clone();
		spDtmp->change_storage_type(PST_UNIQUE);
		pD = spDtmp.get();
		pMat = &m_A;
	}
#endif

	//	Check, if initialized.
	if(!m_bInit)
	{
		init_in_first_step(*pMat, *pC);
		m_bInit = true;

	}

	// Obtain defects.
	SmartPtr<vector_type> ff(m_slicing.template slice_clone<vector_type>(*pD, UZAWA_CMP_DEFAULT) );
	SmartPtr<vector_type> gg(m_slicing.template slice_clone<vector_type>(*pD, UZAWA_CMP_SCHUR) );

	// Clear corrections.
	SmartPtr<vector_type> cRegular(m_slicing.template slice_clone_without_values<vector_type>(*pC, UZAWA_CMP_DEFAULT));
	SmartPtr<vector_type> cSchur(m_slicing.template slice_clone_without_values<vector_type>(*pC, UZAWA_CMP_SCHUR));

	pC->set(0.0);
	cRegular->set(0.0);
	cSchur->set(0.0);

	// Step 1: Solve problem \tilde A11 \tilde u = f
	my_write_debug(*pC, "UZAWA_Correction0");
	my_write_debug(*(GridFunction<TDomain, TAlgebra>*) pD, "UZAWA_Defect0");
	if (m_spForwardInverse.valid())
	{
		// Solve problem, also compute  g:=(f - A11 \tilde u)
		// WARNING: Should use exact A11.
	   // m_spForwardInverse->apply_update_defect(*cRegular, *ff);
		m_spForwardInverse->apply(*cRegular, *ff);

		m_slicing.template set_vector_slice<vector_type>(*cRegular, *pC, UZAWA_CMP_DEFAULT);
		my_write_debug(*pC, "UZAWA_Correction1");

		// Update g:= (g - B21 \tilde u^k)
		m_auxMat[B21]->apply_sub(*gg, *cRegular);
	}

	// Step 2: Solve problem: (\tilde S22) pDelta = (g - B21 \tilde u)
	if (m_SpSchurComplementInverse.valid()) {

		//m_auxMat[AUX_C22]->apply_sub(*gg, *cRegular);
		m_SpSchurComplementInverse->apply(*cSchur, *gg);
		m_slicing.template set_vector_slice<vector_type>(*cSchur, *pC, UZAWA_CMP_SCHUR);
		my_write_debug(*pC, "UZAWA_Correction2");

	    // update defect
		m_auxMat[B12]->apply_sub(*ff, *cSchur);
	}

	// Step 3: Solve problem \tilde A11 uDelta^{k+1} = (f - \tilde A11 u^k) - B12 p^k
	if (m_spBackwardInverse.valid()) {

		// ff has been updated
		m_spBackwardInverse->apply(*cRegular, *ff);
		// m_slicing.template add_vector_slice<vector_type>(*cRegular, *pC, UZAWA_CMP_DEFAULT);
		m_slicing.template set_vector_slice<vector_type>(*cRegular, *pC, UZAWA_CMP_DEFAULT);
		my_write_debug(*pC, "UZAWA_Correction3");
	}


#ifdef UG_PARALLEL
	 pC->set_storage_type(PST_UNIQUE);
#endif

	 // UG_ASSERT(0, "STOP HERE!")
	return true;
}

template <typename TDomain, typename TAlgebra>
void UzawaBase<TDomain, TAlgebra>::
extract_sub_matrices(const matrix_type& K, const TGridFunction& c)
{
	// Copy matrices.
	m_slicing.get_matrix(K, UZAWA_CMP_DEFAULT, UZAWA_CMP_DEFAULT, *(m_auxMat[AUX_A11].template cast_static<matrix_type>()) );
	m_slicing.get_matrix(K, UZAWA_CMP_DEFAULT, UZAWA_CMP_SCHUR,   *(m_auxMat[B12].template cast_static<matrix_type>()) );
	m_slicing.get_matrix(K, UZAWA_CMP_SCHUR,   UZAWA_CMP_DEFAULT, *(m_auxMat[B21].template cast_static<matrix_type>()) );
	m_slicing.get_matrix(K, UZAWA_CMP_SCHUR,   UZAWA_CMP_SCHUR,   *(m_auxMat[AUX_C22].template cast_static<matrix_type>()) );

	UG_DLOG(SchurDebug, 4, "A11 =" << *m_auxMat[AUX_A11] << ", ");
	UG_DLOG(SchurDebug, 4, "B12 =" << *m_auxMat[B12] << ", ");
	UG_DLOG(SchurDebug, 4, "B21 =" << *m_auxMat[B21] << ", ");
	UG_DLOG(SchurDebug, 4, "C22 =" << *m_auxMat[AUX_C22] << std::endl);

#ifdef UG_PARALLEL
	// Copy storage mask.
	uint mask_K = K.get_storage_mask();
	m_auxMat[AUX_A11]->set_storage_type(mask_K);
	m_auxMat[B12]->set_storage_type(mask_K);
	m_auxMat[B21]->set_storage_type(mask_K);
	m_auxMat[AUX_C22]->set_storage_type(mask_K);

	// Extract and set layouts.
	SmartPtr<AlgebraLayouts> defaultLayouts = m_slicing.create_slice_layouts(K.layouts(), UZAWA_CMP_DEFAULT);
	SmartPtr<AlgebraLayouts> schurLayouts   = m_slicing.create_slice_layouts(K.layouts(), UZAWA_CMP_SCHUR);

	m_auxMat[AUX_A11]->set_layouts(defaultLayouts);
	m_auxMat[AUX_C22]->set_layouts(schurLayouts);

#endif


	// Write matrices for debug, if applicable.
	if (m_spSliceDebugWriter[UZAWA_CMP_DEFAULT].valid()) {
		m_spSliceDebugWriter[UZAWA_CMP_DEFAULT]->write_matrix(*m_auxMat[AUX_A11], "UZAWA_init_A11_AfterExtract.mat");
	}
	//write_debug(*m_auxMat[B12], "Uzawa_init_B12_AfterExtract");
	//write_debug(*m_auxMat[B21], "Uzawa_init_B21_AfterExtract");
	if (m_spSliceDebugWriter[UZAWA_CMP_SCHUR].valid()) {
		m_spSliceDebugWriter[UZAWA_CMP_SCHUR]->write_matrix(*m_auxMat[AUX_C22], "UZAWA_init_C22_AfterExtract.mat");
	}
}


template <typename TDomain, typename TAlgebra>
bool UzawaBase<TDomain, TAlgebra>::
postprocess()
{
	postprocess_block_iterations();
	m_bInit = false;

	return true;
}



template <typename TDomain, typename TAlgebra>
SmartPtr<ILinearIterator<typename TAlgebra::vector_type> >
UzawaBase<TDomain, TAlgebra>::clone()
{
	SmartPtr<UzawaBase<TDomain, TAlgebra> >
					newInst(new UzawaBase<TDomain, TAlgebra>(m_vSchurCmp));

	newInst->set_debug(debug_writer());

	// clone approximate inverses
	newInst->m_spForwardInverse = (m_spForwardInverse.valid()) ? m_spForwardInverse->clone().template cast_dynamic<base_type>() : nullptr;
	newInst->m_SpSchurComplementInverse = (m_SpSchurComplementInverse.valid()) ? m_SpSchurComplementInverse->clone().template cast_dynamic<base_type>() : nullptr;
	newInst->m_spBackwardInverse = (m_spBackwardInverse.valid()) ? m_spBackwardInverse->clone().template cast_dynamic<base_type>() : nullptr;

	// todo: need to clone here?
	newInst->m_spSchurUpdateOp = m_spSchurUpdateOp;
	newInst->m_dSchurUpdateWeight = m_dSchurUpdateWeight;

	return newInst;
}



template <typename TDomain, typename TAlgebra>
void UzawaBase<TDomain, TAlgebra>::
create_slice_debug_writers()
{
	std::string basePath = debug_writer()->get_base_dir();

	m_spSliceDebugWriter[UZAWA_CMP_DEFAULT] = make_sp(new AlgebraDebugWriter<algebra_type>);
	m_spSliceDebugWriter[UZAWA_CMP_DEFAULT] ->set_base_dir(basePath.c_str());

	m_spSliceDebugWriter[UZAWA_CMP_SCHUR] = make_sp(new AlgebraDebugWriter<algebra_type>);
	m_spSliceDebugWriter[UZAWA_CMP_SCHUR] ->set_base_dir(basePath.c_str());
}

template <typename TDomain, typename TAlgebra>
template<int dim>
void UzawaBase<TDomain, TAlgebra>::
reset_slice_debug_writers()
{
	UG_LOG("reset_slice_debug_writers"<< std::endl);

	if (m_spSliceDebugWriter[UZAWA_CMP_DEFAULT].invalid() ||
		m_spSliceDebugWriter[UZAWA_CMP_SCHUR].invalid()) return;

	const std::vector<MathVector<dim> > &positions = (m_spGridFunctionDebugWriter.valid()) ? m_spGridFunctionDebugWriter->template get_positions<dim>() : debug_writer()->template get_positions<dim>();
	std::vector<MathVector<dim> > cmpPositions[2];

	m_slicing.get_vector_slice(positions, UZAWA_CMP_DEFAULT, cmpPositions[UZAWA_CMP_DEFAULT]);
	m_spSliceDebugWriter[UZAWA_CMP_DEFAULT]->template set_positions<dim>(cmpPositions[UZAWA_CMP_DEFAULT]);

	m_slicing.get_vector_slice(positions, UZAWA_CMP_SCHUR, cmpPositions[UZAWA_CMP_SCHUR]);
	m_spSliceDebugWriter[UZAWA_CMP_SCHUR]->template set_positions<dim>(cmpPositions[UZAWA_CMP_SCHUR]);

}

template <typename TDomain, typename TAlgebra>
void UzawaBase<TDomain, TAlgebra>::set_debug(SmartPtr<IDebugWriter<algebra_type> > spDebugWriter)
{
	base_type::set_debug(spDebugWriter);
	m_spGridFunctionDebugWriter = debug_writer().template cast_dynamic<GridFunctionDebugWriter<TDomain, TAlgebra> >();

	if (m_spGridFunctionDebugWriter.invalid()) return;

	create_slice_debug_writers();

	debug_writer()->update_positions();
	switch(debug_writer()->get_dim())
	{
		case 1:	reset_slice_debug_writers<1>(); break;
		case 2:	reset_slice_debug_writers<2>(); break;
		case 3:	reset_slice_debug_writers<3>(); break;
		default: UG_LOG("debug dim ???");
	}

}


} // namespace ug
#endif
