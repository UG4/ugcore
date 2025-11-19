/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#ifndef __H__UG__LIB_DISC__FUNCTION_SPACE__GRID_FUNCTION_USER_DATA__
#define __H__UG__LIB_DISC__FUNCTION_SPACE__GRID_FUNCTION_USER_DATA__

#include "common/common.h"

#include "lib_grid/tools/subset_group.h"

#include "lib_disc/common/function_group.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/quadrature/quadrature.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/spatial_disc/user_data/std_user_data.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"


namespace ug{

/**
 * Class for StdDependentUserData interface to scalar values from GridFunction objects.
 *
 * This class provides a StdDependentUserData interface to scalar values from GridFunction
 * objects so that they can be used in import parameters. Note that the GridFunction object
 * must have the same FunctionPattern as the owner of the import parameter.
 *
 * REMARK:
 * This class is used to reuse the unknowns of variables of discretizations. It provides
 * "derivatives w.r.t. the unknowns", too. Note that in this case, the GridFunction object
 * must be the same as the solution. An attempt to use a grid function based on a different
 * ApproximationSpace would lead to errors. Consider ExplicitGridFunctionValue
 * for this purpose instead.
 */
template <typename TGridFunction>
class GridFunctionNumberData
: public StdDependentUserData<GridFunctionNumberData<TGridFunction>,
  number, TGridFunction::dim>
{
	public:
		//	world dimension of grid function
		static constexpr int dim = TGridFunction::dim;

		private:
		// grid function
		SmartPtr<TGridFunction> m_spGridFct;

		//	component of function
		size_t m_fct;

		//	local finite element id
		LFEID m_lfeID;

	public:
		/// constructor
		GridFunctionNumberData(SmartPtr<TGridFunction> spGridFct, const char* cmp)
		{
			assign(spGridFct, cmp);
		}

		GridFunctionNumberData(){}

		void assign(SmartPtr<TGridFunction> spGridFct, const char* cmp)
		{
			m_spGridFct = spGridFct;
			this->set_functions(cmp);

			//	get function id of name
			m_fct = spGridFct->fct_id_by_name(cmp);

			//	check that function exists
			if(m_fct >= spGridFct->num_fct())
				UG_THROW("GridFunctionNumberData: Function space does not contain"
						" a function with name " << cmp << ".");

			//	local finite element id
			m_lfeID = spGridFct->local_finite_element_id(m_fct);
		}

		virtual bool continuous() const
		{
			return LocalFiniteElementProvider::continuous(m_lfeID);
		}

		template <int refDim>
		void eval_and_deriv(number vValue[],
		                    const MathVector<dim> vGlobIP[],
		                    number time, int si,
		                    GridObject* elem,
		                    const MathVector<dim> vCornerCoords[],
		                    const MathVector<refDim> vLocIP[],
		                    const size_t nip,
		                    LocalVector* u,
		                    bool bDeriv,
		                    int s,
		                    std::vector<std::vector<number> > vvvDeriv[],
		                    const MathMatrix<refDim, dim>* vJT = nullptr) const
		{
			//	reference object id
			const ReferenceObjectID roid = elem->reference_object_id();

			//	get trial space
			try{
				const LocalShapeFunctionSet<refDim>& rTrialSpace =
						LocalFiniteElementProvider::get<refDim>(roid, m_lfeID);

				//	memory for shapes
				std::vector<number> vShape;
				//	memory for indices
				std::vector<DoFIndex> ind;

				//	loop ips
				for(size_t ip = 0; ip < nip; ++ip)
				{
					//	evaluate at shapes at ip
					rTrialSpace.shapes(vShape, vLocIP[ip]);

					//	get multiindices of element
					m_spGridFct->dof_indices(elem, m_fct, ind);

					// 	compute solution at integration point
					vValue[ip] = 0.0;
					for(size_t sh = 0; sh < vShape.size(); ++sh)
					{
						const number valSH = DoFRef(*m_spGridFct, ind[sh]);
						vValue[ip] += valSH * vShape[sh];
					}
				}

				if(bDeriv){
					for(size_t ip = 0; ip < nip; ++ip){
						//	evaluate at shapes at ip
						rTrialSpace.shapes(vShape, vLocIP[ip]);

						for(size_t sh = 0; sh < vShape.size(); ++sh)
							vvvDeriv[ip][0][sh] = vShape[sh];
					}
				}
			}
			UG_CATCH_THROW("GridFunctionNumberData: Shape Function Set missing for"
					" Reference Object: "<<roid<<", Trial Space: "
					<<m_lfeID<<", refDim="<<refDim);

		}
	
		virtual void operator() (number& value,
									const MathVector<dim>& globIP,
									number time, int si,
									Vertex* vrt) const
		{
			std::vector<DoFIndex> ind;
			if (m_spGridFct->dof_indices(vrt, m_fct, ind) != 1)
				UG_THROW ("GridFunctionNumberData: Values at vertices of the grid function are not uniquely defined.");
			value = DoFRef(*m_spGridFct, ind[0]);
		}
};

/**
 * Class for StdDependentUserData interface to vectors from GridFunction objects.
 *
 * This class provides a StdDependentUserData interface to vectors from GridFunction
 * objects so that they can be used in import parameters. Note that the GridFunction object
 * must have the same FunctionPattern as the owner of the import parameter.
 *
 * REMARK:
 * This class is used to reuse the unknowns of variables of discretizations. It provides
 * "derivatives w.r.t. the unknowns", too. Note that in this case the GridFunction object
 * must be the same as the solution. An attempt to use a grid function based on a different
 * ApproximationSpace would lead to errors. Consider ExplicitGridFunctionVector
 * for this purpose instead.
 */
template <typename TGridFunction>
class GridFunctionVectorData
: public StdDependentUserData<GridFunctionVectorData<TGridFunction>,
  MathVector<TGridFunction::dim>, TGridFunction::dim>
{
	public:
	//	world dimension of grid function
	static constexpr int dim = TGridFunction::dim;

	private:
	// grid function
	SmartPtr<TGridFunction> m_spGridFct;

	//	component of function
	size_t m_vfct[dim];

	//	local finite element id
	LFEID m_vlfeID[dim];

	public:
	/// constructor
	GridFunctionVectorData(SmartPtr<TGridFunction> spGridFct, const char* cmp)
	: m_spGridFct(spGridFct)
	{
		this->set_functions(cmp);

		//	create function group of this elem disc
		try{
			//	get strings
			std::string fctString = std::string(cmp);

			//	tokenize strings and select functions
			std::vector<std::string> tokens;
			TokenizeString(fctString, tokens, ',');

			for(size_t i = 0; i < tokens.size(); ++i)
				RemoveWhitespaceFromString(tokens[i]);

			if((int)tokens.size() != dim)
				UG_THROW("GridFunctionVectorData: Needed "<<dim<<" components "
				         "in symbolic function names, but given: "<<cmp);

			//	get function id of name
			for(int i = 0; i < dim; ++i){
				m_vfct[i] = spGridFct->fct_id_by_name(tokens[i].c_str());
				m_vlfeID[i] = spGridFct->local_finite_element_id(m_vfct[i]);
			}

		}UG_CATCH_THROW("GridFunctionVectorData: Cannot find  some symbolic function "
				"name in '"<<cmp<<"'.");
	};

	GridFunctionVectorData(SmartPtr<TGridFunction> spGridFct, std::vector<std::string> vCmp)
	: m_spGridFct(spGridFct)
	{
		this->set_functions(vCmp);

		//	create function group of this elem disc
		try
		{
			if ((int)vCmp.size() != dim)
				UG_THROW("GridFunctionVectorData: Needed "<<dim<<" components "
						 "in symbolic function names, but given: "<< (int) vCmp.size());

			//	get function id of name
			for (size_t i = 0; i < (size_t) dim; ++i)
			{
				m_vfct[i] = spGridFct->fct_id_by_name(vCmp[i].c_str());
				m_vlfeID[i] = spGridFct->local_finite_element_id(m_vfct[i]);
			}

		}UG_CATCH_THROW("GridFunctionVectorData: Cannot find some symbolic function "
				"name in component vector.");
	}

	virtual bool continuous() const
	{
		for(int i = 0; i < dim; ++i)
			if(!LocalFiniteElementProvider::continuous(m_vlfeID[i]))
				return false;
		return true;
	}

	template <int refDim>
	void eval_and_deriv(MathVector<dim> vValue[],
	                    const MathVector<dim> vGlobIP[],
	                    number time, int si,
	                    GridObject* elem,
	                    const MathVector<dim> vCornerCoords[],
	                    const MathVector<refDim> vLocIP[],
	                    const size_t nip,
	                    LocalVector* u,
	                    bool bDeriv,
	                    int s,
	                    std::vector<std::vector<MathVector<dim> > > vvvDeriv[],
	                    const MathMatrix<refDim, dim>* vJT = nullptr) const
	{
		//	reference object id
		const ReferenceObjectID roid = elem->reference_object_id();

		//	memory for shapes
		std::vector<number> vShape;
		//	memory for indices
		std::vector<DoFIndex> ind;

		//	loop ips
		try{
			for(size_t ip = 0; ip < nip; ++ip)
			{
				for(int d = 0; d < dim; ++d)
				{
					const LocalShapeFunctionSet<refDim>& rTrialSpace =
							LocalFiniteElementProvider::get<refDim>(roid, m_vlfeID[d]);

					//	evaluate at shapes at ip
					rTrialSpace.shapes(vShape, vLocIP[ip]);

					//	get multiindices of element
					m_spGridFct->dof_indices(elem, m_vfct[d], ind);

					// 	compute solution at integration point
					vValue[ip][d] = 0.0;
					for(size_t sh = 0; sh < vShape.size(); ++sh)
					{
						const number valSH = DoFRef( *m_spGridFct, ind[sh]);
						vValue[ip][d] += valSH * vShape[sh];
					}
				}
			}


			if(bDeriv){
				for(size_t ip = 0; ip < nip; ++ip){
					for(int d = 0; d < dim; ++d)
					{
						const LocalShapeFunctionSet<refDim>& rTrialSpace =
								LocalFiniteElementProvider::get<refDim>(roid, m_vlfeID[d]);
						//	evaluate at shapes at ip
						rTrialSpace.shapes(vShape, vLocIP[ip]);

						for(size_t sh = 0; sh < vShape.size(); ++sh)
							vvvDeriv[ip][d][sh][d] = vShape[sh];
					}
				}
			}


		}UG_CATCH_THROW("GridFunctionVectorData: Reference Object: "
				<<roid<<", refDim="<<refDim);
	}
	
	virtual void operator() (MathVector<dim>& value,
										const MathVector<dim>& globIP,
										number time, int si,
										Vertex* vrt) const
	{
		std::vector<DoFIndex> ind;
		for(int d = 0; d < dim; ++d)
		{
			if (m_spGridFct->dof_indices(vrt, m_vfct[d], ind) != 1)
				UG_THROW ("GridFunctionVectorData: Values at vertices of the grid function are not uniquely defined.");
			value[d] = DoFRef(*m_spGridFct, ind[0]);
		}
	}
};


/**
 * Class for StdDependentUserData interface to the gradient of scalar values from GridFunction objects.
 *
 * This class provides a StdDependentUserData interface to the gradient of scalar values
 * from GridFunction objects so that they can be used in import parameters. Note that the
 * GridFunction object must have the same FunctionPattern as the owner of the import parameter.
 *
 * REMARK:
 * This class is used to reuse the unknowns of variables of discretizations. It provides
 * "derivatives w.r.t. the unknowns", too. Note that in this case the GridFunction object
 * must be the same as the solution. An attempt to use a grid function based on a different
 * ApproximationSpace would lead to errors. Consider ExplicitGridFunctionGradient
 * for this purpose instead.
 */
template <typename TGridFunction>
class GridFunctionGradientData
: public StdDependentUserData<GridFunctionGradientData<TGridFunction> ,
  MathVector<TGridFunction::dim>, TGridFunction::dim>
{
	public:
	//	world dimension of grid function
	static constexpr int dim = TGridFunction::dim;

	private:
	// grid function
	SmartPtr<TGridFunction> m_spGridFct;

	//	component of function
	size_t m_fct;

	//	local finite element id
	LFEID m_lfeID;

	public:
	/// constructor
	GridFunctionGradientData(SmartPtr<TGridFunction> spGridFct, const char* cmp)
	: m_spGridFct(spGridFct)
	{
		this->set_functions(cmp);

		//	get function id of name
		m_fct = spGridFct->fct_id_by_name(cmp);

		//	check that function exists
		if(m_fct >= spGridFct->num_fct())
			UG_THROW("GridFunctionGradientData: Function space does not contain"
					" a function with name " << cmp << ".");

		//	local finite element id
		m_lfeID = spGridFct->local_finite_element_id(m_fct);
	};

	virtual bool continuous() const
	{
		return false;
	}

	template <int refDim>
	void eval_and_deriv(MathVector<dim> vValue[],
	                    const MathVector<dim> vGlobIP[],
	                    number time, int si,
	                    GridObject* elem,
	                    const MathVector<dim> vCornerCoords[],
	                    const MathVector<refDim> vLocIP[],
	                    const size_t nip,
	                    LocalVector* u,
	                    bool bDeriv,
	                    int s,
	                    std::vector<std::vector<MathVector<dim> > > vvvDeriv[],
	                    const MathMatrix<refDim, dim>* vJT = nullptr) const
	{
		//	reference object id
		const ReferenceObjectID roid = elem->reference_object_id();

		//	get reference element mapping by reference object id
		std::vector<MathMatrix<refDim, dim> > vJTTmp(nip);
		if(vJT == nullptr){
			try{
				DimReferenceMapping<refDim, dim>& mapping
				= ReferenceMappingProvider::get<refDim, dim>(roid, vCornerCoords);

				//	compute transformation matrices
				mapping.jacobian_transposed(&(vJTTmp[0]), vLocIP, nip);

				//	store tmp Gradient
				vJT = &(vJTTmp[0]);
			}UG_CATCH_THROW("GridFunctionGradientData: failed.");
		}

		//	get trial space
		try{
			const LocalShapeFunctionSet<refDim>& rTrialSpace =
					LocalFiniteElementProvider::get<refDim>(roid, m_lfeID);

			//	storage for shape function at ip
			std::vector<MathVector<refDim> > vLocGrad;
			MathVector<refDim> locGrad;

			//	Reference Mapping
			MathMatrix<dim, refDim> JTInv;

			//	loop ips
			for(size_t ip = 0; ip < nip; ++ip)
			{
				//	evaluate at shapes at ip
				rTrialSpace.grads(vLocGrad, vLocIP[ip]);

				//	get multiindices of element
				std::vector<DoFIndex > ind;
				m_spGridFct->dof_indices(elem, m_fct, ind);

				//	compute grad at ip
				VecSet(locGrad, 0.0);
				for(size_t sh = 0; sh < vLocGrad.size(); ++sh)
				{
					const number valSH = DoFRef( *m_spGridFct, ind[sh]);
					VecScaleAppend(locGrad, valSH, vLocGrad[sh]);
				}

				RightInverse (JTInv, vJT[ip]);
				MatVecMult(vValue[ip], JTInv, locGrad);
			}
		}
		UG_CATCH_THROW("GridFunctionNumberData: Shape Function Set missing for"
				" Reference Object: "<<roid<<", Trial Space: "
				<<m_lfeID<<", refDim="<<refDim);

		if(bDeriv)
			UG_THROW("Not implemented.");
	}
};


/**
 * \brief Retrieve component of gradient of GridFunction
 * \details Helper construct to retrieve specific vector element of the gradient 
 *   of a GridFunction inside loops (e.g. integrals) over that GridFunction.
 * \code{.lua}
   -- integrate the second component of the gradient of the GridFunction u
   Integrate( GridFunctionGradientComponentData(u, "c", 2), u )
 * \endcode
 */
template <typename TGridFunction>
class GridFunctionGradientComponentData
: public StdDependentUserData<GridFunctionGradientComponentData<TGridFunction>,
  number, TGridFunction::dim>
{
	public:
	///	World dimension of GridFunction m_spGridFct
	static constexpr int dim = TGridFunction::dim;

	private:
	///	GridFunction to loop over
	SmartPtr<TGridFunction> m_spGridFct;

	///	Component ID of function to be used
	size_t m_fct;

	///	Component index of gradient to return (0-based)
	size_t m_component;

	///	Local Finite Element ID
	LFEID m_lfeID;

	public:
	/**
	 * \brief Constructor
	 * \param[in] spGridFct GridFunction to loop over
	 * \param[in] cmp Name of the GridFunction's function to calculate gradient of
	 * \param[in] component Index of gradient vector component to return (1-based)
	 */
	GridFunctionGradientComponentData( SmartPtr<TGridFunction> spGridFct,
	                                   const char* cmp,
	                                   size_t component /* 1-based */ )
	: m_spGridFct( spGridFct )
	{
		this->set_functions(cmp);

		//	check validity of component index
		if ( component > static_cast<size_t>(dim) && component > 0 ) {
			UG_THROW( "GridFunctionGradientComponentData: Requested component index "
					<< component << " out of bounds [1," << dim << "]." );
		}
		m_component = component - 1;

		//	get function id of name
		m_fct = spGridFct->fct_id_by_name( cmp );

		//	check that function exists
		if( m_fct >= spGridFct->num_fct() ) {
			UG_THROW( "GridFunctionGradientComponentData: Function space does not contain"
					" a function with name " << cmp << "." );
		}

		//	local finite element id
		m_lfeID = spGridFct->local_finite_element_id( m_fct );
	};

	virtual bool continuous() const
	{
		return false;
	}

	/**
	 * \param[out] vValue Array of the <tt>nip</tt> gradient components
	 */
	template <int refDim>
	void eval_and_deriv(number vValue[],
	                    const MathVector<dim> vGlobIP[],
	                    number time, int si,
	                    GridObject* elem,
	                    const MathVector<dim> vCornerCoords[],
	                    const MathVector<refDim> vLocIP[],
	                    const size_t nip,
	                    LocalVector* u,
	                    bool bDeriv,
	                    int s,
	                    std::vector<std::vector<number> > vvvDeriv[],
	                    const MathMatrix<refDim, dim>* vJT = nullptr) const
	{
		//	reference object id
		const ReferenceObjectID roid = elem->reference_object_id();

		//	get reference element mapping by reference object id
		std::vector<MathMatrix<refDim, dim> > vJTTmp( nip );
		if( vJT == nullptr ) {
			try{
				DimReferenceMapping<refDim, dim>& mapping
				= ReferenceMappingProvider::get< refDim, dim >( roid, vCornerCoords );

				//	compute transformation matrices
				mapping.jacobian_transposed( &(vJTTmp[0]), vLocIP, nip );

				//	store tmp Gradient
				vJT = &(vJTTmp[0]);
			} UG_CATCH_THROW( "GridFunctionGradientComponentData: failed.");
		}

		//	get trial space
		try {
			const LocalShapeFunctionSet<refDim>& rTrialSpace =
					LocalFiniteElementProvider::get<refDim>( roid, m_lfeID );

			//	storage for shape function at ip
			std::vector<MathVector<refDim> > vLocGrad;
			MathVector<refDim> locGrad;
			std::vector<MathVector<dim> > vValueVec;
			vValueVec.resize( nip );

			//	Reference Mapping
			MathMatrix<dim, refDim> JTInv;

			//	loop ips
			for( size_t ip = 0; ip < nip; ++ip ) {
				//	evaluate at shapes at ip
				rTrialSpace.grads( vLocGrad, vLocIP[ip] );

				//	get multiindices of element
				std::vector<DoFIndex> ind;
				m_spGridFct->dof_indices( elem, m_fct, ind );

				//	compute grad at ip
				VecSet( locGrad, 0.0 );
				for( size_t sh = 0; sh < vLocGrad.size(); ++sh ) {
					const number valSH = DoFRef( *m_spGridFct, ind[sh] );
					VecScaleAppend( locGrad, valSH, vLocGrad[sh] );
				}

				RightInverse( JTInv, vJT[ip] );
				MatVecMult( vValueVec[ip], JTInv, locGrad );

				vValue[ip] = vValueVec[ip][m_component];

				if(bDeriv){
					for( size_t ip = 0; ip < nip; ++ip ) {
						UG_ASSERT(vvvDeriv[ip].size() == 1,
						          "Single component expected, but "<<vvvDeriv[ip].size())
						UG_ASSERT(vvvDeriv[ip][0].size() == vLocGrad.size(),
						          "Wrong number sh: "<<vvvDeriv[ip][0].size()<<", but expected: "<<vLocGrad.size())

						for( size_t sh = 0; sh < vLocGrad.size(); ++sh ) {
							MatVecMult( vValueVec[ip], JTInv, vLocGrad[sh] );
							vvvDeriv[ip][0][sh] = vValueVec[ip][m_component];
						}
					}
				}
			}
		}
		UG_CATCH_THROW("GridFunctionNumberData: Shape Function Set missing for"
				" Reference Object: "<<roid<<", Trial Space: "
				<<m_lfeID<<", refDim="<<refDim);
	}
};

} // end namespace ug

#endif