/*
 * grid_function_user_data.h
 *
 *  Created on: 03.07.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__FUNCTION_SPACE__GRID_FUNCTION_USER_DATA__
#define __H__UG__LIB_DISC__FUNCTION_SPACE__GRID_FUNCTION_USER_DATA__

#include "common/common.h"

#include "lib_disc/common/subset_group.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/quadrature/quadrature.h"
#include "lib_disc/local_finite_element/local_shape_function_set.h"
#include "lib_disc/spatial_disc/user_data/std/std_user_data.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"


namespace ug{

template <typename TGridFunction>
class GridFunctionNumberData
	: public StdDependentUserData<GridFunctionNumberData<TGridFunction>,
	  	  	  	  	  	  	  	  number, TGridFunction::dim>
{
	public:
	//	world dimension of grid function
		static const int dim = TGridFunction::dim;

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
			: m_spGridFct(spGridFct)
		{
			this->set_functions(cmp);

		//	get function id of name
			m_fct = spGridFct->fct_id_by_name(cmp);

		//	check that function exists
			if(m_fct >= spGridFct->num_fct())
				UG_THROW("GridFunctionNumberData: Function space does not contain"
						" a function with name " << cmp << ".");

		//	local finite element id
			m_lfeID = spGridFct->local_finite_element_id(m_fct);
		};

		virtual bool continuous() const
		{
			if(m_lfeID.type() == LFEID::LAGRANGE) return true;
			else return false;
		}

		template <int refDim>
		inline void evaluate(number vValue[],
		                     const MathVector<dim> vGlobIP[],
		                     number time, int si,
		                     LocalVector& u,
		                     GeometricObject* elem,
		                     const MathVector<dim> vCornerCoords[],
		                     const MathVector<refDim> vLocIP[],
		                     const size_t nip,
		                     const MathMatrix<refDim, dim>* vJT = NULL) const
		{
		//	reference object id
			const ReferenceObjectID roid = elem->reference_object_id();

		//	get trial space
			try{
			const LocalShapeFunctionSet<refDim>& rTrialSpace =
					LocalShapeFunctionSetProvider::get<refDim>(roid, m_lfeID);

		//	memory for shapes
			std::vector<number> vShape;

		//	loop ips
			for(size_t ip = 0; ip < nip; ++ip)
			{
			//	evaluate at shapes at ip
				rTrialSpace.shapes(vShape, vLocIP[ip]);

			//	get multiindices of element
				std::vector<MultiIndex<2> > ind;
				m_spGridFct->multi_indices(elem, m_fct, ind);

			// 	compute solution at integration point
				vValue[ip] = 0.0;
				for(size_t sh = 0; sh < vShape.size(); ++sh)
				{
					const number valSH = DoFRef(*m_spGridFct, ind[sh]);
					vValue[ip] += valSH * vShape[sh];
				}
			}

			}
			UG_CATCH_THROW("GridFunctionNumberData: Shape Function Set missing for"
							" Reference Object: "<<roid<<", Trial Space: "
							<<m_lfeID<<", refDim="<<refDim);
		}
		
		template <int refDim>
		void eval_deriv(GeometricObject* elem){
			//	reference object id
			const ReferenceObjectID roid = elem->reference_object_id();

			//	get trial space
			try{
				const LocalShapeFunctionSet<refDim>& rTrialSpace =
						LocalShapeFunctionSetProvider::get<refDim>(roid, m_lfeID);

				//	memory for shapes
				std::vector<number> vShape;

				for(size_t s = 0; s < this->num_series(); ++s){
					for(size_t ip = 0; ip < this->num_ip(s); ++ip){
						//	evaluate at shapes at ip
						rTrialSpace.shapes(vShape, this->template local_ip<refDim>(s,ip));

						for(size_t sh = 0; sh < vShape.size(); ++sh)
							this->deriv(s, ip, 0, sh) = vShape[sh];
					}
				}
			}
			UG_CATCH_THROW("GridFunctionNumberData: Shape Function Set missing for"
					" Reference Object: "<<roid<<", Trial Space: "
					<<m_lfeID<<", refDim="<<refDim);

		}

		virtual void compute(LocalVector* u, GeometricObject* elem,
							 const MathVector<dim> vCornerCoords[], bool bDeriv = false){
			const number t = this->time();
			const int si = this->subset();
			for(size_t s = 0; s < this->num_series(); ++s)
				evaluate<dim>(this->values(s), this->ips(s), t, si,
                  *u, elem, NULL, this->template local_ips<dim>(s),
                  this->num_ip(s));

			if(!bDeriv) return;

			switch(elem->base_object_id()){
				case EDGE: eval_deriv<1>(elem); break;
				case FACE: eval_deriv<2>(elem); break;
				case VOLUME: eval_deriv<3>(elem); break;
				default: UG_THROW("Not impl.");
			}
		}
};

template <typename TGridFunction>
class GridFunctionVectorData
	: public StdDependentUserData<GridFunctionVectorData<TGridFunction>,
	  	  	  	  	  	  	  	  MathVector<TGridFunction::dim>, TGridFunction::dim>
{
	public:
	//	world dimension of grid function
		static const int dim = TGridFunction::dim;

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

		virtual bool continuous() const
		{
			for(int i = 0; i < dim; ++i)
				if(m_vlfeID[i].type() != LFEID::LAGRANGE)
					return false;
			return true;
		}

		template <int refDim>
		inline void evaluate(MathVector<dim> vValue[],
							 const MathVector<dim> vGlobIP[],
							 number time, int si,
							 LocalVector& u,
							 GeometricObject* elem,
							 const MathVector<dim> vCornerCoords[],
							 const MathVector<refDim> vLocIP[],
							 const size_t nip,
							 const MathMatrix<refDim, dim>* vJT = NULL) const
		{
		//	reference object id
			const ReferenceObjectID roid = elem->reference_object_id();

		//	memory for shapes
			std::vector<number> vShape;

		//	loop ips
			try{
			for(size_t ip = 0; ip < nip; ++ip)
			{
				for(int d = 0; d < dim; ++d)
				{
					const LocalShapeFunctionSet<refDim>& rTrialSpace =
							LocalShapeFunctionSetProvider::get<refDim>(roid, m_vlfeID[d]);

				//	evaluate at shapes at ip
					rTrialSpace.shapes(vShape, vLocIP[ip]);

				//	get multiindices of element
					std::vector<MultiIndex<2> > ind;
					m_spGridFct->multi_indices(elem, m_vfct[d], ind);

				// 	compute solution at integration point
					vValue[ip][d] = 0.0;
					for(size_t sh = 0; sh < vShape.size(); ++sh)
					{
						const number valSH = DoFRef( *m_spGridFct, ind[sh]);
						vValue[ip][d] += valSH * vShape[sh];
					}
				}
			}

			}UG_CATCH_THROW("GridFunctionVectorData: Reference Object: "
						 <<roid<<", refDim="<<refDim);
		}

		virtual void compute(LocalVector* u, GeometricObject* elem,
							 const MathVector<dim> vCornerCoords[], bool bDeriv = false){
			const number t = this->time();
			const int si = this->subset();
			for(size_t s = 0; s < this->num_series(); ++s)
				evaluate<dim>(this->values(s), this->ips(s), t, si,
                  *u, elem, NULL, this->template local_ips<dim>(s),
                  this->num_ip(s));

			if(bDeriv)
				UG_THROW("Not implemented.");
		};
};


template <typename TGridFunction>
class GridFunctionGradientData
	: public StdDependentUserData<GridFunctionGradientData<TGridFunction> ,
	  	  	  	  	  	  	  	  MathVector<TGridFunction::dim>, TGridFunction::dim>
{
	public:
	//	world dimension of grid function
		static const int dim = TGridFunction::dim;

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

		virtual void compute(LocalVector* u, GeometricObject* elem,
							 const MathVector<dim> vCornerCoords[], bool bDeriv = false){
			UG_THROW("Not implemented.");
		}

		template <int refDim>
		inline void evaluate(MathVector<dim> vValue[],
							 const MathVector<dim> vGlobIP[],
							 number time, int si,
							 LocalVector& u,
							 GeometricObject* elem,
							 const MathVector<dim> vCornerCoords[],
							 const MathVector<refDim> vLocIP[],
							 const size_t nip,
							 const MathMatrix<refDim, dim>* vJT = NULL) const
		{
		//	reference object id
			const ReferenceObjectID roid = elem->reference_object_id();

		//	get reference element mapping by reference object id
			std::vector<MathMatrix<refDim, dim> > vJTTmp(nip);
			if(vJT == NULL){
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
					LocalShapeFunctionSetProvider::get<refDim>(roid, m_lfeID);

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
				std::vector<MultiIndex<2> > ind;
				m_spGridFct->multi_indices(elem, m_fct, ind);

			//	compute grad at ip
				VecSet(locGrad, 0.0);
				for(size_t sh = 0; sh < vLocGrad.size(); ++sh)
				{
					const number valSH = DoFRef( *m_spGridFct, ind[sh]);
					VecScaleAppend(locGrad, valSH, vLocGrad[sh]);
				}

				Inverse(JTInv, vJT[ip]);
				MatVecMult(vValue[ip], JTInv, locGrad);
			}
			}
			UG_CATCH_THROW("GridFunctionNumberData: Shape Function Set missing for"
							" Reference Object: "<<roid<<", Trial Space: "
							<<m_lfeID<<", refDim="<<refDim);
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
		static const int dim = TGridFunction::dim;

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

		virtual void compute(LocalVector* u, GeometricObject* elem,
							 const MathVector<dim> vCornerCoords[], bool bDeriv = false){
			UG_THROW("Not implemented.");
		}

		/**
		 * \param[out] vValue Array of the <tt>nip</tt> gradient components
		 */
		template <int refDim>
		inline void evaluate( number vValue[],
		                      const MathVector<dim> vGlobIP[],
		                      number time, int si,
		                      LocalVector& u,
		                      GeometricObject* elem,
		                      const MathVector<dim> vCornerCoords[],
		                      const MathVector<refDim> vLocIP[],
		                      const size_t nip,
		                      const MathMatrix<refDim, dim>* vJT = NULL ) const
		{
			//	reference object id
			const ReferenceObjectID roid = elem->reference_object_id();

			//	get reference element mapping by reference object id
			std::vector<MathMatrix<refDim, dim> > vJTTmp( nip );
			if( vJT == NULL ) {
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
				  LocalShapeFunctionSetProvider::get<refDim>( roid, m_lfeID );

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
					std::vector<MultiIndex<2> > ind;
					m_spGridFct->multi_indices( elem, m_fct, ind );

					//	compute grad at ip
					VecSet( locGrad, 0.0 );
					for( size_t sh = 0; sh < vLocGrad.size(); ++sh ) {
						const number valSH = DoFRef( *m_spGridFct, ind[sh] );
						VecScaleAppend( locGrad, valSH, vLocGrad[sh] );
					}

					Inverse( JTInv, vJT[ip] );
					MatVecMult( vValueVec[ip], JTInv, locGrad );

					vValue[ip] = vValueVec[ip][m_component];
				}
			}
			UG_CATCH_THROW("GridFunctionNumberData: Shape Function Set missing for"
					" Reference Object: "<<roid<<", Trial Space: "
					<<m_lfeID<<", refDim="<<refDim);
		}
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__FUNCTION_SPACE__GRID_FUNCTION_USER_DATA__ */
