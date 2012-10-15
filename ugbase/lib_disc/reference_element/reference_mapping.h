/*
 * reference_element_mapping.h
 *
 *  Created on: 13.05.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__REFERENCE_ELEMENT__REFERENCE_ELEMENT_MAPPING__
#define __H__UG__LIB_DISC__REFERENCE_ELEMENT__REFERENCE_ELEMENT_MAPPING__

#include <cassert>
#include <iostream>
#include <sstream>
#include "common/common.h"
#include "common/math/ugmath.h"
#include "lib_disc/reference_element/reference_element.h"

namespace ug{

/**
 * This class describes the mapping from a reference element into the real
 * (physical) world. The mapping is initialized by the physical positions of
 * the vertices of the real world element. The order of those points must be
 * given as indicated by the corresponding reference element.
 *
 * Let \f$R\f$ be the reference element and \f$T\f$ be the element. Then, the
 * reference mapping is a mapping:
 * \f[
 * 	\phi:	R \mapsto T
 * \f]
 *
 * \tparam	TRefElem		reference element
 * \tparam	TWorldDim		world dimension
 */
template <typename TRefElem, int TWorldDim>
class ReferenceMapping
{
	public:
	///	world dimension (range space dimension)
		static const int worldDim = TWorldDim;

	///	reference dimension (domain space dimension)
		static const int dim = TRefElem::dim;

	///	flag if mapping is linear (i.e. Jacobian does not depend on x)
		static const bool isLinear = false;

	public:
	///	Default Constructor
		ReferenceMapping();

	///	Constructor setting the corners of the element
		ReferenceMapping(const MathVector<worldDim>* vCornerCoord);

	///	Constructor setting the corners of the element
		ReferenceMapping(const std::vector<MathVector<worldDim> >& vCornerCoord);

	///	refresh mapping for new set of corners
		void update(const MathVector<worldDim>* vCornerCoord);

	///	refresh mapping for new set of corners
		void update(const std::vector<MathVector<worldDim> >& vCornerCoord);

	///	returns if mapping is affine
		bool is_linear() const;

	///	map local coordinate to global coordinate
		void local_to_global(MathVector<worldDim>& globPos,
		                     const MathVector<dim>& locPos) const;

	///	map local coordinate to global coordinate for n local positions
		void local_to_global(MathVector<worldDim>* vGlobPos,
							 const MathVector<dim>* vLocPos, size_t n) const;

	///	map local coordinate to global coordinate for a vector of local positions
		void local_to_global(std::vector<MathVector<worldDim> >& vGlobPos,
							 const std::vector<MathVector<dim> >& vLocPos) const;

	///	map global coordinate to local coordinate
		void global_to_local(MathVector<dim>& locPos,
							 const MathVector<worldDim>& globPos,
							 const size_t maxIter = 1000,
							 const number tol = 1e-10) const;

	///	map global coordinate to local coordinate for n local positions
		void global_to_local(MathVector<dim>* vLocPos,
							 const MathVector<worldDim>* vGlobPos, size_t n,
							 const size_t maxIter = 1000,
							 const number tol = 1e-10) const;

	///	map global coordinate to local coordinate for a vector of local positions
		void global_to_local(std::vector<MathVector<dim> >& vLocPos,
							 const std::vector<MathVector<worldDim> >& vGlobPos,
							 const size_t maxIter = 1000,
							 const number tol = 1e-10) const;

	///	returns jacobian
		void jacobian(MathMatrix<worldDim, dim>& J,
		              const MathVector<dim>& locPos) const;

	///	returns jacobian for n local positions
		void jacobian(MathMatrix<worldDim, dim>* vJ,
					  const MathVector<dim>* vLocPos, size_t n) const;

	///	returns jacobian for a vector of local positions
		void jacobian(std::vector<MathMatrix<worldDim, dim> >& J,
					  const std::vector<MathVector<dim> >& vLocPos) const;

	///	returns transposed of jacobian
		void jacobian_transposed(MathMatrix<dim, worldDim>& JT,
		                         const MathVector<dim>& locPos) const;

	///	returns transposed of jacobian for n local positions
		void jacobian_transposed(MathMatrix<dim, worldDim>* vJT,
		                         const MathVector<dim>* vLocPos, size_t n) const;

	///	returns transposed of jacobian for a vector of positions
		void jacobian_transposed(std::vector<MathMatrix<dim, worldDim> >& vJT,
								 const std::vector<MathVector<dim> >& vLocPos) const;

	///	returns transposed of the inverse of the jacobian and sqrt of gram determinante
		number jacobian_transposed_inverse(MathMatrix<worldDim, dim>& JTInv,
		                                   const MathVector<dim>& locPos) const;

	///	returns transposed of the inverse of the jacobian for n local positions
		void jacobian_transposed_inverse(MathMatrix<worldDim, dim>* vJTInv,
		                                 const MathVector<dim>* vLocPos, size_t n) const;

	///	returns transposed of the inverse of the jacobian for n local positions
		void jacobian_transposed_inverse(MathMatrix<worldDim, dim>* vJTInv,
		                                 number* vDet,
										 const MathVector<dim>* vLocPos, size_t n) const;

	///	returns transposed of the inverse of the jacobian for a vector of positions
		void jacobian_transposed_inverse(std::vector<MathMatrix<worldDim, dim> >& vJTInv,
										 const std::vector<MathVector<dim> >& vLocPos) const;

	///	returns transposed of the inverse of the jacobian for a vector of positions
		void jacobian_transposed_inverse(std::vector<MathMatrix<worldDim, dim> >& vJTInv,
										 std::vector<number>& vDet,
										 const std::vector<MathVector<dim> >& vLocPos) const;

	///	returns the determinate of the jacobian
		number sqrt_gram_det(const MathVector<dim>& locPos) const;

	///	returns the determinate of the jacobian for n local positions
		void sqrt_gram_det(number* vDet, const MathVector<dim>* vLocPos, size_t n) const;

	///	returns the determinate of the jacobian for a vector of local positions
		void sqrt_gram_det(std::vector<number> vDet,
						  const std::vector<MathVector<dim> >& vLocPos) const;
};


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Concrete Reference Mappings
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/// Base class for Reference mappings helping to implement interface
template <int dim, int worldDim, bool isLinear, typename TImpl>
class BaseReferenceMapping
{
	public:
	///	returns if mapping is affine
		bool is_linear() const {return isLinear;}

	///	map local coordinate to global coordinate for n local positions
		void local_to_global(MathVector<worldDim>* vGlobPos,
							 const MathVector<dim>* vLocPos, size_t n) const
		{
			for(size_t ip = 0; ip < n; ++ip)
				getImpl().local_to_global(vGlobPos[ip], vLocPos[ip]);
		}

	///	map local coordinate to global coordinate for a vector of local positions
		void local_to_global(std::vector<MathVector<worldDim> >& vGlobPos,
							 const std::vector<MathVector<dim> >& vLocPos) const
		{
			const size_t n = vLocPos.size();
			vGlobPos.resize(n);
			local_to_global(&vGlobPos[0], &vLocPos[0], n);
		}

	///	map global coordinate to local coordinate
		void global_to_local(MathVector<dim>& locPos,
							 const MathVector<worldDim>& globPos,
							 const size_t maxIter = 1000,
							 const number tol = 1e-10) const
		{
			MathMatrix<worldDim, dim> J;
			MathMatrix<dim, worldDim> JInv;
			MathVector<worldDim> dist, compGlobPos;

			for (size_t i = 0; i < maxIter; ++i) {

			//	f(x) := \phi(x) - x_{glob}
				getImpl().local_to_global(compGlobPos, locPos);
				VecSubtract(dist, compGlobPos, globPos);

			//	check if tol reached
				if(fabs(VecTwoNorm(dist)) < tol) return;

			//	compute jacobian df/dx = d \phi(x) / dx =: J
				getImpl().jacobian(J, locPos);

			//	solve x -= J^{-1} f
				LeftInverse(JInv, J);
				MatVecScaleMultAppend(locPos, -1.0, JInv, dist);
			}

			UG_THROW("ReferenceMapping::global_to_local: Newton method did not"
					" reach a tolerance "<<tol<<" after "<<maxIter<<" steps.");
		}

	///	map global coordinate to local coordinate for n local positions
		void global_to_local(MathVector<dim>* vLocPos,
							 const MathVector<worldDim>* vGlobPos, size_t n,
							 const size_t maxIter = 1000,
							 const number tol = 1e-10) const
		{
			if(isLinear){
				if(n == 0) return;

				MathMatrix<worldDim, dim> J;
				MathMatrix<dim, worldDim> JInv;
				MathVector<worldDim> dist, compGlobPos;

			//	compute jacobian df/dx = d \phi(x) / dx =: J
				getImpl().jacobian(J, vLocPos[0]);

			//	solve x -= J^{-1} f
				LeftInverse(JInv, J);

				for(size_t ip = 0; ip < n; ++ip)
				{
				//	f(x) := \phi(x) - x_{glob}
					getImpl().local_to_global(compGlobPos, vLocPos[ip]);
					VecSubtract(dist, compGlobPos, vGlobPos[ip]);

					MatVecScaleMultAppend(vLocPos[ip], -1.0, JInv, dist);
				}
			}
			else{
				for(size_t ip = 0; ip < n; ++ip)
					getImpl().global_to_local(vLocPos[ip], vGlobPos[ip], maxIter, tol);
			}
		}

	///	map global coordinate to local coordinate for a vector of local positions
		void global_to_local(std::vector<MathVector<dim> >& vLocPos,
							 const std::vector<MathVector<worldDim> >& vGlobPos,
							 const size_t maxIter = 1000,
							 const number tol = 1e-10) const
		{
			const size_t n = vGlobPos.size();
			vLocPos.resize(n);
			global_to_local(&vLocPos[0], &vGlobPos[0], n, maxIter, tol);
		}

	///	returns jacobian
		void jacobian(MathMatrix<worldDim, dim>& J,
					  const MathVector<dim>& locPos) const
		{
			MathMatrix<dim, worldDim> JT;
			getImpl().jacobian_transposed(JT, locPos);
			Transpose(J, JT);
		}

	///	returns jacobian for n local positions
		void jacobian(MathMatrix<worldDim, dim>* vJ,
					  const MathVector<dim>* vLocPos, size_t n) const
		{
			if(isLinear){
				if(n == 0) return;
				getImpl().jacobian(vJ[0], vLocPos[0]);
				for(size_t ip = 1; ip < n; ++ip) vJ[ip] = vJ[0];
			}
			else {
				for(size_t ip = 0; ip < n; ++ip)
					getImpl().jacobian(vJ[ip], vLocPos[ip]);
			}
		}

	///	returns jacobian for a vector of local positions
		void jacobian(std::vector<MathMatrix<worldDim, dim> >& vJ,
					  const std::vector<MathVector<dim> >& vLocPos) const
		{
			const size_t n = vLocPos.size();
			vJ.resize(n);
			jacobian(&vJ[0], &vLocPos[0], n);
		}


	///	returns transposed of jacobian for n local positions
		void jacobian_transposed(MathMatrix<dim, worldDim>* vJT,
								 const MathVector<dim>* vLocPos, size_t n) const
		{
			if(isLinear){
				if(n == 0) return;
				getImpl().jacobian_transposed(vJT[0], vLocPos[0]);
				for(size_t ip = 1; ip < n; ++ip) vJT[ip] = vJT[0];
			}
			else {
				for(size_t ip = 0; ip < n; ++ip)
					getImpl().jacobian_transposed(vJT[ip], vLocPos[ip]);
			}
		}

	///	returns transposed of jacobian for a vector of positions
		void jacobian_transposed(std::vector<MathMatrix<dim, worldDim> >& vJT,
								 const std::vector<MathVector<dim> >& vLocPos) const
		{
			const size_t n = vLocPos.size();
			vJT.resize(n);
			jacobian_transposed(&vJT[0], &vLocPos[0], n);
		}

	///	returns transposed of the inverse of the jacobian
		number jacobian_transposed_inverse(MathMatrix<worldDim, dim>& JTInv,
										 const MathVector<dim>& locPos) const
		{
			MathMatrix<dim, worldDim> JT;
			getImpl().jacobian_transposed(JT, locPos);
			return RightInverse(JTInv, JT);
		}

	///	returns transposed of the inverse of the jacobian for n local positions
		void jacobian_transposed_inverse(MathMatrix<worldDim, dim>* vJTInv,
		                                 number* vDet,
										 const MathVector<dim>* vLocPos, size_t n) const
		{
			if(isLinear){
				if(n == 0) return;
				vDet[0] = getImpl().jacobian_transposed_inverse(vJTInv[0], vLocPos[0]);
				for(size_t ip = 1; ip < n; ++ip) vJTInv[ip] = vJTInv[0];
				for(size_t ip = 1; ip < n; ++ip) vDet[ip] = vDet[0];
			}
			else {
				for(size_t ip = 0; ip < n; ++ip)
					vDet[ip] = getImpl().jacobian_transposed_inverse(vJTInv[ip], vLocPos[ip]);
			}
		}

	///	returns transposed of the inverse of the jacobian for n local positions
		void jacobian_transposed_inverse(MathMatrix<worldDim, dim>* vJTInv,
										 const MathVector<dim>* vLocPos, size_t n) const
		{
			if(isLinear){
				if(n == 0) return;
				getImpl().jacobian_transposed_inverse(vJTInv[0], vLocPos[0]);
				for(size_t ip = 1; ip < n; ++ip) vJTInv[ip] = vJTInv[0];
			}
			else {
				for(size_t ip = 0; ip < n; ++ip)
					getImpl().jacobian_transposed_inverse(vJTInv[ip], vLocPos[ip]);
			}
		}

	///	returns transposed of the inverse of the jacobian for a vector of positions
		void jacobian_transposed_inverse(std::vector<MathMatrix<worldDim, dim> >& vJTInv,
		                                 std::vector<number>& vDet,
										 const std::vector<MathVector<dim> >& vLocPos) const
		{
			const size_t n = vLocPos.size();
			vJTInv.resize(n); vDet.resize(n);
			jacobian_transposed_inverse(&vJTInv[0], &vDet[0], &vLocPos[0], n);
		}

	///	returns transposed of the inverse of the jacobian for a vector of positions
		void jacobian_transposed_inverse(std::vector<MathMatrix<worldDim, dim> >& vJTInv,
										 const std::vector<MathVector<dim> >& vLocPos) const
		{
			const size_t n = vLocPos.size();
			vJTInv.resize(n);
			jacobian_transposed_inverse(&vJTInv[0], &vLocPos[0], n);
		}

	///	returns the determinate of the jacobian
		number sqrt_gram_det(const MathVector<dim>& locPos) const
		{
			MathMatrix<dim, worldDim> JT;
			getImpl().jacobian_transposed(JT, locPos);
			return SqrtGramDeterminant(JT);
		}

	///	returns the determinate of the jacobian for n local positions
		void sqrt_gram_det(number* vDet, const MathVector<dim>* vLocPos, size_t n) const
		{
			if(isLinear){
				if(n == 0) return;
				vDet[0] = sqrt_gram_det(vLocPos[0]);
				for(size_t ip = 1; ip < n; ++ip) vDet[ip] = vDet[0];
			}
			else {
				for(size_t ip = 0; ip < n; ++ip)
					vDet[ip] = sqrt_gram_det(vLocPos[ip]);
			}
		}

	///	returns the determinate of the jacobian for a vector of local positions
		void sqrt_gram_det(std::vector<number>& vDet,
		                  const std::vector<MathVector<dim> >& vLocPos) const
		{
			const size_t n = vLocPos.size();
			vDet.resize(n);
			sqrt_gram_det(&vDet[0], &vLocPos[0], n);
		}

	protected:
	///	access to implementation
		TImpl& getImpl() {return static_cast<TImpl&>(*this);}

	///	const access to implementation
		const TImpl& getImpl() const {return static_cast<const TImpl&>(*this);}
};

///////////////////////////////////////////////////////////////////////////////
// Reference Mapping Vertex
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Reference Mapping Edge
///////////////////////////////////////////////////////////////////////////////
template <int TWorldDim>
class ReferenceMapping<ReferenceEdge, TWorldDim>
	: public BaseReferenceMapping<ReferenceEdge::dim, TWorldDim, true,
	  	  	  	  	  	  	  	  ReferenceMapping<ReferenceEdge, TWorldDim> >
{
	public:
	///	world dimension
		static const int worldDim = TWorldDim;

	///	reference dimension
		static const int dim = ReferenceEdge::dim;

	///	flag if mapping is linear (i.e. Jacobian does not depend on x)
		static const bool isLinear = true;

	public:
		typedef BaseReferenceMapping<ReferenceEdge::dim, TWorldDim, true,
	  	  	  	  	  ReferenceMapping<ReferenceEdge, TWorldDim> > base_type;
		using base_type::local_to_global;
		using base_type::jacobian;
		using base_type::jacobian_transposed;
		using base_type::jacobian_transposed_inverse;
		using base_type::sqrt_gram_det;

	public:
	///	Default Constructor
		ReferenceMapping() {}

	///	Constructor setting the corners
	/// \{
		ReferenceMapping(const MathVector<worldDim>* vCornerCoord) {update(vCornerCoord);}
		ReferenceMapping(const std::vector<MathVector<worldDim> >& vCornerCoord) {update(vCornerCoord);}
	/// \}

	///	refresh mapping for new set of corners
		void update(const std::vector<MathVector<worldDim> >& vCornerCoord)
		{
			UG_ASSERT((int)vCornerCoord.size() >= ReferenceEdge::numCorners,
			          "ReferenceMapping: to few Corner Coordinates.");
			update(&vCornerCoord[0]);
		}

	///	update the mapping for a new set of corners
		void update(const MathVector<worldDim>* vCornerCoord)
		{
			co0 = vCornerCoord[0];
			VecSubtract(a10, vCornerCoord[1], vCornerCoord[0]);
		}

	///	map local coordinate to global coordinate
		void local_to_global(MathVector<worldDim>& globPos,
							 const MathVector<dim>& locPos) const
		{
			VecScaleAdd(globPos, 1.0, co0, locPos[0], a10);
		}

	///	returns transposed of jacobian
		void jacobian_transposed(MathMatrix<dim, worldDim>& JT,
								 const MathVector<dim>& locPos) const
		{
			for(int i = 0; i < worldDim; ++i) JT(0,i) = a10[i];
		}

	private:
		MathVector<worldDim> co0, a10;
};

///////////////////////////////////////////////////////////////////////////////
// Reference Mapping Triangle
///////////////////////////////////////////////////////////////////////////////
template <int TWorldDim>
class ReferenceMapping<ReferenceTriangle, TWorldDim>
	: public BaseReferenceMapping<ReferenceTriangle::dim, TWorldDim, true,
	  	  	  	  	  	  	  	  ReferenceMapping<ReferenceTriangle, TWorldDim> >
{
	public:
	///	world dimension
		static const int worldDim = TWorldDim;

	///	reference dimension
		static const int dim = ReferenceTriangle::dim;

	///	flag if mapping is linear (i.e. Jacobian does not depend on x)
		static const bool isLinear = true;

	public:
		typedef BaseReferenceMapping<ReferenceTriangle::dim, TWorldDim, true,
	  	  	  	  ReferenceMapping<ReferenceTriangle, TWorldDim> > base_type;
		using base_type::local_to_global;
		using base_type::jacobian;
		using base_type::jacobian_transposed;
		using base_type::jacobian_transposed_inverse;
		using base_type::sqrt_gram_det;

	public:
	///	Default Constructor
		ReferenceMapping() {}

	///	Constructor setting the corners
	/// \{
		ReferenceMapping(const MathVector<worldDim>* vCornerCoord) {update(vCornerCoord);}
		ReferenceMapping(const std::vector<MathVector<worldDim> >& vCornerCoord) {update(vCornerCoord);}
	/// \}

	///	refresh mapping for new set of corners
		void update(const std::vector<MathVector<worldDim> >& vCornerCoord)
		{
			UG_ASSERT((int)vCornerCoord.size() >= ReferenceTriangle::numCorners,
			          "ReferenceMapping: to few Corner Coordinates.");
			update(&vCornerCoord[0]);
		}

	///	update the mapping for a new set of corners
		void update(const MathVector<worldDim>* vCornerCoord)
		{
			co0 = vCornerCoord[0];
			VecSubtract(a10, vCornerCoord[1], vCornerCoord[0]);
			VecSubtract(a20, vCornerCoord[2], vCornerCoord[0]);
		}

	///	map local coordinate to global coordinate
		void local_to_global(MathVector<worldDim>& globPos,
							 const MathVector<dim>& locPos) const
		{
			VecScaleAdd(globPos, 1.0, co0, locPos[0], a10, locPos[1], a20);
		}

	///	returns transposed of jacobian
		void jacobian_transposed(MathMatrix<dim, worldDim>& JT,
								 const MathVector<dim>& locPos) const
		{
			for(int i = 0; i < worldDim; ++i) JT(0, i) = a10[i];
			for(int i = 0; i < worldDim; ++i) JT(1, i) = a20[i];
		}

	private:
		MathVector<worldDim> co0, a10, a20;

};
///////////////////////////////////////////////////////////////////////////////
// Reference Mapping Quadrilateral
///////////////////////////////////////////////////////////////////////////////


template <int TWorldDim>
class ReferenceMapping<ReferenceQuadrilateral, TWorldDim>
	: public BaseReferenceMapping<ReferenceQuadrilateral::dim, TWorldDim, false,
								  ReferenceMapping<ReferenceQuadrilateral, TWorldDim> >
{
	public:
	///	world dimension
		static const int worldDim = TWorldDim;

	///	reference dimension
		static const int dim = ReferenceQuadrilateral::dim;

	///	flag if mapping is linear (i.e. Jacobian does not depend on x)
		static const bool isLinear = false;

	public:
		typedef BaseReferenceMapping<ReferenceQuadrilateral::dim, TWorldDim, false,
				  ReferenceMapping<ReferenceQuadrilateral, TWorldDim> > base_type;
		using base_type::local_to_global;
		using base_type::jacobian;
		using base_type::jacobian_transposed;
		using base_type::jacobian_transposed_inverse;
		using base_type::sqrt_gram_det;

	public:
	///	Default Constructor
		ReferenceMapping() {}

	///	Constructor setting the corners
	/// \{
		ReferenceMapping(const MathVector<worldDim>* vCornerCoord) {update(vCornerCoord);}
		ReferenceMapping(const std::vector<MathVector<worldDim> >& vCornerCoord) {update(vCornerCoord);}
	/// \}

	///	refresh mapping for new set of corners
		void update(const std::vector<MathVector<worldDim> >& vCornerCoord)
		{
			UG_ASSERT((int)vCornerCoord.size() >= ReferenceQuadrilateral::numCorners,
			          "ReferenceMapping: to few Corner Coordinates.");
			update(&vCornerCoord[0]);
		}

	///	update the mapping for a new set of corners
		void update(const MathVector<worldDim>* vCornerCoord)
		{
			for(int co = 0; co < ReferenceQuadrilateral::numCorners; ++co)
				x[co] = vCornerCoord[co];
		}

	///	map local coordinate to global coordinate
		void local_to_global(MathVector<worldDim>& globPos,
							 const MathVector<dim>& locPos) const
		{
			const number a = (1.-locPos[0]);
			const number b = (1.-locPos[1]);

			VecScaleAdd(globPos,        		a*b, x[0],
			            				locPos[0]*b, x[1],
								locPos[0]*locPos[1], x[2],
										a*locPos[1], x[3]);
		}

	///	returns transposed of jacobian
		void jacobian_transposed(MathMatrix<dim, worldDim>& JT,
								 const MathVector<dim>& locPos) const
		{
			const number a = 1. - locPos[1];
			const number b = 1. - locPos[0];

			for(int i = 0; i < worldDim; ++i)
				JT(0, i) = a*(x[1][i] - x[0][i]) + locPos[1]*(x[2][i] - x[3][i]);

			for(int i = 0; i < worldDim; ++i)
				JT(1, i) = b*(x[3][i] - x[0][i]) + locPos[0]*(x[2][i] - x[1][i]);
		}

	private:
		MathVector<worldDim> x[ReferenceQuadrilateral::numCorners];
};

///////////////////////////////////////////////////////////////////////////////
// Reference Mapping Tetrahedron
///////////////////////////////////////////////////////////////////////////////

template <int TWorldDim>
class ReferenceMapping<ReferenceTetrahedron, TWorldDim>
	: public BaseReferenceMapping<ReferenceTetrahedron::dim, TWorldDim, true,
								  ReferenceMapping<ReferenceTetrahedron, TWorldDim> >
{
	public:
	///	world dimension
		static const int worldDim = TWorldDim;

	///	reference dimension
		static const int dim = ReferenceTetrahedron::dim;

	///	flag if mapping is linear (i.e. Jacobian does not depend on x)
		static const bool isLinear = true;

	public:
		typedef BaseReferenceMapping<ReferenceTetrahedron::dim, TWorldDim, true,
				  ReferenceMapping<ReferenceTetrahedron, TWorldDim> > base_type;
		using base_type::local_to_global;
		using base_type::jacobian;
		using base_type::jacobian_transposed;
		using base_type::jacobian_transposed_inverse;
		using base_type::sqrt_gram_det;

	public:
	///	Default Constructor
		ReferenceMapping() {}

	///	Constructor setting the corners
	/// \{
		ReferenceMapping(const MathVector<worldDim>* vCornerCoord) {update(vCornerCoord);}
		ReferenceMapping(const std::vector<MathVector<worldDim> >& vCornerCoord) {update(vCornerCoord);}
	/// \}

	///	refresh mapping for new set of corners
		void update(const std::vector<MathVector<worldDim> >& vCornerCoord)
		{
			UG_ASSERT((int)vCornerCoord.size() >= ReferenceTetrahedron::numCorners,
			          "ReferenceMapping: to few Corner Coordinates.");
			update(&vCornerCoord[0]);
		}

	///	update the mapping for a new set of corners
		void update(const MathVector<worldDim>* vCornerCoord)
		{
			co0 = vCornerCoord[0];
			VecSubtract(a10, vCornerCoord[1], vCornerCoord[0]);
			VecSubtract(a20, vCornerCoord[2], vCornerCoord[0]);
			VecSubtract(a30, vCornerCoord[3], vCornerCoord[0]);
		}

	///	map local coordinate to global coordinate
		void local_to_global(MathVector<worldDim>& globPos,
							 const MathVector<dim>& locPos) const
		{
			VecScaleAdd(globPos,       1.0, co0,
			            		 locPos[0], a10,
								 locPos[1], a20,
								 locPos[2], a30);
		}

	///	returns transposed of jacobian
		void jacobian_transposed(MathMatrix<dim, worldDim>& JT,
								 const MathVector<dim>& locPos) const
		{
			for(int i = 0; i < worldDim; ++i) JT[0][i] = a10[i];
			for(int i = 0; i < worldDim; ++i) JT[1][i] = a20[i];
			for(int i = 0; i < worldDim; ++i) JT[2][i] = a30[i];
		}

	private:
		MathVector<worldDim> co0, a10, a20, a30;
};

///////////////////////////////////////////////////////////////////////////////
// Reference Mapping Pyramid
///////////////////////////////////////////////////////////////////////////////
template <int TWorldDim>
class ReferenceMapping<ReferencePyramid, TWorldDim>
	: public BaseReferenceMapping<ReferencePyramid::dim, TWorldDim, false,
								  ReferenceMapping<ReferencePyramid, TWorldDim> >
{
	public:
	///	world dimension
		static const int worldDim = TWorldDim;

	///	reference dimension
		static const int dim = ReferencePyramid::dim;

	///	flag if mapping is linear (i.e. Jacobian does not depend on x)
		static const bool isLinear = false;

	public:
		typedef BaseReferenceMapping<ReferencePyramid::dim, TWorldDim, false,
				  ReferenceMapping<ReferencePyramid, TWorldDim> > base_type;
		using base_type::local_to_global;
		using base_type::jacobian;
		using base_type::jacobian_transposed;
		using base_type::jacobian_transposed_inverse;
		using base_type::sqrt_gram_det;

	public:
	///	Default Constructor
		ReferenceMapping() {}

	///	Constructor setting the corners
	/// \{
		ReferenceMapping(const MathVector<worldDim>* vCornerCoord) {update(vCornerCoord);}
		ReferenceMapping(const std::vector<MathVector<worldDim> >& vCornerCoord) {update(vCornerCoord);}
	/// \}

	///	refresh mapping for new set of corners
		void update(const std::vector<MathVector<worldDim> >& vCornerCoord)
		{
			UG_ASSERT((int)vCornerCoord.size() >= ReferencePyramid::numCorners,
			          "ReferenceMapping: to few Corner Coordinates.");
			update(&vCornerCoord[0]);
		}

	///	update the mapping for a new set of corners
		void update(const MathVector<worldDim>* vCornerCoord)
		{
			for(int co = 0; co < ReferencePyramid::numCorners; ++co)
				x[co] = vCornerCoord[co];
		}

	///	map local coordinate to global coordinate
		void local_to_global(MathVector<worldDim>& globPos,
							 const MathVector<dim>& locPos) const
		{
			const number a = 1.0 - locPos[0];
			const number b = 1.0 - locPos[1];
			if (locPos[0] > locPos[1])
			{
				const number a0 = a * b - locPos[2] * b;
				const number a1 = locPos[0] * b - locPos[2]*locPos[1];
				const number a2 = locPos[0] * locPos[1] + locPos[2]*locPos[1];
				const number a3 = a * locPos[1] - locPos[2] * locPos[1];

				for(int d = 0; d < worldDim; ++d)
					globPos[d] = a0*x[0][d]+a1*x[1][d]+a2*x[2][d]
								+a3*x[3][d]+locPos[2]*x[4][d];
			}
			else
			{
				const number a0 = a * b - locPos[2] * a;
				const number a1 = locPos[0] * b - locPos[2]*locPos[0];
				const number a2 = locPos[0] * locPos[1] + locPos[2]*locPos[0];
				const number a3 = a * locPos[1] - locPos[2] * locPos[0];
				for(int d = 0; d < worldDim; ++d)
					globPos[d] = a0*x[0][d]+a1*x[1][d]+
								 a2*x[2][d]+a3*x[3][d]+locPos[2]*x[4][d];
			}
		}

	///	returns transposed of jacobian
		void jacobian_transposed(MathMatrix<dim, worldDim>& JT,
								 const MathVector<dim>& locPos) const
	   {
			number a[3];
			for(int d = 0; d < worldDim; ++d)
				a[d] = x[0][d]-x[1][d]+x[2][d]-x[3][d];

			if (locPos[0] > locPos[1])
			{
				for(int d = 0; d < worldDim; ++d)
					JT(0,d) = x[1][d]-x[0][d]+locPos[1]*a[d];

				for(int d = 0; d < worldDim; ++d)
					JT(1,d) = x[3][d]-x[0][d]+(locPos[0]+locPos[2])*a[d];

				for(int d = 0; d < worldDim; ++d)
					JT(2,d) = x[4][d]-x[0][d]+locPos[1]*a[d];
			}
			else
			{
				for(int d = 0; d < worldDim; ++d)
					JT(0,d) = x[1][d]-x[0][d]+(locPos[1]+locPos[2])*a[d];

				for(int d = 0; d < worldDim; ++d)
					JT(1,d) = x[3][d]-x[0][d]+locPos[0]*a[d];

				for(int d = 0; d < worldDim; ++d)
					JT(2,d) = x[4][d]-x[0][d]+locPos[0]*a[d];
			}
		}

	private:
		MathVector<worldDim> x[ReferencePyramid::numCorners];
};

///////////////////////////////////////////////////////////////////////////////
// Reference Mapping Prism
///////////////////////////////////////////////////////////////////////////////

template <int TWorldDim>
class ReferenceMapping<ReferencePrism, TWorldDim>
	: public BaseReferenceMapping<ReferencePrism::dim, TWorldDim, false,
								  ReferenceMapping<ReferencePrism, TWorldDim> >
{
	public:
	///	world dimension
		static const int worldDim = TWorldDim;

	///	reference dimension
		static const int dim = ReferencePrism::dim;

	///	flag if mapping is linear (i.e. Jacobian does not depend on x)
		static const bool isLinear = false;

	public:
		typedef BaseReferenceMapping<ReferencePrism::dim, TWorldDim, false,
				  ReferenceMapping<ReferencePrism, TWorldDim> > base_type;
		using base_type::local_to_global;
		using base_type::jacobian;
		using base_type::jacobian_transposed;
		using base_type::jacobian_transposed_inverse;
		using base_type::sqrt_gram_det;

	public:
	///	Default Constructor
		ReferenceMapping() {}

	///	Constructor setting the corners
	/// \{
		ReferenceMapping(const MathVector<worldDim>* vCornerCoord) {update(vCornerCoord);}
		ReferenceMapping(const std::vector<MathVector<worldDim> >& vCornerCoord) {update(vCornerCoord);}
	/// \}

	///	refresh mapping for new set of corners
		void update(const std::vector<MathVector<worldDim> >& vCornerCoord)
		{
			UG_ASSERT((int)vCornerCoord.size() >= ReferencePrism::numCorners,
						"ReferenceMapping: to few Corner Coordinates.");
			update(&vCornerCoord[0]);
		}

	///	update the mapping for a new set of corners
		void update(const MathVector<worldDim>* vCornerCoord)
		{
			for(int co = 0; co < ReferencePrism::numCorners; ++co)
				x[co] = vCornerCoord[co];
		}

	///	map local coordinate to global coordinate
		void local_to_global(MathVector<worldDim>& globPos,
							 const MathVector<dim>& locPos) const
		{
			const number a = 1.0 - locPos[0] - locPos[1];
			const number b = 1.0 - locPos[2];
			const number a0 = a * b;
			const number a1 = locPos[0] * b;
			const number a2 = locPos[1] * b;
			const number a3 = a * locPos[2];
			const number a4 = locPos[0] * locPos[2];
			const number a5 = locPos[1] * locPos[2];

			for(int d = 0; d < worldDim; ++d)
				globPos[d] = a0*x[0][d]+a1*x[1][d]+a2*x[2][d]+a3*x[3][d]+a4*x[4][d]+a5*x[5][d];
		}

	///	returns transposed of jacobian
		void jacobian_transposed(MathMatrix<dim, worldDim>& JT,
								 const MathVector<dim>& locPos) const
	   {
	        number a[worldDim], b[worldDim];

			for(int d = 0; d < worldDim; ++d)
	          a[d] = x[0][d]-x[1][d]-x[3][d]+x[4][d];

			for(int d = 0; d < worldDim; ++d)
			  b[d] = x[0][d]-x[2][d]-x[3][d]+x[5][d];

			for(int d = 0; d < worldDim; ++d){
	          JT(0,d) = x[1][d]-x[0][d]+locPos[2]*a[d];
	          JT(1,d) = x[2][d]-x[0][d]+locPos[2]*b[d];
	          JT(2,d) = x[3][d]-x[0][d]+locPos[0]*a[d]+locPos[1]*b[d];
			}
		}

	private:
		MathVector<worldDim> x[ReferencePrism::numCorners];
};


///////////////////////////////////////////////////////////////////////////////
// Reference Mapping Hexahedron
///////////////////////////////////////////////////////////////////////////////

template <int TWorldDim>
class ReferenceMapping<ReferenceHexahedron, TWorldDim>
	: public BaseReferenceMapping<ReferenceHexahedron::dim, TWorldDim, false,
								  ReferenceMapping<ReferenceHexahedron, TWorldDim> >
{
	public:
	///	world dimension
		static const int worldDim = TWorldDim;

	///	reference dimension
		static const int dim = ReferenceHexahedron::dim;

	///	flag if mapping is linear (i.e. Jacobian does not depend on x)
		static const bool isLinear = false;

	public:
		typedef BaseReferenceMapping<ReferenceHexahedron::dim, TWorldDim, false,
				  ReferenceMapping<ReferenceHexahedron, TWorldDim> > base_type;
		using base_type::local_to_global;
		using base_type::jacobian;
		using base_type::jacobian_transposed;
		using base_type::jacobian_transposed_inverse;
		using base_type::sqrt_gram_det;

	public:
	///	Default Constructor
		ReferenceMapping() {}

	///	Constructor setting the corners
	/// \{
		ReferenceMapping(const MathVector<worldDim>* vCornerCoord) {update(vCornerCoord);}
		ReferenceMapping(const std::vector<MathVector<worldDim> >& vCornerCoord) {update(vCornerCoord);}
	/// \}

	///	refresh mapping for new set of corners
		void update(const std::vector<MathVector<worldDim> >& vCornerCoord)
		{
			UG_ASSERT((int)vCornerCoord.size() >= ReferenceHexahedron::numCorners,
					 "ReferenceMapping: to few Corner Coordinates.");
			update(&vCornerCoord[0]);
		}

	///	update the mapping for a new set of corners
		void update(const MathVector<worldDim>* vCornerCoord)
		{
			for(int co = 0; co < ReferenceHexahedron::numCorners; ++co)
				x[co] = vCornerCoord[co];
		}

	///	map local coordinate to global coordinate
		void local_to_global(MathVector<worldDim>& globPos,
							 const MathVector<dim>& locPos) const
		{
			const number a = 1.0 - locPos[0];
			const number b = 1.0 - locPos[1];
			const number c = 1.0 - locPos[2];
			const number a0 = a * b * c;
			const number a1 = locPos[0] * b * c;
			const number a2 = locPos[0] * locPos[1] * c;
			const number a3 = a * locPos[1] * c;
			const number a4 = a * b * locPos[2];
			const number a5 = locPos[0] * b * locPos[2];
			const number a6 = locPos[0] * locPos[1] * locPos[2];
			const number a7 = a * locPos[1] * locPos[2];

			for(int d = 0; d < worldDim; ++d){
				globPos[d] = a0*x[0][d]+a1*x[1][d]+a2*x[2][d]+a3*x[3][d]+
							 a4*x[4][d]+a5*x[5][d]+a6*x[6][d]+a7*x[7][d];
			}
		}

	///	returns transposed of jacobian
		void jacobian_transposed(MathMatrix<dim, worldDim>& JT,
								 const MathVector<dim>& locPos) const
	   {
			number a0,a1,a2,a3;
			const number a = 1.0 - locPos[0];
			const number b = 1.0 - locPos[1];
			const number c = 1.0 - locPos[2];
			a0 = b * c;
			a1 = locPos[1] * c;
			a2 = locPos[1] * locPos[2];
			a3 = b * locPos[2];
			for(int d = 0; d < worldDim; ++d)
				JT(0,d) = a0*(x[1][d]-x[0][d])+a1*(x[2][d]-x[3][d])
						+ a2*(x[6][d]-x[7][d])+a3*(x[5][d]-x[4][d]);

	        a0 = a * c;
	        a1 = locPos[0] * c;
	        a2 = locPos[0] * locPos[2];
	        a3 = a * locPos[2];
			for(int d = 0; d < worldDim; ++d)
				JT(1,d) = a0*(x[3][d]-x[0][d])+a1*(x[2][d]-x[1][d])
						+ a2*(x[6][d]-x[5][d])+a3*(x[7][d]-x[4][d]);

	        a0 = a * b;
	        a1 = locPos[0] * b;
	        a2 = locPos[0] * locPos[1];
	        a3 = a * locPos[1];
			for(int d = 0; d < worldDim; ++d)
				JT(2,d) = a0*(x[4][d]-x[0][d])+a1*(x[5][d]-x[1][d])
						+ a2*(x[6][d]-x[2][d])+a3*(x[7][d]-x[3][d]);
		}

	private:
		MathVector<worldDim> x[ReferenceHexahedron::numCorners];
};


} // end namespace ug

#endif /* __H__UG__LIB_DISC__REFERENCE_ELEMENT__REFERENCE_ELEMENT_MAPPING__ */
