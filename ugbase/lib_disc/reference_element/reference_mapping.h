/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__REFERENCE_ELEMENT__REFERENCE_ELEMENT_MAPPING__
#define __H__UG__LIB_DISC__REFERENCE_ELEMENT__REFERENCE_ELEMENT_MAPPING__

#include <cassert>
#include <iostream>
#include <sstream>
#include "common/common.h"
#include "common/math/ugmath.h"
#include "lib_disc/reference_element/reference_element.h"

namespace ug{

extern DebugID DID_REFERENCE_MAPPING;
extern DebugID DID_REFERENCE_MAPPING_GLOB_TO_LOC;

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
		static constexpr int worldDim = TWorldDim;

	///	reference dimension (domain space dimension)
		static constexpr int dim = TRefElem::dim;

	///	flag if mapping is linear (i.e. Jacobian does not depend on x)
		static constexpr bool isLinear = false;

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
							 size_t maxIter = 1000,
							 number tol = 1e-10) const;

	///	map global coordinate to local coordinate for n local positions
		void global_to_local(MathVector<dim>* vLocPos,
							 const MathVector<worldDim>* vGlobPos, size_t n,
							 size_t maxIter = 1000,
							 number tol = 1e-10) const;

	///	map global coordinate to local coordinate for a vector of local positions
		void global_to_local(std::vector<MathVector<dim> >& vLocPos,
							 const std::vector<MathVector<worldDim> >& vGlobPos,
							 size_t maxIter = 1000,
							 number tol = 1e-10) const;

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
		void global_to_local
		(
			MathVector<dim>& locPos, ///< (o) for the computed local coordinates (i) specify the initial guess!
			const MathVector<worldDim>& globPos, ///< (i) the global coordinates to transform
			const size_t maxIter = 1000, ///< (i) maximum number of the Newton iterations
			const number tol = 1e-10 ///< (i) tolerance (smalles possible correction in the Newton iterations)
		) const
		{
			// We apply the Newton's method for the transformation. We assume here that
			// the Newton's method converges without the linesearch, and the Jacobian is
			// non-singular in the whole iteration process. This is true in particular for
			// the bilinear transformation of convex quadrilaterals if the initial guess
			// is inside of the element.
			MathMatrix<worldDim, dim> J;
			MathMatrix<dim, worldDim> JInv;
			MathVector<worldDim> dist, compGlobPos, minDist;
			MathVector<dim> corr;
			
			VecScale(minDist, globPos, tol);

			for (size_t i = 0; i < maxIter; ++i) {

			//	f(x) := \phi(x) - x_{glob}
				getImpl().local_to_global(compGlobPos, locPos);
				VecSubtract(dist, compGlobPos, globPos);
				
			//	Floating-point cancellation protection:
				if(VecAbsIsLess(dist, minDist))
					return;
				// REMARK: Note that the tolerance is specified not only to provide the
				// sufficient accuracy for the solution but also to protect the iteration
				// from the cancellation phenomena in the floating-point arithmetics.
				// Computation of the distance involves subtraction which is a reason for
				// the loss of precision phenomena. If a small element is located very far
				// from the coordinate origin, the accuracy of the local->global transform
				// is restricted and this cannot be overcome. After we reach this bound,
				// the iterates will oscillate and the defect will not be reduced. We use
				// globPos for the reference values.
				//      Note that this check may not be used as the only stopping criterion:
				// globPos may be 0 by specification. 

				UG_DLOG(DID_REFERENCE_MAPPING_GLOB_TO_LOC, 2,
						"reference_mapping.h: global_to_local() Newton iteration: Iter# "
						<< i << "; fabs(VecTwoNorm(dist)) = " << fabs(VecTwoNorm(dist)) <<
						"; dist = " << dist << "; locPos: " << locPos << std::endl);

			//	compute jacobian df/dx = d \phi(x) / dx =: J
				getImpl().jacobian(J, locPos);

			//	solve c -= J^{-1} f
				LeftInverse(JInv, J);
				MatVecMult(corr, JInv, dist);
				
			//	check if tol reached
				if(VecAbsIsLess(corr, tol))
					return;
				// REMARK: Note that using the Euclidean or maximum norm directly to dist
				// would need tuning of the tolerance for every particular grid: For big
				// elements with the diameter of order 1, the tolerance 1e-10 is good accuracy,
				// but for a refined grid with the elements of the size of order 1e-5 it
				// can be pour. But ||corr|| = ||J^{-1} dist|| is also a norm of dist, and
				// it is rescaled according to the element dimensions.
				//      Besides that, ||corr|| measures whether we can get any further progress
				// in the iterations.

			//	apply the correction
				VecSubtract(locPos, locPos, corr);
			}

			// compiler does not know that maxIter will never be 0
			// therefore it warns that dist may be uninit'ed;
			// as we will throw here anyway, we may as well check that
			UG_COND_THROW(!maxIter, "Without a single iteration, local-to-global "
						  "mapping can never converge.");

			UG_DLOG(DID_REFERENCE_MAPPING_GLOB_TO_LOC, 2, "Last JInv:" << std::endl);
			for(int i = 0; i < 3; ++i)
			{
				UG_DLOG(DID_REFERENCE_MAPPING_GLOB_TO_LOC, 2,
						JInv(i, 0) << "; " << JInv(i, 1) << "; " << JInv(i, 2) << std::endl);
			}

			UG_THROW("ReferenceMapping::global_to_local: Newton method did not"
					" reach a tolerance "<<tol<<" after "<<maxIter<<" steps. "
					"Global Pos: "<<globPos<<", dim: "<<dim<<", worldDim: "<<
					worldDim<<", last newton defect: "<<fabs(VecTwoNorm(dist)));
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
			UG_DLOG(DID_REFERENCE_MAPPING, 2, ">>OCT_DISC_DEBUG:: " << "reference_mapping.h: " << "jacobian_transposed_inverse(): JT " << std::endl);
			for(int i = 0; i < 3; ++i)
				UG_DLOG(DID_REFERENCE_MAPPING, 2, "	JT:" << JT(i, 0) << "; " << JT(i, 1) << "; " << JT(i, 2) << std::endl);
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
// Reference Mapping RegularVertex
///////////////////////////////////////////////////////////////////////////////

template <int TWorldDim>
class ReferenceMapping<ReferenceVertex, TWorldDim>
	: public BaseReferenceMapping<ReferenceVertex::dim, TWorldDim, true,
	  	  	  	  	  	  	  	  ReferenceMapping<ReferenceVertex, TWorldDim> >
{
	public:
	///	world dimension
		static constexpr int worldDim = TWorldDim;

	///	reference dimension
		static constexpr int dim = ReferenceVertex::dim;

	///	flag if mapping is linear (i.e. Jacobian does not depend on x)
		static constexpr bool isLinear = true;

	public:
		using base_type = BaseReferenceMapping<ReferenceVertex::dim, TWorldDim, true, ReferenceMapping >;
		using base_type::local_to_global;
		using base_type::jacobian;
		using base_type::jacobian_transposed;
		using base_type::jacobian_transposed_inverse;
		using base_type::sqrt_gram_det;

	public:
	///	Default Constructor
		ReferenceMapping() = default;

	///	Constructor setting the corners
	/// \{
		explicit ReferenceMapping(const MathVector<worldDim>* vCornerCoord) {update(vCornerCoord);}
		explicit ReferenceMapping(const std::vector<MathVector<worldDim> >& vCornerCoord) {update(vCornerCoord);}
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
		}

	///	map local coordinate to global coordinate
		void local_to_global(MathVector<worldDim>& globPos,
							 const MathVector<dim>& locPos) const
		{
			globPos = co0;
		}

	///	returns transposed of jacobian
		void jacobian_transposed(MathMatrix<dim, worldDim>& JT,
								 const MathVector<dim>& locPos) const
		{
			//for(int i = 0; i < worldDim; ++i) JT(0,i) = 1;
		}

	private:
		MathVector<worldDim> co0;
};

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
		static constexpr int worldDim = TWorldDim;

	///	reference dimension
		static constexpr int dim = ReferenceEdge::dim;

	///	flag if mapping is linear (i.e. Jacobian does not depend on x)
		static constexpr bool isLinear = true;

	public:
		using base_type = BaseReferenceMapping<ReferenceEdge::dim, TWorldDim, true, ReferenceMapping >;
		using base_type::local_to_global;
		using base_type::jacobian;
		using base_type::jacobian_transposed;
		using base_type::jacobian_transposed_inverse;
		using base_type::sqrt_gram_det;

	public:
	///	Default Constructor
		ReferenceMapping() = default;

	///	Constructor setting the corners
	/// \{
		explicit ReferenceMapping(const MathVector<worldDim>* vCornerCoord) {update(vCornerCoord);}
		explicit ReferenceMapping(const std::vector<MathVector<worldDim> >& vCornerCoord) {update(vCornerCoord);}
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
		static constexpr int worldDim = TWorldDim;

	///	reference dimension
		static constexpr int dim = ReferenceTriangle::dim;

	///	flag if mapping is linear (i.e. Jacobian does not depend on x)
		static constexpr bool isLinear = true;

	public:
		using base_type = BaseReferenceMapping<ReferenceTriangle::dim, TWorldDim, true, ReferenceMapping >;
		using base_type::local_to_global;
		using base_type::jacobian;
		using base_type::jacobian_transposed;
		using base_type::jacobian_transposed_inverse;
		using base_type::sqrt_gram_det;

	public:
	///	Default Constructor
		ReferenceMapping() = default;

	///	Constructor setting the corners
	/// \{
		explicit ReferenceMapping(const MathVector<worldDim>* vCornerCoord) {update(vCornerCoord);}
		explicit ReferenceMapping(const std::vector<MathVector<worldDim> >& vCornerCoord) {update(vCornerCoord);}
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
		static constexpr int worldDim = TWorldDim;

	///	reference dimension
		static constexpr int dim = ReferenceQuadrilateral::dim;

	///	flag if mapping is linear (i.e. Jacobian does not depend on x)
		static constexpr bool isLinear = false;

	public:
		using base_type = BaseReferenceMapping<ReferenceQuadrilateral::dim, TWorldDim, false, ReferenceMapping >;
		using base_type::local_to_global;
		using base_type::jacobian;
		using base_type::jacobian_transposed;
		using base_type::jacobian_transposed_inverse;
		using base_type::sqrt_gram_det;

	public:
	///	Default Constructor
		ReferenceMapping() = default;

	///	Constructor setting the corners
	/// \{
		explicit ReferenceMapping(const MathVector<worldDim>* vCornerCoord) {update(vCornerCoord);}
		explicit ReferenceMapping(const std::vector<MathVector<worldDim> >& vCornerCoord) {update(vCornerCoord);}
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
		static constexpr int worldDim = TWorldDim;

	///	reference dimension
		static constexpr int dim = ReferenceTetrahedron::dim;

	///	flag if mapping is linear (i.e. Jacobian does not depend on x)
		static constexpr bool isLinear = true;

	public:
		using base_type = BaseReferenceMapping<ReferenceTetrahedron::dim, TWorldDim, true,
			ReferenceMapping<ReferenceTetrahedron, TWorldDim> >;
		using base_type::local_to_global;
		using base_type::jacobian;
		using base_type::jacobian_transposed;
		using base_type::jacobian_transposed_inverse;
		using base_type::sqrt_gram_det;

	public:
	///	Default Constructor
		ReferenceMapping() = default;

	///	Constructor setting the corners
	/// \{
		explicit ReferenceMapping(const MathVector<worldDim>* vCornerCoord) {update(vCornerCoord);}
		explicit ReferenceMapping(const std::vector<MathVector<worldDim> >& vCornerCoord) {update(vCornerCoord);}
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
		static constexpr int worldDim = TWorldDim;

	///	reference dimension
		static constexpr int dim = ReferencePyramid::dim;

	///	flag if mapping is linear (i.e. Jacobian does not depend on x)
		static constexpr bool isLinear = false;

	public:
		using base_type = BaseReferenceMapping<ReferencePyramid::dim, TWorldDim, false, ReferenceMapping >;
		using base_type::local_to_global;
		using base_type::jacobian;
		using base_type::jacobian_transposed;
		using base_type::jacobian_transposed_inverse;
		using base_type::sqrt_gram_det;

	public:
	///	Default Constructor
		ReferenceMapping() = default;

	///	Constructor setting the corners
	/// \{
		explicit ReferenceMapping(const MathVector<worldDim>* vCornerCoord) {update(vCornerCoord);}
		explicit ReferenceMapping(const std::vector<MathVector<worldDim> >& vCornerCoord) {update(vCornerCoord);}
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
		static constexpr int worldDim = TWorldDim;

	///	reference dimension
		static constexpr int dim = ReferencePrism::dim;

	///	flag if mapping is linear (i.e. Jacobian does not depend on x)
		static constexpr bool isLinear = false;

	public:
		using base_type = BaseReferenceMapping<ReferencePrism::dim, TWorldDim, false, ReferenceMapping >;
		using base_type::local_to_global;
		using base_type::jacobian;
		using base_type::jacobian_transposed;
		using base_type::jacobian_transposed_inverse;
		using base_type::sqrt_gram_det;

	public:
	///	Default Constructor
		ReferenceMapping() = default;

	///	Constructor setting the corners
	/// \{
		explicit ReferenceMapping(const MathVector<worldDim>* vCornerCoord) {update(vCornerCoord);}
		explicit ReferenceMapping(const std::vector<MathVector<worldDim> >& vCornerCoord) {update(vCornerCoord);}
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
		static constexpr int worldDim = TWorldDim;

	///	reference dimension
		static constexpr int dim = ReferenceHexahedron::dim;

	///	flag if mapping is linear (i.e. Jacobian does not depend on x)
		static constexpr bool isLinear = false;

	public:
		using base_type = BaseReferenceMapping<ReferenceHexahedron::dim, TWorldDim, false, ReferenceMapping >;
		using base_type::local_to_global;
		using base_type::jacobian;
		using base_type::jacobian_transposed;
		using base_type::jacobian_transposed_inverse;
		using base_type::sqrt_gram_det;

	public:
	///	Default Constructor
		ReferenceMapping() = default;

	///	Constructor setting the corners
	/// \{
		explicit ReferenceMapping(const MathVector<worldDim>* vCornerCoord) {update(vCornerCoord);}
		explicit ReferenceMapping(const std::vector<MathVector<worldDim> >& vCornerCoord) {update(vCornerCoord);}
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
		    const number a = 1.0 - locPos[0];
			const number b = 1.0 - locPos[1];
			const number c = 1.0 - locPos[2];
			number a0 = b * c;
			number a1 = locPos[1] * c;
			number a2 = locPos[1] * locPos[2];
			number a3 = b * locPos[2];
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

///////////////////////////////////////////////////////////////////////////////
// Reference Mapping Octahedron
///////////////////////////////////////////////////////////////////////////////
template <int TWorldDim>
class ReferenceMapping<ReferenceOctahedron, TWorldDim>
	: public BaseReferenceMapping<ReferenceOctahedron::dim, TWorldDim, false,
								  ReferenceMapping<ReferenceOctahedron, TWorldDim> >
{
	public:
	///	world dimension
		static constexpr int worldDim = TWorldDim;

	///	reference dimension
		static constexpr int dim = ReferenceOctahedron::dim;

	///	flag if mapping is linear (i.e. Jacobian does not depend on x)
		static constexpr bool isLinear = false;

	public:
		using base_type = BaseReferenceMapping<ReferenceOctahedron::dim, TWorldDim, false, ReferenceMapping >;
		using base_type::local_to_global;
		using base_type::jacobian;
		using base_type::jacobian_transposed;
		using base_type::jacobian_transposed_inverse;
		using base_type::sqrt_gram_det;

	public:
	///	Default Constructor
		ReferenceMapping() = default;

	///	Constructor setting the corners
	/// \{
		explicit ReferenceMapping(const MathVector<worldDim>* vCornerCoord) {update(vCornerCoord);}
		explicit ReferenceMapping(const std::vector<MathVector<worldDim> >& vCornerCoord) {update(vCornerCoord);}
	/// \}

	///	refresh mapping for new set of corners
		void update(const std::vector<MathVector<worldDim> >& vCornerCoord)
		{
			UG_ASSERT((int)vCornerCoord.size() >= ReferenceOctahedron::numCorners,
			          "ReferenceMapping: to few Corner Coordinates.");
			update(&vCornerCoord[0]);
		}

	///	update the mapping for a new set of corners
		void update(const MathVector<worldDim>* vCornerCoord)
		{
			for(int co = 0; co < ReferenceOctahedron::numCorners; ++co)
				x[co] = vCornerCoord[co];
		}

	///	map local coordinate to global coordinate
		void local_to_global(MathVector<worldDim>& globPos,
							 const MathVector<dim>& locPos) const
		{
		//	the octahedral shape functions correspond to the tetrahedral ones,
		//	but locally piecewise linear on a subdivision of the octahedron
		//	into 4 sub-tetrahedrons.
			const number z_sgn 	= (locPos[2] < 0) ? -1.0 : 1.0;
			const number a0	 	= (locPos[2] < 0) ? -locPos[2] : 0.0;
			const number a5	 	= (locPos[2] < 0) ?  0.0 : locPos[2];

			if (locPos[0] > locPos[1])
			{
				const number a1 = 1.0 - locPos[0] - z_sgn * locPos[2];
				const number a2 = locPos[0] - locPos[1];
				const number a3 = locPos[1];
				const number a4 = 0.0;

				for(int d = 0; d < worldDim; ++d)
					globPos[d] = a0*x[0][d]+a1*x[1][d]+a2*x[2][d]
								+a3*x[3][d]+a4*x[4][d]+a5*x[5][d];
			}
			else
			{
				const number a1 = 1.0 - locPos[1] - z_sgn * locPos[2];
				const number a2 = 0.0;
				const number a3 = locPos[0];
				const number a4 = locPos[1] - locPos[0];

				for(int d = 0; d < worldDim; ++d)
					globPos[d] = a0*x[0][d]+a1*x[1][d]+a2*x[2][d]
								+a3*x[3][d]+a4*x[4][d]+a5*x[5][d];
			}
		}

	///	returns transposed of jacobian
		void jacobian_transposed(MathMatrix<dim, worldDim>& JT,
								 const MathVector<dim>& locPos) const
	   {
		//	the octahedral shape functions correspond to the tetrahedral ones,
		//	but locally piecewise linear on a subdivision of the octahedron
		//	into 4 sub-tetrahedrons.
			const number z_sgn 	= (locPos[2] < 0) ? -1.0 : 1.0;
			const number Da0	= (locPos[2] < 0) ? -1.0 : 0.0;
			const number Da5	= (locPos[2] < 0) ?  0.0 : 1.0;

			if (locPos[0] > locPos[1])
			{
				for(int d = 0; d < worldDim; ++d)
					JT(0,d) = -x[1][d]+x[2][d];

				for(int d = 0; d < worldDim; ++d)
					JT(1,d) = -x[2][d]+x[3][d];

				for(int d = 0; d < worldDim; ++d)
					JT(2,d) = Da0*x[0][d] - z_sgn*x[1][d] + Da5*x[5][d];
			}
			else
			{
				for(int d = 0; d < worldDim; ++d)
					JT(0,d) = x[3][d]-x[4][d];

				for(int d = 0; d < worldDim; ++d)
					JT(1,d) = -x[1][d]+x[4][d];

				for(int d = 0; d < worldDim; ++d)
					JT(2,d) = Da0*x[0][d] - z_sgn*x[1][d] + Da5*x[5][d];
			}
		}

	private:
		MathVector<worldDim> x[ReferenceOctahedron::numCorners];
};


} // end namespace ug

#endif