/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Raphael Prohl
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

#ifndef __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__OBSTACLE_IN_NORMAL_DIR_IMPL__
#define __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__OBSTACLE_IN_NORMAL_DIR_IMPL__

#include "obstacle_in_normal_dir.h"

#include "lib_disc/domain_util.h"
#include "lib_disc/common/geometry_util.h"
#include "lib_disc/spatial_disc/disc_util/fe_geom.h"

namespace ug{

template <typename TDomain, typename TAlgebra>
void
ObstacleInNormalDir<TDomain,TAlgebra>::preprocess()
{
	//	for debugging
	MathVector<dim> normal;
	normal[0] = 2.0;
	if (dim > 1) normal[1] = 1.5;
	if (dim > 2) normal[2] = -3.0;

	MathVector<dim> transformedONB[dim];
	transform_eulerian_coord_sys(transformedONB, normal);

}

template <typename TDomain, typename TAlgebra>
void
ObstacleInNormalDir<TDomain,TAlgebra>::transform_eulerian_coord_sys(MathVector<dim> transformedONB[],
		const MathVector<dim>& firstTransformedBaseVec)
{
	//	create first unit eulerian base vector
	MathVector<dim> unityX;
	unityX[0] = 1.0;
	if (dim > 1) unityX[1] = 0.0;
	if (dim > 2) unityX[2] = 0.0;

	//	normalize first base vector of the transformed system
	MathVector<dim> normalizedFirstBaseVec;
	const number normOfFirstBaseVec = VecLength(firstTransformedBaseVec);
	VecScale(normalizedFirstBaseVec, firstTransformedBaseVec, 1.0/normOfFirstBaseVec);

	//	compute vector, which is orthogonal to the householder hypersphere
	MathVector<dim> orthoVec;
	VecScaleAdd(orthoVec, 0.5, unityX, -0.5, normalizedFirstBaseVec);

	//	compute householder matrix
	MathMatrix<dim,dim> hMat;
	MatHouseholder(hMat, orthoVec);

	//	get and store transformed orthonormal base vectors
	MathVector<dim> transBaseVec;
	for(size_t i = 0; i < (size_t)dim; ++i)
	{
		for(size_t j = 0; j < (size_t)dim; ++j){
			transBaseVec[j] = hMat(i,j);
		}
		transformedONB[i] = transBaseVec;
	}
}

template <int dim> struct face_type_traits
{
	using face_type0 = void;
    using face_type1 = void;
    using DimFEGeo = void;
};

template <> struct face_type_traits<1>
{
	using face_type0 = ReferenceVertex;
    using face_type1 = ReferenceVertex;
    using DimFEGeo = DimFEGeometry<1, 1>;
};

template <> struct face_type_traits<2>
{
	using face_type0 = ReferenceEdge;
    using face_type1 = ReferenceEdge;
    using DimFEGeo = DimFEGeometry<2, 1>;
};

template <> struct face_type_traits<3>
{
	using face_type0 = ReferenceTriangle;
    using face_type1 = ReferenceQuadrilateral;
    using DimFEGeo = DimFEGeometry<3, 2>;
};

template <typename TDomain, typename TAlgebra>
template <typename TElem, typename TIterator>
void
ObstacleInNormalDir<TDomain,TAlgebra>::
adjust_sol_and_cor_elem(TIterator iterBegin,
		TIterator iterEnd, value_type& sol_i, value_type& c_i, bool& dofIsActive,
		const DoFIndex& dof)
{
	using face_type0 = typename face_type_traits<dim>::face_type0;
	using face_type1 = typename face_type_traits<dim>::face_type1;

	//	storage for corner coordinates
	vector<MathVector<dim> > vCorner;
	vector<MathVector<dim> > vCoPos;

	//	outer normal vector
	MathVector<dim> normal;

	for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		TElem* elem = *iter;

	//	reference object type
		ReferenceObjectID roid = elem->reference_object_id();

		const DimReferenceElement<dim-1>& rRefElem
				= ReferenceElementProvider::get<dim-1>(roid);

	//	get corners of element
		CollectCornerCoordinates(vCorner, *elem, *m_spDomain);

	//	here the ordering of the corners in the reference element is exploited
	//	in order to compute the outer normal later on
		int nCorner = (int)vCorner.size();
		for (int i = 0; i < nCorner; ++i)
			vCoPos.push_back(vCorner[rRefElem.id(dim-1, 0, 0, i)]);

		if ((int)vCorner.size() == dim)
			ElementNormal<face_type0, dim>(normal, &vCoPos[0]);
		else
			ElementNormal<face_type1, dim>(normal, &vCoPos[0]);

		/*for (int i = 0; i < (int)vCorner.size(); ++i)
			UG_LOG("coPos: " << coPos[i] << "\n");
		UG_LOG("normal: " << normal << "\n");

		//const number normOfNormal = VecLength(normal);
		//const value_type tmpVSol = sol_i + c_i;
		//const number uTimesNormal = tmpSol / normOfNormal;
		UG_LOG("\n");*/
	}

	const size_t comp = dof[1];

	//	tmpSol := u_{s-1/2} = u_{s-1} + c
	const number tmpSol = BlockRef(sol_i, comp) + BlockRef(c_i, comp);

	//	get obstacle value corresponding to the dof
	const number obsVal = m_mObstacleValues[dof];

	//	check, if dof is active (tmpSol >= obsVal)
	//	TODO: check if u * n > g, i.e. tmpSol * n > g!
	if (!(tmpSol < obsVal))
	{
		//	is active DoF
		m_vActiveDofs.push_back(dof);

		//	adjust correction & set solution to obstacle-value
		BlockRef(c_i, comp) = obsVal - BlockRef(sol_i, comp);
		BlockRef(sol_i, comp) = obsVal;
		dofIsActive = true;
	}
}

template <typename TDomain, typename TAlgebra>
void
ObstacleInNormalDir<TDomain,TAlgebra>::
adjust_sol_and_cor(value_type& sol_i, value_type& c_i, bool& dofIsActive,
		const DoFIndex& dof)
{
	//	loop over all obstacle subsets
	for (vector<int>::iterator obsSI = m_vObsSubsets.begin();
			obsSI != m_vObsSubsets.end(); ++obsSI)
	{
		const int si = *obsSI;
		UG_LOG("si: " <<si<<"\n");
		const int subsetDim = DimensionOfSubset(*m_spDD->subset_handler(), si);

		switch(subsetDim)
		{
		case 0:
			break;
		case 1:
			adjust_sol_and_cor_elem<RegularEdge>
				(m_spDD->template begin<RegularEdge>(si), m_spDD->template end<RegularEdge>(si),
						sol_i, c_i, dofIsActive, dof);
			break;
		case 2:
			adjust_sol_and_cor_elem<Triangle>
				(m_spDD->template begin<Triangle>(si), m_spDD->template end<Triangle>(si),
						sol_i, c_i, dofIsActive, dof);
			adjust_sol_and_cor_elem<Quadrilateral>
				(m_spDD->template begin<Quadrilateral>(si), m_spDD->template end<Quadrilateral>(si),
						sol_i, c_i, dofIsActive, dof);
			break;
		default:
			UG_THROW("ObstacleInNormalDir::adjust_sol_and_cor:"
				"SubsetDimension "<< subsetDim <<" (subset="<< si <<") not supported.");
		}
	}
}

template <typename TDomain, typename TAlgebra>
void
ObstacleInNormalDir<TDomain,TAlgebra>::
adjust_defect_to_constraint(vector_type& d)
{
	for (vector<MultiIndex<2> >::iterator itActiveInd = m_vActiveDofs.begin();
			itActiveInd < m_vActiveDofs.end(); ++itActiveInd)
	{
		//	check, if Ax >= b. For that case the new defect is set to zero,
		//	since all equations/constraints are fulfilled
		number defect = BlockRef(d[(*itActiveInd)[0]], (*itActiveInd)[1]);
		if (defect > 0.0)
			BlockRef(d[(*itActiveInd)[0]], (*itActiveInd)[1]) = 0.0;
	}
}

template <typename TDomain, typename TAlgebra>
void
ObstacleInNormalDir<TDomain,TAlgebra>::
restrict_obs_values()
{}

} // end namespace ug

#endif