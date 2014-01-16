/*
 *	obstacle_in_normal_dir_impl.h
 *
 *  Created on: 26.11.2013
 *      Author: raphaelprohl
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
    typedef void face_type0;
	typedef void face_type1;
	typedef void DimFEGeo;
};

template <> struct face_type_traits<1>
{
    typedef ReferenceVertex face_type0;
	typedef ReferenceVertex face_type1;
	typedef DimFEGeometry<1, 1> DimFEGeo;
};

template <> struct face_type_traits<2>
{
    typedef ReferenceEdge face_type0;
	typedef ReferenceEdge face_type1;
	typedef DimFEGeometry<2, 1> DimFEGeo;
};

template <> struct face_type_traits<3>
{
    typedef ReferenceTriangle face_type0;
	typedef ReferenceQuadrilateral face_type1;
	typedef DimFEGeometry<3, 2> DimFEGeo;
};

template <typename TDomain, typename TAlgebra>
template <typename TElem, typename TIterator>
void
ObstacleInNormalDir<TDomain,TAlgebra>::
adjust_sol_and_cor_elem(TIterator iterBegin,
		TIterator iterEnd, value_type& sol_i, value_type& c_i, bool& dofIsActive,
		const DoFIndex& dof)
{
	typedef typename face_type_traits<dim>::face_type0 face_type0;
	typedef typename face_type_traits<dim>::face_type1 face_type1;

	//	storage for corner coordinates
	vector<MathVector<dim> > vCorner;

	//	outer normal vector
	MathVector<dim> coPos[dim+1], normal;

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
		for (int i = 0; i < (int)vCorner.size(); ++i)
			coPos[i] = vCorner[rRefElem.id(dim-1, 0, 0, i)];

		if ((int)vCorner.size() == dim)
		{
			ElementNormal<face_type0, dim>(normal, coPos);
			//UG_LOG("face_type0 \n");
		}
		else
		{
			ElementNormal<face_type1, dim>(normal, coPos);
			//UG_LOG("face_type1 \n");
		}

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
			adjust_sol_and_cor_elem<Edge>
				(m_spDD->template begin<Edge>(si), m_spDD->template end<Edge>(si),
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
adjust_defect(vector_type& d)
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

#endif /* __H__UG__LIB_ALGEBRA__OPERATOR__PRECONDITIONER__PROJECTED_GAUSS_SEIDEL__OBSTACLE_IN_NORMAL_DIR_IMPL__ */

