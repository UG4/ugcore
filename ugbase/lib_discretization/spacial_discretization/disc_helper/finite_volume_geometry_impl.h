/*
 * finite_volume_geometry_impl.h
 *
 *  Created on: 04.09.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DISC_HELPER__FINITE_VOLUME_GEOMETRY_IMPL__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DISC_HELPER__FINITE_VOLUME_GEOMETRY_IMPL__

#include "./finite_volume_geometry.h"

namespace ug{


template <typename TElem, int TWorldDim>
FV1Geometry<TElem, TWorldDim>::
FV1Geometry() : m_pElem(NULL)
{
	// set corners of element as local centers of nodes
	for(size_t i = 0; i < m_rRefElem.num_obj(0); ++i)
		m_locMid[0][i] = m_rRefElem.corner(i);

	// compute local midpoints for all geometric objects with  0 < d <= dim
	for(int d = 1; d <= dim; ++d)
	{
		// loop geometric objects of dimension d
		for(size_t i = 0; i < m_rRefElem.num_obj(d); ++i)
		{
			// set first node
			const size_t coID0 = m_rRefElem.id(d, i, 0, 0);
			m_locMid[d][i] = m_locMid[0][coID0];

			// add corner coordinates of the corners of the geometric object
			for(size_t j = 1; j < m_rRefElem.num_obj_of_obj(d, i, 0); ++j)
			{
				const size_t coID = m_rRefElem.id(d, i, 0, j);
				m_locMid[d][i] += m_locMid[0][coID];
			}

			// scale for correct averaging
			m_locMid[d][i] *= 1./(m_rRefElem.num_obj_of_obj(d, i, 0));
		}
	}

	// set up local informations for SubControlVolumeFaces (scvf)
	// each scvf is associated to one edge of the element
	for(size_t i = 0; i < num_scvf(); ++i)
	{
		m_vSCVF[i].m_from = m_rRefElem.id(1, i, 0, 0);
		m_vSCVF[i].m_to = m_rRefElem.id(1, i, 0, 1);

		// set mid ids
		{
			// start at edge midpoint
			m_vSCVF[i].midId[0] = MidID(1,i);

			// loop up dimension
			if(dim == 2)
			{
				m_vSCVF[i].midId[1] = MidID(dim, 0); // center of element
			}
			else if (dim == 3)
			{
				m_vSCVF[i].midId[1] = MidID(2, m_rRefElem.id(1, i, 2, 0)); // side 0
				m_vSCVF[i].midId[2] = MidID(dim, 0); // center of element
				m_vSCVF[i].midId[3] = MidID(2, m_rRefElem.id(1, i, 2, 1)); // side 1
			}
		}

		// copy local corners of scvf
		copy_local_corners(m_vSCVF[i]);

		// integration point
		AveragePositions(m_vSCVF[i].localIP, m_vSCVF[i].m_vLocPos, SCVF::m_numCorners);
	}

	// set up local informations for SubControlVolumes (scv)
	// each scv is associated to one corner of the element
	for(size_t i = 0; i < num_scv(); ++i)
	{
		m_vSCV[i].nodeId = i;

		if(dim == 1)
		{
			m_vSCV[i].midId[0] = MidID(0, i); // set node as corner of scv
			m_vSCV[i].midId[1] = MidID(dim, 0);	// center of element
		}
		else if(dim == 2)
		{
			m_vSCV[i].midId[0] = MidID(0, i); // set node as corner of scv
			m_vSCV[i].midId[1] = MidID(1, m_rRefElem.id(0, i, 1, 0)); // edge 1
			m_vSCV[i].midId[2] = MidID(dim, 0);	// center of element
			m_vSCV[i].midId[3] = MidID(1, m_rRefElem.id(0, i, 1, 1)); // edge 2
		}
		else if(dim == 3 && (ref_elem_type::REFERENCE_OBJECT_ID != ROID_PYRAMID || i != num_scv()-1))
		{
			m_vSCV[i].midId[0] = MidID(0, i); // set node as corner of scv
			m_vSCV[i].midId[1] = MidID(1, m_rRefElem.id(0, i, 1, 1)); // edge 1
			m_vSCV[i].midId[2] = MidID(2, m_rRefElem.id(0, i, 2, 0)); // face 0
			m_vSCV[i].midId[3] = MidID(1, m_rRefElem.id(0, i, 1, 0)); // edge 0
			m_vSCV[i].midId[4] = MidID(1, m_rRefElem.id(0, i, 1, 2)); // edge 2
			m_vSCV[i].midId[5] = MidID(2, m_rRefElem.id(0, i, 2, 2)); // face 2
			m_vSCV[i].midId[6] = MidID(dim, 0);	// center of element
			m_vSCV[i].midId[7] = MidID(2, m_rRefElem.id(0, i, 2, 1)); // face 1
		}
		// TODO: Implement last ControlVolume for Pyramid
		else if(dim == 3 && ref_elem_type::REFERENCE_OBJECT_ID == ROID_PYRAMID && i == num_scv()-1)
		{
			// this scv has 10 corners
			m_vSCV[i].m_numCorners = 10;
			//UG_ASSERT(0, "Last SCV for Pyramid must be implemented");
		}
		else {UG_ASSERT(0, "Dimension higher that 3 not implemented.");}

		// copy local corners of scv
		copy_local_corners(m_vSCV[i]);
	}

	/////////////////////////
	// Shapes and Derivatives
	/////////////////////////
	for(size_t i = 0; i < num_scvf(); ++i)
	{
		const LocalShapeFunctionSet<ref_elem_type>& TrialSpace =
				LocalShapeFunctionSetFactory::inst().get_local_shape_function_set<ref_elem_type>(LSFS_LAGRANGEP1);

		const size_t num_sh = m_numSCV;
		m_vSCVF[i].vShape.resize(num_sh);
		m_vSCVF[i].localGrad.resize(num_sh);
		m_vSCVF[i].globalGrad.resize(num_sh);
		for(size_t sh = 0 ; sh < num_sh; ++sh)
		{
			if(!TrialSpace.evaluate(sh, m_vSCVF[i].localIP, (m_vSCVF[i].vShape)[sh]))
				{UG_LOG("Cannot evaluate local shape.\n"); UG_ASSERT(0, "Error in Constructor.");}
			if(!TrialSpace.evaluate_grad(sh, m_vSCVF[i].localIP, (m_vSCVF[i].localGrad)[sh]))
				{UG_LOG("Cannot evaluate local grad.\n"); UG_ASSERT(0, "Error in Constructor.");}
		}
	}
}


/// update data for given element
template <typename TElem, int TWorldDim>
bool
FV1Geometry<TElem, TWorldDim>::
update(TElem* elem, const Grid& grid, const MathVector<world_dim>* vCornerCoords)
{
	// If already update for this element, do nothing
	if(m_pElem == elem) return true;
	else m_pElem = elem;

	// remember global position of nodes
	for(size_t i = 0; i < m_rRefElem.num_obj(0); ++i)
		m_gloMid[0][i] = vCornerCoords[i];

	// compute global midpoints for all geometric objects with  0 < d <= dim
	for(int d = 1; d <= dim; ++d)
	{
		// loop geometric objects of dimension d
		for(size_t i = 0; i < m_rRefElem.num_obj(d); ++i)
		{
			// set first node
			const size_t coID0 = m_rRefElem.id(d, i, 0, 0);
			m_gloMid[d][i] = m_gloMid[0][coID0];

			// add corner coordinates of the corners of the geometric object
			for(size_t j = 1; j < m_rRefElem.num_obj_of_obj(d, i, 0); ++j)
			{
				const size_t coID = m_rRefElem.id(d, i, 0, j);
				m_gloMid[d][i] += m_gloMid[0][coID];
			}

			// scale for correct averaging
			m_gloMid[d][i] *= 1./(m_rRefElem.num_obj_of_obj(d, i, 0));
		}
	}

	// compute global informations for scvf
	for(size_t i = 0; i < num_scvf(); ++i)
	{
		// copy local corners of scvf
		copy_global_corners(m_vSCVF[i]);

		// integration point
		AveragePositions(m_vSCVF[i].globalIP, m_vSCVF[i].m_vGloPos, SCVF::m_numCorners);

		// normal on scvf
		NormalOnSCVF<ref_elem_type, world_dim>(m_vSCVF[i].Normal, m_vSCVF[i].m_vGloPos);
	}

	// compute size of scv
	for(size_t i = 0; i < num_scv(); ++i)
	{
		// copy global corners
		copy_global_corners(m_vSCV[i]);

		// compute volume of scv
		if(m_vSCV[i].m_numCorners != 10)
		{
			m_vSCV[i].vol = ElementSize<scv_type, world_dim>(m_vSCV[i].m_vGloPos);
		}
		else
		{
			// special case for pyramid, last scv
		}
	}

	/////////////////////////
	// Shapes and Derivatives
	/////////////////////////
	m_rMapping.update(vCornerCoords);

	for(size_t i = 0; i < num_scvf(); ++i)
	{
		if(!m_rMapping.jacobian_transposed_inverse(m_vSCVF[i].localIP, m_vSCVF[i].JtInv))
			{UG_LOG("Cannot compute jacobian transposed.\n"); return false;}
		if(!m_rMapping.jacobian_det(m_vSCVF[i].localIP, m_vSCVF[i].detj))
			{UG_LOG("Cannot compute jacobian determinate.\n"); return false;}

		for(size_t sh = 0 ; sh < num_scv(); ++sh)
			MatVecMult((m_vSCVF[i].globalGrad)[sh], m_vSCVF[i].JtInv, (m_vSCVF[i].localGrad)[sh]);
	}

	//print();
	return true;
}





} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DISC_HELPER__FINITE_VOLUME_GEOMETRY_IMPL__ */
