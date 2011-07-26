/*
 * finite_volume_geometry.cpp
 *
 *  Created on: 04.09.2010
 *      Author: andreasvogel
 */


#include "common/util/provider.h"
#include "finite_volume_geometry.h"
#include "lib_discretization/reference_element/reference_element.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// FV1 Geometry for Reference Element Type
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <typename TElem, int TWorldDim>
FV1Geometry<TElem, TWorldDim>::
FV1Geometry()
	: m_pElem(NULL), m_rRefElem(Provider::get<ref_elem_type>())
{
// 	set corners of element as local centers of nodes
	for(size_t i = 0; i < m_rRefElem.num(0); ++i)
		m_locMid[0][i] = m_rRefElem.corner(i);

// 	compute local midpoints for all geometric objects with  0 < d <= dim
	for(int d = 1; d <= dim; ++d)
	{
	// 	loop geometric objects of dimension d
		for(size_t i = 0; i < m_rRefElem.num(d); ++i)
		{
		// 	set first node
			const size_t coID0 = m_rRefElem.id(d, i, 0, 0);
			m_locMid[d][i] = m_locMid[0][coID0];

		// 	add corner coordinates of the corners of the geometric object
			for(size_t j = 1; j < m_rRefElem.num(d, i, 0); ++j)
			{
				const size_t coID = m_rRefElem.id(d, i, 0, j);
				m_locMid[d][i] += m_locMid[0][coID];
			}

		// 	scale for correct averaging
			m_locMid[d][i] *= 1./(m_rRefElem.num(d, i, 0));
		}
	}

// 	set up local informations for SubControlVolumeFaces (scvf)
// 	each scvf is associated to one edge of the element
	for(size_t i = 0; i < num_scvf(); ++i)
	{
		m_vSCVF[i].m_from = m_rRefElem.id(1, i, 0, 0);
		m_vSCVF[i].m_to = m_rRefElem.id(1, i, 0, 1);

	//	set mid ids
		{
		// 	start at edge midpoint
			m_vSCVF[i].midId[0] = MidID(1,i);

		// 	loop up dimension
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

	// 	copy local corners of scvf
		copy_local_corners(m_vSCVF[i]);

	// 	integration point
		if(dim != 1)
			AveragePositions(m_vSCVF[i].localIP, m_vSCVF[i].m_vLocPos, SCVF::numCorners);
		else
			m_vSCVF[i].localIP = m_locMid[1][0];
	}

// 	set up local informations for SubControlVolumes (scv)
// 	each scv is associated to one corner of the element
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
			throw(UGFatalError("Last SCV for Pyramid must be implemented"));
		}
		else {throw(UGFatalError("Dimension higher that 3 not implemented."));}

	// 	copy local corners of scv
		copy_local_corners(m_vSCV[i]);
	}

	/////////////////////////
	// Shapes and Derivatives
	/////////////////////////

	for(size_t i = 0; i < num_scvf(); ++i)
	{
		const LocalShapeFunctionSet<ref_elem_type>& TrialSpace =
				LocalShapeFunctionSetProvider::
					get<ref_elem_type>
					(LFEID(LFEID::LAGRANGE, 1));

		TrialSpace.shapes(&(m_vSCVF[i].vShape[0]), m_vSCVF[i].localIP);
		TrialSpace.grads(&(m_vSCVF[i].localGrad[0]), m_vSCVF[i].localIP);

	}

	///////////////////////////
	// Copy ip pos in list
	///////////////////////////

// 	loop Sub Control Volumes Faces (SCVF)
	for(size_t i = 0; i < num_scvf(); ++i)
	{
	//	get current SCVF
		const SCVF& rSCVF = scvf(i);

	// 	loop ips
		m_vLocSCVF_IP[i] = rSCVF.local_ip();
	}
}


/// update data for given element
template <typename TElem, int TWorldDim>
bool
FV1Geometry<TElem, TWorldDim>::
update(TElem* elem, const ISubsetHandler& ish, const MathVector<worldDim>* vCornerCoords)
{
// 	If already update for this element, do nothing
	if(m_pElem == elem) return true;
	else m_pElem = elem;

// 	remember global position of nodes
	for(size_t i = 0; i < m_rRefElem.num(0); ++i)
		m_gloMid[0][i] = vCornerCoords[i];

// 	compute global midpoints for all geometric objects with  0 < d <= dim
	for(int d = 1; d <= dim; ++d)
	{
	// 	loop geometric objects of dimension d
		for(size_t i = 0; i < m_rRefElem.num(d); ++i)
		{
		// 	set first node
			const size_t coID0 = m_rRefElem.id(d, i, 0, 0);
			m_gloMid[d][i] = m_gloMid[0][coID0];

		// 	add corner coordinates of the corners of the geometric object
			for(size_t j = 1; j < m_rRefElem.num(d, i, 0); ++j)
			{
				const size_t coID = m_rRefElem.id(d, i, 0, j);
				m_gloMid[d][i] += m_gloMid[0][coID];
			}

		// 	scale for correct averaging
			m_gloMid[d][i] *= 1./(m_rRefElem.num(d, i, 0));
		}
	}

// 	compute global informations for scvf
	for(size_t i = 0; i < num_scvf(); ++i)
	{
	// 	copy local corners of scvf
		copy_global_corners(m_vSCVF[i]);

	// 	integration point
		if(dim != 1)
		AveragePositions(m_vSCVF[i].globalIP, m_vSCVF[i].m_vGloPos, SCVF::numCorners);
		else
			m_vSCVF[i].globalIP = m_gloMid[1][0];

	// 	normal on scvf
		if(ref_elem_type::dim == worldDim)
		{
			NormalOnSCVF<ref_elem_type, worldDim>(m_vSCVF[i].Normal, m_vSCVF[i].m_vGloPos);
		}
		else
		{
			if(ref_elem_type::dim == 1)
				NormalOnSCVF<ref_elem_type, worldDim>(m_vSCVF[i].Normal, vCornerCoords);
			else
				throw(UGFatalError("Not Implemented"));
		}
	}

// 	compute size of scv
	for(size_t i = 0; i < num_scv(); ++i)
	{
	// 	copy global corners
		copy_global_corners(m_vSCV[i]);

	// 	compute volume of scv
		if(m_vSCV[i].m_numCorners != 10)
		{
			m_vSCV[i].vol = ElementSize<scv_type, worldDim>(m_vSCV[i].m_vGloPos);
		}
		else
		{
			// special case for pyramid, last scv
			throw(UGFatalError("Pyramid Not Implemented"));
		}
	}

	/////////////////////////
	// Shapes and Derivatives
	/////////////////////////


	m_rMapping.update(vCornerCoords);

//	if mapping is linear, compute jacobian only once and copy
	if(ReferenceMapping<ref_elem_type, worldDim>::isLinear)
	{
		m_rMapping.jacobian_transposed_inverse(m_vSCVF[0].JtInv, m_vSCVF[0].localIP);
		m_vSCVF[0].detj = m_rMapping.jacobian_det(m_vSCVF[0].localIP);
		for(size_t i = 1; i < num_scvf(); ++i)
		{
			m_vSCVF[i].JtInv = m_vSCVF[0].JtInv;
			m_vSCVF[i].detj = m_vSCVF[0].detj;
		}

		for(size_t i = 0; i < num_scvf(); ++i)
			for(size_t sh = 0 ; sh < num_scv(); ++sh)
				MatVecMult((m_vSCVF[i].globalGrad)[sh], m_vSCVF[0].JtInv, (m_vSCVF[i].localGrad)[sh]);

	}
//	else compute jacobian for each integration point
	else
	{
		for(size_t i = 0; i < num_scvf(); ++i)
		{
			m_rMapping.jacobian_transposed_inverse(m_vSCVF[i].JtInv, m_vSCVF[i].localIP);
			m_vSCVF[i].detj = m_rMapping.jacobian_det(m_vSCVF[i].localIP);

			for(size_t sh = 0 ; sh < num_scv(); ++sh)
				MatVecMult((m_vSCVF[i].globalGrad)[sh], m_vSCVF[i].JtInv, (m_vSCVF[i].localGrad)[sh]);
		}
	}

	///////////////////////////
	// Copy ip pos in list
	///////////////////////////

// 	loop Sub Control Volumes Faces (SCVF)
	for(size_t i = 0; i < numSCVF; ++i)
	{
	//	get current SCVF
		const SCVF& rSCVF = scvf(i);

		m_vGlobSCVF_IP[i] = rSCVF.global_ip();
	}

	/////////////////////////
	// Boundary Faces
	/////////////////////////

//	if no boundary subsets required, return
	if(num_boundary_subsets() == 0) return true;

//	get grid
	Grid& grid = *ish.get_assigned_grid();

//	vector of subset indices of side
	std::vector<int> vSubsetIndex;

//	get subset indices for sides (i.e. edge in 2d, faces in 3d)
	if(dim == 2)
	{
		std::vector<EdgeBase*> vEdges;
		CollectEdgesSorted(vEdges, grid, elem);

		vSubsetIndex.resize(vEdges.size());
		for(size_t i = 0; i < vEdges.size(); ++i)
		{
			vSubsetIndex[i] = ish.get_subset_index(vEdges[i]);
		}
	}

	if(dim == 3)
	{
		std::vector<Face*> vFaces;
		CollectFacesSorted(vFaces, grid, elem);

		vSubsetIndex.resize(vFaces.size());
		for(size_t i = 0; i < vFaces.size(); ++i)
		{
			vSubsetIndex[i] = ish.get_subset_index(vFaces[i]);
		}
	}

//	loop requested subset
	typename std::map<int, std::vector<BF> >::iterator it;
	for (it=m_mapVectorBF.begin() ; it != m_mapVectorBF.end(); ++it)
	{
	//	get subset index
		const int bndIndex = (*it).first;

	//	get vector of BF for element
		std::vector<BF>& vBF = (*it).second;

	//	clear vector
		vBF.clear();

	//	current number of bf
		size_t curr_bf = 0;

	//	loop sides of element
		for(size_t side = 0; side < vSubsetIndex.size(); ++side)
		{
		//	skip non boundary sides
			if(vSubsetIndex[side] != bndIndex) continue;

		///////////////////////////
		//	add Boundary faces
		///////////////////////////

		//	number of corners of side
			const int coOfSide = m_rRefElem.num(dim-1, side, 0);

		//	resize vector
			vBF.resize(curr_bf + coOfSide);

		//	loop corners
			for(int co = 0; co < coOfSide; ++co)
			{
			//	get current bf
				BF& bf = vBF[curr_bf];

			//	set node id == scv this bf belongs to
				bf.m_nodeId = m_rRefElem.id(dim-1, side, 0, co);

			// 	set mid ids
				if(dim == 2)
				{
					bf.midId[co%2] = MidID(0, m_rRefElem.id(1, side, 0, co)); // corner of side
					bf.midId[(co+1)%2] = MidID(1, side); // side midpoint
				}
				else if (dim == 3)
				{
					bf.midId[0] = MidID(0, m_rRefElem.id(2, side, 0, co)); // corner of side
					bf.midId[1] = MidID(1, m_rRefElem.id(2, side, 1, co)); // edge co
					bf.midId[2] = MidID(2, side); // side midpoint
					bf.midId[3] = MidID(1, m_rRefElem.id(2, side, 1, (co -1 + coOfSide)%coOfSide)); // edge co-1
				}

			// 	copy corners of bf
				copy_local_corners(bf);
				copy_global_corners(bf);

			// 	integration point
				AveragePositions(bf.localIP, bf.m_vLocPos, SCVF::numCorners);
				AveragePositions(bf.globalIP, bf.m_vGloPos, SCVF::numCorners);

			// 	normal on scvf
				NormalOnSCVF<ref_elem_type, worldDim>(bf.Normal, bf.m_vGloPos);

			//	compute volume
				bf.m_volume = VecTwoNorm(bf.Normal);

			/////////////////////////
			// Shapes and Derivatives
			/////////////////////////
				const LocalShapeFunctionSet<ref_elem_type>& TrialSpace =
						LocalShapeFunctionSetProvider::
							get<ref_elem_type>
								(LFEID(LFEID::LAGRANGE, 1));

				TrialSpace.shapes(&(bf.vShape[0]), bf.localIP);
				TrialSpace.grads(&(bf.localGrad[0]), bf.localIP);

				m_rMapping.jacobian_transposed_inverse(bf.JtInv, bf.localIP);
				bf.detj = m_rMapping.jacobian_det(bf.localIP);

				for(size_t sh = 0 ; sh < num_scv(); ++sh)
					MatVecMult((bf.globalGrad)[sh], bf.JtInv, (bf.localGrad)[sh]);

			//	increase curr_bf
				++curr_bf;
			}
		}
	}

	return true;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Dim-dependent Finite Volume Geometry
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/// update data for given element
template <int TDim, int TWorldDim>
bool
DimFV1Geometry<TDim, TWorldDim>::
update(GeometricObject* pElem, const ISubsetHandler& ish, const MathVector<worldDim>* vCornerCoords)
{
// 	If already update for this element, do nothing
	if(m_pElem == pElem) return true;
	else m_pElem = pElem;

	//////////////////////////////////
	// COMPUTATION OF LOCAL DATA
	//////////////////////////////////

//	refresh local data, if different roid given
	if(m_roid != pElem->reference_object_id())
	{
	//	remember new roid
		m_roid = (ReferenceObjectID) pElem->reference_object_id();

	//	get reference element
		try{
		const DimReferenceElement<dim>& rRefElem
			= ReferenceElementProvider::get<dim>(m_roid);

	// 	set corners of element as local centers of nodes
		for(size_t i = 0; i < rRefElem.num(0); ++i)
			m_locMid[0][i] = rRefElem.corner(i);

	// 	compute local midpoints for all geometric objects with  0 < d <= dim
		for(int d = 1; d <= dim; ++d)
		{
		// 	loop geometric objects of dimension d
			for(size_t i = 0; i < rRefElem.num(d); ++i)
			{
			// 	set first node
				const size_t coID0 = rRefElem.id(d, i, 0, 0);
				m_locMid[d][i] = m_locMid[0][coID0];

			// 	add corner coordinates of the corners of the geometric object
				for(size_t j = 1; j < rRefElem.num(d, i, 0); ++j)
				{
					const size_t coID = rRefElem.id(d, i, 0, j);
					m_locMid[d][i] += m_locMid[0][coID];
				}

			// 	scale for correct averaging
				m_locMid[d][i] *= 1./(rRefElem.num(d, i, 0));
			}
		}

	//	set number of scvf / scv of this roid
		m_numSCV = rRefElem.num(0); // number of corners
		m_numSCVF = rRefElem.num(1); // number of edges

	// 	set up local informations for SubControlVolumeFaces (scvf)
	// 	each scvf is associated to one edge of the element
		for(size_t i = 0; i < num_scvf(); ++i)
		{
			m_vSCVF[i].m_from = rRefElem.id(1, i, 0, 0);
			m_vSCVF[i].m_to = rRefElem.id(1, i, 0, 1);

		//	set mid ids
			{
			// 	start at edge midpoint
				m_vSCVF[i].midId[0] = MidID(1,i);

			// 	loop up dimension
				if(dim == 2)
				{
					m_vSCVF[i].midId[1] = MidID(dim, 0); // center of element
				}
				else if (dim == 3)
				{
					m_vSCVF[i].midId[1] = MidID(2, rRefElem.id(1, i, 2, 0)); // side 0
					m_vSCVF[i].midId[2] = MidID(dim, 0); // center of element
					m_vSCVF[i].midId[3] = MidID(2, rRefElem.id(1, i, 2, 1)); // side 1
				}
			}

		// 	copy local corners of scvf
			copy_local_corners(m_vSCVF[i]);

		// 	integration point
			if(dim != 1)
				AveragePositions(m_vSCVF[i].localIP, m_vSCVF[i].m_vLocPos, SCVF::numCorners);
			else
				m_vSCVF[i].localIP = m_locMid[1][0];
		}

	// 	set up local informations for SubControlVolumes (scv)
	// 	each scv is associated to one corner of the element
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
				m_vSCV[i].midId[1] = MidID(1, rRefElem.id(0, i, 1, 0)); // edge 1
				m_vSCV[i].midId[2] = MidID(dim, 0);	// center of element
				m_vSCV[i].midId[3] = MidID(1, rRefElem.id(0, i, 1, 1)); // edge 2
			}
			else if(dim == 3 && (m_roid != ROID_PYRAMID || i != num_scv()-1))
			{
				m_vSCV[i].midId[0] = MidID(0, i); // set node as corner of scv
				m_vSCV[i].midId[1] = MidID(1, rRefElem.id(0, i, 1, 1)); // edge 1
				m_vSCV[i].midId[2] = MidID(2, rRefElem.id(0, i, 2, 0)); // face 0
				m_vSCV[i].midId[3] = MidID(1, rRefElem.id(0, i, 1, 0)); // edge 0
				m_vSCV[i].midId[4] = MidID(1, rRefElem.id(0, i, 1, 2)); // edge 2
				m_vSCV[i].midId[5] = MidID(2, rRefElem.id(0, i, 2, 2)); // face 2
				m_vSCV[i].midId[6] = MidID(dim, 0);	// center of element
				m_vSCV[i].midId[7] = MidID(2, rRefElem.id(0, i, 2, 1)); // face 1
			}
			// TODO: Implement last ControlVolume for Pyramid
			else if(dim == 3 && m_roid == ROID_PYRAMID && i == num_scv()-1)
			{
				// this scv has 10 corners
				m_vSCV[i].m_numCorners = 10;
				throw(UGFatalError("Last SCV for Pyramid must be implemented"));
			}
			else {throw(UGFatalError("Dimension higher that 3 not implemented."));}

		// 	copy local corners of scv
			copy_local_corners(m_vSCV[i]);
		}

		/////////////////////////
		// Shapes and Derivatives
		/////////////////////////

		try{
		const DimLocalShapeFunctionSet<dim>& TrialSpace =
			LocalShapeFunctionSetProvider::get<dim>(m_roid, LFEID(LFEID::LAGRANGE, 1));

		for(size_t i = 0; i < num_scvf(); ++i)
		{
		//	remember number of shape functions
			m_vSCVF[i].m_numSH = TrialSpace.num_sh();

			TrialSpace.shapes(&(m_vSCVF[i].vShape[0]), m_vSCVF[i].localIP);
			TrialSpace.grads(&(m_vSCVF[i].localGrad[0]), m_vSCVF[i].localIP);
		}

		for(size_t i = 0; i < num_scv(); ++i)
		{
		//	remember number of shape functions
			m_vSCV[i].m_numSH = TrialSpace.num_sh();

			TrialSpace.shapes(&(m_vSCV[i].vShape[0]), m_vSCV[i].m_vLocPos[0]);
			TrialSpace.grads(&(m_vSCV[i].localGrad[0]), m_vSCV[i].m_vLocPos[0]);
		}

		}catch(UG_ERROR_LocalShapeFunctionSetNotRegistered& ex)
		{
			UG_LOG("ERROR in 'DimFV1Geometry::update': "<<ex.get_msg()<<"\n");
			return false;
		}

		///////////////////////////
		// Copy ip pos in list
		///////////////////////////

	// 	loop Sub Control Volumes Faces (SCVF)
		for(size_t i = 0; i < num_scvf(); ++i)
		{
		//	get current SCVF
			const SCVF& rSCVF = scvf(i);

		// 	loop ips
			m_vLocSCVF_IP[i] = rSCVF.local_ip();
		}

		}catch(UG_ERROR_ReferenceElementMissing& ex)
		{
			UG_LOG("ERROR in 'DimFV1Geometry::update': "<<ex.get_msg()<<"\n");
			return false;
		}
	}

	//////////////////////////////////
	// COMPUTATION OF GLOBAL DATA
	//////////////////////////////////

//	get reference element
	try{
	const DimReferenceElement<dim>& rRefElem
		= ReferenceElementProvider::get<dim>(m_roid);

// 	remember global position of nodes
	for(size_t i = 0; i < rRefElem.num(0); ++i)
		m_gloMid[0][i] = vCornerCoords[i];

// 	compute global midpoints for all geometric objects with  0 < d <= dim
	for(int d = 1; d <= dim; ++d)
	{
	// 	loop geometric objects of dimension d
		for(size_t i = 0; i < rRefElem.num(d); ++i)
		{
		// 	set first node
			const size_t coID0 = rRefElem.id(d, i, 0, 0);
			m_gloMid[d][i] = m_gloMid[0][coID0];

		// 	add corner coordinates of the corners of the geometric object
			for(size_t j = 1; j < rRefElem.num(d, i, 0); ++j)
			{
				const size_t coID = rRefElem.id(d, i, 0, j);
				m_gloMid[d][i] += m_gloMid[0][coID];
			}

		// 	scale for correct averaging
			m_gloMid[d][i] *= 1./(rRefElem.num(d, i, 0));
		}
	}

// 	compute global informations for scvf
	for(size_t i = 0; i < num_scvf(); ++i)
	{
	// 	copy local corners of scvf
		copy_global_corners(m_vSCVF[i]);

	// 	integration point
		if(dim != 1)
		AveragePositions(m_vSCVF[i].globalIP, m_vSCVF[i].m_vGloPos, SCVF::numCorners);
		else
			m_vSCVF[i].globalIP = m_gloMid[1][0];

	// 	normal on scvf
		if(dim == worldDim)
		{
			fv1_dim_traits<dim, worldDim>::NormalOnSCVF(m_vSCVF[i].Normal, m_vSCVF[i].m_vGloPos);
		}
		else
		{
			if(dim == 1)
				fv1_dim_traits<dim, worldDim>::NormalOnSCVF(m_vSCVF[i].Normal, vCornerCoords);
			else
				throw(UGFatalError("Not Implemented"));
		}
	}

// 	compute size of scv
	for(size_t i = 0; i < num_scv(); ++i)
	{
	// 	copy global corners
		copy_global_corners(m_vSCV[i]);

	// 	compute volume of scv
		if(m_vSCV[i].m_numCorners != 10)
		{
			m_vSCV[i].vol = ElementSize<scv_type, worldDim>(m_vSCV[i].m_vGloPos);
		}
		else
		{
			// special case for pyramid, last scv
			throw(UGFatalError("Pyramid Not Implemented"));
		}
	}

	/////////////////////////
	// Shapes and Derivatives
	/////////////////////////

//	get reference mapping
	try{
	DimReferenceMapping<dim, worldDim>& rMapping = ReferenceMappingProvider::get<dim, worldDim>(m_roid);

	rMapping.update(vCornerCoords);

	//\todo compute with on virt. call
//	compute jacobian for lin
	if(rMapping.is_linear())
	{
		rMapping.jacobian_transposed_inverse(m_vSCVF[0].JtInv, m_vSCVF[0].localIP);
		rMapping.jacobian_transposed_inverse(m_vSCV[0].JtInv, m_vSCV[0].m_vLocPos[0]);
		m_vSCVF[0].detj = rMapping.jacobian_det(m_vSCVF[0].localIP);
		m_vSCV[0].detj = rMapping.jacobian_det(m_vSCV[0].m_vLocPos[0]);
		for(size_t i = 1; i < num_scvf(); ++i)
		{
			m_vSCVF[i].JtInv = m_vSCVF[0].JtInv;
			m_vSCVF[i].detj = m_vSCVF[0].detj;
		}
		for(size_t i = 1; i < num_scv(); ++i)
		{
			m_vSCV[i].JtInv = m_vSCV[0].JtInv;
			m_vSCV[i].detj = m_vSCV[0].detj;
		}
	}
//	else compute jacobian for each integration point
	else
	{
		for(size_t i = 0; i < num_scvf(); ++i)
		{
			rMapping.jacobian_transposed_inverse(m_vSCVF[i].JtInv, m_vSCVF[i].localIP);
			m_vSCVF[i].detj = rMapping.jacobian_det(m_vSCVF[i].localIP);
		}
		for(size_t i = 0; i < num_scv(); ++i)
		{
			rMapping.jacobian_transposed_inverse(m_vSCV[i].JtInv, m_vSCV[i].m_vLocPos[0]);
			m_vSCV[i].detj = rMapping.jacobian_det(m_vSCV[i].m_vLocPos[0]);
		}
	}

	}catch(UG_ERROR_ReferenceElementMissing& ex)
	{
		UG_LOG("ERROR in 'DimFV1Geometry::update': "<<ex.get_msg()<<"\n");
		return false;
	}

//	compute global gradients
	for(size_t i = 0; i < num_scvf(); ++i)
		for(size_t sh = 0 ; sh < num_scv(); ++sh)
			MatVecMult((m_vSCVF[i].globalGrad)[sh], m_vSCVF[i].JtInv, (m_vSCVF[i].localGrad)[sh]);

	for(size_t i = 0; i < num_scv(); ++i)
		for(size_t sh = 0 ; sh < num_scv(); ++sh)
			MatVecMult((m_vSCV[i].globalGrad)[sh], m_vSCV[i].JtInv, (m_vSCV[i].localGrad)[sh]);

	///////////////////////////
	// Copy ip pos in list
	///////////////////////////

// 	loop Sub Control Volumes Faces (SCVF)
	for(size_t i = 0; i < num_scvf(); ++i)
	{
	//	get current SCVF
		const SCVF& rSCVF = scvf(i);

		m_vGlobSCVF_IP[i] = rSCVF.global_ip();
	}

	/////////////////////////
	// Boundary Faces
	/////////////////////////

//	if no boundary subsets required, return
	if(num_boundary_subsets() == 0) return true;

//	get grid
	Grid& grid = *ish.get_assigned_grid();

//	vector of subset indices of side
	std::vector<int> vSubsetIndex;

//	get subset indices for sides (i.e. edge in 2d, faces in 3d)
	if(dim == 2)
	{
		std::vector<EdgeBase*> vEdges;
		CollectEdgesSorted(vEdges, grid, pElem);

		vSubsetIndex.resize(vEdges.size());
		for(size_t i = 0; i < vEdges.size(); ++i)
		{
			vSubsetIndex[i] = ish.get_subset_index(vEdges[i]);
		}
	}

	if(dim == 3)
	{
		std::vector<Face*> vFaces;
		CollectFacesSorted(vFaces, grid, pElem);

		vSubsetIndex.resize(vFaces.size());
		for(size_t i = 0; i < vFaces.size(); ++i)
		{
			vSubsetIndex[i] = ish.get_subset_index(vFaces[i]);
		}
	}

//	loop requested subset
	typename std::map<int, std::vector<BF> >::iterator it;
	for (it=m_mapVectorBF.begin() ; it != m_mapVectorBF.end(); ++it)
	{
	//	get subset index
		const int bndIndex = (*it).first;

	//	get vector of BF for element
		std::vector<BF>& vBF = (*it).second;

	//	clear vector
		vBF.clear();

	//	current number of bf
		size_t curr_bf = 0;

	//	loop sides of element
		for(size_t side = 0; side < vSubsetIndex.size(); ++side)
		{
		//	skip non boundary sides
			if(vSubsetIndex[side] != bndIndex) continue;

		///////////////////////////
		//	add Boundary faces
		///////////////////////////

		//	number of corners of side
			const int coOfSide = rRefElem.num(dim-1, side, 0);

		//	resize vector
			vBF.resize(curr_bf + coOfSide);

		//	loop corners
			for(int co = 0; co < coOfSide; ++co)
			{
			//	get current bf
				BF& bf = vBF[curr_bf];

			//	set node id == scv this bf belongs to
				bf.m_nodeId = rRefElem.id(dim-1, side, 0, co);

			// 	set mid ids
				if(dim == 2)
				{
					bf.midId[co%2] = MidID(0, rRefElem.id(1, side, 0, co)); // corner of side
					bf.midId[(co+1)%2] = MidID(1, side); // side midpoint
				}
				else if (dim == 3)
				{
					bf.midId[0] = MidID(0, rRefElem.id(2, side, 0, co)); // corner of side
					bf.midId[1] = MidID(1, rRefElem.id(2, side, 1, co)); // edge co
					bf.midId[2] = MidID(2, side); // side midpoint
					bf.midId[3] = MidID(1, rRefElem.id(2, side, 1, (co -1 + coOfSide)%coOfSide)); // edge co-1
				}

			// 	copy corners of bf
				copy_local_corners(bf);
				copy_global_corners(bf);

			// 	integration point
				AveragePositions(bf.localIP, bf.m_vLocPos, SCVF::numCorners);
				AveragePositions(bf.globalIP, bf.m_vGloPos, SCVF::numCorners);

			// 	normal on scvf
				fv1_dim_traits<dim, worldDim>::NormalOnSCVF(bf.Normal, bf.m_vGloPos);

			//	compute volume
				bf.m_volume = VecTwoNorm(bf.Normal);

			/////////////////////////
			// Shapes and Derivatives
			/////////////////////////

				//\todo: perform not here
				try{
				const DimLocalShapeFunctionSet<dim>& TrialSpace =
					LocalShapeFunctionSetProvider::get<dim>(m_roid, LFEID(LFEID::LAGRANGE, 1));

			//	remember number of shape functions
				bf.m_numSH = TrialSpace.num_sh();

				TrialSpace.shapes(&(bf.vShape[0]), bf.localIP);
				TrialSpace.grads(&(bf.localGrad[0]), bf.localIP);

				}catch(UG_ERROR_LocalShapeFunctionSetNotRegistered& ex)
				{
					UG_LOG("ERROR in 'DimFV1Geometry::update': "<<ex.get_msg()<<"\n");
					return false;
				}

			//	get reference mapping
				try{
				DimReferenceMapping<dim, worldDim>& rMapping = ReferenceMappingProvider::get<dim, worldDim>(m_roid);
				rMapping.update(vCornerCoords);

				rMapping.jacobian_transposed_inverse(bf.JtInv, bf.localIP);
				bf.detj = rMapping.jacobian_det(bf.localIP);

				}catch(UG_ERROR_ReferenceElementMissing& ex)
				{
					UG_LOG("ERROR in 'DimFV1Geometry::update': "<<ex.get_msg()<<"\n");
					return false;
				}

			//	compute global gradients
				for(size_t sh = 0 ; sh < num_scv(); ++sh)
					MatVecMult((bf.globalGrad)[sh], bf.JtInv, (bf.localGrad)[sh]);

			//	increase curr_bf
				++curr_bf;
			}
		}
	}

	}catch(UG_ERROR_ReferenceElementMissing& ex)
	{
		UG_LOG("ERROR in 'DimFV1Geometry::update'"<<ex.get_msg()<<"\n");
		return false;
	}

//	done
	return true;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// FV1ManifoldBoundary
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <typename TElem, int TWorldDim>
FV1ManifoldBoundary<TElem, TWorldDim>::
FV1ManifoldBoundary() : m_pElem(NULL), m_rRefElem(Provider::get<ref_elem_type>())
{
	// set corners of element as local centers of nodes
	for (size_t i = 0; i < m_rRefElem.num(0); ++i)
		m_locMid[0][i] = m_rRefElem.corner(i);

	// compute local midpoints for all geometric objects with  0 < d <= dim
	for (int d = 1; d <= dim; ++d)
	{
		// loop geometric objects of dimension d
		for(size_t i = 0; i < m_rRefElem.num(d); ++i)
		{
			// set first node
			const size_t coID0 = m_rRefElem.id(d, i, 0, 0);
			m_locMid[d][i] = m_locMid[0][coID0];

			// add corner coordinates of the corners of the geometric object
			for(size_t j = 1; j < m_rRefElem.num(d, i, 0); ++j)
			{
				const size_t coID = m_rRefElem.id(d, i, 0, j);
				m_locMid[d][i] += m_locMid[0][coID];
			}

			// scale for correct averaging
			m_locMid[d][i] *= 1./(m_rRefElem.num(d, i, 0));
		}
	}

	// set up local information for Boundary Faces (bf)
	// each bf is associated to one corner of the element
	for (size_t i = 0; i < num_bf(); ++i)
	{
		m_vBF[i].nodeId = i;

		if (dim == 1) // Edge
		{
			m_vBF[i].midId[0] = MidID(0, i);	// set node as corner of bf
			m_vBF[i].midId[1] = MidID(dim, 0);	// center of bnd element
			
			// copy local corners of bf
			copy_local_corners(m_vBF[i]);
			
			// local integration point
			AveragePositions(m_vBF[i].localIP, m_vBF[i].m_vLocPos, 2);
		}
		else if (dim == 2)	// Quadrilateral
		{
			m_vBF[i].midId[0] = MidID(0, i); // set node as corner of bf
			m_vBF[i].midId[1] = MidID(1, m_rRefElem.id(0, i, 1, 0)); // edge 1
			m_vBF[i].midId[2] = MidID(dim, 0);	// center of bnd element
			m_vBF[i].midId[3] = MidID(1, m_rRefElem.id(0, i, 1, 1)); // edge 2
			
			// copy local corners of bf
			copy_local_corners(m_vBF[i]);
			
			// local integration point
			AveragePositions(m_vBF[i].localIP, m_vBF[i].m_vLocPos, 4);
		}
		else {UG_ASSERT(0, "Dimension higher that 2 not implemented.");}
	}

	/////////////
	// Shapes
	/////////////
	// A word of warning: This is only meaningful,
	// if the trial space is piecewise linear on tetrahedrons/triangles!
	for (size_t i = 0; i < num_bf(); ++i)
	{
		const LocalShapeFunctionSet<ref_elem_type>& TrialSpace =
				LocalShapeFunctionSetProvider::
					get<ref_elem_type>
					(LFEID(LFEID::LAGRANGE, 1));

		const size_t num_sh = m_numBF;
		m_vBF[i].vShape.resize(num_sh);

		TrialSpace.shapes(&(m_vBF[i].vShape[0]), m_vBF[i].localIP);
	}

	///////////////////////////
	// Copy ip pos in list
	///////////////////////////

	// 	loop Boundary Faces (BF)
	m_vLocBFIP.clear();
	for (size_t i = 0; i < num_bf(); ++i)
	{
	//	get current BF
		const BF& rBF = bf(i);

	// 	loop ips
		for (size_t ip = 0; ip < rBF.num_ip(); ++ip)
		{
			m_vLocBFIP.push_back(rBF.local_ip(ip));
		}
	}
}





/// update data for given element
template <typename TElem, int TWorldDim>
bool
FV1ManifoldBoundary<TElem, TWorldDim>::
update(TElem* elem, const ISubsetHandler& ish, const MathVector<worldDim>* vCornerCoords)
{
	// 	If already update for this element, do nothing
	if (m_pElem == elem) return true;
	else m_pElem = elem;

	// 	remember global position of nodes
	for (size_t i = 0; i < m_rRefElem.num(0); ++i)
		m_gloMid[0][i] = vCornerCoords[i];

	// 	compute global midpoints for all the other geometric objects (with  0 < d <= dim)
	for (int d = 1; d <= dim; ++d)
	{
		// 	loop geometric objects of dimension d
		for (size_t i = 0; i < m_rRefElem.num(d); ++i)
		{
			// set first node
			const size_t coID0 = m_rRefElem.id(d, i, 0, 0);
			m_gloMid[d][i] = m_gloMid[0][coID0];

		// 	add corner coordinates of the corners of the geometric object
			for (size_t j = 1; j < m_rRefElem.num(d, i, 0); ++j)
			{
				const size_t coID = m_rRefElem.id(d, i, 0, j);
				m_gloMid[d][i] += m_gloMid[0][coID];
			}

		// 	scale for correct averaging
			m_gloMid[d][i] *= 1./(m_rRefElem.num(d, i, 0));
		}
	}
	
	// set local integration points
	for (size_t i = 0; i < num_bf(); ++i)
	{
		// copy global corners of bf
		copy_global_corners(m_vBF[i]);
		
		if (dim == 1) // Edge
			{AveragePositions(m_vBF[i].localIP, m_vBF[i].m_vLocPos, 2);}
		else if (dim == 2)	// Quadrilateral
			{AveragePositions(m_vBF[i].localIP, m_vBF[i].m_vLocPos, 4);}
		else {UG_ASSERT(0, "Dimension higher that 2 not implemented.");}
	}
	
	// 	compute size of BFs
	for (size_t i = 0; i < num_bf(); ++i)
	{
	// 	copy global corners
		copy_global_corners(m_vBF[i]);

	// 	compute volume of bf
		m_vBF[i].vol = ElementSize<bf_type, worldDim>(m_vBF[i].m_vGloPos);
	}
	
	///////////////////////////
	// Copy ip pos in list
	///////////////////////////

	// 	loop Boundary Faces (BF)
	m_vGlobBFIP.clear();
	for (size_t i = 0; i < num_bf(); ++i)
	{
	//	get current BF
		const BF& rBF = bf(i);

	// 	loop ips
		for (size_t ip = 0; ip < rBF.num_ip(); ++ip)
		{
			m_vGlobBFIP.push_back(rBF.global_ip(ip));
		}
	}

	return true;
}


template class FV1Geometry<Edge, 1>;
template class FV1Geometry<Edge, 2>;
template class FV1Geometry<Edge, 3>;

template class FV1Geometry<Triangle, 2>;
template class FV1Geometry<Triangle, 3>;

template class FV1Geometry<Quadrilateral, 2>;
template class FV1Geometry<Quadrilateral, 3>;

template class FV1Geometry<Tetrahedron, 3>;
template class FV1Geometry<Prism, 3>;
template class FV1Geometry<Pyramid, 3>;
template class FV1Geometry<Hexahedron, 3>;


template class DimFV1Geometry<1, 1>;
template class DimFV1Geometry<1, 2>;
template class DimFV1Geometry<1, 3>;

template class DimFV1Geometry<2, 2>;
template class DimFV1Geometry<2, 3>;

template class DimFV1Geometry<3, 3>;


template class FV1ManifoldBoundary<Edge, 1>;
template class FV1ManifoldBoundary<Edge, 2>;
template class FV1ManifoldBoundary<Triangle, 2>;
template class FV1ManifoldBoundary<Quadrilateral, 2>;
template class FV1ManifoldBoundary<Tetrahedron, 3>;
template class FV1ManifoldBoundary<Prism, 3>;
template class FV1ManifoldBoundary<Pyramid, 3>;
template class FV1ManifoldBoundary<Hexahedron, 3>;

} // end namespace ug
