/*
 * interface_handler_local_tools.h
 *
 *  Created on: 19.01.2015
 *      Author: suze
 */

#ifndef INTERFACE_HANDLER_LOCAL_PARTICLE_TOOLS_H_
#define INTERFACE_HANDLER_LOCAL_PARTICLE_TOOLS_H_

namespace ug{
 

template <int TWorldDim>
bool InterfaceHandlerLocalParticle<TWorldDim>::
is_boundary_face_for2(const size_t sideID)
{
// REMARK: not finally tested!
    UG_THROW("InterfaceHandlerLocalParticle::is_boundary_face_for2(): not finally tested!!!!\n");

	if ( dim == 3 )
		UG_THROW("in 'is_boundary_face_for2()': only implemented for 2d!!\n");

// get data
	const DimReferenceElement<dim>& rRefElem
		= ReferenceElementProvider::get<dim>(this->m_roid);

	//	number of corners of side (special case bottom side pyramid)
		const int coOfSide = (this->m_roid != ROID_PYRAMID || sideID != 0)
							? rRefElem.num(dim-1, sideID, 0) : rRefElem.num(dim-1, sideID, 0) + 2;

	if ( coOfSide > 2 )
		UG_THROW("in 'is_boundary_face_for2()': hmm...coOfSide > 2 ??? \n coOfSide = " << coOfSide << "\n");

	std::vector<int> prtIndex(2);

	for(int co = 0; co < coOfSide; ++co)
	{
		size_t cornerID;
		if (this->m_roid != ROID_PYRAMID || sideID != 0)
			cornerID = rRefElem.id(dim-1, sideID, 0, co);
		else
			cornerID = rRefElem.id(dim-1, sideID, 0, (co % 3) + (co>3 ? 1 : 0));

		if ( !this->lies_onInterface(cornerID) )
			return false;
		else
			prtIndex[co] = getPrtIndex(co);

	}

  
	if (  prtIndex[0] ==  prtIndex[1] )
	{
        UG_THROW("========= is_boundary_face = TRUE\n");
		return true;
	}
	UG_LOG("========= is_boundary_face = FALSE\n");

	return false;
}


template<int TWorldDim>
void InterfaceHandlerLocalParticle<TWorldDim>::
update_inner_boundary_radial(const MathVector<TWorldDim>* vCornerCoords)
{
// some check
	if (dim == 2) {
		if (this->m_vBF.size() != 2 && this->m_vBF.size() != 0)
			UG_THROW("Error in 'update_inner_boundary_radial()': this->m_vBF.size() != 2:"
                     << this->m_vBF.size() << "\n");
	}
    const int prtIndex = get_prtIndex();
    if (prtIndex == -1)
        UG_THROW("InterfaceHandlerLocalParticle::update_inner_boundary_radial(): prtIndex not set!\n");

// get data
// we choose the modusIP = 0 and useResized = true
//  --> see implemenation and cases below
	const int modusIP = 0;
	const bool useResized = true;
	const MathVector<dim> center = m_spCutElementHandler->get_center(prtIndex);

// reset data
	m_vRadialAtCo.clear();
	m_vRadialAtCo.resize(this->numCo(), 0.0);
	m_vRadialAtIP.clear();
	m_vRadialAtIP.resize(this->numCo(), 0.0);

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// REMARK: for elements with only ONE outside corner ON the interface, bf.size() = 0!
// 			=> for assembling connections u_fluid -> u_prt, although 'm_vRadialAtCo' is needed!
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// A. set 'm_vRadialAtCo' --> independent of boundary faces!
	for (size_t i = 0; i < this->numCo(); ++i) {
		if (useResized) {
			VecSubtract(m_vRadialAtCo[i], this->corner(i), center);
        } else { UG_THROW("in 'update_inner_boundary_radial()': case not implemented anymore!\n");}

	}

// B. set 'm_vRadialAtIP' --> for boundary faces!
	for (size_t i = 0; i < this->m_vBF.size(); ++i) {
		interfaceBF bf = this->m_vBF[i];

		size_t nodeID = bf.node_id();
		if (!useResized)
			nodeID = this->corner_orig(nodeID);

        typedef DimFV1FTGeometry<dim, dim, InterfaceHandlerLocalParticle<dim> > TFVGeom;
		TFVGeom& geo = GeomProvider<TFVGeom>::get(LFEID(LFEID::LAGRANGE, dim, 1), 1);

		int ip = (bf.node_id() - 1);
		if (ip < 0)
			ip = ip + 3;

		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
		const MathVector<dim>& SCVFip_Pos = scvf.global_ip();

		switch (modusIP) {
		case 0: 	// AtIP == AtCo
			m_vRadialAtIP[nodeID] = m_vRadialAtCo[nodeID];
			break;
		case 1: 	// AtIP == bf.global_ip()
			VecSubtract(m_vRadialAtIP[nodeID], bf.global_ip(), center);
			//		VecNormalize(m_vRadialAtIP[nodeID], m_vRadialAtIP[nodeID]);
			//		VecScale(m_vRadialAtIP[nodeID], m_vRadialAtIP[nodeID], radius);
			break;
		case 2: 	// AtIP == bf.normal()
            UG_THROW("in 'update_inner_boundary_radial()': case not implemented anymore!\n");
            break;
		case 3:
			VecSubtract(m_vRadialAtIP[nodeID], SCVFip_Pos, center);
			//	VecNormalize(m_vRadialAtIP[ip], m_vRadialAtIP[ip]);
			//	VecScale(m_vRadialAtIP[ip], m_vRadialAtIP[ip], VecLength(bf.normal()));
			break;
		default:
			throw(UGError(
					"Error in IInterfaceMapper::update_inner_boundary_radial()!"));
		}

	}

}

template<int TWorldDim>
void InterfaceHandlerLocalParticle<TWorldDim>::
update_inner_boundary_radial_for2_StdFV()
{
// REMARK: not finally tested!
    UG_THROW("InterfaceHandlerLocalParticle::update_inner_boundary_radial_for2_StdFV(): not finally tested!!!!\n");

    const int prtIndex = get_prtIndex();

// some check
	if (prtIndex == -1)
		UG_THROW("InterfaceHandlerLocalParticle::update_inner_boundary_radial_for2_StdFV(): prtIndex not set!\n");

	const MathVector<dim> center1 = m_spCutElementHandler->get_center(0);
 	const MathVector<dim> center2 = m_spCutElementHandler->get_center(1);
 
	UG_LOG("center1: " << center1 << "\n");
	UG_LOG("center2: " << center2 << "\n");

// reset data
	m_vRadialAtCo.clear();
	m_vRadialAtCo.resize(this->numCo(), 0.0);
	m_vRadialAtIP.clear();
	m_vRadialAtIP.resize(this->numCo(), 0.0);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// REMARK: for elements with only ONE outside corner ON the interface, bf.size() = 0!
	// 			=> for assembling connections u_fluid -> u_prt, although 'm_vRadialAtCo' is needed!
	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// A. set 'm_vRadialAtCo' --> independent of boundary faces!
	// fill m_vRadialAtCO for center1:

	MathVector<dim> RadialIPRef1;
	MathVector<dim> RadialIPRef2;
	RadialIPRef1 = 0.0;
	RadialIPRef2 = 0.0;
	size_t numOutCo1 = 0;
	size_t numOutCo2 = 0;

	if (dim == 2 && this->numCo() != 3)
		UG_THROW(
				"' in update_inner_boundary_radial_for2_StdFV(): wrong number of corners! -> numCo = " << this->numCo() << "\n");

	for (size_t i = 0; i < this->numCo(); ++i) {
		const int prtIndex = getPrtIndex(i);

    // set 'm_vRadialAtCo'
		if (prtIndex == 0) {
			VecSubtract(m_vRadialAtCo[i], this->corner(i), center1);
			RadialIPRef1 += this->corner(i);
			numOutCo1++;

			if (!this->lies_onInterface(i))
				UG_THROW(
						"inconsistent in 'update_inner_boundary_radial_for2_StdFV()'!\n");
		} else if (prtIndex == 1) {
			VecSubtract(m_vRadialAtCo[i], this->corner(i), center2);
			RadialIPRef2 += this->corner(i);
			numOutCo2++;

			if (!this->lies_onInterface(i))
				UG_THROW(
						"inconsistent in 'update_inner_boundary_radial_for2_StdFV()'!\n");
		}

	}

// set 'm_vRadialAtIP'
	VecScale(RadialIPRef1, RadialIPRef1, 1.0 / numOutCo1);
	VecScale(RadialIPRef2, RadialIPRef2, 1.0 / numOutCo2);
	UG_LOG("RadialIPRef1 = " << RadialIPRef1 << "\n");
	UG_LOG("RadialIPRef2 = " << RadialIPRef2 << "\n");

	for (size_t i = 0; i < this->numCo(); ++i) {
		const int prtIndex = getPrtIndex(i);

		if (prtIndex == 0) {
			VecSubtract(m_vRadialAtIP[i], RadialIPRef1, center1);
 		}
        else if (prtIndex == 1) {
			VecSubtract(m_vRadialAtIP[i], RadialIPRef2, center2);
 		}

	}

}

// VORSICHT: update_inner_boundary_radial() wird in fv1FT_geom_impl.h aufgerufen, BEVOR remap_BF() stattfindet!
//	=> bf.node_id() immer mit interfaceID, nicht OriginalCornerID!
template<int TWorldDim>
void InterfaceHandlerLocalParticle<TWorldDim>::
update_inner_boundary_radial_for2()
{
// REMARK: not finally tested!
    UG_THROW("InterfaceHandlerLocalParticle::update_inner_boundary_radial_for2(): not finally tested!!!!\n");

// get data:
	std::vector<MathVector<dim> > vCornerCoords;
	for (size_t i = 0; i < 4; ++i)
		vCornerCoords.push_back(this->m_vQuadriCorners_for2[i]);

// we choose the modusIP = 0 and useResized = true
//  --> see implemenation and cases below
	const int modusIP = 1;
	const bool useResized = true;

    const int prtIndex = get_prtIndex();

// some Check
	if (prtIndex == -1)
		UG_THROW("InterfaceHandlerLocalParticle::update_inner_boundary_radial(): prtIndex not set!\n");

	if (this->m_vBF.size() != 4)
		UG_THROW("in 'update_inner_boundary_radial_for2()': this->m_vBF.size() should be 4, but is "
                 << this->m_vBF.size() << "\n");

	const MathVector<dim> center1 = m_spCutElementHandler->get_center(0);
 	const MathVector<dim> center2 = m_spCutElementHandler->get_center(1);
 
	UG_LOG("center1: " << center1 << "\n");
	UG_LOG("center2: " << center2 << "\n");

// reset data
	m_vRadialAtCo.clear();
	m_vRadialAtCo.resize(4, 0.0);
	m_vRadialAtIP.clear();
	m_vRadialAtIP.resize(4, 0.0);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// REMARK: for elements with only ONE outside corner ON the interface, bf.size() = 0!
	// 			=> for assembling connections u_fluid -> u_prt, although 'm_vRadialAtCo' is needed!
	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// A. set 'm_vRadialAtCo' --> independent of boundary faces!
	// fill m_vRadialAtCO for center1:
	for (size_t i = 0; i < 2; ++i) {
		if (useResized)
			VecSubtract(m_vRadialAtCo[i], vCornerCoords[i], center1);
		else
			UG_THROW("in 'update_inner_boundary_radial_for2()': not implemented in the for2-case!\n");

 	}
	// fill m_vRadial for center2:
	for (size_t i = 2; i < 4; ++i) {
  		if (useResized)
			VecSubtract(m_vRadialAtCo[i], vCornerCoords[i], center2);
		else
			UG_THROW("in 'update_inner_boundary_radial_for2()': not implemented in the for2-case!\n");

		UG_LOG("-> vCornerCoords[" << i << "] = " << vCornerCoords[i] << "\n");
		UG_LOG("-> m_vRadialAtCo[" << i << "] = " << m_vRadialAtCo[i] << "\n");
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// B. set 'm_vRadialAtIP' --> for boundary faces!
	// fill m_vRadialAtCO for center1:
	for (size_t i = 0; i < 2; ++i) {
		interfaceBF bf = this->m_vBF[i];

		size_t nodeID = bf.node_id(); // ok, since bf.nodeID = 0,1,2,3
		if (!useResized)
			UG_THROW("in 'update_inner_boundary_radial_for2()': case '!useResized' not implemented!\n");

		int ip = (bf.node_id() - 1);
		if (ip < 0)
			ip = ip + 3;

		switch (modusIP) {
		case 0: 	// AtIP == AtCo
			m_vRadialAtIP[nodeID] = m_vRadialAtCo[nodeID];
			break;
		case 1: 	// AtIP == bf.global_ip()
			VecSubtract(m_vRadialAtIP[nodeID], bf.global_ip(), center1);
			break;
		default:
			throw(UGError(
					"Error in InterfaceHandlerLocalParticle::update_inner_boundary_radial_for2()!"));
		}

	}

	// fill m_vRadialAtCO for center2:
	for (size_t i = 2; i < 4; ++i) {
		interfaceBF bf = this->m_vBF[i];

		size_t nodeID = bf.node_id(); // ok, since bf.nodeID = 0,1,2,3
		if (!useResized)
			UG_THROW("in 'update_inner_boundary_radial_for2()': case '!useResized' not implemented!\n");

		int ip = (bf.node_id() - 1);
		if (ip < 0)
			ip = ip + 3;

		switch (modusIP) {
		case 0: 	// AtIP == AtCo
			m_vRadialAtIP[nodeID] = m_vRadialAtCo[nodeID];
			break;
		case 1: 	// AtIP == bf.global_ip()
			VecSubtract(m_vRadialAtIP[nodeID], bf.global_ip(), center2);
			break;
		default:
			throw(UGError(
					"Error in InterfaceHandlerLocalParticle::update_inner_boundary_radial_for2()!"));
		}

	}
 
}

template<int TWorldDim>
void InterfaceHandlerLocalParticle<TWorldDim>::
update_inner_boundary_radial_StdFV(const MathVector<TWorldDim>* vCornerCoords)
{
    const int prtIndex = get_prtIndex();

	if (prtIndex == -1)
		UG_THROW("InterfaceHandlerLocalParticle::update_inner_boundary_radial(): prtIndex not set!\n");

	const MathVector<dim> center = m_spCutElementHandler->get_center(prtIndex);

	/*
	 typedef DimFV1FTGeometry<dim, dim, InterfaceHandlerLocalParticle<dim> > TFVGeom;
	 TFVGeom& geo = GeomProvider<TFVGeom>::get(LFEID(LFEID::LAGRANGE, dim, 1), 1);
	 */
// reset data
	m_vRadialAtCo.clear();
	m_vRadialAtCo.resize(this->numCo(), 0.0);
	m_vRadialAtIP.clear();
	m_vRadialAtIP.resize(this->numCo(), 0.0);

//  Loop the boundary faces to assemble impulse equations
	MathVector<dim> RadialIPRef;
	RadialIPRef = 0.0;
	size_t numOutCo = 0;
	for (size_t i = 0; i < this->numCo(); ++i) {
		// set 'm_vRadialAtCo'
		VecSubtract(m_vRadialAtCo[i], this->corner(i), center);

		// set 'm_vRadialAtIP'
		/*		const typename TFVGeom::SCVF& scvf = geo.scvf(i);
		 const MathVector<dim>& SCVFip_Pos = scvf.global_ip();
		 VecSubtract(m_vRadialAtIP[i], SCVFip_Pos, center);
		 */
		//VecSubtract(m_vRadialAtIP[i], this->corner(i), center);
		// compute radialAtIP for current element:
		if (this->lies_onInterface(i)) {
			RadialIPRef += this->corner(i);
			numOutCo++;
		}

	}

// set 'm_vRadialAtIP'
	VecScale(RadialIPRef, RadialIPRef, 1.0 / numOutCo);
	for (size_t i = 0; i < this->numCo(); ++i)
		VecSubtract(m_vRadialAtIP[i], RadialIPRef, center);

}

} // end namespace ug



#endif /* INTERFACE_HANDLER_LOCAL_PARTICLE_TOOLS_H_ */
