/*
 * support3D.h
 *
 *  Created on: 31.10.2024
 *      Author: Markus Knodel
 */

#ifndef UGCORE_UGBASE_LIB_GRID_ALGORITHMS_EXTRUSION_SUPPORT3D_H_
#define UGCORE_UGBASE_LIB_GRID_ALGORITHMS_EXTRUSION_SUPPORT3D_H_

#include <iostream>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <assert.h>
#include <string>
#include <sstream>
#include <utility>
#include <vector>
#include <type_traits>

#include "support.h"

namespace ug
{

namespace support
{
template<typename ELEMTYP, typename INDEX_TYP>
class ElemInfo
{
public:

	ElemInfo( ELEMTYP const & elem, INDEX_TYP sudo )
	: m_elem(elem), m_sudo(sudo)
	{
	}



private:

	ELEMTYP m_elem;
	INDEX_TYP m_sudo;

	ElemInfo() {};


};

template<
typename FACETYP,
typename NORMALTYP,
typename VOLUMETYP,
typename EDGETYP,
typename INDEX_TYP
>
class VertexFractureQuadrupel
{
public:

	using ElemInfoFac = ElemInfo<FACETYP,INDEX_TYP>;

	//face, normal, volume, edge

//	VertexFractureQuadrupel()
//	{};


	VertexFractureQuadrupel( ElemInfoFac const & fracFaceInfo,
							 VOLUMETYP const & attVolume,
							 NORMALTYP const & normal,
							 std::pair<EDGETYP,EDGETYP> const & volCutEdges,
							 std::pair<ElemInfoFac,ElemInfoFac> const & volCutEdgeFaces )
	: m_fracFaceInfo(fracFaceInfo),
	  m_attVolume(attVolume),
	  m_normal(normal),
	  m_volCutEdges(volCutEdges),
	  m_volCutEdgeFaces(volCutEdgeFaces)
	{
	}

	// todo fixme getter und ggf auch setter, aber vermutlich nur getter implementieren!!!

private:

//	FACETYP const getFace() const { return m_face; }
//	NORMALTYP const getNormal() const { return m_normal; }
//	VOLUMETYP const getVolume() const { return m_volume; }
//	EDGETYP const getEdge() const { return m_edge; }

	ElemInfo<FACETYP,INDEX_TYP> m_fracFaceInfo;
	VOLUMETYP m_attVolume;
	NORMALTYP m_normal;
	std::pair<EDGETYP,EDGETYP> m_volCutEdges;
	std::pair<ElemInfoFac,ElemInfoFac> m_volCutEdgeFaces;

//private:
//
//	FACETYP m_face;
//	NORMALTYP m_normal;
//	VOLUMETYP m_volume;
//	EDGETYP m_edge;

	VertexFractureQuadrupel()
	{};
};

}

}

#endif /* UGCORE_UGBASE_LIB_GRID_ALGORITHMS_EXTRUSION_SUPPORT3D_H_ */
