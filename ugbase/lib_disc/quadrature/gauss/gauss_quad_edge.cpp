//  This file is parsed from UG 3.9.
//  It provides the Gauss Quadratures for a reference edge.


#include "../quadrature.h"
#include "gauss_quad_edge.h"
#include "common/util/provider.h"

namespace ug{

template <>
number GaussQuadBase<GaussQuadrature<ReferenceEdge, 1>, 1, 1, 1>::m_vWeight[1] =
{
 1.00000000000000000000000000000000
};

template <>
MathVector<1> GaussQuadBase<GaussQuadrature<ReferenceEdge, 1>, 1, 1, 1>::m_vPoint[1] =
{
MathVector<1>(0.50000000000000000000000000000000)
};

template <>
number GaussQuadBase<GaussQuadrature<ReferenceEdge, 3>, 1, 3, 2>::m_vWeight[2] =
{
 0.50000000000000000000000000000000,
 0.50000000000000000000000000000000
};

template <>
MathVector<1> GaussQuadBase<GaussQuadrature<ReferenceEdge, 3>, 1, 3, 2>::m_vPoint[2] =
{
MathVector<1>(0.21132486540518711774542560974902),
MathVector<1>(0.78867513459481288225457439025098)
};

template <>
number GaussQuadBase<GaussQuadrature<ReferenceEdge, 5>, 1, 5, 3>::m_vWeight[3] =
{
 0.27777777777777777777777777777778,
 0.44444444444444444444444444444444,
 0.27777777777777777777777777777778
};

template <>
MathVector<1> GaussQuadBase<GaussQuadrature<ReferenceEdge, 5>, 1, 5, 3>::m_vPoint[3] =
{
MathVector<1>(0.11270166537925831148207346002176),
MathVector<1>(0.50000000000000000000000000000000),
MathVector<1>(0.88729833462074168851792653997824)
};

template <>
number GaussQuadBase<GaussQuadrature<ReferenceEdge, 7>, 1, 7, 4>::m_vWeight[4] =
{
 0.17392742256872692868653197461100,
 0.32607257743127307131346802538900,
 0.32607257743127307131346802538900,
 0.17392742256872692868653197461100
};

template <>
MathVector<1> GaussQuadBase<GaussQuadrature<ReferenceEdge, 7>, 1, 7, 4>::m_vPoint[4] =
{
MathVector<1>(0.06943184420297371238802675555360),
MathVector<1>(0.33000947820757186759866712044838),
MathVector<1>(0.66999052179242813240133287955162),
MathVector<1>(0.93056815579702628761197324444640)
};

template <>
number GaussQuadBase<GaussQuadrature<ReferenceEdge, 9>, 1, 9, 5>::m_vWeight[5] =
{
 0.11846344252809454375713202035996,
 0.23931433524968323402064575741782,
 0.28444444444444444444444444444444,
 0.23931433524968323402064575741782,
 0.11846344252809454375713202035996
};

template <>
MathVector<1> GaussQuadBase<GaussQuadrature<ReferenceEdge, 9>, 1, 9, 5>::m_vPoint[5] =
{
MathVector<1>(0.04691007703066800360118656085030),
MathVector<1>(0.23076534494715845448184278964990),
MathVector<1>(0.50000000000000000000000000000000),
MathVector<1>(0.76923465505284154551815721035010),
MathVector<1>(0.95308992296933199639881343914970)
};

template <>
number GaussQuadBase<GaussQuadrature<ReferenceEdge, 11>, 1, 11, 6>::m_vWeight[6] =
{
 0.08566224618958517252014807108637,
 0.18038078652406930378491675691886,
 0.23395696728634552369493517199478,
 0.23395696728634552369493517199478,
 0.18038078652406930378491675691886,
 0.08566224618958517252014807108637
};

template <>
MathVector<1> GaussQuadBase<GaussQuadrature<ReferenceEdge, 11>, 1, 11, 6>::m_vPoint[6] =
{
MathVector<1>(0.03376524289842398609384922275300),
MathVector<1>(0.16939530676686774316930020249005),
MathVector<1>(0.38069040695840154568474913915964),
MathVector<1>(0.61930959304159845431525086084036),
MathVector<1>(0.83060469323313225683069979750995),
MathVector<1>(0.96623475710157601390615077724700)
};

template <>
number GaussQuadBase<GaussQuadrature<ReferenceEdge, 13>, 1, 13, 7>::m_vWeight[7] =
{
 0.06474248308443484663530571633954,
 0.13985269574463833395073388571189,
 0.19091502525255947247518488774449,
 0.20897959183673469387755102040816,
 0.19091502525255947247518488774449,
 0.13985269574463833395073388571189,
 0.06474248308443484663530571633954
};

template <>
MathVector<1> GaussQuadBase<GaussQuadrature<ReferenceEdge, 13>, 1, 13, 7>::m_vPoint[7] =
{
MathVector<1>(0.02544604382862073773690515797607),
MathVector<1>(0.12923440720030278006806761335961),
MathVector<1>(0.29707742431130141654669679396152),
MathVector<1>(0.50000000000000000000000000000000),
MathVector<1>(0.70292257568869858345330320603848),
MathVector<1>(0.87076559279969721993193238664039),
MathVector<1>(0.97455395617137926226309484202393)
};

template <>
number GaussQuadBase<GaussQuadrature<ReferenceEdge, 15>, 1, 15, 8>::m_vWeight[8] =
{
 0.05061426814518812957626567715498,
 0.11119051722668723527217799721312,
 0.15685332293894364366898110099330,
 0.18134189168918099148257522463860,
 0.18134189168918099148257522463860,
 0.15685332293894364366898110099330,
 0.11119051722668723527217799721312,
 0.05061426814518812957626567715498
};

template <>
MathVector<1> GaussQuadBase<GaussQuadrature<ReferenceEdge, 15>, 1, 15, 8>::m_vPoint[8] =
{
MathVector<1>(0.01985507175123188415821956571526),
MathVector<1>(0.10166676129318663020422303176208),
MathVector<1>(0.23723379504183550709113047540538),
MathVector<1>(0.40828267875217509753026192881991),
MathVector<1>(0.59171732124782490246973807118009),
MathVector<1>(0.76276620495816449290886952459462),
MathVector<1>(0.89833323870681336979577696823792),
MathVector<1>(0.98014492824876811584178043428474)
};

template <>
number GaussQuadBase<GaussQuadrature<ReferenceEdge, 17>, 1, 17, 9>::m_vWeight[9] =
{
 0.04063719418078720598594607905526,
 0.09032408034742870202923601562146,
 0.13030534820146773115937143470932,
 0.15617353852000142003431520329222,
 0.16511967750062988158226253464349,
 0.15617353852000142003431520329222,
 0.13030534820146773115937143470932,
 0.09032408034742870202923601562146,
 0.04063719418078720598594607905526
};

template <>
MathVector<1> GaussQuadBase<GaussQuadrature<ReferenceEdge, 17>, 1, 17, 9>::m_vPoint[9] =
{
MathVector<1>(0.01591988024618695508221189854816),
MathVector<1>(0.08198444633668210285028510596513),
MathVector<1>(0.19331428364970480134564898032926),
MathVector<1>(0.33787328829809553548073099267833),
MathVector<1>(0.50000000000000000000000000000000),
MathVector<1>(0.66212671170190446451926900732167),
MathVector<1>(0.80668571635029519865435101967074),
MathVector<1>(0.91801555366331789714971489403487),
MathVector<1>(0.98408011975381304491778810145184)
};

template <>
number GaussQuadBase<GaussQuadrature<ReferenceEdge, 19>, 1, 19, 10>::m_vWeight[10] =
{
 0.03333567215434406879678440494667,
 0.07472567457529029657288816982885,
 0.10954318125799102199776746711408,
 0.13463335965499817754561346078473,
 0.14776211235737643508694649732567,
 0.14776211235737643508694649732567,
 0.13463335965499817754561346078473,
 0.10954318125799102199776746711408,
 0.07472567457529029657288816982885,
 0.03333567215434406879678440494667
};

template <>
MathVector<1> GaussQuadBase<GaussQuadrature<ReferenceEdge, 19>, 1, 19, 10>::m_vPoint[10] =
{
MathVector<1>(0.01304673574141413996101799395777),
MathVector<1>(0.06746831665550774463395165578825),
MathVector<1>(0.16029521585048779688283631744256),
MathVector<1>(0.28330230293537640460036702841711),
MathVector<1>(0.42556283050918439455758699943514),
MathVector<1>(0.57443716949081560544241300056486),
MathVector<1>(0.71669769706462359539963297158289),
MathVector<1>(0.83970478414951220311716368255744),
MathVector<1>(0.93253168334449225536604834421175),
MathVector<1>(0.98695326425858586003898200604223)
};




template <>
FlexGaussQuadrature<ReferenceEdge>::FlexGaussQuadrature(int order)
{
	switch(order)
	{
	case 0:
	case 1:
		const static GaussQuadrature<ReferenceEdge, 1>& q1 
			= Provider<GaussQuadrature<ReferenceEdge, 1> >::get();

		m_order = q1.order();
		m_numPoints = q1.size();
		m_pvPoint = q1.points();
		m_pvWeight = q1.weights();
		break;

	case 2:
	case 3:
		const static GaussQuadrature<ReferenceEdge, 3>& q3 
			= Provider<GaussQuadrature<ReferenceEdge, 3> >::get();

		m_order = q3.order();
		m_numPoints = q3.size();
		m_pvPoint = q3.points();
		m_pvWeight = q3.weights();
		break;

	case 4:
	case 5:
		const static GaussQuadrature<ReferenceEdge, 5>& q5 
			= Provider<GaussQuadrature<ReferenceEdge, 5> >::get();

		m_order = q5.order();
		m_numPoints = q5.size();
		m_pvPoint = q5.points();
		m_pvWeight = q5.weights();
		break;

	case 6:
	case 7:
		const static GaussQuadrature<ReferenceEdge, 7>& q7 
			= Provider<GaussQuadrature<ReferenceEdge, 7> >::get();

		m_order = q7.order();
		m_numPoints = q7.size();
		m_pvPoint = q7.points();
		m_pvWeight = q7.weights();
		break;

	case 8:
	case 9:
		const static GaussQuadrature<ReferenceEdge, 9>& q9 
			= Provider<GaussQuadrature<ReferenceEdge, 9> >::get();

		m_order = q9.order();
		m_numPoints = q9.size();
		m_pvPoint = q9.points();
		m_pvWeight = q9.weights();
		break;

	case 10:
	case 11:
		const static GaussQuadrature<ReferenceEdge, 11>& q11 
			= Provider<GaussQuadrature<ReferenceEdge, 11> >::get();

		m_order = q11.order();
		m_numPoints = q11.size();
		m_pvPoint = q11.points();
		m_pvWeight = q11.weights();
		break;

	case 12:
	case 13:
		const static GaussQuadrature<ReferenceEdge, 13>& q13 
			= Provider<GaussQuadrature<ReferenceEdge, 13> >::get();

		m_order = q13.order();
		m_numPoints = q13.size();
		m_pvPoint = q13.points();
		m_pvWeight = q13.weights();
		break;

	case 14:
	case 15:
		const static GaussQuadrature<ReferenceEdge, 15>& q15 
			= Provider<GaussQuadrature<ReferenceEdge, 15> >::get();

		m_order = q15.order();
		m_numPoints = q15.size();
		m_pvPoint = q15.points();
		m_pvWeight = q15.weights();
		break;

	case 16:
	case 17:
		const static GaussQuadrature<ReferenceEdge, 17>& q17 
			= Provider<GaussQuadrature<ReferenceEdge, 17> >::get();

		m_order = q17.order();
		m_numPoints = q17.size();
		m_pvPoint = q17.points();
		m_pvWeight = q17.weights();
		break;

	case 18:
	case 19:
		const static GaussQuadrature<ReferenceEdge, 19>& q19 
			= Provider<GaussQuadrature<ReferenceEdge, 19> >::get();

		m_order = q19.order();
		m_numPoints = q19.size();
		m_pvPoint = q19.points();
		m_pvWeight = q19.weights();
		break;

	default: UG_THROW("Order "<<order<<" not available for GaussQuadrature of edge.");
	}
}
}; // namespace ug

