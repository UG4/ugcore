//  This file is parsed from UG 3.9.
//  It provides the Gauss Quadratures for a reference edge.


#include "gauss_quad_edge.h"

namespace ug{

GaussQuadrature<ReferenceEdge, 1>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.500000000000000000;

	m_vWeight[0] =  1.000000000000000000;
}

GaussQuadrature<ReferenceEdge, 3>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.211324865405187118;

	m_vPoint[1][0] = 0.788675134594812882;

	m_vWeight[0] =  0.500000000000000000;
	m_vWeight[1] =  0.500000000000000000;
}

GaussQuadrature<ReferenceEdge, 5>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.112701665379258311;

	m_vPoint[1][0] = 0.500000000000000000;

	m_vPoint[2][0] = 0.887298334620741689;

	m_vWeight[0] =  0.277777777777777778;
	m_vWeight[1] =  0.444444444444444444;
	m_vWeight[2] =  0.277777777777777778;
}

GaussQuadrature<ReferenceEdge, 7>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.069431844202973712;

	m_vPoint[1][0] = 0.330009478207571868;

	m_vPoint[2][0] = 0.669990521792428132;

	m_vPoint[3][0] = 0.930568155797026288;

	m_vWeight[0] =  0.173927422568726929;
	m_vWeight[1] =  0.326072577431273071;
	m_vWeight[2] =  0.326072577431273071;
	m_vWeight[3] =  0.173927422568726929;
}

GaussQuadrature<ReferenceEdge, 9>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.046910077030668004;

	m_vPoint[1][0] = 0.230765344947158454;

	m_vPoint[2][0] = 0.500000000000000000;

	m_vPoint[3][0] = 0.769234655052841546;

	m_vPoint[4][0] = 0.953089922969331996;

	m_vWeight[0] =  0.118463442528094544;
	m_vWeight[1] =  0.239314335249683234;
	m_vWeight[2] =  0.284444444444444444;
	m_vWeight[3] =  0.239314335249683234;
	m_vWeight[4] =  0.118463442528094544;
}

GaussQuadrature<ReferenceEdge, 11>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.033765242898423986;

	m_vPoint[1][0] = 0.169395306766867743;

	m_vPoint[2][0] = 0.380690406958401546;

	m_vPoint[3][0] = 0.619309593041598454;

	m_vPoint[4][0] = 0.830604693233132257;

	m_vPoint[5][0] = 0.966234757101576014;

	m_vWeight[0] =  0.085662246189585173;
	m_vWeight[1] =  0.180380786524069304;
	m_vWeight[2] =  0.233956967286345524;
	m_vWeight[3] =  0.233956967286345524;
	m_vWeight[4] =  0.180380786524069304;
	m_vWeight[5] =  0.085662246189585173;
}

GaussQuadrature<ReferenceEdge, 13>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.025446043828620738;

	m_vPoint[1][0] = 0.129234407200302780;

	m_vPoint[2][0] = 0.297077424311301417;

	m_vPoint[3][0] = 0.500000000000000000;

	m_vPoint[4][0] = 0.702922575688698583;

	m_vPoint[5][0] = 0.870765592799697220;

	m_vPoint[6][0] = 0.974553956171379262;

	m_vWeight[0] =  0.064742483084434847;
	m_vWeight[1] =  0.139852695744638334;
	m_vWeight[2] =  0.190915025252559472;
	m_vWeight[3] =  0.208979591836734694;
	m_vWeight[4] =  0.190915025252559472;
	m_vWeight[5] =  0.139852695744638334;
	m_vWeight[6] =  0.064742483084434847;
}

GaussQuadrature<ReferenceEdge, 15>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.019855071751231884;

	m_vPoint[1][0] = 0.101666761293186630;

	m_vPoint[2][0] = 0.237233795041835507;

	m_vPoint[3][0] = 0.408282678752175098;

	m_vPoint[4][0] = 0.591717321247824902;

	m_vPoint[5][0] = 0.762766204958164493;

	m_vPoint[6][0] = 0.898333238706813370;

	m_vPoint[7][0] = 0.980144928248768116;

	m_vWeight[0] =  0.050614268145188130;
	m_vWeight[1] =  0.111190517226687235;
	m_vWeight[2] =  0.156853322938943644;
	m_vWeight[3] =  0.181341891689180991;
	m_vWeight[4] =  0.181341891689180991;
	m_vWeight[5] =  0.156853322938943644;
	m_vWeight[6] =  0.111190517226687235;
	m_vWeight[7] =  0.050614268145188130;
}

GaussQuadrature<ReferenceEdge, 17>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.015919880246186955;

	m_vPoint[1][0] = 0.081984446336682103;

	m_vPoint[2][0] = 0.193314283649704801;

	m_vPoint[3][0] = 0.337873288298095535;

	m_vPoint[4][0] = 0.500000000000000000;

	m_vPoint[5][0] = 0.662126711701904465;

	m_vPoint[6][0] = 0.806685716350295199;

	m_vPoint[7][0] = 0.918015553663317897;

	m_vPoint[8][0] = 0.984080119753813045;

	m_vWeight[0] =  0.040637194180787206;
	m_vWeight[1] =  0.090324080347428702;
	m_vWeight[2] =  0.130305348201467731;
	m_vWeight[3] =  0.156173538520001420;
	m_vWeight[4] =  0.165119677500629882;
	m_vWeight[5] =  0.156173538520001420;
	m_vWeight[6] =  0.130305348201467731;
	m_vWeight[7] =  0.090324080347428702;
	m_vWeight[8] =  0.040637194180787206;
}

GaussQuadrature<ReferenceEdge, 19>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.013046735741414140;

	m_vPoint[1][0] = 0.067468316655507745;

	m_vPoint[2][0] = 0.160295215850487797;

	m_vPoint[3][0] = 0.283302302935376405;

	m_vPoint[4][0] = 0.425562830509184395;

	m_vPoint[5][0] = 0.574437169490815605;

	m_vPoint[6][0] = 0.716697697064623595;

	m_vPoint[7][0] = 0.839704784149512203;

	m_vPoint[8][0] = 0.932531683344492255;

	m_vPoint[9][0] = 0.986953264258585860;

	m_vWeight[0] =  0.033335672154344069;
	m_vWeight[1] =  0.074725674575290297;
	m_vWeight[2] =  0.109543181257991022;
	m_vWeight[3] =  0.134633359654998178;
	m_vWeight[4] =  0.147762112357376435;
	m_vWeight[5] =  0.147762112357376435;
	m_vWeight[6] =  0.134633359654998178;
	m_vWeight[7] =  0.109543181257991022;
	m_vWeight[8] =  0.074725674575290297;
	m_vWeight[9] =  0.033335672154344069;
}




template <>
FlexGaussQuadrature<ReferenceEdge>::FlexGaussQuadrature(int order)
{
	switch(order)
	{
	case 1:
		const static GaussQuadrature<ReferenceEdge, 1>& q1 
			= QuadRuleProvider::get<GaussQuadrature<ReferenceEdge, 1> >();

		m_order = q1.order();
		m_numPoints = q1.size();
		m_pvPoint = q1.points();
		m_pvWeight = q1.weights();
		break;

	case 3:
		const static GaussQuadrature<ReferenceEdge, 3>& q3 
			= QuadRuleProvider::get<GaussQuadrature<ReferenceEdge, 3> >();

		m_order = q3.order();
		m_numPoints = q3.size();
		m_pvPoint = q3.points();
		m_pvWeight = q3.weights();
		break;

	case 5:
		const static GaussQuadrature<ReferenceEdge, 5>& q5 
			= QuadRuleProvider::get<GaussQuadrature<ReferenceEdge, 5> >();

		m_order = q5.order();
		m_numPoints = q5.size();
		m_pvPoint = q5.points();
		m_pvWeight = q5.weights();
		break;

	case 7:
		const static GaussQuadrature<ReferenceEdge, 7>& q7 
			= QuadRuleProvider::get<GaussQuadrature<ReferenceEdge, 7> >();

		m_order = q7.order();
		m_numPoints = q7.size();
		m_pvPoint = q7.points();
		m_pvWeight = q7.weights();
		break;

	case 9:
		const static GaussQuadrature<ReferenceEdge, 9>& q9 
			= QuadRuleProvider::get<GaussQuadrature<ReferenceEdge, 9> >();

		m_order = q9.order();
		m_numPoints = q9.size();
		m_pvPoint = q9.points();
		m_pvWeight = q9.weights();
		break;

	case 11:
		const static GaussQuadrature<ReferenceEdge, 11>& q11 
			= QuadRuleProvider::get<GaussQuadrature<ReferenceEdge, 11> >();

		m_order = q11.order();
		m_numPoints = q11.size();
		m_pvPoint = q11.points();
		m_pvWeight = q11.weights();
		break;

	case 13:
		const static GaussQuadrature<ReferenceEdge, 13>& q13 
			= QuadRuleProvider::get<GaussQuadrature<ReferenceEdge, 13> >();

		m_order = q13.order();
		m_numPoints = q13.size();
		m_pvPoint = q13.points();
		m_pvWeight = q13.weights();
		break;

	case 15:
		const static GaussQuadrature<ReferenceEdge, 15>& q15 
			= QuadRuleProvider::get<GaussQuadrature<ReferenceEdge, 15> >();

		m_order = q15.order();
		m_numPoints = q15.size();
		m_pvPoint = q15.points();
		m_pvWeight = q15.weights();
		break;

	case 17:
		const static GaussQuadrature<ReferenceEdge, 17>& q17 
			= QuadRuleProvider::get<GaussQuadrature<ReferenceEdge, 17> >();

		m_order = q17.order();
		m_numPoints = q17.size();
		m_pvPoint = q17.points();
		m_pvWeight = q17.weights();
		break;

	case 19:
		const static GaussQuadrature<ReferenceEdge, 19>& q19 
			= QuadRuleProvider::get<GaussQuadrature<ReferenceEdge, 19> >();

		m_order = q19.order();
		m_numPoints = q19.size();
		m_pvPoint = q19.points();
		m_pvWeight = q19.weights();
		break;

	default: UG_ASSERT(0, "Order not availabile. Can not construct GaussQuadrature.\n");
	}
}



// register rules
template <>
bool RegisterQuadratureRule(QuadratureRuleProvider<ReferenceEdge>& factory)
{
	static FlexGaussQuadrature<ReferenceEdge> gaussQuadratureReferenceEdge_1(1);
	static FlexGaussQuadrature<ReferenceEdge> gaussQuadratureReferenceEdge_3(3);
	static FlexGaussQuadrature<ReferenceEdge> gaussQuadratureReferenceEdge_5(5);
	static FlexGaussQuadrature<ReferenceEdge> gaussQuadratureReferenceEdge_7(7);
	static FlexGaussQuadrature<ReferenceEdge> gaussQuadratureReferenceEdge_9(9);
	static FlexGaussQuadrature<ReferenceEdge> gaussQuadratureReferenceEdge_11(11);
	static FlexGaussQuadrature<ReferenceEdge> gaussQuadratureReferenceEdge_13(13);
	static FlexGaussQuadrature<ReferenceEdge> gaussQuadratureReferenceEdge_15(15);
	static FlexGaussQuadrature<ReferenceEdge> gaussQuadratureReferenceEdge_17(17);
	static FlexGaussQuadrature<ReferenceEdge> gaussQuadratureReferenceEdge_19(19);

	bool success = true;
	success &= factory.register_rule(gaussQuadratureReferenceEdge_1);
	success &= factory.register_rule(gaussQuadratureReferenceEdge_3);
	success &= factory.register_rule(gaussQuadratureReferenceEdge_5);
	success &= factory.register_rule(gaussQuadratureReferenceEdge_7);
	success &= factory.register_rule(gaussQuadratureReferenceEdge_9);
	success &= factory.register_rule(gaussQuadratureReferenceEdge_11);
	success &= factory.register_rule(gaussQuadratureReferenceEdge_13);
	success &= factory.register_rule(gaussQuadratureReferenceEdge_15);
	success &= factory.register_rule(gaussQuadratureReferenceEdge_17);
	success &= factory.register_rule(gaussQuadratureReferenceEdge_19);

	return success;
};

}; // namespace ug

