//  This file is parsed from UG 3.9.
//  It provides the Gauss Quadratures for a reference triangle.


#include "../quadrature.h"
#include "gauss_quad_triangle.h"
#include "common/util/provider.h"

namespace ug{

GaussQuadrature<ReferenceTriangle, 1>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.333333333333333333;
	m_vPoint[0][1] = 0.333333333333333333;

	m_vWeight[0] = 0.5 * 1.000000000000000000;
}

GaussQuadrature<ReferenceTriangle, 2>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.166666666666666667;
	m_vPoint[0][1] = 0.166666666666666667;

	m_vPoint[1][0] = 0.166666666666666667;
	m_vPoint[1][1] = 0.666666666666666667;

	m_vPoint[2][0] = 0.666666666666666667;
	m_vPoint[2][1] = 0.166666666666666667;

	m_vWeight[0] = 0.5 * 0.333333333333333333;
	m_vWeight[1] = 0.5 * 0.333333333333333333;
	m_vWeight[2] = 0.5 * 0.333333333333333333;
}

GaussQuadrature<ReferenceTriangle, 3>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.333333333333333333;
	m_vPoint[0][1] = 0.333333333333333333;

	m_vPoint[1][0] = 0.200000000000000000;
	m_vPoint[1][1] = 0.200000000000000000;

	m_vPoint[2][0] = 0.200000000000000000;
	m_vPoint[2][1] = 0.600000000000000000;

	m_vPoint[3][0] = 0.600000000000000000;
	m_vPoint[3][1] = 0.200000000000000000;

	m_vWeight[0] = 0.5 * -0.562500000000000000;
	m_vWeight[1] = 0.5 * 0.520833333333333333;
	m_vWeight[2] = 0.5 * 0.520833333333333333;
	m_vWeight[3] = 0.5 * 0.520833333333333333;
}

GaussQuadrature<ReferenceTriangle, 4>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.091576213509770743;
	m_vPoint[0][1] = 0.091576213509770743;

	m_vPoint[1][0] = 0.091576213509770743;
	m_vPoint[1][1] = 0.816847572980458513;

	m_vPoint[2][0] = 0.816847572980458513;
	m_vPoint[2][1] = 0.091576213509770743;

	m_vPoint[3][0] = 0.445948490915964886;
	m_vPoint[3][1] = 0.445948490915964886;

	m_vPoint[4][0] = 0.445948490915964886;
	m_vPoint[4][1] = 0.108103018168070227;

	m_vPoint[5][0] = 0.108103018168070227;
	m_vPoint[5][1] = 0.445948490915964886;

	m_vWeight[0] = 0.5 * 0.109951743655321868;
	m_vWeight[1] = 0.5 * 0.109951743655321868;
	m_vWeight[2] = 0.5 * 0.109951743655321868;
	m_vWeight[3] = 0.5 * 0.223381589678011466;
	m_vWeight[4] = 0.5 * 0.223381589678011466;
	m_vWeight[5] = 0.5 * 0.223381589678011466;
}

GaussQuadrature<ReferenceTriangle, 5>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.333333333333333333;
	m_vPoint[0][1] = 0.333333333333333333;

	m_vPoint[1][0] = 0.101286507323456339;
	m_vPoint[1][1] = 0.101286507323456339;

	m_vPoint[2][0] = 0.101286507323456339;
	m_vPoint[2][1] = 0.797426985353087322;

	m_vPoint[3][0] = 0.797426985353087322;
	m_vPoint[3][1] = 0.101286507323456339;

	m_vPoint[4][0] = 0.470142064105115090;
	m_vPoint[4][1] = 0.470142064105115090;

	m_vPoint[5][0] = 0.470142064105115090;
	m_vPoint[5][1] = 0.059715871789769820;

	m_vPoint[6][0] = 0.059715871789769820;
	m_vPoint[6][1] = 0.470142064105115090;

	m_vWeight[0] = 0.5 * 0.225000000000000000;
	m_vWeight[1] = 0.5 * 0.125939180544827153;
	m_vWeight[2] = 0.5 * 0.125939180544827153;
	m_vWeight[3] = 0.5 * 0.125939180544827153;
	m_vWeight[4] = 0.5 * 0.132394152788506181;
	m_vWeight[5] = 0.5 * 0.132394152788506181;
	m_vWeight[6] = 0.5 * 0.132394152788506181;
}

GaussQuadrature<ReferenceTriangle, 6>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.063089014491502228;
	m_vPoint[0][1] = 0.063089014491502228;

	m_vPoint[1][0] = 0.063089014491502228;
	m_vPoint[1][1] = 0.873821971016995543;

	m_vPoint[2][0] = 0.873821971016995543;
	m_vPoint[2][1] = 0.063089014491502228;

	m_vPoint[3][0] = 0.249286745170910421;
	m_vPoint[3][1] = 0.249286745170910421;

	m_vPoint[4][0] = 0.249286745170910421;
	m_vPoint[4][1] = 0.501426509658179157;

	m_vPoint[5][0] = 0.501426509658179157;
	m_vPoint[5][1] = 0.249286745170910421;

	m_vPoint[6][0] = 0.053145049844816947;
	m_vPoint[6][1] = 0.310352451033784405;

	m_vPoint[7][0] = 0.053145049844816947;
	m_vPoint[7][1] = 0.636502499121398647;

	m_vPoint[8][0] = 0.310352451033784405;
	m_vPoint[8][1] = 0.053145049844816947;

	m_vPoint[9][0] = 0.310352451033784405;
	m_vPoint[9][1] = 0.636502499121398647;

	m_vPoint[10][0] = 0.636502499121398647;
	m_vPoint[10][1] = 0.053145049844816947;

	m_vPoint[11][0] = 0.636502499121398647;
	m_vPoint[11][1] = 0.310352451033784405;

	m_vWeight[0] = 0.5 * 0.050844906370206817;
	m_vWeight[1] = 0.5 * 0.050844906370206817;
	m_vWeight[2] = 0.5 * 0.050844906370206817;
	m_vWeight[3] = 0.5 * 0.116786275726379366;
	m_vWeight[4] = 0.5 * 0.116786275726379366;
	m_vWeight[5] = 0.5 * 0.116786275726379366;
	m_vWeight[6] = 0.5 * 0.082851075618373575;
	m_vWeight[7] = 0.5 * 0.082851075618373575;
	m_vWeight[8] = 0.5 * 0.082851075618373575;
	m_vWeight[9] = 0.5 * 0.082851075618373575;
	m_vWeight[10] = 0.5 * 0.082851075618373575;
	m_vWeight[11] = 0.5 * 0.082851075618373575;
}

GaussQuadrature<ReferenceTriangle, 7>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.062382265094402118;
	m_vPoint[0][1] = 0.067517867073916085;

	m_vPoint[1][0] = 0.067517867073916085;
	m_vPoint[1][1] = 0.870099867831681796;

	m_vPoint[2][0] = 0.870099867831681796;
	m_vPoint[2][1] = 0.062382265094402118;

	m_vPoint[3][0] = 0.055225456656926612;
	m_vPoint[3][1] = 0.321502493851981823;

	m_vPoint[4][0] = 0.321502493851981823;
	m_vPoint[4][1] = 0.623272049491091566;

	m_vPoint[5][0] = 0.623272049491091566;
	m_vPoint[5][1] = 0.055225456656926612;

	m_vPoint[6][0] = 0.034324302945097146;
	m_vPoint[6][1] = 0.660949196186735658;

	m_vPoint[7][0] = 0.660949196186735658;
	m_vPoint[7][1] = 0.304726500868167196;

	m_vPoint[8][0] = 0.304726500868167196;
	m_vPoint[8][1] = 0.034324302945097146;

	m_vPoint[9][0] = 0.515842334353591779;
	m_vPoint[9][1] = 0.277716166976391783;

	m_vPoint[10][0] = 0.277716166976391783;
	m_vPoint[10][1] = 0.206441498670016438;

	m_vPoint[11][0] = 0.206441498670016438;
	m_vPoint[11][1] = 0.515842334353591779;

	m_vWeight[0] = 0.5 * 0.053034056314872503;
	m_vWeight[1] = 0.5 * 0.053034056314872503;
	m_vWeight[2] = 0.5 * 0.053034056314872503;
	m_vWeight[3] = 0.5 * 0.087762817428892110;
	m_vWeight[4] = 0.5 * 0.087762817428892110;
	m_vWeight[5] = 0.5 * 0.087762817428892110;
	m_vWeight[6] = 0.5 * 0.057550085569963171;
	m_vWeight[7] = 0.5 * 0.057550085569963171;
	m_vWeight[8] = 0.5 * 0.057550085569963171;
	m_vWeight[9] = 0.5 * 0.134986374019605549;
	m_vWeight[10] = 0.5 * 0.134986374019605549;
	m_vWeight[11] = 0.5 * 0.134986374019605549;
}

GaussQuadrature<ReferenceTriangle, 8>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.333333333333333333;
	m_vPoint[0][1] = 0.333333333333333333;

	m_vPoint[1][0] = 0.170569307751760207;
	m_vPoint[1][1] = 0.170569307751760207;

	m_vPoint[2][0] = 0.170569307751760207;
	m_vPoint[2][1] = 0.658861384496479587;

	m_vPoint[3][0] = 0.658861384496479587;
	m_vPoint[3][1] = 0.170569307751760207;

	m_vPoint[4][0] = 0.050547228317030975;
	m_vPoint[4][1] = 0.050547228317030975;

	m_vPoint[5][0] = 0.050547228317030975;
	m_vPoint[5][1] = 0.898905543365938049;

	m_vPoint[6][0] = 0.898905543365938049;
	m_vPoint[6][1] = 0.050547228317030975;

	m_vPoint[7][0] = 0.459292588292723156;
	m_vPoint[7][1] = 0.459292588292723156;

	m_vPoint[8][0] = 0.459292588292723156;
	m_vPoint[8][1] = 0.081414823414553688;

	m_vPoint[9][0] = 0.081414823414553688;
	m_vPoint[9][1] = 0.459292588292723156;

	m_vPoint[10][0] = 0.728492392955404281;
	m_vPoint[10][1] = 0.263112829634638113;

	m_vPoint[11][0] = 0.728492392955404281;
	m_vPoint[11][1] = 0.008394777409957605;

	m_vPoint[12][0] = 0.263112829634638113;
	m_vPoint[12][1] = 0.728492392955404281;

	m_vPoint[13][0] = 0.263112829634638113;
	m_vPoint[13][1] = 0.008394777409957605;

	m_vPoint[14][0] = 0.008394777409957605;
	m_vPoint[14][1] = 0.728492392955404281;

	m_vPoint[15][0] = 0.008394777409957605;
	m_vPoint[15][1] = 0.263112829634638113;

	m_vWeight[0] = 0.5 * 0.144315607677787168;
	m_vWeight[1] = 0.5 * 0.103217370534718250;
	m_vWeight[2] = 0.5 * 0.103217370534718250;
	m_vWeight[3] = 0.5 * 0.103217370534718250;
	m_vWeight[4] = 0.5 * 0.032458497623198080;
	m_vWeight[5] = 0.5 * 0.032458497623198080;
	m_vWeight[6] = 0.5 * 0.032458497623198080;
	m_vWeight[7] = 0.5 * 0.095091634267284625;
	m_vWeight[8] = 0.5 * 0.095091634267284625;
	m_vWeight[9] = 0.5 * 0.095091634267284625;
	m_vWeight[10] = 0.5 * 0.027230314174434994;
	m_vWeight[11] = 0.5 * 0.027230314174434994;
	m_vWeight[12] = 0.5 * 0.027230314174434994;
	m_vWeight[13] = 0.5 * 0.027230314174434994;
	m_vWeight[14] = 0.5 * 0.027230314174434994;
	m_vWeight[15] = 0.5 * 0.027230314174434994;
}

GaussQuadrature<ReferenceTriangle, 9>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.333333333333333333;
	m_vPoint[0][1] = 0.333333333333333333;

	m_vPoint[1][0] = 0.489682519198737628;
	m_vPoint[1][1] = 0.489682519198737628;

	m_vPoint[2][0] = 0.489682519198737628;
	m_vPoint[2][1] = 0.020634961602524744;

	m_vPoint[3][0] = 0.020634961602524744;
	m_vPoint[3][1] = 0.489682519198737628;

	m_vPoint[4][0] = 0.437089591492936637;
	m_vPoint[4][1] = 0.437089591492936637;

	m_vPoint[5][0] = 0.437089591492936637;
	m_vPoint[5][1] = 0.125820817014126725;

	m_vPoint[6][0] = 0.125820817014126725;
	m_vPoint[6][1] = 0.437089591492936637;

	m_vPoint[7][0] = 0.188203535619032730;
	m_vPoint[7][1] = 0.188203535619032730;

	m_vPoint[8][0] = 0.188203535619032730;
	m_vPoint[8][1] = 0.623592928761934540;

	m_vPoint[9][0] = 0.623592928761934540;
	m_vPoint[9][1] = 0.188203535619032730;

	m_vPoint[10][0] = 0.044729513394452710;
	m_vPoint[10][1] = 0.044729513394452710;

	m_vPoint[11][0] = 0.044729513394452710;
	m_vPoint[11][1] = 0.910540973211094580;

	m_vPoint[12][0] = 0.910540973211094580;
	m_vPoint[12][1] = 0.044729513394452710;

	m_vPoint[13][0] = 0.741198598784498021;
	m_vPoint[13][1] = 0.036838412054736284;

	m_vPoint[14][0] = 0.741198598784498021;
	m_vPoint[14][1] = 0.221962989160765696;

	m_vPoint[15][0] = 0.036838412054736284;
	m_vPoint[15][1] = 0.741198598784498021;

	m_vPoint[16][0] = 0.036838412054736284;
	m_vPoint[16][1] = 0.221962989160765696;

	m_vPoint[17][0] = 0.221962989160765696;
	m_vPoint[17][1] = 0.741198598784498021;

	m_vPoint[18][0] = 0.221962989160765696;
	m_vPoint[18][1] = 0.036838412054736284;

	m_vWeight[0] = 0.5 * 0.097135796282798834;
	m_vWeight[1] = 0.5 * 0.031334700227139071;
	m_vWeight[2] = 0.5 * 0.031334700227139071;
	m_vWeight[3] = 0.5 * 0.031334700227139071;
	m_vWeight[4] = 0.5 * 0.077827541004774279;
	m_vWeight[5] = 0.5 * 0.077827541004774279;
	m_vWeight[6] = 0.5 * 0.077827541004774279;
	m_vWeight[7] = 0.5 * 0.079647738927210253;
	m_vWeight[8] = 0.5 * 0.079647738927210253;
	m_vWeight[9] = 0.5 * 0.079647738927210253;
	m_vWeight[10] = 0.5 * 0.025577675658698031;
	m_vWeight[11] = 0.5 * 0.025577675658698031;
	m_vWeight[12] = 0.5 * 0.025577675658698031;
	m_vWeight[13] = 0.5 * 0.043283539377289377;
	m_vWeight[14] = 0.5 * 0.043283539377289377;
	m_vWeight[15] = 0.5 * 0.043283539377289377;
	m_vWeight[16] = 0.5 * 0.043283539377289377;
	m_vWeight[17] = 0.5 * 0.043283539377289377;
	m_vWeight[18] = 0.5 * 0.043283539377289377;
}

GaussQuadrature<ReferenceTriangle, 10>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.333333333333333333;
	m_vPoint[0][1] = 0.333333333333333333;

	m_vPoint[1][0] = 0.425086210602090573;
	m_vPoint[1][1] = 0.425086210602090573;

	m_vPoint[2][0] = 0.425086210602090573;
	m_vPoint[2][1] = 0.149827578795818854;

	m_vPoint[3][0] = 0.149827578795818854;
	m_vPoint[3][1] = 0.425086210602090573;

	m_vPoint[4][0] = 0.023308867510000191;
	m_vPoint[4][1] = 0.023308867510000191;

	m_vPoint[5][0] = 0.023308867510000191;
	m_vPoint[5][1] = 0.953382264979999619;

	m_vPoint[6][0] = 0.953382264979999619;
	m_vPoint[6][1] = 0.023308867510000191;

	m_vPoint[7][0] = 0.628307400213492556;
	m_vPoint[7][1] = 0.223766973576973006;

	m_vPoint[8][0] = 0.628307400213492556;
	m_vPoint[8][1] = 0.147925626209534437;

	m_vPoint[9][0] = 0.223766973576973006;
	m_vPoint[9][1] = 0.628307400213492556;

	m_vPoint[10][0] = 0.223766973576973006;
	m_vPoint[10][1] = 0.147925626209534437;

	m_vPoint[11][0] = 0.147925626209534437;
	m_vPoint[11][1] = 0.628307400213492556;

	m_vPoint[12][0] = 0.147925626209534437;
	m_vPoint[12][1] = 0.223766973576973006;

	m_vPoint[13][0] = 0.611313826181397649;
	m_vPoint[13][1] = 0.358740141864431465;

	m_vPoint[14][0] = 0.611313826181397649;
	m_vPoint[14][1] = 0.029946031954170887;

	m_vPoint[15][0] = 0.358740141864431465;
	m_vPoint[15][1] = 0.611313826181397649;

	m_vPoint[16][0] = 0.358740141864431465;
	m_vPoint[16][1] = 0.029946031954170887;

	m_vPoint[17][0] = 0.029946031954170887;
	m_vPoint[17][1] = 0.611313826181397649;

	m_vPoint[18][0] = 0.029946031954170887;
	m_vPoint[18][1] = 0.358740141864431465;

	m_vPoint[19][0] = 0.821072069985629373;
	m_vPoint[19][1] = 0.143295370426867145;

	m_vPoint[20][0] = 0.821072069985629373;
	m_vPoint[20][1] = 0.035632559587503481;

	m_vPoint[21][0] = 0.143295370426867145;
	m_vPoint[21][1] = 0.821072069985629373;

	m_vPoint[22][0] = 0.143295370426867145;
	m_vPoint[22][1] = 0.035632559587503481;

	m_vPoint[23][0] = 0.035632559587503481;
	m_vPoint[23][1] = 0.821072069985629373;

	m_vPoint[24][0] = 0.035632559587503481;
	m_vPoint[24][1] = 0.143295370426867145;

	m_vWeight[0] = 0.5 * 0.079894504741239708;
	m_vWeight[1] = 0.5 * 0.071123802232377335;
	m_vWeight[2] = 0.5 * 0.071123802232377335;
	m_vWeight[3] = 0.5 * 0.071123802232377335;
	m_vWeight[4] = 0.5 * 0.008223818690464196;
	m_vWeight[5] = 0.5 * 0.008223818690464196;
	m_vWeight[6] = 0.5 * 0.008223818690464196;
	m_vWeight[7] = 0.5 * 0.045430592296170018;
	m_vWeight[8] = 0.5 * 0.045430592296170018;
	m_vWeight[9] = 0.5 * 0.045430592296170018;
	m_vWeight[10] = 0.5 * 0.045430592296170018;
	m_vWeight[11] = 0.5 * 0.045430592296170018;
	m_vWeight[12] = 0.5 * 0.045430592296170018;
	m_vWeight[13] = 0.5 * 0.037359856234305277;
	m_vWeight[14] = 0.5 * 0.037359856234305277;
	m_vWeight[15] = 0.5 * 0.037359856234305277;
	m_vWeight[16] = 0.5 * 0.037359856234305277;
	m_vWeight[17] = 0.5 * 0.037359856234305277;
	m_vWeight[18] = 0.5 * 0.037359856234305277;
	m_vWeight[19] = 0.5 * 0.030886656884563989;
	m_vWeight[20] = 0.5 * 0.030886656884563989;
	m_vWeight[21] = 0.5 * 0.030886656884563989;
	m_vWeight[22] = 0.5 * 0.030886656884563989;
	m_vWeight[23] = 0.5 * 0.030886656884563989;
	m_vWeight[24] = 0.5 * 0.030886656884563989;
}

GaussQuadrature<ReferenceTriangle, 11>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.858870281282636704;
	m_vPoint[0][1] = 0.141129718717363296;

	m_vPoint[1][0] = 0.858870281282636704;
	m_vPoint[1][1] = 0.000000000000000000;

	m_vPoint[2][0] = 0.141129718717363296;
	m_vPoint[2][1] = 0.858870281282636704;

	m_vPoint[3][0] = 0.141129718717363296;
	m_vPoint[3][1] = 0.000000000000000000;

	m_vPoint[4][0] = 0.000000000000000000;
	m_vPoint[4][1] = 0.858870281282636704;

	m_vPoint[5][0] = 0.000000000000000000;
	m_vPoint[5][1] = 0.141129718717363296;

	m_vPoint[6][0] = 0.333333333333333333;
	m_vPoint[6][1] = 0.333333333333333333;

	m_vPoint[7][0] = 0.025989140928287395;
	m_vPoint[7][1] = 0.025989140928287395;

	m_vPoint[8][0] = 0.025989140928287395;
	m_vPoint[8][1] = 0.948021718143425209;

	m_vPoint[9][0] = 0.948021718143425209;
	m_vPoint[9][1] = 0.025989140928287395;

	m_vPoint[10][0] = 0.094287502647922496;
	m_vPoint[10][1] = 0.094287502647922496;

	m_vPoint[11][0] = 0.094287502647922496;
	m_vPoint[11][1] = 0.811424994704155009;

	m_vPoint[12][0] = 0.811424994704155009;
	m_vPoint[12][1] = 0.094287502647922496;

	m_vPoint[13][0] = 0.494636775017213814;
	m_vPoint[13][1] = 0.494636775017213814;

	m_vPoint[14][0] = 0.494636775017213814;
	m_vPoint[14][1] = 0.010726449965572373;

	m_vPoint[15][0] = 0.010726449965572373;
	m_vPoint[15][1] = 0.494636775017213814;

	m_vPoint[16][0] = 0.207343382614511333;
	m_vPoint[16][1] = 0.207343382614511333;

	m_vPoint[17][0] = 0.207343382614511333;
	m_vPoint[17][1] = 0.585313234770977333;

	m_vPoint[18][0] = 0.585313234770977333;
	m_vPoint[18][1] = 0.207343382614511333;

	m_vPoint[19][0] = 0.438907805700492095;
	m_vPoint[19][1] = 0.438907805700492095;

	m_vPoint[20][0] = 0.438907805700492095;
	m_vPoint[20][1] = 0.122184388599015810;

	m_vPoint[21][0] = 0.122184388599015810;
	m_vPoint[21][1] = 0.438907805700492095;

	m_vPoint[22][0] = 0.677937654882590402;
	m_vPoint[22][1] = 0.044841677589130443;

	m_vPoint[23][0] = 0.677937654882590402;
	m_vPoint[23][1] = 0.277220667528279155;

	m_vPoint[24][0] = 0.044841677589130443;
	m_vPoint[24][1] = 0.677937654882590402;

	m_vPoint[25][0] = 0.044841677589130443;
	m_vPoint[25][1] = 0.277220667528279155;

	m_vPoint[26][0] = 0.277220667528279155;
	m_vPoint[26][1] = 0.677937654882590402;

	m_vPoint[27][0] = 0.277220667528279155;
	m_vPoint[27][1] = 0.044841677589130443;

	m_vWeight[0] = 0.5 * 0.007362383783300554;
	m_vWeight[1] = 0.5 * 0.007362383783300554;
	m_vWeight[2] = 0.5 * 0.007362383783300554;
	m_vWeight[3] = 0.5 * 0.007362383783300554;
	m_vWeight[4] = 0.5 * 0.007362383783300554;
	m_vWeight[5] = 0.5 * 0.007362383783300554;
	m_vWeight[6] = 0.5 * 0.087977301162232239;
	m_vWeight[7] = 0.5 * 0.008744311553736023;
	m_vWeight[8] = 0.5 * 0.008744311553736023;
	m_vWeight[9] = 0.5 * 0.008744311553736023;
	m_vWeight[10] = 0.5 * 0.038081571993934938;
	m_vWeight[11] = 0.5 * 0.038081571993934938;
	m_vWeight[12] = 0.5 * 0.038081571993934938;
	m_vWeight[13] = 0.5 * 0.018855448056131292;
	m_vWeight[14] = 0.5 * 0.018855448056131292;
	m_vWeight[15] = 0.5 * 0.018855448056131292;
	m_vWeight[16] = 0.5 * 0.072159697544739526;
	m_vWeight[17] = 0.5 * 0.072159697544739526;
	m_vWeight[18] = 0.5 * 0.072159697544739526;
	m_vWeight[19] = 0.5 * 0.069329138705535900;
	m_vWeight[20] = 0.5 * 0.069329138705535900;
	m_vWeight[21] = 0.5 * 0.069329138705535900;
	m_vWeight[22] = 0.5 * 0.041056315429288567;
	m_vWeight[23] = 0.5 * 0.041056315429288567;
	m_vWeight[24] = 0.5 * 0.041056315429288567;
	m_vWeight[25] = 0.5 * 0.041056315429288567;
	m_vWeight[26] = 0.5 * 0.041056315429288567;
	m_vWeight[27] = 0.5 * 0.041056315429288567;
}

GaussQuadrature<ReferenceTriangle, 12>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.023565220452390000;
	m_vPoint[0][1] = 0.488217389773805000;

	m_vPoint[1][0] = 0.488217389773805000;
	m_vPoint[1][1] = 0.023565220452390000;

	m_vPoint[2][0] = 0.488217389773805000;
	m_vPoint[2][1] = 0.488217389773805000;

	m_vPoint[3][0] = 0.439724392294460000;
	m_vPoint[3][1] = 0.439724392294460000;

	m_vPoint[4][0] = 0.439724392294460000;
	m_vPoint[4][1] = 0.120551215411079000;

	m_vPoint[5][0] = 0.120551215411079000;
	m_vPoint[5][1] = 0.439724392294460000;

	m_vPoint[6][0] = 0.271210385012116000;
	m_vPoint[6][1] = 0.271210385012116000;

	m_vPoint[7][0] = 0.271210385012116000;
	m_vPoint[7][1] = 0.457579229975768000;

	m_vPoint[8][0] = 0.457579229975768000;
	m_vPoint[8][1] = 0.271210385012116000;

	m_vPoint[9][0] = 0.127576145541586000;
	m_vPoint[9][1] = 0.127576145541586000;

	m_vPoint[10][0] = 0.127576145541586000;
	m_vPoint[10][1] = 0.744847708916827900;

	m_vPoint[11][0] = 0.744847708916827900;
	m_vPoint[11][1] = 0.127576145541586000;

	m_vPoint[12][0] = 0.021317350453210000;
	m_vPoint[12][1] = 0.021317350453210000;

	m_vPoint[13][0] = 0.021317350453210000;
	m_vPoint[13][1] = 0.957365299093579900;

	m_vPoint[14][0] = 0.957365299093579900;
	m_vPoint[14][1] = 0.021317350453210000;

	m_vPoint[15][0] = 0.115343494534698000;
	m_vPoint[15][1] = 0.275713269685514000;

	m_vPoint[16][0] = 0.115343494534698000;
	m_vPoint[16][1] = 0.608943235779787900;

	m_vPoint[17][0] = 0.275713269685514000;
	m_vPoint[17][1] = 0.115343494534698000;

	m_vPoint[18][0] = 0.275713269685514000;
	m_vPoint[18][1] = 0.608943235779787900;

	m_vPoint[19][0] = 0.608943235779787900;
	m_vPoint[19][1] = 0.115343494534698000;

	m_vPoint[20][0] = 0.608943235779787900;
	m_vPoint[20][1] = 0.275713269685514000;

	m_vPoint[21][0] = 0.022838332222257000;
	m_vPoint[21][1] = 0.281325580989940000;

	m_vPoint[22][0] = 0.022838332222257000;
	m_vPoint[22][1] = 0.695836086787803100;

	m_vPoint[23][0] = 0.281325580989940000;
	m_vPoint[23][1] = 0.022838332222257000;

	m_vPoint[24][0] = 0.281325580989940000;
	m_vPoint[24][1] = 0.695836086787803100;

	m_vPoint[25][0] = 0.695836086787803100;
	m_vPoint[25][1] = 0.022838332222257000;

	m_vPoint[26][0] = 0.695836086787803100;
	m_vPoint[26][1] = 0.281325580989940000;

	m_vPoint[27][0] = 0.025734050548330000;
	m_vPoint[27][1] = 0.116251915907597000;

	m_vPoint[28][0] = 0.025734050548330000;
	m_vPoint[28][1] = 0.858014033544073000;

	m_vPoint[29][0] = 0.116251915907597000;
	m_vPoint[29][1] = 0.025734050548330000;

	m_vPoint[30][0] = 0.116251915907597000;
	m_vPoint[30][1] = 0.858014033544073000;

	m_vPoint[31][0] = 0.858014033544073000;
	m_vPoint[31][1] = 0.025734050548330000;

	m_vPoint[32][0] = 0.858014033544073000;
	m_vPoint[32][1] = 0.116251915907597000;

	m_vWeight[0] = 0.5 * 0.025731066440455000;
	m_vWeight[1] = 0.5 * 0.025731066440455000;
	m_vWeight[2] = 0.5 * 0.025731066440455000;
	m_vWeight[3] = 0.5 * 0.043692544538038000;
	m_vWeight[4] = 0.5 * 0.043692544538038000;
	m_vWeight[5] = 0.5 * 0.043692544538038000;
	m_vWeight[6] = 0.5 * 0.062858224217885000;
	m_vWeight[7] = 0.5 * 0.062858224217885000;
	m_vWeight[8] = 0.5 * 0.062858224217885000;
	m_vWeight[9] = 0.5 * 0.034796112930709000;
	m_vWeight[10] = 0.5 * 0.034796112930709000;
	m_vWeight[11] = 0.5 * 0.034796112930709000;
	m_vWeight[12] = 0.5 * 0.006166261051559000;
	m_vWeight[13] = 0.5 * 0.006166261051559000;
	m_vWeight[14] = 0.5 * 0.006166261051559000;
	m_vWeight[15] = 0.5 * 0.040371557766381000;
	m_vWeight[16] = 0.5 * 0.040371557766381000;
	m_vWeight[17] = 0.5 * 0.040371557766381000;
	m_vWeight[18] = 0.5 * 0.040371557766381000;
	m_vWeight[19] = 0.5 * 0.040371557766381000;
	m_vWeight[20] = 0.5 * 0.040371557766381000;
	m_vWeight[21] = 0.5 * 0.022356773202303000;
	m_vWeight[22] = 0.5 * 0.022356773202303000;
	m_vWeight[23] = 0.5 * 0.022356773202303000;
	m_vWeight[24] = 0.5 * 0.022356773202303000;
	m_vWeight[25] = 0.5 * 0.022356773202303000;
	m_vWeight[26] = 0.5 * 0.022356773202303000;
	m_vWeight[27] = 0.5 * 0.017316231108659000;
	m_vWeight[28] = 0.5 * 0.017316231108659000;
	m_vWeight[29] = 0.5 * 0.017316231108659000;
	m_vWeight[30] = 0.5 * 0.017316231108659000;
	m_vWeight[31] = 0.5 * 0.017316231108659000;
	m_vWeight[32] = 0.5 * 0.017316231108659000;
}




template <>
FlexGaussQuadrature<ReferenceTriangle>::FlexGaussQuadrature(int order)
{
	switch(order)
	{
	case 1:{
		const static GaussQuadrature<ReferenceTriangle, 1>& q1 
			= Provider<GaussQuadrature<ReferenceTriangle, 1> >::get();

		m_order = q1.order();
		m_numPoints = q1.size();
		m_pvPoint = q1.points();
		m_pvWeight = q1.weights();
		}break;

	case 2:{
		const static GaussQuadrature<ReferenceTriangle, 2>& q2 
			= Provider<GaussQuadrature<ReferenceTriangle, 2> >::get();

		m_order = q2.order();
		m_numPoints = q2.size();
		m_pvPoint = q2.points();
		m_pvWeight = q2.weights();
		}break;

	case 3:{
		const static GaussQuadrature<ReferenceTriangle, 3>& q3 
			= Provider<GaussQuadrature<ReferenceTriangle, 3> >::get();

		m_order = q3.order();
		m_numPoints = q3.size();
		m_pvPoint = q3.points();
		m_pvWeight = q3.weights();
		}break;

	case 4:{
		const static GaussQuadrature<ReferenceTriangle, 4>& q4 
			= Provider<GaussQuadrature<ReferenceTriangle, 4> >::get();

		m_order = q4.order();
		m_numPoints = q4.size();
		m_pvPoint = q4.points();
		m_pvWeight = q4.weights();
		}break;

	case 5:{
		const static GaussQuadrature<ReferenceTriangle, 5>& q5 
			= Provider<GaussQuadrature<ReferenceTriangle, 5> >::get();

		m_order = q5.order();
		m_numPoints = q5.size();
		m_pvPoint = q5.points();
		m_pvWeight = q5.weights();
		}break;

	case 6:{
		const static GaussQuadrature<ReferenceTriangle, 6>& q6 
			= Provider<GaussQuadrature<ReferenceTriangle, 6> >::get();

		m_order = q6.order();
		m_numPoints = q6.size();
		m_pvPoint = q6.points();
		m_pvWeight = q6.weights();
		}break;

	case 7:{
		const static GaussQuadrature<ReferenceTriangle, 7>& q7 
			= Provider<GaussQuadrature<ReferenceTriangle, 7> >::get();

		m_order = q7.order();
		m_numPoints = q7.size();
		m_pvPoint = q7.points();
		m_pvWeight = q7.weights();
		}break;

	case 8:{
		const static GaussQuadrature<ReferenceTriangle, 8>& q8 
			= Provider<GaussQuadrature<ReferenceTriangle, 8> >::get();

		m_order = q8.order();
		m_numPoints = q8.size();
		m_pvPoint = q8.points();
		m_pvWeight = q8.weights();
		}break;

	case 9:{
		const static GaussQuadrature<ReferenceTriangle, 9>& q9 
			= Provider<GaussQuadrature<ReferenceTriangle, 9> >::get();

		m_order = q9.order();
		m_numPoints = q9.size();
		m_pvPoint = q9.points();
		m_pvWeight = q9.weights();
		}break;

	case 10:{
		const static GaussQuadrature<ReferenceTriangle, 10>& q10 
			= Provider<GaussQuadrature<ReferenceTriangle, 10> >::get();

		m_order = q10.order();
		m_numPoints = q10.size();
		m_pvPoint = q10.points();
		m_pvWeight = q10.weights();
		}break;

	case 11:{
		const static GaussQuadrature<ReferenceTriangle, 11>& q11 
			= Provider<GaussQuadrature<ReferenceTriangle, 11> >::get();

		m_order = q11.order();
		m_numPoints = q11.size();
		m_pvPoint = q11.points();
		m_pvWeight = q11.weights();
		}break;

	case 12:{
		const static GaussQuadrature<ReferenceTriangle, 12>& q12 
			= Provider<GaussQuadrature<ReferenceTriangle, 12> >::get();

		m_order = q12.order();
		m_numPoints = q12.size();
		m_pvPoint = q12.points();
		m_pvWeight = q12.weights();
		}break;

	default: UG_ASSERT(0, "Order not availabile. Can not construct GaussQuadrature.\n");
	}
}



// register rules
template <>
bool RegisterGaussQuadRule<ReferenceTriangle>(QuadratureRuleProvider<ReferenceTriangle::dim>& factory)
{
	static FlexGaussQuadrature<ReferenceTriangle> gaussQuadratureReferenceTriangle_1(1);
	static FlexGaussQuadrature<ReferenceTriangle> gaussQuadratureReferenceTriangle_2(2);
	static FlexGaussQuadrature<ReferenceTriangle> gaussQuadratureReferenceTriangle_3(3);
	static FlexGaussQuadrature<ReferenceTriangle> gaussQuadratureReferenceTriangle_4(4);
	static FlexGaussQuadrature<ReferenceTriangle> gaussQuadratureReferenceTriangle_5(5);
	static FlexGaussQuadrature<ReferenceTriangle> gaussQuadratureReferenceTriangle_6(6);
	static FlexGaussQuadrature<ReferenceTriangle> gaussQuadratureReferenceTriangle_7(7);
	static FlexGaussQuadrature<ReferenceTriangle> gaussQuadratureReferenceTriangle_8(8);
	static FlexGaussQuadrature<ReferenceTriangle> gaussQuadratureReferenceTriangle_9(9);
	static FlexGaussQuadrature<ReferenceTriangle> gaussQuadratureReferenceTriangle_10(10);
	static FlexGaussQuadrature<ReferenceTriangle> gaussQuadratureReferenceTriangle_11(11);
	static FlexGaussQuadrature<ReferenceTriangle> gaussQuadratureReferenceTriangle_12(12);

	bool success = true;
	factory.register_rule<ReferenceTriangle>(gaussQuadratureReferenceTriangle_1);
	factory.register_rule<ReferenceTriangle>(gaussQuadratureReferenceTriangle_2);
	factory.register_rule<ReferenceTriangle>(gaussQuadratureReferenceTriangle_3);
	factory.register_rule<ReferenceTriangle>(gaussQuadratureReferenceTriangle_4);
	factory.register_rule<ReferenceTriangle>(gaussQuadratureReferenceTriangle_5);
	factory.register_rule<ReferenceTriangle>(gaussQuadratureReferenceTriangle_6);
	factory.register_rule<ReferenceTriangle>(gaussQuadratureReferenceTriangle_7);
	factory.register_rule<ReferenceTriangle>(gaussQuadratureReferenceTriangle_8);
	factory.register_rule<ReferenceTriangle>(gaussQuadratureReferenceTriangle_9);
	factory.register_rule<ReferenceTriangle>(gaussQuadratureReferenceTriangle_10);
	factory.register_rule<ReferenceTriangle>(gaussQuadratureReferenceTriangle_11);
	factory.register_rule<ReferenceTriangle>(gaussQuadratureReferenceTriangle_12);

	return success;
};

}; // namespace ug

