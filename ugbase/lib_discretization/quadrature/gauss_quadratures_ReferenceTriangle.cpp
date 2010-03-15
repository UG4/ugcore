//This file is parsed from UG 3.9.


#include "quadrature.h"

namespace ug{

template <>
GaussQuadrature<ReferenceTriangle>::GaussQuadrature(int order)
{
	switch(order)
	{
	case 1:
		m_order = 1;
		m_num_points = 1;

		if(allocate_memory(m_num_points) != true)
		{
			assert(0 && "Error while constructing Quadrature Rule.");
		}

		m_points[0][0] = 0.333333333333333333;
		m_points[0][1] = 0.333333333333333333;

		m_weights[0] = 0.5 *  1.000000000000000000;
		break;

	case 2:
		m_order = 2;
		m_num_points = 3;

		if(allocate_memory(m_num_points) != true)
		{
			assert(0 && "Error while constructing Quadrature Rule.");
		}

		m_points[0][0] = 0.166666666666666667;
		m_points[0][1] = 0.166666666666666667;
		m_points[1][0] = 0.166666666666666667;
		m_points[1][1] = 0.666666666666666667;
		m_points[2][0] = 0.666666666666666667;
		m_points[2][1] = 0.166666666666666667;

		m_weights[0] = 0.5 *  0.333333333333333333;
		m_weights[1] = 0.5 *  0.333333333333333333;
		m_weights[2] = 0.5 *  0.333333333333333333;
		break;

	case 3:
		m_order = 3;
		m_num_points = 4;

		if(allocate_memory(m_num_points) != true)
		{
			assert(0 && "Error while constructing Quadrature Rule.");
		}

		m_points[0][0] = 0.333333333333333333;
		m_points[0][1] = 0.333333333333333333;
		m_points[1][0] = 0.200000000000000000;
		m_points[1][1] = 0.200000000000000000;
		m_points[2][0] = 0.200000000000000000;
		m_points[2][1] = 0.600000000000000000;
		m_points[3][0] = 0.600000000000000000;
		m_points[3][1] = 0.200000000000000000;

		m_weights[0] = 0.5 *  -0.562500000000000000;
		m_weights[1] = 0.5 *  0.520833333333333333;
		m_weights[2] = 0.5 *  0.520833333333333333;
		m_weights[3] = 0.5 *  0.520833333333333333;
		break;

	case 4:
		m_order = 4;
		m_num_points = 6;

		if(allocate_memory(m_num_points) != true)
		{
			assert(0 && "Error while constructing Quadrature Rule.");
		}

		m_points[0][0] = 0.091576213509770743;
		m_points[0][1] = 0.091576213509770743;
		m_points[1][0] = 0.091576213509770743;
		m_points[1][1] = 0.816847572980458513;
		m_points[2][0] = 0.816847572980458513;
		m_points[2][1] = 0.091576213509770743;
		m_points[3][0] = 0.445948490915964886;
		m_points[3][1] = 0.445948490915964886;
		m_points[4][0] = 0.445948490915964886;
		m_points[4][1] = 0.108103018168070227;
		m_points[5][0] = 0.108103018168070227;
		m_points[5][1] = 0.445948490915964886;

		m_weights[0] = 0.5 *  0.109951743655321868;
		m_weights[1] = 0.5 *  0.109951743655321868;
		m_weights[2] = 0.5 *  0.109951743655321868;
		m_weights[3] = 0.5 *  0.223381589678011466;
		m_weights[4] = 0.5 *  0.223381589678011466;
		m_weights[5] = 0.5 *  0.223381589678011466;
		break;

	case 5:
		m_order = 5;
		m_num_points = 7;

		if(allocate_memory(m_num_points) != true)
		{
			assert(0 && "Error while constructing Quadrature Rule.");
		}

		m_points[0][0] = 0.333333333333333333;
		m_points[0][1] = 0.333333333333333333;
		m_points[1][0] = 0.101286507323456339;
		m_points[1][1] = 0.101286507323456339;
		m_points[2][0] = 0.101286507323456339;
		m_points[2][1] = 0.797426985353087322;
		m_points[3][0] = 0.797426985353087322;
		m_points[3][1] = 0.101286507323456339;
		m_points[4][0] = 0.470142064105115090;
		m_points[4][1] = 0.470142064105115090;
		m_points[5][0] = 0.470142064105115090;
		m_points[5][1] = 0.059715871789769820;
		m_points[6][0] = 0.059715871789769820;
		m_points[6][1] = 0.470142064105115090;

		m_weights[0] = 0.5 *  0.225000000000000000;
		m_weights[1] = 0.5 *  0.125939180544827153;
		m_weights[2] = 0.5 *  0.125939180544827153;
		m_weights[3] = 0.5 *  0.125939180544827153;
		m_weights[4] = 0.5 *  0.132394152788506181;
		m_weights[5] = 0.5 *  0.132394152788506181;
		m_weights[6] = 0.5 *  0.132394152788506181;
		break;

	case 6:
		m_order = 6;
		m_num_points = 12;

		if(allocate_memory(m_num_points) != true)
		{
			assert(0 && "Error while constructing Quadrature Rule.");
		}

		m_points[0][0] = 0.063089014491502228;
		m_points[0][1] = 0.063089014491502228;
		m_points[1][0] = 0.063089014491502228;
		m_points[1][1] = 0.873821971016995543;
		m_points[2][0] = 0.873821971016995543;
		m_points[2][1] = 0.063089014491502228;
		m_points[3][0] = 0.249286745170910421;
		m_points[3][1] = 0.249286745170910421;
		m_points[4][0] = 0.249286745170910421;
		m_points[4][1] = 0.501426509658179157;
		m_points[5][0] = 0.501426509658179157;
		m_points[5][1] = 0.249286745170910421;
		m_points[6][0] = 0.053145049844816947;
		m_points[6][1] = 0.310352451033784405;
		m_points[7][0] = 0.053145049844816947;
		m_points[7][1] = 0.636502499121398647;
		m_points[8][0] = 0.310352451033784405;
		m_points[8][1] = 0.053145049844816947;
		m_points[9][0] = 0.310352451033784405;
		m_points[9][1] = 0.636502499121398647;
		m_points[10][0] = 0.636502499121398647;
		m_points[10][1] = 0.053145049844816947;
		m_points[11][0] = 0.636502499121398647;
		m_points[11][1] = 0.310352451033784405;

		m_weights[0] = 0.5 *  0.050844906370206817;
		m_weights[1] = 0.5 *  0.050844906370206817;
		m_weights[2] = 0.5 *  0.050844906370206817;
		m_weights[3] = 0.5 *  0.116786275726379366;
		m_weights[4] = 0.5 *  0.116786275726379366;
		m_weights[5] = 0.5 *  0.116786275726379366;
		m_weights[6] = 0.5 *  0.082851075618373575;
		m_weights[7] = 0.5 *  0.082851075618373575;
		m_weights[8] = 0.5 *  0.082851075618373575;
		m_weights[9] = 0.5 *  0.082851075618373575;
		m_weights[10] = 0.5 *  0.082851075618373575;
		m_weights[11] = 0.5 *  0.082851075618373575;
		break;

	case 7:
		m_order = 7;
		m_num_points = 12;

		if(allocate_memory(m_num_points) != true)
		{
			assert(0 && "Error while constructing Quadrature Rule.");
		}

		m_points[0][0] = 0.062382265094402118;
		m_points[0][1] = 0.067517867073916085;
		m_points[1][0] = 0.067517867073916085;
		m_points[1][1] = 0.870099867831681796;
		m_points[2][0] = 0.870099867831681796;
		m_points[2][1] = 0.062382265094402118;
		m_points[3][0] = 0.055225456656926612;
		m_points[3][1] = 0.321502493851981823;
		m_points[4][0] = 0.321502493851981823;
		m_points[4][1] = 0.623272049491091566;
		m_points[5][0] = 0.623272049491091566;
		m_points[5][1] = 0.055225456656926612;
		m_points[6][0] = 0.034324302945097146;
		m_points[6][1] = 0.660949196186735658;
		m_points[7][0] = 0.660949196186735658;
		m_points[7][1] = 0.304726500868167196;
		m_points[8][0] = 0.304726500868167196;
		m_points[8][1] = 0.034324302945097146;
		m_points[9][0] = 0.515842334353591779;
		m_points[9][1] = 0.277716166976391783;
		m_points[10][0] = 0.277716166976391783;
		m_points[10][1] = 0.206441498670016438;
		m_points[11][0] = 0.206441498670016438;
		m_points[11][1] = 0.515842334353591779;

		m_weights[0] = 0.5 *  0.053034056314872503;
		m_weights[1] = 0.5 *  0.053034056314872503;
		m_weights[2] = 0.5 *  0.053034056314872503;
		m_weights[3] = 0.5 *  0.087762817428892110;
		m_weights[4] = 0.5 *  0.087762817428892110;
		m_weights[5] = 0.5 *  0.087762817428892110;
		m_weights[6] = 0.5 *  0.057550085569963171;
		m_weights[7] = 0.5 *  0.057550085569963171;
		m_weights[8] = 0.5 *  0.057550085569963171;
		m_weights[9] = 0.5 *  0.134986374019605549;
		m_weights[10] = 0.5 *  0.134986374019605549;
		m_weights[11] = 0.5 *  0.134986374019605549;
		break;

	case 8:
		m_order = 8;
		m_num_points = 16;

		if(allocate_memory(m_num_points) != true)
		{
			assert(0 && "Error while constructing Quadrature Rule.");
		}

		m_points[0][0] = 0.333333333333333333;
		m_points[0][1] = 0.333333333333333333;
		m_points[1][0] = 0.170569307751760207;
		m_points[1][1] = 0.170569307751760207;
		m_points[2][0] = 0.170569307751760207;
		m_points[2][1] = 0.658861384496479587;
		m_points[3][0] = 0.658861384496479587;
		m_points[3][1] = 0.170569307751760207;
		m_points[4][0] = 0.050547228317030975;
		m_points[4][1] = 0.050547228317030975;
		m_points[5][0] = 0.050547228317030975;
		m_points[5][1] = 0.898905543365938049;
		m_points[6][0] = 0.898905543365938049;
		m_points[6][1] = 0.050547228317030975;
		m_points[7][0] = 0.459292588292723156;
		m_points[7][1] = 0.459292588292723156;
		m_points[8][0] = 0.459292588292723156;
		m_points[8][1] = 0.081414823414553688;
		m_points[9][0] = 0.081414823414553688;
		m_points[9][1] = 0.459292588292723156;
		m_points[10][0] = 0.728492392955404281;
		m_points[10][1] = 0.263112829634638113;
		m_points[11][0] = 0.728492392955404281;
		m_points[11][1] = 0.008394777409957605;
		m_points[12][0] = 0.263112829634638113;
		m_points[12][1] = 0.728492392955404281;
		m_points[13][0] = 0.263112829634638113;
		m_points[13][1] = 0.008394777409957605;
		m_points[14][0] = 0.008394777409957605;
		m_points[14][1] = 0.728492392955404281;
		m_points[15][0] = 0.008394777409957605;
		m_points[15][1] = 0.263112829634638113;

		m_weights[0] = 0.5 *  0.144315607677787168;
		m_weights[1] = 0.5 *  0.103217370534718250;
		m_weights[2] = 0.5 *  0.103217370534718250;
		m_weights[3] = 0.5 *  0.103217370534718250;
		m_weights[4] = 0.5 *  0.032458497623198080;
		m_weights[5] = 0.5 *  0.032458497623198080;
		m_weights[6] = 0.5 *  0.032458497623198080;
		m_weights[7] = 0.5 *  0.095091634267284625;
		m_weights[8] = 0.5 *  0.095091634267284625;
		m_weights[9] = 0.5 *  0.095091634267284625;
		m_weights[10] = 0.5 *  0.027230314174434994;
		m_weights[11] = 0.5 *  0.027230314174434994;
		m_weights[12] = 0.5 *  0.027230314174434994;
		m_weights[13] = 0.5 *  0.027230314174434994;
		m_weights[14] = 0.5 *  0.027230314174434994;
		m_weights[15] = 0.5 *  0.027230314174434994;
		break;

	case 9:
		m_order = 9;
		m_num_points = 19;

		if(allocate_memory(m_num_points) != true)
		{
			assert(0 && "Error while constructing Quadrature Rule.");
		}

		m_points[0][0] = 0.333333333333333333;
		m_points[0][1] = 0.333333333333333333;
		m_points[1][0] = 0.489682519198737628;
		m_points[1][1] = 0.489682519198737628;
		m_points[2][0] = 0.489682519198737628;
		m_points[2][1] = 0.020634961602524744;
		m_points[3][0] = 0.020634961602524744;
		m_points[3][1] = 0.489682519198737628;
		m_points[4][0] = 0.437089591492936637;
		m_points[4][1] = 0.437089591492936637;
		m_points[5][0] = 0.437089591492936637;
		m_points[5][1] = 0.125820817014126725;
		m_points[6][0] = 0.125820817014126725;
		m_points[6][1] = 0.437089591492936637;
		m_points[7][0] = 0.188203535619032730;
		m_points[7][1] = 0.188203535619032730;
		m_points[8][0] = 0.188203535619032730;
		m_points[8][1] = 0.623592928761934540;
		m_points[9][0] = 0.623592928761934540;
		m_points[9][1] = 0.188203535619032730;
		m_points[10][0] = 0.044729513394452710;
		m_points[10][1] = 0.044729513394452710;
		m_points[11][0] = 0.044729513394452710;
		m_points[11][1] = 0.910540973211094580;
		m_points[12][0] = 0.910540973211094580;
		m_points[12][1] = 0.044729513394452710;
		m_points[13][0] = 0.741198598784498021;
		m_points[13][1] = 0.036838412054736284;
		m_points[14][0] = 0.741198598784498021;
		m_points[14][1] = 0.221962989160765696;
		m_points[15][0] = 0.036838412054736284;
		m_points[15][1] = 0.741198598784498021;
		m_points[16][0] = 0.036838412054736284;
		m_points[16][1] = 0.221962989160765696;
		m_points[17][0] = 0.221962989160765696;
		m_points[17][1] = 0.741198598784498021;
		m_points[18][0] = 0.221962989160765696;
		m_points[18][1] = 0.036838412054736284;

		m_weights[0] = 0.5 *  0.097135796282798834;
		m_weights[1] = 0.5 *  0.031334700227139071;
		m_weights[2] = 0.5 *  0.031334700227139071;
		m_weights[3] = 0.5 *  0.031334700227139071;
		m_weights[4] = 0.5 *  0.077827541004774279;
		m_weights[5] = 0.5 *  0.077827541004774279;
		m_weights[6] = 0.5 *  0.077827541004774279;
		m_weights[7] = 0.5 *  0.079647738927210253;
		m_weights[8] = 0.5 *  0.079647738927210253;
		m_weights[9] = 0.5 *  0.079647738927210253;
		m_weights[10] = 0.5 *  0.025577675658698031;
		m_weights[11] = 0.5 *  0.025577675658698031;
		m_weights[12] = 0.5 *  0.025577675658698031;
		m_weights[13] = 0.5 *  0.043283539377289377;
		m_weights[14] = 0.5 *  0.043283539377289377;
		m_weights[15] = 0.5 *  0.043283539377289377;
		m_weights[16] = 0.5 *  0.043283539377289377;
		m_weights[17] = 0.5 *  0.043283539377289377;
		m_weights[18] = 0.5 *  0.043283539377289377;
		break;

	case 10:
		m_order = 10;
		m_num_points = 25;

		if(allocate_memory(m_num_points) != true)
		{
			assert(0 && "Error while constructing Quadrature Rule.");
		}

		m_points[0][0] = 0.333333333333333333;
		m_points[0][1] = 0.333333333333333333;
		m_points[1][0] = 0.425086210602090573;
		m_points[1][1] = 0.425086210602090573;
		m_points[2][0] = 0.425086210602090573;
		m_points[2][1] = 0.149827578795818854;
		m_points[3][0] = 0.149827578795818854;
		m_points[3][1] = 0.425086210602090573;
		m_points[4][0] = 0.023308867510000191;
		m_points[4][1] = 0.023308867510000191;
		m_points[5][0] = 0.023308867510000191;
		m_points[5][1] = 0.953382264979999619;
		m_points[6][0] = 0.953382264979999619;
		m_points[6][1] = 0.023308867510000191;
		m_points[7][0] = 0.628307400213492556;
		m_points[7][1] = 0.223766973576973006;
		m_points[8][0] = 0.628307400213492556;
		m_points[8][1] = 0.147925626209534437;
		m_points[9][0] = 0.223766973576973006;
		m_points[9][1] = 0.628307400213492556;
		m_points[10][0] = 0.223766973576973006;
		m_points[10][1] = 0.147925626209534437;
		m_points[11][0] = 0.147925626209534437;
		m_points[11][1] = 0.628307400213492556;
		m_points[12][0] = 0.147925626209534437;
		m_points[12][1] = 0.223766973576973006;
		m_points[13][0] = 0.611313826181397649;
		m_points[13][1] = 0.358740141864431465;
		m_points[14][0] = 0.611313826181397649;
		m_points[14][1] = 0.029946031954170887;
		m_points[15][0] = 0.358740141864431465;
		m_points[15][1] = 0.611313826181397649;
		m_points[16][0] = 0.358740141864431465;
		m_points[16][1] = 0.029946031954170887;
		m_points[17][0] = 0.029946031954170887;
		m_points[17][1] = 0.611313826181397649;
		m_points[18][0] = 0.029946031954170887;
		m_points[18][1] = 0.358740141864431465;
		m_points[19][0] = 0.821072069985629373;
		m_points[19][1] = 0.143295370426867145;
		m_points[20][0] = 0.821072069985629373;
		m_points[20][1] = 0.035632559587503481;
		m_points[21][0] = 0.143295370426867145;
		m_points[21][1] = 0.821072069985629373;
		m_points[22][0] = 0.143295370426867145;
		m_points[22][1] = 0.035632559587503481;
		m_points[23][0] = 0.035632559587503481;
		m_points[23][1] = 0.821072069985629373;
		m_points[24][0] = 0.035632559587503481;
		m_points[24][1] = 0.143295370426867145;

		m_weights[0] = 0.5 *  0.079894504741239708;
		m_weights[1] = 0.5 *  0.071123802232377335;
		m_weights[2] = 0.5 *  0.071123802232377335;
		m_weights[3] = 0.5 *  0.071123802232377335;
		m_weights[4] = 0.5 *  0.008223818690464196;
		m_weights[5] = 0.5 *  0.008223818690464196;
		m_weights[6] = 0.5 *  0.008223818690464196;
		m_weights[7] = 0.5 *  0.045430592296170018;
		m_weights[8] = 0.5 *  0.045430592296170018;
		m_weights[9] = 0.5 *  0.045430592296170018;
		m_weights[10] = 0.5 *  0.045430592296170018;
		m_weights[11] = 0.5 *  0.045430592296170018;
		m_weights[12] = 0.5 *  0.045430592296170018;
		m_weights[13] = 0.5 *  0.037359856234305277;
		m_weights[14] = 0.5 *  0.037359856234305277;
		m_weights[15] = 0.5 *  0.037359856234305277;
		m_weights[16] = 0.5 *  0.037359856234305277;
		m_weights[17] = 0.5 *  0.037359856234305277;
		m_weights[18] = 0.5 *  0.037359856234305277;
		m_weights[19] = 0.5 *  0.030886656884563989;
		m_weights[20] = 0.5 *  0.030886656884563989;
		m_weights[21] = 0.5 *  0.030886656884563989;
		m_weights[22] = 0.5 *  0.030886656884563989;
		m_weights[23] = 0.5 *  0.030886656884563989;
		m_weights[24] = 0.5 *  0.030886656884563989;
		break;

	case 11:
		m_order = 11;
		m_num_points = 28;

		if(allocate_memory(m_num_points) != true)
		{
			assert(0 && "Error while constructing Quadrature Rule.");
		}

		m_points[0][0] = 0.858870281282636704;
		m_points[0][1] = 0.141129718717363296;
		m_points[1][0] = 0.858870281282636704;
		m_points[1][1] = 0.000000000000000000;
		m_points[2][0] = 0.141129718717363296;
		m_points[2][1] = 0.858870281282636704;
		m_points[3][0] = 0.141129718717363296;
		m_points[3][1] = 0.000000000000000000;
		m_points[4][0] = 0.000000000000000000;
		m_points[4][1] = 0.858870281282636704;
		m_points[5][0] = 0.000000000000000000;
		m_points[5][1] = 0.141129718717363296;
		m_points[6][0] = 0.333333333333333333;
		m_points[6][1] = 0.333333333333333333;
		m_points[7][0] = 0.025989140928287395;
		m_points[7][1] = 0.025989140928287395;
		m_points[8][0] = 0.025989140928287395;
		m_points[8][1] = 0.948021718143425209;
		m_points[9][0] = 0.948021718143425209;
		m_points[9][1] = 0.025989140928287395;
		m_points[10][0] = 0.094287502647922496;
		m_points[10][1] = 0.094287502647922496;
		m_points[11][0] = 0.094287502647922496;
		m_points[11][1] = 0.811424994704155009;
		m_points[12][0] = 0.811424994704155009;
		m_points[12][1] = 0.094287502647922496;
		m_points[13][0] = 0.494636775017213814;
		m_points[13][1] = 0.494636775017213814;
		m_points[14][0] = 0.494636775017213814;
		m_points[14][1] = 0.010726449965572373;
		m_points[15][0] = 0.010726449965572373;
		m_points[15][1] = 0.494636775017213814;
		m_points[16][0] = 0.207343382614511333;
		m_points[16][1] = 0.207343382614511333;
		m_points[17][0] = 0.207343382614511333;
		m_points[17][1] = 0.585313234770977333;
		m_points[18][0] = 0.585313234770977333;
		m_points[18][1] = 0.207343382614511333;
		m_points[19][0] = 0.438907805700492095;
		m_points[19][1] = 0.438907805700492095;
		m_points[20][0] = 0.438907805700492095;
		m_points[20][1] = 0.122184388599015810;
		m_points[21][0] = 0.122184388599015810;
		m_points[21][1] = 0.438907805700492095;
		m_points[22][0] = 0.677937654882590402;
		m_points[22][1] = 0.044841677589130443;
		m_points[23][0] = 0.677937654882590402;
		m_points[23][1] = 0.277220667528279155;
		m_points[24][0] = 0.044841677589130443;
		m_points[24][1] = 0.677937654882590402;
		m_points[25][0] = 0.044841677589130443;
		m_points[25][1] = 0.277220667528279155;
		m_points[26][0] = 0.277220667528279155;
		m_points[26][1] = 0.677937654882590402;
		m_points[27][0] = 0.277220667528279155;
		m_points[27][1] = 0.044841677589130443;

		m_weights[0] = 0.5 *  0.007362383783300554;
		m_weights[1] = 0.5 *  0.007362383783300554;
		m_weights[2] = 0.5 *  0.007362383783300554;
		m_weights[3] = 0.5 *  0.007362383783300554;
		m_weights[4] = 0.5 *  0.007362383783300554;
		m_weights[5] = 0.5 *  0.007362383783300554;
		m_weights[6] = 0.5 *  0.087977301162232239;
		m_weights[7] = 0.5 *  0.008744311553736023;
		m_weights[8] = 0.5 *  0.008744311553736023;
		m_weights[9] = 0.5 *  0.008744311553736023;
		m_weights[10] = 0.5 *  0.038081571993934938;
		m_weights[11] = 0.5 *  0.038081571993934938;
		m_weights[12] = 0.5 *  0.038081571993934938;
		m_weights[13] = 0.5 *  0.018855448056131292;
		m_weights[14] = 0.5 *  0.018855448056131292;
		m_weights[15] = 0.5 *  0.018855448056131292;
		m_weights[16] = 0.5 *  0.072159697544739526;
		m_weights[17] = 0.5 *  0.072159697544739526;
		m_weights[18] = 0.5 *  0.072159697544739526;
		m_weights[19] = 0.5 *  0.069329138705535900;
		m_weights[20] = 0.5 *  0.069329138705535900;
		m_weights[21] = 0.5 *  0.069329138705535900;
		m_weights[22] = 0.5 *  0.041056315429288567;
		m_weights[23] = 0.5 *  0.041056315429288567;
		m_weights[24] = 0.5 *  0.041056315429288567;
		m_weights[25] = 0.5 *  0.041056315429288567;
		m_weights[26] = 0.5 *  0.041056315429288567;
		m_weights[27] = 0.5 *  0.041056315429288567;
		break;

	case 12:
		m_order = 12;
		m_num_points = 33;

		if(allocate_memory(m_num_points) != true)
		{
			assert(0 && "Error while constructing Quadrature Rule.");
		}

		m_points[0][0] = 0.023565220452390000;
		m_points[0][1] = 0.488217389773805000;
		m_points[1][0] = 0.488217389773805000;
		m_points[1][1] = 0.023565220452390000;
		m_points[2][0] = 0.488217389773805000;
		m_points[2][1] = 0.488217389773805000;
		m_points[3][0] = 0.439724392294460000;
		m_points[3][1] = 0.439724392294460000;
		m_points[4][0] = 0.439724392294460000;
		m_points[4][1] = 0.120551215411079000;
		m_points[5][0] = 0.120551215411079000;
		m_points[5][1] = 0.439724392294460000;
		m_points[6][0] = 0.271210385012116000;
		m_points[6][1] = 0.271210385012116000;
		m_points[7][0] = 0.271210385012116000;
		m_points[7][1] = 0.457579229975768000;
		m_points[8][0] = 0.457579229975768000;
		m_points[8][1] = 0.271210385012116000;
		m_points[9][0] = 0.127576145541586000;
		m_points[9][1] = 0.127576145541586000;
		m_points[10][0] = 0.127576145541586000;
		m_points[10][1] = 0.744847708916827900;
		m_points[11][0] = 0.744847708916827900;
		m_points[11][1] = 0.127576145541586000;
		m_points[12][0] = 0.021317350453210000;
		m_points[12][1] = 0.021317350453210000;
		m_points[13][0] = 0.021317350453210000;
		m_points[13][1] = 0.957365299093579900;
		m_points[14][0] = 0.957365299093579900;
		m_points[14][1] = 0.021317350453210000;
		m_points[15][0] = 0.115343494534698000;
		m_points[15][1] = 0.275713269685514000;
		m_points[16][0] = 0.115343494534698000;
		m_points[16][1] = 0.608943235779787900;
		m_points[17][0] = 0.275713269685514000;
		m_points[17][1] = 0.115343494534698000;
		m_points[18][0] = 0.275713269685514000;
		m_points[18][1] = 0.608943235779787900;
		m_points[19][0] = 0.608943235779787900;
		m_points[19][1] = 0.115343494534698000;
		m_points[20][0] = 0.608943235779787900;
		m_points[20][1] = 0.275713269685514000;
		m_points[21][0] = 0.022838332222257000;
		m_points[21][1] = 0.281325580989940000;
		m_points[22][0] = 0.022838332222257000;
		m_points[22][1] = 0.695836086787803100;
		m_points[23][0] = 0.281325580989940000;
		m_points[23][1] = 0.022838332222257000;
		m_points[24][0] = 0.281325580989940000;
		m_points[24][1] = 0.695836086787803100;
		m_points[25][0] = 0.695836086787803100;
		m_points[25][1] = 0.022838332222257000;
		m_points[26][0] = 0.695836086787803100;
		m_points[26][1] = 0.281325580989940000;
		m_points[27][0] = 0.025734050548330000;
		m_points[27][1] = 0.116251915907597000;
		m_points[28][0] = 0.025734050548330000;
		m_points[28][1] = 0.858014033544073000;
		m_points[29][0] = 0.116251915907597000;
		m_points[29][1] = 0.025734050548330000;
		m_points[30][0] = 0.116251915907597000;
		m_points[30][1] = 0.858014033544073000;
		m_points[31][0] = 0.858014033544073000;
		m_points[31][1] = 0.025734050548330000;
		m_points[32][0] = 0.858014033544073000;
		m_points[32][1] = 0.116251915907597000;

		m_weights[0] = 0.5 *  0.025731066440455000;
		m_weights[1] = 0.5 *  0.025731066440455000;
		m_weights[2] = 0.5 *  0.025731066440455000;
		m_weights[3] = 0.5 *  0.043692544538038000;
		m_weights[4] = 0.5 *  0.043692544538038000;
		m_weights[5] = 0.5 *  0.043692544538038000;
		m_weights[6] = 0.5 *  0.062858224217885000;
		m_weights[7] = 0.5 *  0.062858224217885000;
		m_weights[8] = 0.5 *  0.062858224217885000;
		m_weights[9] = 0.5 *  0.034796112930709000;
		m_weights[10] = 0.5 *  0.034796112930709000;
		m_weights[11] = 0.5 *  0.034796112930709000;
		m_weights[12] = 0.5 *  0.006166261051559000;
		m_weights[13] = 0.5 *  0.006166261051559000;
		m_weights[14] = 0.5 *  0.006166261051559000;
		m_weights[15] = 0.5 *  0.040371557766381000;
		m_weights[16] = 0.5 *  0.040371557766381000;
		m_weights[17] = 0.5 *  0.040371557766381000;
		m_weights[18] = 0.5 *  0.040371557766381000;
		m_weights[19] = 0.5 *  0.040371557766381000;
		m_weights[20] = 0.5 *  0.040371557766381000;
		m_weights[21] = 0.5 *  0.022356773202303000;
		m_weights[22] = 0.5 *  0.022356773202303000;
		m_weights[23] = 0.5 *  0.022356773202303000;
		m_weights[24] = 0.5 *  0.022356773202303000;
		m_weights[25] = 0.5 *  0.022356773202303000;
		m_weights[26] = 0.5 *  0.022356773202303000;
		m_weights[27] = 0.5 *  0.017316231108659000;
		m_weights[28] = 0.5 *  0.017316231108659000;
		m_weights[29] = 0.5 *  0.017316231108659000;
		m_weights[30] = 0.5 *  0.017316231108659000;
		m_weights[31] = 0.5 *  0.017316231108659000;
		m_weights[32] = 0.5 *  0.017316231108659000;
		break;

	default: assert(0 && "Order not availabile. Can not construct GaussQuadrature.\n");
	}
}
}; // namespace ug

 // register quadratures at factory
namespace {
using namespace ug;

template <>
std::vector<const QuadratureRule<ReferenceTriangle>* > QuadratureRuleFactory<ReferenceTriangle>::m_rules =
	std::vector<const QuadratureRule<ReferenceTriangle>* >();

GaussQuadrature<ReferenceTriangle> gaussQuadratureReferenceTriangle_1(1);
GaussQuadrature<ReferenceTriangle> gaussQuadratureReferenceTriangle_2(2);
GaussQuadrature<ReferenceTriangle> gaussQuadratureReferenceTriangle_3(3);
GaussQuadrature<ReferenceTriangle> gaussQuadratureReferenceTriangle_4(4);
GaussQuadrature<ReferenceTriangle> gaussQuadratureReferenceTriangle_5(5);
GaussQuadrature<ReferenceTriangle> gaussQuadratureReferenceTriangle_6(6);
GaussQuadrature<ReferenceTriangle> gaussQuadratureReferenceTriangle_7(7);
GaussQuadrature<ReferenceTriangle> gaussQuadratureReferenceTriangle_8(8);
GaussQuadrature<ReferenceTriangle> gaussQuadratureReferenceTriangle_9(9);
GaussQuadrature<ReferenceTriangle> gaussQuadratureReferenceTriangle_10(10);
GaussQuadrature<ReferenceTriangle> gaussQuadratureReferenceTriangle_11(11);
GaussQuadrature<ReferenceTriangle> gaussQuadratureReferenceTriangle_12(12);

static const bool registered_1 = QuadratureRuleFactory<ReferenceTriangle>::instance().register_rule(gaussQuadratureReferenceTriangle_1);
static const bool registered_2 = QuadratureRuleFactory<ReferenceTriangle>::instance().register_rule(gaussQuadratureReferenceTriangle_2);
static const bool registered_3 = QuadratureRuleFactory<ReferenceTriangle>::instance().register_rule(gaussQuadratureReferenceTriangle_3);
static const bool registered_4 = QuadratureRuleFactory<ReferenceTriangle>::instance().register_rule(gaussQuadratureReferenceTriangle_4);
static const bool registered_5 = QuadratureRuleFactory<ReferenceTriangle>::instance().register_rule(gaussQuadratureReferenceTriangle_5);
static const bool registered_6 = QuadratureRuleFactory<ReferenceTriangle>::instance().register_rule(gaussQuadratureReferenceTriangle_6);
static const bool registered_7 = QuadratureRuleFactory<ReferenceTriangle>::instance().register_rule(gaussQuadratureReferenceTriangle_7);
static const bool registered_8 = QuadratureRuleFactory<ReferenceTriangle>::instance().register_rule(gaussQuadratureReferenceTriangle_8);
static const bool registered_9 = QuadratureRuleFactory<ReferenceTriangle>::instance().register_rule(gaussQuadratureReferenceTriangle_9);
static const bool registered_10 = QuadratureRuleFactory<ReferenceTriangle>::instance().register_rule(gaussQuadratureReferenceTriangle_10);
static const bool registered_11 = QuadratureRuleFactory<ReferenceTriangle>::instance().register_rule(gaussQuadratureReferenceTriangle_11);
static const bool registered_12 = QuadratureRuleFactory<ReferenceTriangle>::instance().register_rule(gaussQuadratureReferenceTriangle_12);

};
