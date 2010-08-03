//This file is parsed from UG 3.9.


#include "quadrature.h"

namespace ug{

template <>
GaussQuadrature<ReferenceQuadrilateral>::GaussQuadrature(int order)
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

		m_points[0][0] = 0.500000000000000000;
		m_points[0][1] = 0.500000000000000000;

		m_weights[0] =  1.000000000000000000;
		break;

	case 2:
		m_order = 2;
		m_num_points = 4;

		if(allocate_memory(m_num_points) != true)
		{
			assert(0 && "Error while constructing Quadrature Rule.");
		}

		m_points[0][0] = 0.211324865405187000;
		m_points[0][1] = 0.211324865405187000;
		m_points[1][0] = 0.788675134594813000;
		m_points[1][1] = 0.211324865405187000;
		m_points[2][0] = 0.211324865405187000;
		m_points[2][1] = 0.788675134594813000;
		m_points[3][0] = 0.788675134594813000;
		m_points[3][1] = 0.788675134594813000;

		m_weights[0] =  0.250000000000000000;
		m_weights[1] =  0.250000000000000000;
		m_weights[2] =  0.250000000000000000;
		m_weights[3] =  0.250000000000000000;
		break;

	case 3:
		m_order = 3;
		m_num_points = 4;

		if(allocate_memory(m_num_points) != true)
		{
			assert(0 && "Error while constructing Quadrature Rule.");
		}

		m_points[0][0] = 0.091751709536136984;
		m_points[0][1] = 0.500000000000000000;
		m_points[1][0] = 0.500000000000000000;
		m_points[1][1] = 0.091751709536136984;
		m_points[2][0] = 0.500000000000000000;
		m_points[2][1] = 0.908248290463863016;
		m_points[3][0] = 0.908248290463863016;
		m_points[3][1] = 0.500000000000000000;

		m_weights[0] =  0.250000000000000000;
		m_weights[1] =  0.250000000000000000;
		m_weights[2] =  0.250000000000000000;
		m_weights[3] =  0.250000000000000000;
		break;

	case 4:
		m_order = 4;
		m_num_points = 6;

		if(allocate_memory(m_num_points) != true)
		{
			assert(0 && "Error while constructing Quadrature Rule.");
		}

		m_points[0][0] = 0.500000000000000000;
		m_points[0][1] = 0.500000000000000000;
		m_points[1][0] = 0.983045891539647952;
		m_points[1][1] = 0.500000000000000000;
		m_points[2][0] = 0.727801863918096421;
		m_points[2][1] = 0.074042673347699754;
		m_points[3][0] = 0.727801863918096421;
		m_points[3][1] = 0.925957326652300246;
		m_points[4][0] = 0.134185024213432735;
		m_points[4][1] = 0.184543605511622987;
		m_points[5][0] = 0.134185024213432735;
		m_points[5][1] = 0.815456394488377013;

		m_weights[0] =  0.285714285714285714;
		m_weights[1] =  0.109890109890109890;
		m_weights[2] =  0.141518051751883026;
		m_weights[3] =  0.141518051751883026;
		m_weights[4] =  0.160679750445919171;
		m_weights[5] =  0.160679750445919171;
		break;

	case 5:
		m_order = 5;
		m_num_points = 7;

		if(allocate_memory(m_num_points) != true)
		{
			assert(0 && "Error while constructing Quadrature Rule.");
		}

		m_points[0][0] = 0.500000000000000000;
		m_points[0][1] = 0.500000000000000000;
		m_points[1][0] = 0.016954108460352048;
		m_points[1][1] = 0.500000000000000000;
		m_points[2][0] = 0.983045891539647952;
		m_points[2][1] = 0.500000000000000000;
		m_points[3][0] = 0.211324865405187118;
		m_points[3][1] = 0.112701665379258311;
		m_points[4][0] = 0.211324865405187118;
		m_points[4][1] = 0.887298334620741689;
		m_points[5][0] = 0.788675134594812882;
		m_points[5][1] = 0.112701665379258311;
		m_points[6][0] = 0.788675134594812882;
		m_points[6][1] = 0.887298334620741689;

		m_weights[0] =  0.285714285714285714;
		m_weights[1] =  0.079365079365079365;
		m_weights[2] =  0.079365079365079365;
		m_weights[3] =  0.138888888888888889;
		m_weights[4] =  0.138888888888888889;
		m_weights[5] =  0.138888888888888889;
		m_weights[6] =  0.138888888888888889;
		break;

	case 6:
		m_order = 6;
		m_num_points = 10;

		if(allocate_memory(m_num_points) != true)
		{
			assert(0 && "Error while constructing Quadrature Rule.");
		}

		m_points[0][0] = 0.918202816848812802;
		m_points[0][1] = 0.500000000000000000;
		m_points[1][0] = 0.321269917304346408;
		m_points[1][1] = 0.500000000000000000;
		m_points[2][0] = 0.936050765596565295;
		m_points[2][1] = 0.055617992672617733;
		m_points[3][0] = 0.936050765596565295;
		m_points[3][1] = 0.944382007327382267;
		m_points[4][0] = 0.652992581077713331;
		m_points[4][1] = 0.197571180267657486;
		m_points[5][0] = 0.652992581077713331;
		m_points[5][1] = 0.802428819732342514;
		m_points[6][0] = 0.294864550266671081;
		m_points[6][1] = 0.022276246679468126;
		m_points[7][0] = 0.294864550266671081;
		m_points[7][1] = 0.977723753320531874;
		m_points[8][0] = 0.063565344421560330;
		m_points[8][1] = 0.217270003280622977;
		m_points[9][0] = 0.063565344421560330;
		m_points[9][1] = 0.782729996719377023;

		m_weights[0] =  0.113835811428543360;
		m_weights[1] =  0.206848993300741373;
		m_weights[2] =  0.036000221149911346;
		m_weights[3] =  0.036000221149911346;
		m_weights[4] =  0.167064776065666287;
		m_weights[5] =  0.167064776065666287;
		m_weights[6] =  0.056368501222669839;
		m_weights[7] =  0.056368501222669839;
		m_weights[8] =  0.080224099197110161;
		m_weights[9] =  0.080224099197110161;
		break;

	case 7:
		m_order = 7;
		m_num_points = 12;

		if(allocate_memory(m_num_points) != true)
		{
			assert(0 && "Error while constructing Quadrature Rule.");
		}

		m_points[0][0] = 0.037089950113724269;
		m_points[0][1] = 0.500000000000000000;
		m_points[1][0] = 0.500000000000000000;
		m_points[1][1] = 0.037089950113724269;
		m_points[2][0] = 0.500000000000000000;
		m_points[2][1] = 0.962910049886275731;
		m_points[3][0] = 0.962910049886275731;
		m_points[3][1] = 0.500000000000000000;
		m_points[4][0] = 0.309722783395842172;
		m_points[4][1] = 0.309722783395842172;
		m_points[5][0] = 0.309722783395842172;
		m_points[5][1] = 0.690277216604157828;
		m_points[6][0] = 0.690277216604157828;
		m_points[6][1] = 0.309722783395842172;
		m_points[7][0] = 0.690277216604157828;
		m_points[7][1] = 0.690277216604157828;
		m_points[8][0] = 0.097010108540700628;
		m_points[8][1] = 0.097010108540700628;
		m_points[9][0] = 0.097010108540700628;
		m_points[9][1] = 0.902989891459299372;
		m_points[10][0] = 0.902989891459299372;
		m_points[10][1] = 0.097010108540700628;
		m_points[11][0] = 0.902989891459299372;
		m_points[11][1] = 0.902989891459299372;

		m_weights[0] =  0.060493827160493827;
		m_weights[1] =  0.060493827160493827;
		m_weights[2] =  0.060493827160493827;
		m_weights[3] =  0.060493827160493827;
		m_weights[4] =  0.130148229166848614;
		m_weights[5] =  0.130148229166848614;
		m_weights[6] =  0.130148229166848614;
		m_weights[7] =  0.130148229166848614;
		m_weights[8] =  0.059357943672657559;
		m_weights[9] =  0.059357943672657559;
		m_weights[10] =  0.059357943672657559;
		m_weights[11] =  0.059357943672657559;
		break;

	case 8:
		m_order = 8;
		m_num_points = 16;

		if(allocate_memory(m_num_points) != true)
		{
			assert(0 && "Error while constructing Quadrature Rule.");
		}

		m_points[0][0] = 0.500000000000000000;
		m_points[0][1] = 0.500000000000000000;
		m_points[1][0] = 0.878814588830252724;
		m_points[1][1] = 0.500000000000000000;
		m_points[2][0] = 0.381564078872149053;
		m_points[2][1] = 0.500000000000000000;
		m_points[3][0] = 0.005141035477736659;
		m_points[3][1] = 0.500000000000000000;
		m_points[4][0] = 0.975260477822833442;
		m_points[4][1] = 0.180454347549815173;
		m_points[5][0] = 0.975260477822833442;
		m_points[5][1] = 0.819545652450184827;
		m_points[6][0] = 0.831941368442816589;
		m_points[6][1] = 0.031465461537504778;
		m_points[7][0] = 0.831941368442816589;
		m_points[7][1] = 0.968534538462495222;
		m_points[8][0] = 0.652105340862052255;
		m_points[8][1] = 0.231458234729253161;
		m_points[9][0] = 0.652105340862052255;
		m_points[9][1] = 0.768541765270746839;
		m_points[10][0] = 0.381751640731940260;
		m_points[10][1] = 0.056405746775187680;
		m_points[11][0] = 0.381751640731940260;
		m_points[11][1] = 0.943594253224812320;
		m_points[12][0] = 0.150523261956718220;
		m_points[12][1] = 0.252650589664901834;
		m_points[13][0] = 0.150523261956718220;
		m_points[13][1] = 0.747349410335098166;
		m_points[14][0] = 0.049804612894210133;
		m_points[14][1] = 0.051252090860116231;
		m_points[15][0] = 0.049804612894210133;
		m_points[15][1] = 0.948747909139883769;

		m_weights[0] =  0.013841176405359943;
		m_weights[1] =  0.101097342181518853;
		m_weights[2] =  0.133386651238158614;
		m_weights[3] =  0.029263547196684802;
		m_weights[4] =  0.031403604403436701;
		m_weights[5] =  0.031403604403436701;
		m_weights[6] =  0.034136146183397117;
		m_weights[7] =  0.034136146183397117;
		m_weights[8] =  0.120852119802814238;
		m_weights[9] =  0.120852119802814238;
		m_weights[10] =  0.063132126607385923;
		m_weights[11] =  0.063132126607385923;
		m_weights[12] =  0.090315580970543142;
		m_weights[13] =  0.090315580970543142;
		m_weights[14] =  0.021366063521561773;
		m_weights[15] =  0.021366063521561773;
		break;

	case 9:
		m_order = 9;
		m_num_points = 17;

		if(allocate_memory(m_num_points) != true)
		{
			assert(0 && "Error while constructing Quadrature Rule.");
		}

		m_points[0][0] = 0.500000000000000000;
		m_points[0][1] = 0.500000000000000000;
		m_points[1][0] = 0.015575016819011140;
		m_points[1][1] = 0.184659940134165573;
		m_points[2][0] = 0.184659940134165573;
		m_points[2][1] = 0.984424983180988860;
		m_points[3][0] = 0.815340059865834427;
		m_points[3][1] = 0.015575016819011140;
		m_points[4][0] = 0.984424983180988860;
		m_points[4][1] = 0.815340059865834427;
		m_points[5][0] = 0.036019177020215166;
		m_points[5][1] = 0.875138549989450267;
		m_points[6][0] = 0.124861450010549733;
		m_points[6][1] = 0.036019177020215166;
		m_points[7][0] = 0.875138549989450267;
		m_points[7][1] = 0.963980822979784834;
		m_points[8][0] = 0.963980822979784834;
		m_points[8][1] = 0.124861450010549733;
		m_points[9][0] = 0.238132089892785332;
		m_points[9][1] = 0.273330089432176405;
		m_points[10][0] = 0.273330089432176405;
		m_points[10][1] = 0.761867910107214668;
		m_points[11][0] = 0.726669910567823595;
		m_points[11][1] = 0.238132089892785332;
		m_points[12][0] = 0.761867910107214668;
		m_points[12][1] = 0.726669910567823595;
		m_points[13][0] = 0.073692135333168846;
		m_points[13][1] = 0.538104164096308587;
		m_points[14][0] = 0.461895835903691413;
		m_points[14][1] = 0.073692135333168846;
		m_points[15][0] = 0.538104164096308587;
		m_points[15][1] = 0.926307864666831154;
		m_points[16][0] = 0.926307864666831154;
		m_points[16][1] = 0.461895835903691413;

		m_weights[0] =  0.131687242798353909;
		m_weights[1] =  0.022219844542549677;
		m_weights[2] =  0.022219844542549677;
		m_weights[3] =  0.022219844542549677;
		m_weights[4] =  0.022219844542549677;
		m_weights[5] =  0.028024900532399121;
		m_weights[6] =  0.028024900532399121;
		m_weights[7] =  0.028024900532399121;
		m_weights[8] =  0.028024900532399121;
		m_weights[9] =  0.099570609815517524;
		m_weights[10] =  0.099570609815517524;
		m_weights[11] =  0.099570609815517524;
		m_weights[12] =  0.099570609815517524;
		m_weights[13] =  0.067262834409945201;
		m_weights[14] =  0.067262834409945201;
		m_weights[15] =  0.067262834409945201;
		m_weights[16] =  0.067262834409945201;
		break;

	case 11:
		m_order = 11;
		m_num_points = 24;

		if(allocate_memory(m_num_points) != true)
		{
			assert(0 && "Error while constructing Quadrature Rule.");
		}

		m_points[0][0] = 0.008680388229572264;
		m_points[0][1] = 0.150961947725216218;
		m_points[1][0] = 0.150961947725216218;
		m_points[1][1] = 0.991319611770427736;
		m_points[2][0] = 0.849038052274783782;
		m_points[2][1] = 0.008680388229572264;
		m_points[3][0] = 0.991319611770427736;
		m_points[3][1] = 0.849038052274783782;
		m_points[4][0] = 0.030256808591631546;
		m_points[4][1] = 0.912887917951481969;
		m_points[5][0] = 0.087112082048518031;
		m_points[5][1] = 0.030256808591631546;
		m_points[6][0] = 0.912887917951481969;
		m_points[6][1] = 0.969743191408368454;
		m_points[7][0] = 0.969743191408368454;
		m_points[7][1] = 0.087112082048518031;
		m_points[8][0] = 0.023230235899233992;
		m_points[8][1] = 0.594293069359320977;
		m_points[9][0] = 0.405706930640679023;
		m_points[9][1] = 0.023230235899233992;
		m_points[10][0] = 0.594293069359320977;
		m_points[10][1] = 0.976769764100766008;
		m_points[11][0] = 0.976769764100766008;
		m_points[11][1] = 0.405706930640679023;
		m_points[12][0] = 0.093739725847593450;
		m_points[12][1] = 0.342188283542372902;
		m_points[13][0] = 0.342188283542372902;
		m_points[13][1] = 0.906260274152406550;
		m_points[14][0] = 0.657811716457627098;
		m_points[14][1] = 0.093739725847593450;
		m_points[15][0] = 0.906260274152406550;
		m_points[15][1] = 0.657811716457627098;
		m_points[16][0] = 0.143999043462331847;
		m_points[16][1] = 0.762660125182273881;
		m_points[17][0] = 0.237339874817726119;
		m_points[17][1] = 0.143999043462331847;
		m_points[18][0] = 0.762660125182273881;
		m_points[18][1] = 0.856000956537668153;
		m_points[19][0] = 0.856000956537668153;
		m_points[19][1] = 0.237339874817726119;
		m_points[20][0] = 0.287576375575665375;
		m_points[20][1] = 0.520829035956011184;
		m_points[21][0] = 0.479170964043988816;
		m_points[21][1] = 0.287576375575665375;
		m_points[22][0] = 0.520829035956011184;
		m_points[22][1] = 0.712423624424334625;
		m_points[23][0] = 0.712423624424334625;
		m_points[23][1] = 0.479170964043988816;

		m_weights[0] =  0.012005190837680954;
		m_weights[1] =  0.012005190837680954;
		m_weights[2] =  0.012005190837680954;
		m_weights[3] =  0.012005190837680954;
		m_weights[4] =  0.016517832291137649;
		m_weights[5] =  0.016517832291137649;
		m_weights[6] =  0.016517832291137649;
		m_weights[7] =  0.016517832291137649;
		m_weights[8] =  0.024346694339667041;
		m_weights[9] =  0.024346694339667041;
		m_weights[10] =  0.024346694339667041;
		m_weights[11] =  0.024346694339667041;
		m_weights[12] =  0.052934087499737150;
		m_weights[13] =  0.052934087499737150;
		m_weights[14] =  0.052934087499737150;
		m_weights[15] =  0.052934087499737150;
		m_weights[16] =  0.056406515432215847;
		m_weights[17] =  0.056406515432215847;
		m_weights[18] =  0.056406515432215847;
		m_weights[19] =  0.056406515432215847;
		m_weights[20] =  0.087789679599561359;
		m_weights[21] =  0.087789679599561359;
		m_weights[22] =  0.087789679599561359;
		m_weights[23] =  0.087789679599561359;
		break;

	case 13:
		m_order = 13;
		m_num_points = 33;

		if(allocate_memory(m_num_points) != true)
		{
			assert(0 && "Error while constructing Quadrature Rule.");
		}

		m_points[0][0] = 0.500000000000000000;
		m_points[0][1] = 0.500000000000000000;
		m_points[1][0] = 0.008256658780063868;
		m_points[1][1] = 0.889404855777209711;
		m_points[2][0] = 0.110595144222790289;
		m_points[2][1] = 0.008256658780063868;
		m_points[3][0] = 0.889404855777209711;
		m_points[3][1] = 0.991743341219936132;
		m_points[4][0] = 0.991743341219936132;
		m_points[4][1] = 0.110595144222790289;
		m_points[5][0] = 0.021351150106846317;
		m_points[5][1] = 0.070221997179180536;
		m_points[6][0] = 0.070221997179180536;
		m_points[6][1] = 0.978648849893153683;
		m_points[7][0] = 0.929778002820819464;
		m_points[7][1] = 0.021351150106846317;
		m_points[8][0] = 0.978648849893153683;
		m_points[8][1] = 0.929778002820819464;
		m_points[9][0] = 0.020537414856232571;
		m_points[9][1] = 0.569091729931232677;
		m_points[10][0] = 0.430908270068767323;
		m_points[10][1] = 0.020537414856232571;
		m_points[11][0] = 0.569091729931232677;
		m_points[11][1] = 0.979462585143767429;
		m_points[12][0] = 0.979462585143767429;
		m_points[12][1] = 0.430908270068767323;
		m_points[13][0] = 0.029336387063537382;
		m_points[13][1] = 0.304631891935269500;
		m_points[14][0] = 0.304631891935269500;
		m_points[14][1] = 0.970663612936462618;
		m_points[15][0] = 0.695368108064730500;
		m_points[15][1] = 0.029336387063537382;
		m_points[16][0] = 0.970663612936462618;
		m_points[16][1] = 0.695368108064730500;
		m_points[17][0] = 0.074961663150125712;
		m_points[17][1] = 0.737904312609137953;
		m_points[18][0] = 0.262095687390862047;
		m_points[18][1] = 0.074961663150125712;
		m_points[19][0] = 0.737904312609137953;
		m_points[19][1] = 0.925038336849874288;
		m_points[20][0] = 0.925038336849874288;
		m_points[20][1] = 0.262095687390862047;
		m_points[21][0] = 0.122097321713959282;
		m_points[21][1] = 0.176089181406494634;
		m_points[22][0] = 0.176089181406494634;
		m_points[22][1] = 0.877902678286040718;
		m_points[23][0] = 0.823910818593505366;
		m_points[23][1] = 0.122097321713959282;
		m_points[24][0] = 0.877902678286040718;
		m_points[24][1] = 0.823910818593505366;
		m_points[25][0] = 0.151874960754125293;
		m_points[25][1] = 0.464629245501777532;
		m_points[26][0] = 0.464629245501777532;
		m_points[26][1] = 0.848125039245874707;
		m_points[27][0] = 0.535370754498222468;
		m_points[27][1] = 0.151874960754125293;
		m_points[28][0] = 0.848125039245874707;
		m_points[28][1] = 0.535370754498222468;
		m_points[29][0] = 0.295347719152980578;
		m_points[29][1] = 0.671358278020203395;
		m_points[30][0] = 0.328641721979796605;
		m_points[30][1] = 0.295347719152980578;
		m_points[31][0] = 0.671358278020203395;
		m_points[31][1] = 0.704652280847019422;
		m_points[32][0] = 0.704652280847019422;
		m_points[32][1] = 0.328641721979796605;

		m_weights[0] =  0.075095528857806340;
		m_weights[1] =  0.007497959716124783;
		m_weights[2] =  0.007497959716124783;
		m_weights[3] =  0.007497959716124783;
		m_weights[4] =  0.007497959716124783;
		m_weights[5] =  0.009543605329270917;
		m_weights[6] =  0.009543605329270917;
		m_weights[7] =  0.009543605329270917;
		m_weights[8] =  0.009543605329270917;
		m_weights[9] =  0.015106230954437495;
		m_weights[10] =  0.015106230954437495;
		m_weights[11] =  0.015106230954437495;
		m_weights[12] =  0.015106230954437495;
		m_weights[13] =  0.019373184633276335;
		m_weights[14] =  0.019373184633276335;
		m_weights[15] =  0.019373184633276335;
		m_weights[16] =  0.019373184633276335;
		m_weights[17] =  0.029711166825148900;
		m_weights[18] =  0.029711166825148900;
		m_weights[19] =  0.029711166825148900;
		m_weights[20] =  0.029711166825148900;
		m_weights[21] =  0.032440887592500678;
		m_weights[22] =  0.032440887592500678;
		m_weights[23] =  0.032440887592500678;
		m_weights[24] =  0.032440887592500678;
		m_weights[25] =  0.053335395364297347;
		m_weights[26] =  0.053335395364297347;
		m_weights[27] =  0.053335395364297347;
		m_weights[28] =  0.053335395364297347;
		m_weights[29] =  0.064217687370491959;
		m_weights[30] =  0.064217687370491959;
		m_weights[31] =  0.064217687370491959;
		m_weights[32] =  0.064217687370491959;
		break;

	default: assert(0 && "Order not availabile. Can not construct GaussQuadrature.\n");
	}
}

template <>
bool RegisterQuadratureRule(QuadratureRuleFactory<ReferenceQuadrilateral>& factory)
{
	static GaussQuadrature<ReferenceQuadrilateral> gaussQuadratureReferenceQuadrilateral_1(1);
	static GaussQuadrature<ReferenceQuadrilateral> gaussQuadratureReferenceQuadrilateral_2(2);
	static GaussQuadrature<ReferenceQuadrilateral> gaussQuadratureReferenceQuadrilateral_3(3);
	static GaussQuadrature<ReferenceQuadrilateral> gaussQuadratureReferenceQuadrilateral_4(4);
	static GaussQuadrature<ReferenceQuadrilateral> gaussQuadratureReferenceQuadrilateral_5(5);
	static GaussQuadrature<ReferenceQuadrilateral> gaussQuadratureReferenceQuadrilateral_6(6);
	static GaussQuadrature<ReferenceQuadrilateral> gaussQuadratureReferenceQuadrilateral_7(7);
	static GaussQuadrature<ReferenceQuadrilateral> gaussQuadratureReferenceQuadrilateral_8(8);
	static GaussQuadrature<ReferenceQuadrilateral> gaussQuadratureReferenceQuadrilateral_9(9);
	static GaussQuadrature<ReferenceQuadrilateral> gaussQuadratureReferenceQuadrilateral_11(11);
	static GaussQuadrature<ReferenceQuadrilateral> gaussQuadratureReferenceQuadrilateral_13(13);

	bool success = true;
	success &= factory.register_rule(gaussQuadratureReferenceQuadrilateral_1);
	success &= factory.register_rule(gaussQuadratureReferenceQuadrilateral_2);
	success &= factory.register_rule(gaussQuadratureReferenceQuadrilateral_3);
	success &= factory.register_rule(gaussQuadratureReferenceQuadrilateral_4);
	success &= factory.register_rule(gaussQuadratureReferenceQuadrilateral_5);
	success &= factory.register_rule(gaussQuadratureReferenceQuadrilateral_6);
	success &= factory.register_rule(gaussQuadratureReferenceQuadrilateral_7);
	success &= factory.register_rule(gaussQuadratureReferenceQuadrilateral_8);
	success &= factory.register_rule(gaussQuadratureReferenceQuadrilateral_9);
	success &= factory.register_rule(gaussQuadratureReferenceQuadrilateral_11);
	success &= factory.register_rule(gaussQuadratureReferenceQuadrilateral_13);
	return success;
}


}; // namespace ug

 // register quadratures at factory
namespace {
using namespace ug;

template <>
std::vector<const QuadratureRule<ReferenceQuadrilateral>* > QuadratureRuleFactory<ReferenceQuadrilateral>::m_rules =
	std::vector<const QuadratureRule<ReferenceQuadrilateral>* >();

};
