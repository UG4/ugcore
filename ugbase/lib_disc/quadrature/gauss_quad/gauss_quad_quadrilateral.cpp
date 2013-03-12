//  This file is parsed from UG 3.9.
//  It provides the Gauss Quadratures for a reference quadrilateral.


#include "../quadrature.h"
#include "gauss_quad_quadrilateral.h"
#include "common/util/provider.h"

namespace ug{

GaussQuadrature<ReferenceQuadrilateral, 1>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.500000000000000000;
	m_vPoint[0][1] = 0.500000000000000000;

	m_vWeight[0] =  1.000000000000000000;
}

GaussQuadrature<ReferenceQuadrilateral, 2>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.211324865405187000;
	m_vPoint[0][1] = 0.211324865405187000;

	m_vPoint[1][0] = 0.788675134594813000;
	m_vPoint[1][1] = 0.211324865405187000;

	m_vPoint[2][0] = 0.211324865405187000;
	m_vPoint[2][1] = 0.788675134594813000;

	m_vPoint[3][0] = 0.788675134594813000;
	m_vPoint[3][1] = 0.788675134594813000;

	m_vWeight[0] =  0.250000000000000000;
	m_vWeight[1] =  0.250000000000000000;
	m_vWeight[2] =  0.250000000000000000;
	m_vWeight[3] =  0.250000000000000000;
}

GaussQuadrature<ReferenceQuadrilateral, 3>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.091751709536136984;
	m_vPoint[0][1] = 0.500000000000000000;

	m_vPoint[1][0] = 0.500000000000000000;
	m_vPoint[1][1] = 0.091751709536136984;

	m_vPoint[2][0] = 0.500000000000000000;
	m_vPoint[2][1] = 0.908248290463863016;

	m_vPoint[3][0] = 0.908248290463863016;
	m_vPoint[3][1] = 0.500000000000000000;

	m_vWeight[0] =  0.250000000000000000;
	m_vWeight[1] =  0.250000000000000000;
	m_vWeight[2] =  0.250000000000000000;
	m_vWeight[3] =  0.250000000000000000;
}

GaussQuadrature<ReferenceQuadrilateral, 4>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.500000000000000000;
	m_vPoint[0][1] = 0.500000000000000000;

	m_vPoint[1][0] = 0.983045891539647952;
	m_vPoint[1][1] = 0.500000000000000000;

	m_vPoint[2][0] = 0.727801863918096421;
	m_vPoint[2][1] = 0.074042673347699754;

	m_vPoint[3][0] = 0.727801863918096421;
	m_vPoint[3][1] = 0.925957326652300246;

	m_vPoint[4][0] = 0.134185024213432735;
	m_vPoint[4][1] = 0.184543605511622987;

	m_vPoint[5][0] = 0.134185024213432735;
	m_vPoint[5][1] = 0.815456394488377013;

	m_vWeight[0] =  0.285714285714285714;
	m_vWeight[1] =  0.109890109890109890;
	m_vWeight[2] =  0.141518051751883026;
	m_vWeight[3] =  0.141518051751883026;
	m_vWeight[4] =  0.160679750445919171;
	m_vWeight[5] =  0.160679750445919171;
}

GaussQuadrature<ReferenceQuadrilateral, 5>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.500000000000000000;
	m_vPoint[0][1] = 0.500000000000000000;

	m_vPoint[1][0] = 0.016954108460352048;
	m_vPoint[1][1] = 0.500000000000000000;

	m_vPoint[2][0] = 0.983045891539647952;
	m_vPoint[2][1] = 0.500000000000000000;

	m_vPoint[3][0] = 0.211324865405187118;
	m_vPoint[3][1] = 0.112701665379258311;

	m_vPoint[4][0] = 0.211324865405187118;
	m_vPoint[4][1] = 0.887298334620741689;

	m_vPoint[5][0] = 0.788675134594812882;
	m_vPoint[5][1] = 0.112701665379258311;

	m_vPoint[6][0] = 0.788675134594812882;
	m_vPoint[6][1] = 0.887298334620741689;

	m_vWeight[0] =  0.285714285714285714;
	m_vWeight[1] =  0.079365079365079365;
	m_vWeight[2] =  0.079365079365079365;
	m_vWeight[3] =  0.138888888888888889;
	m_vWeight[4] =  0.138888888888888889;
	m_vWeight[5] =  0.138888888888888889;
	m_vWeight[6] =  0.138888888888888889;
}

GaussQuadrature<ReferenceQuadrilateral, 6>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.918202816848812802;
	m_vPoint[0][1] = 0.500000000000000000;

	m_vPoint[1][0] = 0.321269917304346408;
	m_vPoint[1][1] = 0.500000000000000000;

	m_vPoint[2][0] = 0.936050765596565295;
	m_vPoint[2][1] = 0.055617992672617733;

	m_vPoint[3][0] = 0.936050765596565295;
	m_vPoint[3][1] = 0.944382007327382267;

	m_vPoint[4][0] = 0.652992581077713331;
	m_vPoint[4][1] = 0.197571180267657486;

	m_vPoint[5][0] = 0.652992581077713331;
	m_vPoint[5][1] = 0.802428819732342514;

	m_vPoint[6][0] = 0.294864550266671081;
	m_vPoint[6][1] = 0.022276246679468126;

	m_vPoint[7][0] = 0.294864550266671081;
	m_vPoint[7][1] = 0.977723753320531874;

	m_vPoint[8][0] = 0.063565344421560330;
	m_vPoint[8][1] = 0.217270003280622977;

	m_vPoint[9][0] = 0.063565344421560330;
	m_vPoint[9][1] = 0.782729996719377023;

	m_vWeight[0] =  0.113835811428543360;
	m_vWeight[1] =  0.206848993300741373;
	m_vWeight[2] =  0.036000221149911346;
	m_vWeight[3] =  0.036000221149911346;
	m_vWeight[4] =  0.167064776065666287;
	m_vWeight[5] =  0.167064776065666287;
	m_vWeight[6] =  0.056368501222669839;
	m_vWeight[7] =  0.056368501222669839;
	m_vWeight[8] =  0.080224099197110161;
	m_vWeight[9] =  0.080224099197110161;
}

GaussQuadrature<ReferenceQuadrilateral, 7>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.037089950113724269;
	m_vPoint[0][1] = 0.500000000000000000;

	m_vPoint[1][0] = 0.500000000000000000;
	m_vPoint[1][1] = 0.037089950113724269;

	m_vPoint[2][0] = 0.500000000000000000;
	m_vPoint[2][1] = 0.962910049886275731;

	m_vPoint[3][0] = 0.962910049886275731;
	m_vPoint[3][1] = 0.500000000000000000;

	m_vPoint[4][0] = 0.309722783395842172;
	m_vPoint[4][1] = 0.309722783395842172;

	m_vPoint[5][0] = 0.309722783395842172;
	m_vPoint[5][1] = 0.690277216604157828;

	m_vPoint[6][0] = 0.690277216604157828;
	m_vPoint[6][1] = 0.309722783395842172;

	m_vPoint[7][0] = 0.690277216604157828;
	m_vPoint[7][1] = 0.690277216604157828;

	m_vPoint[8][0] = 0.097010108540700628;
	m_vPoint[8][1] = 0.097010108540700628;

	m_vPoint[9][0] = 0.097010108540700628;
	m_vPoint[9][1] = 0.902989891459299372;

	m_vPoint[10][0] = 0.902989891459299372;
	m_vPoint[10][1] = 0.097010108540700628;

	m_vPoint[11][0] = 0.902989891459299372;
	m_vPoint[11][1] = 0.902989891459299372;

	m_vWeight[0] =  0.060493827160493827;
	m_vWeight[1] =  0.060493827160493827;
	m_vWeight[2] =  0.060493827160493827;
	m_vWeight[3] =  0.060493827160493827;
	m_vWeight[4] =  0.130148229166848614;
	m_vWeight[5] =  0.130148229166848614;
	m_vWeight[6] =  0.130148229166848614;
	m_vWeight[7] =  0.130148229166848614;
	m_vWeight[8] =  0.059357943672657559;
	m_vWeight[9] =  0.059357943672657559;
	m_vWeight[10] =  0.059357943672657559;
	m_vWeight[11] =  0.059357943672657559;
}

GaussQuadrature<ReferenceQuadrilateral, 8>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.500000000000000000;
	m_vPoint[0][1] = 0.500000000000000000;

	m_vPoint[1][0] = 0.878814588830252724;
	m_vPoint[1][1] = 0.500000000000000000;

	m_vPoint[2][0] = 0.381564078872149053;
	m_vPoint[2][1] = 0.500000000000000000;

	m_vPoint[3][0] = 0.005141035477736659;
	m_vPoint[3][1] = 0.500000000000000000;

	m_vPoint[4][0] = 0.975260477822833442;
	m_vPoint[4][1] = 0.180454347549815173;

	m_vPoint[5][0] = 0.975260477822833442;
	m_vPoint[5][1] = 0.819545652450184827;

	m_vPoint[6][0] = 0.831941368442816589;
	m_vPoint[6][1] = 0.031465461537504778;

	m_vPoint[7][0] = 0.831941368442816589;
	m_vPoint[7][1] = 0.968534538462495222;

	m_vPoint[8][0] = 0.652105340862052255;
	m_vPoint[8][1] = 0.231458234729253161;

	m_vPoint[9][0] = 0.652105340862052255;
	m_vPoint[9][1] = 0.768541765270746839;

	m_vPoint[10][0] = 0.381751640731940260;
	m_vPoint[10][1] = 0.056405746775187680;

	m_vPoint[11][0] = 0.381751640731940260;
	m_vPoint[11][1] = 0.943594253224812320;

	m_vPoint[12][0] = 0.150523261956718220;
	m_vPoint[12][1] = 0.252650589664901834;

	m_vPoint[13][0] = 0.150523261956718220;
	m_vPoint[13][1] = 0.747349410335098166;

	m_vPoint[14][0] = 0.049804612894210133;
	m_vPoint[14][1] = 0.051252090860116231;

	m_vPoint[15][0] = 0.049804612894210133;
	m_vPoint[15][1] = 0.948747909139883769;

	m_vWeight[0] =  0.013841176405359943;
	m_vWeight[1] =  0.101097342181518853;
	m_vWeight[2] =  0.133386651238158614;
	m_vWeight[3] =  0.029263547196684802;
	m_vWeight[4] =  0.031403604403436701;
	m_vWeight[5] =  0.031403604403436701;
	m_vWeight[6] =  0.034136146183397117;
	m_vWeight[7] =  0.034136146183397117;
	m_vWeight[8] =  0.120852119802814238;
	m_vWeight[9] =  0.120852119802814238;
	m_vWeight[10] =  0.063132126607385923;
	m_vWeight[11] =  0.063132126607385923;
	m_vWeight[12] =  0.090315580970543142;
	m_vWeight[13] =  0.090315580970543142;
	m_vWeight[14] =  0.021366063521561773;
	m_vWeight[15] =  0.021366063521561773;
}

GaussQuadrature<ReferenceQuadrilateral, 9>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.500000000000000000;
	m_vPoint[0][1] = 0.500000000000000000;

	m_vPoint[1][0] = 0.015575016819011140;
	m_vPoint[1][1] = 0.184659940134165573;

	m_vPoint[2][0] = 0.184659940134165573;
	m_vPoint[2][1] = 0.984424983180988860;

	m_vPoint[3][0] = 0.815340059865834427;
	m_vPoint[3][1] = 0.015575016819011140;

	m_vPoint[4][0] = 0.984424983180988860;
	m_vPoint[4][1] = 0.815340059865834427;

	m_vPoint[5][0] = 0.036019177020215166;
	m_vPoint[5][1] = 0.875138549989450267;

	m_vPoint[6][0] = 0.124861450010549733;
	m_vPoint[6][1] = 0.036019177020215166;

	m_vPoint[7][0] = 0.875138549989450267;
	m_vPoint[7][1] = 0.963980822979784834;

	m_vPoint[8][0] = 0.963980822979784834;
	m_vPoint[8][1] = 0.124861450010549733;

	m_vPoint[9][0] = 0.238132089892785332;
	m_vPoint[9][1] = 0.273330089432176405;

	m_vPoint[10][0] = 0.273330089432176405;
	m_vPoint[10][1] = 0.761867910107214668;

	m_vPoint[11][0] = 0.726669910567823595;
	m_vPoint[11][1] = 0.238132089892785332;

	m_vPoint[12][0] = 0.761867910107214668;
	m_vPoint[12][1] = 0.726669910567823595;

	m_vPoint[13][0] = 0.073692135333168846;
	m_vPoint[13][1] = 0.538104164096308587;

	m_vPoint[14][0] = 0.461895835903691413;
	m_vPoint[14][1] = 0.073692135333168846;

	m_vPoint[15][0] = 0.538104164096308587;
	m_vPoint[15][1] = 0.926307864666831154;

	m_vPoint[16][0] = 0.926307864666831154;
	m_vPoint[16][1] = 0.461895835903691413;

	m_vWeight[0] =  0.131687242798353909;
	m_vWeight[1] =  0.022219844542549677;
	m_vWeight[2] =  0.022219844542549677;
	m_vWeight[3] =  0.022219844542549677;
	m_vWeight[4] =  0.022219844542549677;
	m_vWeight[5] =  0.028024900532399121;
	m_vWeight[6] =  0.028024900532399121;
	m_vWeight[7] =  0.028024900532399121;
	m_vWeight[8] =  0.028024900532399121;
	m_vWeight[9] =  0.099570609815517524;
	m_vWeight[10] =  0.099570609815517524;
	m_vWeight[11] =  0.099570609815517524;
	m_vWeight[12] =  0.099570609815517524;
	m_vWeight[13] =  0.067262834409945201;
	m_vWeight[14] =  0.067262834409945201;
	m_vWeight[15] =  0.067262834409945201;
	m_vWeight[16] =  0.067262834409945201;
}

GaussQuadrature<ReferenceQuadrilateral, 11>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.008680388229572264;
	m_vPoint[0][1] = 0.150961947725216218;

	m_vPoint[1][0] = 0.150961947725216218;
	m_vPoint[1][1] = 0.991319611770427736;

	m_vPoint[2][0] = 0.849038052274783782;
	m_vPoint[2][1] = 0.008680388229572264;

	m_vPoint[3][0] = 0.991319611770427736;
	m_vPoint[3][1] = 0.849038052274783782;

	m_vPoint[4][0] = 0.030256808591631546;
	m_vPoint[4][1] = 0.912887917951481969;

	m_vPoint[5][0] = 0.087112082048518031;
	m_vPoint[5][1] = 0.030256808591631546;

	m_vPoint[6][0] = 0.912887917951481969;
	m_vPoint[6][1] = 0.969743191408368454;

	m_vPoint[7][0] = 0.969743191408368454;
	m_vPoint[7][1] = 0.087112082048518031;

	m_vPoint[8][0] = 0.023230235899233992;
	m_vPoint[8][1] = 0.594293069359320977;

	m_vPoint[9][0] = 0.405706930640679023;
	m_vPoint[9][1] = 0.023230235899233992;

	m_vPoint[10][0] = 0.594293069359320977;
	m_vPoint[10][1] = 0.976769764100766008;

	m_vPoint[11][0] = 0.976769764100766008;
	m_vPoint[11][1] = 0.405706930640679023;

	m_vPoint[12][0] = 0.093739725847593450;
	m_vPoint[12][1] = 0.342188283542372902;

	m_vPoint[13][0] = 0.342188283542372902;
	m_vPoint[13][1] = 0.906260274152406550;

	m_vPoint[14][0] = 0.657811716457627098;
	m_vPoint[14][1] = 0.093739725847593450;

	m_vPoint[15][0] = 0.906260274152406550;
	m_vPoint[15][1] = 0.657811716457627098;

	m_vPoint[16][0] = 0.143999043462331847;
	m_vPoint[16][1] = 0.762660125182273881;

	m_vPoint[17][0] = 0.237339874817726119;
	m_vPoint[17][1] = 0.143999043462331847;

	m_vPoint[18][0] = 0.762660125182273881;
	m_vPoint[18][1] = 0.856000956537668153;

	m_vPoint[19][0] = 0.856000956537668153;
	m_vPoint[19][1] = 0.237339874817726119;

	m_vPoint[20][0] = 0.287576375575665375;
	m_vPoint[20][1] = 0.520829035956011184;

	m_vPoint[21][0] = 0.479170964043988816;
	m_vPoint[21][1] = 0.287576375575665375;

	m_vPoint[22][0] = 0.520829035956011184;
	m_vPoint[22][1] = 0.712423624424334625;

	m_vPoint[23][0] = 0.712423624424334625;
	m_vPoint[23][1] = 0.479170964043988816;

	m_vWeight[0] =  0.012005190837680954;
	m_vWeight[1] =  0.012005190837680954;
	m_vWeight[2] =  0.012005190837680954;
	m_vWeight[3] =  0.012005190837680954;
	m_vWeight[4] =  0.016517832291137649;
	m_vWeight[5] =  0.016517832291137649;
	m_vWeight[6] =  0.016517832291137649;
	m_vWeight[7] =  0.016517832291137649;
	m_vWeight[8] =  0.024346694339667041;
	m_vWeight[9] =  0.024346694339667041;
	m_vWeight[10] =  0.024346694339667041;
	m_vWeight[11] =  0.024346694339667041;
	m_vWeight[12] =  0.052934087499737150;
	m_vWeight[13] =  0.052934087499737150;
	m_vWeight[14] =  0.052934087499737150;
	m_vWeight[15] =  0.052934087499737150;
	m_vWeight[16] =  0.056406515432215847;
	m_vWeight[17] =  0.056406515432215847;
	m_vWeight[18] =  0.056406515432215847;
	m_vWeight[19] =  0.056406515432215847;
	m_vWeight[20] =  0.087789679599561359;
	m_vWeight[21] =  0.087789679599561359;
	m_vWeight[22] =  0.087789679599561359;
	m_vWeight[23] =  0.087789679599561359;
}

GaussQuadrature<ReferenceQuadrilateral, 13>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.500000000000000000;
	m_vPoint[0][1] = 0.500000000000000000;

	m_vPoint[1][0] = 0.008256658780063868;
	m_vPoint[1][1] = 0.889404855777209711;

	m_vPoint[2][0] = 0.110595144222790289;
	m_vPoint[2][1] = 0.008256658780063868;

	m_vPoint[3][0] = 0.889404855777209711;
	m_vPoint[3][1] = 0.991743341219936132;

	m_vPoint[4][0] = 0.991743341219936132;
	m_vPoint[4][1] = 0.110595144222790289;

	m_vPoint[5][0] = 0.021351150106846317;
	m_vPoint[5][1] = 0.070221997179180536;

	m_vPoint[6][0] = 0.070221997179180536;
	m_vPoint[6][1] = 0.978648849893153683;

	m_vPoint[7][0] = 0.929778002820819464;
	m_vPoint[7][1] = 0.021351150106846317;

	m_vPoint[8][0] = 0.978648849893153683;
	m_vPoint[8][1] = 0.929778002820819464;

	m_vPoint[9][0] = 0.020537414856232571;
	m_vPoint[9][1] = 0.569091729931232677;

	m_vPoint[10][0] = 0.430908270068767323;
	m_vPoint[10][1] = 0.020537414856232571;

	m_vPoint[11][0] = 0.569091729931232677;
	m_vPoint[11][1] = 0.979462585143767429;

	m_vPoint[12][0] = 0.979462585143767429;
	m_vPoint[12][1] = 0.430908270068767323;

	m_vPoint[13][0] = 0.029336387063537382;
	m_vPoint[13][1] = 0.304631891935269500;

	m_vPoint[14][0] = 0.304631891935269500;
	m_vPoint[14][1] = 0.970663612936462618;

	m_vPoint[15][0] = 0.695368108064730500;
	m_vPoint[15][1] = 0.029336387063537382;

	m_vPoint[16][0] = 0.970663612936462618;
	m_vPoint[16][1] = 0.695368108064730500;

	m_vPoint[17][0] = 0.074961663150125712;
	m_vPoint[17][1] = 0.737904312609137953;

	m_vPoint[18][0] = 0.262095687390862047;
	m_vPoint[18][1] = 0.074961663150125712;

	m_vPoint[19][0] = 0.737904312609137953;
	m_vPoint[19][1] = 0.925038336849874288;

	m_vPoint[20][0] = 0.925038336849874288;
	m_vPoint[20][1] = 0.262095687390862047;

	m_vPoint[21][0] = 0.122097321713959282;
	m_vPoint[21][1] = 0.176089181406494634;

	m_vPoint[22][0] = 0.176089181406494634;
	m_vPoint[22][1] = 0.877902678286040718;

	m_vPoint[23][0] = 0.823910818593505366;
	m_vPoint[23][1] = 0.122097321713959282;

	m_vPoint[24][0] = 0.877902678286040718;
	m_vPoint[24][1] = 0.823910818593505366;

	m_vPoint[25][0] = 0.151874960754125293;
	m_vPoint[25][1] = 0.464629245501777532;

	m_vPoint[26][0] = 0.464629245501777532;
	m_vPoint[26][1] = 0.848125039245874707;

	m_vPoint[27][0] = 0.535370754498222468;
	m_vPoint[27][1] = 0.151874960754125293;

	m_vPoint[28][0] = 0.848125039245874707;
	m_vPoint[28][1] = 0.535370754498222468;

	m_vPoint[29][0] = 0.295347719152980578;
	m_vPoint[29][1] = 0.671358278020203395;

	m_vPoint[30][0] = 0.328641721979796605;
	m_vPoint[30][1] = 0.295347719152980578;

	m_vPoint[31][0] = 0.671358278020203395;
	m_vPoint[31][1] = 0.704652280847019422;

	m_vPoint[32][0] = 0.704652280847019422;
	m_vPoint[32][1] = 0.328641721979796605;

	m_vWeight[0] =  0.075095528857806340;
	m_vWeight[1] =  0.007497959716124783;
	m_vWeight[2] =  0.007497959716124783;
	m_vWeight[3] =  0.007497959716124783;
	m_vWeight[4] =  0.007497959716124783;
	m_vWeight[5] =  0.009543605329270917;
	m_vWeight[6] =  0.009543605329270917;
	m_vWeight[7] =  0.009543605329270917;
	m_vWeight[8] =  0.009543605329270917;
	m_vWeight[9] =  0.015106230954437495;
	m_vWeight[10] =  0.015106230954437495;
	m_vWeight[11] =  0.015106230954437495;
	m_vWeight[12] =  0.015106230954437495;
	m_vWeight[13] =  0.019373184633276335;
	m_vWeight[14] =  0.019373184633276335;
	m_vWeight[15] =  0.019373184633276335;
	m_vWeight[16] =  0.019373184633276335;
	m_vWeight[17] =  0.029711166825148900;
	m_vWeight[18] =  0.029711166825148900;
	m_vWeight[19] =  0.029711166825148900;
	m_vWeight[20] =  0.029711166825148900;
	m_vWeight[21] =  0.032440887592500678;
	m_vWeight[22] =  0.032440887592500678;
	m_vWeight[23] =  0.032440887592500678;
	m_vWeight[24] =  0.032440887592500678;
	m_vWeight[25] =  0.053335395364297347;
	m_vWeight[26] =  0.053335395364297347;
	m_vWeight[27] =  0.053335395364297347;
	m_vWeight[28] =  0.053335395364297347;
	m_vWeight[29] =  0.064217687370491959;
	m_vWeight[30] =  0.064217687370491959;
	m_vWeight[31] =  0.064217687370491959;
	m_vWeight[32] =  0.064217687370491959;
}




template <>
FlexGaussQuadrature<ReferenceQuadrilateral>::FlexGaussQuadrature(int order)
{
	switch(order)
	{
	case 1:{
		const static GaussQuadrature<ReferenceQuadrilateral, 1>& q1 
			= Provider<GaussQuadrature<ReferenceQuadrilateral, 1> >::get();

		m_order = q1.order();
		m_numPoints = q1.size();
		m_pvPoint = q1.points();
		m_pvWeight = q1.weights();
		}break;

	case 2:{
		const static GaussQuadrature<ReferenceQuadrilateral, 2>& q2 
			= Provider<GaussQuadrature<ReferenceQuadrilateral, 2> >::get();

		m_order = q2.order();
		m_numPoints = q2.size();
		m_pvPoint = q2.points();
		m_pvWeight = q2.weights();
		}break;

	case 3:{
		const static GaussQuadrature<ReferenceQuadrilateral, 3>& q3 
			= Provider<GaussQuadrature<ReferenceQuadrilateral, 3> >::get();

		m_order = q3.order();
		m_numPoints = q3.size();
		m_pvPoint = q3.points();
		m_pvWeight = q3.weights();
		}break;

	case 4:{
		const static GaussQuadrature<ReferenceQuadrilateral, 4>& q4 
			= Provider<GaussQuadrature<ReferenceQuadrilateral, 4> >::get();

		m_order = q4.order();
		m_numPoints = q4.size();
		m_pvPoint = q4.points();
		m_pvWeight = q4.weights();
		}break;

	case 5:{
		const static GaussQuadrature<ReferenceQuadrilateral, 5>& q5 
			= Provider<GaussQuadrature<ReferenceQuadrilateral, 5> >::get();

		m_order = q5.order();
		m_numPoints = q5.size();
		m_pvPoint = q5.points();
		m_pvWeight = q5.weights();
		}break;

	case 6:{
		const static GaussQuadrature<ReferenceQuadrilateral, 6>& q6 
			= Provider<GaussQuadrature<ReferenceQuadrilateral, 6> >::get();

		m_order = q6.order();
		m_numPoints = q6.size();
		m_pvPoint = q6.points();
		m_pvWeight = q6.weights();
		}break;

	case 7:{
		const static GaussQuadrature<ReferenceQuadrilateral, 7>& q7 
			= Provider<GaussQuadrature<ReferenceQuadrilateral, 7> >::get();

		m_order = q7.order();
		m_numPoints = q7.size();
		m_pvPoint = q7.points();
		m_pvWeight = q7.weights();
		}break;

	case 8:{
		const static GaussQuadrature<ReferenceQuadrilateral, 8>& q8 
			= Provider<GaussQuadrature<ReferenceQuadrilateral, 8> >::get();

		m_order = q8.order();
		m_numPoints = q8.size();
		m_pvPoint = q8.points();
		m_pvWeight = q8.weights();
		}break;

	case 9:{
		const static GaussQuadrature<ReferenceQuadrilateral, 9>& q9 
			= Provider<GaussQuadrature<ReferenceQuadrilateral, 9> >::get();

		m_order = q9.order();
		m_numPoints = q9.size();
		m_pvPoint = q9.points();
		m_pvWeight = q9.weights();
		}break;

	case 11:{
		const static GaussQuadrature<ReferenceQuadrilateral, 11>& q11 
			= Provider<GaussQuadrature<ReferenceQuadrilateral, 11> >::get();

		m_order = q11.order();
		m_numPoints = q11.size();
		m_pvPoint = q11.points();
		m_pvWeight = q11.weights();
		}break;

	case 13:{
		const static GaussQuadrature<ReferenceQuadrilateral, 13>& q13 
			= Provider<GaussQuadrature<ReferenceQuadrilateral, 13> >::get();

		m_order = q13.order();
		m_numPoints = q13.size();
		m_pvPoint = q13.points();
		m_pvWeight = q13.weights();
		}break;

	default: UG_ASSERT(0, "Order not availabile. Can not construct GaussQuadrature.\n");
	}
}



// register rules
template <>
bool RegisterGaussQuadRule<ReferenceQuadrilateral>(QuadratureRuleProvider<ReferenceQuadrilateral::dim>& factory)
{
	static FlexGaussQuadrature<ReferenceQuadrilateral> gaussQuadratureReferenceQuadrilateral_1(1);
	static FlexGaussQuadrature<ReferenceQuadrilateral> gaussQuadratureReferenceQuadrilateral_2(2);
	static FlexGaussQuadrature<ReferenceQuadrilateral> gaussQuadratureReferenceQuadrilateral_3(3);
	static FlexGaussQuadrature<ReferenceQuadrilateral> gaussQuadratureReferenceQuadrilateral_4(4);
	static FlexGaussQuadrature<ReferenceQuadrilateral> gaussQuadratureReferenceQuadrilateral_5(5);
	static FlexGaussQuadrature<ReferenceQuadrilateral> gaussQuadratureReferenceQuadrilateral_6(6);
	static FlexGaussQuadrature<ReferenceQuadrilateral> gaussQuadratureReferenceQuadrilateral_7(7);
	static FlexGaussQuadrature<ReferenceQuadrilateral> gaussQuadratureReferenceQuadrilateral_8(8);
	static FlexGaussQuadrature<ReferenceQuadrilateral> gaussQuadratureReferenceQuadrilateral_9(9);
	static FlexGaussQuadrature<ReferenceQuadrilateral> gaussQuadratureReferenceQuadrilateral_11(11);
	static FlexGaussQuadrature<ReferenceQuadrilateral> gaussQuadratureReferenceQuadrilateral_13(13);

	bool success = true;
	factory.register_rule<ReferenceQuadrilateral>(gaussQuadratureReferenceQuadrilateral_1);
	factory.register_rule<ReferenceQuadrilateral>(gaussQuadratureReferenceQuadrilateral_2);
	factory.register_rule<ReferenceQuadrilateral>(gaussQuadratureReferenceQuadrilateral_3);
	factory.register_rule<ReferenceQuadrilateral>(gaussQuadratureReferenceQuadrilateral_4);
	factory.register_rule<ReferenceQuadrilateral>(gaussQuadratureReferenceQuadrilateral_5);
	factory.register_rule<ReferenceQuadrilateral>(gaussQuadratureReferenceQuadrilateral_6);
	factory.register_rule<ReferenceQuadrilateral>(gaussQuadratureReferenceQuadrilateral_7);
	factory.register_rule<ReferenceQuadrilateral>(gaussQuadratureReferenceQuadrilateral_8);
	factory.register_rule<ReferenceQuadrilateral>(gaussQuadratureReferenceQuadrilateral_9);
	factory.register_rule<ReferenceQuadrilateral>(gaussQuadratureReferenceQuadrilateral_11);
	factory.register_rule<ReferenceQuadrilateral>(gaussQuadratureReferenceQuadrilateral_13);

	return success;
};

}; // namespace ug

