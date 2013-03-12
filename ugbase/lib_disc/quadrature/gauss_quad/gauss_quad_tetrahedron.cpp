//  This file is parsed from UG 3.9.
//  It provides the Gauss Quadratures for a reference tetrahedron.


#include "../quadrature.h"
#include "gauss_quad_tetrahedron.h"
#include "common/util/provider.h"

namespace ug{

GaussQuadrature<ReferenceTetrahedron, 0>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.250000000000000000;
	m_vPoint[0][1] = 0.250000000000000000;
	m_vPoint[0][2] = 0.250000000000000000;

	m_vPoint[1][0] = 0.100000000000000000;
	m_vPoint[1][1] = 0.100000000000000000;
	m_vPoint[1][2] = 0.400000000000000000;

	m_vPoint[2][0] = 0.100000000000000000;
	m_vPoint[2][1] = 0.400000000000000000;
	m_vPoint[2][2] = 0.100000000000000000;

	m_vPoint[3][0] = 0.400000000000000000;
	m_vPoint[3][1] = 0.100000000000000000;
	m_vPoint[3][2] = 0.100000000000000000;

	m_vPoint[4][0] = 0.100000000000000000;
	m_vPoint[4][1] = 0.400000000000000000;
	m_vPoint[4][2] = 0.400000000000000000;

	m_vPoint[5][0] = 0.400000000000000000;
	m_vPoint[5][1] = 0.100000000000000000;
	m_vPoint[5][2] = 0.400000000000000000;

	m_vPoint[6][0] = 0.400000000000000000;
	m_vPoint[6][1] = 0.400000000000000000;
	m_vPoint[6][2] = 0.100000000000000000;

	m_vWeight[0] = 1./6. * 0.142857142857143000;
	m_vWeight[1] = 1./6. * 0.142857142857143000;
	m_vWeight[2] = 1./6. * 0.142857142857143000;
	m_vWeight[3] = 1./6. * 0.142857142857143000;
	m_vWeight[4] = 1./6. * 0.142857142857143000;
	m_vWeight[5] = 1./6. * 0.142857142857143000;
	m_vWeight[6] = 1./6. * 0.142857142857143000;
}

GaussQuadrature<ReferenceTetrahedron, 1>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.250000000000000000;
	m_vPoint[0][1] = 0.250000000000000000;
	m_vPoint[0][2] = 0.250000000000000000;

	m_vWeight[0] = 1./6. * 1.000000000000000000;
}

GaussQuadrature<ReferenceTetrahedron, 2>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.585410200000000000;
	m_vPoint[0][1] = 0.138196600000000000;
	m_vPoint[0][2] = 0.138196600000000000;

	m_vPoint[1][0] = 0.138196600000000000;
	m_vPoint[1][1] = 0.585410200000000000;
	m_vPoint[1][2] = 0.138196600000000000;

	m_vPoint[2][0] = 0.138196600000000000;
	m_vPoint[2][1] = 0.138196600000000000;
	m_vPoint[2][2] = 0.585410200000000000;

	m_vPoint[3][0] = 0.138196600000000000;
	m_vPoint[3][1] = 0.138196600000000000;
	m_vPoint[3][2] = 0.138196600000000000;

	m_vWeight[0] = 1./6. * 0.250000000000000000;
	m_vWeight[1] = 1./6. * 0.250000000000000000;
	m_vWeight[2] = 1./6. * 0.250000000000000000;
	m_vWeight[3] = 1./6. * 0.250000000000000000;
}

GaussQuadrature<ReferenceTetrahedron, 3>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.000000000000000000;
	m_vPoint[0][1] = 0.000000000000000000;
	m_vPoint[0][2] = 0.000000000000000000;

	m_vPoint[1][0] = 1.000000000000000000;
	m_vPoint[1][1] = 0.000000000000000000;
	m_vPoint[1][2] = 0.000000000000000000;

	m_vPoint[2][0] = 0.000000000000000000;
	m_vPoint[2][1] = 1.000000000000000000;
	m_vPoint[2][2] = 0.000000000000000000;

	m_vPoint[3][0] = 0.000000000000000000;
	m_vPoint[3][1] = 0.000000000000000000;
	m_vPoint[3][2] = 1.000000000000000000;

	m_vPoint[4][0] = 0.333333333333000000;
	m_vPoint[4][1] = 0.333333333333000000;
	m_vPoint[4][2] = 0.000000000000000000;

	m_vPoint[5][0] = 0.333333333333000000;
	m_vPoint[5][1] = 0.000000000000000000;
	m_vPoint[5][2] = 0.333333333333000000;

	m_vPoint[6][0] = 0.000000000000000000;
	m_vPoint[6][1] = 0.333333333333000000;
	m_vPoint[6][2] = 0.333333333333000000;

	m_vPoint[7][0] = 0.333333333333000000;
	m_vPoint[7][1] = 0.333333333333000000;
	m_vPoint[7][2] = 0.333333333333000000;

	m_vWeight[0] = 1./6. * 0.025000000000000000;
	m_vWeight[1] = 1./6. * 0.025000000000000000;
	m_vWeight[2] = 1./6. * 0.025000000000000000;
	m_vWeight[3] = 1./6. * 0.025000000000000000;
	m_vWeight[4] = 1./6. * 0.225000000000000000;
	m_vWeight[5] = 1./6. * 0.225000000000000000;
	m_vWeight[6] = 1./6. * 0.225000000000000000;
	m_vWeight[7] = 1./6. * 0.225000000000000000;
}

GaussQuadrature<ReferenceTetrahedron, 5>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.250000000000000000;
	m_vPoint[0][1] = 0.250000000000000000;
	m_vPoint[0][2] = 0.250000000000000000;

	m_vPoint[1][0] = 0.091971078000000000;
	m_vPoint[1][1] = 0.091971078000000000;
	m_vPoint[1][2] = 0.091971078000000000;

	m_vPoint[2][0] = 0.724086770000000000;
	m_vPoint[2][1] = 0.091971078000000000;
	m_vPoint[2][2] = 0.091971078000000000;

	m_vPoint[3][0] = 0.091971078000000000;
	m_vPoint[3][1] = 0.724086770000000000;
	m_vPoint[3][2] = 0.091971078000000000;

	m_vPoint[4][0] = 0.091971078000000000;
	m_vPoint[4][1] = 0.091971078000000000;
	m_vPoint[4][2] = 0.724086770000000000;

	m_vPoint[5][0] = 0.319793630000000000;
	m_vPoint[5][1] = 0.319793630000000000;
	m_vPoint[5][2] = 0.319793630000000000;

	m_vPoint[6][0] = 0.040619117000000000;
	m_vPoint[6][1] = 0.319793630000000000;
	m_vPoint[6][2] = 0.319793630000000000;

	m_vPoint[7][0] = 0.319793630000000000;
	m_vPoint[7][1] = 0.040619117000000000;
	m_vPoint[7][2] = 0.319793630000000000;

	m_vPoint[8][0] = 0.319793630000000000;
	m_vPoint[8][1] = 0.319793630000000000;
	m_vPoint[8][2] = 0.040619117000000000;

	m_vPoint[9][0] = 0.443649170000000000;
	m_vPoint[9][1] = 0.056350833000000000;
	m_vPoint[9][2] = 0.056350833000000000;

	m_vPoint[10][0] = 0.056350833000000000;
	m_vPoint[10][1] = 0.443649170000000000;
	m_vPoint[10][2] = 0.056350833000000000;

	m_vPoint[11][0] = 0.056350833000000000;
	m_vPoint[11][1] = 0.056350833000000000;
	m_vPoint[11][2] = 0.443649170000000000;

	m_vPoint[12][0] = 0.443649170000000000;
	m_vPoint[12][1] = 0.443649170000000000;
	m_vPoint[12][2] = 0.056350833000000000;

	m_vPoint[13][0] = 0.443649170000000000;
	m_vPoint[13][1] = 0.056350833000000000;
	m_vPoint[13][2] = 0.443649170000000000;

	m_vPoint[14][0] = 0.056350833000000000;
	m_vPoint[14][1] = 0.443649170000000000;
	m_vPoint[14][2] = 0.443649170000000000;

	m_vWeight[0] = 1./6. * 0.118518520000000000;
	m_vWeight[1] = 1./6. * 0.071937084000000000;
	m_vWeight[2] = 1./6. * 0.071937084000000000;
	m_vWeight[3] = 1./6. * 0.071937084000000000;
	m_vWeight[4] = 1./6. * 0.071937084000000000;
	m_vWeight[5] = 1./6. * 0.069068207000000000;
	m_vWeight[6] = 1./6. * 0.069068207000000000;
	m_vWeight[7] = 1./6. * 0.069068207000000000;
	m_vWeight[8] = 1./6. * 0.069068207000000000;
	m_vWeight[9] = 1./6. * 0.052910053000000000;
	m_vWeight[10] = 1./6. * 0.052910053000000000;
	m_vWeight[11] = 1./6. * 0.052910053000000000;
	m_vWeight[12] = 1./6. * 0.052910053000000000;
	m_vWeight[13] = 1./6. * 0.052910053000000000;
	m_vWeight[14] = 1./6. * 0.052910053000000000;
}

GaussQuadrature<ReferenceTetrahedron, 6>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.214602871259152029;
	m_vPoint[0][1] = 0.214602871259152029;
	m_vPoint[0][2] = 0.214602871259152029;

	m_vPoint[1][0] = 0.356191386222543912;
	m_vPoint[1][1] = 0.214602871259152029;
	m_vPoint[1][2] = 0.214602871259152029;

	m_vPoint[2][0] = 0.214602871259152029;
	m_vPoint[2][1] = 0.356191386222543912;
	m_vPoint[2][2] = 0.214602871259152029;

	m_vPoint[3][0] = 0.214602871259152029;
	m_vPoint[3][1] = 0.214602871259152029;
	m_vPoint[3][2] = 0.356191386222543912;

	m_vPoint[4][0] = 0.040673958534611353;
	m_vPoint[4][1] = 0.040673958534611353;
	m_vPoint[4][2] = 0.040673958534611353;

	m_vPoint[5][0] = 0.877978124396165941;
	m_vPoint[5][1] = 0.040673958534611353;
	m_vPoint[5][2] = 0.040673958534611353;

	m_vPoint[6][0] = 0.040673958534611353;
	m_vPoint[6][1] = 0.877978124396165941;
	m_vPoint[6][2] = 0.040673958534611353;

	m_vPoint[7][0] = 0.040673958534611353;
	m_vPoint[7][1] = 0.040673958534611353;
	m_vPoint[7][2] = 0.877978124396165941;

	m_vPoint[8][0] = 0.322337890142275510;
	m_vPoint[8][1] = 0.322337890142275510;
	m_vPoint[8][2] = 0.322337890142275510;

	m_vPoint[9][0] = 0.032986329573173469;
	m_vPoint[9][1] = 0.322337890142275510;
	m_vPoint[9][2] = 0.322337890142275510;

	m_vPoint[10][0] = 0.322337890142275510;
	m_vPoint[10][1] = 0.032986329573173469;
	m_vPoint[10][2] = 0.322337890142275510;

	m_vPoint[11][0] = 0.322337890142275510;
	m_vPoint[11][1] = 0.322337890142275510;
	m_vPoint[11][2] = 0.032986329573173469;

	m_vPoint[12][0] = 0.063661001875017525;
	m_vPoint[12][1] = 0.063661001875017525;
	m_vPoint[12][2] = 0.269672331458315808;

	m_vPoint[13][0] = 0.063661001875017525;
	m_vPoint[13][1] = 0.269672331458315808;
	m_vPoint[13][2] = 0.063661001875017525;

	m_vPoint[14][0] = 0.269672331458315808;
	m_vPoint[14][1] = 0.063661001875017525;
	m_vPoint[14][2] = 0.063661001875017525;

	m_vPoint[15][0] = 0.063661001875017525;
	m_vPoint[15][1] = 0.063661001875017525;
	m_vPoint[15][2] = 0.603005664791649141;

	m_vPoint[16][0] = 0.063661001875017525;
	m_vPoint[16][1] = 0.603005664791649141;
	m_vPoint[16][2] = 0.063661001875017525;

	m_vPoint[17][0] = 0.603005664791649141;
	m_vPoint[17][1] = 0.063661001875017525;
	m_vPoint[17][2] = 0.063661001875017525;

	m_vPoint[18][0] = 0.063661001875017525;
	m_vPoint[18][1] = 0.269672331458315808;
	m_vPoint[18][2] = 0.603005664791649141;

	m_vPoint[19][0] = 0.063661001875017525;
	m_vPoint[19][1] = 0.603005664791649141;
	m_vPoint[19][2] = 0.269672331458315808;

	m_vPoint[20][0] = 0.269672331458315808;
	m_vPoint[20][1] = 0.063661001875017525;
	m_vPoint[20][2] = 0.603005664791649141;

	m_vPoint[21][0] = 0.603005664791649141;
	m_vPoint[21][1] = 0.063661001875017525;
	m_vPoint[21][2] = 0.269672331458315808;

	m_vPoint[22][0] = 0.269672331458315808;
	m_vPoint[22][1] = 0.603005664791649141;
	m_vPoint[22][2] = 0.063661001875017525;

	m_vPoint[23][0] = 0.603005664791649141;
	m_vPoint[23][1] = 0.269672331458315808;
	m_vPoint[23][2] = 0.063661001875017525;

	m_vWeight[0] = 1./6. * 0.039922750258167492;
	m_vWeight[1] = 1./6. * 0.039922750258167492;
	m_vWeight[2] = 1./6. * 0.039922750258167492;
	m_vWeight[3] = 1./6. * 0.039922750258167492;
	m_vWeight[4] = 1./6. * 0.010077211055320643;
	m_vWeight[5] = 1./6. * 0.010077211055320643;
	m_vWeight[6] = 1./6. * 0.010077211055320643;
	m_vWeight[7] = 1./6. * 0.010077211055320643;
	m_vWeight[8] = 1./6. * 0.055357181543654722;
	m_vWeight[9] = 1./6. * 0.055357181543654722;
	m_vWeight[10] = 1./6. * 0.055357181543654722;
	m_vWeight[11] = 1./6. * 0.055357181543654722;
	m_vWeight[12] = 1./6. * 0.048214285714285714;
	m_vWeight[13] = 1./6. * 0.048214285714285714;
	m_vWeight[14] = 1./6. * 0.048214285714285714;
	m_vWeight[15] = 1./6. * 0.048214285714285714;
	m_vWeight[16] = 1./6. * 0.048214285714285714;
	m_vWeight[17] = 1./6. * 0.048214285714285714;
	m_vWeight[18] = 1./6. * 0.048214285714285714;
	m_vWeight[19] = 1./6. * 0.048214285714285714;
	m_vWeight[20] = 1./6. * 0.048214285714285714;
	m_vWeight[21] = 1./6. * 0.048214285714285714;
	m_vWeight[22] = 1./6. * 0.048214285714285714;
	m_vWeight[23] = 1./6. * 0.048214285714285714;
}

GaussQuadrature<ReferenceTetrahedron, 7>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.500000000000000000;
	m_vPoint[0][1] = 0.500000000000000000;
	m_vPoint[0][2] = 0.000000000000000000;

	m_vPoint[1][0] = 0.500000000000000000;
	m_vPoint[1][1] = 0.000000000000000000;
	m_vPoint[1][2] = 0.500000000000000000;

	m_vPoint[2][0] = 0.000000000000000000;
	m_vPoint[2][1] = 0.500000000000000000;
	m_vPoint[2][2] = 0.500000000000000000;

	m_vPoint[3][0] = 0.000000000000000000;
	m_vPoint[3][1] = 0.000000000000000000;
	m_vPoint[3][2] = 0.500000000000000000;

	m_vPoint[4][0] = 0.000000000000000000;
	m_vPoint[4][1] = 0.500000000000000000;
	m_vPoint[4][2] = 0.000000000000000000;

	m_vPoint[5][0] = 0.500000000000000000;
	m_vPoint[5][1] = 0.000000000000000000;
	m_vPoint[5][2] = 0.000000000000000000;

	m_vPoint[6][0] = 0.250000000000000000;
	m_vPoint[6][1] = 0.250000000000000000;
	m_vPoint[6][2] = 0.250000000000000000;

	m_vPoint[7][0] = 0.078213192330318064;
	m_vPoint[7][1] = 0.078213192330318064;
	m_vPoint[7][2] = 0.078213192330318064;

	m_vPoint[8][0] = 0.765360423009045807;
	m_vPoint[8][1] = 0.078213192330318064;
	m_vPoint[8][2] = 0.078213192330318064;

	m_vPoint[9][0] = 0.078213192330318064;
	m_vPoint[9][1] = 0.765360423009045807;
	m_vPoint[9][2] = 0.078213192330318064;

	m_vPoint[10][0] = 0.078213192330318064;
	m_vPoint[10][1] = 0.078213192330318064;
	m_vPoint[10][2] = 0.765360423009045807;

	m_vPoint[11][0] = 0.121843216663905175;
	m_vPoint[11][1] = 0.121843216663905175;
	m_vPoint[11][2] = 0.121843216663905175;

	m_vPoint[12][0] = 0.634470350008284476;
	m_vPoint[12][1] = 0.121843216663905175;
	m_vPoint[12][2] = 0.121843216663905175;

	m_vPoint[13][0] = 0.121843216663905175;
	m_vPoint[13][1] = 0.634470350008284476;
	m_vPoint[13][2] = 0.121843216663905175;

	m_vPoint[14][0] = 0.121843216663905175;
	m_vPoint[14][1] = 0.121843216663905175;
	m_vPoint[14][2] = 0.634470350008284476;

	m_vPoint[15][0] = 0.332539164446420624;
	m_vPoint[15][1] = 0.332539164446420624;
	m_vPoint[15][2] = 0.332539164446420624;

	m_vPoint[16][0] = 0.002382506660738128;
	m_vPoint[16][1] = 0.332539164446420624;
	m_vPoint[16][2] = 0.332539164446420624;

	m_vPoint[17][0] = 0.332539164446420624;
	m_vPoint[17][1] = 0.002382506660738128;
	m_vPoint[17][2] = 0.332539164446420624;

	m_vPoint[18][0] = 0.332539164446420624;
	m_vPoint[18][1] = 0.332539164446420624;
	m_vPoint[18][2] = 0.002382506660738128;

	m_vPoint[19][0] = 0.100000000000000000;
	m_vPoint[19][1] = 0.100000000000000000;
	m_vPoint[19][2] = 0.200000000000000000;

	m_vPoint[20][0] = 0.100000000000000000;
	m_vPoint[20][1] = 0.200000000000000000;
	m_vPoint[20][2] = 0.100000000000000000;

	m_vPoint[21][0] = 0.200000000000000000;
	m_vPoint[21][1] = 0.100000000000000000;
	m_vPoint[21][2] = 0.100000000000000000;

	m_vPoint[22][0] = 0.100000000000000000;
	m_vPoint[22][1] = 0.100000000000000000;
	m_vPoint[22][2] = 0.600000000000000000;

	m_vPoint[23][0] = 0.100000000000000000;
	m_vPoint[23][1] = 0.600000000000000000;
	m_vPoint[23][2] = 0.100000000000000000;

	m_vPoint[24][0] = 0.600000000000000000;
	m_vPoint[24][1] = 0.100000000000000000;
	m_vPoint[24][2] = 0.100000000000000000;

	m_vPoint[25][0] = 0.100000000000000000;
	m_vPoint[25][1] = 0.200000000000000000;
	m_vPoint[25][2] = 0.600000000000000000;

	m_vPoint[26][0] = 0.100000000000000000;
	m_vPoint[26][1] = 0.600000000000000000;
	m_vPoint[26][2] = 0.200000000000000000;

	m_vPoint[27][0] = 0.600000000000000000;
	m_vPoint[27][1] = 0.100000000000000000;
	m_vPoint[27][2] = 0.200000000000000000;

	m_vPoint[28][0] = 0.200000000000000000;
	m_vPoint[28][1] = 0.100000000000000000;
	m_vPoint[28][2] = 0.600000000000000000;

	m_vPoint[29][0] = 0.200000000000000000;
	m_vPoint[29][1] = 0.600000000000000000;
	m_vPoint[29][2] = 0.100000000000000000;

	m_vPoint[30][0] = 0.600000000000000000;
	m_vPoint[30][1] = 0.200000000000000000;
	m_vPoint[30][2] = 0.100000000000000000;

	m_vWeight[0] = 1./6. * 0.005820105820105820;
	m_vWeight[1] = 1./6. * 0.005820105820105820;
	m_vWeight[2] = 1./6. * 0.005820105820105820;
	m_vWeight[3] = 1./6. * 0.005820105820105820;
	m_vWeight[4] = 1./6. * 0.005820105820105820;
	m_vWeight[5] = 1./6. * 0.005820105820105820;
	m_vWeight[6] = 1./6. * 0.109585340796652922;
	m_vWeight[7] = 1./6. * 0.063599649146482121;
	m_vWeight[8] = 1./6. * 0.063599649146482121;
	m_vWeight[9] = 1./6. * 0.063599649146482121;
	m_vWeight[10] = 1./6. * 0.063599649146482121;
	m_vWeight[11] = 1./6. * -0.375106440685991110;
	m_vWeight[12] = 1./6. * -0.375106440685991110;
	m_vWeight[13] = 1./6. * -0.375106440685991110;
	m_vWeight[14] = 1./6. * -0.375106440685991110;
	m_vWeight[15] = 1./6. * 0.029348551578440996;
	m_vWeight[16] = 1./6. * 0.029348551578440996;
	m_vWeight[17] = 1./6. * 0.029348551578440996;
	m_vWeight[18] = 1./6. * 0.029348551578440996;
	m_vWeight[19] = 1./6. * 0.165343915343915344;
	m_vWeight[20] = 1./6. * 0.165343915343915344;
	m_vWeight[21] = 1./6. * 0.165343915343915344;
	m_vWeight[22] = 1./6. * 0.165343915343915344;
	m_vWeight[23] = 1./6. * 0.165343915343915344;
	m_vWeight[24] = 1./6. * 0.165343915343915344;
	m_vWeight[25] = 1./6. * 0.165343915343915344;
	m_vWeight[26] = 1./6. * 0.165343915343915344;
	m_vWeight[27] = 1./6. * 0.165343915343915344;
	m_vWeight[28] = 1./6. * 0.165343915343915344;
	m_vWeight[29] = 1./6. * 0.165343915343915344;
	m_vWeight[30] = 1./6. * 0.165343915343915344;
}

GaussQuadrature<ReferenceTetrahedron, 8>::GaussQuadrature()
{
	m_vPoint[0][0] = 0.250000000000000000;
	m_vPoint[0][1] = 0.250000000000000000;
	m_vPoint[0][2] = 0.250000000000000000;

	m_vPoint[1][0] = 0.206829931610673204;
	m_vPoint[1][1] = 0.206829931610673204;
	m_vPoint[1][2] = 0.206829931610673204;

	m_vPoint[2][0] = 0.379510205167980388;
	m_vPoint[2][1] = 0.206829931610673204;
	m_vPoint[2][2] = 0.206829931610673204;

	m_vPoint[3][0] = 0.206829931610673204;
	m_vPoint[3][1] = 0.379510205167980388;
	m_vPoint[3][2] = 0.206829931610673204;

	m_vPoint[4][0] = 0.206829931610673204;
	m_vPoint[4][1] = 0.206829931610673204;
	m_vPoint[4][2] = 0.379510205167980388;

	m_vPoint[5][0] = 0.082103588310546723;
	m_vPoint[5][1] = 0.082103588310546723;
	m_vPoint[5][2] = 0.082103588310546723;

	m_vPoint[6][0] = 0.753689235068359831;
	m_vPoint[6][1] = 0.082103588310546723;
	m_vPoint[6][2] = 0.082103588310546723;

	m_vPoint[7][0] = 0.082103588310546723;
	m_vPoint[7][1] = 0.753689235068359831;
	m_vPoint[7][2] = 0.082103588310546723;

	m_vPoint[8][0] = 0.082103588310546723;
	m_vPoint[8][1] = 0.082103588310546723;
	m_vPoint[8][2] = 0.753689235068359831;

	m_vPoint[9][0] = 0.005781950505197997;
	m_vPoint[9][1] = 0.005781950505197997;
	m_vPoint[9][2] = 0.005781950505197997;

	m_vPoint[10][0] = 0.982654148484406008;
	m_vPoint[10][1] = 0.005781950505197997;
	m_vPoint[10][2] = 0.005781950505197997;

	m_vPoint[11][0] = 0.005781950505197997;
	m_vPoint[11][1] = 0.982654148484406008;
	m_vPoint[11][2] = 0.005781950505197997;

	m_vPoint[12][0] = 0.005781950505197997;
	m_vPoint[12][1] = 0.005781950505197997;
	m_vPoint[12][2] = 0.982654148484406008;

	m_vPoint[13][0] = 0.050532740018894224;
	m_vPoint[13][1] = 0.050532740018894224;
	m_vPoint[13][2] = 0.449467259981105776;

	m_vPoint[14][0] = 0.050532740018894224;
	m_vPoint[14][1] = 0.449467259981105776;
	m_vPoint[14][2] = 0.050532740018894224;

	m_vPoint[15][0] = 0.449467259981105776;
	m_vPoint[15][1] = 0.050532740018894224;
	m_vPoint[15][2] = 0.050532740018894224;

	m_vPoint[16][0] = 0.050532740018894224;
	m_vPoint[16][1] = 0.449467259981105776;
	m_vPoint[16][2] = 0.449467259981105776;

	m_vPoint[17][0] = 0.449467259981105776;
	m_vPoint[17][1] = 0.050532740018894224;
	m_vPoint[17][2] = 0.449467259981105776;

	m_vPoint[18][0] = 0.449467259981105776;
	m_vPoint[18][1] = 0.449467259981105776;
	m_vPoint[18][2] = 0.050532740018894224;

	m_vPoint[19][0] = 0.229066536116811140;
	m_vPoint[19][1] = 0.229066536116811140;
	m_vPoint[19][2] = 0.035639582788534044;

	m_vPoint[20][0] = 0.229066536116811140;
	m_vPoint[20][1] = 0.035639582788534044;
	m_vPoint[20][2] = 0.229066536116811140;

	m_vPoint[21][0] = 0.035639582788534044;
	m_vPoint[21][1] = 0.229066536116811140;
	m_vPoint[21][2] = 0.229066536116811140;

	m_vPoint[22][0] = 0.229066536116811140;
	m_vPoint[22][1] = 0.229066536116811140;
	m_vPoint[22][2] = 0.506227344977843677;

	m_vPoint[23][0] = 0.229066536116811140;
	m_vPoint[23][1] = 0.506227344977843677;
	m_vPoint[23][2] = 0.229066536116811140;

	m_vPoint[24][0] = 0.506227344977843677;
	m_vPoint[24][1] = 0.229066536116811140;
	m_vPoint[24][2] = 0.229066536116811140;

	m_vPoint[25][0] = 0.229066536116811140;
	m_vPoint[25][1] = 0.035639582788534044;
	m_vPoint[25][2] = 0.506227344977843677;

	m_vPoint[26][0] = 0.229066536116811140;
	m_vPoint[26][1] = 0.506227344977843677;
	m_vPoint[26][2] = 0.035639582788534044;

	m_vPoint[27][0] = 0.506227344977843677;
	m_vPoint[27][1] = 0.229066536116811140;
	m_vPoint[27][2] = 0.035639582788534044;

	m_vPoint[28][0] = 0.035639582788534044;
	m_vPoint[28][1] = 0.229066536116811140;
	m_vPoint[28][2] = 0.506227344977843677;

	m_vPoint[29][0] = 0.035639582788534044;
	m_vPoint[29][1] = 0.506227344977843677;
	m_vPoint[29][2] = 0.229066536116811140;

	m_vPoint[30][0] = 0.506227344977843677;
	m_vPoint[30][1] = 0.035639582788534044;
	m_vPoint[30][2] = 0.229066536116811140;

	m_vPoint[31][0] = 0.036607749553197424;
	m_vPoint[31][1] = 0.036607749553197424;
	m_vPoint[31][2] = 0.190486041934633456;

	m_vPoint[32][0] = 0.036607749553197424;
	m_vPoint[32][1] = 0.190486041934633456;
	m_vPoint[32][2] = 0.036607749553197424;

	m_vPoint[33][0] = 0.190486041934633456;
	m_vPoint[33][1] = 0.036607749553197424;
	m_vPoint[33][2] = 0.036607749553197424;

	m_vPoint[34][0] = 0.036607749553197424;
	m_vPoint[34][1] = 0.036607749553197424;
	m_vPoint[34][2] = 0.736298458958971697;

	m_vPoint[35][0] = 0.036607749553197424;
	m_vPoint[35][1] = 0.736298458958971697;
	m_vPoint[35][2] = 0.036607749553197424;

	m_vPoint[36][0] = 0.736298458958971697;
	m_vPoint[36][1] = 0.036607749553197424;
	m_vPoint[36][2] = 0.036607749553197424;

	m_vPoint[37][0] = 0.036607749553197424;
	m_vPoint[37][1] = 0.190486041934633456;
	m_vPoint[37][2] = 0.736298458958971697;

	m_vPoint[38][0] = 0.036607749553197424;
	m_vPoint[38][1] = 0.736298458958971697;
	m_vPoint[38][2] = 0.190486041934633456;

	m_vPoint[39][0] = 0.736298458958971697;
	m_vPoint[39][1] = 0.036607749553197424;
	m_vPoint[39][2] = 0.190486041934633456;

	m_vPoint[40][0] = 0.190486041934633456;
	m_vPoint[40][1] = 0.036607749553197424;
	m_vPoint[40][2] = 0.736298458958971697;

	m_vPoint[41][0] = 0.190486041934633456;
	m_vPoint[41][1] = 0.736298458958971697;
	m_vPoint[41][2] = 0.036607749553197424;

	m_vPoint[42][0] = 0.736298458958971697;
	m_vPoint[42][1] = 0.190486041934633456;
	m_vPoint[42][2] = 0.036607749553197424;

	m_vWeight[0] = 1./6. * -0.123001131951839495;
	m_vWeight[1] = 1./6. * 0.085501834937201407;
	m_vWeight[2] = 1./6. * 0.085501834937201407;
	m_vWeight[3] = 1./6. * 0.085501834937201407;
	m_vWeight[4] = 1./6. * 0.085501834937201407;
	m_vWeight[5] = 1./6. * 0.011802199878803406;
	m_vWeight[6] = 1./6. * 0.011802199878803406;
	m_vWeight[7] = 1./6. * 0.011802199878803406;
	m_vWeight[8] = 1./6. * 0.011802199878803406;
	m_vWeight[9] = 1./6. * 0.001019004654557324;
	m_vWeight[10] = 1./6. * 0.001019004654557324;
	m_vWeight[11] = 1./6. * 0.001019004654557324;
	m_vWeight[12] = 1./6. * 0.001019004654557324;
	m_vWeight[13] = 1./6. * 0.027478102946803691;
	m_vWeight[14] = 1./6. * 0.027478102946803691;
	m_vWeight[15] = 1./6. * 0.027478102946803691;
	m_vWeight[16] = 1./6. * 0.027478102946803691;
	m_vWeight[17] = 1./6. * 0.027478102946803691;
	m_vWeight[18] = 1./6. * 0.027478102946803691;
	m_vWeight[19] = 1./6. * 0.034226914852091511;
	m_vWeight[20] = 1./6. * 0.034226914852091511;
	m_vWeight[21] = 1./6. * 0.034226914852091511;
	m_vWeight[22] = 1./6. * 0.034226914852091511;
	m_vWeight[23] = 1./6. * 0.034226914852091511;
	m_vWeight[24] = 1./6. * 0.034226914852091511;
	m_vWeight[25] = 1./6. * 0.034226914852091511;
	m_vWeight[26] = 1./6. * 0.034226914852091511;
	m_vWeight[27] = 1./6. * 0.034226914852091511;
	m_vWeight[28] = 1./6. * 0.034226914852091511;
	m_vWeight[29] = 1./6. * 0.034226914852091511;
	m_vWeight[30] = 1./6. * 0.034226914852091511;
	m_vWeight[31] = 1./6. * 0.012843114846972556;
	m_vWeight[32] = 1./6. * 0.012843114846972556;
	m_vWeight[33] = 1./6. * 0.012843114846972556;
	m_vWeight[34] = 1./6. * 0.012843114846972556;
	m_vWeight[35] = 1./6. * 0.012843114846972556;
	m_vWeight[36] = 1./6. * 0.012843114846972556;
	m_vWeight[37] = 1./6. * 0.012843114846972556;
	m_vWeight[38] = 1./6. * 0.012843114846972556;
	m_vWeight[39] = 1./6. * 0.012843114846972556;
	m_vWeight[40] = 1./6. * 0.012843114846972556;
	m_vWeight[41] = 1./6. * 0.012843114846972556;
	m_vWeight[42] = 1./6. * 0.012843114846972556;
}




template <>
FlexGaussQuadrature<ReferenceTetrahedron>::FlexGaussQuadrature(int order)
{
	switch(order)
	{
	case 0:{
		const static GaussQuadrature<ReferenceTetrahedron, 0>& q0 
			= Provider<GaussQuadrature<ReferenceTetrahedron, 0> >::get();

		m_order = q0.order();
		m_numPoints = q0.size();
		m_pvPoint = q0.points();
		m_pvWeight = q0.weights();
		}break;

	case 1:{
		const static GaussQuadrature<ReferenceTetrahedron, 1>& q1 
			= Provider<GaussQuadrature<ReferenceTetrahedron, 1> >::get();

		m_order = q1.order();
		m_numPoints = q1.size();
		m_pvPoint = q1.points();
		m_pvWeight = q1.weights();
		}break;

	case 2:{
		const static GaussQuadrature<ReferenceTetrahedron, 2>& q2 
			= Provider<GaussQuadrature<ReferenceTetrahedron, 2> >::get();

		m_order = q2.order();
		m_numPoints = q2.size();
		m_pvPoint = q2.points();
		m_pvWeight = q2.weights();
		}break;

	case 3:{
		const static GaussQuadrature<ReferenceTetrahedron, 3>& q3 
			= Provider<GaussQuadrature<ReferenceTetrahedron, 3> >::get();

		m_order = q3.order();
		m_numPoints = q3.size();
		m_pvPoint = q3.points();
		m_pvWeight = q3.weights();
		}break;

	case 5:{
		const static GaussQuadrature<ReferenceTetrahedron, 5>& q5 
			= Provider<GaussQuadrature<ReferenceTetrahedron, 5> >::get();

		m_order = q5.order();
		m_numPoints = q5.size();
		m_pvPoint = q5.points();
		m_pvWeight = q5.weights();
		}break;

	case 6:{
		const static GaussQuadrature<ReferenceTetrahedron, 6>& q6 
			= Provider<GaussQuadrature<ReferenceTetrahedron, 6> >::get();

		m_order = q6.order();
		m_numPoints = q6.size();
		m_pvPoint = q6.points();
		m_pvWeight = q6.weights();
		}break;

	case 7:{
		const static GaussQuadrature<ReferenceTetrahedron, 7>& q7 
			= Provider<GaussQuadrature<ReferenceTetrahedron, 7> >::get();

		m_order = q7.order();
		m_numPoints = q7.size();
		m_pvPoint = q7.points();
		m_pvWeight = q7.weights();
		}break;

	case 8:{
		const static GaussQuadrature<ReferenceTetrahedron, 8>& q8 
			= Provider<GaussQuadrature<ReferenceTetrahedron, 8> >::get();

		m_order = q8.order();
		m_numPoints = q8.size();
		m_pvPoint = q8.points();
		m_pvWeight = q8.weights();
		}break;

	default: UG_ASSERT(0, "Order not availabile. Can not construct GaussQuadrature.\n");
	}
}



// register rules
template <>
bool RegisterGaussQuadRule<ReferenceTetrahedron>(QuadratureRuleProvider<ReferenceTetrahedron::dim>& factory)
{
	static FlexGaussQuadrature<ReferenceTetrahedron> gaussQuadratureReferenceTetrahedron_0(0);
	static FlexGaussQuadrature<ReferenceTetrahedron> gaussQuadratureReferenceTetrahedron_1(1);
	static FlexGaussQuadrature<ReferenceTetrahedron> gaussQuadratureReferenceTetrahedron_2(2);
	static FlexGaussQuadrature<ReferenceTetrahedron> gaussQuadratureReferenceTetrahedron_3(3);
	static FlexGaussQuadrature<ReferenceTetrahedron> gaussQuadratureReferenceTetrahedron_5(5);
	static FlexGaussQuadrature<ReferenceTetrahedron> gaussQuadratureReferenceTetrahedron_6(6);
	static FlexGaussQuadrature<ReferenceTetrahedron> gaussQuadratureReferenceTetrahedron_7(7);
	static FlexGaussQuadrature<ReferenceTetrahedron> gaussQuadratureReferenceTetrahedron_8(8);

	bool success = true;
	factory.register_rule<ReferenceTetrahedron>(gaussQuadratureReferenceTetrahedron_0);
	factory.register_rule<ReferenceTetrahedron>(gaussQuadratureReferenceTetrahedron_1);
	factory.register_rule<ReferenceTetrahedron>(gaussQuadratureReferenceTetrahedron_2);
	factory.register_rule<ReferenceTetrahedron>(gaussQuadratureReferenceTetrahedron_3);
	factory.register_rule<ReferenceTetrahedron>(gaussQuadratureReferenceTetrahedron_5);
	factory.register_rule<ReferenceTetrahedron>(gaussQuadratureReferenceTetrahedron_6);
	factory.register_rule<ReferenceTetrahedron>(gaussQuadratureReferenceTetrahedron_7);
	factory.register_rule<ReferenceTetrahedron>(gaussQuadratureReferenceTetrahedron_8);

	return success;
};

}; // namespace ug

