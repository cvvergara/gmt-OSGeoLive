/*--------------------------------------------------------------------
 *
 *	grdmath.h [Generated by make_math.sh]
 *
 *	Copyright (c) 1991-2011 by P. Wessel and W. H. F. Smith
 *	See LICENSE.TXT file for copying and redistribution conditions.
 *
 *	This program is free software; you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation; version 2 or any later version.
 *
 *	This program is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 *
 *	Contact info: gmt.soest.hawaii.edu
 *--------------------------------------------------------------------*/
/*	grdmath.h is automatically generated by make_math.sh;
 *	Do NOT edit manually!
 */

void grdmath_init (PFV ops[], GMT_LONG n_args[], GMT_LONG n_out[])
{

	/* Operator function		# of operands  		# of outputs */

	ops[0] = grd_ABS;		n_args[0] = 1;		n_out[0] = 1;
	ops[1] = grd_ACOS;		n_args[1] = 1;		n_out[1] = 1;
	ops[2] = grd_ACOSH;		n_args[2] = 1;		n_out[2] = 1;
	ops[3] = grd_ACOT;		n_args[3] = 1;		n_out[3] = 1;
	ops[4] = grd_ACSC;		n_args[4] = 1;		n_out[4] = 1;
	ops[5] = grd_ADD;		n_args[5] = 2;		n_out[5] = 1;
	ops[6] = grd_AND;		n_args[6] = 2;		n_out[6] = 1;
	ops[7] = grd_ASEC;		n_args[7] = 1;		n_out[7] = 1;
	ops[8] = grd_ASIN;		n_args[8] = 1;		n_out[8] = 1;
	ops[9] = grd_ASINH;		n_args[9] = 1;		n_out[9] = 1;
	ops[10] = grd_ATAN;		n_args[10] = 1;		n_out[10] = 1;
	ops[11] = grd_ATAN2;		n_args[11] = 2;		n_out[11] = 1;
	ops[12] = grd_ATANH;		n_args[12] = 1;		n_out[12] = 1;
	ops[13] = grd_BEI;		n_args[13] = 1;		n_out[13] = 1;
	ops[14] = grd_BER;		n_args[14] = 1;		n_out[14] = 1;
	ops[15] = grd_CAZ;		n_args[15] = 2;		n_out[15] = 1;
	ops[16] = grd_CBAZ;		n_args[16] = 2;		n_out[16] = 1;
	ops[17] = grd_CDIST;		n_args[17] = 2;		n_out[17] = 1;
	ops[18] = grd_CEIL;		n_args[18] = 1;		n_out[18] = 1;
	ops[19] = grd_CHICRIT;		n_args[19] = 2;		n_out[19] = 1;
	ops[20] = grd_CHIDIST;		n_args[20] = 2;		n_out[20] = 1;
	ops[21] = grd_CORRCOEFF;		n_args[21] = 2;		n_out[21] = 1;
	ops[22] = grd_COS;		n_args[22] = 1;		n_out[22] = 1;
	ops[23] = grd_COSD;		n_args[23] = 1;		n_out[23] = 1;
	ops[24] = grd_COSH;		n_args[24] = 1;		n_out[24] = 1;
	ops[25] = grd_COT;		n_args[25] = 1;		n_out[25] = 1;
	ops[26] = grd_COTD;		n_args[26] = 1;		n_out[26] = 1;
	ops[27] = grd_CPOISS;		n_args[27] = 2;		n_out[27] = 1;
	ops[28] = grd_CSC;		n_args[28] = 1;		n_out[28] = 1;
	ops[29] = grd_CSCD;		n_args[29] = 1;		n_out[29] = 1;
	ops[30] = grd_CURV;		n_args[30] = 1;		n_out[30] = 1;
	ops[31] = grd_D2DX2;		n_args[31] = 1;		n_out[31] = 1;
	ops[32] = grd_D2DXY;		n_args[32] = 1;		n_out[32] = 1;
	ops[33] = grd_D2DY2;		n_args[33] = 1;		n_out[33] = 1;
	ops[34] = grd_D2R;		n_args[34] = 1;		n_out[34] = 1;
	ops[35] = grd_DDX;		n_args[35] = 1;		n_out[35] = 1;
	ops[36] = grd_DDY;		n_args[36] = 1;		n_out[36] = 1;
	ops[37] = grd_DILOG;		n_args[37] = 1;		n_out[37] = 1;
	ops[38] = grd_DIV;		n_args[38] = 2;		n_out[38] = 1;
	ops[39] = grd_DUP;		n_args[39] = 1;		n_out[39] = 2;
	ops[40] = grd_EQ;		n_args[40] = 2;		n_out[40] = 1;
	ops[41] = grd_ERF;		n_args[41] = 1;		n_out[41] = 1;
	ops[42] = grd_ERFC;		n_args[42] = 1;		n_out[42] = 1;
	ops[43] = grd_ERFINV;		n_args[43] = 1;		n_out[43] = 1;
	ops[44] = grd_EXCH;		n_args[44] = 2;		n_out[44] = 2;
	ops[45] = grd_EXP;		n_args[45] = 1;		n_out[45] = 1;
	ops[46] = grd_EXTREMA;		n_args[46] = 1;		n_out[46] = 1;
	ops[47] = grd_FACT;		n_args[47] = 1;		n_out[47] = 1;
	ops[48] = grd_FCRIT;		n_args[48] = 3;		n_out[48] = 1;
	ops[49] = grd_FDIST;		n_args[49] = 3;		n_out[49] = 1;
	ops[50] = grd_FLIPLR;		n_args[50] = 1;		n_out[50] = 1;
	ops[51] = grd_FLIPUD;		n_args[51] = 1;		n_out[51] = 1;
	ops[52] = grd_FLOOR;		n_args[52] = 1;		n_out[52] = 1;
	ops[53] = grd_FMOD;		n_args[53] = 2;		n_out[53] = 1;
	ops[54] = grd_GE;		n_args[54] = 2;		n_out[54] = 1;
	ops[55] = grd_GT;		n_args[55] = 2;		n_out[55] = 1;
	ops[56] = grd_HYPOT;		n_args[56] = 2;		n_out[56] = 1;
	ops[57] = grd_I0;		n_args[57] = 1;		n_out[57] = 1;
	ops[58] = grd_I1;		n_args[58] = 1;		n_out[58] = 1;
	ops[59] = grd_IN;		n_args[59] = 2;		n_out[59] = 1;
	ops[60] = grd_INRANGE;		n_args[60] = 3;		n_out[60] = 1;
	ops[61] = grd_INSIDE;		n_args[61] = 1;		n_out[61] = 1;
	ops[62] = grd_INV;		n_args[62] = 1;		n_out[62] = 1;
	ops[63] = grd_ISNAN;		n_args[63] = 1;		n_out[63] = 1;
	ops[64] = grd_J0;		n_args[64] = 1;		n_out[64] = 1;
	ops[65] = grd_J1;		n_args[65] = 1;		n_out[65] = 1;
	ops[66] = grd_JN;		n_args[66] = 2;		n_out[66] = 1;
	ops[67] = grd_K0;		n_args[67] = 1;		n_out[67] = 1;
	ops[68] = grd_K1;		n_args[68] = 1;		n_out[68] = 1;
	ops[69] = grd_KEI;		n_args[69] = 1;		n_out[69] = 1;
	ops[70] = grd_KER;		n_args[70] = 1;		n_out[70] = 1;
	ops[71] = grd_KN;		n_args[71] = 2;		n_out[71] = 1;
	ops[72] = grd_KURT;		n_args[72] = 1;		n_out[72] = 1;
	ops[73] = grd_LDIST;		n_args[73] = 1;		n_out[73] = 1;
	ops[74] = grd_LE;		n_args[74] = 2;		n_out[74] = 1;
	ops[75] = grd_LMSSCL;		n_args[75] = 1;		n_out[75] = 1;
	ops[76] = grd_LOG;		n_args[76] = 1;		n_out[76] = 1;
	ops[77] = grd_LOG10;		n_args[77] = 1;		n_out[77] = 1;
	ops[78] = grd_LOG1P;		n_args[78] = 1;		n_out[78] = 1;
	ops[79] = grd_LOG2;		n_args[79] = 1;		n_out[79] = 1;
	ops[80] = grd_LOWER;		n_args[80] = 1;		n_out[80] = 1;
	ops[81] = grd_LRAND;		n_args[81] = 2;		n_out[81] = 1;
	ops[82] = grd_LT;		n_args[82] = 2;		n_out[82] = 1;
	ops[83] = grd_MAD;		n_args[83] = 1;		n_out[83] = 1;
	ops[84] = grd_MAX;		n_args[84] = 2;		n_out[84] = 1;
	ops[85] = grd_MEAN;		n_args[85] = 1;		n_out[85] = 1;
	ops[86] = grd_MED;		n_args[86] = 1;		n_out[86] = 1;
	ops[87] = grd_MIN;		n_args[87] = 2;		n_out[87] = 1;
	ops[88] = grd_MOD;		n_args[88] = 2;		n_out[88] = 1;
	ops[89] = grd_MODE;		n_args[89] = 1;		n_out[89] = 1;
	ops[90] = grd_MUL;		n_args[90] = 2;		n_out[90] = 1;
	ops[91] = grd_NAN;		n_args[91] = 2;		n_out[91] = 1;
	ops[92] = grd_NEG;		n_args[92] = 1;		n_out[92] = 1;
	ops[93] = grd_NEQ;		n_args[93] = 2;		n_out[93] = 1;
	ops[94] = grd_NOT;		n_args[94] = 1;		n_out[94] = 1;
	ops[95] = grd_NRAND;		n_args[95] = 2;		n_out[95] = 1;
	ops[96] = grd_OR;		n_args[96] = 2;		n_out[96] = 1;
	ops[97] = grd_PDIST;		n_args[97] = 1;		n_out[97] = 1;
	ops[98] = grd_PLM;		n_args[98] = 3;		n_out[98] = 1;
	ops[99] = grd_PLMg;		n_args[99] = 3;		n_out[99] = 1;
	ops[100] = grd_POP;		n_args[100] = 1;		n_out[100] = 0;
	ops[101] = grd_POW;		n_args[101] = 2;		n_out[101] = 1;
	ops[102] = grd_PQUANT;		n_args[102] = 2;		n_out[102] = 1;
	ops[103] = grd_PSI;		n_args[103] = 1;		n_out[103] = 1;
	ops[104] = grd_PV;		n_args[104] = 3;		n_out[104] = 1;
	ops[105] = grd_QV;		n_args[105] = 3;		n_out[105] = 1;
	ops[106] = grd_R2;		n_args[106] = 2;		n_out[106] = 1;
	ops[107] = grd_R2D;		n_args[107] = 1;		n_out[107] = 1;
	ops[108] = grd_RAND;		n_args[108] = 2;		n_out[108] = 1;
	ops[109] = grd_RINT;		n_args[109] = 1;		n_out[109] = 1;
	ops[110] = grd_ROTX;		n_args[110] = 2;		n_out[110] = 1;
	ops[111] = grd_ROTY;		n_args[111] = 2;		n_out[111] = 1;
	ops[112] = grd_SAZ;		n_args[112] = 2;		n_out[112] = 1;
	ops[113] = grd_SBAZ;		n_args[113] = 2;		n_out[113] = 1;
	ops[114] = grd_SDIST;		n_args[114] = 2;		n_out[114] = 1;
	ops[115] = grd_SEC;		n_args[115] = 1;		n_out[115] = 1;
	ops[116] = grd_SECD;		n_args[116] = 1;		n_out[116] = 1;
	ops[117] = grd_SIGN;		n_args[117] = 1;		n_out[117] = 1;
	ops[118] = grd_SIN;		n_args[118] = 1;		n_out[118] = 1;
	ops[119] = grd_SINC;		n_args[119] = 1;		n_out[119] = 1;
	ops[120] = grd_SIND;		n_args[120] = 1;		n_out[120] = 1;
	ops[121] = grd_SINH;		n_args[121] = 1;		n_out[121] = 1;
	ops[122] = grd_SKEW;		n_args[122] = 1;		n_out[122] = 1;
	ops[123] = grd_SQR;		n_args[123] = 1;		n_out[123] = 1;
	ops[124] = grd_SQRT;		n_args[124] = 1;		n_out[124] = 1;
	ops[125] = grd_STD;		n_args[125] = 1;		n_out[125] = 1;
	ops[126] = grd_STEP;		n_args[126] = 1;		n_out[126] = 1;
	ops[127] = grd_STEPX;		n_args[127] = 1;		n_out[127] = 1;
	ops[128] = grd_STEPY;		n_args[128] = 1;		n_out[128] = 1;
	ops[129] = grd_SUB;		n_args[129] = 2;		n_out[129] = 1;
	ops[130] = grd_TAN;		n_args[130] = 1;		n_out[130] = 1;
	ops[131] = grd_TAND;		n_args[131] = 1;		n_out[131] = 1;
	ops[132] = grd_TANH;		n_args[132] = 1;		n_out[132] = 1;
	ops[133] = grd_TCRIT;		n_args[133] = 2;		n_out[133] = 1;
	ops[134] = grd_TDIST;		n_args[134] = 2;		n_out[134] = 1;
	ops[135] = grd_TN;		n_args[135] = 2;		n_out[135] = 1;
	ops[136] = grd_UPPER;		n_args[136] = 1;		n_out[136] = 1;
	ops[137] = grd_XOR;		n_args[137] = 2;		n_out[137] = 1;
	ops[138] = grd_Y0;		n_args[138] = 1;		n_out[138] = 1;
	ops[139] = grd_Y1;		n_args[139] = 1;		n_out[139] = 1;
	ops[140] = grd_YLM;		n_args[140] = 2;		n_out[140] = 2;
	ops[141] = grd_YLMg;		n_args[141] = 2;		n_out[141] = 2;
	ops[142] = grd_YN;		n_args[142] = 2;		n_out[142] = 1;
	ops[143] = grd_ZCRIT;		n_args[143] = 1;		n_out[143] = 1;
	ops[144] = grd_ZDIST;		n_args[144] = 1;		n_out[144] = 1;
}