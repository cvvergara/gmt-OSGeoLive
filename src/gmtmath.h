/*--------------------------------------------------------------------
 *
 *	gmtmath.h [Generated by make_math.sh]
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
/*	gmtmath.h is automatically generated by make_math.sh;
 *	Do NOT edit manually!
 */

void gmtmath_init (PFV ops[], GMT_LONG n_args[], GMT_LONG n_out[])
{

	/* Operator function		# of operands  		# of outputs */

	ops[0] = table_ABS;		n_args[0] = 1;		n_out[0] = 1;
	ops[1] = table_ACOS;		n_args[1] = 1;		n_out[1] = 1;
	ops[2] = table_ACOSH;		n_args[2] = 1;		n_out[2] = 1;
	ops[3] = table_ACOT;		n_args[3] = 1;		n_out[3] = 1;
	ops[4] = table_ACSC;		n_args[4] = 1;		n_out[4] = 1;
	ops[5] = table_ADD;		n_args[5] = 2;		n_out[5] = 1;
	ops[6] = table_AND;		n_args[6] = 2;		n_out[6] = 1;
	ops[7] = table_ASEC;		n_args[7] = 1;		n_out[7] = 1;
	ops[8] = table_ASIN;		n_args[8] = 1;		n_out[8] = 1;
	ops[9] = table_ASINH;		n_args[9] = 1;		n_out[9] = 1;
	ops[10] = table_ATAN;		n_args[10] = 1;		n_out[10] = 1;
	ops[11] = table_ATAN2;		n_args[11] = 2;		n_out[11] = 1;
	ops[12] = table_ATANH;		n_args[12] = 1;		n_out[12] = 1;
	ops[13] = table_BEI;		n_args[13] = 1;		n_out[13] = 1;
	ops[14] = table_BER;		n_args[14] = 1;		n_out[14] = 1;
	ops[15] = table_CEIL;		n_args[15] = 1;		n_out[15] = 1;
	ops[16] = table_CHICRIT;		n_args[16] = 2;		n_out[16] = 1;
	ops[17] = table_CHIDIST;		n_args[17] = 2;		n_out[17] = 1;
	ops[18] = table_COL;		n_args[18] = 1;		n_out[18] = 1;
	ops[19] = table_CORRCOEFF;		n_args[19] = 2;		n_out[19] = 1;
	ops[20] = table_COS;		n_args[20] = 1;		n_out[20] = 1;
	ops[21] = table_COSD;		n_args[21] = 1;		n_out[21] = 1;
	ops[22] = table_COSH;		n_args[22] = 1;		n_out[22] = 1;
	ops[23] = table_COT;		n_args[23] = 1;		n_out[23] = 1;
	ops[24] = table_COTD;		n_args[24] = 1;		n_out[24] = 1;
	ops[25] = table_CPOISS;		n_args[25] = 2;		n_out[25] = 1;
	ops[26] = table_CSC;		n_args[26] = 1;		n_out[26] = 1;
	ops[27] = table_CSCD;		n_args[27] = 1;		n_out[27] = 1;
	ops[28] = table_D2DT2;		n_args[28] = 1;		n_out[28] = 1;
	ops[29] = table_D2R;		n_args[29] = 1;		n_out[29] = 1;
	ops[30] = table_DDT;		n_args[30] = 1;		n_out[30] = 1;
	ops[31] = table_DILOG;		n_args[31] = 1;		n_out[31] = 1;
	ops[32] = table_DIV;		n_args[32] = 2;		n_out[32] = 1;
	ops[33] = table_DUP;		n_args[33] = 1;		n_out[33] = 2;
	ops[34] = table_EQ;		n_args[34] = 2;		n_out[34] = 1;
	ops[35] = table_ERF;		n_args[35] = 1;		n_out[35] = 1;
	ops[36] = table_ERFC;		n_args[36] = 1;		n_out[36] = 1;
	ops[37] = table_ERFINV;		n_args[37] = 1;		n_out[37] = 1;
	ops[38] = table_EXCH;		n_args[38] = 2;		n_out[38] = 2;
	ops[39] = table_EXP;		n_args[39] = 1;		n_out[39] = 1;
	ops[40] = table_FACT;		n_args[40] = 1;		n_out[40] = 1;
	ops[41] = table_FCRIT;		n_args[41] = 3;		n_out[41] = 1;
	ops[42] = table_FDIST;		n_args[42] = 3;		n_out[42] = 1;
	ops[43] = table_FLIPUD;		n_args[43] = 1;		n_out[43] = 1;
	ops[44] = table_FLOOR;		n_args[44] = 1;		n_out[44] = 1;
	ops[45] = table_FMOD;		n_args[45] = 2;		n_out[45] = 1;
	ops[46] = table_GE;		n_args[46] = 2;		n_out[46] = 1;
	ops[47] = table_GT;		n_args[47] = 2;		n_out[47] = 1;
	ops[48] = table_HYPOT;		n_args[48] = 2;		n_out[48] = 1;
	ops[49] = table_I0;		n_args[49] = 1;		n_out[49] = 1;
	ops[50] = table_I1;		n_args[50] = 1;		n_out[50] = 1;
	ops[51] = table_IN;		n_args[51] = 2;		n_out[51] = 1;
	ops[52] = table_INRANGE;		n_args[52] = 3;		n_out[52] = 1;
	ops[53] = table_INT;		n_args[53] = 1;		n_out[53] = 1;
	ops[54] = table_INV;		n_args[54] = 1;		n_out[54] = 1;
	ops[55] = table_ISNAN;		n_args[55] = 1;		n_out[55] = 1;
	ops[56] = table_J0;		n_args[56] = 1;		n_out[56] = 1;
	ops[57] = table_J1;		n_args[57] = 1;		n_out[57] = 1;
	ops[58] = table_JN;		n_args[58] = 2;		n_out[58] = 1;
	ops[59] = table_K0;		n_args[59] = 1;		n_out[59] = 1;
	ops[60] = table_K1;		n_args[60] = 1;		n_out[60] = 1;
	ops[61] = table_KEI;		n_args[61] = 1;		n_out[61] = 1;
	ops[62] = table_KER;		n_args[62] = 1;		n_out[62] = 1;
	ops[63] = table_KN;		n_args[63] = 2;		n_out[63] = 1;
	ops[64] = table_KURT;		n_args[64] = 1;		n_out[64] = 1;
	ops[65] = table_LE;		n_args[65] = 2;		n_out[65] = 1;
	ops[66] = table_LMSSCL;		n_args[66] = 1;		n_out[66] = 1;
	ops[67] = table_LOG;		n_args[67] = 1;		n_out[67] = 1;
	ops[68] = table_LOG10;		n_args[68] = 1;		n_out[68] = 1;
	ops[69] = table_LOG1P;		n_args[69] = 1;		n_out[69] = 1;
	ops[70] = table_LOG2;		n_args[70] = 1;		n_out[70] = 1;
	ops[71] = table_LOWER;		n_args[71] = 1;		n_out[71] = 1;
	ops[72] = table_LRAND;		n_args[72] = 2;		n_out[72] = 1;
	ops[73] = table_LSQFIT;		n_args[73] = 1;		n_out[73] = 0;
	ops[74] = table_LT;		n_args[74] = 2;		n_out[74] = 1;
	ops[75] = table_MAD;		n_args[75] = 1;		n_out[75] = 1;
	ops[76] = table_MAX;		n_args[76] = 2;		n_out[76] = 1;
	ops[77] = table_MEAN;		n_args[77] = 1;		n_out[77] = 1;
	ops[78] = table_MED;		n_args[78] = 1;		n_out[78] = 1;
	ops[79] = table_MIN;		n_args[79] = 2;		n_out[79] = 1;
	ops[80] = table_MOD;		n_args[80] = 2;		n_out[80] = 1;
	ops[81] = table_MODE;		n_args[81] = 1;		n_out[81] = 1;
	ops[82] = table_MUL;		n_args[82] = 2;		n_out[82] = 1;
	ops[83] = table_NAN;		n_args[83] = 2;		n_out[83] = 1;
	ops[84] = table_NEG;		n_args[84] = 1;		n_out[84] = 1;
	ops[85] = table_NEQ;		n_args[85] = 2;		n_out[85] = 1;
	ops[86] = table_NOT;		n_args[86] = 1;		n_out[86] = 1;
	ops[87] = table_NRAND;		n_args[87] = 2;		n_out[87] = 1;
	ops[88] = table_OR;		n_args[88] = 2;		n_out[88] = 1;
	ops[89] = table_PLM;		n_args[89] = 3;		n_out[89] = 1;
	ops[90] = table_PLMg;		n_args[90] = 3;		n_out[90] = 1;
	ops[91] = table_POP;		n_args[91] = 1;		n_out[91] = 0;
	ops[92] = table_POW;		n_args[92] = 2;		n_out[92] = 1;
	ops[93] = table_PQUANT;		n_args[93] = 2;		n_out[93] = 1;
	ops[94] = table_PSI;		n_args[94] = 1;		n_out[94] = 1;
	ops[95] = table_PV;		n_args[95] = 3;		n_out[95] = 1;
	ops[96] = table_QV;		n_args[96] = 3;		n_out[96] = 1;
	ops[97] = table_R2;		n_args[97] = 2;		n_out[97] = 1;
	ops[98] = table_R2D;		n_args[98] = 1;		n_out[98] = 1;
	ops[99] = table_RAND;		n_args[99] = 2;		n_out[99] = 1;
	ops[100] = table_RINT;		n_args[100] = 1;		n_out[100] = 1;
	ops[101] = table_ROOTS;		n_args[101] = 2;		n_out[101] = 1;
	ops[102] = table_ROTT;		n_args[102] = 2;		n_out[102] = 1;
	ops[103] = table_SEC;		n_args[103] = 1;		n_out[103] = 1;
	ops[104] = table_SECD;		n_args[104] = 1;		n_out[104] = 1;
	ops[105] = table_SIGN;		n_args[105] = 1;		n_out[105] = 1;
	ops[106] = table_SIN;		n_args[106] = 1;		n_out[106] = 1;
	ops[107] = table_SINC;		n_args[107] = 1;		n_out[107] = 1;
	ops[108] = table_SIND;		n_args[108] = 1;		n_out[108] = 1;
	ops[109] = table_SINH;		n_args[109] = 1;		n_out[109] = 1;
	ops[110] = table_SKEW;		n_args[110] = 1;		n_out[110] = 1;
	ops[111] = table_SQR;		n_args[111] = 1;		n_out[111] = 1;
	ops[112] = table_SQRT;		n_args[112] = 1;		n_out[112] = 1;
	ops[113] = table_STD;		n_args[113] = 1;		n_out[113] = 1;
	ops[114] = table_STEP;		n_args[114] = 1;		n_out[114] = 1;
	ops[115] = table_STEPT;		n_args[115] = 1;		n_out[115] = 1;
	ops[116] = table_SUB;		n_args[116] = 2;		n_out[116] = 1;
	ops[117] = table_SUM;		n_args[117] = 1;		n_out[117] = 1;
	ops[118] = table_TAN;		n_args[118] = 1;		n_out[118] = 1;
	ops[119] = table_TAND;		n_args[119] = 1;		n_out[119] = 1;
	ops[120] = table_TANH;		n_args[120] = 1;		n_out[120] = 1;
	ops[121] = table_TCRIT;		n_args[121] = 2;		n_out[121] = 1;
	ops[122] = table_TDIST;		n_args[122] = 2;		n_out[122] = 1;
	ops[123] = table_TN;		n_args[123] = 2;		n_out[123] = 1;
	ops[124] = table_UPPER;		n_args[124] = 1;		n_out[124] = 1;
	ops[125] = table_XOR;		n_args[125] = 2;		n_out[125] = 1;
	ops[126] = table_Y0;		n_args[126] = 1;		n_out[126] = 1;
	ops[127] = table_Y1;		n_args[127] = 1;		n_out[127] = 1;
	ops[128] = table_YN;		n_args[128] = 2;		n_out[128] = 1;
	ops[129] = table_ZCRIT;		n_args[129] = 1;		n_out[129] = 1;
	ops[130] = table_ZDIST;		n_args[130] = 1;		n_out[130] = 1;
}