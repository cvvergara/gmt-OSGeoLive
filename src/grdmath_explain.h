/*--------------------------------------------------------------------
 *
 *	grdmath_explain.h [Generated by make_math.sh]
 *
 *	Copyright (c) 1991-2013 by P. Wessel and W. H. F. Smith
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
/*	grdmath_explain.h is automatically generated by make_math.sh;
 *	Do NOT edit manually!
 *
 */
		fprintf (stderr, "	ABS	1 1	abs (A).\n");
		fprintf (stderr, "	ACOS	1 1	acos (A).\n");
		fprintf (stderr, "	ACOSH	1 1	acosh (A).\n");
		fprintf (stderr, "	ACOT	1 1	acot (A).\n");
		fprintf (stderr, "	ACSC	1 1	acsc (A).\n");
		fprintf (stderr, "	ADD	2 1	A + B.\n");
		fprintf (stderr, "	AND	2 1	NaN if A and B == NaN, B if A == NaN, else A.\n");
		fprintf (stderr, "	ASEC	1 1	asec (A).\n");
		fprintf (stderr, "	ASIN	1 1	asin (A).\n");
		fprintf (stderr, "	ASINH	1 1	asinh (A).\n");
		fprintf (stderr, "	ATAN	1 1	atan (A).\n");
		fprintf (stderr, "	ATAN2	2 1	atan2 (A, B).\n");
		fprintf (stderr, "	ATANH	1 1	atanh (A).\n");
		fprintf (stderr, "	BEI	1 1	bei (A).\n");
		fprintf (stderr, "	BER	1 1	ber (A).\n");
		fprintf (stderr, "	CAZ	2 1	Cartesian azimuth from grid nodes to stack x,y.\n");
		fprintf (stderr, "	CBAZ	2 1	Cartesian backazimuth from grid nodes to stack x,y.\n");
		fprintf (stderr, "	CDIST	2 1	Cartesian distance between grid nodes and stack x,y.\n");
		fprintf (stderr, "	CEIL	1 1	ceil (A) (smallest integer >= A).\n");
		fprintf (stderr, "	CHICRIT	2 1	Critical value for chi-squared-distribution, with alpha = A and n = B.\n");
		fprintf (stderr, "	CHIDIST	2 1	chi-squared-distribution P(chi2,n), with chi2 = A and n = B.\n");
		fprintf (stderr, "	CORRCOEFF	2 1	Correlation coefficient r(A, B).\n");
		fprintf (stderr, "	COS	1 1	cos (A) (A in radians).\n");
		fprintf (stderr, "	COSD	1 1	cos (A) (A in degrees).\n");
		fprintf (stderr, "	COSH	1 1	cosh (A).\n");
		fprintf (stderr, "	COT	1 1	cot (A) (A in radians).\n");
		fprintf (stderr, "	COTD	1 1	cot (A) (A in degrees).\n");
		fprintf (stderr, "	CPOISS	2 1	Cumulative Poisson distribution F(x,lambda), with x = A and lambda = B.\n");
		fprintf (stderr, "	CSC	1 1	csc (A) (A in radians).\n");
		fprintf (stderr, "	CSCD	1 1	csc (A) (A in degrees).\n");
		fprintf (stderr, "	CURV	1 1	Curvature of A (Laplacian).\n");
		fprintf (stderr, "	D2DX2	1 1	d^2(A)/dx^2 2nd derivative.\n");
		fprintf (stderr, "	D2DXY	1 1	d^2(A)/dxdy 2nd derivative.\n");
		fprintf (stderr, "	D2DY2	1 1	d^2(A)/dy^2 2nd derivative.\n");
		fprintf (stderr, "	D2R	1 1	Converts Degrees to Radians.\n");
		fprintf (stderr, "	DDX	1 1	d(A)/dx Central 1st derivative.\n");
		fprintf (stderr, "	DDY	1 1	d(A)/dy Central 1st derivative.\n");
		fprintf (stderr, "	DILOG	1 1	dilog (A).\n");
		fprintf (stderr, "	DIV	2 1	A / B.\n");
		fprintf (stderr, "	DUP	1 2	Places duplicate of A on the stack.\n");
		fprintf (stderr, "	EQ	2 1	1 if A == B, else 0.\n");
		fprintf (stderr, "	ERF	1 1	Error function erf (A).\n");
		fprintf (stderr, "	ERFC	1 1	Complementary Error function erfc (A).\n");
		fprintf (stderr, "	ERFINV	1 1	Inverse error function of A.\n");
		fprintf (stderr, "	EXCH	2 2	Exchanges A and B on the stack.\n");
		fprintf (stderr, "	EXP	1 1	exp (A).\n");
		fprintf (stderr, "	EXTREMA	1 1	Local Extrema: +2/-2 is max/min, +1/-1 is saddle with max/min in x, 0 elsewhere.\n");
		fprintf (stderr, "	FACT	1 1	A! (A factorial).\n");
		fprintf (stderr, "	FCRIT	3 1	Critical value for F-distribution, with alpha = A, n1 = B, and n2 = C.\n");
		fprintf (stderr, "	FDIST	3 1	F-distribution Q(F,n1,n2), with F = A, n1 = B, and n2 = C.\n");
		fprintf (stderr, "	FLIPLR	1 1	Reverse order of values in each row.\n");
		fprintf (stderr, "	FLIPUD	1 1	Reverse order of values in each column.\n");
		fprintf (stderr, "	FLOOR	1 1	floor (A) (greatest integer <= A).\n");
		fprintf (stderr, "	FMOD	2 1	A %% B (remainder after truncated division).\n");
		fprintf (stderr, "	GE	2 1	1 if A >= B, else 0.\n");
		fprintf (stderr, "	GT	2 1	1 if A > B, else 0.\n");
		fprintf (stderr, "	HYPOT	2 1	hypot (A, B) = sqrt (A*A + B*B).\n");
		fprintf (stderr, "	I0	1 1	Modified Bessel function of A (1st kind, order 0).\n");
		fprintf (stderr, "	I1	1 1	Modified Bessel function of A (1st kind, order 1).\n");
		fprintf (stderr, "	IN	2 1	Modified Bessel function of A (1st kind, order B).\n");
		fprintf (stderr, "	INRANGE	3 1	1 if B <= A <= C, else 0.\n");
		fprintf (stderr, "	INSIDE	1 1	1 when inside or on polygon(s) in A, else 0.\n");
		fprintf (stderr, "	INV	1 1	1 / A.\n");
		fprintf (stderr, "	ISNAN	1 1	1 if A == NaN, else 0.\n");
		fprintf (stderr, "	J0	1 1	Bessel function of A (1st kind, order 0).\n");
		fprintf (stderr, "	J1	1 1	Bessel function of A (1st kind, order 1).\n");
		fprintf (stderr, "	JN	2 1	Bessel function of A (1st kind, order B).\n");
		fprintf (stderr, "	K0	1 1	Modified Kelvin function of A (2nd kind, order 0).\n");
		fprintf (stderr, "	K1	1 1	Modified Bessel function of A (2nd kind, order 1).\n");
		fprintf (stderr, "	KEI	1 1	kei (A).\n");
		fprintf (stderr, "	KER	1 1	ker (A).\n");
		fprintf (stderr, "	KN	2 1	Modified Bessel function of A (2nd kind, order B).\n");
		fprintf (stderr, "	KURT	1 1	Kurtosis of A.\n");
		fprintf (stderr, "	LDIST	1 1	Compute distance from lines in multi-segment ASCII file A.\n");
		fprintf (stderr, "	LE	2 1	1 if A <= B, else 0.\n");
		fprintf (stderr, "	LMSSCL	1 1	LMS scale estimate (LMS STD) of A.\n");
		fprintf (stderr, "	LOG	1 1	log (A) (natural log).\n");
		fprintf (stderr, "	LOG10	1 1	log10 (A) (base 10).\n");
		fprintf (stderr, "	LOG1P	1 1	log (1+A) (accurate for small A).\n");
		fprintf (stderr, "	LOG2	1 1	log2 (A) (base 2).\n");
		fprintf (stderr, "	LOWER	1 1	The lowest (minimum) value of A.\n");
		fprintf (stderr, "	LRAND	2 1	Laplace random noise with mean A and std. deviation B.\n");
		fprintf (stderr, "	LT	2 1	1 if A < B, else 0.\n");
		fprintf (stderr, "	MAD	1 1	Median Absolute Deviation (L1 STD) of A.\n");
		fprintf (stderr, "	MAX	2 1	Maximum of A and B.\n");
		fprintf (stderr, "	MEAN	1 1	Mean value of A.\n");
		fprintf (stderr, "	MED	1 1	Median value of A.\n");
		fprintf (stderr, "	MIN	2 1	Minimum of A and B.\n");
		fprintf (stderr, "	MOD	2 1	A mod B (remainder after floored division).\n");
		fprintf (stderr, "	MODE	1 1	Mode value (Least Median of Squares) of A.\n");
		fprintf (stderr, "	MUL	2 1	A * B.\n");
		fprintf (stderr, "	NAN	2 1	NaN if A == B, else A.\n");
		fprintf (stderr, "	NEG	1 1	-A.\n");
		fprintf (stderr, "	NEQ	2 1	1 if A != B, else 0.\n");
		fprintf (stderr, "	NOT	1 1	NaN if A == NaN, 1 if A == 0, else 0.\n");
		fprintf (stderr, "	NRAND	2 1	Normal, random values with mean A and std. deviation B.\n");
		fprintf (stderr, "	OR	2 1	NaN if A or B == NaN, else A.\n");
		fprintf (stderr, "	PDIST	1 1	Compute distance from points in ASCII file A.\n");
		fprintf (stderr, "	PLM	3 1	Associated Legendre polynomial P(A) degree B order C.\n");
		fprintf (stderr, "	PLMg	3 1	Normalized associated Legendre polynomial P(A) degree B order C (geophysical convention).\n");
		fprintf (stderr, "	POP	1 0	Delete top element from the stack.\n");
		fprintf (stderr, "	POW	2 1	A ^ B.\n");
		fprintf (stderr, "	PQUANT	2 1	The B'th Quantile (0-100%%) of A.\n");
		fprintf (stderr, "	PSI	1 1	Psi (or Digamma) of A.\n");
		fprintf (stderr, "	PV	3 1	Legendre function Pv(A) of degree v = real(B) + imag(C).\n");
		fprintf (stderr, "	QV	3 1	Legendre function Qv(A) of degree v = real(B) + imag(C).\n");
		fprintf (stderr, "	R2	2 1	R2 = A^2 + B^2.\n");
		fprintf (stderr, "	R2D	1 1	Convert Radians to Degrees.\n");
		fprintf (stderr, "	RAND	2 1	Uniform random values between A and B.\n");
		fprintf (stderr, "	RINT	1 1	rint (A) (nearest integer).\n");
		fprintf (stderr, "	ROTX	2 1	Rotate A by the (constant) shift B in x-direction.\n");
		fprintf (stderr, "	ROTY	2 1	Rotate A by the (constant) shift B in y-direction.\n");
		fprintf (stderr, "	SAZ	2 1	Spherical azimuth from grid nodes to stack x,y.\n");
		fprintf (stderr, "	SBAZ	2 1	Spherical backazimuth from grid nodes to stack x,y.\n");
		fprintf (stderr, "	SDIST	2 1	Spherical (Great circle) distance (in degrees) between grid nodes and stack lon,lat (A, B).\n");
		fprintf (stderr, "	SEC	1 1	sec (A) (A in radians).\n");
		fprintf (stderr, "	SECD	1 1	sec (A) (A in degrees).\n");
		fprintf (stderr, "	SIGN	1 1	sign (+1 or -1) of A.\n");
		fprintf (stderr, "	SIN	1 1	sin (A) (A in radians).\n");
		fprintf (stderr, "	SINC	1 1	sinc (A) (sin (pi*A)/(pi*A)).\n");
		fprintf (stderr, "	SIND	1 1	sin (A) (A in degrees).\n");
		fprintf (stderr, "	SINH	1 1	sinh (A).\n");
		fprintf (stderr, "	SKEW	1 1	Skewness of A.\n");
		fprintf (stderr, "	SQR	1 1	A^2.\n");
		fprintf (stderr, "	SQRT	1 1	sqrt (A).\n");
		fprintf (stderr, "	STD	1 1	Standard deviation of A.\n");
		fprintf (stderr, "	STEP	1 1	Heaviside step function: H(A).\n");
		fprintf (stderr, "	STEPX	1 1	Heaviside step function in x: H(x-A).\n");
		fprintf (stderr, "	STEPY	1 1	Heaviside step function in y: H(y-A).\n");
		fprintf (stderr, "	SUB	2 1	A - B.\n");
		fprintf (stderr, "	TAN	1 1	tan (A) (A in radians).\n");
		fprintf (stderr, "	TAND	1 1	tan (A) (A in degrees).\n");
		fprintf (stderr, "	TANH	1 1	tanh (A).\n");
		fprintf (stderr, "	TCRIT	2 1	Critical value for Student's t-distribution, with alpha = A and n = B.\n");
		fprintf (stderr, "	TDIST	2 1	Student's t-distribution A(t,n), with t = A, and n = B.\n");
		fprintf (stderr, "	TN	2 1	Chebyshev polynomial Tn(-1<t<+1,n), with t = A, and n = B.\n");
		fprintf (stderr, "	UPPER	1 1	The highest (maximum) value of A.\n");
		fprintf (stderr, "	XOR	2 1	B if A == NaN, else A.\n");
		fprintf (stderr, "	Y0	1 1	Bessel function of A (2nd kind, order 0).\n");
		fprintf (stderr, "	Y1	1 1	Bessel function of A (2nd kind, order 1).\n");
		fprintf (stderr, "	YLM	2 2	Re and Im orthonormalized spherical harmonics degree A order B.\n");
		fprintf (stderr, "	YLMg	2 2	Cos and Sin normalized spherical harmonics degree A order B (geophysical convention).\n");
		fprintf (stderr, "	YN	2 1	Bessel function of A (2nd kind, order B).\n");
		fprintf (stderr, "	ZCRIT	1 1	Critical value for the normal-distribution, with alpha = A.\n");
		fprintf (stderr, "	ZDIST	1 1	Cumulative normal-distribution C(x), with x = A.\n");
