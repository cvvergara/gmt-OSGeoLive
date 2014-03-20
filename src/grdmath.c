/*--------------------------------------------------------------------
 *	$Id: grdmath.c 10217 2014-02-27 19:47:18Z pwessel $
 *
 *	Copyright (c) 1991-2014 by P. Wessel and W. H. F. Smith
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
/*
 * grdmath.c is a reverse polish calculator that operates on grid files
 * (and constants) and perform basic mathematical operations
 * on them like add, multiply, etc.
 * Some operators only work on one operand (e.g., log, exp)
 *
 * Author:	Paul Wessel
 * Date:	22-AUG-1988
 * Revised:	27-JUL-2000 PW: Added EXTREMA operator
 *		18-AUG-2000 PW: Added LT, LE, EQ, GE, GT, NAN, ISNAN, XOR, MODE, MAD, LMSSCL
 *		20-AUG-2000 PW: Added FDIST, TDIST, CHI2DIST
 *		23-AUG-2000 PW: Added LOWER and UPPER
 *		28-NOV-2001 PW: Added Critical values for Chi2, F, T, and Z distributions
 *		12-NOV-2002 PW: Added Cartesian and Spherical Azimuth operators
 *		27-JAN-2004 PW: Added SINC and NEQ
 *		28-MAR-2004 PW: Added FLIPUD, FLIPLR, ROTX, ROTY, INRANGE, LDIST, PDIST, INSIDE
 *		01-JUL-2004 PW: Added LRAND
 *		17-JUL-2004 PW: Added LOG2
 *		27-JUL-2005 Pw: Added TN (Chebyshev)
 *		10-AUG-2005 PW: Added lower case x and y for normalized -1/+1 coordinates
 *		07-SEP-2005 PW: Added CORRCOEFF
 *		12-JAN-2006 PW: Use size_t for counters to allow > 2 Gb arrays
 *		22-MAR-2006 PW: Added CPOISS
 *		24-MAR-2006 PW: Added ZDIST
 *		06-JUL-2007 PW: Added PSI, PV, QV, COT, COTD, ACOT, SEC, SECD, ASEC, CSC, CSCD, ACSC
 *		21-SEP-2007 PW: Added KURT, SKEW, PQUANT, EULER
 *		28-SEP-2007 RS: Added PLMg, YLMg and set Imaginary component of YL0 = 0.
 *		07-DEC-2007 PW: Added XMIN, XMAX, XINC, NX and same for Y.
 *		13-AUG-2008 PW: Added NOT, fixed BCs for 2nd derivatives when some nodes are NaN
 *		10-JAN-2014 PW: Fixed SDIST to yield km like the other ?DIST operators, and added
 *			         KM2DEG and DEG2KM to help with getting degrees from km, mostly.
 * Version:	4.5
 *
 */

#define _XOPEN_SOURCE
#include "gmt.h"

#define GRDMATH_ARG_IS_OPERATOR		 0
#define GRDMATH_ARG_IS_FILE		-1
#define GRDMATH_ARG_IS_NUMBER		-2
#define GRDMATH_ARG_IS_PI		-3
#define GRDMATH_ARG_IS_E		-4
#define GRDMATH_ARG_IS_EULER		-5
#define GRDMATH_ARG_IS_XMIN		-6
#define GRDMATH_ARG_IS_XMAX		-7
#define GRDMATH_ARG_IS_XINC		-8
#define GRDMATH_ARG_IS_NX		-9
#define GRDMATH_ARG_IS_YMIN		-10
#define GRDMATH_ARG_IS_YMAX		-11
#define GRDMATH_ARG_IS_YINC		-12
#define GRDMATH_ARG_IS_NY		-13
#define GRDMATH_ARG_IS_X_MATRIX		-14
#define GRDMATH_ARG_IS_x_MATRIX		-15
#define GRDMATH_ARG_IS_Y_MATRIX		-16
#define GRDMATH_ARG_IS_y_MATRIX		-17
#define GRDMATH_ARG_IS_ASCIIFILE	-18
#define GRDMATH_ARG_IS_SAVE		-19
#define GRDMATH_N_SKIP_ARGS		 8

#define GRDMATH_STACK_SIZE		100

struct GRDMATH_CTRL {	/* All control options for this program (except common args) */
	/* active is TRUE if the option has been activated */
	struct F {	/* -F */
		GMT_LONG active;
	} F;
	struct I {	/* -Idx[/dy] */
		GMT_LONG active;
		double xinc, yinc;
	} I;
	struct M {	/* -M */
		GMT_LONG active;
	} M;
	struct N {	/* -N */
		GMT_LONG active;
	} N;
};

struct GRDMATH_INFO {
	GMT_LONG nm;
	char ASCII_file[BUFSIZ];
	GMT_LONG convert;		/* Reflects -M */
	float *grd_x, *grd_y;
	float *grd_xn, *grd_yn;
	float *dx, dy;		/* In flat-Earth m if -M is set */
	struct GRD_HEADER header;
};

#include "grdmath_def.h"

/* Helper functions */
void grd_AZ_sub (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG reverse);
void grd_PVQV (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG kind);
void grd_YLM_sub (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG ortho);

int main (int argc, char **argv)
{
	GMT_LONG i, j, k, i64, arg, op = 0, nstack = 0, new_stack = -1, ok = 1, n_items = 0, this_stack;
	GMT_LONG consumed_operands[GRDMATH_N_OPERATORS], produced_operands[GRDMATH_N_OPERATORS];

	GMT_LONG constant[GRDMATH_STACK_SIZE], error = FALSE, skip_arg, set_r = FALSE;

	float *stack[GRDMATH_STACK_SIZE];

	double factor[GRDMATH_STACK_SIZE], value, x_noise, y_noise, off, scale, special_symbol[GRDMATH_ARG_IS_PI-GRDMATH_ARG_IS_NY+1];

	char *outfile = CNULL;
	char *skip_args[GRDMATH_N_SKIP_ARGS] = {"-R", "-I", "-F", "-M", "-N", "-V", "-f", "-b"};	/* 2nd time we want to skip these */

	struct GRD_HEADER grd[GRDMATH_STACK_SIZE];
	struct GMT_HASH *p = NULL, *current = NULL, localhashnode[GRDMATH_N_OPERATORS];
	struct GRDMATH_INFO info;
	struct GRDMATH_CTRL *Ctrl = NULL;

	PFV call_operator[GRDMATH_N_OPERATORS];

	GMT_LONG decode_argument(char *txt, double *value, struct GMT_HASH *H);
	void grdmath_init (PFV ops[], GMT_LONG n_args[], GMT_LONG n_out[]);
	void *New_grdmath_Ctrl (), Free_grdmath_Ctrl (struct GRDMATH_CTRL *C);

	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct GRDMATH_CTRL *) New_grdmath_Ctrl ();		/* Allocate and initialize defaults in a new control structure */

	memset ((void *)&info, 0, sizeof (struct GRDMATH_INFO));

	if (argc == 2 && !strcmp (argv[1], "-")) error = GMT_give_synopsis_and_exit = TRUE;

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "grdmath %s - Reverse Polish Notation (RPN) calculator for grid files (element by element)\n\n", GMT_VERSION);
		fprintf (stderr, "usage: grdmath [%s] [-F] [%s] [-M] [-N] [-V] %s %s\n", GMT_Rgeo_OPT, GMT_I_OPT, GMT_bi_OPT, GMT_f_OPT);
		fprintf (stderr, "	A B op C op D op ... = outfile\n\n");

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\tA, B, etc are grid files, constants, or symbols (see below).\n");
		fprintf (stderr, "\tThe stack can hold up to %d entries (given enough memory).\n", GRDMATH_STACK_SIZE);
		fprintf (stderr, "\tTrigonometric operators expect radians.\n");
		fprintf (stderr, "\tThe operators and number of input and output arguments are:\n\n");
		fprintf (stderr, "\tName    #args   Returns\n");
		fprintf (stderr, "\t-----------------------\n");
#include "grdmath_explain.h"
		fprintf (stderr, "\n\tThe special symbols are:\n\n");
		fprintf (stderr, "\t  PI	= 3.1415926...\n");
		fprintf (stderr, "\t  E	= 2.7182818...\n");
		fprintf (stderr, "\t  EULER	= 0.5772156...\n");
		fprintf (stderr, "\t  XMIN, XMAX, XINC or NX	= the corresponding constants\n");
		fprintf (stderr, "\t  YMIN, YMAX, YINC or NY	= the corresponding constants\n");
		fprintf (stderr, "\t  X	= grid with x-coordinates\n");
		fprintf (stderr, "\t  Y	= grid with y-coordinates\n");
		fprintf (stderr, "\t  Xn	= grid with normalized [-1|+1] x-coordinates\n");
		fprintf (stderr, "\t  Yn	= grid with normalized [-1|+1] y-coordinates\n");
		fprintf (stderr, "\n\tOPTIONS: (only used if no grid files are passed as arguments)\n\n");
		fprintf (stderr, "\t-F Set pixel grid registration [Default is gridline orientation]. Requires -R and -I.\n");
		GMT_inc_syntax ('I', 0);
		fprintf (stderr, "\t-M Handle map units in derivatives.  In this case, dx,dy of grid\n");
		fprintf (stderr, "\t   will be converted from degrees lon,lat into meters (Flat-earth approximation).\n");
		fprintf (stderr, "\t   Default computes derivatives in units of data/grid_distance.\n");
		fprintf (stderr, "\t-N Do not perform strict domain check if several grids are involved\n");
		fprintf (stderr, "\t   [Default checks that domain is within %g * [xinc or yinc] of each other].\n", GMT_SMALL);
		GMT_explain_option ('R');
		GMT_explain_option ('V');
		GMT_explain_option ('i');
		fprintf (stderr, "\t   (Only applies to the input files for operators LDIST, PDIST, and INSIDE).\n");
		GMT_explain_option ('f');
		exit (EXIT_FAILURE);
	}

	for (i = 1, error = TRUE; error && i < argc; i++) if (argv[i][0] == '=' && argv[i][1] == '\0') error = FALSE;
	if (error) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Usage is <operations> = outfile\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	GMT_hash_init (localhashnode, operator, GRDMATH_N_OPERATORS, GRDMATH_N_OPERATORS);

	for (i = 0; i < GRDMATH_STACK_SIZE; i++) {
		constant[i] = FALSE;
		factor[i] = 0.0;
		grd[i].nx = grd[i].ny = 0;
		stack[i] = (float *)NULL;
	}

	/* Get header from one file so we can allocate space */

	GMT_grd_init (&info.header, argc, argv, FALSE);

	for (arg = 1; info.nm == 0 && arg < argc; arg++) {

		if (argv[arg][0] == '-') continue;	/* Command line option or a negative number */
		if (decode_argument (argv[arg], &value, localhashnode) == GRDMATH_ARG_IS_SAVE) {	/* Output case = */
			arg++;	/* Skip the output file name since the file does not exist yet */
			while (argv[arg][0] == '-' && arg < argc) arg++;	/* In case someone placed an option bewteen = and output gridname */
			continue;
		}
		if (decode_argument (argv[arg], &value, localhashnode) != GRDMATH_ARG_IS_FILE) continue;
		if (arg < (argc - 1) && !(strncmp (argv[arg+1], "LDIST", (size_t)5) && strncmp (argv[arg+1], "PDIST", (size_t)5) && strncmp (argv[arg+1], "INSIDE", (size_t)6))) continue;	/* Not a grid file */

		GMT_err_fail (GMT_read_grd_info (argv[arg], &info.header), argv[arg]);

		info.nm = GMT_get_nm (info.header.nx, info.header.ny);
	}

	/* Scan command line for -R, -I, -F, -M, -N,  -f, -b, -V */

	for (arg = 1; arg < argc; arg++) {
		if (argv[arg][0] == '-') {

			switch (argv[arg][1]) {

				case 'R':
					set_r = TRUE;
				case 'V':
				case 'b':
				case 'f':
					error += GMT_parse_common_options (argv[arg], &info.header.x_min, &info.header.x_max, &info.header.y_min, &info.header.y_max);
					break;

				case 'F':
					Ctrl->F.active = TRUE;
					break;

				case 'I':
					Ctrl->I.active = TRUE;
					if (GMT_getinc (&argv[arg][2], &Ctrl->I.xinc, &Ctrl->I.yinc)) {
						GMT_inc_syntax ('I', 1);
						error = TRUE;
					}
					break;
				case 'M':
					Ctrl->M.active = TRUE;
					break;
				case 'N':
					Ctrl->N.active = TRUE;
					break;

			}
		}
	}

	GMT_check_lattice (&Ctrl->I.xinc, &Ctrl->I.yinc, &Ctrl->F.active, &Ctrl->I.active);

	if (info.nm && (set_r || Ctrl->I.active || Ctrl->F.active)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Cannot use -R, -I, [-F] when grid files are specified\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if ((set_r + Ctrl->I.active)%2) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  -R and -I must both be specified\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (Ctrl->I.active && (Ctrl->I.xinc <= 0.0 || Ctrl->I.yinc <= 0.0)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -I option.  Must specify positive increment(s)\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (set_r && Ctrl->I.active) {
		info.header.node_offset = (int)Ctrl->F.active;
		info.header.x_inc = Ctrl->I.xinc;
		info.header.y_inc = Ctrl->I.yinc;
		GMT_RI_prepare (&info.header);	/* Ensure -R -I consistency and set nx, ny */
		GMT_err_fail (GMT_grd_RI_verify (&info.header, 1), "");
		info.nm = GMT_get_nm (info.header.nx, info.header.ny);
		info.header.xy_off = 0.5 * info.header.node_offset;
	}

	if (info.nm == 0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Expression must contain at least one grid file or -R, -I\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	stack[0] = (float *) GMT_memory (VNULL, info.nm, sizeof (float), GMT_program);

	/* Get x and y vectors */

	info.grd_x = (float *) GMT_memory (VNULL, (size_t)info.header.nx, sizeof (float), GMT_program);
	info.grd_y = (float *) GMT_memory (VNULL, (size_t)info.header.ny, sizeof (float), GMT_program);
	info.grd_xn = (float *) GMT_memory (VNULL, (size_t)info.header.nx, sizeof (float), GMT_program);
	info.grd_yn = (float *) GMT_memory (VNULL, (size_t)info.header.ny, sizeof (float), GMT_program);
	for (j = 0; j < info.header.ny; j++) info.grd_y[j] = (float)GMT_j_to_y (j, info.header.y_min, info.header.y_max, info.header.y_inc, info.header.xy_off, info.header.ny);
	for (i = 0; i < info.header.nx; i++) info.grd_x[i] = (float)GMT_i_to_x (i, info.header.x_min, info.header.x_max, info.header.x_inc, info.header.xy_off, info.header.nx);
	off = 0.5 * (info.header.x_max + info.header.x_min);
	scale = 2.0 / (info.header.x_max - info.header.x_min);
	for (i = 0; i < info.header.nx; i++) info.grd_xn[i] = (float)((info.grd_x[i] - off) * scale);
	off = 0.5 * (info.header.y_max + info.header.y_min);
	scale = 2.0 / (info.header.y_max - info.header.y_min);
	for (j = 0; j < info.header.ny; j++) info.grd_yn[j] = (float)((info.grd_y[j] - off) * scale);
	x_noise = GMT_SMALL * info.header.x_inc;	y_noise = GMT_SMALL * info.header.y_inc;
	info.dx = (float *) GMT_memory (VNULL, (size_t)info.header.ny, sizeof (float), GMT_program);

	if (Ctrl->M.active) {	/* Use flat earth distances for gradients */
		for (j = 0; j < info.header.ny; j++) info.dx[j] = (float) (project_info.DIST_M_PR_DEG * info.header.x_inc * cosd (info.grd_y[j]));
		info.dy = (float)(project_info.DIST_M_PR_DEG * info.header.y_inc);
		info.convert = TRUE;
	}
	else {	/* Constant increments in user units */
		for (j = 0; j < info.header.ny; j++) info.dx[j] = (float)info.header.x_inc;
		info.dy = (float)info.header.y_inc;
	}

	grdmath_init (call_operator, consumed_operands, produced_operands);

	nstack = 0;

	special_symbol[GRDMATH_ARG_IS_PI-GRDMATH_ARG_IS_PI] = M_PI;
	special_symbol[GRDMATH_ARG_IS_PI-GRDMATH_ARG_IS_E] = M_E;
	special_symbol[GRDMATH_ARG_IS_PI-GRDMATH_ARG_IS_EULER] = M_EULER;
	special_symbol[GRDMATH_ARG_IS_PI-GRDMATH_ARG_IS_XMIN] = info.header.x_min;
	special_symbol[GRDMATH_ARG_IS_PI-GRDMATH_ARG_IS_XMAX] = info.header.x_max;
	special_symbol[GRDMATH_ARG_IS_PI-GRDMATH_ARG_IS_XINC] = info.header.x_inc;
	special_symbol[GRDMATH_ARG_IS_PI-GRDMATH_ARG_IS_NX] = info.header.nx;
	special_symbol[GRDMATH_ARG_IS_PI-GRDMATH_ARG_IS_YMIN] = info.header.y_min;
	special_symbol[GRDMATH_ARG_IS_PI-GRDMATH_ARG_IS_YMAX] = info.header.y_max;
	special_symbol[GRDMATH_ARG_IS_PI-GRDMATH_ARG_IS_YINC] = info.header.y_inc;
	special_symbol[GRDMATH_ARG_IS_PI-GRDMATH_ARG_IS_NY] = info.header.ny;
	if (gmtdefs.verbose) fprintf (stderr, "%s: ", GMT_program);
	for (arg = 1; !error && arg < argc; arg++) {

		/* First check if we should skip optional arguments */

		for (k = 0, skip_arg = FALSE; !skip_arg && k < GRDMATH_N_SKIP_ARGS; k++) skip_arg = !strncmp (argv[arg], skip_args[k], (size_t)2);
		if (skip_arg) continue;

		op = decode_argument (argv[arg], &value, localhashnode);

		if (op != GRDMATH_ARG_IS_FILE && !GMT_access(argv[arg], R_OK)) fprintf (stderr, "%s Warning: The number or operator %s may be confused with an existing file %s!\n", GMT_program, argv[arg], argv[arg]);

		if (op == GRDMATH_ARG_IS_SAVE) {	/* Time to save the current stack to output and pop the stack */
			if (nstack <= 0 && n_items) {
				fprintf (stderr, "%s: GMT SYNTAX ERROR:  No items left on the stack for output!\n", GMT_program);
				exit (EXIT_FAILURE);
			}
			arg++;	/* Skip to the file name */
			outfile = argv[arg];

			if (gmtdefs.verbose) fprintf (stderr, "= %s", outfile);

			GMT_grd_init (&info.header, argc, argv, TRUE);

			if (n_items && new_stack < 0 && constant[nstack-1]) {	/* Only a constant provided, set grid accordingly */
				for (i64 = 0; i64 < info.nm; i64++) stack[nstack-1][i64] = (float)factor[nstack-1];
			}
			this_stack = (n_items) ? nstack - 1 : 0;	/* Since giving no args just go get an empty grid is OK */
			GMT_err_fail (GMT_write_grd (outfile, &info.header, stack[this_stack], 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), outfile);
			if (n_items) nstack--;	/* Pop off the current stack if there is one */
			new_stack = nstack;
			continue;
		}
		if (op < GRDMATH_ARG_IS_OPERATOR) {	/* File name or factor */

			if (op == GRDMATH_ARG_IS_FILE && !(strncmp (argv[arg+1], "LDIST", (size_t)5) && strncmp (argv[arg+1], "PDIST", (size_t)5) && strncmp (argv[arg+1], "INSIDE", (size_t)6))) op = GRDMATH_ARG_IS_ASCIIFILE;

			if (nstack == GRDMATH_STACK_SIZE) {	/* Stack overflow */
				error = TRUE;
				continue;
			}
			n_items++;
			if (op == GRDMATH_ARG_IS_NUMBER) {
				constant[nstack] = TRUE;
				factor[nstack] = value;
				error = FALSE;
				if (gmtdefs.verbose) fprintf (stderr, "%g ", factor[nstack]);
				nstack++;
				continue;
			}
			else if (op <= GRDMATH_ARG_IS_PI && op >= GRDMATH_ARG_IS_NY) {
				constant[nstack] = TRUE;
				factor[nstack] = special_symbol[GRDMATH_ARG_IS_PI-op];
				if (gmtdefs.verbose) fprintf (stderr, "%g ", factor[nstack]);
				nstack++;
				continue;
			}

			/* Here we need a matrix */

			GMT_grd_init (&grd[nstack], argc, argv, TRUE);
			if (!stack[nstack]) stack[nstack] = (float *) GMT_memory (VNULL, info.nm, sizeof (float), GMT_program);
			constant[nstack] = FALSE;

			if (op == GRDMATH_ARG_IS_X_MATRIX) {		/* Need to set up matrix of x-values */
				if (gmtdefs.verbose) fprintf (stderr, "X ");
				for (j = k = 0; j < info.header.ny; j++, k += info.header.nx)
					memcpy ((void *)&stack[nstack][k], (void *)info.grd_x, (size_t)(info.header.nx * sizeof (float)));
			}
			else if (op == GRDMATH_ARG_IS_x_MATRIX) {		/* Need to set up matrix of normalized x-values */
				if (gmtdefs.verbose) fprintf (stderr, "Xn ");
				for (j = k = 0; j < info.header.ny; j++, k += info.header.nx)
					memcpy ((void *)&stack[nstack][k], (void *)info.grd_xn, (size_t)(info.header.nx * sizeof (float)));
			}
			else if (op == GRDMATH_ARG_IS_Y_MATRIX) {	/* Need to set up matrix of y-values */
				if (gmtdefs.verbose) fprintf (stderr, "Y ");
				for (j = k = 0; j < info.header.ny; j++) for (i = 0; i < info.header.nx; i++, k++)
					stack[nstack][k] = info.grd_y[j];
			}
			else if (op == GRDMATH_ARG_IS_y_MATRIX) {	/* Need to set up matrix of normalized y-values */
				if (gmtdefs.verbose) fprintf (stderr, "Yn ");
				for (j = k = 0; j < info.header.ny; j++) for (i = 0; i < info.header.nx; i++, k++)
					stack[nstack][k] = info.grd_yn[j];
			}
			else if (op == GRDMATH_ARG_IS_ASCIIFILE) {
				strcpy (info.ASCII_file, argv[arg]);
				if (gmtdefs.verbose) fprintf (stderr, "(%s) ", argv[arg]);
			}
			else if (op == GRDMATH_ARG_IS_FILE) {		/* Filename given */
				if (gmtdefs.verbose) fprintf (stderr, "%s ", argv[arg]);
				GMT_err_fail (GMT_read_grd_info (argv[arg], &grd[nstack]), argv[arg]);
				if (grd[nstack].nx != info.header.nx || grd[nstack].ny != info.header.ny) {
					fprintf (stderr, "%s: grid files not of same size!\n", GMT_program);
					exit (EXIT_FAILURE);
				}
				else if (!Ctrl->N.active && (fabs (grd[nstack].x_min - info.header.x_min) > x_noise || fabs (grd[nstack].x_max - info.header.x_max) > x_noise ||
					fabs (grd[nstack].y_min - info.header.y_min) > y_noise || fabs (grd[nstack].y_max - info.header.y_max) > y_noise)) {
					fprintf (stderr, "%s: grid files do not cover the same area!\n", GMT_program);
					exit (EXIT_FAILURE);
				}
				GMT_err_fail (GMT_read_grd (argv[arg], &grd[nstack], stack[nstack], 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), argv[arg]);
			}
			nstack++;
			continue;
		}

		/* Here we have an operator */

		if ((new_stack = nstack - consumed_operands[op] + produced_operands[op]) >= GRDMATH_STACK_SIZE) {
			error = TRUE;
			continue;
		}

		if (nstack < consumed_operands[op]) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR:  Operation \"%s\" requires %ld operands\n", GMT_program, operator[op], consumed_operands[op]);
			exit (EXIT_FAILURE);
		}

		n_items++;
		if (gmtdefs.verbose) fprintf (stderr, "%s ", operator[op]);

		for (i = produced_operands[op] - consumed_operands[op]; i > 0; i--) {
			if (stack[nstack+i-1])	continue;

			/* Must make space for more */

			stack[nstack+i-1] = (float *) GMT_memory (VNULL, info.nm, sizeof (float), GMT_program);
		}

		/* If operators operates on constants only we may have to make space as well */

		for (j = 0, i = nstack - consumed_operands[op]; j < produced_operands[op]; j++, i++) {
			if (constant[i] && !stack[i]) stack[i] = (float *) GMT_memory (VNULL, info.nm, sizeof (float), GMT_program);
		}

		(*call_operator[op]) (&info, stack, constant, factor, nstack - 1);	/* Do it */

		nstack = new_stack;

		for (i = 1; i <= produced_operands[op]; i++)
			constant[nstack-i] = FALSE;	/* Now filled with grid */
	}

	if (error && !ok) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Unable to decode constant %s (File not found?)\n", GMT_program, argv[i-1]);
		exit (EXIT_FAILURE);
	}

	if (error) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Stack overflow (%s)\n", GMT_program, argv[i-1]);
		exit (EXIT_FAILURE);
	}

	for (i = 0; i < GRDMATH_STACK_SIZE; i++) if (stack[i]) GMT_free ((void *)stack[i]);
	GMT_free ((void *)info.grd_x);
	GMT_free ((void *)info.grd_y);
	GMT_free ((void *)info.grd_xn);
	GMT_free ((void *)info.grd_yn);
	GMT_free ((void *)info.dx);
	for (i = 0; i < GRDMATH_N_OPERATORS; i++) {
		p = localhashnode[i].next;
		while ((current = p)) {
			p = p->next;
			GMT_free ((void *)current);
		}
	}

	if (gmtdefs.verbose) fprintf (stderr, "\n");

	if (nstack > 0) fprintf (stderr, "%s: Warning: %ld more operands left on the stack!\n", GMT_program, nstack);

	Free_grdmath_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

/* -----------------------------------------------------------------
 *              Definitions of all operator functions
 * -----------------------------------------------------------------*/

void grd_ABS (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: ABS 1 1 abs (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand == 0!\n", GMT_program);
	if (constant[last]) a = fabs (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : fabs ((double)stack[last][i]));
}

void grd_ACOS (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: ACOS 1 1 acos (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && fabs (factor[last]) > 1.0) fprintf (stderr, "%s: Warning, |operand| > 1 for ACOS!\n", GMT_program);
	if (constant[last]) a = d_acos (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : d_acos ((double)stack[last][i]));
}

void grd_ACOSH (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: ACOSH 1 1 acosh (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && fabs (factor[last]) > 1.0) fprintf (stderr, "%s: Warning, operand < 1 for ACOSH!\n", GMT_program);
	if (constant[last]) a = acosh (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : acosh ((double)stack[last][i]));
}

void grd_ACOT (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: ACOT 1 1 acot (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && fabs (factor[last]) > 1.0) fprintf (stderr, "%s: Warning, |operand| > 1 for ACOT!\n", GMT_program);
	if (constant[last]) a = atan (1.0 / factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : atan ((double)(1.0 / stack[last][i])));
}

void grd_ACSC (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: ACSC 1 1 acsc (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && fabs (factor[last]) > 1.0) fprintf (stderr, "%s: Warning, |operand| > 1 for ACSC!\n", GMT_program);
	if (constant[last]) a = d_asin (1.0 / factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : d_asin ((double)(1.0 / stack[last][i])));
}

void grd_ADD (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: ADD 2 1 A + B.  */
{
	GMT_LONG prev, i;
	double a, b;

	prev = last - 1;
	if (gmtdefs.verbose && constant[prev] && factor[prev] == 0.0) fprintf (stderr, "%s: Warning, operand one == 0!\n", GMT_program);
	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand two == 0!\n", GMT_program);
	for (i = 0; i < info->nm; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][i];
		b = (constant[last]) ? factor[last] : stack[last][i];
		stack[prev][i] = (float)(a + b);
	}
}

void grd_AND (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: AND 2 1 NaN if A and B == NaN, B if A == NaN, else A.  */
{
	GMT_LONG prev, i;
	double a, b;

	prev = last - 1;
	for (i = 0; i < info->nm; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][i];
		b = (constant[last]) ? factor[last] : stack[last][i];
		stack[prev][i] = (float)((GMT_is_dnan (a)) ? b : a);
	}
}

void grd_ASEC (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: ASEC 1 1 asec (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && fabs (factor[last]) > 1.0) fprintf (stderr, "%s: Warning, |operand| > 1 for ASEC!\n", GMT_program);
	if (constant[last]) a = d_acos (1.0 / factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : d_acos ((double)(1.0 / stack[last][i])));
}

void grd_ASIN (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: ASIN 1 1 asin (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && fabs (factor[last]) > 1.0) fprintf (stderr, "%s: Warning, |operand| > 1 for ASIN!\n", GMT_program);
	if (constant[last]) a = d_asin (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : d_asin ((double)stack[last][i]));
}

void grd_ASINH (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: ASINH 1 1 asinh (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = asinh (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : asinh ((double)stack[last][i]));
}

void grd_ATAN (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: ATAN 1 1 atan (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = atan (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : atan ((double)stack[last][i]));
}

void grd_ATAN2 (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: ATAN2 2 1 atan2 (A, B).  */
{
	GMT_LONG prev, i;
	double a, b;

	prev = last - 1;
	if (gmtdefs.verbose && constant[prev] && factor[prev] == 0.0) fprintf (stderr, "%s: Warning, operand one == 0 for ATAN2!\n", GMT_program);
	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand two == 0 for ATAN2!\n", GMT_program);
	for (i = 0; i < info->nm; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][i];
		b = (constant[last]) ? factor[last] : stack[last][i];
		stack[prev][i] = (float)d_atan2 (a, b);
	}
}

void grd_ATANH (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: ATANH 1 1 atanh (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && fabs (factor[last]) >= 1.0) fprintf (stderr, "%s: Warning, |operand| >= 1 for ATANH!\n", GMT_program);
	if (constant[last]) a = atanh (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : atanh ((double)stack[last][i]));
}

void grd_BEI (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: BEI 1 1 bei (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = GMT_bei (fabs (factor[last]));
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : GMT_bei (fabs((double)stack[last][i])));
}

void grd_BER (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: BER 1 1 ber (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = GMT_ber (fabs (factor[last]));
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : GMT_ber (fabs ((double)stack[last][i])));
}

void grd_CAZ (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: CAZ 2 1 Cartesian azimuth from grid nodes to stack x,y.  */
{
	GMT_LONG i, j, prev, k;
	double x, y, az;

	prev = last - 1;
	for (k = j = i = 0; k < info->nm; k++) {
		x = (constant[prev]) ? factor[prev] : stack[prev][k];
		y = (constant[last]) ? factor[last] : stack[last][k];
		az = 90.0 - atan2d (y - (double)info->grd_y[j], x - (double)info->grd_x[i]);
		while (az < -180.0) az += 360.0;
		while (az > +180.0) az -= 360.0;
		stack[prev][k] = (float)az;
		i++;
		if (i == info->header.nx) i = 0, j++;
	}
}

void grd_CBAZ (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: CBAZ 2 1 Cartesian backazimuth from grid nodes to stack x,y.  */
{
	GMT_LONG i, j, prev, k;
	double x, y, az;

	prev = last - 1;
	for (k = j = i = 0; k < info->nm; k++) {
		x = (constant[prev]) ? factor[prev] : stack[prev][k];
		y = (constant[last]) ? factor[last] : stack[last][k];
		az = 270.0 - atan2d (y - (double)info->grd_y[j], x - (double)info->grd_x[i]);
		while (az < -180.0) az += 360.0;
		while (az > +180.0) az -= 360.0;
		stack[prev][k] = (float)az;
		i++;
		if (i == info->header.nx) i = 0, j++;
	}
}

void grd_CDIST (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: CDIST 2 1 Cartesian distance between grid nodes and stack x,y.  */
{
	GMT_LONG i, j, prev, k;
	double a, b;

	prev = last - 1;
	for (k = j = i = 0; k < info->nm; k++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][k];
		b = (constant[last]) ? factor[last] : stack[last][k];
		stack[prev][k] = (float)hypot (a - (double)info->grd_x[i], b - (double)info->grd_y[j]);
		i++;
		if (i == info->header.nx) i = 0, j++;
	}
}

void grd_CEIL (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: CEIL 1 1 ceil (A) (smallest integer >= A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = ceil (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : ceil ((double)stack[last][i]));
}

void grd_CHICRIT (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: CHICRIT 2 1 Critical value for chi-squared-distribution, with alpha = A and n = B.  */
{
	GMT_LONG prev, i;
	double a, b;

	prev = last - 1;
	if (gmtdefs.verbose && constant[prev] && factor[prev] == 0.0) fprintf (stderr, "%s: Warning, operand one == 0 for CHICRIT!\n", GMT_program);
	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand two == 0 for CHICRIT!\n", GMT_program);
	for (i = 0; i < info->nm; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][i];
		b = (constant[last]) ? factor[last] : stack[last][i];
		stack[prev][i] = (float)GMT_chi2crit (a, b);
	}
}

void grd_CHIDIST (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: CHIDIST 2 1 chi-squared-distribution P(chi2,n), with chi2 = A and n = B.  */
{
	GMT_LONG prev, i;
	double a, b, prob;

	prev = last - 1;
	if (gmtdefs.verbose && constant[prev] && factor[prev] == 0.0) fprintf (stderr, "%s: Warning, operand one == 0 for CHIDIST!\n", GMT_program);
	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand two == 0 for CHIDIST!\n", GMT_program);
	for (i = 0; i < info->nm; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][i];
		b = (constant[last]) ? factor[last] : stack[last][i];
		GMT_chi2 (a, b, &prob);
		stack[prev][i] = (float)prob;
	}
}

void grd_CORRCOEFF (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: CORRCOEFF 2 1 Correlation coefficient r(A, B).  */
{
	GMT_LONG prev, i;
	double coeff;

	prev = last - 1;
	if (constant[prev] || constant[last]) {
		fprintf (stderr, "%s: Error, cannot have constant operands for CORRCOEFF!\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	coeff = GMT_corrcoeff_f (stack[prev], stack[last], (GMT_LONG)info->nm, 0);
	for (i = 0; i < info->nm; i++) stack[prev][i] = (float)coeff;
}

void grd_COS (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: COS 1 1 cos (A) (A in radians).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = cos (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : cos ((double)stack[last][i]));
}

void grd_COSD (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: COSD 1 1 cos (A) (A in degrees).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = cosd (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : cosd ((double)stack[last][i]));
}

void grd_COSH (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: COSH 1 1 cosh (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = cosh (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : cosh ((double)stack[last][i]));
}

void grd_COT (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: COT 1 1 cot (A) (A in radians).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = (1.0 / tan (factor[last]));
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : (1.0 / tan ((double)stack[last][i])));
}

void grd_COTD (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: COTD 1 1 cot (A) (A in degrees).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = (1.0 / tand (factor[last]));
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : (1.0 / tand ((double)stack[last][i])));
}

void grd_CPOISS (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: CPOISS 2 1 Cumulative Poisson distribution F(x,lambda), with x = A and lambda = B.  */
{
	GMT_LONG prev, i;
	double a, b, prob;

	prev = last - 1;
	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand two == 0 for CHIDIST!\n", GMT_program);
	for (i = 0; i < info->nm; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][i];
		b = (constant[last]) ? factor[last] : stack[last][i];
		GMT_cumpoisson (a, b, &prob);
		stack[prev][i] = (float)prob;
	}
}

void grd_CSC (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: CSC 1 1 csc (A) (A in radians).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = (1.0 / sin (factor[last]));
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : (1.0 / sin ((double)stack[last][i])));
}

void grd_CSCD (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: CSCD 1 1 csc (A) (A in degrees).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = (1.0 / sind (factor[last]));
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : (1.0 / sind ((double)stack[last][i])));
}

void grd_CURV (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: CURV 1 1 Curvature of A (Laplacian).  */
{
	GMT_LONG i, j, k, nx;
	double cy;
	float *z, *cx;

	/* Curvature (Laplacian) */

	if (constant[last]) {
		if (gmtdefs.verbose) fprintf (stderr, "%s: Warning, operand to CURV is constant!\n", GMT_program);
		memset ((void *)stack[last], 0, info->nm * sizeof (float));
		return;
	}

	z = (float *) (float *) GMT_memory (VNULL, info->nm, sizeof (float), GMT_program);
	cx = (float *) (float *) GMT_memory (VNULL, (size_t)info->header.ny, sizeof (float), GMT_program);
	for (j = 0; j < info->header.ny; j++) cx[j] = (float)1.0 / (info->dx[j] * info->dx[j]);

	nx = (size_t)info->header.nx;
	cy = 1.0 / (info->dy * info->dy);

	/* First left/right (here d2/dx2 == 0 by way of BC) */

	for (j = 1, k = nx; j < info->header.ny-1; j++, k += nx)
		z[k] = (float)(cy * (stack[last][k+nx] - 2.0 * stack[last][k] + stack[last][k-nx]));
	for (j = 1, k = 2*nx-1; j < info->header.ny-1; j++, k += nx)
		z[k] = (float)(cy * (stack[last][k+nx] - 2.0 * stack[last][k] + stack[last][k-nx]));

	/* Then top/bottom (here d2/dy2 == 0 by way of BC) */

	for (i = k = 1, j = 0; i < info->header.nx - 1; i++, k++)
		z[k] = (float)(cx[j] * (stack[last][k+1] - 2.0 * stack[last][k] + stack[last][k-1]));
	for (i = 1, k = info->nm - nx + 1, j = info->header.ny-1; i < info->header.nx - 1; i++, k++)
		z[k] = (float)(cx[j] * (stack[last][k+1] - 2.0 * stack[last][k] + stack[last][k-1]));

	/* Then inside */

	for (j = 1, k = nx; j < info->header.ny-1; j++) {
		k++;
		for (i = 1; i < info->header.nx-1; i++, k++) {
			z[k] = (float)(cx[j] * (stack[last][k+1] - 2.0 * stack[last][k] + stack[last][k-1]) + cy * (stack[last][k+nx] - 2.0 * stack[last][k] + stack[last][k-nx]));
		}
		k++;
	}

	/* The 4 corners have d2/dxy = 0 as BC (as initialized by GMT_memory).
	   However, if neighbors are NaN then the corner should be NaN as well */
	
	/* Top Left: k = 0*/
	k = 0;
	if (GMT_is_fnan (stack[last][k]) || GMT_is_fnan (stack[last][k+1]) || GMT_is_fnan (stack[last][k+nx]) || GMT_is_fnan (stack[last][k+nx+1])) z[k] = GMT_f_NaN;
	/* Top Right: k = nx-1 */
	k = nx-1;
	if (GMT_is_fnan (stack[last][k]) || GMT_is_fnan (stack[last][k-1]) || GMT_is_fnan (stack[last][k+nx]) || GMT_is_fnan (stack[last][k+nx-1])) z[k] = GMT_f_NaN;
	/* Bottom Left: k = nx-1 */
	k = (info->header.ny-1) * nx;
	if (GMT_is_fnan (stack[last][k]) || GMT_is_fnan (stack[last][k+1]) || GMT_is_fnan (stack[last][k-nx]) || GMT_is_fnan (stack[last][k-nx+1])) z[k] = GMT_f_NaN;
	/* Top Right: k = nx-1 */
	k = info->header.ny * nx - 1;
	if (GMT_is_fnan (stack[last][k]) || GMT_is_fnan (stack[last][k-1]) || GMT_is_fnan (stack[last][k-nx]) || GMT_is_fnan (stack[last][k-nx-1])) z[k] = GMT_f_NaN;

	memcpy ((void *)stack[last], (void *)z, info->nm * sizeof (float));
	GMT_free ((void *)z);
	GMT_free ((void *)cx);
}

void grd_D2DX2 (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: D2DX2 1 1 d^2(A)/dx^2 2nd derivative.  */
{
	GMT_LONG i, j, k;
	double c, left, next_left;

	/* Central 2nd difference in x */

	if (constant[last]) {
		if (gmtdefs.verbose) fprintf (stderr, "%s: Warning, operand to D2DX2 is constant!\n", GMT_program);
		memset ((void *)stack[last], 0, info->nm * sizeof (float));
		return;
	}

	for (j = k = 0; j < info->header.ny; j++) {
		c = 1.0 / (info->dx[j] * info->dx[j]);
		next_left = stack[last][k];
		stack[last][k] = (GMT_is_fnan (stack[last][k]) || GMT_is_fnan (stack[last][k+1])) ? GMT_f_NaN : (float)0.0;	/* Natural BC has curvature = 0 at edges unless we have NaNs */
		k++;
		for (i = 1; i < info->header.nx-1; i++, k++) {
			left = next_left;
			next_left = stack[last][k];
			stack[last][k] = (float)(c * (stack[last][k+1] - 2.0 * stack[last][k] + left));
		}
		stack[last][k] = (GMT_is_fnan (stack[last][k]) || GMT_is_fnan (stack[last][k-1])) ? GMT_f_NaN : (float)0.0;	/* Natural BC has curvature = 0 at edges unless we have NaNs */
		k++;
	}
}

void grd_D2DY2 (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: D2DY2 1 1 d^2(A)/dy^2 2nd derivative.  */
{
	GMT_LONG i, j, k, nx;
	double c, bottom, next_bottom;

	/* Central 2nd difference in y */

	if (constant[last]) {
		if (gmtdefs.verbose) fprintf (stderr, "%s: Warning, operand to D2DY2 is constant!\n", GMT_program);
		memset ((void *)stack[last], 0, info->nm * sizeof (float));
		return;
	}

	c = 1.0 / (info->dy * info->dy);
	nx = (size_t)info->header.nx;
	for (i = 0; i < info->header.nx; i++) {
		k = i;
		next_bottom = stack[last][k];
		stack[last][k] = (GMT_is_fnan (stack[last][k]) || GMT_is_fnan (stack[last][k+nx])) ? GMT_f_NaN : (float)0.0;	/* Natural BC has curvature = 0 at edges unless we have NaNs */
		k += nx;
		for (j = 1; j < info->header.ny-1; j++, k += nx) {
			bottom = next_bottom;
			next_bottom = stack[last][k];
			stack[last][k] = (float)(c * (stack[last][k+nx] - 2 * stack[last][k] + bottom));
		}
		stack[last][k] = (GMT_is_fnan (stack[last][k]) || GMT_is_fnan (stack[last][k-nx])) ? GMT_f_NaN : (float)0.0;	/* Natural BC has curvature = 0 at edges unless we have NaNs */
	}
}

void grd_D2DXY (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: D2DXY 1 1 d^2(A)/dxdy 2nd derivative.  */
{
	GMT_LONG i, j, k, nx;
	double *cx, cy, cxy;
	float *z;

	/* Cross derivative d2/dxy = d2/dyx  */

	if (constant[last]) {
		if (gmtdefs.verbose) fprintf (stderr, "%s: Warning, operand to D2DXY is constant!\n", GMT_program);
		memset ((void *)stack[last], 0, info->nm * sizeof (float));
		return;
	}

	z = (float *) (float *) GMT_memory (VNULL, info->nm, sizeof (float), GMT_program);
	cx = (double *) (double *) GMT_memory (VNULL, (size_t)info->header.ny, sizeof (double), GMT_program);
	for (j = 0; j < info->header.ny; j++) cx[j] = 0.5 / (info->dx[j] * info->dx[j]);

	nx = (size_t)info->header.nx;
	cy = 0.5 / (info->dy * info->dy);

	/* First left/right using the BC that d2/dx2 = 0 */

	for (j = 1, k = nx; j < info->header.ny-1; j++, k += nx)
		z[k] = (float)(2.0 * cx[j]*cy * (stack[last][k-nx+1] - stack[last][k-nx] + stack[last][k+nx] - stack[last][k+nx+1]));
	for (j = 1, k = 2*nx-1; j < info->header.ny-1; j++, k += nx)
		z[k] = (float)(2.0 * cx[j]*cy * (stack[last][k-nx] - stack[last][k-nx-1] + stack[last][k+nx-1] - stack[last][k+nx]));

	/* Then top/bottom  using the BC that d2/dy2 = 0 */

	cxy = cx[0] * cy;
	for (i = k = 1; i < info->header.nx - 1; i++, k++)
		z[k] = (float)(2.0 * cxy * (stack[last][k+1] - stack[last][k+nx+1] + stack[last][k+nx-1] - stack[last][k-1]));
	cxy = cx[info->header.ny-1] * cy;
	for (i = 1, k = info->nm - nx + 1; i < info->header.nx - 1; i++, k++)
		z[k] = (float)(2.0 * cxy * (stack[last][k+1-nx] - stack[last][k-nx-1] + stack[last][k-1] - stack[last][k+1]));

	/* The 4 corners have d2/dxy = 0 as BC (as initialized by GMT_memory).
	   However, if neighbors are NaN then the corner should be NaN as well */
	
	/* Top Left: k = 0*/
	k = 0;
	if (GMT_is_fnan (stack[last][k]) || GMT_is_fnan (stack[last][k+1]) || GMT_is_fnan (stack[last][k+nx]) || GMT_is_fnan (stack[last][k+nx+1])) z[k] = GMT_f_NaN;
	/* Top Right: k = nx-1 */
	k = nx-1;
	if (GMT_is_fnan (stack[last][k]) || GMT_is_fnan (stack[last][k-1]) || GMT_is_fnan (stack[last][k+nx]) || GMT_is_fnan (stack[last][k+nx-1])) z[k] = GMT_f_NaN;
	/* Bottom Left: k = nx-1 */
	k = (info->header.ny-1) * nx;
	if (GMT_is_fnan (stack[last][k]) || GMT_is_fnan (stack[last][k+1]) || GMT_is_fnan (stack[last][k-nx]) || GMT_is_fnan (stack[last][k-nx+1])) z[k] = GMT_f_NaN;
	/* Top Right: k = nx-1 */
	k = info->header.ny * nx - 1;
	if (GMT_is_fnan (stack[last][k]) || GMT_is_fnan (stack[last][k-1]) || GMT_is_fnan (stack[last][k-nx]) || GMT_is_fnan (stack[last][k-nx-1])) z[k] = GMT_f_NaN;
	
	/* Then inside */

	for (j = 1, k = nx; j < info->header.ny-1; j++) {
		k++;
		for (i = 1; i < info->header.nx-1; i++, k++) {
			z[k] = (float)(cx[j]*cy * (stack[last][k-nx+1] - stack[last][k-nx-1] + stack[last][k+nx-1] - stack[last][k+nx+1]));
		}
		k++;
	}

	memcpy ((void *)stack[last], (void *)z, info->nm * sizeof (float));
	GMT_free ((void *)z);
	GMT_free ((void *)cx);
}

void grd_D2R (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: D2R 1 1 Converts Degrees to Radians.  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = factor[last] * D2R;
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : (stack[last][i] * D2R));
}

void grd_DDX (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: DDX 1 1 d(A)/dx Central 1st derivative.  */
{
	GMT_LONG i, j, k;
	double c, left, next_left;

	/* Central 1st difference in x */

	if (constant[last]) {
		if (gmtdefs.verbose) fprintf (stderr, "%s: Warning, operand to DDX is constant!\n", GMT_program);
		memset ((void *)stack[last], 0, info->nm * sizeof (float));
		return;
	}

	for (j = k = 0; j < info->header.ny; j++) {
		c = 0.5 / info->dx[j];
		next_left = 2.0 * stack[last][k] - stack[last][k+1];
		for (i = 0; i < info->header.nx-1; i++, k++) {
			left = next_left;
			next_left = stack[last][k];
			stack[last][k] = (float)(c * (stack[last][k+1] - left));
		}
		stack[last][k] = (float)(2.0 * c * (stack[last][k] - next_left));
		k++;
	}
}

void grd_DDY (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: DDY 1 1 d(A)/dy Central 1st derivative.  */
{
	GMT_LONG i, j, k, nx;
	double c, bottom, next_bottom;

	/* Central 1st difference in y */

	if (constant[last]) {
		if (gmtdefs.verbose) fprintf (stderr, "%s: Warning, operand to DDY is constant!\n", GMT_program);
		memset ((void *)stack[last], 0, info->nm * sizeof (float));
		return;
	}

	c = -0.5 / info->dy;	/* Because the loop over j below goes from ymax to ymin we compensate with a minus sign here */
	nx = (size_t)info->header.nx;
	for (i = 0; i < info->header.nx; i++) {
		k = i;
		next_bottom = 2.0 * stack[last][k] - stack[last][k+nx];
		for (j = 0; j < info->header.ny - 1; j++, k += nx) {
			bottom = next_bottom;
			next_bottom = stack[last][k];
			stack[last][k] = (float)(c * (stack[last][k+nx] - bottom));
		}
		stack[last][k] = (float)(2.0 * c * (stack[last][k] - next_bottom));
	}
}

void grd_DEG2KM (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: DEG2KM 1 1 Converts Spherical Degrees to Kilometers.  */
{
	GMT_LONG i;
	double a = 0.0;

	if (GMT_io.in_col_type[0] & GMT_IS_GEO) {	/* Geographic data */
		if (!GMT_IS_SPHERICAL) fprintf (stderr, "%s: Warning, DEG2KM is only exact when ELLIPSOID == sphere\n", GMT_program);
	}
	else
		fprintf (stderr, "%s: Warning, DEG2KM used with Cartesian data\n", GMT_program);
	if (constant[last]) a = factor[last] * project_info.DIST_KM_PR_DEG;
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : (stack[last][i] * project_info.DIST_KM_PR_DEG));
}

void grd_DILOG (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: DILOG 1 1 dilog (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = GMT_dilog (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : GMT_dilog (stack[last][i]));
}

void grd_DIV (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: DIV 2 1 A / B.  */
{
	GMT_LONG prev, i;
	double a, b;

	prev = last - 1;

	if (constant[last] && factor[last] == 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Cannot divide by zero\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (constant[last]) {	/* Turn divide into multiply */
		a = factor[last];	/* Save original factor */
		factor[last] = 1.0 / factor[last];
		grd_MUL (info, stack, constant, factor, last);
		factor[last] = a;	/* Restore factor */
		return;
	}

	for (i = 0; i < info->nm; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][i];
		b = (constant[last]) ? factor[last] : stack[last][i];
		stack[prev][i] = (float)(a / b);
	}
}

void grd_DUP (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: DUP 1 2 Places duplicate of A on the stack.  */
{
	GMT_LONG next, i;

	next = last + 1;
	constant[next] = constant[last];
	factor[next] = factor[last];
	if (constant[last]) {	/* Time to fess up */
		for (i = 0; i < info->nm; i++) stack[last][i] = (float)factor[last];
	}

	memcpy ((void *)stack[next], (void *)stack[last], info->nm * sizeof (float));
}

void grd_ERF (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: ERF 1 1 Error function erf (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = erf (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : erf ((double)stack[last][i]));
}

void grd_ERFC (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: ERFC 1 1 Complementary Error function erfc (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = erfc (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : erfc ((double)stack[last][i]));
}

void grd_EQ (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: EQ 2 1 1 if A == B, else 0.  */
{
	GMT_LONG prev, i;
	double a, b;

	prev = last - 1;
	for (i = 0; i < info->nm; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][i];
		b = (constant[last]) ? factor[last] : stack[last][i];
		stack[prev][i] = (GMT_is_dnan (a) || GMT_is_dnan (b)) ? GMT_f_NaN : (float)(a == b);
	}
}

void grd_ERFINV (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: ERFINV 1 1 Inverse error function of A.  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = GMT_erfinv (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : GMT_erfinv ((double)stack[last][i]));
}

void grd_EXCH (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: EXCH 2 2 Exchanges A and B on the stack.  */
{
	GMT_LONG prev, i;

	prev = last - 1;
	for (i = 0; i < info->nm; i++) {
		if (constant[prev]) stack[prev][i] = (float)factor[prev];
		if (constant[last]) stack[last][i] = (float)factor[last];
		f_swap (stack[last][i], stack[prev][i]);
	}
	d_swap (factor[last], factor[prev]);
	l_swap (constant[last], constant[prev]);
}

void grd_EXP (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: EXP 1 1 exp (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = exp (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : exp ((double)stack[last][i]));
}

void grd_FACT (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: FACT 1 1 A! (A factorial).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = GMT_factorial ((GMT_LONG)irint(factor[last]));
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : GMT_factorial ((GMT_LONG)irint(stack[last][i])));
}

void grd_EXTREMA (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: EXTREMA 1 1 Local Extrema: +2/-2 is max/min, +1/-1 is saddle with max/min in x, 0 elsewhere.  */
{
	GMT_LONG i, j, nx1, ny1, dx, dy, diag, k;
	GMT_LONG do_derivative (float *z, GMT_LONG this_node, GMT_LONG this_pos, GMT_LONG last_i, GMT_LONG last_j, GMT_LONG type);
	GMT_LONG do_diagonal (float *z, GMT_LONG this_node, GMT_LONG this_i, GMT_LONG this_j, GMT_LONG last_i, GMT_LONG last_j, GMT_LONG type);
	float *z;

	/* Find local extrema in grid */

	if (constant[last]) {
		if (gmtdefs.verbose) fprintf (stderr, "%s: Warning, operand to EXTREMA is constant!\n", GMT_program);
		memset ((void *)stack[last], 0, info->nm * sizeof (float));
		return;
	}

	z = (float *) GMT_memory (VNULL, info->nm, sizeof (float), GMT_program);

	nx1 = info->header.nx - 1;
	ny1 = info->header.ny - 1;

	/* We will visit each node on the grid and determine if there are extrema.  We do this
	 * by looking at the along-x and along-y profiles separately.  If both of them shows the
	 * central node in a min or max location with respect to its two neighbors (one if at the
	 * edge of the grid or in the presence of NaNs) we must do further checking.  If the two
	 * extrema are of different sign we call it a saddle point and we are done.  If both are
	 * min or max then we look at the two diagonal lines to see if they can confirm this.
	 * Here, NaNs are not held against you - it takes real values to overrule the dx,dy values.
	 *
	 * Min is given -2, Max is given +2, Saddle points are +1 (if max in x) or -1 (if max in y)
	 * Default is 0 which means no extrema.
	 */

	for (j = k = 0; j < info->header.ny; j++) {	/* k is middle (current) node */
		for (i = 0; i < info->header.nx; i++, k++) {

			if (GMT_is_fnan (stack[last][k])) continue;	/* No extrema if point is NaN */

			if ((dx = do_derivative (stack[last], k, i, nx1, info->header.nx, 0)) == -2) continue;	/* Too many NaNs or flat x-line */
			if ((dy = do_derivative (stack[last], k, j, ny1, info->header.nx, 1)) == -2) continue;	/* Too many NaNs or flat y-line */

			if ((dx * dy) == 0) continue;	/* No min or max possible */
			if ((dx * dy) < 0) {	/* Saddle point - don't need to check diagonals */
				z[k] = (float)((dx > 0) ? +1 : -1);
				continue;
			}

			/* Need to examine diagonal trends to verify min or max */

			if ((diag = do_diagonal (stack[last], k, i, j, nx1, ny1, 0)) == -2) continue;	/* Sorry, no extrema along diagonal N45E */
			if (diag != 0 && diag != dx) continue;						/* Sorry, extrema of opposite sign along diagonal N45E  */
			if ((diag = do_diagonal (stack[last], k, i, j, nx1, ny1, 1)) == -2) continue;	/* Sorry, no extrema along diagonal N135E */
			if (diag != 0 && diag != dx) continue;						/* Sorry, extrema of opposite sign along diagonal N135E  */

			/* OK, we have a min or max point; just use dx to check which kind */

			z[k] = (float)((dx > 0) ? +2 : -2);
		}
	}

	memcpy ((void *)stack[last], (void *)z, info->nm * sizeof (float));
	GMT_free ((void *)z);
}

/* Subroutines for grd_EXTREMA */

GMT_LONG do_derivative (float *z, GMT_LONG this_node, GMT_LONG this_pos, GMT_LONG last, GMT_LONG nx, GMT_LONG type)
{
	/* z is the data matrix */
	/* this_node is the index that gives us the center node */
	/* this_pos is either the x index (i) or the y index (j), depending on type */
	/* last_i is the last x index (nx-1) */
	/* last_j is the last y index (ny-1) */
	/* type = 0 means x-derivative, type = 1 means y-derivative */

	GMT_LONG off, next_node, prev_node;

	off = (type == 0) ? 1 : nx;

	if (this_pos == 0) {	/* Node at left column (x) or top row (y) - hence only one neighbor  - hence only one neighbor to the right (x) or below (y) */
		next_node = this_node + off;
		if (GMT_is_fnan (z[next_node])) return (-2);		/* One of the two points is a NaN */
		if (z[this_node] == z[next_node]) return (-2);		/* Flat line, no extrema possible */
		if (z[this_node] < z[next_node]) return (-1);		/* A local minimum */
		return (+1);						/* Else it must be a local maximum */
	}
	else if (this_pos == last) {	/* Node at right column (x) or bottom row (y) - hence only one neighbor to the left (x) or above (y) */
		prev_node = this_node - off;
		if (GMT_is_fnan (z[prev_node])) return (-2);		/* One of the two points is a NaN */
		if (z[this_node] == z[prev_node]) return (-2);		/* Flat line, no extrema possible */
		if (z[this_node] < z[prev_node]) return (-1);		/* A local minimum */
		return (+1);						/* Else it must be a local maximum */
	}
	else {	/* Node has neighbors on either side */
		prev_node = this_node - off;	/* Node to the left */
		next_node = this_node + off;	/* Node to the right */
		if (GMT_is_fnan (z[prev_node])) {			/* At least one of the two neighbor points is a NaN */
			if (GMT_is_fnan (z[next_node])) return (-2);	/* Both points are NaN */
			if (z[this_node] == z[next_node]) return (-2);	/* Flat line, no extrema possible */
			if (z[this_node] < z[next_node]) return (-1);	/* A local minimum */
			return (+1);					/* Else it must be a local maximum */
		}
		else if (GMT_is_fnan (z[next_node])) {			/* One of the two neighbor points is a NaN */
			if (z[this_node] == z[prev_node]) return (-2);	/* Flat line, no extrema possible */
			if (z[this_node] < z[prev_node]) return (-1);	/* A local minimum */
			return (+1);					/* Else it must be a local maximum */
		}
		else if (z[this_node] == z[prev_node] && z[this_node] == z[next_node])
			return (-2);					/* Flat line, no extrema possible */
		else if (z[this_node] < z[prev_node] && z[this_node] < z[next_node])
			return (-1);					/* A local minimum */
		else if (z[this_node] > z[prev_node] && z[this_node] > z[next_node])
			return (+1);					/* A local maximum */
		else
			return (0);					/* Nuthin' */
	}
}

GMT_LONG do_diagonal (float *z, GMT_LONG this_node, GMT_LONG this_i, GMT_LONG this_j, GMT_LONG last_i, GMT_LONG last_j, GMT_LONG type)
{
	/* Need to take d/dn in either the N45E direction (type = 0) or N135E (type = 1) direction */

	GMT_LONG next_node, prev_node, off;

	off = last_i + 1;
	if (type == 0) off = -off;

	/* diagonal check is only fatal if there is actual data to rule out a min or max */
	/* If there are NaNs we do not rule out the min/max */

	if (this_j == 0) {	/* Node is somewhere along top row so limited checking only */

		if (this_i == 0) {	/* Node is upper left corner - can only do limited checking on one half-diagonal */
			next_node = this_node + last_i + 2;
			if (GMT_is_fnan (z[next_node])) return (0);		/* The only neighbor points is a NaN */
			if (z[this_node] == z[next_node]) return (-2);		/* Flat line, no extrema possible */
			if (z[this_node] < z[next_node]) return (-1);		/* A local minimum */
			return (+1);						/* Else it must be a local maximum */
		}
		else if (this_i == last_i) {	/* Node is upper right corner - can only do limited checking on one half-diagonal */
			prev_node = this_node + last_i;
			if (GMT_is_fnan (z[prev_node])) return (0);		/* The only neighbor points is a NaN */
			if (z[this_node] == z[prev_node]) return (-2);		/* Flat line, no extrema possible */
			if (z[this_node] < z[prev_node]) return (-1);		/* A local minimum */
			return (+1);						/* Else it must be a local maximum */
		}
		else {	/* Node is along top row but not at the corners - must check one of the two half-diagonals below */
			next_node = this_node + ((type) ? last_i + 2 : last_i);
			if (GMT_is_fnan (z[next_node])) return (0);		/* One of the two neighbor points is a NaN */
			if (z[this_node] == z[next_node]) return (-2);		/* Flat line, no extrema possible */
			if (z[this_node] < z[next_node]) return (-1);		/* A local minimum */
			return (+1);						/* Else it must be a local maximum */
		}
	}
	else if (this_j == last_j) {	/* Node is somewhere along bottom row so limited checking only */

		if (this_i == 0) {	/* Node is lower left corner - can only do limited checking on one half-diagonal */
			next_node = this_node - last_i;
			if (GMT_is_fnan (z[next_node])) return (0);		/* The only neighbor points is a NaN */
			if (z[this_node] == z[next_node]) return (-2);		/* Flat line, no extrema possible */
			if (z[this_node] < z[next_node]) return (-1);		/* A local minimum */
			return (+1);						/* Else it must be a local maximum */
		}
		else if (this_i == last_i) {	/* Node is lower right corner - can only do limited checking on one half-diagonal */
			prev_node = this_node - last_i - 2;
			if (GMT_is_fnan (z[prev_node])) return (0);		/* The only neighbor points is a NaN */
			if (z[this_node] == z[prev_node]) return (-2);		/* Flat line, no extrema possible */
			if (z[this_node] < z[prev_node]) return (-1);		/* A local minimum */
			return (+1);						/* Else it must be a local maximum */
		}
		else {	/* Node is along bottom row but not at the corners - must check one of the two half-diagonals below */
			next_node = this_node - ((type) ? last_i + 2 : last_i);
			if (GMT_is_fnan (z[next_node])) return (0);		/* One of the two neighbor points is a NaN */
			if (z[this_node] == z[next_node]) return (-2);		/* Flat line, no extrema possible */
			if (z[this_node] < z[next_node]) return (-1);		/* A local minimum */
			return (+1);						/* Else it must be a local maximum */
		}
	}
	else {	/* Node is along one of the intermediate rows - checking is limited by x-position only */

		if (this_i == 0) {	/* Node is on left row - can only do limited checking on two half-diagonals */
			next_node = this_node + off + 1;
			if (GMT_is_fnan (z[next_node])) return (0);		/* One of the two neighbor points is a NaN */
			if (z[this_node] == z[next_node]) return (-2);		/* Flat line, no extrema possible */
			if (z[this_node] < z[next_node]) return (-1);		/* A local minimum */
			return (+1);						/* Else it must be a local maximum */
		}
		else if (this_i == last_i) {	/* Node is on right row - can only do limited checking on two half-diagonals */
			next_node = this_node + off - 1;
			if (GMT_is_fnan (z[next_node])) return (0);		/* One of the two neighbor points is a NaN */
			if (z[this_node] == z[next_node]) return (-2);		/* Flat line, no extrema possible */
			if (z[this_node] < z[next_node]) return (-1);		/* A local minimum */
			return (+1);						/* Else it must be a local maximum */
		}
		else {	/* Node has all required neighbors */
			next_node = this_node + off + 1;
			prev_node = this_node - off - 1;
			if (GMT_is_fnan (z[prev_node])) {			/* At least one of the two neighbor points is a NaN */
				if (GMT_is_fnan (z[next_node])) return (0);	/* Both points are NaN */
				if (z[this_node] == z[next_node]) return (-2);	/* Flat line, no extrema possible */
				if (z[this_node] < z[next_node]) return (-1);	/* A local minimum */
				return (+1);					/* Else it must be a local maximum */
			}
			else if (GMT_is_fnan (z[next_node])) {			/* One of the two neighbor points is a NaN */
				if (z[this_node] == z[prev_node]) return (-2);	/* Flat line, no extrema possible */
				if (z[this_node] < z[prev_node]) return (-1);	/* A local minimum */
				return (+1);					/* Else it must be a local maximum */
			}
			else if (z[this_node] == z[prev_node] && z[this_node] == z[next_node])
				return (-2);					/* Flat line, no extrema possible */
			else if (z[this_node] < z[prev_node] && z[this_node] < z[next_node])
				return (-1);					/* A local minimum */
			else if (z[this_node] > z[prev_node] && z[this_node] > z[next_node])
				return (+1);					/* A local maximum */
			else
				return (0);					/* Nuthin' */
		}
	}
}

void grd_FCRIT (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: FCRIT 3 1 Critical value for F-distribution, with alpha = A, n1 = B, and n2 = C.  */
{
	GMT_LONG i, nu1, nu2, prev1, prev2;
	double alpha;

	prev1 = last - 1;
	prev2 = last - 2;
	if (gmtdefs.verbose && constant[prev2] && factor[prev2] == 0.0) fprintf (stderr, "%s: Warning, operand one == 0 for FCRIT!\n", GMT_program);
	if (gmtdefs.verbose && constant[prev1] && factor[prev1] == 0.0) fprintf (stderr, "%s: Warning, operand two == 0 for FCRIT!\n", GMT_program);
	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand three == 0 for FCRIT!\n", GMT_program);
	for (i = 0; i < info->nm; i++) {
		alpha = (constant[prev2]) ? factor[prev2] : stack[prev2][i];
		nu1 = irint ((double)((constant[prev1]) ? factor[prev1] : stack[prev1][i]));
		nu2= irint ((double)((constant[last]) ? factor[last] : stack[last][i]));
		stack[prev2][i] = (float)GMT_Fcrit (alpha, (double)nu1, (double)nu2);
	}
}

void grd_FDIST (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: FDIST 3 1 F-distribution Q(F,n1,n2), with F = A, n1 = B, and n2 = C.  */
{
	GMT_LONG i, nu1, nu2, prev1, prev2;
	double F, chisq1, chisq2 = 1.0, prob;

	prev1 = last - 1;
	prev2 = last - 2;
	if (gmtdefs.verbose && constant[prev1] && factor[prev1] == 0.0) fprintf (stderr, "%s: Warning, operand two == 0 for FDIST!\n", GMT_program);
	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand three == 0 for FDIST!\n", GMT_program);
	for (i = 0; i < info->nm; i++) {
		F = (constant[prev2]) ? factor[prev2] : stack[prev2][i];
		nu1 = (GMT_LONG)(irint ((double)((constant[prev1]) ? factor[prev1] : stack[prev1][i])));
		nu2 = (GMT_LONG)(irint ((double)((constant[last]) ? factor[last] : stack[last][i])));
		/* Since GMT_f_q needs chisq1 and chisq2, we set chisq2 = 1 and solve for chisq1 */
		chisq1 = F * nu1 / nu2;
		(void) GMT_f_q (chisq1, nu1, chisq2, nu2, &prob);
		stack[prev2][i] = (float)prob;
	}
}

void grd_FLIPLR (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: FLIPLR 1 1 Reverse order of values in each row.  */
{
	GMT_LONG k, k0, nx1, i, j, nx_half;

	/* Reverse order of all rows */

	if (constant[last]) {
		if (gmtdefs.verbose) fprintf (stderr, "%s: Warning, operand to FLIPLR is constant!\n", GMT_program);
		return;
	}

	nx_half = info->header.nx / 2;
	nx1 = (size_t)info->header.nx - 1;
	for (j = k0 = 0; j < info->header.ny; j++, k0 += (size_t)info->header.nx) {	/* Do this to all rows */
		for (i = 0, k = nx1; i < nx_half; i++, k--) {
			f_swap (stack[last][k0+i], stack[last][k0+k]);
		}
	}
}

void grd_FLIPUD (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: FLIPUD 1 1 Reverse order of values in each column.  */
{
	GMT_LONG  k, ny1, nx, i, j, ny_half;

	/* Reverse order of all columns */

	if (constant[last]) {
		if (gmtdefs.verbose) fprintf (stderr, "%s: Warning, operand to FLIPLR is constant!\n", GMT_program);
		return;
	}

	ny_half = info->header.ny / 2;
	ny1 = (size_t)info->header.ny - 1;
	nx = (size_t)info->header.nx;
	for (i = 0; i < info->header.nx; i++) {
		for (j = 0, k = ny1; j < ny_half; j++, k--) {	/* Do this to all rows */
			f_swap (stack[last][j*nx+i], stack[last][k*nx+i]);
		}
	}
}

void grd_FLOOR (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: FLOOR 1 1 floor (A) (greatest integer <= A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = floor (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : floor ((double)stack[last][i]));
}

void grd_FMOD (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: FMOD 2 1 A % B (remainder after truncated division).  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;
	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, using FMOD 0!\n", GMT_program);
	for (i = 0; i < info->nm; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][i];
		b = (constant[last]) ? factor[last] : stack[last][i];
		stack[prev][i] = (float)fmod (a, b);
	}
}

void grd_GE (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: GE 2 1 1 if A >= B, else 0.  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;
	for (i = 0; i < info->nm; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][i];
		b = (constant[last]) ? factor[last] : stack[last][i];
		stack[prev][i] = (GMT_is_dnan (a) || GMT_is_dnan (b)) ? GMT_f_NaN : (float)(a >= b);
	}
}

void grd_GT (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: GT 2 1 1 if A > B, else 0.  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;
	for (i = 0; i < info->nm; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][i];
		b = (constant[last]) ? factor[last] : stack[last][i];
		stack[prev][i] = (GMT_is_dnan (a) || GMT_is_dnan (b)) ? GMT_f_NaN : (float)(a > b);
	}
}

void grd_HYPOT (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: HYPOT 2 1 hypot (A, B) = sqrt (A*A + B*B).  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;
	if (gmtdefs.verbose && constant[prev] && factor[prev] == 0.0) fprintf (stderr, "%s: Warning, operand one == 0!\n", GMT_program);
	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand two == 0!\n", GMT_program);
	for (i = 0; i < info->nm; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][i];
		b = (constant[last]) ? factor[last] : stack[last][i];
		stack[prev][i] = (float)hypot (a, b);
	}
}

void grd_I0 (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: I0 1 1 Modified Bessel function of A (1st kind, order 0).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = GMT_i0 (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : GMT_i0 ((double)stack[last][i]));
}

void grd_I1 (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: I1 1 1 Modified Bessel function of A (1st kind, order 1).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = GMT_i1 (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : GMT_i1 ((double)stack[last][i]));
}

void grd_IN (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: IN 2 1 Modified Bessel function of A (1st kind, order B).  */
{
	GMT_LONG i, prev, order = 0;
	GMT_LONG simple = FALSE;
	double b = 0.0;

	prev = last - 1;
	if (constant[last]) {
		if (gmtdefs.verbose && factor[last] < 0.0) fprintf (stderr, "%s: Warning, order < 0 for IN!\n", GMT_program);
		if (gmtdefs.verbose && fabs (rint(factor[last]) - factor[last]) > GMT_SMALL) fprintf (stderr, "%s: Warning, order not an integer for IN!\n", GMT_program);
		order = irint (fabs (factor[last]));
		if (constant[prev]) {
			b = GMT_in (order, fabs (factor[prev]));
			simple = TRUE;
		}
	}
	for (i = 0; i < info->nm; i++) {
		if (simple)
			stack[prev][i] = (float)b;
		else {
			if (!constant[last]) order = irint (fabs ((double)stack[last][i]));
			stack[last][i] = (float)GMT_in (order, fabs ((double)stack[prev][i]));
		}
	}
}

void grd_INRANGE (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: INRANGE 3 1 1 if B <= A <= C, else 0.  */
{
	GMT_LONG i, prev1, prev2;
	GMT_LONG inrange;
	float a = 0.0, b = 0.0, c = 0.0;

	/* last is C */
	prev1 = last - 1;	/* This is B */
	prev2 = last - 2;	/* This is A */

	/* Set to 1 where B <= A <= C, 0 elsewhere, except where
	 * A, B, or C = NaN, in which case we set answer to NaN */

	if (constant[prev2]) a = (float)factor[prev2];
	if (constant[prev1]) b = (float)factor[prev1];
	if (constant[last])  c = (float)factor[last];

	for (i = 0; i < info->nm; i++) {
		if (!constant[prev2]) a = stack[prev2][i];
		if (!constant[prev1]) b = stack[prev1][i];
		if (!constant[last])  c = stack[last][i];

		if (GMT_is_fnan (a) || GMT_is_fnan (b) || GMT_is_fnan (c)) {
			stack[prev2][i] = GMT_f_NaN;
			continue;
		}

		inrange = (b <= a && a <= c);
		stack[prev2][i] = (float)inrange;
	}
}

void grd_INSIDE (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: INSIDE 1 1 1 when inside or on polygon(s) in A, else 0.  */
{
	void grd_INSIDE_GEO (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last);
	void grd_INSIDE_XY (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last);

	if (GMT_io.in_col_type[0] & GMT_IS_GEO)	/* Geographic data */
		grd_INSIDE_GEO (info, stack, constant, factor, last);
	else					/* Cartesian data */
		grd_INSIDE_XY (info, stack, constant, factor, last);
}

void grd_INSIDE_GEO (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
{	/* Suitable for geographic (lon, lat) data and polygons */
	GMT_LONG i, j, P, inside, k;
	GMT_LONG greenwich;
	struct GMT_TABLE *polygon;

	greenwich = (info->header.x_min < 0.0 && info->header.x_max > 0.0);

	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = 2;	/* Default is 2 input cols for binary data */
	GMT_io.skip_duplicates = TRUE;	/* The inonout algorithm assumes no duplicates */
	GMT_import_table ((void *)(info->ASCII_file), GMT_IS_FILE, &polygon, 0.0, greenwich, TRUE, TRUE);	/* Closed polygons */

	for (i = j = k = 0; k < info->nm; k++) {	/* Visit each node */
		for (P = inside = 0; !inside && P < polygon->n_segments; P++) {
			inside = GMT_inonout_sphpol ((double)info->grd_x[i], (double)info->grd_y[j], polygon->segment[P]);
		}
		stack[last][k] = (float)((inside) ? 1.0 : 0.0);
		i++;
		if (i == info->header.nx) {	/* Go to next row */
			i = 0, j++;
		}
	}
	GMT_io.skip_duplicates = FALSE;	/* Reset to FALSE */

	/* Free memory used for pol */

	GMT_free_table (polygon);
}

void grd_INSIDE_XY (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
{	/* Version for Cartesian data */
	GMT_LONG i, j, P, inside, k;
	struct GMT_TABLE *polygon;

	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = 2;	/* Default is 2 input cols for binary data */
	GMT_import_table ((void *)(info->ASCII_file), GMT_IS_FILE, &polygon, 0.0, FALSE, TRUE, TRUE);	/* Closed polygons */

	for (i = j = k = 0; k < info->nm; k++) {	/* Visit each node */
		for (P = inside = 0; !inside && P < polygon->n_segments; P++) {
			if (info->grd_y[j] < polygon->segment[P]->min[GMT_Y] || info->grd_y[j] > polygon->segment[P]->max[GMT_Y]) continue;	/* Outside y range */
			if (info->grd_x[j] < polygon->segment[P]->min[GMT_X] || info->grd_x[j] > polygon->segment[P]->max[GMT_X]) continue;	/* Outside x range */
			inside = GMT_non_zero_winding ((double)info->grd_x[i], (double)info->grd_y[j], polygon->segment[P]->coord[GMT_X], polygon->segment[P]->coord[GMT_Y], polygon->segment[P]->n_rows);	/* Must test */
		}
		stack[last][k] = (float)((inside) ? 1.0 : 0.0);
		i++;
		if (i == info->header.nx) {	/* Go to next row */
			i = 0, j++;
		}
	}

	/* Free memory used for pol */

	GMT_free_table (polygon);
}

void grd_INV (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: INV 1 1 1 / A.  */
{
	GMT_LONG i;
	double a;

	if (constant[last] && factor[last] == 0.0) {
		fprintf (stderr, "%s: Error, Cannot take inverse of zero!\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (constant[last]) factor[last] = 1.0 / factor[last];
	for (i = 0; i < info->nm; i++) {
		a = (constant[last]) ? factor[last] : 1.0 / stack[last][i];
		stack[last][i] = (float)a;
	}
}

void grd_ISNAN (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: ISNAN 1 1 1 if A == NaN, else 0.  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = GMT_is_dnan (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : GMT_is_fnan (stack[last][i]));
}

void grd_J0 (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: J0 1 1 Bessel function of A (1st kind, order 0).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = j0 (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : j0 ((double)stack[last][i]));
}

void grd_J1 (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: J1 1 1 Bessel function of A (1st kind, order 1).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = j1 (fabs (factor[last]));
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : j1 (fabs ((double)stack[last][i])));
}

void grd_JN (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: JN 2 1 Bessel function of A (1st kind, order B).  */
{
	GMT_LONG i, prev;
	int order = 0;
	GMT_LONG simple = FALSE;
	double b = 0.0;

	prev = last - 1;
	if (constant[last]) {
		if (gmtdefs.verbose && factor[last] < 0.0) fprintf (stderr, "%s: Warning, order < 0 for JN!\n", GMT_program);
		if (gmtdefs.verbose && fabs (rint(factor[last]) - factor[last]) > GMT_SMALL) fprintf (stderr, "%s: Warning, order not an integer for JN!\n", GMT_program);
		order = irint (fabs (factor[last]));
		if (constant[prev]) {
			b = jn (order, fabs (factor[prev]));
			simple = TRUE;
		}
	}
	for (i = 0; i < info->nm; i++) {
		if (simple)
			stack[prev][i] = (float)b;
		else {
			if (!constant[last]) order = irint (fabs ((double)stack[last][i]));
			stack[last][i] = (float)jn (order, fabs ((double)stack[prev][i]));
		}
	}
}

void grd_K0 (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: K0 1 1 Modified Kelvin function of A (2nd kind, order 0).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = GMT_k0 (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : GMT_k0 ((double)stack[last][i]));
}

void grd_K1 (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: K1 1 1 Modified Bessel function of A (2nd kind, order 1).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = GMT_k1 (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : GMT_k1 ((double)stack[last][i]));
}

void grd_KEI (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: KEI 1 1 kei (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = GMT_kei (fabs (factor[last]));
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : GMT_kei (fabs ((double)stack[last][i])));
}

void grd_KER (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: KER 1 1 ker (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = GMT_ker (fabs (factor[last]));
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : GMT_ker (fabs ((double)stack[last][i])));
}

void grd_KM2DEG (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: KM2DEG 1 1 Converts Kilometers to Spherical Degrees.  */
{
	GMT_LONG i;
	double a = 0.0, f = 1.0 / project_info.DIST_KM_PR_DEG;

	if (GMT_io.in_col_type[0] & GMT_IS_GEO) {	/* Geographic data */
		if (!GMT_IS_SPHERICAL) fprintf (stderr, "%s: Warning, KM2DEG is only exact when ELLIPSOID == sphere\n", GMT_program);
	}
	else
		fprintf (stderr, "%s: Warning, KM2DEG used with Cartesian data\n", GMT_program);
	if (constant[last]) a = factor[last] * f;
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : (stack[last][i] * project_info.DIST_KM_PR_DEG));
}

void grd_KN (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: KN 2 1 Modified Bessel function of A (2nd kind, order B).  */
{
	GMT_LONG i, prev, order = 0;
	GMT_LONG simple = FALSE;
	double b = 0.0;

	prev = last - 1;
	if (constant[last]) {
		if (gmtdefs.verbose && factor[last] < 0.0) fprintf (stderr, "%s: Warning, order < 0 for KN!\n", GMT_program);
		if (gmtdefs.verbose && fabs (rint(factor[last]) - factor[last]) > GMT_SMALL) fprintf (stderr, "%s: Warning, order not an integer for KN!\n", GMT_program);
		order = irint (fabs (factor[last]));
		if (constant[prev]) {
			b = GMT_kn (order, fabs (factor[prev]));
			simple = TRUE;
		}
	}
	for (i = 0; i < info->nm; i++) {
		if (simple)
			stack[prev][i] = (float)b;
		else {
			if (!constant[last]) order = irint (fabs ((double)stack[last][i]));
			stack[last][i] = (float)GMT_kn (order, fabs ((double)stack[prev][i]));
		}
	}
}

void grd_KURT (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: KURT 1 1 Kurtosis of A.  */
{
	GMT_LONG i, n = 0;
	double mean = 0.0, sum2 = 0.0, kurt = 0.0, delta;
	float f_kurt;

	if (constant[last]) {	/* Trivial case */
		for (i = 0; i < info->nm; i++) stack[last][i] = GMT_f_NaN;
		return;
	}

	/* Use Welford (1962) algorithm to compute mean and corrected sum of squares */
	for (i = 0; i < info->nm; i++) {
		if (GMT_is_fnan (stack[last][i])) continue;
		n++;
		delta = (double)stack[last][i] - mean;
		mean += delta / n;
		sum2 += delta * ((double)stack[last][i] - mean);
	}
	if (n > 1) {
		for (i = 0; i < info->nm; i++) {
			if (GMT_is_fnan (stack[last][i])) continue;
			delta = (double)stack[last][i] - mean;
			kurt += pow (delta, 4.0);
		}
		sum2 /= (n - 1);
		kurt = kurt / (n * sum2 * sum2) - 3.0;
		f_kurt = (float)kurt;
	}
	else
		f_kurt = GMT_f_NaN;
	for (i = 0; i < info->nm; i++) stack[last][i] = f_kurt;
}

void grd_LDIST (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: LDIST 1 1 Compute distance from lines in multi-segment ASCII file A.  */
{
	GMT_LONG i, j, k;
	double d;
	GMT_LONG greenwich;
	PFL near_a_line;
	struct GMT_TABLE *line;

	if (GMT_io.in_col_type[0] & GMT_IS_GEO) {	/* Geographic data */
		GMT_distance_func = (PFD) ((GMT_IS_SPHERICAL) ? GMT_great_circle_dist_km : GMT_geodesic_dist_km);
		near_a_line = (PFL) GMT_near_a_line_spherical;
		greenwich = (info->header.x_min < 0.0 && info->header.x_max > 0.0);
	}
	else {
		GMT_distance_func = (PFD) GMT_cartesian_dist;
		near_a_line = (PFL) GMT_near_a_line_cartesian;
		greenwich = FALSE;
	}

	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = 2;	/* Default is 2 input cols for binary data */
	GMT_import_table ((void *)(info->ASCII_file), GMT_IS_FILE, &line, 0.0, greenwich, FALSE, TRUE);

	for (i = j = k = 0; k < info->nm; k++) {	/* Visit each node */
		(void) near_a_line ((double)info->grd_x[i], (double)info->grd_y[j], line, TRUE, &d, NULL, NULL);
		stack[last][k] = (float)d;
		i++;
		if (i == info->header.nx) {
			if (gmtdefs.verbose == 2) fprintf (stderr, "Line %ld\n", j);
			i = 0, j++;
		}
	}

	/* Free memory used for line */

	GMT_free_table (line);
}

void grd_LE (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: LE 2 1 1 if A <= B, else 0.  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;
	for (i = 0; i < info->nm; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][i];
		b = (constant[last]) ? factor[last] : stack[last][i];
		stack[prev][i] = (GMT_is_dnan (a) || GMT_is_dnan (b)) ? GMT_f_NaN : (float)(a <= b);
	}
}

void grd_LOG (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: LOG 1 1 log (A) (natural log).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, argument to log = 0\n", GMT_program);

	if (constant[last]) a = d_log (fabs (factor[last]));
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : d_log (fabs ((double)stack[last][i])));
}

void grd_LOG10 (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: LOG10 1 1 log10 (A) (base 10).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, argument to log10 = 0\n", GMT_program);

	if (constant[last]) a = d_log10 (fabs (factor[last]));
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : d_log10 (fabs ((double)stack[last][i])));
}

void grd_LOG1P (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: LOG1P 1 1 log (1+A) (accurate for small A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && factor[last] < 0.0) fprintf (stderr, "%s: Warning, argument to log1p < 0\n", GMT_program);

	if (constant[last]) a = d_log1p (fabs (factor[last]));
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : d_log1p (fabs ((double)stack[last][i])));
}

void grd_LOG2 (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: LOG2 1 1 log2 (A) (base 2).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, argument to log2 = 0\n", GMT_program);

	if (constant[last]) a = d_log (fabs (factor[last])) * M_LN2_INV;
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)(((constant[last]) ? a : d_log (fabs ((double)stack[last][i]))) * M_LN2_INV);
}

void grd_LMSSCL (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: LMSSCL 1 1 LMS scale estimate (LMS STD) of A.  */
{
	GMT_LONG GMT_mode_selection = 0, GMT_n_multiples = 0, i;
	double mode, lmsscl;
	float lmsscl_f;

	if (constant[last]) {	/* Trivial case */
		memset ((void *)stack[last], 0, info->nm * sizeof (float));
		return;
	}

	/* Sort will put any NaNs to the end - we then count to find the real data */

	qsort ((void *)stack[last], info->nm, sizeof (float), GMT_comp_float_asc);
	for (i = info->nm; GMT_is_fnan (stack[last][i-1]) && i > 1; i--);
	if (i) {
		GMT_mode_f (stack[last], (GMT_LONG)i, (GMT_LONG)(i/2), 0, GMT_mode_selection, &GMT_n_multiples, &mode);
		GMT_getmad_f (stack[last], (GMT_LONG)i, mode, &lmsscl);
		lmsscl_f = (float)lmsscl;
	}
	else
		lmsscl_f = GMT_f_NaN;

	for (i = 0; i < info->nm; i++) stack[last][i] = lmsscl_f;
	if (GMT_n_multiples > 0) fprintf (stderr, "%s: WARNING: %ld Multiple modes found\n", GMT_program, GMT_n_multiples);
}

void grd_LOWER (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: LOWER 1 1 The lowest (minimum) value of A.  */
{
	GMT_LONG i;
	float low;

	if (constant[last]) {	/* Trivial case */
		for (i = 0; i < info->nm; i++) stack[last][i] = (float)factor[last];
		return;
	}

	for (i = 0, low = FLT_MAX; i < info->nm; i++) {
		if (GMT_is_fnan (stack[last][i])) continue;
		if (stack[last][i] < low) low = stack[last][i];
	}
	for (i = 0; i < info->nm; i++) if (!GMT_is_fnan (stack[last][i])) stack[last][i] = low;
}

void grd_LRAND (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: LRAND 2 1 Laplace random noise with mean A and std. deviation B.  */
{
	GMT_LONG i, prev;
	double a = 0.0, b = 0.0;

	prev = last - 1;
	if (constant[prev]) a = factor[prev];
	if (constant[last]) b = factor[last];
	for (i = 0; i < info->nm; i++) {
		if (!constant[prev]) a = (double)stack[prev][i];
		if (!constant[last]) b = (double)stack[last][i];
		stack[prev][i] = (float)(a + b * GMT_lrand ());
	}
}

void grd_LT (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: LT 2 1 1 if A < B, else 0.  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;
	for (i = 0; i < info->nm; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][i];
		b = (constant[last]) ? factor[last] : stack[last][i];
		stack[prev][i] = (GMT_is_dnan (a) || GMT_is_dnan (b)) ? GMT_f_NaN : (float)(a < b);
	}
}

void grd_MAD (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: MAD 1 1 Median Absolute Deviation (L1 STD) of A.  */
{
	GMT_LONG i;
	double mad, med;
	float mad_f;

	if (constant[last]) {	/* Trivial case */
		memset ((void *)stack[last], 0, info->nm * sizeof (float));
		return;
	}

	/* Sort will put any NaNs to the end - we then count to find the real data */

	qsort ((void *)stack[last], info->nm, sizeof (float), GMT_comp_float_asc);
	for (i = info->nm; GMT_is_fnan (stack[last][i-1]) && i > 1; i--);
	if (i) {
		med = (i%2) ? stack[last][i/2] : (float)(0.5 * (stack[last][(i-1)/2] + stack[last][i/2]));
		GMT_getmad_f (stack[last], (GMT_LONG)i, med, &mad);
		mad_f = (float)mad;
	}
	else
		mad_f = GMT_f_NaN;

	for (i = 0; i < info->nm; i++) stack[last][i] = mad_f;
}

void grd_MAX (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: MAX 2 1 Maximum of A and B.  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;
	for (i = 0; i < info->nm; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][i];
		b = (constant[last]) ? factor[last] : stack[last][i];
		stack[prev][i] = (GMT_is_dnan (a) || GMT_is_dnan (b)) ? GMT_f_NaN : (float)MAX (a, b);
	}
}

void grd_MEAN (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: MEAN 1 1 Mean value of A.  */
{
	GMT_LONG i, n_a = 0;
	double sum_a = 0.0;

	if (constant[last]) {	/* Trivial case */
		for (i = 0; i < info->nm; i++) stack[last][i] = (float)factor[last];
		return;
	}

	for (i = 0; i < info->nm; i++) {
		if (GMT_is_fnan (stack[last][i])) continue;
		sum_a += stack[last][i];
		n_a++;
	}
	sum_a = (n_a) ? sum_a / (double)n_a : 0.0;
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)sum_a;
}

void grd_MED (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: MED 1 1 Median value of A.  */
{
	GMT_LONG i;
	float med;

	if (constant[last]) {	/* Trivial case */
		for (i = 0; i < info->nm; i++) stack[last][i] = (float)factor[last];
		return;
	}

	qsort ((void *)stack[last], info->nm, sizeof (float), GMT_comp_float_asc);
	for (i = info->nm; GMT_is_fnan (stack[last][i-1]) && i > 1; i--);
	if (i)
		med = (i%2) ? stack[last][i/2] : (float)(0.5 * (stack[last][(i-1)/2] + stack[last][i/2]));
	else
		med = GMT_f_NaN;

	for (i = 0; i < info->nm; i++) stack[last][i] = med;
}

void grd_MIN (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: MIN 2 1 Minimum of A and B.  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;
	for (i = 0; i < info->nm; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][i];
		b = (constant[last]) ? factor[last] : stack[last][i];
		stack[prev][i] = (GMT_is_dnan (a) || GMT_is_dnan (b)) ? GMT_f_NaN : (float)MIN (a, b);
	}
}

void grd_MOD (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: MOD 2 1 A mod B (remainder after floored division).  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;
	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, using MOD 0!\n", GMT_program);
	for (i = 0; i < info->nm; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][i];
		b = (constant[last]) ? factor[last] : stack[last][i];
		stack[prev][i] = (float)MOD (a, b);
	}
}

void grd_MODE (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: MODE 1 1 Mode value (Least Median of Squares) of A.  */
{
	GMT_LONG GMT_mode_selection = 0, GMT_n_multiples = 0, i;
	double mode = 0.0;

	if (constant[last]) {	/* Trivial case */
		for (i = 0; i < info->nm; i++) stack[last][i] = (float)factor[last];
		return;
	}

	qsort ((void *)stack[last], info->nm, sizeof (float), GMT_comp_float_asc);
	for (i = info->nm; GMT_is_fnan (stack[last][i-1]) && i > 1; i--);
	if (i)
		GMT_mode_f (stack[last], (GMT_LONG)i, (GMT_LONG)(i/2), 0, GMT_mode_selection, &GMT_n_multiples, &mode);
	else
		mode = GMT_f_NaN;

	for (i = 0; i < info->nm; i++) stack[last][i] = (float)mode;
	if (GMT_n_multiples > 0) fprintf (stderr, "%s: WARNING: %ld Multiple modes found\n", GMT_program, GMT_n_multiples);
}

void grd_MUL (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: MUL 2 1 A * B.  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;
	if (gmtdefs.verbose && constant[prev] && factor[prev] == 0.0) fprintf (stderr, "%s: Warning, operand one == 0!\n", GMT_program);
	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand two == 0!\n", GMT_program);
	for (i = 0; i < info->nm; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][i];
		b = (constant[last]) ? factor[last] : stack[last][i];
		stack[prev][i] = (float)(a * b);
	}
}

void grd_NAN (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: NAN 2 1 NaN if A == B, else A.  */
{
	GMT_LONG i, prev;
	double a = 0.0, b = 0.0;

	prev = last - 1;
	if (constant[prev]) a = factor[prev];
	if (constant[last]) b = factor[last];
	for (i = 0; i < info->nm; i++) {
		if (!constant[prev]) a = stack[prev][i];
		if (!constant[last]) b = stack[last][i];
		stack[prev][i] = ((a == b) ? GMT_f_NaN : (float)a);
	}
}

void grd_NEG (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: NEG 1 1 -A.  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand == 0!\n", GMT_program);
	if (constant[last]) a = -factor[last];
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : -stack[last][i]);
}

void grd_NEQ (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: NEQ 2 1 1 if A != B, else 0.  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;
	for (i = 0; i < info->nm; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][i];
		b = (constant[last]) ? factor[last] : stack[last][i];
		stack[prev][i] = (float)(a != b);
	}
}

void grd_NOT (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: NOT 1 1 NaN if A == NaN, 1 if A == 0, else 0.  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand == 0!\n", GMT_program);
	if (constant[last]) a = (fabs (factor[last]) > GMT_CONV_LIMIT) ? 0.0 : 1.0;
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : ((fabs (stack[last][i]) > GMT_CONV_LIMIT) ? 0.0 : 1.0));
}

void grd_NRAND (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: NRAND 2 1 Normal, random values with mean A and std. deviation B.  */
{
	GMT_LONG i, prev;
	double a = 0.0, b = 0.0;

	prev = last - 1;
	if (constant[prev]) a = factor[prev];
	if (constant[last]) b = factor[last];
	for (i = 0; i < info->nm; i++) {
		if (!constant[prev]) a = (double)stack[prev][i];
		if (!constant[last]) b = (double)stack[last][i];
		stack[prev][i] = (float)(a + b * GMT_nrand ());
	}
}

void grd_OR (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: OR 2 1 NaN if A or B == NaN, else A.  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;
	for (i = 0; i < info->nm; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][i];
		b = (constant[last]) ? factor[last] : stack[last][i];
		stack[prev][i] = (float)((GMT_is_dnan (a) || GMT_is_dnan (b)) ? GMT_f_NaN : a);
	}
}

void grd_PDIST (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: PDIST 1 1 Compute distance from points in ASCII file A.  */
{
	GMT_LONG i, j, dummy[2], k;
	struct GMT_TABLE *T;

	if (GMT_io.in_col_type[0] & GMT_IS_GEO) {	/* Geographic data */
		GMT_distance_func = (PFD) ((GMT_IS_SPHERICAL) ? GMT_great_circle_dist_km : GMT_geodesic_dist_km);
	}
	else {
		GMT_distance_func = (PFD) GMT_cartesian_dist;
	}

	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = 2;	/* Default is 2 input cols for binary data */
	GMT_import_table ((void *)info->ASCII_file, GMT_IS_FILE, &T, 0.0, FALSE, FALSE, TRUE);

	for (i = j = k = 0; k < info->nm; k++) {	/* Visit each node */
		stack[last][k] = (float)GMT_dist_to_point ((double)info->grd_x[i], (double)info->grd_y[j], T, dummy);
		i++;
		if (i == info->header.nx) i = 0, j++;
	}

	/* Free memory used for points */

	GMT_free_table (T);
}

void grd_POP (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: POP 1 0 Delete top element from the stack.  */
{

	/* Dummy routine that does nothing but consume the top element of stack */
}

void grd_PLM (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: PLM 3 1 Associated Legendre polynomial P(A) degree B order C.  */
{
	GMT_LONG prev, first, L, M, i;
	double a = 0.0;
				/* last holds the order M */
	prev  = last - 1;	/* prev holds the degree L */
	first = prev - 1;	/* first holds the argument x = cos(colat) */

	if (!(constant[prev] && constant[last])) {
		fprintf (stderr, "%s: L and M must be constants in PLM!\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	L = irint (factor[prev]);
	M = irint (factor[last]);

	if (constant[first]) a = GMT_plm (L, M, factor[first]);
	for (i = 0; i < info->nm; i++) stack[first][i] = (float)((constant[first]) ? a : GMT_plm (L, M, stack[first][i]));
}


void grd_PLMg (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: PLMg 3 1 Normalized associated Legendre polynomial P(A) degree B order C (geophysical convention).  */
{
	GMT_LONG prev, first, L, M, i;
	double a = 0.0;
				/* last holds the order M */
	prev  = last - 1;	/* prev holds the degree L */
	first = prev - 1;	/* first holds the argument x = cos(colat) */

	if (!(constant[prev] && constant[last])) {
		fprintf (stderr, "%s: L and M must be constants in PLMg!\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	L = irint (factor[prev]);
	M = irint (factor[last]);

	if (constant[first]) a = GMT_plm_bar (L, M, factor[first], FALSE);
	for (i = 0; i < info->nm; i++) stack[first][i] = (float)((constant[first]) ? a : GMT_plm_bar (L, M, stack[first][i], FALSE));
}

void grd_POW (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: POW 2 1 A ^ B.  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;

	if (gmtdefs.verbose && constant[prev] && factor[prev] == 0.0) fprintf (stderr, "%s: Warning, operand one == 0!\n", GMT_program);
	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand two == 0!\n", GMT_program);
	for (i = 0; i < info->nm; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][i];
		b = (constant[last]) ? factor[last] : stack[last][i];
		stack[prev][i] = (float)pow (a, b);
	}
}

void grd_PQUANT (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: PQUANT 2 1 The B'th Quantile (0-100%) of A.  */
{
	GMT_LONG i, prev;
	float p;

	prev  = last - 1;	/* last holds the selected quantile (0-100), prev the data % */
	if (!constant[last]) {
		fprintf (stderr, "%s: Error: PQUANT must be given a constant quantile\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (factor[last] < 0.0 || factor[last] > 100.0) {
		fprintf (stderr, "%s: Error: PQUANT must be given a constant quantile between 0-100%%\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (constant[prev]) {	/* Trivial case */
		fprintf (stderr, "%s: Warning: PQUANT of a constant is set to NaN\n", GMT_program);
		p = GMT_f_NaN;
	}
	else {
		qsort ((void *)stack[prev], (size_t)info->nm, sizeof (float), GMT_comp_float_asc);
		p = (float) GMT_quantile_f (stack[prev], factor[last], (GMT_LONG)info->nm);
	}

	for (i = 0; i < info->nm; i++) stack[prev][i] = p;
}

void grd_PSI (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: PSI 1 1 Psi (or Digamma) of A.  */
{
	GMT_LONG i;
	double a = 0.0, x[2];

	x[1] = 0.0;	/* No imaginary part */
	if (constant[last]) {
		x[0] = factor[last];
		a = GMT_psi (x, (double *)NULL);
	}

	for (i = 0; i < info->nm; i++) {
		if (!constant[last]) {
			x[0] = (double)stack[last][i];
			a = GMT_psi (x, (double *)NULL);
		}
		stack[last][i] = (float)a;
	}
}

void grd_PV (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: PV 3 1 Legendre function Pv(A) of degree v = real(B) + imag(C).  */
{
	grd_PVQV (info, stack, constant, factor, last, 0);
}

void grd_QV (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: QV 3 1 Legendre function Qv(A) of degree v = real(B) + imag(C).  */
{
	grd_PVQV (info, stack, constant, factor, last, 1);
}

void grd_PVQV (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG kind)
{
	GMT_LONG prev, first, n, i;
	double a = 0.0, x = 0.0, nu[2], pq[4];
	static char *name[2] = {"PV", "QV"};
	GMT_LONG calc;
				/* last holds the imaginary order vi */
	prev  = last - 1;	/* prev holds the real order vr */
	first = prev - 1;	/* first holds the argument x = cos(colat) */

	calc = !(constant[prev] && constant[last] && constant[first]);	/* Only constant it all args are constant */
	if (!calc) {	/* All constants */
		nu[0] = factor[prev];
		nu[1] = factor[last];
		if (gmtdefs.verbose && (factor[first] < -1.0 || factor[first] > 1.0)) fprintf (stderr, "%s: Warning, argument to %s outside domain!\n", GMT_program, name[kind]);
		GMT_PvQv (factor[first], nu, pq, &n);
		a = pq[2*kind];
	}
	if (constant[prev]) nu[0] = factor[prev];
	if (constant[last]) nu[1] = factor[last];
	if (constant[first])    x = factor[first];
	kind *= 2;
	for (i = 0; i < info->nm; i++) {
		if (calc){
			if (!constant[prev]) nu[0] = (double)stack[prev][i];
			if (!constant[last]) nu[1] = (double)stack[last][i];
			if (!constant[first])    x = (double)stack[first][i];
			GMT_PvQv (x, nu, pq, &n);
			a = pq[kind];
		}
		stack[first][i] = (float)a;
	}
}

void grd_R2 (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: R2 2 1 R2 = A^2 + B^2.  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;
	if (gmtdefs.verbose && constant[prev] && factor[prev] == 0.0) fprintf (stderr, "%s: Warning, operand one == 0!\n", GMT_program);
	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand two == 0!\n", GMT_program);
	if (constant[prev]) factor[prev] *= factor[prev];
	if (constant[last]) factor[last] *= factor[last];
	for (i = 0; i < info->nm; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][i] * stack[prev][i];
		b = (constant[last]) ? factor[last] : stack[last][i] * stack[last][i];
		stack[prev][i] = (float)(a + b);
	}
}

void grd_R2D (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: R2D 1 1 Convert Radians to Degrees.  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = R2D * factor[last];
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : R2D * stack[last][i]);
}

void grd_RAND (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: RAND 2 1 Uniform random values between A and B.  */
{
	GMT_LONG i, prev;
	double a = 0.0, b = 0.0;

	prev = last - 1;
	if (constant[prev]) a = factor[prev];
	if (constant[last]) b = factor[last];
	for (i = 0; i < info->nm; i++) {
		if (!constant[prev]) a = (double)stack[prev][i];
		if (!constant[last]) b = (double)stack[last][i];
		stack[prev][i] = (float)(a + GMT_rand () * (b - a));
	}
}

void grd_RINT (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: RINT 1 1 rint (A) (nearest integer).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = rint (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : rint ((double)stack[last][i]));
}

void grd_ROTX (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: ROTX 2 1 Rotate A by the (constant) shift B in x-direction.  */
{
	GMT_LONG i, j, prev, shift, *node, k, nx;
	float *z;

	/* Shift grid A by the x-shift B.  B must be a constant */

	prev = last - 1;
	if (!constant[last]) {
		fprintf (stderr, "%s: DX shift (B) must be a constant in ROTX!\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	shift = irint (factor[last] / info->header.x_inc);	/* Shift of nodes */

	if (constant[prev] || !shift) return;	/* Trivial since A is a constant or shift is zero */
	if (shift < 0) shift += info->header.nx;	/* Same thing */
	nx = (size_t) info->header.nx;
	/* Set up permutation vector */

	node = (GMT_LONG *) GMT_memory (VNULL, nx, sizeof (GMT_LONG), GMT_program);
	z =  (float *) GMT_memory (VNULL, nx, sizeof (float), GMT_program);
	for (i = 0; i < info->header.nx; i++) node[i] = (i + shift) % info->header.nx;	/* Move by shift but rotate around */
	for (j = k = 0; j < info->header.ny; j++, k += nx) {	/* For each row */
		for (i = 0; i < info->header.nx; i++) z[node[i]] = stack[prev][k+i];
		memcpy ((void *)&stack[prev][k], (void *)z, nx * sizeof (float));
	}
	GMT_free ((void *)z);
	GMT_free ((void *)node);
}

void grd_ROTY (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: ROTY 2 1 Rotate A by the (constant) shift B in y-direction.  */
{
	GMT_LONG i, j, prev, shift, *node, nx;
	float *z;

	/* Shift grid A by the y-shift B.  B must be a constant */

	prev = last - 1;
	if (!constant[last]) {
		fprintf (stderr, "%s: DY shift (B) must be a constant in ROTY!\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	shift = irint (factor[last] / info->header.y_inc);	/* Shift of nodes */

	if (constant[prev] || !shift) return;	/* Trivial since A is a constant or shift is zero */
	if (shift < 0) shift += info->header.ny;	/* Same thing */
	nx = (size_t) info->header.nx;
	/* Set up permutation vector */

	node = (GMT_LONG *) GMT_memory (VNULL, (size_t)info->header.ny, sizeof (GMT_LONG), GMT_program);
	z =  (float *) GMT_memory (VNULL, (size_t)info->header.ny, sizeof (float), GMT_program);
	for (j = 0; j < info->header.ny; j++) node[j] = (j + info->header.ny - shift) % info->header.ny;	/* Move by shift but rotate around */
	for (i = 0; i < info->header.nx; i++) {	/* For each column */
		for (j = 0; j < info->header.ny; j++) z[node[j]] = stack[prev][j*nx+i];
		for (j = 0; j < info->header.ny; j++) stack[prev][j*nx+i] = z[j];
	}
	GMT_free ((void *)z);
	GMT_free ((void *)node);
}

void grd_SDIST (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: SDIST 2 1 Spherical (Great circle|geodesic) distance (in km) between grid nodes and stack lon,lat (A, B).  */
{
	GMT_LONG i, j, prev, k;
	double a, b;

	GMT_distance_func = (PFD) ((GMT_IS_SPHERICAL) ? GMT_great_circle_dist_km : GMT_geodesic_dist_km);
	prev = last - 1;
	for (k = j = i = 0; k < info->nm; k++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][k];
		b = (constant[last]) ? factor[last] : stack[last][k];
		stack[prev][k] = (float)(GMT_distance_func) (a, b, (double)info->grd_x[i], (double)info->grd_y[j]);
		i++;
		if (i == info->header.nx) i = 0, j++;
	}
}

void grd_SAZ (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: SAZ 2 1 Spherical azimuth from grid nodes to stack x,y.  */
/* Azimuth from grid ones to stack point */
{
	grd_AZ_sub (info, stack, constant, factor, last, FALSE);
}

void grd_SBAZ (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: SBAZ 2 1 Spherical backazimuth from grid nodes to stack x,y.  */
/* Azimuth from stack point to grid ones (back azimuth) */
{
	grd_AZ_sub (info, stack, constant, factor, last, TRUE);
}

void grd_AZ_sub (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG reverse)
{
	GMT_LONG i, j, prev, k;
	double x0 = 0.0, y0 = 0.0, az;
	PFD azimuth_func;

	azimuth_func = (PFD) ((GMT_IS_SPHERICAL) ? GMT_az_backaz_sphere : GMT_az_backaz_geodesic);
	prev = last - 1;
	if (constant[prev]) x0 = factor[prev];
	if (constant[last]) y0 = factor[last];
	for (k = j = i = 0; k < info->nm; k++) {
		if (!constant[prev]) x0 = (double)stack[prev][i];
		if (!constant[last]) y0 = (double)stack[last][i];
		az = (*azimuth_func) (info->grd_x[i], info->grd_y[j], x0, y0, reverse);
		while (az < -180.0) az += 360.0;
		while (az > +180.0) az -= 360.0;
		stack[prev][k] = (float)az;
		i++;
		if (i == info->header.nx) i = 0, j++;
	}
}

void grd_SEC (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: SEC 1 1 sec (A) (A in radians).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = (1.0 / cos (factor[last]));
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : (1.0 / cos ((double)stack[last][i])));
}

void grd_SECD (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: SECD 1 1 sec (A) (A in degrees).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = (1.0 / cosd (factor[last]));
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : (1.0 / cosd ((double)stack[last][i])));
}

void grd_SIGN (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: SIGN 1 1 sign (+1 or -1) of A.  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand == 0!\n", GMT_program);
	if (constant[last]) a = copysign (1.0, factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : copysign (1.0, (double)stack[last][i]));
}

void grd_SIN (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: SIN 1 1 sin (A) (A in radians).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = sin (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : sin ((double)stack[last][i]));
}

void grd_SINC (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: SINC 1 1 sinc (A) (sin (pi*A)/(pi*A)).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = GMT_sinc (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : GMT_sinc ((double)stack[last][i]));
}

void grd_SIND (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: SIND 1 1 sin (A) (A in degrees).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = sind (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : sind ((double)stack[last][i]));
}

void grd_SINH (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: SINH 1 1 sinh (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = sinh (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : sinh ((double)stack[last][i]));
}

void grd_SKEW (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: SKEW 1 1 Skewness of A.  */
{
	GMT_LONG i, n = 0;
	double mean = 0.0, sum2 = 0.0, skew = 0.0, delta;
	float f_skew;

	if (constant[last]) {	/* Trivial case */
		for (i = 0; i < info->nm; i++) stack[last][i] = GMT_f_NaN;
		return;
	}

	/* Use Welford (1962) algorithm to compute mean and corrected sum of squares */
	for (i = 0; i < info->nm; i++) {
		if (GMT_is_fnan (stack[last][i])) continue;
		n++;
		delta = (double)stack[last][i] - mean;
		mean += delta / n;
		sum2 += delta * ((double)stack[last][i] - mean);
	}
	if (n > 1) {
		for (i = 0; i < info->nm; i++) {
			if (GMT_is_fnan (stack[last][i])) continue;
			delta = (double)stack[last][i] - mean;
			skew += pow (delta, 3.0);
		}
		sum2 /= (n - 1);
		skew /= n * pow (sum2, 1.5);
		f_skew = (float)skew;
	}
	else
		f_skew = GMT_f_NaN;
	for (i = 0; i < info->nm; i++) stack[last][i] = f_skew;
}

void grd_SQR (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: SQR 1 1 A^2.  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = factor[last] * factor[last];
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : stack[last][i] * stack[last][i]);
}

void grd_SQRT (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: SQRT 1 1 sqrt (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && factor[last] < 0.0) fprintf (stderr, "%s: Warning, operand one < 0!\n", GMT_program);
	if (constant[last]) a = sqrt (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : sqrt ((double)stack[last][i]));
}

void grd_STD (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: STD 1 1 Standard deviation of A.  */
{
	GMT_LONG i, n = 0;
	double mean = 0.0, sum2 = 0.0, delta;

	if (constant[last]) {	/* Trivial case */
		memset ((void *)stack[last], 0, info->nm * sizeof (float));
		return;
	}

	/* Use Welford (1962) algorithm to compute mean and corrected sum of squares */
	for (i = 0; i < info->nm; i++) {
		if (GMT_is_fnan (stack[last][i])) continue;
		n++;
		delta = (double)stack[last][i] - mean;
		mean += delta / n;
		sum2 += delta * ((double)stack[last][i] - mean);
	}
	sum2 = (n > 1) ? sqrt (sum2 / (n - 1)) : 0.0;
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)sum2;
}

void grd_STEP (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: STEP 1 1 Heaviside step function: H(A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = factor[last];
	for (i = 0; i < info->nm; i++) {
		if (!constant[last]) a = (double)stack[last][i];
		if (a == 0.0)
			stack[last][i] = (float)0.5;
		else
			stack[last][i] = (float)((a < 0.0) ? 0.0 : 1.0);
	}
}

void grd_STEPX (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: STEPX 1 1 Heaviside step function in x: H(x-A).  */
{
	GMT_LONG i;
	double a;

	for (i = 0; i < info->nm; i++) {
		a = info->grd_x[i%info->header.nx] - ((constant[last]) ? factor[last] : stack[last][i]);
		if (a == 0.0)
			stack[last][i] = (float)0.5;
		else
			stack[last][i] = (float)((a < 0.0) ? 0.0 : 1.0);
	}
}

void grd_STEPY (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: STEPY 1 1 Heaviside step function in y: H(y-A).  */
{
	GMT_LONG i;
	double a;

	for (i = 0; i < info->nm; i++) {
		a = info->grd_y[i/info->header.nx] - ((constant[last]) ? factor[last] : stack[last][i]);
		if (a == 0.0)
			stack[last][i] = (float)0.5;
		else
			stack[last][i] = (float)((a < 0.0) ? 0.0 : 1.0);
	}
}

void grd_SUB (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: SUB 2 1 A - B.  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;
	if (gmtdefs.verbose && constant[prev] && factor[prev] == 0.0) fprintf (stderr, "%s: Warning, operand one == 0!\n", GMT_program);
	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand two == 0!\n", GMT_program);
	for (i = 0; i < info->nm; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][i];
		b = (constant[last]) ? factor[last] : stack[last][i];
		stack[prev][i] = (float)(a - b);
	}
}

void grd_TAN (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: TAN 1 1 tan (A) (A in radians).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = tan (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : tan ((double)stack[last][i]));
}

void grd_TAND (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: TAND 1 1 tan (A) (A in degrees).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = tand (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : tand ((double)stack[last][i]));
}

void grd_TANH (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: TANH 1 1 tanh (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = tanh (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : tanh ((double)stack[last][i]));
}

void grd_TN (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: TN 2 1 Chebyshev polynomial Tn(-1<t<+1,n), with t = A, and n = B.  */
{
	GMT_LONG i, prev, n;
	double a = 0.0, t;

	prev = last - 1;
	for (i = 0; i < info->nm; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][i];
		n = irint ((double)((constant[last]) ? factor[last] : stack[last][i]));
		GMT_chebyshev (a, n, &t);
		stack[prev][i] = (float)t;
	}
}

void grd_TCRIT (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: TCRIT 2 1 Critical value for Student's t-distribution, with alpha = A and n = B.  */
{
	GMT_LONG i, b, prev;
	double a;

	prev = last - 1;
	if (gmtdefs.verbose && constant[prev] && factor[prev] == 0.0) fprintf (stderr, "%s: Warning, operand one == 0 for TCRIT!\n", GMT_program);
	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand two == 0 for TCRIT!\n", GMT_program);
	for (i = 0; i < info->nm; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][i];
		b = irint ((double)((constant[last]) ? factor[last] : stack[last][i]));
		stack[prev][i] = (float)GMT_tcrit (a, (double)b);
	}
}

void grd_TDIST (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: TDIST 2 1 Student's t-distribution A(t,n), with t = A, and n = B.  */
{
	GMT_LONG i, b, prev;
	double a, prob;

	prev = last - 1;
	if (gmtdefs.verbose && constant[prev] && factor[prev] == 0.0) fprintf (stderr, "%s: Warning, operand one == 0 for TDIST!\n", GMT_program);
	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand two == 0 for TDIST!\n", GMT_program);
	for (i = 0; i < info->nm; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][i];
		b = irint ((double)((constant[last]) ? factor[last] : stack[last][i]));
		(void) GMT_student_t_a (a, (GMT_LONG)b, &prob);
		stack[prev][i] = (float)prob;
	}
}

void grd_UPPER (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: UPPER 1 1 The highest (maximum) value of A.  */
{
	GMT_LONG i;
	float high;

	if (constant[last]) {	/* Trivial case */
		for (i = 0; i < info->nm; i++) stack[last][i] = (float)factor[last];
		return;
	}

	for (i = 0, high = -FLT_MAX; i < info->nm; i++) {
		if (GMT_is_fnan (stack[last][i])) continue;
		if (stack[last][i] > high) high = stack[last][i];
	}
	for (i = 0; i < info->nm; i++) if (!GMT_is_fnan (stack[last][i])) stack[last][i] = high;
}

void grd_XOR (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: XOR 2 1 0 if A == NaN and B == NaN, NaN if B == NaN, else A.  */
{
	GMT_LONG i, prev;
	double a = 0.0, b = 0.0;

	prev = last - 1;
	if (constant[prev]) a = factor[prev];
	if (constant[last]) b = factor[last];
	for (i = 0; i < info->nm; i++) {
		if (!constant[prev]) a = stack[prev][i];
		if (!constant[last]) b = stack[last][i];
		stack[prev][i] = (float)((GMT_is_fnan (a) && GMT_is_fnan (b)) ? 0.0 : (GMT_is_fnan (b) ? GMT_f_NaN : a));
	}
}

void grd_Y0 (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: Y0 1 1 Bessel function of A (2nd kind, order 0).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = y0 (fabs (factor[last]));
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : y0 (fabs ((double)stack[last][i])));
}

void grd_Y1 (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: Y1 1 1 Bessel function of A (2nd kind, order 1).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = y1 (fabs (factor[last]));
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : y1 (fabs ((double)stack[last][i])));
}

void grd_YLM (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: YLM 2 2 Re and Im orthonormalized spherical harmonics degree A order B.  */
{
	grd_YLM_sub (info, stack, constant, factor, last, TRUE);
}

void grd_YLMg (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: YLMg 2 2 Cos and Sin normalized spherical harmonics degree A order B (geophysical convention).  */
{
	grd_YLM_sub (info, stack, constant, factor, last, FALSE);
}

void grd_YLM_sub (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG ortho)
{
	/* Returns geophysical normalization, unless M < 0, then orthonormalized form */
	GMT_LONG i, j, prev, L, M, k;
	double x, z, P, C, S;

	prev = last - 1;
	if (!(constant[prev] && constant[last])) {
		fprintf (stderr, "%s: L and M must be constants in YLM[g]!\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	L = irint (factor[prev]);
	M = irint (factor[last]);
	z = GMT_abs (M) * D2R;	/* abs() just in case routine is called with -M to add (-1)^M */

	for (k = j = 0; j < info->header.ny; j++) {	/* For each latitude */

		x = sind (info->grd_y[j]);	/* Plm takes cos(colatitude) = sin(latitude) */
		P = GMT_plm_bar (L, M, x, ortho);
		if (M == 0) {
			for (i = 0; i < info->header.nx; i++, k++) {
				stack[prev][k] = (float)P;
				stack[last][k] = 0.0;
			}
		}
		else {
			for (i = 0; i < info->header.nx; i++, k++) {
				sincos (z * info->grd_x[i], &S, &C);
				stack[prev][k] = (float)(P * C);
				stack[last][k] = (float)(P * S);
			}
		}
	}
}

void grd_YN (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: YN 2 1 Bessel function of A (2nd kind, order B).  */
{
	GMT_LONG i, prev;
	int order = 0;
	double b = 0.0;
	GMT_LONG simple = FALSE;

	prev = last - 1;
	if (gmtdefs.verbose && constant[prev] && factor[prev] == 0.0) fprintf (stderr, "%s: Warning, argument = 0 for YN!\n", GMT_program);
	if (constant[last]) {
		if (gmtdefs.verbose && factor[last] < 0.0) fprintf (stderr, "%s: Warning, order < 0 for YN!\n", GMT_program);
		if (gmtdefs.verbose && (rint(factor[last]) != factor[last])) fprintf (stderr, "%s: Warning, order not an integer for YN!\n", GMT_program);
		order = irint (fabs (factor[last]));
		if (constant[prev]) {
			b = yn (order, fabs (factor[prev]));
			simple = TRUE;
		}
	}
	for (i = 0; i < info->nm; i++) {
		if (simple)
			stack[prev][i] = (float)b;
		else {
			if (!constant[last]) order = irint (fabs ((double)stack[last][i]));
			stack[last][i] = (float)yn (order, fabs ((double)stack[prev][i]));
		}
	}
}

void grd_ZCRIT (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: ZCRIT 1 1 Critical value for the normal-distribution, with alpha = A.  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = GMT_zcrit (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : GMT_zcrit (a));
}

void grd_ZDIST (struct GRDMATH_INFO *info, float *stack[], GMT_LONG *constant, double *factor, GMT_LONG last)
/*OPERATOR: ZDIST 1 1 Cumulative normal-distribution C(x), with x = A.  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = GMT_zdist (factor[last]);
	for (i = 0; i < info->nm; i++) stack[last][i] = (float)((constant[last]) ? a : GMT_zdist (a));
}

/* ---------------------- end operator functions --------------------- */

GMT_LONG decode_argument (char *txt, double *value, struct GMT_HASH *H)
{
	GMT_LONG i, expect, check = GMT_IS_NAN;
	GMT_LONG possible_number = FALSE;
	char copy[GMT_LONG_TEXT];
	double tmp = 0.0;
	GMT_LONG get_operator (char *choice, struct GMT_HASH *H);

	if (!strcmp (txt, "=")) return GRDMATH_ARG_IS_SAVE;	/* Time to save stack; next arg is filename */

	/* Check if argument is operator */

	if ((i = get_operator (txt, H)) >= GRDMATH_ARG_IS_OPERATOR) return (i);

	/* Next look for symbols with special meaning */

	if (!(strcmp (txt, "PI") && strcmp (txt, "pi"))) return GRDMATH_ARG_IS_PI;
	if (!(strcmp (txt, "E") && strcmp (txt, "e"))) return GRDMATH_ARG_IS_E;
	if (!strcmp (txt, "EULER")) return GRDMATH_ARG_IS_EULER;
	if (!strcmp (txt, "XMIN")) return GRDMATH_ARG_IS_XMIN;
	if (!strcmp (txt, "XMAX")) return GRDMATH_ARG_IS_XMAX;
	if (!strcmp (txt, "XINC")) return GRDMATH_ARG_IS_XINC;
	if (!strcmp (txt, "NX")) return GRDMATH_ARG_IS_NX;
	if (!strcmp (txt, "YMIN")) return GRDMATH_ARG_IS_YMIN;
	if (!strcmp (txt, "YMAX")) return GRDMATH_ARG_IS_YMAX;
	if (!strcmp (txt, "YINC")) return GRDMATH_ARG_IS_YINC;
	if (!strcmp (txt, "NY")) return GRDMATH_ARG_IS_NY;
	if (!strcmp (txt, "X")) return GRDMATH_ARG_IS_X_MATRIX;
	if (!strcmp (txt, "Xn")) return GRDMATH_ARG_IS_x_MATRIX;
	if (!strcmp (txt, "Y")) return GRDMATH_ARG_IS_Y_MATRIX;
	if (!strcmp (txt, "Yn")) return GRDMATH_ARG_IS_y_MATRIX;

	/* Preliminary test-conversion to a number */

	sscanf (txt, "%[^=?]", copy);	/* Exclude netcdf 3/-D grid extensions */
	if (!GMT_not_numeric (copy)) {	/* Only check if we are not sure this is NOT a number */
		expect = (strchr (copy, 'T')) ? GMT_IS_ABSTIME : GMT_IS_UNKNOWN;	/* Watch out for dateTclock-strings */
		check = GMT_scanf (copy, expect, &tmp);
		possible_number = TRUE;
	}

	/* Determine if argument is file. But first strip off suffix */

	if (!GMT_access (copy, R_OK)) {	/* Yes it is */
		if (check != GMT_IS_NAN && possible_number) fprintf (stderr, "%s: WARNING: Your argument %s is both a file and a number.  File is selected\n", GMT_program, txt);
		return GRDMATH_ARG_IS_FILE;
	}

	if (check != GMT_IS_NAN) {	/* OK it is a number */
		*value = tmp;
		return GRDMATH_ARG_IS_NUMBER;
	}

	if (txt[0] == '-') {	/* Probably a bad commandline option */
		fprintf (stderr, "%s: ERROR: Option %s not recognized\n", GMT_program, txt);
		exit (EXIT_FAILURE);
	}

	fprintf (stderr, "%s: GMT SYNTAX ERROR: %s is not a number, operator or file name\n", GMT_program, txt);
	exit (EXIT_FAILURE);
	return (0);	/* Dummy return to satisfy some compilers */
}

GMT_LONG get_operator (char *choice, struct GMT_HASH *H) {
	GMT_LONG op;
	/* Returns -1 if not a registered operator */

	if (!strncmp (choice, "GDIST", (size_t)5)) choice[0] = 'S';	/* Changing GDIST to SDIST for backwards compatibility */

	op = GMT_hash_lookup (choice, H, GRDMATH_N_OPERATORS, GRDMATH_N_OPERATORS);

	if (op < 0 && strlen (choice) == 1) {	/* Check for old-style operators */

		switch (choice[0]) {
			case '+':
				op = ADD;
				break;
			case '-':
				op = SUB;
				break;
			case 'x':
				op = MUL;
				break;
			case '/':
				op = DIV;
				break;
			case '^':
				op = POW;
				break;
		}
	}

	return (op);
}

void *New_grdmath_Ctrl () {	/* Allocate and initialize a new control structure */
	struct GRDMATH_CTRL *C;

	C = (struct GRDMATH_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GRDMATH_CTRL), "New_grdmath_Ctrl");

	/* Initialize values whose defaults are not 0/FALSE/NULL */

	return ((void *)C);
}

void Free_grdmath_Ctrl (struct GRDMATH_CTRL *C) {	/* Deallocate control structure */
	GMT_free ((void *)C);
}

#include "grdmath.h"
