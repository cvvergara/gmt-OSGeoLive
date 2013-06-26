/*--------------------------------------------------------------------
 *	$Id: gmtmath.c 9923 2012-12-18 20:45:53Z pwessel $
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
/*
 * gmtmath.c is a reverse polish calculator that operates on table files
 * (and constants) and perform basic mathematical operations
 * on them like add, multiply, etc.
 * Some operators only work on one operand (e.g., log, exp)
 *
 * Author:	Paul Wessel
 * Date:	10-NOV-1998
 * Version:	1.0 based on 3.1 grdmath and sample1d
 *		3.1.2 PW 03/26/99 Added -H capability
 *		3.1.2 PW 04/07/99 Added -Q for quick scalar calculator
 *		3.3.2 PW 09/10/99 Added erfinv
 *		3.3.3 PW 12/10/99 Added RAND and NRAND
 *		3.3.4 PW 03/14/00 Fix problems with EXCH and DUP for constants
 *		3.3.5 PW 07/14/00 Changed STEP to STEPT and added new STEP
 *		3.3.5 PW 07/31/00 Empty -T means there are no time columns.
 *		3.3.6 PW 08/16/00 Add LT, LE, EQ, GE, GT, NAN, CHIDIST, FDIST, TDIST
 *		3.3.6 PW 08/18/00 Add INT, ISNAN, XOR, MODE, MAD, LMSSCL, SUM
 *		3.3.6 PW 08/18/00 Added -S to just return first row
 *		3.3.6 PW 08/23/00 Added LOWER and UPPER
 *		3.4   PW 03/01/01
 *		4.0   PW 11/28/01 Added Critical values for Chi2, F, T, and Z distributions
 *		      PW 11/30/01 Added LSQFIT to solve a general least squares system
 *		      PW 12/08/01 Allow for -T<filename> with irregular time coordinates
 *		      PW 01/27/04 Added SINC and NEQ
 *		      PW 03/24/04 Added ROOTS
 *		      PW 03/28/04 Added FLIPUD, ROTT
 *		      PW 07/01/04 Added LRAND
 *		      PW 07/17/04 Added LOG2
 *		      PW 07/27/05 Added TN (Chebyshev)
 *		      PW 08/05/05 Added -I to output descending times (reverse t)
 *		      PW 08/10/05 Added Tn for normalized T coordinates [-1 | +1 ]
 *		      PW 09/07/05 Added CORRCOEFF
 *		      PW 02/16/06 If STDIN is given, read <stdin> and put it on the stack
 *				  Also added -F to select which columns should be output [all]
 *		      PW 03/22/06 Added CPOISS
 *		      PW 03/25/06 Removed use of global variables, added ZDIST
 *		      PW 07/06/07 Added PSI, PV, QV, COT, COTD, ACOT, SEC, SECD, ASEC, CSC, CSCD, ACSC
 *		      PW 09/21/07 Added KURT, SKEW, PQUANT, EULER
 *		25-SEP-2007 RS: Added PLMg.
 *		07-DEC-2007 PW: Added TMIN, TMAX, TINC, N as special constants.
 *		13-AUG-2008 PW: Added NOT, fixed BCs for D2DT2 when some nodes are NaN. Also tolerate NaNs in data columns.
 *		08-OCT-2008 RS: Added INRANGE (which was added to grdmath on 28-MAR-2004, but not to gmtmath)
 */

#include "gmt.h"

#define GMTMATH_ARG_IS_OPERATOR	 0
#define GMTMATH_ARG_IS_FILE	-1
#define GMTMATH_ARG_IS_NUMBER	-2
#define GMTMATH_ARG_IS_PI	-3
#define GMTMATH_ARG_IS_E	-4
#define GMTMATH_ARG_IS_EULER	-5
#define GMTMATH_ARG_IS_TMIN	-6
#define GMTMATH_ARG_IS_TMAX	-7
#define GMTMATH_ARG_IS_TINC	-8
#define GMTMATH_ARG_IS_N	-9
#define GMTMATH_ARG_IS_T_MATRIX	-10
#define GMTMATH_ARG_IS_t_MATRIX	-11

#define GMTMATH_STACK_SIZE	100

struct GMTMATH_CTRL {	/* All control options for this program (except common args) */
	/* active is TRUE if the option has been activated */
	struct A {	/* -A<t_f(t).d> */
		GMT_LONG active;
		char *file;
	} A;
	struct C {	/* -C<cols> */
		GMT_LONG active;
		GMT_LONG *cols;
	} C;
	struct F {	/* -F<cols> */
		GMT_LONG active;
		GMT_LONG *cols;
	} F;
	struct I {	/* -I */
		GMT_LONG active;
	} I;
	struct N {	/* -N<n_col>/<t_col> */
		GMT_LONG active;
		GMT_LONG ncol, tcol;
	} N;
	struct Q {	/* -Q */
		GMT_LONG active;
	} Q;
	struct S {	/* -S[f|l] */
		GMT_LONG active;
		GMT_LONG mode;
	} S;
	struct T {	/* -T[<tmin/tmax/t_inc>] | -T<file> */
		GMT_LONG active;
		GMT_LONG notime;
		GMT_LONG mode;	/* = 1 if t_inc really is number of desired nodes */
		double min, max, inc;
		char *file;
	} T;
};

struct TABLE_HEADER {
	GMT_LONG n_row;		/* Number of time-nodes (rows) */
	GMT_LONG n_col;		/* Number of columns */
	double t_min;		/* Minimum t value */
	double t_max;		/* Maximum t value */
	double t_inc;		/* t increment */
};

struct GMTMATH_INFO {
	GMT_LONG irregular;	/* TRUE if t_inc varies */
	GMT_LONG roots_found;	/* TRUE if roots have been solved for */
	GMT_LONG very_first;	/* TRUE the very first time */
	GMT_LONG *skip_row;	/* TRUE for each row to be skipped */
	GMT_LONG n_roots;	/* Number of roots found */
	GMT_LONG r_col;		/* The column used to find roots */
	double *t_coordinates;	/* Array with t values */
	double *tn_coordinates;	/* Array with t normalized values [-1,+1] */
	struct TABLE_HEADER header;
	char **segment_header;	/* List of segment headers to output when -m is in effect */
	char head_record[BUFSIZ];
};

#include "gmtmath_def.h"

/* Helper functions */
void new_table (double ***s, GMT_LONG n_col, GMT_LONG n);
void table_PVQV (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row, GMT_LONG kind);

int main (int argc, char **argv)
{
	GMT_LONG i, j, arg, op = 0, nstack = 0, new_stack = -1, last_arg, ok = 1;
	GMT_LONG use_t_col = 0, first_last_all = 0, nm = 0;
	GMT_LONG consumed_operands[GMTMATH_N_OPERATORS], produced_operands[GMTMATH_N_OPERATORS];

	GMT_LONG constant[GMTMATH_STACK_SIZE], error = FALSE, set_t = FALSE, got_t_from_file = FALSE;
	GMT_LONG set_q = FALSE, read_stdin = FALSE, t_check_required = TRUE;

	double **stack[GMTMATH_STACK_SIZE], **rhs = NULL, **tmp_stack = NULL, **stdin_stack = NULL;

	double factor[GMTMATH_STACK_SIZE], t_noise, value, off, scale, special_symbol[GMTMATH_ARG_IS_PI-GMTMATH_ARG_IS_N+1];

	char *outfile = CNULL, file[BUFSIZ];

	struct TABLE_HEADER tbl[GMTMATH_STACK_SIZE], stdin_header, rhs_header;

	struct GMT_HASH *p = NULL, *current = NULL, localhashnode[GMTMATH_N_OPERATORS];
	struct GMTMATH_INFO info;
	struct GMTMATH_CTRL *Ctrl = NULL;

	FILE *fp = NULL;

	PFV call_operator[GMTMATH_N_OPERATORS];

	GMT_LONG decode_argument (char *txt, double *value, struct GMT_HASH *H);
	void gmtmath_init(PFV ops[], GMT_LONG n_args[], GMT_LONG n_out[]);
	void GMT_read_table (struct GMTMATH_INFO *info, char *file, struct TABLE_HEADER *h, double ***p, GMT_LONG t_col, GMT_LONG init);
	void GMT_write_table (struct GMTMATH_INFO *info, char *file, struct TABLE_HEADER *h, double **p, GMT_LONG first_last_all, GMT_LONG cols[]);
	void decode_columns (char *txt, GMT_LONG *skip, GMT_LONG n_col, GMT_LONG t_col, GMT_LONG mode);
	void solve_LSQFIT (struct GMTMATH_INFO *info, double **stack[], GMT_LONG last, GMT_LONG n_col, GMT_LONG n_row, GMT_LONG skip[], char *file);
	void *New_gmtmath_Ctrl (), Free_gmtmath_Ctrl (struct GMTMATH_CTRL *C);

	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct GMTMATH_CTRL *) New_gmtmath_Ctrl ();		/* Allocate and initialize defaults in a new control structure */

	memset ((void *)&info, 0, sizeof (struct GMTMATH_INFO));
	info.very_first = TRUE;

	if (argc == 2 && !strcmp (argv[1], "-")) error = GMT_give_synopsis_and_exit = TRUE;

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "gmtmath %s - Reverse Polish Notation (RPN) calculator for table data\n\n", GMT_VERSION);
		fprintf (stderr, "usage: gmtmath [-A<t_f(t).d>] [-C<cols>] [-F<cols>] [%s] [-I] [-N<n_col>/<t_col>] [-Q]\n", GMT_H_OPT);
		fprintf (stderr, "\t[-S[f|l]] [-T[<tmin/tmax/t_inc>[+]]] [-V] [-%s]\n\t[%s] [%s] A B op C op ... = [outfile]\n\n", GMT_b_OPT, GMT_f_OPT, GMT_m_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\tA, B, etc are table files, constants, or symbols (see below).\n");
		fprintf (stderr, "\tTo read stdin give filename as STDIN (which can appear more than once).\n");
		fprintf (stderr, "\tThe stack can hold up to %d entries (given enough memory).\n", GMTMATH_STACK_SIZE);
		fprintf (stderr, "\tTrigonometric operators expect radians.\n");
		fprintf (stderr, "\tThe operators and number of input and output arguments:\n\n");
		fprintf (stderr, "\tName    #args   Returns\n");
		fprintf (stderr, "\t-----------------------\n");
#include "gmtmath_explain.h"
		fprintf (stderr, "\n\tThe special symbols are:\n\n");
		fprintf (stderr, "\t  PI	= 3.1415926...\n");
		fprintf (stderr, "\t  E	= 2.7182818...\n");
		fprintf (stderr, "\t  EULER	= 0.5772156...\n");
		fprintf (stderr, "\t  TMIN, TMAX, or TINC	= the corresponding constant\n");
		fprintf (stderr, "\t  N	= number of records\n");
		fprintf (stderr, "\t  T	= table with t-coordinates\n");
		fprintf (stderr, "\t  Tn	= table with normalized [-1 to +1] t-coordinates\n");
		fprintf (stderr, "\n\tOPTIONS:\n\n");
		fprintf (stderr, "\t-A Requires -N and will initialize table with file containing t and f(t) only.\n");
		fprintf (stderr, "\t   t goes into column <t_col> while f(t) goes into column <n_col> - 1.\n");
		fprintf (stderr, "\t-C change which columns to operate on [Default is all except time].\n");
		fprintf (stderr, "\t   -C reverts to the default, -Cr toggles current settings, and -Ca selects all columns.\n");
		fprintf(stderr,"\t-F Give comma-separated list of desired columns or ranges to output (0 is first column) [Default is all].\n");
		GMT_explain_option ('H');
		fprintf (stderr, "\t-I Reverses the output sequence into descending order [ascending].\n");
		fprintf (stderr, "\t-N sets the number of columns and the id of the time column (0 is first) [2/0].\n");
		fprintf (stderr, "\t-Q quick scalar calculator. Shorthand for -Ca -N1/0 -T0/0/1.\n");
		fprintf (stderr, "\t-S Only write first row upon completion of calculations [write all rows].\n");
		fprintf (stderr, "\t   Optionally, append l for last row or f for first row [Default].\n");
		fprintf (stderr, "\t-T Set domain from t_min to t_max in steps of t_inc.\n");
		fprintf (stderr, "\t   Append + to t_inc to indicate the number of points instead.\n");
		fprintf (stderr, "\t   If a filename is given instead we read t coordinates from first column.\n");
		fprintf (stderr, "\t   If no domain is given we assume no time, i.e., only data columns are present.\n");
		fprintf (stderr, "\t   This choice also implies -Ca.\n");
		GMT_explain_option ('V');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		GMT_explain_option ('o');
		GMT_explain_option ('n');
		GMT_explain_option ('f');
		GMT_explain_option ('m');
		exit (EXIT_FAILURE);
	}

	for (i = 1, error = TRUE; error && i < argc; i++) if (argv[i][0] == '=' && argv[i][1] == '\0') error = FALSE;
	if (error) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Usage is <operations> = [outfile]\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	GMT_hash_init (localhashnode, operator, GMTMATH_N_OPERATORS, GMTMATH_N_OPERATORS);

	for (i = 0; i < GMTMATH_STACK_SIZE; i++) {
		constant[i] = FALSE;
		factor[i] = 0.0;
		tbl[i].n_col = tbl[i].n_row = 0;
		stack[i] = (double **)NULL;
	}

	last_arg = 1;
	while (last_arg < argc && !(argv[last_arg][0] == '=' && argv[last_arg][1] == '\0')) last_arg++;	/* Find position of the = argument */
	i = last_arg + 1;	/* This is normally the position of the output file name (if one was given) */
	while (i < argc && argv[i][0] == '-') i++;	/* Skip past any options between = and output file */
	outfile = (i < argc) ? argv[i] : NULL;
	GMT_io.skip_if_NaN[GMT_X] = GMT_io.skip_if_NaN[GMT_Y] = FALSE;	/* Turn off default GMT NaN-handling of x/y (e.g. lon/lat columns) */

	/* Must first scan command line for -b, -f, -H, -m|M, and -N switches before reading any file */

	for (arg = 1; arg < last_arg; arg++) {
		if (!strncmp (argv[arg], "-H", (size_t)2)) error += GMT_parse_common_options (argv[arg], NULL, NULL, NULL, NULL);
		if (!strncmp (argv[arg], "-M", (size_t)2)) GMT_parse_m_option (&argv[arg][2]);
		if (!strncmp (argv[arg], "-m", (size_t)2)) GMT_parse_m_option (&argv[arg][2]);
		if (!strncmp (argv[arg], "-m", (size_t)2)) GMT_parse_m_option (&argv[arg][2]);
		if (!strncmp (argv[arg], "-b", (size_t)2)) error += GMT_parse_common_options (argv[arg], NULL, NULL, NULL, NULL);
		if (!strncmp (argv[arg], "-f", (size_t)2)) error += GMT_parse_common_options (argv[arg], NULL, NULL, NULL, NULL);
		if (!strcmp (argv[arg], "-T")) t_check_required = FALSE;	/* Turn off default GMT NaN-handling in t column */
		if (!strncmp (argv[arg], "-N", (size_t)2)) {	/* Correctly determine which column is time, needed if -T is set */
			Ctrl->N.active = TRUE;
			sscanf (&argv[arg][2], "%" GMT_LL "d/%" GMT_LL "d", &Ctrl->N.ncol, &Ctrl->N.tcol);
		}
	}
	GMT_io.skip_if_NaN[Ctrl->N.tcol] = t_check_required;	/* Determines if the t-column may have NaNs */
	
	/* Get header from one file so we can allocate space */

	for (arg = 1; nm == 0 && arg < last_arg; arg++) {

		if (argv[arg][0] == '-' && argv[arg][1] != 0) continue;	/* Command line option */
		if (decode_argument (argv[arg], &value, localhashnode) != GMTMATH_ARG_IS_FILE) continue;

		strcpy (file, argv[arg]);
		if (!strcmp (file, "STDIN")) {
			GMT_read_table (&info, argv[arg], &stdin_header, &stdin_stack, Ctrl->N.tcol, TRUE);
			memcpy ((void *)&info.header, (void *)&stdin_header, sizeof (struct TABLE_HEADER));
			read_stdin = TRUE;
		}
		else
			GMT_read_table (&info, argv[arg], &info.header, &tmp_stack, Ctrl->N.tcol, TRUE);

		nm = info.header.n_row * info.header.n_col;
		Ctrl->N.ncol = info.header.n_col;
		got_t_from_file = 1;
		use_t_col = Ctrl->N.tcol;
	}

	/* Scan command line for -A, -I, -T, -Q, -S, -V */

	for (arg = 1; arg < last_arg; arg++) {
		if (argv[arg][0] == '-') {

			switch (argv[arg][1]) {

				case 'V':
					error += GMT_parse_common_options (argv[arg], NULL, NULL, NULL, NULL);
					break;

				case 'A':	/* y(x) table for LSQFIT operations */
					Ctrl->A.active = TRUE;
					Ctrl->A.file = strdup (&argv[arg][2]);
					break;
				case 'F':
					decode_columns (&argv[arg][2], Ctrl->F.cols, GMT_MAX_COLUMNS, 0, 1);
					Ctrl->F.active = TRUE;
					break;
				case 'I':
					Ctrl->I.active = TRUE;
					break;
				case 'Q':	/* Quick for -Ca -N1/0 -T0/0/1 */
					Ctrl->Q.active = TRUE;
					break;
				case 'S':	/* Only want one row (first or last) */
					Ctrl->S.active = TRUE;
					if (!argv[arg][2] || argv[arg][2] == 'F' || argv[arg][2] == 'f')
						Ctrl->S.mode = -1;
					else if (argv[arg][2] == 'L' || argv[arg][2] == 'l')
						Ctrl->S.mode = +1;
					else {
						fprintf (stderr, "%s: GMT SYNTAX ERROR: Syntax is -S[f|l]\n", GMT_program);
						exit (EXIT_FAILURE);
					}
					break;
				case 'T':	/* Either get a file with time coordinate or a min/max/dt setting */
					Ctrl->T.active = TRUE;
					if (argv[arg][2] && !(GMT_access (&argv[arg][2], R_OK)))	/* Argument given and file can be opened */
						Ctrl->T.file = strdup (&argv[arg][2]);
					else {
						if (sscanf (&argv[arg][2], "%lf/%lf/%lf", &info.header.t_min, &info.header.t_max, &info.header.t_inc) != 3) Ctrl->T.notime = TRUE;
						if (argv[arg][strlen(argv[arg])-1] == '+') Ctrl->T.mode = 1;
					}
					break;
			}
		}
	}
	if (Ctrl->Q.active && (Ctrl->T.active || Ctrl->N.active || Ctrl->C.active)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Cannot use -T, -N, or -C when -Q has been set\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
	if (Ctrl->N.active && (Ctrl->N.ncol <= 0 || Ctrl->N.tcol < 0 || Ctrl->N.tcol >= Ctrl->N.ncol)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  -N must have positive n_cols and 0 <= t_col < n_col\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	set_t = (Ctrl->T.active && !Ctrl->T.file && !Ctrl->T.notime);
	if (nm && set_t) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Cannot use -T when data files are specified\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (Ctrl->A.active) {
		GMT_read_table (&info, Ctrl->A.file, &rhs_header, &rhs, 0, TRUE);	/* Always store as t and f(t) in cols 0 and 1 */
		if (rhs_header.n_col != 2) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR:  -A must take a file with 2 (t,f(t)) columns\n", GMT_program);
			exit (EXIT_FAILURE);
		}
	}
	if (Ctrl->N.active) {
		GMT_io.skip_if_NaN[GMT_X] = GMT_io.skip_if_NaN[GMT_Y] = FALSE;
		GMT_io.skip_if_NaN[Ctrl->N.tcol] = TRUE;
	}
	if (Ctrl->Q.active) {
		Ctrl->N.ncol = 1;
		Ctrl->N.tcol = 0;
		Ctrl->N.active = set_t = set_q = TRUE;
		info.header.t_min = info.header.t_max = 0;
		info.header.t_inc = 1.0;
	}
	if (Ctrl->T.active) {
		if (Ctrl->T.file) {	/* Got a filename */
			if (got_t_from_file) {
				fprintf (stderr, "%s: GMT SYNTAX ERROR:  Cannot use -T when data files are specified\n", GMT_program);
				exit (EXIT_FAILURE);
			}
			GMT_read_table (&info, Ctrl->T.file, &info.header, &tmp_stack, 0, TRUE);
			use_t_col = 0;
			got_t_from_file = 2;
		}
	}
	if (set_t && !set_q) {
		if (Ctrl->T.mode == 1) {	/* Got n, now set t_inc */
			info.header.t_inc = (info.header.t_max - info.header.t_min) / (info.header.t_inc - 1.0);
		}
		switch (GMT_minmaxinc_verify (info.header.t_min, info.header.t_max, info.header.t_inc, GMT_SMALL)) {
			case 1:
				fprintf (stderr, "%s: GMT SYNTAX ERROR -T:  (max - min) is not a whole multiple of inc\n", GMT_program);
				exit (EXIT_FAILURE);
				break;
			case 2:
				if (info.header.t_inc != 1.0) {	/* Allow for somebody explicitly saying -T0/0/1 */
					fprintf (stderr, "%s: GMT SYNTAX ERROR -T:  (max - min) is <= 0\n", GMT_program);
					exit (EXIT_FAILURE);
				}
				break;
			case 3:
				fprintf (stderr, "%s: GMT SYNTAX ERROR -T:  inc is <= 0\n", GMT_program);
				exit (EXIT_FAILURE);
				break;
			default:	/* OK */
				break;
		}

		info.header.n_row = irint ((info.header.t_max - info.header.t_min) / info.header.t_inc) + 1;
		info.header.n_col = Ctrl->N.ncol;
		nm = info.header.n_row * info.header.n_col;
	}
	first_last_all = Ctrl->S.mode;

	if (Ctrl->A.active) {	/* Get number of rows and time from the file, but not n_cols (that takes -N, which defaults to 2) */
		info.header.n_row = rhs_header.n_row;
		info.header.n_col = Ctrl->N.ncol;
		info.header.t_min = rhs_header.t_min;
		info.header.t_max = rhs_header.t_max;
		info.header.t_inc = rhs_header.t_inc;
		nm = info.header.n_row * info.header.n_col;
	}
	if (set_q) info.header.n_row = info.header.n_col = nm = 1;
	if (Ctrl->T.file) {
		info.header.n_col = Ctrl->N.ncol;
		nm = info.header.n_row * info.header.n_col;
	}
	if (nm == 0) {	/* Neither a file nor -T given; must read data from stdin */
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Expression must contain at least one table file or -T [and -N]\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	new_table (&stack[0], info.header.n_col, info.header.n_row);

	if (!Ctrl->T.notime && info.header.n_col > 1) Ctrl->C.cols[Ctrl->N.tcol] = (set_q) ? FALSE : TRUE;
	Ctrl->F.cols = (GMT_LONG *) GMT_memory ((void *)Ctrl->F.cols, (size_t)info.header.n_col, sizeof (GMT_LONG), GMT_program);
	if (!Ctrl->F.active) for (i = 0; i < info.header.n_col; i++) Ctrl->F.cols[i] = TRUE;

	/* Get t vector */

	info.t_coordinates = (double *) GMT_memory (VNULL, (size_t)info.header.n_row, sizeof (double), GMT_program);
	info.tn_coordinates = (double *) GMT_memory (VNULL, (size_t)info.header.n_row, sizeof (double), GMT_program);
	if (read_stdin) {
		memcpy ((void *)info.t_coordinates, (void *)stdin_stack[use_t_col], (size_t)(info.header.n_row * sizeof (double)));
		for (i = 1; i < info.header.n_row && (info.skip_row[i] || info.skip_row[i-1]); i++);	/* Find the first real two records in a row */
		info.header.t_inc = (i == info.header.n_row) ? GMT_d_NaN : stdin_stack[use_t_col][i] - stdin_stack[use_t_col][i-1];
		t_noise = fabs (GMT_SMALL * info.header.t_inc);
		for (i = 1; i < info.header.n_row && !info.irregular; i++) if (fabs (fabs (info.t_coordinates[i] - info.t_coordinates[i-1]) - fabs (info.header.t_inc)) > t_noise && !(info.skip_row[i] || info.skip_row[i-1])) info.irregular = TRUE;
	}
	else if (got_t_from_file) {
		memcpy ((void *)info.t_coordinates, (void *)tmp_stack[use_t_col], (size_t)(info.header.n_row * sizeof (double)));
		for (i = 1; i < info.header.n_row && (info.skip_row[i] || info.skip_row[i-1]); i++);	/* Find the first real two records in a row */
		info.header.t_inc = (i == info.header.n_row) ? GMT_d_NaN : tmp_stack[use_t_col][i] - tmp_stack[use_t_col][i-1];
		t_noise = fabs (GMT_SMALL * info.header.t_inc);
		for (i = 1; i < info.header.n_row && !info.irregular; i++) if (fabs (fabs (info.t_coordinates[i] - info.t_coordinates[i-1]) - fabs (info.header.t_inc)) > t_noise && !(info.skip_row[i] || info.skip_row[i-1])) info.irregular = TRUE;
		j = (got_t_from_file == 1) ? info.header.n_col : 1;
		for (i = 0; i < j; i++) GMT_free ((void *)tmp_stack[i]);
		GMT_free ((void *)tmp_stack);
	}
	else {
		for (i = 0; i < info.header.n_row; i++) info.t_coordinates[i] = (i == (info.header.n_row-1)) ? info.header.t_max: info.header.t_min + i * info.header.t_inc;
		t_noise = fabs (GMT_SMALL * info.header.t_inc);
	}
	off = 0.5 * (info.t_coordinates[info.header.n_row-1] + info.t_coordinates[0]);
	scale = 2.0 / (info.t_coordinates[info.header.n_row-1] - info.t_coordinates[0]);
	if (Ctrl->I.active) for (i = 0; i < info.header.n_row/2; i++) d_swap (info.t_coordinates[i], info.t_coordinates[info.header.n_row-1-i]);	/* Reverse time-series */
	for (i = 0; i < info.header.n_row; i++) info.tn_coordinates[i] = (info.t_coordinates[i] - off) * scale;
	if (!read_stdin) memcpy ((void *)stack[0][Ctrl->N.tcol], (void *)info.t_coordinates, (size_t)(info.header.n_row * sizeof (double)));
	if (Ctrl->A.active) {
		memcpy ((void *)stack[0][Ctrl->N.ncol-1], (void *)rhs[1], (size_t)(info.header.n_row * sizeof (double)));
		GMT_free ((void *)rhs[0]);
		GMT_free ((void *)rhs[1]);
		nstack = 1;
	}
	else
		nstack = 0;

	special_symbol[GMTMATH_ARG_IS_PI-GMTMATH_ARG_IS_PI] = M_PI;
	special_symbol[GMTMATH_ARG_IS_PI-GMTMATH_ARG_IS_E] = M_E;
	special_symbol[GMTMATH_ARG_IS_PI-GMTMATH_ARG_IS_EULER] = M_EULER;
	special_symbol[GMTMATH_ARG_IS_PI-GMTMATH_ARG_IS_TMIN] = info.header.t_min;
	special_symbol[GMTMATH_ARG_IS_PI-GMTMATH_ARG_IS_TMAX] = info.header.t_max;
	special_symbol[GMTMATH_ARG_IS_PI-GMTMATH_ARG_IS_TINC] = info.header.t_inc;
	special_symbol[GMTMATH_ARG_IS_PI-GMTMATH_ARG_IS_N] = (double)info.header.n_row;
	
	if (!info.skip_row) info.skip_row = (GMT_LONG *) GMT_memory (VNULL, (size_t)info.header.n_row, sizeof (GMT_LONG), GMT_program);	/* All FALSE if generated by -T */

	gmtmath_init (call_operator, consumed_operands, produced_operands);


	for (arg = 1; !error && arg < last_arg; arg++) {

		/* First check if we should skip optional arguments */

		if (!(strncmp (argv[arg], "-T", (size_t)2) && strncmp (argv[arg], "-b", (size_t)2) && strncmp (argv[arg], "-f", (size_t)2) && strncmp (argv[arg], "-N", (size_t)2))) continue;
		if (!(strncmp (argv[arg], "-H", (size_t)2) && strncmp (argv[arg], "-Q", (size_t)2) && strncmp (argv[arg], "-S", (size_t)2) && strncmp (argv[arg], "-V", (size_t)2))) continue;
		if (!(strncmp (argv[arg], "-A", (size_t)2) && strncmp (argv[arg], "-I", (size_t)2) && strncmp (argv[arg], "-F", (size_t)2) && strncmp (argv[arg], "-M", (size_t)2) && strncmp (argv[arg], "-m", (size_t)2))) continue;

		if (!strncmp (argv[arg], "-C", (size_t)2)) {	/* Change affected columns */
			decode_columns (&argv[arg][2], Ctrl->C.cols, Ctrl->N.ncol, Ctrl->N.tcol, 0);
			continue;
		}

		op = decode_argument (argv[arg], &value, localhashnode);

		if (op != GMTMATH_ARG_IS_FILE && !GMT_access(argv[arg], R_OK)) fprintf (stderr, "%s Warning: The number or operator %s may be confused with an existing file %s!\n", GMT_program, argv[arg], argv[arg]);

		if (op < GMTMATH_ARG_IS_OPERATOR) {	/* File name or factor */

			if (nstack == GMTMATH_STACK_SIZE) {	/* Stack overflow */
				error = TRUE;
				continue;
			}

			if (op == GMTMATH_ARG_IS_NUMBER) {
				constant[nstack] = TRUE;
				factor[nstack] = value;
				error = FALSE;
				if (gmtdefs.verbose) fprintf (stderr, "%g ", factor[nstack]);
				nstack++;
				continue;
			}
			else if (op <= GMTMATH_ARG_IS_PI && op >= GMTMATH_ARG_IS_N) {
				constant[nstack] = TRUE;
				factor[nstack] = special_symbol[GMTMATH_ARG_IS_PI-op];
				if (gmtdefs.verbose) fprintf (stderr, "%g ", factor[nstack]);
				nstack++;
				continue;
			}

			/* Here we need a matrix */

			if (!stack[nstack]) new_table (&stack[nstack], info.header.n_col, info.header.n_row);

			constant[nstack] = FALSE;

			if (op == GMTMATH_ARG_IS_T_MATRIX) {	/* Need to set up matrix of t-values */
				if (Ctrl->T.notime) {
					fprintf (stderr, "%s: T is not defined for plain data files!\n", GMT_program);
					exit (EXIT_FAILURE);
				}
				if (gmtdefs.verbose) fprintf (stderr, "T ");
				for (j = 0; j < info.header.n_col; j++) memcpy ((void *)stack[nstack][j], (void *)info.t_coordinates, (size_t)(info.header.n_row * sizeof (double)));
			}
			else if (op == GMTMATH_ARG_IS_t_MATRIX) {	/* Need to set up matrix of normalized t-values */
				if (Ctrl->T.notime) {
					fprintf (stderr, "%s: Tn is not defined for plain data files!\n", GMT_program);
					exit (EXIT_FAILURE);
				}
				if (gmtdefs.verbose) fprintf (stderr, "Tn ");
				for (j = 0; j < info.header.n_col; j++) if (j != Ctrl->N.tcol) memcpy ((void *)stack[nstack][j], (void *)info.tn_coordinates, (size_t)(info.header.n_row * sizeof (double)));
			}
			else if (op == GMTMATH_ARG_IS_FILE) {		/* Filename given */
				if (!strcmp (argv[arg], "STDIN")) {	/* stdin file */
					if (gmtdefs.verbose) fprintf (stderr, "<stdin> ");
					memcpy ((void *)&tbl[nstack], (void *)&stdin_header, sizeof (struct TABLE_HEADER));
					for (j = 0; j < info.header.n_col; j++) memcpy ((void *)stack[nstack][j], (void *)stdin_stack[j], (size_t)(info.header.n_row * sizeof (double)));
				}
				else {
					if (gmtdefs.verbose) fprintf (stderr, "%s ", argv[arg]);
					GMT_read_table (&info, argv[arg], &tbl[nstack], &stack[nstack], Ctrl->N.tcol, FALSE);
				}
				if (tbl[nstack].n_row != info.header.n_row || tbl[nstack].n_col != info.header.n_col) {
					fprintf (stderr, "%s: tables not of same size!\n", GMT_program);
					exit (EXIT_FAILURE);
				}
				else if (!Ctrl->T.notime && (fabs (tbl[nstack].t_min - info.header.t_min) > t_noise || fabs (tbl[nstack].t_max - info.header.t_max) > t_noise)) {
					fprintf (stderr, "%s: tables do not cover the same domain!\n", GMT_program);
					exit (EXIT_FAILURE);
				}
			}
			nstack++;
			continue;
		}

		/* Here we have an operator */

		if (!strncmp (argv[arg], "ROOTS", (size_t)5) && !((arg+1) == last_arg && argv[arg+1][0] == '=')) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR:  Only = may follow operator ROOTS\n", GMT_program);
			exit (EXIT_FAILURE);
		}

		if ((new_stack = nstack - consumed_operands[op] + produced_operands[op]) >= GMTMATH_STACK_SIZE) {
			error = TRUE;
			continue;
		}

		if (nstack < consumed_operands[op]) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR:  Operation \"%s\" requires %ld operands\n", GMT_program, operator[op], consumed_operands[op]);
			exit (EXIT_FAILURE);
		}

		if (gmtdefs.verbose) fprintf (stderr, "%s ", operator[op]);

		for (i = produced_operands[op] - consumed_operands[op]; i > 0; i--) {
			if (stack[nstack+i-1])	continue;

			/* Must make space for more */

			new_table (&stack[nstack+i-1], info.header.n_col, info.header.n_row);
		}

		/* If operators operates on constants only we may have to make space as well */

		for (j = 0, i = nstack - consumed_operands[op]; j < produced_operands[op]; j++, i++) {
			if (constant[i] && !stack[i]) new_table (&stack[i], info.header.n_col, info.header.n_row);
		}

		if (!strcmp (operator[op], "LSQFIT")) {	/* Special case, solve LSQ system and exit */
			solve_LSQFIT (&info, stack, nstack - 1, Ctrl->N.ncol, info.header.n_row, Ctrl->C.cols, outfile);
			exit (EXIT_SUCCESS);
		}

		for (j = 0; j < Ctrl->N.ncol; j++) {
			if (Ctrl->C.cols[j]) continue;
			(*call_operator[op]) (&info, stack, constant, factor, nstack - 1, j, info.header.n_row);	/* Do it */
		}

		nstack = new_stack;

		for (i = 1; i <= produced_operands[op]; i++)
			constant[nstack-i] = FALSE;	/* Now filled with table */
	}

	if (error && !ok) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Unable to decode constant %s (File not found?)\n", GMT_program, argv[i-1]);
		exit (EXIT_FAILURE);
	}

	if (error) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Stack overflow (%s)\n", GMT_program, argv[i-1]);
		exit (EXIT_FAILURE);
	}

	if (gmtdefs.verbose) {
		(outfile) ? fprintf (stderr, "= %s", outfile) : fprintf (stderr, "= <stdout>");
	}

	if (new_stack < 0 && constant[0]) {	/* Only a constant provided, set table accordingly */
		for (j = 0; j < info.header.n_col; j++) {
			if (Ctrl->C.cols[j]) continue;
			for (i = 0; i < info.header.n_row; i++) stack[0][j][i] = factor[0];
		}
	}

	if (gmtdefs.verbose) fprintf (stderr, "\n");

	if (info.roots_found) {	/* Special treatment of root finding */
		if (outfile) {
			if ((fp = GMT_fopen (outfile, "w")) == NULL) {
				fprintf (stderr, "%s: GMT ERROR:  Could not create file (%s)\n", GMT_program, outfile);
				exit (EXIT_FAILURE);
			}
		}
		else
			fp = GMT_stdout;
		for (i = 0; i < info.n_roots; i++) GMT_output (fp, 1, &stack[0][info.r_col][i]);
		if (fp != GMT_stdout) GMT_fclose (fp);
	}
	else
		GMT_write_table (&info, outfile, &info.header, stack[0], first_last_all, Ctrl->F.cols);

	for (i = 0; i < GMTMATH_STACK_SIZE; i++) if (stack[i]) {
		for (j = 0; j < info.header.n_col; j++) GMT_free ((void *)stack[i][j]);
		GMT_free ((void *)stack[i]);
	}
	if (read_stdin) {
		for (j = 0; j < info.header.n_col; j++) GMT_free ((void *)stdin_stack[j]);
		GMT_free ((void *)stdin_stack);
	}
	GMT_free ((void *)info.t_coordinates);
	GMT_free ((void *)info.tn_coordinates);
	if (info.segment_header) GMT_free ((void *)info.segment_header);
	if (info.skip_row) GMT_free ((void *)info.skip_row);
	for (i = 0; i < GMTMATH_N_OPERATORS; i++) {
		p = localhashnode[i].next;
		while ((current = p)) {
			p = p->next;
			GMT_free ((void *)current);
		}
	}

	if (nstack > 1) fprintf (stderr, "%s: Warning: %ld more operands left on the stack!\n", GMT_program, nstack-1);

	Free_gmtmath_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void GMT_read_table (struct GMTMATH_INFO *info, char *file, struct TABLE_HEADER *h, double ***p, GMT_LONG t_col, GMT_LONG init)
{	/* init is TRUE when we use GMT_read_table to determine size of rows, columns etc and may have to allocate more
	 * memory on the fly.  Once size is known then future calls pass FALSE since *p already has the right size allocated.
	 * Since all files read by gmtmath MUST have the same format, rows, and cols, we check that this file matches the
	 * information obtained by the very first file (when init was TRUE) and exit otherwise. */
	GMT_LONG init_skip_row;
	GMT_LONG j, n_expected_fields, n_fields;
	GMT_LONG n_alloc_d = 0, n_alloc_h = 0;
	GMT_LONG n = 0;
	double *in, **table = NULL;
	char buffer[BUFSIZ];
	FILE *fp;

	n_expected_fields = (GMT_io.binary[GMT_IN]) ? GMT_io.ncol[GMT_IN] : GMT_MAX_COLUMNS;

	if (!strcmp (file, "STDIN")) {
		fp = GMT_stdin;
#ifdef SET_IO_MODE
		GMT_setmode (GMT_IN);
#endif
	}
	else if ((fp = GMT_fopen (file, GMT_io.r_mode)) == NULL) {
		fprintf (stderr, "%s: Error opening file %s\n", GMT_program, file);
		exit (EXIT_FAILURE);
	}

	for (j = 0; GMT_io.io_header[GMT_IN] && j < GMT_io.n_header_recs; j++) {
		GMT_fgets (buffer, BUFSIZ, fp);
		if (info->very_first && j == 0) strcpy (info->head_record, buffer);
		info->very_first = FALSE;
	}

	GMT_input (fp, &n_expected_fields, &in);
	if (!init && n_expected_fields != info->header.n_col) {
		fprintf (stderr, "%s: ERROR: Input file has different number of columns (%ld) than expected (%ld)\n", GMT_program, n_expected_fields, info->header.n_col);
		exit (EXIT_FAILURE);
	}
	if (!init) memcpy ((void *)h, (void *)&(info->header), sizeof (struct TABLE_HEADER));
	
	if (GMT_io.status & GMT_IO_EOF) {
		fprintf (stderr, "%s: Error reading 1st record of file %s\n", GMT_program, file);
		exit (EXIT_FAILURE);
	}

	/* If working on a multiple-segment data file there might potentially be tons of empty headers
	 * that we need to keep track of.  The record number (n) below refers to original records in the
	 * file; thus the table records corresponding to a multisegment header will have NaNs; this is
	 * set later given the skip_row array which is TRUE for records that are headers and FALSE otherwise.
	 * We maintain two separate n_alloc counters for the header information array and the data table
	 * since a devious data file might contain more than GMT_CHUNK of multiple segment headers before
	 * getting to the first data records.  Once a set of headers have been processed the data table is
	 * automatically updated to have at least as much memory allocated to it as the headers have.
	 */

	if (!init) table = *p;	/* *p has the correct memory preallocated */
	init_skip_row = (init && info->segment_header == NULL);	/* TRUE the first time which is the only time we need to allocate space */
	if (init_skip_row) {	/* Get memory for record information arrays */
		n_alloc_h = GMT_CHUNK;
		info->skip_row = (GMT_LONG *) GMT_memory (VNULL, n_alloc_h, sizeof (GMT_LONG), GMT_program);	/* All FALSE by default */
		if (!GMT_io.binary[GMT_IN]) info->segment_header = (char **) GMT_memory (VNULL, n_alloc_h, sizeof (char *), GMT_program);
	}
	do {	/* Process all input records, headers and all */
		while ((GMT_io.status & GMT_IO_SEGMENT_HEADER) && !(GMT_io.status & GMT_IO_EOF)) {	/* Encountered a multisegment header */
			if (init_skip_row) {
				info->skip_row[n] = TRUE;	/* Flag as a header record */
				if (!GMT_io.binary[GMT_IN]) {	/* Only ascii files have useful header info - store in array */
					info->segment_header[n] = strdup (GMT_io.segment_header);	/* Save this header */
				}
			}
			n_fields = GMT_input (fp, &n_expected_fields, &in);
			n++;
			if (init_skip_row && n == n_alloc_h) {	/* Need to allocate more header information memory */
				n_alloc_h <<= 1;
				info->skip_row = (GMT_LONG *) GMT_memory ((void *)info->skip_row, n_alloc_h, sizeof (GMT_LONG), GMT_program);
				if (!GMT_io.binary[GMT_IN]) info->segment_header = (char **) GMT_memory ((void *)info->segment_header, n_alloc_h, sizeof (char *), GMT_program);
			}
		}
		if (init && (n_alloc_d == 0 || (n_alloc_d < n_alloc_h) || (n_alloc_d < n))) {	/* Either first time we have read a data record or if we have read more headers that we have allocated data space for so far */
			n_alloc_d = MAX(GMT_CHUNK, n_alloc_h);		/* Since we might have read a lot of multisegment headers before getting here */
			while (n >= n_alloc_d) n_alloc_d <<= 1;	/* Since when init_skip_row is FALSE then n_alloc_h = 0 but n might be large */
			h->n_col = n_expected_fields;
			new_table (&table, n_expected_fields, (GMT_LONG)n_alloc_d);
		}
		for (j = 0; j < h->n_col; j++) table[j][n] = in[j];	/* Copy current record values to the current row of the table */
		if (init)
			info->skip_row[n] = FALSE;	/* Flag as a data record */
		else if (info->skip_row[n]) {	/* Does not match previous files */
			fprintf (stderr, "%s: ERROR: Input files have segment headers in different places!\n", GMT_program);
			exit (EXIT_FAILURE);
		}
		n++;
		if (init && n == n_alloc_d) {	/* Time to allocate more data and header record arrays */
			n_alloc_d <<= 1;
			if (init_skip_row) {
				n_alloc_h = (n_alloc_h == 0) ? GMT_CHUNK : (n_alloc_h << 1);
				info->skip_row = (GMT_LONG *) GMT_memory ((void *)info->skip_row, n_alloc_h, sizeof (GMT_LONG), GMT_program);
				if (!GMT_io.binary[GMT_IN]) info->segment_header = (char **) GMT_memory ((void *)info->segment_header, n_alloc_h, sizeof (char *), GMT_program);
			}
			for (j = 0; j < h->n_col; j++) table[j] = (double *) GMT_memory ((void *)table[j], (size_t)n_alloc_d, sizeof (double), GMT_program);
		}
	} while ((n_fields = GMT_input (fp, &n_expected_fields, &in)) >= 0 && !(GMT_io.status & GMT_IO_EOF));

	if (init ){
		h->t_min = table[t_col][0];
		h->t_max = table[t_col][n-1];
		h->t_inc = (h->t_max - h->t_min) / (n-1);
		*p = table;
		h->n_row = n;
		if (init_skip_row) {	/* Finalize size of header arrays */
			info->skip_row = (GMT_LONG *) GMT_memory ((void *)info->skip_row, (size_t)n, sizeof (GMT_LONG), GMT_program);
			if (!GMT_io.binary[GMT_IN]) info->segment_header = (char **) GMT_memory ((void *)info->segment_header, (size_t)n, sizeof (char *), GMT_program);
		}
	}
	else if (h->n_row != n) {
		fprintf (stderr, "%s: ERROR: Input file has different number of rows (%ld) than expected (%ld)\n", GMT_program, n, info->header.n_row);
		exit (EXIT_FAILURE);
	}
	/* Fill in NaNs for multisegment headers which are flagged in the skip_row array */
	for (n = 0; n < h->n_row; n++) {
		if (!info->skip_row[n]) continue;
		for (j = 0; j < h->n_col; j++) table[j][n] = GMT_d_NaN;
	}
	if (fp != GMT_stdin) GMT_fclose (fp);
}

void GMT_write_table (struct GMTMATH_INFO *info, char *file, struct TABLE_HEADER *h, double **p, GMT_LONG first_last_all, GMT_LONG cols[])
{	/* first_last_all will write first [-1], last [+1], or all [0] records */
	GMT_LONG i, j, k, start, stop;
	double *out;
	FILE *fp;

	if (!file) {
		fp = GMT_stdout;
#ifdef SET_IO_MODE
		GMT_setmode (GMT_OUT);
#endif
	}
	else if ((fp = GMT_fopen (file, GMT_io.w_mode)) == NULL) {
		fprintf (stderr, "%s: Error creating file %s\n", GMT_program, file);
		exit (EXIT_FAILURE);
	}

	out = (double *) GMT_memory (VNULL, (size_t)h->n_col, sizeof (double), GMT_program);

	if (GMT_io.io_header[GMT_OUT] && !GMT_io.binary[GMT_OUT]) {
		GMT_fputs (info->head_record, fp);
		for (i = 1; i < GMT_io.n_header_recs; i++) GMT_fputs ("# gmtmath header record\n", fp);
	}

	start = (first_last_all == +1) ? h->n_row - 1 : 0;
	stop = (first_last_all == -1) ? 1 : h->n_row;
	for (i = start; i < stop; i++) {
		if (info->skip_row[i]) {
			if (!(GMT_io.binary[GMT_IN] || GMT_io.binary[GMT_OUT])) strcpy (GMT_io.segment_header, info->segment_header[i]);
			GMT_write_segmentheader (fp, h->n_col);
		}
		else {
			for (j = k = 0; j < h->n_col; j++) if (cols[j]) out[k++] = p[j][i];
			GMT_output (fp, k, out);
		}
	}

	if (fp != GMT_stdout) GMT_fclose (fp);
	GMT_free ((void *)out);
}

void new_table (double ***s, GMT_LONG n_col, GMT_LONG n)
{	/* First time it is called for a table the the s pointer is NULL */
	GMT_LONG j;
	double **p;
	p = (double **) GMT_memory ((void *)(*s), (size_t)n_col, sizeof (double *), GMT_program);
	for (j = 0; j < n_col; j++) p[j] = (double *) GMT_memory ((void *)p[j], (size_t)n, sizeof (double), GMT_program);
	*s = p;
}

void decode_columns (char *txt, GMT_LONG *skip, GMT_LONG n_col, GMT_LONG t_col, GMT_LONG mode)
{
	GMT_LONG i, start, stop, pos, T, F;
	char p[BUFSIZ];

	if (mode == 0) {	/* For marking columns to skip */
		T = TRUE;
		F = FALSE;
	}
	else {			/* For marking columns to output */
		T = FALSE;
		F = TRUE;
	}
	if (mode == 0 && !txt[0]) {	/* Reset to default */
		for (i = 0; i < n_col; i++) skip[i] = F;
		skip[t_col] = T;
	}
	else if (mode == 0 && txt[0] == 'r' && txt[1] == '\0') {	/* Reverse all settings */
		for (i = 0; i < n_col; i++) skip[i] = !skip[i];
	}
	else if (mode == 0 && txt[0] == 'a') {	/* Select all columns */
		for (i = 0; i < n_col; i++) skip[i] = F;
	}
	else {	/* Set the selected columns */
		for (i = 0; i < n_col; i++) skip[i] = T;
		pos = 0;
		while ((GMT_strtok (txt, ",", &pos, p))) {
			if (strchr (p, '-'))
				sscanf (p, "%" GMT_LL "d-%" GMT_LL "d", &start, &stop);
			else {
				sscanf (p, "%" GMT_LL "d", &start);
				stop = start;
			}
			stop = MIN (stop, n_col-1);
			for (i = start; i <= stop; i++) skip[i] = F;
		}
	}
}


/* -----------------------------------------------------------------
 *              Definitions of all operator functions
 * -----------------------------------------------------------------*/

void table_ABS (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: ABS 1 1 abs (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand == 0!\n", GMT_program);
	if (constant[last]) a = fabs (factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : fabs (stack[last][col][i]);
}

void table_ACOS (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: ACOS 1 1 acos (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && fabs (factor[last]) > 1.0) fprintf (stderr, "%s: Warning, |operand| > 1 for ACOS!\n", GMT_program);
	if (constant[last]) a = d_acos (factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : d_acos (stack[last][col][i]);
}

void table_ACOSH (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: ACOSH 1 1 acosh (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && fabs (factor[last]) > 1.0) fprintf (stderr, "%s: Warning, operand < 1 for ACOSH!\n", GMT_program);
	if (constant[last]) a = acosh (factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : acosh (stack[last][col][i]);
}

void table_ACSC (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: ACSC 1 1 acsc (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && fabs (factor[last]) > 1.0) fprintf (stderr, "%s: Warning, |operand| > 1 for ACSC!\n", GMT_program);
	if (constant[last]) a = d_asin (1.0 / factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : d_asin (1.0 / stack[last][col][i]);
}

void table_ACOT (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: ACOT 1 1 acot (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && fabs (factor[last]) > 1.0) fprintf (stderr, "%s: Warning, |operand| > 1 for ACOT!\n", GMT_program);
	if (constant[last]) a = atan (1.0 / factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : atan (1.0 / stack[last][col][i]);
}

void table_ADD (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: ADD 2 1 A + B.  */
{
	GMT_LONG i;
	GMT_LONG prev;
	double a, b;

	prev = last - 1;
	if (gmtdefs.verbose && constant[prev] && factor[prev] == 0.0) fprintf (stderr, "%s: Warning, operand one == 0!\n", GMT_program);
	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand two == 0!\n", GMT_program);
	for (i = 0; i < n_row; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][col][i];
		b = (constant[last]) ? factor[last] : stack[last][col][i];
		stack[prev][col][i] = a + b;
	}
}

void table_AND (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: AND 2 1 NaN if A and B == NaN, B if A == NaN, else A.  */
{
	GMT_LONG i;
	GMT_LONG prev;
	double a, b;

	prev = last - 1;
	for (i = 0; i < n_row; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][col][i];
		b = (constant[last]) ? factor[last] : stack[last][col][i];
		stack[prev][col][i] = (GMT_is_dnan (a)) ? b : a;
	}
}

void table_ASEC (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: ASEC 1 1 asec (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && fabs (factor[last]) > 1.0) fprintf (stderr, "%s: Warning, |operand| > 1 for ASEC!\n", GMT_program);
	if (constant[last]) a = d_acos (1.0 / factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : d_acos (1.0 / stack[last][col][i]);
}

void table_ASIN (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: ASIN 1 1 asin (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && fabs (factor[last]) > 1.0) fprintf (stderr, "%s: Warning, |operand| > 1 for ASIN!\n", GMT_program);
	if (constant[last]) a = d_asin (factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : d_asin (stack[last][col][i]);
}

void table_ASINH (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: ASINH 1 1 asinh (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = asinh (factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : asinh (stack[last][col][i]);
}

void table_ATAN (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: ATAN 1 1 atan (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = atan (factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : atan (stack[last][col][i]);
}

void table_ATAN2 (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: ATAN2 2 1 atan2 (A, B).  */
{
	GMT_LONG i;
	GMT_LONG prev;
	double a, b;

	prev = last - 1;
	if (gmtdefs.verbose && constant[prev] && factor[prev] == 0.0) fprintf (stderr, "%s: Warning, operand one == 0 for ATAN2!\n", GMT_program);
	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand two == 0 for ATAN2!\n", GMT_program);
	for (i = 0; i < n_row; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][col][i];
		b = (constant[last]) ? factor[last] : stack[last][col][i];
		stack[prev][col][i] = d_atan2 (a, b);
	}
}

void table_ATANH (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: ATANH 1 1 atanh (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && fabs (factor[last]) >= 1.0) fprintf (stderr, "%s: Warning, |operand| >= 1 for ATANH!\n", GMT_program);
	if (constant[last]) a = atanh (factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : atanh (stack[last][col][i]);
}

void table_BEI (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: BEI 1 1 bei (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = GMT_bei (fabs (factor[last]));
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : GMT_bei (fabs (stack[last][col][i]));
}

void table_BER (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: BER 1 1 ber (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = GMT_ber (fabs (factor[last]));
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : GMT_ber (fabs (stack[last][col][i]));
}

void table_CEIL (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: CEIL 1 1 ceil (A) (smallest integer >= A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = ceil (factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : ceil (stack[last][col][i]);
}

void table_CHICRIT (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: CHICRIT 2 1 Critical value for chi-squared-distribution, with alpha = A and n = B.  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;
	if (gmtdefs.verbose && constant[prev] && factor[prev] == 0.0) fprintf (stderr, "%s: Warning, operand one == 0 for CHICRIT!\n", GMT_program);
	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand two == 0 for CHICRIT!\n", GMT_program);
	for (i = 0; i < n_row; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][col][i];
		b = (constant[last]) ? factor[last] : stack[last][col][i];
		stack[prev][col][i] = GMT_chi2crit (a, b);
	}
}

void table_CHIDIST (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: CHIDIST 2 1 chi-squared-distribution P(chi2,n), with chi2 = A and n = B.  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;
	if (gmtdefs.verbose && constant[prev] && factor[prev] == 0.0) fprintf (stderr, "%s: Warning, operand one == 0 for CHIDIST!\n", GMT_program);
	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand two == 0 for CHIDIST!\n", GMT_program);
	for (i = 0; i < n_row; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][col][i];
		b = (constant[last]) ? factor[last] : stack[last][col][i];
		GMT_chi2 (a, b, &stack[prev][col][i]);
	}
}

void table_COL (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: COL 1 1 Places column A on the stack.  */
{
	GMT_LONG i, k, prev;

	if (!constant[last]) {
		fprintf (stderr, "%s: Error, argument to COL must be a constant column number (0 <= k < n_col)!\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	prev = last - 1;
	k = irint (factor[last]);
	for (i = 0; i < n_row; i++) {
		stack[last][col][i] = stack[prev][k][i];
	}
}

void table_CORRCOEFF (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: CORRCOEFF 2 1 Correlation coefficient r(A, B).  */
{
	GMT_LONG i;
	GMT_LONG prev;
	double *a, *b, coeff;

	prev = last - 1;
	if (constant[prev] && constant[last]) {	/* Correlation is undefined */
		for (i = 0; i < n_row; i++) stack[prev][col][i] = GMT_d_NaN;
		return;
	}

	if (constant[prev]) {		/* Must create the missing (constant) column */
		a = GMT_memory (VNULL, (size_t)n_row, sizeof (double), GMT_program);
		for (i = 0; i < n_row; i++) a[i] = factor[prev];
		b = stack[last][col];
	}
	else if (constant[last]) {	/* Must create the missing (constant) column */
		a = stack[prev][col];
		b = GMT_memory (VNULL, (size_t)n_row, sizeof (double), GMT_program);
		for (i = 0; i < n_row; i++) b[i] = factor[last];
	}
	else {
		a = stack[prev][col];
		b = stack[last][col];
	}
	coeff = GMT_corrcoeff (a, b, n_row, 0);
	for (i = 0; i < n_row; i++) stack[prev][col][i] = coeff;
	if (constant[prev]) GMT_free ((void *)a);
	if (constant[last]) GMT_free ((void *)b);
}


void table_COS (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: COS 1 1 cos (A) (A in radians).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = cos (factor[last]);
	for (i = 0; i < n_row; i++) {
		stack[last][col][i] = (constant[last]) ? a : cos (stack[last][col][i]);
	}
}

void table_COSD (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: COSD 1 1 cos (A) (A in degrees).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = cosd (factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : cosd (stack[last][col][i]);
}

void table_COSH (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: COSH 1 1 cosh (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = cosh (factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : cosh (stack[last][col][i]);
}

void table_COT (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: COT 1 1 cot (A) (A in radians).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = (1.0 / tan (factor[last]));
	for (i = 0; i < n_row; i++) {
		stack[last][col][i] = (constant[last]) ? a : (1.0 / tan (stack[last][col][i]));
	}
}


void table_COTD (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: COTD 1 1 cot (A) (A in degrees).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = (1.0 / tand (factor[last]));
	for (i = 0; i < n_row; i++) {
		stack[last][col][i] = (constant[last]) ? a : (1.0 / tand (stack[last][col][i]));
	}
}

void table_CSC (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: CSC 1 1 csc (A) (A in radians).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = (1.0 / sin (factor[last]));
	for (i = 0; i < n_row; i++) {
		stack[last][col][i] = (constant[last]) ? a : (1.0 / sin (stack[last][col][i]));
	}
}

void table_CSCD (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: CSCD 1 1 csc (A) (A in degrees).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = (1.0 / sind (factor[last]));
	for (i = 0; i < n_row; i++) {
		stack[last][col][i] = (constant[last]) ? a : (1.0 / sind (stack[last][col][i]));
	}
}

void table_CPOISS (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: CPOISS 2 1 Cumulative Poisson distribution F(x,lambda), with x = A and lambda = B.  */
{
	GMT_LONG i;
	GMT_LONG prev;
	double a, b;

	prev = last - 1;
	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand two == 0 for CPOISS!\n", GMT_program);
	for (i = 0; i < n_row; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][col][i];
		b = (constant[last]) ? factor[last] : stack[last][col][i];
		GMT_cumpoisson (a, b, &stack[prev][col][i]);
	}
}

void table_DDT (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: DDT 1 1 d(A)/dt Central 1st derivative.  */
{
	GMT_LONG i;
	double c, left, next_left;

	/* Central 1st difference in t */

	if (info->irregular) fprintf (stderr, "%s: Warning, DDT called on irregularly spaced data (not supported)!\n", GMT_program);
	if (gmtdefs.verbose && constant[last]) fprintf (stderr, "%s: Warning, operand to DDT is constant!\n", GMT_program);

	c = 0.5 / info->header.t_inc;
	i = 0;
	while (info->skip_row[i] && i < n_row) i++;	/* Start of first segment */
	while (i < n_row) {	/* Process each segment */
		next_left = 2.0 * stack[last][col][i] - stack[last][col][i+1];
		while (i < n_row - 1 && !info->skip_row[i+1]) {
			left = next_left;
			next_left = stack[last][col][i];
			stack[last][col][i] = (constant[last]) ? 0.0 : c * (stack[last][col][i+1] - left);
			i++;
		}
		stack[last][col][i] = (constant[last]) ? 0.0 : 2.0 * c * (stack[last][col][i] - next_left);
		i++;
		while (info->skip_row[i] && i < n_row) i++;	/* Start of next segment */
	}

}

void table_D2DT2 (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: D2DT2 1 1 d^2(A)/dt^2 2nd derivative.  */
{
	GMT_LONG i;
	double c, left, next_left;

	/* Central 2nd difference in t */

	if (info->irregular) fprintf (stderr, "%s: Warning, D2DT2 called on irregularly spaced data (not supported)!\n", GMT_program);
	if (gmtdefs.verbose && constant[last]) fprintf (stderr, "%s: Warning, operand to D2DT2 is constant!\n", GMT_program);

	c = 1.0 / (info->header.t_inc * info->header.t_inc);
	i = 0;
	while (info->skip_row[i] && i < n_row) i++;	/* Start of first segment */
	while (i < n_row) {	/* Process each segment */
		next_left = stack[last][col][i];
		stack[last][col][i] = (GMT_is_dnan (stack[last][col][i]) || GMT_is_dnan (stack[last][col][i+1])) ? GMT_d_NaN : 0.0;
		i++;
		while (i < n_row - 1 && !info->skip_row[i+1]) {
			left = next_left;
			next_left = stack[last][col][i];
			stack[last][col][i] = (constant[last]) ? 0.0 : c * (stack[last][col][i+1] - 2 * stack[last][col][i] + left);
			i++;
		}
		stack[last][col][i] = (GMT_is_dnan (stack[last][col][i]) || GMT_is_dnan (stack[last][col][i-1])) ? GMT_d_NaN : 0.0;
		i++;
		while (info->skip_row[i] && i < n_row) i++;	/* Start of next segment */
	}
}

void table_D2R (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: D2R 1 1 Converts Degrees to Radians.  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = factor[last] * D2R;
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : stack[last][col][i] * D2R;
}

void table_DILOG (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: DILOG 1 1 dilog (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = GMT_dilog (factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : GMT_dilog (stack[last][col][i]);
}

void table_DIV (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: DIV 2 1 A / B.  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;

	if (constant[last] && factor[last] == 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Cannot divide by zero\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (constant[last]) {	/* Turn divide into multiply */
		a = factor[last];	/* Save old factor */
		factor[last] = 1.0 / factor[last];
		table_MUL (info, stack, constant, factor, last, col, n_row);
		factor[last] = a;	/* Restore factor to original value */
		return;
	}

	for (i = 0; i < n_row; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][col][i];
		b = (constant[last]) ? factor[last] : stack[last][col][i];
		stack[prev][col][i] = a / b;
	}
}

void table_DUP (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: DUP 1 2 Places duplicate of A on the stack.  */
{
	GMT_LONG next, i;

	next = last + 1;
	factor[next] = factor[last];
	constant[next] = constant[last];
	if (constant[last])
		for (i = 0; i < n_row; i++) stack[next][col][i] = stack[last][col][i] = factor[next];
	else
		memcpy ((void *)stack[next][col], (void *)stack[last][col], (size_t)(n_row * sizeof (double)));
}

void table_ERF (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: ERF 1 1 Error function erf (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = erf (factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : erf (stack[last][col][i]);
}

void table_ERFC (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: ERFC 1 1 Complementary Error function erfc (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = erfc (factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : erfc (stack[last][col][i]);
}

void table_ERFINV (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: ERFINV 1 1 Inverse error function of A.  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = GMT_erfinv (factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : GMT_erfinv (stack[last][col][i]);
}

void table_EQ (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: EQ 2 1 1 if A == B, else 0.  */
{
	GMT_LONG i;
	GMT_LONG prev;
	double a, b;

	prev = last - 1;
	for (i = 0; i < n_row; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][col][i];
		b = (constant[last]) ? factor[last] : stack[last][col][i];
		stack[prev][col][i] = (double)(a == b);
	}
}

void table_EXCH (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: EXCH 2 2 Exchanges A and B on the stack.  */
{
	GMT_LONG i;
	GMT_LONG prev;

	prev = last - 1;

	for (i = 0; i < n_row; i++) {
		if (constant[last]) stack[last][col][i] = factor[last];
		if (constant[prev]) stack[prev][col][i] = factor[prev];
		d_swap (stack[last][col][i], stack[prev][col][i]);
	}
	l_swap (constant[last], constant[prev]);
	d_swap (factor[last], factor[prev]);
}

void table_EXP (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: EXP 1 1 exp (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = exp (factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : exp (stack[last][col][i]);
}

void table_FACT (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: FACT 1 1 A! (A factorial).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = GMT_factorial ((GMT_LONG)irint(factor[last]));
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : GMT_factorial ((GMT_LONG)irint(stack[last][col][i]));
}

void table_FCRIT (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: FCRIT 3 1 Critical value for F-distribution, with alpha = A, n1 = B, and n2 = C.  */
{
	GMT_LONG i, nu1, nu2, prev1, prev2;
	double alpha;

	prev1 = last - 1;
	prev2 = last - 2;
	if (gmtdefs.verbose && constant[prev2] && factor[prev2] == 0.0) fprintf (stderr, "%s: Warning, operand one == 0 for FCRIT!\n", GMT_program);
	if (gmtdefs.verbose && constant[prev1] && factor[prev1] == 0.0) fprintf (stderr, "%s: Warning, operand two == 0 for FCRIT!\n", GMT_program);
	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand three == 0 for FCRIT!\n", GMT_program);
	for (i = 0; i < n_row; i++) {
		alpha = (constant[prev2]) ? factor[prev2] : stack[prev2][col][i];
		nu1 = irint ((double)((constant[prev1]) ? factor[prev1] : stack[prev1][col][i]));
		nu2 = irint ((double)((constant[last]) ? factor[last] : stack[last][col][i]));
		stack[prev2][col][i] = GMT_Fcrit (alpha, (double)nu1, (double)nu2);
	}
}

void table_FDIST (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: FDIST 3 1 F-distribution Q(F,n1,n2), with F = A, n1 = B, and n2 = C.  */
{
	GMT_LONG i, nu1, nu2;
	GMT_LONG prev1, prev2;
	double F, chisq1, chisq2 = 1.0;

	prev1 = last - 1;
	prev2 = last - 2;
	if (gmtdefs.verbose && constant[prev1] && factor[prev1] == 0.0) fprintf (stderr, "%s: Warning, operand two == 0 for FDIST!\n", GMT_program);
	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand three == 0 for FDIST!\n", GMT_program);
	for (i = 0; i < n_row; i++) {
		F = (constant[prev2]) ? factor[prev2] : stack[prev2][col][i];
		nu1 = irint ((double)((constant[prev1]) ? factor[prev1] : stack[prev1][col][i]));
		nu2 = irint ((double)((constant[last]) ? factor[last] : stack[last][col][i]));
		/* Since GMT_f_q needs chisq1 and chisq2, we set chisq2 = 1 and solve for chisq1 */
		chisq1 = F * nu1 / nu2;
		(void) GMT_f_q (chisq1, nu1, chisq2, nu2, &stack[prev2][col][i]);
	}
}

void table_FLIPUD (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: FLIPUD 1 1 Reverse order of each column.  */
{
	GMT_LONG i, k;
	/* Reverse the order of points in a column */
	if (constant[last]) return;
	if (gmtdefs.verbose && GMT_io.multi_segments[GMT_IN]) {
		fprintf (stderr, "%s: Warning, FLIPUD on multisegment file not supported!\n", GMT_program);
		return;
	}
	for (i = 0, k = n_row-1; i < n_row/2; i++, k--) d_swap (stack[last][col][i], stack[last][col][k]);
}

void table_FLOOR (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: FLOOR 1 1 floor (A) (greatest integer <= A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = floor (factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : floor (stack[last][col][i]);
}

void table_FMOD (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: FMOD 2 1 A % B (remainder after truncated division).  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;
	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, using FMOD 0!\n", GMT_program);
	for (i = 0; i < n_row; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][col][i];
		b = (constant[last]) ? factor[last] : stack[last][col][i];
		stack[prev][col][i] = fmod (a, b);
	}
}

void table_GE (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: GE 2 1 1 if A >= B, else 0.  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;
	for (i = 0; i < n_row; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][col][i];
		b = (constant[last]) ? factor[last] : stack[last][col][i];
		stack[prev][col][i] = (double)(a >= b);
	}
}

void table_GT (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: GT 2 1 1 if A > B, else 0.  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;
	for (i = 0; i < n_row; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][col][i];
		b = (constant[last]) ? factor[last] : stack[last][col][i];
		stack[prev][col][i] = (double)(a > b);
	}
}

void table_HYPOT (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: HYPOT 2 1 hypot (A, B) = sqrt (A*A + B*B).  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;
	if (gmtdefs.verbose && constant[prev] && factor[prev] == 0.0) fprintf (stderr, "%s: Warning, operand one == 0!\n", GMT_program);
	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand two == 0!\n", GMT_program);
	for (i = 0; i < n_row; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][col][i];
		b = (constant[last]) ? factor[last] : stack[last][col][i];
		stack[prev][col][i] = hypot (a, b);
	}
}

void table_I0 (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: I0 1 1 Modified Bessel function of A (1st kind, order 0).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = GMT_i0 (factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : GMT_i0 (stack[last][col][i]);
}

void table_I1 (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: I1 1 1 Modified Bessel function of A (1st kind, order 1).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = GMT_i1 (factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : GMT_i1 (stack[last][col][i]);
}

void table_IN (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
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
	for (i = 0; i < n_row; i++) {
		if (simple)
			stack[prev][col][i] = b;
		else {
			if (!constant[last]) order = irint (fabs (stack[last][col][i]));
			stack[prev][col][i] = GMT_in (order, fabs (stack[prev][col][i]));
		}
	}
}

void table_INRANGE (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: INRANGE 3 1 1 if B <= A <= C, else 0.  */
{
	GMT_LONG i, prev1, prev2;
	GMT_LONG inrange;
	double a = 0.0, b = 0.0, c = 0.0;

	/* last is C */
	prev1 = last - 1;	/* This is B */
	prev2 = last - 2;	/* This is A */

	/* Set to 1 where B <= A <= C, 0 elsewhere, except where
	 * A, B, or C = NaN, in which case we set answer to NaN */

	if (constant[prev2]) a = (double)factor[prev2];
	if (constant[prev1]) b = (double)factor[prev1];
	if (constant[last])  c = (double)factor[last];

	for (i = 0; i < n_row; i++) {
		if (!constant[prev2]) a = stack[prev2][col][i];
		if (!constant[prev1]) b = stack[prev1][col][i];
		if (!constant[last])  c = stack[last][col][i];

		if (GMT_is_dnan (a) || GMT_is_dnan (b) || GMT_is_dnan (c)) {
			stack[prev2][col][i] = GMT_d_NaN;
			continue;
		}

		inrange = (b <= a && a <= c);
		stack[prev2][col][i] = (double)inrange;
	}
}

void table_INT (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: INT 1 1 Numerically integrate A.  */
{
	GMT_LONG i;
	double f = 0.0, left, right, sum;

	if (constant[last]) {	/* Trivial case */
		sum = factor[last] * info->header.t_inc;
		for (i = 0; i < n_row; i++) stack[last][col][i] = i * sum;
		return;
	}

	/* We use dumb trapezoidal rule - one day we will replace with more sophisticated rules */

	sum = 0.0;
	if (!info->irregular) f = 0.5 * info->header.t_inc;
	i = 0;
	while (info->skip_row[i] && i < n_row) i++;	/* Wind to first segment */
	while (i < n_row) {
		left = stack[last][col][i];
		stack[last][col][i] = sum;
		i++;
		while (i < n_row && !info->skip_row[i]) {	/* Dumb trapezoidal rule */
			if (info->irregular) f = 0.5 * (info->t_coordinates[i] - info->t_coordinates[i-1]);
			right = stack[last][col][i];
			sum += f * (left + right);
			stack[last][col][i] = sum;
			left = right;
			i++;
		}
		while (info->skip_row[i] && i < n_row) i++;	/* Wind to first segment */
	}
}

void table_INV (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: INV 1 1 1 / A.  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last] && factor[last] == 0.0) {
		fprintf (stderr, "%s: Error, Cannot take inverse of zero!\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (constant[last]) a = 1.0 / factor[last];
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : 1.0 / stack[last][col][i];
}

void table_ISNAN (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: ISNAN 1 1 1 if A == NaN, else 0.  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = (double)GMT_is_dnan (factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : (double)GMT_is_dnan (stack[last][col][i]);
}

void table_J0 (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: J0 1 1 Bessel function of A (1st kind, order 0).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = j0 (factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : j0 (stack[last][col][i]);
}

void table_J1 (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: J1 1 1 Bessel function of A (1st kind, order 1).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = j1 (fabs (factor[last]));
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : j1 (fabs (stack[last][col][i]));
}

void table_JN (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
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
	for (i = 0; i < n_row; i++) {
		if (simple)
			stack[prev][col][i] = b;
		else {
			if (!constant[last]) order = irint (fabs (stack[last][col][i]));
			stack[prev][col][i] = jn (order, fabs (stack[prev][col][i]));
		}
	}
}

void table_K0 (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: K0 1 1 Modified Kelvin function of A (2nd kind, order 0).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = GMT_k0 (factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : GMT_k0 (stack[last][col][i]);
}

void table_K1 (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: K1 1 1 Modified Bessel function of A (2nd kind, order 1).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = GMT_k1 (factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : GMT_k1 (stack[last][col][i]);
}

void table_KN (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
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
	for (i = 0; i < n_row; i++) {
		if (simple)
			stack[prev][col][i] = b;
		else {
			if (!constant[last]) order = irint (fabs (stack[last][col][i]));
			stack[prev][col][i] = GMT_kn (order, fabs (stack[prev][col][i]));
		}
	}
}

void table_KEI (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: KEI 1 1 kei (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = GMT_kei (fabs (factor[last]));
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : GMT_kei (fabs (stack[last][col][i]));
}

void table_KER (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: KER 1 1 ker (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = GMT_ker (fabs (factor[last]));
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : GMT_ker (fabs (stack[last][col][i]));
}

void table_KURT (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: KURT 1 1 Kurtosis of A.  */
{
	GMT_LONG i, n = 0;
	double mean = 0.0, sum2 = 0.0, kurt = 0.0, delta;

	if (constant[last]) {	/* Trivial case */
		for (i = 0; i < n_row; i++) stack[last][col][i] = GMT_d_NaN;
		return;
	}

	/* Use Welford (1962) algorithm to compute mean and corrected sum of squares */
	for (i = 0; i < n_row; i++) {
		if (GMT_is_dnan (stack[last][col][i])) continue;
		n++;
		delta = stack[last][col][i] - mean;
		mean += delta / n;
		sum2 += delta * (stack[last][col][i] - mean);
	}
	if (n > 1) {
		for (i = 0; i < n_row; i++) {
			if (GMT_is_dnan (stack[last][col][i])) continue;
			delta = stack[last][col][i] - mean;
			kurt += pow (delta, 4.0);
		}
		sum2 /= (n - 1);
		kurt = kurt / (n * sum2 * sum2) - 3.0;
	}
	else
		kurt = GMT_d_NaN;
	for (i = 0; i < n_row; i++) stack[last][col][i] = kurt;
}

void table_LE (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: LE 2 1 1 if A <= B, else 0.  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;
	for (i = 0; i < n_row; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][col][i];
		b = (constant[last]) ? factor[last] : stack[last][col][i];
		stack[prev][col][i] = (double)(a <= b);
	}
}

void table_LMSSCL (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: LMSSCL 1 1 LMS scale estimate (LMS STD) of A.  */
{
	GMT_LONG i, GMT_mode_selection = 0, GMT_n_multiples = 0;
	double lmsscl, mode;

	if (constant[last]) {	/* Trivial case */
		memset ((void *)stack[last][col], 0, (size_t)(n_row * sizeof (double)));
		return;
	}

	qsort ((void *)stack[last][col], (size_t)n_row, sizeof (double), GMT_comp_double_asc);
	for (i = n_row; GMT_is_dnan (stack[last][col][i-1]) && i > 1; i--);
	if (i) {
		GMT_mode (stack[last][col], i, i/2, 0, GMT_mode_selection, &GMT_n_multiples, &mode);
		GMT_getmad (stack[last][col], i, mode, &lmsscl);
	}
	else
		lmsscl = GMT_d_NaN;

	for (i = 0; i < n_row; i++) stack[last][col][i] = lmsscl;
	if (GMT_n_multiples > 0) fprintf (stderr, "%s: WARNING: %ld Multiple modes found\n", GMT_program, GMT_n_multiples);
}

void table_LOG (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: LOG 1 1 log (A) (natural log).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, argument to log = 0!\n", GMT_program);

	if (constant[last]) a = d_log (fabs (factor[last]));
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : d_log (fabs (stack[last][col][i]));
}

void table_LOG10 (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: LOG10 1 1 log10 (A) (base 10).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, argument to log10 = 0!\n", GMT_program);

	if (constant[last]) a = d_log10 (fabs (factor[last]));
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : d_log10 (fabs (stack[last][col][i]));
}

void table_LOG1P (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: LOG1P 1 1 log (1+A) (accurate for small A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && factor[last] < 0.0) fprintf (stderr, "%s: Warning, argument to log1p < 0!\n", GMT_program);

	if (constant[last]) a = d_log1p (fabs (factor[last]));
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : d_log1p (fabs (stack[last][col][i]));
}

void table_LOG2 (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: LOG2 1 1 log2 (A) (base 2).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, argument to log2 = 0!\n", GMT_program);

	if (constant[last]) a = d_log (fabs (factor[last])) * M_LN2_INV;
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : d_log (fabs (stack[last][col][i])) * M_LN2_INV;
}

void table_LOWER (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: LOWER 1 1 The lowest (minimum) value of A.  */
{
	GMT_LONG i;
	double low;

	if (constant[last]) {	/* Trivial case */
		for (i = 0; i < n_row; i++) stack[last][col][i] = factor[last];
		return;
	}

	for (i = 0, low = DBL_MAX; i < n_row; i++) {
		if (GMT_is_dnan (stack[last][col][i])) continue;
		if (stack[last][col][i] < low) low = stack[last][col][i];
	}
	if (low == DBL_MAX) low = GMT_d_NaN;	/* All rows were NaN */
	for (i = 0; i < n_row; i++) stack[last][col][i] = low;
}

void table_LRAND (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: LRAND 2 1 Laplace random noise with mean A and std. deviation B.  */
{
	GMT_LONG i, prev;
	double a = 0.0, b = 0.0;

	prev = last - 1;
	if (constant[prev]) a = factor[prev];
	if (constant[last]) b = factor[last];
	for (i = 0; i < n_row; i++) {
		if (!constant[prev]) a = stack[prev][col][i];
		if (!constant[last]) b = stack[last][col][i];
		stack[prev][col][i] = a + b * GMT_lrand ();
	}
}

void table_LSQFIT (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: LSQFIT 1 0 Let current table be [A | b]; return least squares solution x = A \\ b.  */
{
	/* Dummy routine needed since the automatically generated include file will have table_LSQFIT
	 * with these parameters just like any other function.  However, when we find LSQFIT we will
	 * instead call solve_LSQFIT which can be found at the end of these functions */
}

void table_LT (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: LT 2 1 1 if A < B, else 0.  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;
	for (i = 0; i < n_row; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][col][i];
		b = (constant[last]) ? factor[last] : stack[last][col][i];
		stack[prev][col][i] = (double)(a < b);
	}
}

void table_MAD (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: MAD 1 1 Median Absolute Deviation (L1 STD) of A.  */
{
	GMT_LONG i;
	double mad, med;

	if (constant[last]) {	/* Trivial case */
		memset ((void *)stack[last][col], 0, (size_t)(n_row * sizeof (double)));
		return;
	}

	qsort ((void *)stack[last][col], (size_t)n_row, sizeof (double), GMT_comp_double_asc);
	for (i = n_row; GMT_is_dnan (stack[last][col][i-1]) && i > 1; i--);
	if (i) {
		med = (i%2) ? stack[last][col][i/2] : 0.5 * (stack[last][col][(i-1)/2] + stack[last][col][i/2]);
		GMT_getmad (stack[last][col], i, med, &mad);
	}
	else
		mad = GMT_d_NaN;

	for (i = 0; i < n_row; i++) stack[last][col][i] = mad;
}

void table_MAX (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: MAX 2 1 Maximum of A and B.  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;
	for (i = 0; i < n_row; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][col][i];
		b = (constant[last]) ? factor[last] : stack[last][col][i];
		stack[prev][col][i] = MAX (a, b);
	}
}

void table_MEAN (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: MEAN 1 1 Mean value of A.  */
{
	GMT_LONG i, n_a = 0;
	double sum_a = 0.0;

	if (constant[last]) {	/* Trivial case */
		for (i = 0; i < n_row; i++) stack[last][col][i] = factor[last];
		return;
	}

	for (i = 0; i < n_row; i++) {
		if (GMT_is_dnan (stack[last][col][i])) continue;
		sum_a += stack[last][col][i];
		n_a++;
	}
	sum_a = (n_a) ? sum_a / n_a : GMT_d_NaN;
	for (i = 0; i < n_row; i++) stack[last][col][i] = sum_a;
}

void table_MED (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: MED 1 1 Median value of A.  */
{
	GMT_LONG i;
	double med;

	if (constant[last]) {	/* Trivial case */
		for (i = 0; i < n_row; i++) stack[last][col][i] = factor[last];
		return;
	}

	qsort ((void *)stack[last][col], (size_t)n_row, sizeof (double), GMT_comp_double_asc);
	for (i = n_row; GMT_is_dnan (stack[last][col][i-1]) && i > 1; i--);
	if (i)
		med = (i%2) ? stack[last][col][i/2] : 0.5 * (stack[last][col][(i-1)/2] + stack[last][col][i/2]);
	else
		med = GMT_d_NaN;

	for (i = 0; i < n_row; i++) stack[last][col][i] = med;
}

void table_MIN (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: MIN 2 1 Minimum of A and B.  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;
	for (i = 0; i < n_row; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][col][i];
		b = (constant[last]) ? factor[last] : stack[last][col][i];
		stack[prev][col][i] = MIN (a, b);
	}
}

void table_MOD (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: MOD 2 1 A mod B (remainder after floored division).  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;
	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, using MOD 0!\n", GMT_program);
	for (i = 0; i < n_row; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][col][i];
		b = (constant[last]) ? factor[last] : stack[last][col][i];
		stack[prev][col][i] = MOD (a, b);
	}
}

void table_MODE (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: MODE 1 1 Mode value (Least Median of Squares) of A.  */
{
	GMT_LONG i, GMT_mode_selection = 0, GMT_n_multiples = 0;
	double mode;

	if (constant[last]) {	/* Trivial case */
		for (i = 0; i < n_row; i++) stack[last][col][i] = factor[last];
		return;
	}

	qsort ((void *)stack[last][col], (size_t)n_row, sizeof (double), GMT_comp_double_asc);
	for (i = n_row; GMT_is_dnan (stack[last][col][i-1]) && i > 1; i--);
	if (i)
		GMT_mode (stack[last][col], i, i/2, 0, GMT_mode_selection, &GMT_n_multiples, &mode);
	else
		mode = GMT_d_NaN;

	for (i = 0; i < n_row; i++) stack[last][col][i] = mode;
	if (GMT_n_multiples > 0) fprintf (stderr, "%s: WARNING: %ld Multiple modes found\n", GMT_program, GMT_n_multiples);
}

void table_MUL (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: MUL 2 1 A * B.  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;
	if (gmtdefs.verbose && constant[prev] && factor[prev] == 0.0) fprintf (stderr, "%s: Warning, operand one == 0!\n", GMT_program);
	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand two == 0!\n", GMT_program);
	for (i = 0; i < n_row; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][col][i];
		b = (constant[last]) ? factor[last] : stack[last][col][i];
		stack[prev][col][i] = a * b;
	}
}

void table_NAN (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: NAN 2 1 NaN if A == B, else A.  */
{
	GMT_LONG i, prev;
	double a = 0.0, b = 0.0;

	prev = last - 1;
	if (constant[prev]) a = factor[prev];
	if (constant[last]) b = factor[last];
	for (i = 0; i < n_row; i++) {
		if (!constant[prev]) a = stack[prev][col][i];
		if (!constant[last]) b = stack[last][col][i];
		stack[prev][col][i] = ((a == b) ? GMT_d_NaN : a);
	}
}

void table_NEG (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: NEG 1 1 -A.  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand == 0!\n", GMT_program);
	if (constant[last]) a = -factor[last];
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : -stack[last][col][i];
}

void table_NEQ (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: NEQ 2 1 1 if A != B, else 0.  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;
	for (i = 0; i < n_row; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][col][i];
		b = (constant[last]) ? factor[last] : stack[last][col][i];
		stack[prev][col][i] = (double)(a != b);
	}
}

void table_NOT (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: NOT 1 1 NaN if A == NaN, 1 if A == 0, else 0.  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand == 0!\n", GMT_program);
	if (constant[last]) a = (fabs (factor[last]) > GMT_CONV_LIMIT) ? 0.0 : 1.0;
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : ((fabs (stack[last][col][i]) > GMT_CONV_LIMIT) ? 0.0 : 1.0);
}

void table_NRAND (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: NRAND 2 1 Normal, random values with mean A and std. deviation B.  */
{
	GMT_LONG i, prev;
	double a = 0.0, b = 0.0;

	prev = last - 1;
	if (constant[prev]) a = factor[prev];
	if (constant[last]) b = factor[last];
	for (i = 0; i < n_row; i++) {
		if (!constant[prev]) a = stack[prev][col][i];
		if (!constant[last]) b = stack[last][col][i];
		stack[prev][col][i] = a + b * GMT_nrand ();
	}
}

void table_OR (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: OR 2 1 NaN if A or B == NaN, else A.  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;
	for (i = 0; i < n_row; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][col][i];
		b = (constant[last]) ? factor[last] : stack[last][col][i];
		stack[prev][col][i] = (GMT_is_dnan (a) || GMT_is_dnan (b)) ? GMT_d_NaN : a;
	}
}

void table_PLM (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: PLM 3 1 Associated Legendre polynomial P(A) degree B order C.  */
{
	GMT_LONG i, prev, first, L, M;
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
	for (i = 0; i < n_row; i++) stack[first][col][i] = (constant[first]) ? a : GMT_plm (L, M, stack[first][col][i]);
}

void table_PLMg (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: PLMg 3 1 Normalized associated Legendre polynomial P(A) degree B order C (geophysical convention).  */
{
	GMT_LONG i, prev, first, L, M;
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
	for (i = 0; i < n_row; i++) stack[first][col][i] = (constant[first]) ? a : GMT_plm_bar (L, M, stack[first][col][i], FALSE);
}

void table_POP (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: POP 1 0 Delete top element from the stack.  */
{
	/* Dummy routine that does nothing but consume the top element of stack */
}

void table_POW (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: POW 2 1 A ^ B.  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;

	if (gmtdefs.verbose && constant[prev] && factor[prev] == 0.0) fprintf (stderr, "%s: Warning, operand one == 0!\n", GMT_program);
	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand two == 0!\n", GMT_program);
	for (i = 0; i < n_row; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][col][i];
		b = (constant[last]) ? factor[last] : stack[last][col][i];
		stack[prev][col][i] = pow (a, b);
	}
}

void table_PQUANT (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: PQUANT 2 1 The B'th Quantile (0-100%) of A.  */
{
	GMT_LONG i, prev;
	double p;

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
		p = GMT_d_NaN;
	}
	else {
		qsort ((void *)stack[prev][col], (size_t)n_row, sizeof (double), GMT_comp_double_asc);
		p = GMT_quantile (stack[prev][col], factor[last], n_row);
	}

	for (i = 0; i < n_row; i++) stack[prev][col][i] = p;
}

void table_PSI (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: PSI 1 1 Psi (or Digamma) of A.  */
{
	GMT_LONG i;
	double a = 0.0, x[2];

	x[1] = 0.0;	/* No imaginary part */
	if (constant[last]) {
		x[0] = factor[last];
		a = GMT_psi (x, (double *)NULL);
	}
	for (i = 0; i < n_row; i++) {
		if (!constant[last]) {
			x[0] = stack[last][col][i];
			a = GMT_psi (x, (double *)NULL);
		}
		stack[last][col][i] = a;
	}
}

void table_PV (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: PV 3 1 Legendre function Pv(A) of degree v = real(B) + imag(C).  */
{
	table_PVQV (info, stack, constant, factor, last, col, n_row, 0);
}

void table_QV (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: QV 3 1 Legendre function Qv(A) of degree v = real(B) + imag(C).  */
{
	table_PVQV (info, stack, constant, factor, last, col, n_row, 1);
}

void table_PVQV (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row, GMT_LONG kind)
{	/* kind: 0 = Pv, 1 = Qv */
	GMT_LONG i, prev, first, n;
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
	for (i = 0; i < n_row; i++) {
		if (calc){
			if (!constant[prev]) nu[0] = stack[prev][col][i];
			if (!constant[last]) nu[1] = stack[last][col][i];
			if (!constant[first])    x = stack[first][col][i];
			GMT_PvQv (x, nu, pq, &n);
			a = pq[kind];
		}
		stack[first][col][i] = a;
	}
}

void table_R2 (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: R2 2 1 R2 = A^2 + B^2.  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;
	if (gmtdefs.verbose && constant[prev] && factor[prev] == 0.0) fprintf (stderr, "%s: Warning, operand one == 0!\n", GMT_program);
	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand two == 0!\n", GMT_program);
	if (constant[prev]) factor[prev] *= factor[prev];
	if (constant[last]) factor[last] *= factor[last];
	for (i = 0; i < n_row; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][col][i] * stack[prev][col][i];
		b = (constant[last]) ? factor[last] : stack[last][col][i] * stack[last][col][i];
		stack[prev][col][i] = a + b;
	}
}

void table_R2D (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: R2D 1 1 Convert Radians to Degrees.  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = factor[last] * R2D;
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : stack[last][col][i] * R2D;
}

void table_RAND (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: RAND 2 1 Uniform random values between A and B.  */
{
	GMT_LONG i, prev;
	double a = 0.0, b = 0.0;

	prev = last - 1;
	if (constant[prev]) a = factor[prev];
	if (constant[last]) b = factor[last];
	for (i = 0; i < n_row; i++) {
		if (!constant[prev]) a = stack[prev][col][i];
		if (!constant[last]) b = stack[last][col][i];
		stack[prev][col][i] = a + GMT_rand () * (b - a);
	}
}

void table_RINT (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: RINT 1 1 rint (A) (nearest integer).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = rint (factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : rint (stack[last][col][i]);
}

void table_ROTT (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: ROTT 2 1 Rotate A by the (constant) shift B in the t-direction.  */
{
	GMT_LONG i, shift, prev;
	double *z;

	if (!constant[last]) {
		fprintf (stderr, "%s: T-shift must be a constant in ROTT\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	prev = last - 1;
	shift = irint (factor[last] / info->header.t_inc);
	if (constant[prev] || !shift) return;	/* Easy, constant or no shift */
	if (shift < 0) shift += n_row;		/* Same thing */

	z = (double *) GMT_memory (VNULL, (size_t)(n_row), sizeof (double), GMT_program);

	for (i = 0; i < n_row; i++) z[(i+shift)%n_row] = stack[prev][col][i];
	memcpy ((void *)stack[prev][col], (void *)z, (size_t)(n_row * sizeof (double)));
	GMT_free ((void *)z);
}

void table_SEC (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: SEC 1 1 sec (A) (A in radians).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = (1.0 / cos (factor[last]));
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : (1.0 / cos (stack[last][col][i]));
}

void table_SECD (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: SECD 1 1 sec (A) (A in degrees).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = (1.0 / cosd (factor[last]));
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : (1.0 / cosd (stack[last][col][i]));
}

void table_SIGN (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: SIGN 1 1 sign (+1 or -1) of A.  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand == 0!\n", GMT_program);
	if (constant[last]) a = copysign (1.0, factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : copysign (1.0, stack[last][col][i]);
}

void table_SIN (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: SIN 1 1 sin (A) (A in radians).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = sin (factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : sin (stack[last][col][i]);
}

void table_SINC (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: SINC 1 1 sinc (A) (sin (pi*A)/(pi*A)).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = GMT_sinc (factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : GMT_sinc (stack[last][col][i]);
}

void table_SIND (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: SIND 1 1 sin (A) (A in degrees).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = sind (factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : sind (stack[last][col][i]);
}

void table_SINH (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: SINH 1 1 sinh (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = sinh (factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : sinh (stack[last][col][i]);
}

void table_SKEW (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: SKEW 1 1 Skewness of A.  */
{
	GMT_LONG i, n = 0;
	double mean = 0.0, sum2 = 0.0, skew = 0.0, delta;

	if (constant[last]) {	/* Trivial case */
		for (i = 0; i < n_row; i++) stack[last][col][i] = GMT_d_NaN;
		return;
	}

	/* Use Welford (1962) algorithm to compute mean and corrected sum of squares */
	for (i = 0; i < n_row; i++) {
		if (GMT_is_dnan (stack[last][col][i])) continue;
		n++;
		delta = stack[last][col][i] - mean;
		mean += delta / n;
		sum2 += delta * (stack[last][col][i] - mean);
	}
	if (n > 1) {
		for (i = 0; i < n_row; i++) {
			if (GMT_is_dnan (stack[last][col][i])) continue;
			delta = stack[last][col][i] - mean;
			skew += pow (delta, 3.0);
		}
		sum2 /= (n - 1);
		skew /= n * pow (sum2, 1.5);
	}
	else
		skew = GMT_d_NaN;
	for (i = 0; i < n_row; i++) stack[last][col][i] = skew;
}

void table_SQRT (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: SQRT 1 1 sqrt (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && factor[last] < 0.0) fprintf (stderr, "%s: Warning, operand one < 0!\n", GMT_program);
	if (constant[last]) a = sqrt (factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : sqrt (stack[last][col][i]);
}

void table_SQR (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: SQR 1 1 A^2.  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = factor[last] * factor[last];
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : stack[last][col][i] * stack[last][col][i];
}

void table_STD (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: STD 1 1 Standard deviation of A.  */
{
	GMT_LONG i, n = 0;
	double mean = 0.0, sum2 = 0.0, delta;

	if (constant[last]) {	/* Trivial case */
		memset ((void *)stack[last][col], 0, (size_t)(n_row * sizeof (double)));
		return;
	}

	/* Use Welford (1962) algorithm to compute mean and corrected sum of squares */
	for (i = 0; i < n_row; i++) {
		if (GMT_is_dnan (stack[last][col][i])) continue;
		n++;
		delta = stack[last][col][i] - mean;
		mean += delta / n;
		sum2 += delta * (stack[last][col][i] - mean);
	}
	sum2 = (n > 1) ? sqrt (sum2 / (n - 1)) : GMT_d_NaN;
	for (i = 0; i < n_row; i++) stack[last][col][i] = sum2;
}

void table_STEP (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: STEP 1 1 Heaviside step function H(A).  */
{
	GMT_LONG i;
	double a;

	for (i = 0; i < n_row; i++) {
		a = (constant[last]) ? factor[last] : stack[last][col][i];
		if (a == 0.0)
			stack[last][col][i] = 0.5;
		else
			stack[last][col][i] = (a < 0.0) ? 0.0 : 1.0;
	}
}

void table_STEPT (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: STEPT 1 1 Heaviside step function H(t-A).  */
{
	GMT_LONG i;
	double a;

	for (i = 0; i < n_row; i++) {
		a = info->t_coordinates[i] - ((constant[last]) ? factor[last] : stack[last][col][i]);
		if (a == 0.0)
			stack[last][col][i] = 0.5;
		else
			stack[last][col][i] = (a < 0.0) ? 0.0 : 1.0;
	}
}

void table_SUB (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: SUB 2 1 A - B.  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;
	if (gmtdefs.verbose && constant[prev] && factor[prev] == 0.0) fprintf (stderr, "%s: Warning, operand one == 0!\n", GMT_program);
	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand two == 0!\n", GMT_program);
	for (i = 0; i < n_row; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][col][i];
		b = (constant[last]) ? factor[last] : stack[last][col][i];
		stack[prev][col][i] = a - b;
	}
}

void table_SUM (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: SUM 1 1 Cumulative sum of A.  */
{
	GMT_LONG i;
	double a = 0.0, sum = 0.0;

	if (constant[last]) a = factor[last];
	for (i = 0; i < n_row; i++) {
		if (!constant[last]) a = stack[last][col][i];
		if (!GMT_is_dnan (a)) sum += a;
		stack[last][col][i] = sum;
	}
}

void table_TAN (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: TAN 1 1 tan (A) (A in radians).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = tan (factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : tan (stack[last][col][i]);
}

void table_TAND (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: TAND 1 1 tan (A) (A in degrees).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = tand (factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : tand (stack[last][col][i]);
}

void table_TANH (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: TANH 1 1 tanh (A).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = tanh (factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : tanh (stack[last][col][i]);
}

void table_TN (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: TN 2 1 Chebyshev polynomial Tn(-1<A<+1) of degree B.  */
{
	GMT_LONG i, prev, n;
	double a;

	prev = last - 1;
	for (i = 0; i < n_row; i++) {
		n = irint ((constant[last]) ? factor[last] : stack[last][col][i]);
		a = (constant[prev]) ? factor[prev] : stack[prev][col][i];
		GMT_chebyshev (a, n, &stack[prev][col][i]);
	}
}

void table_TCRIT (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: TCRIT 2 1 Critical value for Student's t-distribution, with alpha = A and n = B.  */
{
	GMT_LONG i, prev;
	double a, b;

	prev = last - 1;
	if (gmtdefs.verbose && constant[prev] && factor[prev] == 0.0) fprintf (stderr, "%s: Warning, operand one == 0 for TCRIT!\n", GMT_program);
	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand two == 0 for TCRIT!\n", GMT_program);
	for (i = 0; i < n_row; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][col][i];
		b = (constant[last]) ? factor[last] : stack[last][col][i];
		stack[prev][col][i] = GMT_tcrit (a, b);
	}
}

void table_TDIST (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: TDIST 2 1 Student's t-distribution A(t,n), with t = A, and n = B.  */
{
	GMT_LONG i, b, prev;
	double a;

	prev = last - 1;
	if (gmtdefs.verbose && constant[prev] && factor[prev] == 0.0) fprintf (stderr, "%s: Warning, operand one == 0 for TDIST!\n", GMT_program);
	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand two == 0 for TDIST!\n", GMT_program);
	for (i = 0; i < n_row; i++) {
		a = (constant[prev]) ? factor[prev] : stack[prev][col][i];
		b = irint ((constant[last]) ? factor[last] : stack[last][col][i]);
		(void) GMT_student_t_a (a, (GMT_LONG)b, &stack[prev][col][i]);
	}
}

void table_UPPER (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: UPPER 1 1 The highest (maximum) value of A.  */
{
	GMT_LONG i;
	double high;

	if (constant[last]) {	/* Trivial case */
		for (i = 0; i < n_row; i++) stack[last][col][i] = factor[last];
		return;
	}

	for (i = 0, high = -DBL_MAX; i < n_row; i++) {
		if (GMT_is_dnan (stack[last][col][i])) continue;
		if (stack[last][col][i] > high) high = stack[last][col][i];
	}
	if (high == -DBL_MAX) high = GMT_d_NaN;	/* All rows were NaN */
	for (i = 0; i < n_row; i++) stack[last][col][i] = high;
}

void table_XOR (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: XOR 2 1 B if A == NaN, else A.  */
{
	GMT_LONG i, prev;
	double a = 0.0, b = 0.0;

	prev = last - 1;
	if (constant[prev]) a = factor[prev];
	if (constant[last]) b = factor[last];
	for (i = 0; i < n_row; i++) {
		if (!constant[prev]) a = stack[prev][col][i];
		if (!constant[last]) b = stack[last][col][i];
		stack[prev][col][i] = (GMT_is_dnan (a)) ? b : a;
	}
}

void table_Y0 (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: Y0 1 1 Bessel function of A (2nd kind, order 0).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand = 0 for Y0!\n", GMT_program);
	if (constant[last]) a = y0 (fabs (factor[last]));
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : y0 (fabs (stack[last][col][i]));
}

void table_Y1 (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: Y1 1 1 Bessel function of A (2nd kind, order 1).  */
{
	GMT_LONG i;
	double a = 0.0;

	if (gmtdefs.verbose && constant[last] && factor[last] == 0.0) fprintf (stderr, "%s: Warning, operand = 0 for Y1!\n", GMT_program);
	if (constant[last]) a = y1 (fabs (factor[last]));
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : y1 (fabs (stack[last][col][i]));
}

void table_YN (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: YN 2 1 Bessel function of A (2nd kind, order B).  */
{
	GMT_LONG i, prev;
	int order = 0;
	GMT_LONG simple = FALSE;
	double b = 0.0;

	prev = last - 1;
	if (gmtdefs.verbose && constant[last] && factor[last] < 0.0) fprintf (stderr, "%s: Warning, order < 0 for YN!\n", GMT_program);
	if (gmtdefs.verbose && constant[last] && fabs (rint(factor[last]) - factor[last]) > GMT_SMALL) fprintf (stderr, "%s: Warning, order not an integer for YN!\n", GMT_program);
	if (gmtdefs.verbose && constant[prev] && factor[prev] == 0.0) fprintf (stderr, "%s: Warning, argument = 0 for YN!\n", GMT_program);
	if (constant[last]) order = irint (fabs (factor[last]));
	if (constant[last] && constant[prev]) {
		b = yn (order, fabs (factor[prev]));
		simple = TRUE;
	}
	for (i = 0; i < n_row; i++) {
		if (simple)
			stack[prev][col][i] = b;
		else {
			if (!constant[last]) order = irint (fabs (stack[last][col][i]));
			stack[prev][col][i] = yn (order, fabs (stack[prev][col][i]));
		}
	}
}

void table_ZCRIT (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: ZCRIT 1 1 Critical value for the normal-distribution, with alpha = A.  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = GMT_zcrit (factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : GMT_zcrit (stack[last][col][i]);
}

void table_ZDIST (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: ZDIST 1 1 Cumulative normal-distribution C(x), with x = A.  */
{
	GMT_LONG i;
	double a = 0.0;

	if (constant[last]) a = GMT_zdist (factor[last]);
	for (i = 0; i < n_row; i++) stack[last][col][i] = (constant[last]) ? a : GMT_zdist (stack[last][col][i]);
}

void table_ROOTS (struct GMTMATH_INFO *info, double **stack[], GMT_LONG *constant, double *factor, GMT_LONG last, GMT_LONG col, GMT_LONG n_row)
/*OPERATOR: ROOTS 2 1 Treats col A as f(t) = 0 and returns its roots.  */
{
	GMT_LONG i, prev;
	double *roots;

	/* Treats the chosen column (at there is only one) as f(t) and solves for t that makes f(t) == 0.
	 * For now we only solve using a linear spline but in the future this should depend on the users
	 * choice of INTERPOLANT. */

	if (info->roots_found) return;	/* Already been here */
	if (!constant[last]) {
		fprintf (stderr, "%s: Argument to operator ROOTS must be a constant: the column number. Reset to 0\n", GMT_program);
		info->r_col = 0;
	}
	else
		info->r_col = irint (factor[last]);
	if (info->r_col < 0 || info->r_col >= info->header.n_col) {
		fprintf (stderr, "%s: Argument to operator ROOTS must be a column number 0 < col < %ld. Reset to 0\n", GMT_program, info->header.n_col);
		info->r_col = 0;
	}
	roots = (double *) GMT_memory (VNULL, (size_t)(n_row), sizeof (double), GMT_program);
	info->n_roots = 0;
	prev = last - 1;
	if (stack[prev][info->r_col][0] == 0.0) roots[info->n_roots++] = info->t_coordinates[0];
	for (i = 1; i < n_row; i++) {
		if (stack[prev][info->r_col][i] == 0.0) {
			roots[info->n_roots++] = info->t_coordinates[i];
			continue;
		}

		if ((stack[prev][info->r_col][i-1] * stack[prev][info->r_col][i]) < 0.0) {	/* Crossing 0 */
			roots[info->n_roots] = info->t_coordinates[i-1] - stack[prev][info->r_col][i-1] * (info->t_coordinates[i] - info->t_coordinates[i-1]) / (stack[prev][info->r_col][i] - stack[prev][info->r_col][i-1]);
			info->n_roots++;
		}
	}
	for (i = 0; i < info->n_roots; i++) stack[prev][info->r_col][i] = roots[i];
	GMT_free ((void *)roots);
	info->roots_found = TRUE;
}

/* ---------------------- end operator functions --------------------- */

void solve_LSQFIT (struct GMTMATH_INFO *info, double **stack[], GMT_LONG last, GMT_LONG n_col, GMT_LONG n_row, GMT_LONG skip[], char *file)
{
	/* Consider the current table the augmented matrix [A | b], making up the linear system Ax = b.
	 * We will set up the normal equations, solve for x, and output the solution before quitting.
	 * This function is special since it operates across columns and returns n_col scalars.
	 * We try to solve this positive definite & symmetric matrix with Cholsky methods; if that fails
	 * we do a full SVD decomposition and set small eigenvalues to zero, yielding an approximate solution.
	 */

	GMT_LONG i, j, k0, i2, j2, rhs, k, n, ier;
	double cond, *N, *B, *d, *x, *b, *z, *v, *lambda;
	FILE *fp;

	for (i = n = 0; i < n_col; i++) if (!skip[i]) n++;	/* Need to find how many active columns we have */
	if (n < 2) {
		fprintf (stderr, "%s: Error, LSQFIT requires at least 2 active columns!\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	rhs = n_col - 1;
	while (skip[rhs] && rhs > 0) rhs--;	/* Get last active col number as the rhs vector b */
	n--;					/* Account for b, the rhs vector, to get row & col dimensions of normal matrix N */

	N = (double *) GMT_memory (VNULL, (size_t)(n*n), sizeof (double), GMT_program);
	B = (double *) GMT_memory (VNULL, (size_t)n_row, sizeof (double), GMT_program);

	/* Do the row & col dot products, skipping inactive columns as we go along */
	for (j = j2 = 0; j < n; j2++) {	/* j2 is table column, j is row in N matrix */
		if (skip[j2]) continue;
		for (i = i2 = 0; i < n; i2++) {	/* i2 is table column, i is column in N matrix */
			if (skip[i2]) continue;
			k0 = j * n + i;
			N[k0] = 0.0;
			for (k = 0; k < n_row; k++) N[k0] += stack[last][j2][k] * stack[last][i2][k];
			i++;
		}
		B[j] = 0.0;
		for (k = 0; k < n_row; k++) B[j] += stack[last][j2][k] * stack[last][rhs][k];
		j++;
	}

	d = (double *) GMT_memory (VNULL, (size_t)n, sizeof(double), GMT_program);
	x = (double *) GMT_memory (VNULL, (size_t)n, sizeof(double), GMT_program);
	if ( (ier = GMT_chol_dcmp (N, d, &cond, n, n) ) != 0) {	/* Decomposition failed, use SVD method */
		GMT_LONG nrots;
		GMT_chol_recover (N, d, n, n, ier, TRUE);		/* Restore to former matrix N */
		/* Solve instead using GMT_jacobi */
		lambda = (double *) GMT_memory (VNULL, (size_t)n, sizeof(double), GMT_program);
		b = (double *) GMT_memory (VNULL, (size_t)n, sizeof(double), GMT_program);
		z = (double *) GMT_memory (VNULL, (size_t)n, sizeof(double), GMT_program);
		v = (double *) GMT_memory (VNULL, (size_t)n*n, sizeof(double), GMT_program);

		if (GMT_jacobi (N, &n, &n, lambda, v, b, z, &nrots)) {
			fprintf (stderr, "%s: Eigenvalue routine failed to converge in 50 sweeps.\n", GMT_program);
			fprintf (stderr, "%s: The reported L2 positions might be garbage.\n", GMT_program);
		}
		/* Solution x = v * lambda^-1 * v' * B */

		/* First do d = V' * B, so x = v * lambda^-1 * d */
		for (j = 0; j < n; j++) for (k = 0, d[j] = 0.0; k < n; k++) d[j] += v[k*n+j] * B[k];
		/* Then do d = lambda^-1 * d by setting small lambda's to zero */
		for (j = k = 0; j < n; j++) {
			if (lambda[j] < 1.0e7) {
				d[j] = 0.0;
				k++;
			}
			else
				d[j] /= lambda[j];
		}
		if (k) fprintf (stderr,"%s: %ld eigenvalues < 1.0e-7 set to zero to yield a stable solution\n", GMT_program, k);

		/* Finally do x = v * d */
		for (j = 0; j < n; j++) for (k = 0; k < n; k++) x[j] += v[j*n+k] * d[k];

		GMT_free ((void *)b);
		GMT_free ((void *)z);
		GMT_free ((void *)v);
	}
	else {	/* Decomposition worked, now solve system */
		GMT_chol_solv (N, x, B, n, n);
	}

	if (!file) {
		fp = GMT_stdout;
#ifdef SET_IO_MODE
		GMT_setmode (GMT_OUT);
#endif
	}
	else if ((fp = GMT_fopen (file, GMT_io.w_mode)) == NULL) {
		fprintf (stderr, "%s: Error creating file %s\n", GMT_program, file);
		exit (EXIT_FAILURE);
	}

	GMT_output (fp, n, x);
	if (fp != GMT_stdout) GMT_fclose (fp);

	GMT_free ((void *)x);
	GMT_free ((void *)d);
	GMT_free ((void *)N);
	GMT_free ((void *)B);
}

/* ---------------------- start convenience functions --------------------- */

GMT_LONG decode_argument (char *txt, double *value, struct GMT_HASH *H) {
	GMT_LONG expect, i, check = GMT_IS_NAN;
	GMT_LONG possible_number = FALSE;
	char copy[GMT_LONG_TEXT];
	char *mark;
	double tmp = 0.0;
	GMT_LONG get_operator (char *choice, struct GMT_HASH *H);

	/* Check if argument is operator */

	if ((i = get_operator (txt, H)) >= GMTMATH_ARG_IS_OPERATOR) return (i);

	/* Next look for symbols with special meaning */

	if (!(strcmp (txt, "STDIN"))) return GMTMATH_ARG_IS_FILE;	/* read from stdin */
	if (!(strcmp (txt, "PI") && strcmp (txt, "pi"))) return GMTMATH_ARG_IS_PI;
	if (!(strcmp (txt, "E") && strcmp (txt, "e"))) return GMTMATH_ARG_IS_E;
	if (!strcmp (txt, "EULER")) return GMTMATH_ARG_IS_EULER;
	if (!strcmp (txt, "TMIN")) return GMTMATH_ARG_IS_TMIN;
	if (!strcmp (txt, "TMAX")) return GMTMATH_ARG_IS_TMAX;
	if (!strcmp (txt, "TINC")) return GMTMATH_ARG_IS_TINC;
	if (!strcmp (txt, "N")) return GMTMATH_ARG_IS_N;
	if (!(strcmp (txt, "T") && strcmp (txt, "t"))) return GMTMATH_ARG_IS_T_MATRIX;
	if (!(strcmp (txt, "Tn") && strcmp (txt, "tn"))) return GMTMATH_ARG_IS_t_MATRIX;

	/* Preliminary test-conversion to a number */

	strcpy (copy, txt);
	if (!GMT_not_numeric (copy)) {	/* Only check if we are not sure this is NOT a number */
		expect = (strchr (copy, 'T')) ? GMT_IS_ABSTIME : GMT_IS_UNKNOWN;	/* Watch out for dateTclock-strings */
		check = GMT_scanf (copy, expect, &tmp);
		possible_number = TRUE;
	}

	/* Determine if argument is file. Remove possible question mark. */

	mark = strchr (copy, '?');
	if (mark) *mark = '\0';
	if (!GMT_access (copy, R_OK)) {	/* Yes it is */
		if (check != GMT_IS_NAN && possible_number) fprintf (stderr, "%s: WARNING: Your argument %s is both a file and a number.  File is selected\n", GMT_program, txt);
		return GMTMATH_ARG_IS_FILE;
	}

	if (check != GMT_IS_NAN) {	/* OK it is a number */
		*value = tmp;
		return GMTMATH_ARG_IS_NUMBER;
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

	op = GMT_hash_lookup (choice, H, GMTMATH_N_OPERATORS, GMTMATH_N_OPERATORS);

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

void *New_gmtmath_Ctrl () {	/* Allocate and initialize a new control structure */
	struct GMTMATH_CTRL *C;

	C = (struct GMTMATH_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GMTMATH_CTRL), "New_gmtmath_Ctrl");

	/* Initialize values whose defaults are not 0/FALSE/NULL */

	C->C.cols = (GMT_LONG *) GMT_memory (VNULL, (size_t)GMT_MAX_COLUMNS, sizeof (GMT_LONG), GMT_program);
	C->F.cols = (GMT_LONG *) GMT_memory (VNULL, (size_t)GMT_MAX_COLUMNS, sizeof (GMT_LONG), GMT_program);
	C->N.ncol = 2;

	return ((void *)C);
}

void Free_gmtmath_Ctrl (struct GMTMATH_CTRL *C) {	/* Deallocate control structure */
	if (C->A.file) free ((void *)C->A.file);
	GMT_free ((void *)C->C.cols);
	GMT_free ((void *)C->F.cols);
	if (C->T.file) free ((void *)C->T.file);
	GMT_free ((void *)C);
}

#include "gmtmath.h"
