/*-----------------------------------------------------------------
 *	$Id: x2sys_solve.c 10173 2014-01-01 09:52:34Z pwessel $
 *
 *      Copyright (c) 1999-2014 by P. Wessel
 *      See LICENSE.TXT file for copying and redistribution conditions.
 *
 *      This program is free software; you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation; version 2 or any later version.
 *
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *
 *      Contact info: gmt.soest.hawaii.edu
 *--------------------------------------------------------------------*/
/* x2sys_solve will read the crossover data base and determine the least-
 * squares correction coefficients for the specified field.  The correction
 * for each track is a sum of basis functions scaled by unknown coefficients
 * which we will solve for using least squares.  The normal equation matrix
 * is built directly and we solve the square linear system using Gauss-
 * Jordan elimination.
 *
 * Author:	Paul Wessel
 * Date:	18-SEPT-2008
 * Version:	1.0, based on the spirit of the old x_system code x_solve_dc_drift
 *		but completely rewritten from the ground up.
 *
 *
 */

/* #define DEBUGX */	/* Uncomment for testing */

#include "x2sys.h"

#define N_COE_PARS	12	/* Total number of items that might be known at each crossover */
#define COL_COE		0	/* The crossover value in whatever field we are studying */
#define COL_XX		1	/* lon or x at crossover */
#define COL_YY		1	/* lat or y at crossover */
#define COL_T1		3	/* Time along track one (in seconds) */
#define COL_T2		4	/* Time along track two (in seconds) */
#define COL_D1		5	/* Distance along track one */
#define COL_D2		6	/* Distance along track two */
#define COL_H1		7	/* Heading along track one (in degrees) */
#define COL_H2		8	/* Heading along track two (in degrees) */
#define COL_Z1		9	/* Observation value along track one */
#define COL_Z2		10	/* Observation value along track two */
#define COL_WW		11	/* Composite weight at crossover */

#define N_BASIS		7	/* Number of basis functions currently encoded */

/* The 6 different kinds of corrections coded so far */
#define F_IS_CONSTANT	1	/* Subtract a constant from each track */
#define F_IS_DRIFT_D	2	/* Subtract a trend with distance from each track */
#define F_IS_HEADING	3	/* Subtract a magnetic heading correction from each track */
#define F_IS_GRAV1930	4	/* Subtract a trend with latitude from each track */
#define F_IS_SCALE	5	/* Apply a scale to the observations for each track */
#define F_IS_DRIFT_T	6	/* Subtract a trend with time from each track */

/* The data matrix holds the COE information in these predefined column entries.
 * Not all columns are necessarily used and allocated:
 * Col 0:	COE value [Always used]
 * Col 1:	lon or x
 * Col 2:	lat or y
 * Col 3:	time_1
 * Col 4:	time_2
 * Col 5:	dist_1
 * Col 6:	dist_2
 * Col 7:	head_1
 * Col 8:	head_2
 * Col 9:	z_1
 * Col 10:	z_2
 * The array active_col[N_COE_PARS] will be TRUE for those columns actually used
 */

/* Available basis functions.  To add more basis functions:
 * 1. Add another one based on the examples below.
 * 2. Increase N_BASIS by one.
 * 3. If new data columns are needed you need to modify N_COE_PARS and comment above
 * The arguments to each basis functions are:
 * P:	  The 2-D array that holds all the parameters for all crossovers.  Each basis function
 *	  only uses those columns it requires.
 * which: Either 0 or 1 to indicate the first or second track in the crossover
 * col:	  The crossover number we are working on.
 */

double basis_constant (double **P, int which, int col) {	/* Basis function f for a constant c*f = c*1 : == 1 */
	return (1.0);	/* which, col are not used here */
}

double basis_tdrift (double **P, int which, int col) {	/* Basis function f for a linear drift rate in time c*f = c*t : t */
	return (P[COL_T1+which][col]);	/* When this is called the time columns have been corrected for t0 (start of track) */
}

double basis_ddrift (double **P, int which, int col) {	/* Basis function f for a linear drift rate in dist c*f = c*d : d */
	return (P[COL_D1+which][col]);
}

double basis_cosh (double **P, int which, int col) {	/* Basis function f for a dependence on cos(h)  c*f = c*cos(h) : cos(h) */
	return (cosd(P[COL_H1+which][col]));
}

double basis_cos2h (double **P, int which, int col) {	/* Basis function f for a dependence on cos(2*h)  c*f = c*cos(2*h) : cos(2*h) */
	return (cosd(2.0*P[COL_H1+which][col]));
}

double basis_sinh (double **P, int which, int col) {	/* Basis function f for a dependence on sin(h)  c*f = c*sin(h) : sin(h) */
	return (sind(P[COL_H1+which][col]));
}

double basis_sin2h (double **P, int which, int col) {	/* Basis function f for a dependence on sin(2*h)  c*f = c*sin(2*h) : sin(2*h) */
	return (sind(2.0*P[COL_H1+which][col]));
}

double basis_siny2 (double **P, int which, int col) {	/* Basis function f for a dependence on sin^2(y)  c*f = c*sin^2(y) : sin^2(y) */
	return (pow (sind(P[COL_YY][col]), 2.0));	/* which not used since y is common to both tracks */
}

double basis_z (double **P, int which, int col) {	/* Basis function f for a dependence on value c*f = c*z : z */
	return (P[COL_Z1+which][col]);
}

int main (int argc, char **argv)
{
	char *TAG = CNULL, *column = CNULL, **trk_list = NULL, *dbase = CNULL;
	char  trk[2][GMT_TEXT_LEN], t_txt[2][GMT_TEXT_LEN], z_txt[GMT_TEXT_LEN], w_txt[GMT_TEXT_LEN], line[BUFSIZ];
	GMT_LONG error = FALSE, grow_list = FALSE, normalize = FALSE, use_weights = FALSE, active_col[N_COE_PARS], reset_zero = FALSE;
	GMT_LONG write_mat_eq = FALSE, scale_range = FALSE;
	int *ID[2] = {NULL, NULL};
	int	it, n_iterations = 1;
	GMT_LONG n_par = 0, n, m, t, n_tracks = 0, mode = 0, n_active;
	GMT_LONG i, p, j, k, r, s, off, row, n_COE = 0, n_alloc = GMT_CHUNK, n_alloc_t = GMT_CHUNK, ierror;
	double *N = NULL, *a = NULL, *b = NULL, *data[N_COE_PARS], sgn, zero_test = 1.0e-08, old_mean, new_mean, sw2;
	double old_stdev, new_stdev, e_k, Sw, Sx, Sxx, range = 0.0, *start = NULL;
	double sc_range = 1;
	struct X2SYS_INFO *S = NULL;
	struct X2SYS_BIX B;
	FILE *fp = NULL;
	PFD basis[N_BASIS];

	argc = (int)GMT_begin (argc, argv);
	
	for (i = strlen(argv[0]); i >= 0 && argv[0][i] != '/'; i--);
	X2SYS_program = &argv[0][i+1];	/* Name without full path */
	memset ((void *)active_col, FALSE, N_COE_PARS * sizeof (GMT_LONG));

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				case 'V':
				case 'b':
				case '\0':
					error += GMT_parse_common_options (argv[i], NULL, NULL, NULL, NULL);
					break;
				case 'C':	/* Needed to report correctly */
					column = &argv[i][2];
					break;
				case 'E':	/* Which model to fit */
					switch (argv[i][2]) {
						case 'c':
							mode = F_IS_CONSTANT;
							break;
						case 'd':
							mode = F_IS_DRIFT_D;
							break;
						case 'g':
							mode = F_IS_GRAV1930;
							break;
						case 'h':
							mode = F_IS_HEADING;
							break;
						case 's':
							mode = F_IS_SCALE;
							break;
						case 't':
							mode = F_IS_DRIFT_T;
							break;
					}
					break;
				case 'M':
					write_mat_eq = TRUE;
					break;
				case 'S':
					scale_range = TRUE;
					sc_range = atof(&argv[i][2]);
					break;
				case 'T':
					TAG = &argv[i][2];
					break;
				case 'W':
					use_weights = TRUE;
					break;
				case 'Z':
					reset_zero = TRUE;
					break;
				default:
					error = TRUE;
					fprintf (stderr, "%s: Unrecognized option -%c\n", GMT_program, argv[i][1]);
					break;
			}
		}
		else
			dbase = argv[i];
	}

	if (argc == 1 || error || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "x2sys_solve %s - Determine least-squares systematic correction from crossovers\n\n", X2SYS_VERSION);
		fprintf (stderr, "usage: x2sys_solve -C<column> -E<flag> -T<TAG> [<coedata>] [-V] [-W] [Z] [%s]\n\n", GMT_bi_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\t-C specifies the column name to process (e.g., faa, mag)\n");
		fprintf (stderr, "\t-E Equation to fit: specify <flag> as c (constant), d (drift over distance),\n");
		fprintf (stderr, "\t     g (latitude), h (heading), s (scale with data), or t (drift over time) [c].\n");
		fprintf (stderr, "\t-T <TAG> is the x2sys tag for the data set.\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t<coedata> is the ASCII data output file from x2sys_list [or we read stdin]\n");
		GMT_explain_option ('V');
		fprintf (stderr, "\t-W weights are present in last column for weighted fit [no weights]\n");
		fprintf (stderr, "\t-Z Use the earliest time and shortest distance as local origin [use 0].\n");
		GMT_explain_option ('i');
		exit (EXIT_FAILURE);
	}

	if (mode < 0) {
		fprintf (stderr, "%s: ERROR -E: Choose among c, d, g, h, s and t\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	
	/* Initialize system via the tag */
	
	x2sys_err_fail (x2sys_set_system (TAG, &S, &B, &GMT_io), TAG);

	/* Verify that the chosen column is known to the system */
	
	if (column) x2sys_err_fail (x2sys_pick_fields (column, S), "-C");
	if (S->n_out_columns != 1) {
		fprintf (stderr, "%s: ERROR: -C must specify a single column name\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	active_col[COL_COE] = TRUE;	/* Always used */
	switch (mode) {	/* Set up pointers to basis functions and assign constants */
		case F_IS_CONSTANT:
			n_par = 1;
			basis[0] = basis_constant;
			break;
		case F_IS_DRIFT_T:
			active_col[COL_T1] = active_col[COL_T2] = TRUE;
			n_par = 2;
			basis[0] = basis_constant;
			basis[1] = basis_tdrift;
			break;
		case F_IS_DRIFT_D:
			active_col[COL_D1] = active_col[COL_D2] = TRUE;
			n_par = 2;
			basis[0] = basis_constant;
			basis[1] = basis_ddrift;
			break;
		case F_IS_GRAV1930:
			active_col[COL_YY] = TRUE;
			n_par = 2;
			basis[0] = basis_constant;
			basis[1] = basis_siny2;
			break;
		case F_IS_HEADING:
			active_col[COL_H1] = active_col[COL_H2] = TRUE;
			n_par = 5;
			basis[0] = basis_constant;
			basis[1] = basis_cosh;
			basis[2] = basis_cos2h;
			basis[3] = basis_sinh;
			basis[4] = basis_sin2h;
			break;
		case F_IS_SCALE:
			active_col[COL_Z1] = active_col[COL_Z2] = TRUE;
			n_par = 1;
			basis[0] = basis_z;
			break;
	}
	
	/* Allocate memory for COE data */
	
	for (i = n_active = 0; i < N_COE_PARS; i++) if (active_col[i]) {
		data[i] = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (double), GMT_program);
		n_active++;
	}
	data[COL_WW] = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (double), GMT_program);
	for (i = 0; i < 2; i++) ID[i] = (int *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (int), GMT_program);
	
	/* Read the crossover info */

	fp = GMT_stdin;
	if (dbase && (fp = GMT_fopen (dbase, GMT_io.r_mode)) == NULL) {
		fprintf (stderr, "%s: ERROR: Cannot open file %s\n", GMT_program, dbase);
		exit (EXIT_FAILURE);	
	}
	
	if (n_tracks == 0)	{	/* Create track list on the go */
		grow_list = TRUE;
		trk_list = (char **) GMT_memory (VNULL, (size_t)n_alloc_t, sizeof (char *), GMT_program);
	}
	
	n_COE = 0;
	if (GMT_io.binary[GMT_IN]) {	/* Binary input */
		/* Here, first two cols have track IDs and we do not write track names */
		int min_ID, max_ID;
		GMT_LONG n_fields, n_tracks2, n_expected_fields;
		double *in;
		char *check;
		
		min_ID = INT_MAX;	max_ID = -INT_MAX;
		n_expected_fields = n_active + use_weights;
		while ((n_fields = GMT_input (fp, &n_expected_fields, &in)) >= 0 && !(GMT_io.status & GMT_IO_EOF)) {	/* Not yet EOF */
			for (i = 0; i < 2; i++) {	/* Get IDs and keept track of min/max values */
				ID[i][n_COE] = irint (in[i]);
				if (ID[i][n_COE] < min_ID) min_ID = ID[i][n_COE];
				if (ID[i][n_COE] > max_ID) max_ID = ID[i][n_COE];
			}
			switch (mode) {	/* Handle input differently depending on what is expected */
				case F_IS_CONSTANT:
					data[COL_COE][n_COE] = in[2];
					break;
				case F_IS_DRIFT_T:
					data[COL_T1][n_COE] = in[2];
					data[COL_T2][n_COE] = in[3];
					data[COL_COE][n_COE] = in[4];
					break;
				case F_IS_DRIFT_D:
					data[COL_D1][n_COE] = in[2];
					data[COL_D2][n_COE] = in[3];
					data[COL_COE][n_COE] = in[4];
					break;
				case F_IS_GRAV1930:
					data[COL_YY][n_COE] = in[2];
					data[COL_COE][n_COE] = in[3];
					break;
				case F_IS_HEADING:
					data[COL_H1][n_COE] = in[2];
					data[COL_H2][n_COE] = in[3];
					data[COL_COE][n_COE] = in[4];
					break;
				case F_IS_SCALE:
					data[COL_Z1][n_COE] = in[2];
					data[COL_Z2][n_COE] = in[3];
					data[COL_COE][n_COE] = in[2] - in[3];
					break;
			}
			data[COL_WW][n_COE] = (use_weights) ? in[n_active] : 1.0;	/* Weight */
			if (GMT_is_dnan (data[COL_COE][n_COE])) {
				fprintf (stderr, "%s: Warning: COE == NaN skipped during reading\n", GMT_program);
				continue;
			}
			n_COE++;
			if (n_COE == n_alloc) {
				n_alloc <<= 1;
				for (i = 0; i < N_COE_PARS; i++) if (active_col[i]) data[i] = (double *) GMT_memory ((void *)data[i], (size_t)n_alloc, sizeof (double), GMT_program);
				data[COL_WW] = (double *) GMT_memory ((void *)data[COL_WW], (size_t)n_alloc, sizeof (double), GMT_program);
			}
		}
		/* Check that IDs are all contained within 0 <= ID < n_tracks and that there are no gaps */
		n_tracks2 = max_ID - min_ID + 1;
		if (n_tracks && n_tracks2 != n_tracks) {
			fprintf (stderr, "%s: ERROR: The ID numbers in the binary file %s are not compatible with the <trklist> length\n", GMT_program, dbase);
			error = TRUE;	
		}
		else {	/* Either no tracks read before or the two numbers did match properly */
			/* Look for ID gaps */
			n_tracks = n_tracks2;
			check = (char *) GMT_memory (VNULL, (size_t)n_tracks, sizeof (char), GMT_program);
			for (k = 0; k < n_COE; k++) for (i = 0; i < 2; i++) check[ID[i][k]] = TRUE;
			for (k = 0; k < n_tracks && check[k]; k++);
			GMT_free ((void *)check);
			if (k < n_tracks) {
				fprintf (stderr, "%s: ERROR: The ID numbers in the binary file %s to not completely cover the range 0 <= ID < n_tracks!\n", GMT_program, dbase);
				error = TRUE;
			}
		}
		if (error) {	/* Delayed the cleanup until here */
			for (i = 0; i < N_COE_PARS; i++) if (active_col[i]) GMT_free ((void *)data[i]);
			GMT_free ((void *)data[COL_WW]);
			for (i = 0; i < 2; i++) GMT_free ((void *)ID[i]);
			exit (EXIT_FAILURE);
		}
	}
	else {	/* Ascii input with track names */
		char file_TAG[GMT_TEXT_LEN], file_column[GMT_TEXT_LEN], *not_used = NULL;
		not_used = GMT_fgets (line, BUFSIZ, fp);	/* Read first line with TAG and column */
		sscanf (&line[7], "%s %s", file_TAG, file_column);
		if (strcmp (TAG, file_TAG) && strcmp (column, file_column)) {
			fprintf (stderr, "%s: ERROR: The TAG and column info in the ASCII file %s are not compatible with the -C -T options\n", GMT_program, dbase);
			exit (EXIT_FAILURE);	
		}
		while (GMT_fgets (line, BUFSIZ, fp)) {    /* Not yet EOF */
			GMT_chop (line);					/* Rid the world of CR/LF */
			if (line[0] == '#') continue;	/* Skip other comments */
			switch (mode) {	/* Handle input differently depending on what is expected */
				case F_IS_CONSTANT:
					sscanf (line, "%s %s %s %s", trk[0], trk[1], z_txt, w_txt);
					if (GMT_scanf (z_txt, GMT_IS_FLOAT, &data[COL_COE][n_COE]) == GMT_IS_NAN) 
						data[COL_COE][n_COE] = GMT_d_NaN;
					break;
				case F_IS_DRIFT_T:
					sscanf (line, "%s %s %s %s %s %s", trk[0], trk[1], t_txt[0], t_txt[1], z_txt, w_txt);
					if (GMT_scanf (t_txt[0], GMT_IS_ABSTIME, &data[COL_T1][n_COE]) == GMT_IS_NAN) 
						data[COL_T1][n_COE] = GMT_d_NaN;
					if (GMT_scanf (t_txt[1], GMT_IS_ABSTIME, &data[COL_T2][n_COE]) == GMT_IS_NAN) 
						data[COL_T2][n_COE] = GMT_d_NaN;
					if (GMT_scanf (z_txt, GMT_IS_FLOAT, &data[COL_COE][n_COE]) == GMT_IS_NAN) 
						data[COL_COE][n_COE] = GMT_d_NaN;
					break;
				case F_IS_DRIFT_D:
					sscanf (line, "%s %s %s %s %s %s", trk[0], trk[1], t_txt[0], t_txt[1], z_txt, w_txt);
					if (GMT_scanf (t_txt[0], GMT_IS_FLOAT, &data[COL_D1][n_COE]) == GMT_IS_NAN) 
						data[COL_D1][n_COE] = GMT_d_NaN;
					if (GMT_scanf (t_txt[1], GMT_IS_FLOAT, &data[COL_D2][n_COE]) == GMT_IS_NAN) 
						data[COL_D2][n_COE] = GMT_d_NaN;
					if (GMT_scanf (z_txt, GMT_IS_FLOAT, &data[COL_COE][n_COE]) == GMT_IS_NAN) 
						data[COL_COE][n_COE] = GMT_d_NaN;
					break;
				case F_IS_GRAV1930:
					sscanf (line, "%s %s %s %s %s", trk[0], trk[1], t_txt[0], z_txt, w_txt);
					if (GMT_scanf (t_txt[0], GMT_IS_LAT, &data[COL_YY][n_COE]) == GMT_IS_NAN) 
						data[COL_YY][n_COE] = GMT_d_NaN;
					if (GMT_scanf (z_txt, GMT_IS_FLOAT, &data[COL_COE][n_COE]) == GMT_IS_NAN) 
						data[COL_COE][n_COE] = GMT_d_NaN;
					break;
				case F_IS_HEADING:
					sscanf (line, "%s %s %s %s %s %s", trk[0], trk[1], t_txt[0], t_txt[1], z_txt, w_txt);
					if (GMT_scanf (t_txt[0], GMT_IS_FLOAT, &data[COL_H1][n_COE]) == GMT_IS_NAN) 
						data[COL_H1][n_COE] = GMT_d_NaN;
					if (GMT_scanf (t_txt[1], GMT_IS_FLOAT, &data[COL_H2][n_COE]) == GMT_IS_NAN) 
						data[COL_H2][n_COE] = GMT_d_NaN;
					if (GMT_scanf (z_txt, GMT_IS_FLOAT, &data[COL_COE][n_COE]) == GMT_IS_NAN) 
						data[COL_COE][n_COE] = GMT_d_NaN;
					break;
				case F_IS_SCALE:
					sscanf (line, "%s %s %s %s %s", trk[0], trk[1], t_txt[0], t_txt[1], w_txt);
					if (GMT_scanf (t_txt[0], GMT_IS_FLOAT, &data[COL_Z1][n_COE]) == GMT_IS_NAN) 
						data[COL_Z1][n_COE] = GMT_d_NaN;
					if (GMT_scanf (t_txt[1], GMT_IS_FLOAT, &data[COL_Z2][n_COE]) == GMT_IS_NAN) 
						data[COL_Z2][n_COE] = GMT_d_NaN;
					break;
			}
			if (use_weights) {
				if (GMT_scanf (w_txt, GMT_IS_FLOAT, &data[COL_WW][n_COE]) == GMT_IS_NAN) data[COL_WW][n_COE] = GMT_d_NaN;
			}
			else
				data[COL_WW][n_COE] = 1.0;
			if (GMT_is_dnan (data[COL_COE][n_COE])) {
				fprintf (stderr, "%s: Warning: COE == NaN skipped during reading\n", GMT_program);
				continue;
			}
			
			for (i = 0; i < 2; i++) {	/* Look up track IDs */
				ID[i][n_COE] = x2sys_find_track (trk[i], trk_list, (int)n_tracks);	/* Return track id # for this leg */
				if (ID[i][n_COE] == -1) {	/* Leg not in the data base yet */
					if (grow_list) {	/* Add it */
						trk_list[n_tracks] = strdup (trk[i]);
						ID[i][n_COE] = (int)n_tracks++;
						if (n_tracks == n_alloc_t) {
							n_alloc_t <<= 1;
							trk_list = (char **) GMT_memory ((void *)trk_list, (size_t)n_alloc_t, sizeof (char *), GMT_program);
						}
					}
				}
			}
		
			n_COE++;
			if (n_COE == n_alloc) {
				n_alloc <<= 1;
				for (i = 0; i < N_COE_PARS; i++) 
					if (active_col[i]) 
						data[i] = (double *) GMT_memory ((void *)data[i], (size_t)n_alloc, sizeof (double), GMT_program);
				data[COL_WW] = (double *) GMT_memory ((void *)data[COL_WW], (size_t)n_alloc, sizeof (double), GMT_program);
				for (i = 0; i < 2; i++) 
					ID[i] = (int *) GMT_memory ((void *)ID[i], (size_t)n_alloc, sizeof (int), GMT_program);
			}
		}
	}
	if (fp != GMT_stdin) GMT_fclose(fp);
	if (gmtdefs.verbose) fprintf (stderr, "%s: Found %ld COE records\n", GMT_program, n_COE);
	for (i = 0; i < N_COE_PARS; i++) 
		if (active_col[i]) data[i] = (double *) GMT_memory ((void *)data[i], (size_t)n_COE, sizeof (double), GMT_program);
	data[COL_WW] = (double *) GMT_memory ((void *)data[COL_WW], (size_t)n_COE, sizeof (double), GMT_program);
	
	normalize = (mode == F_IS_DRIFT_T || mode == F_IS_DRIFT_D);
	if (normalize) {	/* For numerical stability, normalize distances or times to fall in 0-1 range */
		double min_extent, max_extent;
		start = (double *) GMT_memory (NULL, (size_t)n_tracks, sizeof (double), GMT_program);
		j = (mode == F_IS_DRIFT_T) ? COL_T1 : COL_D1;	/* Which variable we are working on */
		
		if (reset_zero) {	/* Must determine smallest t or d for all tracks */
			for (k = 0; k < n_tracks; k++) start[k] = DBL_MAX;
			for (k = 0; k < n_COE; k++) {
				for (i = 0; i < 2; i++) {
					if (data[j+i][k] < start[ID[i][k]]) start[ID[i][k]] = data[j+i][k];
				}
			}
		}
		/* Now we know the start time/distance [or all 0].  Remove it to reduce mangitudes of values */
		min_extent = DBL_MAX;	max_extent = -DBL_MAX;
		for (k = 0; k < n_COE; k++) {
			for (i = 0; i < 2; i++) {
				if (reset_zero) data[j+i][k] -= start[ID[i][k]];
				if (data[j+i][k] < min_extent) min_extent = data[j+i][k];
				if (data[j+i][k] > max_extent) max_extent = data[j+i][k];
			}
		}
		range = max_extent - min_extent;
		if (scale_range) range /= sc_range;
		for (k = 0; k < n_COE; k++) for (i = 0; i < 2; i++) data[j+i][k] /= range;
	}
	
	/* Estimate old weighted mean and std.dev */
	
	for (k = 0, Sw = Sx = Sxx = 0.0; k < n_COE; k++) {	/* For each crossover */
		Sw += data[COL_WW][k];
		Sx += (data[COL_WW][k] * data[COL_COE][k]);
		Sxx += (data[COL_WW][k] * data[COL_COE][k] * data[COL_COE][k]);
	}
	old_mean = Sx / Sw;
	old_stdev = sqrt ((n_COE * Sxx - Sx * Sx) / (Sw*Sw*(n_COE - 1.0)/n_COE));
	
	/* Set up matrix and column vectors */
	
	n = n_tracks * n_par;	/* Total number of unknowns */
	m = (mode == F_IS_SCALE) ? n : n + 1;	/* Need extra row/column to handle Lagrange's multiplier for unknown absolute level */
	N = (double *) GMT_memory (VNULL, (size_t)(m*m), sizeof (double), GMT_program);
	b = (double *) GMT_memory (VNULL, (size_t)m, sizeof (double), GMT_program);

	/* Build A'A and A'b ==> N*x = b normal equation matrices directly since A may be too big to do A'A by multiplication.
	 * For all corrections that involve a constant shift we must add the constraint that such shifts sum to zero; this
	 * is implemented by adding an extra row/column with the appropriate 0s and 1s and a Lagrange multiplier. */

	for (it = 0; it < n_iterations; it++) {
		if (it > 0) {
			memset((void *)N, 0, m*m*sizeof(double));
			memset((void *)b, 0, m*sizeof(double));
		}

	for (p = row = 0; p < n_tracks; p++) {	/* For each track */
		for (s = 0; s < n_par; s++, row++) {	/* For each track's unknown parameter  */
			for (k = 0; k < n_COE; k++) {	/* For each crossover */
				i = ID[0][k];	/* Get track # 1 ID */
				j = ID[1][k];	/* Get track # 2 ID */
				if (i == p) {
					sgn = +1.0;	t = 0;
				} else if (p == j) {
					sgn = -1.0;	t = 1;
				} else continue;
				sw2 = sgn * data[COL_WW][k] * data[COL_WW][k];
				for (r = 0, off = m * row; r < n_par; r++) {	/* For each track's parameter in f(p)  */
					N[off+i*n_par+r] += sw2 * (basis[r](data,0,k) * basis[s](data,t,k));
					N[off+j*n_par+r] -= sw2 * (basis[r](data,1,k) * basis[s](data,t,k));
				}
				b[row] += sw2 * (data[COL_COE][k] * basis[s](data,t,k));
			}
			if (mode != F_IS_SCALE && s == 0) N[m*row+m-1] = 1.0;	/* Augmented column entry for Lagrange multiplier */
		}
	}

	if (mode != F_IS_SCALE) {	/* Augmented row for Lagrange multiplier for constants */
		for (i = 0, off = m*n; i < n; i += n_par) N[off+i] = 1.0;
	}
	
#ifdef DEBUGX	
	fprintf (stderr, "Matrix equation N * a = b: (N = %ld)\n", m);
	for (i = 0; i < m; i++) {
		for (j = 0; j < m; j++) fprintf (stderr, "%8.2f\t", N[i*m+j]);
		fprintf (stderr, "\t%8.2f\n", b[i]);
	}
#endif

	if (write_mat_eq) {
		char *Na, *B;
		FILE *fp;
		Na = "Na.dat";		B = "b.dat";
		fp = fopen(Na, "w");
		for (i = 0; i < m; i++) {
			for (j = 0; j < m; j++) fprintf (fp, "%f%s", N[i*m+j], gmtdefs.field_delimiter);
			fprintf (fp, "\n");
		}
		fclose(fp);

		fp = fopen(B, "w");
		for (i = 0; i < m; i++)
			fprintf (fp, "%f\n", b[i]);
		fclose(fp);
	}

	/* Get LS solution */

	GMT_gauss (N, b, m, m, zero_test, &ierror, 1);
	if (it == n_iterations-1) GMT_free ((void *)N);
	a = b;	/* Convenience since the solution is called a in the notes */

	/* Estimate new st.dev. */
	
	for (k = 0, Sw = Sx = Sxx = 0.0; k < n_COE; k++) {	/* For each crossover */
		i = ID[0][k];	/* Get track # 1 ID */
		j = ID[1][k];	/* Get track # 2 ID */
		e_k = data[COL_COE][k];
		for (r = 0; r < n_par; r++) {	/* Correct crossover for track adjustments  */
			e_k += a[j*n_par+r]*basis[r](data,1,k) - a[i*n_par+r]*basis[r](data,0,k);
		}
		Sw += data[COL_WW][k];
		Sx += (data[COL_WW][k] * e_k);
		Sxx += (data[COL_WW][k] * e_k * e_k);
#ifdef DEBUGX	
		fprintf (stderr, "COE # %ld: Was %g Is %g\n", k, data[COL_COE][k], e_k);
#endif
		if (it > 0) {
			data[COL_COE][k] = e_k;
		}
	}
	new_mean = Sx / Sw;
	new_stdev = sqrt ((n_COE * Sxx - Sx * Sx) / (Sw*Sw*(n_COE - 1.0)/n_COE));
	
		if (gmtdefs.verbose)
			fprintf (stderr, "%s: Iteration %d:\tOld mean and st.dev.: %g %g New mean and st.dev.: %g %g\n", 
				GMT_program, it+1, old_mean, old_stdev, new_mean, new_stdev);

		old_mean  = new_mean;
		old_stdev = new_stdev;

		if (it > 0) {
			for (p = 0; p < n_tracks; p++)		/* For each track */ 
				for (r = 0; r < n_par; r++) 	/* For each track's parameter in f(p)  */
					a[p*n_par + r] += a[p*n_par + r];
		}
	}	/* end iteration loop */
	
	/* Write correction table */
	
	for (p = 0; p < n_tracks; p++) {
		if (normalize) a[p*n_par+1] /= range;	/* Unnormalize slopes */
		(GMT_io.binary[GMT_IN]) ? printf ("%ld", p) : printf ("%s", trk_list[p]);
		printf ("\t%s", column);
		switch (mode) {	/* Set up pointers to basis functions and assign constants */
			case F_IS_CONSTANT:
				printf ("\t%g\n", a[p]);
				break;
			case F_IS_DRIFT_T:
				if (reset_zero)
					printf ("\t%g\t%g*((time-T))\n", a[p*n_par], a[p*n_par+1]);
				else
					printf ("\t%g\t%g*((time-%g))\n", a[p*n_par], a[p*n_par+1], start[p]);
				break;
			case F_IS_DRIFT_D:
				if (reset_zero)
					printf ("\t%g\t%g*((dist-%g))\n", a[p*n_par], a[p*n_par+1], start[p]);
				else
					printf ("\t%g\t%g*((dist))\n", a[p*n_par], a[p*n_par+1]);
				break;
			case F_IS_GRAV1930:
				printf ("\t%g\t%g*sin((lat))^2\n", a[p*n_par], a[p*n_par+1]);
				break;
			case F_IS_HEADING:
				printf ("\t%g\t%g*cos((azim))\t%g*cos(2*(azim))\t%g*sin((azim))\t%g*sin(2*(azim))\n", a[p*n_par], a[p*n_par+1], a[p*n_par+2], a[p*n_par+3], a[p*n_par+4]);
				break;
			case F_IS_SCALE:
				printf ("\t%g*((%s))\n", 1.0 - a[p], column);
				break;
		}
	}
	
	/* Free up memory */
	
	for (i = 0; i < N_COE_PARS; i++) if (active_col[i]) GMT_free ((void *)data[i]);
	GMT_free ((void *)data[COL_WW]);
	for (i = 0; i < 2; i++) GMT_free ((void *)ID[i]);
	GMT_free ((void *)b);
	if (normalize) GMT_free ((void *)start);
	x2sys_free_list (trk_list, (int)n_tracks);

	x2sys_end (S);
	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}
