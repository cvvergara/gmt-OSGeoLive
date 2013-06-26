/*--------------------------------------------------------------------
 *    $Id: filter1d.c 9923 2012-12-18 20:45:53Z pwessel $
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
 * filter1d will read N columns of data from file/GMT_stdin and return
 * filtered output at user-selected positions.  The time variable
 * can be in any specified column of the input.  Several filters
 * are available, with robustness as an option.
 *
 * Filters:
 * Convolutions: Boxcar, Cosine Arch, Gaussian, Median, or Mode.
 * Geospatial:   Median, Mode, Extreme values.
 * Robust:	Option replaces outliers with medians in filter.
 * Output:	At input times, or from t_start to t_stop by t_int.
 * Lack:	Option checks for data gaps in input series.
 * Symmetry:	Option checks for asymmetry in filter window.
 * Quality:	Option checks for low mean weight in window.
 *
 * Author:	Walter H. F. Smith
 * Date:	13 January, 1989
 * Version:	4.1.2 Free of global variables.
 */

#include "gmt.h"

struct FILTER1D_CTRL {
	struct D {	/* -D<inc> */
		GMT_LONG active;
		double inc;
	} D;
	struct E {	/* -E */
		GMT_LONG active;
	} E;
	struct F {	/* -F<type><width>[<mode>] */
		GMT_LONG active;
		char filter;	/* Character codes for the filter */
		double width;
		GMT_LONG mode;
		char *file;	/* Character codes for the filter */
	} F;
	struct I {	/* -I<ignoreval> */
		GMT_LONG active;	/* TRUE when input values that equal value should be discarded */
		double value;
	} I;
	struct L {	/* -L<lackwidth> */
		GMT_LONG active;
		double value;
	} L;
	struct N {	/* -N<n_cols>/<t_col> */
		GMT_LONG active;
		GMT_LONG col, ncols;
	} N;
	struct Q {	/* -Q<factor> */
		GMT_LONG active;
		double value;
	} Q;
	struct S {	/* -S<symmetry> */
		GMT_LONG active;
		double value;
	} S;
	struct T {	/* -T[<tmin/tmax/t_inc>] */
		GMT_LONG active;
		double min, max, inc;
	} T;
};

#define FILTER1D_BOXCAR 0
#define FILTER1D_COS_ARCH 1
#define FILTER1D_GAUSSIAN 2
#define FILTER1D_MEDIAN 3
#define FILTER1D_MODE 4
#define FILTER1D_LOWER_ALL 5
#define FILTER1D_LOWER_POS 6
#define FILTER1D_UPPER_ALL 7
#define FILTER1D_UPPER_NEG 8
#define FILTER1D_N_FILTERS 9
#define FILTER1D_CONVOLVE 2		/* If filter_type > FILTER1D_CONVOLVE then a FILTER1D_MEDIAN, FILTER1D_MODE, or EXTREME filter is selected  */

struct FILTER1D_INFO {	/* Control structure for all aspects of the filter setup */
	GMT_LONG	use_ends;	/* True to start/stop at ends of series instead of 1/2 width inside  */
	GMT_LONG check_asym;	/* TRUE to test whether the data are asymmetric about the output time  */
	GMT_LONG check_lack;	/* TRUE to test for lack of data (gap) in the filter window */
	GMT_LONG check_q;	/* TRUE to test average weight or N in median */
	GMT_LONG robust;		/* Look for outliers in data when TRUE */
	GMT_LONG equidist;	/* Data is evenly sampled in t */
	GMT_LONG out_at_time;	/* TRUE when output is required at evenly spaced intervals */
	GMT_LONG f_operator;	/* TRUE if custom weights coefficients sum to zero */

	GMT_LONG	*n_this_col;		/* Pointer to array of counters [one per column]  */
	GMT_LONG	*n_left;		/* Pointer to array of counters [one per column]  */
	GMT_LONG	*n_right;		/* Pointer to array of counters [one per column]  */
	GMT_LONG	n_rows;			/* Number of rows of input  */
	GMT_LONG	n_cols;			/* Number of columns of input  */
	GMT_LONG	t_col;			/* Column of time abscissae (independent variable)  */
	GMT_LONG	n_f_wts;		/* Number of filter weights  */
	GMT_LONG	half_n_f_wts;		/* Half the number of filter weights  */
	GMT_LONG	n_row_alloc;		/* Number of rows of data to allocate  */
	GMT_LONG	n_work_alloc;		/* Number of rows of workspace to allocate  */
	GMT_LONG	filter_type;		/* Flag indicating desired filter type  */
	GMT_LONG	kind;			/* -1 skip +ve, +1 skip -ve, else use all  [for the l|L|u|U filter] */
	GMT_LONG	way;			/* -1 find minimum, +1 find maximum  [for the l|L|u|U filter] */
	GMT_LONG	mode_selection;
	GMT_LONG	n_multiples;
	
	double	*f_wt;			/* Pointer for array of filter coefficients  */
	double	*min_loc;		/* Pointer for array of values, one per [column]  */
	double	*max_loc;
	double	*last_loc;
	double	*this_loc;
	double	*min_scl;
	double	*max_scl;
	double	*last_scl;
	double	*this_scl;
	double	**work;			/* Pointer to array of pointers to doubles for work  */
	double	**data;			/* Pointer to array of pointers to doubles for data  */
	double	dt;			/* Delta time resolution for filter computation  */
	double	q_factor;		/* Quality level for mean weights or n in median  */
	double filter_width;		/* Full width of filter in user's units */
	double half_width;
	double t_start;			/* x-value of first output point */
	double t_stop;			/* x-value of last output point */
	double t_int;			/* Output interval */
	double sym_coeff;		/* Symmetry coefficient  */
	double lack_width;		/* Lack of data width  */
	double extreme;			/* Extreme value [for the l|L|u|U filter] */
	
	FILE	*fp_wt;		/* File pointer to custom weight coefficients (optional) */
};

void allocate_more_work_space (struct FILTER1D_INFO *F);	/* Does just what it says  */

int main (int argc, char **argv)
{
	GMT_LONG error = FALSE;
	
	char	c, txt_a[GMT_TEXT_LEN], txt_b[GMT_TEXT_LEN];
	char	buffer[BUFSIZ];		/* Used to scan input lines  */
	
	GMT_LONG	i, n_fields, n_expected_fields, n_read;
	
	double	last_time, new_time, *in = NULL;
	
	FILE	*fp_in = NULL;		/* File pointer to data file or GMT_stdin */
	
	struct FILTER1D_INFO F;
	struct FILTER1D_CTRL *Ctrl = NULL;

	GMT_LONG	load_data_and_check_extrema (struct FILTER1D_INFO *F);	/* Does just what it says  */
	GMT_LONG	set_up_filter (struct FILTER1D_INFO *F);		/* Creates filter weights or reads them from file and sets start & stop  */
	GMT_LONG	do_the_filter (struct FILTER1D_INFO *F);		/* Does the job  */
	GMT_LONG	allocate_space (struct FILTER1D_INFO *F);		/* Does just what it says  */
	void	free_space (struct FILTER1D_INFO *F);			/* Does just what it says  */
	void	allocate_more_data_space (struct FILTER1D_INFO *F);	/* Does just what it says  */
	void load_parameters (struct FILTER1D_INFO *C, struct FILTER1D_CTRL *Ctrl);
	void *New_filter1d_Ctrl (), Free_filter1d_Ctrl (struct FILTER1D_CTRL *C);
	
	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct FILTER1D_CTRL *)New_filter1d_Ctrl ();	/* Allocate and initialize a new control structure */

	memset ((void *)&F, 0, sizeof (struct FILTER1D_INFO));	/* Init control structure to NULL */
	F.n_row_alloc = GMT_CHUNK;
	F.n_work_alloc = GMT_CHUNK;
	F.equidist = TRUE;
	
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {

				/* Common parameters */

				case 'H':
				case 'V':
				case ':':
				case 'b':
				case 'f':
				case '\0':
					error += GMT_parse_common_options (argv[i], 0, 0, 0, 0);
					break;

				/* Supplemental parameters */

				case 'F':	/* Filter selection  */
					if (strchr ("BbCcGgMmPpLlUuFf", argv[i][2])) {	/* OK filter code */
						Ctrl->F.active = TRUE;
						Ctrl->F.filter = argv[i][2];
						Ctrl->F.width = atof (&argv[i][3]);
						switch (Ctrl->F.filter) {	/* Get some futher info from some filters */
							case 'P':
							case 'p':
								c = argv[i][strlen(argv[i]-1)];
								if (c == '-') Ctrl->F.mode = -1;
								if (c == '+') Ctrl->F.mode = +1;
								break;
							case 'F':
							case 'f':
								Ctrl->F.width = DBL_MAX;	/* To avoid range test errors before reading coefficients */
								if (argv[i][3] && !GMT_access (&argv[i][3], R_OK))
									Ctrl->F.file = strdup (&argv[i][3]);
								else
									error++;
								break;
						}
					}
					else
						error++;
					if (error) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -F option.  Correct syntax: -FX<width>, X one of BbCcGgMmPpFflLuU\n", GMT_program);
						error++;
					}
					break;

				case 'D':
					Ctrl->D.inc = atof (&argv[i][2]);
					Ctrl->D.active = TRUE;
					break;
				case 'E':
					Ctrl->E.active = TRUE;
					break;
				case 'I':
					Ctrl->I.value = atof (&argv[i][2]);
					Ctrl->I.active = TRUE;
					break;
				case 'L':
					Ctrl->L.active = TRUE;
					Ctrl->L.value = atof (&argv[i][2]);
					break;
				case 'N':
					if (sscanf (&argv[i][2], "%" GMT_LL "d/%" GMT_LL "d", &Ctrl->N.ncols, &Ctrl->N.col) != 2) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -N option:  Syntax is -N<ncol>/<tcol>\n", GMT_program);
						error++;
					}
					else
						Ctrl->N.active = TRUE;
					break;
				case 'Q':
					Ctrl->Q.value = atof(&argv[i][2]);
					Ctrl->Q.active = TRUE;
					break;
				case 'S':
					Ctrl->S.active = TRUE;
					Ctrl->S.value = atof (&argv[i][2]);
					break;
				case 'T':
					if (sscanf (&argv[i][2], "%[^/]/%[^/]/%lf", txt_a, txt_b, &Ctrl->T.inc) != 3) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -T option:  Syntax is -T<start>/<stop>/<inc>\n", GMT_program);
						error++;
					}
					else {
						GMT_scanf_arg (txt_a, GMT_IS_UNKNOWN, &Ctrl->T.min);
						GMT_scanf_arg (txt_b, GMT_IS_UNKNOWN, &Ctrl->T.max);
						Ctrl->T.active = TRUE;
					}
					break;
				default:
					error++;
                                        GMT_default_error (argv[i][1]);
					break;
			}
		}
		else {
			if ((fp_in = GMT_fopen(argv[i], GMT_io.r_mode)) == NULL) {
				fprintf (stderr, "%s: Cannot open file %s\n", GMT_program, argv[i]);
				exit (EXIT_FAILURE);
			}
		}
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "filter1d %s - Time domain filtering of 1-D time series\n\n", GMT_VERSION);
		fprintf (stderr, "usage: filter1d [infile] -F<type><width>[<mode>] [-D<increment>] [-E] [%s]\n", GMT_H_OPT);
		fprintf (stderr, "\t[-I<ignore_val>] [-L<lack_width>] [-N<n_cols>/<t_col>] [-Q<q_factor>] [-S<symmetry>]\n");
		fprintf (stderr, "\t[-T<start>/<stop>/<int>] [-V] [%s] [%s]\n\n", GMT_b_OPT, GMT_f_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\t-F sets Filtertype.  Choose from convolution and non-convolution filters\n");
		fprintf (stderr, "\t   and append full filter <width> in same units as time column.\n");
		fprintf (stderr, "\t   Convolution filters:\n");
		fprintf (stderr, "\t     b: Boxcar : Weights are equal.\n");
		fprintf (stderr, "\t     c: Cosine arch : Weights given by cosine arch.\n");
		fprintf (stderr, "\t     g: Gaussian : Weights given by Gaussian function.\n");
		fprintf (stderr, "\t     f<name>: Custom : Weights given in one-column file <name>.\n");
		fprintf (stderr, "\t   Non-convolution filters:\n");
		fprintf (stderr, "\t     m: Median : Return the median value.\n");
		fprintf (stderr, "\t     p: Maximum likelihood probability (mode) estimator : Return the mode.\n");
		fprintf (stderr, "\t        By default, we return the average mode if more than one is found.\n");
		fprintf (stderr, "\t        Append - or + to the width to return the smallest or largest mode instead.\n");
		fprintf (stderr, "\t     l: Lower : Return minimum of all points.\n");
		fprintf (stderr, "\t     L: Lower+ : Return minimum of all positive points.\n");
		fprintf (stderr, "\t     u: Upper : Return maximum of all points.\n");
		fprintf (stderr, "\t     U: Upper- : Return maximum of all negative points.\n");
		fprintf (stderr, "\t   Upper case type B, C, G, M, P, F will use robust filter versions,\n");
		fprintf (stderr, "\t   i.e., replace outliers (2.5 L1 scale off median) with median during filtering.\n");
		
		fprintf (stderr, "\n\tOPTIONS:\n");

		fprintf (stderr, "\t-D used when series is NOT equidistantly sampled.\n");
		fprintf (stderr, "\t   Then <increment> will be the abscissae resolution, i.e. all abscissae\n");
		fprintf (stderr, "\t   will be rounded off to a multiple of <increment>.\n");
		fprintf (stderr, "\t-E include Ends of time series in output.  Default loses half_width at each end.\n");
		GMT_explain_option ('H');
		fprintf (stderr, "\t-I to ignore values; If an input value == <ignore_val> it will be set to NaN.\n");
		fprintf (stderr, "\t-L checks for Lack of data condition.  If input data has a gap exceeding\n");
		fprintf (stderr, "\t   <width> then no output will be given at that point.  Default does not check Lack.\n");
		fprintf (stderr, "\t-N sets # columns in input and which column contains the independent\n");
		fprintf (stderr, "\t   variable (time). The left-most column is # 0, the right-most is # (<n_cols> - 1).\n");
		fprintf (stderr, "\t   Default is <n_cols> = 2, <t_col> = 0; i.e. file has t, f(t) pairs.\n");
		fprintf (stderr, "\t-Q assess Quality of output value by checking mean weight in convolution.\n");
		fprintf (stderr, "\t   Enter <q_factor> between 0 and 1.  If mean weight < q_factor, output is suppressed\n");
		fprintf (stderr, "\t   at this point.  Default does not check Quality.\n");
		fprintf (stderr, "\t-S checks symmetry of data about window center.  Enter a factor\n");
		fprintf (stderr, "\t   between 0 and 1.  If ( (abs(n_left - F->n_right)) / (n_left + F->n_right) ) > factor,\n");
		fprintf (stderr, "\t   then no output will be given at this point.  Default does not check Symmetry.\n");
		fprintf (stderr, "\t-T make evenly spaced timesteps from <start> to <stop> by <int>. Default uses input times.\n");
		GMT_explain_option ('V');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		GMT_explain_option ('o');
		GMT_explain_option ('n');
		GMT_explain_option ('f');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	/* Check arguments */

	if (!Ctrl->F.active) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  -F is required\n", GMT_program);
		error++;
	}
	if (Ctrl->F.active && Ctrl->F.width <= 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -F option: Filterwidth must be positive\n", GMT_program);
		error++;
	}
	if (Ctrl->T.active && (Ctrl->T.max - Ctrl->T.min) < Ctrl->F.width) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -T option:  Output interval < filterwidth\n", GMT_program);
		error++;
	}
	if (Ctrl->L.active && (Ctrl->L.value < 0.0 || Ctrl->L.value > Ctrl->F.width) ) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -L option:  Unreasonable lack-of-data interval\n", GMT_program);
		error++;
	}
	if (Ctrl->S.active && (Ctrl->S.value < 0.0 || Ctrl->S.value > 1.0) ) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -S option:  Enter a factor between 0 and 1\n", GMT_program);
		error++;
	}
	if (Ctrl->N.col >= Ctrl->N.ncols) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -N option.  time column exceeds number of columns\n", GMT_program);
		error++;
	}
	if (Ctrl->Q.active && Ctrl->Q.value < 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -Q option:  Enter a factor between 0 and 1\n", GMT_program);
		error++;
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] < F.n_cols) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data must have at least %ld fields\n", GMT_program, F.n_cols);
		error++;
	}
	if (GMT_io.binary[GMT_IN] && !Ctrl->N.active && GMT_io.ncol[GMT_IN] == 0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Must specify number of columns in binary input data (-bi)\n", GMT_program);
		error++;
	}

	if (error) exit (EXIT_FAILURE);

	if (GMT_io.binary[GMT_IN] && gmtdefs.verbose) {
		char *type[2] = {"double", "single"};
		fprintf (stderr, "%s: Expects %ld-column %s-precision binary data\n", GMT_program, GMT_io.ncol[GMT_IN], type[GMT_io.single_precision[GMT_IN]]);
	}

	load_parameters (&F, Ctrl);	/* Pass parameters from Control structure to Filter structure */
	
	if (strchr ("BCGMPF", Ctrl->F.filter)) {	/* First deal with robustness request */
		F.robust = TRUE;
		Ctrl->F.filter = (char)tolower ((int)Ctrl->F.filter);
	}
	switch (Ctrl->F.filter) {	/* Set filter parameters */
		case 'b':
			F.filter_type = FILTER1D_BOXCAR;
			break;
		case 'c':
			F.filter_type = FILTER1D_COS_ARCH;
			break;
		case 'g':
			F.filter_type = FILTER1D_GAUSSIAN;
			break;
		case 'm':
			F.filter_type = FILTER1D_MEDIAN;
			break;
		case 'p':
			F.filter_type = FILTER1D_MODE;
			F.mode_selection = Ctrl->F.mode;
			break;
		case 'l':
			F.filter_type = FILTER1D_LOWER_ALL;
			F.way = -1;
			F.extreme = DBL_MAX;
			break;
		case 'L':
			F.filter_type = FILTER1D_LOWER_POS;
			F.way = -1;
			F.kind = +1;
			break;
		case 'u':
			F.filter_type = FILTER1D_UPPER_ALL;
			F.way = +1;
			F.extreme = -DBL_MAX;
			break;
		case 'U':
			F.filter_type = FILTER1D_UPPER_NEG;
			F.way = +1;
			F.kind = -1;
			break;
		case 'f':
			if ((F.fp_wt = GMT_fopen (Ctrl->F.file,"r")) == NULL) {
				fprintf (stderr, "%s: Cannot open file %s\n", GMT_program, Ctrl->F.file);
				exit (EXIT_FAILURE);
			}
			break;
	}

	GMT_io.skip_if_NaN[GMT_X] = GMT_io.skip_if_NaN[GMT_Y] = FALSE;	/* Turn off default GMT NaN-handling */
	GMT_io.skip_if_NaN[F.t_col] = TRUE;			/* ... But disallow NaN in "time" column */

#ifdef SET_IO_MODE
	GMT_setmode (GMT_OUT);
#endif
	if (F.filter_type > FILTER1D_CONVOLVE) F.robust = FALSE;

	if (allocate_space (&F)) {
		fprintf (stderr, "%s: fatal error memory allocation.\n", GMT_program);
		free_space (&F);
		exit (EXIT_FAILURE);
	}

	if (fp_in == NULL) {
		fp_in = GMT_stdin;
#ifdef SET_IO_MODE
		GMT_setmode (GMT_IN);
#endif
	}

	if (GMT_io.io_header[GMT_IN]) {
		for (i = 0; i < GMT_io.n_header_recs; i++) {
			GMT_fgets (buffer, BUFSIZ, fp_in);
			if (GMT_io.io_header[GMT_OUT]) GMT_fputs (buffer, GMT_stdout);
		}
	}

	F.n_rows = n_read = 0;
	if (F.robust || (F.filter_type == FILTER1D_MEDIAN) ) {
		for (i = 0; i < F.n_cols; i++) {
			F.min_loc[i] = DBL_MAX;
			F.max_loc[i] = -DBL_MAX;
		}
	}
	last_time = -DBL_MAX;

	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] > 0)
		n_expected_fields = GMT_io.ncol[GMT_IN];
	else
		n_expected_fields = F.n_cols;

	while ((n_fields = GMT_input (fp_in, &n_expected_fields, &in)) >= 0 && !(GMT_io.status & GMT_IO_EOF)) {	/* Not yet EOF */

		n_read++;

		if (GMT_io.status & GMT_IO_MISMATCH) {
			fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld (skipped)\n", GMT_program, n_fields,  n_expected_fields, n_read);
			continue;
		}

		for (i = 0; i < F.n_cols; i++) {
			if (Ctrl->I.active && in[i] == Ctrl->I.value)
				F.data[i][F.n_rows] = GMT_d_NaN;
			else {
				F.data[i][F.n_rows] = in[i];
				if (F.robust || (F.filter_type == FILTER1D_MEDIAN) ) {
					if (in[i] > F.max_loc[i]) F.max_loc[i] = in[i];
					if (in[i] < F.min_loc[i]) F.min_loc[i] = in[i];
				}
			}
		}

		if (!GMT_is_dnan (F.data[F.t_col][F.n_rows]) )
		{	/* Don't n_rows++ if time is bad this row  */
			new_time = F.data[F.t_col][F.n_rows];
			F.n_rows++;
			if (new_time < last_time) {
				fprintf (stderr, "%s: Error! Time decreases at line # %ld\n", GMT_program, F.n_rows);
				fprintf (stderr, "\tUse UNIX utility sort and then try again.\n");
				exit (EXIT_FAILURE);
			}
			last_time = new_time;
		}

		if (F.n_rows == F.n_row_alloc) {	/* Need more memory */
			F.n_row_alloc <<= 1;
			allocate_more_data_space (&F);
		}
	}

	if (fp_in != GMT_stdin) GMT_fclose (fp_in);

	for (i = 0; i < F.n_cols; i++) F.data[i] = (double *) GMT_memory ((void *)F.data[i], (size_t)F.n_rows, sizeof (double), GMT_program);
	if (gmtdefs.verbose) fprintf (stderr, "%s: Read %ld data lines\n", GMT_program, F.n_rows);

	/* Initialize scale parameters and last_loc based on min and max of data  */

	if (F.robust || (F.filter_type == FILTER1D_MEDIAN) ) {
		for (i = 0; i < F.n_cols; i++) {
			F.min_scl[i] = 0.0;
			F.max_scl[i] = 0.5 * (F.max_loc[i] - F.min_loc[i]);
			F.last_scl[i] = 0.5 * F.max_scl[i];
			F.last_loc[i] = 0.5 * (F.max_loc[i] + F.min_loc[i]);
		}
	}

	if (set_up_filter (&F) ) {
		fprintf (stderr, "%s: fatal error during coefficient setup.\n", GMT_program);
		free_space (&F);
		exit (EXIT_FAILURE);
	}

	if (do_the_filter (&F) ) {
		fprintf (stderr, "%s: fatal error in filtering routine.\n", GMT_program);
		free_space (&F);
		exit (EXIT_FAILURE);
	}

	free_space (&F);

	if (F.n_multiples > 0 && gmtdefs.verbose) fprintf (stderr, "%s: WARNING: %ld multiple modes found\n", GMT_program, F.n_multiples);

	Free_filter1d_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

GMT_LONG set_up_filter (struct FILTER1D_INFO *F)
{
	GMT_LONG	i, i1, i2;
	double	t_0, t_1, time, w_sum;
	char buffer[BUFSIZ];
	GMT_LONG	normalize = FALSE;
	PFD get_weight[3];		/* Selects desired weight function.  */

	double	boxcar_weight (double radius, double half_width);	/* Returns the weight for a given delta_time  */
	double	cosine_weight (double radius, double half_width);
	double	gaussian_weight (double radius, double half_width);

	t_0 = F->data[F->t_col][0];
	t_1 = F->data[F->t_col][F->n_rows-1];
	if (F->equidist) F->dt = (t_1 - t_0) / (F->n_rows - 1);

	if (F->fp_wt != NULL) {
		F->f_wt = (double *) GMT_memory (VNULL, (size_t)F->n_work_alloc, sizeof(double), GMT_program);
		F->n_f_wts = 0;
		w_sum = 0.0;
		while (GMT_fgets (buffer, BUFSIZ, F->fp_wt)) {
			if (GMT_is_a_blank_line (buffer)) continue;	/* Skip blank lines or # comments */
			GMT_chop (buffer);
			if (sscanf (buffer, "%lf", &F->f_wt[F->n_f_wts]) != 1) {
				fprintf (stderr, "%s: Error decoding filter weights from file\n", GMT_program);
				exit (EXIT_FAILURE);
			}
			w_sum += F->f_wt[F->n_f_wts];
			F->n_f_wts++;
			if (F->n_f_wts == F->n_work_alloc) {	/* Need more memory */
				F->n_work_alloc <<= 1;
				allocate_more_work_space (F);
			}
		}
		GMT_fclose (F->fp_wt);
		F->f_operator = (w_sum == 0.0);	/* If weights sum to zero it is an operator like {-1 1] or [1 -2 1] */
		F->half_n_f_wts = F->n_f_wts / 2;
		F->half_width = F->half_n_f_wts * F->dt;
		F->filter_width = 2.0 * F->half_width;
		F->f_wt = (double *) GMT_memory ((void *)F->f_wt, (size_t)F->n_f_wts, sizeof (double), GMT_program);
		if (gmtdefs.verbose) fprintf (stderr, "%s: Read %ld filter weights from file.\n", GMT_program, F->n_f_wts);
	}
	else if (F->filter_type <= FILTER1D_CONVOLVE) {
		get_weight[FILTER1D_BOXCAR] = (PFD)boxcar_weight;
		get_weight[FILTER1D_COS_ARCH] = (PFD)cosine_weight;
		get_weight[FILTER1D_GAUSSIAN] = (PFD)gaussian_weight;
		F->half_width = 0.5 * F->filter_width;
		F->half_n_f_wts = (GMT_LONG)floor (F->half_width / F->dt);
		F->n_f_wts = 2 * F->half_n_f_wts + 1;

		F->f_wt = (double *) GMT_memory (VNULL, (size_t)F->n_f_wts, sizeof(double), GMT_program);
		for (i = 0; i <= F->half_n_f_wts; i++) {
			time = i * F->dt;
			i1 = F->half_n_f_wts - i;
			i2 = F->half_n_f_wts + i;
			F->f_wt[i1] = F->f_wt[i2] = ( *get_weight[F->filter_type]) (time, F->half_width);
		}
		if (normalize) {
			w_sum = 0.0;
			for (i = 0; i < F->n_f_wts; i++) w_sum += F->f_wt[i];
			for (i = 0; i < F->n_f_wts; i++) F->f_wt[i] /= w_sum;
		}
	}
	else {
		F->half_width = 0.5 * F->filter_width;
	}

	/* Initialize start/stop time */

	if (F->out_at_time) {
		if (F->use_ends) {
			while (F->t_start < t_0) F->t_start += F->t_int;
			while (F->t_stop > t_1) F->t_stop -= F->t_int;
		}
		else {
			while ( (F->t_start - F->half_width) < t_0) F->t_start += F->t_int;
			while ( (F->t_stop + F->half_width) > t_1) F->t_stop -= F->t_int;
		}
	}
	else {
		if (F->use_ends) {
			F->t_start = t_0;
			F->t_stop = t_1;
		}
		else {
			for (i = 0; (F->data[F->t_col][i] - t_0) < F->half_width; i++);
			F->t_start = F->data[F->t_col][i];
			for (i = F->n_rows - 1; (t_1 - F->data[F->t_col][i]) < F->half_width; i--);
			F->t_stop = F->data[F->t_col][i];
		}
	}

	if (gmtdefs.verbose) fprintf (stderr, "F width: %g Resolution: %g Start: %g Stop: %g\n",
		F->filter_width, F->dt, F->t_start, F->t_stop);

	return(0);
}

GMT_LONG do_the_filter (struct FILTER1D_INFO *F)
{
	GMT_LONG	i_row, left, right, n_l, n_r, iq, i_f_wt;
	GMT_LONG	i_t_output = 0, n_in_filter, n_for_call, n_good_ones;
	GMT_LONG i_col;
	double	time, delta_time, *outval, wt, val, med, scl, small;
	GMT_LONG	*good_one;		/* Pointer to array of logicals [one per column]  */
	double	*wt_sum;		/* Pointer for array of weight sums [each column]  */
	double	*data_sum;		/* Pointer for array of data * weight sums [columns]  */

	GMT_LONG	lack_check (struct FILTER1D_INFO *F, GMT_LONG i_col, GMT_LONG left, GMT_LONG right);	/* Tests for lack of data condition  (optional)  */
	void	get_robust_estimates (struct FILTER1D_INFO *F, GMT_LONG j, GMT_LONG n, GMT_LONG both);	/* Does just what it says  */

	outval = (double *)GMT_memory (VNULL, (size_t)F->n_cols, sizeof (double), GMT_program);
	good_one = (GMT_LONG *) GMT_memory (VNULL, (size_t)F->n_cols, sizeof(GMT_LONG), GMT_program);
	wt_sum = (double *) GMT_memory (VNULL, (size_t)F->n_cols, sizeof(double), GMT_program);
	data_sum = (double *) GMT_memory (VNULL, (size_t)F->n_cols, sizeof(double), GMT_program);

	if(!F->out_at_time) {	/* Position i_t_output at first output time  */
		for(i_t_output = 0; F->data[F->t_col][i_t_output] < F->t_start; i_t_output++);
		small = (F->data[F->t_col][1] - F->data[F->t_col][0]);
	}
	else
		small = F->t_int;

	small *= GMT_CONV_LIMIT;
	time = F->t_start;
	left = right = 0;		/* Left/right end of filter window */

	iq = irint (F->q_factor);

	while (time <= (F->t_stop + small)) {
		while ((time - F->data[F->t_col][left] - small) > F->half_width) left++;
		while (right < F->n_rows && (F->data[F->t_col][right] - time - small) <= F->half_width) right++;
		n_in_filter = right - left;
		if ( (!(n_in_filter)) || (F->check_lack && ( (F->filter_width / n_in_filter) > F->lack_width) ) ) {
			if (F->out_at_time)
				time += F->t_int;
			else {
				i_t_output++;
				time = (i_t_output < F->n_rows) ? F->data[F->t_col][i_t_output] : F->t_stop + 1.0;
			}
			continue;
		}

		for (i_col = 0; i_col < F->n_cols; i_col++) {
			F->n_this_col[i_col] = 0;
			wt_sum[i_col] = 0.0;
			data_sum[i_col] = 0.0;
			if (i_col == F->t_col) {
				good_one[i_col] = FALSE;
			}
			else if (F->check_lack) {
				good_one[i_col] = !(lack_check(F, i_col, left, right));
			}
			else {
				good_one[i_col] = TRUE;
			}
			if (F->check_asym) {
				F->n_left[i_col] = 0;
				F->n_right[i_col] = 0;
			}
		}

		if (F->robust || F->filter_type > FILTER1D_CONVOLVE) {
			if (n_in_filter > F->n_work_alloc) {
				F->n_work_alloc = n_in_filter;
				allocate_more_work_space (F);
			}
			for (i_row = left; i_row < right; i_row++) {
				for (i_col = 0; i_col < F->n_cols; i_col++) {
					if (!(good_one[i_col])) continue;
					if (!GMT_is_dnan (F->data[i_col][i_row])) {
						F->work[i_col][F->n_this_col[i_col]] = F->data[i_col][i_row];
						F->n_this_col[i_col]++;
						if (F->check_asym) {
							if (F->data[F->t_col][i_row] < time) F->n_left[i_col]++;
							if (F->data[F->t_col][i_row] > time) F->n_right[i_col]++;
						}
					}
				}
			}
			if (F->check_asym) {
				for (i_col = 0; i_col < F->n_cols; i_col++) {
					if (!(good_one[i_col])) continue;
					n_l = F->n_left[i_col];
					n_r = F->n_right[i_col];
					if ((((double)GMT_abs(n_l - n_r))/(n_l + n_r)) > F->sym_coeff) good_one[i_col] = FALSE;
				}
			}
			if ( (F->filter_type > FILTER1D_CONVOLVE) && F->check_q) {
				for (i_col = 0; i_col < F->n_cols; i_col++) {
					if (F->n_this_col[i_col] < iq) good_one[i_col] = FALSE;
				}
			}

			for (i_col = 0; i_col < F->n_cols; i_col++) {
				if (good_one[i_col]) {
					n_for_call = F->n_this_col[i_col];
					get_robust_estimates (F, i_col, n_for_call, F->robust);
				}
			}

		}	/* That's it for the robust work  */

		if (F->filter_type > FILTER1D_CONVOLVE) {

			/* Need to count how many good ones; use data_sum area  */

			n_good_ones = 0;
			for (i_col = 0; i_col < F->n_cols; i_col++) {
				if (i_col == F->t_col) {
					data_sum[i_col] = time;
				}
				else if (good_one[i_col]) {
					data_sum[i_col] = F->this_loc[i_col];
					n_good_ones++;
				}
				else {
					data_sum[i_col] = GMT_d_NaN;
				}
			}
			if (n_good_ones) GMT_output (GMT_stdout, F->n_cols, data_sum);
		}
		else {
			if (F->robust) for (i_col = 0; i_col < F->n_cols; i_col++) F->n_this_col[i_col] = 0;

			for (i_row = left; i_row < right; i_row++) {
				delta_time = time - F->data[F->t_col][i_row];
				i_f_wt = F->half_n_f_wts + (GMT_LONG)floor(0.5 + delta_time/F->dt);
				if ( (i_f_wt < 0) || (i_f_wt >= F->n_f_wts) ) continue;

				for(i_col = 0; i_col < F->n_cols; i_col++) {
					if ( !(good_one[i_col]) ) continue;
					if (!GMT_is_dnan (F->data[i_col][i_row])) {
						wt = F->f_wt[i_f_wt];
						val = F->data[i_col][i_row];
						if (F->robust) {
							med = F->this_loc[i_col];
							scl = F->this_scl[i_col];
							val = ((fabs(val-med)) > (2.5 * scl)) ? med : val;
						}
						else if (F->check_asym) {	/* This wasn't already done  */
							if (F->data[F->t_col][i_row] < time) F->n_left[i_col]++;
							if (F->data[F->t_col][i_row] > time) F->n_right[i_col]++;
						}
						wt_sum[i_col] += wt;
						data_sum[i_col] += (wt * val);
						F->n_this_col[i_col]++;
					}
				}
			}
			n_good_ones = 0;
			for (i_col = 0; i_col < F->n_cols; i_col++) {
				if ( !(good_one[i_col]) ) continue;
				if ( !(F->n_this_col[i_col]) ) {
					good_one[i_col] = FALSE;
					continue;
				}
				if (F->check_asym && !(F->robust) ) {
					n_l = F->n_left[i_col];
					n_r = F->n_right[i_col];
					if ((((double)GMT_abs(n_l - n_r))/(n_l + n_r)) > F->sym_coeff) {
						good_one[i_col] = FALSE;
						continue;
					}
				}
				if (F->check_q && ((wt_sum[i_col] / F->n_this_col[i_col]) < F->q_factor)) {
					good_one[i_col] = FALSE;
					continue;
				}
				n_good_ones++;
			}
			if (n_good_ones) {
				for (i_col = 0; i_col < F->n_cols; i_col++) {
					if (i_col == F->t_col)
						outval[i_col] = time;
					else if (good_one[i_col])
						outval[i_col] = (F->f_operator) ? data_sum[i_col] : data_sum[i_col] / wt_sum[i_col];
					else
						outval[i_col] = GMT_d_NaN;
				}
				GMT_output (GMT_stdout, F->n_cols, outval);
			}
		}

		/* Go to next output time */

		if (F->out_at_time)
			time += F->t_int;
		else {
			i_t_output++;
			time = (i_t_output < F->n_rows) ? F->data[F->t_col][i_t_output] : F->t_stop + 1.0;
		}
	}

	GMT_free ((void *)outval);
	GMT_free ((void *)good_one);
	GMT_free ((void *)wt_sum);
	GMT_free ((void *)data_sum);

	return (0);
}

GMT_LONG allocate_space (struct FILTER1D_INFO *F)
{
	GMT_LONG	i;

	F->data = (double **) GMT_memory (VNULL, (size_t)F->n_cols, sizeof(double *), GMT_program);
	for (i = 0; i < F->n_cols; i++) F->data[i] = (double *) GMT_memory (VNULL, (size_t)F->n_row_alloc, sizeof(double), GMT_program);

	F->n_this_col = (GMT_LONG *) GMT_memory (VNULL, (size_t)F->n_cols, sizeof(GMT_LONG), GMT_program);

	if (F->check_asym) F->n_left = (GMT_LONG *) GMT_memory (VNULL, (size_t)F->n_cols, sizeof(GMT_LONG), GMT_program);
	if (F->check_asym) F->n_right = (GMT_LONG *) GMT_memory (VNULL, (size_t)F->n_cols, sizeof(GMT_LONG), GMT_program);

	if (F->robust || (F->filter_type > FILTER1D_CONVOLVE) ) {	/* Then we need workspace  */

		F->work = (double **) GMT_memory (VNULL, (size_t)F->n_cols, sizeof(double *), GMT_program);
		for (i = 0; i < F->n_cols; i++) F->work[i] = (double *) GMT_memory (VNULL, (size_t)F->n_work_alloc, sizeof(double), GMT_program);
		F->min_loc = (double *) GMT_memory (VNULL, (size_t)F->n_cols, sizeof(double), GMT_program);
		F->max_loc = (double *) GMT_memory (VNULL, (size_t)F->n_cols, sizeof(double), GMT_program);
		F->last_loc = (double *) GMT_memory (VNULL, (size_t)F->n_cols, sizeof(double), GMT_program);
		F->this_loc = (double *) GMT_memory (VNULL, (size_t)F->n_cols, sizeof(double), GMT_program);
		F->min_scl = (double *) GMT_memory (VNULL, (size_t)F->n_cols, sizeof(double), GMT_program);
		F->max_scl = (double *) GMT_memory (VNULL, (size_t)F->n_cols, sizeof(double), GMT_program);
		F->this_scl = (double *) GMT_memory (VNULL, (size_t)F->n_cols, sizeof(double), GMT_program);
		F->last_scl = (double *) GMT_memory (VNULL, (size_t)F->n_cols, sizeof(double), GMT_program);
	}
	return (0);
}

void allocate_more_data_space (struct FILTER1D_INFO *F)
{
	GMT_LONG	i;

	for (i = 0; i < F->n_cols; i++) F->data[i] = (double *) GMT_memory ((void *)F->data[i], (size_t)F->n_row_alloc, sizeof(double), GMT_program);
}

void allocate_more_work_space (struct FILTER1D_INFO *F)
{
	GMT_LONG	i;

	for (i = 0; i < F->n_cols; i++) F->work[i] = (double *) GMT_memory ((void *)F->work[i], (size_t)F->n_work_alloc, sizeof(double), GMT_program);
}


void free_space (struct FILTER1D_INFO *F)
{
	GMT_LONG i;
	if (F->robust || (F->filter_type > FILTER1D_CONVOLVE) ) {	
		for (i = 0; i < F->n_cols; i++)	GMT_free ((void *)F->work[i]);
		GMT_free ((void *)F->work);
	}
	for (i = 0; i < F->n_cols; i++)	GMT_free ((void *)F->data[i]);
	GMT_free ((void *)F->data);
	GMT_free ((void *)F->n_this_col);
	if (F->check_asym) GMT_free ((void *)F->n_left);
	if (F->check_asym) GMT_free ((void *)F->n_right);
	if (F->min_loc) GMT_free ((void *)F->min_loc);
	if (F->max_loc) GMT_free ((void *)F->max_loc);
	if (F->last_loc) GMT_free ((void *)F->last_loc);
	if (F->this_loc) GMT_free ((void *)F->this_loc);
	if (F->min_scl) GMT_free ((void *)F->min_scl);
	if (F->max_scl) GMT_free ((void *)F->max_scl);
	if (F->last_scl) GMT_free ((void *)F->last_scl);
	if (F->this_scl) GMT_free ((void *)F->this_scl);
	if (F->n_f_wts) GMT_free ((void *)F->f_wt);
}

GMT_LONG	lack_check (struct FILTER1D_INFO *F, GMT_LONG i_col, GMT_LONG left, GMT_LONG right)
{
	GMT_LONG	last_row, this_row;
	GMT_LONG	lacking = FALSE;
	double	last_t;

	last_row = left;
	while (!(GMT_is_dnan (F->data[i_col][last_row])) && last_row < (right - 1)) last_row++;

	last_t = F->data[F->t_col][last_row];
	this_row = last_row + 1;
	while ( !(lacking) && this_row < (right - 1)) {
		while (!(GMT_is_dnan (F->data[i_col][this_row])) && this_row < (right - 1)) this_row++;

		if ( (F->data[F->t_col][this_row] - last_t) > F->lack_width)
			lacking = TRUE;
		else {
			last_t = F->data[F->t_col][this_row];
			last_row = this_row;
			this_row++;
		}
	}
	return (lacking);
}

void get_robust_estimates (struct FILTER1D_INFO *F, GMT_LONG j, GMT_LONG n, GMT_LONG both)
{
	GMT_LONG i, n_smooth;
	GMT_LONG sort_me = TRUE;
	double	low, high, last, temp;

	if (F->filter_type > FILTER1D_MODE) {
		temp = GMT_extreme (F->work[j], n, F->extreme, F->kind, F->way);
	}
	else if (F->filter_type == FILTER1D_MODE) {

		n_smooth = n / 2;

		GMT_mode (F->work[j], n, n_smooth, sort_me, F->mode_selection, &F->n_multiples, &temp);
	}
	else {
		low = F->min_loc[j];
		high = F->max_loc[j];
		last = F->last_loc[j];

		GMT_median (F->work[j], n, low, high, last, &temp);
	}

	F->last_loc[j] = F->this_loc[j] = temp;

	if (both) {
		for (i = 0; i < n; i++)
			F->work[j][i] = fabs(F->work[j][i] - F->this_loc[j]);
		low = F->min_scl[j];
		high = F->max_scl[j];
		last = F->last_scl[j];
		GMT_median(F->work[j], n, low, high, last, &temp);
		F->last_scl[j] = F->this_scl[j] = temp;
	} 
}

double	boxcar_weight (double radius, double half_width)
{
	double weight;

	if (radius > half_width)
		weight = 0.0;
	else
		weight = 1.0;

	return (weight);
}

double	cosine_weight (double radius, double half_width)
{
	double weight, cosine_constant;
	cosine_constant = M_PI / half_width;
	if (radius > half_width)
		weight = 0.0;
	else
		weight = 1.0 + cos (radius * cosine_constant);

	return (weight);
}

double	gaussian_weight (double radius, double half_width)
{
	double weight, gauss_constant;
	gauss_constant = -4.5 / (half_width * half_width);
	if (radius > half_width)
		weight = 0.0;
	else
		weight = exp (radius * radius * gauss_constant);

	return (weight);
}

void load_parameters (struct FILTER1D_INFO *F, struct FILTER1D_CTRL *Ctrl)
{
	F->filter_width = Ctrl->F.width;
	F->dt = Ctrl->D.inc;
	F->equidist = !Ctrl->D.active;
	F->use_ends = Ctrl->E.active;
	F->check_lack = Ctrl->L.active;
	F->lack_width = Ctrl->L.value;
	F->n_cols = Ctrl->N.ncols;
	F->t_col = Ctrl->N.col;
	F->q_factor = Ctrl->Q.value;
	F->check_q = Ctrl->Q.active;
	F->check_asym = Ctrl->S.active;
	F->sym_coeff = Ctrl->S.value;
	F->t_start = Ctrl->T.min;
	F->t_stop = Ctrl->T.max;
	F->t_int =Ctrl->T.inc;
	F->out_at_time = Ctrl->T.active;
}

void *New_filter1d_Ctrl () {	/* Allocate and initialize a new control structure */
	struct FILTER1D_CTRL *C;

	C = (struct FILTER1D_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct FILTER1D_CTRL), "New_filter1d_Ctrl");

	/* Initialize values whose defaults are not 0/FALSE/NULL */
	C->N.ncols = 2;						/* Expected # of columns */

	return ((void *)C);
}

void Free_filter1d_Ctrl (struct FILTER1D_CTRL *C) {	/* Deallocate control structure */
	if (C->F.file) free ((void *)C->F.file);	
	GMT_free ((void *)C);	
}
