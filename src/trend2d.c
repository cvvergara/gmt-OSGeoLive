/*--------------------------------------------------------------------
 *	$Id: trend2d.c 10173 2014-01-01 09:52:34Z pwessel $
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
 * trend2d [<xyz[w]file>] -F<output_flags> -N<n_m_parameters>[r] 
 *	[-C<condition_#>] [-I[<confid>]] [-V] [-W]
 *
 * where:
 *	[<xyz[w]file>] is an ascii file with x y z in first 3 columns
 *		[or x y z w in first 4 columns].  Default reads from GMT_stdin.
 *	-F<output_flags> is a string of at least one, up to five, in
 *		and order, from the set {x y z m r w}.  x,y,z = input,
 *		m = model, r = residual = z-m, and w= weight used.
 *	-N<n_m_parameters>[r]
 *		If iterative Robust fitting desired, use -N<#>r, else -N.
 *		[Max] Number of terms in the model is <n_m_parameters>.
 *		Example:  Robust bilinear surface:  -N4r.  Max n = 10.
 *	[-C<condition_#>] Cut off eigenvalue spectrum; use only eigen-
 *		values such that (lambda_max / lambda[i]) < condition_#.
 *	[-I[<confid>]] Iteratively Increment the number of model parameters,
 *		searching for the significant model size, up to a maximum
 *		set by <n_m_parameters>.  We start with a 1 parameter
 *		model and then iteratively increase the number of
 *		model parameters, m, while m <= <n_m_parameters> &&
 *		reduction in variance from i to i+1 is significant
 *		at the <confid> level according to F test.  If user sets
 *		-I without giving <confid> then <confid> = 0.95.
 *	[-V]	Verbose operation.
 *	[-W]	Weighted data are input.  Read 4 cols and use 4th as weight.
 *
 *
 * Read GMT_stdin or file of x y z triples, or weighted data, x y z w.  Fit 
 * a regression model z = f(x,y) + e, where e are error misfits and f(x,y)
 * has some user-prescribed functional form.  The user may choose the number
 * of terms in the model to fit, whether to seek iterative refinement robust
 * w.r.t. outliers, and whether to seek automatic discovery of the significant
 * number of model parameters.
 *
 * Adapted from trend1d by w. h. f. smith.
 *
 *
 * During model fitting the data x,y coordinates are normalized into the domain
 * [-1, 1] for Chebyshev Polynomial fitting.  Before writing out the data the
 * coordinates are rescaled to match the original input values.
 *
 *
 * Author:	W. H. F. Smith
 * Date:	17 June 1991-2000.
 * Revised:	12-JUN-1998 for GMT 3.1 (PW)
 *		10-JUL-2000 for GMT 3.3.5 (PW) Added -L option
 * Version:	4
 */

#include "gmt.h"

#define TREND2D_N_OUTPUT_CHOICES 6

struct TREND2D_CTRL {
	struct C {	/* -C<condition_#> */
		GMT_LONG active;
		double value;
	} C;
	struct F {	/* -F<xymrw> */
		GMT_LONG active;
		char col[TREND2D_N_OUTPUT_CHOICES];	/* Character codes for desired output in the right order */
	} F;
	struct I {	/* -I[<confidence>] */
		GMT_LONG active;
		double value;
	} I;
	struct N {	/* -N<n_model>[r] */
		GMT_LONG active;
		GMT_LONG robust;
		GMT_LONG value;
	} N;
	struct W {	/* -W */
		GMT_LONG active;
	} W;
};

struct	TREND2D_DATA {
	double	x;
	double	y;
	double	z;
	double	m;
	double	r;
	double	w;
};

int main(int argc, char **argv)
{
	GMT_LONG	i, j, k, n_outputs, n_model, significant, rank, n_req, np;
	GMT_LONG	n_data;

	GMT_LONG	error = FALSE, weighted_output = FALSE;

	double	*gtg = NULL, *v = NULL, *gtd = NULL, *lambda = NULL, *workb = NULL, *workz = NULL, *c_model = NULL, *o_model = NULL, *w_model = NULL, *work = NULL;	/* Arrays  */
	double	xmin, xmax, ymin, ymax, c_chisq, o_chisq = 0.0, w_chisq, scale = 1.0, prob;
	double	get_chisq(struct TREND2D_DATA *data, GMT_LONG n_data, GMT_LONG n_model);

	char format[BUFSIZ];

	FILE	*fp = NULL;

	struct	TREND2D_DATA *data = NULL;
	struct TREND2D_CTRL *Ctrl = NULL;

	void read_data(struct TREND2D_DATA **data, GMT_LONG *n_data, double *xmin, double *xmax, double *ymin, double *ymax, GMT_LONG weighted_input, double **work, FILE *fp);
	void write_output(struct TREND2D_DATA *data, GMT_LONG n_data, char *output_choice, GMT_LONG n_outputs);
	void transform_x(struct TREND2D_DATA *data, GMT_LONG n_data, double xmin, double xmax, double ymin, double ymax);
	void untransform_x(struct TREND2D_DATA *data, GMT_LONG n_data, double xmin, double xmax, double ymin, double ymax);
	void recompute_weights(struct TREND2D_DATA *data, GMT_LONG n_data, double *work, double *scale);
	void allocate_array_space(GMT_LONG np, double **gtg, double **v, double **gtd, double **lambda, double **workb, double **workz, double **c_model, double **o_model, double **w_model);
	void free_the_memory(double *gtg, double *v, double *gtd, double *lambda, double *workb, double *workz, double *c_model, double *o_model, double *w_model, struct TREND2D_DATA *data, double *work);
	void calc_m_and_r(struct TREND2D_DATA *data, GMT_LONG n_data, double *model, GMT_LONG n_model, double *grow);
	void move_model_a_to_b(double *model_a, double *model_b, GMT_LONG n_model, double *chisq_a, double *chisq_b);
	void load_gtg_and_gtd(struct TREND2D_DATA *data, GMT_LONG n_data, double *gtg, double *gtd, double *grow, GMT_LONG n_model, GMT_LONG mp);
	void solve_system(double *gtg, double *gtd, double *model, GMT_LONG n_model, GMT_LONG mp, double *lambda, double *v, double *b, double *z, double c_no, GMT_LONG *ir);
	void *New_trend2d_Ctrl (), Free_trend2d_Ctrl (struct TREND2D_CTRL *C);
	
	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct TREND2D_CTRL *)New_trend2d_Ctrl ();	/* Allocate and initialize a new control structure */

	n_outputs = 0;
	sprintf (format, "%s\t", gmtdefs.d_format);

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
					error += (GMT_LONG)GMT_parse_common_options (argv[i], 0, 0, 0, 0);
					break;

				/* Supplemental parameters */


				case 'C':
					Ctrl->C.active = TRUE;
					Ctrl->C.value = atof(&argv[i][2]);
					break;
				case 'F':
					Ctrl->F.active = TRUE;
					for (j = 2, k = 0; argv[i][j]; j++, k++) {
						if (k < TREND2D_N_OUTPUT_CHOICES)
							Ctrl->F.col[k] = argv[i][j];
						else {
							error++;
							fprintf (stderr, "%s: GMT SYNTAX ERROR -F option: Too many output columns selected\n", GMT_program);
							fprintf (stderr, "%s: GMT SYNTAX ERROR -F option: Choose from -Fxyzmrw\n", GMT_program);
						}
					}
					break;
				case 'I':
					Ctrl->I.active = TRUE;
					Ctrl->I.value = (argv[i][2]) ? atof(&argv[i][2]) : 0.51;
					break;
				case 'L':
					GMT_io.in_col_type[0] = GMT_io.out_col_type[0] = GMT_IS_LON;
					GMT_io.in_col_type[1] = GMT_io.out_col_type[1] = GMT_IS_LAT;
					fprintf (stderr, "%s: Option -L is obsolete (but is processed correctly).  Please use -f instead\n", GMT_program);
					break;
				case 'N':
					Ctrl->N.active = TRUE;
					if (strchr (argv[i], 'r')) Ctrl->N.robust = TRUE;
					j = (argv[i][2] == 'r') ? 3 : 2;
					Ctrl->N.value = (argv[i][j]) ? atoi(&argv[i][j]) : 0;
					break;
				case 'W':
					Ctrl->W.active = TRUE;
					break;
				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else {
			if ((fp = GMT_fopen(argv[i], GMT_io.r_mode)) == NULL) {
				fprintf (stderr, "%s:  Could not open file %s\n", GMT_program, argv[i]);
				error = TRUE;
			}
		}
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf(stderr,"trend2d %s - Fit a [weighted] [robust] polynomial for z = f(x,y) to ascii xyz[w]\n\n", GMT_VERSION);
		fprintf(stderr,"usage:  trend2d -F<xyzmrw> -N<n_model>[r] [<xyz[w]file>] [-C<condition_#>] [%s] [-I[<confidence>]]\n\n", GMT_H_OPT);
		fprintf(stderr,"	[-V] [-W] [%s] [%s] [%s]\n\n", GMT_t_OPT, GMT_b_OPT, GMT_f_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf(stderr,"\t-F Choose at least 1, up to 6, any order, of xyzmrw for ascii output to stdout.\n");
		fprintf(stderr,"\t-N fit a [robust] model with <n_model> terms.  <n_model> in [1,10].  E.g., robust planar = -N3r.\n");
		fprintf(stderr,"\t   Model parameters order is given as follows:\n");
		fprintf(stderr,"\t   z = m1 + m2*x + m3*y + m4*x*y + m5*x^2 + m6*y^2 + m7*x^3 + m8*x^2*y + m9*x*y^2 + m10*y^3.\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf(stderr,"\t[<xyz[w]file>] name of ascii file, first 3 cols = x y z [4 cols = x y z w]; [Default reads stdin].\n");
		fprintf(stderr,"\t   x=x, y=y, z=z, m=model, r=residual=z-m, w=weight.  w determined iteratively if robust fit used.\n");
		fprintf(stderr,"\t-C Truncate eigenvalue spectrum so matrix has <condition_#> [Default = 1.0e06].\n");
		GMT_explain_option ('H');
		fprintf(stderr,"\t-I Iteratively Increase # model parameters, to a max of <n_model> so long as the\n");
		fprintf(stderr,"\t   reduction in variance is significant at the <confidence> level.\n");
		fprintf(stderr,"\t   Give -I without a number to default to 0.51 confidence level.\n");
		GMT_explain_option ('V');
		fprintf(stderr,"\t-W Weighted input given, weights in 4th column.  [Default is unweighted].\n");
		GMT_explain_option (':');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf(stderr,"\t   Default is 3 (or 4 if -W is set) columns.\n");
		GMT_explain_option ('o');
		GMT_explain_option ('n');
		GMT_explain_option ('f');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	if (Ctrl->C.value <= 1.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -C option.  Condition number must be larger than unity\n", GMT_program);
		error++;
	}
	if (Ctrl->I.value < 0.0 || Ctrl->I.value > 1.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -C option.  Give 0 < confidence level < 1.0\n", GMT_program);
		error++;
	}
	if (Ctrl->N.value <= 0 || Ctrl->N.value > 10) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -N option.  Must request 1-10 parameters\n", GMT_program);
		error++;
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
	n_req = (Ctrl->W.active) ? 4 : 3;
	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = n_req;
	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] < n_req) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data (-bi) must have at least %ld columns\n", GMT_program, n_req);
		error++;
	}

	for (k = n_outputs = 0; k < TREND2D_N_OUTPUT_CHOICES && Ctrl->F.col[k]; k++) {
		if (!strchr ("xyzmrw", Ctrl->F.col[k])) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR -F option.  Unrecognized output choice %c\n", GMT_program, Ctrl->F.col[k]);
			error++;
		}
		else if (Ctrl->F.col[k] == 'w')
			weighted_output = TRUE;

		n_outputs++;
	}
        if (n_outputs == 0) {
                fprintf (stderr, "%s: GMT SYNTAX ERROR -F option.  Must specify at least one output column\n", GMT_program);
                error++;
        }
	if (n_outputs > TREND2D_N_OUTPUT_CHOICES) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -F option.  Too many output columns specified (%ld)\n", GMT_program, n_outputs);
		error++;
	}
	if (error) exit (EXIT_FAILURE);

	if (GMT_io.binary[GMT_IN] && gmtdefs.verbose) {
		char *type[2] = {"double", "single"};
		fprintf (stderr, "%s: Expects %ld-column %s-precision binary data\n", GMT_program, GMT_io.ncol[GMT_IN], type[GMT_io.single_precision[GMT_IN]]);
	}

#ifdef SET_IO_MODE
	GMT_setmode (GMT_OUT);
#endif
	project_info.w = 0;	project_info.e = 360.0;	/* For -L not to cause trouble in GMT_input */
	np = Ctrl->N.value;	/* Row dimension for matrices gtg and v  */
	allocate_array_space((GMT_LONG)np, &gtg, &v, &gtd, &lambda, &workb, &workz, &c_model, &o_model, &w_model);

	read_data(&data, &n_data, &xmin, &xmax, &ymin, &ymax, Ctrl->W.active, &work, fp);

	if (xmin == xmax || ymin == ymax) {
		fprintf(stderr,"%s:  Fatal error in input data.  X min = X max.\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (n_data == 0) {
		fprintf(stderr,"%s:  Fatal error.  Could not read any data.\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (n_data < Ctrl->N.value) {
		fprintf(stderr,"%s:  Warning.  Ill-posed problem.  n_data < n_model_max.\n", GMT_program);
	}

	transform_x(data, n_data, xmin, xmax, ymin, ymax);	/* Set domain to [-1, 1] or [-pi, pi]  */

	if (gmtdefs.verbose) {
		fprintf(stderr,"%s:  Read %ld data with X values from %.8g to %.8g\n", GMT_program, n_data, xmin, xmax);
		fprintf(stderr,"N_model\tRank\tChi_Squared\tSignificance\n");
	}

	sprintf (format, "%%ld\t%%ld\t%s\t%s\n", gmtdefs.d_format, gmtdefs.d_format);

	if (Ctrl->I.active) {
		n_model = 1;

		/* Fit first model  */
		load_gtg_and_gtd(data, n_data, gtg, gtd, workb, n_model, np);
		solve_system(gtg, gtd, c_model, n_model, np, lambda, v, workb, workz, Ctrl->C.value, &rank);
		calc_m_and_r(data, n_data, c_model, n_model, workb);
		c_chisq = get_chisq(data, n_data, n_model);
		if (gmtdefs.verbose) fprintf(stderr, format, n_model, rank, c_chisq, 1.0);
		if (Ctrl->N.robust) {
			do {
				recompute_weights(data, n_data, work, &scale);
				move_model_a_to_b(c_model, w_model, n_model, &c_chisq, &w_chisq);
				load_gtg_and_gtd(data, n_data, gtg, gtd, workb, n_model, np);
				solve_system(gtg, gtd, c_model, n_model, np, lambda, v, workb, workz, Ctrl->C.value, &rank);
				calc_m_and_r(data, n_data, c_model, n_model, workb);
				c_chisq = get_chisq(data, n_data, n_model);
				significant = GMT_sig_f(c_chisq, n_data-n_model, w_chisq, n_data-n_model, Ctrl->I.value, &prob);
				if (gmtdefs.verbose) fprintf(stderr, format, n_model, rank, c_chisq, prob);
			} while (significant);
			/* Go back to previous model only if w_chisq < c_chisq  */
			if (w_chisq < c_chisq) {
				move_model_a_to_b(w_model, c_model, n_model, &w_chisq, &c_chisq);
				calc_m_and_r(data, n_data, c_model, n_model, workb);
				if (weighted_output && n_model == Ctrl->N.value) recompute_weights(data, n_data, work, &scale);
			}
		}
		/* First [robust] model has been found  */

		significant = TRUE;
		while(n_model < Ctrl->N.value && significant) {
			move_model_a_to_b(c_model, o_model, n_model, &c_chisq, &o_chisq);
			n_model++;

			/* Fit next model  */
			load_gtg_and_gtd(data, n_data, gtg, gtd, workb, n_model, np);
			solve_system(gtg, gtd, c_model, n_model, np, lambda, v, workb, workz, Ctrl->C.value, &rank);
			calc_m_and_r(data, n_data, c_model, n_model, workb);
			c_chisq = get_chisq(data, n_data, n_model);
			if (gmtdefs.verbose) fprintf(stderr, format, n_model, rank, c_chisq, 1.0);
			if (Ctrl->N.robust) {
				do {
					recompute_weights(data, n_data, work, &scale);
					move_model_a_to_b(c_model, w_model, n_model, &c_chisq, &w_chisq);
					load_gtg_and_gtd(data, n_data, gtg, gtd, workb, n_model, np);
					solve_system(gtg, gtd, c_model, n_model, np, lambda, v, workb, workz, Ctrl->C.value, &rank);
					calc_m_and_r(data, n_data, c_model, n_model, workb);
					c_chisq = get_chisq(data, n_data, n_model);
					significant = GMT_sig_f(c_chisq, n_data-n_model, w_chisq, n_data-n_model, Ctrl->I.value, &prob);
					if (gmtdefs.verbose) fprintf(stderr, format, n_model, rank, c_chisq, prob);
				} while (significant);
				/* Go back to previous model only if w_chisq < c_chisq  */
				if (w_chisq < c_chisq) {
					move_model_a_to_b(w_model, c_model, n_model, &w_chisq, &c_chisq);
					calc_m_and_r(data, n_data, c_model, n_model, workb);
					if (weighted_output && n_model == Ctrl->N.value) recompute_weights(data, n_data, work, &scale);
				}
			}
			/* Next [robust] model has been found  */
			significant = GMT_sig_f(c_chisq, n_data-n_model, o_chisq, n_data-n_model-1, Ctrl->I.value, &prob);
		}

		if (!(significant) ) {	/* Go back to previous [robust] model, stored in o_model  */
			n_model--;
			rank--;
			move_model_a_to_b(o_model, c_model, n_model, &o_chisq, &c_chisq);
			calc_m_and_r(data, n_data, c_model, n_model, workb);
			if (Ctrl->N.robust && weighted_output) recompute_weights(data, n_data, work, &scale);
		}
	}
	else {
		n_model = Ctrl->N.value;
		load_gtg_and_gtd(data, n_data, gtg, gtd, workb, n_model, np);
		solve_system(gtg, gtd, c_model, n_model, np, lambda, v, workb, workz, Ctrl->C.value, &rank);
		calc_m_and_r(data, n_data, c_model, n_model, workb);
		c_chisq = get_chisq(data, n_data, n_model);
		if (gmtdefs.verbose) fprintf(stderr, format, n_model, rank, c_chisq, 1.0);
		if (Ctrl->N.robust) {
			do {
				recompute_weights(data, n_data, work, &scale);
				move_model_a_to_b(c_model, w_model, n_model, &c_chisq, &w_chisq);
				load_gtg_and_gtd(data, n_data, gtg, gtd, workb, n_model, np);
				solve_system(gtg, gtd, c_model, n_model, np, lambda, v, workb, workz, Ctrl->C.value, &rank);
				calc_m_and_r(data, n_data, c_model, n_model, workb);
				c_chisq = get_chisq(data, n_data, n_model);
				significant = GMT_sig_f(c_chisq, n_data-n_model, w_chisq, n_data-n_model, Ctrl->I.value, &prob);
				if (gmtdefs.verbose) fprintf(stderr, format, n_model, rank, c_chisq, prob);
			} while (significant);
			/* Go back to previous model only if w_chisq < c_chisq  */
			if (w_chisq < c_chisq) {
				move_model_a_to_b(w_model, c_model, n_model, &w_chisq, &c_chisq);
				calc_m_and_r(data, n_data, c_model, n_model, workb);
				if (weighted_output && n_model == Ctrl->N.value) recompute_weights(data, n_data, work, &scale);
			}
		}
	}

	if (gmtdefs.verbose) {
		sprintf (format, "%%s: Final model stats:  N model parameters %%d.  Rank %%d.  Chi-Squared:  %s\n", gmtdefs.d_format);
		fprintf(stderr, format, GMT_program, n_model, rank, c_chisq);
		fprintf(stderr,"Model Coefficients: ");
		sprintf (format, "%s\t", gmtdefs.d_format);
		for (i = 0; i < n_model; i++) {
			fprintf(stderr, format, c_model[i]);
		}
		fprintf(stderr,"\n");
	}

	untransform_x(data, n_data, xmin, xmax, ymin, ymax);

	write_output(data, n_data, Ctrl->F.col, n_outputs);

	free_the_memory(gtg, v, gtd, lambda, workb, workz, c_model, o_model, w_model, data, work);

	Free_trend2d_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void read_data (struct TREND2D_DATA **data, GMT_LONG *n_data, double *xmin, double *xmax, double *ymin, double *ymax, GMT_LONG weighted_input, double **work, FILE *fp)
{

	GMT_LONG	n_read, i, n_alloc = GMT_CHUNK;
	GMT_LONG n_expected_fields, n_fields;
	char	buffer[BUFSIZ];
	double	*in;

	if (fp == NULL) {
		fp = GMT_stdin;
#ifdef SET_IO_MODE
		GMT_setmode (GMT_IN);
#endif
	}
	(*data) = (struct TREND2D_DATA *) GMT_memory (VNULL, (size_t)n_alloc, sizeof(struct TREND2D_DATA), GMT_program);

	if (GMT_io.io_header[GMT_IN]) for (i = 0; i < GMT_io.n_header_recs; i++) GMT_fgets (buffer, BUFSIZ, fp);
	i = n_read = 0;
	n_expected_fields = (GMT_io.binary[GMT_IN]) ? GMT_io.ncol[GMT_IN] : 3 + weighted_input;

	while ((n_fields = GMT_input (fp, &n_expected_fields, &in)) >= 0 && !(GMT_io.status & GMT_IO_EOF)) {

		n_read++;
		if (GMT_io.status & GMT_IO_MISMATCH) {
			fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld (skipped)\n", GMT_program, n_fields, n_expected_fields, n_read);
			continue;
		}

		(*data)[i].x = in[GMT_X];
		(*data)[i].y = in[GMT_Y];
		(*data)[i].z = in[GMT_Z];
		(*data)[i].w = (weighted_input) ? in[3] : 1.0;

		if (i) {
			if (*xmin > (*data)[i].x) *xmin = (*data)[i].x;
			if (*xmax < (*data)[i].x) *xmax = (*data)[i].x;
			if (*ymin > (*data)[i].y) *ymin = (*data)[i].y;
			if (*ymax < (*data)[i].y) *ymax = (*data)[i].y;
		}
		else {
			*xmin = (*data)[i].x;
			*xmax = (*data)[i].x;
			*ymin = (*data)[i].y;
			*ymax = (*data)[i].y;
		}
		i++;

		if (i == n_alloc) {
			n_alloc <<= 1;
			*data = (struct TREND2D_DATA *) GMT_memory ((void *)*data, (size_t)n_alloc, sizeof(struct TREND2D_DATA), GMT_program);
		}
		if (i == INT_MAX) {
			fprintf (stderr, "%s: ERROR: Cannot process more than %d data points\n", GMT_program, INT_MAX);
			GMT_free ((void *)data);
			exit (EXIT_FAILURE);
		}
	}
	if (fp != GMT_stdin) GMT_fclose(fp);
	*data = (struct TREND2D_DATA *) GMT_memory ((void *)*data, (size_t)i, sizeof(struct TREND2D_DATA), GMT_program);
	*work = (double *) GMT_memory (VNULL, (size_t)i, sizeof(double), GMT_program);
	*n_data = i;
}

void allocate_array_space (GMT_LONG np, double **gtg, double **v, double **gtd, double **lambda, double **workb, double **workz, double **c_model, double **o_model, double **w_model)
{
	*gtg = (double *) GMT_memory (VNULL, (size_t)(np*np), sizeof(double), GMT_program);
	*v = (double *) GMT_memory (VNULL, (size_t)(np*np), sizeof(double), GMT_program);
	*gtd = (double *) GMT_memory (VNULL, (size_t)np, sizeof(double), GMT_program);
	*lambda = (double *) GMT_memory (VNULL, (size_t)np, sizeof(double), GMT_program);
	*workb = (double *) GMT_memory (VNULL, (size_t)np, sizeof(double), GMT_program);
	*workz = (double *) GMT_memory (VNULL, (size_t)np, sizeof(double), GMT_program);
	*c_model = (double *) GMT_memory (VNULL, (size_t)np, sizeof(double), GMT_program);
	*o_model = (double *) GMT_memory (VNULL, (size_t)np, sizeof(double), GMT_program);
	*w_model = (double *) GMT_memory (VNULL, (size_t)np, sizeof(double), GMT_program);
}

void write_output (struct TREND2D_DATA *data, GMT_LONG n_data, char *output_choice, GMT_LONG n_outputs)
{
	GMT_LONG	i;
	GMT_LONG j;
	double out[6];

	for (i = 0; i < n_data; i++) {
		for (j = 0; j < n_outputs; j++) {
			switch (output_choice[j]) {
				case 'x':
					out[j] = data[i].x;
					break;
				case 'y':
					out[j] = data[i].y;
					break;
				case 'z':
					out[j] = data[i].z;
					break;
				case 'm':
					out[j] = data[i].m;
					break;
				case 'r':
					out[j] = data[i].r;
					break;
				case 'w':
					out[j] = data[i].w;
					break;
			}
		}
		GMT_output (GMT_stdout, n_outputs, out);
	}
}

void free_the_memory (double *gtg, double *v, double *gtd, double *lambda, double *workb, double *workz, double *c_model, double *o_model, double *w_model, struct TREND2D_DATA *data, double *work)
{
	GMT_free ((void *)work);
	GMT_free ((void *)data);
	GMT_free ((void *)w_model);
	GMT_free ((void *)o_model);
	GMT_free ((void *)c_model);
	GMT_free ((void *)workz);
	GMT_free ((void *)workb);
	GMT_free ((void *)lambda);
	GMT_free ((void *)gtd);
	GMT_free ((void *)v);
	GMT_free ((void *)gtg);
}

void transform_x (struct TREND2D_DATA *data, GMT_LONG n_data, double xmin, double xmax, double ymin, double ymax)
{
	GMT_LONG	i;
	double	offsetx, scalex;
	double	offsety, scaley;

	offsetx = 0.5 * (xmin + xmax);	/* Mid Range  */
	offsety = 0.5 * (ymin + ymax);
	scalex = 2.0 / (xmax - xmin);	/* 1 / (1/2 Range)  */
	scaley = 2.0 / (ymax - ymin);

	for (i = 0; i < n_data; i++) {
		data[i].x = (data[i].x - offsetx) * scalex;
		data[i].y = (data[i].y - offsety) * scaley;
	}
}

void untransform_x (struct TREND2D_DATA *data, GMT_LONG n_data, double xmin, double xmax, double ymin, double ymax)
{
	GMT_LONG	i;
	double	offsetx, scalex;
	double	offsety, scaley;

	offsetx = 0.5 * (xmin + xmax);	/* Mid Range  */
	offsety = 0.5 * (ymin + ymax);
	scalex = 0.5 * (xmax - xmin);	/* 1/2 Range  */
	scaley = 0.5 * (ymax - ymin);

	for (i = 0; i < n_data; i++) {
		data[i].x = (data[i].x * scalex) + offsetx;
		data[i].y = (data[i].y * scaley) + offsety;
	}
}

double get_chisq (struct TREND2D_DATA *data, GMT_LONG n_data, GMT_LONG n_model)
{
	GMT_LONG	i, nu;
	double	chi = 0.0;


	for (i = 0; i < n_data; i++) {	/* Weight is already squared  */
		if (data[i].w == 1.0) {
			chi += (data[i].r * data[i].r);
		}
		else {
			chi += (data[i].r * data[i].r * data[i].w);
		}
	}
	nu = n_data - n_model;
	if (nu > 1) return(chi/nu);
	return(chi);
}

void recompute_weights (struct TREND2D_DATA *data, GMT_LONG n_data, double *work, double *scale)
{
	GMT_LONG	i;
	double	k, ksq, rr;

	/* First find median { fabs(data[].r) },
		estimate scale from this,
		and compute chisq based on this.  */ 

	for (i = 0; i < n_data; i++) {
		work[i] = fabs(data[i].r);
	}
	qsort((void *)work, (size_t)n_data, sizeof(double), GMT_comp_double_asc);

	if (n_data%2) {
		*scale = 1.4826 * work[n_data/2];
	}
	else {
		*scale = 0.7413 * (work[n_data/2 - 1] + work[n_data/2]);
	}

	k = 1.5 * (*scale);	/*  Huber[1964] weight; 95% efficient for Normal data  */
	ksq = k * k;

	for (i = 0; i < n_data; i++) {
		rr = fabs(data[i].r);
		if (rr <= k) {
			data[i].w = 1.0;
		}
		else {
			data[i].w = (2*k/rr) - (ksq/(rr*rr) );	/* This is really w-squared  */
		}
	}
}

void load_g_row (double x, double y, GMT_LONG n, double *gr)
      	     	/* Current data position, appropriately normalized.  */
   	  	/* Number of model parameters, and elements of gr[]  */
      	     	/* Elements of row of G matrix.  */
{
	/* Routine computes the elements gr[j] in the ith row of the
		G matrix (Menke notation), where x,y is the ith datum's
		location.  */

	GMT_LONG	j;

	j = 0;
	while (j < n) {
		switch (j) {
			case 0:
				gr[j] = 1.0;
				break;
			case 1:
				gr[j] = x;
				break;
			case 2:
				gr[j] = y;
				break;
			case 3:
				gr[j] = x*y;
				break;
			case 4:
				gr[j] = 2 * x * gr[1] - gr[0];
				break;
			case 5:
				gr[j] = 2 * y * gr[2] - gr[0];
				break;
			case 6:
				gr[j] = 2 * x * gr[4] - gr[1];
				break;
			case 7:
				gr[j] = gr[4] * gr[2];
				break;
			case 8:
				gr[j] = gr[5] * gr[1];
				break;
			case 9:
				gr[j] = 2 * y * gr[5] - gr[2];
				break;
		}
		j++;
	}
}

void calc_m_and_r (struct TREND2D_DATA *data, GMT_LONG n_data, double *model, GMT_LONG n_model, double *grow)
{
	/*	model[n_model] holds solved coefficients of m_type model.
		grow[n_model] is a vector for a row of G matrix.  */

	GMT_LONG	i, j;
	for (i = 0; i < n_data; i++) {
		load_g_row(data[i].x, data[i].y, n_model, grow);
		data[i].m = 0.0;
		for (j = 0; j < n_model; j++) {
			data[i].m += model[j]*grow[j];
		}
		data[i].r = data[i].z - data[i].m;
	}
}

void move_model_a_to_b (double *model_a, double *model_b, GMT_LONG n_model, double *chisq_a, double *chisq_b)
{
	GMT_LONG	i;
	for(i = 0; i<  n_model; i++) {
		model_b[i] = model_a[i];
	}
	*chisq_b = *chisq_a;
}

void load_gtg_and_gtd (struct TREND2D_DATA *data, GMT_LONG n_data, double *gtg, double *gtd, double *grow, GMT_LONG n_model, GMT_LONG mp)
      	    	      
      			  
   			    	/* mp is row dimension of gtg  */
{

	GMT_LONG	i;
	GMT_LONG j, k;
	double	wz;

	/* First zero the contents for summing:  */

	for (j = 0; j < n_model; j++) {
		for (k = 0; k < n_model; k++) {
			gtg[j + k*mp] = 0.0;
		}
		gtd[j] = 0.0;
	}

	/* Sum over all data  */
	for (i = 0; i < n_data; i++) {
		load_g_row(data[i].x, data[i].y, n_model, grow);
		if (data[i].w != 1.0) {
			wz = data[i].w * data[i].z;
			for (j = 0; j < n_model; j++) {
				for (k = 0; k < n_model; k++) {
					gtg[j + k*mp] += (data[i].w * grow[j] * grow[k]);
				}
				gtd[j] += (wz * grow[j]);
			}
		}
		else {
			for (j = 0; j < n_model; j++) {
				for (k = 0; k < n_model; k++) {
					gtg[j + k*mp] += (grow[j] * grow[k]);
				}
				gtd[j] += (data[i].z * grow[j]);
			}
		}
	}
}

void solve_system (double *gtg, double *gtd, double *model, GMT_LONG n_model, GMT_LONG mp, double *lambda, double *v, double *b, double *z, double c_no, GMT_LONG *ir)
{

	GMT_LONG	i, j, k, rank = 0, nrots;
	GMT_LONG n, m;
	double	c_test, temp_inverse_ij;

	if (n_model == 1) {
		model[0] = gtd[0] / gtg[0];
		*ir = 1;
	}
	else {
		n = n_model;
		m = mp;
		if(GMT_jacobi(gtg, &n, &m, lambda, v, b, z, &nrots)) {
			fprintf(stderr,"trend2d:  Warning:  Matrix Solver Convergence Failure.\n");
		}
		c_test = fabs(lambda[0])/c_no;
		while(rank < n_model && lambda[rank] > 0.0 && lambda[rank] > c_test) rank++;
		for (i = 0; i < n_model; i++) {
			model[i] = 0.0;
			for (j = 0; j < n_model; j++) {
				temp_inverse_ij = 0.0;
				for (k = 0; k <  rank; k++) {
					temp_inverse_ij += (v[i + k*mp] * v[j + k*mp] / lambda[k]);
				}
				model[i] += (temp_inverse_ij * gtd[j]);
			}
		}
		*ir = rank;
	}
}

void *New_trend2d_Ctrl () {	/* Allocate and initialize a new control structure */
	struct TREND2D_CTRL *C;
	
	C = (struct TREND2D_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct TREND2D_CTRL), "New_trend2d_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	C->C.value = 1.0e06;		/* Condition number for matrix solution  */	
	C->I.value = 0.51;		/* Confidence interval for significance test  */
	return ((void *)C);
}

void Free_trend2d_Ctrl (struct TREND2D_CTRL *C) {	/* Deallocate control structure */
	GMT_free ((void *)C);	
}
