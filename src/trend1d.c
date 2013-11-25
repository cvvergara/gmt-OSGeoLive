/*--------------------------------------------------------------------
 *	$Id: trend1d.c 10008 2013-04-08 20:08:47Z pwessel $
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
 * trend1d [<xy[w]file>] -F<output_flags> -N[f]<n_m_parameters>[r] 
 *	[-C<condition_#>] [-I[<confid>]] [-V] [-W]
 *
 * where:
 *	[<xy[w]file>] is an ascii file with x y in first 2 columns [or
 *		x y w in first 3 columns].  Default reads from GMT_stdin.
 *	-F<output_flags> is a string of at least one, up to five, in
 *		and order, from the set {x y m r w}.  x,y = input,
 *		m = model, r = residual = y-m, and w= weight used.
 *	-N[f]<n_m_parameters>[r]
 *		If iterative Robust fitting desired, use append r.
 *		To fit a Fourier model, use -Nf.
 *		Number of terms in the model is <n_m_parameters>.
 *		Example:  Robust quadratic polynomial:  -N2r.
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
 *	[-W]	Weighted data are input.  Read 3 cols and use 3rd as weight.
 *
 *
 * Read GMT_stdin or file of x y pairs, or weighted pairs as x,y w data.  Fit 
 * a regression model y = f(x) + e, where e are error misfits and f(x) has
 * some user-prescribed functional form.  Presently available models are
 * polynomials and Fourier series.  The user may choose the number of terms
 * in the model to fit, whether to seek iterative refinement robust w.r.t.
 * outliers, and whether to seek automatic discovery of the significant
 * number of model parameters.
 *
 *
 * In trend1d I chose to construct the polynomial model using Chebyshev 
 * Polynomials so that the user may easily compare the sizes of the
 * coefficients (and compare with a Fourier series as well).  Tn(x)
 * is an n-degree polynomial with n zero-crossings in [-1,1] and n+1
 * extrema, at which the value of Tn(x) is +/- 1.  It is this property
 * which makes it easy to compare the size of the coefficients.
 *
 * During model fitting the data x coordinate is normalized into the domain
 * [-1, 1] for Chebyshev Polynomial fitting, or into the domain [-pi, pi]
 * for Fourier series fitting.  Before writing out the data the coordinate
 * is rescaled to match the original input values.
 *
 * An n degree polynomial can be written with terms of the form a0 + a1*x
 * + a2*x*x + ...  But it can also be written using other polynomial 
 * basis functions, such as a0*P0 + a1*P1 + a2*P2..., the Legendre
 * polynomials, and a0*T0 + a1*T1 + a2*T2..., the Chebyshev polynomials.
 * (The domain of the x values has to be in [-1, 1] in order to use P or T.)
 * It is well known that the ordinary polynomial basis 1, x, x*x, ... gives 
 * terribly ill- conditioned matrices.  The Ps and Ts do much better.
 * This is because the ordinary basis is far from orthogonal.  The Ps
 * are orthogonal on [-1,1] and the Ts are orthogonal on [-1,1] under a
 * simple weight function.
 * Because the Ps have ordinary orthogonality on [-1,1], I expected them
 * to be the best basis for a regression model; best meaning that they 
 * would lead to the most balanced G'G (matrix of normal equations) with
 * the smallest condition number and the most nearly diagonal model
 * parameter covariance matrix ((G'G)inverse).  It turns out, however, that
 * the G'G obtained from the Ts is very similar and usually has a smaller 
 * condition number than the Ps G'G.  Both of these are vastly superior to
 * the usual polynomials 1, x, x*x.  In a test with 1000 equally spaced
 * data and 8 model parameters, the Chebyshev system had a condition # = 10.6,
 * Legendre = 14.8, and traditional = 54722.7.  For 1000 randomly spaced data
 * and 8 model parameters, the results were C = 13.1, L = 15.6, and P = 54916.6.
 * As the number of model parameters approaches the number of data, the 
 * situation still holds, although all matrices get ill-conditioned; for 8 
 * random data and 8 model parameters, C = 1.8e+05, L = 2.6e+05, P = 1.1e+08.
 * I expected the Legendre polynomials to have a covariance matrix more nearly 
 * diagonal than that of the Chebyshev polynomials, but on this criterion also
 * the Chebyshev turned out to do better.  Only as ndata -> n_model_parameters
 * does the Legendre covariance matrix do better than the Chebyshev.   So for
 * all these reasons I use Chebyshev polynomials.
 *
 * Author:	W. H. F. Smith
 * Date:	25 February 1991-2000.
 * Revised:	11 June, 1991-2000 for v2.0 of GMT-SYSTEM.
 *		13-JUN-1998, for GMT 3.1 (PW)
 *		13-JUL-2000, for GMT 3.3.5 (PW)
 *		10-MAY-2001, PW: Use Numerical Recipes scheme to also output polynomial coefficients
 * Version:	3.4 18-APR-2001
 * Version:	4.1.x
 */

#include "gmt.h"

#define TREND1D_N_OUTPUT_CHOICES 5

struct TREND1D_CTRL {
	struct C {	/* -C<condition_#> */
		GMT_LONG active;
		double value;
	} C;
	struct F {	/* -F<xymrw> */
		GMT_LONG active;
		char col[TREND1D_N_OUTPUT_CHOICES];	/* Character codes for desired output in the right order */
	} F;
	struct I {	/* -I[<confidence>] */
		GMT_LONG active;
		double value;
	} I;
	struct N {	/* -N[f]<n_model>[r] */
		GMT_LONG active;
		GMT_LONG robust;
		GMT_LONG mode;
		GMT_LONG value;
	} N;
	struct W {	/* -W */
		GMT_LONG active;
	} W;
};

#define TREND1D_POLYNOMIAL 0
#define TREND1D_FOURIER 1

struct	TREND1D_DATA {
	double	x;
	double	y;
	double	m;
	double	r;
	double	w;
};

int main (int argc, char **argv)
{
	GMT_LONG i, j, k,n_outputs, n_model, significant, rank, n_req;
	GMT_LONG  n_data, np;
	GMT_LONG	error = FALSE, weighted_output = FALSE;

	double	*gtg = NULL, *v = NULL, *gtd = NULL, *lambda = NULL, *workb = NULL, *workz = NULL, *c_model = NULL, *o_model = NULL, *w_model = NULL, *work = NULL;	/* Arrays  */
	double	xmin, xmax, c_chisq, o_chisq, w_chisq, scale = 1.0, prob;
	double	get_chisq(struct TREND1D_DATA *data, GMT_LONG n_data, GMT_LONG n_model);

	char	format[BUFSIZ];

	FILE	*fp = NULL;

	struct	TREND1D_DATA *data = NULL;
	struct TREND1D_CTRL *Ctrl = NULL;

	void read_data(struct TREND1D_DATA **data, GMT_LONG *n_data, double *xmin, double *xmax, GMT_LONG weighted_input, double **work, FILE *fp);
	void write_output(struct TREND1D_DATA *data, GMT_LONG n_data, char *output_choice, GMT_LONG n_outputs), transform_x(struct TREND1D_DATA *data, GMT_LONG n_data, GMT_LONG model_type, double xmin, double xmax);
	void untransform_x(struct TREND1D_DATA *data, GMT_LONG n_data, GMT_LONG model_type, double xmin, double xmax);
	void recompute_weights(struct TREND1D_DATA *data, GMT_LONG n_data, double *work, double *scale);
	void allocate_array_space(GMT_LONG np, double **gtg, double **v, double **gtd, double **lambda, double **workb, double **workz, double **c_model, double **o_model, double **w_model);
	void free_the_memory(double *gtg, double *v, double *gtd, double *lambda, double *workb, double *workz, double *c_model, double *o_model, double *w_model, struct TREND1D_DATA *data, double *work);
	void calc_m_and_r(struct TREND1D_DATA *data, GMT_LONG n_data, double *model, GMT_LONG n_model, GMT_LONG m_type, double *grow);
	void move_model_a_to_b(double *model_a, double *model_b, GMT_LONG n_model, double *chisq_a, double *chisq_b);
	void load_gtg_and_gtd(struct TREND1D_DATA *data, GMT_LONG n_data, double *gtg, double *gtd, double *grow, GMT_LONG n_model, GMT_LONG mp, GMT_LONG m_type);
	void solve_system(double *gtg, double *gtd, double *model, GMT_LONG n_model, GMT_LONG mp, double *lambda, double *v, double *b, double *z, double c_no, GMT_LONG *ir);
	void GMT_cheb_to_pol (double c[], GMT_LONG n, double a, double b);
	void *New_trend1d_Ctrl (), Free_trend1d_Ctrl (struct TREND1D_CTRL *C);
	
	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct TREND1D_CTRL *)New_trend1d_Ctrl ();	/* Allocate and initialize a new control structure */

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
						if (k < TREND1D_N_OUTPUT_CHOICES)
							Ctrl->F.col[k] = argv[i][j];
						else {
							error++;
							fprintf (stderr, "%s: GMT SYNTAX ERROR -F option: Too many output columns selected\n", GMT_program);
							fprintf (stderr, "%s: GMT SYNTAX ERROR -F option: Choose from -Fxymrw\n", GMT_program);
						}
					}
					break;
				case 'I':
					Ctrl->I.active = TRUE;
					Ctrl->I.value = (argv[i][2]) ? atof(&argv[i][2]) : 0.51;
					break;
				case 'N':
					Ctrl->N.active = TRUE;
					if (strchr (argv[i], 'r')) Ctrl->N.robust = TRUE;
					j = (argv[i][2] == 'r') ? 3 : 2;
					if (argv[i][j] == 'F' || argv[i][j] == 'f') {
						Ctrl->N.mode = TREND1D_FOURIER;
						j++;
					}
					else if (argv[i][j] == 'P' || argv[i][j] == 'p') {
						Ctrl->N.mode = TREND1D_POLYNOMIAL;
						j++;
					}
					if (argv[i][j])
						Ctrl->N.value = atoi(&argv[i][j]);
					else {
						error = TRUE;
 						fprintf (stderr, "%s: GMT SYNTAX ERROR -N option.  No model specified\n", GMT_program);
					}
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
		fprintf(stderr,"trend1d %s - Fit a [weighted] [robust] polynomial [or Fourier] model for y = f(x) to ascii xy[w]\n\n", GMT_VERSION);
		fprintf(stderr,"usage:  trend1d -F<xymrw> -N[f]<n_model>[r] [<xy[w]file>] [-C<condition_#>]\n");
		fprintf(stderr,"\t[%s] [-I[<confidence>]] [-V] [-W] [%s] [%s] [%s]\n\n", GMT_H_OPT, GMT_t_OPT, GMT_b_OPT, GMT_f_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf(stderr,"\t-F Choose at least 1, up to 5, any order, of xymrw for ascii output to stdout.\n");
		fprintf(stderr,"\t   x=x, y=y, m=model, r=residual=y-m, w=weight.  w determined iteratively if robust fit used.\n");
		fprintf(stderr,"\t-N fit a Polynomial [Default] or Fourier (-Nf) model with <n_model> terms.\n");
		fprintf(stderr,"\t   Append r for robust model.  E.g., robust quadratic = -N3r.\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf(stderr,"\t[<xy[w]file>] name of ascii file, first 2 cols = x y [3 cols = x y w]; [Default reads stdin].\n");
		fprintf(stderr,"\t-C Truncate eigenvalue spectrum so matrix has <condition_#> [Default = 1.0e06].\n");
		GMT_explain_option ('H');
		fprintf(stderr,"\t-I Iteratively Increase # model parameters, to a max of <n_model> so long as the\n");
		fprintf(stderr,"\t   reduction in variance is significant at the <confidence> level.\n");
		fprintf(stderr,"\t   Give -I without a number to default to 0.51 confidence level.\n");
		GMT_explain_option ('V');
		fprintf(stderr,"\t-W Weighted input given, weights in 3rd column.  [Default is unweighted].\n");
		GMT_explain_option (':');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf(stderr,"\t   Default is 2 (or 3 if -W is set) input columns.\n");
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
	if (Ctrl->N.value <= 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -N option.  A positive number of terms must be specified\n", GMT_program);
		error++;
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
	n_req = (Ctrl->W.active) ? 3 : 2;
	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = n_req;
	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] < n_req) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data (-bi) must have at least %ld columns\n", GMT_program, n_req);
		error++;
	}
	for (k = n_outputs = 0; k < TREND1D_N_OUTPUT_CHOICES && Ctrl->F.col[k]; k++) {
		if (!strchr ("xymrw", Ctrl->F.col[k])) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR -F option.  Unrecognized output choice %c\n", GMT_program, Ctrl->F.col[k]);
			error++;
		}
		else if (Ctrl->F.col[k] == 'w')
			weighted_output = TRUE;

		n_outputs++;
	}
        if (n_outputs == 0) {
                fprintf (stderr, "%s: GMT SYNTAX ERROR -F option.  Must specify at least one output columns \n", GMT_program);
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

	np = Ctrl->N.value;	/* Row dimension for matrices gtg and v  */
	allocate_array_space(np, &gtg, &v, &gtd, &lambda, &workb, &workz, &c_model, &o_model, &w_model);

	read_data(&data, &n_data, &xmin, &xmax, Ctrl->W.active, &work, fp);

	if (xmin == xmax) {
		fprintf(stderr,"%s:  Fatal error in input data.  X min = X max.\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (n_data == 0) {
		fprintf(stderr,"%s:  Fatal error.  Could not read any data.\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (n_data < Ctrl->N.value) fprintf(stderr,"%s: Warning. Ill-posed problem.  n_data < n_model_max.\n", GMT_program);

	transform_x(data, n_data, Ctrl->N.mode, xmin, xmax);	/* Set domain to [-1, 1] or [-pi, pi]  */

	if (gmtdefs.verbose) {
		sprintf(format,"%%s:  Read %%ld data with X values from %s to %s\n", gmtdefs.d_format, gmtdefs.d_format);
		fprintf(stderr, format, GMT_program, n_data, xmin, xmax);
		fprintf(stderr,"N_model\tRank\tChi_Squared\tSignificance\n");
	}

	sprintf (format, "%%ld\t%%ld\t%s\t%s\n", gmtdefs.d_format, gmtdefs.d_format);

	if (Ctrl->I.active) {
		n_model = 1;

		/* Fit first model  */
		load_gtg_and_gtd(data, n_data, gtg, gtd, workb, n_model, np, Ctrl->N.mode);
		solve_system(gtg, gtd, c_model, n_model, np, lambda, v, workb, workz, Ctrl->C.value, &rank);
		calc_m_and_r(data, n_data, c_model, n_model, Ctrl->N.mode, workb);
		c_chisq = get_chisq(data, n_data, n_model);
		if (gmtdefs.verbose) fprintf(stderr, format, n_model, rank, c_chisq, 1.0);
		if (Ctrl->N.robust) {
			do {
				recompute_weights(data, n_data, work, &scale);
				move_model_a_to_b(c_model, w_model, n_model, &c_chisq, &w_chisq);
				load_gtg_and_gtd(data, n_data, gtg, gtd, workb, n_model, np, Ctrl->N.mode);
				solve_system(gtg, gtd, c_model, n_model, np, lambda, v, workb, workz, Ctrl->C.value, &rank);
				calc_m_and_r(data, n_data, c_model, n_model, Ctrl->N.mode, workb);
				c_chisq = get_chisq(data, n_data, n_model);
				significant = GMT_sig_f(c_chisq, n_data-n_model, w_chisq, n_data-n_model, Ctrl->I.value, &prob);
				if (gmtdefs.verbose) fprintf(stderr, format, n_model, rank, c_chisq, prob);
			} while (significant);
			/* Go back to previous model only if w_chisq < c_chisq  */
			if (w_chisq < c_chisq) {
				move_model_a_to_b(w_model, c_model, n_model, &w_chisq, &c_chisq);
				calc_m_and_r(data, n_data, c_model, n_model, Ctrl->N.mode, workb);
				if (weighted_output && n_model == Ctrl->N.value) recompute_weights(data, n_data, work, &scale);
			}
		}
		/* First [robust] model has been found  */

		significant = TRUE;
		while(n_model < Ctrl->N.value && significant) {
			move_model_a_to_b(c_model, o_model, n_model, &c_chisq, &o_chisq);
			n_model++;

			/* Fit next model  */
			load_gtg_and_gtd(data, n_data, gtg, gtd, workb, n_model, np, Ctrl->N.mode);
			solve_system(gtg, gtd, c_model, n_model, np, lambda, v, workb, workz, Ctrl->C.value, &rank);
			calc_m_and_r(data, n_data, c_model, n_model, Ctrl->N.mode, workb);
			c_chisq = get_chisq(data, n_data, n_model);
			if (gmtdefs.verbose) fprintf(stderr, format, n_model, rank, c_chisq, 1.00);
			if (Ctrl->N.robust) {
				do {
					recompute_weights(data, n_data, work, &scale);
					move_model_a_to_b(c_model, w_model, n_model, &c_chisq, &w_chisq);
					load_gtg_and_gtd(data, n_data, gtg, gtd, workb, n_model, np, Ctrl->N.mode);
					solve_system(gtg, gtd, c_model, n_model, np, lambda, v, workb, workz, Ctrl->C.value, &rank);
					calc_m_and_r(data, n_data, c_model, n_model, Ctrl->N.mode, workb);
					c_chisq = get_chisq(data, n_data, n_model);
					significant = GMT_sig_f(c_chisq, n_data-n_model, w_chisq, n_data-n_model, Ctrl->I.value, &prob);
					if (gmtdefs.verbose) fprintf(stderr, format, n_model, rank, c_chisq, prob);
				} while (significant);
				/* Go back to previous model only if w_chisq < c_chisq  */
				if (w_chisq < c_chisq) {
					move_model_a_to_b(w_model, c_model, n_model, &w_chisq, &c_chisq);
					calc_m_and_r(data, n_data, c_model, n_model, Ctrl->N.mode, workb);
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
			calc_m_and_r(data, n_data, c_model, n_model, Ctrl->N.mode, workb);
			if (Ctrl->N.robust && weighted_output) recompute_weights(data, n_data, work, &scale);
		}
	}
	else {
		n_model = Ctrl->N.value;
		load_gtg_and_gtd(data, n_data, gtg, gtd, workb, n_model, np, Ctrl->N.mode);
		solve_system(gtg, gtd, c_model, n_model, np, lambda, v, workb, workz, Ctrl->C.value, &rank);
		calc_m_and_r(data, n_data, c_model, n_model, Ctrl->N.mode, workb);
		c_chisq = get_chisq(data, n_data, n_model);
		if (gmtdefs.verbose) fprintf(stderr, format, n_model, rank, c_chisq, 1.00);
		if (Ctrl->N.robust) {
			do {
				recompute_weights(data, n_data, work, &scale);
				move_model_a_to_b(c_model, w_model, n_model, &c_chisq, &w_chisq);
				load_gtg_and_gtd(data, n_data, gtg, gtd, workb, n_model, np, Ctrl->N.mode);
				solve_system(gtg, gtd, c_model, n_model, np, lambda, v, workb, workz, Ctrl->C.value, &rank);
				calc_m_and_r(data, n_data, c_model, n_model, Ctrl->N.mode, workb);
				c_chisq = get_chisq(data, n_data, n_model);
				significant = GMT_sig_f(c_chisq, n_data-n_model, w_chisq, n_data-n_model, Ctrl->I.value, &prob);
				if (gmtdefs.verbose) fprintf(stderr, format, n_model, rank, c_chisq, prob);
			} while (significant);
			/* Go back to previous model only if w_chisq < c_chisq  */
			if (w_chisq < c_chisq) {
				move_model_a_to_b(w_model, c_model, n_model, &w_chisq, &c_chisq);
				calc_m_and_r(data, n_data, c_model, n_model, Ctrl->N.mode, workb);
				if (weighted_output && n_model == Ctrl->N.value) recompute_weights(data, n_data, work, &scale);
			}
		}
	}

	if (gmtdefs.verbose) {
		sprintf (format, "%%s: Final model stats:  N model parameters %%ld.  Rank %%ld.  Chi-Squared:  %s\n", gmtdefs.d_format);
		fprintf(stderr, format, GMT_program, n_model, rank, c_chisq);
		fprintf(stderr,"%s: Model Coefficients  (Chebyshev): ", GMT_program);
		sprintf (format, "\t%s", gmtdefs.d_format);
		for (i = 0; i < n_model; i++) fprintf (stderr, format, c_model[i]);
		fprintf(stderr,"\n");
		GMT_cheb_to_pol (c_model, n_model, xmin, xmax);
		fprintf(stderr,"%s: Model Coefficients (Polynomial): ", GMT_program);
		for (i = 0; i < n_model; i++) fprintf (stderr, format, c_model[i]);
		fprintf(stderr,"\n");
	}

	untransform_x(data, n_data, Ctrl->N.mode, xmin, xmax);

	write_output(data, n_data, Ctrl->F.col, n_outputs);

	free_the_memory(gtg, v, gtd, lambda, workb, workz, c_model, o_model, w_model, data, work);

	Free_trend1d_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void read_data (struct TREND1D_DATA **data, GMT_LONG *n_data, double *xmin, double *xmax, GMT_LONG weighted_input, double **work, FILE *fp)
{
	GMT_LONG	n_alloc = GMT_CHUNK, n_expected_fields, n_fields;
	GMT_LONG i, n_read;
	double	*in;
	char	buffer[BUFSIZ];

	if (fp == NULL) {
		fp = GMT_stdin;
#ifdef SET_IO_MODE
		GMT_setmode (GMT_IN);
#endif
	}
	(*data) = (struct TREND1D_DATA *) GMT_memory (VNULL, (size_t)n_alloc, sizeof(struct TREND1D_DATA), GMT_program);

	if (GMT_io.io_header[GMT_IN]) for (i = 0; i < GMT_io.n_header_recs; i++) GMT_fgets (buffer, BUFSIZ, fp);
	i = n_read = 0;
	n_expected_fields = (GMT_io.binary[GMT_IN]) ? GMT_io.ncol[GMT_IN] : 2 + weighted_input;

	while ((n_fields = GMT_input (fp, &n_expected_fields, &in)) >= 0 && !(GMT_io.status & GMT_IO_EOF)) {

		n_read++;
		if (GMT_io.status & GMT_IO_MISMATCH) {
			fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld (skipped)\n", GMT_program, n_fields, n_expected_fields, n_read);
			continue;
		}

		(*data)[i].x = in[GMT_X];
		(*data)[i].y = in[GMT_Y];
		(*data)[i].w = (weighted_input) ? in[GMT_Z] : 1.0;

		if (i) {
			if (*xmin > (*data)[i].x) *xmin = (*data)[i].x;
			if (*xmax < (*data)[i].x) *xmax = (*data)[i].x;
		}
		else {
			*xmin = (*data)[i].x;
			*xmax = (*data)[i].x;
		}
		i++;

		if (i == n_alloc) {
			n_alloc <<= 1;
			*data = (struct TREND1D_DATA *) GMT_memory ((void *)*data, (size_t)n_alloc, sizeof(struct TREND1D_DATA), GMT_program);
		}
		if (i == INT_MAX) {
			fprintf (stderr, "%s: ERROR: Cannot process more than %d data points\n", GMT_program, INT_MAX);
			GMT_free ((void *)data);
			exit (EXIT_FAILURE);
		}
	}
	if (fp != GMT_stdin) GMT_fclose(fp);

	*data = (struct TREND1D_DATA *) GMT_memory ((void *)*data, (size_t)i, sizeof(struct TREND1D_DATA), GMT_program);
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

void write_output (struct TREND1D_DATA *data, GMT_LONG n_data, char *output_choice, GMT_LONG n_outputs)
{
	GMT_LONG	i;
	GMT_LONG j;
	double out[5];

	for (i = 0; i < n_data; i++) {
		for (j = 0; j < n_outputs; j++) {
			switch (output_choice[j]) {
				case 'x':
					out[j] = data[i].x;
					break;
				case 'y':
					out[j] = data[i].y;
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

void free_the_memory (double *gtg, double *v, double *gtd, double *lambda, double *workb, double *workz, double *c_model, double *o_model, double *w_model, struct TREND1D_DATA *data, double *work)
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

void transform_x (struct TREND1D_DATA *data, GMT_LONG n_data, GMT_LONG model_type, double xmin, double xmax)
{
	GMT_LONG	i;
	double	offset, scale;

	offset = 0.5 * (xmin + xmax);	/* Mid Range  */
	scale = 2.0 / (xmax - xmin);	/* 1 / (1/2 Range)  */

	if (model_type == TREND1D_FOURIER) {	/* Set Range to 1 period  */
		scale *= M_PI;
	}

	for (i = 0; i < n_data; i++) {
		data[i].x = (data[i].x - offset) * scale;
	}
}

void untransform_x (struct TREND1D_DATA *data, GMT_LONG n_data, GMT_LONG model_type, double xmin, double xmax)
{
	GMT_LONG	i;
	double	offset, scale;

	offset = 0.5 * (xmin + xmax);	/* Mid Range  */
	scale = 0.5 * (xmax - xmin);	/* 1/2 Range  */

	if (model_type == TREND1D_FOURIER) {
		scale /= M_PI;
	}

	for (i = 0; i < n_data; i++) {
		data[i].x = (data[i].x * scale) + offset;
	}
}

double get_chisq (struct TREND1D_DATA *data, GMT_LONG n_data, GMT_LONG n_model)
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

void recompute_weights (struct TREND1D_DATA *data, GMT_LONG n_data, double *work, double *scale)
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

void load_g_row (double x, GMT_LONG n, double *gr, GMT_LONG m)
      	  	/* Current data position, appropriately normalized.  */
   	  	/* Number of model parameters, and elements of gr[]  */
      	     	/* Elements of row of G matrix.  */
   	  	/* Parameter indicating model type  */
{
	/* Routine computes the elements gr[j] in the ith row of the
		G matrix (Menke notation), where x is the ith datum's
		abscissa.  */

	GMT_LONG	j, k;

	if (n) {

		gr[0] = 1.0;

		switch (m) {

			case TREND1D_POLYNOMIAL:
				/* Create Chebyshev polynomials  */
				if (n > 1) gr[1] = x;
				for (j = 2; j < n; j++) {
					gr[j] = 2 * x * gr[j-1] - gr[j-2];
				}
				break;

			case TREND1D_FOURIER:
				for (j = 1; j < n; j++) {
					k = (j + 1)/2;
					if (k > 1) {
						if (j%2) {
							gr[j] = cos(k*x);
						}
						else {
							gr[j] = sin(k*x);
						}
					}
					else {
						if (j%2) {
							gr[j] = cos(x);
						}
						else {
							gr[j] = sin(x);
						}
					}
				}
				break;
		}
	}
}

void calc_m_and_r (struct TREND1D_DATA *data, GMT_LONG n_data, double *model, GMT_LONG n_model, GMT_LONG m_type, double *grow)
{
	/*	model[n_model] holds solved coefficients of m_type model.
		grow[n_model] is a vector for a row of G matrix.  */

	GMT_LONG	i;
	GMT_LONG j;
	for (i = 0; i < n_data; i++) {
		load_g_row(data[i].x, n_model, grow, m_type);
		data[i].m = 0.0;
		for (j = 0; j < n_model; j++) {
			data[i].m += model[j]*grow[j];
		}
		data[i].r = data[i].y - data[i].m;
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

void load_gtg_and_gtd (struct TREND1D_DATA *data, GMT_LONG n_data, double *gtg, double *gtd, double *grow, GMT_LONG n_model, GMT_LONG mp, GMT_LONG m_type)
{
   	/* mp is row dimension of gtg  */

	GMT_LONG	i;
	GMT_LONG j, k;
	double	wy;

	/* First zero the contents for summing:  */

	for (j = 0; j < n_model; j++) {
		for (k = 0; k < n_model; k++) {
			gtg[j + k*mp] = 0.0;
		}
		gtd[j] = 0.0;
	}

	/* Sum over all data  */
	for (i = 0; i < n_data; i++) {
		load_g_row(data[i].x, n_model, grow, m_type);
		if (data[i].w != 1.0) {
			wy = data[i].w * data[i].y;
			for (j = 0; j < n_model; j++) {
				for (k = 0; k < n_model; k++) {
					gtg[j + k*mp] += (data[i].w * grow[j] * grow[k]);
				}
				gtd[j] += (wy * grow[j]);
			}
		}
		else {
			for (j = 0; j < n_model; j++) {
				for (k = 0; k < n_model; k++) {
					gtg[j + k*mp] += (grow[j] * grow[k]);
				}
				gtd[j] += (data[i].y * grow[j]);
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
		n = (GMT_LONG)n_model;
		m = (GMT_LONG)mp;
		if(GMT_jacobi(gtg, &n, &m, lambda, v, b, z, &nrots)) {
			fprintf(stderr,"%s:  Warning:  Matrix Solver Convergence Failure.\n", GMT_program);
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

void GMT_cheb_to_pol (double c[], GMT_LONG n, double a, double b)
{
	/* Convert from Chebyshev coefficients used on a t =  [-1,+1] interval
	 * to polynomial coefficients on the original x = [a b] interval.
	 * Modified from Numerical Miracles, ...eh Recipes */
	 
	 GMT_LONG j, k;
	 double sv, cnst, fac;
	 double *d, *dd;
	 
	 d  = GMT_memory (VNULL, (size_t)n, sizeof (double), GMT_program);
	 dd = GMT_memory (VNULL, (size_t)n, sizeof (double), GMT_program);
	 
	 /* First we generate coefficients for a polynomial in t */
	 
	 d[0] = c[n-1];
	 for (j = n - 2; j >= 1; j--) {
	 	for (k = n - j; k >= 1; k--) {
			sv = d[k];
			d[k] = 2.0 * d[k-1] - dd[k];
			dd[k] = sv;
		}
		sv = d[0];
		d[0] = -dd[0] + c[j];
		dd[0] = sv;
	}
	for (j = n - 1; j >= 1; j--) d[j] = d[j-1] - dd[j];
	/* d[0] = -dd[0] + 0.5 * c[0]; */	/* This is what Num. Rec. says, but we do not do the approx with 0.5 * c[0] */
	d[0] = -dd[0] + c[0];

	/* Next step is to undo the scaling so we can use coefficients with x */

	cnst = fac = 2.0 / (b - a);
	for (j = 1; j < n; j++) {
		d[j] *= fac;
		fac *= cnst;
	}
	cnst = 0.5 * (a + b);
	for (j = 0; j <= n - 2; j++) for (k = n - 2; k >= j; k--) d[k] -= cnst * d[k+1];

	/* Return the new coefficients via c */

	memcpy ((void *)c, (void *)d, (size_t)(n * sizeof (double)));

	GMT_free ((void *) d);
	GMT_free ((void *) dd);
}

void *New_trend1d_Ctrl () {	/* Allocate and initialize a new control structure */
	struct TREND1D_CTRL *C;
	
	C = (struct TREND1D_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct TREND1D_CTRL), "New_trend1d_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	C->C.value = 1.0e06;		/* Condition number for matrix solution  */	
	C->I.value = 0.51;		/* Confidence interval for significance test  */
	C->N.mode = TREND1D_POLYNOMIAL;
	return ((void *)C);
}

void Free_trend1d_Ctrl (struct TREND1D_CTRL *C) {	/* Deallocate control structure */
	GMT_free ((void *)C);	
}
