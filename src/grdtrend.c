/*--------------------------------------------------------------------
 *	$Id: grdtrend.c 10173 2014-01-01 09:52:34Z pwessel $
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
/* grdtrend <input.grd> -N[r]<n_model> [-T<trend.grd>] [-V]
	[-W<weight.grd] [-D<differences.grd]

Reads a grid file and fits a trend surface.  Trend surface
is defined by:

m1 +m2*x + m3*y + m4*xy + m5*x*x + m6*y*y + m7*x*x*x
	+ m8*x*x*y + m9*x*y*y + m10*y*y*y.

n_model is set by the user to be an integer in [1,10]
which sets the number of model coefficients to fit.
Thus:
n_model = 1 gives the mean value of the surface,
n_model = 3 fits a plane,
n_model = 4 fits a bilinear surface,
n_model = 6 fits a biquadratic,
n_model = 10 fits a bicubic surface.

The user may write out grid files of the fitted surface
[-T<trend.grd>] and / or of the residuals (input data
minus fitted trend) [-D<differences.grd] and / or of
the weights used in iterative fitting [-W<weight.grd].
This last option applies only when the surface is fit
iteratively [-N<n>[r]].

A robust fit may be achieved by iterative fitting of
a weighted least squares problem, where the weights
are set according to a scale length based on the 
Median absolute deviation (MAD: Huber, 1982).  The
-N<n>r option achieves this.

Author:		W. H. F. Smith
Date:		21 May, 1991.
Version:	4
Calls:		uses the QR solution of the Normal
		equations furnished by Wm. Menke's
		C routine "gauss".  We gratefully
		acknowledge this contribution, now
		as GMT_gauss in gmt_vector.c
Revised:	12-JUN-1998 PW, for GMT 3.1

Remarks:

We adopt a translation and scaling of the x,y coordinates.
We choose x,y such that they are in [-1,1] over the range
of the grid file.  If the problem is unweighted, all input
values are filled (no "holes" or NaNs in the input grid file),
and n_model <= 4 (bilinear or simpler), then the normal
equations matrix (G'G in Menke notation) is diagonal under
this change of coordinates, and the solution is trivial.
In this case, it would be dangerous to try to accumulate
the sums which are the elements of the normal equations;
while they analytically cancel to zero, the addition errors
would likely prevent this.  Therefore we have written a
routine, grd_trivial_model(), to handle this case.

If the problem is more complex than the above trivial case,
(missing values, weighted problem, or n_model > 4), then
G'G is not trivial and we just naively accumulate sums in
the G'G matrix.  We hope that the changed coordinates in
[-1,1] will help the accuracy of the problem.  We also use
Legendre polynomials in this case so that the matrix elements
are conveniently sized near 0 or 1.

*/

#include "gmt.h"

struct GRDTREND_CTRL {	/* All control options for this program (except common args) */
	/* ctive is TRUE if the option has been activated */
	struct D {	/* -D<diffgrid> */
		GMT_LONG active;
		char *file;
	} D;
	struct N {	/* -N[r]<n_model> */
		GMT_LONG active;
		GMT_LONG robust;
		GMT_LONG value;
	} N;
	struct T {	/* -T<trend.grd> */
		GMT_LONG active;
		char *file;
	} T;
	struct W {	/* -W<weight.grd> */
		GMT_LONG active;
		char *file;
	} W;
};

int main(int argc, char **argv)
{
	GMT_LONG	error = FALSE, trivial, weighted;

	char	*i_filename = NULL, format[BUFSIZ];

	GMT_LONG	j, k, ierror = 0, iterations;
	GMT_LONG	i, nxy;

	float	*data = NULL;		/* Pointer for array from input grid file  */
	float	*trend = NULL;		/* Pointer for array containing fitted surface  */
	float	*resid = NULL;		/* Pointer for array containing residual surface  */
	float	*weight = VNULL;/* Pointer for array containing data weights  */

	double	chisq, old_chisq, zero_test = 1.0e-08, scale = 1.0;
	double	*xval = NULL;		/* Pointer for array of change of variable:  x[i]  */
	double	*yval = NULL;		/* Pointer for array of change of variable:  y[j]  */
	double	*gtg = NULL;		/* Pointer for array for matrix G'G normal equations  */
	double	*gtd = NULL;		/* Pointer for array for vector G'd normal equations  */
	double	*old = NULL;		/* Pointer for array for old model, used for robust sol'n  */
	double	*pstuff = NULL;	/* Pointer for array for Legendre polynomials of x[i],y[j]  */

	struct GRD_HEADER head_d, head_w;
	struct GRDTREND_CTRL *Ctrl = NULL;

	void	set_up_vals(double *val, GMT_LONG nval, double vmin, double vmax, double dv, GMT_LONG pixel_reg);		/* Store x[i], y[j] once for all to save time  */
	void	grd_trivial_model(float *data, GMT_LONG nx, GMT_LONG ny, double *xval, double *yval, double *gtd, GMT_LONG n_model);	/* Fit trivial models.  See Remarks above.  */
	void	load_pstuff(double *pstuff, GMT_LONG n_model, double x, double y, GMT_LONG newx, GMT_LONG newy);		/* Compute Legendre polynomials of x[i],y[j] as needed  */
	void	compute_trend(float *trend, GMT_LONG nx, GMT_LONG ny, double *xval, double *yval, double *gtd, GMT_LONG n_model, double *pstuff);	/* Find trend from a model  */
	void	compute_resid(float *data, float *trend, float *resid, GMT_LONG nxy);	/* Find residuals from a trend  */
	void	compute_chisq(float *resid, float *weight, GMT_LONG nxy, double *chisq, double scale);	/* Find Chi-Squared from weighted residuals  */
	void	compute_robust_weight(float *resid, float *weight, GMT_LONG nxy, double *scale);	/* Find weights from residuals  */
	void	write_model_parameters(double *gtd, GMT_LONG n_model);	/* Do reports if gmtdefs.verbose == TRUE  */
	void	load_gtg_and_gtd(float *data, GMT_LONG nx, GMT_LONG ny, double *xval, double *yval, double *pstuff, double *gtg, double *gtd, GMT_LONG n_model, float *weight, GMT_LONG weighted);		/* Fill normal equations matrices  */
	void	*New_grdtrend_Ctrl (), Free_grdtrend_Ctrl (struct GRDTREND_CTRL *C);

/* Execution begins here with loop over arguments:  */

	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct GRDTREND_CTRL *) New_grdtrend_Ctrl ();		/* Allocate and initialize defaults in a new control structure */

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {

				/* Common parameters */

				case 'V':
				case '\0':
					error += GMT_parse_common_options (argv[i], 0, 0, 0, 0);
					break;

				/* Supplemental parameters */

				case 'D':
					Ctrl->D.active = TRUE;
					if (argv[i][2])
						Ctrl->D.file = strdup (&argv[i][2]);
					else {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -D option:  Must specify file name\n", GMT_program);
						error = TRUE;
					}
					break;
				case 'N':
					/* Must check for both -N[r]<n_model> and -N<n_model>[r] due to confusion */
					Ctrl->N.active = TRUE;
					if (strchr (argv[i], 'r')) Ctrl->N.robust = TRUE;
					j = (argv[i][2] == 'r') ? 3 : 2;
					if (argv[i][j]) Ctrl->N.value = atoi(&argv[i][j]);
					break;
				case 'T':
					Ctrl->T.active = TRUE;
					if (argv[i][2])
						Ctrl->T.file = strdup (&argv[i][2]);
					else {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -T option:  Must specify file name\n", GMT_program);
						error = TRUE;
					}
					break;
				case 'W':
					Ctrl->W.active = TRUE;
					if (argv[i][2])
						Ctrl->W.file = strdup (&argv[i][2]);
					else {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -W option:  Must specify file name\n", GMT_program);
						error = TRUE;
					}
					/* OK if this file doesn't exist:  */
					break;
				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else
			i_filename = argv[i];
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "grdtrend %s - Fit trend surface to gridded data\n\n", GMT_VERSION);
		fprintf (stderr,"usage:  grdtrend <input.grd> -N<n_model>[r] [-D<diff.grd>]\n");
		fprintf (stderr,"\t[-T<trend.grd>] [-V] [-W<weight.grd>]\n\n");

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr,"\t<input.grd> is name of grid file to fit trend to.\n");
		fprintf(stderr,"\t-N fit a [robust] model with <n_model> terms.  <n_model> in [1,10].  E.g., robust planar = -N3r.\n");
		fprintf(stderr,"\t   Model parameters order is given as follows:\n");
		fprintf(stderr,"\t   z = m1 + m2*x + m3*y + m4*x*y + m5*x^2 + m6*y^2 + m7*x^3 + m8*x^2*y + m9*x*y^2 + m10*y^3.\n");
		fprintf (stderr,"\n\tOPTIONS:\n");
		fprintf (stderr,"\t-D Supply filename to write grid file of differences (input - trend).\n");
		fprintf (stderr,"\t-T Supply filename to write grid file of trend.\n");
		GMT_explain_option ('V');
		fprintf (stderr,"\t-W Supply filename if you want to [read and] write grid file of weights.\n");
		fprintf (stderr,"\t   If <weight.grd> can be read at run, and if robust = FALSE, weighted problem will be solved.\n");
		fprintf (stderr,"\t   If robust = TRUE, weights used for robust fit will be written to <weight.grd>.\n");
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	if (!i_filename) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify input file\n", GMT_program);
		error++;
	}
	if (Ctrl->N.value <= 0 || Ctrl->N.value > 10) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -N option:  Specify 1-10 model parameters\n", GMT_program);
		error++;
	}

	if (error) exit (EXIT_FAILURE);

/* End of argument parsing.  */

	weighted = (Ctrl->N.robust || Ctrl->W.active);
	trivial = (Ctrl->N.value < 5 && !weighted);

/* Read the input file:  */

	GMT_err_fail (GMT_read_grd_info (i_filename, &head_d), i_filename);

	GMT_grd_init (&head_d, argc, argv, TRUE);
	nxy = GMT_get_nm (head_d.nx, head_d.ny);
	data = (float *) GMT_memory (VNULL, (size_t)nxy, sizeof (float), GMT_program);
	GMT_err_fail (GMT_read_grd (i_filename, &head_d, data, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), i_filename);

	/* Check for NaNs:  */
	i = 0;
	while (trivial && i < nxy) {
		if (GMT_is_fnan (data[i])) trivial = FALSE;
		i++;
	}

/* End input read section.  */

/* Allocate other required arrays:  */

	trend = (float *) GMT_memory (VNULL, (size_t)nxy, sizeof (float), GMT_program);
	resid = (float *) GMT_memory (VNULL, (size_t)nxy, sizeof (float), GMT_program);
	xval = (double *) GMT_memory (VNULL, (size_t)head_d.nx, sizeof (double), GMT_program);
	yval = (double *) GMT_memory (VNULL, (size_t)head_d.ny, sizeof (double), GMT_program);
	gtg = (double *) GMT_memory (VNULL, (size_t)(Ctrl->N.value*Ctrl->N.value), sizeof (double), GMT_program);
	gtd = (double *) GMT_memory (VNULL, (size_t)Ctrl->N.value, sizeof (double), GMT_program);
	old = (double *) GMT_memory (VNULL, (size_t)Ctrl->N.value, sizeof (double), GMT_program);
	pstuff = (double *) GMT_memory (VNULL, (size_t)Ctrl->N.value, sizeof (double), GMT_program);
	pstuff[0] = 1.0; /* This is P0(x) = 1, which is not altered in this program. */

/* If a weight array is needed, get one:  */

	if (weighted) {
		weight = (float *) GMT_memory (VNULL, (size_t)nxy, sizeof (float), GMT_program);
		if (!GMT_access (Ctrl->W.file, R_OK)) {	/* We have weights on input  */
			GMT_grd_init (&head_w, argc, argv, FALSE);
			GMT_err_fail (GMT_read_grd_info (Ctrl->W.file, &head_w), Ctrl->W.file);
			if (head_w.nx != head_d.nx || head_w.ny != head_d.ny) {
				fprintf (stderr,"%s:  Input weight file does not match input data file.  Ignoring.\n", GMT_program);
				for (i = 0; i < nxy; i++) weight[i] = 1.0;
			}
			else
				GMT_err_fail (GMT_read_grd (Ctrl->W.file, &head_w, weight, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), Ctrl->W.file);
		}
		else
			for (i = 0; i < nxy; i++) weight[i] = 1.0;
	}

/* End of weight set up.  */

/* Set up xval and yval lookup tables:  */

	set_up_vals(xval, head_d.nx, head_d.x_min, head_d.x_max, head_d.x_inc,
		 head_d.node_offset);
	set_up_vals(yval, head_d.ny, head_d.y_min, head_d.y_max, head_d.y_inc,
		 head_d.node_offset);

/* End of set up of lookup values.  */

/* Do the problem:  */

	if (trivial) {
		grd_trivial_model(data, head_d.nx, head_d.ny, xval, yval, gtd, Ctrl->N.value);
		compute_trend(trend, head_d.nx, head_d.ny, xval, yval, gtd, Ctrl->N.value, pstuff);
		compute_resid(data, trend, resid, nxy);
	}
	else {	/* Problem is not trivial  !!  */

		load_gtg_and_gtd(data, head_d.nx, head_d.ny, xval, yval, pstuff, gtg, gtd, Ctrl->N.value, weight, weighted);
		GMT_gauss (gtg, gtd, Ctrl->N.value, Ctrl->N.value, zero_test, &ierror, 1);
		if (ierror) {
			fprintf (stderr,"%s:  Gauss returns error code %ld\n", GMT_program, ierror);
			exit (EXIT_FAILURE);
		}
		compute_trend(trend, head_d.nx, head_d.ny, xval, yval, gtd, Ctrl->N.value, pstuff);
		compute_resid(data, trend, resid, nxy);

		if (Ctrl->N.robust) {
			compute_chisq(resid, weight, nxy, &chisq, scale);
			iterations = 1;
			sprintf(format, "%%s Robust iteration %%d:  Old Chi Squared:  %s  New Chi Squared %s\n", gmtdefs.d_format, gmtdefs.d_format);
			do {
				old_chisq = chisq;
				for (k = 0; k < Ctrl->N.value; k++) old[k] = gtd[k];
				compute_robust_weight(resid, weight, nxy, &scale);
				load_gtg_and_gtd(data, head_d.nx, head_d.ny, xval, yval, pstuff, gtg, gtd, Ctrl->N.value, weight, weighted);
				GMT_gauss (gtg, gtd, Ctrl->N.value, Ctrl->N.value, zero_test, &ierror, 1);
				if (ierror) {
					fprintf (stderr,"%s:  Gauss returns error code %ld\n", GMT_program, ierror);
					exit (EXIT_FAILURE);
				}
				compute_trend(trend, head_d.nx, head_d.ny, xval, yval, gtd, Ctrl->N.value, pstuff);
				compute_resid(data, trend, resid, nxy);
				compute_chisq(resid, weight, nxy, &chisq, scale);
				if (gmtdefs.verbose) fprintf (stderr, format, GMT_program, iterations, old_chisq, chisq);
				iterations++;
			} while (old_chisq / chisq > 1.0001);

			/* Get here when new model not significantly better; use old one:  */

			for (k = 0; k < Ctrl->N.value; k++) gtd[k] = old[k];
			compute_trend(trend, head_d.nx, head_d.ny, xval, yval, gtd, Ctrl->N.value, pstuff);
			compute_resid(data, trend, resid, nxy);
		}
	}

/* End of do the problem section.  */

/* Get here when ready to do output:  */

	if (gmtdefs.verbose) write_model_parameters(gtd, Ctrl->N.value);
	if (Ctrl->T.file) {
		strcpy (head_d.title, "trend surface");
		GMT_err_fail (GMT_write_grd (Ctrl->T.file, &head_d, trend, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), Ctrl->T.file);
	}
	if (Ctrl->D.file) {
		strcpy (head_d.title, "trend residuals");
		GMT_err_fail (GMT_write_grd (Ctrl->D.file, &head_d, resid, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), Ctrl->D.file);
	}
	if (Ctrl->W.file && Ctrl->N.robust) {
		strcpy (head_d.title, "trend weights");
		GMT_err_fail (GMT_write_grd (Ctrl->W.file, &head_d, weight, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), Ctrl->W.file);
	}

/* That's all, folks!  */

	if (weighted) GMT_free ((void *)weight);
	GMT_free ((void *)pstuff);
	GMT_free ((void *)gtd);
	GMT_free ((void *)gtg);
	GMT_free ((void *)old);
	GMT_free ((void *)yval);
	GMT_free ((void *)xval);
	GMT_free ((void *)resid);
	GMT_free ((void *)trend);
	GMT_free ((void *)data);

	Free_grdtrend_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void set_up_vals (double *val, GMT_LONG nval, double vmin, double vmax, double dv, GMT_LONG pixel_reg)
{
	GMT_LONG	i;
	double  v, middle, drange, true_min, true_max;

	true_min = (pixel_reg) ? vmin + 0.5 * dv : vmin;
	true_max = (pixel_reg) ? vmax - 0.5 * dv : vmax;

	middle = 0.5 * (true_min + true_max);
	drange = 2.0 / (true_max - true_min);
	for (i = 0; i < nval; i++) {
		v = true_min + i * dv;
		val[i] = (v - middle) * drange;
	}
	/* Just to be sure no rounding outside:  */
	val[0] = -1.0;
	val[nval - 1] = 1.0;
	return;
}

void load_pstuff (double *pstuff, GMT_LONG n_model, double x, double y, GMT_LONG newx, GMT_LONG newy)
{
	/* If either x or y has changed, compute new Legendre polynomials as needed  */

	if (newx) {
		if (n_model >= 2) pstuff[1] = x;
		if (n_model >= 5) pstuff[4] = 0.5*(3.0*pstuff[1]*pstuff[1] - 1.0);
		if (n_model >= 7) pstuff[6] = (5.0*pstuff[1]*pstuff[4] - 2.0*pstuff[1])/3.0;
	}
	if (newy) {
		if (n_model >= 3) pstuff[2] = y;
		if (n_model >= 6) pstuff[5] = 0.5*(3.0*pstuff[2]*pstuff[2] - 1.0);
		if (n_model >= 10) pstuff[9] = (5.0*pstuff[2]*pstuff[5] - 2.0*pstuff[2])/3.0;
	}
	/* In either case, refresh cross terms:  */

	if (n_model >= 4) pstuff[3] = pstuff[1]*pstuff[2];
	if (n_model >= 8) pstuff[7] = pstuff[4]*pstuff[2];
	if (n_model >= 9) pstuff[8] = pstuff[1]*pstuff[5];

	return;
}

void compute_trend (float *trend, GMT_LONG nx, GMT_LONG ny, double *xval, double *yval, double *gtd, GMT_LONG n_model, double *pstuff)
{
	GMT_LONG	i, j, k;
	GMT_LONG	ij;

	for (ij = 0, j = 0; j < ny; j++) {
		for (i = 0; i < nx; i++, ij++) {
			load_pstuff(pstuff, n_model, xval[i], yval[j], 1, (!(i)));
			trend[ij] = 0.0;
			for (k = 0; k < n_model; k++) trend[ij] += (float)(pstuff[k]*gtd[k]);
		}
	}
}

void compute_resid (float *data, float *trend, float *resid, GMT_LONG nxy)
{
	GMT_LONG	i;

	for (i = 0; i < nxy; i++) resid[i] = data[i] - trend[i];
	return;
}

void grd_trivial_model (float *data, GMT_LONG nx, GMT_LONG ny, double *xval, double *yval, double *gtd, GMT_LONG n_model)
{
	/* Routine to fit up elementary polynomial model of grd data, 
	model = gtd[0] + gtd[1]*x + gtd[2]*y + gtd[3] * x * y,
	where x,y are normalized to range [-1,1] and there are no
	NaNs in grid file, and problem is unweighted least squares.  */

	GMT_LONG	i, j;
	GMT_LONG	ij;
	double	x2, y2, sumx2 = 0.0, sumy2 = 0.0, sumx2y2 = 0.0;

	/* First zero the model parameters to use for sums:  */

	for (i = 0; i < n_model; i++) gtd[i] = 0.0;

	/* Now accumulate sums:  */

	for (ij = 0, j = 0; j < ny; j++) {
		y2 = yval[j] * yval[j];
		for (i = 0; i < nx; i++, ij++) {
			x2 = xval[i] * xval[i];
			sumx2 += x2;
			sumy2 += y2;
			sumx2y2 += (x2 * y2);
			gtd[0] += data[ij];
			if (n_model >= 2) gtd[1] += data[ij] * xval[i];
			if (n_model >= 3) gtd[2] += data[ij] * yval[j];
			if (n_model == 4) gtd[3] += data[ij] * xval[i] * yval[j];
		}
	}

	/* See how trivial it is?  */

	gtd[0] /= (nx * ny);
	if (n_model >= 2) gtd[1] /= sumx2;
	if (n_model >= 3) gtd[2] /= sumy2;
	if (n_model == 4) gtd[3] /= sumx2y2;

	return;
}

void compute_chisq (float *resid, float *weight, GMT_LONG nxy, double *chisq, double scale)
{
	GMT_LONG	i;
	double	tmp;

	*chisq = 0.0;
	for (i = 0; i < nxy; i++) {
		if (GMT_is_fnan (resid[i])) continue;
		tmp = resid[i];
		if (scale != 1.0) tmp /= scale;
		tmp *= tmp;
		if (weight[i] != 1.0) tmp *= weight[i];
		*chisq += tmp;
	}
	return;
}

void compute_robust_weight (float *resid, float *weight, GMT_LONG nxy, double *scale)
{
	GMT_LONG	i, j, j2;
	double	r, mad;

	for (i = j = 0; i < nxy; i++) {
		if (GMT_is_fnan (resid[i]))continue;
		weight[j] = (float)fabs((double)resid[i]);
		j++;
	}

	qsort ((void *)weight, (size_t)j, sizeof(float), GMT_comp_float_asc);

	j2 = j / 2;
	if (j%2)
		mad = weight[j2];
	else
		mad = 0.5 *(weight[j2] + weight[j2 - 1]);

	/* Adjust mad to equal Gaussian sigma:  */

	*scale = 1.4826 * mad;

	/* Use weight according to Huber (1981), but squared:  */

	for (i = 0; i < nxy; i++) {
		if (GMT_is_fnan (resid[i])) {
			weight[i] = resid[i];
			continue;
		}
		r = fabs(resid[i]) / (*scale);

		weight[i] = (float)((r <= 1.5) ? 1.0 : (3.0 - 2.25/r) / r);
	}
	return;
}

void write_model_parameters (double *gtd, GMT_LONG n_model)
{
	GMT_LONG	i;
	char	pbasis[10][16];
	char format[BUFSIZ];

	sprintf(pbasis[0], "Mean");
	sprintf(pbasis[1], "X");
	sprintf(pbasis[2], "Y");
	sprintf(pbasis[3], "X*Y");
	sprintf(pbasis[4], "P2(x)");
	sprintf(pbasis[5], "P2(y)");
	sprintf(pbasis[6], "P3(x)");
	sprintf(pbasis[7], "P2(x)*P1(y)");
	sprintf(pbasis[8], "P1(x)*P2(y)");
	sprintf(pbasis[9], "P3(y)");

	sprintf(format, "Coefficient fit to %%s:  %s\n", gmtdefs.d_format);
	for (i = 0; i < n_model; i++) fprintf (stderr, format, pbasis[i], gtd[i]);

	return;
}

void load_gtg_and_gtd (float *data, GMT_LONG nx, GMT_LONG ny, double *xval, double *yval, double *pstuff, double *gtg, double *gtd, GMT_LONG n_model, float *weight, GMT_LONG weighted)
{
	/* Routine to load the matrix G'G (gtg) and vector G'd (gtd)
	for the normal equations.  Routine uses indices i,j to refer
	to the grid file of data, and k,l to refer to the k_row, l_col
	of the normal equations matrix.  We need sums of [weighted]
	data and model functions in gtg and gtd.  We save time by
	loading only lower triangular part of gtg and then filling
	by symmetry after i,j loop.  */

	GMT_LONG	i, j, k, l, n_used;
	GMT_LONG	ij;

/*	First zero things out to start:  */

	n_used = 0;
	for (k = 0; k < n_model; k++) {
		gtd[k] = 0.0;
		for (l = 0; l < n_model; l++) gtg[k*n_model+l] = 0.0;
	}

/*  Now get going.  Have to load_pstuff separately in i and j,
	because it is possible that we skip data when i = 0.
	Loop over all data:  */

	for (ij = 0, j = 0; j < ny; j++ ) {
		load_pstuff(pstuff, n_model, xval[0], yval[j], 0, 1);
		for (i = 0; i < nx; i++, ij++) {

			if (GMT_is_fnan (data[ij]))continue;

			n_used++;
			load_pstuff(pstuff, n_model, xval[i], yval[j], 1, 0);

/* If weighted  */	if (weighted) {
				/* Loop over all gtg and gtd elements:  */
				gtd[0] += (data[ij] * weight[ij]);
				gtg[0] += (weight[ij]);
				for (k = 1; k < n_model; k++) {
					gtd[k] += (data[ij] * weight[ij] * pstuff[k]);
					gtg[k] += (weight[ij] * pstuff[k]);
					for (l = k; l < n_model; l++) gtg[k + l*n_model] += (pstuff[k]*pstuff[l]*weight[ij]);
				}
			}
/* If !weighted  */	else {
				/* Loop over all gtg and gtd elements:  */
				gtd[0] += data[ij];
				for (k = 1; k < n_model; k++) {
					gtd[k] += (data[ij] * pstuff[k]);
					gtg[k] += pstuff[k];
					for (l = k; l < n_model; l++) gtg[k + l*n_model] += (pstuff[k]*pstuff[l]);
				}
/* End if  */		}
		}
	}
/* End of loop over data i,j  */

/* Now if !weighted, use more accurate sum for gtg[0], and set symmetry:  */

	if (!weighted) gtg[0] = (double)n_used;

	for (k = 0; k < n_model; k++) {
		for (l = 0; l < k; l++) gtg[l + k*n_model] = gtg[k + l*n_model];
	}
/* That is all there is to it!  */

	return;
}

void *New_grdtrend_Ctrl () {	/* Allocate and initialize a new control structure */
	struct GRDTREND_CTRL *C;
	
	C = (struct GRDTREND_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GRDTREND_CTRL), "New_grdtrend_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
		
	return ((void *)C);
}

void Free_grdtrend_Ctrl (struct GRDTREND_CTRL *C) {	/* Deallocate control structure */
	if (C->D.file) free ((void *)C->D.file);	
	if (C->T.file) free ((void *)C->T.file);	
	if (C->W.file) free ((void *)C->W.file);	
	GMT_free ((void *)C);	
}
