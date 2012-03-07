/*--------------------------------------------------------------------
 *    $Id: fitcircle.c,v 1.53 2011/07/08 21:27:05 guru Exp $
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
/*
 * fitcircle <lonlatfile> [-L1] [-L2] [-S[<lat>]]
 *
 * Read lon,lat pairs from GMT_stdin[file].  Find mean position and pole
 * of best-fit circle through these points.  By default, fit great
 * circle.  If -S, fit small circle.  In this case, fit great circle
 * first, and then search for minimum small circle by bisection.
 *
 * Formally, we want to minimize some norm on the distance between
 * each point and the circle, measured perpendicular to the circle.
 * For both L1 and L2 norms this is a rather intractable problem.
 * (L2 is non-linear, and in L1 it is not clear how to proceed).
 * However, some approximations exist which work well and are simple
 * to compute.  We create a list of x,y,z vectors on the unit sphere,
 * representing the original data.  To find a great circle, do this:
 * For L1:
 * 	Find the Fisher mean of these data, call it mean position.
 *	Find the (Fisher) mean of all cross-products between data and
 *		the mean position; call this the pole to the great circle.
 *	Note that the cross-products are proportional to the distance
 *		between datum and mean; hence above average gives data far
 *		from mean larger weight in determining pole.  This is
 *		analogous to fitting line in plane, where data far from
 *		average abscissa have large leverage in determining slope.
 * For L2:
 *	Create 3 x 3 matrix of sums of products of data vector elements.
 *	Find eigenvectors and eigenvalues of this matrix.
 *	Find mean as eigenvector corresponding to max eigenvalue.
 *	Find pole as eigenvector corresponding to min eigenvalue.
 *	Eigenvalue-eigenvector decomposition performed by Jacobi's iterative
 *		method of successive Givens rotations.  Trials suggest that
 *		this converges extremely rapidly (3 sweeps, 9 rotations).
 *
 * To find a small circle, first find the great circle pole and the mean
 * position.  Suppose the small circle pole to lie in the plane containing
 * the mean and great circle pole, and narrow down its location by bisection.
 * Alternatively, specify the latitude of the small circle.
 *
 * Author:	W. H. F. Smith
 * Date:	16 September 1991.
 * Version:	4
 *
 */

#include "gmt.h"

struct FITCIRCLE_CTRL {	/* All control options for this program (except common args) */
	/* active is TRUE if the option has been activated */
	struct L {	/* -L[<n>] */
		GMT_LONG active;
		GMT_LONG norm;	/* 1, 2, or 3 (both) */
	} L;
	struct S {	/* -S */
		GMT_LONG active;
		GMT_LONG mode;	/* 0 = find latitude, 1 = use specified latitude */
		double lat;	/* 0 for great circle */
	} S;
};

struct	FITCIRCLE_DATA {
	double	x[3];
};

int main (int argc, char **argv)
{
	GMT_LONG	error = FALSE, greenwich = FALSE, allocate;

	char	format[BUFSIZ], *not_used = NULL;

	GMT_LONG	imin, imax, nrots, n_fields, n_expected_fields, n_read;
	GMT_LONG	i, j, k, n_alloc, n_data, n, np;

	double	lonsum, latsum, *in = NULL;
	double	meanv[3], cross[3], cross_sum[3], gcpole[3], scpole[3];		/* Extra vectors  */
	double	*a = NULL, *lambda = NULL, *v = NULL, *b = NULL, *z = NULL;	/* Matrix stuff */
	double	get_small_circle (struct FITCIRCLE_DATA *data, GMT_LONG ndata, double *center, double *gcpole, double *scpole, GMT_LONG norm, double *work, GMT_LONG mode, double slat);
	double	rad, *work = VNULL;

	FILE	*fp = NULL;

	struct	FITCIRCLE_DATA *data = NULL;
	struct FITCIRCLE_CTRL *Ctrl = NULL;

	void *New_fitcircle_Ctrl (), Free_fitcircle_Ctrl (struct FITCIRCLE_CTRL *C);

	argc = (int)GMT_begin (argc, argv);

 	Ctrl = (struct FITCIRCLE_CTRL *) New_fitcircle_Ctrl ();		/* Allocate and initialize defaults in a new control structure */

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

				case 'L':
					Ctrl->L.active = TRUE;
					Ctrl->L.norm = 3;
					if (argv[i][2]) Ctrl->L.norm = atoi(&argv[i][2]);
					break;
				case 'S':
					Ctrl->S.active = TRUE;
					if (argv[i][2]) {
						Ctrl->S.lat = atof (&argv[i][2]);
						Ctrl->S.mode = 1;
					}
                                        break;
				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
                        }
                }
		else {
			if ((fp = GMT_fopen(argv[i], GMT_io.r_mode)) == NULL) {
				fprintf (stderr, "%s: Could not open file %s\n", GMT_program, argv[i]);
				error = TRUE;
			}
		}
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "fitcircle %s - Find best-fitting great circle to points on sphere\n\n", GMT_VERSION);
		fprintf(stderr,"usage:  fitcircle [<input_file>] -L[<n>] [%s] [-S[<lat>]] [-V] [%s]\n", GMT_H_OPT, GMT_t_OPT);
		fprintf(stderr,"\t[%s] [%s]\n\n", GMT_bi_OPT, GMT_f_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf(stderr,"\tReads from input_file or standard input\n");
		fprintf(stderr,"\t-L specify norm as -L1 or -L2; or use -L or -L3 to give both.\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		GMT_explain_option ('H');
		fprintf(stderr,"\t-S will attempt to fit a small circle rather than a great circle.\n");
		fprintf(stderr,"\t   Optionally append the latitude of the small circle you want to fit\n");
		GMT_explain_option ('V');
		GMT_explain_option (':');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		GMT_explain_option ('f');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

 	if (Ctrl->L.norm < 1 || Ctrl->L.norm > 3) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -L option:  Choose between 1, 2, or 3\n", GMT_program);
		error++;
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
        if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = 2;
        if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] < 2) {
                fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data (-bi) must have at least 2 columns\n", GMT_program);
		error++;
	}
	if (error) exit (EXIT_FAILURE);

	if (GMT_io.binary[GMT_IN] && gmtdefs.verbose) {
		char *type[2] = {"double", "single"};
		fprintf (stderr, "%s: Expects %ld-column %s-precision binary data\n", GMT_program, GMT_io.ncol[GMT_IN], type[GMT_io.single_precision[GMT_IN]]);
	}


	if (fp == NULL) {
		fp = GMT_stdin;
#ifdef SET_IO_MODE
		GMT_setmode (GMT_IN);
#endif
	}
	n_alloc = GMT_CHUNK;
	n_data = n_read = 0;
	lonsum = latsum = 0.0;

	data = (struct FITCIRCLE_DATA *) GMT_memory (VNULL, (size_t)n_alloc, sizeof(struct FITCIRCLE_DATA), GMT_program);

	if (GMT_io.io_header[GMT_IN]) for (i = 0; i < GMT_io.n_header_recs; i++) not_used = GMT_fgets (format, BUFSIZ, fp);

	sprintf (format, "%s\t%s", gmtdefs.d_format, gmtdefs.d_format);
	n_expected_fields = (GMT_io.ncol[GMT_IN]) ? GMT_io.ncol[GMT_IN] : 2;

	while ((n_fields = GMT_input (fp, &n_expected_fields, &in)) >= 0 && !(GMT_io.status & GMT_IO_EOF)) {	/* Not yet EOF */

		n_read++;

		if (GMT_io.status & GMT_IO_MISMATCH) {
			fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld (skipped)\n", GMT_program, n_fields,  n_expected_fields, n_data);
			continue;
		}

		lonsum += in[GMT_X];
		latsum += in[GMT_Y];
		GMT_geo_to_cart (in[GMT_Y], in[GMT_X], data[n_data].x, TRUE);
		n_data++;

		if (n_data == n_alloc) {
			n_alloc <<= 1;
      			data = (struct FITCIRCLE_DATA *) GMT_memory ((void *)data, (size_t)n_alloc, sizeof(struct FITCIRCLE_DATA), GMT_program);
	        }
        }
        if (fp != GMT_stdin) GMT_fclose(fp);
        
      	data = (struct FITCIRCLE_DATA *) GMT_memory ((void *)data, (size_t)n_data, sizeof(struct FITCIRCLE_DATA), GMT_program);
	allocate = (Ctrl->S.active && (Ctrl->L.norm%2)) ;
	if (allocate) work = (double *) GMT_memory (VNULL, (size_t)n_data, sizeof(double), GMT_program);

	lonsum /= n_data;
	latsum /= n_data;
	if (gmtdefs.verbose) {
		fprintf (stderr, "%s: %ld points read, Average Position (Flat Earth): ", GMT_program, n_data);
		fprintf (stderr, format, lonsum, latsum);
		fprintf (stderr, "\n");
	}

	if (lonsum > 180.0) greenwich = TRUE;

	/* Get Fisher mean in any case, in order to set L2 mean correctly, if needed.  */


	meanv[0] = meanv[1] = meanv[2] = 0.0;

	for (i = 0; i < n_data; i++) for (j = 0; j < 3; j++) meanv[j] += data[i].x[j];

	GMT_normalize3v (meanv);

	if (Ctrl->L.norm%2) {
		GMT_cart_to_geo (&latsum, &lonsum, meanv, TRUE);
		if (greenwich && lonsum < 0.0) lonsum += 360.0;
		fprintf (stdout, format, lonsum, latsum);
		fprintf (stdout, "\tL1 Average Position (Fisher's Method)\n");

		cross_sum[0] = cross_sum[1] = cross_sum[2] = 0.0;
		for (i = 0; i < n_data; i++) {
			GMT_cross3v (&data[i].x[0], meanv, cross);
			if (cross[2] < 0.0) {
				cross_sum[0] -= cross[0];
				cross_sum[1] -= cross[1];
				cross_sum[2] -= cross[2];
			}
			else {
				cross_sum[0] += cross[0];
				cross_sum[1] += cross[1];
				cross_sum[2] += cross[2];
			}
		}
		GMT_normalize3v (cross_sum);
		if (Ctrl->S.active) for (i = 0; i < 3; i++) gcpole[i] = cross_sum[i];

		GMT_cart_to_geo (&latsum, &lonsum, cross_sum, TRUE);
		if (greenwich && lonsum < 0.0) lonsum += 360.0;
		fprintf (stdout, format, lonsum, latsum);
		fprintf (stdout, "\tL1 N Hemisphere Great Circle Pole (Cross-Averaged)\n");
		latsum = -latsum;
		lonsum = d_atan2d(-cross_sum[1], -cross_sum[0]);
		if (greenwich && lonsum < 0.0) lonsum += 360.0;
		fprintf (stdout, format, lonsum, latsum);
		fprintf (stdout, "\tL1 S Hemisphere Great Circle Pole (Cross-Averaged)\n");
		if (Ctrl->S.active) {
			if (gmtdefs.verbose) fprintf (stderr,"Fitting small circle using L1 norm.\n");
			rad = get_small_circle (data, n_data, meanv, gcpole, scpole, 1, work, Ctrl->S.mode, Ctrl->S.lat);
			if (rad >= 0.0) {
				GMT_cart_to_geo (&latsum, &lonsum, scpole, TRUE);
				if (greenwich && lonsum < 0.0) lonsum += 360.0;
				fprintf (stdout, format, lonsum, latsum);
				fprintf (stdout, "\tL1 Small Circle Pole.  ");
				sprintf(format, "Distance from Pole to L1 Small Circle (degrees):     %s\n", gmtdefs.d_format);
				fprintf (stdout, format, rad);
			}
		}
	}

	sprintf (format, "%s\t%s", gmtdefs.d_format, gmtdefs.d_format);
	if (Ctrl->L.norm/2) {

		n = 3;
		np = n;

		a = (double *) GMT_memory (VNULL, (size_t)np*np, sizeof(double), GMT_program);
		lambda = (double *) GMT_memory (VNULL, (size_t)np, sizeof(double), GMT_program);
		b = (double *) GMT_memory (VNULL, (size_t)np, sizeof(double), GMT_program);
		z = (double *) GMT_memory (VNULL, (size_t)np, sizeof(double), GMT_program);
		v = (double *) GMT_memory (VNULL, (size_t)np*np, sizeof(double), GMT_program);

		for (i = 0; i < n_data; i++) for (j = 0; j < n; j++) for (k = 0; k < n; k++)
			a[j + k*np] += (data[i].x[j]*data[i].x[k]);

		if (GMT_jacobi (a, &n, &np, lambda, v, b, z, &nrots)) {
			fprintf(stderr,"%s: Eigenvalue routine failed to converge in 50 sweeps.\n", GMT_program);
			fprintf(stderr,"%s: The reported L2 positions might be garbage.\n", GMT_program);
		}
		if (gmtdefs.verbose) fprintf(stderr,"%s: Eigenvalue routine converged in %ld rotations.\n", GMT_program, nrots);
		imax = 0;
		imin = 2;
		if (d_acos (GMT_dot3v (v, meanv)) > M_PI_2)
			for (i = 0; i < 3; i++) meanv[i] = -v[imax*np+i];
		else
			for (i = 0; i < 3; i++) meanv[i] = v[imax*np+i];
		GMT_cart_to_geo (&latsum, &lonsum, meanv, TRUE);
		if (greenwich && lonsum < 0.0) lonsum += 360.0;
		fprintf (stdout, format, lonsum, latsum);
		fprintf (stdout, "\tL2 Average Position (Eigenval Method)\n");

		if (v[imin*np+2] < 0.0)	/* Eigvec is in S Hemisphere  */
			for (i = 0; i < 3; i++) gcpole[i] = -v[imin*np+i];
		else
			for (i = 0; i < 3; i++) gcpole[i] = v[imin*np+i];

		GMT_cart_to_geo (&latsum, &lonsum, gcpole, TRUE);
		if (greenwich && lonsum < 0.0) lonsum += 360.0;
		fprintf (stdout, format, lonsum, latsum);
		fprintf (stdout, "\tL2 N Hemisphere Great Circle Pole (Eigenval Method)\n");
		latsum = -latsum;
		lonsum = d_atan2d(-gcpole[1], -gcpole[0]);
		if (greenwich && lonsum < 0.0) lonsum += 360.0;
		fprintf (stdout, format, lonsum, latsum);
		fprintf (stdout, "\tL2 S Hemisphere Great Circle Pole (Eigenval Method)\n");

		GMT_free ((void *)v);
		GMT_free ((void *)z);
		GMT_free ((void *)b);
		GMT_free ((void *)lambda);
		GMT_free ((void *)a);
		if (Ctrl->S.active) {
			if (gmtdefs.verbose) fprintf(stderr,"%s: Fitting small circle using L2 norm.\n", GMT_program);
			rad = get_small_circle(data, n_data, meanv, gcpole, scpole, 2, work, Ctrl->S.mode, Ctrl->S.lat);
			if (rad >= 0.0) {
				/* True when small circle fits better than great circle */
				GMT_cart_to_geo (&latsum, &lonsum, scpole, TRUE);
				if (greenwich && lonsum < 0.0) lonsum += 360.0;
				fprintf (stdout, format, lonsum, latsum);
				fprintf (stdout, "\tL2 Small Circle Pole.  ");
				sprintf(format, "Distance from Pole to L2 Small Circle (degrees):     %s\n", gmtdefs.d_format);
				fprintf (stdout, format, rad);
			}
		}
	}
	if (allocate) GMT_free ((void *)work);
	GMT_free ((void *)data);

	Free_fitcircle_Ctrl (Ctrl);	/* Deallocate control structure */
	
	GMT_end (argc, argv);
	
	exit (EXIT_SUCCESS);
}

double	get_small_circle (struct FITCIRCLE_DATA *data, GMT_LONG ndata, double *center, double *gcpole, double *scpole, GMT_LONG norm, double *work, GMT_LONG mode, double slat)
{

	/* Find scpole, the pole to the best-fit small circle, 
		by L(norm) iterative search along arc between
		center and +/- gcpole, the pole to the best fit
		great circle.  */

	GMT_LONG	i, j;
	double	temppole[3], a[3], b[3], oldpole[3];
	double	trypos, tryneg, afit, bfit, afactor, bfactor, fit, oldfit;
	double	length_ab, length_aold, length_bold, circle_misfit(struct FITCIRCLE_DATA *data, GMT_LONG ndata, double *pole, GMT_LONG norm, double *work, double *circle_distance), circle_distance;

	/* First find out if solution is between center and gcpole,
		or center and -gcpole:  */

	for (i = 0; i < 3; i++) temppole[i] = (center[i] + gcpole[i]);
	GMT_normalize3v (temppole);
	trypos = circle_misfit(data, ndata, temppole, norm, work, &circle_distance);

	for (i = 0; i < 3; i++) temppole[i] = (center[i] - gcpole[i]);
	GMT_normalize3v (temppole);
	tryneg = circle_misfit(data, ndata, temppole, norm, work, &circle_distance);

	if (tryneg < trypos) {
		for (i = 0; i < 3; i++) a[i] = center[i];
		for (i = 0; i < 3; i++) b[i] = -gcpole[i];
	}
	else {
		for (i = 0; i < 3; i++) a[i] = center[i];
		for (i = 0; i < 3; i++) b[i] = gcpole[i];
	}

	/* Now a is at center and b is at pole on correct side. */
	
	if (mode) {	/* Want a specified latitude */
		afactor = sind (slat);
		bfactor = cosd (slat);
		for (i = 0; i < 3; i++) scpole[i] = (afactor * a[i] + bfactor * b[i]);
		GMT_normalize3v (scpole);
		fit = circle_misfit (data, ndata, scpole, norm, work, &circle_distance);
		return (90.0-slat);
	}
	
	/* Try to bracket a minimum.  Move from b toward
	a in 1 degree steps:  */

	afit = circle_misfit(data, ndata, a, norm, work, &circle_distance);
	bfit = circle_misfit(data, ndata, b, norm, work, &circle_distance);
	j = 1;
	do {
		afactor = sind(j);
		bfactor = cosd(j);
		for (i = 0; i < 3; i++) temppole[i] = (afactor * a[i] + bfactor * b[i]);
		GMT_normalize3v (temppole);
		fit = circle_misfit(data, ndata, temppole, norm, work, &circle_distance);
		j++;
	} while (j < 90 && fit > bfit && fit > afit);

	if (j == 90) {
		/* Bad news.  There isn't a better fitting pole anywhere.  */
		fprintf(stderr,"%s: Sorry.  Cannot find small circle fitting better than great circle.\n", GMT_program);
		for (i = 0; i < 3; i++) scpole[i] = gcpole[i];
		return(-1.0);
	}
	/* Get here when temppole points to a minimum bracketed by a and b.  */

	for (i = 0; i < 3; i++) oldpole[i] = temppole[i];
	oldfit = fit;

	/* Now, while not converged, take golden section of wider interval.  */
	length_ab = d_acos (GMT_dot3v (a, b));
	length_aold = d_acos (GMT_dot3v (a, oldpole));
	length_bold = d_acos (GMT_dot3v (b, oldpole));
	do {
		if (length_aold > length_bold) {
			/* Section a_old  */
			for (i = 0; i < 3; i++) temppole[i] = (0.38197*a[i] + 0.61803*oldpole[i]);
			GMT_normalize3v (temppole);
			fit = circle_misfit(data, ndata, temppole, norm, work, &circle_distance);
			if (fit < oldfit) {
				/* Improvement.  b = oldpole, oldpole = temppole  */
				for (i = 0; i < 3; i++) {
					b[i] = oldpole[i];
					oldpole[i] = temppole[i];
				}
				oldfit = fit;
			}
			else {
				/* Not improved.  a = temppole  */
				for (i = 0; i < 3; i++) a[i] = temppole[i];
			}
		}
		else {
			/* Section b_old  */
			for (i = 0; i < 3; i++) temppole[i] = (0.38197*b[i] + 0.61803*oldpole[i]);
			GMT_normalize3v (temppole);
			fit = circle_misfit(data, ndata, temppole, norm, work, &circle_distance);
			if (fit < oldfit) {
				/* Improvement.  a = oldpole, oldpole = temppole  */
				for (i = 0; i < 3; i++) {
					a[i] = oldpole[i];
					oldpole[i] = temppole[i];
				}
				oldfit = fit;
			}
			else {
				/* Not improved.  b = temppole  */
				for (i = 0; i < 3; i++) b[i] = temppole[i];
			}
		}
		length_ab = d_acos (GMT_dot3v (a, b));
		length_aold = d_acos (GMT_dot3v (a, oldpole));
		length_bold = d_acos (GMT_dot3v (b, oldpole));
	} while (length_ab > 0.0001);	/* 1 milliradian = 0.05 degree  */

	for (i = 0; i < 3; i++) scpole[i] = oldpole[i];
	return (R2D * circle_distance);
}

double	circle_misfit(struct FITCIRCLE_DATA *data, GMT_LONG ndata, double *pole, GMT_LONG norm, double *work, double *circle_distance)
{
	/* Find the L(norm) misfit between a small circle through
		center with pole pole.  Return misfit in radians.  */

	double	distance, delta_distance, misfit = 0.0;
	GMT_LONG	i;

	/* At first, I thought we could use the center to define
		circle_dist = distance between pole and center.
		Then sum over data {dist[i] - circle_dist}.
		But it turns out that if the data are tightly
		curved, so that they are on a small circle 
		within a few degrees of the pole, then the
		center point is not on the small circle, and
		we cannot use it.  So, we first have to fit
		the circle_dist correctly:  */

	if (norm == 1) {
		for (i = 0; i < ndata; i++) work[i] = d_acos (GMT_dot3v (&data[i].x[0], pole));
		qsort((void *)work, (size_t)ndata, sizeof(double), GMT_comp_double_asc);
		if (ndata%2)
			*circle_distance = work[ndata/2];
		else
			*circle_distance = 0.5 * (work[(ndata/2)-1] + work[ndata/2]);
	}
	else {
		*circle_distance = 0.0;
		for (i = 0; i < ndata; i++) *circle_distance += d_acos (GMT_dot3v (&data[i].x[0], pole));
		*circle_distance /= ndata;
	}

	/* Now do each data point:  */

	for (i = 0; i < ndata; i++) {
		distance = d_acos (GMT_dot3v (&data[i].x[0], pole));
		delta_distance = fabs(*circle_distance - distance);
		misfit += ((norm == 1) ? delta_distance : delta_distance * delta_distance);
	}
	return (norm == 1) ? misfit : sqrt (misfit);
}

void *New_fitcircle_Ctrl () {	/* Allocate and initialize a new control structure */
	struct FITCIRCLE_CTRL *C;
	
	C = (struct FITCIRCLE_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct FITCIRCLE_CTRL), "New_fitcircle_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	
	C->L.norm = -1;
	
	return ((void *)C);
}

void Free_fitcircle_Ctrl (struct FITCIRCLE_CTRL *C) {	/* Deallocate control structure */
	GMT_free ((void *)C);	
}
