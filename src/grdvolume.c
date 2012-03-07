/*--------------------------------------------------------------------
 *	$Id: grdvolume.c,v 1.56 2011/07/08 21:27:06 guru Exp $
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
 * grdvolume reads a 2d binary gridded grid file, and calculates the volume
 * under the surface using exact integration of the bilinear interpolating
 * surface.  As an option, the user may supply a contour value; then the
 * volume is only integrated inside the chosen contour.
 *
 * Author:	Paul Wessel
 * Date:	23-SEP-1997
 * Revised:	02-JUN-1999
 * Version:	4
 */

#include "gmt.h"

struct GRDVOLUME_CTRL {
	struct C {	/* -C */
		GMT_LONG active;
		double low, high, inc;
	} C;
	struct L {	/* -L<base> */
		GMT_LONG active;
		double value;
	} L;
	struct S {	/* -S */
		GMT_LONG active;
		char unit;
	} S;
	struct T {	/* -T */
		GMT_LONG active;
	} T;
	struct Z {	/* Z<fact>[/<shift>] */
		GMT_LONG active;
		double scale, offset;
	} Z;
};

double vol_prism_frac_x (float *z, GMT_LONG ij, GMT_LONG nx, double x0, double x1, double a, double b, double c, double d);
double vol_prism_frac_y (float *z, GMT_LONG ij, GMT_LONG nx, double y0, double y1, double a, double b, double c, double d);

int main (int argc, char **argv)
{
	GMT_LONG error = FALSE, full = FALSE, bad, cut[4];

	char *grdfile = CNULL, format[BUFSIZ];

	GMT_LONG i, j, n = 0, c, k, pos, neg, nc, n_contours, nx, ny, mode = 0, nz = 0; 
	GMT_LONG ij, nm, ij_inc[4];
	
	float *f = NULL;

	double take_out, west, east, south, north, dv, da, cval = 0.0, cellsize, fact, dist_pr_deg;
	double *area = NULL, *vol = NULL, *height = NULL, this_base, small;

	struct GRDVOLUME_CTRL *Ctrl = NULL;

	struct GRD_HEADER grd;

	void SW_triangle (float f[], GMT_LONG ij, GMT_LONG nx, GMT_LONG triangle, double *dv, double *da);
	void NE_triangle (float f[], GMT_LONG ij, GMT_LONG nx, GMT_LONG triangle, double *dv, double *da);
	void SE_triangle (float f[], GMT_LONG ij, GMT_LONG nx, GMT_LONG triangle, double *dv, double *da);
	void NW_triangle (float f[], GMT_LONG ij, GMT_LONG nx, GMT_LONG triangle, double *dv, double *da);
	void NS_trapezoid (float f[], GMT_LONG ij, GMT_LONG nx, GMT_LONG right, double *dv, double *da);
	void EW_trapezoid (float f[], GMT_LONG ij, GMT_LONG nx, GMT_LONG top, double *dv, double *da);
	GMT_LONG ors_find_kink (double y[], GMT_LONG n, GMT_LONG mode);
	void *New_grdvolume_Ctrl (), Free_grdvolume_Ctrl (struct GRDVOLUME_CTRL *C);

	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct GRDVOLUME_CTRL *) New_grdvolume_Ctrl ();	/* Allocate and initialize a new control structure */
	
	west = east = south = north = 0.0;

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {

				/* Common parameters */

				case 'V':
					if (argv[i][2] == 'L' || argv[i][2] == 'l') full = TRUE;
				case 'R':
				case ':':
				case 'f':
				case '\0':
					error += GMT_parse_common_options (argv[i], &west, &east, &south, &north);
					break;


				/* Supplemental parameters */

				case 'C':
					Ctrl->C.active = TRUE;
					n = sscanf (&argv[i][2], "%lf/%lf/%lf", &Ctrl->C.low, &Ctrl->C.high, &Ctrl->C.inc);
					if (n == 3) {
						if (Ctrl->C.low >= Ctrl->C.high || Ctrl->C.inc <= 0.0) {
							fprintf (stderr, "%s: GMT SYNTAX ERROR -C:  high must exceed low and delta must be positive\n", GMT_program);
							error++;
						}
					}
					else
						Ctrl->C.high = Ctrl->C.low, Ctrl->C.inc = 1.0;	/* So calculation of ncontours will yield 1 */
					break;
				case 'L':
					Ctrl->L.active = TRUE;
					Ctrl->L.value = (argv[i][2]) ? atof (&argv[i][2]) : GMT_d_NaN;
					break;
				case 'S':
					Ctrl->S.active = TRUE;
                                        if (argv[i][2]) Ctrl->S.unit = argv[i][2];
					break;
				case 'T':
 					Ctrl->T.active = TRUE;
					break;
				case 'Z':
 					Ctrl->Z.active = TRUE;
					nz = (argv[i][2]) ? sscanf (&argv[i][2], "%lf/%lf", &Ctrl->Z.scale, &Ctrl->Z.offset) : -1;
					if (nz < 0 || nz > 2) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR option -Z: Must specify <fact> and optionally <shift>\n", GMT_program);
						error++;
					}
					break;
				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else
			grdfile = argv[i];
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr,"grdvolume %s - Calculating volume under a surface within a contour\n\n", GMT_VERSION);
		fprintf (stderr, "usage: grdvolume <grdfile> [-C<cval> or -C<low/high/delta>] [-L<base>] [-S[k]] [-T]\n\t[%s] [-V] [-Z<fact>[/<shift>]] [%s]\n", GMT_Rgeo_OPT, GMT_f_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\t<grdfile> is the name of the 2-D binary data set\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-C find area and volume inside the <cval> contour\n");
		fprintf (stderr, "\t   OR search using all contours from low to high\n");
		fprintf (stderr, "\t   [Default returns entire area and volume of grid]\n");
		fprintf (stderr, "\t-L Add volume from <base> up to contour [Default is from contour and up only]\n");
		fprintf (stderr, "\t-S Convert degrees to m, append k for km [Default is Cartesian]\n");
		fprintf (stderr, "\t-T Use curvature rather than maximum to find best contour value\n");
		GMT_explain_option ('R');
		GMT_explain_option ('V');
		fprintf (stderr, "\t   Append l for listing of all results (when contour search is selected)\n");
		fprintf (stderr, "\t-Z Subtract <shift> and then multiply data by <fact> before processing [1/0]\n");
		GMT_explain_option ('f');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	if (!grdfile) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify input grid file\n", GMT_program);
		error++;
	}
	if (Ctrl->C.active && !(n == 1 || n == 3)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR option -C: Must specify 1 or 3 arguments\n", GMT_program);
		error++;
	}
	if (Ctrl->S.active && !(Ctrl->S.unit == '\0' || Ctrl->S.unit == 'k')) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR option -S: May append k only\n", GMT_program);
		error++;
	}
	if (Ctrl->L.active && GMT_is_dnan (Ctrl->L.value)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR option -L: Must specify base\n", GMT_program);
		error++;
	}

	if (error) exit (EXIT_FAILURE);

	GMT_grd_init (&grd, argc, argv, FALSE);
	GMT_err_fail (GMT_read_grd_info (grdfile, &grd), grdfile);

	if (!project_info.region_supplied) {	/* No subset asked for */
		west = grd.x_min;
		east = grd.x_max;
		south = grd.y_min;
		north = grd.y_max;
	}
	else if (!project_info.region)	/* Got w/s/e/n, make into w/e/s/n */
		d_swap (south, east);

	nx = GMT_get_n (west, east, grd.x_inc, grd.node_offset);
	ny = GMT_get_n (south, north, grd.y_inc, grd.node_offset);
	nm = GMT_get_nm (nx, ny);

	f = (float *) GMT_memory (VNULL, nm, sizeof (float), GMT_program);

	GMT_err_fail (GMT_read_grd (grdfile, &grd, f, west, east, south, north, GMT_pad, FALSE), grdfile);

	ij_inc[0] = 0;	ij_inc[1] = 1;	ij_inc[2] = 1 - nx;	ij_inc[3] = -(long)nx;
	mode = (Ctrl->T.active) ? 1 : 0;
	cellsize = grd.x_inc * grd.y_inc;
	if (Ctrl->S.active) {
		dist_pr_deg = project_info.DIST_M_PR_DEG;
		if (Ctrl->S.unit == 'k') dist_pr_deg *= 0.001;	/* Use km instead */
		cellsize *= dist_pr_deg * dist_pr_deg;
	}

	n_contours = (Ctrl->C.active) ? irint ((Ctrl->C.high - Ctrl->C.low) / Ctrl->C.inc) + 1 : 1;

	height = (double *) GMT_memory (VNULL, n_contours, sizeof (double), GMT_program);
	vol = (double *) GMT_memory (VNULL, n_contours, sizeof (double), GMT_program);
	area = (double *) GMT_memory (VNULL, n_contours, sizeof (double), GMT_program);

	if (!(Ctrl->Z.scale == 1.0 && Ctrl->Z.offset == 0.0)) {
		if (gmtdefs.verbose) fprintf (stderr, "%s: Subtracting %g and multiplying by %g\n", GMT_program, Ctrl->Z.offset, Ctrl->Z.scale);
		for (ij = 0; ij < nm; ij++) f[ij] = (float)((f[ij] - Ctrl->Z.offset) * Ctrl->Z.scale);
		grd.z_min = (grd.z_min - Ctrl->Z.offset) * Ctrl->Z.scale;
		grd.z_max = (grd.z_max - Ctrl->Z.offset) * Ctrl->Z.scale;
		if (Ctrl->Z.scale < 0.0) d_swap (grd.z_min, grd.z_max);
	}

	this_base = (Ctrl->L.active) ? Ctrl->L.value : 0.0;
	small = Ctrl->C.inc * 1.0e-6;

	for (c = 0; Ctrl->C.active && c < n_contours; c++) {	/* Trace contour, only count volumes inside contours */

		cval = Ctrl->C.low + c * Ctrl->C.inc;
		take_out = (c == 0) ? cval : Ctrl->C.inc;	/* Take out start contour the first time and just the increment subsequent times */

		for (ij = 0; ij < nm; ij++) {
			f[ij] -= (float)take_out;		/* Take out the zero value */
			if (f[ij] == 0.0) f[ij] = (float)small;	/* But we dont want exactly zero, just + or - */
		}
		if (Ctrl->L.active) this_base -= take_out;

		if (Ctrl->L.active && this_base >= 0.0) {
			fprintf (stderr, "%s: Base is > than contour - exiting\n", GMT_program);
			exit (EXIT_FAILURE);
		}

		for (j = 1, ij = nx; j < (int)ny; j++) {

			dv = da = 0.0;	/* Reset these for each row */

			for (i = 0; i < (int)(nx-1); i++, ij++) {

				/* Find if a contour goes through this bin */

				for (k = neg = pos = 0, bad = FALSE; !bad && k < 4; k++) {
					(f[ij+ij_inc[k]] <= (float)small) ? neg++ : pos++;
					if (GMT_is_fnan (f[ij+ij_inc[k]])) bad = TRUE;
				}

                                if (bad || neg == 4) continue;	/* Contour not crossing, go to next bin */

				if (pos == 4) {	/* Need entire prism */
					dv += 0.25 * (f[ij] + f[ij+1] + f[ij-nx] + f[ij-nx+1]);
					da += 1.0;
				}
				else {	/* Need partial prisms */

					for (k = nc = 0; k < 4; k++) cut[k] = FALSE;
					if ((f[ij+1] * f[ij]) < 0.0)       nc++, cut[0] = TRUE;	/* Crossing the S border */
					if ((f[ij+1] * f[ij+1-nx]) < 0.0)  nc++, cut[1] = TRUE;	/* Crossing the E border */
					if ((f[ij-nx] * f[ij+1-nx]) < 0.0) nc++, cut[2] = TRUE;	/* Crossing the N border */
					if ((f[ij-nx] * f[ij]) < 0.0)      nc++, cut[3] = TRUE;	/* Crossing the W border */

					if (nc < 2) continue;	/* Can happen if some nodes were 0 and then reset to small, thus passing the test */

					if (nc == 4) {	/* Saddle scenario */
						if (f[ij] > 0) {	/* Need both SW and NE triangles */
							SW_triangle (f, ij, nx, TRUE, &dv, &da);
							NE_triangle (f, ij, nx, TRUE, &dv, &da);
						}
						else {			/* Need both SE and NW corners */
							SE_triangle (f, ij, nx, TRUE, &dv, &da);
							NW_triangle (f, ij, nx, TRUE, &dv, &da);
						}

					}
					else if (cut[0]) {	/* Contour enters at S border ... */
						if (cut[1])	/* and exits at E border */
							SE_triangle (f, ij, nx, (f[ij+1] > 0.0), &dv, &da);
						else if (cut[2])	/* or exits at N border */
							NS_trapezoid (f, ij, nx, f[ij] < 0.0, &dv, &da);
						else			/* or exits at W border */
							SW_triangle (f, ij, nx, (f[ij] > 0.0), &dv, &da);
					}
					else if (cut[1]) {	/* Contour enters at E border */
						if (cut[2])	/* exits at N border */
							NE_triangle (f, ij, nx, (f[ij+1-nx] > 0.0), &dv, &da);
						else			/* or exits at W border */
							EW_trapezoid (f, ij, nx, f[ij] < 0.0, &dv, &da);
					}
					else			/* Contours enters at N border and exits at W */
						NW_triangle (f, ij, nx, (f[ij-nx] > 0.0), &dv, &da);
				}
			}
			ij++;

			fact = cellsize;
			/* Allow for shrinking of longitudes with latitude */
			if (Ctrl->S.active) fact *= cosd (grd.y_max - (j+0.5) * grd.y_inc);

			vol[c]  += dv * fact;
			area[c] += da * fact;
		}

		/* Adjust for lower starting base */
		if (Ctrl->L.active) vol[c] -= area[c] * this_base;
	}
	if (!Ctrl->C.active) {	/* Columns with bilinear tops */
		if (Ctrl->L.active && Ctrl->L.value >= grd.z_min) {
			fprintf (stderr, "%s: Base is > than minimum z - exiting\n", GMT_program);
			exit (EXIT_FAILURE);
		}
		for (j = 0, ij = 0; j < grd.ny; j++) {
			dv = da = 0.0;
			for (i = 0; i < grd.nx; i++, ij++) {
				if (GMT_is_fnan (f[ij])) continue;

				/* Half the leftmost and rightmost cell */
				if (!grd.node_offset && (i == 0 || i == grd.nx-1)) {
					dv += 0.5 * f[ij];
					da += 0.5;
				}
				else {
					dv += f[ij];
					da += 1.0;
				}
			}

			fact = cellsize;
			/* Allow for shrinking of longitudes with latitude */
			if (Ctrl->S.active) fact *= cosd (grd.y_max - j * grd.y_inc);
			/* Half the top and bottom row */
			if (!grd.node_offset && (j == 0 || j == grd.ny-1)) fact *= 0.5;

			vol[0]  += dv * fact;
			area[0] += da * fact;
		}

		/* Adjust for lower starting base */
		if (Ctrl->L.active) vol[0] -= area[0] * this_base;
	}

	/* Compute average heights */

	for (c = 0; c < n_contours; c++) height[c] = (area[c] > 0.0) ? vol[c] / area[c] : GMT_d_NaN;

	/* Find the best contour that gives largest height */

	c = (Ctrl->C.active) ? ors_find_kink (height, n_contours, mode) : 0;

	/* Print out final estimates */

        sprintf (format, "%s\t%s\t%s\t%s\n", gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format);

	if (full) {
		sprintf (format, "%s\t%s\t%s\t%s\n", gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format);
		for (c = 0; c < n_contours; c++) fprintf (stdout, format, Ctrl->C.low + c * Ctrl->C.inc, area[c], vol[c], height[c]);
	}
	else if (Ctrl->C.active) {
        	sprintf (format, "%s\t%s\t%s\t%s\n", gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format);
		fprintf (stdout, format, Ctrl->C.low + c * Ctrl->C.inc, area[c], vol[c], height[c]);
	}
	else {
        	sprintf (format, "%s\t%s\t%s\n", gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format);
		fprintf (stdout, format, area[c], vol[c], height[c]);
	}

	GMT_free ((void *)f);
	GMT_free ((void *)area);
	GMT_free ((void *)vol);
	GMT_free ((void *)height);

	Free_grdvolume_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void SW_triangle (float f[], GMT_LONG ij, GMT_LONG nx, GMT_LONG triangle, double *dv, double *da)
{	/* Calculates area of a SW-corner triangle */
	/* triangle = TRUE gets triangle, FALSE gives the complementary area */
	double x1, y0, frac;

	x1 = f[ij] / (f[ij] - f[ij+1]);
	y0 = f[ij] / (f[ij] - f[ij-nx]);
	frac = (x1 == 0.0) ? 0.0 : vol_prism_frac_x (f, ij, nx, 0.0, x1, 0.0, 0.0, -y0 / x1, y0);
	if (triangle) {
		*dv += frac;
		*da += 0.5 * x1 * y0;
	}
	else {
		*dv += 0.25 * (f[ij] + f[ij+1] + f[ij-nx] + f[ij-nx+1]) - frac;
		*da += 1.0 - 0.5 * x1 * y0;
	}
}

void NE_triangle (float f[], GMT_LONG ij, GMT_LONG nx, GMT_LONG triangle, double *dv, double *da)
{	/* Calculates area of a NE-corner triangle */
	/* triangle = TRUE gets triangle, FALSE gives the complementary area */
	double x0, y1, a, x0_1, y1_1, frac = 0.0;

	x0 = f[ij-nx] / (f[ij-nx] - f[ij+1-nx]);
	y1 = f[ij+1] / (f[ij+1] - f[ij+1-nx]);
	x0_1 = 1.0 - x0;
	y1_1 = y1 - 1.0;
	if (x0_1 != 0.0) {
		a = y1_1 / x0_1;
		frac = vol_prism_frac_x (f, ij, nx, x0, 1.0, a, 1.0 - a * x0, 0.0, 0.0);
	}
	if (triangle) {
		*dv += frac;
		*da -= 0.5 * x0_1 * y1_1;	/* -ve because we need 1 - y1 */
	}
	else {
		*dv += 0.25 * (f[ij] + f[ij+1] + f[ij-nx] + f[ij-nx+1]) - frac;
		*da += 1.0 + 0.5 * x0_1 * y1_1;	/* +ve because we need 1 - y1 */
	}
}

void SE_triangle (float f[], GMT_LONG ij, GMT_LONG nx, GMT_LONG triangle, double *dv, double *da)
{	/* Calculates area of a SE-corner triangle */
	/* triangle = TRUE gets triangle, FALSE gives the complementary area */
	double x0, y1, c, x0_1, frac = 0.0;

	x0 = f[ij] / (f[ij] - f[ij+1]);
	y1 = f[ij+1] / (f[ij+1] - f[ij+1-nx]);
	x0_1 = 1.0 - x0;
	if (x0_1 != 0.0) {
		c = y1 / x0_1;
		frac = vol_prism_frac_x (f, ij, nx, x0, 1.0, 0.0, 0.0, c, -c * x0);
	}
	if (triangle) {
		*dv += frac;
		*da += 0.5 * x0_1 * y1;
	}
	else {
		*dv += 0.25 * (f[ij] + f[ij+1] + f[ij-nx] + f[ij-nx+1]) - frac;
		*da += 1.0 - 0.5 * x0_1 * y1;
	}
}

void NW_triangle (float f[], GMT_LONG ij, GMT_LONG nx, GMT_LONG triangle, double *dv, double *da)
{	/* Calculates area of a NW-corner triangle */
	/* triangle = TRUE gets triangle, FALSE gives the complementary area */
	double x1, y0, y0_1, frac;

	x1 = f[ij-nx] / (f[ij-nx] - f[ij+1-nx]);
	y0 = f[ij] / (f[ij] - f[ij-nx]);
	y0_1 = 1.0 - y0;
	frac = (x1 == 0.0) ? 0.0 : vol_prism_frac_x (f, ij, nx, 0.0, x1, y0_1 / x1, y0, 0.0, 1.0);
	if (triangle) {
		*dv += frac;
		*da += 0.5 * x1 * y0_1;
	}
	else {
		*dv += 0.25 * (f[ij] + f[ij+1] + f[ij-nx] + f[ij-nx+1]) - frac;
		*da += 1.0 - 0.5 * x1 * y0_1;
	}
}

void NS_trapezoid (float f[], GMT_LONG ij, GMT_LONG nx, GMT_LONG right, double *dv, double *da)
{	/* Calculates area of a NS trapezoid */
	/* right = TRUE gets the right trapezoid, FALSE gets the left */
	double x0, x1;

	x0 = f[ij] / (f[ij] - f[ij+1]);
	x1 = f[ij-nx] / (f[ij-nx] - f[ij+1-nx]);
	if (right) {	/* Need right piece */
		*dv += vol_prism_frac_y (f, ij, nx, 0.0, 1.0, x1 - x0, x0, 0.0, 1.0);
		*da += 0.5 * (2.0 - x0 - x1);
	}
	else {
		*dv += vol_prism_frac_y (f, ij, nx, 0.0, 1.0, 0.0, 0.0, x1 - x0, x0);
		*da += 0.5 * (x0 + x1);
	}
}

void EW_trapezoid (float f[], GMT_LONG ij, GMT_LONG nx, GMT_LONG top, double *dv, double *da)
{	/* Calculates area of a EW trapezoid */
	/* top = TRUE gets the top trapezoid, FALSE gets the bottom */
	double y0, y1;

	y0 = f[ij] / (f[ij] - f[ij-nx]);
	y1 = f[ij+1] / (f[ij+1] - f[ij+1-nx]);
	if (top) {	/* Need top piece */
		*dv += vol_prism_frac_x (f, ij, nx, 0.0, 1.0, y1 - y0, y0, 0.0, 1.0);
		*da += 0.5 * (2.0 - y0 - y1);
	}
	else {
		*dv += vol_prism_frac_x (f, ij, nx, 0.0, 1.0, 0.0, 0.0, y1 - y0, y0);
		*da += 0.5 * (y0 + y1);
	}
}

/* This function returns the volume bounded by a trapezoid based on two vertical
 * lines x0 and x1 and two horizontal lines y0 = ax +b and y1 = cx + d
 */

double vol_prism_frac_x (float *z, GMT_LONG ij, GMT_LONG nx, double x0, double x1, double a, double b, double c, double d)
{
	double dzdx, dzdy, dzdxy, ca, db, c2a2, d2b2, cdab, v, x02, x12, x03, x04, x13, x14;

	dzdx = (z[ij+1] - z[ij]);
	dzdy = (z[ij-nx] - z[ij]);
	dzdxy = (z[ij-nx+1] + z[ij] - z[ij+1] - z[ij-nx]);

	ca = c - a;
	db = d - b;
	c2a2 = c * c - a * a;
	d2b2 = d * d - b * b;
	cdab = c * d - a * b;
	x02 = x0 * x0;	x03 = x02 * x0;	x04 = x02 * x02;
	x12 = x1 * x1;	x13 = x12 * x1;	x14 = x12 * x12;

	v = (3.0 * dzdxy * c2a2 * (x14 - x04) +
	     4.0 * (2.0 * dzdx * ca + dzdy * c2a2 + 2.0 * dzdxy * cdab) * (x13 - x03) +
	     6.0 * (2.0 * z[ij] * ca + 2.0 * dzdx * db + 2.0 * dzdy * cdab + dzdxy * d2b2) * (x12 - x02) +
	     12.0 * (2.0 * z[ij] * db + dzdy * d2b2) * (x1 - x0)) / 24.0;

	return (v);
}

/* This function returns the volume bounded by a trapezoid based on two horizontal
 * lines y0 and y1 and two vertical lines x0 = ay +b and x1 = cy + d
 */

double vol_prism_frac_y (float *z, GMT_LONG ij, GMT_LONG nx, double y0, double y1, double a, double b, double c, double d)
{
	double dzdx, dzdy, dzdxy, ca, db, c2a2, d2b2, cdab, v, y02, y03, y04, y12, y13, y14;

	dzdx = (z[ij+1] - z[ij]);
	dzdy = (z[ij-nx] - z[ij]);
	dzdxy = (z[ij-nx+1] + z[ij] - z[ij+1] - z[ij-nx]);

	ca = c - a;
	db = d - b;
	c2a2 = c * c - a * a;
	d2b2 = d * d - b * b;
	cdab = c * d - a * b;
	y02 = y0 * y0;	y03 = y02 * y0;	y04 = y02 * y02;
	y12 = y1 * y1;	y13 = y12 * y1;	y14 = y12 * y12;

	v = (3.0 * dzdxy * c2a2 * (y14 - y04) +
	     4.0 * (2.0 * dzdy * ca + dzdx * c2a2 + 2.0 * dzdxy * cdab) * (y13 - y03) +
	     6.0 * (2.0 * z[ij] * ca + 2.0 * dzdy * db + 2.0 * dzdx * cdab + dzdxy * d2b2) * (y12 - y02) +
	     12.0 * (2.0 * z[ij] * db + dzdx * d2b2) * (y1 - y0)) / 24.0;

	return (v);
}

GMT_LONG ors_find_kink (double y[], GMT_LONG n, GMT_LONG mode)
{	/* mode: 0 = find maximum, 1 = find curvature kink */
	GMT_LONG i, ic, im;
	double *c, *f;
	double median3 (double x[]);

	if (mode == 0) {	/* Find maximum value */

		for (i = im = 0; i < n; i++) if (y[i] > y[im]) im = i;
		return (im);
	}

	/* Calculate curvatures */

	c = (double *) GMT_memory (VNULL, n, sizeof (double), GMT_program);

	for (i = 1; i < (n-1); i++) c[i] = y[i+1] - 2.0 * y[i] + y[i-1];
	c[0] = c[1];
	if (n > 1) c[n-1] = c[n-2];

	/* Apply 3-point median filter to curvatures  */

	f = (double *) GMT_memory (VNULL, n, sizeof (double), GMT_program);
	for (i = 1; i < (n-1); i++) f[i] = median3 (&c[i-1]);

	/* Find maximum negative filtered curvature */

	for (i = ic = 1; i < (n-1); i++) if (f[i] < f[ic]) ic = i;

	GMT_free ((void *)c);
	GMT_free ((void *)f);

	return (ic);
}

double median3 (double x[])
{

	if (x[0] < x[1]) {
		if (x[2] > x[1]) return (x[1]);
		if (x[2] > x[0]) return (x[2]);
		return (x[0]);
	}
	else {
		if (x[2] > x[0]) return (x[0]);
		if (x[2] < x[1]) return (x[1]);
		return (x[2]);
	}
}

void *New_grdvolume_Ctrl () {	/* Allocate and initialize a new control structure */
	struct GRDVOLUME_CTRL *C;
	
	C = (struct GRDVOLUME_CTRL *) GMT_memory (VNULL, 1, sizeof (struct GRDVOLUME_CTRL), "New_grdvolume_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	C->L.value = GMT_d_NaN;
	C->Z.scale = 1.0;
	return ((void *)C);
}

void  Free_grdvolume_Ctrl (struct GRDVOLUME_CTRL *C) {	/* Deallocate control structure */
	GMT_free ((void *)C);	
}
