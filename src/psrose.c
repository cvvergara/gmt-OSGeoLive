/*--------------------------------------------------------------------
 *	$Id: psrose.c 10173 2014-01-01 09:52:34Z pwessel $
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
 * psrose reads a file [or standard input] with azimuth and length information
 * and draws a sector or rose diagram.  Several options for plot layout are available.
 * 2 diagrams are possible: Full circle (360) or half circle (180).  In the
 * latter case azimuths > 180 are reversed (-= 180).  psrose are fully
 * compatible with other gmtsystem v2.0 plot files and overlays.
 *
 * To be compatible with GMT, I assume radial distance to be "x"
 * and azimuth to be "y".  Hence, west = 0.0 and east = max_radius
 * south/north is -90,90 for halfcircle and 0,360 for full circle
 *
 * Author:	Paul Wessel
 * Date:	20-JUN-2000
 * Version:	4
 */
 
#include "gmt.h"
#include "pslib.h"

struct PSROSE_CTRL {	/* All control options for this program (except common args) */
	/* active is TRUE if the option has been activated */
	struct A {	/* -A<sector_angle>[r] */
		GMT_LONG active;
		GMT_LONG rose;
		double inc;
	} A;
	struct C {	/* -C[<modefile>] */
		GMT_LONG active;
		char *file;
	} C;
	struct D {	/* -D */
		GMT_LONG active;
	} D;
	struct E {	/* -Eazimuth/elevation */
		GMT_LONG active;
		double azimuth, elevation;
	} E;
	struct F {	/* -F */
		GMT_LONG active;
	} F;
	struct G {	/* -G<fill> */
		GMT_LONG active;
		struct GMT_FILL fill;
	} G;
	struct I {	/* -I */
		GMT_LONG active;
	} I;
	struct L {	/* -L */
		GMT_LONG active;
		char w[GMT_LONG_TEXT];
		char e[GMT_LONG_TEXT];
		char s[GMT_LONG_TEXT];
		char n[GMT_LONG_TEXT];
	} L;
	struct M {	/* -M<params> */
		GMT_LONG active;
		double width, length, thickness;
		int rgb[3];
	} M;
	struct N {	/* -N */
		GMT_LONG active;
	} N;
	struct S {	/* -Sscale[n] */
		GMT_LONG active;
		GMT_LONG normalize;
		double scale;
	} S;
	struct T {	/* -T */
		GMT_LONG active;
	} T;
	struct W {	/* -W<pen> */
		GMT_LONG active;
		struct GMT_PEN pen;
	} W;
	struct Z {	/* -Zscale */
		GMT_LONG active;
		double scale;
	} Z;
};

int main (int argc, char **argv)
{
	GMT_LONG error = FALSE, find_mean = FALSE, half_only = FALSE, windrose = TRUE;
	GMT_LONG automatic = FALSE, sector_plot = FALSE;

	char text[BUFSIZ], *file = CNULL, format[BUFSIZ], *choice[2] = {"OFF", "ON"};
	char txt_a[GMT_LONG_TEXT], txt_b[GMT_LONG_TEXT], txt_c[GMT_LONG_TEXT], txt_d[GMT_LONG_TEXT];

	GMT_LONG n_bins, n_annot;
	GMT_LONG status, n_alpha, n_modes, n_read, n_fields, n_expected_fields;
	GMT_LONG i, bin, n = 0, n_alloc = GMT_CHUNK;

	double max = 0.0, radius, az, x_origin, y_origin, tmp, one_or_two = 1.0, s, c, half_bin_width;
	double angle1, angle2, angle, x, y, mean_theta, mean_radius, xr = 0.0, yr = 0.0;
	double x1, x2, y1, y2, total = 0.0, total_arc, off, max_radius, az_offset, *in = NULL;
	double west = 0.0, east = 0.0, south = 0.0, north = 0.0, asize, lsize, this_az, mean_vector, mean_resultant;
	double *xx = NULL, *yy = NULL, *sum = NULL, *azimuth = NULL, *length = NULL, *mode_direction = NULL, *mode_length = NULL, dim[3];

	FILE *fp = NULL, *fpm = NULL;

	struct PSROSE_CTRL *Ctrl = NULL;
	struct GMT_FILL f;

	void plot3d (double x, double y, GMT_LONG pen);
	void arc3d (double x, double y, double r, double a1, double a2, GMT_LONG status);
	void *New_psrose_Ctrl (), Free_psrose_Ctrl (struct PSROSE_CTRL *C);

	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct PSROSE_CTRL *) New_psrose_Ctrl ();	/* Allocate and initialize a new control structure */

	asize = gmtdefs.annot_font_size[0] * GMT_u2u[GMT_PT][GMT_INCH];
	lsize = gmtdefs.annot_font_size[0] * GMT_u2u[GMT_PT][GMT_INCH];
	mode_direction = (double *) GMT_memory (VNULL, (size_t)1, sizeof (double), "psrose");
	mode_length = (double *) GMT_memory (VNULL, (size_t)1, sizeof (double), "psrose");

	/* Check and interpret the command line arguments */

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch(argv[i][1]) {

				/* Common parameters */

				case 'B':
				case 'H':
				case 'K':
				case 'O':
				case 'P':
				case 'R':
				case 'U':
				case 'V':
				case 'X':
				case 'x':
				case 'Y':
				case 'y':
				case 'b':
				case 'c':
				case ':':
				case '\0':
					error += GMT_parse_common_options (argv[i], &west, &east, &south, &north);
					west = 0.0;
					break;

				/* Supplemental options */

				case 'A':		/* Get Sector angle in degrees */
					Ctrl->A.active = TRUE;
					Ctrl->A.inc = atof (&argv[i][2]);
					if (argv[i][strlen (argv[i])-1] == 'r') Ctrl->A.rose = TRUE;
					break;
				case 'C':		/* Read mode file and plot directions */
					Ctrl->C.active = TRUE;
					if (argv[i][2]) Ctrl->C.file = strdup (&argv[i][2]);
					break;
				case 'D':		/* Center the bins */
					Ctrl->D.active = TRUE;
					break;
				case 'E':
					Ctrl->E.active = TRUE;
					sscanf (&argv[i][2], "%lf/%lf", &Ctrl->E.azimuth, &Ctrl->E.elevation);
					break;
				case 'F':		/* Disable scalebar plotting */
					Ctrl->F.active = TRUE;
					break;
				case 'G':		/* Set Gray shade */
					Ctrl->G.active = TRUE;
					if (GMT_getfill (&argv[i][2], &Ctrl->G.fill)) {
						GMT_fill_syntax ('G', " ");
						error++;
					}
					break;
				case 'I':		/* Compute statistics only - no plot */
					Ctrl->I.active = TRUE;
					break;
				case 'L':		/* Overwride default labeling */
					Ctrl->L.active = TRUE;
					if (argv[i][2]) {
						n = sscanf (&argv[i][2], "%[^/]/%[^/]/%[^/]/%s", Ctrl->L.w, Ctrl->L.e, Ctrl->L.s, Ctrl->L.n);
						if (n != 4) {
							fprintf (stderr, "%s: GMT SYNTAX ERROR -L option.  Correct syntax:\n", GMT_program);
							fprintf (stderr, "\t-L<westlabel/eastlabel/southlabel/<northlabel>>\n");
							error = TRUE;
						}
					}
					else
						Ctrl->L.w[0] = Ctrl->L.e[0] = Ctrl->L.s[0] = Ctrl->L.n[0] = '\0';
					break;
				case 'M':		/* Get arrow parameters */
					Ctrl->M.active = TRUE;
					if (argv[i][2]) {
						n = sscanf (&argv[i][2], "%[^/]/%[^/]/%[^/]/%s", txt_a, txt_b, txt_c, txt_d);
						if (n != 4 || GMT_getrgb (txt_d, Ctrl->M.rgb)) {
							fprintf (stderr, "%s: GMT SYNTAX ERROR -M option.  Correct syntax:\n", GMT_program);
							fprintf (stderr, "\t-M<tailwidth/headlength/headwidth/<color>>\n");
							error = TRUE;
						}
						else {
							Ctrl->M.thickness = GMT_convert_units (txt_a, GMT_INCH);
							Ctrl->M.length = GMT_convert_units (txt_b, GMT_INCH);
							Ctrl->M.width = GMT_convert_units (txt_c, GMT_INCH);
						}
					}
					break;
				case 'N':		/* Make sectors area be proportional to frequency instead of radius */
					Ctrl->N.active = TRUE;
					break;
				case 'S':		/* Get radius of unit circle in inches */
					Ctrl->S.active = TRUE;
					n = strlen (argv[i]) - 1;
					if (argv[i][n] == 'n') {
						Ctrl->S.normalize = TRUE;
						argv[i][n] = 0;
					}
					Ctrl->S.scale = GMT_convert_units (&argv[i][2], GMT_INCH);
					break;
				case 'T':		/* Oriented instead of directed data */
					Ctrl->T.active = TRUE;
					break;
				case 'W':		/* Get pen width for outline */
					Ctrl->W.active = TRUE;
					if (GMT_getpen (&argv[i][2], &Ctrl->W.pen)) {
						GMT_pen_syntax ('W', " ");
						error++;
					}
					break;
				case 'Z':		/* Scale radii before using data */
					Ctrl->Z.active = TRUE;
					Ctrl->Z.scale = atof (&argv[i][2]);
					break;
				default:		/* Options not recognized */
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else
			file = argv[i];
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr,"psrose %s - Polar histogram (rose diagram) plotter\n\n", GMT_VERSION);
		fprintf (stderr, "usage: psrose <infile> [-A<sector_angle>[r]] [%s]\n", GMT_B_OPT);
		fprintf (stderr, "\t[-C[<modes>]] [-D] [-E<azim>/<elev>] [-G<fill>] [%s] [-I] [-K] [-L[<wlab/elab/slab/nlab>]]\n", GMT_H_OPT);
		fprintf (stderr, "\t[-M<parameters>] [-N] [-O] [-P] [-R<r0/r1/theta0/theta1>] [-S<scale>[n]] [-T] [%s] [-V] [%s]\n", GMT_U_OPT, GMT_t_OPT);
		fprintf (stderr, "\t[-W<pen>] [%s] [%s] [-Z<scale>] [%s] [%s]\n\n", GMT_X_OPT, GMT_Y_OPT, GMT_c_OPT, GMT_bi_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\t<infile> (in ASCII or binary) has (length,azimuth) pairs.  If not given read standard input.\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-A sector width in degrees for sector diagram [Default is windrose].\n");
		fprintf (stderr, "\t   append r to get rose diagram.\n");
		GMT_explain_option ('B');
		fprintf (stderr, "\t   (Remember: radial is x-direction, azimuthal is y-direction).\n");
		fprintf (stderr, "\t-C plot vectors listed in the <modes> file.  If no file, use mean direction.\n");
		fprintf (stderr, "\t-D will center the sectors.\n");
		fprintf (stderr, "\t-E set azimuth and elevation of viewpoint for 3-D perspective [180/90].\n");
		fprintf (stderr, "\t-F Do not draw the scale length bar [Default plots scale in lower right corner].\n");
		GMT_fill_syntax ('G', "Specifies color for diagram [Default is no fill].");
		GMT_explain_option ('H');
		fprintf (stderr, "\t-I for inquire.  Only compute statistics - no plot is created.\n");
		GMT_explain_option ('K');
		fprintf (stderr, "\t-L Override default labels [Default is WEST/EAST/SOUTH/NORTH for full circle and 90W/90E/-/0 for half-circle].\n");
		fprintf (stderr, "\t   If no argument is given then labels will be disabled.  Give - to disable an individual label.\n");
		fprintf (stderr, "\t-M Append tailwidth/headlength/headwidth/r/g/b to set arrow attributes [0.03i/0.12i/0.1i/0/0/0].\n");
		fprintf (stderr, "\t-N normalizes rose plots for area, i.e., takes sqrt(r) before plotting [FALSE].\n");
		fprintf (stderr, "\t   Only applicable if normalization has been specified with -S<radius>n.\n");
		GMT_explain_option ('O');
		GMT_explain_option ('P');
		fprintf (stderr, "\t-R specifies the region.  (r0 = 0, r1 = max_radius.  For azimuth:\n");
		fprintf (stderr, "\t   Specify theta0/theta1 = -90/90 (half-circle) or 0/360 only).\n");
		fprintf (stderr, "\t   If r0 = r1 = 0, psrose will compute a reasonable r1 value.\n");
		fprintf (stderr, "\t-S specifies the radius of the unit circle in %s [%g]. Normalize r if n is appended.\n", GMT_unit_names[gmtdefs.measure_unit], Ctrl->S.scale);
		fprintf (stderr, "\t-T indicates that the vectors are oriented (two-headed), not directed [Default].\n");
		GMT_explain_option ('U');
		GMT_explain_option ('V');
		GMT_pen_syntax ('W', "sets pen attributes for outline of rose. [Default is no outline].");
		GMT_explain_option ('X');
		fprintf (stderr, "\t-Z multiply the radii by scale before plotting.\n");
		GMT_explain_option ('c');
		fprintf (stderr, "\t-: Expect (azimuth,radius) input rather than (radius,azimuth) [%s].\n", choice[gmtdefs.xy_toggle[GMT_IN]]);
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf (stderr, "\t   Default is 2 input columns.\n");
		fprintf (stderr, "\t(See man page for gmtdefaults to set these and other defaults to ON or OFF).\n");
		exit (EXIT_FAILURE);
	}

	/* Check that the options selected are mutually consistent */

	if (Ctrl->E.elevation <= 0.0 || Ctrl->E.elevation > 90.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -E option:  Elevation must be in 0-90 range\n", GMT_program);
		error++;
	}
	if (Ctrl->S.scale <= 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -S option:  radius must be nonzero\n", GMT_program);
		error++;
	}
	if (GMT_IS_ZERO (Ctrl->Z.scale)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -Z option:  factor must be nonzero\n", GMT_program);
		error++;
	}
	if (Ctrl->A.inc < 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -A option:  sector width must be positive\n", GMT_program);
		error++;
	}
	if (!((south == -90.0 && north == 90.0) || (south == 0.0 && north == 360.0))) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -R option:  theta0/theta1 must be either -90/90 or 0/360\n", GMT_program);
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

	z_project.view_azimuth = Ctrl->E.azimuth;
	z_project.view_elevation = Ctrl->E.elevation;
	if (Ctrl->C.active) {	/* Either calculate or read mean vector to plot */
		if (!Ctrl->C.file)	/* Not given, calculate mean direction */
			find_mean = TRUE;
		else if ((fpm = GMT_fopen (Ctrl->C.file, "r")) == NULL) {	/* oops */
			fprintf (stderr, "%s: GMT SYNTAX ERROR -C.  Cannot open file %s\n", GMT_program, Ctrl->C.file);
			exit (EXIT_FAILURE);
		}
	}
		
	if (file && (fp = GMT_fopen (file, GMT_io.r_mode)) == NULL) {
		fprintf (stderr, "%s:  Cannot open file %s\n", GMT_program, file);
		exit (EXIT_FAILURE);
	}

	if (fp == NULL) {
		fp = GMT_stdin;
#ifdef SET_IO_MODE
		GMT_setmode (GMT_IN);
#endif
	}

	max_radius = east;
	half_only = GMT_IS_ZERO (south + 90.0);
	if (Ctrl->A.rose) windrose = FALSE;
	sector_plot = (Ctrl->A.inc > 0.0);
	if (sector_plot) windrose = FALSE;	/* Draw rose diagram instead of sector diagram */
	if (!Ctrl->S.normalize) Ctrl->N.active = FALSE;	/* Only do this is data is normalized for length also */
	if (!Ctrl->I.active && !project_info.region_supplied) automatic = TRUE;
	if (Ctrl->T.active) one_or_two = 2.0;
	half_bin_width = Ctrl->D.active * Ctrl->A.inc * 0.5;
	if (half_only) {
		total_arc = 180.0;
		az_offset = 90.0;
	}
	else {
		total_arc = 360.0;
		az_offset = 0.0;
	}
	n_bins = (Ctrl->A.inc <= 0.0) ? 1 : irint (total_arc / Ctrl->A.inc);

	sum = (double *) GMT_memory (VNULL, (size_t)n_bins, sizeof (double), GMT_program);
	xx = (double *) GMT_memory (VNULL, (size_t)(n_bins+2), sizeof (double), GMT_program);
	yy = (double *) GMT_memory (VNULL, (size_t)(n_bins+2), sizeof (double), GMT_program);
	azimuth = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (double), GMT_program);
	length = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (double), GMT_program);

	if (GMT_io.io_header[GMT_IN]) for (i = 0; i < GMT_io.n_header_recs; i++) GMT_fgets (text, BUFSIZ, fp);

	/* Read data and do some stats */

	n = n_read = 0;
	n_expected_fields = (GMT_io.ncol[GMT_IN]) ? GMT_io.ncol[GMT_IN] : 2;

	n_fields = GMT_input (fp,  &n_expected_fields, &in);

	while (! (GMT_io.status & GMT_IO_EOF)) {	/* Not yet EOF */

		n_read++;

		if (GMT_io.status & GMT_IO_MISMATCH) {
			fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld\n", GMT_program, n_fields,  n_expected_fields, n_read);
			continue;
		}

		length[n] = in[GMT_X];
		azimuth[n] = in[GMT_Y];

		if (Ctrl->Z.scale != 1.0) length[n] *= Ctrl->Z.scale;

		/* Make sure azimuth is in 0 <= az < 360 range */

		while (azimuth[n] < 0.0) azimuth[n] += 360.0;
		while (azimuth[n] >= 360.0) azimuth[n] -= 360.0;

		if (half_only) {	/* Flip azimuths about E-W line i.e. -90 < az <= 90 */
			if (azimuth[n] > 90.0 && azimuth[n] <= 270.0) azimuth[n] -= 180.0;
			if (azimuth[n] > 270.0) azimuth[n] -= 360.0;
		}
		else if (Ctrl->T.active) {
			azimuth[n] = 0.5 * fmod (2.0 * azimuth[n], 360.0);
		}

		/* Double angle to find mean azimuth */

		sincosd (one_or_two * azimuth[n], &s, &c);
		xr += length[n] * c;
		yr += length[n] * s;

		total += length[n];

		n++;
		if (n == n_alloc) {	/* Get more memory */
			n_alloc <<= 1;
			azimuth = (double *) GMT_memory ((void *)azimuth, (size_t)n_alloc, sizeof (double), GMT_program);
			length = (double *) GMT_memory ((void *)length, (size_t)n_alloc, sizeof (double), GMT_program);
		}
		n_fields = GMT_input (fp, &n_expected_fields, &in);
	}
	if (fp != GMT_stdin) GMT_fclose (fp);

	if (Ctrl->A.inc > 0.0) {	/* Sum up sector diagram info */
		for (i = 0; i < n; i++) {
			if (Ctrl->D.active) {	/* Center bin by removing half bin width here */
				this_az = azimuth[i] - half_bin_width;
				if (!half_only && this_az < 0.0)   this_az += 360.0;
				if (half_only  && this_az < -90.0) this_az += 180.0;
			}
			else
				this_az = azimuth[i];
			bin = (int) ((this_az + az_offset) / Ctrl->A.inc);
			sum[bin] += length[i];
			if (Ctrl->T.active) {	/* Also count its other end */
				this_az += 180.0;	if (this_az >= 360.0) this_az -= 360.0;
				bin = (int) ((this_az + az_offset) / Ctrl->A.inc);
				sum[bin] += length[i];
			}
		}
	}

	mean_theta = d_atan2d (yr, xr) / one_or_two;
	if (mean_theta < 0.0) mean_theta += 360.0;
	mean_vector = hypot (xr, yr) / n;
	mean_resultant = mean_radius = hypot (xr, yr) / total;
	if (!Ctrl->S.normalize) mean_radius *= max_radius;

	if (Ctrl->A.inc > 0.0) {	/* Find max of the bins */
		for (bin = 0; bin < n_bins; bin++) if (sum[bin] > max) max = sum[bin];
		if (Ctrl->S.normalize) for (bin = 0; bin < n_bins; bin++) sum[bin] /= max;
		if (Ctrl->N.active) for (bin = 0; bin < n_bins; bin++) sum[bin] = sqrt (sum[bin]);
	}
	else {
		for (i = 0; i < n; i++) if (length[i] > max) max = length[i];
		if (Ctrl->S.normalize) {
			max = 1.0 / max;
			for (i = 0; i < n; i++) length[i] *= max;
			max = 1.0 / max;
		}
	}

	if (Ctrl->I.active || gmtdefs.verbose) {
		char *kind[2] = {"r", "bin sum"};
		if (file) strcpy (text, file); else strcpy (text, "<stdin>");
		sprintf (format, "%%s: Info for %%s: n = %%ld mean az = %s mean r = %s mean resultant length = %s max %s = %s scaled mean r = %s linear length sum = %s\n",
			gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format, kind[Ctrl->A.active], gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format);
		fprintf (stderr, format, GMT_program, text, n, mean_theta, mean_vector, mean_resultant, max, mean_radius, total);
		if (Ctrl->I.active) {
			GMT_free ((void *)sum);
			GMT_free ((void *)xx);
			GMT_free ((void *)yy);
			GMT_free ((void *)azimuth);
			GMT_free ((void *)length);
			exit (EXIT_FAILURE);
		}
	}

	if (automatic) {
		if (GMT_IS_ZERO (frame_info.axis[0].item[GMT_ANNOT_UPPER].interval)) {
			tmp = pow (10.0, floor (d_log10 (max)));
			if ((max / tmp) < 3.0) tmp *= 0.5;
		}
		else
			tmp = frame_info.axis[0].item[GMT_ANNOT_UPPER].interval;
		max_radius = ceil (max / tmp) * tmp;
		if (GMT_IS_ZERO (frame_info.axis[0].item[GMT_ANNOT_UPPER].interval) || GMT_IS_ZERO (frame_info.axis[0].item[GMT_GRID_UPPER].interval)) {	/* Tickmarks not set */
			frame_info.axis[0].item[GMT_ANNOT_UPPER].interval = frame_info.axis[0].item[GMT_GRID_UPPER].interval = tmp;
			frame_info.plot = TRUE;
		}
	}

	if (frame_info.plot && GMT_IS_ZERO (frame_info.axis[1].item[GMT_ANNOT_UPPER].interval)) frame_info.axis[1].item[GMT_ANNOT_UPPER].interval = frame_info.axis[1].item[GMT_GRID_UPPER].interval = 30.0;

	/* Ready to plot.  So set up GMT projections (not used by psrose), we set region to actual plot width and scale to 1 */

	project_info.projection = GMT_LINEAR;
	project_info.region = TRUE;
	project_info.pars[0] = project_info.pars[1] = 1.0;
	GMT_err_fail (GMT_map_setup (-Ctrl->S.scale, Ctrl->S.scale, -Ctrl->S.scale, Ctrl->S.scale), "");

	GMT_plotinit (argc, argv);

	x_origin = Ctrl->S.scale;	y_origin = ((half_only) ? 0.0 : Ctrl->S.scale);
	ps_transrotate (x_origin, y_origin, 0.0);
	if (project_info.three_D) ps_transrotate (-z_project.xmin, -z_project.ymin - y_origin, 0.0);
	if (!Ctrl->S.normalize) Ctrl->S.scale /= max_radius;

	GMT_setpen (&Ctrl->W.pen);
	if (windrose) {
		for (i = 0; i < n; i++) {
			sincosd (90.0 - azimuth[i], &s, &c);
			radius = length[i] * Ctrl->S.scale;
			x = radius * c;
			y = radius * s;
			if (Ctrl->T.active)
				plot3d (-x, -y, 3);
			else
				plot3d (0.0, 0.0, 3);
			plot3d (x, y, -2);
		}
	}

	if (sector_plot && !Ctrl->A.rose && Ctrl->G.fill.rgb[0] >= 0) {	/* Draw pie slices for sector plot if fill is requested */

		for (bin = 0; bin < n_bins; bin++) {
			az = bin * Ctrl->A.inc - az_offset + half_bin_width;
			dim[1] = (90.0 - az - Ctrl->A.inc);
			dim[2] = dim[1] + Ctrl->A.inc;
			dim[0] = sum[bin] * Ctrl->S.scale;
			GMT_pie (0.0, 0.0, project_info.z_level, dim, &(Ctrl->G.fill), FALSE);
		}
	}
	else if (Ctrl->A.rose) {	/* Draw rose diagram */

		for (bin = i = 0; bin < n_bins; bin++, i++) {
			az = (bin + 0.5) * Ctrl->A.inc - az_offset - half_bin_width;
			sincosd (90.0 - az, &s, &c);
			xx[i] = Ctrl->S.scale * sum[bin] * c;
			yy[i] = Ctrl->S.scale * sum[bin] * s;
		}
		if (half_only) {
			xx[i] = Ctrl->S.scale * 0.5 * (sum[0] + sum[n_bins-1]);
			yy[i] = 0.0;
			i++;
			xx[i] = -xx[i-1];
			yy[i] = 0.0;
			i++;
		}
		if (project_info.three_D) GMT_2D_to_3D (xx, yy, project_info.z_level, i);
		ps_polygon (xx, yy, i, Ctrl->G.fill.rgb, Ctrl->W.active);
	}

	if (sector_plot && Ctrl->W.active && !Ctrl->A.rose) {	/* Draw a line outlining the pie slices */
		angle1 = ((half_only) ? 180.0 : 90.0) - half_bin_width;
		angle2 = ((half_only) ?   0.0 : 90.0) - half_bin_width;
		sincosd (angle1, &s, &c);
		x1 = (sum[0] * Ctrl->S.scale) * c;
		y1 = (sum[0] * Ctrl->S.scale) * s;
		sincosd (angle2, &s, &c);
		x2 = (sum[n_bins-1] * Ctrl->S.scale) * c;
		y2 = (sum[n_bins-1] * Ctrl->S.scale) * s;
		plot3d (x1, y1, 3);
		plot3d (x2, y2, 2);
		for (bin = n_bins-1; bin >= 0; bin--) {
			status = (bin == 0) ? 2 : 0;
			az = bin * Ctrl->A.inc - az_offset + half_bin_width;
			angle1 = 90.0 - az - Ctrl->A.inc;
			angle2 = angle1 + Ctrl->A.inc;
			arc3d (0.0, 0.0, sum[bin] * Ctrl->S.scale, angle1, angle2, status);
		}
	}

	if (Ctrl->C.active) {

		if (find_mean) {	/* Use the mean direction only */
			n_modes = 1;
			mode_direction[0] = mean_theta;
			mode_length[0] = mean_radius;
		}
		else {	/* Get mode parameters from separate file */
			n_modes = 0;
			n_alloc = 1;
			while (GMT_fgets (text, BUFSIZ, fpm)) {
				if (GMT_is_a_blank_line (text)) continue;	/* Skip blank lines or # comments */
				sscanf (text, "%lf %lf", &mode_direction[n_modes], &mode_length[n_modes]);
				n_modes++;
				if (n_modes == n_alloc) {
					n_alloc <<= 1;
					mode_direction = (double *) GMT_memory ((void *)mode_direction, (size_t)n_alloc, sizeof (double), "psrose");
					mode_length = (double *) GMT_memory ((void *)mode_length, (size_t)n_alloc, sizeof (double), "psrose");
				}
			}
			GMT_fclose (fpm);
			mode_direction = (double *) GMT_memory ((void *)mode_direction, (size_t)n_modes, sizeof (double), "psrose");
			mode_length = (double *) GMT_memory ((void *)mode_length, (size_t)n_modes, sizeof (double), "psrose");
		}
		for (i = 0; i < n_modes; i++) {
			if (Ctrl->N.active) mode_length[i] = sqrt (mode_length[i]);
			if (half_only && mode_direction[i] > 90.0 && mode_direction[i] <= 270.0) mode_direction[i] -= 180.0;
			angle = 90.0 - mode_direction[i];
			sincosd (angle, &s, &c);
			xr = Ctrl->S.scale * mode_length[i] * c;
			yr = Ctrl->S.scale * mode_length[i] * s;
			GMT_init_fill (&f, Ctrl->M.rgb[0], Ctrl->M.rgb[1], Ctrl->M.rgb[2]);       /* Initialize fill structure */
			GMT_vector (0.0, 0.0, xr, yr, project_info.z_level, Ctrl->M.thickness, Ctrl->M.length, Ctrl->M.width, gmtdefs.vector_shape, &f, TRUE);
		}

	}

	if (Ctrl->L.active) {	/* Deactivate those with - */
		if (Ctrl->L.w[0] == '-' && Ctrl->L.w[1] == '\0') Ctrl->L.w[0] = '\0';
		if (Ctrl->L.e[0] == '-' && Ctrl->L.e[1] == '\0') Ctrl->L.e[0] = '\0';
		if (Ctrl->L.s[0] == '-' && Ctrl->L.s[1] == '\0') Ctrl->L.s[0] = '\0';
		if (Ctrl->L.n[0] == '-' && Ctrl->L.n[1] == '\0') Ctrl->L.n[0] = '\0';
	}
	else {	/* Use default labels */
		strcpy (Ctrl->L.w, "WEST");	strcpy (Ctrl->L.e, "EAST");	strcpy (Ctrl->L.s, "SOUTH");	strcpy (Ctrl->L.n, "NORTH");
	}
	if (frame_info.plot) {	/* Draw grid lines etc */
		GMT_setpen (&gmtdefs.grid_pen[0]);
		off = max_radius * Ctrl->S.scale;
		n_alpha = (frame_info.axis[1].item[GMT_GRID_UPPER].interval > 0.0) ? irint (total_arc / frame_info.axis[1].item[GMT_GRID_UPPER].interval) : -1;
		for (i = 0; i <= n_alpha; i++) {
			angle = i * frame_info.axis[1].item[GMT_GRID_UPPER].interval;
			sincosd (angle, &s, &c);
			x = max_radius * Ctrl->S.scale * c;
			y = max_radius * Ctrl->S.scale * s;
			plot3d (0.0, 0.0, 3);
			plot3d (x, y, -2);
		}

		n_bins = (frame_info.axis[0].item[GMT_GRID_UPPER].interval > 0.0) ? irint (max_radius / frame_info.axis[0].item[GMT_GRID_UPPER].interval) : -1;
		for (i = 1; i <= n_bins; i++)
			arc3d (0.0, 0.0, i * frame_info.axis[0].item[GMT_GRID_UPPER].interval * Ctrl->S.scale, 0.0, total_arc, 3);
		ps_setpaint (gmtdefs.basemap_frame_rgb);
		y = lsize + 6.0 * gmtdefs.annot_offset[0];
		GMT_text3D (0.0, off + y, project_info.z_level, gmtdefs.header_font_size, gmtdefs.header_font, frame_info.header, 0.0, 2, 0);

		GMT_get_format (frame_info.axis[0].item[GMT_ANNOT_UPPER].interval, frame_info.axis[0].unit, frame_info.axis[0].prefix, format);

		if (half_only) {
			GMT_LONG which, no_degree;
			which = (gmtdefs.degree_format >= 100) ? gmt_ring : gmt_degree;
			no_degree = (gmtdefs.degree_format >= 1000);	/* No, we don't want the degree symbol at all */
			if (!Ctrl->L.active) {	/* Use default labels */
				if (no_degree) {
					sprintf (Ctrl->L.w, "90W");
					sprintf (Ctrl->L.e, "90E");
					sprintf (Ctrl->L.n, "0");
				}
				else {
					sprintf (Ctrl->L.w, "90%cW", (int)gmtdefs.encoding.code[which]);
					sprintf (Ctrl->L.e, "90%cE", (int)gmtdefs.encoding.code[which]);
					sprintf (Ctrl->L.n, "0%c",   (int)gmtdefs.encoding.code[which]);
				}
			}
			y = -(3.0 * gmtdefs.annot_offset[0] + GMT_font[gmtdefs.annot_font[0]].height * asize);
			if (frame_info.axis[0].label[0]) GMT_text3D (0.0, y, project_info.z_level, gmtdefs.label_font_size, gmtdefs.label_font, frame_info.axis[0].label, 0.0, 10, 0);
			y = -(5.0 * gmtdefs.annot_offset[0] + GMT_font[gmtdefs.annot_font[0]].height * lsize + GMT_font[gmtdefs.label_font].height * lsize);
			if (frame_info.axis[1].label[0]) GMT_text3D (0.0, y, project_info.z_level, gmtdefs.label_font_size, gmtdefs.label_font, frame_info.axis[1].label, 0.0, 10, 0);
			GMT_text3D (0.0, -gmtdefs.annot_offset[0], project_info.z_level, gmtdefs.annot_font_size[0], gmtdefs.annot_font[0], "0", 0.0, 10, 0);
			n_annot = (frame_info.axis[0].item[GMT_ANNOT_UPPER].interval > 0.0) ? irint (max_radius / frame_info.axis[0].item[GMT_ANNOT_UPPER].interval) : -1;
			for (i = 1; i <= n_annot; i++) {
				x = i * frame_info.axis[0].item[GMT_ANNOT_UPPER].interval;
				sprintf (text, format, x);
				x *= Ctrl->S.scale;
				GMT_text3D (x, -gmtdefs.annot_offset[0], project_info.z_level, gmtdefs.annot_font_size[0], gmtdefs.annot_font[0], text, 0.0, 10, 0);
				GMT_text3D (-x, -gmtdefs.annot_offset[0], project_info.z_level, gmtdefs.annot_font_size[0], gmtdefs.annot_font[0], text, 0.0, 10, 0);
			}
		}
		else {
			GMT_text3D (0.0, -off - 2.0 * gmtdefs.annot_offset[0], project_info.z_level, gmtdefs.label_font_size, gmtdefs.label_font, Ctrl->L.s, 0.0, 10, 0);
			if (!Ctrl->F.active) {	/* Draw scale bar */
				plot3d (off, -off, 3);
				plot3d ((max_radius - frame_info.axis[0].item[GMT_GRID_UPPER].interval) * Ctrl->S.scale, -off, 2);
				plot3d (off, -off, 3);
				plot3d (off, gmtdefs.tick_length - off, 2);
				plot3d ((max_radius - frame_info.axis[0].item[GMT_GRID_UPPER].interval) * Ctrl->S.scale, -off, 3);
				plot3d ((max_radius - frame_info.axis[0].item[GMT_GRID_UPPER].interval) * Ctrl->S.scale, gmtdefs.tick_length - off, -2);
				if (frame_info.axis[0].label[0]) {
					strcat (format, " %s");
					sprintf (text, format, frame_info.axis[0].item[GMT_GRID_UPPER].interval, frame_info.axis[0].label);
				}
				else
					sprintf (text, format, frame_info.axis[0].item[GMT_GRID_UPPER].interval);
				GMT_text3D ((max_radius - 0.5 * frame_info.axis[0].item[GMT_GRID_UPPER].interval) * Ctrl->S.scale, -(off + gmtdefs.annot_offset[0]), project_info.z_level, gmtdefs.annot_font_size[0], gmtdefs.annot_font[0], text, 0.0, 10, 0);
			}
			y = -(off + 5.0 * gmtdefs.annot_offset[0] + GMT_font[gmtdefs.annot_font[0]].height * lsize + GMT_font[gmtdefs.label_font].height * lsize);
			if (frame_info.axis[1].label[0]) GMT_text3D (0.0, y, project_info.z_level, gmtdefs.label_font_size, gmtdefs.label_font, frame_info.axis[1].label, 0.0, 10, 0);
		}
		GMT_text3D (off + 2.0 * gmtdefs.annot_offset[0], 0.0, project_info.z_level, gmtdefs.label_font_size, gmtdefs.label_font, Ctrl->L.e, 0.0, 5, 0);
		GMT_text3D (-off - 2.0 * gmtdefs.annot_offset[0], 0.0, project_info.z_level, gmtdefs.label_font_size, gmtdefs.label_font, Ctrl->L.w, 0.0, 7, 0);
		GMT_text3D (0.0, off + 2.0 * gmtdefs.annot_offset[0], project_info.z_level, gmtdefs.label_font_size, gmtdefs.label_font, Ctrl->L.n, 0.0, 2, 0);
		ps_setpaint (gmtdefs.background_rgb);
	}
	if (project_info.three_D) ps_rotatetrans (z_project.xmin, z_project.ymin + y_origin, 0.0);
	ps_rotatetrans (-x_origin, -y_origin, 0.0);

	GMT_plotend ();

	GMT_free ((void *)sum);
	GMT_free ((void *)xx);
	GMT_free ((void *)yy);
	GMT_free ((void *)azimuth);
	GMT_free ((void *)length);
	GMT_free ((void *)mode_length);
	GMT_free ((void *)mode_direction);

	Free_psrose_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void arc3d (double x, double y, double r, double a1, double a2, GMT_LONG status)
{
	if (project_info.three_D) {	/* Must project and draw line segments */
		GMT_LONG n;
		double *xx, *yy;

		n = GMT_get_arc (x, y, r, a1, a2, &xx, &yy);
		GMT_2D_to_3D (xx, yy, project_info.z_level, (GMT_LONG)n);
		ps_line (xx, yy, (GMT_LONG)n, status, FALSE);
		GMT_free ((void *)xx);
		GMT_free ((void *)yy);
	}
	else
		ps_arc (x, y, r, a1, a2, status);
}

void plot3d (double x, double y, GMT_LONG pen)
{
	if (project_info.three_D) GMT_xyz_to_xy (x, y, project_info.z_level, &x, &y);
	ps_plot (x, y, (int)pen);
}

void *New_psrose_Ctrl () {	/* Allocate and initialize a new control structure */
	struct PSROSE_CTRL *C;
	
	C = (struct PSROSE_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct PSROSE_CTRL), "New_psrose_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	C->E.azimuth = 180.0;
	C->E.elevation = 90.0;
	C->M.rgb[0] = C->M.rgb[1] = C->M.rgb[2] = 0;
	GMT_init_fill (&C->G.fill, -1, -1, -1);
	GMT_init_pen (&C->W.pen, GMT_PENWIDTH);
	C->M.thickness = 0.03;
	C->M.length = 0.12;
	C->M.width = 0.1;
	C->S.scale = 3.0;
	C->Z.scale = 1.0;
	if (gmtdefs.measure_unit == GMT_CM) {
		C->S.scale = 7.5 / 2.54;
		C->M.thickness = 0.075 / 2.54;
		C->M.length = 0.3 / 2.54;
		C->M.width = 0.25 / 2.54;
	}
	return ((void *)C);
}

void Free_psrose_Ctrl (struct PSROSE_CTRL *C) {	/* Deallocate control structure */
	if (C->C.file) free ((void *)C->C.file);	
	GMT_free ((void *)C);	
}
