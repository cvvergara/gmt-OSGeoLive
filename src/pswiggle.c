/*--------------------------------------------------------------------
 *	$Id: pswiggle.c 9923 2012-12-18 20:45:53Z pwessel $
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
 * pswiggle reads x,y,z from GMT_stdin and plots a wiggleplot using the
 * specified map projection. If the distance between 2 consecutive points
 * exceeds the value of dist_gap, a data gap is assumed.  The user may
 * select a preferred direction where s/he wants the positive anomalies
 * to be pointing towards.  Positive normal vectors to the trackline
 * will then always be within +- 90 degrees of this direction.
 * Separate colors may be specified for positive anomaly, outline of
 * anomaly, and trackline.  Plotting of the outline and track is optional.
 *
 * Author:	Paul Wessel
 * Date:	20-JUN-2000
 * Version:	4
 *
 */
 
#include "gmt.h"
#include "pslib.h"

struct PSWIGGLE_CTRL {
	struct A {	/* -A<azimuth> */
		GMT_LONG active;
		double value;
	} A;
	struct C {	/* -C<center> */
		GMT_LONG active;
		double value;
	} C;
	struct D {	/* -D[x]<gap> */
		GMT_LONG active;
		GMT_LONG mode;	/* 0 is km, 1 is projected inches */
		double value;
	} D;
	struct E {	/* -Eazim/elev */
		GMT_LONG active;
		double azimuth, elevation;
	} E;
	struct G {	/* -G<fill> */
		GMT_LONG active;
		struct GMT_FILL fill;
	} G;
	struct I {	/* -I<azimuth> */
		GMT_LONG active;
		double value;
	} I;
	struct N {	/* -N */
		GMT_LONG active;
	} N;
	struct S {	/* -S[x]<lon0>/<lat0>/<length>/<units> */
		GMT_LONG active;
		GMT_LONG cartesian;
		double lon, lat, length;
		char *label;
	} S;
	struct T {	/* -T<pen> */
		GMT_LONG active;
		struct GMT_PEN pen;
	} T;
	struct W {	/* -W<pen> */
		GMT_LONG active;
		struct GMT_PEN pen;
	} W;
	struct Z {	/* -Z<scale> */
		GMT_LONG active;
		double scale;
		char unit;
	} Z;
};

#define PSWIGGLE_MAX_POINTS 1000

int main (int argc, char **argv)
{
	GMT_LONG error = FALSE, nofile = TRUE, paint_wiggle, done;

	char line[BUFSIZ], *units = CNULL, txt_a[GMT_LONG_TEXT], txt_b[GMT_LONG_TEXT];

	GMT_LONG i, j, n, n_alloc = GMT_CHUNK;
	GMT_LONG k, n_files = 0, fno, n_args, n_expected_fields, n_fields, n_read, wantx, wanty;

	double west = 0.0, east = 360.0, south = 0.0, north = 0.0, x_2, y_2, d2m = 1.0;
	double dx, dy, dz, ds, *in = NULL, *xx = NULL, *yy = NULL, *zz = NULL, *lon = NULL, *lat = NULL, *z = NULL, start_az, stop_az, fix_az;

	FILE *fp = NULL;

	struct PSWIGGLE_CTRL *Ctrl = NULL;

	void GMT_draw_z_scale (double x0, double y0, double length, double zscale, GMT_LONG gave_xy, char *units);
	void plot_wiggle (double *x, double *y, double *z, GMT_LONG np, double zscale, double start_az, double stop_az, GMT_LONG fixed, double fix_az, struct GMT_FILL *fill, struct GMT_PEN *pen_o, struct GMT_PEN *pen_t, GMT_LONG paint_wiggle, GMT_LONG negative, GMT_LONG outline, GMT_LONG track);
	void *New_pswiggle_Ctrl (), Free_pswiggle_Ctrl (struct PSWIGGLE_CTRL *C);
	
	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct PSWIGGLE_CTRL *)New_pswiggle_Ctrl ();	/* Allocate and initialize a new control structure */

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch(argv[i][1]) {

				/* Common parameters */
	
				case 'B':
				case 'H':
				case 'J':
				case 'K':
				case 'M':
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
				case 'f':
				case 'm':
				case ':':
				case '\0':
					error += GMT_parse_common_options (argv[i], &west, &east, &south, &north);
					break;
		
				/* Supplemental parameters */
	
				case 'A':
					Ctrl->A.active = TRUE;
					Ctrl->A.value = atof (&argv[i][2]);
					break;
				case 'C':
					Ctrl->C.active = TRUE;
					Ctrl->C.value = atof (&argv[i][2]);
					break;
				case 'D':
					Ctrl->D.active = TRUE;
					if (argv[i][2] == 'x') {
						Ctrl->D.mode = 1;	/* Projected distances (e.g., inches) */
						Ctrl->D.value = GMT_convert_units (&argv[i][3], GMT_INCH);
					}
					else	
						Ctrl->D.value = atof (&argv[i][2]);	/* In km */
					break;
				case 'E':
					Ctrl->E.active = TRUE;
					error += GMT_get_proj3D (&argv[i][2], &Ctrl->E.azimuth, &Ctrl->E.elevation);
					break;
				case 'G':
					Ctrl->G.active = TRUE;
					if (GMT_getfill (&argv[i][2], &Ctrl->G.fill)) {
						GMT_fill_syntax ('G', " ");
						error++;
					}
					break;
				case 'I':
					Ctrl->I.value = atof (&argv[i][2]);
					Ctrl->I.active = TRUE;
					break;
				case 'N':
					Ctrl->N.active = TRUE;
					break;
				case 'T':
					Ctrl->T.active = TRUE;
					if (GMT_getpen (&argv[i][2], &Ctrl->T.pen)) {
						GMT_pen_syntax ('T', " ");
						error++;
					}
					break;
				case 'W':
					Ctrl->W.active = TRUE;
					if (GMT_getpen (&argv[i][2], &Ctrl->W.pen)) {
						GMT_pen_syntax ('W', " ");
						error++;
					}
					break;
				case 'S':
					Ctrl->S.active = TRUE;
					j = 0;
					if (argv[i][2] == 'x') Ctrl->S.cartesian = TRUE, j = 1;
					k = sscanf (&argv[i][2+j], "%[^/]/%[^/]/%lf", txt_a, txt_b, &Ctrl->S.length);
					wantx = (Ctrl->S.cartesian) ? GMT_IS_FLOAT : GMT_IS_LON;
					wanty = (Ctrl->S.cartesian) ? GMT_IS_FLOAT : GMT_IS_LAT;
					error += GMT_verify_expectations (wantx, GMT_scanf_arg (txt_a, wantx, &Ctrl->S.lon), txt_a);
					error += GMT_verify_expectations (wanty, GMT_scanf_arg (txt_b, wanty, &Ctrl->S.lat), txt_b);
					if ((units = strrchr (argv[i], '/'))) {
						units++;
						Ctrl->S.label = strdup (units);
					}
					if (k != 3) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -S option:  Correct syntax\n", GMT_program);
						fprintf (stderr, "	-S[x]<x0>/<y0>/<length>[/<units>]\n");
						error++;
					}
					break;
			
				case 'Z':
					Ctrl->Z.active = TRUE;
					j = strlen (argv[i]) - 1;
					if (strchr ("cimpCIMP", (int)argv[i][j])) Ctrl->Z.unit = argv[i][j];
					Ctrl->Z.scale = atof (&argv[i][2]);
					break;
			
				/* Illegal options */
	
				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else
			n_files++;
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr,"pswiggle %s - Plot xyz-series along tracks\n\n", GMT_VERSION);
		fprintf (stderr, "usage: pswiggle <xyz-files> %s %s -Z<scale>\n", GMT_J_OPT, GMT_Rgeo_OPT);
		fprintf (stderr, "\t[-A<azimuth>] [%s] [-C<center>] [-D[x]<gap>] [%s] [-G<fill>] [%s]\n", GMT_B_OPT, GMT_E_OPT, GMT_Ho_OPT);
		fprintf (stderr, "\t[-I<az>] -K [-N] [-O] [-P] [-S[x]<lon0>/<lat0>/<length>/<units>] [-T<trackpen>]\n");
		fprintf (stderr, "\t[%s] [-V] [-W<outlinepen>] [%s] [%s]\n\t[%s] [%s] [%s] [%s]] [%s]\n", GMT_U_OPT, GMT_X_OPT, GMT_Y_OPT, GMT_c_OPT, GMT_t_OPT, GMT_bi_OPT, GMT_f_OPT, GMT_mo_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\t<xyz_files> is one or more files.  If none, read standard input.\n");
		GMT_explain_option ('j');
		GMT_explain_option ('R');
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-A set azimuth for preferred positive wiggle orientation [0.0 (north)].\n");
		GMT_explain_option ('b');
		fprintf (stderr, "\t-C sets center value to be removed from z before plotting [0].\n");
		fprintf (stderr, "\t-D means there is a datagap if 2 points are more than <gap> distance units apart.\n");
		fprintf (stderr, "\t   For map projections, <gap> is assumed to be in km.\n");
		fprintf (stderr, "\t   Use -Dx to specify gaps in projected distances (append unit c|i|m|p).\n");
		fprintf (stderr, "\t   Otherwise, user-units are assumed.\n");
		GMT_explain_option ('E');
		GMT_fill_syntax ('G', "Specify color/pattern for positive areas.");
		GMT_explain_option ('H');
		fprintf (stderr, "\t-I set fixed projection azimuths for wiggles.\n");
		GMT_explain_option ('K');
		fprintf (stderr, "\t-N Fill negative wiggles instead [Default is positive].\n");
		GMT_explain_option ('O');
		GMT_explain_option ('P');
		fprintf (stderr, "\t-S draws a simple vertical scale centered on <lon0>/<lat0>.  Use -Sx to specify cartesian coordinates instead.\n");
		fprintf (stderr, "\t   <length> is in z-units, append unit name for labeling.\n");
		fprintf (stderr, "\t-T specifies track pen attributes. [Default is no track].\n");
		GMT_explain_option ('U');
		GMT_explain_option ('V');
		GMT_pen_syntax ('W', "specifies outline pen attributes [Default is no outline].");
		GMT_explain_option ('X');
		fprintf (stderr, "\t-Z gives the wiggle scale in data-units per %s.\n", GMT_unit_names[gmtdefs.measure_unit]);
		GMT_explain_option ('c');
		GMT_explain_option (':');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf (stderr, "\t   Default is 3 input columns.\n");
		GMT_explain_option ('f');
		GMT_explain_option ('m');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	if (!project_info.region_supplied) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify -R option\n", GMT_program);
		error++;
	}
	if (!(Ctrl->W.active || Ctrl->G.active)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify at least one of -G, -W\n", GMT_program);
		error++;
	}
	if (Ctrl->E.elevation <= 0.0 || Ctrl->E.elevation > 90.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -E option:  Elevation must be in 0-90 range\n", GMT_program);
		error++;
	}
	if (Ctrl->Z.scale == 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -Z option:  scale must be nonzero\n", GMT_program);
		error++;
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
        if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = 3;
        if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] < 3) {
                fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data (-bi) must have at least 3 columns\n", GMT_program);
		error++;
	}
	      
	if (error) exit (EXIT_FAILURE);

	if (GMT_io.binary[GMT_IN] && gmtdefs.verbose) {
		char *type[2] = {"double", "single"};
		fprintf (stderr, "%s: Expects %ld-column %s-precision binary data\n", GMT_program, GMT_io.ncol[GMT_IN], type[GMT_io.single_precision[GMT_IN]]);
	}

	z_project.view_azimuth = Ctrl->E.azimuth;
	z_project.view_elevation = Ctrl->E.elevation;

	GMT_err_fail (GMT_map_setup (west, east, south, north), "");

	if (GMT_IS_MAPPING) {
		Ctrl->D.value *= 1000.0;    /* Distance is now in meters */
		d2m = TWO_PI * project_info.EQ_RAD / 360.0;
	}

	GMT_plotinit (argc, argv);

	if (project_info.three_D) ps_transrotate (-z_project.xmin, -z_project.ymin, 0.0);

	Ctrl->Z.scale = 1.0 / Ctrl->Z.scale;

	switch (Ctrl->Z.unit) {	/* Adjust for possible unit selection */
		case 'C':
		case 'c':
			Ctrl->Z.scale *= GMT_u2u[GMT_CM][GMT_INCH];
			break;
		case 'I':
		case 'i':
			Ctrl->Z.scale *= GMT_u2u[GMT_INCH][GMT_INCH];
			break;
		case 'M':
		case 'm':
			Ctrl->Z.scale *= GMT_u2u[GMT_M][GMT_INCH];
			break;
		case 'P':
		case 'p':
			Ctrl->Z.scale *= GMT_u2u[GMT_PT][GMT_INCH];
			break;
		default:
			Ctrl->Z.scale *= GMT_u2u[gmtdefs.measure_unit][GMT_INCH];
			break;
	}

	/* Now convert angles to radians and keep it that way */
	start_az = (Ctrl->A.value - 90.0) * D2R;
	stop_az  = (Ctrl->A.value + 90.0) * D2R;
	fix_az = Ctrl->I.value * D2R;

	GMT_map_clip_on (GMT_no_rgb, 3);

	/* Allocate memory */

	lon = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (double), GMT_program);
	lat = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (double), GMT_program);
	z   = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (double), GMT_program);
	xx  = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (double), GMT_program);
	yy  = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (double), GMT_program);
	zz  = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (double), GMT_program);

	if (n_files != 0) nofile = FALSE;
	done = FALSE;
	n_args = (argc > 1) ? argc : 2;

	for (fno = 1; !done && fno < n_args; fno++) {	/* Loop over all the files */
		if (!nofile && argv[fno][0] == '-') continue;
		if (nofile) {
			fp = GMT_stdin;
			done = TRUE;
			if (gmtdefs.verbose) fprintf (stderr, "%s: Reading data from standard input\n", GMT_program);
#ifdef SET_IO_MODE
			GMT_setmode (GMT_IN);
#endif
		}
		else if ((fp = GMT_fopen (argv[fno], GMT_io.r_mode)) == NULL) {
			fprintf (stderr, "%s: Cannot open file %s\n", GMT_program, argv[fno]);
			continue;
		}

		if (!nofile && gmtdefs.verbose) {
			fprintf (stderr, "%s: Working on file %s\n", GMT_program, argv[fno]);
			sprintf (line, "File %s", argv[fno]);
			ps_comment (line);
		}

		if (GMT_io.io_header[GMT_IN]) for (i = 0; i < GMT_io.n_header_recs; i++) GMT_fgets (line, BUFSIZ, fp);
		if (GMT_io.multi_segments[GMT_IN]) {
			GMT_fgets (line, BUFSIZ, fp);
			if (gmtdefs.verbose) ps_comment (line);
		}
		n_expected_fields = (GMT_io.ncol[GMT_IN]) ? GMT_io.ncol[GMT_IN] : 3;

		n_fields = GMT_input (fp, &n_expected_fields, &in);
		while (! (GMT_io.status & GMT_IO_EOF)) {	/* Not yet EOF */

			while (GMT_io.status & GMT_IO_SEGMENT_HEADER) {
				if (gmtdefs.verbose) ps_comment (GMT_io.segment_header);
				n_fields = GMT_input (fp, &n_expected_fields, &in);
			}
			if (GMT_io.status & GMT_IO_EOF) continue;	/* At EOF */

			n = n_read = 0;
			while (! (GMT_io.status & (GMT_IO_SEGMENT_HEADER | GMT_IO_EOF))) {	/* Keep going until FALSE or = 2 segment header */
				n_read++;

				if (GMT_io.status & GMT_IO_MISMATCH) {
					fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld\n", GMT_program, n_fields, n_expected_fields, n_read);
					continue;
				}

				lon[n] = in[GMT_X];
				lat[n] = in[GMT_Y];
				z[n] = in[GMT_Z] - Ctrl->C.value;
				n++;
				if (n >= n_alloc) {
					n_alloc <<= 1;
					lon = (double *) GMT_memory ((void *)lon, (size_t)n_alloc, sizeof (double), GMT_program);
					lat = (double *) GMT_memory ((void *)lat, (size_t)n_alloc, sizeof (double), GMT_program);
					z = (double *) GMT_memory ((void *)z, (size_t)n_alloc, sizeof (double), GMT_program);
					xx = (double *) GMT_memory ((void *)xx, (size_t)n_alloc, sizeof (double), GMT_program);
					yy = (double *) GMT_memory ((void *)yy, (size_t)n_alloc, sizeof (double), GMT_program);
					zz = (double *) GMT_memory ((void *)zz, (size_t)n_alloc, sizeof (double), GMT_program);
				}

				n_fields = GMT_input (fp, &n_expected_fields, &in);
			}
			
			GMT_geo_to_xy (lon[0], lat[0], &xx[0], &yy[0]);
			zz[0] = z[0];
			paint_wiggle = (Ctrl->G.active && ((Ctrl->N.active && z[0] <= 0.0) || (!Ctrl->N.active && z[0] >= 0.0)));
			j = 1;
			for (i = 1; i < n; i++) {	/* Convert to inches/cm and get distance increments */
				
				GMT_geo_to_xy (lon[i], lat[i], &x_2, &y_2);

				if (Ctrl->D.mode)	/* Check for gap using projected units */
					ds = hypot (x_2 - xx[j-1], y_2 - yy[j-1]);
				else {	/* Distance in meters */
					dx = lon[i] - lon[i-1];
					dy = lat[i] - lat[i-1];
					if (GMT_IS_MAPPING) {	/* Distance in meters */
						dx *= cosd (0.5 * (lat[i] + lat[i-1]));
						ds = d2m * hypot (dx, dy);
					}
					else	/* Plain Cartesian distances of the users coordinates */
						ds = hypot (dx, dy);
				}
			
				if (j > 0 && (ds > Ctrl->D.value || GMT_is_dnan (z[i]))) {	/* Data gap, plot what we have */
					paint_wiggle = (Ctrl->G.active && ((Ctrl->N.active && zz[j-1] <= 0.0) || (!Ctrl->N.active && zz[j-1] >= 0.0)));
					plot_wiggle (xx, yy, zz, j, Ctrl->Z.scale, start_az, stop_az, Ctrl->I.active, fix_az, &Ctrl->G.fill, &Ctrl->W.pen, &Ctrl->T.pen, paint_wiggle, Ctrl->N.active, Ctrl->W.active, Ctrl->T.active);
					j = 0;
				}
				else if (j >= PSWIGGLE_MAX_POINTS) {
					paint_wiggle = (Ctrl->G.active && ((Ctrl->N.active && zz[j-1] <= 0.0) || (!Ctrl->N.active && zz[j-1] >= 0.0)));
					plot_wiggle (xx, yy, zz, j, Ctrl->Z.scale, start_az, stop_az, Ctrl->I.active, fix_az, &Ctrl->G.fill, &Ctrl->W.pen, &Ctrl->T.pen, paint_wiggle, Ctrl->N.active, Ctrl->W.active, Ctrl->T.active);
					xx[0] = xx[j-1];
					yy[0] = yy[j-1];
					zz[0] = zz[j-1];
					j = 1;
				}
				else if (!GMT_is_dnan (z[i-1]) && (z[i]*z[i-1] < 0.0 || z[i] == 0.0)) {	/* Crossed 0, add new point and plot */
					dz = z[i] - z[i-1];
					xx[j] = (dz == 0.0) ? xx[j-1] : xx[j-1] + fabs (z[i-1] / dz) * (x_2 - xx[j-1]);
					yy[j] = (dz == 0.0) ? yy[j-1] : yy[j-1] + fabs (z[i-1] / dz) * (y_2 - yy[j-1]);
					zz[j] = 0.0;
					j++;
					paint_wiggle = (Ctrl->G.active && ((Ctrl->N.active && zz[j-2] <= 0.0) || (!Ctrl->N.active && zz[j-2] >= 0.0)));
					plot_wiggle (xx, yy, zz, j, Ctrl->Z.scale, start_az, stop_az, Ctrl->I.active, fix_az, &Ctrl->G.fill, &Ctrl->W.pen, &Ctrl->T.pen, paint_wiggle, Ctrl->N.active, Ctrl->W.active, Ctrl->T.active);
					xx[0] = xx[j-1];
					yy[0] = yy[j-1];
					zz[0] = zz[j-1];
					j = 1;
				}
				xx[j] = x_2;
				yy[j] = y_2;
				zz[j] = z[i];
				if (!GMT_is_dnan (z[i])) j++;
			}
	
			if (j > 1) {
				paint_wiggle = (Ctrl->G.active && ((Ctrl->N.active && zz[j-1] <= 0.0) || (!Ctrl->N.active && zz[j-1] >= 0.0)));
				plot_wiggle (xx, yy, zz, j, Ctrl->Z.scale, start_az, stop_az, Ctrl->I.active, fix_az, &Ctrl->G.fill, &Ctrl->W.pen, &Ctrl->T.pen, paint_wiggle, Ctrl->N.active, Ctrl->W.active, Ctrl->T.active);
			}
		}
		if (fp != GMT_stdin) GMT_fclose (fp);
	}
	
	GMT_map_clip_off();
	GMT_map_basemap ();

	if (Ctrl->S.active) GMT_draw_z_scale (Ctrl->S.lon, Ctrl->S.lat, Ctrl->S.length, Ctrl->Z.scale, Ctrl->S.cartesian, Ctrl->S.label);

	if (project_info.three_D) ps_rotatetrans (z_project.xmin, z_project.ymin, 0.0);
	GMT_plotend ();

	GMT_free ((void *)lon);
	GMT_free ((void *)lat);
	GMT_free ((void *)z);
	GMT_free ((void *)xx);
	GMT_free ((void *)yy);
	GMT_free ((void *)zz);

	Free_pswiggle_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void plot_wiggle (double *x, double *y, double *z, GMT_LONG np, double zscale, double start_az, double stop_az, GMT_LONG fixed, double fix_az, struct GMT_FILL *fill, struct GMT_PEN *pen_o, struct GMT_PEN *pen_t, GMT_LONG paint_wiggle, GMT_LONG negative, GMT_LONG outline, GMT_LONG track)
{
	double dx, dy, len, az = 0.0, s = 0.0, c = 0.0, x_inc, y_inc;
	GMT_LONG n = 0, i;
	static char *name[2] = {"Positive", "Negative"}, comment[GMT_LONG_TEXT];

	if (fixed) {
		az = fix_az;
		sincos (az, &s, &c);
	}

	if (paint_wiggle || outline) {
		if (!GMT_n_alloc) GMT_get_plot_array ();
		GMT_x_plot[0] = x[0];
		GMT_y_plot[0] = y[0];
		n = 1;
		for (i = 0; i < np; i++) {
			if (!fixed && i < np-1) {
				dx = x[i+1] - x[i];
				dy = y[i+1] - y[i];
				if (!(dx == 0.0 && dy == 0.0)) az = -d_atan2 (dy, dx) - TWO_PI;	/* Azimuth of normal to track */
				while (az < start_az) az += TWO_PI;
				if (az > stop_az) az -= M_PI;
			}
			if (fabs (z[i]) > 0.0) {
				if (!fixed) sincos (az, &s, &c);
				len = zscale * z[i];
				x_inc = len * s;
				y_inc = len * c;
			}
			else
				x_inc = y_inc = 0.0;
		
			GMT_x_plot[n] = x[i] + x_inc;
			GMT_y_plot[n] = y[i] + y_inc;
			n++;
			if (n == GMT_n_alloc) GMT_get_plot_array ();
		}
		GMT_x_plot[n] = x[np-1];
		GMT_y_plot[n] = y[np-1];
		n++;

		if (paint_wiggle) {
			for (i = np - 2; i >= 0; i--, n++) {	/* Go back to 1st point along track */
				if (n == GMT_n_alloc) GMT_get_plot_array ();
				GMT_x_plot[n] = x[i];
				GMT_y_plot[n] = y[i];
			}
		}
		if (project_info.three_D) GMT_2D_to_3D (GMT_x_plot, GMT_y_plot, project_info.z_level, n);
	}


	if (paint_wiggle) { /* First shade wiggles */
		sprintf (comment, "%s wiggle", name[negative]);
		ps_comment (comment);
		GMT_fill (GMT_x_plot, GMT_y_plot, n, fill, FALSE);
	}

	if (outline) { /* Then draw wiggle outline */
		sprintf (comment, "Wiggle line");
		ps_comment (comment);
		GMT_setpen (pen_o);
		ps_line (&GMT_x_plot[1], &GMT_y_plot[1], np, 3, FALSE);
	}

	if (track) {	/* Finally draw track line */
		if (project_info.three_D) GMT_2D_to_3D (x, y, project_info.z_level, np);
		sprintf (comment, "Track line");
		ps_comment (comment);
		GMT_setpen (pen_t);
		ps_line (x, y, np, 3, FALSE);
	}
}

void GMT_draw_z_scale (double x0, double y0, double length, double zscale, GMT_LONG gave_xy, char *units)
{
	double dy, off, xx[4], yy[4];
	char txt[GMT_LONG_TEXT];

	GMT_setpen (&gmtdefs.tick_pen);

	if (!gave_xy) {
		GMT_geo_to_xy (x0, y0, &xx[0], &yy[0]);
		x0 = xx[0];	y0 = yy[0];
	}

	if (units)
		sprintf (txt, "%g %s", length, units);
	else
		sprintf (txt, "%g", length);
	dy = 0.5 * length * zscale;
	GMT_xyz_to_xy (x0 + gmtdefs.map_scale_height, y0 - dy, 0.0, &xx[0], &yy[0]);
	GMT_xyz_to_xy (x0, y0 - dy, 0.0, &xx[1], &yy[1]);
	GMT_xyz_to_xy (x0, y0 + dy, 0.0, &xx[2], &yy[2]);
	GMT_xyz_to_xy (x0 + gmtdefs.map_scale_height, y0 + dy, 0.0, &xx[3], &yy[3]);
	ps_line (xx, yy, (GMT_LONG)4, 3, FALSE);
	off = ((gmtdefs.map_scale_height > 0.0) ? gmtdefs.tick_length : 0.0) + gmtdefs.annot_offset[0];
	GMT_text3D (x0 + off, y0, project_info.z_level, gmtdefs.annot_font_size[0], gmtdefs.annot_font[0], txt, 0.0, 5, 0);
}

void *New_pswiggle_Ctrl () {	/* Allocate and initialize a new control structure */
	struct PSWIGGLE_CTRL *C;
	
	C = (struct PSWIGGLE_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct PSWIGGLE_CTRL), "New_pswiggle_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	C->D.value = DBL_MAX;
	C->E.azimuth = 180.0;
	C->E.elevation = 90.0;
	GMT_init_pen (&C->W.pen, GMT_PENWIDTH);
	GMT_init_pen (&C->T.pen, GMT_PENWIDTH);
	GMT_init_fill (&C->G.fill, gmtdefs.basemap_frame_rgb[0], gmtdefs.basemap_frame_rgb[1], gmtdefs.basemap_frame_rgb[2]);
		
	return ((void *)C);
}

void Free_pswiggle_Ctrl (struct PSWIGGLE_CTRL *C) {	/* Deallocate control structure */
	if (C->S.label) free ((void *)C->S.label);	
	GMT_free ((void *)C);	
}
