/*--------------------------------------------------------------------
 *	$Id: grdvector.c 10216 2014-02-27 19:24:21Z pwessel $
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
   grdvector reads 2 grid files that contains the 2 components of a vector
   field (cartesian or polar) and plots vectors at the grid positions.
   This is basically a short-hand for using grd2xyz | psxy -SV and is
   more convenient for such plots on a grid.

   Author:	Paul Wessel
   Date:	12-JUN-1995
   Revised:	15-FEB-2000
   Version:	4
 
 */


#include "gmt.h"
#include "pslib.h"

struct GRDVECTOR_CTRL {
	struct A {	/* -A */
		GMT_LONG active;
	} A;
	struct C {	/* -C<cpt> */
		GMT_LONG active;
		char *file;
	} C;
	struct E {	/* -E */
		GMT_LONG active;
	} E;
	struct G {	/* -G<fill> */
		GMT_LONG active;
		struct GMT_FILL fill;
	} G;
	struct I {	/* -Idx[/dy] */
		GMT_LONG active;
		double xinc, yinc;
	} I;
	struct N {	/* -N */
		GMT_LONG active;
	} N;
	struct Q {	/* -Q<params> */
		GMT_LONG active;
		double width, length, thickness, norm;
	} Q;
	struct S {	/* -S[l]<scale>][<unit>] */
		GMT_LONG active;
		GMT_LONG constant;
		char unit;
		double factor;
	} S;
	struct T {	/* -T */
		GMT_LONG active;
	} T;
	struct W {	/* -W<pen> */
		GMT_LONG active;
		struct GMT_PEN pen;
	} W;
	struct Z {	/* -Z */
		GMT_LONG active;
	} Z;
};

int main (int argc, char **argv)
{
	GMT_LONG	i, j, n = 0, nx, ny, i0, j0, di, dj;

	GMT_LONG ij, nm;
	
	GMT_LONG shrink_properties = FALSE, error = FALSE;

	char *file[2], txt_a[GMT_LONG_TEXT], txt_b[GMT_LONG_TEXT], txt_c[GMT_LONG_TEXT];

	float *r = NULL, *theta = NULL;

	double v_w, h_l, h_w, v_shrink = 1.0, tmp, x, y, plot_x, plot_y, x_off, y_off;
	double west, east, south, north, x2, y2;
	double data_west, data_east, data_south, data_north, value, c, s;

	struct GRD_HEADER h[2];
	struct GRDVECTOR_CTRL *Ctrl = NULL;

	void *New_grdvector_Ctrl (), Free_grdvector_Ctrl (struct GRDVECTOR_CTRL *C);

	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct GRDVECTOR_CTRL *)New_grdvector_Ctrl ();	/* Allocate and initialize a new control structure */
	
	west = east = south = north = 0.0;
	di = dj = 1;
	i0 = j0 = 0;

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				/* Common parameters */

				case 'B':
				case 'J':
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
				case 'c':
				case 'f':
				case '\0':
					error += GMT_parse_common_options (argv[i], &west, &east, &south, &north);
					break;

				/* Supplemental parameters */

				case 'A':
					Ctrl->A.active = TRUE;
					break;
				case 'C':	/* Vary symbol color with z */
					Ctrl->C.active = TRUE;
					Ctrl->C.file = strdup (&argv[i][2]);
					break;
				case 'E':
					Ctrl->E.active = TRUE;
					break;
				case 'G':		/* Set Gray shade for polygon */
					Ctrl->G.active = TRUE;
					if (GMT_getfill (&argv[i][2], &Ctrl->G.fill)) {
						GMT_fill_syntax ('G', " ");
						error++;
					}
					break;
				case 'I':	/* Only use gridnodes Ctrl->I.xinc,Ctrl->I.yinc apart */
					Ctrl->I.active = TRUE;
					if (GMT_getinc (&argv[i][2], &Ctrl->I.xinc, &Ctrl->I.yinc)) {
						GMT_inc_syntax ('I', 1);
						error = TRUE;
					}
					break;
				case 'N':	/* Do not clip at border */
					Ctrl->N.active = TRUE;
					break;
				case 'Q':
					Ctrl->Q.active = TRUE;
					for (j = 2; argv[i][j] && argv[i][j] != 'n'; j++);
					if (argv[i][j]) {	/* Normalize option used */
						Ctrl->Q.norm = GMT_convert_units (&argv[i][j+1], GMT_INCH);
						if (Ctrl->Q.norm <= 0.0) {
							fprintf (stderr, "%s: GMT SYNTAX ERROR -Qn option:  No reference length given\n", GMT_program);
							error++;
						}
						argv[i][j] = '\0';	/* Temporarily chop of the n<norm> string */
					}
					if (argv[i][2] && argv[i][3] != 'n') {	/* We specified the three parameters */
						if (sscanf (&argv[i][2], "%[^/]/%[^/]/%s", txt_a, txt_b, txt_c) != 3) {
							fprintf (stderr, "%s: GMT SYNTAX ERROR -Q option:  Could not decode arrowwidth/headlength/headwidth\n", GMT_program);
							error++;
						}
						else {
							Ctrl->Q.thickness  = GMT_convert_units (txt_a, GMT_INCH);
							Ctrl->Q.length = GMT_convert_units (txt_b, GMT_INCH);
							Ctrl->Q.width  = GMT_convert_units (txt_c, GMT_INCH);
						}
					}
					if (Ctrl->Q.norm > 0.0) argv[i][j] = 'n';	/* Restore the n<norm> string */
					break;
				case 'S':
					Ctrl->S.active = TRUE;
					j = (int)strlen (argv[i]) - 1;
					if (strchr ("cimpCIMP", (int)argv[i][j]))	/* Recognized unit character */
						Ctrl->S.unit = argv[i][j];
					else if (! (argv[i][j] == '.' || isdigit ((int)argv[i][j]))) {	/* Not decimal point or digit means trouble */
						fprintf (stderr, "%s: GMT SYNTAX ERROR -S option:  Unrecognized unit %c\n", GMT_program, argv[i][j]);
						error++;
					}
					if (argv[i][2] == 'l' || argv[i][2] == 'L') {
						Ctrl->S.constant = TRUE;
						Ctrl->S.factor = atof (&argv[i][3]);
					}
					else
						Ctrl->S.factor = atof (&argv[i][2]);
					break;
				case 'T':
					Ctrl->T.active = TRUE;
					break;
				case 'W':		/* Set line attributes */
					Ctrl->W.active = TRUE;
					if (argv[i][2] && GMT_getpen (&argv[i][2], &Ctrl->W.pen)) {
						GMT_pen_syntax ('W', " ");
						error++;
					}
					break;
				case 'Z':
					Ctrl->Z.active = TRUE;
					break;
				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else if (n < 2)
			file[n++] = argv[i];
		else
			n++;
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "grdvector %s - Plot vector fields from grid files\n\n", GMT_VERSION);
		fprintf (stderr, "usage: grdvector compx.grd compy.grd %s %s [-A]\n", GMT_J_OPT, GMT_Rgeo_OPT);
		fprintf (stderr, "\t[%s] [-C<cpt>] [-E] [-G<fill>] [-I<dx/dy>] [-K] [-N] [-O] [-P] [-Q<params>] [-S[l]<scale>[<unit>]]\n", GMT_B_OPT);
		fprintf (stderr, "\t[-T] [%s] [-V] [-W<pen>] [%s] [%s] [-Z] [%s]\n\n", GMT_U_OPT, GMT_X_OPT, GMT_Y_OPT, GMT_c_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\tcompx & compy are grid files with the 2 vector components.\n");
		GMT_explain_option ('j');
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-A means grids have polar (r, theta) components [Default is Cartesian].\n");
		GMT_explain_option ('b');
		fprintf (stderr, "\t-C Use cpt-file to assign colors based on vector length.\n");
		fprintf (stderr, "\t-E cEnter vectors on grid nodes [Default draws from grid node].\n");
		GMT_fill_syntax ('G', "Select vector fill [Default is outlines only].");
		fprintf (stderr, "\t-I plots only those nodes that are <dx/dy> apart [Default is all nodes].\n");
		GMT_explain_option ('K');
		fprintf (stderr, "\t-N Do Not clip vectors that exceed the map boundaries [Default will clip].\n");
		GMT_explain_option ('O');
		GMT_explain_option ('P');
		fprintf (stderr, "\t-Q Select vector plot [Default is stick-plot].\n");
		fprintf (stderr, "\t   Optionally, specify vector parameters\n");
		fprintf (stderr, "\t   <params> are arrowwidth/headlength/headwidth [Default is 0.03i/0.12i/0.09i].\n");
		fprintf (stderr, "\t   Append n<size>[unit] which will cause vectors shorter than <size> to be\n");
		fprintf (stderr, "\t     scaled down.\n");
		GMT_explain_option ('R');
		fprintf (stderr, "\t-S sets scale for vector length in data units per %s [1].\n", GMT_unit_names[gmtdefs.measure_unit]);
		fprintf (stderr, "\t   Append c, i, m, or p to indicate cm, inch, m, or points as the distance unit.\n");
		fprintf (stderr, "\t   Alternatively, prepend l to indicate a fixed length for all vectors.\n");
		fprintf (stderr, "\t-T means azimuth should be converted to angles based on map projection.\n");
		GMT_explain_option ('U');
		GMT_explain_option ('V');
		GMT_pen_syntax ('W', "sets pen attributes.");
		fprintf (stderr, "\t   Default pen attributes [width = %gp, color = (%d/%d/%d), solid line].\n", 
			Ctrl->W.pen.width, Ctrl->W.pen.rgb[0], Ctrl->W.pen.rgb[1], Ctrl->W.pen.rgb[2]);
		GMT_explain_option ('X');
		fprintf (stderr, "\t-Z means the angles provided are azimuths rather than direction.\n");
		GMT_explain_option ('c');
		GMT_explain_option ('f');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	GMT_check_lattice (&Ctrl->I.xinc, &Ctrl->I.yinc, NULL, &Ctrl->I.active);

	if (Ctrl->I.active && (Ctrl->I.xinc <= 0.0 || Ctrl->I.yinc <= 0.0)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -I option:  Must specify positive increments\n", GMT_program);
		error++;
	}
	if (Ctrl->S.factor == 0.0 && !Ctrl->S.constant) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -S option:  Scale must be nonzero\n", GMT_program);
		error++;
	}
	if (Ctrl->S.factor <= 0.0 && Ctrl->S.constant) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -Sl option:  Length must be positive\n", GMT_program);
		error++;
	}
	if (Ctrl->Q.active && Ctrl->Q.norm > 0.0) {
		v_shrink = 1.0 / Ctrl->Q.norm;
		shrink_properties = TRUE;
	}
	if (Ctrl->S.constant && shrink_properties) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -Sl, -Q options:  Cannot use -Q..n<size> with -Sl\n", GMT_program);
		error++;
	}
	if (Ctrl->Z.active && !Ctrl->A.active) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -Z option:  Azimuths not valid input for Cartesian data\n", GMT_program);
		error++;
	}
	if (Ctrl->C.active && !Ctrl->C.file) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -C option:  Must specify a color palette table\n", GMT_program);
		error++;
	}
	if (!(Ctrl->G.active || Ctrl->W.active || Ctrl->C.active)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify at least one of -G, -W, -C\n", GMT_program);
		error++;
	}
	if (n != 2) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify two input grid files\n", GMT_program);
		error++;
	}

	if (error) exit (EXIT_FAILURE);

	if (Ctrl->C.active) GMT_read_cpt (Ctrl->C.file);

	if (!(strcmp (file[0], "=") || strcmp (file[1], "="))) {
		fprintf (stderr, "%s: Piping of grid files not supported!\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	for (i = 0; i < 2; i++) {
		GMT_grd_init (&h[i], argc, argv, FALSE);
		GMT_err_fail (GMT_read_grd_info (file[i], &h[i]), file[i]);
	}

	if (!(h[0].nx == h[1].nx && h[0].ny == h[1].ny && h[0].x_min == h[1].x_min && h[0].y_min == h[1].y_min 
		&& h[0].x_inc == h[1].x_inc && h[0].y_inc == h[1].y_inc)) {
		fprintf (stderr, "%s: files %s and %s does not match!\n", GMT_program, file[0], file[1]);
		exit (EXIT_FAILURE);
	}

	/* Determine what wesn to pass to map_setup */

	if (!project_info.region_supplied) {
		west = h[0].x_min;
		east = h[0].x_max;
		south = h[0].y_min;
		north = h[0].y_max;
	}

	GMT_err_fail (GMT_map_setup (west, east, south, north), "");

	/* Determine the wesn to be used to read the grid file */

	if (GMT_grd_setregion (&h[0], &data_west, &data_east, &data_south, &data_north, BCR_BILINEAR)) {
		/* No grid to plot; just do empty map and exit */
		GMT_plotinit (argc, argv);
		GMT_map_basemap ();
		GMT_plotend ();
		GMT_end (argc, argv);
		exit (EXIT_SUCCESS);
	}

	/* Read data */

	nx = GMT_get_n (data_west, data_east, h[0].x_inc, h[0].node_offset);
	ny = GMT_get_n (data_south, data_north, h[0].y_inc, h[0].node_offset);
	nm = GMT_get_nm (nx, ny);
	r = (float *) GMT_memory (VNULL, (size_t)nm, sizeof (float), GMT_program);
	theta = (float *) GMT_memory (VNULL, (size_t)nm, sizeof (float), GMT_program);
	GMT_err_fail (GMT_read_grd (file[0], &h[0], r, data_west, data_east, data_south, data_north, GMT_pad, FALSE), file[0]);
	GMT_err_fail (GMT_read_grd (file[1], &h[1], theta, data_west, data_east, data_south, data_north, GMT_pad, FALSE), file[1]);

	if (!Ctrl->S.constant) Ctrl->S.factor = 1.0 / Ctrl->S.factor;

	switch (Ctrl->S.unit) {	/* Adjust for possible unit selection */
		case 'C':
		case 'c':
			Ctrl->S.factor *= GMT_u2u[GMT_CM][GMT_INCH];
			break;
		case 'I':
		case 'i':
			Ctrl->S.factor *= GMT_u2u[GMT_INCH][GMT_INCH];
			break;
		case 'M':
		case 'm':
			Ctrl->S.factor *= GMT_u2u[GMT_M][GMT_INCH];
			break;
		case 'P':
		case 'p':
			Ctrl->S.factor *= GMT_u2u[GMT_PT][GMT_INCH];
			break;
		default:
			Ctrl->S.factor *= GMT_u2u[gmtdefs.measure_unit][GMT_INCH];
			break;
	}

	GMT_plotinit (argc, argv);

	GMT_setpen (&Ctrl->W.pen);

        if (!Ctrl->N.active) GMT_map_clip_on (GMT_no_rgb, 3);

	if (Ctrl->I.xinc != 0.0 && Ctrl->I.yinc != 0.0) {	/* Gave a coarser grid spacing, we hope */
		double val = Ctrl->I.yinc / h[0].y_inc;
		dj = irint (val);
		if (dj == 0 || fabs (val - dj) > GMT_CONV_LIMIT) {
			fprintf (stderr, "%s: Error: New y grid increment (%g) is not a multiple of actual grid increment (%g).\n", GMT_program, Ctrl->I.xinc, h[0].x_inc);
			exit (EXIT_FAILURE);
		}
		val = Ctrl->I.xinc / h[0].x_inc;
		di = irint (val);
		if (di == 0 || fabs (val - di) > GMT_CONV_LIMIT) {
			fprintf (stderr, "%s: Error: New x grid increment (%g) is not a multiple of actual grid increment (%g).\n", GMT_program, Ctrl->I.xinc, h[0].x_inc);
			exit (EXIT_FAILURE);
		}
		/* Determine starting point for straddled access */
		tmp = ceil (h[0].y_max / Ctrl->I.yinc) * Ctrl->I.yinc;
		if (tmp > h[0].y_max) tmp -= Ctrl->I.yinc;
		j0 = irint ((h[0].y_max - tmp) / h[0].y_inc);
		tmp = floor (h[0].x_min / Ctrl->I.xinc) * Ctrl->I.xinc;
		if (tmp < h[0].x_min) tmp += Ctrl->I.xinc;
		i0 = irint ((tmp - h[0].x_min) / h[0].x_inc);
	}

	for (j = j0; j < h[1].ny; j += dj) {
		y = GMT_j_to_y (j, h[0].y_min, h[0].y_max, h[0].y_inc, h[0].xy_off, h[0].ny);
		for (i = i0; i < h[1].nx; i += di) {

			ij = GMT_IJ (j, i, h[0].nx);
			if (GMT_is_fnan (r[ij]) || GMT_is_fnan (theta[ij])) continue;	/* Cannot plot NaN-vectors */

			value = r[ij];

			if (!Ctrl->A.active) {
				value = hypot (theta[ij], r[ij]);
				if (value == 0.0) continue;
				theta[ij] = (float)(atan2d (theta[ij], r[ij]));
				r[ij] = (float)value;
			}
			else if (r[ij] < 0.0) {
				r[ij] = -r[ij];
				theta[ij] += 180.0;
			}
			else if (r[ij] == 0.0) continue;

			if (Ctrl->C.active) GMT_get_rgb_from_z (value, Ctrl->G.fill.rgb);

			x = GMT_i_to_x (i, h[0].x_min, h[0].x_max, h[0].x_inc, h[0].xy_off, h[0].nx);
			if (!Ctrl->N.active) {
				GMT_map_outside (x, y);
				if (GMT_abs (GMT_x_status_new) > 1 || GMT_abs (GMT_y_status_new) > 1) continue;
			}
			GMT_geo_to_xy (x, y, &plot_x, &plot_y);

			if (Ctrl->T.active) {
				if (!Ctrl->Z.active) theta[ij] = (float)90.0 - theta[ij];
				GMT_azim_to_angle (x, y, 0.1, (double)theta[ij], &tmp);
				theta[ij] = (float)(tmp * D2R);
			}
			else
				theta[ij] *= (float)D2R;

			if (Ctrl->S.constant)
				r[ij] = (float)Ctrl->S.factor;
			else
				r[ij] *= (float)Ctrl->S.factor;

			if (project_info.projection == GMT_LINEAR) {	/* Must check if negative scales were used */
				if (!project_info.xyz_pos[0]) {	/* Negative x scale */
					if (!project_info.xyz_pos[1])	/* Negative y-scale too */
						theta[ij] += 180.0f;
					else
						theta[ij] = 180.0f - theta[ij];
				}
				else if (!project_info.xyz_pos[1])	/* Negative y-scale only */
					theta[ij] = - theta[ij];
			}
			GMT_flip_angle_f (&theta[ij]);
			sincos (theta[ij], &s, &c);
			x2 = plot_x + r[ij] * c;
			y2 = plot_y + r[ij] * s;

			if (Ctrl->E.active) {
				x_off = 0.5 * (x2 - plot_x);
				y_off = 0.5 * (y2 - plot_y);
				plot_x -= x_off;
				plot_y -= y_off;
				x2 -= x_off;
				y2 -= y_off;
			}

			if (!Ctrl->Q.active) {
				if (Ctrl->C.active) ps_setpaint (Ctrl->G.fill.rgb);
				ps_segment (plot_x, plot_y, x2, y2);
				continue;
			}

			if (shrink_properties && r[ij] < Ctrl->Q.norm) {	/* Scale arrow attributes down with length */
				v_w = Ctrl->Q.thickness * r[ij] * v_shrink;
				h_l = Ctrl->Q.length * r[ij] * v_shrink;
				h_w = Ctrl->Q.width * r[ij] * v_shrink;
				ps_vector (plot_x, plot_y, x2, y2, v_w, h_l, h_w, gmtdefs.vector_shape, Ctrl->G.fill.rgb, Ctrl->W.active);
			}
			else	/* Leave as specified */
				ps_vector (plot_x, plot_y, x2, y2, Ctrl->Q.thickness, Ctrl->Q.length, Ctrl->Q.width, gmtdefs.vector_shape, Ctrl->G.fill.rgb, Ctrl->W.active);
		}
	}

        if (!Ctrl->N.active) GMT_map_clip_off ();

	GMT_map_basemap ();

	GMT_plotend ();

	GMT_free ((void *)r);
	GMT_free ((void *)theta);

	Free_grdvector_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void *New_grdvector_Ctrl () {	/* Allocate and initialize a new control structure */
	struct GRDVECTOR_CTRL *C;
	
	C = (struct GRDVECTOR_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GRDVECTOR_CTRL), "New_grdvector_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	GMT_init_fill (&C->G.fill, -1, -1, -1);
	GMT_init_pen (&C->W.pen, GMT_PENWIDTH);
	C->Q.thickness = 0.03;
	C->Q.length = 0.12;
	C->Q.width = 0.1;
	C->S.factor = 1.0;
	return ((void *)C);
}

void Free_grdvector_Ctrl (struct GRDVECTOR_CTRL *C) {	/* Deallocate control structure */
	if (C->C.file) free ((void *)C->C.file);	
	GMT_free ((void *)C);	
}
