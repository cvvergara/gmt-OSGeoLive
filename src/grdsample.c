/*--------------------------------------------------------------------
 *	$Id: grdsample.c,v 1.70 2011/07/08 21:27:06 guru Exp $
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
 * grdsample reads a grid file and evaluates the grid at new grid positions
 * specified by new dx/dy values using a 2-D Taylor expansion of order 3.
 * In order to evaluate derivatives along the edges of the surface, I assume 
 * natural bicubic spline conditions, i.e. both the second and third normal 
 * derivatives are zero, and that the dxdy derivative in the corners are zero, too.
 *
 * Author:	Paul Wessel
 * Date:	19-JUL-1989
 * Revised:	6-JAN-1990	PW: Updated to v.2.0
 * Revised:	16-JUN-1998	PW: Updated to v.3.1
 * Version:	4
 */

#include "gmt.h"

struct GRDSAMPLE_CTRL {
	struct F {	/* -F */
		GMT_LONG active;
	} F;
	struct G {	/* -G<grdfile> */
		GMT_LONG active;
		char *file;
	} G;
	struct I {	/* -Idx[/dy] */
		GMT_LONG active;
		double xinc, yinc;
	} I;
	struct L {	/* -L<flag> */
		GMT_LONG active;
		char mode[4];
	} L;
	struct Q {	/* -Q[b|c|l|n][[/]<threshold>] */
		GMT_LONG active;
		GMT_LONG interpolant;
		double threshold;
	} Q;
	struct T {	/* -T */
		GMT_LONG active;
	} T;
};

int main (int argc, char **argv)
{
	GMT_LONG error = FALSE;

	char *infile = CNULL, format[BUFSIZ];

	GMT_LONG ij, nm, i, j, ii = 0, jj = 0;
	
	float *a = NULL, *b = NULL;

	double *lon, lat;

	struct GRD_HEADER grd_a, grd_b;
	struct GMT_EDGEINFO edgeinfo;
	struct GMT_BCR bcr;
	struct GRDSAMPLE_CTRL *Ctrl = NULL;

	void *New_grdsample_Ctrl (), Free_grdsample_Ctrl (struct GRDSAMPLE_CTRL *C);

	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct GRDSAMPLE_CTRL *)New_grdsample_Ctrl ();	/* Allocate and initialize a new control structure */
	
	GMT_grd_init (&grd_b, argc, argv, FALSE);

	GMT_boundcond_init (&edgeinfo);

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {

				/* Common parameters */

				case 'R':
				case 'V':
				case 'f':
				case '\0':
					error += (int)GMT_parse_common_options (argv[i], &grd_b.x_min, &grd_b.x_max, &grd_b.y_min, &grd_b.y_max);
					break;

				/* Supplemental parameters */

				case 'F':
					Ctrl->F.active = TRUE;
					break;
				case 'G':
					Ctrl->G.active = TRUE;
					Ctrl->G.file = strdup (&argv[i][2]);
					break;
				case 'N':	/* Backwards compatible.  nx/ny can now be set with -I */
					Ctrl->I.active = TRUE;
					sscanf (&argv[i][2], "%" GMT_LL "d/%" GMT_LL "d", &ii, &jj);
					if (jj == 0) jj = ii;
					sprintf (format, "%" GMT_LL "d+/%" GMT_LL "d+", ii, jj);
					GMT_getinc (format, &Ctrl->I.xinc, &Ctrl->I.yinc);
					break;
				case 'I':
					Ctrl->I.active = TRUE;
					if (GMT_getinc (&argv[i][2], &Ctrl->I.xinc, &Ctrl->I.yinc)) {
						GMT_inc_syntax ('I', 1);
						error = TRUE;
					}
					break;
				case 'L':
					Ctrl->L.active = TRUE;
					strncpy (Ctrl->L.mode, &argv[i][2], (size_t)4);
					break;
				case 'Q':
					Ctrl->Q.active = TRUE;
					Ctrl->Q.interpolant = BCR_BILINEAR;
					for (j = 2; j < 5 && argv[i][j]; j++) {
						switch (argv[i][j]) {
							case 'n':
								Ctrl->Q.interpolant = BCR_NEARNEIGHBOR; break;
							case 'l':
								Ctrl->Q.interpolant = BCR_BILINEAR; break;
							case 'b':
								Ctrl->Q.interpolant = BCR_BSPLINE; break;
							case 'c':
								Ctrl->Q.interpolant = BCR_BICUBIC; break;
							case '/':
							default:
								Ctrl->Q.threshold = atof (&argv[i][j]);
								if (j == 2 && Ctrl->Q.threshold < GMT_SMALL) {
									Ctrl->Q.interpolant = BCR_NEARNEIGHBOR;
									fprintf (stderr, "%s: Warning: Option -Q0 deprecated. Use -Qn instead.\n", GMT_program);
								}
								j = 5; break;
						}
					}
					break;
				case 'T':	/* Convert from pixel file <-> gridfile */
					Ctrl->T.active = TRUE;
					break;
				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else 
			infile = argv[i];
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "grdsample %s - Resample a grid file onto a new grid\n\n", GMT_VERSION);
		fprintf (stderr, "usage: grdsample <old_grdfile> -G<new_grdfile> [-F] [%s] [-L<flag>]\n", GMT_I_OPT);
		fprintf (stderr, "\t[-Q[<value>]] [%s] [-T] [-V] [%s]\n", GMT_Rgeo_OPT, GMT_f_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\t<old_grdfile> is data set to be resampled\n");
		fprintf (stderr, "\t-G sets the name of the interpolated output grid file\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-F Force pixel registration  [Default is same as input]\n");
		GMT_inc_syntax ('I', 0);
		fprintf (stderr, "\t   When omitted: grid spacing is copied from input grid.\n");
		fprintf (stderr, "\t-L sets boundary conditions.  <flag> can be either\n");
		fprintf (stderr, "\t   g for geographic boundary conditions\n");
		fprintf (stderr, "\t   or one or both of\n");
		fprintf (stderr, "\t   x for periodic boundary conditions on x\n");
		fprintf (stderr, "\t   y for periodic boundary conditions on y\n");
		fprintf (stderr, "\t-Q Quick mode, use bilinear rather than bicubic [Default] interpolation.\n");
		fprintf (stderr, "\t   Alternatively, select interpolation mode by adding b = B-spline, c = bicubic,\n");
		fprintf (stderr, "\t   l = bilinear, or n = nearest-neighbor.\n");
		fprintf (stderr, "\t   Optionally, append <threshold> in the range [0,1]. [Default = 1 requires all\n");
		fprintf (stderr, "\t   4 or 16 nodes to be non-NaN.], <threshold> = 0.5 will interpolate about 1/2 way\n");
		fprintf (stderr, "\t   from a non-NaN to a NaN node, while 0.1 will go about 90%% of the way, etc.\n");
		fprintf (stderr, "\t   -Q0 will return the value of the nearest node instead of interpolating (Same as -Qn).\n");
		fprintf (stderr, "\t-R specifies a subregion [Default is old region]\n");
		fprintf (stderr, "\t-T Toggles between grid registration and pixel registration\n");
		GMT_explain_option ('V');
		GMT_explain_option ('f');

		exit (EXIT_FAILURE);
	}

	GMT_check_lattice (&Ctrl->I.xinc, &Ctrl->I.yinc, &Ctrl->F.active, &Ctrl->I.active);
	
	if (Ctrl->Q.active && (Ctrl->Q.threshold < 0.0 || Ctrl->Q.threshold > 1.0)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -Q:  threshold must be in [0,1] range\n", GMT_program);
		error++;
	}
	if (!infile) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify input file\n", GMT_program);
		error++;
	}
	if (!Ctrl->G.file) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -G:  Must specify output file\n", GMT_program);
		error++;
	}
	if (Ctrl->F.active && Ctrl->T.active) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Only one of -F, -T may be specified\n", GMT_program);
		error++;
	}
	if (Ctrl->I.active && (Ctrl->I.xinc <= 0.0 || Ctrl->I.yinc <= 0.0)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -I:  Must specify positive increments\n", GMT_program);
		error++;
	}
	if (Ctrl->L.active && GMT_boundcond_parse (&edgeinfo, Ctrl->L.mode)) error++;
	if (error) exit (EXIT_FAILURE);
	
	GMT_err_fail (GMT_read_grd_info (infile, &grd_a), infile);

	if (!project_info.region_supplied) {
		grd_b.x_min = grd_a.x_min;
		grd_b.x_max = grd_a.x_max;
		grd_b.y_min = grd_a.y_min;
		grd_b.y_max = grd_a.y_max;
	}

	if (Ctrl->I.active) {
		grd_b.x_inc = Ctrl->I.xinc;
		grd_b.y_inc = Ctrl->I.yinc;
	}
	else {
		grd_b.x_inc = grd_a.x_inc;
		grd_b.y_inc = grd_a.y_inc;
	}

	if (Ctrl->T.active)
		grd_b.node_offset = !grd_a.node_offset;
	else if (Ctrl->F.active)
		grd_b.node_offset = TRUE;
	else
		grd_b.node_offset = grd_a.node_offset;

	grd_b.xy_off = 0.5 * grd_b.node_offset;

	GMT_RI_prepare (&grd_b);	/* Ensure -R -I consistency and set nx, ny */

	GMT_boundcond_param_prep (&grd_a, &edgeinfo);

	nm = GMT_get_nm (4 + grd_a.nx, 4 + grd_a.ny);	/* Reserve space for padding */
	a = (float *) GMT_memory (VNULL, (size_t)nm, sizeof(float), GMT_program);

	if (project_info.region_supplied) {
		if (GMT_io.in_col_type[0] == GMT_IS_LON) {
			grd_b.x_min -= 720.0;
			grd_b.x_max -= 720.0;
			while (grd_b.x_max < grd_a.x_min) grd_b.x_min += 360.0, grd_b.x_max += 360.0;
		}
		if (!edgeinfo.nxp && ((grd_b.x_min + GMT_CONV_LIMIT) < grd_a.x_min || (grd_b.x_max - GMT_CONV_LIMIT) > grd_a.x_max)) {
			fprintf (stderr, "%s:  Selected region exceeds the X-boundaries of the grid file!\n", GMT_program);
			exit (EXIT_FAILURE);
		}
		else if (!edgeinfo.nyp && ((grd_b.y_min + GMT_CONV_LIMIT) < grd_a.y_min || (grd_b.y_max - GMT_CONV_LIMIT) > grd_a.y_max)) {
			fprintf (stderr, "%s:  Selected region exceeds the Y-boundaries of the grid file!\n", GMT_program);
			exit (EXIT_FAILURE);
		}
	}

	if (!Ctrl->I.active) {
		grd_b.x_inc = GMT_get_inc (grd_b.x_min, grd_b.x_max, grd_b.nx, grd_b.node_offset);
		grd_b.y_inc = GMT_get_inc (grd_b.y_min, grd_b.y_max, grd_b.ny, grd_b.node_offset);
	}

	GMT_err_fail (GMT_grd_RI_verify (&grd_b, 1), Ctrl->G.file);

	nm = GMT_get_nm (grd_b.nx, grd_b.ny);
	b = (float *) GMT_memory (VNULL, (size_t)nm, sizeof(float), GMT_program);

	sprintf (format, "%%s: New grid (%s/%s/%s/%s) nx = %%d ny = %%d dx = %s dy = %s node_offset = %%d\n", gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format);
	if (gmtdefs.verbose) fprintf (stderr, format, GMT_program, grd_b.x_min, grd_b.x_max, grd_b.y_min, grd_b.y_max, grd_b.nx, grd_b.ny, grd_b.x_inc, grd_b.y_inc, grd_b.node_offset);

	GMT_pad[0] = GMT_pad[1] = GMT_pad[2] = GMT_pad[3] = 2;	/* Leave room for 2 empty boundary rows/cols */

	GMT_err_fail (GMT_read_grd (infile, &grd_a, a, grd_a.x_min, grd_a.x_max, grd_a.y_min, grd_a.y_max, GMT_pad, FALSE), infile);

	/* Initialize bcr structure:  */

	GMT_bcr_init (&grd_a, GMT_pad, Ctrl->Q.interpolant, Ctrl->Q.threshold, &bcr);

	/* Set boundary conditions  */

	GMT_boundcond_set (&grd_a, &edgeinfo, GMT_pad, a);

	/* Precalculate longitudes */

	lon = (double *) GMT_memory (VNULL, (size_t)grd_b.nx, sizeof (double), GMT_program);
	for (i = 0; i < grd_b.nx; i++) {
		lon[i] = GMT_i_to_x (i, grd_b.x_min, grd_b.x_max, grd_b.x_inc, grd_b.xy_off, grd_b.nx);
		if (!edgeinfo.nxp)
			/* Nothing */;
		else if (lon[i] > grd_a.x_max)
			lon[i] -= grd_a.x_inc * edgeinfo.nxp;
		else if (lon[i] < grd_a.x_min)
			lon[i] += grd_a.x_inc * edgeinfo.nxp;
	}

	for (j = ij = 0; j < grd_b.ny; j++) {
		lat = GMT_j_to_y (j, grd_b.y_min, grd_b.y_max, grd_b.y_inc, grd_b.xy_off, grd_b.ny);
		if (!edgeinfo.nyp)
			/* Nothing */;
		else if (lat > grd_a.y_max)
			lat -= grd_a.y_inc * edgeinfo.nyp;
		else if (lat < grd_a.y_min)
			lat += grd_a.y_inc * edgeinfo.nyp;
		for (i = 0; i < grd_b.nx; i++, ij++) {
			b[ij] = (float)GMT_get_bcr_z (&grd_a, lon[i], lat, a, &edgeinfo, &bcr);
		}
	}

	GMT_pad[0] = GMT_pad[1] = GMT_pad[2] = GMT_pad[3] = 0;	/* No boundary rows/cols on output */
	GMT_err_fail (GMT_write_grd (Ctrl->G.file, &grd_b, b, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), Ctrl->G.file);

	GMT_free ((void *)a);
	GMT_free ((void *)b);
	GMT_free ((void *)lon);

	Free_grdsample_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void *New_grdsample_Ctrl () {	/* Allocate and initialize a new control structure */
	struct GRDSAMPLE_CTRL *C;
	
	C = (struct GRDSAMPLE_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GRDSAMPLE_CTRL), "New_grdsample_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	C->Q.interpolant = BCR_BICUBIC; C->Q.threshold = 1.0;
	return ((void *)C);
}

void Free_grdsample_Ctrl (struct GRDSAMPLE_CTRL *C) {	/* Deallocate control structure */
	if (C->G.file) free ((void *)C->G.file);	
	GMT_free ((void *)C);	
}
