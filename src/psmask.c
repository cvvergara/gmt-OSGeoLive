/*--------------------------------------------------------------------
 *	$Id: psmask.c,v 1.96 2011/07/08 21:27:06 guru Exp $
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
 * psmask tries to achieve masking using one of two different approaches:
 * The default way of operation is: Instead of painting tiles where there's no
 * data [or where there is data] psmask uses contouring to find the polygons
 * that contain the data [or contain the regions with no data].  For many types
 * of data coverage, this results in a manageable path instead of thousands of
 * tiles.  So, instead of painting polygons, psmask sets up actual clip paths.
 * As an option, the user may specify a rgb combination to fill the clipped areas.
 * To avoid having to deal with the problems that arise if the data distribution
 * is such that part of the map boundary should be part of the clip paths, we
 * internally enlarge the grid by one gridsize unit so no nodes along the edges
 * have data.  Before using the clippaths we move points outside the real region
 * onto the boundary.
 * Should the paths be too long for PostScript to handle, the user may override
 * the default with the -T switch.  Then, masking is achieved using tiling
 *
 * Author:	Paul Wessel
 * Date:	01-JUL-2000
 * Version:	4
 *
 */

#include "gmt.h"
#include "pslib.h"

struct PSMASK_CTRL {
	struct C {	/* -C */
		GMT_LONG active;
	} C;
	struct D {	/* -D<dumpfile> */
		GMT_LONG active;
		char *file;
		GMT_LONG min_pts;
	} D;
	struct E {	/* -Eazim/elev */
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
	struct I {	/* -Idx[/dy] */
		GMT_LONG active;
		double xinc, yinc;
	} I;
	struct N {	/* -N */
		GMT_LONG active;
	} N;
	struct S {	/* -S<radius>[m|c|k|K] */
		GMT_LONG active;
		double radius;
		char unit;
	} S;
	struct T {	/* -T */
		GMT_LONG active;
	} T;
};

struct PSMASK_INFO {
	GMT_LONG first_dump;
	GMT_LONG p[5], i_off[5], j_off[5], k_off[5], offset;
	unsigned int bit[32];
};

int main (int argc, char **argv)
{
	GMT_LONG ij, n, nm, n_read, n_alloc;
	GMT_LONG i, j, n_edges, di, dj, ii, jj, n_expected_fields, *edge = NULL;
	GMT_LONG section, n_fields, distance_flag = 0, n_plus = 0, k;
	GMT_LONG error = FALSE, first = TRUE, node_only;

	char line[BUFSIZ], *grd = NULL, *not_used = NULL;

	double *in = NULL, distance, x0, y0, x1, y1, shrink = 1.0;
	double xinc2, yinc2, *x = NULL, *y = NULL;

	FILE *fp = NULL;

	struct GRD_HEADER h;
	struct PSMASK_INFO info;
	struct PSMASK_CTRL *Ctrl = NULL;

#ifdef DEBUG
	GMT_LONG debug = FALSE;
#endif
	void draw_clip_contours (double *xx, double *yy, GMT_LONG nn, int rgb[], GMT_LONG id, GMT_LONG flag);
	void dump_clip_contours (struct PSMASK_INFO *info, double *xx, double *yy, GMT_LONG nn, GMT_LONG id, char *file);
	void shrink_clip_contours (double *x, double *y, GMT_LONG n, double w, double e);
	GMT_LONG clip_contours (struct PSMASK_INFO *info, char *grd, struct GRD_HEADER *h, double xinc2, double yinc2, GMT_LONG *edge, GMT_LONG first, double **x, double **y, GMT_LONG *max);
	void *New_psmask_Ctrl (), Free_psmask_Ctrl (struct PSMASK_CTRL *C);
	
	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct PSMASK_CTRL *)New_psmask_Ctrl ();	/* Allocate and initialize a new control structure */

	memset ((void *)&info, 0, sizeof (struct PSMASK_INFO));
	info.first_dump = TRUE;

	GMT_grd_init (&h, argc, argv, FALSE);

	GMT_init_fill (&Ctrl->G.fill, -1, -1, -1);
	xinc2 = yinc2 = 0.0;

	Ctrl->D.min_pts = 0;

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
				case 'b':
				case 'c':
				case ':':
				case '\0':
					error += GMT_parse_common_options (argv[i], &h.x_min, &h.x_max, &h.y_min, &h.y_max);
					break;

				/* Supplemental parameters */

				case 'C':
					Ctrl->C.active = TRUE;
					break;
				case 'D':	/* Dump the polygons to files */
					Ctrl->D.active = TRUE;
					free ((void *)Ctrl->D.file);

					for (n_plus = 0, k = 2; argv[i][k]; k++) {
						if (argv[i][k] == '+' && argv[i][k+1] == 'n') {
							Ctrl->D.min_pts = atoi(&argv[i][k + 2]);
							n_plus = 1;
							break;
						}
					}

					if (n_plus)	/* If extra option rip it before check if there is a prefix */
						argv[i][k] = '\0';
					Ctrl->D.file = strdup (&argv[i][2]);
					break;
				case 'E':
					Ctrl->E.active = TRUE;
					error += GMT_get_proj3D (&argv[i][2], &Ctrl->E.azimuth, &Ctrl->E.elevation);
					break;
				case 'F':
					Ctrl->F.active = TRUE;
					break;
				case 'G':
					Ctrl->G.active = TRUE;
					if (GMT_getfill (&argv[i][2], &Ctrl->G.fill)) {
						GMT_fill_syntax ('G', " ");
						error++;
					}
					break;
				case 'I':
					Ctrl->I.active = TRUE;
					if (GMT_getinc (&argv[i][2], &Ctrl->I.xinc, &Ctrl->I.yinc)) {
						GMT_inc_syntax ('I', 1);
						error = TRUE;
					}
					break;
				case 'N':
					Ctrl->N.active = TRUE;
					break;
				case 'S':	/* Radius of influence */
					Ctrl->S.active = TRUE;
					Ctrl->S.radius = GMT_getradius (&argv[i][2]);
					Ctrl->S.unit = argv[i][strlen(argv[i])-1];
					break;
				case 'T':
					Ctrl->T.active = TRUE;
					break;
#ifdef DEBUG
				case 'd':
					debug = TRUE;
					break;
#endif
				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else {
			if ((fp = GMT_fopen(argv[i], GMT_io.r_mode)) == NULL) {
				fprintf (stderr, "%s: Could not open file %s\n", GMT_program, argv[i]);
				exit (EXIT_FAILURE);
			}
		}
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr,"psmask %s - Masking or clipping of 2-D data sets\n\n", GMT_VERSION);
		fprintf (stderr, "usage: psmask <xyz-file> %s %s %s\n", GMT_I_OPT, GMT_J_OPT, GMT_Rgeo_OPT);
		fprintf (stderr, "\t[%s] [-C] -D[-]<file>[+n<n_pts>] [%s] [-F] [-G<fill>] [%s] [-K]\n", GMT_B_OPT, GMT_E_OPT, GMT_H_OPT);
		fprintf (stderr, "\t[-N] [-O] [-P] [-S<radius>[k]] [-T] [%s] [-V] [%s]\n", GMT_U_OPT, GMT_X_OPT);
		fprintf (stderr, "\t[%s] [%s] [%s] [%s] [%s] [-n]\n\n", GMT_Y_OPT, GMT_c_OPT, GMT_t_OPT, GMT_b_OPT, GMT_mo_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\t<xyz-file> is the datafile.  If not given, read standard input\n");
		GMT_inc_syntax ('I', 0);
		GMT_explain_option ('j');
		GMT_explain_option ('R');
		fprintf (stderr, "\n\tOPTIONS:\n");
		GMT_explain_option ('b');
		fprintf (stderr, "\t-C means stop existing clip-path.  No other options required\n");
		fprintf (stderr, "\t-D dumps the clip-paths to files using the prefix <file>_ [mask_]\n");
		fprintf (stderr, "\t   Ignored if -T is specified.  If prefix begins with - we will instead\n");
		fprintf (stderr, "\t   use the rest of the name for a single multi-segment file.\n");
		fprintf (stderr, "\t   Append +n<n_pts> to limit the number of points in files to a minimum of n_pts.\n");
		GMT_explain_option ('E');
		fprintf (stderr, "\t-F Force pixel registration [Default is gridline registration]\n");
		GMT_fill_syntax ('G', "Select fill color/pattern [Default is no fill].");
		GMT_explain_option ('H');
		GMT_explain_option ('K');
		fprintf (stderr, "\t-N will invert the sense of the clipping [or tiling]\n");
		GMT_explain_option ('O');
		GMT_explain_option ('P');
		fprintf (stderr, "\t-S sets search radius in -R, -I units; append m or c for minutes or seconds.\n");
		fprintf (stderr, "\t   This means nodes inside circles of <radius> centered on\n");
		fprintf (stderr, "\t   the input data points are considered to be reliable estimates of the surface\n");
		fprintf (stderr, "\t   Default is -S0, i.e., only the nearest node is considered reliable\n");
		fprintf (stderr, "\t   Append k for km (implies -R,-I in degrees), use flat Earth approximation.\n");
		fprintf (stderr, "\t   Append K for km (implies -R,-I in degrees), use exact geodesic distances.\n");
		fprintf (stderr, "\t   If the current ELLIPSOID is spherical then great circle distances are used.\n");
		fprintf (stderr, "\t-T will paint tiles.  [Default will trace data outline]\n");
		fprintf (stderr, "\t   If set you must also specify a color/fill with -G\n");
		GMT_explain_option ('U');
		GMT_explain_option ('V');
		GMT_explain_option ('X');
		GMT_explain_option ('c');
		GMT_explain_option (':');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf (stderr, "\t   Default is 2 input columns\n");
		GMT_explain_option ('o');
		GMT_explain_option ('n');
		GMT_explain_option ('m');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	GMT_check_lattice (&Ctrl->I.xinc, &Ctrl->I.yinc, &Ctrl->F.active, &Ctrl->I.active);

	if (!Ctrl->C.active) {
		if (Ctrl->E.elevation <= 0.0 || Ctrl->E.elevation > 90.0) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR -E option:  Elevation must be in 0-90 range\n", GMT_program);
			error++;
		}
		if (!project_info.region_supplied) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify -R option\n", GMT_program);
			error++;
		}
		if (Ctrl->T.active && !GMT_IS_RECT_GRATICULE) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR -T option:  Only available with Linear, Mercator, or basic cylindrical projections\n", GMT_program);
			error++;
		}
		if (Ctrl->T.active && !(Ctrl->G.fill.rgb[0] >= 0 || Ctrl->G.fill.use_pattern)) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR -T option:  Must also specify a tile color with -G\n", GMT_program);
			error++;
		}
		if (!Ctrl->I.active) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify -I option\n", GMT_program);
			error++;
		}
		else if (Ctrl->I.xinc <= 0.0 || Ctrl->I.yinc <= 0.0) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR -I option:  Must specify positive increments\n", GMT_program);
			error++;
		}
		if (Ctrl->S.active && (Ctrl->S.radius <= 0.0 || GMT_is_dnan (Ctrl->S.radius))) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR.  Radius is NaN or negative\n", GMT_program);
			error++;
		}
		if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = 2;
		if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] < 2) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data (-bi) must have at least 2 columns\n", GMT_program);
			error++;
		}
		if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
			error++;
		}
		h.node_offset = (int)Ctrl->F.active;
		h.x_inc = Ctrl->I.xinc;
		h.y_inc = Ctrl->I.yinc;
		h.xy_off = 0.5 * h.node_offset;
		xinc2 = 0.5 * h.x_inc;
		yinc2 = 0.5 * h.y_inc;
		GMT_RI_prepare (&h);	/* Ensure -R -I consistency and set nx, ny */
	}

	if (error) exit (EXIT_FAILURE);

	/* Now test also if we must provide a default -J */
	if (Ctrl->D.active && project_info.projection < 0) {		/* Is this the right way of testing it? */
		GMT_parse_J_option ("x1d");
	}

	if (GMT_io.binary[GMT_IN] && gmtdefs.verbose) {
		char *type[2] = {"double", "single"};
		fprintf (stderr, "%s: Expects %ld-column %s-precision binary data\n", GMT_program, GMT_io.ncol[GMT_IN], type[GMT_io.single_precision[GMT_IN]]);
	}
	if (Ctrl->S.unit == 'k') distance_flag = 1;
	if (Ctrl->S.unit == 'K') distance_flag = 2;
	if (distance_flag) {
		shrink = cosd (0.5 * (h.y_min + h.y_max));
		di = (GMT_LONG)ceil (Ctrl->S.radius / (project_info.DIST_KM_PR_DEG * h.x_inc * shrink) + GMT_CONV_LIMIT);
		dj = (GMT_LONG)ceil (Ctrl->S.radius / (project_info.DIST_KM_PR_DEG * h.y_inc) + GMT_CONV_LIMIT);
	}
	else {
		di = irint (0.5 * Ctrl->S.radius / h.x_inc + GMT_CONV_LIMIT);
		dj = irint (0.5 * Ctrl->S.radius / h.y_inc + GMT_CONV_LIMIT);
	}
	z_project.view_azimuth = Ctrl->E.azimuth;
	z_project.view_elevation = Ctrl->E.elevation;
	if (!Ctrl->C.active) GMT_err_fail (GMT_map_setup (h.x_min, h.x_max, h.y_min, h.y_max), "");

	if (!project_info.x_off_supplied && GMT_ps.overlay) GMT_ps.x_origin = 0.0;
	if (!project_info.y_off_supplied && GMT_ps.overlay) GMT_ps.y_origin = 0.0;

	if (Ctrl->C.active)
		GMT_ps.clip = -1;	/* Signal that this program terminates clipping that initiated prior to this process */
	else if (!Ctrl->T.active)
		GMT_ps.clip = +1;	/* Signal that this program initiates clipping that wil outlive this process */
		
	GMT_plotinit (argc, argv);

	if (!Ctrl->C.active) {	/* Start new clip_path */

		if (project_info.three_D) ps_transrotate (-z_project.xmin, -z_project.ymin, 0.0);

		if (gmtdefs.verbose) fprintf (stderr, "%s: Allocate memory, read and process data file\n", GMT_program);

		/* Enlarge region by 1 row/column */

		h.x_min -= h.x_inc;	h.x_max += h.x_inc;	h.y_min -= h.y_inc;	h.y_max += h.y_inc;

		h.nx = GMT_get_n (h.x_min, h.x_max, h.x_inc, h.node_offset);
		h.ny = GMT_get_n (h.y_min, h.y_max, h.y_inc, h.node_offset);

		nm = h.nx * h.ny;
		grd = (char *) GMT_memory (VNULL, (size_t)nm, sizeof (char), GMT_program);

		/* Add GMT_CONV_LIMIT to ensure that special case radius = inc --> irint(0.5) actually rounds to 1 */
		
		if (distance_flag == 2 && !GMT_IS_SPHERICAL) distance_flag = 3;	/* Use geodesics */
		switch (distance_flag) {	/* Take different action depending on how we want distances calculated */

			case 0:		/* Cartesian distance */
				GMT_distance_func = GMT_cartesian_dist;
				break;
			case 1:		/* Flat Earth Approximation */
				GMT_distance_func = GMT_flatearth_dist_km;
				break;
			case 2:		/* Full spherical calculation */
				GMT_distance_func = GMT_great_circle_dist_km;
				break;
			case 3:		/* Full Ellipsoidal calculation */
				GMT_distance_func = GMT_geodesic_dist_km;
			break;
		}
		node_only = (di == 0 && dj == 0);
		if (node_only && Ctrl->S.radius > 0.0) {
			fprintf (stderr, "%s: Warning: Your search radius is too small to have any effect and is ignored.\n", GMT_program);
		}
		if (fp == NULL) {
			fp = GMT_stdin;
#ifdef SET_IO_MODE
			GMT_setmode (GMT_IN);
#endif
		}

		if (GMT_io.io_header[GMT_IN]) for (i = 0; i < GMT_io.n_header_recs; i++) not_used = GMT_fgets (line, BUFSIZ, fp);

		n_expected_fields = (GMT_io.ncol[GMT_IN]) ? GMT_io.ncol[GMT_IN] : 2;
		n_read = 0;
		while ((n_fields = GMT_input (fp, &n_expected_fields, &in)) >= 0 && !(GMT_io.status & GMT_IO_EOF)) {

			while ((GMT_io.status & GMT_IO_SEGMENT_HEADER) && !(GMT_io.status & GMT_IO_EOF)) { /* Skip headers */
				n_fields = GMT_input (fp, &n_expected_fields, &in);
			}
			if ((GMT_io.status & GMT_IO_EOF)) continue;	/* At EOF */

			n_read++;

			if (GMT_io.status & GMT_IO_MISMATCH) {
				fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld (skipped)\n", GMT_program, n_fields, n_expected_fields, n_read);
				continue;
			}

			if (GMT_y_is_outside (in[GMT_Y], h.y_min, h.y_max)) continue;		/* Outside y-range */
			if (GMT_x_is_outside (&in[GMT_X], h.x_min, h.x_max)) continue;	/* Outside x-range (or longitude) */

			/* Determine the node closest to the data point */

			i = GMT_x_to_i (in[GMT_X], h.x_min, h.x_inc, h.xy_off, h.nx);
			if (i < 0 || i >= h.nx) continue;
			j = GMT_y_to_j (in[GMT_Y], h.y_min, h.y_inc, h.xy_off, h.ny);
			if (j < 0 || j >= h.ny) continue;

			if (node_only) {
				grd[j*h.nx+i] = 1;
			}
			else {

				/* Set coordinate of this node */

				x0 = GMT_i_to_x (i, h.x_min, h.x_max, h.x_inc, h.xy_off, h.nx);
				y0 = GMT_j_to_y (j, h.y_min, h.y_max, h.y_inc, h.xy_off, h.ny);

				/* Set this and all nodes within radius distance to 1 */

				for (ii = i - di; ii <= i + di; ii++) {
					if (ii < 0 || ii >= h.nx) continue;
					x1 = GMT_i_to_x (ii, h.x_min, h.x_max, h.x_inc, h.xy_off, h.nx);
					for (jj = j - dj; jj <= j + dj; jj++) {
						if (jj < 0 || jj >= h.ny) continue;
						y1 = GMT_j_to_y (jj, h.y_min, h.y_max, h.y_inc, h.xy_off, h.ny);
						distance = (GMT_distance_func) (x1, y1, x0, y0);
						if (distance > Ctrl->S.radius) continue;
						grd[jj*h.nx+ii] = 1;
					}
				}
			}
		}

		if (fp != GMT_stdin) GMT_fclose (fp);
		if (gmtdefs.verbose) fprintf (stderr, "%s: Read %ld data points\n", GMT_program, n_read);

		if (Ctrl->N.active) for (i = 0; i < nm; i++) grd[i] = 1 - grd[i];	/* Reverse sense of test */

		/* Force perimeter nodes to be FALSE */

		for (i = 0, ij = (h.ny-1) * h.nx; i < h.nx; i++) grd[i] = grd[i+ij] = FALSE;
		for (j = 0; j < h.ny; j++) grd[j*h.nx] = grd[(j+1)*h.nx-1] = FALSE;

#ifdef DEBUG
		if (debug) {
			float *z;
			z = (float *) GMT_memory (VNULL, (size_t)nm, sizeof (float), GMT_program);
			for (i = 0; i < nm; i++) z[i] = (float)grd[i];
			GMT_write_grd ("psmask.grd", &h, z, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE);
			GMT_free ((void *)z);
		}
#endif
		if (!Ctrl->T.active) {	/* Must trace the outline of ON/OFF values in grd */
			/* Arrays holding the contour xy values */
			x = (double *) GMT_memory (VNULL, (size_t)GMT_CHUNK, sizeof (double), GMT_program);
			y = (double *) GMT_memory (VNULL, (size_t)GMT_CHUNK, sizeof (double), GMT_program);
			n_alloc = GMT_CHUNK;

			n_edges = h.ny * (GMT_LONG )ceil (h.nx / 16.0);
			edge = (GMT_LONG *) GMT_memory (VNULL, (size_t)n_edges, sizeof (GMT_LONG), GMT_program);

			GMT_map_basemap ();

			if (gmtdefs.verbose) fprintf (stderr, "%s: Tracing the clip path\n", GMT_program);

			section = 0;
			first = TRUE;
			while ((n = clip_contours (&info, grd, &h, xinc2, yinc2, edge, first, &x, &y, &n_alloc)) > 0) {
				shrink_clip_contours (x, y, n, h.x_min, h.x_max);
				if (Ctrl->D.active && n > Ctrl->D.min_pts) dump_clip_contours (&info, x, y, n, section, Ctrl->D.file);
				draw_clip_contours (x, y, n, Ctrl->G.fill.rgb, section, first);
				first = FALSE;
				section++;
			}

			draw_clip_contours (x, y, (GMT_LONG)0, Ctrl->G.fill.rgb, section, 2);	/* Activate clip-path */

			GMT_free ((void *)edge);
			GMT_free ((void *)x);
			GMT_free ((void *)y);
		}
		else {	/* Just paint tiles */
			GMT_LONG start, n_use, np, plot_n;
			double *grd_x, grd_y, y_bot, y_top, *xx, *yy, *xp, *yp;
			if (gmtdefs.verbose) fprintf (stderr, "%s: Tiling...\n", GMT_program);
			grd_x = (double *) GMT_memory (VNULL, (size_t)h.nx, sizeof (double), GMT_program);
			for (i = 0; i < h.nx; i++) grd_x[i] = GMT_i_to_x (i, h.x_min, h.x_max, h.x_inc, h.xy_off, h.nx);

			for (j = 0; j < h.ny; j++) {
				grd_y = GMT_j_to_y (j, h.y_min, h.y_max, h.y_inc, h.xy_off, h.ny);
				y_bot = grd_y - yinc2;
				y_top = grd_y + yinc2;
				ij = GMT_IJ (j, 1, h.nx);
				for (i = 0; i < h.nx; i++, ij++) {
					if (((GMT_LONG)grd[ij]) == 0) continue;

					np = GMT_graticule_path (&xx, &yy, 1, grd_x[i] - xinc2, grd_x[i] + xinc2, y_bot, y_top);
					plot_n = GMT_clip_to_map (xx, yy, np, &xp, &yp);
					GMT_free ((void *)xx);
					GMT_free ((void *)yy);
					if (plot_n == 0) continue;	/* Outside */
					
					if ((*GMT_will_it_wrap) (xp, yp, plot_n, &start)) {	/* Polygon wraps */

						/* First truncate against left border */

						GMT_n_plot = GMT_truncate (xp, yp, plot_n, start, -1);
						n_use = GMT_compact_line (GMT_x_plot, GMT_y_plot, GMT_n_plot, FALSE, 0);
						GMT_fill (GMT_x_plot, GMT_y_plot, n_use, &Ctrl->G.fill, FALSE);

						/* Then truncate against right border */

						GMT_n_plot = GMT_truncate (xp, yp, plot_n, start, +1);
						n_use = GMT_compact_line (GMT_x_plot, GMT_y_plot, GMT_n_plot, FALSE, 0);
						GMT_fill (GMT_x_plot, GMT_y_plot, n_use, &Ctrl->G.fill, FALSE);
					}
					else
						GMT_fill (xp, yp, plot_n, &Ctrl->G.fill, FALSE);
					GMT_free ((void *)xp);
					GMT_free ((void *)yp);
				}
			}

			GMT_free ((void *)grd_x);

			GMT_map_basemap ();
		}

		if (project_info.three_D) ps_rotatetrans (z_project.xmin, z_project.ymin, 0.0);

		GMT_free ((void *)grd);
		if (!Ctrl->T.active && gmtdefs.verbose) fprintf (stderr, "%s: clipping on!\n", GMT_program);
	}
	else {	/* Just undo previous clip-path */
		ps_clipoff ();
		GMT_map_basemap ();
		if (gmtdefs.verbose) fprintf (stderr, "%s: clipping off!\n", GMT_program);
	}

	ps_setpaint (gmtdefs.background_rgb);

	GMT_plotend ();

	Free_psmask_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void draw_clip_contours (double *xx, double *yy, GMT_LONG nn, int rgb[], GMT_LONG id, GMT_LONG flag)
{
	GMT_LONG i;
	double x, y;
	char comment[GMT_TEXT_LEN];

	if (nn < 2 && flag < 2) return;

	for (i = 0; i < nn; i++) {
		x = xx[i];	y = yy[i];
		GMT_geo_to_xy (x, y, &xx[i], &yy[i]);
	}
	nn = GMT_compact_line (xx, yy, nn, FALSE, 0);

	if (project_info.three_D) GMT_2D_to_3D (xx, yy, project_info.z_level, nn);

	if (nn > 0) {
		sprintf (comment, "Start of clip path sub-segment %ld", id);
		ps_comment (comment);
	}
	ps_clipon (xx, yy, nn, rgb, flag);
	if (nn > 0) {
		sprintf (comment, "End of clip path sub-segment %ld", id);
		ps_comment (comment);
	}
}

void dump_clip_contours (struct PSMASK_INFO *info, double *xx, double *yy, GMT_LONG nn, GMT_LONG id, char *file)
{
	GMT_LONG i;
	double out[2];
	char fname[BUFSIZ];
	FILE *fp;

	if (nn < 2) return;
	if (nn == 2 && GMT_IS_ZERO (xx[1] - xx[0]) && GMT_IS_ZERO (yy[1] - yy[0])) return;

	if (file[0] == '-') {	/* Want a single multi-segment file */
		sprintf (fname, "%s", &file[1]);
		if (info->first_dump && (fp = GMT_fopen (fname, GMT_io.w_mode)) == NULL) {
			fprintf (stderr, "%s: Unable to create file %s - exiting\n", GMT_program, fname);
			exit (EXIT_FAILURE);
		}
		else if ((fp = GMT_fopen (fname, GMT_io.a_mode)) == NULL) {
			fprintf (stderr, "%s: Unable to append to file %s - exiting\n", GMT_program, fname);
			exit (EXIT_FAILURE);
		}
		info->first_dump = FALSE;
		GMT_write_segmentheader (fp, 2);
	}
	else {
		if (GMT_io.binary[GMT_OUT])
			sprintf (fname, "%s_%ld.b", file, id);
		else
			sprintf (fname, "%s_%ld.xy", file, id);
		if ((fp = GMT_fopen (fname, GMT_io.w_mode)) == NULL) {
			fprintf (stderr, "%s: Unable to create file %s - exiting\n", GMT_program, fname);
			exit (EXIT_FAILURE);
		}
	}
	for (i = 0; i < nn; i++) {
		out[GMT_X] = xx[i];	out[GMT_Y] = yy[i];
		GMT_output (fp, 2, out);
	}
	GMT_fclose (fp);
}

GMT_LONG clip_contours (struct PSMASK_INFO *info, char *grd, struct GRD_HEADER *h, double xinc2, double yinc2, GMT_LONG *edge, GMT_LONG first, double **x, double **y, GMT_LONG *max)
{
	/* The routine finds the zero-contour in the grd dataset.  it assumes that
	 * no node has a value exactly == 0.0.  If more than max points are found
	 * trace_clip_contours will try to allocate more memory in blocks of GMT_CHUNK points
	 */
	 
	static GMT_LONG i0, j0, side;
	GMT_LONG ij;
	GMT_LONG i, j, n = 0, n_edges, edge_word, edge_bit;
	GMT_LONG go_on = TRUE;
	GMT_LONG trace_clip_contours (struct PSMASK_INFO *info, char *grd, GMT_LONG *edge, struct GRD_HEADER *h, double xinc2, double yinc2, double **xx, double **yy, GMT_LONG i, GMT_LONG j, GMT_LONG kk, GMT_LONG *max);
	 
	 
	n_edges = h->ny * (GMT_LONG) ceil (h->nx / 16.0);
	 
	 /* Reset edge-flags to zero, if necessary */
	 if (first) {
		info->offset = n_edges / 2;
	 	i0 = 0;	/* Begin with upper left bin which is i = 0 and j = 1 */
	 	j0 = 1;
		side = 4;	/* Vertical interior gridlines */
		info->p[0] = info->p[4] = 0;	info->p[1] = 1;	info->p[2] = 1 - h->nx;	info->p[3] = -h->nx;
		info->i_off[0] = info->i_off[2] = info->i_off[3] = info->i_off[4] = 0;	info->i_off[1] =  1;
		info->j_off[0] = info->j_off[1] = info->j_off[3] = info->j_off[4] = 0;	info->j_off[2] = -1;
		info->k_off[0] = info->k_off[2] = info->k_off[4] = 0;	info->k_off[1] = info->k_off[3] = 1;
		for (i = 1, info->bit[0] = 1; i < 32; i++) info->bit[i] = info->bit[i-1] << 1;
	 }

	/* Loop over interior boxes */

	if (side == 4) {
		for (j = j0; go_on && j < h->ny; j++) {
			ij = GMT_IJ (j, i0, h->nx);
			for (i = i0; go_on && i < h->nx-1; i++, ij++) {	/* nx-1 since the last bin starts at nx-2 and ends at nx-1 */
				edge_word = (GMT_LONG)(ij / 32 + info->offset);
				edge_bit = (GMT_LONG)(ij % 32);
				if (!(edge[edge_word] & info->bit[edge_bit]) && ((grd[ij]+grd[ij-h->nx]) == 1)) { /* Start tracing contour */
					*x[0] = GMT_i_to_x (i, h->x_min, h->x_max, h->x_inc, h->xy_off, h->nx);
					*y[0] = GMT_j_to_y (j, h->y_min, h->y_max, h->y_inc, h->xy_off, h->ny);
					edge[edge_word] |= info->bit[edge_bit];
					n = trace_clip_contours (info, grd, edge, h, xinc2, yinc2, x, y, i, j, 3, max);
					go_on = FALSE;
					i0 = i + 1;
					j0 = j;	/* Return to finish this row later */
				}
			}
			if (go_on) i0 = 0;	/* Go to start of next row unless we found something */
		}
		if (n == 0) {
			side = 5;
			i0 = 0;
			j0 = 1;
		}
	}
	if (n == 0 && side == 5) {
		for (j = j0; go_on && j < h->ny; j++) {
			ij = GMT_IJ (j, i0, h->nx);
			for (i = i0; go_on && i < h->nx-1; i++, ij++) {
				edge_word = (GMT_LONG)(ij / 32 + info->offset);
				edge_bit = (GMT_LONG)(ij % 32);
				if (!(edge[edge_word] & info->bit[edge_bit]) && ((grd[ij]+grd[ij+1]) == 1)) { /* Start tracing contour */
					*x[0] = GMT_i_to_x (i, h->x_min, h->x_max, h->x_inc, h->xy_off, h->nx);
					*y[0] = GMT_j_to_y (j, h->y_min, h->y_max, h->y_inc, h->xy_off, h->ny);
					edge[edge_word] |= info->bit[edge_bit];
					n = trace_clip_contours (info, grd, edge, h, xinc2, yinc2, x, y, i, j, 2, max);
					go_on = FALSE;
					i0 = i + 1;
					j0 = j;	/* Return to finish this row later */
				}
				if (go_on) i0 = 1;
			}
		}
	}	

	return (n);
}

GMT_LONG trace_clip_contours (struct PSMASK_INFO *info, char *grd, GMT_LONG *edge, struct GRD_HEADER *h, double xinc2, double yinc2, double **xx, double **yy, GMT_LONG i, GMT_LONG j, GMT_LONG kk, GMT_LONG *max)
{
	GMT_LONG n = 1, k, k0, n_cuts, kk_opposite, first_k, more;
	GMT_LONG edge_word, edge_bit;
	GMT_LONG ij, ij0, m;
	double xk[4], yk[4], x0, y0;

	m = *max - 2;
	
	more = TRUE;
	do {
		ij = GMT_IJ (j, i, h->nx);
		x0 = GMT_i_to_x (i, h->x_min, h->x_max, h->x_inc, h->xy_off, h->nx);
		y0 = GMT_j_to_y (j, h->y_min, h->y_max, h->y_inc, h->xy_off, h->ny);
		n_cuts = 0;
		k0 = kk;

		for (k = 0; k < 4; k++) {	/* Loop over box sides */

			/* Skip where we already have a cut (k == k0) */

			if (k == k0) continue;

			/* Skip edge already has been used */

			ij0 = GMT_IJ (j + info->j_off[k], i + info->i_off[k], h->nx);
			edge_word = (GMT_LONG)(ij0 / 32 + info->k_off[k] * info->offset);
			edge_bit = (GMT_LONG)(ij0 % 32);
			if (edge[edge_word] & info->bit[edge_bit]) continue;

			/* Skip if no zero-crossing on this edge */

			if ((grd[ij+info->p[k+1]] + grd[ij+info->p[k]]) != 1) continue;

			/* Here we have a crossing */

			if (k%2) {	/* Cutting a S-N line */
				if (k == 1) {
					xk[1] = x0 + h->x_inc;
					yk[1] = y0 + 0.5*h->y_inc;
				}
				else {
					xk[3] = x0;
					yk[3] = y0 + 0.5*h->y_inc;
				}
			}
			else {	/* Cutting a E-W line */
				if (k == 0) {
					xk[0] = x0 + 0.5*h->x_inc;
					yk[0] = y0;
				}
				else {
					xk[2] = x0 + 0.5*h->x_inc;
					yk[2] = y0 + h->y_inc;
				}
			}
			kk = k;
			n_cuts++;
		}

		if (n > m) {	/* Must try to allocate more memory */
			*max = (*max == 0) ? GMT_CHUNK : ((*max) << 1);
			m = (m == 0) ? GMT_CHUNK : (m << 1);
			*xx = (double *) GMT_memory ((void *)*xx, (size_t)(*max), sizeof (double), "trace_clip_contours");
			*yy = (double *) GMT_memory ((void *)*yy, (size_t)(*max), sizeof (double), "trace_clip_contours");
		}
		if (n_cuts == 0) {	/* Close interior contour and return */
			/* if (fmod ((*xx[0] - xinc2), h->x_inc) == 0.0) */	/* On side 1 or 3 */
			if (GMT_IS_ZERO (fmod ((*xx[0] - xinc2), h->x_inc)))	/* On side 1 or 3 */
				/* first_k = ((*xx[0] - x0) == 0.0) ? 3 : 1; */
				first_k = GMT_IS_ZERO (*xx[0] - x0) ? 3 : 1;
			else 	/* On side 0 or 2 */
				/* first_k = ((*yy[0] - y0) == 0.0) ? 0 : 2; */
				first_k = GMT_IS_ZERO (*yy[0] - y0) ? 0 : 2;
			kk_opposite = (first_k + 2) % 4;
			if (k0 != kk_opposite) {
				(*xx)[n] = x0 + 0.5*h->x_inc;
				(*yy)[n] = y0 + 0.5*h->y_inc;
				n++;
			}
			(*xx)[n] = (*xx)[0];
			(*yy)[n] = (*yy)[0];
			n++;
			more = FALSE;
		}
		else if (n_cuts == 1) {	/* Draw a line to this point and keep tracing */
			/* Add center of box if this and previous cut NOT on opposite edges */
			kk_opposite = (k0 + 2) % 4;
			if (kk != kk_opposite) {
				(*xx)[n] = x0 + 0.5*h->x_inc;
				(*yy)[n] = y0 + 0.5*h->y_inc;
				n++;
			}
			(*xx)[n] = xk[kk];
			(*yy)[n] = yk[kk];
			n++;
		}
		else {	/* Saddle point, we decide to connect to the point nearest previous point */
			kk = (k0 + 1)%4;	/* Pick next edge since it is arbitrarily where we go */
			/* First add center of box */
			(*xx)[n] = x0 + 0.5*h->x_inc;
			(*yy)[n] = y0 + 0.5*h->y_inc;
			n++;
			(*xx)[n] = xk[kk];
			(*yy)[n] = yk[kk];
			n++;
		}
		if (more) {	/* Mark this edge as used */
			ij0 = GMT_IJ (j + info->j_off[kk], i + info->i_off[kk], h->nx);
			edge_word = (GMT_LONG)(ij0 / 32 + info->k_off[kk] * info->offset);
			edge_bit = (GMT_LONG)(ij0 % 32);
			edge[edge_word] |= info->bit[edge_bit];
		}

		/* Get next box (i,j,kk) */

		i -= (kk-2)%2;
		j -= (kk-1)%2;
		kk = (kk+2)%4;

	} while (more);
	return (n);
}

void shrink_clip_contours (double *x, double *y, GMT_LONG n, double w, double e)
{
	/* Moves outside points to boundary */
	GMT_LONG i;

	for (i = 0; i < n; i++) {
		if (x[i] < w)
			x[i] = w;
		if (x[i] > e)
			x[i] = e;
		if (y[i] < project_info.s)
			y[i] = project_info.s;
		if (y[i] > project_info.n)
			y[i] = project_info.n;
	}
}

void *New_psmask_Ctrl () {	/* Allocate and initialize a new control structure */
	struct PSMASK_CTRL *C;
	
	C = (struct PSMASK_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct PSMASK_CTRL), "New_psmask_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
		
	C->D.file = strdup ("mask");
	C->E.azimuth = 180.0;
	C->E.elevation = 90.0;
	return ((void *)C);
}

void Free_psmask_Ctrl (struct PSMASK_CTRL *C) {	/* Deallocate control structure */
	if (C->D.file) free ((void *)C->D.file);	
	GMT_free ((void *)C);	
}
