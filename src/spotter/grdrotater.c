/*--------------------------------------------------------------------
 *	$Id: grdrotater.c,v 1.45 2011/07/11 19:22:06 guru Exp $
 *
 *   Copyright (c) 1999-2011 by P. Wessel
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; version 2 or any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   Contact info: www.soest.hawaii.edu/pwessel
 *--------------------------------------------------------------------*/
/*
 * grdrotater will read a grid file, apply a finite rotation to the grid
 * coordinates, and then interpolate the old grid at the new coordinates.
 *
 * Author:	Paul Wessel
 * Date:	27-JAN-2006
 * Ver:		1
 */

#include "spotter.h"

#define PAD 3

int main (int argc, char **argv)
{
	GMT_LONG i, j, k, ij, ii, jj, nx, ny, mx = 0, my, nm, ix, iy;
	GMT_LONG iii, jjj, ii_o, jj_o;
	GMT_LONG error = FALSE, interpolant = BCR_BICUBIC, user_polygon_given = FALSE;
	GMT_LONG inside, node = FALSE, outline = TRUE, skip_grid = FALSE, global = FALSE;
	
	float *original_grid = VNULL, *rotated_grid = NULL;

	char txt_a[GMT_TEXT_LEN], txt_b[GMT_TEXT_LEN], txt_c[GMT_TEXT_LEN];
	char *grdfile = NULL, *outfile = NULL, *pol_file = NULL;

	double  xx, yy, lon, x_rot, y_rot, w_rot = 0.0, P_original[3], P_rotated[3], R[3][3];
	double west, east, south, north, data_west, data_east, data_south, data_north;
	double w180, w360, e180, e360, threshold = 1.0, out[2], *grd_x = NULL, *grd_y = NULL, *grd_yc = NULL;
	double wp180, wp360, ep180, ep360, spol, npol;

	struct GRD_HEADER g_head, r_head;
	struct GMT_TABLE *pol = NULL;
	struct GMT_EDGEINFO edgeinfo;
	struct GMT_BCR bcr;

	void get_grid_path (struct GRD_HEADER *h, struct GMT_TABLE **p);

	argc = (int)GMT_begin (argc, argv);

	grdfile = outfile = pol_file = CNULL;
	west = east = south = north = 0.0;
	GMT_io.in_col_type[GMT_X] = GMT_IS_LON;
	GMT_io.in_col_type[GMT_Y] = GMT_IS_LAT;

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				/* Common parameters */

				case 'H':
				case 'M':
				case 'R':
				case 'V':
				case ':':
				case 'b':
				case 'm':
				case '\0':
					error += GMT_parse_common_options (argv[i], &west, &east, &south, &north);
					break;

				/* Supplemental parameters */

				case 'F':
					user_polygon_given = TRUE;
					pol_file = &argv[i][2];
					break;
				case 'T':	/* Finite rotation parameters */
					sscanf (&argv[i][2], "%[^/]/%[^/]/%s", txt_a, txt_b, txt_c);
					error += GMT_verify_expectations (GMT_io.in_col_type[GMT_X], GMT_scanf (txt_a, GMT_io.in_col_type[GMT_X], &x_rot), txt_a);
					error += GMT_verify_expectations (GMT_io.in_col_type[GMT_Y], GMT_scanf (txt_b, GMT_io.in_col_type[GMT_Y], &y_rot), txt_b);
					w_rot = atof (txt_c);
					break;
				case 'G':
					outfile = &argv[i][2];
					break;
				case 'N':
					outline = FALSE;
					break;
				case 'Q':
					interpolant = BCR_BILINEAR;
					for (j = 2; j < 5 && argv[i][j]; j++) {
						switch (argv[i][j]) {
							case 'n':
								interpolant = BCR_NEARNEIGHBOR; break;
							case 'l':
								interpolant = BCR_BILINEAR; break;
							case 'b':
								interpolant = BCR_BSPLINE; break;
							case 'c':
								interpolant = BCR_BICUBIC; break;
							case '/':
							default:
								threshold = atof (&argv[i][j]);
								if (j == 2 && threshold < GMT_SMALL) interpolant = BCR_NEARNEIGHBOR;
								j = 5; break;
						}
					}
					break;
				case 'S':
					skip_grid = TRUE;
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
		fprintf (stderr,"grdrotater %s - Finite rotation reconstruction of geographic grid\n\n", GMT_VERSION);
		fprintf (stderr, "usage: grdrotater <grdfile> -T<plon>/<plat>/<pomega> -G<rotgrid> [-F<polygonfile>]\n");
		fprintf (stderr, "\t[%s] [-N] [-Q[b|c|l|n][[/]<threshold>]] [%s] [-S] [-V] [%s]\n", GMT_H_OPT, GMT_Rgeo_OPT, GMT_t_OPT);
		fprintf (stderr, "\t[%s] [%s] > projpol\n\n", GMT_b_OPT, GMT_m_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr,"\t<grdfile> is a gridded data file in geographic coordinates to be rotated\n");
		fprintf (stderr,"\t-G is the output filename of the new, rotated grid.  The boundary of the\n");
		fprintf (stderr,"\t   original grid (or a subset; see -F) after rotation is written to stdout\n");
		fprintf (stderr,"\t   unless the grid is global.\n");
		fprintf (stderr,"\t-T sets the rotation pole and opening angle for the finite rotation (all in degrees)\n");
		fprintf (stderr,"\n\tOPTIONS:\n");
		fprintf (stderr,"\t-F specifies a multi-segment closed polygon file that describes the area of the grid\n");
		fprintf (stderr,"\t   that should be projected [Default projects entire grid]\n");
		GMT_explain_option ('H');
		fprintf (stderr, "\t-N Do not output the rotated polygon or grid outline\n");
		fprintf (stderr, "\t-Q Quick mode, use bilinear rather than bicubic [Default] interpolation.\n");
		fprintf (stderr, "\t   Alternatively, select interpolation mode by adding b = B-spline, c = bicubic,\n");
		fprintf (stderr, "\t   l = bilinear, or n = nearest-neighbor.\n");
		fprintf (stderr, "\t   Optionally, append <threshold> in the range [0,1]. [Default = 1 requires all\n");
		fprintf (stderr, "\t   4 or 16 nodes to be non-NaN.], <threshold> = 0.5 will interpolate about 1/2 way\n");
		fprintf (stderr, "\t   from a non-NaN to a NaN node, while 0.1 will go about 90%% of the way, etc.\n");
		fprintf (stderr, "\t   -Q0 will return the value of the nearest node instead of interpolating (Same as -Qn).\n");
		GMT_explain_option ('R');
		fprintf (stderr, "\t-S Do NOT rotate the grid - just produce the rotated outline\n");
		GMT_explain_option ('V');
		GMT_explain_option (':');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf (stderr, "\t   Default is 2 input columns\n");
		GMT_explain_option ('o');
		GMT_explain_option ('m');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	if (skip_grid) {
		if (outfile) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR:  No output grid file allowed with -S.\n", GMT_program);
			error++;
		}
	}
	else {
		if (!grdfile) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify input file.\n", GMT_program);
			error++;
		}
		if (!outfile) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify output file.\n", GMT_program);
			error++;
		}
	}
	if (skip_grid && !outline) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR: -N and -S cannot both be given.\n", GMT_program);
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

	GMT_lat_swap_init ();	/* Initialize auxiliary latitude machinery */

	/* Check limits and get data file */

	if (grdfile) {

		GMT_grd_init (&g_head, argc, argv, FALSE);

		GMT_err_fail (GMT_read_grd_info (grdfile, &g_head), grdfile);
	
		/* Determine what wesn to pass to map_setup */

		if (!project_info.region_supplied) {
			west = g_head.x_min;
			east = g_head.x_max;
			south = g_head.y_min;
			north = g_head.y_max;
		}
	
		project_info.w = west;	project_info.e = east;
		project_info.s = south;	project_info.n = north;

		/* Determine the wesn to be used to read the grdfile; or exit if file is outside -R */

		if (GMT_grd_setregion (&g_head, &data_west, &data_east, &data_south, &data_north, interpolant)) {
			fprintf (stderr, "%s: No grid values inside selected region - aborting\n", GMT_program);
			exit (EXIT_FAILURE);
		}
		global = (GMT_IS_ZERO (data_east - data_west - 360.0) && GMT_IS_ZERO (data_north - data_south - 180.0));
	}
	if (!skip_grid) {
		/* Read data */

		if (gmtdefs.verbose) fprintf (stderr, "%s: Allocates memory and read data file\n", GMT_program);
		nx = irint ( (data_east - data_west) / g_head.x_inc) + !g_head.node_offset;
		ny = irint ( (data_north - data_south) / g_head.y_inc) + !g_head.node_offset;
		mx = nx + 4;	my = ny + 4;
		GMT_pad[0] = GMT_pad[1] = GMT_pad[2] = GMT_pad[3] = 2;
		nm = GMT_get_nm (mx, my);

		original_grid = (float *) GMT_memory (VNULL, (size_t)nm, sizeof (float), GMT_program);
		GMT_err_fail (GMT_read_grd (grdfile, &g_head, original_grid, data_west, data_east, data_south, data_north, GMT_pad, FALSE), grdfile);

		GMT_boundcond_init (&edgeinfo);
		if (GMT_360_RANGE (g_head.x_max, g_head.x_min)) GMT_boundcond_parse (&edgeinfo, "g");
	
		GMT_boundcond_param_prep (&g_head, &edgeinfo);

		/* Initialize bcr structure:  */

		GMT_bcr_init (&g_head, GMT_pad, interpolant, threshold, &bcr);
		if (threshold == 0.0) node = TRUE;

		/* Set boundary conditions  */

		GMT_boundcond_set (&g_head, &edgeinfo, GMT_pad, original_grid);
	}
	
	if (user_polygon_given)	/* Read the user's polygon file */
		GMT_import_table ((void *)pol_file, GMT_IS_FILE, &pol, 0.0, (west < 0.0 && east > 0.0), TRUE, TRUE);
	else if (!global)				/* Make a single grid-outline polygon */
		get_grid_path (&g_head, &pol);

	spotter_make_rot_matrix (x_rot, y_rot, w_rot, R);	/* Make rotation matrix from rotation parameters */
	
	if (gmtdefs.verbose) fprintf (stderr, "%s: Reconstruct polygon outline\n", GMT_program);
	
	/* First reconstruct the polygon outline */
	
	ix = (gmtdefs.xy_toggle[1]);	iy = 1 - ix;	/* These are used for output purposes only */
	west = w360 = w180 = DBL_MAX;	east = e360 = e180 = -DBL_MAX;	south = DBL_MAX;	north = -DBL_MAX;
	for (i = 0; !global && i < pol->n_segments; i++) {
		if (outline) GMT_write_segmentheader (GMT_stdout, 2);
		wp360 = wp180 = DBL_MAX;	ep360 = ep180 = -DBL_MAX;	spol = DBL_MAX;	npol = -DBL_MAX;
		for (j = 0; j < pol->segment[i]->n_rows; j++) {
			pol->segment[i]->coord[GMT_Y][j] = GMT_lat_swap (pol->segment[i]->coord[GMT_Y][j], GMT_LATSWAP_G2O);	/* Convert to geocentric */
			GMT_geo_to_cart (pol->segment[i]->coord[GMT_Y][j], pol->segment[i]->coord[GMT_X][j], P_original, TRUE);	/* Convert to a Cartesian x,y,z vector; TRUE since we have degrees */
			spotter_matrix_vect_mult (R, P_original, P_rotated);				/* Rotate the vector */
			GMT_cart_to_geo (&pol->segment[i]->coord[GMT_Y][j], &pol->segment[i]->coord[GMT_X][j], P_rotated, TRUE);	/* Recover lon lat representation; TRUE to get degrees */
			pol->segment[i]->coord[GMT_Y][j] = GMT_lat_swap (pol->segment[i]->coord[GMT_Y][j], GMT_LATSWAP_O2G);	/* Convert back to geodetic */
			if (pol->segment[i]->coord[GMT_Y][j] < south) south = pol->segment[i]->coord[GMT_Y][j];
			if (pol->segment[i]->coord[GMT_Y][j] > north) north = pol->segment[i]->coord[GMT_Y][j];
			if (pol->segment[i]->coord[GMT_Y][j] < spol) spol = pol->segment[i]->coord[GMT_Y][j];
			if (pol->segment[i]->coord[GMT_Y][j] > npol) npol = pol->segment[i]->coord[GMT_Y][j];
			lon = pol->segment[i]->coord[GMT_X][j];
			GMT_lon_range_adjust (0, &lon);	/* 0 <= lon < 360 */
			if (lon < w360) w360 = lon;
			if (lon > e360) e360 = lon;
			if (lon < wp360) wp360 = lon;
			if (lon > ep360) ep360 = lon;
			GMT_lon_range_adjust (2, &lon);	/* -180 <= lon < 180 */
			if (lon < w180) w180 = lon;
			if (lon > e180) e180 = lon;
			if (lon < wp180) wp180 = lon;
			if (lon > ep180) ep180 = lon;
			out[ix] = pol->segment[i]->coord[GMT_X][j];	out[iy] = pol->segment[i]->coord[GMT_Y][j];
			if (outline) GMT_output (GMT_stdout, 2, out);
		}
		pol->segment[i]->pole = 0;
		if (fabs (ep360 - wp360) < 180.0) {
			pol->segment[i]->min[GMT_X] = w360;
			pol->segment[i]->max[GMT_X] = e360;
		}
		else if (fabs (ep180 - wp180) < 180.0) {
			pol->segment[i]->min[GMT_X] = w180;
			pol->segment[i]->max[GMT_X] = e180;
		}
		else {
			pol->segment[i]->min[GMT_X] = 0.0;
			pol->segment[i]->max[GMT_X] = 360.0;
			pol->segment[i]->pole = irint (copysign (1.0, npol));
		}
		pol->segment[i]->min[GMT_Y] = spol;	pol->segment[i]->max[GMT_Y] = npol;
	}

	if (skip_grid) {
		if (!global) GMT_free_table (pol);
	
		if (gmtdefs.verbose) fprintf (stderr, "%s: Done!\n", GMT_program);
		GMT_end (argc, argv);
	
		exit (EXIT_SUCCESS);
	}
	
	/* Then, find min/max of reconstructed outline */
	
	if (fabs (e360 - w360) < 180.0) {
		west = w360;
		east = e360;
	}
	else if (fabs (e180 - w180) < 180.0) {
		west = w180;
		east = e180;
	}
	else {
		west = 0.0;
		east = 360.0;
	}
	if (global) {
		west = g_head.x_min;	east = g_head.x_max;
		south = g_head.y_min;	north = g_head.y_max;
	}
	/* Adjust longitude range, as indicated by OUTPUT_DEGREE_FORMAT */
	GMT_lon_range_adjust (GMT_io.geo.range, &west);
	GMT_lon_range_adjust (GMT_io.geo.range, &east);
	if (west >= east) east += 360.0;
	GMT_grd_init (&r_head, argc, argv, FALSE);
	r_head.node_offset = g_head.node_offset;
	r_head.x_inc = g_head.x_inc;
	r_head.y_inc = g_head.y_inc;
	r_head.x_min = floor  (west / r_head.x_inc) * r_head.x_inc;
	r_head.x_max = ceil   (east / r_head.x_inc) * r_head.x_inc;
	r_head.y_min = floor (south / r_head.y_inc) * r_head.y_inc;
	r_head.y_max = ceil  (north / r_head.y_inc) * r_head.y_inc;
	r_head.nx = irint ((r_head.x_max - r_head.x_min) / r_head.x_inc) + !r_head.node_offset;
	r_head.ny = irint ((r_head.y_max - r_head.y_min) / r_head.y_inc) + !r_head.node_offset;
	rotated_grid = (float *) GMT_memory (VNULL, (size_t)(r_head.nx * r_head.ny), sizeof (float), GMT_program);
	grd_x = (double *) GMT_memory (VNULL, (size_t)r_head.nx, sizeof (double), GMT_program);
	grd_y = (double *) GMT_memory (VNULL, (size_t)r_head.ny, sizeof (double), GMT_program);
	grd_yc = (double *) GMT_memory (VNULL, (size_t)r_head.ny, sizeof (double), GMT_program);
	r_head.xy_off = g_head.xy_off;
	/* Precalculate node coordinates in both degrees and radians */
	for (i = 0; i < r_head.nx; i++) grd_x[i] = GMT_i_to_x (i, r_head.x_min, r_head.x_max, r_head.x_inc, r_head.xy_off, r_head.nx);
	for (j = 0; j < r_head.ny; j++) grd_y[j] = GMT_j_to_y (j, r_head.y_min, r_head.y_max, r_head.y_inc, r_head.xy_off, r_head.ny);
	for (j = 0; j < r_head.ny; j++) grd_yc[j] = GMT_lat_swap (grd_y[j], GMT_LATSWAP_G2O);

	/* Loop over all nodes in the new rotated grid and find those inside the reconstructed polygon */
	
	if (gmtdefs.verbose) fprintf (stderr, "%s: Interpolate reconstructed grid\n", GMT_program);

	spotter_make_rot_matrix (x_rot, y_rot, -w_rot, R);	/* Make inverse rotation using negative angle */
	for (j = ij = 0; j < r_head.ny; j++) {
		for (i = 0; i < r_head.nx; i++, ij++) {
			rotated_grid[ij] = GMT_f_NaN;
			if (!global) {
				inside = FALSE;
				k = 0;
				while (k < pol->n_segments && !inside) {	/* Use degrees since function expects it */
					inside = (GMT_inonout_sphpol (grd_x[i], grd_y[j], pol->segment[k]) > 0);
					k++;
				}
				if (!inside) continue;	/* Outside the polygon(s) */
			}
			
			/* Here we are inside; get the coordinates and rotate back to original grid coordinates */
			
			GMT_geo_to_cart (grd_yc[j], grd_x[i], P_rotated, TRUE);	/* Convert degree lon,lat to a Cartesian x,y,z vector */
			spotter_matrix_vect_mult (R, P_rotated, P_original);	/* Rotate the vector */
			GMT_cart_to_geo (&yy, &xx, P_original, TRUE);		/* Recover degree lon lat representation */
			yy = GMT_lat_swap (yy, GMT_LATSWAP_O2G);		/* Convert back to geodetic */
			xx -= 360.0;
			while (xx < g_head.x_min) xx += 360.0;	/* Make sure we deal with 360 issues */
			if (node) {
				ii = GMT_x_to_i (xx, g_head.x_min, g_head.x_inc, g_head.xy_off, g_head.nx);
				jj = GMT_y_to_j (yy, g_head.y_min, g_head.y_inc, g_head.xy_off, g_head.ny);
				rotated_grid[ij] = original_grid[(jj+GMT_pad[3])*mx+ii+GMT_pad[0]];
			}
			else
				rotated_grid[ij] = (float)GMT_get_bcr_z(&g_head, xx, yy, original_grid, &edgeinfo, &bcr);
		}
	}	
	
	/* Also loop over original node locations to make sure the nearest nodes are set */

	for (i = 0; !global && i < pol->n_segments; i++) {
		for (j = 0; j < pol->segment[i]->n_rows; j++) {
			lon = pol->segment[i]->coord[GMT_X][j];
			while (lon < r_head.x_min) lon += 360.0;
			ii = GMT_x_to_i (lon, r_head.x_min, r_head.x_inc, r_head.xy_off, r_head.nx);
			jj = GMT_y_to_j (pol->segment[i]->coord[GMT_Y][j], r_head.y_min, r_head.y_inc, r_head.xy_off, r_head.ny);
			/* Visit the PAD * PAD number of cells centered on ii, jj */
			for (jjj = (jj - PAD); jjj <= (jj + PAD); jjj++) {
				if (jjj < 0 || jjj >= r_head.ny) continue;
				for (iii = (ii - PAD); iii <= (ii + PAD); iii++) {
					if (iii < 0 || iii >= r_head.nx) continue;
					ij = jjj * r_head.nx + iii;
					if (!GMT_is_fnan (rotated_grid[ij])) continue;	/* Already done this */
					GMT_geo_to_cart (grd_yc[jjj], grd_x[iii], P_rotated, TRUE);	/* Convert degree lon,lat to a Cartesian x,y,z vector */
					spotter_matrix_vect_mult (R, P_rotated, P_original);	/* Rotate the vector */
					GMT_cart_to_geo (&xx, &yy, P_original, TRUE);		/* Recover degree lon lat representation */
					yy = GMT_lat_swap (yy, GMT_LATSWAP_O2G);		/* Convert back to geodetic */
					ii_o = GMT_x_to_i (xx, g_head.x_min, g_head.x_inc, g_head.xy_off, g_head.nx);
					if (ii_o < 0 || ii_o >= g_head.nx) continue;
					jj_o = GMT_y_to_j (yy, g_head.y_min, g_head.y_inc, g_head.xy_off, g_head.ny);
					if (jj_o < 0 || jj_o >= g_head.ny) continue;
					rotated_grid[ij] = original_grid[(jj_o+GMT_pad[3])*mx+ii_o+GMT_pad[0]];
				}
			}
		}
	}

	/* Now write rotated grid */
	
	if (gmtdefs.verbose) fprintf (stderr, "%s: Write reconstructed grid\n", GMT_program);

	GMT_pad[0] = GMT_pad[1] = GMT_pad[2] = GMT_pad[3] = 0;	/* No pad for output grid */
	strcpy (r_head.x_units, "degree");
	strcpy (r_head.y_units, "degree");
	sprintf (r_head.remark, "Grid rotated using lon lat omega = %g %g %g", x_rot, y_rot, w_rot);
	GMT_err_fail (GMT_write_grd (outfile, &r_head, rotated_grid, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), outfile);

	GMT_free ((void *)original_grid);
	GMT_free ((void *)rotated_grid);
	GMT_free ((void *)grd_x);
	GMT_free ((void *)grd_y);
	GMT_free ((void *)grd_yc);
	
	if (!global) GMT_free_table (pol);

	if (gmtdefs.verbose) fprintf (stderr, "%s: Done!\n", GMT_program);
	GMT_end (argc, argv);
	
	exit (EXIT_SUCCESS);
}

void get_grid_path (struct GRD_HEADER *h, struct GMT_TABLE **p)
{
	/* Return a single polygon that encloses this geographic grid exactly.
	 * It is used in the case when no particular clip polygon has been given.
	 * Note that the path is the same for pixel or grid-registered grids.
	 */

	GMT_LONG np = 0, add, i, j;
	struct GMT_TABLE *pol;
	
	pol = (struct GMT_TABLE *)GMT_memory (VNULL, (size_t)1, sizeof (struct GMT_TABLE), GMT_program);
	pol->segment = (struct GMT_LINE_SEGMENT **) GMT_memory (VNULL, (size_t)1, sizeof (struct GMT_LINE_SEGMENT *), GMT_program);
	pol->segment[0] = (struct GMT_LINE_SEGMENT *) GMT_memory (VNULL, (size_t)1, sizeof (struct GMT_LINE_SEGMENT), GMT_program);
	pol->segment[0]->coord = (double **)GMT_memory (VNULL, (size_t)2, sizeof (double *), GMT_program);
	pol->segment[0]->min = (double *)GMT_memory (VNULL, (size_t)2, sizeof (double), GMT_program);
	pol->segment[0]->max = (double *)GMT_memory (VNULL,(size_t) 2, sizeof (double), GMT_program);
		
	/* Add south border w->e */
	if (h->y_min == -90.0) {	/* If at the S pole we just add it twice for end longitudes */
		add = 2;
		pol->segment[0]->coord[GMT_X] = (double *)GMT_memory (VNULL, (size_t)add, sizeof (double), GMT_program);
		pol->segment[0]->coord[GMT_Y] = (double *)GMT_memory (VNULL, (size_t)add, sizeof (double), GMT_program);
		pol->segment[0]->coord[GMT_X][0] = h->x_min;	pol->segment[0]->coord[GMT_X][1] = h->x_max;
		pol->segment[0]->coord[GMT_Y][0] = pol->segment[0]->coord[GMT_Y][1] = h->y_min;
	}
	else {				/* Loop along south border from west to east */
		add = h->nx - !h->node_offset;
		pol->segment[0]->coord[GMT_X] = (double *)GMT_memory (VNULL, (size_t)add, sizeof (double), GMT_program);
		pol->segment[0]->coord[GMT_Y] = (double *)GMT_memory (VNULL, (size_t)add, sizeof (double), GMT_program);
		for (i = 0; i < add; i++) {
			pol->segment[0]->coord[GMT_X][i] = GMT_i_to_x (i, h->x_min, h->x_max, h->x_inc, 0.0, h->nx);
			pol->segment[0]->coord[GMT_Y][i] = h->y_min;
		}
	}
	np += add;
	/* Add east border s->n */
	add = h->ny - !h->node_offset;
	pol->segment[0]->coord[GMT_X] = (double *)GMT_memory ((void *)pol->segment[0]->coord[GMT_X], (size_t)(add + np), sizeof (double), GMT_program);
	pol->segment[0]->coord[GMT_Y] = (double *)GMT_memory ((void *)pol->segment[0]->coord[GMT_Y], (size_t)(add + np), sizeof (double), GMT_program);
	for (j = 0; j < add; j++) {	/* Loop along east border from south to north */
		pol->segment[0]->coord[GMT_X][np+j] = h->x_max;
		pol->segment[0]->coord[GMT_Y][np+j] = GMT_j_to_y (h->ny - 1 - j, h->y_min, h->y_max, h->y_inc, 0.0, h->ny);
	}
	np += add;
	/* Add north border e->w */
	if (h->y_max == 90.0) {	/* If at the N pole we just add it twice for end longitudes */
		add = 2;
		pol->segment[0]->coord[GMT_X] = (double *)GMT_memory ((void *)pol->segment[0]->coord[GMT_X], (size_t)(add + np), sizeof (double), GMT_program);
		pol->segment[0]->coord[GMT_Y] = (double *)GMT_memory ((void *)pol->segment[0]->coord[GMT_Y], (size_t)(add + np), sizeof (double), GMT_program);
		pol->segment[0]->coord[GMT_X][np] = h->x_max;	pol->segment[0]->coord[GMT_X][np+1] = h->x_min;
		pol->segment[0]->coord[GMT_Y][np] = pol->segment[0]->coord[GMT_Y][np+1] = h->y_max;
	}
	else {			/* Loop along north border from east to west */
		add = h->nx - !h->node_offset;
		pol->segment[0]->coord[GMT_X] = (double *)GMT_memory ((void *)pol->segment[0]->coord[GMT_X], (size_t)(add + np), sizeof (double), GMT_program);
		pol->segment[0]->coord[GMT_Y] = (double *)GMT_memory ((void *)pol->segment[0]->coord[GMT_Y], (size_t)(add + np), sizeof (double), GMT_program);
		for (i = 0; i < add; i++) {
			pol->segment[0]->coord[GMT_X][np+i] = GMT_i_to_x (h->nx - 1 - i, h->x_min, h->x_max, h->x_inc, 0.0, h->nx);
			pol->segment[0]->coord[GMT_Y][np+i] = h->y_max;
		}
	}
	np += add;
	/* Add west border n->s */
	add = h->ny - !h->node_offset;
	pol->segment[0]->coord[GMT_X] = (double *)GMT_memory ((void *)pol->segment[0]->coord[GMT_X], (size_t)(add + np + 1), sizeof (double), GMT_program);
	pol->segment[0]->coord[GMT_Y] = (double *)GMT_memory ((void *)pol->segment[0]->coord[GMT_Y], (size_t)(add + np + 1), sizeof (double), GMT_program);
	for (j = 0; j < add; j++) {	/* Loop along west border from north to south */
		pol->segment[0]->coord[GMT_X][np+j] = h->x_min;
		pol->segment[0]->coord[GMT_Y][np+j] = GMT_j_to_y (j, h->y_min, h->y_max, h->y_inc, 0.0, h->ny);
	}
	np += add;
	pol->segment[0]->coord[GMT_X][np] = pol->segment[0]->coord[GMT_X][0];	/* Close polygon explicitly */
	pol->segment[0]->coord[GMT_Y][np] = pol->segment[0]->coord[GMT_Y][0];
	np++;
	pol->segment[0]->n_rows = np;
	pol->segment[0]->min[GMT_X] = h->x_min;	pol->segment[0]->max[GMT_X] = h->x_max;
	pol->segment[0]->min[GMT_Y] = h->y_min;	pol->segment[0]->max[GMT_Y] = h->y_max;
	pol->segment[0]->pole = 0;
	pol->n_segments = 1;
	
	*p = pol;
}
