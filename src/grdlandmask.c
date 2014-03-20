/*--------------------------------------------------------------------
 *	$Id: grdlandmask.c 10173 2014-01-01 09:52:34Z pwessel $
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
 * grdlandmask defines a grid based on region and xinc/yinc values,
 * reads a shoreline data base, and sets the grid nodes inside, on the
 * boundary, and outside of the polygons to the user-defined values
 * <in>, <on>, and <out>.  These may be any number, including NaN.
 *
 * Author:	P. Wessel
 * Date:	23-Sep-1994
 * Version:	3.0
 * Modified:	24-JUN-1998, for GMT 3.1
 *		18-AUG-1999, for GMT 3.3.2
 *		13-JUL-2000, for GMT 3.3.5
 * Version:	4
 */
 
#include "gmt.h"

#define GRDLANDMASK_N_CLASSES	(GMT_MAX_GSHHS_LEVEL + 1)	/* Number of bands separated by the levels */

struct GRDLANDMASK_CTRL {	/* All control options for this program (except common args) */
	/* ctive is TRUE if the option has been activated */
	struct A {	/* -A<min_area>[/<min_level>/<max_level>] */
		GMT_LONG active;
		struct GMT_SHORE_SELECT info;
	} A;
	struct D {	/* -D<resolution> */
		GMT_LONG active;
		GMT_LONG force;	/* if TRUE, select next highest level if current set is not avaialble */
		char set;	/* One of f, h, i, l, c */
	} D;
	struct F {	/* -F */
		GMT_LONG active;
	} F;
	struct G {	/* -G<maskfile> */
		GMT_LONG active;
		char *file;
	} G;
	struct I {	/* -Idx[/dy] */
		GMT_LONG active;
		double xinc, yinc;
	} I;
	struct N {	/* -N<maskvalues>[o] */
		GMT_LONG active;
		GMT_LONG edge;	/* TRUE if edges are considere outside */
		double mask[GRDLANDMASK_N_CLASSES];	/* values for each level */
	} N;
};

int main (int argc, char **argv)
{
	GMT_LONG	error = FALSE, dry_wet_only = FALSE, greenwich = FALSE;
	GMT_LONG temp_shift = FALSE, wrap, used_polygons;

	char line[GMT_LONG_TEXT], ptr[BUFSIZ];
	char *shore_resolution[5] = {"full", "high", "intermediate", "low", "crude"};

	GMT_LONG	i, j, k, ii, bin, ind, np, side, i_min, i_max, j_min, j_max, nx1, ny1, np_new, pos;
	GMT_LONG	base = 3, direction, is_inside = 1, err, ij, nm;

	float *data;

	double	*x = NULL, *y = NULL, xmin, xmax, ymin, ymax, west_border, east_border, i_dx_inch, i_dy_inch;
	double i_dx, i_dy, del_off, dummy;

	struct GMT_SHORE c;
	struct GRD_HEADER header;
	struct GMT_GSHHS_POL *p = NULL;
	struct GRDLANDMASK_CTRL *Ctrl = NULL;

	void *New_grdlandmask_Ctrl (), Free_grdlandmask_Ctrl (struct GRDLANDMASK_CTRL *C);

	argc = (int)GMT_begin (argc, argv);
	GMT_io.out_col_type[0] = GMT_IS_LON;
	GMT_io.out_col_type[1] = GMT_IS_LAT;

	Ctrl = (struct GRDLANDMASK_CTRL *) New_grdlandmask_Ctrl ();		/* Allocate and initialize defaults in a new control structure */

	GMT_grd_init (&header, argc, argv, FALSE);

	/* Check command line arguments */

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {

				/* Common parameters */

				case 'R':
				case 'V':
				case ':':
				case '\0':
					error += GMT_parse_common_options (argv[i], &header.x_min, &header.x_max, &header.y_min, &header.y_max);
					break;

				/* Supplemental parameters */

				case 'A':
					Ctrl->A.active = TRUE;
					GMT_set_levels (&argv[i][2], &Ctrl->A.info);
					break;
				case 'D':
					Ctrl->D.active = TRUE;
					Ctrl->D.set = argv[i][2];
					if (argv[i][3] == '+') Ctrl->D.force = TRUE;
					break;
				case 'N':
					Ctrl->N.active = TRUE;
					strcpy (line, &argv[i][2]);
					if (line[strlen(line)-1] == 'o') { /* Edge is considered outside */
						Ctrl->N.edge = TRUE;
						line[strlen(line)-1] = 0;
					}
					j = pos = 0;
					while (j < 5 && (GMT_strtok (line, "/", &pos, ptr))) {
						Ctrl->N.mask[j] = (ptr[0] == 'N' || ptr[0] == 'n') ? GMT_f_NaN : (float)atof (ptr);
						j++;
					}
					if (!(j == 2 || j == 5)) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -N option:  Specify 2 or 5 arguments\n", GMT_program);
						exit (EXIT_FAILURE);
					}
					dry_wet_only = (j == 2);
					break;
				case 'F':
					Ctrl->F.active = TRUE;
					break;
				case 'G':
					Ctrl->G.active = TRUE;
					Ctrl->G.file = strdup (&argv[i][2]);
					break;
				case 'I':
					if (GMT_getinc (&argv[i][2], &Ctrl->I.xinc, &Ctrl->I.yinc)) {
						GMT_inc_syntax ('I', 1);
						error = TRUE;
					}
					Ctrl->I.active = TRUE;
					break;
				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "grdlandmask %s - Create \"wet-dry\" mask grid file from shoreline data base\n\n", GMT_VERSION);
		fprintf (stderr, "usage: grdlandmask -G<mask_grd_file> %s %s\n", GMT_I_OPT, GMT_Rgeo_OPT);
		fprintf (stderr, "\t[%s] [-D<resolution>][+] [-F] [-N<maskvalues>[o]] [-V]\n\n", GMT_A_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\t-G Specify file name for output mask grid file.\n");
		GMT_inc_syntax ('I', 0);
		GMT_explain_option ('R');
		fprintf (stderr, "\n\tOPTIONS:\n");
		GMT_GSHHS_syntax ('A', "Place limits on coastline features from the GSHHS data base.");
		fprintf (stderr, "\t-D Choose one of the following resolutions:\n");
		fprintf (stderr, "\t   f - full resolution (may be very slow for large regions).\n");
		fprintf (stderr, "\t   h - high resolution (may be slow for large regions).\n");
		fprintf (stderr, "\t   i - intermediate resolution.\n");
		fprintf (stderr, "\t   l - low resolution [Default].\n");
		fprintf (stderr, "\t   c - crude resolution, for tasks that need crude continent outlines only.\n");
		fprintf (stderr, "\t   Append + to use a lower resolution should the chosen one not be available [abort].\n");
		fprintf (stderr, "\t-F Force pixel registration for output grid [Default is gridline orientation].\n");
		fprintf (stderr, "\t-N gives values to use if a node is outside or inside a feature.\n");
		fprintf (stderr, "\t   Append o to let feature boundary be considered outside [Default is inside].\n");
		fprintf (stderr, "\t   Specify this information using 1 of 2 formats:\n");
		fprintf (stderr, "\t   -N<wet>/<dry>.\n");
		fprintf (stderr, "\t   -N<ocean>/<land>/<lake>/<island>/<pond>.\n");
		fprintf (stderr, "\t   NaN is a valid entry.  Default values are 0/1/0/1/0 (i.e., 0/1).\n");
		GMT_explain_option ('V');
		exit (EXIT_FAILURE);
	}

	GMT_check_lattice (&Ctrl->I.xinc, &Ctrl->I.yinc, &Ctrl->F.active, &Ctrl->I.active);

	if (!project_info.region_supplied) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify -R option\n", GMT_program);
		error++;
	}
	if (Ctrl->I.xinc <= 0.0 || Ctrl->I.yinc <= 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -I option.  Must specify positive increment(s)\n", GMT_program);
		error = TRUE;
	}
	if (!Ctrl->G.file) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -G:  Must specify output file\n", GMT_program);
		error = TRUE;
	}

	if (error) exit (EXIT_FAILURE);

	header.x_inc = Ctrl->I.xinc;
	header.y_inc = Ctrl->I.yinc;
	header.node_offset = (int)Ctrl->F.active;
	
	GMT_RI_prepare (&header);	/* Ensure -R -I consistency and set nx, ny */

	if (header.x_min < 0.0 && header.x_max < 0.0) {	/* Shift longitudes */
		temp_shift = TRUE;
		header.x_min += 360.0;
		header.x_max += 360.0;
	}

	GMT_err_fail (GMT_grd_RI_verify (&header, 1), Ctrl->G.file);

	if (Ctrl->D.force) Ctrl->D.set = GMT_shore_adjust_res (Ctrl->D.set);
	base = GMT_set_resolution (&Ctrl->D.set, 'D');
	
	is_inside = (Ctrl->N.edge) ? 2 : 1;		/* Whether of not being exactly on an edge is outside */
	if (dry_wet_only) {
		Ctrl->N.mask[3] = Ctrl->N.mask[1];
		Ctrl->N.mask[2] = Ctrl->N.mask[4] = Ctrl->N.mask[0];
	}

	if (GMT_init_shore (Ctrl->D.set, &c, header.x_min, header.x_max, header.y_min, header.y_max, &Ctrl->A.info)) {
		fprintf (stderr, "%s: %s resolution shoreline data base not installed\n", GMT_program, shore_resolution[base]);
		exit (EXIT_FAILURE);
	}
	if (gmtdefs.verbose == 2) fprintf (stderr, "GSHHS version %s\n%s\n%s\n", c.version, c.title, c.source);

	sprintf (line, "%s\n", gmtdefs.d_format);
	if (gmtdefs.verbose && dry_wet_only) {
		fprintf (stderr, "%s: Nodes in water will be set to ", GMT_program);
		(GMT_is_fnan (Ctrl->N.mask[0])) ? fprintf (stderr, "NaN\n") : fprintf (stderr, line, Ctrl->N.mask[0]);
		fprintf (stderr, "%s: Nodes on land will be set to ", GMT_program);
		(GMT_is_fnan (Ctrl->N.mask[1])) ? fprintf (stderr, "NaN\n") : fprintf (stderr, line, Ctrl->N.mask[1]);
	}
	else if (gmtdefs.verbose && !dry_wet_only) {
		fprintf (stderr, "%s: Nodes in the oceans will be set to ", GMT_program);
		(GMT_is_fnan (Ctrl->N.mask[0])) ? fprintf (stderr, "NaN\n") : fprintf (stderr, line, Ctrl->N.mask[0]);
		fprintf (stderr, "%s: Nodes on land will be set to ", GMT_program);
		(GMT_is_fnan (Ctrl->N.mask[1])) ? fprintf (stderr, "NaN\n") : fprintf (stderr, line, Ctrl->N.mask[1]);
		fprintf (stderr, "%s: Nodes in lakes will be set to ", GMT_program);
		(GMT_is_fnan (Ctrl->N.mask[2])) ? fprintf (stderr, "NaN\n") : fprintf (stderr, line, Ctrl->N.mask[2]);
		fprintf (stderr, "%s: Nodes in islands will be set to ", GMT_program);
		(GMT_is_fnan (Ctrl->N.mask[3])) ? fprintf (stderr, "NaN\n") : fprintf (stderr, line, Ctrl->N.mask[3]);
		fprintf (stderr, "%s: Nodes in ponds will be set to ", GMT_program);
		(GMT_is_fnan (Ctrl->N.mask[4])) ? fprintf (stderr, "NaN\n") : fprintf (stderr, line, Ctrl->N.mask[4]);
	}

	i_dx = 1.0 / header.x_inc;
	i_dy = 1.0 / header.y_inc;
	del_off = (header.node_offset) ? 0.5 : 0.0;
	nm = GMT_get_nm (header.nx, header.ny);
	data = (float *) GMT_memory (VNULL, (size_t)nm, sizeof(float), GMT_program);
	/* All data nodes are thus initialized to 0 */
	x = (double *) GMT_memory (VNULL, (size_t)header.nx, sizeof(double), GMT_program);
	y = (double *) GMT_memory (VNULL, (size_t)header.ny, sizeof(double), GMT_program);

	nx1 = header.nx - 1;	ny1 = header.ny - 1;

	if (header.x_min < 0.0 && header.x_max > 0.0) {	/* Must shift longitudes */
		greenwich = TRUE;
	}

	GMT_parse_J_option ("x1d");	/* Fake linear projection */
	GMT_err_fail (GMT_map_setup (header.x_min, header.x_max, header.y_min, header.y_max), "");
	GMT_parallel_straight = GMT_meridian_straight = 2;	/* No resampling along bin boundaries */
	wrap = GMT_world_map = GMT_grd_is_global (&header);
	
	/* Fill out gridnode coordinates and apply the implicit linear projection */

	for (i = 0; i < header.nx; i++) GMT_geo_to_xy (GMT_i_to_x (i, header.x_min, header.x_max, header.x_inc, header.xy_off, header.nx), 0.0, &x[i], &dummy);
	for (j = 0; j < header.ny; j++) GMT_geo_to_xy (0.0, GMT_j_to_y (j, header.y_min, header.y_max, header.y_inc, header.xy_off, header.ny), &dummy, &y[j]);
	i_dx_inch = 1.0 / fabs (x[1] - x[0]);
	i_dy_inch = 1.0 / fabs (y[1] - y[0]);

        west_border = floor (project_info.w / c.bsize) * c.bsize;
        east_border = ceil (project_info.e / c.bsize) * c.bsize;
	for (ind = 0; ind < c.nb; ind++) {	/* Loop over necessary bins only */

		bin = c.bins[ind];
		if (gmtdefs.verbose) fprintf (stderr, "%s: Working on block # %5ld\r", GMT_program, bin);

		if ((err = GMT_get_shore_bin (ind, &c))) {
			fprintf (stderr, "%s: %s [%s resolution shoreline]\n", GMT_program, GMT_strerror(err), shore_resolution[base]);
			exit (EXIT_FAILURE);
		}

		/* Use polygons, if any.  Go in both directions to cover both land and sea */

		used_polygons = FALSE;

		for (direction = -1; c.ns > 0 && direction < 2; direction += 2) {

			/* Assemble one or more segments into polygons */

			np = GMT_assemble_shore (&c, direction, TRUE, greenwich, west_border, east_border, &p);

			/* Get clipped polygons in x,y inches that can be processed */

			np_new = GMT_prep_polygons (&p, np, FALSE, 0.0, -1);

			for (k = 0; k < np_new; k++) {

				if (p[k].n == 0) continue;

				used_polygons = TRUE;	/* At least som points made it to here */

				/* Find min/max of polygon in inches */

				xmin = xmax = p[k].lon[0];
				ymin = ymax = p[k].lat[0];
				for (i = 1; i < p[k].n; i++) {
					if (p[k].lon[i] < xmin) xmin = p[k].lon[i];
					if (p[k].lon[i] > xmax) xmax = p[k].lon[i];
					if (p[k].lat[i] < ymin) ymin = p[k].lat[i];
					if (p[k].lat[i] > ymax) ymax = p[k].lat[i];
				}
				i_min = (GMT_LONG)MAX (0, ceil (xmin * i_dx_inch - del_off - GMT_CONV_LIMIT));
				if (i_min > nx1) i_min = 0;
				i_max = (GMT_LONG)MIN (nx1, floor (xmax * i_dx_inch - del_off + GMT_CONV_LIMIT));
				if (i_max <= 0 || i_max < i_min) i_max = nx1;
				j_min = (GMT_LONG)MAX (0, ceil ((project_info.ymax - ymax) * i_dy_inch - del_off - GMT_CONV_LIMIT));
				j_max = (GMT_LONG)MIN (ny1, floor ((project_info.ymax - ymin) * i_dy_inch - del_off + GMT_CONV_LIMIT));

				for (j = j_min; j <= j_max; j++) {
					for (i = i_min; i <= i_max; i++) {

						if ((side = GMT_non_zero_winding (x[i], y[j], p[k].lon, p[k].lat, p[k].n)) < is_inside) continue;	/* Outside */

						/* Here, point is inside, we must assign value */

						ij = GMT_IJ (j, i, header.nx);
						if (p[k].level > data[ij]) data[ij] = (float)p[k].level;
					}
				}
			}

			GMT_free_polygons (p, np_new);
			GMT_free ((void *)p);
		}

		if (!used_polygons) {	/* Lack of polygons or clipping etc resulted in no polygons after all, must deal with background */
			k = INT_MAX;	/* Initialize to outside range of levels (4 is highest) */
			/* Visit each of the 4 nodes, test if it is inside -R, and if so update lowest level found so far */

			if (!GMT_map_outside (c.lon_sw, c.lat_sw)) k = MIN (k, c.node_level[0]);			/* SW */
			if (!GMT_map_outside (c.lon_sw + c.bsize, c.lat_sw)) k = MIN (k, c.node_level[1]);		/* SE */
			if (!GMT_map_outside (c.lon_sw + c.bsize, c.lat_sw - c.bsize)) k = MIN (k, c.node_level[2]);	/* NE */
			if (!GMT_map_outside (c.lon_sw, c.lat_sw - c.bsize)) k = MIN (k, c.node_level[3]);		/* NW */

			/* If k is still INT_MAX we must assume this patch should have the min level of the bin */

			if (k == INT_MAX) k = MIN (MIN (c.node_level[0], c.node_level[1]) , MIN (c.node_level[2], c.node_level[3]));

			/* Determine nodes to initialize */

			j_min = (GMT_LONG)MAX (0, ceil ((header.y_max - c.lat_sw - c.bsize) * i_dy - del_off));
			j_max = (GMT_LONG)MIN (ny1, floor ((header.y_max - c.lat_sw) * i_dy - del_off));
			if (wrap) {
				i_min = (GMT_LONG)ceil (fmod (c.lon_sw - header.x_min, 360.0) * i_dx - del_off);
				i_max = (GMT_LONG)floor (fmod (c.lon_sw + c.bsize - header.x_min, 360.0) * i_dx - del_off);
				if (i_max < i_min) i_max += header.nx;
			}
			else {	/* Make sure we are inside our grid */
				double lon_w, lon_e;
				lon_w = c.lon_sw - header.x_min;	lon_e = c.lon_sw + c.bsize - header.x_min;
				if (lon_w < header.x_min && (lon_w+360.0) < header.x_max) {
					lon_w += 360.0;	lon_e += 360.0;
				}
				else if (lon_e > header.x_max && (lon_e-360.0) > header.x_min) {
					lon_w -= 360.0;	lon_e -= 360.0;
				}
				i_min = (GMT_LONG)ceil (lon_w * i_dx - del_off);
				i_max = (GMT_LONG)floor (lon_e * i_dx - del_off);
				if (i_min < 0) i_min = 0;
				if (i_max > nx1) i_max = nx1;
			}
			for (j = j_min; j <= j_max; j++) {
				for (i = i_min; i <= i_max; i++) {
					ii = (wrap) ? i % header.nx : i;
					if (ii < 0 || ii > nx1) continue;
					ij = GMT_IJ (j, ii, header.nx);
					data[ij] = (float)k;
				}
			}
		}

		GMT_free_shore (&c);
	}

	GMT_shore_cleanup (&c);
	if (gmtdefs.verbose) fprintf (stderr, "\n");

	for (ij = 0; ij < nm; ij++) {
		k = irint (data[ij]);
		data[ij] = (float)Ctrl->N.mask[k];
	}

	if (wrap && header.node_offset == 0) { /* Copy over values to the repeating right column */
		for (j = 0, ij = 0; j < header.ny; j++, ij += header.nx) data[ij+nx1] = data[ij];
	}
	
	if (temp_shift) {
		header.x_min -= 360.0;
		header.x_max -= 360.0;
	}

	GMT_io.out_col_type[0] = GMT_IS_LON;	GMT_io.out_col_type[1] = GMT_IS_LAT;	/* Since -Jx1d will say output is Cartesian */
	GMT_err_fail (GMT_write_grd (Ctrl->G.file, &header, data, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), Ctrl->G.file);

	if (gmtdefs.verbose) fprintf (stderr, "%s: Done!\n", GMT_program);

	GMT_free ((void *)data);
	GMT_free ((void *)x);
	GMT_free ((void *)y);

	Free_grdlandmask_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void *New_grdlandmask_Ctrl () {	/* Allocate and initialize a new control structure */
	struct GRDLANDMASK_CTRL *C;
	
	C = (struct GRDLANDMASK_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GRDLANDMASK_CTRL), "New_grdlandmask_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	
	C->A.info.high = GMT_MAX_GSHHS_LEVEL;				/* Include all GSHHS levels */
	C->D.set = 'l';							/* Low-resolution coastline data */
	memset ((void *)C->N.mask, 0, (size_t)(GRDLANDMASK_N_CLASSES * sizeof (float)));	/* Default "wet" value = 0 */
	C->N.mask[1] = C->N.mask[3] = 1.0;				/* Default for "dry" areas = 1 (inside) */
	
	return ((void *)C);
}

void Free_grdlandmask_Ctrl (struct GRDLANDMASK_CTRL *C) {	/* Deallocate control structure */
	if (C->G.file) free ((void *)C->G.file);	
	GMT_free ((void *)C);	
}
