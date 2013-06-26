/*--------------------------------------------------------------------
 *	$Id: grdmask.c 9923 2012-12-18 20:45:53Z pwessel $
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
 * grdmask defines a grid based on region and xinc/yinc values,
 * reads xy polygon files, and sets the grid nodes inside, on the
 * boundary, and outside of the polygons to the user-defined values
 * <in>, <on>, and <out>.  These may be any number, including NaN.
 *
 * Author:	Walter. H. F. Smith
 * Date:	23-May-1991
 * Version:	2.0
 * Modified:	PW: 12-MAY-1998, to GMT 3.1
 * 		PW: 18-FEB-2000, to GMT 3.3.4: Handle polarcap polygons
 *				and for -S check if points are inside -R.
 *				Also added -L for geographic data
 * 		PW: 13-JUL-2000, 3.3.5
 *		PW: 20-AUG-2002, Added -A and do fix_up path as in psxy
 * Version:	4
 */
 
#include "gmt.h"

#define GRDMASK_N_CLASSES	3	/* outside, on edge, and inside */
#define GRDMASK_OUTSIDE		0
#define GRDMASK_ONEDGE		1
#define GRDMASK_INSIDE		2

struct GRDMASK_CTRL {
	struct A {	/* -A[m|p|step] */
		GMT_LONG active;
		GMT_LONG mode;
		double step;
	} A;
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
	struct N {	/* -N<maskvalues> */
		GMT_LONG active;
		double mask[GRDMASK_N_CLASSES];	/* values for each level */
	} N;
	struct S {	/* -S<radius>[m|c|k|K] */
		GMT_LONG active;
		double radius;
		char unit;
	} S;
};

int main (int argc, char **argv)
{
	GMT_LONG	error = FALSE, done, nofile = TRUE, periodic = FALSE, resample = FALSE;

	char line[BUFSIZ], ptr[BUFSIZ];

	GMT_LONG	i, j, side, fno, n_files = 0, n_args;
	GMT_LONG	di, dj, i0, j0, n_fields, n_expected_fields, pos;
	GMT_LONG 	distance_flag = 0, n_pol = 0;
	GMT_LONG	k, ij, nm, n_path, n_read, n_alloc = GMT_CHUNK;

	float *data;

	double	*x = NULL, *y = NULL, xx, yy, xmin, xmax, ymin, ymax, lon_sum, lat_sum, dlon;
	double distance, shrink = 1.0, idx, idy, x0, y0, *in = NULL;
	double xmin1, xmin2, xmax1, xmax2, lon_w;

	FILE *fp = NULL;

	struct GRD_HEADER header;
	struct GMT_LINE_SEGMENT P;
	struct GRDMASK_CTRL *Ctrl = NULL;
	
	void *New_grdmask_Ctrl (), Free_grdmask_Ctrl (struct GRDMASK_CTRL *C);
	
	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct GRDMASK_CTRL *)New_grdmask_Ctrl ();	/* Allocate and initialize a new control structure */

	GMT_grd_init (&header, argc, argv, FALSE);

	/* Check command line arguments */

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
				case 'f':
				case 'm':
				case '\0':
					error += GMT_parse_common_options (argv[i], &header.x_min, &header.x_max, &header.y_min, &header.y_max);
					break;

				/* Supplemental parameters */

				case 'A':	/* Turn off draw_arc mode */
					if (argv[i][2] == 'm')
						Ctrl->A.mode = 1;
					else if (argv[i][2] == 'p')
						Ctrl->A.mode = 2;
					else if (argv[i][2])
						Ctrl->A.step = atof (&argv[i][2]);	/* Undocumented test feature */
					else
						Ctrl->A.active = TRUE;
					break;
				case 'N':
					Ctrl->N.active = TRUE;
					j = pos = 0;
					while (j < GRDMASK_N_CLASSES && (GMT_strtok (&argv[i][2], "/", &pos, ptr))) {
						Ctrl->N.mask[j] = (ptr[0] == 'N' || ptr[0] == 'n') ? GMT_f_NaN : (float)atof (ptr);
						j++;
					}
					break;
				case 'F':
					Ctrl->F.active = TRUE;
					break;
				case 'G':
					Ctrl->G.active = TRUE;
					Ctrl->G.file = strdup (&argv[i][2]);
					break;
				case 'I':
					Ctrl->I.active = TRUE;
					if (GMT_getinc (&argv[i][2], &Ctrl->I.xinc, &Ctrl->I.yinc)) {
						GMT_inc_syntax ('I', 1);
						error = TRUE;
					}
					break;
				case 'L':	/* Obsolete, but backward compatibility prevails [use -f instead] */
					GMT_io.in_col_type[GMT_X] = GMT_io.out_col_type[GMT_X] = GMT_IS_LON;
					GMT_io.in_col_type[GMT_Y] = GMT_io.out_col_type[GMT_Y] = GMT_IS_LAT;
					fprintf (stderr, "%s: Option -L is obsolete (but is processed correctly).  Please use -f instead\n", GMT_program);
					break;
				case 'S':
					Ctrl->S.active = TRUE;
					Ctrl->S.radius = GMT_getradius (&argv[i][2]);
					Ctrl->S.unit = argv[i][strlen(argv[i])-1];
					break;
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
		fprintf (stderr, "grdmask %s - Create mask grid file from polygons or point coverage\n\n", GMT_VERSION);
		fprintf (stderr, "usage: grdmask [<xyfiles>] -G<mask_grd_file> %s\n", GMT_I_OPT);
		fprintf (stderr, "\t%s [-A[m|p]] [-F] [%s] [-N<out>/<edge>/<in>]\n", GMT_Rgeo_OPT, GMT_H_OPT);
		fprintf (stderr, "\t[-S<radius>[m|c|k|K] [-V] [%s] [%s] [%s] [%s]\n\n", GMT_t_OPT, GMT_bi_OPT, GMT_f_OPT, GMT_mo_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\txyfiles is one or more polygon [or point] files\n");
		fprintf (stderr, "\t-G Specify file name for output mask grid file.\n");
		GMT_inc_syntax ('I', 0);
		GMT_explain_option ('R');
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-A Suppress connecting points using great circle arcs, i.e., connect by straight lines\n");
		fprintf (stderr, "\t   unless m or p is appended to first follow meridian then parallel, or vice versa.\n");
		fprintf (stderr, "\t-F Force pixel registration for output grid [Default is gridline orientation].\n");
		GMT_explain_option ('H');
		fprintf (stderr, "\t-N sets values to use if point is outside, on the path, or inside.\n");
		fprintf (stderr, "\t   NaN is a valid entry.  Default values are 0/0/1.\n");
		fprintf (stderr, "\t-S sets search radius in -R, -I units; append m or c for minutes or seconds.\n");
		fprintf (stderr, "\t   This means input data are points and the mask nodes are set to <in> or <out> depending on\n");
		fprintf (stderr, "\t   whether they are inside a circle of specified radius [0] from the nearest data point.\n");
		fprintf (stderr, "\t   Append k for km (implies -R,-I in degrees), use flat Earth approximation.\n");
		fprintf (stderr, "\t   Append K for km (implies -R,-I in degrees), use exact geodesic distances.\n");
		fprintf (stderr, "\t   If the current ELLIPSOID is spherical then great circle distances are used.\n");
		fprintf (stderr, "\t    [Default is to treat xyfiles as polygons and use inside/outside searching].\n");
		GMT_explain_option ('V');
		GMT_explain_option (':');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf (stderr, "\t    [Default is 2 input columns].\n");
		GMT_explain_option ('f');
		GMT_explain_option ('m');
		GMT_explain_option ('.');
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
	if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
        if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = 2;
        if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] < 2) {
                fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data (-bi) must have at least 2 columns\n", GMT_program);
		error++;
	}
	if (Ctrl->S.active && (Ctrl->S.radius < 0.0 || GMT_is_dnan (Ctrl->S.radius))) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Radius is NaN or negative\n", GMT_program);
		error++;
	}
	if (error) exit (EXIT_FAILURE);

	if (GMT_io.binary[GMT_IN] && gmtdefs.verbose) {
		char *type[2] = {"double", "single"};
		fprintf (stderr, "%s: Expects %ld-column %s-precision binary data\n", GMT_program, GMT_io.ncol[GMT_IN], type[GMT_io.single_precision[GMT_IN]]);
	}

	if (Ctrl->S.unit == 'k') distance_flag = 1;
	if (Ctrl->S.unit == 'K') distance_flag = 2;

	header.node_offset = (int)Ctrl->F.active;
	header.x_inc = Ctrl->I.xinc;
	header.y_inc = Ctrl->I.yinc;
	GMT_RI_prepare (&header);	/* Ensure -R -I consistency and set nx, ny */
	GMT_err_fail (GMT_grd_RI_verify (&header, 1), Ctrl->G.file);

	header.xy_off = 0.5 * header.node_offset;
	idx = 1.0 / header.x_inc;
	idy = 1.0 / header.y_inc;
	nm = GMT_get_nm (header.nx, header.ny);
	data = (float *) GMT_memory (VNULL, (size_t)nm, sizeof(float), GMT_program);
	x = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof(double), GMT_program);
	y = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof(double), GMT_program);

	sprintf (line, "%s\n", gmtdefs.d_format);
	if (gmtdefs.verbose) {
		fprintf (stderr, "%s: Nodes completely outside the polygons will be set to ", GMT_program);
		(GMT_is_fnan (Ctrl->N.mask[GRDMASK_OUTSIDE])) ? fprintf (stderr, "NaN\n") : fprintf (stderr, line, Ctrl->N.mask[GRDMASK_OUTSIDE]);
		fprintf (stderr, "%s: Nodes completely inside the polygons will be set to ", GMT_program);
		(GMT_is_fnan (Ctrl->N.mask[GRDMASK_INSIDE])) ? fprintf (stderr, "NaN\n") : fprintf (stderr, line, Ctrl->N.mask[GRDMASK_INSIDE]);
		fprintf (stderr, "%s: Nodes on the polygons boundary will be set to ", GMT_program);
		(GMT_is_fnan (Ctrl->N.mask[GRDMASK_ONEDGE])) ? fprintf (stderr, "NaN\n") : fprintf (stderr, line, Ctrl->N.mask[GRDMASK_ONEDGE]);
	}

	if (distance_flag > 0) {
		double width;
		shrink = cosd (MAX(fabs(header.y_min), fabs(header.y_max)));
		width = (shrink < GMT_CONV_LIMIT) ? DBL_MAX : ceil (Ctrl->S.radius / (project_info.DIST_KM_PR_DEG * header.x_inc * shrink));
		di = (width < (double)header.nx) ? (GMT_LONG) width : header.nx;
		dj = (GMT_LONG)ceil (Ctrl->S.radius / (project_info.DIST_KM_PR_DEG * header.y_inc));
	}
	else {
		di = (GMT_LONG)ceil (Ctrl->S.radius * idx);
		dj = (GMT_LONG)ceil (Ctrl->S.radius * idy);
	}
	periodic = (GMT_io.in_col_type[GMT_X] == GMT_IS_LON);	/* Dealing with geographic coordinates */
	resample = ((!Ctrl->A.active || Ctrl->A.mode) && periodic);
	if (distance_flag == 2 && !GMT_IS_SPHERICAL) distance_flag = 3;	/* Use geodesics */

	if (gmtdefs.verbose) {
		char *type[2] = {"Cartesian", "spherical"};
		fprintf (stderr, "%s: You have chosen to perform %s calculations\n", GMT_program, type[periodic]);
	}
	
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

	memset ((void *)&P, 0, sizeof (struct GMT_LINE_SEGMENT));
	P.coord = (double **) GMT_memory (VNULL, (size_t)2, sizeof (double *), GMT_program);	/* Needed as pointers below */
	P.min = (double *) GMT_memory (VNULL, (size_t)2, (size_t)sizeof (double), GMT_program);		/* Needed to hold min lon/lat */
	P.max = (double *) GMT_memory (VNULL, (size_t)2, (size_t)sizeof (double), GMT_program);		/* Needed to hold max lon/lat */
	
	/* Initialize all nodes to the 'outside' value */

	for (ij = 0; ij < nm; ij++) data[ij] = (float)Ctrl->N.mask[GRDMASK_OUTSIDE];

	if (n_files > 0)
		nofile = FALSE;
	else
		n_files = 1;
	n_args = (argc > 1) ? argc : 2;

	done = FALSE;
	n_expected_fields = (GMT_io.ncol[GMT_IN]) ? GMT_io.ncol[GMT_IN] : 2;
	GMT_io.skip_duplicates = TRUE;	/* The inonout algorithm assumes no duplicates */

	for (fno = 1; !done && fno < n_args; fno++) {	/* Loop over all input files */
		if (!nofile && argv[fno][0] == '-') continue;
		if (nofile) {
			fp = GMT_stdin;
			done = TRUE;
#ifdef SET_IO_MODE
			GMT_setmode (GMT_IN);
#endif
		}
		else if ((fp = GMT_fopen (argv[fno], GMT_io.r_mode)) == NULL) {
			fprintf (stderr, "%s: Cannot open file %s\n", GMT_program, argv[fno]);
			continue;
		}

		if (!nofile && gmtdefs.verbose) fprintf (stderr, "%s: Working on file %s\n", GMT_program, argv[fno]);

		if (GMT_io.io_header[GMT_IN]) for (i = 0; i < GMT_io.n_header_recs; i++) GMT_fgets (line, BUFSIZ, fp);

		n_read = 0;
		n_fields = GMT_input (fp,  &n_expected_fields, &in);
		while (! (GMT_io.status & GMT_IO_EOF)) {	/* Not yet EOF */

			while (GMT_io.status & GMT_IO_SEGMENT_HEADER && !(GMT_io.status & GMT_IO_EOF)) n_fields = GMT_input (fp, &n_expected_fields, &in);
			if ((GMT_io.status & GMT_IO_EOF)) continue;	/* At EOF */

			n_path = 0;
			xmin = ymin = +DBL_MAX;
			xmax = ymax = -DBL_MAX;
			xmin1 = xmin2 = 360.0;
			xmax1 = xmax2 = -360.0;
			lon_sum = lat_sum = 0.0;

			while (!(GMT_io.status & (GMT_IO_SEGMENT_HEADER | GMT_IO_EOF))) {	/* Keep going until FALSE  */
				n_read++;
				if (GMT_io.status & GMT_IO_MISMATCH) {
					fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld (skipped)\n", GMT_program, n_fields,  n_expected_fields, n_read);
					continue;
				}
				x[n_path] = in[GMT_X];
				y[n_path] = in[GMT_Y];

				if (y[n_path] > ymax) ymax = y[n_path];
				if (y[n_path] < ymin) ymin = y[n_path];
				if (periodic) {	/* Longitudes */
					lon_w = x[n_path];
					if (x[n_path] < 0.0) x[n_path] += 360.0;	/* Start off with everything in 0-360 range */
					if (x[n_path] < 180.0) {
						xmin1 = MIN (x[n_path], xmin1);
						xmax1 = MAX (x[n_path], xmax1);
					}
					else {
						xmin2 = MIN (lon_w, xmin2);
						xmax2 = MAX (lon_w, xmax2);
					}
					if (n_path > 0) {	/* Do longitude-difference check sum */
						dlon = x[n_path] - x[n_path-1];
						if (fabs (dlon) > 180.0) dlon = copysign (360.0 - fabs (dlon), -dlon);
						lon_sum += dlon;
					}
					lat_sum += y[n_path];
				}
				else {	/* Cartesian x */
					if (x[n_path] > xmax) xmax = x[n_path];
					if (x[n_path] < xmin) xmin = x[n_path];
				}
				n_path++;

				if (n_path == (n_alloc-1)) {	/* n_alloc-1 since we may need 1 more to close polygon */
					n_alloc <<= 1;
					x = (double *) GMT_memory ((void *)x, (size_t)n_alloc, sizeof(double), GMT_program);
					y = (double *) GMT_memory ((void *)y, (size_t)n_alloc, sizeof(double), GMT_program);
				}
				n_fields = GMT_input (fp, &n_expected_fields, &in);
			}
			n_pol++;
			if (periodic) {	/* Longitudes */
				if (xmin1 == 360.0) {		/* Only negative longitudes found */
					xmin = xmin2;
					xmax = xmax2;
				}
				else if (xmin2 == 360.0) {	/* Only positive longitudes found */
					xmin = xmin1;
					xmax = xmax1;
				}
				else if ((xmin1 - xmax2) < 90.0) {	/* Crossed Greenwich */
					xmin = xmin2;
					xmax = xmax1;
				}
				else {					/* Crossed Dateline */
					xmin = xmin1;
					xmax = xmax2;
				}
				if (xmin > xmax) xmin -= 360.0;
			}
			if (Ctrl->S.active) {	/* Assign 'inside' to nodes within given distance of data constrains */
				for (k = 0; k < n_path; k++) {

					if (GMT_y_is_outside (y[k], header.y_min, header.y_max)) continue;	/* Outside y-range */
					if (GMT_x_is_outside (&x[k], header.x_min, header.x_max)) continue;	/* Outside x-range (or longitude) */

					/* OK, this point is within bounds, but may be exactly on the border */

					i0 = GMT_x_to_i (x[k], header.x_min, header.x_inc, header.xy_off, header.nx);
					if (i0 == header.nx) i0--;	/* Was exactly on the xmax edge */
					j0 = GMT_y_to_j (y[k], header.y_min, header.y_inc, header.xy_off, header.ny);
					if (j0 == header.ny) j0--;	/* Was exactly on the ymin edge */
					data[j0 * header.nx + i0] = (float)Ctrl->N.mask[GRDMASK_INSIDE];	/* This is the nearest node */
					if (Ctrl->S.radius == 0.0) continue;

					for (j = j0 - dj; j <= (j0 + dj); j++) {
						if (j < 0 || j >= header.ny) continue;
						for (i = i0 - di; i <= (i0 + di); i++) {
							if (i < 0 || i >= header.nx) continue;
							ij = GMT_IJ (j, i, header.nx);
							x0 = GMT_i_to_x (i, header.x_min, header.x_max, header.x_inc, header.xy_off, header.nx);
							y0 = GMT_j_to_y (j, header.y_min, header.y_max, header.y_inc, header.xy_off, header.ny);
							distance = (GMT_distance_func) (x[k], y[k], x0, y0);
							if (distance > Ctrl->S.radius) continue;
							data[ij] = (float)Ctrl->N.mask[GRDMASK_INSIDE];	/* The in value */
						}
					}
				}
			}
			else if (n_path) {	/* assign 'inside' to nodes if they are inside given polygon */
				dlon = x[n_path-1] - x[0];
				if (periodic && fabs(dlon) == 360.0) dlon = 0.0;
				if (!(y[n_path-1] == y[0] && dlon == 0.0)) {	/* Explicitly close the polygon */
					x[n_path] = x[0];
					y[n_path] = y[0];
					dlon = x[n_path] - x[n_path-1];
					n_path++;
					if (periodic) {
						if (fabs (dlon) > 180.0) dlon = copysign (360.0 - fabs (dlon), -dlon);
						lon_sum += dlon;
					}
				}
				if (n_path < 3) {	/* Cannot be a polygon; skip */
					if (gmtdefs.verbose) fprintf (stderr, "%s: Polygon %ld only had %ld points; skipped\r", GMT_program, n_pol, n_path);
					continue;
				}
				if (resample) {
					n_path = GMT_fix_up_path (&x, &y, n_path, Ctrl->A.step, Ctrl->A.mode);
					n_alloc = n_path;	/* Since it got reset */
				}

				if (periodic) {	/* Make polygon structure so we can use spherical polygon machinery */
					P.coord[GMT_X] = x;
					P.coord[GMT_Y] = y;
					P.n_rows = n_path;
					P.min[GMT_Y] = ymin;	P.max[GMT_Y] = ymax;
					P.min[GMT_X] = xmin;	P.max[GMT_X] = xmax;
					if (GMT_360_RANGE (lon_sum, 0.0)) {	/* Contains a pole, convert to polar x,y */
						P.pole = (lat_sum < 0.0) ? -1 : +1;	/* S or N pole */
					}
					else
						P.pole = 0;
				}

				for (j = 0; j < header.ny; j++) {

					yy = GMT_j_to_y (j, header.y_min, header.y_max, header.y_inc, header.xy_off, header.ny);
					
					/* First check if point is outside, then there is no need to assign value */
					
					if (periodic) {	/* Containing annulus test */
						if (P.pole != +1 && yy > ymax) continue;	/* No N polar cap and beyond north */
						if (P.pole != -1 && yy < ymin) continue;	/* No S polar cap and beyond south */
					}
					else if (yy < ymin || yy > ymax)	/* Cartesian case */
						continue;

					for (i = 0; i < header.nx; i++) {
						xx = GMT_i_to_x (i, header.x_min, header.x_max, header.x_inc, header.xy_off, header.nx);
						if (periodic) {
							if (P.pole)	/* 360-degree polar cap, must check fully */
								side = GMT_inonout_sphpol (xx, yy, &P);
							else {	/* See if we are outside range of longitudes for polygon */
								while (xx > xmin) xx -= 360.0;	/* Wind clear of west */
								while (xx < xmin) xx += 360.0;	/* Wind east until inside or beyond east */
								if (xx > xmax) continue;	/* Point outside, no need to assign value */
								side = GMT_inonout_sphpol (xx, yy, &P);
							}
						}
						else {
							if (xx < xmin || xx > xmax) continue;	/* Point outside, no need to assign value */
							side = GMT_non_zero_winding (xx, yy, x, y, n_path);
						}

						if (side == 0) continue;	/* Outside */

						/* Here, point is inside or on edge, we must assign value */

						ij = GMT_IJ (j, i, header.nx);
						data[ij] = (float)Ctrl->N.mask[side];
					}
					if (gmtdefs.verbose) fprintf (stderr, "%s: Polygon %ld scanning row %5.5ld\r", GMT_program, n_pol, j);
				}
			}
		}
		if (gmtdefs.verbose) fprintf (stderr, "\n");

		if (fp != GMT_stdin) GMT_fclose (fp);
	}
	GMT_io.skip_duplicates = FALSE;	/* Reset to FALSE */

	GMT_err_fail (GMT_write_grd (Ctrl->G.file, &header, data, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), Ctrl->G.file);

	GMT_free ((void *)data);
	GMT_free ((void *)x);
	GMT_free ((void *)y);
	GMT_free ((void *)P.coord);
	GMT_free ((void *)P.min);
	GMT_free ((void *)P.max);

	Free_grdmask_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void *New_grdmask_Ctrl () {	/* Allocate and initialize a new control structure */
	struct GRDMASK_CTRL *C;
	
	C = (struct GRDMASK_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GRDMASK_CTRL), "New_grdmask_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	C->A.step = 0.1;	/* In degrees */
	C->N.mask[GRDMASK_INSIDE] = 1.0;			/* Default inside value */
	return ((void *)C);
}

void Free_grdmask_Ctrl (struct GRDMASK_CTRL *C) {	/* Deallocate control structure */
	if (C->G.file) free ((void *)C->G.file);	
	GMT_free ((void *)C);	
}
