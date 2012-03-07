/*--------------------------------------------------------------------
 *	$Id: gmtselect.c,v 1.119 2011/07/08 21:27:06 guru Exp $
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
 * gmtselect is a general-purpose spatial filter.  Data pass or fail based
 * on one or more conditions.  Six conditions may be set:
 *
 *	1. Only data inside a rectangular area may pass
 *	2. Only data within a certain distance from given points may pass
 *	3. Only data within a certain distance from given lines may pass
 *	4. Only data within given polygons may pass
 *	5. Only data within the coastline may pass
 *	6. Only data with z-values within specified range may pass
 *
 * Distances are calculated in the users units using Euclidian geometry
 * unless a map projection and region (-R -J) are used.  Then, distances
 * are calculated using spherical geometry and converted to km, and any
 * distances given in options or via headers are assumed to be in km.
 *
 * Any one of these conditions may be negated for the opposite result
 * Both binary and ASCII data files are accommodated
 *
 * Author:	Paul Wessel
 * Date:	25-AUG-1998
 * Version:	3.1
 *		3.2, 15-MAR-1999
 *		3.3, 13-APR-1999.	Added quick check for inside a square before checking
 *					for inside a circle (suggested by Joaquim Luis)
 *					Also added option -Cf for faster, approximate distances
 *		3.3.2 8-SEPT-1999.	Sorting -C points in x to speed up search as
 *				 	suggested by Joaquim Luis.
 *		3.3.4 17-FEB-2000.	Polygons containing either S or N pole will work correctly
 *		3.3.5 10-JUL-2000	Added plain -L for periodicity
 * 		3.4
 * Version:	4.1.2 06-APR-2006	No longer declares global variables
 * Version:	4.2.1 23-APR-2007	Extended -L to -Lp to limit eligible points
 */
 
#include "gmt.h"

#define GMTSELECT_N_TESTS	6				/* Number of specific tests available */
#define GMTSELECT_N_CLASSES	(GMT_MAX_GSHHS_LEVEL + 1)	/* Number of bands separated by the levels */

struct GMTSELECT_DATA {	/* Used for temporary storage when sorting data on x coordinate */
	double   x;
	double   y;
	double   d;
};

struct GMTSELECT_CTRL {	/* All control options for this program (except common args) */
	/* active is TRUE if the option has been activated */
	struct A {	/* -A<min_area>[/<min_level>/<max_level>] */
		GMT_LONG active;
		struct GMT_SHORE_SELECT info;
	} A;
	struct C {	/* [-C[f]<dist>/<ptfile>] */
		GMT_LONG active;
		GMT_LONG fast;	/* FALSE for normal distance calculation, TRUE for flat-earth fast calculation */
		double dist;	/* Radius of influence for each point */
		char *file;	/* Name of file with points */
	} C;
	struct D {	/* -D<resolution> */
		GMT_LONG active;
		GMT_LONG force;	/* if TRUE, select next highest level if current set is not avaialble */
		char set;	/* One of f, h, i, l, c */
	} D;
	struct L {	/* -L[p][<dist>/<lfile>] */
		GMT_LONG active;
		GMT_LONG mode;	/* Controls what happens beyond segment endpoints */
		double dist;	/* Distance of influence for each line */
		char *file;	/* Name of file with lines */
	} L;
	struct F {	/* -F<polygon> */
		GMT_LONG active;
		char *file;	/* Name of file with polygons */
	} F;
	struct I {	/* -Icflrsz */
		GMT_LONG active;
		GMT_LONG pass[GMTSELECT_N_TESTS];	/* One flag for each setting */
	} I;
	struct N {	/* -N<maskvalues>[o] */
		GMT_LONG active;
		GMT_LONG edge;	/* TRUE if edges are considere outside */
		GMT_LONG mask[GMTSELECT_N_CLASSES];	/* Mask for each level */
	} N;
	struct Z {	/* -Z<min>/<max> */
		GMT_LONG active;
		double min;	/* Smallest z-value to pass through */
		double max;	/* Largest z-value to pass through */
	} Z;
};

int main (int argc, char **argv)
{
	GMT_LONG i, j, k, fno, n_files = 0, n_args, err, n_minimum = 2;
	GMT_LONG n_fields, n_expected_fields, ind, bin, last_bin = -1, n_output = -1;
	GMT_LONG np[2], base = 3, wd[2], id, this_node, side, is_inside, row, col, pos;
	GMT_LONG n_read = 0, n_pass = 0, no_resample = 0;

	GMT_LONG error = FALSE, dry_wet_only = FALSE, long_verbose = FALSE, need_header, shuffle;
	GMT_LONG cartesian = FALSE, nofile = TRUE, done = FALSE, first = TRUE, inside, plain_xy;
	GMT_LONG greenwich = FALSE, output_header = FALSE, do_project = FALSE, just_copy_record = FALSE;

	double west = 0.0, east = 0.0, south = 0.0, north = 0.0, xx, yy;
	double *in, west_border = 0.0, east_border = 0.0;
	double xmin, xmax, ymin, ymax, lon, step = 0.0;

	char buffer[BUFSIZ], ptr[BUFSIZ], za[GMT_TEXT_LEN], zb[GMT_TEXT_LEN], *not_used = NULL;
	char *shore_resolution[5] = {"full", "high", "intermediate", "low", "crude"};

	FILE *fp = NULL;

	struct GMT_TABLE *pol = NULL, *line = NULL, *point = NULL;
	struct GMT_GSHHS_POL *p[2] = {NULL, NULL};
	struct GMT_SHORE c;
	struct GMTSELECT_CTRL *Ctrl = NULL;

	PFL near_a_line, near_a_point;

	int compare_x (const void *point_1, const void *point_2);
	void *New_gmtselect_Ctrl (), Free_gmtselect_Ctrl (struct GMTSELECT_CTRL *C);

	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct GMTSELECT_CTRL *) New_gmtselect_Ctrl ();		/* Allocate and initialize defaults in a new control structure */

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {

				/* Common parameters */

				case 'V':
					if (argv[i][2] == 'l') long_verbose = TRUE;
				case 'H':
				case 'M':
				case 'J':
				case 'R':
				case ':':
				case 'b':
				case 'f':
				case 'm':
				case '\0':
					error += GMT_parse_common_options (argv[i], &west, &east, &south, &north);
					break;

				/* Supplemental parameters */

				case 'A':
					Ctrl->A.active = TRUE;
					GMT_set_levels (&argv[i][2], &Ctrl->A.info);
					break;
				case 'C':
					Ctrl->C.active = TRUE;
					k = 2;
					if (argv[i][2] == 'f') Ctrl->C.fast = TRUE, k = 3;
					for (j = k; argv[i][j] && argv[i][j] != '/'; j++);
					if (!argv[i][j]) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -C:  Expects -C[f]dist/file\n", GMT_program);
						error++;
					}
					else {
						Ctrl->C.file = strdup (&argv[i][j+1]);
						Ctrl->C.dist = atof (&argv[i][k]);
					}
					break;
				case 'D':
					Ctrl->D.active = TRUE;
					Ctrl->D.set = argv[i][2];
					if (argv[i][3] == '+') Ctrl->D.force = TRUE;
					break;
				case 'L':
					if (argv[i][2]) {	/* Set line options */
						Ctrl->L.active = TRUE;
						k = 2;
						if (argv[i][k] == 'p') {	/* Disallow points beyond endpoints */
							Ctrl->L.mode = 10;
							k++;
						}
						for (j = k; argv[i][j] && argv[i][j] != '/'; j++);
						if (!argv[i][j]) {
							fprintf (stderr, "%s: GMT SYNTAX ERROR -L:  Expects -L[p]dist/file\n", GMT_program);
							error++;
						}
						else {
							Ctrl->L.file = strdup (&argv[i][j+1]);
							Ctrl->L.dist = atof (&argv[i][k]);
						}
					}
					else {	/* Obsolete flag to process geographic data (use -f instead)  */
						GMT_io.in_col_type[0] =  GMT_io.out_col_type[0] = GMT_IS_LON;
						GMT_io.in_col_type[1] =  GMT_io.out_col_type[1] = GMT_IS_LAT;
					}
					break;
				case 'F':
					Ctrl->F.active = TRUE;
					Ctrl->F.file = strdup (&argv[i][2]);
					break;
				case 'I':
					Ctrl->I.active = TRUE;
					for (j = 2; argv[i][j]; j++) {
						switch (argv[i][j]) {
							case 'r':
								Ctrl->I.pass[0] = FALSE;
								break;
							case 'c':
								Ctrl->I.pass[1] = FALSE;
								break;
							case 'l':
								Ctrl->I.pass[2] = FALSE;
								break;
							case 'f':
								Ctrl->I.pass[3] = FALSE;
								break;
							case 's':
								Ctrl->I.pass[4] = FALSE;
								break;
							case 'z':
								Ctrl->I.pass[5] = FALSE;
								break;
							default:
								fprintf (stderr, "%s: GMT SYNTAX ERROR -I:  Expects -Icflrsz\n", GMT_program);
								error++;
								break;
						}
					}
					break;
				case 'N':
					Ctrl->N.active = TRUE;
					strcpy (buffer, &argv[i][2]);
					if (buffer[strlen(buffer)-1] == 'o') { /* Edge is considered outside */
						Ctrl->N.edge = TRUE;
						buffer[strlen(buffer)-1] = 0;
					}
					j = pos = 0;
					while (j < GMTSELECT_N_CLASSES && (GMT_strtok (buffer, "/", &pos, ptr))) {
						switch (ptr[0]) {
							case 's':	/* Skip points in this level */
								Ctrl->N.mask[j] = 0;
								break;
							case 'k':
								Ctrl->N.mask[j] = 1;
								break;
							default:
								fprintf (stderr, "%s: GMT SYNTAX ERROR -N option:  Bad modifier (use s or k)\n", GMT_program);
								error++;
						}
						j++;
					}
					if (!(j == 2 || j == GMTSELECT_N_CLASSES)) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -N option:  Specify 2 or 5 arguments\n", GMT_program);
						exit (EXIT_FAILURE);
					}
					dry_wet_only = (j == 2);
					break;
				case 'Z':
					Ctrl->Z.active = TRUE;
					j = sscanf (&argv[i][2], "%[^/]/%s", za, zb);
					if (j != 2) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -Z option:  Specify z_min and z_max\n", GMT_program);
						exit (EXIT_FAILURE);
					}
					if (!(za[0] == '-' && za[1] == '\0')) error += GMT_verify_expectations (GMT_io.in_col_type[GMT_Z], GMT_scanf_arg (za, GMT_io.in_col_type[GMT_Z], &Ctrl->Z.min), za);
					if (!(zb[0] == '-' && zb[1] == '\0')) error += GMT_verify_expectations (GMT_io.in_col_type[GMT_Z], GMT_scanf_arg (zb, GMT_io.in_col_type[GMT_Z], &Ctrl->Z.max), zb);
					
					break;
				case '+':	/* Undocumented option to increase path-fix resolution */
					step = atof (&argv[i][2]);
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
		fprintf (stderr, "gmtselect %s - Select data subsets based on multiple spatial criteria\n\n", GMT_VERSION);
		fprintf (stderr, "usage: gmtselect <infiles> [%s] [-C[f]<dist>/<ptfile>]\n", GMT_A_OPT);
		fprintf (stderr, "\t[-D<resolution>][+] [-F<polygon>] [%s] [-L[p]<dist>/<lfile>] [%s] [-I[cflrsz]\n", GMT_J_OPT, GMT_H_OPT);
		fprintf (stderr, "\t[-N<maskvalues>[o]] [%s] [-V[l]]\n\t[%s] [-Z<min>/<max>] [-%s] [%s] [%s]\n\n", GMT_Rgeo_OPT, GMT_t_OPT, GMT_b_OPT, GMT_f_OPT, GMT_m_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\tinfiles (in ASCII, binary, netCDF) have 2 or more columns with (x,y) or (y,x) in first columns.\n");
		fprintf (stderr, "\t  If no file(s) is given, standard input is read.\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		GMT_GSHHS_syntax ('A', "Place limits on coastline features from the GSHHS data base (ignored  unless -N is set).");
		fprintf (stderr, "\t-C pass locations that are within <dist> of any point in the ASCII <ptfile>\n");
		fprintf (stderr, "\t   Give distance as 0 if 3rd column of <ptfile> has individual distances.\n");
		fprintf (stderr, "\t   Distances are Cartesian and in user units [or spherical in km if -fg is used].\n");
		fprintf (stderr, "\t   Use -Cf for approximate (flat earth) rather than exact geodesic distances.\n");
		fprintf (stderr, "\t   If current ELLIPSOID is spherical then great circle distances replace geodesics.\n");
		fprintf (stderr, "\t   Use -R -J to compute mapped Cartesian distances in cm, inch, m, or points [%s]\n", GMT_unit_names[gmtdefs.measure_unit]);
		fprintf (stderr, "\t-D Choose one of the following resolutions: (Ignored unless -N is set)\n");
		fprintf (stderr, "\t   f - full resolution (may be very slow for large regions)\n");
		fprintf (stderr, "\t   h - high resolution (may be slow for large regions)\n");
		fprintf (stderr, "\t   i - intermediate resolution\n");
		fprintf (stderr, "\t   l - low resolution [Default]\n");
		fprintf (stderr, "\t   c - crude resolution, for tasks that need crude continent outlines only\n");
		fprintf (stderr, "\t   Append + to use a lower resolution should the chosen one not be available [abort].\n");
		fprintf (stderr, "\t-L Pass locations that are within <dist> of any line in ASCII <linefile>\n");
		fprintf (stderr, "\t   Give distance as 0 if 2nd column of segment headers have individual distances.\n");
		fprintf (stderr, "\t   Distances are Cartesian and in user units [or spherical in km if -fg is used].\n");
		fprintf (stderr, "\t   If current ELLIPSOID is spherical then great circle distances replace geodesics.\n");
		fprintf (stderr, "\t   Use -R -J to compute mapped Cartesian distances in cm, inch, m, or points [%s]\n", GMT_unit_names[gmtdefs.measure_unit]);
		fprintf (stderr, "\t   Optionally, use -Lp to exclude points projecting beyond a line's endpoints.\n");
		fprintf (stderr, "\t-F pass locations that are inside the polygons in the ASCII <polygon> file\n");
		GMT_explain_option ('H');
		fprintf (stderr, "\t-I Used to reverse the tests, i.e. pass locations outside the region\n");
		fprintf (stderr, "\t   Supply a combination of cflrz where each flag means:\n");
		fprintf (stderr, "\t   c will pass locations beyond the minimum distance to the points in -C\n");
		fprintf (stderr, "\t   f will pass locations outside the polygons in -F\n");
		fprintf (stderr, "\t   l will pass locations beyond the minimum distance to the lines in -L\n");
		fprintf (stderr, "\t   r will pass locations outside the region given in -R [and -J]\n");
		fprintf (stderr, "\t   s will pass locations that otherwise would be skipped in -N\n");
		fprintf (stderr, "\t   z will pass locations outside the range given in -Z\n");
		GMT_explain_option ('J');
		fprintf (stderr, "\t-N set if a point outside or inside a geographic feature should be s(kipped) or k(ept).\n");
		fprintf (stderr, "\t   Append o to let feature boundary be considered outside [Default is inside].\n");
		fprintf (stderr, "\t   Specify this information with s or k using 1 of 2 formats:\n");
		fprintf (stderr, "\t   -N<wet>/<dry>.\n");
		fprintf (stderr, "\t   -N<ocean>/<land>/<lake>/<island>/<pond>.\n");
		fprintf (stderr, "\t   k means keep and s means skip [Default is s/k/s/k/s (i.e., s/k)]\n");
		GMT_explain_option ('R');
		GMT_explain_option ('V');
		fprintf (stderr, "\t   Append l for long verbose, reporting every 1000 points.\n");
		fprintf (stderr, "\t-Z assumes the 3rd data column contains z-values and we want to keep records with\n");
		fprintf (stderr, "\t   <min> <= z <= <max>.  Use - for <min> or <max> if there is no lower/upper limit.\n");
		GMT_explain_option (':');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf (stderr, "\t   Default is 2 input columns (3 if -Z is used)\n");
		GMT_explain_option ('o');
		GMT_explain_option ('n');
		GMT_explain_option ('f');
		GMT_explain_option ('m');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	if (!Ctrl->N.active && (Ctrl->A.active || Ctrl->D.active)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  -A and -D requires -N!\n", GMT_program);
		error++;
	}
	if (Ctrl->Z.active && Ctrl->Z.max <= Ctrl->Z.min) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  -Z must have zmax > zmin!\n", GMT_program);
		error++;
	}
	if ((GMT_io.binary[GMT_IN] || GMT_io.netcdf[GMT_IN]) && GMT_io.io_header[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary and netCDF input data cannot have header -H\n", GMT_program);
		error++;
	}
        if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = 2 + Ctrl->Z.active;
        if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] < (2 + Ctrl->Z.active)) {
                fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data (-bi) must have at least %ld columns\n", GMT_program, 2 + Ctrl->Z.active);
		error++;
	}
	if (error) exit (EXIT_FAILURE);

	if (Ctrl->C.active && !GMT_IS_MAPPING && !(GMT_io.in_col_type[0] == GMT_IS_LON && GMT_io.in_col_type[1] == GMT_IS_LAT)) cartesian = TRUE;

	if (GMT_io.binary[GMT_IN] && gmtdefs.verbose) {
		char *type[2] = {"double", "single"};
		fprintf (stderr, "%s: Expects %ld-column %s-precision binary data\n", GMT_program, GMT_io.ncol[GMT_IN], type[GMT_io.single_precision[GMT_IN]]);
	}

#ifdef SET_IO_MODE
	GMT_setmode (GMT_OUT);
#endif

	shuffle = (gmtdefs.xy_toggle[GMT_IN] != gmtdefs.xy_toggle[GMT_OUT]);	/* Must rewrite output record */
	n_minimum = (Ctrl->Z.active) ? 3 : 2;	/* Minimum number of columns in ASCII input */
	
	if (!project_info.region_supplied && Ctrl->N.active) {	/* If we use coastline data or used -fg but didnt give -R we imply -Rg */
		project_info.region_supplied = TRUE;
		west = 0.0;	east = 360.0;	south = -90.0;	north = +90.0;
		GMT_io.in_col_type[0] = GMT_IS_LON;
		GMT_io.in_col_type[1] = GMT_IS_LAT;
	}
	if (project_info.region_supplied) {	/* -R was set directly or indirectly; hence must set -J if not supplied */
		if (project_info.projection == GMT_NO_PROJ) {	/* -J not specified, set one implicitly */
			if (GMT_io.in_col_type[0] == GMT_IS_LON) {	/* We know we have geographic data (either via -f or -N) */
				project_info.degree[0] = project_info.degree[1] = no_resample = TRUE;
			}
			/* Supply dummy linear proj */
			project_info.projection = project_info.xyz_projection[0] = project_info.xyz_projection[1] = GMT_LINEAR;
			project_info.pars[0] = project_info.pars[1] = 1.0;
		}
		else
			do_project = TRUE;	/* Only TRUE when the USER selected -J, not when we supply dummy -Jx1d */
		step = 0.01;
		if (GMT_io.in_col_type[0] == GMT_IS_LON) {
			if (west < 0.0 && east < 0.0) {
				west += 360.0;
				east += 360.0;
			}
			greenwich = (west < 0.0 && east > 0.0);
		}
		GMT_err_fail (GMT_map_setup (west, east, south, north), "");
		if (no_resample) GMT_parallel_straight = GMT_meridian_straight = 2;	/* No resampling along bin boundaries */
	}

	if (do_project && gmtdefs.verbose) fprintf (stderr, "%s: Warning: -J means all data will be projected before tests are applied\n", GMT_program);
	 
	if (Ctrl->N.active) {
		if (Ctrl->D.force) Ctrl->D.set = GMT_shore_adjust_res (Ctrl->D.set);
		if (Ctrl->D.active) base = GMT_set_resolution (&Ctrl->D.set, 'D');
		if (dry_wet_only) {	/* Post-process -N choice */
			Ctrl->N.mask[3] = Ctrl->N.mask[1];
			Ctrl->N.mask[2] = Ctrl->N.mask[4] = Ctrl->N.mask[0];
		}
		if (GMT_init_shore (Ctrl->D.set, &c, west, east, south, north, &Ctrl->A.info)) {
			fprintf (stderr, "%s: %s resolution shoreline data base not installed\n", GMT_program, shore_resolution[base]);
			exit (EXIT_FAILURE);
		}
		if (gmtdefs.verbose == 2) fprintf (stderr, "GSHHS version %s\n%s\n%s\n", c.version, c.title, c.source);
		west_border = floor (project_info.w / c.bsize) * c.bsize;
		east_border = ceil (project_info.e / c.bsize) * c.bsize;
		wd[0] = 1;	wd[1] = -1;
		np[0] = np[1] = 0;
	}

	just_copy_record = !(GMT_io.binary[GMT_IN] || GMT_io.netcdf[GMT_IN] || GMT_io.binary[GMT_OUT] || GMT_io.netcdf[GMT_OUT] || shuffle);

	/* Initiate pointer to distance calculation function */
	if (GMT_io.in_col_type[0] & GMT_IS_GEO && !do_project) {	/* Geographic data and no -R -J conversion */
		if (Ctrl->C.fast)
			GMT_distance_func = GMT_flatearth_dist_km;
		else
			GMT_distance_func = (PFD) ((GMT_IS_SPHERICAL) ? GMT_geodesic_dist_km : GMT_great_circle_dist_km);
		near_a_line  = (PFL) GMT_near_a_line_spherical;
		near_a_point = (PFL) GMT_near_a_point_spherical;
	}
	else {	/* Cartesian data (or lon/lat projected via -R -J) */
		GMT_distance_func = (PFD) GMT_cartesian_dist;
		near_a_line  = (PFL) GMT_near_a_line_cartesian;
		near_a_point = (PFL) GMT_near_a_point_cartesian;
	}

	if (Ctrl->C.active) { 	/* Initialize point structure used in test for proximity to points */
		GMT_import_table ((void *)Ctrl->C.file, GMT_IS_FILE, &point, Ctrl->C.dist, greenwich, FALSE, FALSE);
		if (point->segment[0]->n_columns < 2) {	/* Trouble */
			fprintf (stderr, "%s: GMT SYNTAX ERROR -C:  %s does not have at least 2 columns with coordinates\n", GMT_program, Ctrl->C.file);
			exit (EXIT_FAILURE);
		}
		if (Ctrl->C.dist == 0.0 && point->segment[0]->n_columns <= 2) {	/* Trouble */
			fprintf (stderr, "%s: GMT SYNTAX ERROR -C:  %s does not have a 3rd column with distances yet -C0/<file> was given\n", GMT_program, Ctrl->C.file);
			exit (EXIT_FAILURE);
		}
		if (do_project) {	/* Convert all the points using the map projection */
			for (i = 0; i < point->n_segments; i++) {
				for (j = 0; j < point->segment[i]->n_rows; j++) {
					GMT_geo_to_xy (point->segment[i]->coord[GMT_X][j], point->segment[i]->coord[GMT_Y][j], &xx, &yy);
					point->segment[i]->coord[GMT_X][j] = xx;
					point->segment[i]->coord[GMT_Y][j] = yy;
				}
			}
			cartesian = TRUE;	/* Well, now it is */
		}
		if (cartesian) {	/* Speed up testing by sorting points on the x-coordinate first */
			struct GMTSELECT_DATA *data;	/* Used for temporary storage when sorting data on x coordinate */

			/* Copy xp into struct data, sort, and copy back */

			data = (struct GMTSELECT_DATA *) GMT_memory (VNULL, (size_t)point->n_records, sizeof (struct GMTSELECT_DATA), GMT_program);

			for (i = k = 0; i < point->n_segments; i++) {
				for (j = 0; j < point->segment[i]->n_rows; j++, k++) {
					data[k].x = point->segment[i]->coord[GMT_X][j];
					data[k].y = point->segment[i]->coord[GMT_Y][j];
					data[k].d = (Ctrl->C.dist == 0.0) ? point->segment[i]->coord[GMT_Z][j] : Ctrl->C.dist;
				}
			}
			
			/* Sort on x to speed up inside testing */
			qsort ((void *)data, (size_t)point->n_records, sizeof (struct GMTSELECT_DATA), compare_x);
			
			for (i = k = 0; i < point->n_segments; i++) {	/* Put back the new order */
				for (j = 0; j < point->segment[i]->n_rows; j++, k++) {
					point->segment[i]->coord[GMT_X][j] = data[k].x;
					point->segment[i]->coord[GMT_Y][j] = data[k].y;
					if (Ctrl->C.dist == 0.0) point->segment[i]->coord[GMT_Z][j] = data[k].d ;
				}
			}
			GMT_free ((void *)data);
		}
	}

	if (Ctrl->L.active) {	/* Initialize lines structure used in test for proximity to lines */
		GMT_import_table ((void *)Ctrl->L.file, GMT_IS_FILE, &line, Ctrl->L.dist, greenwich, FALSE, FALSE);
		if (line->segment[0]->n_columns < 2) {	/* Trouble */
			fprintf (stderr, "%s: GMT SYNTAX ERROR -L:  %s does not have at least 2 columns with coordinates\n", GMT_program, Ctrl->L.file);
			exit (EXIT_FAILURE);
		}
		if (do_project) {	/* Convert all the line points using the map projection */
			for (i = 0; i < line->n_segments; i++) {
				for (j = 0; j < line->segment[i]->n_rows; j++) {
					GMT_geo_to_xy (line->segment[i]->coord[GMT_X][j], line->segment[i]->coord[GMT_Y][j], &xx, &yy);
					line->segment[i]->coord[GMT_X][j] = xx;
					line->segment[i]->coord[GMT_Y][j] = yy;
				}
			}
		}
	}
	if (Ctrl->F.active) {	/* Initialize polygon structure used in test for polygon in/out test */
		GMT_io.skip_duplicates = TRUE;
		GMT_import_table ((void *)Ctrl->F.file, GMT_IS_FILE, &pol, 0.0, greenwich, TRUE, FALSE);
		if (pol->segment[0]->n_columns < 2) {	/* Trouble */
			fprintf (stderr, "%s: GMT SYNTAX ERROR -F:  %s does not have at least 2 columns with coordinates\n", GMT_program, Ctrl->F.file);
			exit (EXIT_FAILURE);
		}
		if (do_project) {	/* Convert all the polygons points using the map projection */
			for (i = 0; i < pol->n_segments; i++) {
				for (j = 0; j < pol->segment[i]->n_rows; j++) {
					GMT_geo_to_xy (pol->segment[i]->coord[GMT_X][j], pol->segment[i]->coord[GMT_Y][j], &xx, &yy);
					pol->segment[i]->coord[GMT_X][j] = xx;
					pol->segment[i]->coord[GMT_Y][j] = yy;
				}
			}
		}
		GMT_io.skip_duplicates = FALSE;	/* Reset to FALSE */
	}
	plain_xy = (do_project || GMT_io.in_col_type[GMT_X] != GMT_IS_LON);
	
	need_header = GMT_io.multi_segments[GMT_OUT];	/* Only need to break up segments */
	is_inside = (Ctrl->N.edge) ? 2 : 1;		/* Whether of not being exactly on an edge is outside */
	/* Now we are ready to take on some input values */

	if (n_files > 0)
		nofile = FALSE;
	else
		n_files = 1;
	n_args = (argc > 1) ? argc : 2;

	for (fno = 1; !done && fno < n_args; fno++) {	/* Loop over input files, if any */
		if (!nofile && argv[fno][0] == '-') continue;

		if (nofile) {	/* Just read standard input */
			fp = GMT_stdin;
			done = TRUE;
			if (gmtdefs.verbose) fprintf (stderr, "%s: Reading from standard input\n", GMT_program);
#ifdef SET_IO_MODE
			GMT_setmode (GMT_IN);
#endif
		}
		else if ((fp = GMT_fopen (argv[fno], GMT_io.r_mode)) == NULL) {
			fprintf (stderr, "%s: Cannot open file %s\n", GMT_program, argv[fno]);
			continue;
		}

		if (!nofile && gmtdefs.verbose) fprintf (stderr, "%s: Working on file %s\n", GMT_program, argv[fno]);

		if (GMT_io.io_header[GMT_IN]) {
			for (i = 0; i < GMT_io.n_header_recs; i++) {
				not_used = GMT_fgets (buffer, BUFSIZ, fp);
				if (first && GMT_io.io_header[GMT_OUT]) GMT_fputs (buffer, GMT_stdout);
			}
			first = FALSE;
		}
		n_expected_fields = (GMT_io.ncol[GMT_IN]) ? GMT_io.ncol[GMT_IN] : GMT_MAX_COLUMNS;

		while ((n_fields = GMT_input (fp, &n_expected_fields, &in)) >= 0 && ! (GMT_io.status & GMT_IO_EOF)) {	/* Not yet EOF */

			while (GMT_io.status & GMT_IO_SEGMENT_HEADER && !(GMT_io.status & GMT_IO_EOF)) {
				output_header = TRUE;
				n_fields = GMT_input (fp, &n_expected_fields, &in);
			}
			if ((GMT_io.status & GMT_IO_EOF)) continue;	/* At EOF */

			n_read++;

			if (!just_copy_record && (GMT_io.status & GMT_IO_MISMATCH)) {
				fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld (skipped)\n", GMT_program, n_fields,  n_expected_fields, n_read);
				continue;
			}

			if (long_verbose && n_read%1000 == 0) fprintf (stderr, "%s: Read %ld records, passed %ld records\r", GMT_program, n_read, n_pass);

			if (n_fields < n_minimum) {
				if (Ctrl->Z.active)
					fprintf (stderr, "%s: -Z requires a data file with at least 3 columns; this file only has %ld near line %ld (skipped)\n", GMT_program, n_fields, n_read);
				else
					fprintf (stderr, "%s: Data file must have at least 2 columns; this file only has %ld near line %ld (skipped)\n", GMT_program, n_fields, n_read);
				continue;
			}

			if (Ctrl->Z.active) {
				if (GMT_is_dnan (in[GMT_Z])) { output_header = need_header; continue;}	/* cannot keep when no z */
				inside = (in[GMT_Z] >= Ctrl->Z.min && in[GMT_Z] <= Ctrl->Z.max); 
				if (inside != Ctrl->I.pass[5]) { output_header = need_header; continue;}
			}

			lon = in[GMT_X];
			if (project_info.region_supplied) {
				inside = !GMT_map_outside (lon, in[GMT_Y]);
				if (inside != Ctrl->I.pass[0]) { output_header = need_header; continue;}
			}

			if (do_project)
				GMT_geo_to_xy (lon, in[GMT_Y], &xx, &yy);
			else {
				xx = lon;
				yy = in[GMT_Y];
			}
			
			if (Ctrl->C.active) {	/* Check for distance to points */
				inside = near_a_point (xx, yy, point, Ctrl->C.dist); 
				if (inside != Ctrl->I.pass[1]) { output_header = need_header; continue;}
			}

			if (Ctrl->L.active) {
				inside = near_a_line (xx, yy, line, Ctrl->L.mode, NULL, NULL, NULL);
				if (inside != Ctrl->I.pass[2]) { output_header = need_header; continue;}
			}

			if (Ctrl->F.active) {
				inside = FALSE;
				for (i = 0; i < pol->n_segments && !inside; i++) {	/* Check each polygon until we find that our point is inside */
					if (plain_xy) {	/* Non-geographic (or projected lon/lat) */
						if (yy < pol->segment[i]->min[GMT_Y] || yy > pol->segment[i]->max[GMT_Y]) continue;	/* Outside polygon's y-range */
						if (xx < pol->segment[i]->min[GMT_X] || xx > pol->segment[i]->max[GMT_X]) continue;	/* Outside polygon's x-range */
						/* Here we must do the full Cartesian polygon check */
						inside = (GMT_non_zero_winding (xx, yy, pol->segment[i]->coord[GMT_X], pol->segment[i]->coord[GMT_Y], pol->segment[i]->n_rows) >= is_inside);
					}
					else {	/* Lon-lat data */
						if (pol->segment[i]->pole) {	/* Special testing for polar caps */
							if (pol->segment[i]->pole == +1 && yy < pol->segment[i]->min[GMT_Y]) continue;	/* Below a N-polar cap */
							if (pol->segment[i]->pole == -1 && yy > pol->segment[i]->max[GMT_Y]) continue;	/* Above a S-polar cap */
						}
						else {
							if (yy < pol->segment[i]->min[GMT_Y] || yy > pol->segment[i]->max[GMT_Y]) continue;	/* Outside polygon's y-range */
							xx = lon - 360.0;
							while (xx < pol->segment[i]->min[GMT_X]) xx += 360.0;	/* Wind to the east of the west boundary */
							if (xx > pol->segment[i]->max[GMT_X]) continue;		/* Outside polygon's longitude-range */
						}
						/* Here we must do the full spherical polygon check */
						inside = (GMT_inonout_sphpol (xx, yy, pol->segment[i]) >= is_inside);
					}
				}
				if (inside != Ctrl->I.pass[3]) { output_header = need_header; continue;}
			}

			if (Ctrl->N.active) {
				xx = lon;
				while (xx < 0.0) xx += 360.0;
				row = ((int)floor ((90.0 - in[GMT_Y]) / c.bsize));
				if (row >= c.bin_ny) row = c.bin_ny - 1;	/* Presumably only kicks in for south pole */
				col = (int)floor (xx / c.bsize);
				bin = row * c.bin_nx + col;
				if (bin != last_bin) {	/* Do this upon entering new bin */
					ind = 0;
					while (ind < c.nb && c.bins[ind] != bin) ind++;	/* Set ind to right bin */
					if (ind == c.nb) continue;			/* Bin not among the chosen ones */
					last_bin = bin;
					GMT_free_shore (&c);	/* Free previously allocated arrays */
					if ((err = GMT_get_shore_bin (ind, &c))) {
						fprintf (stderr, "%s: %s [%s resolution shoreline]\n", GMT_program, GMT_strerror(err), shore_resolution[base]);
						exit (EXIT_FAILURE);
					}

					/* Must use polygons.  Go in both directions to cover both land and sea */
					for (id = 0; id < 2; id++) {
						GMT_free_polygons (p[id], np[id]);
						if (np[id]) GMT_free ((void *)p[id]);
						np[id] = GMT_assemble_shore (&c, wd[id], TRUE, greenwich, west_border, east_border, &p[id]);
						np[id] = GMT_prep_polygons (&p[id], np[id], !no_resample, step, -1);
					}
				}

				if (c.ns == 0) {	/* No lines go through, check node level */
					this_node = MIN (MIN (c.node_level[0], c.node_level[1]) , MIN (c.node_level[2], c.node_level[3]));
				}
				else {
					this_node = 0;
					GMT_geo_to_xy (lon, in[GMT_Y], &xx, &yy);
					for (id = 0; id < 2; id++) {

						for (k = 0; k < np[id]; k++) {

							if (p[id][k].n == 0) continue;

							/* Find min/max of polygon */

							xmin = xmax = p[id][k].lon[0];
							ymin = ymax = p[id][k].lat[0];

							for (i = 1; i < p[id][k].n; i++) {
								if (p[id][k].lon[i] < xmin) xmin = p[id][k].lon[i];
								if (p[id][k].lon[i] > xmax) xmax = p[id][k].lon[i];
								if (p[id][k].lat[i] < ymin) ymin = p[id][k].lat[i];
								if (p[id][k].lat[i] > ymax) ymax = p[id][k].lat[i];
							}

							if (yy < ymin || yy > ymax) continue;
							if (xx < xmin || xx > xmax) continue;

							/* Must compare with polygon */
							
							if ((side = GMT_non_zero_winding (xx, yy, p[id][k].lon, p[id][k].lat, p[id][k].n)) < is_inside) continue;	/* Outside polygon */

							/* Here, point is inside, we must assign value */

							if (p[id][k].level > this_node) this_node = p[id][k].level;
						}
					}
				}
				inside = Ctrl->N.mask[this_node];
				if (inside != Ctrl->I.pass[4]) { output_header = need_header; continue;}
			}

			/* Here, we have passed all test and the point is output */

			if (output_header) {	/* First output point for this segment - write the header */
				GMT_write_segmentheader (GMT_stdout, n_expected_fields);
				output_header = FALSE;
			}

			if (just_copy_record)
				GMT_fputs (GMT_io.current_record, GMT_stdout);
			else {
				if (n_output < 0) n_output = n_expected_fields;
				GMT_output (GMT_stdout, n_output, in);
			}
			n_pass++;
		}
		if (fp != GMT_stdin) GMT_fclose (fp);
	}

	if (gmtdefs.verbose) fprintf (stderr, "%s: Read %ld records, passed %ld records\n", GMT_program, n_read, n_pass);

	if (Ctrl->C.active) GMT_free_table (point);
	if (Ctrl->L.active) GMT_free_table (line);
	if (Ctrl->F.active) GMT_free_table (pol);
	if (Ctrl->N.active) {
		GMT_free_shore (&c);
		GMT_shore_cleanup (&c);
		for (id = 0; id < 2; id++) {
			GMT_free_polygons (p[id], np[id]);
			if (np[id]) GMT_free ((void *)p[id]);
		}
	}
	
	Free_gmtselect_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

int compare_x (const void *point_1, const void *point_2)
{
	struct GMTSELECT_DATA *p1, *p2;

	p1 = (struct GMTSELECT_DATA *)point_1;
	p2 = (struct GMTSELECT_DATA *)point_2;

	if (p1->x < p2->x)
		return (-1);
	else if (p1->x > p2->x)
		return (1);
	else
		return (0);
}

void *New_gmtselect_Ctrl () {	/* Allocate and initialize a new control structure */
	GMT_LONG i;
	struct GMTSELECT_CTRL *C;
	
	C = (struct GMTSELECT_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GMTSELECT_CTRL), "New_gmtselect_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	
	C->A.info.high = GMT_MAX_GSHHS_LEVEL;				/* Include all GSHHS levels */
	C->D.set = 'l';							/* Low-resolution coastline data */
	for (i = 0; i < GMTSELECT_N_TESTS; i++) C->I.pass[i] = TRUE;	/* Default is to pass if we are inside */
	memset ((void *)C->N.mask, 0, (size_t)(GMTSELECT_N_CLASSES * sizeof (GMT_LONG)));	/* Default for "wet" areas = 0 (outside) */
	C->N.mask[1] = C->N.mask[3] = 1;				/* Default for "dry" areas = 1 (inside) */
	C->Z.min = -DBL_MAX;	C->Z.max = DBL_MAX;			/* No limits on z-range */
	
	return ((void *)C);
}

void Free_gmtselect_Ctrl (struct GMTSELECT_CTRL *C) {	/* Deallocate control structure */
	if (C->C.file) free ((void *)C->C.file);	
	if (C->F.file) free ((void *)C->F.file);	
	if (C->L.file) free ((void *)C->L.file);	
	GMT_free ((void *)C);	
}
