/*-----------------------------------------------------------------
 *	$Id: x2sys_binlist.c 10173 2014-01-01 09:52:34Z pwessel $
 *
 *      Copyright (c) 1999-2014 by P. Wessel
 *      See LICENSE.TXT file for copying and redistribution conditions.
 *
 *      This program is free software; you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation; version 2 or any later version.
 *
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *
 *      Contact info: www.soest.hawaii.edu/pwessel
 *--------------------------------------------------------------------*/
/* x2sys_binlist will read one or several data files and dump their
 * contents to stdout in ascii or binary (double precision) mode.
 * Input data file formats are determined by the definition file
 * given by the -D option.
 *
 * Author:	Paul Wessel
 * Date:	15-JUN-2004
 * Version:	1.1, based on the spirit of the old xsystem code
 *
 */

#include "x2sys.h"

#define EA_LAT "37:04:17.1660757541775"

struct BINCROSS {
	double x, y, d;
};

int main (int argc, char **argv)
{
	char *TAG = CNULL, **trk_name = NULL, buffer[BUFSIZ];

	GMT_LONG i, j, k, n_tracks, k1, ij, ii, jj, nx, bi, bj, start_i, end_i, jump_180, jump_360;
	GMT_LONG this_bin_i;	/* This i node for bin */
	GMT_LONG this_bin_j;	/* This j node for bin */
	GMT_LONG this_bin_ij;/* This bin */
	GMT_LONG last_bin_i;	/* Previous i node for bin */
	GMT_LONG last_bin_j;	/* Previous j node for bin */
	GMT_LONG last_bin_ij;/* Previous bin */
	unsigned int nav_flag;
	size_t nx_alloc = GMT_SMALL_CHUNK;

	GMT_LONG error = FALSE, trackline = FALSE, gap, ea_mode = FALSE, cmdline_files;

	double **data = NULL, *dist_km = NULL, *dist_bin = NULL, dist_scale, x, y, dx, del_x, del_y, y_max = 90.0;

	struct X2SYS_INFO *s = NULL;
	struct X2SYS_FILE_INFO p;		/* File information */
	struct X2SYS_BIX B;
	struct BINCROSS *X = NULL;

	int outside (double x, double y, struct X2SYS_BIX *B, GMT_LONG geo);
	unsigned int get_data_flag (double *data[], GMT_LONG j, struct X2SYS_INFO *s);
	int comp_bincross (const void *p1, const void *p2);

	argc = (int)GMT_begin (argc, argv);
	
	for (i = strlen(argv[0]); i >= 0 && argv[0][i] != '/'; i--);
	X2SYS_program = &argv[0][i+1];	/* Name without full path */

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				/* Common parameters */

				case 'V':
				case '\0':
					error += GMT_parse_common_options (argv[i], NULL, NULL, NULL, NULL);
					break;

				/* Supplemental parameters */

				case 'D':
					trackline = TRUE;
					break;
				case 'E':
					ea_mode = TRUE;
					break;
				case 'T':
					TAG = &argv[i][2];
					break;
				default:
					error = TRUE;
					break;
			}
		}
	}

	if (argc == 1 || error || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "x2sys_binlist %s - Create bin index listing from data files\n\n", X2SYS_VERSION);
		fprintf (stderr, "usage: x2sys_binlist <files> -T<TAG> [-D] [-E] [-V]\n\n");

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\t<files> is one or more datafiles, or give =<files.lis> for a file with a list of datafiles\n");
		fprintf (stderr, "\t-T <TAG> is the system tag for this compilation\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-D Calculate track-lengths per bin (see -C for method and -N for units)\n");
		fprintf (stderr, "\t-E Do binning using equal-area bins (with -D only)\n");
		GMT_explain_option ('V');
		exit (EXIT_FAILURE);
	}

	if (ea_mode && !trackline) {
		fprintf (stderr, "%s: -E requires -D\n", GMT_program);
		exit (EXIT_FAILURE);		
	}
	
	if ((n_tracks = x2sys_get_tracknames (argc, argv, &trk_name, &cmdline_files)) == 0) {
		fprintf (stderr, "%s: No datafiles given!\n", GMT_program);
		exit (EXIT_FAILURE);		
	}

	x2sys_err_fail (x2sys_set_system (TAG, &s, &B, &GMT_io), TAG);

	if (ea_mode && !s->geographic) {
		fprintf (stderr, "%s: -E requires geographic data; your TAG implies Cartesian\n", GMT_program);
		exit (EXIT_FAILURE);		
	}

	if (s->geographic) {
		GMT_io.out_col_type[0] = GMT_IS_LON;
		GMT_io.out_col_type[1] = GMT_IS_LAT;
		GMT_io.geo.range = s->geodetic;
	}
	else
		GMT_io.out_col_type[0] = GMT_io.out_col_type[1] = GMT_IS_FLOAT;
	GMT_io.out_col_type[2] = GMT_IS_FLOAT;
	
	MGD77_Set_Unit (s->unit[X2SYS_DIST_SELECTION], &dist_scale, -1);	/* Gets scale which multiplies meters to chosen distance unit */

	if (ea_mode) {
		double mid;
		char proj[80];
		/* Do the equal area map projection so W = 360 and H = 180 */
		if (!(GMT_IS_ZERO (B.x_max - B.x_min - 360.0) && GMT_IS_ZERO (B.y_max - B.y_min - 180.0))) {
			if (gmtdefs.verbose) fprintf (stderr, "%s: -E requires a global region (-Rg or -Rd)", GMT_program);
				exit (EXIT_FAILURE);
		}
		gmtdefs.ellipsoid = GMT_N_ELLIPSOIDS - 1;	/* Make sure we use a spherical projection */
		mid = 0.5 * (B.x_max + B.x_min);	/* Central longitude to use */
		if (gmtdefs.verbose) fprintf (stderr, "%s: To undo equal-area projection, use -R%g/%g/%g/%g -JY%g/%s/360i\n", GMT_program, B.x_min, B.x_max, B.y_min, B.y_max, mid, EA_LAT);
		sprintf (proj, "Y%g/%s/360", mid, EA_LAT);
		GMT_parse_J_option (proj);
		GMT_err_fail (GMT_map_setup (B.x_min, B.x_max, B.y_min, B.y_max), "");
		GMT_geo_to_xy (B.x_min, B.y_min, &B.x_min, &B.y_min);
		GMT_geo_to_xy (B.x_max, B.y_max, &B.x_max, &B.y_max);
		y_max = B.y_max;
	}
	
	x2sys_bix_init (&B, TRUE);
	nav_flag = (1 << s->x_col[GMT_IN]) + (1 << s->y_col[GMT_IN]);	/* For bins just cut by track but no points inside the bin */
	jump_180 = irint (180.0 / B.bin_x);
	jump_360 = irint (360.0 / B.bin_x);

	X = (struct BINCROSS *)GMT_memory (VNULL, nx_alloc, sizeof (struct BINCROSS), GMT_program);
	
	if (trackline) {
		dist_bin = (double *)GMT_memory (VNULL, (size_t)B.nm_bin, sizeof (double), GMT_program);
	}

	sprintf (buffer, "# %s\n", TAG);
	GMT_fputs (buffer, GMT_stdout);

	for (i = 0; i < n_tracks; i++) {

		if (gmtdefs.verbose) fprintf (stderr, "x2sys_binlist: Reading file %s ", trk_name[i]);

		x2sys_err_fail ((s->read_file) (trk_name[i], &data, s, &p, &GMT_io, &j), trk_name[i]);
		if (gmtdefs.verbose) fprintf (stderr, "[%s]\n", s->path);
		
		if (ea_mode) {	/* Project coordinates */
			for (j = 0; j < p.n_rows; j++) GMT_geo_to_xy (data[s->x_col[GMT_IN]][j], data[s->y_col[GMT_IN]][j], &data[s->x_col[GMT_IN]][j], &data[s->y_col[GMT_IN]][j]);
		}
		
		/* Reset bin flags */

		memset ((void *)B.binflag, 0, (size_t)(B.nm_bin * sizeof (unsigned int)));
		if (trackline) {
			memset ((void *)dist_bin, 0, (size_t)(B.nm_bin * sizeof (double)));
			GMT_err_fail (GMT_distances (data[s->x_col[GMT_IN]], data[s->y_col[GMT_IN]], p.n_rows, dist_scale, -s->dist_flag, &dist_km), "");	/* -ve gives increments */
		}

		last_bin_ij = last_bin_i = last_bin_j = -1;
		for (j = 0; j < p.n_rows; j++) {
			if (outside (data[s->x_col[GMT_IN]][j], data[s->y_col[GMT_IN]][j], &B, s->geographic)) continue;
			 x2sys_err_fail (x2sys_bix_get_ij (data[s->x_col[GMT_IN]][j], data[s->y_col[GMT_IN]][j], &this_bin_i, &this_bin_j, &B, &this_bin_ij), "");

			/* While this may be the same bin as the last bin, the data available may have changed so we keep
			 * turning the data flags on again and again. */
			 
			B.binflag[this_bin_ij] |= get_data_flag (data, j, s);

			if (!trackline) continue;	/* Not worried about trackline lengths */
			
			if (last_bin_ij == -1) last_bin_ij = this_bin_ij;	/* Initialize last bin to this bin the first time */
			
			if (j > 0) { /* Can check for gaps starting with 1st to 2nd point */
				gap = FALSE;
				if (s->t_col[GMT_IN] >= 0) {	/* There is a time column in the data*/
					if (GMT_is_dnan (data[s->t_col[GMT_IN]][j])&& (dist_km[j] - dist_km[j-1]) > B.dist_gap) /* but time = NaN, so test for gaps based on distance */
						gap = TRUE;
			   		else if ((data[s->t_col[GMT_IN]][j] - data[s->t_col[GMT_IN]][j-1]) > B.time_gap)	/* We have a time data gap so we skip this interval */
						gap = TRUE;
				}
				else if ((dist_km[j] - dist_km[j-1]) > B.dist_gap) /* There is no time column, must test for gaps based on distance */
					gap = TRUE;
				
				if (gap) {
					last_bin_ij = this_bin_ij;	/* Update the last point's index info */
					last_bin_i = this_bin_i;
					last_bin_j = this_bin_j;
					continue;
				}
			}

			if (this_bin_ij == last_bin_ij) {	/* Same bin, keep adding up incremental distances (this adds 0 the very first time) */
				dist_bin[this_bin_ij] += dist_km[j];
			}
			else  {	/* Crossed into another bin */
				if (s->geographic && (this_bin_i - last_bin_i) > jump_180) {		/* Jumped from east to west across Greenwich */
					start_i = this_bin_i + 1;
					end_i = last_bin_i + jump_360;
					dx = (data[s->x_col[GMT_IN]][j] - 360.0) - data[s->x_col[GMT_IN]][j-1];
				}
				else if (s->geographic && (this_bin_i - last_bin_i) < -jump_180) {	/* Jumped from west to east across Greenwich */
					start_i = last_bin_i + 1;
					end_i = this_bin_i + jump_360;
					dx = data[s->x_col[GMT_IN]][j] - (data[s->x_col[GMT_IN]][j-1] - 360.0);
				}
				else {								/* Did no such thing */
					start_i = MIN (last_bin_i, this_bin_i) + 1;
					end_i = MAX (last_bin_i, this_bin_i);
					dx = data[s->x_col[GMT_IN]][j] - data[s->x_col[GMT_IN]][j-1];
				}
				
				/* Find all the bin-line intersections */
				
				/* Add the start and stop coordinates to the xc/yc arrays so we can get mid-points of the intervals */
				
				X[0].x = data[s->x_col[GMT_IN]][j-1];	X[0].y = data[s->y_col[GMT_IN]][j-1];	X[0].d = 0.0;
				X[1].x = data[s->x_col[GMT_IN]][j];	X[1].y = data[s->y_col[GMT_IN]][j];	X[1].d = hypot (dx, data[s->y_col[GMT_IN]][j] - data[s->y_col[GMT_IN]][j-1]);
				nx = 2;
				for (bj = MIN (last_bin_j, this_bin_j) + 1; bj <= MAX (last_bin_j, this_bin_j); bj++) {	/* If we go in here we know dy is non-zero */
					y = B.y_min + bj * B.bin_y;
					del_y = y - data[s->y_col[GMT_IN]][j-1];
					del_x = del_y * dx / (data[s->y_col[GMT_IN]][j] - data[s->y_col[GMT_IN]][j-1]);
					x = data[s->x_col[GMT_IN]][j-1] + del_x;
					X[nx].x = x;	X[nx].y = y;	X[nx].d = hypot (del_x , del_y);
					nx++;
					if (nx == (int)nx_alloc) {
						nx_alloc <<= 1;
						X = (struct BINCROSS *)GMT_memory ((void *)X, nx_alloc, sizeof (struct BINCROSS), GMT_program);
					}
				}
				for (bi = start_i; bi <= end_i; bi++) {	/* If we go in here we think dx is non-zero (we do a last-ditch dx check just in case) */
					x = B.x_min + bi * B.bin_x;
					if (s->geographic && x >= 360.0) x -= 360.0;
					del_x = x - data[s->x_col[GMT_IN]][j-1];
					if (fabs (del_x) > 180.0) del_x = copysign (360.0 - fabs (del_x), -del_x);
					del_y = (dx == 0.0) ? 0.5 * (data[s->y_col[GMT_IN]][j] - data[s->y_col[GMT_IN]][j-1]) : del_x * (data[s->y_col[GMT_IN]][j] - data[s->y_col[GMT_IN]][j-1]) / dx;
					y = data[s->y_col[GMT_IN]][j-1] + del_y;
					if (s->geographic && fabs (y) > y_max) {
						y = copysign (y_max, y);
						del_y = y - data[s->y_col[GMT_IN]][j-1];
					}
					X[nx].x = x;	X[nx].y = y;	X[nx].d = hypot (del_x, del_y);
					nx++;
					if (nx == (int)nx_alloc) {
						nx_alloc <<= 1;
						X = (struct BINCROSS *)GMT_memory ((void *)X, nx_alloc, sizeof (struct BINCROSS), GMT_program);
					}
				}
				
				/* Here we have 1 or more intersections */
				
				qsort ((void *)X, (size_t)nx, sizeof (struct BINCROSS), comp_bincross);
				
				for (k = 1, k1 = 0; k < nx; k++, k1++) {	/* Process the intervals, getting mid-points and using that to get bin */
					dx = X[k].x - X[k1].x;
					if (s->geographic && dx < -180.0)
						x = 0.5 * (X[k].x + (X[k1].x - 360.0));
					else if (s->geographic && dx > +180.0)
						x = 0.5 * (X[k].x - 360.0 + X[k1].x);
					else
						x = 0.5 * (X[k].x + X[k1].x);
					y = 0.5 * (X[k].y + X[k1].y);
					if (s->geographic && fabs (y) > y_max) y = copysign (y_max, y);
					x2sys_err_fail (x2sys_bix_get_ij (x, y, &ii, &jj, &B, &ij), "");
					dist_bin[ij] += (GMT_distance_func) (X[k].x, X[k].y, X[k1].x, X[k1].y);
					B.binflag[ij] |= nav_flag;		/* Only update nav flags we have not been here already */
				}
			}
			last_bin_ij = this_bin_ij;
			last_bin_i = this_bin_i;
			last_bin_j = this_bin_j;
		}

		x2sys_free_data (data, s->n_fields, &p);

		/* Time for bin index output */

		sprintf (buffer, "> %s\n", trk_name[i]);
		GMT_fputs (buffer, GMT_stdout);
		for (ij = 0; ij < B.nm_bin; ij++) {
			if (B.binflag[ij] == 0) continue;
			x = B.x_min + ((ij % B.nx_bin) + 0.5) * B.bin_x;
			y = B.y_min + ((ij / B.nx_bin) + 0.5) * B.bin_y;
			GMT_ascii_output_one (GMT_stdout, x, 0);
			GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
			GMT_ascii_output_one (GMT_stdout, y, 1);
			sprintf (buffer, "%s%ld%s%u", gmtdefs.field_delimiter, ij, gmtdefs.field_delimiter, B.binflag[ij]);
			GMT_fputs (buffer, GMT_stdout);
			if (trackline) {
				GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
				GMT_ascii_output_one (GMT_stdout, dist_bin[ij], 2);
			}
			GMT_fputs ("\n", GMT_stdout);
		}

		if (trackline) GMT_free ((void *)dist_km);
	}

	GMT_free ((void *)X);
	x2sys_end (s);
	GMT_free ((void *)B.binflag);
	if (trackline) GMT_free ((void *)dist_bin);
	x2sys_free_list (trk_name, (int)n_tracks);

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

int outside (double x, double y, struct X2SYS_BIX *B, GMT_LONG geo)
{
	if (y < B->y_min || y > B->y_max) return (1);
	if (geo) {	/* Geographic data with periodic longitudes */
		while (x < B->x_min) x += 360.0;
		while (x > B->x_max) x -= 360.0;
		if (x < B->x_min) return (1);
	}
	else {	/* Plain Cartesian test */
		if (x < B->x_min || x > B->x_max) return (1);
	}
	return (0);	/* Inside */
}

unsigned int get_data_flag (double *data[], GMT_LONG j, struct X2SYS_INFO *s)
{
	int i;
	unsigned int bit, flag;
	for (i = flag = 0, bit = 1; i < s->n_fields; i++, bit <<= 1) {
		if (GMT_is_dnan (data[i][j])) continue;	/* NaN, so no data here */
		flag |= bit;
	}
	return (flag);
}

int comp_bincross (const void *p1, const void *p2)
{
	struct BINCROSS *a, *b;

	a = (struct BINCROSS *)p1;
	b = (struct BINCROSS *)p2;
	if (a->d < b->d) return (-1);
	if (a->d > b->d) return (1);
	return (0);
}
