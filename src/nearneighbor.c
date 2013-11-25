/*--------------------------------------------------------------------
 *	$Id: nearneighbor.c 10075 2013-07-23 18:39:24Z pwessel $
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
 * Based on a specified grid size, nearneighbor reads an xyz file and
 * determines the nearest points to each node in sectors.  The default
 * looks for the nearest point for each quadrant.  The points must also
 * be within a maximum search-radius from the node.  For the nodes that
 * have a full set of nearest neighbors, a weighted average value is
 * computed.  New feature is full support for boundary conditions so
 * that geographic or periodic conditions are explicitly dealt with
 * in the sense that a data point may wrap around to serve as a
 * constraint on the other side of the periodic boundary.
 *
 * Author:	Paul Wessel
 * Date:	14-JUL-2000
 * Version:	4
 */
 
#include "gmt.h"

#define NN_DEF_SECTORS	4

struct NEARNEIGHBOR_CTRL {	/* All control options for this program (except common args) */
	/* active is TRUE if the option has been activated */
	struct I {	/* -Idx[/dy] */
		GMT_LONG active;
		double xinc, yinc;
	} I;
	struct E {	/* -E<empty> */
		GMT_LONG active;
		double value;
	} E;
	struct F {	/* -F */
		GMT_LONG active;
	} F;
	struct G {	/* -G<grdfile> */
		GMT_LONG active;
		char *file;
	} G;
	struct L {	/* -L<flag> */
		GMT_LONG active;
		char mode[4];
	} L;
	struct N {	/* -N<sectors> */
		GMT_LONG active;
		GMT_LONG sectors, min_sectors;
	} N;
	struct S {	/* -S<radius>[m|c|k|K] */
		GMT_LONG active;
		double radius;
		char unit;
	} S;
	struct W {	/* -W */
		GMT_LONG active;
	} W;
};

struct NEARNEIGHBOR_NODE {	/* Structure with point id and distance pairs for all sectors */
	float *distance;	/* Distance of nearest datapoint to this node per sector */
	GMT_LONG *datum;		/* Point id of this data point */
};

struct NEARNEIGHBOR_POINT {	/* Structure with input data constraints */
	float x, y, z, w;
};

int main (int argc, char **argv)
{
	GMT_LONG i, j, k, i0, j0, ij, n, nm, n_alloc = GMT_CHUNK, n_set, n_filled, n_almost, n_none;
	GMT_LONG *di = NULL, dj, sector, n_read, n_fields, ii, jj;
	GMT_LONG n_files = 0, n_args, fno, n_expected_fields, distance_flag = 0;
	GMT_LONG max_di, actual_max_di, n_req, x_wrap, y_wrap;

	GMT_LONG error = FALSE, done = FALSE, nofile = TRUE;
	GMT_LONG wrap_180, replicate_x, replicate_y;

	double weight, weight_sum, grd_sum, *x0 = NULL, *y0 = NULL, dx, dy, delta, distance = 0.0, factor;
	double *in = NULL, *shrink = NULL, x_left, x_right, y_top, y_bottom, xinc2, yinc2, idx, idy;
	double half_y_width, y_width, half_x_width, x_width, three_over_radius;

	float *grd = NULL;

	char line[BUFSIZ];

	FILE *fp = NULL;

	struct GRD_HEADER header;
	struct GMT_EDGEINFO edgeinfo;
	struct NEARNEIGHBOR_NODE **grid_node = NULL;
	struct NEARNEIGHBOR_POINT  *point = NULL;
	struct NEARNEIGHBOR_CTRL *Ctrl = NULL;

	void assign_node (struct NEARNEIGHBOR_NODE **node, GMT_LONG n_sector, GMT_LONG sector, double distance, GMT_LONG id);
	void free_node (struct NEARNEIGHBOR_NODE *node);
	void *New_nearneighbor_Ctrl (), Free_nearneighbor_Ctrl (struct NEARNEIGHBOR_CTRL *C);

	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct NEARNEIGHBOR_CTRL *) New_nearneighbor_Ctrl ();	/* Allocate and initialize a new control structure */

	GMT_boundcond_init (&edgeinfo);

	GMT_grd_init (&header, argc, argv, FALSE);

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				/* Common parameters */

				case 'H':
				case 'R':
				case 'V':
				case ':':
				case 'b':	/* Input triplets [quadruplets] are binary, not ascii */
				case 'f':
				case '\0':
					error += GMT_parse_common_options (argv[i], &header.x_min, &header.x_max, &header.y_min, &header.y_max);
					break;

				/* Supplemental parameters */

				case 'I':
					Ctrl->I.active = TRUE;
					if (GMT_getinc (&argv[i][2], &Ctrl->I.xinc, &Ctrl->I.yinc)) {
						GMT_inc_syntax ('I', 1);
						error = TRUE;
					}
					break;
				case 'E':
					Ctrl->E.active = TRUE;
					if (!argv[i][2]) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -E option:  Must specify value or NaN\n", GMT_program);
						error++;
					}
					else
						Ctrl->E.value = (argv[i][2] == 'N' || argv[i][2] == 'n') ? GMT_d_NaN : atof (&argv[i][2]);
					break;
				case 'F':
					Ctrl->F.active = TRUE;
					break;
				case 'G':
					Ctrl->G.active = TRUE;
					Ctrl->G.file = strdup (&argv[i][2]);
					break;
				case 'L':
					if (argv[i][2]) {
						Ctrl->L.active = TRUE;
						strncpy (Ctrl->L.mode, &argv[i][2], (size_t)4);
					}
					else {
						GMT_io.in_col_type[0] = GMT_io.out_col_type[0] = GMT_IS_LON;
						GMT_io.in_col_type[1] = GMT_io.out_col_type[1] = GMT_IS_LAT;
						fprintf (stderr, "%s: Option -L is obsolete (but is processed correctly).  Please use -f instead\n", GMT_program);
					}
					break;
				case 'N':	/* -N[sectors[/minsectors]] */
					Ctrl->N.active = TRUE;
					n = sscanf (&argv[i][2], "%" GMT_LL "d/%" GMT_LL "d", &Ctrl->N.sectors, &Ctrl->N.min_sectors);
					if (n < 1) Ctrl->N.sectors = NN_DEF_SECTORS;	/* Just gave -N with no args means -N4/4 */
					if (n < 2) Ctrl->N.min_sectors = irint (Ctrl->N.sectors / 2.0);	/* Giving just -N<sectors> means -N<sectors>/(<sectors>/2) */
					break;
				case 'S':
					Ctrl->S.active = TRUE;
					Ctrl->S.radius = GMT_getradius (&argv[i][2]);
					Ctrl->S.unit = argv[i][strlen(argv[i])-1];
					break;
				case 'W':
					Ctrl->W.active = TRUE;
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
		fprintf (stderr, "nearneighbor %s - A \"Nearest neighbor\" gridding algorithm\n\n", GMT_VERSION);
		fprintf(stderr, "usage: nearneighbor [xyzfile(s)] -G<out_grdfile> %s\n", GMT_I_OPT);
		fprintf(stderr, "\t-N<sectors>[/<min_sectors>] %s -S<radius>[m|c|k|K] [-E<empty>] [-F]\n", GMT_Rgeo_OPT);
		fprintf(stderr, "\t[%s] [-L<flags>] [-V ] [-W] [%s] [%s] [%s]\n\n", GMT_H_OPT, GMT_t_OPT, GMT_bi_OPT, GMT_f_OPT);
		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);
		fprintf(stderr, "\t-G name of output grid.\n");
		GMT_inc_syntax ('I', 0);
		fprintf(stderr, "\t-N Set number of sectors and the minimum number of sectors with data required for averaging.\n");
		fprintf(stderr, "\t   If <min_sectors> is omitted it defaults to ~50%% of <sectors>.\n");
		fprintf(stderr, "\t   Default is -N%d/%d, i.e., a quadrant search, requiring all sectors to be filled.\n", NN_DEF_SECTORS, NN_DEF_SECTORS);
		GMT_explain_option ('R');
		fprintf(stderr, "\t-S sets search radius in -R, -I units; append m or c for minutes or seconds.\n");
		fprintf(stderr, "\t   Append k for km (implies -R,-I in degrees), use flat Earth approximation.\n");
		fprintf(stderr, "\t   Append K for km (implies -R,-I in degrees), use exact geodesic distances.\n");
		fprintf(stderr, "\t   If the current ELLIPSOID is sperical then great circle distances are used.\n");
		fprintf(stderr, "\n\tOPTIONS:\n");
		fprintf(stderr, "\t-E value to use for empty nodes [Default is NaN].\n");
		fprintf(stderr, "\t-F Force pixel registration [Default is gridline registration].\n");
		GMT_explain_option ('H');
		fprintf(stderr, "\t-L sets boundary conditions.  <flags> can be either\n");
		fprintf(stderr, "\t   g for geographic boundary conditions, or one or both of\n");
		fprintf(stderr, "\t   x for periodic boundary conditions on x,\n");
		fprintf(stderr, "\t   y for periodic boundary conditions on y.\n");
		GMT_explain_option ('V');
		fprintf(stderr, "\t-W input file has observation weights in 4th column.\n");
		GMT_explain_option (':');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf(stderr, "\t   Default is 3 (or 4 if -W is set) columns.\n");
		GMT_explain_option ('f');
		GMT_explain_option ('.');

		exit (EXIT_FAILURE);
	}

	GMT_check_lattice (&Ctrl->I.xinc, &Ctrl->I.yinc, &Ctrl->F.active, &Ctrl->I.active);

	if (!project_info.region_supplied) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify -R option\n", GMT_program);
		error++;
	}
	if (!Ctrl->G.file) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR option -G:  Must specify output file\n", GMT_program);
		error++;
	}
	if (Ctrl->N.sectors <= 0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -N option:  Must specify a positive number of sectors\n", GMT_program);
		error++;
	}
	if (Ctrl->S.radius <= 0.0 || GMT_is_dnan (Ctrl->S.radius)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Radius is NaN, zero or negative\n", GMT_program);
		error++;
	}
	if (Ctrl->I.xinc <= 0.0 || Ctrl->I.yinc <= 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -I option.  Must specify positive increment(s)\n", GMT_program);
		error++;
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
	n_req = (Ctrl->W.active) ? 4 : 3;
	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = n_req;
	if (GMT_io.binary[GMT_IN] && n_req > GMT_io.ncol[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  binary input data must have at least %ld columns\n", GMT_program, n_req);
		error++;
	}
	if (Ctrl->L.active && GMT_boundcond_parse (&edgeinfo, Ctrl->L.mode)) error++;

	if (error) exit (EXIT_FAILURE);

	if (GMT_io.binary[GMT_IN] && gmtdefs.verbose) {
		char *type[2] = {"double", "single"};
		fprintf (stderr, "%s: Expects %ld-column %s-precision binary data\n", GMT_program, GMT_io.ncol[GMT_IN], type[GMT_io.single_precision[GMT_IN]]);
	}
	if (Ctrl->S.unit == 'k') distance_flag = 1;
	if (Ctrl->S.unit == 'K') distance_flag = 2;

	if (Ctrl->N.sectors < Ctrl->N.min_sectors) Ctrl->N.min_sectors = Ctrl->N.sectors;	/* Minimum cannot be larger than desired */
	
	header.x_inc = Ctrl->I.xinc;
	header.y_inc = Ctrl->I.yinc;
	header.node_offset = (int)Ctrl->F.active;
	GMT_RI_prepare (&header);	/* Ensure -R -I consistency and set nx, ny */
	GMT_err_fail (GMT_grd_RI_verify (&header, 1), Ctrl->G.file);

	n_expected_fields = (GMT_io.binary[GMT_IN]) ? GMT_io.ncol[GMT_IN] : 3 + Ctrl->W.active;

        if (n_files > 0)
        	nofile = FALSE;
        else
        	n_files = 1;
        n_args = (argc > 1) ? argc : 2;
        
	header.xy_off = 0.5 * header.node_offset;
	xinc2 = header.xy_off * header.x_inc;
	yinc2 = header.xy_off * header.y_inc;
	idx = 1.0 / header.x_inc;
	idy = 1.0 / header.y_inc;

	GMT_boundcond_param_prep (&header, &edgeinfo);

	if (gmtdefs.verbose) fprintf (stderr, "%s: Grid dimensions are nx = %d, ny = %d\n", GMT_program, header.nx, header.ny);

	nm = GMT_get_nm (header.nx, header.ny);
	grid_node = (struct NEARNEIGHBOR_NODE **) GMT_memory (VNULL, (size_t)nm, sizeof (struct NEARNEIGHBOR_NODE *), GMT_program);
	point = (struct NEARNEIGHBOR_POINT *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (struct NEARNEIGHBOR_POINT), GMT_program);

	di = (GMT_LONG *) GMT_memory (VNULL, header.ny, sizeof (GMT_LONG), GMT_program);
	shrink = (double *) GMT_memory (VNULL, header.ny, sizeof (double), GMT_program);

	x0 = (double *) GMT_memory (VNULL, header.nx, sizeof (double), GMT_program);
	y0 = (double *) GMT_memory (VNULL, header.ny, sizeof (double), GMT_program);
	for (i = 0; i < header.nx; i++) x0[i] = header.x_min + i * header.x_inc + xinc2;
	for (j = 0; j < header.ny; j++) y0[j] = header.y_max - j * header.y_inc - yinc2;
	if (distance_flag) {	/* Input data is geographical */
		if (!GMT_IS_SPHERICAL) distance_flag = 3;	/* Use geodesics */
		max_di = (int) (ceil (header.nx / 2.0) + 0.1);
		actual_max_di = 0;
		for (j = 0; j < header.ny; j++) {
			shrink[j] = cosd (y0[j]);
			di[j] = (fabs (y0[j]) == 90.0) ? max_di : (int)(ceil (Ctrl->S.radius / (project_info.DIST_KM_PR_DEG * header.x_inc * shrink[j])) + 0.1);
			if (di[j] > max_di) di[j] = max_di;
			if (di[j] > actual_max_di) actual_max_di = di[j];
		}
		dj = (int) (ceil (Ctrl->S.radius / (project_info.DIST_KM_PR_DEG * header.y_inc)) + 0.1);
	}
	else {	/* Plain Cartesian data */
		max_di = (int) (ceil (Ctrl->S.radius * idx) + 0.1);
		for (j = 0; j < header.ny; j++) di[j] = max_di;
		dj = (int) (ceil (Ctrl->S.radius * idy) + 0.1);
		actual_max_di = max_di;
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

	factor = Ctrl->N.sectors / (2.0 * M_PI);

	/* To allow data points falling outside -R but within the search radius we extend the data domain in all directions */

	x_left = header.x_min;	x_right = header.x_max;	/* This is what -R says */
	if (GMT_io.in_col_type[0] != GMT_IS_LON || !GMT_grd_is_global (&header)) {
		x_left  -= actual_max_di * header.x_inc;	/* OK to extend x-domain since not a periodic geographic grid */
		x_right += actual_max_di * header.x_inc;
	}
	y_top = header.y_max + dj * header.y_inc;	y_bottom = header.y_min - dj * header.y_inc;
	if (GMT_io.in_col_type[1] == GMT_IS_LAT) {	/* For geographic grids we must ensure the extended y-domain is physically possible */
		if (y_bottom < -90.0) y_bottom = -90.0;
		if (y_top > 90.0) y_top = 90.0;
	}
	x_width = header.x_max - header.x_min;		y_width = header.y_max - header.y_min;
	half_x_width = 0.5 * x_width;			half_y_width = 0.5 * y_width;
	n = n_read = 0;
	replicate_x = (edgeinfo.nxp && !header.node_offset);	/* Gridline registration has duplicate column */
	replicate_y = (edgeinfo.nyp && !header.node_offset);	/* Gridline registration has duplicate row */
	x_wrap = header.nx - 1;			/* Add to node index to go to right column */
	y_wrap = (header.ny - 1) * header.nx;	/* Add to node index to go to bottom row */

#ifdef DEBUG
	GMT_memtrack_off (GMT_mem_keeper);
#endif
	for (fno = 1; !done && fno < n_args; fno++) {	/* Loop over input files, if any */
		if (!nofile && argv[fno][0] == '-') continue;

		if (nofile) {	/* Just read standard input */
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

		while ((n_fields = GMT_input (fp, &n_expected_fields, &in)) >= 0 && !(GMT_io.status & GMT_IO_EOF)) {	/* Not yet EOF */

			n_read++;

			if (GMT_io.status & GMT_IO_MISMATCH) {
				fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld (skipped)\n", GMT_program, n_fields, n_expected_fields, n_read);
				continue;
			}

			if (GMT_is_dnan (in[GMT_Z])) continue;					/* Skip if z = NaN */
			if (GMT_y_is_outside (in[GMT_Y], y_bottom, y_top)) continue;		/* Outside y-range */
			if (GMT_x_is_outside (&in[GMT_X], x_left, x_right)) continue;	/* Outside x-range (or longitude) */

			point[n].x = (float)in[GMT_X];
			point[n].y = (float)in[GMT_Y];
			point[n].z = (float)in[GMT_Z];
			if (Ctrl->W.active) point[n].w = (float)in[3];

			/* Find indices of the node closest to this data point */

			i0 = GMT_x_to_i(in[GMT_X], header.x_min, header.x_inc, header.xy_off, header.nx);
			j0 = GMT_y_to_j(in[GMT_Y], header.y_min, header.y_inc, header.xy_off, header.ny);

			/* Loop over all nodes within radius of this node */

			for (j = j0 - dj; j <= (j0 + dj); j++) {

				jj = (int)j;
				if (GMT_y_out_of_bounds (&jj, &header, &edgeinfo, &wrap_180)) continue;	/* Outside y-range */

				for (i = i0 - di[jj]; i <= (i0 + di[jj]); i++) {

					ii = (int)i;
					if (GMT_x_out_of_bounds (&ii, &header, &edgeinfo, wrap_180)) continue;	/* Outside x-range */ 

					/* Here, (ii,jj) is index of a node (k) inside the grid */

					distance = (GMT_distance_func) (x0[ii], y0[jj], in[GMT_X], in[GMT_Y]);

					if (distance > Ctrl->S.radius) continue;	/* Data constraint is too far from this node */

					k = GMT_IJ (jj, ii, header.nx);
					dx = in[GMT_X] - x0[ii];	dy = in[GMT_Y] - y0[jj];

					/* Check for wrap-around in x or y.  This should only occur if the
					   search radius is larger than 1/2 the grid width/height so that
					   the shortest distance is going through the periodic boundary.
					   For longitudes the dx obviously cannot exceed 180 (half_x_width)
					   since we could then go the other direction instead.
					*/
					if (edgeinfo.nxp && fabs (dx) > half_x_width) dx -= copysign (x_width, dx);
					if (edgeinfo.nyp && fabs (dy) > half_y_width) dy -= copysign (y_width, dy);

					/* OK, this point should constrain this node.  Calculate which sector and assign the value */

					sector = ((int)((d_atan2 (dy, dx) + M_PI) * factor)) % Ctrl->N.sectors;
					assign_node (&grid_node[k], Ctrl->N.sectors, sector, distance, n);

					/* With periodic, gridline-registered grids there are duplicate rows and/or columns
					   so we may have to assign the point to more than one node.  The next section deals
					   with this situation.
					*/

					if (replicate_x) {	/* Must check if we have to replicate a column */
						if (ii == 0) 	/* Must replicate left to right column */
							assign_node (&grid_node[k+x_wrap], Ctrl->N.sectors, sector, distance, n);
						else if (ii == edgeinfo.nxp)	/* Must replicate right to left column */
							assign_node (&grid_node[k-x_wrap], Ctrl->N.sectors, sector, distance, n);
					}
					if (replicate_y) {	/* Must check if we have to replicate a row */
						if (jj == 0)	/* Must replicate top to bottom row */
							assign_node (&grid_node[k+y_wrap], Ctrl->N.sectors, sector, distance, n);
						else if (jj == edgeinfo.nyp)	/* Must replicate bottom to top row */
							assign_node (&grid_node[k-y_wrap], Ctrl->N.sectors, sector, distance, n);
					}
				}
			}
			n++;
			if (gmtdefs.verbose && !(n%1000)) fprintf (stderr, "%s: Processed record %10ld\r", GMT_program, n);
			if (n == n_alloc) {
				n_alloc <<= 1;
				point = (struct NEARNEIGHBOR_POINT *) GMT_memory ((void *)point, (size_t)n_alloc, sizeof (struct NEARNEIGHBOR_POINT), GMT_program);
			}
		}
		if (fp != GMT_stdin) GMT_fclose (fp);
	}
	if (gmtdefs.verbose) fprintf (stderr, "%s: Processed record %10ld\n", GMT_program, n);

#ifdef DEBUG
	GMT_memtrack_on (GMT_mem_keeper);
#endif
	point = (struct NEARNEIGHBOR_POINT *) GMT_memory ((void *)point, (size_t)n, sizeof (struct NEARNEIGHBOR_POINT), GMT_program);
	grd = (float *) GMT_memory (VNULL, (size_t)(header.nx * header.ny), sizeof (float), GMT_program);

	/* Compute weighted averages based on the nearest neighbors */

	n_set = n_almost = n_none = 0;

	if (!Ctrl->E.active) Ctrl->E.value = GMT_d_NaN;
	three_over_radius = 3.0 / Ctrl->S.radius;

#ifdef DEBUG
	GMT_memtrack_off (GMT_mem_keeper);
#endif
	for (j = ij = 0; j < header.ny; j++) {
		for (i = 0; i < header.nx; i++, ij++) {

			if (!grid_node[ij]) {	/* No nearest neighbors, set to empty and goto next node */
				n_none++;
				grd[ij] = (float)Ctrl->E.value;
				continue;
			}

			for (k = 0, n_filled = 0; k < Ctrl->N.sectors; k++) if (grid_node[ij]->datum[k] >= 0) n_filled++;
			if (n_filled < Ctrl->N.min_sectors) { 	/* Not minimum set of neighbors in all sectors, set to empty and goto next node */
				n_almost++;
				grd[ij] = (float)Ctrl->E.value;
				free_node (grid_node[ij]);
				continue;
			}

			/* OK, here we have enough data and need to calculate the weighted value */

			n_set++;
			weight_sum = grd_sum = 0.0;	/* Initialize sums */
			for (k = 0; k < Ctrl->N.sectors; k++) {
				if (grid_node[ij]->datum[k] >= 0) {
					delta = three_over_radius * grid_node[ij]->distance[k];
					weight = 1.0 / (1.0 + delta * delta);	/* This is distance weight */
					if (Ctrl->W.active) weight *= point[grid_node[ij]->datum[k]].w;	/* This is observation weight */
					grd_sum += weight * point[grid_node[ij]->datum[k]].z;
					weight_sum += weight;
				}
			}
			grd[ij] = (float)(grd_sum / weight_sum);
			free_node (grid_node[ij]);
		}
		if (gmtdefs.verbose) fprintf (stderr, "%s: Gridded row %10ld\r", GMT_program, j);
	}
	if (gmtdefs.verbose) fprintf (stderr, "%s: Gridded row %10ld\n", GMT_program, j);
#ifdef DEBUG
	GMT_memtrack_on (GMT_mem_keeper);
#endif

	GMT_err_fail (GMT_write_grd (Ctrl->G.file, &header, grd, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), Ctrl->G.file);

	if (gmtdefs.verbose) {
		sprintf (line, "%s)\n", gmtdefs.d_format);
		fprintf (stderr, "%s: %ld nodes were assigned an average value\n", GMT_program, n_set);
		fprintf (stderr, "%s: %ld nodes failed sector criteria and %ld nodes had no neighbor points (all set to ", GMT_program, n_almost, n_none);
		(GMT_is_dnan (Ctrl->E.value)) ? fprintf (stderr, "NaN)\n") : fprintf (stderr, line, Ctrl->E.value);
	}

	GMT_free ((void *)grd);
	GMT_free ((void *)point);
	GMT_free ((void *)grid_node);
	GMT_free ((void *)shrink);
	GMT_free ((void *)di);
	GMT_free ((void *)x0);
	GMT_free ((void *)y0);

	Free_nearneighbor_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

struct NEARNEIGHBOR_NODE *add_new_node (GMT_LONG n)
{
	struct NEARNEIGHBOR_NODE *new;
	
	new = (struct NEARNEIGHBOR_NODE *) GMT_memory (VNULL, (size_t)1, sizeof (struct NEARNEIGHBOR_NODE), GMT_program);
	new->distance = (float *) GMT_memory (VNULL, (size_t)n, sizeof (float), GMT_program);
	new->datum = (GMT_LONG *) GMT_memory (VNULL, (size_t)n, sizeof (GMT_LONG), GMT_program);
	while (n > 0) new->datum[--n] = -1;

	return (new);
}

void assign_node (struct NEARNEIGHBOR_NODE **node, GMT_LONG n_sector, GMT_LONG sector, double distance, GMT_LONG id)
{
	struct NEARNEIGHBOR_NODE *add_new_node(GMT_LONG n);

	/* Allocates node space if not already used and updates the value if closer to node */

	if (!(*node)) *node = add_new_node (n_sector);
	if ((*node)->datum[sector] == -1 || (*node)->distance[sector] > distance) {
		(*node)->distance[sector] = (float)distance;
		(*node)->datum[sector] = id;
	}
}

void free_node (struct NEARNEIGHBOR_NODE *node)
{
	/* Frees allocated node space */

	if (!node) return;	/* Nothing to do */
	GMT_free ((void *)node->distance);
	GMT_free ((void *)node->datum);
	GMT_free ((void *)node);
}

void *New_nearneighbor_Ctrl () {	/* Allocate and initialize a new control structure */
	struct NEARNEIGHBOR_CTRL *C;

	C = (struct NEARNEIGHBOR_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct NEARNEIGHBOR_CTRL), "New_nearneighbor_Ctrl");

	/* Initialize values whose defaults are not 0/FALSE/NULL */
	C->N.sectors = NN_DEF_SECTORS;
	C->N.min_sectors = NN_DEF_SECTORS;
	return ((void *)C);
}

void Free_nearneighbor_Ctrl (struct NEARNEIGHBOR_CTRL *C) {	/* Deallocate control structure */
	if (C->G.file) free ((void *)C->G.file);
	GMT_free ((void *)C);
}
