/*--------------------------------------------------------------------
 *	$Id: sphdistance.c 10115 2013-11-01 20:41:11Z pwessel $
 *
 *	Copyright (c) 2008-2013 by P. Wessel and W. H. F. Smith
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
 * Spherical nearest distances - via Voronoi polygons.  We read input
 * data, assumed to be things like coastlines, and want to create a grid
 * with distances to the nearest line.  The approach here is to break
 * the data into voronoi polygons and then visit all nodes inside each
 * polygon and use geodesic distance calculation from each node to the
 * unique Voronoi interior data node.
 * Relies on STRIPACK Fortran F77 library (Renka, 1997). Reference:
 * Renka, R, J,, 1997, Algorithm 772: STRIPACK: Delaunay Triangulation
 *  and Voronoi Diagram on the Surface of a Sphere, AMC Trans. Math.
 *  Software, 23 (3), 416-434.
 * We translate to C using f2c -r8 and link with -lf2c
 *
 * Author:      Paul Wessel
 * Date:        16-FEB-2008
 *
 */
 
#include "gmt.h"
#include "sph.h"

struct SPHDISTANCE_CTRL {
	struct C {	/* -C */
		GMT_LONG active;
	} C;
	struct D {	/* -D */
		GMT_LONG active;
	} D;
	struct E {	/* -E */
		GMT_LONG active;
	} E;
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
	struct L {	/* -L<unit>] */
		GMT_LONG active;
		char unit;
	} L;
	struct N {	/* -N */
		GMT_LONG active;
		char *file;
	} N;
	struct Q {	/* -Q */
		GMT_LONG active;
		char *file;
	} Q;
};

void prepare_polygon (struct GMT_LINE_SEGMENT *P);

int main (int argc, char **argv)
{
	GMT_LONG i, j, n = 0, n_dup = 0, duplicate_col;
	GMT_LONG ij, ii, s_row, n_row, w_col, e_col, side, nx1;
	GMT_LONG n_expected_fields, n_args, n_fields, fno, n_files = 0, proj_type = 0;

	GMT_LONG error = FALSE, nofile = TRUE, done = FALSE, first = FALSE;

	double sx, sy, cx, cy, first_x = 0.0, first_y = 0.0, *in = NULL;
	double *xx = NULL, *yy = NULL, *zz = NULL, *lon = NULL, *lat = NULL;
	double d_scale, *grid_lon = NULL, *grid_lat = NULL;
	GMT_LONG node, vertex, node_stop, node_new, vertex_new, node_last, vertex_last;
	size_t n_alloc = GMT_CHUNK, p_alloc = GMT_TINY_CHUNK, nm, n_set = 0;
	float *distgrid = NULL;

	PFD distance_func;

	char line[BUFSIZ], *not_used = NULL;
	FILE *fp = NULL;

	struct SPHDISTANCE_CTRL *Ctrl = NULL;
	struct STRIPACK T;
	struct GRD_HEADER h;
	struct GMT_LINE_SEGMENT *P = NULL;
	struct GMT_TABLE *Table = NULL;
	struct STRIPACK_VORONOI *V = NULL;
	void *New_sphdistance_Ctrl (), Free_sphdistance_Ctrl (struct SPHDISTANCE_CTRL *C);
	
	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct SPHDISTANCE_CTRL *)New_sphdistance_Ctrl ();	/* Allocate and initialize a new control structure */
	GMT_grd_init (&h, argc, argv, FALSE);

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
					error += GMT_parse_common_options (argv[i], &h.x_min, &h.x_max, &h.y_min, &h.y_max);
					break;
        
				/* Supplemental parameters */
        
				case 'C':
					Ctrl->C.active = TRUE;
					break;
				case 'D':
					Ctrl->D.active = TRUE;
					break;
				case 'E':
					Ctrl->E.active = TRUE;
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
				case 'L':
					Ctrl->L.active = TRUE;
					if (!(argv[i][2] && strchr ("ekmndEKMND", argv[i][2]))) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR.  Expected -Le|k|m|n|d]\n", GMT_program);
						error = TRUE;
					}
					else
						Ctrl->L.unit = toupper (argv[i][2]);	/* Make sure we pass upper so Geodesic is possible */
					break;
				case 'N':
					Ctrl->N.active = TRUE;
					Ctrl->N.file = strdup (&argv[i][2]);
					break;
				case 'Q':
					Ctrl->Q.active = TRUE;
					Ctrl->Q.file = strdup (&argv[i][2]);
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
		fprintf (stderr, "sphdistance %s - Make grid of distances to nearest feature on a sphere\n\n", GMT_VERSION);
		fprintf (stderr, "usage: sphdistance [<infiles>] -G<grdfile> %s [-C] [-D] [-E] [-F] [%s] [-L<unit>]\n", GMT_I_OPT, GMT_H_OPT);
		fprintf (stderr, "\t[-N<nodefile>] [-Q<voronoifile>] [-V] [%s] [%s] [%s]\n\n", GMT_t_OPT, GMT_b_OPT, GMT_m_OPT);
                
		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);
		fprintf (stderr, "\t-G Specify file name for output distance grid file.\n");
		GMT_inc_syntax ('I', 0);
                
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\tinfiles (in ASCII) has 2 or more columns.  If no file(s) is given, standard input is read (but see -Q).\n");
		fprintf (stderr, "\t-C Conserve memory (Converts lon/lat <--> x/y/z when needed) [store both in memory]. Not used with -Q.\n");
		fprintf (stderr, "\t-D Used with -M to skip repeated input vertex at the end of a closed segment\n");
		fprintf (stderr, "\t-E Assign to grid nodes the Voronoi polygon ID [Calculate distances]\n");
		fprintf (stderr, "\t-F Force pixel registration for output grid [Default is gridline orientation]\n");
		GMT_explain_option ('H');
		fprintf (stderr, "\t-L Specify distance unit as m(e)ter, (k)m, (m)ile, (n)autical mile, or (d)egree.\n");
		fprintf (stderr, "\t   Calculations uses spherical approximations.  Default unit is meters.\n");
		fprintf (stderr, "\t   Set ELLIPSOID to WGS-84 to get geodesic distances.\n");
		fprintf (stderr, "\t-N Specify node file for the Voronoi polygons (sphtriangulate -N output)\n");
		fprintf (stderr, "\t-Q Specify file with Voronoi polygons in sphtriangulate -Qv format\n");
		fprintf (stderr, "\t   [Default performs Voronoi construction on input data first]\n");
		GMT_explain_option ('R');
		fprintf (stderr, "\t   If no region is specified we default to the entire world [-Rg]\n");
		GMT_explain_option ('V');
		GMT_explain_option (':');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf (stderr, "\t   Default is 2 input columns\n");
		GMT_explain_option ('m');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	GMT_check_lattice (&Ctrl->I.xinc, &Ctrl->I.yinc, &Ctrl->F.active, &Ctrl->I.active);

	if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = 2;
	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] < 2) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data (-bi) must have at least 2 columns\n", GMT_program);
		error++;
	}
	if (Ctrl->Q.active && GMT_io.binary[GMT_IN] && !Ctrl->N.active) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data (-bi) with -Q also requires -N\n", GMT_program);
		error++;
	}
	if (Ctrl->D.active && !GMT_io.multi_segments[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  -D requires the -M[i] option\n", GMT_program);
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
	error += GMT_get_dist_scale (Ctrl->L.unit, &d_scale, &proj_type, &distance_func);
	
	if (error) exit (EXIT_FAILURE);

	if (GMT_io.binary[GMT_IN] && gmtdefs.verbose) {
		char *type[2] = {"double", "single"};
		fprintf (stderr, "%s: Expects %ld-column %s-precision binary data\n", GMT_program, GMT_io.ncol[GMT_IN], type[GMT_io.single_precision[GMT_IN]]);
	}
	if (!project_info.region_supplied) {
		h.x_min = 0.0;	h.x_max = 360.0;	h.y_min = -90.0;	h.y_max = 90.0;
	}

#ifdef SET_IO_MODE
	GMT_setmode (GMT_OUT);
#endif

	/* Now we are ready to take on some input values */

	if (Ctrl->Q.active) {	/* Expect a single file with Voronoi polygons */
		if (gmtdefs.verbose) fprintf (stderr, "%s: Read Volonoi polygons from %s ...", GMT_program, Ctrl->Q.file);
		GMT_import_table ((void *)(Ctrl->Q.file), GMT_IS_FILE, &Table, 0.0, FALSE, TRUE, TRUE);
		if (gmtdefs.verbose) fprintf (stderr, "Found %ld segments\n", Table->n_segments);
	 	lon = (double *) GMT_memory (VNULL, (size_t)Table->n_segments, sizeof (double), GMT_program);
	 	lat = (double *) GMT_memory (VNULL, (size_t)Table->n_segments, sizeof (double), GMT_program);
		if (Ctrl->N.active) {	/* Must get nodes from separate file */
			struct GMT_TABLE *NTable;
			GMT_io.ncol[GMT_IN] = 3;
			if (gmtdefs.verbose) fprintf (stderr, "%s: Read Nodes from %s ...", GMT_program, Ctrl->N.file);
			GMT_import_table ((void *)(Ctrl->N.file), GMT_IS_FILE, &NTable, 0.0, FALSE, TRUE, TRUE);
			for (node = 0; node < NTable->n_records; node++) {
				lon[node] = NTable->segment[0]->coord[GMT_X][node];
				lat[node] = NTable->segment[0]->coord[GMT_Y][node];
			}
			GMT_free_table (NTable);
			if (gmtdefs.verbose) fprintf (stderr, "Found %ld records\n", NTable->n_records);
		}
		else {
			for (node = 0; node < Table->n_segments; node++) {
				sscanf (Table->segment[node]->header, "%*s %*s %*d %lf %lf", &lon[node], &lat[node]);
			}
		}
	}
	else {	/* Must process input point/line data*/
		if (n_files > 0)
			nofile = FALSE;
		else
			n_files = 1;
		n_args = (argc > 1) ? argc : 2;
		n_expected_fields = (GMT_io.ncol[GMT_IN]) ? GMT_io.ncol[GMT_IN] : 2;

		n_alloc = GMT_CHUNK;
		xx = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (double), GMT_program);
		yy = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (double), GMT_program);
		zz = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (double), GMT_program);
		if (!Ctrl->C.active) {
			lon = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (double), GMT_program);
			lat = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (double), GMT_program);
		}
		n = 0;
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

			if (!nofile && gmtdefs.verbose) fprintf (stderr, "%s: Reading file %s\n", GMT_program, argv[fno]);

			if (GMT_io.io_header[GMT_IN]) for (i = 0; i < GMT_io.n_header_recs; i++) not_used = GMT_fgets (line, BUFSIZ, fp);

			while ((n_fields = GMT_input (fp, &n_expected_fields, &in)) >= 0 && !(GMT_io.status & GMT_IO_EOF)) {	/* Not yet EOF */

				if (GMT_io.status & GMT_IO_MISMATCH) {
					fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld\n", GMT_program, n_fields, n_expected_fields, n);
					exit (EXIT_FAILURE);
				}
				if (GMT_io.status & GMT_IO_SEGMENT_HEADER) {	/* Segment header, get next record */
					n_fields = GMT_input (fp, &n_expected_fields, &in);
					first_x = in[0];	first_y = in[1];
					first = TRUE;
				}
				if (Ctrl->D.active && !first) {	/* Look for duplicate point at end of segments that replicate start point */
					if (in[0] == first_x && in[1] == first_y) {	/* If any point after the first matches the first */
						n_dup++;
						continue;
					}
				}
				sincos (D2R * in[1], &sy, &cy);
				sincos (D2R * in[0], &sx, &cx);
				xx[n] = cy * cx;
				yy[n] = cy * sx;
				zz[n] = sy;
				if (!Ctrl->C.active) {
					lon[n] = in[GMT_X];
					lat[n] = in[GMT_Y];
				}
				n++;

				if (n == (int)n_alloc) {	/* Get more memory */
					n_alloc <<= 1;
					xx = (double *) GMT_memory ((void *)xx, (size_t)n_alloc, sizeof (double), GMT_program);
					yy = (double *) GMT_memory ((void *)yy, (size_t)n_alloc, sizeof (double), GMT_program);
					zz = (double *) GMT_memory ((void *)zz, (size_t)n_alloc, sizeof (double), GMT_program);
					if (!Ctrl->C.active) {
						lon = (double *) GMT_memory ((void *)lon, (size_t)n_alloc, sizeof (double), GMT_program);
						lat = (double *) GMT_memory ((void *)lat, (size_t)n_alloc, sizeof (double), GMT_program);
					}
				}
				first = FALSE;
			}

			if (fp != GMT_stdin) GMT_fclose (fp);
		}

		xx = (double *) GMT_memory ((void *)xx, (size_t)n, sizeof (double), GMT_program);
		yy = (double *) GMT_memory ((void *)yy, (size_t)n, sizeof (double), GMT_program);
		zz = (double *) GMT_memory ((void *)zz, (size_t)n, sizeof (double), GMT_program);
		if (!Ctrl->C.active) {
			lon = (double *) GMT_memory ((void *)lon, (size_t)n, sizeof (double), GMT_program);
			lat = (double *) GMT_memory ((void *)lat, (size_t)n, sizeof (double), GMT_program);
		}

		if (Ctrl->D.active && n_dup && gmtdefs.verbose) fprintf (stderr, "%s: Skipped %ld duplicate points in segments\n", GMT_program, n_dup);
		if (gmtdefs.verbose) fprintf (stderr, "%s: Do Voronoi construction using %ld points\n", GMT_program, n);

		memset ((void *)&T, 0, sizeof (struct STRIPACK));
		T.mode = VORONOI;
		stripack_lists (n, xx, yy, zz, gmtdefs.verbose, &T);	/* Do the basic triangulation */
		GMT_free ((void *)T.D.tri);
		if (Ctrl->C.active) {	/* Recompute lon,lat and set pointers */
			cart_to_geo (n, xx, yy, zz, xx, yy);	/* Revert to lon, lat */
			lon = xx;
			lat = yy;
		}
		GMT_free ((void *) zz);
	}
	
	/* OK, time to work on the distance grid */
	
	if (gmtdefs.verbose) fprintf (stderr, "%s: Start processing distance grid\n", GMT_program);
	h.node_offset = (int)Ctrl->F.active;
	h.x_inc = Ctrl->I.xinc;
	h.y_inc = Ctrl->I.yinc;
	GMT_RI_prepare (&h);	/* Ensure -R -I consistency and set nx, ny */
	GMT_err_fail (GMT_grd_RI_verify (&h, 1), Ctrl->G.file);

	h.xy_off = 0.5 * h.node_offset;
	nx1 = (h.node_offset) ? h.nx : h.nx - 1;
	duplicate_col = (GMT_360_RANGE (h.x_min, h.x_max) && h.node_offset == 0);	/* E.g., lon = 0 column should match lon = 360 column */
	nm = GMT_get_nm (h.nx, h.ny);
	distgrid = (float *) GMT_memory (VNULL, (size_t)nm, sizeof(float), GMT_program);
	grid_lon = (double *) GMT_memory (VNULL, h.nx, sizeof (double), GMT_program);
	grid_lat = (double *) GMT_memory (VNULL, h.ny, sizeof (double), GMT_program);
	for (i = 0; i < h.nx; i++) grid_lon[i] = GMT_i_to_x (i, h.x_min, h.x_max, h.x_inc, h.xy_off, h.nx);
	for (j = 0; j < h.ny; j++) grid_lat[j] = GMT_j_to_y (j, h.y_min, h.y_max, h.y_inc, h.xy_off, h.ny);
	
	if (Ctrl->Q.active) {	/* Pre-chewed, just get number of nodes */
		n = Table->n_segments;
	}
	else {	/* Need a single polygon structure that we reuse for each polygon */
		P = (struct GMT_LINE_SEGMENT *) GMT_memory (VNULL, 1, sizeof (struct GMT_LINE_SEGMENT), GMT_program);	/* Needed as pointer below */
		P->coord = (double **) GMT_memory (VNULL, 2, sizeof (double *), GMT_program);	/* Needed as pointers below */
		P->min = (double *) GMT_memory (VNULL, 2, sizeof (double), GMT_program);		/* Needed to hold min lon/lat */
		P->max = (double *) GMT_memory (VNULL, 2, sizeof (double), GMT_program);		/* Needed to hold max lon/lat */

		P->coord[GMT_X] = (double *) GMT_memory (VNULL, p_alloc, sizeof (double), GMT_program);
		P->coord[GMT_Y] = (double *) GMT_memory (VNULL, p_alloc, sizeof (double), GMT_program);
	
		V = &T.V;
	}
	for (node = 0; node < n; node++) {

		if (gmtdefs.verbose) fprintf (stderr, "%s: Processing polygon %7ld\r", GMT_program, node);
		if (Ctrl->Q.active) {	/* Just point to next polygon */
			P = Table->segment[node];
		}
		else {	/* Obtain current polygon from Voronoi listings */
			node_new = node_stop = V->lend[node];
			vertex_new = V->listc[node_new];

			/* Each iteration of this DO walks along one side of the polygon,
			   considering the subtriangle NODE --> VERTEX_LAST --> VERTEX. */

			vertex = 0;
		    do {
				node_last = node_new;
				node_new = V->lptr[node_last];
				vertex_last = vertex_new;
				vertex_new = V->listc[node_new];

				P->coord[GMT_X][vertex] = V->lon[vertex_last];
				P->coord[GMT_Y][vertex] = V->lat[vertex_last];
				if (P->coord[GMT_X][vertex] < 0.0) P->coord[GMT_X][vertex] += 360.0;
				if (P->coord[GMT_X][vertex] == 360.0) P->coord[GMT_X][vertex] = 0.0;
				vertex++;
				if (vertex == (int)p_alloc) {
					p_alloc <<= 1;
					P->coord[GMT_X] = (double *) GMT_memory ((void *)P->coord[GMT_X], p_alloc, sizeof (double), GMT_program);
					P->coord[GMT_Y] = (double *) GMT_memory ((void *)P->coord[GMT_Y], p_alloc, sizeof (double), GMT_program);
				}

				/* When we reach the vertex where we started, we are done with this polygon */
			} while (node_new != node_stop);
			P->coord[GMT_X][vertex] = P->coord[GMT_X][0];
			P->coord[GMT_Y][vertex] = P->coord[GMT_Y][0];
			vertex++;
			if (vertex == (int)p_alloc) {
				p_alloc <<= 1;
				P->coord[GMT_X] = (double *) GMT_memory ((void *)P->coord[GMT_X], p_alloc, sizeof (double), GMT_program);
				P->coord[GMT_Y] = (double *) GMT_memory ((void *)P->coord[GMT_Y], p_alloc, sizeof (double), GMT_program);
			}
			P->n_rows = vertex;
		}
		
		/* Here we have the polygon in P */
		
		prepare_polygon (P);	/* Determine the enclosing sector */

		s_row = GMT_y_to_j (P->min[GMT_Y], h.y_min, h.y_inc, h.xy_off, h.ny);
		n_row = GMT_y_to_j (P->max[GMT_Y], h.y_min, h.y_inc, h.xy_off, h.ny);
		w_col = GMT_x_to_i (P->min[GMT_X], h.x_min, h.x_inc, h.xy_off, h.nx);
		e_col = GMT_x_to_i (P->max[GMT_X], h.x_min, h.x_inc, h.xy_off, h.nx);

		for (j = n_row; j <= s_row; j++) {	/* For each scanline intersecting this polygon */
			for (ii = w_col; ii <= e_col; ii++) {	/* March along the scanline */
				i = (ii >= 0) ? ii : ii + nx1;
				side = GMT_inonout_sphpol (grid_lon[i], grid_lat[j], P);
				if (side == 0) continue;	/* Outside spherical polygon */
				ij = GMT_IJ(j,i,h.nx);
				distgrid[ij] = (Ctrl->E.active) ? (float)node : (float)(d_scale * (*distance_func) (grid_lon[i], grid_lat[j], lon[node], lat[node]));
				n_set++;
				if (duplicate_col) {	/* Duplicate the repeating column on the other side of this one */
					if (i == 0) distgrid[ij+nx1] = distgrid[ij], n_set++;
					else if (i == nx1) distgrid[ij-nx1] = distgrid[ij], n_set++;
				}
			}
		}
	}
	if (gmtdefs.verbose) fprintf (stderr, "%s: Processing polygon %7ld\n", GMT_program, node);
	
	if (!Ctrl->Q.active) {
		GMT_free ((void *)P->coord[GMT_X]);
		GMT_free ((void *)P->coord[GMT_Y]);
		GMT_free ((void *)P->min);
		GMT_free ((void *)P->max);
		GMT_free ((void *)P->coord);
		GMT_free ((void *)P);
		GMT_free ((void *)T.V.lon);
		GMT_free ((void *)T.V.lat);
		GMT_free ((void *)T.V.lend);
		GMT_free ((void *)T.V.listc);
		GMT_free ((void *)T.V.lptr);
		GMT_free ((void *)xx);
		GMT_free ((void *)yy);
	}
	else {
		GMT_free_table (Table);
	}
	GMT_free ((void *)grid_lon);
	GMT_free ((void *)grid_lat);
	if (!Ctrl->C.active) {
		GMT_free ((void *) lon);
		GMT_free ((void *) lat);
	}
	
	if (n_set > nm) n_set = nm;
	GMT_err_fail (GMT_write_grd (Ctrl->G.file, &h, distgrid, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), Ctrl->G.file);
	GMT_free ((void *)distgrid);
	if (gmtdefs.verbose) fprintf (stderr, "%s: Gridding completed, %ld nodes visited (at least once)\n", GMT_program, n_set);
	
	if (gmtdefs.verbose) fprintf (stderr, "%s: Done!\n", GMT_program);

	Free_sphdistance_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void prepare_polygon (struct GMT_LINE_SEGMENT *P)
{
	GMT_LONG i, quad_no, n_quad;
	GMT_LONG quad[4] = {FALSE, FALSE, FALSE, FALSE};
	double lon_sum = 0.0, lat_sum = 0.0, lon, dlon;
	double xmin1 = 360.0, xmin2 = 360.0, xmax1 = -360.0, xmax2 = -360.0;
	
	P->min[GMT_X] = P->max[GMT_X] = P->coord[GMT_X][0];
	P->min[GMT_Y] = P->max[GMT_Y] = P->coord[GMT_Y][0];
	for (i = 1; i < P->n_rows; i++) {
		lon = P->coord[GMT_X][i];
		while (lon < 0.0) lon += 360.0;	/* Start off with everything in 0-360 range */
		xmin1 = MIN (lon, xmin1);
		xmax1 = MAX (lon, xmax1);
		quad_no = (GMT_LONG)floor (lon/90.0);	/* Yields quadrants 0-3 */
		if (quad_no == 4) quad_no = 0;		/* When lon == 360.0 */
		quad[quad_no] = TRUE;
		while (lon > 180.0) lon -= 360.0;	/* Switch to -180+/180 range */
		xmin2 = MIN (lon, xmin2);
		xmax2 = MAX (lon, xmax2);
		dlon = P->coord[GMT_X][i] - P->coord[GMT_X][i-1];
		if (fabs (dlon) > 180.0) dlon = copysign (360.0 - fabs (dlon), -dlon);
		lon_sum += dlon;
		lat_sum += P->coord[GMT_Y][i];
		if (P->coord[GMT_Y][i] < P->min[GMT_Y]) P->min[GMT_Y] = P->coord[GMT_Y][i];
		if (P->coord[GMT_Y][i] > P->max[GMT_Y]) P->max[GMT_Y] = P->coord[GMT_Y][i];
		if (P->coord[GMT_X][i] < P->min[GMT_X]) P->min[GMT_X] = P->coord[GMT_X][i];
		if (P->coord[GMT_X][i] > P->max[GMT_X]) P->max[GMT_X] = P->coord[GMT_X][i];
	}
	if (GMT_360_RANGE (lon_sum, 0.0)) {	/* Contains a pole */
		if (lat_sum < 0.0) { /* S */
			P->pole = -1;
			P->min[GMT_Y] = -90.0;
		}
		else {	/* N */
			P->pole = +1;
			P->max[GMT_Y] = 90.0;
		}
		P->min[GMT_X] = 0.0;
		P->max[GMT_X] = 360.0;
	}
	else {
		P->pole = 0;
		n_quad = quad[0] + quad[1] + quad[2] + quad[3];	/* How many quadrants had data */
		if (quad[0] && quad[3]) {	/* Longitudes on either side of Greenwhich only, must use -180/+180 notation */
			P->min[GMT_X] = xmin2;
			P->max[GMT_X] = xmax2;
		}
		else if (quad[1] && quad[2]) {	/* Longitudes on either side of the date line, must user 0/360 notation */
			P->min[GMT_X] = xmin1;
			P->max[GMT_X] = xmax1;
		}
		else if (n_quad == 2 && ((quad[0] && quad[2]) || (quad[1] && quad[3]))) {	/* Funny quadrant gap, pick shortest longitude extent */
			if ((xmax1 - xmin1) < (xmax2 - xmin2)) {	/* 0/360 more compact */
				P->min[GMT_X] = xmin1;
				P->max[GMT_X] = xmax1;
			}
			else {						/* -180/+180 more compact */
				P->min[GMT_X] = xmin2;
				P->max[GMT_X] = xmax2;
			}
		}
		else {						/* Either will do, use default settings */
			P->min[GMT_X] = xmin1;
			P->max[GMT_X] = xmax1;
		}
		if (P->min[GMT_X] > P->max[GMT_X]) P->min[GMT_X] -= 360.0;
		if (P->min[GMT_X] < 0.0 && P->max[GMT_X] < 0.0) P->min[GMT_X] += 360.0, P->max[GMT_X] += 360.0;
	}
}

void *New_sphdistance_Ctrl () {	/* Allocate and initialize a new control structure */
	struct SPHDISTANCE_CTRL *C;
	
	C = (struct SPHDISTANCE_CTRL *) GMT_memory (VNULL, 1, sizeof (struct SPHDISTANCE_CTRL), "New_sphdistance_Ctrl");
	C->L.unit = 'E';	/* Default is meter distances */
	return ((void *)C);
}

void Free_sphdistance_Ctrl (struct SPHDISTANCE_CTRL *C) {	/* Deallocate control structure */
	if (C->G.file) free ((void *)C->G.file);	
	if (C->Q.file) free ((void *)C->Q.file);	
	GMT_free ((void *)C);	
}
