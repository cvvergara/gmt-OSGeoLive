/*--------------------------------------------------------------------
 *	$Id: sphtriangulate.c,v 1.27 2011/07/11 19:22:06 guru Exp $
 *
 *	Copyright (c) 2008-2011 by P. Wessel and W. H. F. Smith
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
 * Spherical triangulation - with Delaunay or Voronoi output.
 * Relies on STRIPACK Fortran F77 library (Renka, 1997). Reference:
 * Renka, R, J,, 1997, Algorithm 772: STRIPACK: Delaunay Triangulation
 *  and Voronoi Diagram on the Surface of a Sphere, AMC Trans. Math.
 *  Software, 23 (3), 416-434.
 * We translate to C using f2c and link with -lf2c
 *
 * Author:      Paul Wessel
 * Date:        13-FEB-2008
 *
 */
 
#include "gmt.h"
#include "sph.h"

struct SPHTRIANGULATE_CTRL {
	struct A {	/* -A */
		GMT_LONG active;
	} A;
	struct C {	/* -C */
		GMT_LONG active;
	} C;
	struct D {	/* -D */
		GMT_LONG active;
	} D;
	struct G {	/* -G<output_grdfile> */
		GMT_LONG active;
		char *file;
	} G;
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
		GMT_LONG mode;	/* 0 is Delaunay, 1 is Voronoi */
	} Q;
	struct T {	/* -T */
		GMT_LONG active;
	} T;
};

void stripack_voronoi_output (FILE *fp, GMT_LONG n, double *lon, double *lat, struct STRIPACK_VORONOI *V, GMT_LONG get_arcs, GMT_LONG get_area, GMT_LONG nodes, char *fname, double scl, PFD distfunc);
void stripack_delaunay_output (FILE *fp, GMT_LONG n, double *lon, double *lat, struct STRIPACK_DELAUNAY *D, GMT_LONG get_arcs, GMT_LONG get_area, GMT_LONG nodes, char *fname, double scl, PFD distfunc);
char *unit_name (char unit, GMT_LONG arc);

int main (int argc, char **argv)
{
	GMT_LONG i, n = 0, n_alloc, n_dup = 0;
	GMT_LONG n_expected_fields, n_args, n_fields, proj_type = 0, fno, n_files = 0;

	GMT_LONG error = FALSE, nofile = TRUE, done = FALSE, first = FALSE;

	double sx, sy, cx, cy, first_x = 0.0, first_y = 0.0, d_scale, *in = NULL;
	double *xx = NULL, *yy = NULL, *zz = NULL, *lon = NULL, *lat = NULL;

	char line[BUFSIZ], buffer[BUFSIZ], *mode[2] = {"Delaunay", "Voronoi"}, *not_used = NULL;
	FILE *fp = NULL;
	PFD distance_func;

	struct SPHTRIANGULATE_CTRL *Ctrl = NULL;
	struct STRIPACK T;
	void *New_sphtriangulate_Ctrl (), Free_sphtriangulate_Ctrl (struct SPHTRIANGULATE_CTRL *C);
	
	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct SPHTRIANGULATE_CTRL *)New_sphtriangulate_Ctrl ();	/* Allocate and initialize a new control structure */

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
        
				/* Common parameters */
        
				case 'H':
				case 'M':
				case 'V':
				case ':':
				case 'b':
				case 'f':
				case 'm':
				case '\0':
					error += GMT_parse_common_options (argv[i], NULL, NULL, NULL, NULL);
					break;
        
				/* Supplemental parameters */
        
				case 'A':
					Ctrl->A.active = TRUE;
					break;
				case 'C':
					Ctrl->C.active = TRUE;
					break;
				case 'D':
					Ctrl->D.active = TRUE;
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
					Ctrl->Q.mode = (argv[i][2] == 'v') ? VORONOI : DELAUNAY;
					break;
				case 'T':
					Ctrl->T.active = TRUE;
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
		fprintf (stderr, "sphtriangulate %s - Delaunay or Voronoi construction of spherical lon,lat data\n\n", GMT_VERSION);
		fprintf (stderr, "usage: sphtriangulate <infiles> [-A] [-C] [-D] [%s] [-L<unit>]\n", GMT_H_OPT);
		fprintf (stderr, "\t[-N<nodefile>] [-Qd|v] [-T] [-V] [%s] [%s] [%s]\n\n", GMT_t_OPT, GMT_b_OPT, GMT_m_OPT);
                
		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);
                
		fprintf (stderr, "\tinfiles (in ASCII) has 2 or more columns.  If no file(s) is given, standard input is read.\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-A Compute and print triangle or polygon areas in header records (see -L for units)\n");
		fprintf (stderr, "\t   If -T is selected we print arc lengths instead.\n");
		fprintf (stderr, "\t   Cannot be used with the binary output option\n");
		fprintf (stderr, "\t-C Conserve memory (Converts lon/lat <--> x/y/z when needed) [store both in memory]\n");
		fprintf (stderr, "\t-D Used with -M to skip repeated input vertex at the end of a closed segment\n");
		GMT_explain_option ('H');
		fprintf (stderr, "\t-L Specify distance unit as m(e)ter, (k)m, (m)ile, (n)autical mile, or (d)egree.\n");
		fprintf (stderr, "\t   Calculations uses spherical approximations.  Default unit is meters.\n");
		fprintf (stderr, "\t   Set ELLIPSOID to WGS-84 to get geodesic distances.\n");
		fprintf (stderr, "\t-N Output file for Delaunay or Voronoi polygon information [Store in output segment headers]\n");
		fprintf (stderr, "\t   Delaunay: output is the node triplets and area (i, j, k, area)\n");
		fprintf (stderr, "\t   Voronoi: output is the node coordinates and polygon area (lon, lat, area).\n");
		fprintf (stderr, "\t   Cannot be used with -T.\n");
		fprintf (stderr, "\t-Q Append d for Delaunay triangles or v for Voronoi polygons [Delaunay]\n");
		fprintf (stderr, "\t   If -bo is used then -N may be used to specify a separate file where the polygon\n");
		fprintf (stderr, "\t   information normally is written.\n");
		fprintf (stderr, "\t-T Write arcs [Default writes polygons]\n");
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

	if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
    if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = 2;
    if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] < 2) {
            fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data (-bi) must have at least 2 columns\n", GMT_program);
		error++;
	}
	if (GMT_io.binary[GMT_OUT] && Ctrl->A.active && !Ctrl->N.active) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary output does not support storing areas unless -N is used\n", GMT_program);
		error++;
	}
	if (Ctrl->N.active && !Ctrl->N.file) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -N:  Must specify output file\n", GMT_program);
		error = TRUE;
	}
	if (Ctrl->D.active && !GMT_io.multi_segments[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  -D requires the -m[i] option\n", GMT_program);
		error++;
	}
	error += GMT_get_dist_scale (Ctrl->L.unit, &d_scale, &proj_type, &distance_func);
	if (error) exit (EXIT_FAILURE);
	if (Ctrl->L.unit == 'D') d_scale = -d_scale;	/* Flag so we can do steradians */
	if (GMT_io.binary[GMT_IN] && gmtdefs.verbose) {
		char *type[2] = {"double", "single"};
		fprintf (stderr, "%s: Expects %ld-column %s-precision binary data\n", GMT_program, GMT_io.ncol[GMT_IN], type[GMT_io.single_precision[GMT_IN]]);
	}

#ifdef SET_IO_MODE
	GMT_setmode (GMT_OUT);
#endif

	/* Now we are ready to take on some input values */

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
				first_x = in[GMT_X];	first_y = in[GMT_Y];
				first = TRUE;
			}
			if (Ctrl->D.active && !first) {	/* Look for duplicate point at end of segments that replicate start point */
				if (in[GMT_X] == first_x && in[GMT_Y] == first_y) {	/* If any point after the first matches the first */
					n_dup++;
					continue;
				}
				else if (in[GMT_X] == xx[n] && in[GMT_Y] == yy[n]) {	/* Duplicate of previous point, skip*/
					n_dup++;
					continue;
				}
			}
			sincosd (in[GMT_Y], &sy, &cy);
			sincosd (in[GMT_X], &sx, &cx);
			xx[n] = cy * cx;
			yy[n] = cy * sx;
			zz[n] = sy;
			if (!Ctrl->C.active) {
				lon[n] = in[GMT_X];
				lat[n] = in[GMT_Y];
			}
			n++;

			if (n == n_alloc) {	/* Get more memory */
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
	T.mode = (int)Ctrl->Q.mode;
	stripack_lists (n, xx, yy, zz, gmtdefs.verbose, &T);	/* Do the basic triangulation */
	if (Ctrl->C.active) {	/* Recompute lon,lat and set pointers */
		cart_to_geo (n, xx, yy, zz, xx, yy);	/* Revert to lon, lat */
		lon = xx;
		lat = yy;
	}
	GMT_free ((void *) zz);
	if (Ctrl->T.active) GMT_io.multi_segments[GMT_OUT] = TRUE;	/* Must produce multisegment output files */
	if (!GMT_io.binary[GMT_OUT]) {
		sprintf (buffer, "# sphtriangulate %s output via STRPACK.", mode[Ctrl->Q.mode]);
		GMT_fputs (buffer, GMT_stdout);
		if (Ctrl->T.active) {
			sprintf (buffer, "  Arc lengths in %s\n", unit_name (Ctrl->L.unit, TRUE));
			GMT_fputs (buffer, GMT_stdout);
		}
		else {
			sprintf (buffer, "  Areas in %s\n", unit_name (Ctrl->L.unit, FALSE));
			GMT_fputs (buffer, GMT_stdout);
		}
	}
	if (Ctrl->Q.mode == VORONOI) {
		stripack_voronoi_output (GMT_stdout, n, lon, lat, &T.V, Ctrl->T.active, Ctrl->A.active, Ctrl->N.active, Ctrl->N.file, d_scale, distance_func);
		GMT_free ((void *)T.V.lon);
		GMT_free ((void *)T.V.lat);
		GMT_free ((void *)T.V.lend);
		GMT_free ((void *)T.V.listc);
		GMT_free ((void *)T.V.lptr);
	}
	else {
		stripack_delaunay_output (GMT_stdout, n, lon, lat, &T.D, Ctrl->T.active, Ctrl->A.active, Ctrl->N.active, Ctrl->N.file, d_scale, distance_func);
	}
	
	GMT_free ((void *)T.D.tri);
	if (!Ctrl->C.active) {
		GMT_free ((void *) lon);
		GMT_free ((void *) lat);
	}
	GMT_free ((void *) xx);
	GMT_free ((void *) yy);

	if (gmtdefs.verbose) fprintf (stderr, "%s: Done!\n", GMT_program);

	Free_sphtriangulate_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void stripack_delaunay_output (FILE *fp, GMT_LONG n, double *lon, double *lat, struct STRIPACK_DELAUNAY *D, GMT_LONG get_arcs, GMT_LONG get_area, GMT_LONG nodes, char *fname, double scl, PFD distfunc)
{	/* Prints out the Delaunay triangles either as polygons (for filling) or arcs (lines). */
	GMT_LONG i, ij, k;
	double out[4], area_sphere = 0.0, area_triangle = GMT_d_NaN, V[3][3], R2 = 6371007.1810 * 6371007.1810, dist = GMT_d_NaN;
	FILE *fp2 = NULL;

	if (nodes && !get_arcs) {	/* Voronoi node and area information will be written to voronoi_nodes.b */
		if ( (fp2 = GMT_fopen (fname, GMT_io.w_mode)) == NULL) {
			fprintf(stderr,"%s:  Cannot open w %s\n", GMT_program, fname);
			exit (EXIT_FAILURE);
		}
	}

	if (scl < 0.0) /* Steradians */
		R2 = scl = 1.0;
	else
		R2 *= (scl * scl);	/* Get final measure unit for area */
	if (!get_arcs) {	/* Output polygons of three points.  Not explicitly closed (use -L or -G in psxy) */
		if (gmtdefs.verbose) fprintf (stderr, "%s: Output %d unique triangle polygons\n", GMT_program, D->n);
		for (k = ij = 0; k < D->n; k++, ij += TRI_NROW) {	/* For each triangle */
			/* Write multisegment header with triangle # and the three node numbers */
			if (get_area) {
				for (i = 0; i < 3; i++) geo_to_cart (lon[D->tri[ij+i]], lat[D->tri[ij+i]], V[i], TRUE);
				area_triangle = stripack_areas (V[0], V[1], V[2]);
				area_sphere += area_triangle;
			}
			if (nodes) {
				out[GMT_X] = D->tri[ij];
				out[GMT_Y] = D->tri[ij+1];
				out[GMT_Z] = D->tri[ij+2];
				out[3] = area_triangle * R2;
				GMT_output (fp2, 4, out);
				if (!GMT_io.binary[GMT_OUT]) sprintf (GMT_io.segment_header, "%c Triangle: %ld %d-%d-%d Area: %g\n", GMT_io.EOF_flag[GMT_OUT], k,
					D->tri[ij], D->tri[ij+1], D->tri[ij+2], area_triangle * R2);
			}
			else if (!GMT_io.binary[GMT_OUT]) {
				sprintf (GMT_io.segment_header, "%c Triangle: %ld\n", GMT_io.EOF_flag[GMT_OUT], k);
			}
			GMT_write_segmentheader (fp, 2);
			for (i = 0; i < 3; i++) {	/* Write out the three vertices */
				out[GMT_X] = lon[D->tri[ij+i]];
				out[GMT_Y] = lat[D->tri[ij+i]];
				GMT_output (fp, 2, out);
			}
		}
		if (gmtdefs.verbose && get_area) fprintf (stderr, "%s: Total surface area = %g\n", GMT_program, area_sphere * R2);
		if (nodes) GMT_fclose (fp2);
	}
	else {	/* Want just the arcs (to draw then, probably).  This avoids repeating arcs */
		GMT_LONG j, ij1, ij2, ij3, n_arcs;
		struct STRPACK_ARC *arc;
		
		n_arcs = 3 * D->n;
		arc = (struct STRPACK_ARC *) GMT_memory (VNULL, (size_t)n_arcs, sizeof (struct STRPACK_ARC), GMT_program);
		for (k = ij = ij1 = 0, ij2 = 1, ij3 = 2; k < D->n; k++, ij1 += TRI_NROW, ij2 += TRI_NROW, ij3 += TRI_NROW) {	/* For each triangle */
			arc[ij].begin = D->tri[ij1];	arc[ij].end = D->tri[ij2];	ij++;
			arc[ij].begin = D->tri[ij2];	arc[ij].end = D->tri[ij3];	ij++;
			arc[ij].begin = D->tri[ij1];	arc[ij].end = D->tri[ij3];	ij++;
		}
		for (k = 0; k < n_arcs; k++) if (arc[k].begin > arc[k].end) i_swap (arc[k].begin, arc[k].end);

		qsort ((void *)arc, (size_t)n_arcs, sizeof (struct STRPACK_ARC), compare_arc);
		for (i = 1, j = 0; i < n_arcs; i++) {
			if (arc[i].begin != arc[j].begin || arc[i].end != arc[j].end) j++;
			arc[j] = arc[i];
		}
		n_arcs = j + 1;
		if (gmtdefs.verbose) fprintf (stderr, "%s: Output %ld unique triangle arcs\n", GMT_program, n_arcs);

		for (i = 0; i < n_arcs; i++) {
			out[GMT_X] = lon[arc[i].begin];	out[GMT_Y] = lat[arc[i].begin];
			out[2] = lon[arc[i].end];	out[3] = lat[arc[i].end];
			if (!GMT_io.binary[GMT_OUT]) {
				if (get_area) dist = scl * (*distfunc)(out[GMT_X], out[GMT_Y], out[2], out[3]);
				sprintf (GMT_io.segment_header, "%c Arc: %d-%d Length: %g\n", GMT_io.EOF_flag[GMT_OUT], arc[i].begin, arc[i].end, dist);
			}
			GMT_write_segmentheader (fp, 2);
			GMT_output (fp, 2, out);
			GMT_output (fp, 2, &out[2]);
		}
		GMT_free ((void *)arc);
	}
}

void stripack_voronoi_output (FILE *fp, GMT_LONG n, double *lon, double *lat, struct STRIPACK_VORONOI *V, GMT_LONG get_arcs, GMT_LONG get_area, GMT_LONG nodes, char *fname, double scl, PFD distfunc)
{	/* Prints out the Voronoi polygons either as polygons (for filling) or arcs (lines) */
	GMT_LONG i, j, k, node, vertex, node_stop, node_new, vertex_new, node_last, vertex_last, n_arcs = 0;
	GMT_LONG n_alloc = GMT_CHUNK, p_alloc = GMT_TINY_CHUNK;
	double area_sphere = 0.0, area_polygon, area_triangle, area_km2 = GMT_d_NaN, dist = GMT_d_NaN, V1[3], V2[3], V3[3], out[4];
	double *plat, *plon, R2 = 6371007.1810 * 6371007.1810;
	struct STRPACK_ARC *arc = NULL;
	FILE *fp2 = NULL;

	if (nodes && !get_arcs) {	/* Voronoi node and area information will be written to separate node file */
		if ( (fp2 = GMT_fopen (fname, GMT_io.w_mode)) == NULL) {
			fprintf(stderr,"%s:  Cannot open w %s\n", GMT_program, fname);
			exit (EXIT_FAILURE);
		}
	}
		
	if (scl < 0.0) /* Steradians */
		R2 = scl = 1.0;
	else
		R2 *= (scl * scl);	/* Get final measure unit for area */
	if (get_arcs) arc = (struct STRPACK_ARC *) GMT_memory (VNULL, n_alloc, sizeof (struct STRPACK_ARC), GMT_program);
	plon = (double *) GMT_memory (VNULL, p_alloc, sizeof (double), GMT_program);
	plat = (double *) GMT_memory (VNULL, p_alloc, sizeof (double), GMT_program);
	
	for (node = 0; node < n; node++) {

		area_polygon = 0.0;

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

			if (get_arcs) {
				arc[n_arcs].begin = (int)vertex_last;	arc[n_arcs].end = (int)vertex_new;
				n_arcs++;
				if (n_arcs == n_alloc) {
					n_alloc <<= 1;
					arc = (struct STRPACK_ARC *) GMT_memory ((void *)arc, n_alloc, sizeof (struct STRPACK_ARC), GMT_program);
				}
				vertex++;
			}
			else {
				plon[vertex] = V->lon[vertex_last];
				plat[vertex] = V->lat[vertex_last];
				if (get_area) {	/* Convert three corners to Cartesian */
					geo_to_cart (V->lon[node], V->lat[node], V1, TRUE);
					geo_to_cart (V->lon[vertex_last], V->lat[vertex_last], V2, TRUE);
					geo_to_cart (V->lon[vertex_new], V->lat[vertex_new], V3, TRUE);
					area_triangle = stripack_areas (V1, V2, V3);
					area_polygon += area_triangle;
				}
				vertex++;
				if (vertex == p_alloc) {
					p_alloc <<= 1;
					plon = (double *) GMT_memory ((void *)plon, p_alloc, sizeof (double), GMT_program);
					plat = (double *) GMT_memory ((void *)plat, p_alloc, sizeof (double), GMT_program);
				}
			}

			/* When we reach the vertex where we started, we are done with this polygon */
		} while (node_new != node_stop);

		if (!get_arcs) {
			if (get_area) area_km2 = area_polygon * R2;
			if (nodes) {
				out[GMT_X] = lon[node];
				out[GMT_Y] = lat[node];
				out[2] = area_km2;
				GMT_output (fp2, 3, out);
				if (!GMT_io.binary[GMT_OUT]) sprintf (GMT_io.segment_header, "%c Pol: %ld\n", GMT_io.EOF_flag[GMT_OUT], node);
			}
			else if (!GMT_io.binary[GMT_OUT]) {	/* Only ASCII output can have detailed header information */
				sprintf (GMT_io.segment_header, "%c Pol: %ld %g %g Area: %g\n", GMT_io.EOF_flag[GMT_OUT], node, lon[node], lat[node], area_km2);
			}
			
			GMT_write_segmentheader (fp, 2);
			for (i = 0; i < vertex; i++) {
				out[GMT_X] = plon[i];
				out[GMT_Y] = plat[i];
				GMT_output (fp, 2, out);
			}
	
			if (get_area) area_sphere += area_polygon;
		}
	}
	if (get_arcs) {	/* Process arcs */
		for (k = 0; k < n_arcs; k++) if (arc[k].begin > arc[k].end) i_swap (arc[k].begin, arc[k].end);

		qsort ((void *)arc, (size_t)n_arcs, sizeof (struct STRPACK_ARC), compare_arc);
		for (i = 1, j = 0; i < n_arcs; i++) {
			if (arc[i].begin != arc[j].begin || arc[i].end != arc[j].end) j++;
			arc[j] = arc[i];
		}
		n_arcs = j + 1;
		if (gmtdefs.verbose) fprintf (stderr, "%s: Output %ld unique Voronoi arcs\n", GMT_program, n_arcs);

		for (i = 0; i < n_arcs; i++) {
			out[GMT_X] = V->lon[arc[i].end];	out[GMT_Y] = V->lat[arc[i].end];
			out[2] = V->lon[arc[i].begin];	out[3] = V->lat[arc[i].begin];
			if (!GMT_io.binary[GMT_OUT]) {
				if (get_area) dist = scl * (*distfunc)(out[GMT_X], out[GMT_Y], out[2], out[3]);
				sprintf (GMT_io.segment_header, "%c Arc: %d-%d Length: %g\n", GMT_io.EOF_flag[GMT_OUT], arc[i].begin, arc[i].end, dist);
			}
			GMT_write_segmentheader (fp, 2);
			GMT_output (fp, 2, out);
			GMT_output (fp, 2, &out[2]);
		}
		GMT_free ((void *)arc);
	}
	else {
		GMT_free ((void *)plon);
		GMT_free ((void *)plat);
		if (nodes) GMT_fclose (fp2);
		if (gmtdefs.verbose && get_area)
			fprintf (stderr, "%s: Total surface area = %g\n", GMT_program, area_sphere * R2);
	}
}

char *unit_name (char unit, GMT_LONG arc) {
	char *name;
	switch (unit) {
		case 'K':	/* km */
			name = (arc) ? "km" : "km^2";
			break;
		case 'M':	/* Miles */
			name = (arc) ? "miles" : "sq. miles";
			break;
		case 'N':	/* Nautical miles */
			name = (arc) ? "nautical miles" : "sq. nautical miles";
			break;
		case 'D':	/* Degrees */
			name = (arc) ? "degrees" : "steradians";
			break;
		default:
			name = (arc) ? "m" : "m^2";
			break;
	}
	return (name);
}
void *New_sphtriangulate_Ctrl () {	/* Allocate and initialize a new control structure */
	struct SPHTRIANGULATE_CTRL *C;
	
	C = (struct SPHTRIANGULATE_CTRL *) GMT_memory (VNULL, 1, sizeof (struct SPHTRIANGULATE_CTRL), "New_sphtriangulate_Ctrl");
	C->L.unit = 'E';	/* Default is meter distances */
	
	return ((void *)C);
}

void Free_sphtriangulate_Ctrl (struct SPHTRIANGULATE_CTRL *C) {	/* Deallocate control structure */
	if (C->G.file) free ((void *)C->G.file);	
	if (C->N.file) free ((void *)C->N.file);	
	GMT_free ((void *)C);	
}
