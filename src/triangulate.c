/*--------------------------------------------------------------------
 *	$Id: triangulate.c 10173 2014-01-01 09:52:34Z pwessel $
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
 * triangulate reads one or more files (or GMT_stdin) with x,y[,whatever] and
 * outputs the indices of the vertices of the optimal Delaunay triangulation
 * using the method by Watson, D. F., ACORD: Automatic contouring of raw data,
 * Computers & Geosciences, 8, 97-101, 1982.  Optionally, the output may take
 * the form of (1) a multi-segment file with the vertex coordinates needed to
 * draw the triangles, or (2) a grid file based on gridding the plane estimates
 * PS. Instead of Watson's method you may choose to link with the triangulate
 * routine written by Jonathan Shewchuk.  See the file TRIANGLE.HOWTO for
 * details.  That function is far faster than Watson's method.
 *
 * Author:      Paul Wessel
 * Date:        1-JAN-1993
 * Version:     2.0
 * Revised:	13-AUG-1998 for GMT 3.1
 * Revised:	02-MAR-1999 for GMT 3.2
 * Revised:	29-MAR-2000 for GMT 3.3.5
 * Modified:	10 Jul 2000 by PW to add -L
 * Version:	4
 *
 */
 
#include "gmt.h"

struct TRIANGULATE_CTRL {
	struct D {	/* -Dx|y */
		GMT_LONG active;
		GMT_LONG dir;
	} D;
	struct E {	/* -E<Ctrl->E.value> */
		GMT_LONG active;
		double value;
	} E;
	struct F {	/* -F */
		GMT_LONG active;
	} F;
	struct G {	/* -G<output_grdfile> */
		GMT_LONG active;
		char *file;
	} G;
	struct I {	/* -Idx[/dy] */
		GMT_LONG active;
		double xinc, yinc;
	} I;
	struct Q {	/* -Q */
		GMT_LONG active;
	} Q;
};

struct TRIANGULATE_EDGE {
	GMT_LONG begin, end;
};

int main (int argc, char **argv)
{
	int *link, compare_edge(const void *p1, const void *p2);
	GMT_LONG ij, ij1, ij2, ij3, np, n_alloc, n = 0, i, j, k, n_edge;
	GMT_LONG i_min, i_max, j_min, j_max, p, nm;
	GMT_LONG fno, n_files = 0, n_expected_fields, n_args, n_read, n_fields;

	GMT_LONG error = FALSE, map_them = FALSE;
	GMT_LONG nofile = TRUE, done = FALSE, do_grid = FALSE;

	double zj, zk, zl, zlj, zkj, *in;
	double *xx = NULL, *yy = NULL, *zz = VNULL, vx[4], vy[4], xp, yp, a, b, c, f;
	double xkj, xlj, ykj, ylj, out[2], *xe = NULL, *ye = NULL;

	float *grd = NULL;

	char line[BUFSIZ], buffer[BUFSIZ];
#ifdef TRIANGLE_D
	char *algorithm = "Shewchuk";
#else
	char *algorithm = "Watson";
#endif

	FILE *fp = NULL;

	struct GRD_HEADER header;
	struct TRIANGULATE_EDGE *edge = NULL;
	struct TRIANGULATE_CTRL *Ctrl = NULL;

	void *New_triangulate_Ctrl (), Free_triangulate_Ctrl (struct TRIANGULATE_CTRL *C);
	
	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct TRIANGULATE_CTRL *)New_triangulate_Ctrl ();	/* Allocate and initialize a new control structure */

	GMT_grd_init (&header, argc, argv, FALSE);

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
        
				/* Common parameters */
        
				case 'H':
				case 'J':
				case 'M':
				case 'R':
				case 'V':
				case ':':
				case 'b':
				case 'f':
				case 'm':
				case '\0':
					error += (GMT_LONG)GMT_parse_common_options (argv[i], &header.x_min, &header.x_max, &header.y_min, &header.y_max);
					break;
        
				/* Supplemental parameters */
        
				case 'D':
					Ctrl->D.active = TRUE;
					if (argv[i][2] == 'x' || argv[i][2] == 'X')
						Ctrl->D.dir = GMT_X;
					else if (argv[i][2] == 'y' || argv[i][2] == 'Y')
						Ctrl->D.dir = GMT_Y;
					else {
						fprintf (stderr, "%s: GMT SYNTAX ERROR.  Give -Dx or -Dy\n", GMT_program);
						error++;
					}
					break;
				case 'E':
					Ctrl->E.active = TRUE;
					Ctrl->E.value = (argv[i][2] == 'N' || argv[i][2] == 'n') ? GMT_d_NaN : atof (&argv[i][2]);
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
#ifdef TRIANGLE_D
				case 'Q':
					Ctrl->Q.active = TRUE;
					break;
#endif
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
		fprintf (stderr, "triangulate %s - Optimal (Delaunay) triangulation of Cartesian xyz-data [%s]\n\n", GMT_VERSION, algorithm);
		fprintf (stderr, "usage: triangulate <infiles> [-Dx|y] [-E<empty>] [-F] [-G<grdfile>]\n");
		fprintf (stderr, "\t[%s] [%s] [%s]", GMT_H_OPT, GMT_I_OPT, GMT_J_OPT);
#ifdef TRIANGLE_D
		fprintf (stderr, " [-Q]");
#endif
		fprintf (stderr, "\n\t[%s] [-V] [%s] [%s]\n\t[%s] [%s]\n\n", GMT_Rgeo_OPT, GMT_t_OPT, GMT_b_OPT, GMT_f_OPT, GMT_m_OPT);
                
		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);
                
		fprintf (stderr, "\tinfiles (in ASCII) has 2 or more columns.  If no file(s) is given, standard input is read.\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-Dx or -Dy takes derivative in that direction (only with -G) [Default is z value].\n");
		fprintf (stderr, "\t-E value to use for empty nodes [Default is NaN].\n");
		fprintf (stderr, "\t-F Force pixel registration (only with -G) [Default is gridline registration].\n");
		fprintf (stderr, "\t-G Grid data. Give name of output grid file and specify -R -I.\n");
		GMT_explain_option ('H');
		GMT_inc_syntax ('I', 0);
		GMT_explain_option ('J');   
#ifdef TRIANGLE_D
		fprintf (stderr, "\t-Q Compute Voronoi polygon edges instead (requires -m and -R) [Delaunay triangulation].\n");
#endif
		GMT_explain_option ('R');
		GMT_explain_option ('V');
		GMT_explain_option (':');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf (stderr, "\t   Default is 2 input columns.\n");
		fprintf (stderr, "\t-bo writes binary integer index table; ignored if -m is selected [Default is ASCII i/o].\n");
		GMT_explain_option ('f');
		fprintf (stderr, "\t-m output triangle edges as multiple segments separated by a record\n");
		fprintf (stderr, "\t   whose first character is <flag> ['>'].\n");
		fprintf (stderr, "\t   [Default is to output the indices of vertices for each Delaunay triangle].\n");
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

	if (Ctrl->D.active && (Ctrl->D.dir < GMT_X || Ctrl->D.dir > GMT_Y)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -D option.  Must specify x or y\n", GMT_program);
		error = TRUE;
	}
	if (Ctrl->I.active && (Ctrl->I.xinc <= 0.0 || Ctrl->I.yinc <= 0.0)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -I option.  Must specify positive increment(s)\n", GMT_program);
		error = TRUE;
	}
	if (Ctrl->G.active && !Ctrl->G.file) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -G option.  Must specify file name\n", GMT_program);
		error = TRUE;
	}
	if (Ctrl->G.active && (do_grid = (Ctrl->I.active + project_info.region_supplied)) != 2) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Must specify -R, -I, -G for gridding\n", GMT_program);
		error++;
	}
	if (Ctrl->G.active && Ctrl->Q.active) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -G option.  Cannot be used with -Q\n", GMT_program);
		error = TRUE;
	}
	if (Ctrl->Q.active && !GMT_io.multi_segments[GMT_OUT]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -Q option.  Requires -mo\n", GMT_program);
		error = TRUE;
	}
	if (Ctrl->Q.active && !project_info.region_supplied) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -Q option.  Requires -R\n", GMT_program);
		error = TRUE;
	}

	if (error) exit (EXIT_FAILURE);

	if (GMT_io.binary[GMT_IN] && gmtdefs.verbose) {
		char *type[2] = {"double", "single"};
		fprintf (stderr, "%s: Expects %ld-column %s-precision binary data\n", GMT_program, GMT_io.ncol[GMT_IN], type[GMT_io.single_precision[GMT_IN]]);
	}

#ifdef SET_IO_MODE
	GMT_setmode (GMT_OUT);
#endif

	if (do_grid) {
		do_grid = 1;
		header.node_offset = (int)Ctrl->F.active;
		header.x_inc = Ctrl->I.xinc;
		header.y_inc = Ctrl->I.yinc;
		GMT_RI_prepare (&header);	/* Ensure -R -I consistency and set nx, ny */
		GMT_err_fail (GMT_grd_RI_verify (&header, 1), Ctrl->G.file);
	}
	if (project_info.region_supplied && project_info.projection != GMT_NO_PROJ) { /* Gave -R -J */
		map_them = TRUE;
		GMT_err_fail (GMT_map_setup (header.x_min, header.x_max, header.y_min, header.y_max), "");
	}

	/* Now we are ready to take on some input values */

	if (n_files > 0)
		nofile = FALSE;
	else
		n_files = 1;
	n_args = (argc > 1) ? argc : 2;
	n_expected_fields = (GMT_io.ncol[GMT_IN]) ? GMT_io.ncol[GMT_IN] : 2 + do_grid;

	n_alloc = GMT_CHUNK;
	xx = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (double), GMT_program);
	yy = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (double), GMT_program);
	if (do_grid) zz = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (double), GMT_program);

	n = n_read = 0;
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
			n_files--;
			continue;
		}

		if (!nofile && gmtdefs.verbose) fprintf (stderr, "%s: Reading file %s\n", GMT_program, argv[fno]);

		if (GMT_io.io_header[GMT_IN]) for (i = 0; i < GMT_io.n_header_recs; i++) GMT_fgets (line, BUFSIZ, fp);

		while ((n_fields = GMT_input (fp, &n_expected_fields, &in)) >= 0 && !(GMT_io.status & GMT_IO_EOF)) {	/* Not yet EOF */

			n_read++;
			if (GMT_io.status & GMT_IO_MISMATCH) {
				fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld (skipped)\n", GMT_program, n_fields, n_expected_fields, n_read);
				continue;
			}

			xx[n] = in[GMT_X];
			yy[n] = in[GMT_Y];
			if (do_grid) zz[n] = in[GMT_Z];
			n++;

			if (n == n_alloc) {	/* Get more memory */
				n_alloc <<= 1;
				xx = (double *) GMT_memory ((void *)xx, (size_t)n_alloc, sizeof (double), GMT_program);
				yy = (double *) GMT_memory ((void *)yy, (size_t)n_alloc, sizeof (double), GMT_program);
				if (do_grid) zz = (double *) GMT_memory ((void *)zz, (size_t)n_alloc, sizeof (double), GMT_program);
			}
			if (n == INT_MAX) {
				fprintf (stderr, "%s: ERROR: Cannot triangulate more than %d points\n", GMT_program, INT_MAX);
				GMT_free ((void *)xx);
				GMT_free ((void *)yy);
				if (do_grid) GMT_free ((void *)zz);
				exit (EXIT_FAILURE);
			}
		}

		if (fp != GMT_stdin) GMT_fclose (fp);
	}

	if (n_files == 0) {
		fprintf (stderr, "%s: ERROR: No files could be read\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (n_read == 0) {
		fprintf (stderr, "%s: ERROR: No data points read\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	
	xx = (double *) GMT_memory ((void *)xx, (size_t)n, sizeof (double), GMT_program);
	yy = (double *) GMT_memory ((void *)yy, (size_t)n, sizeof (double), GMT_program);
	if (do_grid) zz = (double *) GMT_memory ((void *)zz, (size_t)n, sizeof (double), GMT_program);

	if (map_them) {	/* Must make parallel arrays for projected x/y */

		double *xxp, *yyp;

		xxp = (double *) GMT_memory (VNULL, (size_t)n, sizeof (double), GMT_program);
		yyp = (double *) GMT_memory (VNULL, (size_t)n, sizeof (double), GMT_program);
		for (i = 0; i < n; i++) GMT_geo_to_xy (xx[i], yy[i], &xxp[i], &yyp[i]);

		if (gmtdefs.verbose) fprintf (stderr, "%s: Do Delaunay optimal triangulation on projected coordinates\n", GMT_program);

		if (Ctrl->Q.active) {
			double we[2];
			we[0] = project_info.xmin;	we[1] = project_info.xmax;
			np = GMT_voronoi (xxp, yyp, (int)n, we, &xe, &ye);
		}
		else
			np = GMT_delaunay (xxp, yyp, (int)n, &link);

		GMT_free ((void *)xxp);
		GMT_free ((void *)yyp);
	}
	else {
		if (gmtdefs.verbose) fprintf (stderr, "%s: Do Delaunay optimal triangulation on given coordinates\n", GMT_program);

		if (Ctrl->Q.active) {
			double we[2];
			we[0] = project_info.w;	we[1] = project_info.e;
			np = GMT_voronoi (xx, yy, (int)n, we, &xe, &ye);
		}
		else
			np = GMT_delaunay (xx, yy, (int)n, &link);
	}

	if (gmtdefs.verbose) {
		if (Ctrl->Q.active)
			fprintf (stderr, "%s: %ld Voronoi edges found\n", GMT_program, np);
		else
			fprintf (stderr, "%s: %ld Delaunay triangles found\n", GMT_program, np);
	}
	
	if (do_grid) {

		header.nx = GMT_get_n (header.x_min, header.x_max, header.x_inc, header.node_offset);
		header.ny = GMT_get_n (header.y_min, header.y_max, header.y_inc, header.node_offset);
		header.xy_off = 0.5 * header.node_offset;
		nm = GMT_get_nm (header.nx, header.ny);
		grd = (float *) GMT_memory (VNULL, (size_t)nm, sizeof (float), GMT_program);
		if (!Ctrl->E.active) Ctrl->E.value = GMT_d_NaN;
		for (i = 0; i < nm; i++) grd[i] = (float)Ctrl->E.value;	/* initialize grid */

		for (k = ij = 0; k < np; k++) {

			/* Find equation for the plane as z = ax + by + c */

			vx[0] = vx[3] = xx[link[ij]];	vy[0] = vy[3] = yy[link[ij]];	zj = zz[link[ij++]];
			vx[1] = xx[link[ij]];	vy[1] = yy[link[ij]];	zk = zz[link[ij++]];
			vx[2] = xx[link[ij]];	vy[2] = yy[link[ij]];	zl = zz[link[ij++]];

			xkj = vx[1] - vx[0];	ykj = vy[1] - vy[0];
			zkj = zk - zj;	xlj = vx[2] - vx[0];
			ylj = vy[2] - vy[0];	zlj = zl - zj;

			f = 1.0 / (xkj * ylj - ykj * xlj);
			a = -f * (ykj * zlj - zkj * ylj);
			b = -f * (zkj * xlj - xkj * zlj);
			c = -a * vx[1] - b * vy[1] + zk;

			/* Compute grid indices the current triangle may cover, assuming all triangles are
			   in the -R region (header.x_min/x_max etc.)  Always, i_min <= i_max, j_min <= j_max.
			 */

			xp = MIN (MIN (vx[0], vx[1]), vx[2]);	i_min = GMT_x_to_i (xp, header.x_min, header.x_inc, header.xy_off, header.nx);
			xp = MAX (MAX (vx[0], vx[1]), vx[2]);	i_max = GMT_x_to_i (xp, header.x_min, header.x_inc, header.xy_off, header.nx);
			yp = MAX (MAX (vy[0], vy[1]), vy[2]);	j_min = GMT_y_to_j (yp, header.y_min, header.y_inc, header.xy_off, header.ny);
			yp = MIN (MIN (vy[0], vy[1]), vy[2]);	j_max = GMT_y_to_j (yp, header.y_min, header.y_inc, header.xy_off, header.ny);

			/* Adjustments for triangles outside -R region. */
			/* Triangle to the left or right. */
			if ((i_max < 0) || (i_min >= header.nx)) continue;
			/* Triangle Above or below */
			if ((j_max < 0) || (j_min >= header.ny)) continue;

			/* Triangle covers boundary, left or right. */
			if (i_min < 0) i_min = 0;       if (i_max >= header.nx) i_max = header.nx - 1;
			/* Triangle covers boundary, top or bottom. */
			if (j_min < 0) j_min = 0;       if (j_max >= header.ny) j_max = header.ny - 1;

			for (j = j_min; j <= j_max; j++) {
				yp = GMT_j_to_y (j, header.y_min, header.y_max, header.y_inc, header.xy_off, header.ny);
				p = j * header.nx + i_min;
				for (i = i_min; i <= i_max; i++, p++) {
					xp = GMT_i_to_x (i, header.x_min, header.x_max, header.x_inc, header.xy_off, header.nx);

					if (!GMT_non_zero_winding (xp, yp, vx, vy, (GMT_LONG)4)) continue;	/* Outside */

					if (Ctrl->D.dir == GMT_X)
						grd[p] = (float)a;
					else if (Ctrl->D.dir == GMT_Y)
						grd[p] = (float)b;
					else
						grd[p] = (float)(a * xp + b * yp + c);
				}
			}
		}
		GMT_err_fail (GMT_write_grd (Ctrl->G.file, &header, grd, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), Ctrl->G.file);
		GMT_free ((void *)zz);
		GMT_free ((void *)grd);
	}
	if (GMT_io.multi_segments[GMT_OUT]) {	/* Must find unique edges to output only once */
		if (Ctrl->Q.active) {	/* Voronoi edges */
			for (i = j = 0; i < np; i++) {
				sprintf (buffer, "%c Edge %d\n", GMT_io.EOF_flag[GMT_OUT], (int)i);
				GMT_fputs (buffer, GMT_stdout);
				out[GMT_X] = xe[j];	out[GMT_Y] = ye[j++];
				GMT_output (GMT_stdout, 2, out);
				out[GMT_X] = xe[j];	out[GMT_Y] = ye[j++];
				GMT_output (GMT_stdout, 2, out);
			}
			GMT_free ((void *)xe);
			GMT_free ((void *)ye);
		}
		else {	/* Triangle edges */
			n_edge = 3 * np;
			edge = (struct TRIANGULATE_EDGE *) GMT_memory (VNULL, (size_t)n_edge, sizeof (struct TRIANGULATE_EDGE), GMT_program);
			for (i = ij1 = 0, ij2 = 1, ij3 = 2; i < np; i++, ij1 += 3, ij2 += 3, ij3 += 3) {
				edge[ij1].begin = link[ij1];	edge[ij1].end = link[ij2];
				edge[ij2].begin = link[ij2];	edge[ij2].end = link[ij3];
				edge[ij3].begin = link[ij1];	edge[ij3].end = link[ij3];
			}
			for (i = 0; i < n_edge; i++) if (edge[i].begin > edge[i].end) l_swap (edge[i].begin, edge[i].end);

			qsort ((void *)edge, (size_t)n_edge, sizeof (struct TRIANGULATE_EDGE), compare_edge);
			for (i = 1, j = 0; i < n_edge; i++) {
				if (edge[i].begin != edge[j].begin || edge[i].end != edge[j].end) j++;
				edge[j] = edge[i];
			}
			n_edge = j + 1;

			if (gmtdefs.verbose) fprintf (stderr, "%s: %ld unique triangle edges\n", GMT_program, n_edge);

			for (i = 0; i < n_edge; i++) {
				sprintf (buffer, "%c Edge %ld-%ld\n", GMT_io.EOF_flag[GMT_OUT], edge[i].begin, edge[i].end);
				GMT_fputs (buffer, GMT_stdout);
				out[GMT_X] = xx[edge[i].begin];	out[GMT_Y] = yy[edge[i].begin];
				GMT_output (GMT_stdout, 2, out);
				out[GMT_X] = xx[edge[i].end];	out[GMT_Y] = yy[edge[i].end];
				GMT_output (GMT_stdout, 2, out);
			}
			GMT_free ((void *)edge);
		}
	}
	else if (GMT_io.binary[GMT_OUT])
		i = GMT_fwrite ((void *)link, sizeof (int), (size_t)(3*np), GMT_stdout);
	else
		for (i = ij = 0; i < np; i++, ij += 3) {
			sprintf (buffer, "%d%s%d%s%d\n", link[ij], gmtdefs.field_delimiter, link[ij+1], gmtdefs.field_delimiter, link[ij+2]);
			GMT_fputs (buffer, GMT_stdout);
		}

	GMT_free ((void *) xx);
	GMT_free ((void *) yy);
#ifdef TRIANGLE_D
#ifdef DEBUG
	/* Shewchuk's function allocated the memory separately */
	GMT_memtrack_off (GMT_mem_keeper);
#endif
#endif
	if (!Ctrl->Q.active) GMT_free ((void *) link);
#ifdef TRIANGLE_D
#ifdef DEBUG
	GMT_memtrack_on (GMT_mem_keeper);
#endif
#endif

	if (gmtdefs.verbose) fprintf (stderr, "%s: Done!\n", GMT_program);

	Free_triangulate_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

int compare_edge (const void *p1, const void *p2)
{
	struct TRIANGULATE_EDGE *a, *b;

	a = (struct TRIANGULATE_EDGE *)p1;
	b = (struct TRIANGULATE_EDGE *)p2;
	if (a->begin < b->begin)
		return (-1);
	else if (a->begin > b->begin)
		return (1);
	else {
		if (a->end < b->end)
			return (-1);
		else if (a->end > b->end)
			return (1);
		else
			return (0);
	}
}

void *New_triangulate_Ctrl () {	/* Allocate and initialize a new control structure */
	struct TRIANGULATE_CTRL *C;
	
	C = (struct TRIANGULATE_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct TRIANGULATE_CTRL), "New_triangulate_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	C->D.dir = -1;	/* No derivatives */
	return ((void *)C);
}

void Free_triangulate_Ctrl (struct TRIANGULATE_CTRL *C) {	/* Deallocate control structure */
	if (C->G.file) free ((void *)C->G.file);	
	GMT_free ((void *)C);	
}
