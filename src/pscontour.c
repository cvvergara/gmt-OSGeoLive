/*--------------------------------------------------------------------
 *	$Id: pscontour.c,v 1.126 2011/07/08 21:27:06 guru Exp $
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
/* pscontour will read a file of points in the plane, performs the
 * Delaunay triangulation, and contours these triangles.  As an option
 * the user may provide a file with indices of which vertices constitute
 * the triangles.
 *
 * Author:	Paul Wessel
 * Date:	13-SEP-2001
 * Version:	4
 */

#include "gmt.h"
#include "pslib.h"

struct PSCONTOUR_CTRL {
	struct GMT_CONTOUR contour;
	struct A {	/* -A[-|aint][+a<angle>][+c<dx>[/<dy>]][+f<font>][+g<fill>][+j<just>][+l<label>][+o|O|t][+s<size>][+p<pen>][+u<unit>] */
		GMT_LONG active;
		GMT_LONG off;
		double interval;
	} A;
	struct C {	/* -C<cpt> */
		GMT_LONG active;
		char *file;
	} C;
	struct D {	/* -D<dumpfile> */
		GMT_LONG active;
		char *file;
	} D;
	struct E {	/* -Eazim/elev */
		GMT_LONG active;
		double azimuth, elevation;
	} E;
	struct G {	/* -G[d|f|n|l|L|x|X]<params> */
		GMT_LONG active;
	} G;
	struct I {	/* -I */
		GMT_LONG active;
	} I;
	struct L {	/* -L<pen> */
		GMT_LONG active;
		struct GMT_PEN pen;
	} L;
	struct N {	/* -N */
		GMT_LONG active;
	} N;
	struct S {	/* -S */
		GMT_LONG active;
	} S;
	struct T {	/* -T<indexfile> */
		GMT_LONG active;
		char *file;
	} T;
	struct W {	/* -W[+]<type><pen> */
		GMT_LONG active;
		GMT_LONG color;
		struct GMT_PEN pen;
	} W;
};

/* Returns the id of the node common to the two edges */
#define get_node_index(edge_1, edge_2) (((PSCONTOUR_SUM = (edge_1) + (edge_2)) == 1) ? 1 : ((PSCONTOUR_SUM == 2) ? 0 : 2))
#define get_other_node(node1, node2) (((PSCONTOUR_SUM = (node1 + node2)) == 3) ? 0 : ((PSCONTOUR_SUM == 2) ? 1 : 2))	/* The other node needed */

struct PSCONTOUR_LINE {	/* Beginning and end of straight contour segment */
	double x0, y0;
	double x1, y1;
};

struct PSCONTOUR {
	GMT_LONG n_alloc, nl;
	double val;
	double angle;
	char type;
	struct PSCONTOUR_LINE *L;
};

struct PSCONTOUR_PT {
	double x, y;
	struct PSCONTOUR_PT *next;
};

struct PSCONTOUR_CHAIN {
	struct PSCONTOUR_PT *begin;
	struct PSCONTOUR_PT *end;
	struct PSCONTOUR_CHAIN *next;
};

int main(int argc, char **argv)
{
	GMT_LONG nx, k2, k3, node1, node2, c, PSCONTOUR_SUM;
	GMT_LONG n_alloc, section = 0, *vert = NULL, *cind = NULL, n_read, n_fields, n_contours = 0;
	GMT_LONG add, bad, n_expected_fields, last_entry, last_exit;
	GMT_LONG ij, n, np, k, i, low, high;
	
	int rgb[3], *ind = NULL;
	
	GMT_LONG error = FALSE, more, skip = FALSE, advance = FALSE, close;

	double xx[3], yy[3], zz[3], xout[5], yout[5], *in, west, east, south, north;
	double *xc = NULL, *yc = NULL, *zc = NULL, *x = NULL, *y = NULL, *z = NULL, *xp = NULL, *yp = NULL, current_contour = -DBL_MAX, xyz[2][3];

	char line[BUFSIZ], cont_label[GMT_LONG_TEXT], format[GMT_LONG_TEXT], *not_used = NULL;
#ifdef TRIANGLE_D
	char *algorithm = "Shewchuk";
#else
	char *algorithm = "Watson";
#endif

	FILE *fp = NULL, *fp_d = NULL;

	struct PSCONTOUR *cont = NULL;
	struct PSCONTOUR_CTRL *Ctrl = NULL;

	GMT_LONG get_triangle_crossings (double *x, double *y, double *z, int *ind, double **xc, double **yc, double **zc, GMT_LONG **v, GMT_LONG **cindex);
	void paint_it (double x[], double y[], GMT_LONG n, double z);
	void *New_pscontour_Ctrl (), Free_pscontour_Ctrl (struct PSCONTOUR_CTRL *C);

	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct PSCONTOUR_CTRL *)New_pscontour_Ctrl ();	/* Allocate and initialize a new control structure */

	/* Check and interpret the command line arguments */

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch(argv[i][1]) {

				/* Common parameters */

				case 'B':
				case 'H':
				case 'J':
				case 'K':
				case 'M':
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
				case 'm':
				case ':':
				case '\0':
					error += GMT_parse_common_options (argv[i], &west, &east, &south, &north);
					break;

				/* Supplemental parameters */

				case 'A':
					/* Format can be one of two:
					 * 3.4.x: -A[-][f<fontsize>][a<angle>][/<r/g/b>][t|o]
					 * 4.x:   -A[-][+a<angle>][+c<dx>[/<dy>]][+f<size>][+g[<fill>]][+j<just>][+o][+p[<pen>]][+v]
					 */
					Ctrl->A.active = TRUE;
					if (argv[i][2] == '-') Ctrl->A.off = TRUE;
					bad = GMT_contlabel_specs (&argv[i][2], &Ctrl->contour);
					Ctrl->contour.annot = !Ctrl->A.off;
					if (bad) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -A option.  Correct syntax:\n", GMT_program);
						fprintf (stderr, "\t-A[-][+a<angle>][+c<dx>[/<dy>]][+f<font>][+g[<fill>]][+j<just>][+o][+p[<pen>]][+s<size>][+u<unit>][+v]\n");
						error += bad;
					}
					break;

				case 'C':
					Ctrl->C.active = TRUE;
					Ctrl->C.file = strdup (&argv[i][2]);
					break;
				case 'D':
					Ctrl->D.active = TRUE;
					free ((void *)Ctrl->D.file);
					Ctrl->D.file = strdup (&argv[i][2]);
					break;
				case 'E':
					Ctrl->E.active = TRUE;
					error += GMT_get_proj3D (&argv[i][2], &Ctrl->E.azimuth, &Ctrl->E.elevation);
					break;
				case 'G':
					Ctrl->G.active = TRUE;
					error += GMT_contlabel_info ('G', &argv[i][2], &Ctrl->contour);
					break;
				case 'I':
					Ctrl->I.active = TRUE;
					break;
				case 'L':
					Ctrl->L.active = TRUE;
					if (GMT_getpen (&argv[i][2], &Ctrl->L.pen)) {
						GMT_pen_syntax ('L', " ");
						error++;
					}
					break;
				case 'N':
					Ctrl->N.active = TRUE;
					break;
				case 'S':		/* Skip points outside border */
					Ctrl->S.active = TRUE;
					break;
				case 'T':
					Ctrl->T.file = strdup (&argv[i][2]);
					Ctrl->T.active = TRUE;
					break;
				case 'W':
					Ctrl->W.active = TRUE;
					k = 2;
					if (argv[i][k] == '+') Ctrl->W.color = TRUE, k++;
					if (argv[i][k] && GMT_getpen (&argv[i][k], &Ctrl->W.pen)) {
						GMT_pen_syntax ('W', " ");
						error++;
					}
					break;

				/* Options not recognized */

				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else
			fp = GMT_fopen (argv[i], GMT_io.r_mode);
	}

	if (GMT_give_synopsis_and_exit || argc == 1) {	/* Display usage */
		fprintf (stderr,"pscontour %s - Contour xyz-data by triangulation [%s]\n\n", GMT_VERSION, algorithm);
		fprintf (stderr,"usage: pscontour <xyzfile> -C<cpt_file> %s %s\n", GMT_J_OPT, GMT_Rgeoz_OPT);
		fprintf (stderr, "\t[-A[-|<annot_int>][<labelinfo>] [%s] [-D<dumpfile>] [%s]\n", GMT_B_OPT, GMT_E_OPT);
		fprintf (stderr, "\t[%s] [%s] [-I] [%s] [-K] [-L<pen>] [-N] [-O]\n", GMT_CONTG, GMT_H_OPT, GMT_Jz_OPT);
 		fprintf (stderr, "\t[-P] [-S] [-T<indexfile>] [-U] [-V] [-W[+]<pen>] [%s] [%s]\n", GMT_X_OPT, GMT_Y_OPT);
		fprintf (stderr, "\t[%s] [%s] [%s] [%s]\n\n", GMT_c_OPT, GMT_t_OPT, GMT_b_OPT, GMT_mo_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr,"\t-C Color palette table\n");
		GMT_explain_option ('j');
		GMT_explain_option ('Z');
		GMT_explain_option ('R');
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-A Annotation label information.\n");
		fprintf (stderr, "\t   Give A- to disable all contour annotations implied in -C\n");
		fprintf (stderr, "\t   <labelinfo> controls the specifics of the labels.  Append what you need:\n");
		GMT_label_syntax (5, 0);
		GMT_explain_option ('b');
		fprintf (stderr, "\t-D to Dump contour lines to individual files (but see -m)\n");
		GMT_explain_option ('E');
		fprintf (stderr, "\t-G Controls placement of labels along contours.  Choose among five algorithms:\n");
		GMT_cont_syntax (3, 0);
		GMT_explain_option ('H');
		fprintf (stderr, "\t-I Color triangles using the cpt file\n");
		GMT_explain_option ('K');
		GMT_pen_syntax ('L', "draws the triangular mesh with the specified pen.");
		fprintf (stderr,"\t-N do NOT clip contours/image at the border [Default clips]\n");
		GMT_explain_option ('O');
		GMT_explain_option ('P');
		fprintf (stderr,"\t-S Skip xyz points outside region [Default keeps all]\n");
		fprintf (stderr,"\t-T file with triplets of point indices for each triangle\n");
		fprintf (stderr,"\t   [Default performs the Delaunay triangulation on xyz-data]\n");
		GMT_explain_option ('U');
		GMT_explain_option ('V');
		GMT_pen_syntax ('W', "selects contouring and sets contour pen attributes.");
		fprintf (stderr, "\t   Use + to draw colored contours based on the cpt file\n");
		GMT_explain_option ('X');
		GMT_explain_option ('c');
		GMT_explain_option (':');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf (stderr,"\t   Default is 3 input columns.\n");
		GMT_explain_option ('o');
		GMT_explain_option ('n');
		fprintf (stderr, "\t-m Used with -D.   Create a single multiple segment file where contours are separated by a record\n");
		fprintf (stderr, "\t   whose first character is <flag> ['>'].  This header also has the contour level value\n");
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	/* Check that the options selected are mutually consistent */

	if (!project_info.region_supplied) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify -R option\n", GMT_program);
		error++;
	}
	if (!(Ctrl->W.active || Ctrl->I.active)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify one of -W or -I\n", GMT_program);
		error++;
	}
	if (!Ctrl->C.file) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -C option:  Must specify a color palette table\n", GMT_program);
		error++;
	}
	if (Ctrl->T.active && !Ctrl->T.file) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -T option:  Must specify an index file\n", GMT_program);
		error++;
	}
	if (Ctrl->E.elevation <= 0.0 || Ctrl->E.elevation > 90.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -E option:  Elevation must be in 0-90 range\n", GMT_program);
		error++;
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
        if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = 3;
        if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] < 3) {
                fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data (-bi) must have at least 3 columns\n", GMT_program);
		error++;
	}

	if (error) exit (EXIT_FAILURE);

	if (GMT_io.binary[GMT_IN] && gmtdefs.verbose) {
		char *type[2] = {"double", "single"};
		fprintf (stderr, "%s: Expects %ld-column %s-precision binary data\n", GMT_program, GMT_io.ncol[GMT_IN], type[GMT_io.single_precision[GMT_IN]]);
	}

	if (GMT_io.binary[GMT_OUT] && GMT_io.multi_segments[GMT_OUT]) {
		fprintf (stderr, "%s: GMT Warning.  -m for output ignored with -D, -bo\n", GMT_program);
		error++;
	}
	z_project.view_azimuth = Ctrl->E.azimuth;
	z_project.view_elevation = Ctrl->E.elevation;
	if (Ctrl->T.active) {
		if ((fp_d = GMT_fopen (Ctrl->T.file, GMT_io.r_mode)) == NULL) {
		fprintf (stderr, "%s: Could not open index file %s\n", GMT_program, Ctrl->T.file);
			exit (EXIT_FAILURE);
		}
	}

	if (Ctrl->D.active && Ctrl->D.file[0] == 0) {
		fprintf (stderr, "%s: contours will be written to file contour\n", GMT_program);
		strcpy (Ctrl->D.file,"contour");
	}
	if (GMT_io.multi_segments[GMT_OUT]) GMT_io.multi_segments[GMT_OUT] = 2;

	GMT_read_cpt (Ctrl->C.file);
	if (Ctrl->I.active && GMT_continuous) {
		fprintf (stderr, "%s: -I option requires constant color between contours!\n", GMT_program);
		exit (EXIT_FAILURE);
	}
#ifdef GMT_CPT2	
	if (GMT_categorical) {
		fprintf (stderr, "%s: GMT WARNING:  Categorical data (as implied by CPT file) does not have contours.  Check plot.\n", GMT_program);
		exit (EXIT_FAILURE);
	}
#endif
	GMT_err_fail (GMT_map_setup (west, east, south, north), "");

	GMT_plotinit (argc, argv);

	if (project_info.three_D) ps_transrotate (-z_project.xmin, -z_project.ymin, 0.0);

        if (!Ctrl->N.active) GMT_map_clip_on (GMT_no_rgb, 3);

	if (fp == NULL) {
		fp = GMT_stdin;
#ifdef SET_IO_MODE
		GMT_setmode (GMT_IN);
#endif
	}

	n_alloc = GMT_CHUNK;
	x = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (double), GMT_program);
	y = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (double), GMT_program);
	z = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (double), GMT_program);

	if (GMT_io.io_header[GMT_IN]) for (i = 0; i < GMT_io.n_header_recs; i++) not_used = GMT_fgets (line, BUFSIZ, fp);

	xyz[0][2] = DBL_MAX;	xyz[1][2] = -DBL_MAX;
	n_read = n = 0;
	n_expected_fields = (GMT_io.ncol[GMT_IN]) ? GMT_io.ncol[GMT_IN] : 3;
	n_fields = GMT_input (fp, &n_expected_fields, &in);

	while (! (GMT_io.status & GMT_IO_EOF)) {	/* Not yet EOF */

		n_read++;
		if (GMT_io.status & GMT_IO_MISMATCH) {
			fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld (skipped)\n", GMT_program, n_fields, n_expected_fields, n_read);
			continue;
		}

		if (Ctrl->S.active) {	/* Must check if points are inside plot region */
			GMT_map_outside (in[GMT_X], in[GMT_Y]);
			skip = (GMT_abs (GMT_x_status_new) > 1 || GMT_abs (GMT_y_status_new) > 1);
		}

		if (!(skip || GMT_is_dnan (in[GMT_Z]))) {	/* Unless outside or z = NaN */

			x[n] = in[GMT_X];
			y[n] = in[GMT_Y];
			z[n] = in[GMT_Z];
			if (z[n] < xyz[1][2]) {	/* New minimum */
				xyz[0][0] = x[n];
				xyz[0][1] = y[n];
				xyz[0][2] = z[n];
			}
			if (z[n] > xyz[1][2]) {	/* New maximum */
				xyz[1][0] = x[n];
				xyz[1][1] = y[n];
				xyz[1][2] = z[n];
			}
			n++;

			if (n == n_alloc) {
				n_alloc <<= 1;
				x = (double *) GMT_memory ((void *)x, (size_t)n_alloc, sizeof (double), GMT_program);
				y = (double *) GMT_memory ((void *)y, (size_t)n_alloc, sizeof (double), GMT_program);
				z = (double *) GMT_memory ((void *)z, (size_t)n_alloc, sizeof (double), GMT_program);
			}
			if (n == INT_MAX) {
				fprintf (stderr, "%s: ERROR: Cannot triangulate more than %d points\n", GMT_program, INT_MAX);
				GMT_free ((void *)x);
				GMT_free ((void *)y);
				GMT_free ((void *)z);
				exit (EXIT_FAILURE);
			}	
		}

		n_fields = GMT_input (fp, &n_expected_fields, &in);
	}
	if (fp != GMT_stdin) GMT_fclose (fp);

	x = (double *) GMT_memory ((void *)x, (size_t)n, sizeof (double), GMT_program);
	y = (double *) GMT_memory ((void *)y, (size_t)n, sizeof (double), GMT_program);
	z = (double *) GMT_memory ((void *)z, (size_t)n, sizeof (double), GMT_program);

	if (GMT_contlabel_prep (&Ctrl->contour, xyz)) exit (EXIT_FAILURE);	/* Prep for crossing lines, if any */

	/* Map transform */

	for (i = 0; i < n; i++) GMT_geo_to_xy (x[i], y[i], &x[i], &y[i]);

	if (Ctrl->T.active) {	/* Read precalculated triangulation indices */
		n_alloc = GMT_CHUNK;
		ind = (int *) GMT_memory (VNULL, (size_t)(3*n_alloc), sizeof (int), GMT_program);

		ij = np = 0;

		if (GMT_io.binary[GMT_IN])	/* Binary input */
			more = (GMT_fread ((void *)ind, sizeof (int), (size_t)3, fp_d) == 3);
		else		/* ascii input */
			more = (GMT_fgets (line, BUFSIZ, fp_d) != CNULL);

		while (more) {
			if (!GMT_io.binary[GMT_IN]) {	/* Ascii input, must be more careful */
				advance = FALSE;
				if (!GMT_is_a_blank_line (line)) {	/* Got a record to process */
					i = sscanf (line, "%d %d %d", &ind[ij], &ind[ij+1], &ind[ij+2]);
				 	if (i != 3) {	/* Mangled record? */
						fprintf (stderr, "%s: Error reading triangulation indices near line %ld\n", GMT_program, np);
						exit (EXIT_FAILURE);
					}
					advance = TRUE;
				}
			}
			if (advance) {	/* All is well if we get here */
				ij += 3;
				np++;
			}
			if (np == n_alloc) {
				n_alloc <<= 1;
				ind = (int *) GMT_memory ((void *)ind, (size_t)(3*n_alloc), sizeof (int), GMT_program);
			}
			if (GMT_io.binary[GMT_IN])	/* Binary input */
				more = (GMT_fread ((void *)&ind[ij], sizeof (int), (size_t)3, fp_d) == 3);
			else		/* ascii input */
				more = (GMT_fgets (line, BUFSIZ, fp_d) != CNULL);
		}
		ind = (int *) GMT_memory ((void *)ind, (size_t)ij, sizeof (int), GMT_program);
		GMT_fclose (fp_d);
	}
	else	/* Do Delaunay triangulation */

		np = GMT_delaunay (x, y, (int)n, &ind);

	if (Ctrl->L.active) {	/* Draw triangular mesh */

		GMT_setpen (&Ctrl->L.pen);

		for (i = k = 0; i < np; i++) {	/* For all triangles */

			xx[0] = x[ind[k]];	yy[0] = y[ind[k++]];
			xx[1] = x[ind[k]];	yy[1] = y[ind[k++]];
			xx[2] = x[ind[k]];	yy[2] = y[ind[k++]];

			if (project_info.three_D) for (c = 0; c < 3; c++) GMT_xy_do_z_to_xy (xx[c], yy[c], project_info.z_level, &xx[c], &yy[c]);
			ps_line (xx, yy, (GMT_LONG)3, 3, TRUE);
		}
	}

	/* Get PSCONTOUR structs */

	if (Ctrl->W.active) {
		n_contours = GMT_n_colors + 1;
		cont = (struct PSCONTOUR *) GMT_memory (VNULL, (size_t)n_contours, sizeof (struct PSCONTOUR), GMT_program);

		for (i = 0; i < GMT_n_colors; i++) {
			cont[i].val = GMT_lut[i].z_low;
			cont[i].type = (GMT_lut[i].annot && !Ctrl->A.off) ? 'A' : 'C';
			cont[i].angle = (Ctrl->contour.angle_type == 2) ? Ctrl->contour.label_angle : GMT_d_NaN;
		}
		cont[GMT_n_colors].val = GMT_lut[GMT_n_colors-1].z_high;
		cont[GMT_n_colors].type = ((GMT_lut[GMT_n_colors-1].annot & 2) && !Ctrl->A.off) ? 'A' : 'C';
		cont[GMT_n_colors].angle = (Ctrl->contour.angle_type == 2) ? Ctrl->contour.label_angle : GMT_d_NaN;
		for (i = 0; i < n_contours; i++) {
			cont[i].n_alloc = GMT_SMALL_CHUNK;
			cont[i].L = (struct PSCONTOUR_LINE *) GMT_memory (VNULL, (size_t)GMT_SMALL_CHUNK, sizeof (struct PSCONTOUR_LINE), GMT_program);
		}
	}
	Ctrl->contour.line_pen = Ctrl->W.pen;

	for (i = ij = 0; i < np; i++, ij += 3) {	/* For all triangles */

		k = ij;
		xx[0] = x[ind[k]];	yy[0] = y[ind[k]];	zz[0] = z[ind[k++]];
		xx[1] = x[ind[k]];	yy[1] = y[ind[k]];	zz[1] = z[ind[k++]];
		xx[2] = x[ind[k]];	yy[2] = y[ind[k]];	zz[2] = z[ind[k]];

		nx = get_triangle_crossings (x, y, z, &ind[ij], &xc, &yc, &zc, &vert, &cind);

		if (Ctrl->I.active) {	/* Must color the triangle slices according to cpt file */

			if (nx == 0) {	/* No contours go through - easy, but must check for NaNs */
				int kzz;
				double zzz;
				for (k = kzz = 0, zzz = 0.0; k < 3; k++) {
					if (GMT_is_dnan (zz[k])) continue;
					zzz += zz[k];
					kzz++;
				}
				if (kzz) paint_it (xx, yy, (GMT_LONG)3, zzz / kzz);
			}
			else {	/* Must paint all those slices separately */

				/* Find vertices with the lowest and highest values */

				for (k = 1, low = high = 0; k < 3; k++) {
					if (zz[k] < zz[low])   low = k;
					if (zz[k] > zz[high]) high = k;
				}

				/* Paint the piece delimited by the low node and the first contour */

				xout[0] = xx[low];	yout[0] = yy[low];
				node1 = get_node_index (vert[0], vert[1]);	/* Find single vertex opposing this contour segment */
				if (node1 == low) {	/* Contour and low node make up a triangle */
					xout[1] = xc[0];	yout[1] = yc[0];
					xout[2] = xc[1];	yout[2] = yc[1];
					n = 3;
				}
				else {	/* Need the other two vertices to form a 4-sided polygon */
					node2 = get_other_node (node1, low);	/* The other node needed */
					xout[1] = xx[node2];	yout[1] = yy[node2];
					if (low == vert[0] || node2 == vert[1]) {	/* Add segment in opposite order */
						xout[2] = xc[1];	yout[2] = yc[1];
						xout[3] = xc[0];	yout[3] = yc[0];
					}
					else  {	/* Add in regular order */
						xout[2] = xc[0];	yout[2] = yc[0];
						xout[3] = xc[1];	yout[3] = yc[1];
					}
					n = 4;
				}
				paint_it (xout, yout, n, 0.5 * (zz[low] + zc[1]));	/* z is contour value */

				/* Then loop over contours and paint the part between contours */

				for (k = 1, k2 = 2, k3 = 3; k < nx; k++, k2 += 2, k3 += 2) {
					xout[0] = xc[k2-2];	yout[0] = yc[k2-2];
					xout[1] = xc[k3-2];	yout[1] = yc[k3-2];
					n = 2;
					last_entry = vert[k2-2];
					last_exit  = vert[k3-2];
					if (last_exit == vert[k2]) {
						xout[n] = xc[k2];	yout[n] = yc[k2];	n++;
						xout[n] = xc[k3];	yout[n] = yc[k3];	n++;
						if (vert[k3] != last_entry) {	/* Need to add an intervening corner */
							node1 = get_node_index (last_entry, vert[k3]);	/* Find corner id */
							xout[n] = xx[node1];	yout[n] = yy[node1];	n++;
						}
					}
					else if (last_exit == vert[k3]) {
						xout[n] = xc[k3];	yout[n] = yc[k3];	n++;
						xout[n] = xc[k2];	yout[n] = yc[k2];	n++;
						if (vert[k2] != last_entry) {	/* Need to add an intervening corner */
							node1 = get_node_index (last_entry, vert[k2]);	/* Find corner id */
							xout[n] = xx[node1];	yout[n] = yy[node1];	n++;
						}
					}
					else if (last_entry == vert[k2]) {
						node1 = get_node_index (last_exit, vert[k3]);	/* Find corner id */
						xout[n] = xx[node1];	yout[n] = yy[node1];	n++;
						xout[n] = xc[k3];	yout[n] = yc[k3];	n++;
						xout[n] = xc[k2];	yout[n] = yc[k2];	n++;
					}
					else {
						node1 = get_node_index (last_exit, vert[k2]);	/* Find corner id */
						xout[n] = xx[node1];	yout[n] = yy[node1];	n++;
						xout[n] = xc[k2];	yout[n] = yc[k2];	n++;
						xout[n] = xc[k3];	yout[n] = yc[k3];	n++;
					}
					paint_it (xout, yout, n, 0.5 * (zc[k2]+zc[k2-2]));
				}

				/* Add the last piece between last contour and high node */

				k2 -= 2;	k3 -= 2;
				xout[0] = xx[high];	yout[0] = yy[high];
				node1 = get_node_index (vert[k2], vert[k3]);	/* Find corner id */
				if (node1 == high) {	/* Cut off a triangular piece */
					xout[1] = xc[k2];	yout[1] = yc[k2];
					xout[2] = xc[k3];	yout[2] = yc[k3];
					n = 3;
				}
				else {	/* Need a 4-sided polygon */
					node2 = get_other_node (node1, high);	/* The other node needed */
					xout[1] = xx[node2];	yout[1] = yy[node2];
					if (high == vert[0] || node2 == vert[1]) {	/* On same side, start here */
						xout[2] = xc[k3];	yout[2] = yc[k3];
						xout[3] = xc[k2];	yout[3] = yc[k2];
					}
					else  {	/* On same side, start here */
						xout[2] = xc[k2];	yout[2] = yc[k2];
						xout[3] = xc[k3];	yout[3] = yc[k3];
					}
					n = 4;
				}
				paint_it (xout, yout, n, 0.5 * (zz[high] + zc[k2]));	/* z is contour value */
			}
		}

		if (Ctrl->W.active && nx > 0) {	/* Save contour lines for later */
			for (k = k2 = 0; k < nx; k++) {
				c = cind[k];
				n = cont[c].nl;
				cont[c].L[n].x0 = xc[k2];
				cont[c].L[n].y0 = yc[k2++];
				cont[c].L[n].x1 = xc[k2];
				cont[c].L[n].y1 = yc[k2++];
				n++;
				if (n >= cont[c].n_alloc) {
					cont[c].n_alloc <<= 1;
					cont[c].L = (struct PSCONTOUR_LINE *) GMT_memory ((void *)cont[c].L, (size_t)cont[c].n_alloc, sizeof (struct PSCONTOUR_LINE), GMT_program);
				}
				cont[c].nl = (int)n;
			}
		}

		if (nx > 0) {
			GMT_free ((void *)xc);
			GMT_free ((void *)yc);
			GMT_free ((void *)zc);
			GMT_free ((void *)vert);
			GMT_free ((void *)cind);
		}
	}

	/* Draw contours */

	if (Ctrl->W.active) {

		struct PSCONTOUR_CHAIN *head_c, *last_c, *this_c;
		struct PSCONTOUR_PT *p = VNULL, *q;

		for (c = 0; c < n_contours; c++) {

			if (cont[c].nl == 0) {
				GMT_free ((void *)cont[c].L);
				continue;
			}

			head_c = last_c = (struct PSCONTOUR_CHAIN *) GMT_memory (VNULL, (size_t)1, sizeof (struct PSCONTOUR_CHAIN), GMT_program);

			while (cont[c].nl) {
				this_c = last_c->next = (struct PSCONTOUR_CHAIN *) GMT_memory (VNULL, (size_t)1, sizeof (struct PSCONTOUR_CHAIN), GMT_program);
				k = 0;
				this_c->begin = (struct PSCONTOUR_PT *) GMT_memory (VNULL, (size_t)1, sizeof (struct PSCONTOUR_PT), GMT_program);
				this_c->end = (struct PSCONTOUR_PT *) GMT_memory (VNULL, (size_t)1, sizeof (struct PSCONTOUR_PT), GMT_program);
				this_c->begin->x = cont[c].L[k].x0;
				this_c->begin->y = cont[c].L[k].y0;
				this_c->end->x = cont[c].L[k].x1;
				this_c->end->y = cont[c].L[k].y1;
				this_c->begin->next = this_c->end;
				cont[c].nl--;
				cont[c].L[k] = cont[c].L[cont[c].nl];
				while (k < cont[c].nl) {
					add = 0;
					if (fabs(cont[c].L[k].x0 - this_c->begin->x) < GMT_SMALL && fabs(cont[c].L[k].y0 - this_c->begin->y) < GMT_SMALL) {
						p = (struct PSCONTOUR_PT *) GMT_memory (VNULL, (size_t)1, sizeof (struct PSCONTOUR_PT), GMT_program);
						p->x = cont[c].L[k].x1;
						p->y = cont[c].L[k].y1;
						p->next = this_c->begin;
						add = -1;
					}
					else if (fabs(cont[c].L[k].x1 - this_c->begin->x) < GMT_SMALL && fabs(cont[c].L[k].y1 - this_c->begin->y) < GMT_SMALL) {
						p = (struct PSCONTOUR_PT *) GMT_memory (VNULL, (size_t)1, sizeof (struct PSCONTOUR_PT), GMT_program);
						p->x = cont[c].L[k].x0;
						p->y = cont[c].L[k].y0;
						p->next = this_c->begin;
						add = -1;
					}
					else if (fabs(cont[c].L[k].x0 - this_c->end->x) < GMT_SMALL && fabs(cont[c].L[k].y0 - this_c->end->y) < GMT_SMALL) {
						p = (struct PSCONTOUR_PT *) GMT_memory (VNULL, (size_t)1, sizeof (struct PSCONTOUR_PT), GMT_program);
						p->x = cont[c].L[k].x1;
						p->y = cont[c].L[k].y1;
						this_c->end->next = p;
						add = 1;
					}
					else if (fabs(cont[c].L[k].x1 - this_c->end->x) < GMT_SMALL && fabs(cont[c].L[k].y1 - this_c->end->y) < GMT_SMALL) {
						p = (struct PSCONTOUR_PT *) GMT_memory (VNULL, (size_t)1, sizeof (struct PSCONTOUR_PT), GMT_program);
						p->x = cont[c].L[k].x0;
						p->y = cont[c].L[k].y0;
						this_c->end->next = p;
						add = 1;
					}
					if (add) {	/* Got one */
						if (add == -1)
							this_c->begin = p;
						else if (add == 1)
							this_c->end = p;
						cont[c].nl--;
						cont[c].L[k] = cont[c].L[cont[c].nl];
						k = 0;
					}
					else
						k++;
				}
				last_c = this_c;
			}
			GMT_free ((void *)cont[c].L);

			this_c = head_c->next;
			while (this_c) {
				xp = (double *) GMT_memory (VNULL, (size_t)GMT_SMALL_CHUNK, sizeof (double), GMT_program);
				yp = (double *) GMT_memory (VNULL, (size_t)GMT_SMALL_CHUNK, sizeof (double), GMT_program);
				n_alloc = GMT_SMALL_CHUNK;
				p = this_c->begin;
				n = 0;
				while (p) {
					xp[n] = p->x;
					yp[n++] = p->y;
					q = p;
					p = p->next;
					GMT_free ((void *)q);
					if (n == n_alloc) {
						n_alloc <<= 1;
						xp = (double *) GMT_memory ((void *)xp, (size_t)n_alloc, sizeof (double), GMT_program);
						yp = (double *) GMT_memory ((void *)yp, (size_t)n_alloc, sizeof (double), GMT_program);
					}
				}
				last_c = this_c;
				this_c = this_c->next;
				GMT_free ((void *)last_c);

				close = !GMT_polygon_is_open (xp, yp, n);

				if (Ctrl->W.color) {
					if (current_contour != cont[c].val) {
						GMT_get_rgb_from_z (cont[c].val, rgb);
						ps_setpaint (rgb);
						memcpy ((void *)&Ctrl->contour.line_pen.rgb, (void *)rgb, (size_t)(3 * sizeof (int)));
						if (!Ctrl->contour.got_font_rgb && Ctrl->contour.curved_text) memcpy ((void *)&Ctrl->contour.font_rgb, (void *)rgb, (size_t)(3 * sizeof (int)));
						current_contour = cont[c].val;
					}
				}

				if (cont[c].type == 'A' || cont[c].type == 'a') {	/* Annotated contours */
					GMT_get_format (cont[c].val, Ctrl->contour.unit, CNULL, format);
					sprintf (cont_label, format, cont[c].val);
				}
				if (Ctrl->D.active) {
					int count;
					double *xtmp, *ytmp;
					/* Must first apply inverse map transform */
					xtmp = (double *) GMT_memory (VNULL, (size_t)n, sizeof (double), GMT_program);
					ytmp = (double *) GMT_memory (VNULL, (size_t)n, sizeof (double), GMT_program);
					for (count = 0; count < n; count++) GMT_xy_to_geo (&xtmp[count], &ytmp[count], xp[count], yp[count]);
					GMT_dump_contour (xtmp, ytmp, n, cont[c].val, section++, close, Ctrl->D.file);
					GMT_free ((void *)xtmp);
					GMT_free ((void *)ytmp);
				}

				GMT_hold_contour (&xp, &yp, n, cont[c].val, cont_label, cont[c].type, cont[c].angle, close, &Ctrl->contour);

				GMT_free ((void *)xp);
				GMT_free ((void *)yp);
			}
			GMT_free ((void *)head_c);
		}
		GMT_contlabel_plot (&Ctrl->contour);
		GMT_contlabel_free (&Ctrl->contour);
	}

    if (!Ctrl->N.active) GMT_map_clip_off ();

	GMT_map_basemap ();

	if (project_info.three_D) ps_rotatetrans (z_project.xmin, z_project.ymin, 0.0);

	GMT_plotend ();

	GMT_free ((void *)x);
	GMT_free ((void *)y);
	GMT_free ((void *)z);
	if (Ctrl->W.active) GMT_free ((void *)cont);
	if (Ctrl->T.active)
		GMT_free ((void *)ind);	/* Allocated above by GMT_memory */
	else {
#ifdef TRIANGLE_D
#ifdef DEBUG
	/* Shewchuk's function allocated the memory separately */
		GMT_memtrack_off (GMT_mem_keeper);
#endif
#endif
		GMT_free ((void *) ind);
#ifdef TRIANGLE_D
#ifdef DEBUG
		GMT_memtrack_on (GMT_mem_keeper);
#endif
#endif
	}

	Free_pscontour_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}


GMT_LONG get_triangle_crossings (double *x, double *y, double *z, int *ind, double **xc, double **yc, double **zc, GMT_LONG **v, GMT_LONG **cindex)
{
	/* This routine finds all the contour crossings for this triangle.  Each contour consists of
	 * linesegments made up of two points, with coordinates xc, yc, and contour level zc.
	 */
	 
	GMT_LONG i, j, k, k2, i1, nx, n_alloc, ok, n_ok, *vout, *cind, *ctmp = NULL;
	double xx[3], yy[3], zz[3], zmin, zmax, dz, frac, small = 1.0e-6, *xout, *yout, *zout, *ztmp = NULL;

	xx[0] = x[ind[0]];	yy[0] = y[ind[0]];	zz[0] = z[ind[0]];
	xx[1] = x[ind[1]];	yy[1] = y[ind[1]];	zz[1] = z[ind[1]];
	xx[2] = x[ind[2]];	yy[2] = y[ind[2]];	zz[2] = z[ind[2]];

	if (GMT_is_dnan (zz[0]) || GMT_is_dnan (zz[1]) || GMT_is_dnan (zz[2])) return (0);	/* Cannot have crossings if NaNs are present */
	if (zz[0] == zz[1] && zz[1] == zz[2]) return (0);	/* Cannot have crossings if all nodes are equal */

	zmin = MIN (zz[0], MIN (zz[1], zz[2]));	/* Min z vertex */
	zmax = MAX (zz[0], MAX (zz[1], zz[2]));	/* Max z vertex */

	i = 0;	j = GMT_n_colors - 1;
	while (GMT_lut[i].z_low < zmin && i < GMT_n_colors) i++;
	while (GMT_lut[j].z_high > zmax && j > 0) j--;

	nx = j - i + 2;	/* Total number of contours */

	if (nx <= 0) return (0);

	n_alloc = 2 * nx;
	xout = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (double), GMT_program);
	yout = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (double), GMT_program);
	zout = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (double), GMT_program);
	ztmp = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (double), GMT_program);
	vout = (GMT_LONG *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (GMT_LONG), GMT_program);
	cind = (GMT_LONG *) GMT_memory (VNULL, (size_t)nx, sizeof (GMT_LONG), GMT_program);
	ctmp = (GMT_LONG *) GMT_memory (VNULL, (size_t)nx, sizeof (GMT_LONG), GMT_program);

	/* Fill out array zout which holds the nx contour levels */

	k = k2 = 0;
	while (i <= j) {
		ztmp[k2] = ztmp[k2+1] = GMT_lut[i].z_low;
		ctmp[k++] = i;
		k2 += 2;
		i++;
	}
	ztmp[k2] = ztmp[k2+1] = GMT_lut[j].z_high;
	ctmp[k] = j + 1;

	/* Loop over the contour levels and determine the line segments */

	for (k = k2 = j = n_ok = 0; k < nx; k++, k2 += 2) {
		ok = FALSE;
		for (i = 0; i < 3; i++) if (zz[i] == ztmp[k2]) zz[i] += small;	/* Refuse to go through nodes */
		for (i = 0; i < 3; i++) {	/* Try each side in turn 0-1, 1-2, 2-0 */
			i1 = (i == 2) ? 0 : i + 1;
			if ((ztmp[k2] >= zz[i] && ztmp[k2] < zz[i1]) || (ztmp[k2] <= zz[i] && ztmp[k2] > zz[i1])) {
				dz = zz[i1] - zz[i];
				if (dz == 0.0) {	/* Contour goes along ende */
					xout[j] = xx[i];	yout[j] = yy[i];
				}
				else {
					frac = (ztmp[k2] - zz[i]) / dz;
					xout[j] = xx[i] + frac * (xx[i1] - xx[i]);
					yout[j] = yy[i] + frac * (yy[i1] - yy[i]);
				}
				zout[j] = ztmp[k2];
				vout[j++] = i;	/* Keep track of the side number */
				ok = TRUE;	/* Wish to add this segment */
			}
		}
		if (j%2)
			j--;	/* Contour went through a single vertex only, skip this */
		else if (ok)
			cind[n_ok++] = ctmp[k];
	}

	nx = j / 2;	/* Since j might have changed */
	if (nx) {
		*xc = xout;
		*yc = yout;
		*zc = zout;
		*v = vout;
		*cindex = cind;
	} else {
		GMT_free ((void *)xout);
		GMT_free ((void *)yout);
		GMT_free ((void *)zout);
		GMT_free ((void *)vout);
		GMT_free ((void *)cind);
	}
	GMT_free ((void *)ctmp);
	GMT_free ((void *)ztmp);
	return (nx);
}

void paint_it (double x[], double y[], GMT_LONG n, double z) {
	GMT_LONG index, k;
	int rgb[3];
	struct GMT_FILL *f;

	if (n < 3) return;	/* Need at least 3 points to make a polygon */

	index = GMT_get_rgb_from_z (z, rgb);
	if (GMT_cpt_skip) return;	/* Skip this z-slice */

	/* Now we must paint, with colors or patterns */

	if (project_info.three_D)
		for (k = 0; k < n; k++) GMT_xy_do_z_to_xy (x[k], y[k], project_info.z_level, &x[k], &y[k]);
	if ((index >= 0 && (f = GMT_lut[index].fill)) || (index < 0 && (f = GMT_bfn[index+3].fill))) {
		GMT_fill (x, y, (GMT_LONG)n, f, FALSE);	/* Contours drawn separately later if desired */
	}
	else
		ps_patch (x, y, (GMT_LONG)n, rgb, FALSE);	/* Contours drawn separately later if desired */
}

void *New_pscontour_Ctrl () {	/* Allocate and initialize a new control structure */
	struct PSCONTOUR_CTRL *C;
	
	C = (struct PSCONTOUR_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct PSCONTOUR_CTRL), "New_pscontour_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	GMT_contlabel_init (&C->contour, 1);
	C->D.file = strdup ("contour");
	C->E.azimuth = 180.0;
	C->E.elevation = 90.0;
	GMT_init_pen (&C->L.pen, GMT_PENWIDTH);
	GMT_init_pen (&C->W.pen, GMT_PENWIDTH);
		
	return ((void *)C);
}

void Free_pscontour_Ctrl (struct PSCONTOUR_CTRL *C) {	/* Deallocate control structure */
	if (C->C.file) free ((void *)C->C.file);	
	if (C->D.file) free ((void *)C->D.file);	
	if (C->T.file) free ((void *)C->T.file);	
	GMT_free ((void *)C);	
}
