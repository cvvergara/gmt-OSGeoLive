/*--------------------------------------------------------------------
 *	$Id: grdview.c,v 1.122 2011/07/08 22:39:06 guru Exp $
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
 * grdview will read a topofile and produce a 3-D perspective plot
 * of the surface z = f(x,y) using PostScript. The surface can
 * be represented as:
 *	1) A Mesh plot
 *	2) A shaded (or colored) surface w/wo contourlines and w/wo
 *	   illumination by artificial sun(s).
 *
 * grdview calls contours to find the line segments that make up the
 * contour lines. This allows the user to specify that the contours
 * should be smoothed before plotting. This will make the resulting
 * image smoother, especially if nx and ny are relatively small.
 * As an option, a drape grid file can be specified.  Then, the colors
 * are calculated from that file while the topo file is used for shape.
 * Alternatively, give three drape files (red, green, blue components)
 * to bypass the z -> rgb via the cpt file.
 *
 * Author:	Paul Wessel
 * Date:	17-SEP-2001
 * Version:	4
 */
 
#include "gmt.h"
#include "pslib.h"

/* Declarations needed for binning of smooth contours */

#define GRDVIEW_MESH		1
#define GRDVIEW_SURF		2
#define GRDVIEW_IMAGE		3

struct GRDVIEW_CTRL {
	struct C {	/* -C<cpt> */
		GMT_LONG active;
		char *file;
	} C;
	struct E {	/* -E<azimuth/elevation> */
		GMT_LONG active;
		double azimuth, elevation;
	} E;
	struct G {	/* -G<drapefile> */
		GMT_LONG active;
		GMT_LONG image;
		char *file[3];
	} G;
	struct I {	/* -G<intensfile> */
		GMT_LONG active;
		char *file;
	} I;
	struct L {	/* -L<flag> */
		GMT_LONG active;
		GMT_LONG interpolant;
		char mode[4];
		double threshold;
	} L;
	struct N {	/* -N<level>[/<color>] */
		GMT_LONG active;
		GMT_LONG facade;
		int rgb[3];
		double level;
	} N;
	struct Q {	/* -Q<type>[g] */
		GMT_LONG active, special;
		GMT_LONG outline;
		GMT_LONG mask;
		GMT_LONG monochrome;
		GMT_LONG mode;	/* GRDVIEW_MESH, GRDVIEW_SURF, GRDVIEW_IMAGE */
		GMT_LONG dpi;
		struct GMT_FILL fill;
	} Q;
	struct S {	/* -S<smooth> */
		GMT_LONG active;
		GMT_LONG value;
	} S;
	struct T {	/* -T[s][o[<pen>] */
		GMT_LONG active;
		GMT_LONG skip;
		GMT_LONG outline;
		struct GMT_PEN pen;
	} T;
	struct W {	/* -W[+]<type><pen> */
		GMT_LONG active;
		GMT_LONG contour;
		struct GMT_PEN pen[3];
	} W;
	struct Z {	/* -Z<z_level> */
		GMT_LONG active;
		double level;
	} Z;
};

struct GRDVIEW_BIN {
	struct GRDVIEW_CONT *first_cont;
};

struct GRDVIEW_CONT {
	struct GRDVIEW_POINT *first_point;
	struct GRDVIEW_CONT *next_cont;
	double value;
};

struct GRDVIEW_POINT {
	double x, y;
	struct GRDVIEW_POINT *next_point;
};

int main (int argc, char **argv)
{
	GMT_LONG get_contours, bad, set_z = FALSE, error = FALSE, pen_not_set;
	GMT_LONG first, begin, saddle, subset = FALSE, no_nans, drape_resample = FALSE;

	char *topofile = NULL;
	char *c_method[2] = {"colorimage", "colortiles",};

	GMT_LONG j, n_edges, max, i_bin, j_bin, i_bin_old, j_bin_old, nx_f, ny_f, nx, ny;
	GMT_LONG sw, se, nw, ne, id, n4, nk, c, i_start, i_stop, j_start, j_stop, i_inc, j_inc, ii, jj;
	GMT_LONG bin, i, ij, k, k1, n, nm, nm2, nm_use, mx, my, n_drape = 0, n_out, bin_inc[4], ij_inc[4];
	GMT_LONG PS_colormask_off = 0, PS_colormask, n_commas, q_set = 0, two, way, *edge = NULL;

	int rgb[3];
	
	float *grd[3], *zgrd = NULL, *intensity = NULL, *topo = NULL;

	double cval, x_left, x_right, y_top, y_bottom, small = GMT_SMALL, z_ave, this_intensity = 0.0, *xval, *yval;
	double dx2, dy2, take_out, west = 0.0, east = 0.0, south = 0.0, north = 0.0, new_z_level = 0.0;
	double data_west, data_east, data_south, data_north, delx, dely, z_val, next_up = 0.0, xmesh[4], ymesh[4];
	double x_pixel_size, y_pixel_size, *x_imask = NULL, *y_imask = NULL, x_inc[4], y_inc[4], *x = NULL, *y = NULL, *z = NULL, *v = NULL, *xx = NULL, *yy = NULL;

	struct GRDVIEW_CONT *start_cont = NULL, *this_cont = NULL, *last_cont = NULL;
	struct GRDVIEW_POINT *this_point = VNULL, *last_point = NULL;
	struct GRD_HEADER header, t_head, i_head, d_head[3];
        struct GMT_EDGEINFO edgeinfo;
	struct GMT_BCR t_bcr, i_bcr;
	struct GRDVIEW_BIN *binij = NULL;
	struct GRDVIEW_CTRL *Ctrl = NULL;

	void copy_points_fw (double x[], double y[], double z[], double v[], double xcont[], double ycont[], double zcont[], double vcont[], GMT_LONG cont, GMT_LONG *n);
	void copy_points_bw (double x[], double y[], double z[], double v[], double xcont[], double ycont[], double zcont[], double vcont[], GMT_LONG cont, GMT_LONG *n);
	double get_z_ave (double v[], double next_up, GMT_LONG n);
	void add_node (double x[], double y[], double z[], double v[], GMT_LONG *k, GMT_LONG node, double X_vert[], double Y_vert[], float topo[], float zgrd[], GMT_LONG ij, GMT_LONG bin);
	void paint_it (double x[], double y[], GMT_LONG n, double z, GMT_LONG intens, GMT_LONG monochrome, double intensity, GMT_LONG outline);
	double get_intensity (float *intensity, GMT_LONG k, GMT_LONG nx);
	struct GRDVIEW_CONT *get_cont_struct (GMT_LONG bin, struct GRDVIEW_BIN *binij, double value);
	struct GRDVIEW_POINT *get_point (double x, double y);
#if 0
	void grdview_init_setup (struct GRD_HEADER *header, float *topo, GMT_LONG two, GMT_LONG draw_plane, double plane_level);
#endif
	GMT_LONG pixel_inside (GMT_LONG ip, GMT_LONG jp, GMT_LONG *ix, GMT_LONG *iy, GMT_LONG bin, GMT_LONG bin_inc[]);
	GMT_LONG quick_idist(GMT_LONG x1, GMT_LONG y1, GMT_LONG x2, GMT_LONG y2);
	GMT_LONG get_side (double x, double y, double x_left, double y_bottom, double xinc, double yinc, double dx2, double dy2);
	GMT_LONG GMT_rgb_is_nan_rgb (int rgb[]);
	void *New_grdview_Ctrl (), Free_grdview_Ctrl (struct GRDVIEW_CTRL *C);

	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct GRDVIEW_CTRL *)New_grdview_Ctrl ();	/* Allocate and initialize a new control structure */
	
	topofile = CNULL;
	GMT_3D_mode = 1;	/* Only do background axis first; do foreground at end */
	grd[0] = grd[1] = grd[2] = (float *)NULL;
	GMT_grd_init (&header, argc, argv, FALSE);
	GMT_grd_init (&t_head, argc, argv, FALSE);
	GMT_grd_init (&i_head, argc, argv, FALSE);
	for (i = 0; i < 3; i++) GMT_grd_init (&d_head[i], argc, argv, FALSE);
	
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {

				/* Common parameters */

				case 'B':
				case 'J':
				case 'K':
				case 'O':
				case 'P':
				case 'R':
				case 'U':
				case 'V':
				case 'X':
				case 'x':
				case 'Y':
				case 'y':
				case 'c':
				case '\0':
					error += GMT_parse_common_options (argv[i], &west, &east, &south, &north);
					break;

				/* Supplemental parameters */

				case 'C':
					Ctrl->C.active = TRUE;
					Ctrl->C.file = strdup (&argv[i][2]);
					break;
				case 'E':
					Ctrl->E.active = TRUE;
					error += GMT_get_proj3D (&argv[i][2], &Ctrl->E.azimuth, &Ctrl->E.elevation);
					break;
				case 'G':
					Ctrl->G.active = TRUE;
					for (k = 2, n_commas = 0; argv[i][k]; k++) if (argv[i][k] == ',') n_commas++;
					if (n_commas == 2) {	/* Three r,g,b grids for draping */
						char A[GMT_LONG_TEXT], B[GMT_LONG_TEXT], C[GMT_LONG_TEXT];
						sscanf (&argv[i][2], "%[^,],%[^,],%s", A, B, C);
						Ctrl->G.file[0] = strdup (A);
						Ctrl->G.file[1] = strdup (B);
						Ctrl->G.file[2] = strdup (C);
						Ctrl->G.image = TRUE;
					}
					else if (n_commas == 0) {
						Ctrl->G.file[0] = strdup (&argv[i][2]);
					}
					else {
						fprintf (stderr, "%s: GMT SYNTAX ERROR option -G:  Usage is -G<z.grd> | -G<r.grd>,<g.grd>,<b.grd>\n", GMT_program);
						error = TRUE;
					}
					break;
				case 'I':
					Ctrl->I.active = TRUE;
					Ctrl->I.file = strdup (&argv[i][2]);
					break;
				case 'L':
					if (argv[i][2]) {
						Ctrl->L.active = TRUE;
						strncpy (Ctrl->L.mode, &argv[i][2], (size_t)4);
					}
					else
						Ctrl->L.interpolant = BCR_BILINEAR;
					break;
				case 'N':
					if (argv[i][2]) {
						char colors[GMT_TEXT_LEN];
						Ctrl->N.active = TRUE;
						n = sscanf (&argv[i][2], "%lf/%s", &Ctrl->N.level, colors);
						if (n == 2) {
							if (GMT_getrgb (colors, Ctrl->N.rgb)) {
								fprintf (stderr, "%s: GMT SYNTAX ERROR option -N:  Usage is -N<level>[/<color>]\n", GMT_program);
								error = TRUE;
							}
							Ctrl->N.facade = TRUE;
						}
					}
					else {
						fprintf (stderr, "%s: GMT SYNTAX ERROR option -N:  Usage is -N<level>[/<color>]\n", GMT_program);
						error = TRUE;
					}
					break;
				case 'Q':
					Ctrl->Q.active = TRUE;
					q_set++;
					switch (argv[i][2]) {
						case 'm':	/* Mesh plot */
							Ctrl->Q.mode = GRDVIEW_MESH;
							if (argv[i][3] == '/' && GMT_getfill (&argv[i][4], &Ctrl->Q.fill)) {
								fprintf (stderr, "%s: GMT SYNTAX ERROR -Qm option: To give mesh color, use -Qm/<color>\n", GMT_program);
								error = TRUE;
							}
							break;
						case 's':	/* Color wo/ contours */
							Ctrl->Q.mode = GRDVIEW_SURF;
							if (argv[i][3] == 'm') Ctrl->Q.outline = TRUE;
							break;
						case 't':	/* texture image */
							Ctrl->Q.special = TRUE;
						case 'i':	/* image w/ clipmask */
						case 'I':	/* Backward compatibility, gives -Qi */
							Ctrl->Q.mode = GRDVIEW_IMAGE;
							if (argv[i][3] && isdigit ((int)argv[i][3])) Ctrl->Q.dpi = atoi (&argv[i][3]);
							break;
						case 'c':	/* image w/ colormask */
							if (argv[i][3] && isdigit ((int)argv[i][3])) Ctrl->Q.dpi = atoi (&argv[i][3]);
							Ctrl->Q.mode = GRDVIEW_IMAGE;
							Ctrl->Q.mask = TRUE;
							break;
						default:
							fprintf (stderr, "%s: GMT SYNTAX ERROR:  Unrecognized qualifier (%c) for option -%c\n", GMT_program, argv[i][2], argv[i][1]);
							error = TRUE;
							break;
					}
					Ctrl->Q.monochrome = (argv[i][strlen(argv[i])-1] == 'g');
					break;
				case 'S':
					Ctrl->S.active = TRUE;
					Ctrl->S.value = atoi (&argv[i][2]);
					break;
				case 'T':
					Ctrl->T.active = TRUE;
					k = 2;
					if (argv[i][2] == 's') Ctrl->T.skip = TRUE, k = 3;
					if (argv[i][k] == 'o') {	/* Want tile outline also */
						Ctrl->T.outline = TRUE;
						k++;
						if (argv[i][k] && GMT_getpen (&argv[i][k], &Ctrl->T.pen)) {
							GMT_pen_syntax ('T', " ");
							error++;
						}
					}
					break;
				case 'W':	/* Contour, mesh, or facade pens */
					Ctrl->W.active = TRUE;
					j = (argv[i][2] == 'm' || argv[i][2] == 'c' || argv[i][2] == 'f') ? 3 : 2;
					id = 0;
					if (j == 3) {	/* First check that the m or c is not part of a color name instead */
						char txt_a[GMT_LONG_TEXT];
						n = j+1;
						while (argv[i][n] && argv[i][n] != ',' && argv[i][n] != '/') n++;	/* Wind until end or , or / */
						strncpy (txt_a, &argv[i][2], (size_t)(n-2));	txt_a[n-2] = '\0';
						if (GMT_colorname2index (txt_a) >= 0)	/* Found a colorname; reset j to 2 */
							j = 2;
						else
							id = (argv[i][2] == 'f') ? 2 : ((argv[i][2] == 'm') ? 1 : 0);
					}
					if (GMT_getpen (&argv[i][j], &Ctrl->W.pen[id])) {
						GMT_pen_syntax ('W', " ");
						error++;
					}
					if (j == 2)	/* Copy pen when using just -W */
						Ctrl->W.pen[2] = Ctrl->W.pen[1] = Ctrl->W.pen[0];
					if (id == 0) Ctrl->W.contour = TRUE;	/* Switch contouring ON for -Wc or -W */
					break;
				case 'Z':
					if (argv[i][2]) {
						new_z_level = atof (&argv[i][2]);
						set_z = TRUE;
					}
					break;

				/* Illegal options */

				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else
			topofile = argv[i];
	}

	if (!(Ctrl->Q.mode == GRDVIEW_MESH || Ctrl->Q.mode == GRDVIEW_SURF || Ctrl->Q.mode == GRDVIEW_IMAGE || Ctrl->T.active)) Ctrl->Q.mode = GRDVIEW_MESH;	/* Default is mesh plot */

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "grdview %s - Plot topofiles in 3-D\n\n", GMT_VERSION);
		fprintf (stderr, "usage: grdview <topofile> %s [-B<tickinfo>] [-C<cpt_file>] [%s]\n", GMT_J_OPT, GMT_E_OPT);
		fprintf (stderr, "\t[-G<drapefile> | -G<grd_r>,<grd_g>,<grd_b>] [-I<intensfile>] [%s] [-K] [-L[<flags>]]\n", GMT_Jz_OPT);
		fprintf (stderr, "\t[-N<level>[/<r/g/b>]] [-O] [-P] [-Q<type>[g]] [%s]\n", GMT_Rgeoz_OPT);
		fprintf (stderr, "\t[-S<smooth>] [-T[s][o[<pen>]]] [%s] [-V] [-W<type><pen>]\n\t[%s] [%s] [-Z<z_level>] [%s]\n\n", GMT_U_OPT, GMT_X_OPT, GMT_Y_OPT, GMT_c_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\t<topofile> is data set to be plotted\n");
		GMT_explain_option ('j');
		GMT_explain_option ('Z');
		fprintf (stderr, "\n\tOPTIONS:\n");
		GMT_explain_option ('b');
		fprintf (stderr, "\t-C Color palette file\n");
		GMT_explain_option ('E');
		fprintf (stderr, "\t-G <drapefile> rather than <topofile> is the data set to color-code.\n");
		fprintf (stderr, "\t   Use <topofile> as the relief and \'drape\' the image on top.\n");
		fprintf (stderr, "\t   Note that -Jz and -N always refers to the <topofile>\n");
		fprintf (stderr, "\t   Alternatively, give three grid files with the red, green, and blue components in 0-255 range.\n");
		fprintf (stderr, "\t   If so you must also choose -Qi.\n");
		fprintf (stderr, "\t-I gives name of intensity file and selects illumination\n");
		fprintf (stderr, "\t-L sets boundary conditions when resampling the grid.  <flags> can be either\n");
		fprintf (stderr, "\t   g for geographic boundary conditions\n");
		fprintf (stderr, "\t   or one or both of\n");
		fprintf (stderr, "\t   x for periodic boundary conditions on x\n");
		fprintf (stderr, "\t   y for periodic boundary conditions on y\n");
		fprintf (stderr, "\t   If no <flags> are set, use bilinear rather than bicubic [Default] resampling \n");
		GMT_explain_option ('Z');
		GMT_explain_option ('K');
		fprintf (stderr, "\t-N<level> will draw a plane at z = level.  Append color [/r/g/b] to paint\n");
		fprintf (stderr, "\t   the facade between the plane and the data perimeter.\n");
		GMT_explain_option ('O');
		GMT_explain_option ('P');
		fprintf (stderr, "\t-Q sets plot reQuest. Choose one of the following:\n");
		fprintf (stderr, "\t   -Qm for Mesh plot [Default].  Append /color for mesh paint [%d/%d/%d]\n", Ctrl->Q.fill.rgb[0], Ctrl->Q.fill.rgb[1], Ctrl->Q.fill.rgb[2]);
		fprintf (stderr, "\t   -Qs[m] for colored or shaded Surface.  Append m to draw meshlines on the surface.\n");
		fprintf (stderr, "\t   -Qi for scanline converting polygons to rasterimage.  Append effective dpi [100].\n");
		fprintf (stderr, "\t   -Qc. As -Qi but use PS Level 3 colormasking for nodes with z = NaN.  Append effective dpi [100].\n");
		fprintf (stderr, "\t   To force a monochrome image using the GMT_YIQ transformation, append g\n");
		GMT_explain_option ('R');
		fprintf (stderr, "\t-S will smooth contours first (see grdcontour for <smooth> value info).\n");
		GMT_pen_syntax ('T', "will image the data without interpolation by painting polygonal tiles\n\t   Append s to skip tiles for nodes with z = NaN [Default paints all tiles]\n\t   Append o[<pen>] to draw tile outline [Default uses no outline]");
		fprintf (stderr, "\t   Cannot be used with -Jz|Z as it produces a flat image.\n");
		GMT_explain_option ('U');
		GMT_explain_option ('U');
		GMT_explain_option ('V');
		GMT_pen_syntax ('W', "sets pen attributes for various features in form <type><pen>");
		fprintf (stderr, "\t   <type> can be c for contours, m for mesh, and f for facade\n");
		fprintf (stderr, "\t   c draw scontours on top of surface or mesh.  [Default is no contours]\n");
		fprintf (stderr, "\t     Optionally append pen attributes [width = %gp, color = (%d/%d/%d), solid line].\n", 
			Ctrl->W.pen[0].width, Ctrl->W.pen[0].rgb[0], Ctrl->W.pen[0].rgb[1], Ctrl->W.pen[0].rgb[2]);
		fprintf (stderr, "\t   m sets attributes for mesh lines [[width = %gp, color = (%d/%d/%d), solid line].\n", 
			Ctrl->W.pen[1].width, Ctrl->W.pen[1].rgb[0], Ctrl->W.pen[1].rgb[1], Ctrl->W.pen[1].rgb[2]);
		fprintf (stderr, "\t     Requires -Qm or -Qsm to take effect.\n");
		fprintf (stderr, "\t   f sets attributes for facade outline [[width = %gp, color = (%d/%d/%d), solid line].\n", 
			Ctrl->W.pen[2].width, Ctrl->W.pen[2].rgb[0], Ctrl->W.pen[2].rgb[1], Ctrl->W.pen[2].rgb[2]);
		fprintf (stderr, "\t     Requires -N to take effect.\n");
		GMT_explain_option ('X');
		fprintf (stderr, "\t-Z For 3-D plots: Set the z-level of map [0]\n");
		GMT_explain_option ('c');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	if (q_set > 1) {	/* Gave more than one -Q setting */
		fprintf (stderr, "%s: GMT ERROR:  -Qm, -Qs, -Qc, and -Qi are mutually exclusive options\n", GMT_program);
		error++;
	}
	if (Ctrl->T.active && Ctrl->Q.active) {	/* Gave both -Q and -T */
		fprintf (stderr, "%s: GMT ERROR:  -Q and -T are mutually exclusive options\n", GMT_program);
		error++;
	}
	if (!topofile) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify input file\n", GMT_program);
		error++;
	}
	n_drape = (Ctrl->G.image) ? 3 : 1;
	if (Ctrl->G.active) {
		for (i = 0; i < n_drape; i++) {
			if (!Ctrl->G.file[i][0]) {
				fprintf (stderr, "%s: GMT SYNTAX ERROR option -G:  Must specify drape file\n", GMT_program);
				error++;
			}
		}
		if (n_drape == 3 && !Ctrl->Q.mode == GRDVIEW_IMAGE) {
			fprintf (stderr, "%s: R/G/B drape requires -Qi option\n", GMT_program);
			error++;
		}
	}
	if (Ctrl->I.active && !Ctrl->I.file) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR option -I:  Must specify intensity file\n", GMT_program);
		error++;
	}
	if ((Ctrl->Q.mode == GRDVIEW_SURF || Ctrl->Q.mode == GRDVIEW_IMAGE || Ctrl->W.contour) && !Ctrl->C.file && !Ctrl->G.image) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify color palette table\n", GMT_program);
		error++;
	}
	if (Ctrl->Q.mode == GRDVIEW_IMAGE && Ctrl->Q.dpi <= 0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -Qi option:  Must specify positive dpi\n", GMT_program);
		error++;
	}
	if (Ctrl->T.active && project_info.JZ_set) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -T option:  Cannot specify -JZ|z\n", GMT_program);
		error++;
	}
	if (Ctrl->S.value < 0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -S option:  smooth value must be positive\n", GMT_program);
		error++;
	}
	if (Ctrl->E.elevation <= 0.0 || Ctrl->E.elevation > 90.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -E option:  Elevation must be in 0-90 range\n", GMT_program);
		error++;
	}
	if (Ctrl->L.active && GMT_boundcond_parse (&edgeinfo, Ctrl->L.mode)) error++;

	if (error) exit (EXIT_FAILURE);

	if (Ctrl->C.file) {
		GMT_read_cpt (Ctrl->C.file);
		if (GMT_b_and_w) Ctrl->Q.monochrome = TRUE;
	}
#ifdef GMT_CPT2	
	if (GMT_categorical && Ctrl->W.active) {
		fprintf (stderr, "%s: GMT WARNING -W:  Categorical data (as implied by CPT file) does not have contours.  Check plot.\n", GMT_program);
		exit (EXIT_FAILURE);
	}
#endif	
	z_project.view_azimuth = Ctrl->E.azimuth;
	z_project.view_elevation = Ctrl->E.elevation;
	get_contours = ( (Ctrl->Q.mode == GRDVIEW_MESH && Ctrl->W.contour) || (Ctrl->Q.mode == GRDVIEW_SURF && GMT_n_colors > 1) );
	PS_colormask   = (Ctrl->Q.mask) ? -1 : +1;;

	two = (Ctrl->G.active && get_contours) ? 2 : 0;	/* Must read topofile with 2 boundary columns/rows */

	if (!strcmp (topofile, "=")) {
		fprintf (stderr, "%s: Piping of topofile not supported!\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	GMT_err_fail (GMT_read_grd_info (topofile, &t_head), topofile);
	if (Ctrl->I.active) GMT_err_fail (GMT_read_grd_info (Ctrl->I.file, &i_head), Ctrl->I.file);

	if (Ctrl->I.active && (t_head.nx != i_head.nx || t_head.ny != i_head.ny)) {
		fprintf (stderr, "%s: GMT ERROR: topo and illumination files have not the same size!\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	if (Ctrl->G.active) {
		for (i = 0; i < n_drape; i++) GMT_err_fail (GMT_read_grd_info (Ctrl->G.file[i], &d_head[i]), Ctrl->G.file[i]);
		two = 2; 
	}

	GMT_boundcond_init (&edgeinfo);

	/* Determine what wesn to pass to map_setup */

	if (!project_info.region_supplied) {
		west = t_head.x_min;
		east = t_head.x_max;
		south = t_head.y_min;
		north = t_head.y_max;
	}
	else if (!(west == t_head.x_min && east == t_head.x_max && south == t_head.y_min && north == t_head.y_max))
		subset = TRUE;

	if (project_info.z_bottom == 0.0 && project_info.z_top == 0.0) {
		project_info.z_bottom = t_head.z_min;
		project_info.z_top = t_head.z_max;
		if (Ctrl->N.active && Ctrl->N.level < project_info.z_bottom) project_info.z_bottom = Ctrl->N.level;
		if (Ctrl->N.active && Ctrl->N.level > project_info.z_top) project_info.z_top = Ctrl->N.level;
	}

	GMT_err_fail (GMT_map_setup (west, east, south, north), "");

	/* Determine the wesn to be used to read the grid file */

	if (GMT_grd_setregion (&t_head, &data_west, &data_east, &data_south, &data_north, BCR_BILINEAR)) {
		/* No grid to plot; just do empty map and exit */
		GMT_plotinit (argc, argv);
		if (project_info.three_D) {
			ps_transrotate (-z_project.xmin, -z_project.ymin, 0.0);
			GMT_map_basemap ();
			GMT_vertical_axis (2);	/* Draw background axis */
			ps_rotatetrans (z_project.xmin, z_project.ymin, 0.0);
		}
		else
			GMT_map_basemap ();
		GMT_plotend ();
		GMT_end (argc, argv);
		exit (EXIT_SUCCESS);
	}

	/* Read data */

	if (gmtdefs.verbose) fprintf (stderr, "%s: Processing shape file\n", GMT_program);

	nx = GMT_get_n (data_west, data_east, t_head.x_inc, t_head.node_offset);
	ny = GMT_get_n (data_south, data_north, t_head.y_inc, t_head.node_offset);
	mx = nx + 2 * two;
	my = ny + 2 * two;
	nm = GMT_get_nm (nx, ny);
	nm2 = GMT_get_nm (mx, my);
	topo = (float *) GMT_memory (VNULL, (size_t)nm2, sizeof (float), GMT_program);
	GMT_pad[0] = GMT_pad[1] = GMT_pad[2] = GMT_pad[3] = two;
	GMT_err_fail (GMT_read_grd (topofile, &t_head, topo, data_west, data_east, data_south, data_north, GMT_pad, FALSE), topofile);
	GMT_pad[0] = GMT_pad[1] = GMT_pad[2] = GMT_pad[3] = 0;

	nx_f = t_head.nx;
	ny_f = t_head.ny;

	delx = (t_head.node_offset) ? 0.5 * t_head.x_inc :0.0;
	dely = (t_head.node_offset) ? 0.5 * t_head.y_inc :0.0;

	xval = (double *) GMT_memory (VNULL, (size_t)t_head.nx, sizeof (double), GMT_program);
	yval = (double *) GMT_memory (VNULL, (size_t)t_head.ny, sizeof (double), GMT_program);
	for (i = 0; i < t_head.nx; i++) xval[i] = GMT_i_to_x (i, t_head.x_min, t_head.x_max, t_head.x_inc, t_head.xy_off, t_head.nx);
	for (j = 0; j < t_head.ny; j++) yval[j] = GMT_j_to_y (j, t_head.y_min, t_head.y_max, t_head.y_inc, t_head.xy_off, t_head.ny);

	GMT_boundcond_param_prep (&t_head, &edgeinfo);

	if (Ctrl->G.active) {	/* draping wanted */

		for (i = 0; i < n_drape; i++) {
			if (gmtdefs.verbose) fprintf (stderr, "%s: Processing drape file %s\n", GMT_program, Ctrl->G.file[i]);

			grd[i] = (float *) GMT_memory (VNULL, (size_t)(d_head[i].nx*d_head[i].ny), sizeof (float), GMT_program);
			GMT_err_fail (GMT_read_grd (Ctrl->G.file[i], &d_head[i], grd[i], data_west, data_east, data_south, data_north, GMT_pad, FALSE), Ctrl->G.file[i]);
			if (d_head[i].nx != nx_f || d_head[i].ny != ny_f) {
				drape_resample = TRUE;
				/* fprintf (stderr, "%s: Drape file %s has improper dimensions!\n", GMT_program, Ctrl->G.file[i]);
				exit (EXIT_FAILURE); */
			}
		}
		zgrd = grd[0];
		header = d_head[0];
	}
	else {
		zgrd = topo;
		header = t_head;
	}

	if (!project_info.xyz_pos[2])	/* Negative z-scale, must flip */
		d_swap (project_info.z_bottom, project_info.z_top);

	ij_inc[0] = 0;		ij_inc[1] = 1;	ij_inc[2] = 1 - mx;	ij_inc[3] = -mx;
	nw = two * mx + two;
	ne = nw + t_head.nx - 1;
	sw = (t_head.ny + two - 1) * mx + two;
	se = sw + t_head.nx - 1;

#if 0
	grdview_init_setup (&t_head, topo, two, Ctrl->N.active, Ctrl->N.level);	/* Find projected min/max in y-direction */
#endif
	i_start = (z_project.quadrant == 1 || z_project.quadrant == 2) ? 0 : header.nx - 2;
	i_stop  = (z_project.quadrant == 1 || z_project.quadrant == 2) ? header.nx - 1 : -1;
	i_inc   = (z_project.quadrant == 1 || z_project.quadrant == 2) ? 1 : -1;
	j_start = (z_project.quadrant == 1 || z_project.quadrant == 4) ? header.ny - 1 : 1;
	j_stop  = (z_project.quadrant == 1 || z_project.quadrant == 4) ? 0 : header.ny;
	j_inc   = (z_project.quadrant == 1 || z_project.quadrant == 4) ? -1 : 1;
	bin_inc[0] = 0;		bin_inc[1] = 1;	bin_inc[2] = 1 - header.nx;	bin_inc[3] = -header.nx;
	x_inc[0] = x_inc[3] = 0.0;	x_inc[1] = x_inc[2] = header.x_inc;
	y_inc[0] = y_inc[1] = 0.0;	y_inc[2] = y_inc[3] = header.y_inc;

	if (get_contours) {	/* Need to find contours */
		if (gmtdefs.verbose) fprintf (stderr, "%s: Find contours\n", GMT_program);
		n_edges = header.ny * (GMT_LONG )ceil (header.nx / 16.0);
		edge = (GMT_LONG *) GMT_memory (VNULL, (size_t)n_edges, sizeof (GMT_LONG), GMT_program);
		binij = (struct GRDVIEW_BIN *) GMT_memory (VNULL, (size_t)nm, sizeof (struct GRDVIEW_BIN), GMT_program);
		small = GMT_SMALL * (header.z_max - header.z_min);
		if (gmtdefs.verbose) fprintf (stderr, "%s: Trace and bin contours...\n", GMT_program);
		first = TRUE;
		for (c = 0; c < GMT_n_colors+1; c++) {	/* For each color change */

			/* Reset markers and set up new zero-contour*/

			cval = (c == GMT_n_colors) ? GMT_lut[c-1].z_high : GMT_lut[c].z_low;

			if (cval < header.z_min || cval > header.z_max) continue;

			if (gmtdefs.verbose) fprintf (stderr, "%s: Now tracing contour interval %8g\r", GMT_program, cval);
			take_out = (first) ? cval : cval - GMT_lut[c-1].z_low;
			first = FALSE;
			for (i = 0; i < nm; i++) {
				if (!GMT_is_fnan (zgrd[i])) zgrd[i] -= (float)take_out;
				if (zgrd[i] == 0.0) zgrd[i] += (float)small;
			}

			begin = TRUE;
			while ((n = GMT_contours (zgrd, &header, Ctrl->S.value, gmtdefs.interpolant, 0, edge, &begin, &x, &y)) > 0) {
				i_bin_old = j_bin_old = -1;
				for (i = 1; i < n; i++) {
					i_bin = (GMT_LONG)floor (((0.5 * (x[i-1] + x[i]) - header.x_min) / header.x_inc) - header.xy_off);
					j_bin = (GMT_LONG)floor (((header.y_max - 0.5 * (y[i-1] + y[i])) / header.y_inc) - header.xy_off) + 1;
					if (i_bin != i_bin_old || j_bin != j_bin_old) {	/* Entering new bin */
						bin = j_bin * header.nx + i_bin;
						this_cont = get_cont_struct (bin, binij, cval);
						this_cont->value = cval;
						this_cont->first_point = get_point (x[i-1], y[i-1]);
						this_point = this_cont->first_point;
						i_bin_old = i_bin;
						j_bin_old = j_bin;
					}
					this_point->next_point = get_point (x[i], y[i]);
					this_point = this_point->next_point;
				}
				GMT_free ((void *)x);
				GMT_free ((void *)y);
			}
		}

		/* Remove temporary variables */

		GMT_free ((void *)edge);

		/* Go back to beginning and reread since grd has been destroyed */

		if (Ctrl->G.active) {
			GMT_err_fail (GMT_read_grd_info (Ctrl->G.file[0], &d_head[0]), Ctrl->G.file[0]);
			GMT_err_fail (GMT_read_grd (Ctrl->G.file[0], &d_head[0], grd[0], data_west, data_east, data_south, data_north, GMT_pad, FALSE), Ctrl->G.file[0]);
			header = d_head[0];
		}
		else {
			GMT_err_fail (GMT_read_grd_info (topofile, &t_head), topofile);
			GMT_pad[0] = GMT_pad[1] = GMT_pad[2] = GMT_pad[3] = two;
			GMT_err_fail (GMT_read_grd (topofile, &t_head, topo, data_west, data_east, data_south, data_north, GMT_pad, FALSE), topofile);
			GMT_pad[0] = GMT_pad[1] = GMT_pad[2] = GMT_pad[3] = 0;
			header = t_head;
		}
		if (gmtdefs.verbose) fprintf (stderr, "\n");
	}

	if (Ctrl->I.active) {	/* Illumination wanted */

		if (gmtdefs.verbose) fprintf (stderr, "%s: Processing illumination file\n", GMT_program);

		nm_use = (drape_resample) ? nm2 : nm;
		intensity = (float *) GMT_memory (VNULL, (size_t)nm_use, sizeof (float), GMT_program);
		k = (drape_resample) ? two : 0;
		GMT_pad[0] = GMT_pad[1] = GMT_pad[2] = GMT_pad[3] = k;
		GMT_err_fail (GMT_read_grd (Ctrl->I.file, &i_head, intensity, data_west, data_east, data_south, data_north, GMT_pad, FALSE), Ctrl->I.file);
		if (i_head.nx != nx_f || i_head.ny != ny_f) {
			fprintf (stderr, "%s: Intensity file has improper dimensions!\n", GMT_program);
			exit (EXIT_FAILURE);
		}
		if (drape_resample) GMT_bcr_init (&i_head, GMT_pad, Ctrl->L.interpolant, Ctrl->L.threshold, &i_bcr);
		GMT_pad[0] = GMT_pad[1] = GMT_pad[2] = GMT_pad[3] = 0;
	}

	if (two) {	/* Initialize bcr stuff */
		GMT_pad[0] = GMT_pad[1] = GMT_pad[2] = GMT_pad[3] = two;
		GMT_bcr_init (&t_head, GMT_pad, Ctrl->L.interpolant, Ctrl->L.threshold, &t_bcr);

		/* Set boundary conditions  */

		GMT_boundcond_set (&t_head, &edgeinfo, GMT_pad, topo);
		GMT_pad[0] = GMT_pad[1] = GMT_pad[2] = GMT_pad[3] = 0;
	}

	dx2 = 0.5 * header.x_inc;		dy2 = 0.5 * header.y_inc;

	no_nans = TRUE;

	for (j = 0; j < t_head.ny - 1; j++) {	/* Nodes part of tiles completely outside -R is set to NaN and thus not considered below */
		for (i = 0; i < t_head.nx - 1; i++, bin++, ij++) {
			for (jj = n_out = 0; jj < 2; jj++) {	/* Loop over the 4 nodes making up one tile */
				for (ii = 0; ii < 2; ii++) {
					bin = GMT_IJ (j + jj, i + ii, t_head.nx);
					ij = (two) ? ((j + jj) + 2) * mx + (i + ii) + two : bin;
					if (GMT_is_fnan (topo[ij])) no_nans = FALSE;
					if (GMT_map_outside (xval[i+ii], yval[j+jj])) n_out++;
				}
			}
			if (n_out == 4) topo[ij] = topo[ij-1] = topo[ij-t_head.nx] = topo[ij-t_head.nx-1] = GMT_f_NaN;	/* Entire tile is outside */
		}
	}

	max = 2 * (MAX(1,Ctrl->S.value) * (((header.nx > header.ny) ? header.nx : header.ny) + 2)) + 1;
	x = (double *) GMT_memory (VNULL, (size_t)max, sizeof (double), GMT_program);
	y = (double *) GMT_memory (VNULL, (size_t)max, sizeof (double), GMT_program);
	z = (double *) GMT_memory (VNULL, (size_t)max, sizeof (double), GMT_program);
	v = (double *) GMT_memory (VNULL, (size_t)max, sizeof (double), GMT_program);

	if (gmtdefs.verbose) fprintf (stderr, "%s: Start creating PostScript plot\n", GMT_program);

	GMT_plotinit (argc, argv);

	ps_setformat (3);

	if (project_info.three_D) {
		ps_transrotate (-z_project.xmin, -z_project.ymin, 0.0);
		GMT_map_basemap ();	/* Plot basemap first if 3-D */
	}

	if (project_info.z_pars[0] == 0.0) GMT_map_clip_on (GMT_no_rgb, 3);

	if (set_z) project_info.z_level = new_z_level;

	xx = (double *) GMT_memory (VNULL, (size_t)max, sizeof(double), GMT_program);
	yy = (double *) GMT_memory (VNULL, (size_t)max, sizeof(double), GMT_program);
	if (Ctrl->N.active) {
		ps_comment ("Plot the plane at desired level");
		GMT_setpen (&Ctrl->W.pen[2]);
		if (!z_project.draw[0])	{	/* Southern side */
			if (!project_info.region) {
				GMT_geoz_to_xy (z_project.corner_x[0], z_project.corner_y[0], Ctrl->N.level, &xx[0], &yy[0]);
				GMT_geoz_to_xy (z_project.corner_x[1], z_project.corner_y[1], Ctrl->N.level, &xx[1], &yy[1]);
				ps_line (xx, yy, 2, 3, TRUE);
			}
			else {
				for (i = 0; i < header.nx; i++) GMT_geoz_to_xy (header.x_min + i * header.x_inc + delx, header.y_min + dely, Ctrl->N.level, &xx[i], &yy[i]);
				ps_line (xx, yy, header.nx, 3, TRUE);
			}
		}
		if (!z_project.draw[2])	{	/* Northern side */
			if (!project_info.region) {
				GMT_geoz_to_xy (z_project.corner_x[3], z_project.corner_y[3], Ctrl->N.level, &xx[0], &yy[0]);
				GMT_geoz_to_xy (z_project.corner_x[2], z_project.corner_y[2], Ctrl->N.level, &xx[1], &yy[1]);
				ps_line (xx, yy, 2, 3, TRUE);
			}
			else {
				for (i = 0; i < header.nx; i++) GMT_geoz_to_xy (header.x_min + i * header.x_inc + delx, header.y_max - dely, Ctrl->N.level, &xx[i], &yy[i]);
				ps_line (xx, yy, header.nx, 3, TRUE);
			}
		}
		if (!z_project.draw[3])	{	/* Western side */
			if (!project_info.region) {
				GMT_geoz_to_xy (z_project.corner_x[0], z_project.corner_y[0], Ctrl->N.level, &xx[0], &yy[0]);
				GMT_geoz_to_xy (z_project.corner_x[3], z_project.corner_y[3], Ctrl->N.level, &xx[1], &yy[1]);
				ps_line (xx, yy, 2, 3, TRUE);
			}
			else {
				for (j = 0; j < header.ny; j++) GMT_geoz_to_xy (header.x_min + delx, header.y_max - j * header.y_inc - dely, Ctrl->N.level, &xx[j], &yy[j]);
				ps_line (xx, yy, header.ny, 3, TRUE);
			}
		}
		if (!z_project.draw[1])	{	/* Eastern side */
			if (!project_info.region) {
				GMT_geoz_to_xy (z_project.corner_x[1], z_project.corner_y[1], Ctrl->N.level, &xx[0], &yy[0]);
				GMT_geoz_to_xy (z_project.corner_x[2], z_project.corner_y[2], Ctrl->N.level, &xx[1], &yy[1]);
				ps_line (xx, yy, 2, 3, TRUE);
			}
			else {
				for (j = 0; j < header.ny; j++) GMT_geoz_to_xy (header.x_max - delx, header.y_max - j * header.y_inc - dely, Ctrl->N.level, &xx[j], &yy[j]);
				ps_line (xx, yy, header.ny, 3, TRUE);
			}
		}

		if (project_info.region) {
			GMT_geoz_to_xy (header.x_min + delx, header.y_min + dely, Ctrl->N.level, &xx[0], &yy[0]);
			GMT_geoz_to_xy (header.x_max - delx, header.y_min + dely, Ctrl->N.level, &xx[1], &yy[1]);
			GMT_geoz_to_xy (header.x_max - delx, header.y_max - dely, Ctrl->N.level, &xx[2], &yy[2]);
			GMT_geoz_to_xy (header.x_min + delx, header.y_max - dely, Ctrl->N.level, &xx[3], &yy[3]);
			if (!GMT_is_fnan (topo[nw])) {
				GMT_geoz_to_xy (header.x_min + delx, header.y_max - dely, (double)(topo[nw]), &x_left, &y_top);
				ps_segment (x_left, y_top, xx[3], yy[3]);
			}
			if (!GMT_is_fnan (topo[ne])) {
				GMT_geoz_to_xy (header.x_max - delx, header.y_max - dely, (double)(topo[ne]), &x_right, &y_top);
				ps_segment (x_right, y_top, xx[2], yy[2]);
			}
			if (!GMT_is_fnan (topo[se])) {
				GMT_geoz_to_xy (header.x_max - delx, header.y_min + dely, (double)(topo[se]), &x_right, &y_bottom);
				ps_segment (x_right, y_bottom, xx[1], yy[1]);
			}
			if (!GMT_is_fnan (topo[sw])) {
				GMT_geoz_to_xy (header.x_min + delx, header.y_min + dely, (double)(topo[sw]), &x_left, &y_bottom);
				ps_segment (x_left, y_bottom, xx[0], yy[0]);
			}
		}
	}

	if (Ctrl->T.active) {	/* Plot image as polygonal pieces. Here, -JZ is not set */
		double *xx, *yy;
		struct GMT_FILL fill;
		GMT_init_fill (&fill, -1, -1, -1);	/* Initialize fill structure */

		if (gmtdefs.verbose) fprintf (stderr, "%s: Tiling without interpolation\n", GMT_program);

		if (Ctrl->T.outline) GMT_setpen (&Ctrl->T.pen);
		for (j = k = 0; j < header.ny; j++) {
			for (i = 0; i < header.nx; i++, k++) {	/* Compute rgb for each pixel */
				if (GMT_is_fnan (topo[k]) && Ctrl->T.skip) continue;
				if (Ctrl->I.active && Ctrl->T.skip && GMT_is_fnan (intensity[k])) continue;
				GMT_get_rgb_from_z (topo[k], fill.rgb);
				if (Ctrl->I.active) GMT_illuminate (intensity[k], fill.rgb);
				n = GMT_graticule_path (&xx, &yy, 1, xval[i] - dx2, xval[i] + dx2, yval[j] - dy2, yval[j] + dy2);
				GMT_fill_polygon (xx, yy, 0.0, n, &fill, Ctrl->T.outline);
				GMT_free ((void *)xx);
				GMT_free ((void *)yy);
			}
		}
		GMT_free ((void *) xval);
		GMT_free ((void *) yval);
	}
	
	if (Ctrl->Q.mode == GRDVIEW_IMAGE) {	/* compute image */
		GMT_LONG nx_i, ny_i, kk, ip, jp, min_i, max_i, min_j, max_j, dist, node, nm_i, nm_drape, layers, last_i, last_j, p;
		GMT_LONG *top_jp = VNULL, *bottom_jp = VNULL, *ix, *iy;
		GMT_LONG done;
		double xp, yp, sum_w, w, sum_i, x_width, y_width, value;
		double sum_r, sum_g, sum_b, intval = 0.0, y_drape, *x_drape;
		float *int_drape = VNULL;
		unsigned char *bitimage_24 = (unsigned char *)NULL, *bitimage_8 = (unsigned char *)NULL;

		if (gmtdefs.verbose) {
			if (GMT_cpt_pattern) fprintf (stderr, "%s: Warning: Patterns in cpt file will not work with -Qi\n", GMT_program);
			fprintf (stderr, "%s: get and store projected vertices\n", GMT_program);
		}

		ps_comment ("Plot 3-D surface using scanline conversion of polygons to raster image");

		x_width = z_project.xmax - z_project.xmin;	/* Size of image in inches */
		y_width = z_project.ymax - z_project.ymin;
		nx_i = irint (x_width * Ctrl->Q.dpi);	/* Size of image in pixels */
		ny_i = irint (y_width * Ctrl->Q.dpi);
		last_i = nx_i - 1;	last_j = ny_i - 1;
		if (drape_resample) {
			nm_drape = d_head[0].nx * d_head[0].ny;
			ix = (GMT_LONG *) GMT_memory (VNULL, (size_t)nm_drape, sizeof (GMT_LONG), GMT_program);
			iy = (GMT_LONG *) GMT_memory (VNULL, (size_t)nm_drape, sizeof (GMT_LONG), GMT_program);
			x_drape = (double *) GMT_memory (VNULL, (size_t)d_head[0].nx, sizeof (double), GMT_program);
			if (Ctrl->I.active) int_drape = (float *) GMT_memory (VNULL, (size_t)nm_drape, sizeof (float), GMT_program);
			for (i = 0; i < d_head[0].nx; i++) x_drape[i] = GMT_i_to_x (i, d_head[0].x_min, d_head[0].x_max, d_head[0].x_inc, d_head[0].xy_off, d_head[0].nx);
			for (j = bin = 0; j < d_head[0].ny; j++) {	/* Get projected coordinates converted to pixel locations */
				y_drape = GMT_j_to_y (j, d_head[0].y_min, d_head[0].y_max,d_head[0].y_inc, d_head[0].xy_off, d_head[0].ny);
				for (i = 0; i < d_head[0].nx; i++, bin++) {
					value = GMT_get_bcr_z(&t_head, x_drape[i], y_drape, topo, &edgeinfo, &t_bcr);
					if (GMT_is_fnan (value)) {	/* Outside -R or NaNs not used */
						ix[bin] = iy[bin] = -1;
					}
					else {
						GMT_geoz_to_xy (x_drape[i], y_drape, value, &xp, &yp);
						/* Make sure ix,iy fall in the range (0,nx_i-1), (0,ny_i-1) */
						ix[bin] = MAX(0, MIN((GMT_LONG)floor((xp - z_project.xmin) * Ctrl->Q.dpi), last_i));
						iy[bin] = MAX(0, MIN((GMT_LONG)floor((yp - z_project.ymin) * Ctrl->Q.dpi), last_j));
					}
					if (Ctrl->I.active) int_drape[bin] = (float)GMT_get_bcr_z(&i_head, x_drape[i], y_drape, intensity, &edgeinfo, &i_bcr);
				}
			}
			GMT_free ((void *) x_drape);
			if (Ctrl->I.active) {	/* Reset intensity grid so that we have no boundary row/cols */
				GMT_free ((void *)intensity);
				intensity = int_drape;
			}
		}
		else {
			ix = (GMT_LONG *) GMT_memory (VNULL, (size_t)nm, sizeof (GMT_LONG), GMT_program);
			iy = (GMT_LONG *) GMT_memory (VNULL, (size_t)nm, sizeof (GMT_LONG), GMT_program);
			for (j = bin = ij = 0; j < header.ny; j++) {	/* Get projected coordinates converted to pixel locations */
				ij = (two) ? (j + 2) * mx + two : bin;
				for (i = 0; i < header.nx; i++, bin++, ij++) {
					if (GMT_is_fnan (topo[ij])) {	/* Outside -R or NaNs not used */
						ix[bin] = iy[bin] = -1;
					}
					else {
						GMT_geoz_to_xy (xval[i], yval[j], (double)topo[ij], &xp, &yp);
						/* Make sure ix,iy fall in the range (0,nx_i-1), (0,ny_i-1) */
						ix[bin] = MAX(0, MIN((GMT_LONG)floor((xp - z_project.xmin) * Ctrl->Q.dpi), last_i));
						iy[bin] = MAX(0, MIN((GMT_LONG)floor((yp - z_project.ymin) * Ctrl->Q.dpi), last_j));
					}
				}
			}
		}
		GMT_free ((void *) xval);
		GMT_free ((void *) yval);

		/* Allocate image array and set background to PAGE_COLOR */

		if (Ctrl->Q.monochrome) {
			char gray;

			nm_i = nx_i * ny_i;
			layers = 1;
			bitimage_8 = (unsigned char *) GMT_memory (VNULL, (size_t)nm_i, sizeof (char), GMT_program);
			gray = GMT_YIQ (gmtdefs.page_rgb);
			memset ((void *)bitimage_8, gray, (size_t)nm_i);
		}
		else {
			nm_i = nx_i * ny_i * 3;
			layers = 3;
			if (PS_colormask == -1) PS_colormask_off = 3;
			bitimage_24 = (unsigned char *) GMT_memory (VNULL, (size_t)(nm_i + PS_colormask_off), sizeof (char), GMT_program);
			if (PS_colormask == -1)
				memcpy ((void *)rgb, (void *)GMT_bfn[GMT_NAN].rgb, (size_t)(3 * sizeof (int)));
			else
				memcpy ((void *)rgb, (void *)gmtdefs.page_rgb, (size_t)(3 * sizeof (int)));
			kk = 0;
			while (kk < (nm_i + PS_colormask_off)) {
				bitimage_24[kk++] = (unsigned char)rgb[0];
				bitimage_24[kk++] = (unsigned char)rgb[1];
				bitimage_24[kk++] = (unsigned char)rgb[2];
			}
		}

		if (no_nans) {	/* Set up arrays for staircase clippath and initialize them */
			top_jp = (GMT_LONG *) GMT_memory (VNULL, (size_t)nx_i, sizeof (GMT_LONG), GMT_program);
			bottom_jp = (GMT_LONG *) GMT_memory (VNULL, (size_t)nx_i, sizeof (GMT_LONG), GMT_program);
			for (ip = 0; ip < nx_i; ip++) bottom_jp[ip] = ny_i;
		}

		/* Plot from back to front */

		if (gmtdefs.verbose) fprintf (stderr, "%s: Start rasterization\n", GMT_program);
		for (j = j_start; j != j_stop; j += j_inc) {

			if (gmtdefs.verbose) fprintf (stderr, "%s: Scan line conversion at j-line %.6ld\r", GMT_program, j);

			for (i = i_start; i != i_stop; i += i_inc) {
				bin = j * header.nx + i;
				ij = (two) ? (j + 2) * header.nx + i + 2 : bin;
				for (k = bad = 0; !bad && k < 4; k++) bad = (ix[bin+bin_inc[k]] < 0 || iy[bin+bin_inc[k]] < 0);
				if (bad) continue;

				min_i = max_i = ix[bin];
				min_j = max_j = iy[bin];
				for (k = 1; k < 4; k++) {
					p = bin+bin_inc[k];
					if (ix[p] < min_i) min_i = ix[p];
					if (ix[p] > max_i) max_i = ix[p];
					if (iy[p] < min_j) min_j = iy[p];
					if (iy[p] > max_j) max_j = iy[p];
				}
				for (jp = min_j; jp <= max_j; jp++) {	/* Loop over all the pixels that will make up this tile */
					if (jp < 0 || jp >= ny_i) continue;
					for (ip = min_i; ip <= max_i; ip++) {
						if (ip < 0 || ip >= nx_i) continue;
						if (!pixel_inside (ip, jp, ix, iy, bin, bin_inc)) continue;
						/* These pixels are part of the current tile */
						if (no_nans) {	/* Update clip mask */
							if (jp > top_jp[ip]) top_jp[ip] = jp; 
							if (jp < bottom_jp[ip]) bottom_jp[ip] = jp;
						}

						sum_r = sum_g = sum_b = sum_w = sum_i = 0.0;
						done = FALSE;
						for (k = bad = 0; !done && k < 4; k++) {	/* Loop over the 4 corners of the present tile */
							node = bin + bin_inc[k];
							if (Ctrl->G.image) {	/* Have 3 grids with R,G,B values */
								rgb[0] = irint ((double)grd[0][node]);	if (rgb[0] < 0) rgb[0] = 0; else if (rgb[0] > 255) rgb[0] = 255;
								rgb[1] = irint ((double)grd[1][node]);	if (rgb[1] < 0) rgb[1] = 0; else if (rgb[1] > 255) rgb[1] = 255;
								rgb[2] = irint ((double)grd[2][node]);	if (rgb[2] < 0) rgb[2] = 0; else if (rgb[2] > 255) rgb[2] = 255;
								if (GMT_rgb_is_nan_rgb (rgb)) bad++;	/* watch out for NaN colors */
							}
							else {		/* Use lookup to get color */
								GMT_get_rgb_from_z (zgrd[node], rgb);
								if (GMT_is_fnan (zgrd[node]))	/* watch out for NaNs in the z-data*/
									bad++;
								else if (Ctrl->I.active && GMT_is_fnan (intensity[node]))	/* watch out for NaNs in the intensity data*/
									bad++;
							}
							if (bad) continue;	/* We don't want to blend in the (typically) gray NaN colors with the others. */
							
							dist = quick_idist (ip, jp, ix[node], iy[node]);
							if (dist == 0) {	/* Only need this corner value */
								done = TRUE;
								if (Ctrl->I.active) intval = intensity[node];
							}
							else {	/* Crude weighted average based on 1/distance to the nearest node */
								w = 1.0 / (double)dist;
								sum_r += rgb[0] * w;
								sum_g += rgb[1] * w;
								sum_b += rgb[2] * w;
								if (Ctrl->I.active) sum_i += intensity[node] * w;
								sum_w += w;
							}
						}
						if (!done && bad < 4) {	/* Must get weighted value when more than one non-nan value was found */
							sum_w = 1.0 / sum_w;
							rgb[0] = irint (sum_r * sum_w);
							rgb[1] = irint (sum_g * sum_w);
							rgb[2] = irint (sum_b * sum_w);
							if (Ctrl->I.active) intval = sum_i * sum_w;
						}
						if (Ctrl->Q.special) GMT_get_rgb_from_z (zgrd[bin], rgb);
						if (Ctrl->I.active && bad < 4) GMT_illuminate (intval, rgb);
						kk = layers * ((ny_i-jp-1) * nx_i + ip) + PS_colormask_off;
						if (Ctrl->Q.monochrome) /* GMT_YIQ transformation */
							bitimage_8[kk] = (unsigned char) GMT_YIQ (rgb);
						else {
							bitimage_24[kk++] = (unsigned char) rgb[0];
							bitimage_24[kk++] = (unsigned char) rgb[1];
							bitimage_24[kk] = (unsigned char) rgb[2];
						}
					}
				}
			}
		}
		if (gmtdefs.verbose) fprintf (stderr, "\n");

		if (no_nans) {	/* Must implement the clip path for the image */

			x_pixel_size = x_width / (double)nx_i;
			y_pixel_size = y_width / (double)ny_i;
			n4 = 4 * nx_i;
			x_imask = (double *) GMT_memory (VNULL, (size_t)n4, sizeof (double), GMT_program);
			y_imask = (double *) GMT_memory (VNULL, (size_t)n4, sizeof (double), GMT_program);
			nk = n4 - 1;

			for (ip = k = 0; ip < nx_i; ip++, k+= 2) {
				k1 = k + 1;
				x_imask[k]  = x_imask[nk-k]  = z_project.xmin + ip * x_pixel_size;
				x_imask[k1] = x_imask[nk-k1] = x_imask[k] + x_pixel_size;
				if (top_jp[ip] < bottom_jp[ip]) {	/* No pixels set in this column */
					y_imask[k] = y_imask[k1] = y_imask[nk-k] = y_imask[nk-k1] = z_project.ymin;
				}
				else {	/* Set top of upper pixel and bottom of lower pixel */
					y_imask[k] = y_imask[k1] = z_project.ymin + (top_jp[ip] + 1) * y_pixel_size;
					y_imask[nk-k] = y_imask[nk-k1] = z_project.ymin + bottom_jp[ip] * y_pixel_size;
				}
			}
			ps_clipon (x_imask, y_imask, 4 * nx_i, GMT_no_rgb, 3);
			GMT_free ((void *)x_imask);
			GMT_free ((void *)y_imask);
		}

		if (gmtdefs.verbose) fprintf (stderr, "%s: Creating PostScript image ", GMT_program);
		if (Ctrl->Q.monochrome) {
			if (gmtdefs.verbose) fprintf (stderr, "[B/W image]\n");
			GMT_color_image (z_project.xmin, z_project.ymin, x_width, y_width, bitimage_8, nx_i, ny_i, 8);
			GMT_free ((void *)bitimage_8);
		}
		else {
			if (gmtdefs.verbose) fprintf (stderr, "[%s]\n", c_method[gmtdefs.color_image]);
			GMT_color_image (z_project.xmin, z_project.ymin, x_width, y_width, bitimage_24, nx_i * PS_colormask, ny_i, 24);
			GMT_free ((void *)bitimage_24);
		}

		if (no_nans){
			ps_clipoff ();
			GMT_free ((void *)top_jp);
			GMT_free ((void *)bottom_jp);
		}

		GMT_free ((void *)ix);
		GMT_free ((void *)iy);
	}

	if (Ctrl->Q.mode == GRDVIEW_MESH) {
		ps_comment ("Start of mesh plot");
		GMT_setpen (&Ctrl->W.pen[1]);
		if (Ctrl->Q.monochrome) Ctrl->Q.fill.rgb[0] = Ctrl->Q.fill.rgb[1] = Ctrl->Q.fill.rgb[2] = GMT_YIQ (Ctrl->Q.fill.rgb);	/* Do GMT_YIQ transformation */
		for (j = j_start; j != j_stop; j += j_inc) {
			y_bottom = yval[j];
			y_top = y_bottom + GMT_abs (j_inc) * header.y_inc;
			for (i = i_start; i != i_stop; i += i_inc) {
				bin = j * header.nx + i;
				ij = (two) ? (j + 2) * mx + i + 2 : bin;
				for (k = bad = 0; !bad && k < 4; k++) bad += GMT_is_fnan (topo[ij+ij_inc[k]]);
				if (bad) continue;
				x_left = xval[i];
				x_right = x_left + GMT_abs (i_inc) * header.x_inc;
				GMT_geoz_to_xy (x_left, y_bottom, (double)(topo[ij+ij_inc[0]]), &xx[0], &yy[0]);
				GMT_geoz_to_xy (x_right, y_bottom, (double)(topo[ij+ij_inc[1]]), &xx[1], &yy[1]);
				GMT_geoz_to_xy (x_right, y_top, (double)(topo[ij+ij_inc[2]]), &xx[2], &yy[2]);
				GMT_geoz_to_xy (x_left, y_top, (double)(topo[ij+ij_inc[3]]), &xx[3], &yy[3]);
				ps_patch (xx, yy, 4, Ctrl->Q.fill.rgb, TRUE);
				if (Ctrl->W.contour) {
					pen_not_set = TRUE;
					if (binij[bin].first_cont == NULL) continue;
					for (this_cont = binij[bin].first_cont->next_cont; this_cont; this_cont = this_cont->next_cont) {
						for (k = 0, this_point = this_cont->first_point; this_point; this_point = this_point->next_point) {
							z_val = (Ctrl->G.active) ? GMT_get_bcr_z (&t_head, (double)this_point->x, (double)this_point->y, topo, &edgeinfo, &t_bcr) : this_cont->value;
							if (GMT_is_dnan (z_val)) continue;
							GMT_geoz_to_xy ((double)this_point->x, (double)this_point->y, z_val, &xx[k], &yy[k]);
							k++;
						}
						if (pen_not_set) {
							GMT_setpen (&Ctrl->W.pen[0]);
							pen_not_set = FALSE;
						}
						ps_line (xx, yy, k, 3, FALSE);
					}
					if (!pen_not_set) GMT_setpen (&Ctrl->W.pen[1]);
				}
			}
		}
		GMT_free ((void *) xval);
		GMT_free ((void *) yval);
	}
	else if (Ctrl->Q.mode == GRDVIEW_SURF) {
		GMT_LONG start_side, entry_side, exit_side, next_side, low, ncont, nw_se_diagonal, check;
		GMT_LONG corner[2], bad_side[2][2], p, p1, p2, saddle_sign;
		double *xcont, *ycont, *zcont, *vcont, X_vert[4], Y_vert[4], saddle_small;

		xcont = (double *) GMT_memory (VNULL, (size_t)max, sizeof (double), GMT_program);
		ycont = (double *) GMT_memory (VNULL, (size_t)max, sizeof (double), GMT_program);
		zcont = (double *) GMT_memory (VNULL, (size_t)max, sizeof (double), GMT_program);
		vcont = (double *) GMT_memory (VNULL, (size_t)max, sizeof (double), GMT_program);

		ps_comment ("Start of filled surface");
		if (Ctrl->Q.outline) GMT_setpen (&Ctrl->W.pen[1]);

		for (j = j_start; j != j_stop; j += j_inc) {
			y_bottom = yval[j];
			y_top = y_bottom + header.y_inc;
			for (i = i_start; i != i_stop; i += i_inc) {
				bin = j * header.nx + i;
				ij = (two) ? (j + 2) * mx + i + 2 : bin;
				x_left = xval[i];
				x_right = x_left + header.x_inc;
				for (k = bad = 0; !bad && k < 4; k++) bad += GMT_is_fnan (topo[ij+ij_inc[k]]);
				if (bad) {
					if (GMT_bfn[GMT_NAN].skip || project_info.three_D) continue;

					X_vert[0] = X_vert[3] = x_left;	X_vert[1] = X_vert[2] = x_right;
					Y_vert[0] = Y_vert[1] = y_bottom;	Y_vert[2] = Y_vert[3] = y_top;
					for (k = 0; k < 4; k++) GMT_geoz_to_xy (X_vert[k], Y_vert[k], 0.0, &xmesh[k], &ymesh[k]);
					paint_it (xmesh, ymesh, 4, GMT_d_NaN, FALSE, Ctrl->Q.monochrome, 0.0, Ctrl->Q.outline);
					continue;
				}

				if (Ctrl->I.active) {
					this_intensity = get_intensity (intensity, bin, header.nx);
					if (GMT_is_dnan (this_intensity)) continue;

				}
				/* Get mesh polygon */

				X_vert[0] = X_vert[3] = x_left;	X_vert[1] = X_vert[2] = x_right;
				Y_vert[0] = Y_vert[1] = y_bottom;	Y_vert[2] = Y_vert[3] = y_top;

				if (get_contours && binij[bin].first_cont) {	/* Contours go thru here */

					/* Determine if this bin will give us saddle trouble */

					start_cont = this_cont = binij[bin].first_cont->next_cont;
					saddle = FALSE;
					while (!saddle && this_cont->next_cont) {
						if (this_cont->next_cont->value == this_cont->value)
							saddle = TRUE;
						else
							this_cont = this_cont->next_cont;
					}
					if (saddle) {	/* Must deal with this separately */

						this_point = this_cont->first_point;
						entry_side = get_side (this_point->x, this_point->y, x_left, y_bottom, header.x_inc, header.y_inc, dx2, dy2);
						while (this_point->next_point) this_point = this_point->next_point;	/* Go to end */
						exit_side  = get_side (this_point->x, this_point->y, x_left, y_bottom, header.x_inc, header.y_inc, dx2, dy2);


						if (MIN (zgrd[bin+bin_inc[1]], zgrd[bin+bin_inc[3]]) > MAX (zgrd[bin], zgrd[bin+bin_inc[2]])) {
							saddle_sign = +1;
							check = TRUE;
						}
						else if (MAX (zgrd[bin+bin_inc[1]], zgrd[bin+bin_inc[3]]) < MIN (zgrd[bin], zgrd[bin+bin_inc[2]])) {
							saddle_sign = -1;
							check = TRUE;
						}
						else if (MIN (zgrd[bin], zgrd[bin+bin_inc[2]]) > MAX (zgrd[bin+bin_inc[1]], zgrd[bin+bin_inc[3]])) {
							saddle_sign = +1;
							check = FALSE;
						}
						else {
							saddle_sign = -1;
							check = FALSE;
						}
						nw_se_diagonal = ((entry_side + exit_side) == 3);
						if (nw_se_diagonal != check) saddle_sign = -saddle_sign;
						if (nw_se_diagonal) {	/* Diagonal goes NW - SE */
							corner[0] = 0;	bad_side[0][0] = 1;	bad_side[0][1] = 2;
							corner[1] = 2;	bad_side[1][0] = 0;	bad_side[1][1] = 3;
						}
						else {	/* Diagonal goes NE -SW */
							corner[0] = 1;	bad_side[0][0] = 2;	bad_side[0][1] = 3;
							corner[1] = 3;	bad_side[1][0] = 0;	bad_side[1][1] = 1;
						}
						saddle_small = saddle_sign * small;

						for (p = 0; p < 2; p++) {	/* For each triangular half */

							/* Set this points as the start anchor */

							low = corner[p];
							n = 0;
							add_node (x, y, z, v, &n, low, X_vert, Y_vert, topo, zgrd, ij+ij_inc[low], bin+bin_inc[low]);
							start_side = next_side = low;
							way = 0;

							for (this_cont = start_cont; this_cont; this_cont = this_cont->next_cont) {

								/* First get all the x/y pairs for this contour */

								for (k = 0, this_point = this_cont->first_point; this_point; this_point = this_point->next_point) {
									xcont[k] = this_point->x;
									ycont[k] = this_point->y;
									zcont[k] = (Ctrl->G.active) ? GMT_get_bcr_z (&t_head, xcont[k], ycont[k], topo, &edgeinfo, &t_bcr) : this_cont->value;
									if (GMT_is_dnan (zcont[k])) continue;
									vcont[k] = this_cont->value;
									k++;
								}
								ncont = k;

								entry_side = get_side (xcont[0], ycont[0], x_left, y_bottom, header.x_inc, header.y_inc, dx2, dy2);
								exit_side  = get_side (xcont[ncont-1], ycont[ncont-1], x_left, y_bottom, header.x_inc, header.y_inc, dx2, dy2);

								if (entry_side == bad_side[p][0] || entry_side == bad_side[p][1]) continue;
								if (exit_side == bad_side[p][0] || exit_side == bad_side[p][1]) continue;

								/* OK, got the correct contour */

								next_up = (this_cont->next_cont) ? this_cont->next_cont->value : DBL_MAX;

								exit_side  = get_side (xcont[ncont-1], ycont[ncont-1], x_left, y_bottom, header.x_inc, header.y_inc, dx2, dy2);

								if (way == 0 || next_side == entry_side) {	/* Just hook up */
									copy_points_fw (x, y, z, v, xcont, ycont, zcont, vcont, ncont, &n);
									next_side = exit_side;
								}
								else if (next_side == exit_side) {	/* Just hook up but reverse */
									copy_points_bw (x, y, z, v, xcont, ycont, zcont, vcont, ncont, &n);
									next_side = entry_side;
								}
								/* Compute the xy from the xyz triplets */

								for (k = 0; k < n; k++) GMT_geoz_to_xy (x[k], y[k], z[k], &xx[k], &yy[k]);
								z_ave = (GMT_continuous) ? get_z_ave (v, next_up, n) : this_cont->value;

								/* Now paint the polygon piece */

								paint_it (xx, yy, n, z_ave-saddle_small, Ctrl->I.active, Ctrl->Q.monochrome, this_intensity, FALSE);

								/* Reset the anchor points to previous contour */

								n = 0;
								copy_points_fw (x, y, z, v, xcont, ycont, zcont, vcont, ncont, &n);
								next_side = exit_side;
								start_side = entry_side;
								way = (zgrd[bin+bin_inc[low]] < this_cont->value) ? -1 : 1;
							}

							/* Final contour needs to add diagonal */

							if (corner[p] == 0 || corner[p] == 2) {
								p1 = (next_side < 2) ? 1 : 3;
								p2 = (next_side < 2) ? 3 : 1;
							}
							else {
								p1 = (next_side % 3) ? 2 : 0;
								p2 = (next_side % 3) ? 0 : 2;
							}
							add_node (x, y, z, v, &n, p1, X_vert, Y_vert, topo, zgrd, ij+ij_inc[p1], bin+bin_inc[p1]);
							add_node (x, y, z, v, &n, p2, X_vert, Y_vert, topo, zgrd, ij+ij_inc[p2], bin+bin_inc[p2]);

							/* Compute the xy from the xyz triplets */

							for (k = 0; k < n; k++) GMT_geoz_to_xy (x[k], y[k], z[k], &xx[k], &yy[k]);

							z_ave = (GMT_continuous) ? get_z_ave (v, next_up, n) : v[0];

							/* Now paint the polygon piece */

							paint_it (xx, yy, n, z_ave+saddle_small, Ctrl->I.active, Ctrl->Q.monochrome, this_intensity, FALSE);

						} /* End triangular piece */

					} /* End Saddle section */
					else {
						/* Ok, here we do not have to worry about saddles */

						/* Find lowest corner (id = low) */

						for (k = 1, low = 0; k < 4; k++) if (zgrd[bin+bin_inc[k]] < zgrd[bin+bin_inc[low]]) low = k;

						/* Set this points as the start anchor */

						n = 0;
						add_node (x, y, z, v, &n, low, X_vert, Y_vert, topo, zgrd, ij+ij_inc[low], bin+bin_inc[low]);
						start_side = next_side = low;
						way = 1;

						this_cont = start_cont;
						while (this_cont) {

							next_up = (this_cont->next_cont) ? this_cont->next_cont->value : DBL_MAX;
							/* First get all the x/y pairs for this contour */

							for (k = 0, this_point = this_cont->first_point; this_point; this_point = this_point->next_point) {
								xcont[k] = this_point->x;
								ycont[k] = this_point->y;
								zcont[k] = (Ctrl->G.active) ? GMT_get_bcr_z (&t_head, xcont[k], ycont[k], topo, &edgeinfo, &t_bcr) : this_cont->value;
								if (GMT_is_dnan (zcont[k])) continue;
								vcont[k] = this_cont->value;
								k++;
							}
							ncont = k;

							entry_side = get_side (xcont[0], ycont[0], x_left, y_bottom, header.x_inc, header.y_inc, dx2, dy2);
							exit_side  = get_side (xcont[ncont-1], ycont[ncont-1], x_left, y_bottom, header.x_inc, header.y_inc, dx2, dy2);

							while (!(next_side == entry_side || next_side == exit_side)) {	/* Must add intervening corner */
								if (way == 1) next_side = (next_side + 1) % 4;
								add_node (x, y, z, v, &n, next_side, X_vert, Y_vert, topo, zgrd, ij+ij_inc[next_side], bin+bin_inc[next_side]);
								if (way == -1) next_side = (next_side - 1 + 4) % 4;
							}
							if (next_side == entry_side) {	/* Just hook up */
								copy_points_fw (x, y, z, v, xcont, ycont, zcont, vcont, ncont, &n);
								next_side = exit_side;
							}
							else if (next_side == exit_side) {	/* Just hook up but reverse */
								copy_points_bw (x, y, z, v, xcont, ycont, zcont, vcont, ncont, &n);
								next_side = entry_side;
							}
							/* Now we must complete the polygon if necessary */

							while (!(start_side == next_side)) {	/* Must add intervening corner */
								if (way == 1) next_side = (next_side + 1) % 4;
								add_node (x, y, z, v, &n, next_side, X_vert, Y_vert, topo, zgrd, ij+ij_inc[next_side], bin+bin_inc[next_side]);
								if (way == -1) next_side = (next_side - 1 + 4) % 4;
							}

							/* Compute the xy from the xyz triplets */

							for (k = 0; k < n; k++) GMT_geoz_to_xy (x[k], y[k], z[k], &xx[k], &yy[k]);
							z_ave = (GMT_continuous) ? get_z_ave (v, next_up, n) : this_cont->value;

							/* Now paint the polygon piece */

							paint_it (xx, yy, n, z_ave-small, Ctrl->I.active, Ctrl->Q.monochrome, this_intensity, FALSE);

							/* Reset the anchor points to previous contour */

							n = 0;
							copy_points_fw (x, y, z, v, xcont, ycont, zcont, vcont, ncont, &n);
							next_side = exit_side;
							start_side = entry_side;
							way = (zgrd[bin+bin_inc[start_side]] < this_cont->value) ? -1 : 1;

							this_cont = this_cont->next_cont;	/* Goto next contour */
 						}

						/* Final contour needs to compete with corners only */

						while (!(start_side == next_side)) {	/* Must add intervening corner */
							if (way == 1) next_side = (next_side +1) % 4;
							add_node (x, y, z, v, &n, next_side, X_vert, Y_vert, topo, zgrd, ij+ij_inc[next_side], bin+bin_inc[next_side]);
							if (way == -1) next_side = (next_side - 1 + 4) % 4;
						}

						/* Compute the xy from the xyz triplets */

						for (k = 0; k < n; k++) GMT_geoz_to_xy (x[k], y[k], z[k], &xx[k], &yy[k]);

						z_ave = (GMT_continuous) ? get_z_ave (v, next_up, n) : v[0];

						/* Now paint the polygon piece */

						paint_it (xx, yy, n, z_ave+small, Ctrl->I.active, Ctrl->Q.monochrome, this_intensity, FALSE);

					} /* End non-saddle case */

					/* Draw contour lines if desired */

					pen_not_set = TRUE;
					for (this_cont = start_cont; Ctrl->W.contour && this_cont; this_cont = this_cont->next_cont) {
						for (k = 0, this_point = this_cont->first_point; this_point; this_point = this_point->next_point) {
							z_val = (Ctrl->G.active) ? GMT_get_bcr_z (&t_head, (double)this_point->x, (double)this_point->y, topo, &edgeinfo, &t_bcr) : this_cont->value;
							if (GMT_is_dnan (z_val)) continue;

							GMT_geoz_to_xy ((double)this_point->x, (double)this_point->y, z_val, &xx[k], &yy[k]);
							k++;
						}
						if (pen_not_set) {
							GMT_setpen (&Ctrl->W.pen[0]);
							pen_not_set = FALSE;
						}
						ps_line (xx, yy, k, 3, FALSE);
					}
					if (!pen_not_set) GMT_setpen (&Ctrl->W.pen[1]);
					if (Ctrl->Q.outline) {
						for (k = 0; k < 4; k++) GMT_geoz_to_xy (X_vert[k], Y_vert[k], (double)(topo[ij+ij_inc[k]]), &xmesh[k], &ymesh[k]);
						ps_patch (xmesh, ymesh, 4, GMT_no_rgb, TRUE);
					}
				}
				else {	/* No Contours */

					if (GMT_continuous) {	/* Take the color corresponding to the average value of the four corners */
						for (n = 0, z_ave = 0.0; n < 4; n++) z_ave += zgrd[bin+bin_inc[n]];
						z_ave *= 0.25;
					}
					else	/* Take the value of any corner */
						z_ave = zgrd[bin];
						
					/* Now paint the polygon piece */

					for (k = 0; k < 4; k++) GMT_geoz_to_xy (X_vert[k], Y_vert[k], (double)(topo[ij+ij_inc[k]]), &xmesh[k], &ymesh[k]);
					paint_it (xmesh, ymesh, 4, z_ave, Ctrl->I.active, Ctrl->Q.monochrome, this_intensity, Ctrl->Q.outline);
				}
			}
		}
		GMT_free ((void *) xval);
		GMT_free ((void *) yval);
		GMT_free ((void *) xcont);
		GMT_free ((void *) ycont);
		GMT_free ((void *) zcont);
		GMT_free ((void *) vcont);
	}

	if (Ctrl->W.pen[1].texture || Ctrl->W.pen[0].texture) ps_setdash (CNULL, 0);

	if (project_info.z_pars[0] == 0.0) GMT_map_clip_off();

	if (Ctrl->N.facade) {	/* Cover the two front sides */
		ps_comment ("Painting the frontal facade");
		GMT_setpen (&Ctrl->W.pen[2]);
		if (!z_project.draw[0])	{	/* Southern side */
			for (i = n = 0, ij = sw; i < header.nx; i++, ij++) {
				if (GMT_is_fnan (topo[ij])) continue;
				GMT_geoz_to_xy (header.x_min+i*header.x_inc+delx, header.y_min+dely, (double)(topo[ij]), &xx[n], &yy[n]);
				n++;
			}
			for (i = header.nx - 1; i >= 0; i--, n++) GMT_geoz_to_xy (header.x_min+i*header.x_inc+delx, header.y_min+dely, Ctrl->N.level, &xx[n], &yy[n]);
			ps_polygon (xx, yy, n, Ctrl->N.rgb, TRUE);
		}
		if (!z_project.draw[1]) {	/*	Eastern side */
			for (j = n = 0, ij = ne; j < header.ny; j++, ij += mx) {
				if (GMT_is_fnan (topo[ij])) continue;
				GMT_geoz_to_xy (header.x_max-delx, header.y_max-j*header.y_inc-dely, (double)(topo[ij]), &xx[n], &yy[n]);
				n++;
			}
			for (j = header.ny - 1; j >= 0; j--, n++) GMT_geoz_to_xy (header.x_max-delx, header.y_max-j*header.y_inc-dely, Ctrl->N.level, &xx[n], &yy[n]);
			ps_polygon (xx, yy, n, Ctrl->N.rgb, TRUE);
		}
		if (!z_project.draw[2])	{	/* Northern side */
			for (i = n = 0, ij = nw; i < header.nx; i++, ij++) {
				if (GMT_is_fnan (topo[ij])) continue;
				GMT_geoz_to_xy (header.x_min+i*header.x_inc+delx, header.y_max-dely, (double)(topo[ij]), &xx[n], &yy[n]);
				n++;
			}
			for (i = header.nx - 1; i >= 0; i--, n++) GMT_geoz_to_xy (header.x_min+i*header.x_inc+delx, header.y_max-dely, Ctrl->N.level, &xx[n], &yy[n]);
			ps_polygon (xx, yy, n, Ctrl->N.rgb, TRUE);
		}
		if (!z_project.draw[3]) {	/*	Western side */
			for (j = n = 0, ij = nw; j < header.ny; j++, ij += mx) {
				if (GMT_is_fnan (topo[ij])) continue;
				GMT_geoz_to_xy (header.x_min+delx, header.y_max-j*header.y_inc-dely, (double)(topo[ij]), &xx[n], &yy[n]);
				n++;
			}
			for (j = header.ny - 1; j >= 0; j--, n++) GMT_geoz_to_xy (header.x_min+delx, header.y_max-j*header.y_inc-dely, Ctrl->N.level, &xx[n], &yy[n]);
			ps_polygon (xx, yy, n, Ctrl->N.rgb, TRUE);
		}
	}

	if (project_info.three_D) {
		GMT_vertical_axis (2);	/* Draw background axis */
		ps_rotatetrans (z_project.xmin, z_project.ymin, 0.0);
	}
	else
		GMT_map_basemap ();	/* Plot basemap last if not 3-D */

	GMT_plotend ();

	/* Free memory */

	if (get_contours) {
		for (ij = 0; ij < nm; ij++) {
			if (!binij[ij].first_cont) continue;
			last_cont = binij[ij].first_cont;
			for (this_cont = binij[ij].first_cont->next_cont; this_cont; this_cont = this_cont->next_cont) {
				if (this_cont->first_point) {
					last_point = this_cont->first_point;
					for (this_point = this_cont->first_point->next_point; this_point; this_point = this_point->next_point) {
						GMT_free ((void *)last_point);
						last_point = this_point;
					}
					GMT_free ((void *)last_point);
				}
				GMT_free ((void *)last_cont);
				last_cont = this_cont;
			}
			GMT_free ((void *)last_cont);
		}
		GMT_free ((void *)binij);
	}

	GMT_free ((void *)xx);
	GMT_free ((void *)yy);
	GMT_free ((void *)x);
	GMT_free ((void *)y);
	GMT_free ((void *)z);
	GMT_free ((void *)v);
	GMT_free ((void *)topo);
	if (Ctrl->I.active) GMT_free ((void *)intensity);
	if (Ctrl->G.active) for (i = 0; i < n_drape; i++) GMT_free ((void *)grd[i]);

	if (gmtdefs.verbose) fprintf (stderr, "%s: Done!\n", GMT_program);

	Free_grdview_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

struct GRDVIEW_CONT *get_cont_struct (GMT_LONG bin, struct GRDVIEW_BIN *binij, double value)
{
	struct GRDVIEW_CONT *cont, *new_cont;

	if (!binij[bin].first_cont) binij[bin].first_cont = (struct GRDVIEW_CONT *) GMT_memory (VNULL, (size_t)1, sizeof (struct GRDVIEW_CONT), GMT_program);

	for (cont = binij[bin].first_cont; cont->next_cont && cont->next_cont->value <= value; cont = cont->next_cont);

	new_cont = (struct GRDVIEW_CONT *) GMT_memory (VNULL, (size_t)1, sizeof (struct GRDVIEW_CONT), GMT_program);
	if (cont->next_cont) {	/* Put it in the link */
		new_cont->next_cont = cont->next_cont;
		cont->next_cont = new_cont;
	}
	else	/* End of list */
		cont->next_cont = new_cont;
	return (new_cont);
}

struct GRDVIEW_POINT *get_point (double x, double y)
{
	struct GRDVIEW_POINT *point;

	point = (struct GRDVIEW_POINT *) GMT_memory (VNULL, (size_t)1, sizeof (struct GRDVIEW_POINT), GMT_program);
	point->x = x;
	point->y = y;
	return (point);
}

#if 0
/* Removed this because it yields unpredictable results, making it impossible to line up different 3D plots */

void grdview_init_setup (struct GRD_HEADER *header, float *topo, GMT_LONG two, GMT_LONG draw_plane, double plane_level)
{
	GMT_LONG i, j, ij;
	double xtmp, ytmp, tmp, delx, dely;

	delx = (header->node_offset) ? 0.5 * header->x_inc :0.0;
	dely = (header->node_offset) ? 0.5 * header->y_inc :0.0;
	/* Find projected min/max in y-direction */

	z_project.ymax = z_project.ymin;	/* Reset from whatever it was */

	for (j = 0; j < header->ny; j++) {
		ij = GMT_IJ (j + two, two, header->nx);
		for (i = 0; i < header->nx; i++, ij++) {
			if (GMT_is_fnan (topo[ij])) continue;
			GMT_geoz_to_xy (header->x_min + i * header->x_inc, header->y_max - j * header->y_inc, (double)topo[ij] + project_info.z_bottom, &xtmp, &ytmp);
			z_project.ymin = MIN (z_project.ymin, ytmp);
			z_project.ymax = MAX (z_project.ymax, ytmp);
		}
	}
	if (draw_plane) {	/* plane or facade may exceed the found min/max */
		for (i = 0; i < header->nx; i++) {
			tmp = header->x_min + i * header->x_inc + delx;
			GMT_geoz_to_xy (tmp, header->y_min + dely, plane_level + project_info.z_bottom, &xtmp, &ytmp);
			z_project.ymin = MIN (z_project.ymin, ytmp);
			z_project.ymax = MAX (z_project.ymax, ytmp);
			GMT_geoz_to_xy (tmp, header->y_max-dely, plane_level + project_info.z_bottom, &xtmp, &ytmp);
			z_project.ymin = MIN (z_project.ymin, ytmp);
			z_project.ymax = MAX (z_project.ymax, ytmp);
		}
		for (j = 0; j < header->ny; j++) {
			tmp = header->y_max - j * header->y_inc - dely;
			GMT_geoz_to_xy (header->x_min+delx, tmp, plane_level + project_info.z_bottom, &xtmp, &ytmp);
			z_project.ymin = MIN (z_project.ymin, ytmp);
			z_project.ymax = MAX (z_project.ymax, ytmp);
			GMT_geoz_to_xy (header->x_max-delx, tmp, plane_level + project_info.z_bottom, &xtmp, &ytmp);
			z_project.ymin = MIN (z_project.ymin, ytmp);
			z_project.ymax = MAX (z_project.ymax, ytmp);
		}
	}
}
#endif

double get_intensity (float *intensity, GMT_LONG k, GMT_LONG nx)
{
	/* Finds the average intensity for this polygon */
	return (0.25 * (intensity[k] + intensity[k+1] + intensity[k-nx] + intensity[k-nx+1]));
}

GMT_LONG pixel_inside (GMT_LONG ip, GMT_LONG jp, GMT_LONG *ix, GMT_LONG *iy, GMT_LONG bin, GMT_LONG bin_inc[])
{
	GMT_LONG i, what;
	double x[6], y[6];

	for (i = 0; i < 4; i++) {
		x[i] = (double)ix[bin+bin_inc[i]];
		y[i] = (double)iy[bin+bin_inc[i]];
	}
	x[4] = x[0];	y[4] = y[0];
	what = GMT_non_zero_winding ((double)ip, (double)jp, x, y, 5);
	return (what);
}

GMT_LONG quick_idist (GMT_LONG x1, GMT_LONG y1, GMT_LONG x2, GMT_LONG y2)
{
	if ((x2 -= x1) < 0) x2 = -x2;
	if ((y2 -= y1) < 0) y2 = -y2;
	return (x2 + y2 - (((x2 > y2) ? y2 : x2) >> 1));
}

GMT_LONG get_side (double x, double y, double x_left, double y_bottom, double xinc, double yinc, double dx2, double dy2) {
	/* Figure out on what side this point sites on */

	double del_x, del_y;
	GMT_LONG side;

	del_x = x - x_left;
	if (del_x > dx2) del_x = xinc - del_x;
	del_y = y - y_bottom;
	if (del_y > dy2) del_y = yinc - del_y;
	if (del_x < del_y) /* Cutting N-S gridlines */
		side = ((x-x_left) > dx2) ? 1 : 3;
	else /* Cutting E-W gridlines */
		side = ((y-y_bottom) > dy2) ? 2 : 0;
	return (side);
}

void copy_points_fw (double x[], double y[], double z[], double v[], double xcont[], double ycont[], double zcont[], double vcont[], GMT_LONG ncont, GMT_LONG *n) {
	GMT_LONG k;
	for (k = 0; k < ncont; k++, (*n)++) {
		x[*n] = xcont[k];
		y[*n] = ycont[k];
		z[*n] = zcont[k];
		v[*n] = vcont[k];
	}
}

void copy_points_bw (double x[], double y[], double z[], double v[], double xcont[], double ycont[], double zcont[], double vcont[], GMT_LONG ncont, GMT_LONG *n) {
	GMT_LONG k;
	for (k = ncont - 1; k >= 0; k--, (*n)++) {
		x[*n] = xcont[k];
		y[*n] = ycont[k];
		z[*n] = zcont[k];
		v[*n] = vcont[k];
	}
}

double get_z_ave (double v[], double next_up, GMT_LONG n) {
	GMT_LONG k;
	double z_ave;

	for (k = 0, z_ave = 0.0; k < n; k++) z_ave += MIN (v[k], next_up);
	return (z_ave / n);
}

void add_node (double x[], double y[], double z[], double v[], GMT_LONG *k, GMT_LONG node, double X_vert[], double Y_vert[], float topo[], float zgrd[], GMT_LONG ij, GMT_LONG bin) {
	/* Adds a corner node to list of points and increments counter */
	x[*k] = X_vert[node];
	y[*k] = Y_vert[node];
	z[*k] = topo[ij];
	v[*k] = zgrd[bin];
	(*k)++;
}

void paint_it (double x[], double y[], GMT_LONG n, double z, GMT_LONG intens, GMT_LONG monochrome, double intensity, GMT_LONG outline) {
	GMT_LONG index;
	int rgb[3];
	struct GMT_FILL *f;

	if (n < 3) return;	/* Need at least 3 points to make a polygon */

	index = GMT_get_rgb_from_z (z, rgb);
	if (GMT_cpt_skip) return;	/* Skip this z-slice */

	/* Now we must paint, with colors or patterns */

	if ((index >= 0 && (f = GMT_lut[index].fill)) || (index < 0 && (f = GMT_bfn[index+3].fill))) {	/* Pattern */
		GMT_fill (x, y, n, f, outline);
	}
	else {	/* Solid color/gray */
		if (intens) GMT_illuminate (intensity, rgb);
		if (monochrome) rgb[0] = rgb[1] = rgb[2] = GMT_YIQ (rgb); /* GMT_YIQ transformation */
		ps_patch (x, y, n, rgb, outline);	/* Contours drawn separately (after this call) if desired */
	}
}

GMT_LONG GMT_rgb_is_nan_rgb (int rgb[])
{
	return (rgb[0] == GMT_bfn[GMT_NAN].rgb[0] && rgb[1] == GMT_bfn[GMT_NAN].rgb[1] && rgb[2] == GMT_bfn[GMT_NAN].rgb[2]);
}

void *New_grdview_Ctrl () {	/* Allocate and initialize a new control structure */
	struct GRDVIEW_CTRL *C;
	
	C = (struct GRDVIEW_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GRDVIEW_CTRL), "New_grdview_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	C->E.azimuth = 180.0;
	C->E.elevation = 90.0;
	GMT_init_pen (&C->T.pen, GMT_PENWIDTH);
	GMT_init_pen (&C->W.pen[0], 3.0 * GMT_PENWIDTH);	/* Contour pen */
	GMT_init_pen (&C->W.pen[1], GMT_PENWIDTH);	/* Mesh pen */
	GMT_init_pen (&C->W.pen[2], 3.0 * GMT_PENWIDTH);	/* Facade pen */
	C->Q.dpi = 100;
	GMT_init_fill (&C->Q.fill, 255 - gmtdefs.basemap_frame_rgb[0], 255 - gmtdefs.basemap_frame_rgb[0], 255 - gmtdefs.basemap_frame_rgb[0]);
	C->L.interpolant = BCR_BICUBIC; C->L.threshold = 1.0;
	C->S.value = 1;
	return ((void *)C);
}

void Free_grdview_Ctrl (struct GRDVIEW_CTRL *C) {	/* Deallocate control structure */
	GMT_LONG i;
	if (C->C.file) free ((void *)C->C.file);	
	for (i = 0; i < 3; i++) if (C->G.file[i]) free ((void *)C->G.file[i]);	
	if (C->I.file) free ((void *)C->I.file);	
	GMT_free ((void *)C);	
}
