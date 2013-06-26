/*--------------------------------------------------------------------
 *	$Id: grdcontour.c 9923 2012-12-18 20:45:53Z pwessel $
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
 * grdcontour reads a 2-D grid file and contours it. This algorithm handles
 * cases where the old contourm failed, e.g. when a contour line passes
 * exactly through a data point. This is achieved by adding a tiny
 * number to those data points that would have a contour line going
 * through them. This is cheating, but the human eye cannot tell the
 * difference, so who cares?
 * A multitude of options exist and are explained in the usage message.
 * The 2-D grid file format is outlined in the routines read_grd/write_grd.
 *
 * Author:	Paul Wessel
 * Date:	6-JAN-1991
 * Version:	2.0	Based on old v1.x by Paul Wessel & Walter Smith
 * Revised:	21-JUL-1998 for GMT 3.1
 * Revised:	03-MAR-1999 for GMT 3.2
 * Revised:	25-FEB-2000 for GMT 3.3.4
 * Revised:	25-JUN-2000 for GMT 3.3.5
 * Version:	3.4 Added Kaj Jancke's modification for transparent label boxes
 * Version:	4.0 15-JUL-2004 Implemented label clipping in PostScript. Expanded -G option
 * Version:	4.1.x
 *
 */

#include "gmt.h"
#include "pslib.h"

#define GRDCONTOUR_MIN_LENGTH 0.01

struct GRDCONTOUR_CTRL {
	struct GMT_CONTOUR contour;
	struct A {	/* -A[-|aint][+a<angle>][+c<dx>[/<dy>]][+f<font>][+g<fill>][+j<just>][+l<label>][+o|O|t][+s<size>][+p<pen>][+u<unit>] */
		GMT_LONG active;
		GMT_LONG off;
		double interval;
	} A;
	struct C {	/* -C<cont_int> */
		GMT_LONG active;
		GMT_LONG cpt;
		char *file;
		double interval;
	} C;
	struct D {	/* -D<dumpfile> */
		GMT_LONG active;
		char *file;
	} D;
	struct E {	/* -Eazimuth/elevation */
		GMT_LONG active;
		double azimuth, elevation;
	} E;
	struct F {	/* -F<way> */
		GMT_LONG active;
		GMT_LONG value;
	} F;
	struct G {	/* -G[d|f|n|l|L|x|X]<params> */
		GMT_LONG active;
	} G;
	struct L {	/* -L<Low/high> */
		GMT_LONG active;
		double low, high;
	} L;
	struct Q {	/* -Q<cut> */
		GMT_LONG active;
		GMT_LONG min;
	} Q;
	struct S {	/* -S<smooth> */
		GMT_LONG active;
		GMT_LONG value;
	} S;
	struct T {	/* -T[+|-][<gap>[c|i|m|p]/<length>[c|i|m|p]][:LH] */
		GMT_LONG active;
		GMT_LONG label;
		GMT_LONG low, high;	/* TRUE to tick low and high locals */
		double spacing, length;
		char L_label, H_label;
	} T;
	struct W {	/* -W[+]<type><pen> */
		GMT_LONG active;
		GMT_LONG color;
		struct GMT_PEN pen[2];
	} W;
	struct Z {	/* -Z[<fact>[/shift>]][p] */
		GMT_LONG active;
		GMT_LONG periodic;
		double scale, offset;
	} Z;
};

struct SAVE {
	double *x, *y;
	double cval;
	GMT_LONG n;
	struct GMT_PEN pen;
	GMT_LONG do_it, high;
};

int main (int argc, char **argv)
{
	GMT_LONG j, k, c, n_edges, n_alloc = 0, section, n_contours, id;
	GMT_LONG nx, ny, n_save = 0, got, nm, i, n, nn, *edge = NULL;
	
	int rgb[3];

	GMT_LONG error = FALSE, interior, begin;
#ifdef EXPERIMENT
	GMT_LONG wrap_check;
#endif

	char *grdfile = NULL, line[BUFSIZ], *cont_type = NULL, *cont_do_tick = NULL;
	char txt_a[GMT_LONG_TEXT], txt_b[GMT_LONG_TEXT];
	char cont_label[GMT_LONG_TEXT], format[GMT_LONG_TEXT];

	double aval, west = 0.0, east = 0.0, south = 0.0, north = 0.0, *xp = NULL, *yp = NULL;
	double *contour = NULL, cval, min, max, small, xyz[2][3];
	double small_x, small_y, data_west, data_east, data_south, data_north, tmp, *cont_angle = NULL;
	double *x = NULL, *y = NULL;	/* Arrays holding the contour xy values */

	float *grd_original = NULL, *grd = NULL;

	struct GRD_HEADER header;
	struct SAVE *save = NULL;
	struct GRDCONTOUR_CTRL *Ctrl = NULL;

	void sort_and_plot_ticks (struct SAVE *save, GMT_LONG n, struct GRD_HEADER *h, float *grd, double tick_gap, double tick_length, GMT_LONG tick_low, GMT_LONG tick_high, GMT_LONG tick_label, char Llabel, char Hlabel);
	void GMT_grd_minmax (float *a, struct GRD_HEADER *h, double xyz[2][3]);
	void adjust_hill_label (struct GMT_CONTOUR *G, struct GRD_HEADER *h, float *z);
	void *New_grdcontour_Ctrl (), Free_grdcontour_Ctrl (struct GRDCONTOUR_CTRL *C);
	
	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct GRDCONTOUR_CTRL *)New_grdcontour_Ctrl ();	/* Allocate and initialize a new control structure */
	
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {

				/* Common parameters */

				case 'B':
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
				case 'f':
				case 'm':
				case '\0':
					error += GMT_parse_common_options (argv[i], &west, &east, &south, &north);
					break;

				/* Supplemental parameters */

				case 'A':
					/* Format can be one of two:
					 * 3.4.x: -A[-|aint][f<fontsize>][a<angle>][/<r/g/b>][t|o]
					 * 4.x:   -A[-|aint][+a<angle>][+c<dx>[/<dy>]][+f<font>][+g<fill>][+j<just>][+l<label>][+o|O|t][+s<size>][+p<pen>][+u<unit>]
					 */

					Ctrl->A.active = TRUE;
					error += GMT_contlabel_specs (&argv[i][2], &Ctrl->contour);
					if (argv[i][2] == '-')
						Ctrl->A.off = TRUE;
					else {
						n = sscanf (&argv[i][2], "%lf", &Ctrl->A.interval);
						Ctrl->contour.annot = TRUE;
					}
					if (error) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -A option.  Correct syntax:\n", GMT_program);
						fprintf (stderr, "\t-A[-][aint][+a<angle>][+c<dx>[/<dy>]][+f<font>][+g[<fill>]][+j<just>][+o][+p[<pen>]][+s<size>][+u<unit>][+v]\n");
					}
					break;
				case 'C':
					Ctrl->C.active = TRUE;
					if (!GMT_access (&argv[i][2], R_OK)) {	/* Gave a readable file */
						Ctrl->C.interval = 1.0;
						Ctrl->C.cpt = (!strncmp (&argv[i][strlen(argv[i])-4], ".cpt", (size_t)4)) ? TRUE : FALSE;
						Ctrl->C.file = strdup (&argv[i][2]);
					}
					else
						Ctrl->C.interval = atof (&argv[i][2]);
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
				case 'F':
					Ctrl->F.active = TRUE;
					switch (argv[i][2]) {
						case '\0':
						case 'l':
						case 'L':
							Ctrl->F.value = -1;
							break;
						case 'R':
						case 'r':
							Ctrl->F.value = +1;
							break;
						default:
							fprintf (stderr, "%s: GMT SYNTAX ERROR -F option.  Correct syntax is -F[l|r]:\n", GMT_program);
							break;
					}
					break;
				case 'G':
					Ctrl->G.active = TRUE;
					error += GMT_contlabel_info ('G', &argv[i][2], &Ctrl->contour);
					break;
				case 'L':
					Ctrl->L.active = TRUE;
					sscanf (&argv[i][2], "%lf/%lf", &Ctrl->L.low, &Ctrl->L.high);
					break;
				case 'N':	/* Backwards compatibility - now done in -A instead */
					if (argv[i][2])
						strcpy (Ctrl->contour.unit, &argv[i][2]);
					else
						strcpy (Ctrl->contour.unit, "z");
					break;
				case 'S':
					Ctrl->S.active = TRUE;
					Ctrl->S.value = atoi (&argv[i][2]);
					break;
				case 'Q':
					Ctrl->Q.active = TRUE;
					Ctrl->Q.min = atoi (&argv[i][2]);
					break;
				case 'T':
					Ctrl->T.active = Ctrl->T.high = Ctrl->T.low = TRUE;	/* Default if just -T is given */
					if (argv[i][2]) {	/* But here we gave more options */
						if (argv[i][2] == '+')			/* Only tick local highs */
							Ctrl->T.low = FALSE, j = 1;
						else if (argv[i][2] == '-')		/* Only tick local lows */
							Ctrl->T.high = FALSE, j = 1;
						else
							j = 0;
						n = 0;
						if (strchr (&argv[i][2+j], '/')) {	/* Gave gap/length */
							n = sscanf (&argv[i][2+j], "%[^/]/%[^:]", txt_a, txt_b);
							if (n == 2) {
								Ctrl->T.spacing = GMT_convert_units (txt_a, GMT_INCH);
								Ctrl->T.length = GMT_convert_units (txt_b, GMT_INCH);
							}
						}
						for (j = 2; argv[i][j] && argv[i][j] != ':'; j++);
						if (argv[i][j]) {	/* Gave high/low markers */
							j++;
							if (argv[i][j])
								Ctrl->T.L_label = argv[i][j];
							else
								error = TRUE;
							j++;
							if (argv[i][j])
								Ctrl->T.H_label = argv[i][j];
							else
								error = TRUE;
							Ctrl->T.label = TRUE;
						}
						if (n == 1 || Ctrl->T.spacing <= 0.0 || Ctrl->T.length == 0.0) {
							fprintf (stderr, "%s: GMT SYNTAX ERROR -T option.  Correct syntax:\n", GMT_program);
							fprintf (stderr, "\t-T[+|-][<tick_gap>[c|i|m|p]/<tick_length>[c|i|m|p]][:LH], <tick_gap> must be > 0\n");
							error = TRUE;
						}
					}
					break;
				case 'W':
					Ctrl->W.active = TRUE;
					k = 2;
					if (argv[i][k] == '+') Ctrl->W.color = TRUE, k++;
					j = (argv[i][k] == 'a' || argv[i][k] == 'c') ? k+1 : k;
					if (j == k) {	/* Set both */
						if (GMT_getpen (&argv[i][j], &Ctrl->W.pen[0])) {
							GMT_pen_syntax ('W', " ");
							error++;
						}
						else Ctrl->W.pen[1] = Ctrl->W.pen[0];
					}
					else {	/* Gave a or c.  Because the user may say -Wcyan we must prevent this from being seen as -Wc and color yan! */
						/* Get the argument following a or c and up to first comma, slash (or to the end) */
						n = k+1;
						while (!(argv[i][n] == ',' || argv[i][n] == '/' || argv[i][n] == '\0')) n++;
						strncpy (txt_a, &argv[i][k], (size_t)(n-k));	txt_a[n-k] = '\0';
						if (GMT_colorname2index (txt_a) >= 0) j = k;	/* Found a colorname; wind j back by 1 */
						id = (argv[i][k] == 'a') ? 1 : 0;
						if (GMT_getpen (&argv[i][j], &Ctrl->W.pen[id])) {
							GMT_pen_syntax ('W', " ");
							error++;
						}
						if (j == k) Ctrl->W.pen[1] = Ctrl->W.pen[0];	/* Must copy since it was not -Wc nor -Wa after all */
					}
					break;
				case 'Z':
					Ctrl->Z.active = TRUE;
					if (argv[i][2] && argv[i][2] != 'p') n = sscanf (&argv[i][2], "%lf/%lf", &Ctrl->Z.scale, &Ctrl->Z.offset);
					Ctrl->Z.periodic = (argv[i][strlen(argv[i])-1] == 'p');	/* Phase data */
					break;
				default:
					error++;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else
			grdfile = argv[i];
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr,"grdcontour %s - Contouring of 2-D gridded data sets\n\n", GMT_VERSION);
		fprintf (stderr, "usage: grdcontour <grdfile> -C<cont_int> %s\n", GMT_J_OPT);
		fprintf (stderr, "\t[-A[-|<annot_int>][<labelinfo>] [%s] [-D<dumpfile>] [%s] [-F[l|r]] [%s]\n", GMT_B_OPT, GMT_E_OPT, GMT_CONTG);
		fprintf (stderr, "\t[-K] [-L<low/high>] [-O] [-P] [-Q<cut>] [%s] [-S<smooth>]\n", GMT_Rgeo_OPT);
		fprintf (stderr, "\t[-T[+|-][<gap>[c|i|m|p]/<length>[c|i|m|p]][:LH]] [%s] [-V] [-W[+]<type><pen>]\n", GMT_U_OPT);
		fprintf (stderr, "\t[%s] [%s] [-Z[<fact>[/shift>]][p]] [%s] [%s] [%s]\n\n", GMT_X_OPT, GMT_Y_OPT, GMT_bo_OPT, GMT_c_OPT, GMT_mo_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\t<grdfile> is 2-D netCDF grid file to be contoured.\n");
		fprintf (stderr, "\t-C Contours to be drawn can be specified in one of three ways:\n");
		fprintf (stderr, "\t   1. Fixed contour interval.\n");
		fprintf (stderr, "\t   2. Name of file with contour levels in col 1 and C(ont) or A(nnot) in col 2\n");
		fprintf (stderr, "\t      [and optionally an individual annotation angle in col 3].\n");
		fprintf (stderr, "\t   3. Name of cpt-file.\n");
		fprintf (stderr, "\t   If -T is used, only contours with upper case C or A is ticked\n");
		fprintf (stderr, "\t     [cpt-file contours are set to C unless last column has flags; Use -A to force all to A].\n");
		GMT_explain_option ('j');
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-A Annotation label information. [Default is no annoted contours].\n");
		fprintf (stderr, "\t   Give annotation interval OR - to disable all contour annotations implied in -C.\n");
		fprintf (stderr, "\t   <labelinfo> controls the specifics of the labels.  Append what you need:\n");
		GMT_label_syntax (5, 0);
		GMT_explain_option ('b');
		fprintf (stderr, "\t-D to Dump contour lines to individual files (but see -m).\n");
		fprintf (stderr, "\t   Append file prefix [contour].  Files will be called <dumpfile>_<cont>_#[_i].xyz|b,\n");
		fprintf (stderr, "\t   where <cont> is the contour value and # is a segment counter;\n");
		fprintf (stderr, "\t   _i is inserted for interior (closed) contours, with xyz (ascii) or b (binary) as extension.\n");
		fprintf (stderr, "\t   However, if -D- is given then files are C#_e or C#_i plus extension, where # is a running number.\n");
		GMT_explain_option ('E');
		fprintf (stderr, "\t-F force dumped contours to be oriented so that the higher z-values are to the left (-Fl [Default])\n");
		fprintf (stderr, "\t   or right (-Fr) as we move along the contour [Default is not oriented].\n");
		fprintf (stderr, "\t-G Controls placement of labels along contours.  Choose among five algorithms:\n");
		GMT_cont_syntax (3, 0);
		GMT_explain_option ('K');
		fprintf (stderr, "\t-L only contour inside the specified range.\n");
		GMT_explain_option ('O');
		GMT_explain_option ('P');
		fprintf (stderr, "\t-Q Do not draw closed contours with less than <cut> points [Draw all contours].\n");
		GMT_explain_option ('R');
		fprintf (stderr, "\t   [Default is extent of grid].\n");
		fprintf (stderr, "\t-S will Smooth contours by splining and resampling\n");
		fprintf (stderr, "\t   at approximately gridsize/<smooth> intervals.\n");
		fprintf (stderr, "\t-T will embellish innermost, closed contours with ticks pointing in the downward direction.\n");
		fprintf (stderr, "\t   User may specify to tick only highs (-T+) or lows (-T-) [-T means both].\n");
		fprintf (stderr, "\t   Append spacing/ticklength (append units) to change defaults [%g/%g %s].\n",
			Ctrl->T.spacing*GMT_u2u[GMT_INCH][gmtdefs.measure_unit], Ctrl->T.length*GMT_u2u[GMT_INCH][gmtdefs.measure_unit], GMT_unit_names[gmtdefs.measure_unit]);
		fprintf (stderr, "\t   Append :LH to plot the characters L and H in the center of closed contours\n");
		fprintf (stderr, "\t   for local Lows and Highs (e.g, give :-+ to plot - and + signs).\n");
		GMT_explain_option ('U');
		GMT_explain_option ('V');
		GMT_pen_syntax ('W', "sets pen attributes. Append a<pen> for annotated contours or c<pen> for regular contours [Default].");
		fprintf (stderr, "\t   The default settings are:\n");
		fprintf (stderr, "\t   Contour pen:  width = %gp, color = (%d/%d/%d), texture = solid\n", Ctrl->W.pen[0].width, Ctrl->W.pen[0].rgb[0], Ctrl->W.pen[0].rgb[1], Ctrl->W.pen[0].rgb[2]);
		fprintf (stderr, "\t   Annotate pen: width = %gp, color = (%d/%d/%d), texture = solid\n", Ctrl->W.pen[1].width, Ctrl->W.pen[1].rgb[0], Ctrl->W.pen[1].rgb[1], Ctrl->W.pen[1].rgb[2]);
		fprintf (stderr, "\t   Use + to draw colored contours based on the cpt file.\n");
		GMT_explain_option ('X');
		fprintf (stderr, "\t-Z to subtract <shift> and multiply data by <fact> before contouring [1/0].\n");
		fprintf (stderr, "\t   Append p for z-data that is periodic in 360 (i.e., phase data).\n");
		GMT_explain_option ('o');
		GMT_explain_option ('n');
		GMT_explain_option ('c');
		GMT_explain_option ('f');
		fprintf (stderr, "\t-m Used with -D.   Create a single multiple segment file where contours are separated by a record\n");
		fprintf (stderr, "\t   whose first character is <flag> ['>'].  This header also has the contour level value.\n");
		exit (EXIT_FAILURE);
	}

	if (Ctrl->A.interval > 0.0 && (!Ctrl->C.file && Ctrl->C.interval == 0.0)) Ctrl->C.interval = Ctrl->A.interval;
	if (!Ctrl->C.file && Ctrl->C.interval <= 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -C option:  Must specify contour interval, file name with levels, or cpt-file\n", GMT_program);
		error++;
	}
	if (Ctrl->L.low >= Ctrl->L.high) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -L option: lower limit >= upper!\n", GMT_program);
		error++;
	}
	if (Ctrl->F.active && !Ctrl->D.active) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -F option: Must also specify -D\n", GMT_program);
		error++;
	}
	if (Ctrl->S.value < 0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -S option:  Smooth_factor must be > 0\n", GMT_program);
		error++;
	}
	if (Ctrl->Q.min < 0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -Q option:  Value must be >= 0\n", GMT_program);
		error++;
	}
	if (Ctrl->E.elevation <= 0.0 || Ctrl->E.elevation > 90.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -E option:  Elevation must be in 0-90 range\n", GMT_program);
		error++;
	}
	if (!grdfile) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify input file\n", GMT_program);
		error++;
	}
	if (Ctrl->contour.label_dist_spacing <= 0.0 || Ctrl->contour.half_width <= 0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -G option.  Correct syntax:\n", GMT_program);
		fprintf (stderr, "\t-G<annot_dist>/<npoints>, both values must be > 0\n");
		error++;
	}
	if (Ctrl->Z.scale == 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -Z option:  factor must be nonzero\n", GMT_program);
		error++;
	}
	if (Ctrl->W.color && !Ctrl->C.cpt) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -W option:  + only valid if -C sets a cpt file\n", GMT_program);
		error++;
	}

	if (error) exit (EXIT_FAILURE);

	if (Ctrl->D.active && Ctrl->D.file[0] == 0) {
		fprintf (stderr, "%s: contours will be written to files contour_*.xyz\n", GMT_program);
		strcpy (Ctrl->D.file,"contour");
	}
	if (GMT_io.multi_segments[GMT_OUT]) GMT_io.multi_segments[GMT_OUT] = 2;

	GMT_z_periodic = Ctrl->Z.periodic;	/* Phase data */
	z_project.view_azimuth = Ctrl->E.azimuth;
	z_project.view_elevation = Ctrl->E.elevation;
	if (gmtdefs.verbose) fprintf (stderr, "%s: Allocate memory and read data file\n", GMT_program);

	if (!strcmp (grdfile,  "=")) {
		fprintf (stderr, "%s: Piping of grid file not supported!\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	GMT_err_fail (GMT_read_grd_info (grdfile, &header), grdfile);
	if (!strcmp (Ctrl->contour.unit, "z")) strcpy (Ctrl->contour.unit, header.z_units);

	/* Determine what wesn to pass to map_setup */

	if (!project_info.region_supplied) {
		west = header.x_min;
		east = header.x_max;
		south = header.y_min;
		north = header.y_max;
	}
	GMT_err_fail (GMT_map_setup (west, east, south, north), "");

	/* Determine the wesn to be used to read the grid file */

	if (GMT_grd_setregion (&header, &data_west, &data_east, &data_south, &data_north, BCR_BILINEAR)) {
		/* No grid to plot; just do empty map and exit */
		if (gmtdefs.verbose) fprintf (stderr, "%s: Warning: No data within specified region\n", GMT_program);
		GMT_plotinit (argc, argv);
		GMT_map_basemap ();
		GMT_plotend ();
		Free_grdcontour_Ctrl (Ctrl);	/* Deallocate control structure */
		GMT_end (argc, argv);
		exit (EXIT_SUCCESS);
	}

	nx = GMT_get_n (data_west, data_east, header.x_inc, header.node_offset);
	ny = GMT_get_n (data_south, data_north, header.y_inc, header.node_offset);
	nm = GMT_get_nm (nx, ny);

	/* Read data */
	
	grd_original = (float *) GMT_memory (VNULL, (size_t)nm, sizeof (float), GMT_program);
	GMT_err_fail (GMT_read_grd (grdfile, &header, grd_original, data_west, data_east, data_south, data_north, GMT_pad, FALSE), grdfile);

	if (!(Ctrl->Z.scale == 1.0 && Ctrl->Z.offset == 0.0)) {
		if (gmtdefs.verbose) fprintf (stderr, "%s: Subtracting %g and multiplying grid by %g\n", GMT_program, Ctrl->Z.offset, Ctrl->Z.scale);
		for (i = 0; i < nm; i++) grd_original[i] = (float)((grd_original[i] - Ctrl->Z.offset) * Ctrl->Z.scale);
		header.z_min = (header.z_min - Ctrl->Z.offset) * Ctrl->Z.scale;
		header.z_max = (header.z_max - Ctrl->Z.offset) * Ctrl->Z.scale;
		if (Ctrl->Z.scale < 0.0) d_swap (header.z_min, header.z_max);
	}
	if (Ctrl->L.low > header.z_min) header.z_min = Ctrl->L.low;
	if (Ctrl->L.high < header.z_max) header.z_max = Ctrl->L.high;
	if (Ctrl->L.active && header.z_max < header.z_min) {	/* Specified contour range outside range of data - quit */
		if (gmtdefs.verbose) fprintf (stderr, "%s: Warning: No contours within specified -L range\n", GMT_program);
		GMT_free ((void *)grd_original);
		Free_grdcontour_Ctrl (Ctrl);	/* Deallocate control structure */
		GMT_end (argc, argv);
		exit (EXIT_SUCCESS);
	}
	
	small = Ctrl->C.interval * 1.0e-6;
	if (Ctrl->A.interval == 0.0) Ctrl->A.interval = Ctrl->C.interval;

	small_x = 0.01 * header.x_inc;	small_y = 0.01 * header.y_inc;

	if (Ctrl->contour.annot) {	/* Wants annotated contours */
		aval = floor (header.z_min / Ctrl->A.interval) * Ctrl->A.interval;
		if (aval < header.z_min) aval += Ctrl->A.interval;
	}
	else
		aval = header.z_max + 1.0;

	if (Ctrl->C.cpt) {	/* Presumably got a cpt-file */
		GMT_read_cpt (Ctrl->C.file);
#ifdef GMT_CPT2	
		if (GMT_categorical) {
			fprintf (stderr, "%s: GMT WARNING:  Categorical data (as implied by CPT file) does not have contours.  Check plot.\n", GMT_program);
			exit (EXIT_FAILURE);
		}
#endif
		n_contours = GMT_n_colors + 1;
		(void)GMT_alloc_memory2 ((void **)&contour, (void **)&cont_angle, n_contours, 0, sizeof (double), GMT_program);
		n_alloc = GMT_alloc_memory2 ((void **)&cont_type, (void **)&cont_do_tick, n_contours, 0, sizeof (char), GMT_program);
		for (i = c = 0; i < GMT_n_colors; i++) {
			if (GMT_lut[i].skip) continue;
			contour[c] = GMT_lut[i].z_low;
			if (Ctrl->A.off)
				cont_type[c] = 'C';
			else if (GMT_lut[i].annot)
				cont_type[c] = 'A';
			else
				cont_type[c] = (Ctrl->contour.annot) ? 'A' : 'C';
			cont_angle[c] = (Ctrl->contour.angle_type == 2) ? Ctrl->contour.label_angle : GMT_d_NaN;
			cont_do_tick[c] = (char)Ctrl->T.active;
			c++;
		}
		contour[c] = GMT_lut[GMT_n_colors-1].z_high;
		if (Ctrl->A.off)
			cont_type[c] = 'C';
		else if (GMT_lut[GMT_n_colors-1].annot & 2)
			cont_type[c] = 'A';
		else
			cont_type[c] = (Ctrl->contour.annot) ? 'A' : 'C';
		cont_angle[c] = (Ctrl->contour.angle_type == 2) ? Ctrl->contour.label_angle : GMT_d_NaN;
		cont_do_tick[c] = (char)Ctrl->T.active;
		n_contours = c + 1;
	}
	else if (Ctrl->C.file) {	/* read contour info from file */
		FILE *fpc = NULL;

		n_contours = 0;
		if ((fpc = GMT_fopen (Ctrl->C.file, "r")) == NULL) {
			fprintf (stderr, "%s: Error opening contour info file %s\n", GMT_program, Ctrl->C.file);
			exit (EXIT_FAILURE);
		}
		while (GMT_fgets (line, BUFSIZ, fpc)) {
			GMT_chop (line);
			if (line[0] == '#' || line[0] == '\0') continue;
			if (n_contours == n_alloc) {
				(void)GMT_alloc_memory2 ((void **)&contour, (void **)&cont_angle, n_contours, n_alloc, sizeof (double), GMT_program);
				n_alloc = GMT_alloc_memory2 ((void **)&cont_type, (void **)&cont_do_tick, n_contours, n_alloc, sizeof (char), GMT_program);
			}
			got = sscanf (line, "%lf %c %lf", &contour[n_contours], &cont_type[n_contours], &tmp);
			if (cont_type[n_contours] == 0) cont_type[n_contours] = 'C';
			cont_do_tick[n_contours] = (Ctrl->T.active && ((cont_type[n_contours] == 'C') || (cont_type[n_contours] == 'A'))) ? 1 : 0;
			cont_angle[n_contours] = (got == 3) ? tmp : GMT_d_NaN;
			if (got == 3) Ctrl->contour.angle_type = 2;	/* Must set this directly if angles are provided */
			n_contours++;
		}
		GMT_fclose (fpc);
	}
	else {	/* Set up contour intervals automatically from Ctrl->C.interval and Ctrl->A.interval */
		min = floor (header.z_min / Ctrl->C.interval) * Ctrl->C.interval; if (!GMT_z_periodic && min < header.z_min) min += Ctrl->C.interval;
		max = ceil (header.z_max / Ctrl->C.interval) * Ctrl->C.interval; if (max > header.z_max) max -= Ctrl->C.interval;
		for (c = irint (min/Ctrl->C.interval), n_contours = 0; c <= irint (max/Ctrl->C.interval); c++, n_contours++) {
			if (n_contours == n_alloc) {
				(void)GMT_alloc_memory2 ((void **)&contour, (void **)&cont_angle, n_contours, n_alloc, sizeof (double), GMT_program);
				n_alloc = GMT_alloc_memory2 ((void **)&cont_type, (void **)&cont_do_tick, n_contours, n_alloc, sizeof (char), GMT_program);
			}
			contour[n_contours] = c * Ctrl->C.interval;
			if (Ctrl->contour.annot && (contour[n_contours] - aval) > GMT_SMALL) aval += Ctrl->A.interval;
			cont_type[n_contours] = (fabs (contour[n_contours] - aval) < GMT_SMALL) ? 'A' : 'C';
			cont_angle[n_contours] = (Ctrl->contour.angle_type == 2) ? Ctrl->contour.label_angle : GMT_d_NaN;
			cont_do_tick[n_contours] = (char)Ctrl->T.active;
		}
	}
	if (GMT_z_periodic && n_contours > 1 && fabs (contour[n_contours-1] - contour[0] - 360.0) < GMT_SMALL) {	/* Get rid of redundant contour */
		n_contours--;
	}

	if (n_contours == 0) {	/* No contours within range of data */
		if (gmtdefs.verbose) fprintf (stderr, "%s: Warning: No contours found\n", GMT_program);
		GMT_plotinit (argc, argv);
		GMT_map_basemap ();
		GMT_plotend ();
		GMT_free ((void *)contour);
		GMT_free ((void *)cont_type);
		GMT_free ((void *)cont_angle);
		GMT_free ((void *)cont_do_tick);
		GMT_free ((void *)grd_original);
		Free_grdcontour_Ctrl (Ctrl);	/* Deallocate control structure */
		GMT_end (argc, argv);
		exit (EXIT_SUCCESS);
	}
	
	/* OK, now we know we have contouring to do */
	
	(void)GMT_alloc_memory2 ((void **)&contour, (void **)&cont_angle, 0, n_contours, sizeof (double), GMT_program);
	n_alloc = GMT_alloc_memory2 ((void **)&cont_type, (void **)&cont_do_tick, 0, n_contours, sizeof (char), GMT_program);

	GMT_grd_minmax (grd_original, &header, xyz);
	if (GMT_contlabel_prep (&Ctrl->contour, xyz)) exit (EXIT_FAILURE);	/* Prep for crossing lines, if any */

	for (i = 0; i < 3; i++) GMT_io.out_col_type[i] = GMT_io.in_col_type[i];	/* Used if -D is set */
	
	/* Because we are doing single-precision, we cannot subtract incrementally but must start with the
	 * original grid values and subtract the current contour value. */
	 
	grd = (float *) GMT_memory (VNULL, (size_t)nm, sizeof (float), GMT_program);
	n_edges = header.ny * (int )ceil (header.nx / 16.0);
	edge = (GMT_LONG *) GMT_memory (VNULL, (size_t)n_edges, sizeof (GMT_LONG), GMT_program);

#ifdef EXPERIMENT
	wrap_check = (GMT_io.in_col_type[0] == GMT_IS_LON && (fabs (header.x_max - header.x_min - 360.0) < GMT_SMALL));
#endif
	
	GMT_plotinit (argc, argv);
	if (project_info.three_D) ps_transrotate (-z_project.xmin, -z_project.ymin, 0.0);

	GMT_map_clip_on (GMT_no_rgb, 3);

	for (c = 0; c < n_contours; c++) {	/* For each contour value cval */

		if (Ctrl->L.active && (contour[c] < Ctrl->L.low || contour[c] > Ctrl->L.high)) continue;	/* Outside desired range */

		/* Reset markers and set up new zero-contour*/

		cval = contour[c];
		if (gmtdefs.verbose) fprintf (stderr, "%s: Tracing the %g contour\n", GMT_program, cval);

		/* New approach to avoid round-off */

		for (i = 0; i < nm; i++) {
			grd[i] = grd_original[i] - (float)cval;		/* If there are NaNs they will remain NaNs */
			if (grd[i] == 0.0) grd[i] += (float)small;	/* There will be no actual zero-values, just -ve and +ve values */
		}

		section = 0;
		id = (cont_type[c] == 'A' || cont_type[c] == 'a') ? 1 : 0;

		Ctrl->contour.line_pen = Ctrl->W.pen[id];	/* Load current pen into contour structure */
		if (Ctrl->W.color) {	/* Override pen & label color according to cpt file */
			GMT_get_rgb_from_z (cval, rgb);
			memcpy ((void *)&Ctrl->contour.line_pen.rgb, (void *)rgb, (size_t)(3 * sizeof (int)));
			if (!Ctrl->contour.got_font_rgb && Ctrl->contour.curved_text) memcpy ((void *)&Ctrl->contour.font_rgb, (void *)rgb, (size_t)(3 * sizeof (int)));
		}

		n_alloc = 0;
		begin = TRUE;
		while ((n = GMT_contours (grd, &header, Ctrl->S.value, gmtdefs.interpolant, Ctrl->F.value, edge, &begin, &x, &y)) > 0) {

			if (fabs (x[0] - x[n-1]) < small_x && fabs (y[0] - y[n-1]) < small_y) {
				interior = 3;
				x[n-1] = x[0];	y[n-1] = y[0];	/* Ensure exact closure */
			}
#ifdef EXPERIMENT	
			else if (wrap_check) {
				if (fabs (x[0] - header.x_min) < GMT_SMALL && fabs (x[n-1] - header.x_min) < GMT_SMALL) interior = 2;
				else if (fabs (x[0] - header.x_max) < GMT_SMALL && fabs (x[n-1] - header.x_max) < GMT_SMALL) interior = 1;
			}
#endif
			else
				interior = FALSE;

			if (!interior || n >= Ctrl->Q.min) {
				if (Ctrl->D.active) GMT_dump_contour (x, y, n, cval, section, interior, Ctrl->D.file);
				if ((nn = GMT_clip_to_map (x, y, n, &xp, &yp))) {	/* Lines inside the region */
					/* From here on, xp/yp are map inches */
					if (cont_type[c] == 'A' || cont_type[c] == 'a') {	/* Annotated contours */
						GMT_get_format (cval, Ctrl->contour.unit, CNULL, format);
						sprintf (cont_label, format, cval);
					}
					else
						cont_label[0] = (char)0;
					if (cont_do_tick[c] && interior) {
						if (n_save == n_alloc) n_alloc = GMT_alloc_memory ((void **)&save, n_save, n_alloc, sizeof (struct SAVE), GMT_program);
						(void)GMT_alloc_memory2 ((void **)&save[n_save].x, (void **)&save[n_save].y, nn, 0, sizeof (double), GMT_program);
						memcpy ((void *)save[n_save].x, (void *)xp, (size_t)(nn * sizeof (double)));
						memcpy ((void *)save[n_save].y, (void *)yp, (size_t)(nn * sizeof (double)));
						save[n_save].n = nn;
						memcpy ((void *)&save[n_save].pen, (void *)&Ctrl->W.pen[id], sizeof (struct GMT_PEN));
						save[n_save].do_it = interior;
						save[n_save].cval = cval;
						n_save++;
					}
					GMT_hold_contour (&xp, &yp, nn, cval, cont_label, cont_type[c], cont_angle[c], interior, &Ctrl->contour);
					GMT_free ((void *)xp);
					GMT_free ((void *)yp);
				}
				section++;
			}
			GMT_free ((void *)x);
			GMT_free ((void *)y);
		}
	}

	if (Ctrl->W.pen[0].texture || Ctrl->W.pen[1].texture) ps_setdash (CNULL, 0);

	if (Ctrl->T.active && n_save) {
		(void)GMT_alloc_memory ((void **)&save, 0, n_save, sizeof (struct SAVE), GMT_program);

		sort_and_plot_ticks (save, n_save, &header, grd_original, Ctrl->T.spacing, Ctrl->T.length, Ctrl->T.low, Ctrl->T.high, Ctrl->T.label, Ctrl->T.L_label, Ctrl->T.H_label);
		for (i = 0; i < n_save; i++) {
			GMT_free ((void *)save[i].x);
			GMT_free ((void *)save[i].y);
		}
		GMT_free ((void *)save);
	}

	
	if (Ctrl->contour.hill_label) adjust_hill_label (&Ctrl->contour, &header, grd_original);	/* Must possibly adjust label angles so that label is readable when following contours */
	
	GMT_contlabel_plot (&Ctrl->contour);
	GMT_contlabel_free (&Ctrl->contour);

	GMT_map_clip_off ();

	GMT_map_basemap ();

	ps_setpaint (gmtdefs.background_rgb);

	if (project_info.three_D) ps_rotatetrans (z_project.xmin, z_project.ymin, 0.0);

	GMT_plotend ();

	GMT_free ((void *)grd);
	GMT_free ((void *)grd_original);
	GMT_free ((void *)edge);
	GMT_free ((void *)contour);
	GMT_free ((void *)cont_type);
	GMT_free ((void *)cont_angle);
	GMT_free ((void *)cont_do_tick);

	if (gmtdefs.verbose) fprintf (stderr, "%s: Done!\n", GMT_program);

	Free_grdcontour_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void sort_and_plot_ticks (struct SAVE *save, GMT_LONG n, struct GRD_HEADER *h, float *grd, double tick_gap, double tick_length, GMT_LONG tick_low, GMT_LONG tick_high, GMT_LONG tick_label, char Llabel, char Hlabel)
{
	GMT_LONG np, i, j, inside, ix, iy, done, n_ticks, stop;
	double x0, y0, add, sign, dx, dy, x_back, y_back, x_front, y_front;
	double xmin, xmax, ymin, ymax, x_mean = 0.0, y_mean = 0.0, inc, dist, a, small, this_lon, this_lat, L, sa, ca, *s;
	char txt[2][2];
#ifdef EXPERIMENT	
	GMT_LONG k, m, reverse;
#endif
	/* The x/y coordinates in SAVE are now all projected to map inches */

	small = 0.1 * project_info.xmax / h->nx;
	if (tick_label) {
		txt[0][0] = Llabel;	txt[0][1] = 0;
		txt[1][0] = Hlabel;	txt[1][1] = 0;
	}

	for (i = 0; i < n; i++) {	/* Look for tiny closed contours around a single node == cval */
		for (j = 1, L = 0.0; j < save[i].n; j++) L += hypot (save[i].x[j]-save[i].x[j-1], save[i].y[j]-save[i].y[j-1]);
		if (L < GMT_SMALL) {	/* This is a tiny closed contour; Mark as invalid */
			save[i].do_it = FALSE;
			continue;
		}
	}

#ifdef EXPERIMENT	
	/* Reconnect closed contours split at a periodic boundary */
	
	for (i = 0; i < n; i++) {	/* For all or these "closed" contours */
		if (save[i].do_it == 3) continue;	/* This one is fully closed already */
		for (j = i + 1; j < n; j++) {	/* Compare to all other "closed" contours */
			if (save[j].do_it == 3) continue;	/* This one is also fully closed */
			if (save[i].cval != save[j].cval) continue;	/* Not the same contour level */
			if (!(save[i].do_it + save[j].do_it) == 3) continue;	/* Not west + east split contours */
			/* OK, found a matching pair; now see if they are projected to be next to each other */
			np = save[i].n;
			if (fabs (save[i].x[0] - save[j].x[0]) < SMALL && fabs (save[i].y[0] - save[j].y[0]) < SMALL)
				reverse = TRUE;
			else if (fabs (save[i].x[np-1] - save[j].x[0]) < SMALL && fabs (save[i].y[np-1] - save[j].y[0]) < SMALL)
				reverse = FALSE;
			else
				continue;	/* No match */
			/* Here we append west to east segment */
			save[i].x = (double *) GMT_memory ((void *)save[i].x, (size_t)(np+save[j].n), sizeof (double), GMT_program);
			save[i].y = (double *) GMT_memory ((void *)save[i].y, (size_t)(np+save[j].n), sizeof (double), GMT_program);
			if (reverse) {
				for (k = np, m = save[j].n-1; m >= 0; k++, m--) {
					save[i].x[k] = save[j].x[m];
					save[i].y[k] = save[j].y[m];
				}
			}
			else {
				memcpy ((void *)&save[i].x[np], (void *)save[j].x, save[j].n * sizeof (double));
				memcpy ((void *)&save[i].y[np], (void *)save[j].y, save[j].n * sizeof (double));
			}
			save[i].n += save[j].n;
			save[i].do_it = 3;
			save[j].do_it = FALSE;
		}
	}
#endif
	
	for (i = 0; i < n; i++) {	/* Mark polygons that have other polygons inside them */
		np = save[i].n;
		for (j = 0; save[i].do_it == 3 && j < n; j++) {
			if (!save[j].do_it) continue;	/* No point checking contours that either have others inside or are invalid */
			inside = GMT_non_zero_winding (save[j].x[0], save[j].y[0], save[i].x, save[i].y, (GMT_LONG)np);
			if (inside == 2) save[i].do_it = FALSE;
		}
	}

	/* Here, only the polygons that are innermost (containing the local max/min, will have do_it = TRUE */

	for (i = 0; i < n; i++) {
		if (!save[i].do_it) continue;
		np = save[i].n;

		/* Loop around the contour and get min/max original x,y (longitude, latitude) coordinates */

		GMT_xy_to_geo (&xmax, &ymax, save[i].x[0], save[i].y[0]);
		xmin = xmax;
		ymin = ymax;
		for (j = 1; j < np; j++) {
			GMT_xy_to_geo (&this_lon, &this_lat, save[i].x[j], save[i].y[j]);
			xmin = MIN (xmin, this_lon);
			xmax = MAX (xmax, this_lon);
			ymin = MIN (ymin, this_lat);
			ymax = MAX (ymax, this_lat);
		}

		/* Pick the mid-latitude and march along that line from east to west */
		
		this_lat = 0.5 * (ymax + ymin);	/* Mid latitude, probably not exactly on a y-node */
		iy = MIN (h->ny - 1, GMT_y_to_j (this_lat, h->y_min, h->y_inc, h->xy_off, h->ny));	/* Get closest j-row */
		this_lat = GMT_j_to_y (iy, h->y_min, h->y_max, h->y_inc, h->xy_off, h->ny);		/* Get its matching latitude */
		ix = GMT_x_to_i (xmin, h->x_min, h->x_inc, h->xy_off, h->nx);	/* Westernmost point */
		stop  = GMT_x_to_i (xmax, h->x_min, h->x_inc, h->xy_off, h->nx);	/* Eastermost point */
		done = FALSE;
		while (!done && ix <= stop) {
			this_lon = GMT_i_to_x (ix, h->x_min, h->x_max, h->x_inc, h->xy_off, h->nx);	/* Current longitude */
			GMT_geo_to_xy (this_lon, this_lat, &x0, &y0);
			inside = GMT_non_zero_winding (x0, y0, save[i].x, save[i].y, (GMT_LONG)np);
			if (inside == 2)	/* OK, this point is inside */
				done = TRUE;
			else
				ix++;
		}
		if (!done) continue;	/* Failed to find an inside point */
		save[i].high = (grd[iy*h->nx+ix] > save[i].cval);

		if (save[i].high && !tick_high) continue;	/* Do not tick highs */
		if (!save[i].high && !tick_low) continue;	/* Do not tick lows */

		if (tick_label) {
			x_mean = save[i].x[0];
			y_mean = save[i].y[0];
		}
		s = (double *) GMT_memory (VNULL, (size_t)np, sizeof (double), GMT_program);
		for (j = 1, s[0] = 0.0; j < np; j++) {
			s[j] = s[j-1] + hypot (save[i].x[j]-save[i].x[j-1], save[i].y[j]-save[i].y[j-1]);
			if (tick_label) {
				x_mean += save[i].x[j];
				y_mean += save[i].y[j];
			}
		}
		if (s[np-1] < GRDCONTOUR_MIN_LENGTH) {
			GMT_free ((void *)s);
			continue;
		}
		dx = save[i].x[1] - save[i].x[0];
		dy = save[i].y[1] - save[i].y[0];
		a = atan2 (dy, dx) + M_PI_2;
		sincos (a, &sa, &ca);
		x0 = 0.5 * (save[i].x[0] + save[i].x[1]) + small * ca;
		y0 = 0.5 * (save[i].y[0] + save[i].y[1]) + small * sa;
		inside = GMT_non_zero_winding (x0, y0, save[i].x, save[i].y, (GMT_LONG)np);
		sign = (inside == 2) ? 1.0 : -1.0;

		n_ticks = (int)(s[np-1] / tick_gap);
		if (n_ticks == 0) {
			GMT_free ((void *)s);
			continue;
		}
		GMT_setpen (&save[i].pen);
		inc = s[np-1] / n_ticks;
		if (tick_label) {
			x_mean /= np;
			y_mean /= np;
			GMT_text3D (x_mean, y_mean, project_info.z_level, gmtdefs.annot_font_size[0], gmtdefs.annot_font[0], txt[save[i].high], 0.0, 6, 0);
		}

		x_back = save[i].x[np-1];
		y_back = save[i].y[np-1];
		dist = 0.0;
		j = 0;
		add = sign * ((save[i].high) ? -1.0 : 1.0);
		while (j < np-1) {
			x_front = save[i].x[j+1];
			y_front = save[i].y[j+1];
			if (s[j] >= dist) {	/* Time for tick */
				dx = x_front - x_back;
				dy = y_front - y_back;
				a = atan2 (dy, dx) + add * M_PI_2;
				sincos (a, &sa, &ca);
				if (project_info.three_D) {
					GMT_xy_do_z_to_xy (save[i].x[j], save[i].y[j], project_info.z_level, &x0, &y0);
					ps_plot (x0, y0, PSL_PEN_MOVE);
					GMT_xy_do_z_to_xy (save[i].x[j] + tick_length * ca, save[i].y[j] + tick_length * sa, project_info.z_level, &x0, &y0);
					ps_plot (x0, y0, PSL_PEN_DRAW_AND_STROKE);
				}
				else {
					ps_plot (save[i].x[j], save[i].y[j], PSL_PEN_MOVE);
					ps_plotr (tick_length * ca, tick_length * sa, PSL_PEN_DRAW_AND_STROKE);
				}
				dist += inc;
			}
			x_back = x_front;
			y_back = y_front;
			j++;
		}
		GMT_free ((void *)s);
	}
}

void GMT_grd_minmax (float *a, struct GRD_HEADER *h, double xyz[2][3])
{
	GMT_LONG i, i_min, i_max;
	float z_min, z_max;
	double half;

	z_min = FLT_MAX;	z_max = -FLT_MAX;
	i_min = i_max = 0;
	half = (h->node_offset) ? 0.5 : 0.0;
	for (i = 0; i < h->nx * h->ny; i++) {
		if (GMT_is_fnan (a[i])) continue;
		if (a[i] < z_min) {
			z_min = a[i];
			i_min = i;
		}
		if (a[i] > z_max) {
			z_max = a[i];
			i_max = i;
		}
	}
	xyz[0][0] = h->x_min + (i_min%h->nx + half) * h->x_inc;
	xyz[0][1] = h->y_max - (i_min/h->nx + half) * h->y_inc;
	xyz[1][0] = h->x_min + (i_max%h->nx + half) * h->x_inc;
	xyz[1][1] = h->y_max - (i_max/h->nx + half) * h->y_inc;
	xyz[0][2] = z_min;	xyz[1][2] = z_max;
}

void adjust_hill_label (struct GMT_CONTOUR *G, struct GRD_HEADER *h, float *z)
{	
	GMT_LONG i, k, i0, j0;
	double nx, ny, x_on, y_on, x_node, y_node, x_node_p, y_node_p, dx, dy, dz, dot, angle;
	struct GMT_CONTOUR_LINE *C;
	
	for (i = 0; i < G->n_segments; i++) {
		C = G->segment[i];	/* Pointer to current segment */
		for (k = 0; k < C->n_labels; k++) {
			GMT_xy_to_geo (&x_on, &y_on, C->L[k].x, C->L[k].y);	/* Retrieve original coordinates */
			j0 = GMT_y_to_j (y_on, h->y_min, h->y_inc, h->xy_off, h->ny);
			if (j0 < 0 || j0 >= h->ny) continue;	/* Somehow, outside y range */
			while (GMT_io.in_col_type[0] == GMT_IS_LON && x_on < h->x_min) x_on += 360.0;
			while (GMT_io.in_col_type[0] == GMT_IS_LON && x_on > h->x_max) x_on -= 360.0;
			i0 = GMT_x_to_i (x_on, h->x_min, h->x_inc, h->xy_off, h->nx);
			if (i0 < 0 || i0 >= h->nx) continue;	/* Somehow, outside x range */
			angle = fmod (2.0 * C->L[k].angle, 360.0) * 0.5;	/* 0-180 range */
			if (angle > 90.0) angle -= 180.0;
			sincosd (angle + 90, &ny, &nx);	/* Coordinate of normal to label line */
			x_node = GMT_i_to_x (i0, h->x_min, h->x_max, h->x_inc, h->xy_off, h->nx);
			y_node = GMT_j_to_y (j0, h->y_min, h->y_max, h->y_inc, h->xy_off, h->ny);
			GMT_geo_to_xy (x_node, y_node, &x_node_p, &y_node_p);	/* Projected coordinates of nearest node point */
			dx = x_node_p - C->L[k].x;
			dy = y_node_p - C->L[k].y;
			if (GMT_IS_ZERO (hypot (dx, dy))) {
				fprintf (stderr, "%s: Unable to adjust hill label contour orientation (node point on contour)\n", GMT_program);
				continue;
			}
			dz = z[j0*h->nx+i0] - C->z;
			if (GMT_IS_ZERO (dz)) {
				fprintf (stderr, "%s: Unable to adjust hill label contour orientation (node value = contour value)\n", GMT_program);
				continue;
			}
			dot = dx * nx + dy * ny;	/* Dot product of n and vector from contour to node. +ve if on same side of contour line */
			if (irint (copysign (1.0, dot * dz)) != G->hill_label)
				C->L[k].angle += 180.0;	/* Must turn upside-down */
		}
	}
}
	
void *New_grdcontour_Ctrl () {	/* Allocate and initialize a new control structure */
	struct GRDCONTOUR_CTRL *C;
	
	C = (struct GRDCONTOUR_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GRDCONTOUR_CTRL), "New_grdcontour_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	
	GMT_contlabel_init (&C->contour, 1);
	C->D.file = strdup ("contour");
	C->E.azimuth = 180.0;
	C->E.elevation = 90.0;
	C->L.low = -DBL_MAX;
	C->L.high = DBL_MAX;
	if (gmtdefs.measure_unit == GMT_CM) {	/* Default is 0.5cm and 1 mm */
		C->T.spacing = 0.5 / 2.54;
		C->T.length = 0.1 / 2.54;
	}
	else {					/* Default in inches is 0.2 inches and 0.04 inches */
		C->T.spacing = 0.2;
		C->T.length = 0.04;
	}
	GMT_init_pen (&C->W.pen[0], GMT_PENWIDTH);	/* Normal pen for contours and 3*thicker for annotated contours */
	GMT_init_pen (&C->W.pen[1], 3.0 * GMT_PENWIDTH);
	C->Z.scale = 1.0;
	
	return ((void *)C);
}

void Free_grdcontour_Ctrl (struct GRDCONTOUR_CTRL *C) {	/* Deallocate control structure */
	if (C->C.file) free ((void *)C->C.file);	
	if (C->D.file) free ((void *)C->D.file);	
	GMT_free ((void *)C);	
}
