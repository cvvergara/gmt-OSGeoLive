/*--------------------------------------------------------------------
 *	$Id: psxyz.c 10173 2014-01-01 09:52:34Z pwessel $
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
 * psxyz will read <x,y,z> triplets from GMT_stdin and plot
 * the positions in 3-D using symbols selected by the user. If no
 * symbol size is specified, psxyz expects <x,y,z,size> as input, where
 * size is the symbol size.  Several map projections are supported.
 * For linear projections a 3-D basemap is provided.  All symbols are
 * projected onto the xy plane (so that circles becomes ellipses) except
 * GMT_SYMBOL_COLUMN and GMT_SYMBOL_CUBE which are fully 3-D.
 * PostScript code is then written to GMT_stdout.
 *
 * Author:    Paul Wessel
 * Date:      22-AUG-2000
 * Version:   4
 *
 */

#include "gmt.h"
#include "pslib.h"

struct PSXYZ_CTRL {
	struct A {	/* -A[step] {NOT IMPLEMENTED YET} */
		GMT_LONG active;
		double step;
	} A;
	struct C {	/* -C<cpt> */
		GMT_LONG active;
		char *file;
	} C;
	struct D {	/* -D<dx>/<dy>[/<dz>] */
		GMT_LONG active;
		double dx, dy, dz;
	} D;
	struct E {	/* -E<azimuth/elevation> */
		GMT_LONG active;
		double azimuth, elevation;
	} E;
	struct G {	/* -G<fill> */
		GMT_LONG active;
		struct GMT_FILL fill;
	} G;
	struct I {	/* -I<intensity> */
		GMT_LONG active;
		double value;
	} I;
	struct L {	/* -L */
		GMT_LONG active;
	} L;
	struct N {	/* -N */
		GMT_LONG active;
	} N;
	struct Q {	/* -Q */
		GMT_LONG active;
	} Q;
	struct S {	/* -S */
		GMT_LONG active;
		char *arg;
	} S;
	struct W {	/* -W<pen> */
		GMT_LONG active;
		GMT_LONG mode;	/* 0 = normal, 1 = -C applies to pen color only, 2 = -C applies to symbol fill & pen color */
		struct GMT_PEN pen;
	} W;
	struct Z {	/* -Z<z_level> */
		GMT_LONG active;
		double level;
	} Z;
};

struct PSXYZ_DATA1 {
	char *string;
	GMT_LONG symbol, flag;
	float x_size, y_size, z_size;
	double lon, lat;
	double x, y, z, value, dist;
	struct GMT_FILL f;
	struct GMT_PEN p;
	GMT_LONG outline;
};

struct PSXYZ_DATA2 {
	double x, y, z;
};

int main (int argc, char **argv)
{
	GMT_LONG error = FALSE, nofile = TRUE, polygon = FALSE, done, penset_OK = TRUE;
	GMT_LONG get_rgb = FALSE, read_symbol, clip_set = FALSE, fill_active, bar_clip;
	GMT_LONG default_outline = FALSE, outline_active, save_u = FALSE;

	char line[BUFSIZ], txt_a[GMT_LONG_TEXT], txt_b[GMT_LONG_TEXT], txt_c[GMT_LONG_TEXT];

	GMT_LONG i, n, n_alloc = 0, n_total_read = 0, j, n_files = 0, pos2x, pos2y;
	GMT_LONG fno, n_cols_start = 3, n_args, n_fields, this_outline;
	GMT_LONG three, four, five, n_required = 3, n_expected = 0, change;
	
	int rgb[9];

	double west = 0.0, east = 0.0, south = 0.0, north = 0.0, dim[3], DX, DY;
	double lux[3], tmp, *in = NULL, x_1, x_2, y_1, y_2, xxx, yyy, font_size;
	double x2, y2, v_w, h_l, h_w, dx, dy, xx[2], yy[2], s, c, size;

	FILE *fp = NULL;

	struct GMT_PEN default_pen, current_pen;
	struct GMT_FILL default_fill, current_fill;
	struct GMT_SYMBOL S;
	struct PSXYZ_DATA1 *data1 = NULL;
	struct PSXYZ_DATA2 *data2 = NULL;
	struct PSXYZ_CTRL *Ctrl = NULL;

	void column3D (double x, double y, double z, double base, double x_size, double y_size, int *rgb, GMT_LONG outline);
	void cube3D (double x, double y, double z, double x_size, double y_size, int *rgb, GMT_LONG outline);
	void barx3D (double x, double y, double z, double base, double size, struct GMT_FILL *f, GMT_LONG outline);
	void bary3D (double x, double y, double z, double base, double size, struct GMT_FILL *f, GMT_LONG outline);
	void sort_on_distance (struct PSXYZ_DATA1 *data, GMT_LONG n);
	void *New_psxyz_Ctrl (), Free_psxyz_Ctrl (struct PSXYZ_CTRL *C);
	
	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct PSXYZ_CTRL *)New_psxyz_Ctrl ();	/* Allocate and initialize a new control structure */

	GMT_3D_mode = 1;	/* Only do background axis first; do foreground at end */

	for (i = 0; i < 9; i++) rgb[i] = -1;

	/* Initialize GMT_SYMBOL structure */

	memset ((void *)&S, 0, sizeof (struct GMT_SYMBOL));
	GMT_contlabel_init (&S.G, 0);

	S.v_width = 0.03;
	S.h_length = 0.12;
	S.h_width = 0.1;
	S.base = GMT_d_NaN;
	
	lux[0] = lux[1] = lux[2] = 0.0;	/* Default is no shading */

	if ((S.u = gmtdefs.measure_unit) == GMT_CM) {
		S.v_width = 0.075 / 2.54;	S.h_length = 0.3 / 2.54;
		S.h_width = 0.25 / 2.54;
	}
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
				case 'f':
				case 'g':
				case 'm':
				case ':':
				case '\0':
					error += (GMT_LONG)GMT_parse_common_options (argv[i], &west, &east, &south, &north);
					break;
		
				/* Supplemental parameters */
	
				case 'C':
					Ctrl->C.file = strdup (&argv[i][2]);
					Ctrl->C.active = TRUE;
					break;
				case 'D':
					if ((n = sscanf (&argv[i][2], "%[^/]/%[^/]/%s", txt_a, txt_b, txt_c)) < 2) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -D option: Give x and y [and z] offsets\n", GMT_program);
						error++;
					}
					else {
						Ctrl->D.dx = GMT_convert_units (txt_a, GMT_INCH);
						Ctrl->D.dy = GMT_convert_units (txt_b, GMT_INCH);
						if (n == 3) Ctrl->D.dz = GMT_convert_units (txt_c, GMT_INCH);
						Ctrl->D.active = TRUE;
					}
					break;
				case 'E':
					Ctrl->E.active = TRUE;
					error += (GMT_LONG)GMT_get_proj3D (&argv[i][2], &Ctrl->E.azimuth, &Ctrl->E.elevation);
					break;
				case 'G':		/* Set color for symbol or polygon */
					Ctrl->G.active = TRUE;
					if (!argv[i][2] || (argv[i][2] && GMT_getfill (&argv[i][2], &Ctrl->G.fill))) {
						GMT_fill_syntax ('G', " ");
						error++;
					}
					break;
				case 'I':	/* Adjust symbol color via intensity */
					Ctrl->I.value = atof (&argv[i][2]);
					Ctrl->I.active = TRUE;
					break;
				case 'L':		/* Force closed polygons */
					Ctrl->L.active = TRUE;
					break;
                                case 'N':               /* Do not clip to map */
                                        Ctrl->N.active = TRUE;
                                        break;
				case 'Q':               /* Do not sort symbols based on distance */
                                        Ctrl->Q.active = TRUE;
                                        break;
				case 'S':		/* Get symbol [and size] */
                                        Ctrl->S.active = TRUE;
					Ctrl->S.arg = &argv[i][2];
					break;
				case 'W':		/* Set line attributes */
					Ctrl->W.active = TRUE;
					j = 2;
					if (argv[i][j] == '-') {Ctrl->W.mode = 1; j++;}
					if (argv[i][j] == '+') {Ctrl->W.mode = 2; j++;}
					if (argv[i][j] && GMT_getpen (&argv[i][j], &Ctrl->W.pen)) {
						GMT_pen_syntax ('W', " ");
						error++;
					}
					break;
				case 'Z':
					if (argv[i][2]) {
						Ctrl->Z.level = atof (&argv[i][2]);
						Ctrl->Z.active = TRUE;
					}
					break;

				default:		/* Options not recognized */
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else
			n_files++;
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {	/* Display usage */
		fprintf (stderr,"psxyz %s - Plot lines, polygons, and symbols in 3-D\n\n", GMT_VERSION);
		fprintf(stderr,"usage: psxyz <xyzfiles> %s %s [%s] [%s] [-C<cpt>]\n", GMT_J_OPT, GMT_Rgeoz_OPT, GMT_B_OPT, GMT_Jz_OPT);
		fprintf(stderr, "\t[-D<dx>/<dy>[/<dz>]] [%s] [-G<fill>] [%s] [-I<intens>] [-K] [-L]\n", GMT_E_OPT, GMT_Ho_OPT);
		fprintf(stderr, "\t[-N] [-O] [-P] [-Q] [-S<symbol><size>[/size_y]] [%s] [-V] [-W[+|-][<pen>]] [%s]\n", GMT_U_OPT, GMT_X_OPT);
		fprintf(stderr, "\t[%s] [%s] [%s] [%s] [%s] [%s] [%s]\n", GMT_Y_OPT, GMT_c_OPT, GMT_t_OPT, GMT_bi_OPT, GMT_f_OPT, GMT_g_OPT, GMT_mo_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\t<xyzfiles> is one or more files.  If none, read standard input.\n");
		GMT_explain_option ('j');
		GMT_explain_option ('Z');
		GMT_explain_option ('r');
		fprintf (stderr, "\n\tOPTIONS:\n");
		GMT_explain_option ('b');
		fprintf (stderr, "\t-C Use cpt-file to assign symbol colors based on t-value in 4th column (requires -S).\n");
		fprintf (stderr, "\t   Without -S, psxyz expects -m lines/polygons and looks for -Z<tval> strings in each multiheader.\n");
		fprintf (stderr, "\t-D Offset symbol or line positions by <dx>/<dy>[/<dy>] [no offset].\n");
		GMT_explain_option ('E');
		GMT_fill_syntax ('G', "Specify color or pattern [Default is no fill].");
		fprintf (stderr, "\t   If -G is specified but not -S, psxyz draws a filled polygon.\n");
		GMT_explain_option ('H');
		fprintf (stderr, "\t-I sets a constant intensity used to modulate the fill color (requires -C or -G).\n");
		GMT_explain_option ('K');
		fprintf (stderr, "\t-L Force closed polygons.\n");
		fprintf (stderr, "\t-N Do Not skip symbols that fall outside map border [Default will ignore those outside].\n");
		GMT_explain_option ('O');
		GMT_explain_option ('P');
		fprintf (stderr, "\t-Q Do NOT sort symbols based on distance to viewer before plotting.\n");
		fprintf (stderr, "\t-S to select symbol type and symbol size (in %s).  Choose between\n", GMT_unit_names[gmtdefs.measure_unit]);
		fprintf (stderr, "\t   -(xdash), +(plus), st(a)r, (b|B)ar, (c)ircle, (d)iamond, (e)llipse, (f)ront, octa(g)on, (h)exagon\n");
		fprintf (stderr, "\t   (i)nvtriangle, (j)rotated rectangle, (k)ustom, (%c)etter, pe(n)tagon, (m)athangle, c(o)lumn, (p)oint,\n", 'l');
		fprintf (stderr, "\t   (q)uoted line, (r)ect, (s)quare, (t)riangle, c(u)be, (v)ector, (w)edge, (x)cross, (y)dash, (z)dash.\n");
		fprintf (stderr, "\t   If no size is specified, psxyz expects the 4th column to have sizes.\n");
		fprintf (stderr, "\t   If no symbol is specified, psxyz expects the last column to have symbol code.\n");
		fprintf (stderr, "\t   [Note: if -C is selected then 4th means 5th column, etc.]\n");
		fprintf (stderr, "\t   column and cube are true 3-D objects (give size as xsize/ysize), the rest is 2-D perspective only.\n");
		fprintf (stderr, "\t   By default, column and cube are shaded; use O and U to disable 3-D illumination.\n");
		fprintf (stderr, "\t   Symbols A, C, D, F, H, I, N, S, T are adjusted to have same area of a circle of given size.\n");
		fprintf (stderr, "\t   Bar (or Column): Append b<base> to give the y- (or z-) value of the base [Default = 0 (1 for log-scales)].\n");
		fprintf (stderr, "\t     Use -SB for horizontal bars; then base refers to x.\n");
		fprintf (stderr, "\t   Ellipses: the direction, major, and minor axis must be in input columns 4, 5 and 6.\n");
		fprintf (stderr, "\t     If -SE rather than -Se is selected, psxyz will expect azimuth, and axes in km\n");
		fprintf (stderr, "\t     and convert azimuths based on the chosen map projection.  If projection is linear then\n");
		fprintf (stderr, "\t     we scale the minor/major axes by the map scale.\n");
		fprintf (stderr, "\t   Rotatable Rectangle: Direction, x-dimension, and y-dimension are expected in columns 4, 5, and 6.\n");
		fprintf (stderr, "\t     If -SJ rather than -Sj is selected, psxyz will expect azimuth, and dimensions in km\n");
		fprintf (stderr, "\t     and convert azimuths based on the chosen map projection.  If projection is linear then\n");
		fprintf (stderr, "\t     we scale the rectangle dimension by the map scale.\n");
		fprintf (stderr, "\t   Fronts: Give tickgap/ticklen[dir][type][:offset], where\n");
		fprintf (stderr, "\t     dir    = Plot symbol to the l(eft) or r(ight) of front [Default is centered].\n");
		fprintf (stderr, "\t     type   =  b(ox), c(ircle), f(ault), s(lip), t(riangle) [Default is f].\n");
		fprintf (stderr, "\t       box      : square when centered, half-square otherwise.\n"); 
		fprintf (stderr, "\t       circle   : full when centered, half-circle otherwise.\n"); 
		fprintf (stderr, "\t       fault    : centered cross-tick or tick only in the <dir> direction.\n"); 
		fprintf (stderr, "\t       slip     : left-lateral or right-lateral strike-slip arrows (centered is not defined).\n"); 
		fprintf (stderr, "\t       triangle : diagonal square when centered, directed triangle otherwise.\n"); 
		fprintf (stderr, "\t     offset = start plotting symbols when along-track distance equals offset [0].\n"); 
		fprintf (stderr, "\t   Kustom: append <symbolname> immediately after 'k'; this will look for <symbolname>.def in\n");
		fprintf (stderr, "\t     the current directory, in $GMT_USERDIR, or in $GMT_SHAREDIR (in that order).\n");
		GMT_list_custom_symbols ();
		fprintf (stderr, "\t   Letter: append /<string> after symbol size, and optionally %%<font>.\n");
		fprintf (stderr, "\t   Mathangle: the start and stop directions of the math angle must be in input columns 4, 5.\n");
		fprintf (stderr, "\t     Use -Smf for arrow at first angle, -Sml for last angle, -Smb for both [none].\n");
		fprintf (stderr, "\t   Quoted line: Give [d|f|n|l|x]<info>[:<labelinfo>]. z must be constant for each segment.\n");
		fprintf (stderr, "\t     <code><info> controls the placement of labels along the lines.  Choose among five algorithms:\n");
		GMT_cont_syntax (7, 1);
		fprintf (stderr, "\t   <labelinfo> controls the specifics of the labels.  Append what you need:\n");
		GMT_label_syntax (7, 1);
		fprintf (stderr, "\t   Rectangles: the x- and y-dimension must be in input columns 4 and 5.\n");
		fprintf (stderr, "\t   Vectors: the direction and length must be in input columns 4 and 5.\n");
		fprintf (stderr, "\t     Furthermore, <size> means arrowwidth/headlength/headwidth [Default is %g/%g/%g].\n", S.v_width, S.h_length, S.h_width);
		fprintf (stderr, "\t     If -SV rather than -Sv is selected, psxyz will expect azimuth and length\n");
		fprintf (stderr, "\t     and convert azimuths based on the chosen map projection.\n");
		fprintf (stderr, "\t     Insert h(head), +(center), or t(ail) after -Sv|V to justify vector w.r.t. input (x,y).\n");
		fprintf (stderr, "\t     Insert s(egment) if (x,y) is tail and columns 3 and 4 holds head (x,y).\n");
		fprintf (stderr, "\t     Upper case H, B, T, S will give double-headed vector [Default is t].\n");
		GMT_explain_option ('U');
		fprintf (stderr, "\t   Wedges: the start and stop directions of pie wedge must be in input columns 4 and 5.\n");
		fprintf (stderr, "\t     If -SW rather than -Sw is selected, psxy will expect two azimuths instead.\n");
		GMT_explain_option ('V');
		GMT_pen_syntax ('W', "sets pen attributes:");
		fprintf (stderr, "\t   Implicitly draws symbol outline with this pen.\n");
		fprintf (stderr, "\t   The leading + applies cpt color (-C) to both symbol fill and pen; - applies it to pen only.\n");
		GMT_explain_option ('X');
		GMT_explain_option ('c');
		GMT_explain_option (':');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf (stderr, "\t   Default is the required number of columns.\n");
		GMT_explain_option ('f');
		GMT_explain_option ('g');
		GMT_explain_option ('m');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	/* Check that the options selected are mutually consistent */

	if (Ctrl->S.active && GMT_parse_symbol_option (Ctrl->S.arg, &S, 1, TRUE)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -S option\n", GMT_program);
		error++;
	}
	if (!project_info.region_supplied) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify -R option\n", GMT_program);
		error++;
	}
	if (Ctrl->E.elevation <= 0.0 || Ctrl->E.elevation > 90.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -E option:  Elevation must be in 0-90 range\n", GMT_program);
		error++;
	}
	if (Ctrl->C.active && S.symbol == 0 && !GMT_io.multi_segments[GMT_IN]) {
		error++;
		fprintf (stderr, "%s: GMT SYNTAX ERROR -C option: Must also specify a symbol (see -S) or give polygon file with -m\n", GMT_program);
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
	if (GMT_io.binary[GMT_IN] && S.symbol == GMT_SYMBOL_NOT_SET) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have symbol information\n", GMT_program);
		error++;
	}
	if (Ctrl->W.active && Ctrl->W.mode && !Ctrl->C.active) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  -W option +|-<pen> requires the -C option\n", GMT_program);
		error++;
	}

	if (error) exit (EXIT_FAILURE);

	if (S.symbol > 0 && Ctrl->C.active) get_rgb = TRUE;
	read_symbol = (S.symbol == GMT_SYMBOL_NOT_SET);

	if (S.symbol == 0 && (Ctrl->G.active || Ctrl->L.active)) polygon = TRUE;
	if (S.symbol == GMT_SYMBOL_FRONT || S.symbol == GMT_SYMBOL_QUOTED_LINE) polygon = FALSE;

	z_project.view_azimuth = Ctrl->E.azimuth;
	z_project.view_elevation = Ctrl->E.elevation;
	default_pen = current_pen = Ctrl->W.pen;
	default_fill = current_fill = Ctrl->G.fill;
	default_outline = Ctrl->W.active;
	if (Ctrl->I.active) {
		GMT_illuminate (Ctrl->I.value, current_fill.rgb);
		GMT_illuminate (Ctrl->I.value, default_fill.rgb);
	}

	if (get_rgb) n_cols_start++;

	three = (get_rgb) ? 4 : 3;
	four = (get_rgb) ? 5 : 4;
	five = (get_rgb) ? 6 : 5;
	pos2x = three + gmtdefs.xy_toggle[GMT_IN];	/* Column with a 2nd longitude (for VECTORS with two sets of coordinates) */
	pos2y = four - gmtdefs.xy_toggle[GMT_IN];	/* Column with a 2nd latitude (for VECTORS with two sets of coordinates) */
	n_required = n_cols_start + S.n_required;

	for (j = n_cols_start; j < 10; j++) GMT_io.in_col_type[j] = GMT_IS_DIMENSION;		/* Since these may have units appended */
	for (j = 0; j < S.n_nondim; j++) GMT_io.in_col_type[S.nondim_col[j]+get_rgb] = GMT_IS_FLOAT;	/* Since these are angles, not dimensions */
	n_expected = (GMT_io.binary[GMT_IN]) ? GMT_io.ncol[GMT_IN] : n_required;
	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN]=  n_required;
	if (n_expected < n_required) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data must have at least %ld columns\n", GMT_program, n_required);
		error++;
	}

	if (GMT_io.binary[GMT_IN] && gmtdefs.verbose) {
		char *type[2] = {"double", "single"};
		fprintf (stderr, "%s: Expects %ld-column %s-precision binary data\n", GMT_program, GMT_io.ncol[GMT_IN], type[GMT_io.single_precision[GMT_IN]]);
	}

	if (Ctrl->C.active) GMT_read_cpt (Ctrl->C.file);

	if (n_files > 0)
		nofile = FALSE;
	else
		n_files = 1;
	n_args = (argc > 1) ? argc : 2;
	if (S.symbol == GMT_SYMBOL_QUOTED_LINE) {
		if (GMT_contlabel_prep (&S.G, NULL)) exit (EXIT_FAILURE);
		penset_OK = FALSE;	/* Since it is set in PSL */
	}

	GMT_err_fail (GMT_map_setup (west, east, south, north), "");

	if (S.u_set) {	/* When -Sc<unit> is given we temporarily reset the system unit to these units so conversions will work */
		save_u = gmtdefs.measure_unit;
		gmtdefs.measure_unit = S.u;
	}

	lux[0] = fabs (z_project.sin_az * z_project.cos_el);
	lux[1] = fabs (z_project.cos_az * z_project.cos_el);
	lux[2] = fabs (z_project.sin_el);
	tmp = MAX (lux[0], MAX (lux[1], lux[2]));
	for (i = 0; i < 3; i++) lux[i] = (lux[i] / tmp) - 0.5;
	
	if ((Ctrl->C.active || current_fill.rgb[0]) >= 0 && (S.symbol == GMT_SYMBOL_COLUMN || S.symbol == GMT_SYMBOL_CUBE)) {	/* Modify the color for each facet */
		for (i = 0; i < 3; i++) {
			rgb[3*i] = current_fill.rgb[0];
			rgb[3*i+1] = current_fill.rgb[1];
			rgb[3*i+2] = current_fill.rgb[2];
			if (S.shade3D) GMT_illuminate (lux[i], &rgb[3*i]);
		}
	}

	GMT_plotinit (argc, argv);

	ps_transrotate (-z_project.xmin, -z_project.ymin, 0.0);

	if (Ctrl->Z.active) project_info.z_level = Ctrl->Z.level;

	GMT_map_basemap ();

	if (project_info.z_pars[0] == 0.0 && !Ctrl->N.active) {
		GMT_map_clip_on (GMT_no_rgb, 3);
		clip_set = TRUE;
	}
	if (S.symbol == GMT_SYMBOL_ELLIPSE) Ctrl->N.active = TRUE;

	if (penset_OK) GMT_setpen (&current_pen);

	if (S.symbol == GMT_SYMBOL_TEXT && Ctrl->G.active && !Ctrl->W.active) ps_setpaint (current_fill.rgb);
	if (S.symbol == GMT_SYMBOL_TEXT) {
		font_size = S.size_x * 72.0;	/* To get points */
		ps_setfont (S.font_no);		/* Set the required font */
	}
	if (S.convert_angles && GMT_IS_MAPPING) S.convert_angles = 2;
	if ((S.symbol == GMT_SYMBOL_VECTOR || S.symbol == GMT_SYMBOL_VECTOR) && S.v_just == 3) {
		/* Reading 2nd coordinate so must set column types */
		GMT_io.in_col_type[pos2x] = GMT_io.in_col_type[GMT_X];
		GMT_io.in_col_type[pos2y] = GMT_io.in_col_type[GMT_Y];
	}
	fill_active = (GMT_LONG)Ctrl->G.active;	/* Make copies because we will change the values */
	outline_active =  (GMT_LONG)Ctrl->W.active;
	if (S.symbol > 0 && !outline_active && !fill_active && !get_rgb) outline_active = TRUE;		/* If no fill nor outline for symbols then turn outline on */
	bar_clip = ((S.symbol == GMT_SYMBOL_BARX || S.symbol == GMT_SYMBOL_BARY) && !Ctrl->N.active) ? TRUE : FALSE;
			/* Bars get clipped instead of skipped */
	if (Ctrl->D.active) {
		GMT_xyz_to_xy (Ctrl->D.dx, Ctrl->D.dy, Ctrl->D.dz, &DX, &DY);
		ps_transrotate (DX, DY, 0.0);	/* Shift plot a bit */
	}
	GMT_io.skip_if_NaN[GMT_Z] = TRUE;	/* Extend GMT NaN-handling to the z-coordinate */

	done = FALSE;
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

		GMT_io.file_no++;
		GMT_io.seg_no = -1;
		if (!nofile && gmtdefs.verbose) {
			fprintf (stderr, "%s: Working on file %s\n", GMT_program, argv[fno]);
			sprintf (line, "File: %s", argv[fno]);
			ps_comment (line);
		}
		if (GMT_io.io_header[GMT_IN]) {
			for (i = 0; i < GMT_io.n_header_recs; i++) {
				GMT_fgets (line, BUFSIZ, fp);
				if (i == 0 && S.G.label_type == 2) GMT_extract_label (&line[1], S.G.label);	/* Set initial header as potential label */
			}
		}

		if (S.symbol > 0) {	/* symbol part (not counting GMT_SYMBOL_FRONT and GMT_SYMBOL_QUOTED_LINE) */

			GMT_world_map = (S.symbol == GMT_SYMBOL_ELLIPSE && S.convert_angles) ? FALSE : TRUE;
			if (read_symbol) n_expected = GMT_MAX_COLUMNS;

			n = 0;
			while ((n_fields = GMT_input (fp, &n_expected, &in)) >= 0 && !GMT_REC_IS_EOF) {

				while (GMT_REC_IS_SEG_HEADER && !GMT_REC_IS_EOF) {	/* Skip segment headers */
					if (gmtdefs.verbose) ps_comment (GMT_io.segment_header);
					
					change = GMT_parse_multisegment_header (GMT_io.segment_header, Ctrl->C.active, &fill_active, &current_fill, &default_fill, &outline_active, &current_pen, &default_pen, default_outline);
					if (Ctrl->I.active) {
						GMT_illuminate (Ctrl->I.value, current_fill.rgb);
						GMT_illuminate (Ctrl->I.value, default_fill.rgb);
					}
					if (change & 4) GMT_setpen (&current_pen);

					if (read_symbol) n_expected = GMT_MAX_COLUMNS;
					n_fields = GMT_input (fp, &n_expected, &in);
				}
				if (GMT_REC_IS_EOF) continue;		/* At EOF */

				n_total_read++;
				if (GMT_REC_IS_ERROR && !read_symbol) {
					fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld (skipped)\n", GMT_program, n_fields, n_expected, n_total_read);
					continue;
				}

				if (read_symbol) {
					i = strlen (GMT_io.current_record) - 1;
					GMT_io.current_record[i--] = 0;	/* Chop off the \n */
					while (GMT_io.current_record[i] && !strchr (" ,\t", (int)GMT_io.current_record[i])) i--;
					GMT_parse_symbol_option (&GMT_io.current_record[i+1], &S, 1, FALSE);
					n_required = n_cols_start + S.n_required;
				}

				/* If read_symbol is TRUE we must unto the dimension scaling for items that should have no scaling, e.g., angles */
				
				if (read_symbol) for (j = 0; j < S.n_nondim; j++) in[S.nondim_col[j]+get_rgb] *= GMT_u2u[GMT_INCH][S.u];	/* Undo scales since these are angles, not dimensions */

				if (S.read_size) {
					S.size_x = in[three];
					S.size_y = (S.symbol == GMT_SYMBOL_COLUMN) ? in[four] : S.size_x;
				}

				if (n == n_alloc) n_alloc = GMT_alloc_memory ((void **)&data1, n, n_alloc, sizeof (struct PSXYZ_DATA1), GMT_program);
				data1[n].flag = 0;
				data1[n].lon = in[GMT_X];
				data1[n].lat = in[GMT_Y];
				data1[n].string = NULL;
				if (S.read_vector) {
					data1[n].x_size = (float)in[three];	/* direction */
					data1[n].y_size = (float)in[four];	/* length */
					data1[n].z_size = (float)S.v_shrink;	/* Normalization shrink (if given) */
					data1[n].flag = S.convert_angles;	/* Azim vs Dir. */
					if (data1[n].y_size < S.v_norm) data1[n].flag |= 8;	/* Flag to shrink vector attributes */
				}
				else if (S.symbol == GMT_SYMBOL_ELLIPSE) {
					data1[n].x_size = (float)in[three];	/* direction */
					data1[n].y_size = (float)in[four];	/* major */
					data1[n].z_size = (float)in[five];	/* minor */
					data1[n].flag = S.convert_angles;	/* Azim vs Dir. */
					if (S.convert_angles == 1) {	/* Got axes in user units, change to inches */
						data1[n].y_size *= (float)project_info.x_scale;
						data1[n].z_size *= (float)project_info.x_scale;
						data1[n].x_size = (float)(90.0 - data1[n].x_size);	/* Cartesian azimuth */
					}
				}
				else if (S.symbol == GMT_SYMBOL_ROTATERECT) {
					data1[n].x_size = (float)in[three];	/* direction */
					data1[n].y_size = (float)in[four];	/* width */
					data1[n].z_size = (float)in[five];	/* height */
					data1[n].flag = S.convert_angles;	/* Azim vs Dir. */
					if (S.convert_angles == 1) {	/* Got axes in user units, change to inches */
						data1[n].y_size *= (float)project_info.x_scale;
						data1[n].z_size *= (float)project_info.x_scale;
						data1[n].x_size = (float)(90.0 - data1[n].x_size);	/* Cartesian azimuth */
					}
				}
				else if (S.symbol == GMT_SYMBOL_PIE || S.symbol == GMT_SYMBOL_MANGLE) {
					data1[n].x_size = (float)S.size_x;
					if (S.convert_angles == 1) {	/* Got Cartesian azimuths instead */
						data1[n].z_size = (float)(90.0 - in[three+S.read_size]);	/* Start direction in degrees */
						data1[n].y_size = (float)(90.0 - in[four+S.read_size]);		/* Stop direction in degrees */
					}
					else {
						data1[n].y_size = (float)in[three+S.read_size];			/* Start direction in degrees */
						data1[n].z_size = (float)in[four+S.read_size];			/* Stop direction in degrees */
					}
				}
				else if (S.symbol == GMT_SYMBOL_RECT) {
					data1[n].x_size = (float)in[three+S.read_size];	/* x-dim */
					data1[n].y_size = (float)in[four+S.read_size];	/* y-dim */
				}
				else {
					data1[n].x_size = (float)S.size_x;
					data1[n].y_size = (float)S.size_y;
				}
				if (S.user_unit) data1[n].flag |= 4;

				/* Skip zero-size symbols */
		
				if (S.symbol != GMT_SYMBOL_POINT && S.symbol < GMT_SYMBOL_ELLIPSE && S.size_x <= 0.0) {
					if (gmtdefs.verbose) fprintf (stderr, "%s: Warning: Symbol size <= 0.0 near line %ld (skipped)\n", GMT_program, n_total_read);
					continue;
				}
				if (S.read_vector && S.v_just < 3 && data1[n].y_size <= 0.0) {
					if (gmtdefs.verbose) fprintf (stderr, "%s: Warning: vector length <= 0.0 near line %ld (skipped)\n", GMT_program, n_total_read);
					continue;
				}
				if (S.symbol == GMT_SYMBOL_ELLIPSE && (data1[n].y_size <= 0.0 || data1[n].z_size <= 0.0)) {
					if (gmtdefs.verbose) fprintf (stderr, "%s: Warning: ellipse axes <= 0.0 near line %ld (skipped)\n", GMT_program, n_total_read);
					continue;
				}

				if (bar_clip) {
					if (S.symbol == GMT_SYMBOL_BARX) {
						if (in[GMT_X] < west)
							in[GMT_X] = west;
						else if (in[GMT_X] > east)
							in[GMT_X] = east;
					}
					else if (S.symbol == GMT_SYMBOL_BARY) {
						if (in[GMT_Y] < south)
							in[GMT_Y] = south;
						else if (in[GMT_Y] > north)
							in[GMT_Y] = north;
					}
					if (in[GMT_Z] < project_info.z_bottom || in[GMT_Z] > project_info.z_top) continue;
				}
				else if (!Ctrl->N.active) {
					if (in[GMT_Z] < project_info.z_bottom || in[GMT_Z] > project_info.z_top) continue;
					GMT_map_outside (in[GMT_X], in[GMT_Y]);
					if (GMT_abs (GMT_x_status_new) > 1 || GMT_abs (GMT_y_status_new) > 1) continue;
				}

				GMT_project3D (in[GMT_X], in[GMT_Y], in[GMT_Z], &data1[n].x, &data1[n].y, &data1[n].z);
				if (get_rgb) {	/* Lookup z to get rgb */
					data1[n].value = in[3];
					GMT_get_rgb_from_z (data1[n].value, current_fill.rgb);
					if (Ctrl->I.active) GMT_illuminate (Ctrl->I.value, current_fill.rgb);
				}
				if (S.symbol == GMT_SYMBOL_COLUMN) {
					if (GMT_is_dnan (S.base))
						tmp = 0.0;
					else
						GMT_z_to_zz (S.base, &tmp);
					data1[n].z_size = (float)tmp;
				}
				if (S.symbol == GMT_SYMBOL_BARX) {
					if (GMT_is_dnan (S.base))
						tmp = 0.0;
					else
						GMT_x_to_xx (S.base, &tmp);
					data1[n].z_size = (float)tmp;
				}
				if (S.symbol == GMT_SYMBOL_BARY) {
					if (GMT_is_dnan (S.base))
						tmp = 0.0;
					else
						GMT_y_to_yy (S.base, &tmp);
					data1[n].z_size = (float)tmp;
				}
				if (S.symbol == GMT_SYMBOL_TEXT) data1[n].string = strdup (S.string);
				data1[n].symbol = S.symbol;
				if (Ctrl->W.mode) {
					memcpy ((void *)Ctrl->W.pen.rgb, (void *)current_fill.rgb, 3*sizeof (int));
					current_pen = Ctrl->W.pen;
				}
				if (Ctrl->W.mode & 1) memcpy ((void *)current_fill.rgb, (void *)GMT_no_rgb, 3*sizeof (int));
				data1[n].f = current_fill;
				data1[n].p = current_pen;
				data1[n].outline = outline_active;
				n++;
				if (read_symbol) n_expected = GMT_MAX_COLUMNS;
			}
			n_alloc = GMT_alloc_memory ((void **)&data1, 0, n, sizeof (struct PSXYZ_DATA1), GMT_program);
	
			/* Sort according to distance  from viewer */
	
			if (!Ctrl->Q.active) sort_on_distance (data1, n);

			/* Now plot these symbols one at the time */

			for (i = 0; i < n; i++) {

				if (data1[i].symbol == GMT_SYMBOL_COLUMN || data1[i].symbol == GMT_SYMBOL_CUBE) {
					for (j = 0; j < 3; j++) {
						rgb[3*j] = data1[i].f.rgb[0];
						rgb[3*j+1] = data1[i].f.rgb[1];
						rgb[3*j+2] = data1[i].f.rgb[2];
						if (S.shade3D) GMT_illuminate (lux[j], &rgb[3*j]);
					}
				}
		
				if (GMT_cpt_skip) continue;	/* Skip this because cpt says so */

				if (data1[i].outline) GMT_setpen (&data1[i].p);
		
				switch (data1[i].symbol) {
					case GMT_SYMBOL_NONE:
						break;
					case GMT_SYMBOL_STAR:
						size = (S.equal_area) ? 1.67289326141 * data1[i].x_size : data1[i].x_size;
						GMT_star (data1[i].x, data1[i].y, data1[i].z, &size, &data1[i].f, data1[i].outline);
						break;
					case GMT_SYMBOL_BARX:
						if (data1[i].flag & 4) {
							GMT_geo_to_xy (data1[i].x, data1[i].y-0.5*data1[i].x_size, &x_1, &y_1);
							GMT_geo_to_xy (data1[i].x, data1[i].y+0.5*data1[i].x_size, &x_2, &y_2);
							xxx = 0.5 * hypot (x_1 - x_2, y_1 - y_2);
							barx3D (data1[i].x, data1[i].y, data1[i].z, data1[i].z_size, xxx, &data1[i].f, data1[i].outline);
						}
						else
							barx3D (data1[i].x, data1[i].y, data1[i].z, data1[i].z_size, (double)data1[i].x_size, &data1[i].f, data1[i].outline);
						break;
					case GMT_SYMBOL_BARY:
						if (data1[i].flag & 4) {
							GMT_geo_to_xy (data1[i].x-0.5*data1[i].x_size, data1[i].y, &x_1, &y_1);
							GMT_geo_to_xy (data1[i].x+0.5*data1[i].x_size, data1[i].y, &x_2, &y_2);
							xxx = 0.5 * hypot (x_1 - x_2, y_1 - y_2);
							bary3D (data1[i].x, data1[i].y, data1[i].z, data1[i].z_size, xxx, &data1[i].f, data1[i].outline);
						}
						else
							bary3D (data1[i].x, data1[i].y, data1[i].z, data1[i].z_size, (double)data1[i].x_size, &data1[i].f, data1[i].outline);
						break;
					case GMT_SYMBOL_COLUMN:
						if (data1[i].flag & 4) {
							GMT_geo_to_xy (data1[i].x-data1[i].x_size, data1[i].y, &x_1, &y_1);
							GMT_geo_to_xy (data1[i].x+data1[i].x_size, data1[i].y, &x_2, &y_2);
							xxx = 0.5 * hypot (x_1 - x_2, y_1 - y_2);
							GMT_geo_to_xy (data1[i].x, data1[i].y-data1[i].y_size, &x_1, &y_1);
							GMT_geo_to_xy (data1[i].x, data1[i].y+data1[i].y_size, &x_2, &y_2);
							yyy = 0.5 * hypot (x_1 - x_2, y_1 - y_2);
							column3D (data1[i].x, data1[i].y, data1[i].z, data1[i].z_size, xxx, yyy, rgb, data1[i].outline);
						}
						else
							column3D (data1[i].x, data1[i].y, data1[i].z, data1[i].z_size, (double)data1[i].x_size, (double)data1[i].y_size, rgb, data1[i].outline);
						break;
					case GMT_SYMBOL_CROSS:
						dim[0] = (double)data1[i].x_size;
						GMT_cross (data1[i].x, data1[i].y, data1[i].z, dim, NULL, FALSE);
						break;
					case GMT_SYMBOL_PLUS:
						dim[0] = (double)data1[i].x_size;
						GMT_plus (data1[i].x, data1[i].y, data1[i].z, dim, NULL, FALSE);
						break;
					case GMT_SYMBOL_POINT:
						dim[0] = GMT_POINT_SIZE;
						GMT_plus (data1[i].x, data1[i].y, data1[i].z, dim, NULL, FALSE);
						break;
					case GMT_SYMBOL_CIRCLE:
						size = (double)data1[i].x_size;
						GMT_circle (data1[i].x, data1[i].y, data1[i].z, &size, &data1[i].f, data1[i].outline);
						break;
					case GMT_SYMBOL_SQUARE:
						size = (S.equal_area) ? 1.25331413732 * data1[i].x_size : data1[i].x_size;
						GMT_square (data1[i].x, data1[i].y, data1[i].z, &size, &data1[i].f, data1[i].outline);
						break;
					case GMT_SYMBOL_HEXAGON:
						size = (S.equal_area) ? 1.09963611079 * data1[i].x_size : data1[i].x_size;
						GMT_hexagon (data1[i].x, data1[i].y, data1[i].z, &size, &data1[i].f, data1[i].outline);
						break;
					case GMT_SYMBOL_MANGLE:
						dim[0] = 0.5 * data1[i].x_size; dim[1] = data1[i].y_size; dim[2] = data1[i].z_size;
						GMT_matharc (data1[i].x, data1[i].y, data1[i].z, dim, gmtdefs.vector_shape,  &data1[i].p, S.v_double_heads);
						break;
					case GMT_SYMBOL_PENTAGON:
						size = (S.equal_area) ? 1.14948092619 * data1[i].x_size : data1[i].x_size;
						GMT_pentagon (data1[i].x, data1[i].y, data1[i].z, &size, &data1[i].f, data1[i].outline);
						break;
					case GMT_SYMBOL_OCTAGON:
						size = (S.equal_area) ? 1.05390736526 * data1[i].x_size : data1[i].x_size;
						GMT_octagon (data1[i].x, data1[i].y, data1[i].z, &size, &data1[i].f, data1[i].outline);
						break;
					case GMT_SYMBOL_TRIANGLE:
						size = (S.equal_area) ? 1.55512030156 * data1[i].x_size : data1[i].x_size;
						GMT_triangle (data1[i].x, data1[i].y, data1[i].z, &size, &data1[i].f, data1[i].outline);
						break;
					case GMT_SYMBOL_ITRIANGLE:
						size = (S.equal_area) ? 1.55512030156 * data1[i].x_size : data1[i].x_size;
						GMT_itriangle (data1[i].x, data1[i].y, data1[i].z, &size, &data1[i].f, data1[i].outline);
						break;
					case GMT_SYMBOL_DIAMOND:
						size = (S.equal_area) ? 1.25331413732 * data1[i].x_size : data1[i].x_size;
						GMT_diamond (data1[i].x, data1[i].y, data1[i].z, &size, &data1[i].f, data1[i].outline);
						break;
					case GMT_SYMBOL_CUBE:
						cube3D (data1[i].x, data1[i].y, data1[i].z, (double)data1[i].x_size, (double)data1[i].y_size, rgb, data1[i].outline);
						break;
					case GMT_SYMBOL_ELLIPSE:
						if (data1[i].flag == 2)
							GMT_plot_ellipse (data1[i].lon, data1[i].lat, data1[i].z, data1[i].y_size, data1[i].z_size, data1[i].x_size, current_fill, data1[i].outline);
						else {
							GMT_flip_angle_f (&data1[i].x_size);
							dim[0] = data1[i].x_size; dim[1] = data1[i].y_size; dim[2] = data1[i].z_size;
							GMT_ellipse (data1[i].x, data1[i].y, data1[i].z, dim, &data1[i].f, data1[i].outline);
						}
						break;
					case GMT_SYMBOL_ROTATERECT:
						if (data1[i].flag == 2)
							GMT_plot_rectangle (data1[i].lon, data1[i].lat, data1[i].z, data1[i].y_size, data1[i].z_size, data1[i].x_size, current_fill, data1[i].outline);
						else {
							GMT_flip_angle_f (&data1[i].x_size);
							dim[0] = data1[i].x_size; dim[1] = data1[i].y_size; dim[2] = data1[i].z_size;
							GMT_rotrect (data1[i].x, data1[i].y, data1[i].z, dim, &data1[i].f, data1[i].outline);
						}
						break;
					case GMT_SYMBOL_TEXT:
						font_size = data1[i].x_size * 72.0;	/* To get points */
						if (data1[i].outline && fill_active) {
							ps_setpaint (data1[i].f.rgb);
							GMT_text3D (data1[i].x, data1[i].y, data1[i].z, font_size, S.font_no, data1[i].string, 0.0, 6, FALSE);
							ps_setpaint (data1[i].p.rgb);
							GMT_text3D (data1[i].x, data1[i].y, data1[i].z, font_size, S.font_no, data1[i].string, 0.0, 6, TRUE);
						}
						else if (fill_active)
							GMT_text3D (data1[i].x, data1[i].y, data1[i].z, font_size, S.font_no, data1[i].string, 0.0, 6, FALSE);
						else
							GMT_text3D (data1[i].x, data1[i].y, data1[i].z, font_size, S.font_no, data1[i].string, 0.0, 6, TRUE);
						free ((void *)data1[i].string);
						break;
					case GMT_SYMBOL_VECTOR:
						if (data1[i].flag & 2) {
							GMT_azim_to_angle (data1[i].lon, data1[i].lat, 0.1, data1[i].x_size, &tmp);
							data1[i].x_size = (float)tmp;
						}
						else if (data1[i].flag & 1) {
							data1[i].x_size = (float)(90.0 - data1[i].x_size);
						}
						if (S.v_just == 3) {
							GMT_geo_to_xy (in[pos2x], in[pos2y], &x2, &y2);
							if (GMT_is_dnan (x2) || GMT_is_dnan (y2)) {
								fprintf (stderr, "%s: Warning: Vector head coordinates contain NaNs near line %ld. Skipped\n", GMT_program, n_total_read);
								continue;
							}
						}
						else {
							GMT_flip_angle_f (&data1[i].x_size);
							sincosd (data1[i].x_size, &s, &c);
							x2 = data1[i].x + data1[i].y_size * c;
							y2 = data1[i].y + data1[i].y_size * s;
							if (S.v_just) {
								dx = S.v_just * 0.5 * (x2 - data1[i].x);	dy = S.v_just * 0.5 * (y2 - data1[i].y);
								data1[i].x -= dx;	data1[i].y -= dy;
								x2 -= dx;		y2 -= dy;
							}
						}
						this_outline = (S.v_double_heads) ? data1[i].outline + 8 : data1[i].outline;
						GMT_vector (data1[i].x, data1[i].y, x2, y2, data1[i].z, S.v_width, S.h_length, S.h_width, gmtdefs.vector_shape, &data1[i].f, (GMT_LONG)this_outline);
						break;
					case GMT_SYMBOL_VECTOR2:
						if (data1[i].flag & 2) {
							GMT_azim_to_angle (data1[i].lon, data1[i].lat, 0.1, data1[i].x_size, &tmp);
							data1[i].x_size = (float)tmp;
						}
						else if (data1[i].flag & 1) {
							data1[i].x_size = (float)(90.0 - data1[i].x_size);
						}
						if (S.v_just == 3) {
							GMT_geo_to_xy (in[pos2x], in[pos2y], &x2, &y2);
							if (GMT_is_dnan (x2) || GMT_is_dnan (y2)) {
								fprintf (stderr, "%s: Warning: Vector head coordinates contain NaNs near line %ld. Skipped\n", GMT_program, n_total_read);
								continue;
							}
						}
						else {
							GMT_flip_angle_f (&data1[i].x_size);
							sincosd (data1[i].x_size, &s, &c);
							x2 = data1[i].x + data1[i].y_size * c;
							y2 = data1[i].y + data1[i].y_size * s;
							if (S.v_just) {
								dx = S.v_just * 0.5 * (x2 - data1[i].x);	dy = S.v_just * 0.5 * (y2 - data1[i].y);
								data1[i].x -= dx;	data1[i].y -= dy;
								x2 -= dx;		y2 -= dy;
							}
						}
						if (data1[i].flag & 8) {	/* Scale arrow attributes down with length */
							v_w = S.v_width * data1[i].y_size * data1[i].z_size;
							h_l = S.h_length * data1[i].y_size * data1[i].z_size;
							h_w = S.h_width * data1[i].y_size * data1[i].z_size;
							this_outline = (S.v_double_heads) ? data1[i].outline + 8 : data1[i].outline;
							GMT_vector (data1[i].x, data1[i].y, x2, y2, data1[i].z, v_w, h_l, h_w, gmtdefs.vector_shape, &data1[i].f, (GMT_LONG)this_outline);
						}
						else {	/* Leave as specified */
							this_outline = (S.v_double_heads) ? data1[i].outline + 8 : data1[i].outline;
							GMT_vector (data1[i].x, data1[i].y, x2, y2, data1[i].z, S.v_width, S.h_length, S.h_width, gmtdefs.vector_shape, &data1[i].f, (GMT_LONG)this_outline);
						}
						break;
					case GMT_SYMBOL_PIE:
						if (S.convert_angles == 2) {
							GMT_azim_to_angle (data1[i].lon, data1[i].lat, 0.1, (double)(data1[i].y_size), &tmp);
							data1[i].y_size = (float)(tmp);
							GMT_azim_to_angle (data1[i].lon, data1[i].lat, 0.1, (double)(data1[i].z_size), &tmp);
							data1[i].z_size = (float)tmp;
							f_swap (data1[i].y_size, data1[i].z_size);
						}
						dim[0] = 0.5 * data1[i].x_size; dim[1] = data1[i].y_size; dim[2] = data1[i].z_size;
						GMT_pie (data1[i].x, data1[i].y, data1[i].z, dim, &data1[i].f, data1[i].outline);
						break;
					case GMT_SYMBOL_XDASH:
						GMT_xyz_to_xy (data1[i].x-0.5*(double)data1[i].x_size, data1[i].y, data1[i].z, &xx[0], &yy[0]);
						GMT_xyz_to_xy (data1[i].x+0.5*(double)data1[i].x_size, data1[i].y, data1[i].z, &xx[1], &yy[1]);
						ps_segment (xx[0], yy[0], xx[1], yy[1]);
						break;
					case GMT_SYMBOL_YDASH:
						GMT_xyz_to_xy (data1[i].x, data1[i].y-0.5*(double)data1[i].x_size, data1[i].z, &xx[0], &yy[0]);
						GMT_xyz_to_xy (data1[i].x, data1[i].y+0.5*(double)data1[i].x_size, data1[i].z, &xx[1], &yy[1]);
						ps_segment (xx[0], yy[0], xx[1], yy[1]);
						break;
					case GMT_SYMBOL_ZDASH:
						GMT_xyz_to_xy (data1[i].x, data1[i].y, data1[i].z-0.5*(double)data1[i].x_size, &xx[0], &yy[0]);
						GMT_xyz_to_xy (data1[i].x, data1[i].y, data1[i].z+0.5*(double)data1[i].x_size, &xx[1], &yy[1]);
						ps_segment (xx[0], yy[0], xx[1], yy[1]);
						break;
					case GMT_SYMBOL_RECT:
						dim[0] = data1[i].x_size; dim[1] = data1[i].y_size;
						GMT_rect (data1[i].x-0.5*dim[0], data1[i].y-0.5*dim[1], data1[i].z, dim, &data1[i].f, data1[i].outline);
						break;
					case GMT_SYMBOL_CUSTOM:
						GMT_draw_custom_symbol (data1[i].x, data1[i].y, data1[i].z, (double)data1[i].x_size, S.custom, &current_pen, &current_fill, data1[i].outline);
						break;
				}
			}
		}
		else {	/* Line/polygon part */
			n_required = 3;
			n_fields = GMT_input (fp, &n_expected, &in);
			while (!GMT_REC_IS_EOF) {
				while (GMT_REC_IS_NEW_SEGMENT) {	/* Process multisegment header flags (or NaN gaps) */
					if (GMT_REC_IS_SEG_HEADER) {
						if (gmtdefs.verbose) ps_comment (GMT_io.segment_header);
					
						change = GMT_parse_multisegment_header (GMT_io.segment_header, Ctrl->C.active, &fill_active, &current_fill, &default_fill, &outline_active, &current_pen, &default_pen, default_outline);
						if (Ctrl->I.active) {
							GMT_illuminate (Ctrl->I.value, current_fill.rgb);
							GMT_illuminate (Ctrl->I.value, default_fill.rgb);
						}
						if ((change & 4) && penset_OK) GMT_setpen (&current_pen);
						if ((change & 1)) polygon = TRUE;
						if ((change & 2) && !Ctrl->L.active) {
							polygon = FALSE;
							ps_setpaint (current_fill.rgb);
						}
						if (S.G.label_type == 2) {	/* Update segment header */
							GMT_extract_label (&GMT_io.segment_header[1], S.G.label);
						}
					}
					n_fields = GMT_input (fp, &n_expected, &in);
					n_total_read++;
				}
				if (GMT_REC_IS_EOF) continue;		/* At EOF */
				if (GMT_REC_IS_GAP) GMT_io.status = 0;	/* Done with gap */

				n = 0;
				while (!GMT_REC_IS_LINE_BREAK) {	/* Keep going until FALSE or = 2 segment header or gap */
					if (GMT_REC_IS_ERROR) {
						fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld\n", GMT_program, n_fields, n_expected, n_total_read);
						exit (EXIT_FAILURE);
					}

					if (n == n_alloc) n_alloc = GMT_alloc_memory ((void **)&data2, n, n_alloc, sizeof (struct PSXYZ_DATA2), GMT_program);
					data2[n].x = in[GMT_X];	data2[n].y = in[GMT_Y];
					data2[n].z = in[GMT_Z];
					n++;
					n_fields = GMT_input (fp, &n_expected, &in);
				}
				
				/* Done reading a segment */
				
				if (polygon) {  /* Explicitly close polygon */
					if (n == n_alloc) n_alloc = GMT_alloc_memory ((void **)&data2, n, n_alloc, sizeof (struct PSXYZ_DATA2), GMT_program);
					data2[n].x = data2[0].x;
					data2[n].y = data2[0].y;
					data2[n].z = data2[0].z;
					n++;
				}
				n_alloc = GMT_alloc_memory ((void **)&data2, n, n_alloc, sizeof (struct PSXYZ_DATA2), GMT_program);
		
				if (!GMT_cpt_skip) {
					double *xp, *yp;
					xp = (double *) GMT_memory (VNULL, (size_t)n, sizeof (double), GMT_program);
					yp = (double *) GMT_memory (VNULL, (size_t)n, sizeof (double), GMT_program);
		
					if (polygon) {
						for (i = 0; i < n; i++) GMT_geoz_to_xy (data2[i].x, data2[i].y, data2[i].z, &xp[i], &yp[i]);
						GMT_fill (xp, yp, n, &current_fill, outline_active);
					}
					else if (S.symbol == GMT_SYMBOL_QUOTED_LINE) {
						for (i = 0; i < n; i++) {
							GMT_geo_to_xy (data2[i].x, data2[i].y, &xp[i], &yp[i]);
						}
						S.G.z_level = data2[0].z;	/* Any z will do since we require line to be in x-y plane */
						S.G.line_pen = current_pen;
						GMT_hold_contour (&xp, &yp, n, 0.0, "N/A", 'A', S.G.label_angle, Ctrl->L.active, &S.G);
					}
					else {
						for (i = 0; i < n; i++) GMT_geoz_to_xy (data2[i].x, data2[i].y, data2[i].z, &xp[i], &yp[i]);
						ps_line (xp, yp, n, 3, FALSE);
					}
					if (S.symbol == GMT_SYMBOL_FRONT) {
						for (i = 0; i < n; i++) GMT_geo_to_xy (data2[i].x, data2[i].y, &xp[i], &yp[i]);
						GMT_draw_fence (xp, yp, project_info.z_level, n, &S.f, &current_fill, outline_active); 	/* Must draw fault crossbars */
					}

					GMT_free ((void *)xp);
					GMT_free ((void *)yp);
				}
			}
		}
		if (fp != GMT_stdin) GMT_fclose (fp);
	}
	if (data1) GMT_free ((void *)data1);
	if (data2) GMT_free ((void *)data2);

	if (S.symbol == GMT_SYMBOL_QUOTED_LINE) {
		GMT_contlabel_plot (&S.G);
		GMT_contlabel_free (&S.G);
	}

	if (S.u_set) gmtdefs.measure_unit = save_u;	/* Reset unit */

	if (Ctrl->D.active) ps_transrotate (-DX, -DY, 0.0);	/* Shift plot a bit */

	if (clip_set) GMT_map_clip_off ();

	if (current_pen.texture) ps_setdash (CNULL, 0);
	if (project_info.three_D) GMT_vertical_axis (2);	/* Draw background axis */

	ps_rotatetrans (z_project.xmin, z_project.ymin, 0.0);
	GMT_plotend ();

	Free_psxyz_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void column3D (double x, double y, double z, double base, double x_size, double y_size, int *rgb, GMT_LONG outline)
{
	GMT_LONG i, j, k;
	double zz, xp[4], yp[4], zp[4], plot_x[4], plot_y[4], top, sign;

	x_size *= 0.5;
	y_size *= 0.5;
	top = z;
	if (top < base) d_swap (top, base);

	for (i = 0; i < 3; i++) {
		sign = -1.0;
		zz = base;
		switch (z_project.face[i]) {
			case 0:	/* yz plane positive side */
				sign = 1.0;
			case 1:	/* negative side */
				xp[0] = xp[1] = xp[2] = xp[3] = x + sign * x_size;
				yp[0] = yp[3] = y - y_size;	yp[1] = yp[2] = y + y_size;
				zp[0] = zp[1] = base;	zp[2] = zp[3] = top;
				break;
			case 2:	/* xz plane positive side */
				sign = 1.0;
			case 3:	/* negative side */
				xp[0] = xp[3] = x - x_size;	xp[1] = xp[2] = x + x_size;
				yp[0] = yp[1] = yp[2] = yp[3] = y + sign * y_size;
				zp[0] = zp[1] = base;	zp[2] = zp[3] = top;
				break;
			case 4:	/* xy plane positive side */
				zz = top;
			case 5:	/* negative side */
				xp[0] = xp[3] = x - x_size;	yp[0] = yp[1] = y - y_size;
				xp[1] = xp[2] = x + x_size;	yp[2] = yp[3] = y + y_size;
				zp[0] = zp[1] = zp[2] = zp[3] = zz;
				break;
		}
		k = z_project.face[i] / 2;
		for (j = 0; j < 4; j++) GMT_xyz_to_xy (xp[j], yp[j], zp[j], &plot_x[j], &plot_y[j]);
		ps_patch (plot_x, plot_y, (GMT_LONG)4, &rgb[3*k], outline);
	}
}

void cube3D (double x, double y, double z, double x_size, double y_size, int *rgb, GMT_LONG outline)
{
	GMT_LONG i, j, k;
	double xp[4], yp[4], zp[4], plot_x[4], plot_y[4], sign;

	x_size *= 0.5;
	y_size *= 0.5;
	for (i = 0; i < 3; i++) {
		sign = -1.0;
		switch (z_project.face[i]) {
			case 4:	/* xy plane positive side */
				sign = 1.0;
			case 5:	/* negative side */
				xp[0] = xp[3] = x - x_size;	yp[0] = yp[1] = y - y_size;
				xp[1] = xp[2] = x + x_size;	yp[2] = yp[3] = y + y_size;
				zp[0] = zp[1] = zp[2] = zp[3] = z + sign * x_size;
				break;
			case 2:	/* xz plane positive side */
				sign = 1.0;
			case 3:	/* negative side */
				xp[0] = xp[3] = x - x_size;	xp[1] = xp[2] = x + x_size;
				yp[0] = yp[1] = yp[2] = yp[3] = y + sign * y_size;
				zp[0] = zp[1] = z - x_size;	zp[2] = zp[3] = z + x_size;
				break;
			case 0:	/* yz plane positive side */
				sign = 1.0;
			case 1:	/* negative side */
				xp[0] = xp[1] = xp[2] = xp[3] = x + sign * x_size;
				yp[0] = yp[3] = y - y_size;	yp[1] = yp[2] = y + y_size;
				zp[0] = zp[1] = z - x_size;	zp[2] = zp[3] = z + x_size;
				break;
		}
		k = z_project.face[i] / 2;
		for (j = 0; j < 4; j++) GMT_xyz_to_xy (xp[j], yp[j], zp[j], &plot_x[j], &plot_y[j]);
		ps_patch (plot_x, plot_y, (GMT_LONG)4, &rgb[3*k], outline);
	}
}

void bary3D (double x, double y, double z, double base, double size, struct GMT_FILL *f, GMT_LONG outline)
{
	GMT_LONG i;
	double xp[4], yp[4], plot_x[4], plot_y[4];

	size *= 0.5;
	xp[0] = xp[3] = x - size;	xp[1] = xp[2] = x + size;
	yp[0] = yp[1] = base;	yp[2] = yp[3] = y;
	for (i = 0; i < 4; i++) GMT_xyz_to_xy (xp[i], yp[i], z, &plot_x[i], &plot_y[i]);
	GMT_fill (plot_x, plot_y, (GMT_LONG)4, f, (GMT_LONG)outline);
}

void barx3D (double x, double y, double z, double base, double size, struct GMT_FILL *f, GMT_LONG outline)
{
	GMT_LONG i;
	double xp[4], yp[4], plot_x[4], plot_y[4];

	size *= 0.5;
	yp[0] = yp[3] = y - size;	yp[1] = yp[2] = y + size;
	xp[0] = xp[1] = base;	xp[2] = xp[3] = x;
	for (i = 0; i < 4; i++) GMT_xyz_to_xy (xp[i], yp[i], z, &plot_x[i], &plot_y[i]);
	GMT_fill (plot_x, plot_y, (GMT_LONG)4, f, (GMT_LONG)outline);
}

void sort_on_distance (struct PSXYZ_DATA1 *data, GMT_LONG n)
{
	/* This function sorts the data array such that points farthest away are plotted first */
	GMT_LONG i;

	int dist_compare(const void *a, const void *b);

	for (i = 0; i < n; i++) data[i].dist = - z_project.sin_az * data[i].x - z_project.cos_az * data[i].y;

	qsort ((void *)data, (size_t)n, sizeof (struct PSXYZ_DATA1), dist_compare);
}

int dist_compare (const void *a, const void *b)
{
	double z;

	if ( ((struct PSXYZ_DATA1 *)a)->dist < ((struct PSXYZ_DATA1 *)b)->dist) {
		return (-1);
	}
	else if ( ((struct PSXYZ_DATA1 *)a)->dist > ((struct PSXYZ_DATA1 *)b)->dist) {
		return (1);
	}
	else {
		z = ( ((struct PSXYZ_DATA1 *)a)->z - ((struct PSXYZ_DATA1 *)b)->z) * z_project.sin_el;
		if (z < 0.0)
			return (-1);
		else if (z > 0.0)
			return (1);
		else
			return (0);
	}
}

void *New_psxyz_Ctrl () {	/* Allocate and initialize a new control structure */
	struct PSXYZ_CTRL *C;
	
	C = (struct PSXYZ_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct PSXYZ_CTRL), "New_psxyz_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
		
	C->E.azimuth = 180.0;
	C->E.elevation = 90.0;
	GMT_init_pen (&C->W.pen, GMT_PENWIDTH);
	GMT_init_fill (&C->G.fill, -1, -1, -1);	/* Default is no fill */
	C->A.step = gmtdefs.line_step;
	return ((void *)C);
}

void Free_psxyz_Ctrl (struct PSXYZ_CTRL *C) {	/* Deallocate control structure */
	if (C->C.file) free ((void *)C->C.file);	
	GMT_free ((void *)C);	
}
