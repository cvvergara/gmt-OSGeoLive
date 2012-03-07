/*--------------------------------------------------------------------
 *	$Id: psxy.c,v 1.203 2011/07/08 21:27:06 guru Exp $
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
 * psxy will read <x,y> pairs (or <lon,lat>) from GMT_stdin and plot
 * the lines, polygons, or symbols on a map. A variety of symbols
 * may be specified, some of which require additional columns of
 * data, like symbol-size etc.  Only one symbol may be plotted
 * at the time.  PostScript code is written to GMT_stdout.
 *
 * Author:	Paul Wessel
 * Date:	15-SEP-2001
 * Version:	4
 *
 */

#include "gmt.h"
#include "pslib.h"

struct PSXY_CTRL {
	struct A {	/* -A[m|p|step] */
		GMT_LONG active;
		GMT_LONG mode;
		double step;
	} A;
	struct C {	/* -C<cpt> */
		GMT_LONG active;
		char *file;
	} C;
	struct D {	/* -D<dx>/<dy> */
		GMT_LONG active;
		double dx, dy;
	} D;
	struct E {	/* -E[x|X][y|Y][cap][/[+|-]<pen>] */
		GMT_LONG active;
		GMT_LONG xbar, ybar;	/* 0 = not used, 1 = error bar, 2 = box-whisker, 3 notched box-whisker */
		GMT_LONG mode;	/* 0 = normal, 1 = -C applies to error pen color, 2 = -C applies to symbol fill & error pen color */
		double size;
		struct GMT_PEN pen;
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
	struct S {	/* -S */
		GMT_LONG active;
		char *arg;
	} S;
	struct T {	/* -T */
		GMT_LONG active;
	} T;
	struct W {	/* -W<pen> */
		GMT_LONG active;
		GMT_LONG mode;	/* 0 = normal, 1 = -C applies to pen color only, 2 = -C applies to symbol fill & pen color */
		struct GMT_PEN pen;
	} W;
};

int main (int argc, char **argv)
{
	GMT_LONG error = FALSE, nofile, polygon = FALSE, read_symbol = FALSE, penset_OK = TRUE, outline_active;
	GMT_LONG default_outline = FALSE, done, get_rgb = FALSE, resample, bar_clip, save_u = FALSE;
	GMT_LONG old_GMT_world_map, error_x = FALSE, error_y = FALSE, def_err_xy = FALSE, clip_set = FALSE, fill_active;

	char line[BUFSIZ], txt_a[GMT_LONG_TEXT], txt_b[GMT_LONG_TEXT], *not_used = NULL;

	GMT_LONG i, n, n_alloc = 0, n_total_read = 0, plot_n, pos2x, pos2y;
	GMT_LONG j, j0, n_files = 0, fno, xy_errors[2], two, three, four, error_type[2] = {0,0};
	GMT_LONG n_args, n_fields, n_cols_start = 2, this_outline, change;
	GMT_LONG n_required, n_expected = 0, reset_cap = 0, error_cols[3] = {1,4,5};

	double west = 0.0, east = 0.0, south = 0.0, north = 0.0, s, c, step;
	double plot_x, plot_y, x_1, x_2, y_1, y_2, x0, y0;
	double direction = 0.0, length = 0.0, major = 0.0, minor = 0.0, x2, y2;
	double dir1 = 0.0, dir2 = 0.0, x_len = 0.0, y_len = 0.0, dx, dy, size;
	double v_w, h_l, h_w, tmp, *in = NULL, font_size, dim[3];
	double *xx = NULL, *yy = NULL, *xp = NULL, *yp = NULL;

	FILE *fp = NULL;

	struct GMT_PEN current_pen, default_pen;
	struct GMT_FILL current_fill, default_fill;
	struct GMT_SYMBOL S;
	struct PSXY_CTRL *Ctrl = NULL;

	void plot_x_errorbar (double x, double y, double delta_x, double error_width2, GMT_LONG line);
	void plot_y_errorbar (double x, double y, double delta_y, double error_width2, GMT_LONG line);
	void plot_x_whiskerbar (double x, double y, double hinge[], double error_width2, int rgb[], GMT_LONG line, GMT_LONG notched);
	void plot_y_whiskerbar (double x, double y, double hinge[], double error_width2, int rgb[], GMT_LONG line, GMT_LONG notched);
	void *New_psxy_Ctrl (), Free_psxy_Ctrl (struct PSXY_CTRL *C);
	
	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct PSXY_CTRL *)New_psxy_Ctrl ();	/* Allocate and initialize a new control structure */

	xy_errors[GMT_X] = xy_errors[1] = 0;	/* These will be col # of where to find this info in data */

	/* Initialize GMT_SYMBOL structure */

	memset ((void *)&S, 0, sizeof (struct GMT_SYMBOL));
	GMT_contlabel_init (&S.G, 0);
	
	S.v_width = 0.03;
	S.h_length = 0.12;
	S.h_width = 0.1;
	S.font_no = gmtdefs.annot_font[0];

	if ((S.u = gmtdefs.measure_unit) == GMT_CM) {
		S.v_width = 0.075 / 2.54;	S.h_length = 0.3 / 2.54;
		S.h_width = 0.25 / 2.54; Ctrl->E.size = 0.25 / 2.54;
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
	
				case 'A':	/* Turn off draw_arc mode */
					Ctrl->A.active = TRUE;
					if (argv[i][2] == 'm')
						Ctrl->A.mode = 1;
					else if (argv[i][2] == 'p')
						Ctrl->A.mode = 2;
					else if (argv[i][2])
						Ctrl->A.step = atof (&argv[i][2]);	/* Undocumented test feature */
					break;
			
				case 'C':	/* Vary symbol color with z */
					Ctrl->C.file = strdup (&argv[i][2]);
					Ctrl->C.active = TRUE;
					break;
				case 'D':
					if ((j = sscanf (&argv[i][2], "%[^/]/%s", txt_a, txt_b)) < 1) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -D option: Give x [and y] offsets\n", GMT_program);
						error++;
					}
					else {
						Ctrl->D.dx = GMT_convert_units (txt_a, GMT_INCH);
						Ctrl->D.dy = (j == 2) ? GMT_convert_units (txt_b, GMT_INCH) : Ctrl->D.dx;
						Ctrl->D.active = TRUE;
					}
					break;
				case 'E':	/* Get info for error bars and bow-whisker bars */
					Ctrl->E.active = TRUE;
					j = 2;
					while (argv[i][j] && argv[i][j] != '/') {
						if (argv[i][j] == 'x')		/* Error bar for x */
							Ctrl->E.xbar = 1;
						else if (argv[i][j] == 'X') {	/* Box-whisker instead */
							Ctrl->E.xbar = 2;
							if (argv[i][j+1] == 'n') {Ctrl->E.xbar++; j++;}
						}
						else if (argv[i][j] == 'y')	/* Error bar for y */
							Ctrl->E.ybar = 1;
						else if (argv[i][j] == 'Y') {	/* Box-whisker instead */
							Ctrl->E.ybar = 2;
							if (argv[i][j+1] == 'n') {Ctrl->E.ybar++; j++;}
						}
						else {	/* Get error 'cap' width */
							strcpy (txt_a, &argv[i][j]);
							j0 = 0;
							while (txt_a[j0] && txt_a[j0] != '/') j0++;
							txt_a[j0] = 0;
							Ctrl->E.size = GMT_convert_units (txt_a, GMT_INCH);
							while (argv[i][j] && argv[i][j] != '/') j++;
							j--;
						}
						j++;
					}
					if (argv[i][j] == '/') {
						j++;
						if (argv[i][j] == '-') {Ctrl->E.mode = 1; j++;}
						if (argv[i][j] == '+') {Ctrl->E.mode = 2; j++;}
						if (argv[i][j]) GMT_getpen (&argv[i][j], &Ctrl->E.pen);
					}
					break;
				case 'G':		/* Set fill for symbols or polygon */
					Ctrl->G.active = TRUE;
					if (!argv[i][2] || GMT_getfill (&argv[i][2], &Ctrl->G.fill)) {
						GMT_fill_syntax ('G', " ");
						error++;
					}
					break;
				case 'I':	/* Adjust symbol color via intensity */
					Ctrl->I.value = atof (&argv[i][2]);
					Ctrl->I.active = TRUE;
					break;
				case 'L':		/* Close line segments */
					Ctrl->L.active = TRUE;
					break;
				case 'N':		/* Do not skip points outside border */
					Ctrl->N.active = TRUE;
					break;
				case 'S':		/* Get symbol [and size] */
					Ctrl->S.active = TRUE;
					Ctrl->S.arg = &argv[i][2];
					break;
				case 'T':		/* Skip all input files */
					Ctrl->T.active = TRUE;
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
			
				/* Illegal options */
	
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
		fprintf (stderr,"psxy %s - Plot lines, polygons, and symbols on maps\n\n", GMT_VERSION);
		fprintf (stderr,"usage: psxy <infiles> %s %s [-A[m|p]] [%s] [-C<cpt>] [-D<dx>/<dy>]\n", GMT_J_OPT, GMT_Rgeo_OPT, GMT_B_OPT);
		fprintf (stderr, "\t[-E[x|y|X|Y][n][cap][/[+|-]<pen>]] [-G<fill>] [%s] [-I<intens>] [-K] [-L] [-N] [-O] [-P]\n", GMT_H_OPT);
		fprintf (stderr, "\t[-S[<symbol>][<size>]] [-T] [%s] [-V] [-W[+|-][<pen>]] [%s] [%s]\n", GMT_U_OPT, GMT_X_OPT, GMT_Y_OPT);
		fprintf (stderr, "\t[%s] [%s] [%s] [%s] [%s] [%s]\n\n", GMT_c_OPT, GMT_t_OPT, GMT_bi_OPT, GMT_f_OPT, GMT_g_OPT, GMT_mo_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\t<infiles> is one or more files.  If no, read standard input\n");
		GMT_explain_option ('j');
		GMT_explain_option ('R');
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-A Suppress drawing line segments as great circle arcs, i.e. draw straight lines unless\n");
		fprintf (stderr, "\t   m or p is appended to first follow meridian then parallel, or vice versa.\n");
		GMT_explain_option ('b');
		fprintf (stderr, "\t-C Use cpt-file to assign symbol colors based on z-value in 3rd column (requires -S)\n");
		fprintf (stderr, "\t   Without -S, psxy excepts -m lines/polygons and looks for -Z<val> options in each multiheader\n");
		fprintf (stderr, "\t   Then, color is applied to polygon fill (-L) or polygon pen\n");
		fprintf (stderr, "\t-D Offset symbol or line positions by <dx>/<dy> [no offset].\n");
		fprintf (stderr, "\t-E means draw error bars for x, y, or both.  Add cap-width [%g]. Append pen attributes;\n", Ctrl->E.size);
		fprintf (stderr, "\t   leading + applies cpt color (-C) to both symbol fill and error pen, - applies it to pen only.\n");
		fprintf (stderr, "\t   If X or Y is used instead a box-and-whisker diagram is drawn instead, using data from 4\n");
		fprintf (stderr, "\t   extra columns to get the 0 %%, 25 %%, 75 %%, and 100%% quantiles (point value is assumed to be 50%%)\n");
		fprintf (stderr, "\t   If n is appended after X (or Y) we expect a 5th extra column with the sample size which is needed to draw\n");
		fprintf (stderr, "\t   a notched box-and whisker diagram (notch width represents uncertainty in median).\n");
		fprintf (stderr, "\t   Finally, use -W, -G to control the appearance of the 25-75%% box.\n");
		GMT_fill_syntax ('G', "Specify color or pattern [Default is no fill].");
		fprintf (stderr, "\t   -G option can be present in all subheaders (requires -m and no -S)\n");
		GMT_explain_option ('H');
		fprintf (stderr, "\t-I sets a constant intensity used to modulate the fill color (requires -C or -G)\n");
		GMT_explain_option ('K');
		fprintf (stderr, "\t-L Force closed polygons.\n");
		fprintf (stderr, "\t-N Do Not skip/clip symbols that fall outside map border [Default will ignore those outside]\n");
		GMT_explain_option ('O');
		GMT_explain_option ('P');
		fprintf (stderr, "\t-S to select symbol type and symbol size (in %s).  Choose between\n", GMT_unit_names[gmtdefs.measure_unit]);
		fprintf (stderr, "\t   -(xdash), +(plus), st(a)r, (b|B)ar, (c)ircle, (d)iamond, (e)llipse, (f)ront, octa(g)on, (h)exagon, (i)nvtriangle,\n");
		fprintf (stderr, "\t   (j)rotated rectangle, (k)ustom, (l)etter, (m)athangle, pe(n)tagon, (p)oint, (q)uoted line, (r)ect, (s)quare,\n");
		fprintf (stderr, "\t   (t)riangle, (v)ector, (w)edge, (x)cross, (y)dash\n");
		fprintf (stderr, "\t   If no size is specified, psxy expects the 3rd column to have sizes.\n");
		fprintf (stderr, "\t   If no symbol is specified, psxy expects the last column to have symbol code.\n");
		fprintf (stderr, "\t   [Note: if -C is selected then 3rd means 4th column, etc.]\n");
		fprintf (stderr, "\t   Symbols A, C, D, G, H, I, N, S, T are adjusted to have same area as a circle of given size\n");
		fprintf (stderr, "\t   Bars: Append b<base> to give the y-value of the base [Default = 0]\n");
		fprintf (stderr, "\t      Append u if width is in x-input units [Default is %s]\n", GMT_unit_names[gmtdefs.measure_unit]);
		fprintf (stderr, "\t      Use upper case -SB for horizontal bars (base then refers to x and width may be in y-units [Default is vertical]\n");
		fprintf (stderr, "\t   Ellipses: the direction, major, and minor axis must be in input columns 3, 4 and 5.\n");
		fprintf (stderr, "\t     If -SE rather than -Se is selected, psxy will expect azimuth, and axes in km\n");
		fprintf (stderr, "\t     and convert azimuths based on the chosen map projection.  If projection is linear then\n");
		fprintf (stderr, "\t     we scale the minor/major axes by the map scale.\n");
		fprintf (stderr, "\t   Rotatable Rectangle: Direction, x-dimension, and y-dimension are expected in columns 3, 4, and 5.\n");
		fprintf (stderr, "\t     If -SJ rather than -Sj is selected, psxy will expect azimuth, and dimensions in km\n");
		fprintf (stderr, "\t     and convert azimuths based on the chosen map projection.  If projection is linear then\n");
		fprintf (stderr, "\t     we scale the rectangle dimension by the map scale.\n");
		fprintf (stderr, "\t   Fronts: Give tickgap/ticklen[dir][type][:offset], where\n");
		fprintf (stderr, "\t     dir    = Plot symbol to the l(eft) or r(ight) of front [Default is centered]\n");
		fprintf (stderr, "\t     type   =  b(ox), c(ircle), f(ault), s(lip), t(riangle) [Default is f]\n");
		fprintf (stderr, "\t       box      : square when centered, half-square otherwise.\n"); 
		fprintf (stderr, "\t       circle   : full when centered, half-circle otherwise.\n"); 
		fprintf (stderr, "\t       fault    : centered cross-tick or tick only in the <dir> direction.\n"); 
		fprintf (stderr, "\t       slip     : left-lateral or right-lateral strike-slip arrows (centered is not defined).\n"); 
		fprintf (stderr, "\t       triangle : diagonal square when centered, directed triangle otherwise.\n"); 
		fprintf (stderr, "\t     offset = start plotting symbols when along-track distance equals offset [0].\n"); 
		fprintf (stderr, "\t   Kustom: append <symbolname> immediately after 'k'; this will look for <symbolname>.def in\n");
		fprintf (stderr, "\t     the current directory, in $GMT_USERDIR, or in $GMT_SHAREDIR (in that order).\n");
		GMT_list_custom_symbols ();
		fprintf (stderr, "\t   Letter: append /<string> after symbol size, and optionally %%<font>\n");
		fprintf (stderr, "\t   Mathangle: the start and stop directions of the math angle must be in input columns 3, 4.\n");
		fprintf (stderr, "\t     Use -Smf for arrow at first angle, -Sml for last angle, -Smb for both [none].\n");
		fprintf (stderr, "\t   Quoted line: Give [d|f|n|l|x]<info>[:<labelinfo>]\n");
		fprintf (stderr, "\t     <code><info> controls placement of labels along lines.  Choose among five algorithms:\n");
		GMT_cont_syntax (7, 1);
		fprintf (stderr, "\t     <labelinfo> controls the specifics of the labels.  Append what you need:\n");
		GMT_label_syntax (7, 1);
		fprintf (stderr, "\t   Rectangles: the x- and y-dimension must be in input columns 3 and 4.\n");
		fprintf (stderr, "\t   Vectors: the direction and length must be in input columns 3 and 4.\n");
		fprintf (stderr, "\t     Furthermore, <size> means arrowwidth/headlength/headwith [Default is %g/%g/%g]\n", S.v_width, S.h_length, S.h_width);
		fprintf (stderr, "\t     If -SV rather than -Sv is selected, psxy will expect azimuth and length\n");
		fprintf (stderr, "\t     and convert azimuths based on the chosen map projection\n");
		fprintf (stderr, "\t     Insert h(head), b(balance point), or t(ail) after -Sv|V to justify vector w.r.t. input (x,y).\n");
		fprintf (stderr, "\t     Insert s(egment) if (x,y) is tail and columns 3 and 4 holds head (x,y).\n");
		fprintf (stderr, "\t     Upper case H, B, T, S will give double-headed vector [Default is t]\n");
		fprintf (stderr, "\t   Wedges: the start and stop directions of pie wedge must be in input columns 3, 4.\n");
		fprintf (stderr, "\t     If -SW rather than -Sw is selected, psxy will expect two azimuths instead.\n");
		fprintf (stderr, "\t-T Ignores all input files.\n");
		GMT_explain_option ('U');
		GMT_explain_option ('V');
		GMT_pen_syntax ('W', "sets pen attributes:");
		fprintf (stderr, "\t   The leading + applies cpt color (-C) to both symbol fill and pen; - applies it to pen only.\n");
		GMT_explain_option ('X');
		GMT_explain_option ('c');
		GMT_explain_option (':');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf (stderr, "\t   Default is the required number of columns\n");
		GMT_explain_option ('f');
		GMT_explain_option ('g');
		GMT_explain_option ('m');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	/* Check that the options selected are mutually consistent */

	if (Ctrl->S.active && GMT_parse_symbol_option (Ctrl->S.arg, &S, 0, TRUE)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -S option\n", GMT_program);
		error++;
	}
	if (S.symbol > 0 && Ctrl->C.active) get_rgb = TRUE;	/* Need to read z-vales from input data file */
	if (Ctrl->E.active && (S.read_vector || S.symbol == GMT_SYMBOL_ELLIPSE || S.symbol == GMT_SYMBOL_FRONT || S.symbol == GMT_SYMBOL_QUOTED_LINE || S.symbol == GMT_SYMBOL_ROTATERECT)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -E option: Incompatible with -Se, -Sf, -Sj, -Sq, -Sv\n", GMT_program);
		error++;
	}
	if (Ctrl->C.active && S.symbol == 0 && !GMT_io.multi_segments[GMT_IN]) {
		error++;
		fprintf (stderr, "%s: GMT SYNTAX ERROR -C option: Must also specify a symbol (see -S) or give polygon file with -m\n", GMT_program);
	}
	if (!project_info.region_supplied) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify -R option\n", GMT_program);
		error++;
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
	if (GMT_io.binary[GMT_IN] && S.symbol == GMT_SYMBOL_NOT_SET) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have symbol information\n", GMT_program);
		error++;
	}
	if (Ctrl->E.active && Ctrl->E.mode && !Ctrl->C.active) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  -E option +|-<pen> requires the -C option\n", GMT_program);
		error++;
	}
	if (Ctrl->W.active && Ctrl->W.mode && !Ctrl->C.active) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  -W option +|-<pen> requires the -C option\n", GMT_program);
		error++;
	}
	if ((Ctrl->W.mode + Ctrl->E.mode) == 3) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Conflicting -E and -W options regarding -C option application\n", GMT_program);
		error++;
	}
	if (Ctrl->T.active && n_files) fprintf (stderr, "%s: GMT WARNING.  Option -T ignores all input files\n", GMT_program);

	if (error) exit (EXIT_FAILURE);

	if (Ctrl->E.active) {	/* Set error bar parameters */
		j = 2;	/* Normally, error bar related columns start in position 2 */
		if (Ctrl->E.xbar == 1) {
			xy_errors[GMT_X] = j++;
			error_type[GMT_X] = 0;
		}
		else if (Ctrl->E.xbar == 2) {	/* Box-whisker instead */
			xy_errors[GMT_X] = j++;
			error_type[GMT_X] = 1;
		}
		else if (Ctrl->E.xbar == 3) {	/* Notched Box-whisker instead */
			xy_errors[GMT_X] = j++;
			error_type[GMT_X] = 2;
		}
		if (Ctrl->E.ybar == 1) {
			xy_errors[GMT_Y] = j++;
			error_type[GMT_Y] = 0;
		}
		else if (Ctrl->E.ybar == 2) {	/* Box-whisker instead */
			xy_errors[GMT_Y] = j++;
			error_type[GMT_Y] = 1;
		}
		else if (Ctrl->E.ybar == 3) {	/* Notched Box-whisker instead */
			xy_errors[GMT_Y] = j++;
			error_type[GMT_Y] = 2;
		}
		if (!(xy_errors[GMT_X] || xy_errors[GMT_Y])) {	/* Default is plain error bars for both */
			def_err_xy = TRUE;
			xy_errors[GMT_X] = 2;	/* Assumes xy input, later check for -: */
			xy_errors[GMT_Y] = 3;
			error_type[GMT_X] = error_type[GMT_Y] = 1;
		}
	}

	read_symbol = (S.symbol == GMT_SYMBOL_NOT_SET);	/* Only for ASCII input */

	if (S.symbol == 0 && (Ctrl->G.active || Ctrl->L.active)) polygon = TRUE;
	if (S.symbol == GMT_SYMBOL_FRONT || S.symbol == GMT_SYMBOL_QUOTED_LINE) polygon = FALSE;

	current_pen = default_pen = Ctrl->W.pen;
	current_fill = default_fill = Ctrl->G.fill;
	default_outline = Ctrl->W.active;
	if (Ctrl->I.active) {
		GMT_illuminate (Ctrl->I.value, current_fill.rgb);
		GMT_illuminate (Ctrl->I.value, default_fill.rgb);
	}

	Ctrl->E.size *= 0.5;	/* Since we draw half-way in either direction */

	if (Ctrl->E.active && S.symbol == 0)	/* Assume user only wants error bars */
		S.symbol = GMT_SYMBOL_NONE;

	if (Ctrl->C.active) {
		GMT_read_cpt (Ctrl->C.file);
		if (get_rgb) n_cols_start++;
	}

	/* For most symbols, the data columns beyond two will be dimensions that either have the units appended (e.g., 2c)
	 * or they are assumed to be in the current measure unit (MEASURE_UNIT).  We therefore set the in_col_type to be
	 * GMT_IS_DIMENSION for these so that unit conversions are handled correctly.  However, some symbols also require
	 * angles via the input data file.  S.n_nondim and S.nondim_col are used to reset the in_col_type back to GMT_IS_FLOAT
	 * for those columns are expected to contain angles.  When NO SYMBOL is specified in -S we must parse the ASCII data
	 * record to determine the symbol, and this happens AFTER we have already converted the dimensions.  We must therefore
	 * undo this scaling based on what columns might be angles. */
	
	two   = (get_rgb) ? 3 : 2;
	three = (get_rgb) ? 4 : 3;
	four  = (get_rgb) ? 5 : 4;
	pos2x = two + gmtdefs.xy_toggle[GMT_IN];	/* Column with a 2nd longitude (for VECTORS with two sets of coordinates) */
	pos2y = three - gmtdefs.xy_toggle[GMT_IN];	/* Column with a 2nd latitude (for VECTORS with two sets of coordinates) */
	if (Ctrl->E.active) {
		if (S.read_size) GMT_io.in_col_type[two] = GMT_IS_DIMENSION;	/* Must read symbol size from data record */
		if (xy_errors[GMT_X] && xy_errors[GMT_Y] && error_type[GMT_X] >= 1) xy_errors[GMT_Y] += error_cols[error_type[GMT_X]] - 1;	/* Need 3 or 4 columns for whisker bars */
		if (def_err_xy && gmtdefs.xy_toggle[GMT_IN]) {	/* -E should be -Eyx */
			l_swap (xy_errors[GMT_X], xy_errors[GMT_Y]);
			l_swap (error_type[GMT_X], error_type[GMT_Y]);
		}
		if (xy_errors[GMT_X]) n_cols_start += error_cols[error_type[GMT_X]], error_x = TRUE;
		if (xy_errors[GMT_Y]) n_cols_start += error_cols[error_type[GMT_Y]], error_y = TRUE;
		xy_errors[GMT_X] += (S.read_size + get_rgb);	/* Move 0-2 columns over */
		xy_errors[GMT_Y] += (S.read_size + get_rgb);
	}
	else	/* Here we have the usual x y [z] [size] [other args] [symbol] record */
		for (j = n_cols_start; j < 10; j++) GMT_io.in_col_type[j] = GMT_IS_DIMENSION;		/* Since these may have units appended */
	for (j = 0; j < S.n_nondim; j++) GMT_io.in_col_type[S.nondim_col[j]+get_rgb] = GMT_IS_FLOAT;	/* Since these are angles, not dimensions */
	n_required = n_cols_start + S.n_required;
	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = n_required;
	n_expected = (GMT_io.binary[GMT_IN]) ? GMT_io.ncol[GMT_IN] : n_required;
	if (n_expected < n_required) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data must have at least %ld columns\n", GMT_program, n_required);
		error++;
	}

	if (GMT_io.binary[GMT_IN] && gmtdefs.verbose) {
		char *type[2] = {"double", "single"};
		fprintf (stderr, "%s: Expects %ld-column %s-precision binary data\n", GMT_program, GMT_io.ncol[GMT_IN], type[GMT_io.single_precision[GMT_IN]]);
	}

	GMT_err_fail (GMT_map_setup (west, east, south, north), "");

	if (S.u_set) {	/* When -Sc<unit> is given we temporarily reset the system unit to these units so conversions will work */
		save_u = gmtdefs.measure_unit;
		gmtdefs.measure_unit = S.u;
	}
	
	if (S.symbol == GMT_SYMBOL_QUOTED_LINE) {
		if (GMT_contlabel_prep (&S.G, NULL)) exit (EXIT_FAILURE);
		penset_OK = FALSE;	/* The pen for quoted lines are set within the PSL code itself so we dont do it here in psxy */
	}

	resample = ((!Ctrl->A.active || Ctrl->A.mode) && GMT_IS_MAPPING);
	/* Maximum step size (in degrees) used for interpolation of line segments along great circles (if requested) */
	step = Ctrl->A.step / project_info.x_scale / project_info.DIST_M_PR_DEG;

	GMT_plotinit (argc, argv);

	if (S.symbol == GMT_SYMBOL_FRONT || S.symbol == GMT_SYMBOL_QUOTED_LINE || (S.symbol > 0 && !Ctrl->N.active)) {
		GMT_map_clip_on (GMT_no_rgb, 3);
		clip_set = TRUE;
	}
	if (penset_OK) GMT_setpen (&current_pen);
	if (S.symbol == GMT_SYMBOL_TEXT && Ctrl->G.active && !Ctrl->W.active) ps_setpaint (current_fill.rgb);

	if (S.symbol == GMT_SYMBOL_TEXT) {
		ps_setfont (S.font_no);			/* Set the required font */
	}
	if (S.symbol == GMT_SYMBOL_ELLIPSE) Ctrl->N.active = TRUE;	/* So we can see ellipses that have centers outside -R */
	if (S.symbol == GMT_SYMBOL_POINT && !Ctrl->E.active) {	/* Temporarily change linecap to round (1)*/
		ps_setlinecap (1);
		reset_cap = 1;
	}
	else if ((read_symbol || Ctrl->E.active) && gmtdefs.ps_line_cap != 1) {	/* Must set/rest cap around each ps_point call */
		reset_cap = 2;
	}
	old_GMT_world_map = GMT_world_map;
	if (S.convert_angles && GMT_IS_MAPPING) S.convert_angles = 2;
	if (S.symbol == GMT_SYMBOL_BARX && !S.base_set && project_info.xyz_projection[GMT_X] == GMT_LOG10) S.base = west;	/* Default to west level for horizontal log10 bars */
	if (S.symbol == GMT_SYMBOL_BARY && !S.base_set && project_info.xyz_projection[GMT_Y] == GMT_LOG10) S.base = south;	/* Default to south level for vertical log10 bars */
	if (S.symbol == GMT_SYMBOL_BARX && !S.base_set && project_info.xyz_projection[GMT_X] == GMT_LOG10) S.base = west;	/* Default to west level for horizontal log10 bars */
	bar_clip = ((S.symbol == GMT_SYMBOL_BARX || S.symbol == GMT_SYMBOL_BARY) && !Ctrl->N.active) ? TRUE : FALSE;		/* Bars get clipped instead of skipped */
	if ((S.symbol == GMT_SYMBOL_VECTOR || S.symbol == GMT_SYMBOL_VECTOR2) && S.v_just == 3) {
		/* Reading 2nd coordinate so must set column types */
		GMT_io.in_col_type[pos2x] = GMT_io.in_col_type[GMT_X];
		GMT_io.in_col_type[pos2y] = GMT_io.in_col_type[GMT_Y];
	}
	else if (S.symbol == GMT_SYMBOL_VECTOR || S.symbol == GMT_SYMBOL_VECTOR2) {
		/* Reading 2nd coordinate so must set column types */
		GMT_io.in_col_type[two] = GMT_IS_FLOAT;	/* Direction */
	}
	fill_active = Ctrl->G.active;	/* Make copies because we will change the values */
	outline_active =  Ctrl->W.active;
	if (S.symbol > 0 && !outline_active && !fill_active && !get_rgb) outline_active = TRUE;		/* If no fill nor outline for symbols then turn outline on */

	if (Ctrl->D.active) ps_transrotate (Ctrl->D.dx, Ctrl->D.dy, 0.0);	/* Shift plot a bit */
	done = Ctrl->T.active;	/* When -T is used, skip all input */
	n_args = (argc > 1) ? argc : 2;
	nofile = (n_files == 0);

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
				not_used = GMT_fgets (line, BUFSIZ, fp);
				if (i == 0 && S.G.label_type == 2) GMT_extract_label (&line[1], S.G.label);	/* Set initial header as potential label */
			}
		}
		if (S.symbol > 0) {	/* symbol part (not counting GMT_SYMBOL_FRONT and GMT_SYMBOL_QUOTED_LINE) */

			GMT_world_map = (S.symbol == GMT_SYMBOL_ELLIPSE && S.convert_angles) ? FALSE : TRUE;
			if (read_symbol) n_expected = GMT_MAX_COLUMNS;
			while ((n_fields = GMT_input (fp, &n_expected, &in)) >= 0 && !GMT_REC_IS_EOF) {

				while (GMT_REC_IS_SEG_HEADER && !GMT_REC_IS_EOF) {	/* Process segment headers */
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
	
				if (read_symbol) {	/* The symbol was given as last item on an ASCII data record.  Must parse and set parameters */
					i = strlen (GMT_io.current_record) - 1;
					GMT_io.current_record[i--] = 0;	/* Chop off the \n */
					while (GMT_io.current_record[i] && !strchr (" ,\t", (int)GMT_io.current_record[i])) i--;
					(void ) GMT_parse_symbol_option (&GMT_io.current_record[i+1], &S, 0, FALSE);
					n_required = n_cols_start + S.n_required;
				}
				if (S.read_size) {
					S.size_x = in[two];
					S.size_x2 = 0.5 * S.size_x;
				}

				if (S.symbol == GMT_SYMBOL_BARY) GMT_geo_to_xy (west, S.base, &plot_x, &y0);	/* Zero y level for vertical bars */
				if (S.symbol == GMT_SYMBOL_BARX) GMT_geo_to_xy (S.base, south, &x0, &plot_y);	/* Zero x level for horizontal bars */

				if (get_rgb) {
					GMT_get_rgb_from_z (in[GMT_Z], current_fill.rgb);
					if (Ctrl->I.active) GMT_illuminate (Ctrl->I.value, current_fill.rgb);
					if (S.symbol == GMT_SYMBOL_CROSS || S.symbol == GMT_SYMBOL_PLUS || S.symbol == GMT_SYMBOL_POINT || S.symbol == GMT_SYMBOL_TEXT) ps_setpaint (current_fill.rgb);
				}
				if (GMT_cpt_skip) continue;	/* Chosen cpt file indicates skip for this z */
		
				/* If read_symbol is TRUE we must unto the dimension scaling for items that should have no scaling, e.g., angles, or dimensions given in km which must be projected */
				
				if (read_symbol) for (j = 0; j < S.n_nondim; j++) in[S.nondim_col[j]+get_rgb] *= GMT_u2u[GMT_INCH][S.u];	/* Undo scales since these are angles or km, not dimensions */
				
				if (S.read_vector) {
					direction = in[two];
					length = in[three];
				}
				else if (S.symbol == GMT_SYMBOL_ELLIPSE && S.convert_angles == 2) {	/* Leave axes in km */
					direction = in[two];
					major = in[three];
					minor = in[four];
				}
				else if (S.symbol == GMT_SYMBOL_ELLIPSE && S.convert_angles) {	/* Scale lengths by given scale */
					direction = 90.0 - in[two];	/* Cartesian azimuth */
					major = in[three] * project_info.x_scale;
					minor = in[four] * project_info.x_scale;
				}
				else if (S.symbol == GMT_SYMBOL_ELLIPSE) {
					direction = in[two];
					major = in[three];
					minor = in[four];
				}
				else if (S.symbol == GMT_SYMBOL_ROTATERECT && S.convert_angles == 2) {	/* Leave axes in km */
					direction = in[two];
					x_len = in[three+S.read_size];
					y_len = in[four+S.read_size];
				}
				else if (S.symbol == GMT_SYMBOL_ROTATERECT && S.convert_angles) {	/* Scale lengths by given scale */
					direction = 90.0 - in[two];	/* Cartesian azimuth */
					x_len = in[three+S.read_size] * project_info.x_scale;
					y_len = in[four+S.read_size] * project_info.x_scale;
				}
				else if (S.symbol == GMT_SYMBOL_ROTATERECT) {
					direction = in[two];
					x_len = in[three+S.read_size];
					y_len = in[four+S.read_size];
				}
				else if (S.symbol == GMT_SYMBOL_PIE && S.convert_angles == 1) {	/* Got Cartesian azimuths, get dirs */
					dir2 = 90.0 - in[two+S.read_size];
					dir1 = 90.0 - in[three+S.read_size];
				}
				else if (S.symbol == GMT_SYMBOL_PIE || S.symbol == GMT_SYMBOL_MANGLE) {
					dir1 = in[two+S.read_size];
					dir2 = in[three+S.read_size];
				}
				else if (S.symbol == GMT_SYMBOL_RECT) {
					x_len = in[two+S.read_size];
					y_len = in[three+S.read_size];
				}
		
				/* Skip zero-size symbols */
		
				if (!(S.symbol == GMT_SYMBOL_BARX || S.symbol == GMT_SYMBOL_BARY) && S.symbol < GMT_SYMBOL_ELLIPSE && S.size_x <= 0.0) {
					if (gmtdefs.verbose) fprintf (stderr, "%s: Warning: Symbol size <= 0.0 near line %ld (skipped)\n", GMT_program, n_total_read);
					continue;
				}
		
				if (S.read_vector && S.v_just < 3 && length <= 0.0) {
					if (gmtdefs.verbose) fprintf (stderr, "%s: Warning: Vector length <= 0.0 near line %ld (skipped)\n", GMT_program, n_total_read);
					continue;
				}
				if (S.symbol == GMT_SYMBOL_ELLIPSE && (major <= 0.0 || minor <= 0.0)) {
					if (gmtdefs.verbose) fprintf (stderr, "%s: Warning: Ellipse axes <= 0.0 near line %ld (skipped)\n", GMT_program, n_total_read);
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
				}
				else if (!Ctrl->N.active) {
					GMT_map_outside (in[GMT_X], in[GMT_Y]);
					if (GMT_abs (GMT_x_status_new) > 1 || GMT_abs (GMT_y_status_new) > 1) continue;
				}

				if (GMT_geo_to_xy (in[GMT_X], in[GMT_Y], &plot_x, &plot_y)) continue;	/* NaNs */
				
				if (GMT_is_dnan (plot_x)) {	/* Transformation yielded a NaN (e.g. log (-ve)) */
					fprintf (stderr, "%s: Warning: Data point with x = NaN near line %ld\n", GMT_program, n_total_read);
					continue;
				}
				if (GMT_is_dnan (plot_y)) {
					if (gmtdefs.verbose) fprintf (stderr, "%s: Warning: Data point with y = NaN near line %ld\n", GMT_program, n_total_read);
					continue;
				}
				if (Ctrl->W.mode) memcpy ((void *)Ctrl->W.pen.rgb, (void *)current_fill.rgb, 3*sizeof (int));
				if (Ctrl->E.active) {
					if (Ctrl->E.mode) memcpy ((void *)Ctrl->E.pen.rgb, (void *)current_fill.rgb, 3*sizeof (int));
					if (Ctrl->E.mode & 1) memcpy ((void *)current_fill.rgb, (void *)GMT_no_rgb, 3*sizeof (int));
					GMT_setpen (&Ctrl->E.pen);
					if (error_x) {
						if (error_type[GMT_X] == 0)
							plot_x_errorbar (in[GMT_X], in[GMT_Y], in[xy_errors[GMT_X]], Ctrl->E.size, n_total_read);
						else
							plot_x_whiskerbar (plot_x, in[GMT_Y], &in[xy_errors[GMT_X]], Ctrl->E.size, current_fill.rgb, n_total_read, error_type[GMT_X]);
					}
					if (error_y) {
						if (error_type[GMT_Y] == 0)
							plot_y_errorbar (in[GMT_X], in[GMT_Y], in[xy_errors[GMT_Y]], Ctrl->E.size, n_total_read);
						else
							plot_y_whiskerbar (in[GMT_X], plot_y, &in[xy_errors[GMT_Y]], Ctrl->E.size, current_fill.rgb, n_total_read, error_type[GMT_Y]);
					}
					if (!Ctrl->W.mode) GMT_setpen (&current_pen);
				}
				if (Ctrl->W.mode & 1) memcpy ((void *)current_fill.rgb, (void *)GMT_no_rgb, 3*sizeof (int));
				if (Ctrl->W.mode) GMT_setpen (&Ctrl->W.pen);
				
				switch (S.symbol) {
					case GMT_SYMBOL_NONE:
						break;
					case GMT_SYMBOL_STAR:
						size = (S.equal_area) ? 1.67289326141 * S.size_x : S.size_x;
						GMT_star (plot_x, plot_y, 0.0, &size, &current_fill, outline_active);
						break;
					case GMT_SYMBOL_BARX:
						if (S.user_unit) {	/* Width measured in y units */
							GMT_geo_to_xy (S.base, in[GMT_Y]-S.size_x2, &x_1, &y_1);
							GMT_geo_to_xy (in[GMT_X], in[GMT_Y]+S.size_x2, &x_2, &y_2);
							dim[0] = x_2 - x_1;	dim[1] = y_2 - y_1;
							GMT_rect (x_1, y_1, 0.0, dim, &current_fill, outline_active);
						}
						else {
							dim[1] = S.size_x;	dim[0] = plot_x - x0;
							GMT_rect (x0, plot_y-S.size_x2, 0.0, dim, &current_fill, outline_active);
						}
						break;
					case GMT_SYMBOL_BARY:
						if (S.user_unit) {	/* Width measured in x units */
							GMT_geo_to_xy (in[GMT_X]-S.size_x2, S.base, &x_1, &y_1);
							GMT_geo_to_xy (in[GMT_X]+S.size_x2, in[GMT_Y], &x_2, &y_2);
							dim[0] = x_2 - x_1;	dim[1] = y_2 - y_1;
							GMT_rect (x_1, y_1, 0.0, dim, &current_fill, outline_active);
						}
						else {
							dim[0] = S.size_x;	dim[1] = plot_y - y0;
							GMT_rect (plot_x-S.size_x2, y0, 0.0, dim, &current_fill, outline_active);
						}
						break;
					case GMT_SYMBOL_CROSS:
						ps_cross (plot_x, plot_y, S.size_x);
						break;
					case GMT_SYMBOL_PLUS:
						ps_plus (plot_x, plot_y, S.size_x);
						break;
					case GMT_SYMBOL_POINT:
						if (reset_cap == 2) ps_setlinecap (1);
						ps_point (plot_x, plot_y, S.size_x);
						if (reset_cap == 2) ps_setlinecap (gmtdefs.ps_line_cap);
						break;
					case GMT_SYMBOL_CIRCLE:
						GMT_circle (plot_x, plot_y, 0.0, &S.size_x, &current_fill, outline_active);
						break;
					case GMT_SYMBOL_SQUARE:
						size = (S.equal_area) ? 1.25331413732 * S.size_x : S.size_x;
						GMT_square (plot_x, plot_y, 0.0, &size, &current_fill, outline_active);
						break;
					case GMT_SYMBOL_HEXAGON:
						size = (S.equal_area) ? 1.09963611079 * S.size_x : S.size_x;
						GMT_hexagon (plot_x, plot_y, 0.0, &size, &current_fill, outline_active);
						break;
					case GMT_SYMBOL_MANGLE:
						dim[0] = S.size_x2; dim[1] = dir1; dim[2] = dir2;
						GMT_matharc (plot_x, plot_y, 0.0, dim, gmtdefs.vector_shape, &current_pen, S.v_double_heads);
						break;
					case GMT_SYMBOL_PENTAGON:
						size = (S.equal_area) ? 1.14948092619 * S.size_x : S.size_x;
						GMT_pentagon (plot_x, plot_y, 0.0, &size, &current_fill, outline_active);
						break;
					case GMT_SYMBOL_OCTAGON:
						size = (S.equal_area) ? 1.05390736526 * S.size_x : S.size_x;
						GMT_octagon (plot_x, plot_y, 0.0, &size, &current_fill, outline_active);
						break;
					case GMT_SYMBOL_TRIANGLE:
						size = (S.equal_area) ? 1.55512030156 * S.size_x : S.size_x;
						GMT_triangle (plot_x, plot_y, 0.0, &size, &current_fill, outline_active);
						break;
					case GMT_SYMBOL_ITRIANGLE:
						size = (S.equal_area) ? 1.55512030156 * S.size_x : S.size_x;
						GMT_itriangle (plot_x, plot_y, 0.0, &size, &current_fill, outline_active);
						break;
					case GMT_SYMBOL_DIAMOND:
						size = (S.equal_area) ? 1.25331413732 * S.size_x : S.size_x;
						GMT_diamond (plot_x, plot_y, 0.0, &size, &current_fill, outline_active);
						break;
					case GMT_SYMBOL_TEXT:
						font_size = S.size_x * 72.0;
						if (outline_active && Ctrl->G.active) {
							ps_setpaint (current_fill.rgb);
							ps_text (plot_x, plot_y, font_size, S.string, 0.0, 6, FALSE);
							ps_setpaint (current_pen.rgb);
							ps_text (plot_x, plot_y, font_size, S.string, 0.0, 6, TRUE);
						}
						else if (Ctrl->G.active)
							ps_text (plot_x, plot_y, font_size, S.string, 0.0, 6, FALSE);
						else
							ps_text (plot_x, plot_y, font_size, S.string, 0.0, 6, TRUE);
						break;
					case GMT_SYMBOL_ELLIPSE:
						if (S.convert_angles == 2)
							GMT_plot_ellipse (in[GMT_X], in[GMT_Y], 0.0, major, minor, direction, current_fill, outline_active);
						else {
							dim[0] = direction; dim[1] = major; dim[2] = minor;
							GMT_ellipse (plot_x, plot_y, 0.0, dim, &current_fill, outline_active);
						}
						break;
					case GMT_SYMBOL_ROTATERECT:
						if (S.convert_angles == 2)
							GMT_plot_rectangle (in[GMT_X], in[GMT_Y], 0.0, x_len, y_len, direction, current_fill, outline_active);
						else {
							dim[0] = direction; dim[1] = x_len; dim[2] = y_len;
							GMT_rotrect (plot_x, plot_y, 0.0, dim, &current_fill, outline_active);
						}
						break;
					case GMT_SYMBOL_VECTOR:
						if (S.convert_angles == 2) {
							GMT_azim_to_angle (in[GMT_X], in[GMT_Y], 0.1, direction, &tmp);
							direction = tmp;
						}
						else if (S.convert_angles == 1)	/* Cartesian azimuth */
							direction = 90.0 - direction;
						if (S.v_just == 3) {
							GMT_geo_to_xy (in[pos2x], in[pos2y], &x2, &y2);
							if (GMT_is_dnan (x2) || GMT_is_dnan (y2)) {
								fprintf (stderr, "%s: Warning: Vector head coordinates contain NaNs near line %ld. Skipped\n", GMT_program, n_total_read);
								continue;
							}
						}
						else {
							sincosd (direction, &s, &c);
							x2 = plot_x + length * c;
							y2 = plot_y + length * s;
							if (S.v_just) {
								dx = S.v_just * 0.5 * (x2 - plot_x);	dy = S.v_just * 0.5 * (y2 - plot_y);
								plot_x -= dx;		plot_y -= dy;
								x2 -= dx;		y2 -= dy;
							}
						}
						this_outline = (S.v_double_heads) ? outline_active + 8 : outline_active;
						GMT_vector (plot_x, plot_y, x2, y2, 0.0, S.v_width, S.h_length, S.h_width, gmtdefs.vector_shape, &current_fill, (GMT_LONG)this_outline);
						break;
					case GMT_SYMBOL_VECTOR2:
						if (S.convert_angles == 2) {
							GMT_azim_to_angle (in[GMT_X], in[GMT_Y], 1.0, direction, &tmp);
							direction = tmp;
						}
						else if (S.convert_angles == 1)	/* Cartesian azimuth */
							direction = 90.0 - direction;
						if (S.v_just == 3) {
							GMT_geo_to_xy (in[pos2x], in[pos2y], &x2, &y2);
							if (GMT_is_dnan (x2) || GMT_is_dnan (y2)) {
								fprintf (stderr, "%s: Warning: Vector head coordinates contain NaNs near line %ld. Skipped\n", GMT_program, n_total_read);
								continue;
							}
						}
						else {
							sincosd (direction, &s, &c);
							x2 = plot_x + length * c;
							y2 = plot_y + length * s;
							if (S.v_just) {
								dx = S.v_just * 0.5 * (x2 - plot_x);	dy = S.v_just * 0.5 * (y2 - plot_y);
								plot_x -= dx;		plot_y -= dy;
								x2 -= dx;		y2 -= dy;
							}
						}
						this_outline = (S.v_double_heads) ? outline_active + 8 : outline_active;
						if (length < S.v_norm) {	/* Scale arrow attributes down with length */
							v_w = S.v_width * length * S.v_shrink;
							h_l = S.h_length * length * S.v_shrink;
							h_w = S.h_width * length * S.v_shrink;
							GMT_vector (plot_x, plot_y, x2, y2, 0.0, v_w, h_l, h_w, gmtdefs.vector_shape, &current_fill, (GMT_LONG)this_outline);
						}
						else	/* Leave as specified */
							GMT_vector (plot_x, plot_y, x2, y2, 0.0, S.v_width, S.h_length, S.h_width, gmtdefs.vector_shape, &current_fill, (GMT_LONG)this_outline);
						break;
					case GMT_SYMBOL_PIE:
						if (S.convert_angles == 2) {
							GMT_azim_to_angle (in[GMT_X], in[GMT_Y], 1.0, dir1, &tmp);
							dir1 = tmp;
							GMT_azim_to_angle (in[GMT_X], in[GMT_Y], 1.0, dir2, &tmp);
							dir2 = tmp;
							d_swap (dir1, dir2);
						}
						dim[0] = S.size_x2; dim[1] = dir1; dim[2] = dir2;
						GMT_pie (plot_x, plot_y, 0.0, dim, &current_fill, outline_active);
						break;
					case GMT_SYMBOL_XDASH:
						ps_segment (plot_x - S.size_x2, plot_y, plot_x + S.size_x2, plot_y);
						break;
					case GMT_SYMBOL_YDASH:
						ps_segment (plot_x, plot_y - S.size_x2, plot_x, plot_y + S.size_x2);
						break;
					case GMT_SYMBOL_RECT:
						dim[0] = x_len; dim[1] = y_len;
						GMT_rect (plot_x - 0.5*x_len, plot_y - 0.5*y_len, 0.0, dim, &current_fill, outline_active);
						break;
					case GMT_SYMBOL_CUSTOM:
						GMT_draw_custom_symbol (plot_x, plot_y, 0.0, S.size_x, S.custom, &current_pen, &current_fill, outline_active);
						break;
				}
				if (read_symbol) n_expected = GMT_MAX_COLUMNS;
			}
		}
		else {	/* Line/polygon part */

			n_required = 2;
			n_fields = GMT_input (fp, &n_expected, &in);
			while (!GMT_REC_IS_EOF) {	/* Not yet EOF */
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
					if (n == n_alloc) n_alloc = GMT_alloc_memory2 ((void **)&xx, (void **)&yy, n, n_alloc, sizeof (double), GMT_program);
					xx[n] = in[GMT_X];	yy[n] = in[GMT_Y];	n++;
					n_fields = GMT_input (fp, &n_expected, &in);
				}
				
				/* Done reading a segment */
				
				if (polygon && GMT_polygon_is_open (xx, yy, n)) {	/* Explicitly close polygon so that arc will work */
					if (n == n_alloc) n_alloc = GMT_alloc_memory2 ((void **)&xx, (void **)&yy, n, n_alloc, sizeof (double), GMT_program);
					xx[n] = xx[0];	yy[n] = yy[0];	n++;
				}
				n_alloc = GMT_alloc_memory2 ((void **)&xx, (void **)&yy, 0, n, sizeof (double), GMT_program);
		
				if (GMT_cpt_skip) continue;	/* CPT says skip this polygon */
				
				if (resample) n = n_alloc = GMT_fix_up_path (&xx, &yy, n, step, Ctrl->A.mode);	/* Resample if spacing is too coarse */

				if (polygon)	/* Want a filled polygon */
					GMT_fill_polygon (xx, yy, 0.0, n, &current_fill, outline_active);	/* NO outline since path may have been clipped */
				else {		/* Want a line of some sort */
					if ((plot_n = GMT_geo_to_xy_line (xx, yy, n)) == 0) continue;		/* Nothing there there */
					if (S.symbol == GMT_SYMBOL_QUOTED_LINE) {	/* Labelled lines are dealt with by the contour machinery */
						GMT_LONG start = 0, stop;
						S.G.line_pen = current_pen;
						/* Since the line might exit and enter the region we must isolate the segments before calling GMT_hold_contour */
						while (start < plot_n) {	/* More segments to register */
							while (GMT_pen[start+1] == 3 && start < (plot_n-1)) start++;	/* Skip multiple start points */
							stop = start + 1;
							while (GMT_pen[stop] == 2 && stop < plot_n) stop++;	/* Find the end of this segment */
							n = stop - start;
							xp = (double *) GMT_memory (VNULL, (size_t)n, sizeof (double), GMT_program);
							yp = (double *) GMT_memory (VNULL, (size_t)n, sizeof (double), GMT_program);
							memcpy ((void *)xp, (void *)&GMT_x_plot[start], (size_t)(n*sizeof (double)));
							memcpy ((void *)yp, (void *)&GMT_y_plot[start], (size_t)(n*sizeof (double)));
							GMT_hold_contour (&xp, &yp, n, 0.0, "N/A", 'A', S.G.label_angle, Ctrl->L.active, &S.G);
							start = stop;
							GMT_free ((void *)xp);
							GMT_free ((void *)yp);
						}
					}
					else {
						GMT_plot_line (GMT_x_plot, GMT_y_plot, GMT_pen, plot_n);
						if (S.symbol == GMT_SYMBOL_FRONT) GMT_draw_fence (GMT_x_plot, GMT_y_plot, 0.0, plot_n, &S.f, &current_fill, outline_active); 	/* Must draw front crossbars */
					}
				}
			}
		}
		if (fp != GMT_stdin) GMT_fclose (fp);
	}

	if (S.u_set) gmtdefs.measure_unit = save_u;	/* Reset unit */

	if (reset_cap == 1) ps_setlinecap (gmtdefs.ps_line_cap);	/* Reset linecap to default */

	if (S.symbol == GMT_SYMBOL_QUOTED_LINE) {
		GMT_contlabel_plot (&S.G);
		GMT_contlabel_free (&S.G);
	}

	if (Ctrl->D.active) ps_transrotate (-Ctrl->D.dx, -Ctrl->D.dy, 0.0);	/* Reset shift */

	if (clip_set) GMT_map_clip_off ();

	GMT_world_map = old_GMT_world_map;

	GMT_map_basemap ();

	GMT_plotend ();

	if (S.symbol <= 0) {
		GMT_free ((void *)xx);
		GMT_free ((void *)yy);
	}

	Free_psxy_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void plot_x_errorbar (double x, double y, double delta_x, double error_width2, GMT_LONG line) {
	double x_1, x_2, y_1, y_2;
	GMT_LONG tip1, tip2;

	tip1 = tip2 = (error_width2 > 0.0);
	GMT_geo_to_xy (x - delta_x, y, &x_1, &y_1);
	GMT_geo_to_xy (x + delta_x, y, &x_2, &y_2);
	if (GMT_is_dnan (x_1)) {
		fprintf (stderr, "%s: Warning: X error bar exceeded domain near line %ld. Set to x_min\n", GMT_program, line);
		x_1 = project_info.xmin;
		tip1 = FALSE;
	}
	if (GMT_is_dnan (x_2)) {
		fprintf (stderr, "%s: Warning: X error bar exceeded domain near line %ld. Set to x_max\n", GMT_program, line);
		x_2 = project_info.xmax;
		tip2 = FALSE;
	}
	ps_segment (x_1, y_1, x_2, y_2);
	if (tip1) ps_segment (x_1, y_1 - error_width2, x_1, y_1 + error_width2);
	if (tip2) ps_segment (x_2, y_2 - error_width2, x_2, y_2 + error_width2);
}

void plot_y_errorbar (double x, double y, double delta_y, double error_width2, GMT_LONG line) {
	double x_1, x_2, y_1, y_2;
	GMT_LONG tip1, tip2;

	tip1 = tip2 = (error_width2 > 0.0);
	GMT_geo_to_xy (x, y - delta_y, &x_1, &y_1);
	GMT_geo_to_xy (x, y + delta_y, &x_2, &y_2);
	if (GMT_is_dnan (y_1)) {
		fprintf (stderr, "%s: Warning: Y error bar exceeded domain near line %ld. Set to y_min\n", GMT_program, line);
		y_1 = project_info.ymin;
		tip1 = FALSE;
	}
	if (GMT_is_dnan (y_2)) {
		fprintf (stderr, "%s: Warning: Y error bar exceeded domain near line %ld. Set to y_max\n", GMT_program, line);
		y_2 = project_info.ymax;
		tip2 = FALSE;
	}
	ps_segment (x_1, y_1, x_2, y_2);
	if (tip1) ps_segment (x_1 - error_width2, y_1, x_1 + error_width2, y_1);
	if (tip2) ps_segment (x_2 - error_width2, y_2, x_2 + error_width2, y_2);
}

void plot_x_whiskerbar (double x, double y, double hinge[], double error_width2, int rgb[], GMT_LONG line, GMT_LONG kind) {
	GMT_LONG i;
	static GMT_LONG q[4] = {0, 25, 75, 100};
	double xx[4], yy[4];

	for (i = 0; i < 4; i++) {	/* for 0, 25, 75, 100% hinges */
		GMT_geo_to_xy (hinge[i], y, &xx[i], &yy[i]);
		if (GMT_is_dnan (xx[i])) {
			fprintf (stderr, "%s: Warning: X %ld %% hinge exceeded domain near line %ld\n", GMT_program, q[i], line);
			xx[i] = (i <2 ) ? project_info.xmin :  project_info.xmax;
		}
	}
	yy[1] -= error_width2;
	yy[2] += error_width2;

	ps_segment (xx[0], yy[1], xx[0], yy[2]);		/* Left whisker */
	ps_segment (xx[0], yy[0], xx[1], yy[0]);

	ps_segment (xx[3], yy[1], xx[3], yy[2]);		/* Right whisker */
	ps_segment (xx[3], yy[0], xx[2], yy[0]);

	if (kind == 2) {	/* Notched box-n-whisker plot */
		double xp[10], yp[10], s, p;
		s = 1.7 * ((1.25 * (xx[2] - xx[1])) / (1.35 * hinge[4]));	/* 5th term in hinge has n */
		xp[0] = xp[9] = xx[1];
		xp[1] = xp[8] = ((p = (x - s)) < xp[0]) ? xp[0] : p;
		xp[2] = xp[7] = x;
		xp[4] = xp[5] = xx[2];
		xp[3] = xp[6] = ((p = (x + s)) > xp[4]) ? xp[4] : p;
		yp[0] = yp[1] = yp[3] = yp[4] = yy[1];
		yp[5] = yp[6] = yp[8] = yp[9] = yy[2];
		yp[2] = yy[0] - 0.5 * error_width2;
		yp[7] = yy[0] + 0.5 * error_width2;
		ps_patch (xp, yp, (GMT_LONG)10, rgb, TRUE);
		ps_segment (x, yp[7], x, yp[2]);			/* Median line */
	}
	else {
		ps_rect (xx[1], yy[1], xx[2], yy[2], rgb, TRUE);	/* Main box */
		ps_segment (x, yy[1], x, yy[2]);			/* Median line */
	}
}

void plot_y_whiskerbar (double x, double y, double hinge[], double error_width2, int rgb[], GMT_LONG line, GMT_LONG kind) {
	GMT_LONG i;
	static GMT_LONG q[4] = {0, 25, 75, 100};
	double xx[4], yy[4];

	for (i = 0; i < 4; i++) {	/* for 0, 25, 75, 100% hinges */
		GMT_geo_to_xy (x, hinge[i], &xx[i], &yy[i]);
		if (GMT_is_dnan (yy[i])) {
			fprintf (stderr, "%s: Warning: Y %ld %% hinge exceeded domain near line %ld\n", GMT_program, q[i], line);
			yy[i] = (i <2 ) ? project_info.ymin :  project_info.ymax;
		}
	}
	xx[1] -= error_width2;
	xx[2] += error_width2;

	ps_segment (xx[1], yy[0], xx[2], yy[0]);		/* bottom whisker */
	ps_segment (xx[0], yy[0], xx[0], yy[1]);

	ps_segment (xx[1], yy[3], xx[2], yy[3]);		/* Top whisker */
	ps_segment (xx[0], yy[3], xx[0], yy[2]);

	if (kind == 2) {	/* Notched box-n-whisker plot */
		double xp[10], yp[10], s, p;
		s = 1.7 * ((1.25 * (yy[2] - yy[1])) / (1.35 * hinge[4]));	/* 5th term in hinge has n */
		xp[0] = xp[1] = xp[3] = xp[4] = xx[2];
		xp[5] = xp[6] = xp[8] = xp[9] = xx[1];
		xp[2] = xx[0] + 0.5 * error_width2;
		xp[7] = xx[0] - 0.5 * error_width2;
		yp[0] = yp[9] = yy[1];
		yp[1] = yp[8] = ((p = (y - s)) < yp[0]) ? yp[0] : p;
		yp[2] = yp[7] = y;
		yp[4] = yp[5] = yy[2];
		yp[3] = yp[6] = ((p = (y + s)) > yp[4]) ? yp[4] : p;
		ps_patch (xp, yp, (GMT_LONG)10, rgb, TRUE);
		ps_segment (xp[7], y, xp[2], y);			/* Median line */
	}
	else {
		ps_rect (xx[2], yy[2], xx[1], yy[1], rgb, TRUE);	/* Main box */
		ps_segment (xx[1], y, xx[2], y);			/* Median line */
	}
}

void *New_psxy_Ctrl () {	/* Allocate and initialize a new control structure */
	struct PSXY_CTRL *C;
	
	C = (struct PSXY_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct PSXY_CTRL), "New_psxy_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
		
	GMT_init_pen (&C->W.pen, GMT_PENWIDTH);
	GMT_init_pen (&C->E.pen, GMT_PENWIDTH);
	GMT_init_fill (&C->G.fill, -1, -1, -1);	/* Default is no fill */
	C->A.step = gmtdefs.line_step;
	C->E.size = 0.1;
	return ((void *)C);
}

void Free_psxy_Ctrl (struct PSXY_CTRL *C) {	/* Deallocate control structure */
	if (C->C.file) free ((void *)C->C.file);	
	GMT_free ((void *)C);	
}
