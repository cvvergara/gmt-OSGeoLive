/*--------------------------------------------------------------------
 *	$Id: pstext.c 10071 2013-07-06 02:31:21Z pwessel $
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
 * pstext will read (x, y, angle, size, fontno, justify, text) from GMT_stdin
 * or file and plot the textstrings at (x,y) on a map using the size, font,
 * and justification selected by the user.  Alternatively (with -H), read
 * one or more text paragraphs to be typeset.
 *
 * Author:	Paul Wessel
 * Date:	21-JAN-1991-2000
 * Version:	2.0 based on old v1.x
 * Version:	3.1 based on old 3.0
 * Version:	3.3: PW: Introducing paragraph mode -m
 * Version:	3.3.4: PW: Input x,y can be in dd:mm[:ss] format
 * Version:	4.2.x: Now can handle @_, @:, @; everywhere and @;colorname; is OK
 *
 */

#include "gmt.h"
#include "pslib.h"

struct PSTEXT_CTRL {
	struct A {	/* -A */
		GMT_LONG active;
	} A;
	struct C {	/* -C<dx>/<dy> */
		GMT_LONG active;
		GMT_LONG percent;
		double dx, dy;
	} C;
	struct D {	/* -D[j]<dx>/<dy>[v[<pen>] */
		GMT_LONG active;
		GMT_LONG justify;
		GMT_LONG line;
		double dx, dy;
		struct GMT_PEN pen;
	} D;
	struct E {	/* -Eazim/elev */
		GMT_LONG active;
		double azimuth, elevation;
	} E;
	struct G {	/* -G<fill> */
		GMT_LONG active;
		struct GMT_FILL fill;
	} G;
	struct L {	/* -L */
		GMT_LONG active;
	} L;
	struct N {	/* -N */
		GMT_LONG active;
	} N;
	struct S {	/* -S<pen> */
		GMT_LONG active;
		struct GMT_PEN pen;
	} S;
	struct W {	/* -W[<fill>][o|O|c[<pen>] */
		GMT_LONG active;
		GMT_LONG outline;
		GMT_LONG paint;
		GMT_LONG mode;
		struct GMT_PEN pen;
		struct GMT_FILL fill;
	} W;
	struct Z {	/* -Z<z_level> */
		GMT_LONG active;
		GMT_LONG mode;	/* 0 for common z_level, 1 for optional z-value in 3rd col for each record */
		double level;
	} Z;
};

struct PSTEXT_INFO {
	GMT_LONG text_justify;
	GMT_LONG paragraph_font;
	GMT_LONG block_justify;
	GMT_LONG boxflag;
	GMT_LONG space_flag;
	double x_offset, y_offset;	/* Offset from reference point */
	double line_spacing;
	double paragraph_width;
	double font_size;
	double paragraph_angle;
	double x_space, y_space;	/* Extra spacing between box and text */
	struct GMT_PEN boxpen;
	struct GMT_PEN vecpen;
	struct GMT_PEN txtpen;
	struct GMT_FILL txtfill;
	struct GMT_FILL boxfill;
};

int main(int argc, char **argv)
{
	GMT_LONG i, j, k, last_font, fno, n_files = 0, n_args, nscan;
	GMT_LONG ix, iy, *use_font = NULL, n_paragraphs = 0, pos;
	GMT_LONG n_words = 0, n_read = 0, n_processed = 0, n_alloc = 0;

	double xy[2], west = 0.0, east = 0.0, south = 0.0, north = 0.0, plot_x = 0.0, plot_y = 0.0;
	double z_level = 0.0, tmp, dx, dy, xx[2], yy[2], z_level_orig = 0.0;

	GMT_LONG error = FALSE, nofile = TRUE, master_record = FALSE;
	GMT_LONG do_paragraphs = FALSE, done = FALSE, old_GMT_world_map, add;
	GMT_LONG skip_text_records = FALSE;

	char text[BUFSIZ], line[BUFSIZ], this_font[GMT_LONG_TEXT], just_key[5], txt_a[GMT_LONG_TEXT], txt_b[GMT_LONG_TEXT];
	char txt_x[GMT_LONG_TEXT], txt_y[GMT_LONG_TEXT], txt_z[GMT_LONG_TEXT], p[BUFSIZ];
	char **word = (char **)NULL, pjust;

	FILE *fp = NULL;

	struct PSTEXT_INFO T;
	struct PSTEXT_CTRL *Ctrl = NULL;

	void GMT_putwords (double x, double y, char **text, GMT_LONG n_words, struct PSTEXT_INFO *info);
	void *New_pstext_Ctrl (), Free_pstext_Ctrl (struct PSTEXT_CTRL *C);
	void load_parameters (struct PSTEXT_INFO *T, struct PSTEXT_CTRL *C);
	
	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct PSTEXT_CTRL *)New_pstext_Ctrl ();	/* Allocate and initialize a new control structure */

	memset ((void *)&T, 0, sizeof (struct PSTEXT_INFO));

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
				case 'c':
				case 'f':
				case 'm':
				case ':':
				case '\0':
					error += GMT_parse_common_options (argv[i], &west, &east, &south, &north);
					break;

				/* Supplemental parameters */

				case 'A':	/* Getting azimuths rather than directions, must convert via map projection */
					Ctrl->A.active = TRUE;
					break;
				case 'C':
					Ctrl->C.active = TRUE;
					if (argv[i][2]) {	/* Replace default settings with user settings */
						Ctrl->C.percent = (strchr (argv[i], '%')) ? TRUE : FALSE;
						k = sscanf (&argv[i][2], "%[^/]/%s", txt_a, txt_b);
						for (j = 0; txt_a[j]; j++) if (txt_a[j] == '%') txt_a[j] = '\0';	/* Remove % signs before processing values */
						for (j = 0; k == 2 && txt_b[j]; j++) if (txt_b[j] == '%') txt_b[j] = '\0';
						Ctrl->C.dx = GMT_convert_units (txt_a, GMT_INCH);
						Ctrl->C.dy = (k == 2) ? GMT_convert_units (txt_b, GMT_INCH) : Ctrl->C.dx;
					}
					break;
				case 'D':
					Ctrl->D.active = TRUE;
					k = 2;
					if (argv[i][k] == 'j') Ctrl->D.justify = TRUE, k++;
					for (j = k; argv[i][j] && argv[i][j] != 'v'; j++);
					if (argv[i][j] == 'v') {
						Ctrl->D.line = TRUE;
						if (argv[i][j+1] && GMT_getpen (&argv[i][j+1], &Ctrl->D.pen)) {
							fprintf (stderr, "%s: GMT SYNTAX ERROR -D option: Give pen after c\n", GMT_program);
							error++;
						}
						argv[i][j] = 0;
					}
					j = sscanf (&argv[i][k], "%[^/]/%s", txt_a, txt_b);
					Ctrl->D.dx = GMT_convert_units (txt_a, GMT_INCH);
					Ctrl->D.dy = (j == 2) ? GMT_convert_units (txt_b, GMT_INCH) : Ctrl->D.dx;
					break;
				case 'E':
					Ctrl->E.active = TRUE;
					error += GMT_get_proj3D (&argv[i][2], &Ctrl->E.azimuth, &Ctrl->E.elevation);
					break;
				case 'G':
					Ctrl->G.active = TRUE;
					if (GMT_getfill (&argv[i][2], &Ctrl->G.fill)) {
						GMT_fill_syntax ('G', " ");
						error++;
					}
					break;
				case 'L':
					Ctrl->L.active = TRUE;
					break;
				case 'N':	/* Do not clip at border */
					Ctrl->N.active = TRUE;
					break;
				case 'S':
					Ctrl->S.active = TRUE;
					if (GMT_getpen (&argv[i][2], &Ctrl->S.pen)) {
						GMT_pen_syntax ('S', " ");
						error++;
					}
					break;
				case 'W':
					Ctrl->W.active = TRUE;
					/* Do we have -W<fill>,o|O|c|C<pen> ? */
					for (j = 3; argv[i][j] && (argv[i][j-1] != ',' || !(strchr("cCoO", argv[i][j]))); j++);
					if (argv[i][j]) {	/* Gave comma separated fill and outline pens */
						Ctrl->W.outline = TRUE;
						Ctrl->W.mode = argv[i][j];
						if (argv[i][j+1] && GMT_getpen (&argv[i][j+1], &Ctrl->W.pen)) {
							fprintf (stderr, "%s: GMT SYNTAX ERROR -W option: Bad pen given after %c\n", GMT_program, (int)Ctrl->W.mode);
							error++;
						}
						if (j > 3) {	/* Gave fill */
							char txt_a[GMT_LONG_TEXT];
							strncpy (txt_a, &argv[i][2], (size_t)(j-3));
							txt_a[j-3] = '\0';
							if (GMT_getfill (txt_a, &Ctrl->W.fill)) {
								GMT_fill_syntax ('W', " ");
								error++;
							}
							Ctrl->W.paint = TRUE;
						}
					}
					else {	/* Backwards compatible way of doing things - will fail if things like -Wcyanored is given */
						for (j = 2; argv[i][j] && !(strchr ("cCoO", argv[i][j])); j++);
						if (argv[i][j])  {	/* Found c|C|o|O, but this could be a mistake since -Wcyan or -Worange will trigger this - must check further */
							char txt_a[GMT_LONG_TEXT];
							int n;
							n = 2;
							while (!(argv[i][n] == '/' || argv[i][n] == '\0')) n++;
							strncpy (txt_a, &argv[i][2], (size_t)(n-2));	txt_a[n-2] = '\0';
							if (GMT_colorname2index (txt_a) >= 0) {	/* Found a colorname; thus this was actually a fill statement instead */
								if (GMT_getfill (txt_a, &Ctrl->W.fill)) {
									GMT_fill_syntax ('W', " ");
									error++;
								}
								Ctrl->W.paint = TRUE;
							}
							else {	/* Presumably we did get a pen outline request */
								Ctrl->W.outline = TRUE;	/* Want box outline */
								Ctrl->W.mode = argv[i][j];
								if (argv[i][j+1] && GMT_getpen (&argv[i][j+1], &Ctrl->W.pen)) {
									fprintf (stderr, "%s: GMT SYNTAX ERROR -W option: Bad pen given after %c\n", GMT_program, (int)Ctrl->W.mode);
									error++;
								}
								argv[i][j] = 0;
							}
						}
						if (argv[i][2]) {	/* Also gave fill */
							if (GMT_getfill (&argv[i][2], &Ctrl->W.fill)) {
								GMT_fill_syntax ('W', " ");
								error++;
							}
							Ctrl->W.paint = TRUE;
						}
					}
					break;
				case 'Z':
					Ctrl->Z.active = TRUE;
					if (argv[i][2] == '+' && argv[i][3] == '\0') /* Read level for each item */
						Ctrl->Z.mode = 1;
					else /* Gave specific constant level */
						Ctrl->Z.level = atof (&argv[i][2]);
					break;

				/* Options not recognized */

				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else
			n_files++;
	}

	if (argc == 1 || Ctrl->L.active || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "pstext %s - To plot text on maps\n\n", GMT_VERSION);
		fprintf (stderr, "usage: pstext <txtfile> %s %s\n", GMT_J_OPT, GMT_Rgeoz_OPT);
		fprintf (stderr, "\t[-A] [%s] [-C<dx>/<dy>] [-D[j]<dx>[/<dy>][v[<pen>]] [%s] [-G<color>]\n", GMT_B_OPT, GMT_E_OPT);
		fprintf (stderr, "\t[%s] [-K] [-L] [-N] [-O] [-P] [-S<pen>] [%s]\n", GMT_Ho_OPT, GMT_U_OPT);
		fprintf (stderr, "\t[-V] [-W[<fill>,][o|O|c|C[<pen>]]] [%s] [%s]\n", GMT_X_OPT, GMT_Y_OPT);
		fprintf (stderr, "\t[-Z[<zlevel>|+]] [%s] [%s] [%s] [%s]\n\n", GMT_t_OPT, GMT_c_OPT, GMT_f_OPT, GMT_mo_OPT);
		fprintf (stderr, "\tReads (x,y,size,angle,fontno,justify,text) from <txtfile> [or stdin]\n");
		fprintf (stderr, "\tOR (with -m) one or more text paragraphs with formatting info in the segment header.\n");
		fprintf (stderr, "\tjustify is of the form [T|M|B][L|C|R] (top/middle/bottom/left/center/right).\n");
		fprintf (stderr, "\tBuilt-in escape sequences:\n");
		fprintf (stderr, "\t   @~ toggles between current font and Symbol font.\n");
		fprintf (stderr, "\t   @%%<no>%% switches to font number <no>; @%%%% resets font.\n");
		fprintf (stderr, "\t   @:<size>: switches font size; @:: resets font size.\n");
		fprintf (stderr, "\t   @;<color>; switches font color; @;; resets font color.\n");
		fprintf (stderr, "\t   @+ toggles between normal and superscript mode.\n");
		fprintf (stderr, "\t   @- toggles between normal and subscript mode.\n");
		fprintf (stderr, "\t   @# toggles between normal and Small Caps mode.\n");
		fprintf (stderr, "\t   @_ toggles between normal and underlined text.\n");
		fprintf (stderr, "\t   @!<char1><char2> makes one composite character.\n");
		fprintf (stderr, "\t   @@ prints the @ sign itself.\n");
		fprintf (stderr, "\t   Use @a, @c, @e, @n, @o, @s, @u, @A, @C @E, @N, @O, @U for accented European characters.\n");
		fprintf (stderr, "\t(See manual page for more information).\n");

		if (Ctrl->L.active) {	/* List fonts */
			fprintf (stderr, "\n\tFont #	Font Name\n");
			fprintf (stderr, "\t------------------------------------\n");
			for (i = 0; i < GMT_N_FONTS; i++)
				fprintf (stderr, "\t%3ld\t%s\n", i, GMT_font[i].name);
		}

		if (Ctrl->L.active || GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		GMT_explain_option ('j');
		GMT_explain_option ('R');
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-A Angles given as azimuths; convert to directions using current projection.\n");
		GMT_explain_option ('b');
		fprintf (stderr, "\t-C sets the clearance between characters and surrounding box.  Only used\n");
		fprintf (stderr, "\t   if -W has been set.  Append units {c,i,m,p} or %% of fontsize [15%%].\n");
		fprintf (stderr, "\t-D adds <add_x>,<add_y> to the text origin AFTER projecting with -J [0/0].\n");
		fprintf (stderr, "\t   Use -Dj to move text origin away from point (direction determined by text's justification).\n");
		fprintf (stderr, "\t   Append v[<pen>] to draw line from text to original point.  If <add_y> is not given it equal <add_x>.\n");
		GMT_explain_option ('E');
		fprintf (stderr, "\t-G set the color (red/green/blue (0-255)) for solid text [%d/%d/%d].\n", gmtdefs.basemap_frame_rgb[0], gmtdefs.basemap_frame_rgb[1], gmtdefs.basemap_frame_rgb[2]);
		GMT_explain_option ('H');
		GMT_explain_option ('K');
		fprintf (stderr, "\t-L lists the font-numbers and font-names available, then exits.\n");
		fprintf (stderr, "\t-N Do Not clip text that exceeds the map boundaries [Default will clip].\n");
		GMT_explain_option ('O');
		GMT_explain_option ('P');
		fprintf (stderr, "\t-S draw outline of characters.  Append pen attributes.\n");
		GMT_explain_option ('U');
		GMT_explain_option ('V');
		fprintf (stderr, "\t-W paints and outlines a rectangle underneath the text [Default is no rectangle].\n");
		fprintf (stderr, "\t   Append fill color [none].  To draw outline, append o, O, or c, and optionally a pen [No outline].\n");
		fprintf (stderr, "\t   Separate the fill color and the outline information by a comma if both are present.\n");
		fprintf (stderr, "\t   Lower case o will draw rectangular outline.\n");
		fprintf (stderr, "\t   Upper case O will draw rectangle with rounded corners (-m only).\n");
		fprintf (stderr, "\t   Lower case c will draw concave rectangle (-m only).\n");
		fprintf (stderr, "\t   Upper case C will draw convex rectangle (-m only).\n");
		GMT_explain_option ('X');
		fprintf (stderr, "\t-Z For 3-D plots: Set the z-level of map [0].  If -Z+ is given then we expect\n");
		fprintf (stderr, "\t   records to have an optional z value in the 3rd column (i.e., x y z size ...).\n");
		fprintf (stderr, "\t   Note that -Z+ also sets -N.\n");
		GMT_explain_option ('c');
		GMT_explain_option ('f');
		GMT_explain_option ('m');
		fprintf (stderr, "\t   Expects (x y size angle fontno justify linespace parwidth parjust) in segment header\n");
		fprintf (stderr, "\t   followed by lines with one or more paragraphs of text.\n");
		fprintf (stderr, "\t   parjust is one of (l)eft, (c)enter, (r)ight, or (j)ustified.\n");
		GMT_explain_option (':');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	/* Check that the options selected are mutually consistent */

	if (GMT_io.multi_segments[GMT_IN] || GMT_io.multi_segments[GMT_OUT]) do_paragraphs = TRUE;
	if (Ctrl->C.dx < 0.0 || Ctrl->C.dy < 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -C option:  clearances cannot be negative!\n", GMT_program);
		error++;
	}
	if (Ctrl->C.dx == 0.0 && Ctrl->C.dy == 0.0 && (Ctrl->W.mode == 'O')) {
		fprintf (stderr, "%s: Warning: -WO requires a nonzero -C\n", GMT_program);
		error++;
	}
	if (Ctrl->C.dx == 0.0 && Ctrl->C.dy == 0.0 && (Ctrl->W.mode == 'c' || Ctrl->W.mode == 'C')) {
		fprintf (stderr, "%s: Warning: -Wc|C requires a nonzero -C\n", GMT_program);
		error++;
	}
	if (!(fabs(Ctrl->E.azimuth) == 180.0 && Ctrl->E.elevation == 90.0) && do_paragraphs) {
		fprintf (stderr, "%s: -E not implemented with -m yet\n", GMT_program);
		error++;
	}
	if (Ctrl->D.dx == 0.0 && Ctrl->D.dy == 0.0 && Ctrl->D.line) {
		fprintf (stderr, "%s: Warning: -D<x/y>v requires one nonzero <x/y>\n", GMT_program);
		error++;
	}
	if (Ctrl->E.elevation <= 0.0 || Ctrl->E.elevation > 90.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -E option:  Elevation must be in 0-90 range\n", GMT_program);
		error++;
	}
	if (!project_info.region_supplied) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify -R option\n", GMT_program);
		error++;
	}

	if (error) exit (EXIT_FAILURE);

	load_parameters (&T, Ctrl);	/* Pass info from Ctrl to T */
	
	z_project.view_azimuth = Ctrl->E.azimuth;
	z_project.view_elevation = Ctrl->E.elevation;
	use_font = (GMT_LONG *) GMT_memory (VNULL, (size_t)GMT_N_FONTS, sizeof (GMT_LONG), GMT_program);
	add = !(T.x_offset == 0.0 && T.y_offset == 0.0);
	if (add && Ctrl->D.justify) T.boxflag |= 64;
	if (n_files > 0)
		nofile = FALSE;
	else
		n_files = 1;
	n_args = (argc > 1) ? argc : 2;

	GMT_err_fail (GMT_map_setup (west, east, south, north), "");

	GMT_plotinit (argc, argv);

	if (project_info.three_D) ps_transrotate (-z_project.xmin, -z_project.ymin, 0.0);

	if (Ctrl->Z.mode == 1) z_level_orig = project_info.z_level;
	if (Ctrl->Z.active) {
		project_info.z_level = Ctrl->Z.level;
		GMT_z_to_zz (Ctrl->Z.level, &z_level);
	}
	GMT_setpen (&T.txtpen);

	if (!(Ctrl->N.active || Ctrl->Z.mode == 1)) GMT_map_clip_on (GMT_no_rgb, 3);

	last_font = 0;
	ix = (gmtdefs.xy_toggle[GMT_IN]);	iy = 1 - ix;

	/* if (draw_box && project_info.three_D) draw_box = FALSE; */	/* Not implemented yet */

	/* Mark used fonts */

	if (strcmp (gmtdefs.encoding.name, "Standard") != 0) {
		memset (use_font, 0, GMT_N_FONTS * sizeof (int));
		if (GMT_ps.unix_time) use_font[0] = use_font[1] = TRUE;
		use_font[gmtdefs.annot_font[0]] = TRUE;
		use_font[gmtdefs.annot_font[1]] = TRUE;
		if (frame_info.header[0]) use_font[gmtdefs.header_font] = TRUE;
		if (frame_info.axis[0].label[0] || frame_info.axis[1].label[0] || frame_info.axis[2].label[0]) use_font[gmtdefs.label_font] = TRUE;
	}

	if (do_paragraphs) {
		n_alloc = GMT_CHUNK;
		word = (char **) GMT_memory (VNULL, n_alloc, sizeof (char *), GMT_program);
	}

	old_GMT_world_map = GMT_world_map;
	GMT_world_map = TRUE;

	for (fno = 1; !done && fno < n_args; fno++) {	/* Loop over all input files */
		if (!nofile && argv[fno][0] == '-') continue;
		if (nofile) {
			fp = GMT_stdin;
			done = TRUE;
		}
		else if ((fp = GMT_fopen (argv[fno], "r")) == NULL) {
			fprintf (stderr, "%s: Cannot open file %s\n", GMT_program, argv[fno]);
			continue;
		}

		if (!nofile && gmtdefs.verbose) fprintf (stderr, "%s: Working on file %s\n", GMT_program, argv[fno]);
		if (GMT_io.io_header[GMT_IN]) for (i = 0; i < GMT_io.n_header_recs; i++) GMT_fgets (line, BUFSIZ, fp);

		while (GMT_fgets (line, BUFSIZ, fp)) {

			GMT_enforce_rgb_triplets (line, BUFSIZ);	/* If @; is used, make sure the color information passed on to ps_text is in r/b/g format */
			
			if (do_paragraphs) {	/* Paragraph mode */
				if (line[0] == GMT_io.EOF_flag[GMT_IN]) {

					skip_text_records = FALSE;
					if (n_processed) {	/* Must output what we got */
						GMT_putwords (plot_x, plot_y, word, n_words, &T);
						n_processed = n_words = 0;
						n_paragraphs++;
					}
					
					if (Ctrl->Z.mode == 1) {	/* Expect z in 3rd column */
						nscan = sscanf (&line[2], "%s %s %s %lf %lf %s %s %s %s %c\n", txt_x, txt_y, txt_z, &T.font_size, &T.paragraph_angle, this_font, just_key, txt_a, txt_b, &pjust);
						if ((GMT_scanf (txt_z, GMT_io.in_col_type[2], &Ctrl->Z.level) == GMT_IS_NAN)) {
							fprintf (stderr, "%s: Record %ld had bad z coordinate, must exit)\n", GMT_program, n_read);
							exit (EXIT_FAILURE);
						}
						project_info.z_level = Ctrl->Z.level;
						GMT_z_to_zz (Ctrl->Z.level, &z_level);
					}
					else
						nscan = sscanf (&line[2], "%s %s %lf %lf %s %s %s %s %c\n", txt_x, txt_y, &T.font_size, &T.paragraph_angle, this_font, just_key, txt_a, txt_b, &pjust);
					if ((GMT_scanf (txt_x, GMT_io.in_col_type[0], &xy[ix]) == GMT_IS_NAN) || (GMT_scanf (txt_y, GMT_io.in_col_type[1], &xy[iy]) == GMT_IS_NAN)) {
						fprintf (stderr, "%s: Record %ld had bad x and/or y coordinates, must exit)\n", GMT_program, n_read);
						exit (EXIT_FAILURE);
					}

					if (nscan != (9 + Ctrl->Z.mode)) {
						fprintf (stderr, "%s: Record %ld had incomplete paragraph information, must exit)\n", GMT_program, n_read);
						exit (EXIT_FAILURE);
					}
					GMT_geo_to_xy (xy[0], xy[1], &plot_x, &plot_y);
					if (!Ctrl->N.active) {
						skip_text_records = TRUE;	/* If this record should be skipped we must skip the whole paragraph */
						GMT_map_outside (xy[0], xy[1]);
						if (GMT_abs (GMT_x_status_new) > 1 || GMT_abs (GMT_y_status_new) > 1) continue;
						skip_text_records = FALSE;	/* Since we got here we do not want to skip */
					}
					T.text_justify = (int)pjust;
					T.block_justify = (isdigit ((int)just_key[0])) ? atoi (just_key) : GMT_just_decode (just_key, 12);
					if (T.block_justify == -99) {
						fprintf (stderr, "%s: Record %ld had bad justification info (set to LB)\n", GMT_program, n_read);
						T.block_justify = 1;
					}
					if (Ctrl->A.active) {
						GMT_azim_to_angle (xy[0], xy[1], 0.1, T.paragraph_angle, &tmp);
						T.paragraph_angle = fmod (tmp + 360.0 + 90.0, 180.0) - 90.0;	/* Ensure usable angles for text plotting */
						if (fabs (T.paragraph_angle - tmp) > 179.0) T.block_justify = 4 * (T.block_justify/4) + 2 - (T.block_justify%4 - 2);	/* Flip any L/R code */
					}
					T.paragraph_font = GMT_font_lookup (this_font, GMT_font, GMT_N_FONTS);
					if (T.paragraph_font == GMT_N_FONTS) {
						fprintf (stderr, "%s: Record %ld had bad font (set to %s (0))\n", GMT_program, n_read, GMT_font[0].name);
						T.paragraph_font = 0;
					}
					if (!use_font[T.paragraph_font] && strcmp (gmtdefs.encoding.name, "Standard") != 0) {	/* Must reencode this font */
						ps_encode_font (T.paragraph_font);
						use_font[T.paragraph_font] = TRUE;
					}
					T.line_spacing = GMT_convert_units (txt_a, GMT_INCH);
					T.paragraph_width  = GMT_convert_units (txt_b, GMT_INCH);
					master_record = TRUE;
				}
				else {	/* Text block record */
					if (skip_text_records) continue;	/* Skip all records for this paragraph */
					if (!master_record) {
						fprintf (stderr, "%s: Text record line %ld not preceded by paragraph information, must exit)\n", GMT_program, n_read);
						exit (EXIT_FAILURE);
					}

					line[strlen(line)-1] = 0;

					if (line[0] == 0) {	/* Blank line marked by single NULL character */
						word[n_words] = (char *) calloc ((size_t)(1), sizeof (char));
						n_words++;
						if (n_words == n_alloc) {
							n_alloc <<= 1;
							word = (char **) GMT_memory (VNULL, n_alloc, sizeof (char *), GMT_program);
						}
					}
					pos = 0;
					while ((GMT_strtok (line, " ", &pos, p))) {
						word[n_words] = strdup (p);
						n_words++;
						if (n_words == n_alloc) {
							n_alloc <<= 1;
							word = (char **) GMT_memory (VNULL, n_alloc, sizeof (char *), GMT_program);
						}
					}

					n_processed++;
				}
				n_read++;
			}
			else {	/* Old-style pstext input */
				if (GMT_is_a_blank_line (line)) continue;	/* Skip blank lines or # comments */

				if (Ctrl->Z.mode == 1) {	/* Expect z in 3rd column */
					nscan = sscanf (line, "%s %s %s %lf %lf %s %s %[^\n]\n", txt_x, txt_y, txt_z, &T.font_size, &T.paragraph_angle, this_font, just_key, text);
					if ((GMT_scanf (txt_z, GMT_io.in_col_type[2], &Ctrl->Z.level) == GMT_IS_NAN)) {
						fprintf (stderr, "%s: Record %ld had bad z coordinate, must exit)\n", GMT_program, n_read);
						exit (EXIT_FAILURE);
					}
					project_info.z_level = Ctrl->Z.level;
					GMT_z_to_zz (Ctrl->Z.level, &z_level);
				}
				else
					nscan = sscanf (line, "%s %s %lf %lf %s %s %[^\n]\n", txt_x, txt_y, &T.font_size, &T.paragraph_angle, this_font, just_key, text);
				if ((GMT_scanf (txt_x, GMT_io.in_col_type[0], &xy[ix]) == GMT_IS_NAN) || (GMT_scanf (txt_y, GMT_io.in_col_type[1], &xy[iy]) == GMT_IS_NAN)) {
					fprintf (stderr, "%s: Record %ld had bad x and/or y coordinates, must exit)\n", GMT_program, n_read);
					exit (EXIT_FAILURE);
				}
				if (nscan != (7+Ctrl->Z.mode)) {
					fprintf (stderr, "%s: Record %ld is incomplete (skipped)\n", GMT_program, n_read);
					continue;
				}
				n_read++;
				GMT_geo_to_xy (xy[0], xy[1], &plot_x, &plot_y);
				xx[0] = plot_x;	yy[0] = plot_y;
				if (!Ctrl->N.active) {
					GMT_map_outside (xy[0], xy[1]);
					if (GMT_abs (GMT_x_status_new) > 1 || GMT_abs (GMT_y_status_new) > 1) continue;
				}

				T.block_justify = (isdigit ((int)just_key[0])) ? atoi (just_key) : GMT_just_decode (just_key, 12);
				if (T.block_justify == -99) {
					fprintf (stderr, "%s: Record %ld had bad justification info (skipped)\n", GMT_program, n_read);
					continue;
				}
				if (Ctrl->A.active) {
					GMT_azim_to_angle (xy[0], xy[1], 0.1, T.paragraph_angle, &tmp);
					T.paragraph_angle = fmod (tmp + 360.0 + 90.0, 180.0) - 90.0;	/* Ensure usable angles for text plotting */
					if (fabs (T.paragraph_angle - tmp) > 179.0) T.block_justify = 4 * (T.block_justify/4) + 2 - (T.block_justify%4 - 2);	/* Flip any L/R code */
				}
				if (add) {
					if (Ctrl->D.justify) {	/* Smart offset according to justification (from Dave Huang) */
						GMT_smart_justify (T.block_justify, T.paragraph_angle, T.x_offset, T.y_offset, &plot_x, &plot_y);
					} else {	/* Default hard offset */
						plot_x += T.x_offset;
						plot_y += T.y_offset;
					}
					xx[1] = plot_x;	yy[1] = plot_y;
				}
				T.paragraph_font = GMT_font_lookup (this_font, GMT_font, GMT_N_FONTS);
				if (T.paragraph_font == GMT_N_FONTS) {
					fprintf (stderr, "%s: Record %ld had bad font (set to %s (0))\n", GMT_program, n_read, GMT_font[0].name);
					T.paragraph_font = 0;
				}
				n_paragraphs++;

				if (T.paragraph_font != last_font) {
					ps_setfont (T.paragraph_font);
					last_font = T.paragraph_font;
				}
				if (!use_font[T.paragraph_font] && strcmp (gmtdefs.encoding.name, "Standard") != 0) {	/* Must reencode this font */
					ps_encode_font (T.paragraph_font);
					use_font[T.paragraph_font] = TRUE;
				}
				if (T.boxflag & 32) {	/* Draw line from original point to shifted location */
					GMT_setpen (&T.vecpen);
					if (project_info.three_D) GMT_2D_to_3D (xx, yy, project_info.z_level, (GMT_LONG)2);
					ps_segment (xx[0], yy[0], xx[1], yy[1]);
				}
				if (T.boxflag & 31) {
					GMT_setpen (&T.boxpen);
					if (T.space_flag) {	/* Meant % of fontsize */
						dx = 0.01 * T.x_space * T.font_size / 72.0;
						dy = 0.01 * T.y_space * T.font_size / 72.0;
					}
					else {
						dx = T.x_space;
						dy = T.y_space;
					}
					GMT_textbox3D (plot_x, plot_y, z_level, T.font_size, T.paragraph_font, text, T.paragraph_angle, T.block_justify, (T.boxflag & 5), dx, dy, T.boxfill.rgb);
					GMT_setpen (&T.txtpen);
				}
				if (Ctrl->S.active) {
					GMT_setpen (&T.txtpen);
					GMT_text3D (plot_x, plot_y, z_level, T.font_size, T.paragraph_font, text, T.paragraph_angle, T.block_justify, TRUE);
				}
				if (Ctrl->G.active) {
					ps_setpaint (T.txtfill.rgb);
					GMT_text3D (plot_x, plot_y, z_level, T.font_size, T.paragraph_font, text, T.paragraph_angle, T.block_justify, FALSE);
				}
			}
		} /* Go read next line */

		if (do_paragraphs && n_processed) {	/* Must output what we got */
			GMT_putwords (plot_x, plot_y, word, n_words, &T);
			n_processed = n_words = 0;
			n_paragraphs++;
		}
		if (fp != GMT_stdin) GMT_fclose (fp);
	}

 	if (do_paragraphs) GMT_free ((void *)word);

	if (!(Ctrl->N.active || Ctrl->Z.mode == 1)) GMT_map_clip_off (); 

	GMT_world_map = old_GMT_world_map;

	if (Ctrl->Z.mode == 1) project_info.z_level = z_level_orig;

	GMT_map_basemap ();

	if (project_info.three_D) ps_rotatetrans (z_project.xmin, z_project.ymin, 0.0);

	GMT_plotend ();

	if (gmtdefs.verbose) {
		if (do_paragraphs)
			fprintf (stderr, "pswords: Plotted %ld textblocks\n", n_paragraphs);
		else
			fprintf (stderr, "pstext: Plotted %ld text-strings\n", n_paragraphs);
	}

	GMT_free ((void *)use_font);

	Free_pstext_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void GMT_putwords (double x, double y, char **text, GMT_LONG n_words, struct PSTEXT_INFO *T) {
	char *boxpen_texture = CNULL, *vecpen_texture = CNULL;
	GMT_LONG boxpen_width, boxpen_offset, vecpen_width, vecpen_offset, i;
	int boxpen_rgb[3], vecpen_rgb[3];
	double dx, dy;

	boxpen_texture = GMT_convertpen (&T->boxpen, &boxpen_width, &boxpen_offset, boxpen_rgb);
	vecpen_texture = GMT_convertpen (&T->vecpen, &vecpen_width, &vecpen_offset, vecpen_rgb);

	if (T->space_flag) {	/* Meant % of fontsize */
		dx = 0.01 * T->x_space * T->font_size / 72.0;
		dy = 0.01 * T->y_space * T->font_size / 72.0;
	}
	else {
		dx = T->x_space;
		dy = T->y_space;
	}
	ps_words (x, y, text, n_words, T->line_spacing, T->paragraph_width, T->text_justify, T->paragraph_font, T->font_size,
		T->paragraph_angle, T->txtfill.rgb, T->block_justify, T->boxflag, T->x_offset, T->y_offset, dx, dy,
		boxpen_width, boxpen_texture, boxpen_offset, boxpen_rgb,
		vecpen_width, vecpen_texture, vecpen_offset, vecpen_rgb,
		T->boxfill.rgb);

	if (boxpen_texture) GMT_free ((void *)boxpen_texture);
	if (vecpen_texture) GMT_free ((void *)vecpen_texture);

	for (i = 0; i < n_words; i++) free ((void *)text[i]);
}

void load_parameters (struct PSTEXT_INFO *T, struct PSTEXT_CTRL *C) {
	T->x_space = C->C.dx;
	T->y_space = C->C.dy;
	T->space_flag = (C->C.percent) ? 1 : 0;
	if (C->D.active) {
		T->x_offset = C->D.dx;
		T->y_offset = C->D.dy;
		if (C->D.line) T->boxflag |= 32;
		T->vecpen = C->D.pen;
	}
	if (C->G.active) {
		T->txtfill = C->G.fill;
	}
	if (C->S.active) {
		T->txtpen = C->S.pen;
	}
	if (C->W.active) {
		if (C->W.outline) T->boxflag |= 1;	/* Want box outline */
		if (C->W.mode == 'O') T->boxflag |= 4;	/* Want rounded box outline */
		if (C->W.mode == 'c') T->boxflag |= 8;	/* Want concave box outline */
		if (C->W.mode == 'C') T->boxflag |= 16;	/* Want convex box outline */
		T->boxpen = C->W.pen;
		T->boxfill = C->W.fill;
		if (C->W.paint) T->boxflag |= 2;
		if (!(T->boxflag | 31)) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR -W option: Must give arguments\n", GMT_program);
			exit (EXIT_FAILURE);
		}
	}

}

void *New_pstext_Ctrl () {	/* Allocate and initialize a new control structure */
	struct PSTEXT_CTRL *C;
	
	C = (struct PSTEXT_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct PSTEXT_CTRL), "New_pstext_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
		
	GMT_init_pen (&C->D.pen, GMT_PENWIDTH);
	C->C.dx = C->C.dy = 15.0;	/* 15% of font size is default clearance */
	C->C.percent = TRUE;
	C->E.azimuth = 180.0;
	C->E.elevation = 90.0;
	C->G.active = TRUE;
	GMT_init_fill (&C->G.fill, gmtdefs.basemap_frame_rgb[0], gmtdefs.basemap_frame_rgb[1], gmtdefs.basemap_frame_rgb[2]);
	GMT_init_fill (&C->W.fill, -1, -1, -1);	/* No fill */
	GMT_init_pen (&C->S.pen, GMT_PENWIDTH);
	GMT_init_pen (&C->W.pen, GMT_PENWIDTH);

	return ((void *)C);
}

void Free_pstext_Ctrl (struct PSTEXT_CTRL *C) {	/* Deallocate control structure */
	GMT_free ((void *)C);	
}
