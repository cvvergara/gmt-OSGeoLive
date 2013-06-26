/*--------------------------------------------------------------------
 *	$Id: pslegend.c 9923 2012-12-18 20:45:53Z pwessel $
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
 * pslegend will make map legends from input that specifies what will go
 * into the legend, such as headers, symbols with explanatory text,
 * paragraph text, and empty space and horizontal/vertical lines.
 *
 * Author:	Paul Wessel
 * Date:	18-SEP-2001
 * Version:	4
 *
 */

#include "gmt.h"
#include "pslib.h"

#ifdef WIN32
#include <process.h>
#endif

struct PSLEGEND_CTRL {
	struct C {	/* -C<dx>/<dy> */
		GMT_LONG active;
		double dx, dy;
	} C;
	struct D {	/* -D[x]<x0>/<y0>/w/h/just */
		GMT_LONG active;
		GMT_LONG cartesian;
		double lon, lat, width, height;
		char justify[3];
	} D;
	struct F {	/* -F */
		GMT_LONG active;
	} F;
	struct G {	/* -G<fill> */
		GMT_LONG active;
		struct GMT_FILL fill;
	} G;
	struct L {	/* -L<spacing> */
		GMT_LONG active;
		double spacing;
	} L;
	struct S {	/* -C<script> */
		GMT_LONG active;
		char *file;
	} S;
};

int main (int argc, char **argv)
{
	GMT_LONG flush_paragraph = FALSE, draw_vertical_line = FALSE, gave_label, gave_mapscale_options;
	
#ifdef WIN32
	char *del = "del", *escape = "^", quote = ' ';
#else
	char *del = "rm -f", *escape = "\\", quote = '\'';
#endif
	char txt_a[GMT_LONG_TEXT], txt_b[GMT_LONG_TEXT], txt_c[GMT_LONG_TEXT], txt_d[GMT_LONG_TEXT], txt_e[GMT_LONG_TEXT];
	char txt_f[GMT_LONG_TEXT], key[GMT_LONG_TEXT], sub[GMT_LONG_TEXT];
	char symbol[GMT_LONG_TEXT], text[BUFSIZ], image[BUFSIZ], xx[GMT_LONG_TEXT], yy[GMT_LONG_TEXT];
	char size[GMT_LONG_TEXT], angle[GMT_LONG_TEXT], mapscale[GMT_LONG_TEXT], font[GMT_LONG_TEXT], lspace[GMT_LONG_TEXT];
	char tw[GMT_LONG_TEXT], jj[GMT_LONG_TEXT], line[BUFSIZ], vpen[GMT_LONG_TEXT], script[GMT_LONG_TEXT], tmptxt[GMT_LONG_TEXT];
	char sarg[GMT_LONG_TEXT], sparg[BUFSIZ], txtcolor[GMT_LONG_TEXT], psxy[GMT_LONG_TEXT], pstext[GMT_LONG_TEXT];
	char bar_cpt[GMT_LONG_TEXT], bar_gap[GMT_LONG_TEXT], bar_height[GMT_LONG_TEXT], bar_opts[BUFSIZ], mode, just, *opt = NULL;
	char *barg = CNULL, *jarg = CNULL, *rarg = CNULL, *f = CNULL, *u = CNULL;
	
	GMT_LONG i, k, n, justify = 0, n_columns = 1, error = 0, column_number = 0, n_files = 0, n_scan, ifont;
	
	double x_off, west, east, south, north, x, y, x0, y0, L, off_ss, off_tt, V = 0.0;
	double half_line_spacing, quarter_line_spacing, one_line_spacing, font_size, y_start = 0.0, d_off;
	
	struct imageinfo header;
	struct PSLEGEND_CTRL *Ctrl = NULL;
	
	FILE *fp = NULL, *fpo = NULL;
	
	void *New_pslegend_Ctrl (), Free_pslegend_Ctrl (struct PSLEGEND_CTRL *C);

	/* Define the fraction of the height of the font to the font size */
#define FONT_HEIGHT1 (GMT_font[gmtdefs.annot_font[0]].height)
#define FONT_HEIGHT2 (GMT_font[ifont].height)
#define FONT_HEIGHT3 (GMT_font[gmtdefs.label_font].height)
	/* This used to be all 0.75 */
	
	/* Because pslegend uses system calls we must first make a copy of any arguments that GMT_begin will remove, such
	 * as +gmtdefaults and --PAR=value
	 */
	 
	memset ((void *)sparg, 0, (size_t)BUFSIZ);
	for (i = 1, k = 0; i < argc; i++) {
		if (argv[i][0] == '+' || (argv[i][0] == '-' && argv[i][1] == '-')) {
			if (k) strcat (sparg, " ");
			strcat (sparg, argv[i]);
			k++;
		}
	}

	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct PSLEGEND_CTRL *)New_pslegend_Ctrl ();	/* Allocate and initialize a new control structure */

	/* Check and interpret the command line arguments */

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch(argv[i][1]) {

				/* Common parameters */

				case 'J':
				case 'K':
				case 'O':
				case 'P':
				case 'R':
				case 'V':
				case 'X':
				case 'x':
				case 'Y':
				case 'y':
				case 'c':
				case '\0':
					if (argv[i][1] == 'R') rarg = argv[i];	/* save original -R option */
					if (argv[i][1] == 'J') jarg = argv[i];	/* save original -J option */
					error += GMT_parse_common_options (argv[i], &west, &east, &south, &north);
					break;

				case 'U':	/* Just need a pointer to pass along */
					u = &argv[i][2];
					break;

				/* Supplemental parameters */

				case 'B':
					barg = argv[i];	/* Keep this for later */
					break;
				case 'C':	/* Sets the clearance between frame and internal items */
					Ctrl->C.active = TRUE;
					sscanf (&argv[i][2], "%[^/]/%s", txt_a, txt_b);
					Ctrl->C.dx = GMT_convert_units (txt_a, GMT_INCH);
					Ctrl->C.dy = GMT_convert_units (txt_b, GMT_INCH);
					break;
				case 'D':	/* Sets position and size of legend */
					Ctrl->D.active = TRUE;
					if (argv[i][2] == 'x') {	/* Gave location directly in projected units (inches, cm, etc) */
						Ctrl->D.cartesian = TRUE;
						k = 3;
					}
					else				/* Gave lon, lat */
						k = 2;
					n = sscanf (&argv[i][k], "%[^/]/%[^/]/%[^/]/%[^/]/%s", txt_a, txt_b, txt_c, txt_d, Ctrl->D.justify);
					if (n != 5) {
						fprintf (stderr, "%s ERROR: Syntax is -D[x]<xpos>/<ypos>/<width>/<height>/<justify>\n", GMT_program);
						error++;
					}
					if (argv[i][2] == 'x') {
						Ctrl->D.lon = GMT_convert_units (txt_a, GMT_INCH);
						Ctrl->D.lat = GMT_convert_units (txt_b, GMT_INCH);
					}
					else {	/* Given in user units, likely degrees */
						error += GMT_verify_expectations (GMT_io.in_col_type[0], GMT_scanf (txt_a, GMT_io.in_col_type[0], &Ctrl->D.lon), txt_a);
						error += GMT_verify_expectations (GMT_io.in_col_type[1], GMT_scanf (txt_b, GMT_io.in_col_type[1], &Ctrl->D.lat), txt_b);
					}
					Ctrl->D.width   = GMT_convert_units (txt_c, GMT_INCH);
					Ctrl->D.height  = GMT_convert_units (txt_d, GMT_INCH);
					break;
				case 'F':
					Ctrl->F.active = TRUE;
					break;
				case 'G':	/* Inside legend fill? */
					Ctrl->G.active = TRUE;
					if (GMT_getfill (&argv[i][2], &Ctrl->G.fill)) {	/* We check syntax here */
						GMT_fill_syntax ('G', " ");
						error++;
					}
					f = argv[i];		/* Pointer to fill argument */
					break;
				case 'L':			/* Sets linespacing in units of fontsize [1.1] */
					Ctrl->L.active = TRUE;
					Ctrl->L.spacing = atof (&argv[i][2]);
					break;
				case 'S':
					Ctrl->S.active = TRUE;
					if (argv[i][2]) {	/* Specified script filename */
						Ctrl->S.file = strdup (&argv[i][2]);
					}
					break;

				/* Options not recognized */

				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else if (n_files == 0) {
			if ((fp = GMT_fopen (argv[i], "r")) == NULL) {
				fprintf (stderr, "%s ERROR: Cannot open file %s\n", GMT_program, argv[i]);
				exit (EXIT_FAILURE);
			}
			n_files++;
		}
		else {
			fprintf (stderr, "%s ERROR: Only one file argument allowed\n", GMT_program);
			exit (EXIT_FAILURE);
		}
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "pslegend %s - To plot legends on maps\n\n", GMT_VERSION);
		fprintf (stderr, "usage: pslegend [<infofile>] -D[x]<x0>/<y0>/w/h/just %s %s\n", GMT_J_OPT, GMT_Rgeo_OPT);
		fprintf (stderr, "\t[%s] [-C<dx>/<dy>] [-F] [-G<fill>] [-K] [-L<spacing>] [-O] [-P] [-S[<script>]]\n", GMT_B_OPT);
		fprintf (stderr, "\t[%s] [-V] [%s] [%s] [%s]\n\n", GMT_U_OPT, GMT_X_OPT, GMT_Y_OPT, GMT_c_OPT);
		fprintf (stderr, "\tReads legend layout information from <infofile> [or stdin].\n");
		fprintf (stderr, "\t(See manual page for more information).\n");
		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\t-D sets position and size of legend box.  Prepend x if coordinates are projected.\n");
		fprintf (stderr, "\t   Append the justification of the whole legend box using pstext justification codes.\n");
		GMT_explain_option ('j');
		GMT_explain_option ('R');
		fprintf (stderr, "\n\tOPTIONS:\n");
		GMT_explain_option ('b');
		fprintf (stderr, "\t-C sets the clearance between legend frame and internal items [0.05i/0.05i].\n");
		fprintf (stderr, "\t-F Draw border around the legend (using FRAME_PEN) [Default is no border].\n");
		GMT_fill_syntax ('G', "Set the fill for the legend box [Default is no fill].");
		GMT_explain_option ('K');
		fprintf (stderr, "\t-L Sets the linespacing factor in units of the current annotation font size [1.1].\n");
		GMT_explain_option ('O');
		GMT_explain_option ('P');
		fprintf (stderr, "\t-S Dump legend script to stdout, or optionally to file <script>.\n");
		fprintf (stderr, "\t   [Default is to write PostScript output].\n");
		GMT_explain_option ('U');
		GMT_explain_option ('V');
		GMT_explain_option ('X');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	/* Check that the options selected are mutually consistent */

	if (Ctrl->C.dx < 0.0 || Ctrl->C.dy < 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -C option:  clearances cannot be negative!\n", GMT_program);
		error++;
	}
	if (Ctrl->D.width < 0.0 || Ctrl->D.height < 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -D option:  legend box sizes cannot be negative!\n", GMT_program);
		error++;
	}
	if (!project_info.region_supplied) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify -R option\n", GMT_program);
		error++;
	}

	if (error) exit (EXIT_FAILURE);

	if (Ctrl->S.active) {
		if (Ctrl->S.file) {
			if ((fpo = fopen (Ctrl->S.file, "w")) == NULL) {
				fprintf (stderr, "%s: GMT SYNTAX ERROR -S option:  Cannot create file %s\n", GMT_program, Ctrl->S.file);
				exit (EXIT_FAILURE);
			}
			if (gmtdefs.verbose) fprintf (stderr, "%s: Generate legend script %s\n", GMT_program, Ctrl->S.file);
		}
		else {
			fpo = GMT_stdout;
			if (gmtdefs.verbose) fprintf (stderr, "%s: Generate legend script [stdout]\n", GMT_program);
		}
	}
	else {
		if (GMT_TMPDIR)
			sprintf (script, "%s/pslegend_%ld.bat", GMT_TMPDIR, (GMT_LONG)getpid ());
		else
			sprintf (script, "GMT%ld.bat", (GMT_LONG)getpid ());
		if (gmtdefs.verbose) fprintf (stderr, "%s: Generate temporary legend script %s\n", GMT_program, script);
		fpo = fopen (script, "w");
	}
#ifdef WIN32
	fprintf (fpo, "@ECHO OFF\n");
#endif
	if (!fp) fp = GMT_stdin;

	GMT_err_fail (GMT_map_setup (west, east, south, north), "");

	justify = GMT_just_decode (Ctrl->D.justify, 12);
	
	if (!Ctrl->D.cartesian) {
		GMT_geo_to_xy (Ctrl->D.lon, Ctrl->D.lat, &x, &y);
		Ctrl->D.lon = x;	Ctrl->D.lat = y;
	}

	/* Adjust for -X -Y shifts */

	Ctrl->D.lon += GMT_ps.x_origin;
	Ctrl->D.lat += GMT_ps.y_origin;

	/* Allow for justification */

	Ctrl->D.lon -= 0.5 * ((justify-1)%4) * Ctrl->D.width;
	Ctrl->D.lat -= 0.5 * (justify/4) * Ctrl->D.height;

	/* Initialize psxy and pstext call strings and tmptxt file name with appropriate values */
	
	if (GMT_ps.absolute) {	/* Must pass the -Xa -Ya settings to every psxy and pstext call */
		mode = 'a';
		sprintf (pstext, "pstext -R -JX -O -K -X%c%gi -Y%c%gi %s", mode, Ctrl->D.lon, mode, Ctrl->D.lat, sparg);
		sprintf (psxy, "psxy -R -JX -O -K -X%c%gi -Y%c%gi %s", mode, Ctrl->D.lon, mode, Ctrl->D.lat, sparg);
	}
	else {			/* No need to pass, -X0 -Y0 are set implicitly */
		mode = 'r';
		sprintf (pstext, "pstext -R -JX -O -K %s", sparg);
		sprintf (psxy, "psxy -R -JX -O -K %s", sparg);
	}
	if (GMT_TMPDIR)
		sprintf (tmptxt, "%s/pslegend_%ld.txt", GMT_TMPDIR, (GMT_LONG)getpid ());
	else
		sprintf (tmptxt, "pslegend_%ld.txt", (GMT_LONG)getpid ());
	if (gmtdefs.verbose) fprintf (stderr, "%s: Use temporary input file %s\n", GMT_program, tmptxt);
	
	/* First draw legend frame box.  Note -JX%gi/-%gi which means y is positive down from the top of the box */

	if (Ctrl->F.active || Ctrl->G.active)
		fprintf (fpo, "psbasemap -R0/%g/0/%g -JX%gi/-%gi -X%c%gi -Y%c%gi -K %s", Ctrl->D.width, Ctrl->D.height, Ctrl->D.width, Ctrl->D.height, mode, Ctrl->D.lon, mode, Ctrl->D.lat, sparg);
	else
		fprintf (fpo, "psxy -R0/%g/0/%g -JX%gi/-%gi -T -X%c%gi -Y%c%gi -K %s", Ctrl->D.width, Ctrl->D.height, Ctrl->D.width, Ctrl->D.height, mode, Ctrl->D.lon, mode, Ctrl->D.lat, sparg);
	if (GMT_ps.portrait) fprintf (fpo, " -P");
	if (Ctrl->F.active) fprintf (fpo, " -B0");
	if (GMT_ps.overlay) fprintf (fpo, " -O");
	if (Ctrl->G.active) fprintf (fpo, " %s", f);
	if (u) fprintf (fpo, " -U\"%s\"", u);
	fprintf (fpo, "\n");

	x0 = Ctrl->C.dx;
	y0 = Ctrl->C.dy;	/* y0 follows the bottom of the last printed text box */
	one_line_spacing = Ctrl->L.spacing * gmtdefs.annot_font_size[0] / 72.0;
	half_line_spacing = 0.5 * one_line_spacing;
	quarter_line_spacing = 0.25 * one_line_spacing;
	column_number = 0;

	while (GMT_fgets (line, BUFSIZ, fp)) {

		if (GMT_is_a_blank_line (line)) continue;	/* Skip blank lines or # comments */
		if (line[0] != 'T' && flush_paragraph) {	/* Flush contents of pending paragraph */
			if (flush_paragraph) fprintf (fpo, "%s -m %s\n%s %s\n", pstext, tmptxt, del, tmptxt);
			flush_paragraph = FALSE;
			column_number = 0;
		}

		switch (line[0]) {
			case 'C':	/* Color change */
				sscanf (&line[2], "%[^\n]", txtcolor);
				if (GMT_ps.absolute) {	/* Must pass the -Xa -Ya settings to every psxy and pstext call */
					if (txtcolor[0] == '-') sprintf (pstext, "pstext -R -JX -O -K -X%c%gi -Y%c%gi %s", mode, Ctrl->D.lon, mode, Ctrl->D.lat, sparg);
					else sprintf (pstext, "pstext -R -JX -O -K -G%s -X%c%gi -Y%c%gi %s", txtcolor, mode, Ctrl->D.lon, mode, Ctrl->D.lat, sparg);
				}
				else {
					if (txtcolor[0] == '-') sprintf (pstext, "pstext -R -JX -O -K %s", sparg);
					else sprintf (pstext, "pstext -R -JX -O -K -G%s %s", txtcolor, sparg);
				}
				break;
			case 'B':	/* Color scale Bar */
				bar_opts[0] = '\0';
				sscanf (&line[2], "%s %s %s %[^\n]", bar_cpt, bar_gap, bar_height, bar_opts);
				x_off = GMT_convert_units (bar_gap, GMT_INCH);
				fprintf (fpo, "psscale -C%s -O -K -D%gi/%gi/%gi/%sh %s\n", bar_cpt, 0.5 * Ctrl->D.width, Ctrl->D.height-y0, Ctrl->D.width - 2 * x_off, bar_height, bar_opts);
				y0 += GMT_convert_units (bar_height, GMT_INCH) + gmtdefs.tick_length + gmtdefs.annot_offset[0] + FONT_HEIGHT1 * gmtdefs.annot_font_size[0] / 72.0;
				column_number = 0;
				break;
			case 'H':	/* Header record */
				sscanf (&line[2], "%s %s %[^\n]", size, font, text);
				if (size[0] == '-') sprintf (size, "%g", gmtdefs.header_font_size);
				if (font[0] == '-') sprintf (font, "%ld", gmtdefs.header_font);
				ifont = GMT_font_lookup (font, GMT_font, GMT_N_FONTS);
				font_size = atoi(size);
				d_off = 0.5 * (Ctrl->L.spacing - FONT_HEIGHT2) * font_size / 72.0;	/* To center the text */
				y0 += Ctrl->L.spacing * font_size / 72.0;
				fprintf (fpo, "echo %c%g %g %s 0 %s BC %s%c | %s\n", quote, 0.5 * Ctrl->D.width, y0 - d_off, size, font, text, quote, pstext);
				column_number = 0;
				break;

			case 'I':	/* Image record */
				sscanf (&line[2], "%s %s %s", image, size, key);
				(void) ps_load_image (image, &header);
				justify = GMT_just_decode (key, 12);
				x_off = (justify%4 == 1) ? x0 : ((justify%4 == 3) ? Ctrl->D.width - Ctrl->C.dx : 0.5 * Ctrl->D.width);
				fprintf (fpo, "psimage -O -K %s -W%s -C%gi/%gi/%s\n", image, size, x_off, Ctrl->D.height-y0, key);
				y0 += GMT_convert_units (size, GMT_INCH) * (double)header.height / (double)header.width;
				column_number = 0;
				break;

			case 'L':	/* Label record */
				sscanf (&line[2], "%s %s %s %[^\n]", size, font, key, text);
				if (size[0] == '-') sprintf (size, "%g", gmtdefs.label_font_size);
				if (font[0] == '-') sprintf (font, "%ld", gmtdefs.label_font);
				ifont = GMT_font_lookup (font, GMT_font, GMT_N_FONTS);
				font_size = atoi(size);
				d_off = 0.5 * (Ctrl->L.spacing - FONT_HEIGHT2) * font_size / 72.0;	/* To center the text */
				if (column_number%n_columns == 0) y0 += Ctrl->L.spacing * font_size / 72.0;
				justify = GMT_just_decode (key, 0);
				x_off = (Ctrl->D.width / n_columns) * (column_number%n_columns);
				x_off += (justify%4 == 1) ? x0 : ((justify%4 == 3) ? (Ctrl->D.width / n_columns) - Ctrl->C.dx : 0.5 * (Ctrl->D.width / n_columns));
				fprintf (fpo, "echo %c%g %g %s 0 %s B%s %s%c | %s\n", quote, x_off, y0 - d_off, size, font, key, text, quote, pstext);
				column_number++;
				break;

			case 'M':	/* Map scale record M lon0|- lat0 length[n|m|k][+opts] f|p  [-R -J] */
				n_scan = sscanf (&line[2], "%s %s %s %s %s %s", txt_a, txt_b, txt_c, txt_d, txt_e, txt_f);
				k = (txt_d[0] != 'f') ? 1 : 0;	/* Determines if we start -L with f or not */
				for (i = 0, gave_mapscale_options = FALSE; txt_c[i] && !gave_mapscale_options; i++) if (txt_c[i] == '+') gave_mapscale_options = TRUE;
				/* Default assumes label is added on top */
				just = 't';
				gave_label = TRUE;
				d_off = FONT_HEIGHT3 * gmtdefs.label_font_size / 72.0 + fabs(gmtdefs.label_offset);

				if ((opt = strchr (txt_c, '+'))) {	/* Specified alternate label (could be upper case, hence 0.85) and justification */
					char txt_cpy[BUFSIZ], p[GMT_LONG_TEXT];
					GMT_LONG pos = 0;
					strcpy (txt_cpy, opt);
					while ((GMT_strtok (txt_cpy, "+", &pos, p))) {
						switch (p[0]) {
							case 'u':	/* Label put behind annotation */
								gave_label = FALSE;
								break;
							case 'j':	/* Justification */
								just = p[1];
								break;
							default:	/* Just ignore */
								break;
						}
					}
				}
				if (gave_label && just == 't') y0 += d_off;
				if (!strcmp (txt_a, "-"))	/* No longitude needed */
					sprintf (mapscale, "fx%gi/%gi/%s/%s", 0.5 * Ctrl->D.width, Ctrl->D.height - y0, txt_b, txt_c);
				else				/* Gave both lon and lat for scale */
					sprintf (mapscale, "fx%gi/%gi/%s/%s/%s", 0.5 * Ctrl->D.width, Ctrl->D.height - y0, txt_a, txt_b, txt_c);

				if (n_scan == 6)	/* Gave specific -R -J on M line */
					fprintf (fpo, "psbasemap %s %s -O -K -L%s\n", txt_e, txt_f, &mapscale[k]);
				else	/* Use -R -J supplied to pslegend */
					fprintf (fpo, "psbasemap %s %s -O -K -L%s\n", rarg, jarg, &mapscale[k]);
				/* Reset -R -J by calling a dummy psxy -T */
				fprintf (fpo, "psxy -R0/%g/0/%g -JX%gi/-%gi -O -K -T\n", Ctrl->D.width, Ctrl->D.height, Ctrl->D.width, Ctrl->D.height);
				if (gave_label && just == 'b') y0 += d_off;
				y0 += gmtdefs.map_scale_height + FONT_HEIGHT1 * gmtdefs.annot_font_size[0] / 72.0 + gmtdefs.annot_offset[0];
				column_number = 0;
				break;

			case 'S':	/* Symbol record */
				n_scan = sscanf (&line[2], "%s %s %s %s %s %s %[^\n]", txt_a, symbol, size, txt_c, txt_d, txt_b, text);
				off_ss = GMT_convert_units (txt_a, GMT_INCH);
				off_tt = GMT_convert_units (txt_b, GMT_INCH);
				d_off = 0.5 * (Ctrl->L.spacing - FONT_HEIGHT1) * gmtdefs.annot_font_size[0] / 72.0;	/* To center the text */
				if (column_number%n_columns == 0) y0 += one_line_spacing;
				y0 -= half_line_spacing;	/* Move to center of box */
				x_off = x0 + (Ctrl->D.width / n_columns) * (column_number%n_columns);
				if (symbol[0] == 'f') {	/* Front is different, must plot as a line segment */
					i = 0;
					while (size[i] != '/' && size[i]) i++;
					if (size[i] != '/') {
						fprintf (stderr, "%s: ERROR: -Sf option must have a tick length\n", GMT_program);
						exit (EXIT_FAILURE);
					}
					size[i] = '\0';	/* Temporarily truncate */
					x = 0.5 * GMT_convert_units (size, GMT_INCH);
					size[i] = '/';	/* Undo truncation */
					i++;
					fprintf (fpo, "echo %g %g > %s\necho %g %g >> %s\n", x_off + off_ss-x, y0, tmptxt, x_off + off_ss+x, y0, tmptxt);
					fprintf (fpo, "%s -S%s%s %s", psxy, symbol, &size[i], tmptxt);
					if (txt_c[0] != '-') fprintf (fpo, " -G%s", txt_c);
					if (txt_d[0] != '-') fprintf (fpo, " -W%s", txt_d);
					fprintf (fpo, "\n%s %s\n", del, tmptxt);
				}
				else {	/* Regular symbols */
					if (symbol[0] == 'k')
						sprintf (sub, "%s/%s", symbol, size);
					else
						sprintf (sub, "%s%s", symbol, size);
					if (symbol[0] == 'E' || symbol[0] == 'e') {	/* Ellipse needs more arguments */
						x = GMT_convert_units (size, gmtdefs.measure_unit);
						sprintf (sarg, "%g %g 0 %g %g", x_off + off_ss, y0, x, 0.65*x);
					}
					else if (symbol[0] == 'V' || symbol[0] == 'v') {	/* Vector also need more args */
						i = 0;
						while (size[i] != '/' && size[i]) i++;
						if (size[i] != '/') {	/* The necessary arguments not supplied! */
							sprintf (sub, "vb");
							exit (EXIT_FAILURE);
						}
						else {
							size[i++] = '\0';	/* So GMT_convert_units won'c complain */
							sprintf (sub, "%sb%s", symbol, &size[i]);
						}
						x = GMT_convert_units (size, gmtdefs.measure_unit);
						sprintf (sarg, "%g %g 0 %g", x_off + off_ss, y0, x);
					}
					else if (symbol[0] == 'r') {	/* Rectangle also need more args */
						x = GMT_convert_units (size, GMT_INCH);
						sprintf (sarg, "%g %g %g %g", x_off + off_ss, y0, x, 0.65*x);
					}
					else if (symbol[0] == 'w') {	/* Wedge also need more args */
						x = GMT_convert_units (size, GMT_INCH);
						sprintf (sarg, "%g %g 20 60", x_off + off_ss -0.5*x, y0+0.25*x);
					}
					else
						sprintf (sarg, "%g %g", x_off + off_ss, y0);
					fprintf (fpo, "echo %s | %s -S%s", sarg, psxy, sub);
					if (txt_c[0] != '-') fprintf (fpo, " -G%s", txt_c);
					if (txt_d[0] != '-') fprintf (fpo, " -W%s", txt_d);
					fprintf (fpo, "\n");
				}
				/* Finally, print text; skip when empty */
				y0 += half_line_spacing;	/* Go back to bottom of box */
				if (n_scan == 7) fprintf (fpo, "echo %c%g %g %g 0 %ld BL %s%c | %s\n", quote, x_off + off_tt, y0 - d_off, gmtdefs.annot_font_size[0], gmtdefs.annot_font[0], text, quote, pstext);
				column_number++;
				break;

			case 'D':	/* Delimiter record */
				sscanf (&line[2], "%s %s", txt_a, txt_b);
				L = GMT_convert_units (txt_a, GMT_INCH);
				y0 += quarter_line_spacing;
				fprintf (fpo, "echo %g %g > %s\necho %g %g >> %s\n", L, y0, tmptxt, Ctrl->D.width - L, y0, tmptxt);
				fprintf (fpo, "%s -W%s %s\n%s %s\n", psxy, txt_b, tmptxt, del, tmptxt);
				y0 += quarter_line_spacing;
				column_number = 0;
				break;

			case 'G':	/* Gap record */
				sscanf (&line[2], "%s", txt_a);
				y0 += (txt_a[strlen(txt_a)-1] == 'l') ? atoi (txt_a) * one_line_spacing : GMT_convert_units (txt_a, GMT_INCH);
				column_number = 0;
				break;

			case 'N':	/* n_columns record */
				sscanf (&line[2], "%s", txt_a);
				n_columns = atoi (txt_a);
				column_number = 0;
				break;

			case 'V':	/* Vertical line from here to next V */
				if (draw_vertical_line) {	/* Second time, now draw line */
					fprintf (fpo, "echo # vertical lines > %s\n", tmptxt);
					for (i = 1; i < n_columns; i++) {
						x_off = i * Ctrl->D.width / n_columns;
						fprintf (fpo, "echo %s> bar %ld >> %s\n", escape, i, tmptxt);
						fprintf (fpo, "echo %g %g >> %s\necho %g %g >> %s\n", x_off, y_start+V-quarter_line_spacing, tmptxt, x_off, y0-V+quarter_line_spacing, tmptxt);
					}
					fprintf (fpo, "%s -W%s -H -m %s\n%s %s\n", psxy, vpen, tmptxt, del, tmptxt);
					draw_vertical_line = FALSE;
				}
				else {
					draw_vertical_line = TRUE;
					y_start = y0;
					sscanf (&line[2], "%s %s", txt_a, vpen);
					V = GMT_convert_units (txt_a, GMT_INCH);
				}
				column_number = 0;
				break;

			case '>':	/* Paragraph text header */
				n = sscanf (&line[1], "%s %s %s %s %s %s %s %s %s", xx, yy, size, angle, font, key, lspace, tw, jj);
				if (n < 0) n = 0;	/* Since -1 is returned if no arguments */
				if (!(n == 0 || n == 9)) {
					fprintf (stderr, "%s: ERROR: The > record must have 0 or 9 arguments (only %ld found)\n", GMT_program, n);
					exit (EXIT_FAILURE);
				}
				d_off = 0.5 * (Ctrl->L.spacing - FONT_HEIGHT1) * gmtdefs.annot_font_size[0] / 72.0;
				if (n == 0 || xx[0] == '-') sprintf (xx, "%g", x0);
				if (n == 0 || yy[0] == '-') sprintf (yy, "%g", y0 + d_off);
				if (n == 0 || size[0] == '-') sprintf (size, "%g", gmtdefs.annot_font_size[0]);
				if (n == 0 || angle[0] == '-') sprintf (angle, "0");
				if (n == 0 || font[0] == '-') sprintf (font, "%ld", gmtdefs.annot_font[0]);
				if (n == 0 || key[0] == '-') sprintf (key, "TL");
				if (n == 0 || lspace[0] == '-') sprintf (lspace, "%gi", one_line_spacing);
				if (n == 0 || tw[0] == '-') sprintf (tw, "%gi", Ctrl->D.width - 2.0 * Ctrl->C.dx);
				if (n == 0 || jj[0] == '-') sprintf (jj, "j");
				fprintf (fpo, "echo %s> %s %s %s %s %s %s %s %s %s > %s\n", escape, xx, yy, size, angle, font, key, lspace, tw, jj, tmptxt);
				flush_paragraph = TRUE;
				column_number = 0;
				break;

			case 'T':	/* paragraph text record */
				d_off = 0.5 * (Ctrl->L.spacing - FONT_HEIGHT1) * gmtdefs.annot_font_size[0] / 72.0;
				/* If no previous > record, then use defaults */
				if (!flush_paragraph) fprintf (fpo, "echo %s> %g %g %g 0 %ld TL %gi %gi j > %s\n", escape, x0, y0 + d_off, gmtdefs.annot_font_size[0], gmtdefs.annot_font[0], one_line_spacing, Ctrl->D.width - 2.0 * Ctrl->C.dx, tmptxt);
				sscanf (&line[2], "%[^\n]", text);
				fprintf (fpo, "echo %c%s%c >> %s\n", quote, text, quote, tmptxt);
				flush_paragraph = TRUE;
				column_number = 0;
				break;

			default:
				fprintf (stderr, "%s: ERROR: Unrecognized record (%s)\n", GMT_program, line);
				exit (EXIT_FAILURE);
			break;
		}
	}
	if (fp != GMT_stdin) GMT_fclose (fp);

	if (flush_paragraph) fprintf (fpo, "%s -m %s\n%s %s\n", pstext, tmptxt, del, tmptxt);
	/* Revert to original region (-R), projection and size (-J) and position (-X, -Y) */
	fprintf (fpo, "psxy %s %s -T -X%c%gi -Y%c%gi -O %s", rarg, jarg, mode, -Ctrl->D.lon+GMT_ps.x_origin, mode, -Ctrl->D.lat+GMT_ps.y_origin, sparg);
	if (barg) fprintf (fpo, " %s", barg);
	if (!GMT_ps.last_page) fprintf (fpo, " -K");
	fprintf (fpo, "\n");

	if (!Ctrl->S.active) {	/* Add auto-delete command at the end of the script and then execute it */
		int err = 0;
		if (fpo != stdout) fclose (fpo);
		if (gmtdefs.verbose) fprintf (stderr, "%s: Executing and removing legend script\n", GMT_program);
#ifdef WIN32
		err = system (script);
#else
		sprintf (sub, "sh %s", script);
		err = system (sub);
#endif
		remove (script);
		if (err) fprintf (stderr, "%s: System call return with non-zero status %d\n", GMT_program, err);
	}
	else if (fpo != stdout)
		fclose (fpo);

	if (gmtdefs.verbose) fprintf (stderr, "%s: Done\n", GMT_program);

	Free_pslegend_Ctrl (Ctrl);	/* Deallocate control structure */

 	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void *New_pslegend_Ctrl () {	/* Allocate and initialize a new control structure */
	struct PSLEGEND_CTRL *C;
	
	C = (struct PSLEGEND_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct PSLEGEND_CTRL), "New_pslegend_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	
	C->C.dx = C->C.dy = 0.05;	
	C->D.width = C->D.height = 1.0;
	C->L.spacing = 1.1;
	return ((void *)C);
}

void Free_pslegend_Ctrl (struct PSLEGEND_CTRL *C) {	/* Deallocate control structure */
	if (C->S.file) free ((void *)C->S.file);	
	GMT_free ((void *)C);	
}
