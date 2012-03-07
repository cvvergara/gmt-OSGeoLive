/*--------------------------------------------------------------------
 *	$Id: gmt2kml.c,v 1.37 2011/07/11 19:22:06 guru Exp $
 *
 *	Copyright (c) 2009-2011 by P. Wessel
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
 * gmt2kml is a reformatter that takes GMT tables and converts them to
 * KML files for Google Earth
 *
 * Author:	Paul Wessel
 * Date:	19-MAY-2009
 * Version:	1.0
 */
 
#include "gmt.h"

#define POINT	0
#define EVENT	1
#define SPAN	2
#define LINE	3
#define POLYGON	4
#define KML_GROUND		0
#define KML_GROUND_REL		1
#define KML_ABSOLUTE		2
#define KML_SEAFLOOR_REL	3
#define KML_SEAFLOOR		4
#define KML_DOCUMENT		0
#define KML_FOLDER		1

struct EXT_COL {
	int col;			/* Column in input record */
	char name[GMT_TEXT_LEN];	/* Name of this column */
};

GMT_LONG check_lon_lat (double *lon, double *lat, double west, double east, double north, double south);
void print_altmode (int extrude, int fmode, int altmode);
void place_region_tag (double west, double east, double south, double north, double alt_min, double alt_max, GMT_LONG lod_min, GMT_LONG lod_max, double fade_min, double fade_max);
GMT_LONG ascii_output_one (double x, GMT_LONG col);

int main (int argc, char **argv)
{
	GMT_LONG i, j, k, fno, n_files = 0, n_args, fmode = POINT, altmode = KML_GROUND, n_coord = 0, itransp[5];
	GMT_LONG ix, iy, t1_col, t2_col, pos, t_itransp, n_ext_cols = 0, n_alloc = 0, n_rec = 0;
	GMT_LONG lod_min = 0, lod_max = 0, set_nr = 0, pnt_nr = 0, index = -4;

	GMT_LONG error = FALSE, long_verbose = FALSE, nofile = TRUE, done = FALSE, greenwich = FALSE, first = TRUE;
	GMT_LONG get_label = FALSE, no_label = FALSE, visible = TRUE, open = FALSE, use_folder = FALSE, format_label = FALSE;
	GMT_LONG add_description = FALSE, get_actual_wesn = FALSE, do_pen = TRUE, do_fill = TRUE, get_z = FALSE, get_alt = FALSE, get_rgb = FALSE, extrude = FALSE;
	
	char doc_title[BUFSIZ], folder_name[BUFSIZ], icon[BUFSIZ], description[BUFSIZ];
	char *feature[5] = {"Point", "Point", "Point", "LineString", "Polygon"}, *name[5] = {"Point", "Event", "Timespan", "Line", "Polygon"};
	char *indent[5] = {"\t\t\t", "\t\t\t", "\t\t\t", "\t\t\t", "\t\t\t\t\t"}, *fmt = NULL, buffer[BUFSIZ], label[BUFSIZ];
	char *Document[2] = {"Document", "Folder"}, p[BUFSIZ], C[4][GMT_TEXT_LEN], *c = NULL;

	int rgb[3];

	double west = 0.0, east = 0.0, south = 0.0, north = 0.0, altitude = 0.0, altscale = 1.0, alt_min = 0.0, alt_max = -1.0;
	double transparency[5] = {1.0, 1.0, 1.0, 1.0, 0.75}, t_transp = 1.0, t_scale = 1.0, scale = 1.0, fade_min = 0.0, fade_max = 0.0, out[5];

	FILE *fp = NULL, *fp_h = NULL;

	struct GMT_TABLE *line = NULL;
	struct GMT_FILL fill, t_fill;
	struct GMT_PEN pen;
	struct EXT_COL *D = NULL;

	argc = (int)GMT_begin (argc, argv);

	GMT_io.in_col_type[GMT_X] = GMT_io.out_col_type[GMT_X] = GMT_IS_LON;
	GMT_io.in_col_type[GMT_Y] = GMT_io.out_col_type[GMT_Y] = GMT_IS_LAT;
	GMT_init_fill (&fill, 255, 192, 128);		/* Default fill color */
	GMT_init_fill (&t_fill, 255, 255, 255);		/* Default text color */
	GMT_init_pen (&pen, 1.0);			/* Default pen width */
	doc_title[0] = folder_name[0] = label[0] = '\0';
	strcpy (icon, "http://maps.google.com/mapfiles/kml/pal4/icon57.png");
	
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {

				/* Common parameters */

				case 'V':
					if (argv[i][2] == 'l') long_verbose = TRUE;
				case 'H':
				case 'K':
				case 'M':
				case 'O':
				case ':':
				case 'b':
				case 'f':
				case 'm':
				case '\0':
					error += GMT_parse_common_options (argv[i], &west, &east, &south, &north);
					break;

				/* Supplemental parameters */

				case 'A':	/* Altitude mode */
					switch (argv[i][3]) {
						case 'x':
							altscale = atof (&argv[i][4]);
							get_alt = TRUE;
						case '\0':
							get_alt = TRUE;
							break;
						default:
							altitude = atof (&argv[i][3]);
							break;
					}
					switch (argv[i][2]) {
						case 'a':
							altmode = KML_ABSOLUTE;
							break;
						case 'g':
							altmode = (altitude != 0.0 || get_alt) ? KML_GROUND_REL : KML_GROUND;
							break;
						case 's':
							altmode = (altitude != 0.0 || get_alt) ? KML_SEAFLOOR_REL : KML_SEAFLOOR;
							break;
						default:
							fprintf (stderr, "%s: Bad altitude mode. Use a, g or s.\n", GMT_program);
							error++;
							break;
					}
					break;
				case 'C':	/* Color table */
					if (argv[i][2]) {
						error += GMT_read_cpt (&argv[i][2]);
						get_rgb = TRUE;
						if (GMT_continuous) {
							fprintf (stderr, "%s: Cannot use continuous color palette\n", GMT_program);
							error++;
						}
					}
					else {
						fprintf (stderr, "%s: Need to supply color palette name\n", GMT_program);
						error++;
					}
					break;
				case 'D':	/* Description file */
					if ((fp_h = GMT_fopen (&argv[i][2], GMT_io.r_mode)) == NULL) {
						fprintf (stderr, "%s: Cannot open HTML description file %s\n", GMT_program, &argv[i][2]);
						error++;
					}
					else
						add_description = TRUE;
					break;
					
				case 'E':	/* Extrude feature down to the ground*/
				 	extrude = 1;
					break;
				case 'F':	/* feature type */
					switch (argv[i][2]) {
						case 's':
							fmode = POINT;
							break;
						case 'e':
							fmode = EVENT;
							break;
						case 't':
							fmode = SPAN;
							break;
						case 'l':
							fmode = LINE;
							break;
						case 'p':
							fmode = POLYGON;
							break;
						default:
							fprintf (stderr, "%s: Bad feature type. Use s, e, t, l or p.\n", GMT_program);
							break;
					}
					break;
				case 'G':		/* Set fill for symbols or polygon */
					switch (argv[i][2]) {
						case 'f':	/* Symbol/polygon color fill */
							if (argv[i][3] == '-')
								do_fill = FALSE;
							else if (!argv[i][3] || GMT_getfill (&argv[i][3], &fill)) {
								GMT_fill_syntax ('G', " ");
								error++;
							}
							break;
						case 'n':	/* Label name color */
							if (argv[i][3] == '-')
								t_transp = 0.0;
							else if (!argv[i][3] || GMT_getfill (&argv[i][3], &t_fill)) {
								GMT_fill_syntax ('G', " ");
								error++;
							}
							break;
						default:
							GMT_fill_syntax ('G', " ");
							error++;
							break;
					}
					break;
				case 'I':	/* Custom icon */
					if (argv[i][2] == '+')
						sprintf (icon, "http://maps.google.com/mapfiles/kml/%s", &argv[i][3]);
					else if (argv[i][2])
						strcpy (icon, &argv[i][2]);
					break;
				case 'L':	/* Extended data */
					pos = n_ext_cols = 0;
					while ((GMT_strtok (&argv[i][2], ",", &pos, p))) {
						for (k = 0; p[k] && p[k] != ':'; k++);	/* Find position of colon */
						p[k] = ' ';
						if (n_ext_cols == n_alloc) n_alloc = GMT_alloc_memory ((void **)&D, n_ext_cols, n_alloc, sizeof (struct EXT_COL), GMT_program);
						sscanf (p, "%d %[^:]", &D[n_ext_cols].col, D[n_ext_cols].name);
						n_ext_cols++;
					}
					break;
				case 'N':	/* Feature label */
					if (argv[i][2] == '+') {	/* Special ASCII labelled input file */
						get_label = TRUE;
					}
					else if (!argv[i][2]) {	/* Want no label */
						t_transp = 0.0;
						no_label = TRUE;
					}
					else {
						fmt = &argv[i][2];
						format_label = TRUE;
					}
					break;
				case 'Q':	/* Transparency for symbols, lines, and polygons */
					switch (argv[i][2]) {
						case 'e':
						case 's':
						case 't':
							transparency[POINT] = transparency[EVENT] = transparency[SPAN] = atof (&argv[i][3]);
							break;
						case 'n':
							t_transp = atof (&argv[i][3]);
							break;
						case 'l':
							transparency[LINE] = atof (&argv[i][3]);
							break;
						case 'p':
							transparency[POLYGON] = atof (&argv[i][3]);
							break;
						default:	/* Same for all */
							t_transp = transparency[POINT] = transparency[EVENT] = transparency[SPAN] =  transparency[LINE] = transparency[POLYGON] = atof (&argv[i][2]);
							break;
					}
					break;
				case 'R':
					if (argv[i][2] == 'a')	/* Get args from data domain */
						get_actual_wesn = TRUE;
					else if (argv[i][2])
						error += GMT_parse_common_options (argv[i], &west, &east, &south, &north);
					break;
				case 'S':	/* Scale */
					if (argv[i][2] == 'f')
						scale = atof (&argv[i][3]);
					else if (argv[i][2] == 'n')
						t_scale = atof (&argv[i][3]);
					else {
						fprintf (stderr, "%s: -S requires f or n, then size\n", GMT_program);
						error++;
					}
					break;
				case 'T':	/* Title [and folder] */
					if ((c = strchr (&argv[i][2], '/'))) {	/* Got both title and folder */
						strcpy (folder_name, &c[1]);
						*c = '\0';
						strcpy (doc_title, &argv[i][2]);
						*c = '/';
					}
					else
						strcpy (doc_title, &argv[i][2]);
					break;
				case 'W':
					if (argv[i][2] == '-')
						do_pen = FALSE;
					else if (GMT_getpen (&argv[i][2], &pen)) {
						GMT_pen_syntax ('W', " ");
						error++;
					}
					break;
				case 'Z':
					pos = 0;
					while ((GMT_strtok (&argv[i][3], "+", &pos, p))) {
					switch (p[0]) {
						case 'a':	/* Altitude range */
							if (sscanf (&p[1], "%[^/]/%s", C[0], C[1]) != 2) {
								fprintf (stderr, "%s: -Z+a requires 2 arguments\n", GMT_program);
								error++;
							}
							alt_min = atof (C[0]);	alt_max = atof (C[1]);
							break;
						case 'l':	/* LOD */
							if (sscanf (&p[1], "%[^/]/%s", C[0], C[1]) != 2) {
								fprintf (stderr, "%s: -Z+l requires 2 arguments\n", GMT_program);
								error++;
							}
							lod_min = atoi (C[0]);	lod_max = atoi (C[1]);
							break;
						case 'f':	/* Fading */
							if (sscanf (&p[1], "%[^/]/%s", C[0], C[1]) != 2) {
								fprintf (stderr, "%s: -Z+f requires 2 arguments\n", GMT_program);
								error++;
							}
							fade_min = atoi (C[0]);	fade_max = atoi (C[1]);
							break;
						case 'v':
							visible = FALSE;
							break;
						case 'o':
							open = TRUE;
							break;
						default:
							fprintf (stderr, "%s: -Z unrecognized modifier +%c\n", GMT_program, p[0]);
							error++;
							break;
					}
				}
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
		fprintf (stderr, "gmt2kml %s - Convert GMT data tables to KML files for Google Earth\n\n", GMT_VERSION);
		fprintf (stderr, "usage: gmt2kml <infile> [-Aa|g|s[<altitude>|x<scale>]] [-D<descriptfile>] [-E] [-Fe|s<cpt>|t|l|p]\n");
		fprintf (stderr, "\t[-Gf|n-|[+]<fill>] [%s] [-I<icon>] [-K] [-L<col:name>,col:name>,...] [-N+|<template>|<name>]\n", GMT_H_OPT);
		fprintf (stderr, "\t[-O] [-Q[e|s|t|l|p|n]<transp>] [-Ra|<w/e/s/n>] [-Sc|n<scale>] [-T<title>/[/<foldername>] [-V] [-W-|<pen>] [-Z<opts>]\n");
		fprintf (stderr, "\t[%s] [-%s] [%s] [%s]\n\n", GMT_t_OPT, GMT_bi_OPT, GMT_f_OPT, GMT_m_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\tinfiles (in ASCII or binary) have 2 or more columns with (lon,lat) or (lat,lon) in first columns.\n");
		fprintf (stderr, "\t  If no file(s) is given, standard input is read.\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-A Altitude mode, choose among three modes:\n");
		fprintf (stderr,"\t      a Absolute altitude\n");
		fprintf (stderr,"\t      g Altitude relative to sea surface or ground\n");
		fprintf (stderr,"\t      s Altitude relative to seafloor or ground\n");
		fprintf (stderr,"\t    Optionally, append fixed <altitude>, or x<scale>. [g0: Clamped to sea surface or ground]\n");
		fprintf (stderr, "\t-D File with HTML snippets to use for data description [none]\n");
		fprintf (stderr, "\t-E Extend feature down to the ground [no extrusion]\n");
		fprintf (stderr, "\t-F Feature type; choose from (e)vent, (s)symbol, (t)imespan, (l)ine, or (p)olygon [s]\n");
		fprintf (stderr, "\t   Optionally append color palette name to -Fs to color icons by value. \n");
		fprintf (stderr, "\t   All features expect lon, lat in the first two columns. \n");
		fprintf (stderr, "\t   Value or altitude is given in the third column (see -A and -C)\n");
		fprintf (stderr, "\t   Event requires a timestamp in the next column.\n");
		fprintf (stderr, "\t   Timespan requires begin and end timestamps in the next two columns (use NaN for unlimited).\n");
		GMT_rgb_syntax ('G', "Specify color for symbol/polygon fill (f) [lightorange] or text label (n) [white].");
		fprintf (stderr, "\t   Use -Gf- to turn off polygon fill.\n");
		fprintf (stderr, "\t   Use -Gn- to turn off labels.\n");
		GMT_explain_option ('H');
		fprintf (stderr, "\t-I URL to an alternative icon used for the symbol [Google circle]\n");
		fprintf (stderr, "\t   If URL starts with + we will prepend http://maps.google.com/mapfiles/kml/\n");
		fprintf (stderr, "\t   [Default is a local icon with no directory path].\n");
		fprintf (stderr, "\t-K means allow for more KML code to be appended later [OFF]\n");
		fprintf (stderr, "\t-L Supply extended data informat via <col>:<name> strings [none]\n");
		fprintf (stderr, "\t-N Controls the feature labels.\n");
		fprintf (stderr, "\t   By default, -L\"label\" statements in the multisegment header are used. Alternatively,\n");
		fprintf (stderr, "\t   1. Specify -N+ if the rest of the data record should be used as label.\n");
		fprintf (stderr, "\t   2. Append a string that may contain the format %%d for a running feature count.\n");
		fprintf (stderr, "\t   3. Give no argument to indicate no labels\n");
		fprintf (stderr, "\t-O means append the KML code to an existing document [OFF]\n");
		fprintf (stderr, "\t-Q Set transparency level for selected feature, from 0 (transparent) to 1 (opaque)\n");
		fprintf (stderr, "\t   n is text label transparency.  [1 for points, lines, and labels; 0.75 for polygons]\n");
		fprintf (stderr, "\t-R Issue Region tag.  Append w/e/s/n to set a particular region or append a to use the\n");
		fprintf (stderr, "\t   actual domain of the data (single file only) [no region specified]\n");
		fprintf (stderr, "\t-S Scale for (c)ircle icon size or (n)ame label [1]\n");
		fprintf (stderr, "\t-T Append KML document title name [GMT Data Document]\n");
		fprintf (stderr, "\t   Optionally append /<foldername> to name folder when used with\n");
		fprintf (stderr, "\t   -O and -K to organize features into groups.\n");
		GMT_explain_option ('V');
		GMT_pen_syntax ('W', "Specify pen attributes for lines and polygons [Default is solid line of unit thickness].");
		fprintf (stderr, "\t   Give width in pixels and append p.  Use -W- to turn off polygon outlines.\n");
		fprintf (stderr, "\t-Z Control visibility of features.  Append one or more modifiers:\n");
		fprintf (stderr, "\t   +a<alt_min>/<alt_max> inserts altitude limits [no limit]\n");
		fprintf (stderr, "\t   +l<minLOD>/<maxLOD>] sets Level Of Detail when layer should be active [always active]\n");
		fprintf (stderr, "\t     layer goes inactive when there are fewer than minLOD pixels or more\n");
		fprintf (stderr, "\t     than maxLOD pixels visible.  -1 means never invisible.\n");
		fprintf (stderr, "\t   +f<minfade>/<maxfade>] sets distances over which we fade from opaque to transparent [no fading]\n");
		fprintf (stderr, "\t   +v turns off visibility [feature is visible]\n");
		fprintf (stderr, "\t   +o open document or folder when loaded [closed]\n");
		GMT_explain_option (':');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf (stderr, "\t   Default is 2 input columns.\n");
		GMT_explain_option ('f');
		GMT_explain_option ('m');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	if (get_rgb && fmode >= LINE) {
		fprintf (stderr, "%s: Color palette has no effect on plotting lines or polygons; -C option ignored\n", GMT_program);
		get_rgb = FALSE;
	}

	if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
        if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = 2;
	if (get_actual_wesn && n_files > 1) {	/* Only accept -Ra for single file */
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  -Ra without arguments only accepted for single table\n", GMT_program);
		error++;
	}
	for (i = 0; i < 5; i++) {
		if (transparency[i] < 0.0 || transparency[i] > 1.0) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR.  -Q takes transparencies in range 0-1\n", GMT_program);
			error++;
		}
		itransp[i] = irint (transparency[i] * 255.0);
	}
	t_itransp = irint (t_transp * 255.0);
	if (scale < 0.0 ) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  -S takes scales > 0.0\n", GMT_program);
		error++;
	}
	if (pen.width < 1.0) {	/* Must specify integer pen width */
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  -W given pen width < 1 pixel.  Use integers and append p as unit.\n", GMT_program);
		error++;
	}
	
	if (error) exit (EXIT_FAILURE);

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
	out[GMT_Z] = altitude;
	strcpy (gmtdefs.field_delimiter, ",");	/* Specify comma-separated output */
	GMT_io.geo.range = 2;			/* Want -180/+180 longitude output format */
	strcpy (gmtdefs.d_format, "%.12g");	/* Make sure we use enough decimals */
	ix = gmtdefs.xy_toggle[GMT_IN];	iy = 1 - ix;
	n_coord = (fmode < LINE) ? fmode + 2 : 2;
	get_z = (get_rgb || get_alt);
	if (get_z) n_coord++;
	t1_col = 2 + get_z;
	t2_col = 3 + get_z;
	if (fmode == EVENT || fmode == SPAN) GMT_io.in_col_type[t1_col] = GMT_io.out_col_type[t1_col] = GMT_IS_ABSTIME;
	if (fmode == SPAN)  GMT_io.in_col_type[t2_col] = GMT_io.out_col_type[t2_col] = GMT_IS_ABSTIME;

	if (!doc_title[0]) strcpy (doc_title, "GMT Data Document");
	if (!folder_name[0]) sprintf (folder_name, "%s Features", name[fmode]);
	if (GMT_ps.overlay || !GMT_ps.last_page) use_folder = TRUE;	/* When at least one or -O, -K is used */
	if (GMT_ps.overlay) {
		printf ("<%s>\n\t<name>%s</name>\n", Document[KML_FOLDER], folder_name);
	}
	else {
		/* Create KML header */
		printf ("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
		printf ("<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n");
		printf ("<%s>\n\t<name>%s</name>\n", Document[KML_DOCUMENT], doc_title);
		if (!visible) printf ("\t<visibility>0</visibility>\n");
		if (open) printf ("\t<open>1</open>\n");
		if (use_folder) printf ("<%s>\n\t<name>%s</name>\n", Document[KML_FOLDER], folder_name);
	}
	if (pen.width < 1.0) pen.width = 1.0;	/* Ensure we specify integer pen width */
	
	printf ("\t<Style id=\"GMT-4\">\n");

	/* Set icon style (applies to symbols only */
	printf ("\t\t<IconStyle>\n\t\t\t<scale>%g</scale>\n\t\t\t<color>%2.2x%2.2x%2.2x%2.2x</color>\n\t\t\t<Icon>\n\t\t\t\t<href>%s</href>\n\t\t\t</Icon>\n\t\t</IconStyle>\n",
		scale, (unsigned int)itransp[POINT], (unsigned int)fill.rgb[2], (unsigned int)fill.rgb[1], (unsigned int)fill.rgb[0], icon);

	/* Set shared line and polygon style (also for extrusions) */
	printf ("\t\t<LineStyle>\n\t\t\t<color>%2.2x%2.2x%2.2x%2.2x</color>\n\t\t\t<width>%d</width>\n\t\t</LineStyle>\n",
		(unsigned int)itransp[LINE], (unsigned int)pen.rgb[2], (unsigned int)pen.rgb[1], (unsigned int)pen.rgb[0], irint (pen.width));
	printf ("\t\t<PolyStyle>\n\t\t\t<color>%2.2x%2.2x%2.2x%2.2x</color>\n\t\t\t<fill>%ld</fill>\n\t\t\t<outline>%ld</outline>\n\t\t</PolyStyle>\n",
		(unsigned int)itransp[POLYGON], (unsigned int)fill.rgb[2], (unsigned int)fill.rgb[1], (unsigned int)fill.rgb[0], do_fill, do_pen);

	/* Set style for labels */
	printf ("\t\t<LabelStyle>\n\t\t\t<scale>%g</scale>\n\t\t\t<color>%2.2x%2.2x%2.2x%2.2x</color>\n\t\t</LabelStyle>\n",
		t_scale, (unsigned int)t_itransp, (unsigned int)t_fill.rgb[2], (unsigned int)t_fill.rgb[1], (unsigned int)t_fill.rgb[0]);
	printf ("\t</Style>\n");

	for (i = -3; get_rgb && i < GMT_n_colors; i++) {
		GMT_get_rgb_lookup (i, 0.0, rgb);
		printf ("\t<Style id=\"GMT%ld\">\n", i);
		/* Set icon style (applies to symbols only */
		printf ("\t\t<IconStyle>\n\t\t\t<scale>%g</scale>\n\t\t\t<color>%2.2x%2.2x%2.2x%2.2x</color>\n\t\t\t<Icon>\n\t\t\t\t<href>%s</href>\n\t\t\t</Icon>\n\t\t</IconStyle>\n",
			scale, (unsigned int)itransp[POINT], (unsigned int)rgb[2], (unsigned int)rgb[1], (unsigned int)rgb[0], icon);
		printf ("\t\t<LineStyle>\n\t\t\t<color>%2.2x%2.2x%2.2x%2.2x</color>\n\t\t\t<width>%d</width>\n\t\t</LineStyle>\n",
			(unsigned int)itransp[LINE], (unsigned int)pen.rgb[2], (unsigned int)pen.rgb[1], (unsigned int)pen.rgb[0], irint (pen.width));
		printf ("\t</Style>\n");
	}
	if (add_description) {
		char line[BUFSIZ];
		printf ("\t<description>\n\t\t<![CDATA[\n");
		while (GMT_fgets (line, BUFSIZ, fp_h)) printf ("\t\t\t%s", line);
		GMT_fclose (fp_h);
		printf ("\t\t]]>\n\t</description>\n");
	}
	
	for (fno = 1; !done && fno < n_args; fno++) {	/* Loop over input files, if any */
		if (!nofile && argv[fno][0] == '-') continue;

		if (nofile) {	/* Just read standard input */
			fp = GMT_stdin;
			done = TRUE;
			if (gmtdefs.verbose) fprintf (stderr, "%s: Reading from standard input\n", GMT_program);
#ifdef SET_IO_MODE
			GMT_setmode (GMT_IN);
#endif
		}
		else if ((fp = GMT_fopen (argv[fno], GMT_io.r_mode)) == NULL) {
			fprintf (stderr, "%s: Cannot open file %s\n", GMT_program, argv[fno]);
			continue;
		}
		set_nr = pnt_nr = 0;

		if (!nofile && gmtdefs.verbose) fprintf (stderr, "%s: Working on file %s\n", GMT_program, argv[fno]);
		if (get_label) {	/* Special ascii table processing */
			while (GMT_fgets (buffer, BUFSIZ, fp)) {
				if (buffer[0] == '#') continue;
				switch (n_coord) {
					case 2:	/* Just lon, lat, label */
						sscanf (buffer, "%s %s %[^\n]", C[ix], C[iy], label);
						break;
					case 3:	/* Just lon, lat, a, label */
						sscanf (buffer, "%s %s %s %[^\n]", C[ix], C[iy], C[2], label);
						break;
					case 4:	/* Just lon, lat, a, b, label */
						sscanf (buffer, "%s %s %s %s %[^\n]", C[ix], C[iy], C[2], C[3], label);
						break;
					case 5:	/* Just lon, lat, z, t1, t2, label */
						sscanf (buffer, "%s %s %s %s %s %[^\n]", C[ix], C[iy], C[2], C[3], C[4], label);
						break;
				}
				if (GMT_verify_expectations (GMT_io.in_col_type[GMT_X], GMT_scanf_arg (C[GMT_X], GMT_io.in_col_type[GMT_X], &out[GMT_X]), C[GMT_X])) {
					fprintf (stderr, "%s: ERROR: Could not decode longitude from %s\n", GMT_program, C[GMT_X]);
					exit (EXIT_FAILURE);
				}
				if (GMT_verify_expectations (GMT_io.in_col_type[GMT_Y], GMT_scanf_arg (C[GMT_Y], GMT_io.in_col_type[GMT_Y], &out[GMT_Y]), C[GMT_Y])) {
					fprintf (stderr, "%s: ERROR: Could not decode latitude from %s\n", GMT_program, C[GMT_Y]);
					exit (EXIT_FAILURE);
				}
				if (project_info.region_supplied && check_lon_lat(&out[GMT_X], &out[GMT_Y], west, east, south, north)) continue;
				if (get_z) {
					if (GMT_verify_expectations (GMT_io.in_col_type[GMT_Z], GMT_scanf_arg (C[GMT_Z], GMT_io.in_col_type[GMT_Z], &out[GMT_Z]), C[GMT_Z])) {
						fprintf (stderr, "%s: ERROR: Could not decode altitude from %s\n", GMT_program, C[GMT_Z]);
						exit (EXIT_FAILURE);
					}
					if (get_rgb) index = GMT_get_index (out[GMT_Z]);
					out[GMT_Z] = get_alt ? out[GMT_Z] * altscale : altitude;
				}
				if (fmode == EVENT) {
					if (GMT_verify_expectations (GMT_io.in_col_type[t1_col], GMT_scanf_arg (C[t1_col], GMT_io.in_col_type[t1_col], &out[t1_col]), C[t1_col])) {
						fprintf (stderr, "%s: ERROR: Could not decode time event from %s\n", GMT_program, C[t1_col]);
						exit (EXIT_FAILURE);
					}
				}
				else if (fmode == SPAN) {
					if (!(strcmp (C[t1_col], "NaN")))
						out[t1_col] = GMT_d_NaN;
					else if (GMT_verify_expectations (GMT_io.in_col_type[t1_col], GMT_scanf_arg (C[t1_col], GMT_io.in_col_type[t1_col], &out[t1_col]), C[t1_col])) {
						fprintf (stderr, "%s: ERROR: Could not decode time span beginning from %s\n", GMT_program, C[t1_col]);
						exit (EXIT_FAILURE);
					}
					if (!(strcmp (C[t2_col], "NaN")))
						out[t2_col] = GMT_d_NaN;
					else if (GMT_verify_expectations (GMT_io.in_col_type[t2_col], GMT_scanf_arg (C[t2_col], GMT_io.in_col_type[t2_col], &out[t2_col]), C[t2_col])) {
						fprintf (stderr, "%s: ERROR: Could not decode time span end from %s\n", GMT_program, C[t2_col]);
						exit (EXIT_FAILURE);
					}
				}
				if (project_info.region_supplied && first) {	/* Issue Region tag as given on commmand line*/
					place_region_tag (west, east, south, north, alt_min, alt_max, lod_min, lod_max, fade_min, fade_max);
					first = FALSE;
				}
				printf ("\t<Placemark>\n");
				printf ("\t\t<name>%s</name>\n", label);
				if (fmode == SPAN) {
					printf ("\t\t<TimeSpan>\n");
					if (!GMT_is_dnan(out[t1_col])) printf ("\t\t\t<begin>%s</begin>\n", C[t1_col]);
					if (!GMT_is_dnan(out[t2_col])) printf ("\t\t\t<end>%s</end>\n", C[t2_col]);
					printf ("\t\t</TimeSpan>\n");
				}
				else if (fmode == EVENT) printf ("\t\t<TimeStamp>\n\t\t\t<when>%s</when>\n\t\t</TimeStamp>\n", C[t1_col]);
				printf ("\t\t<styleUrl>#GMT%ld</styleUrl>\n", index);
				printf ("\t\t<%s>\n", feature[fmode]);
				print_altmode ((int)extrude, FALSE, (int)altmode);
				printf ("\t\t\t<coordinates>");
				ascii_output_one (out[GMT_X], GMT_X);	printf (",");
				ascii_output_one (out[GMT_Y], GMT_Y);	printf (",");
				ascii_output_one (out[GMT_Z], GMT_Z);
				printf ("</coordinates>\n");
				printf ("\t\t</%s>\n", feature[fmode]);
				printf ("\t</Placemark>\n");
				n_rec++;
				if (gmtdefs.verbose && !(n_rec%10000)) fprintf (stderr, "%s: Processed %ld points\n", GMT_program, n_rec);
			}
		}
		else {
			if (GMT_import_table ((void *)fp, GMT_IS_STREAM, &line, 0.0, greenwich, FALSE, TRUE) == GMT_IO_EOF) continue;	/* Empty file */
			if (project_info.region_supplied && first) {	/* Issue Region tag as given on commmand line*/
				place_region_tag (west, east, south, north, alt_min, alt_max, lod_min, lod_max, fade_min, fade_max);
				first = FALSE;
			}
			else if (get_actual_wesn) {	/* Issue Region tag */
				place_region_tag (line->min[GMT_X], line->max[GMT_X], line->min[GMT_Y], line->max[GMT_Y], alt_min, alt_max, lod_min, lod_max, fade_min, fade_max);
			}
			for (i = 0; i < line->n_segments; i++) {
				pnt_nr = 0;
				for (j = 0; j < line->segment[i]->n_rows; j++) {
					out[GMT_X] = line->segment[i]->coord[GMT_X][j];
					out[GMT_Y] = line->segment[i]->coord[GMT_Y][j];
					if (project_info.region_supplied && check_lon_lat(&out[GMT_X], &out[GMT_Y], west, east, south, north)) continue;
					if (get_z && line->n_columns > 2) {
						out[GMT_Z] = line->segment[i]->coord[GMT_Z][j];
						if (get_rgb) index = GMT_get_index(out[GMT_Z]);
						out[GMT_Z] = get_alt ? out[GMT_Z] * altscale : altitude;
					}
					if (fmode < LINE && GMT_is_fnan(out[GMT_Z])) continue;	/* Symbols with NaN height are not plotted anyhow */

					/* If first point, produce a segment header */
					if (pnt_nr + set_nr + nofile == 0) printf ("\t<Folder>\n\t\t<name>%s</name>\n", argv[fno]);
					if (pnt_nr > 0) { /* Nothing */ }
					else if (fmode < LINE) {
						printf ("\t<Folder>\n");
						if (line->segment[i]->label)
							printf ("\t\t<name>%s</name>\n", line->segment[i]->label);
						else
							printf ("\t\t<name>%s Set %ld</name>\n", name[fmode], set_nr);
					}
					else {
						printf ("\t<Placemark>\n");
						if (no_label) { /* Nothing */ }
						else if (format_label) {
							printf ("\t\t<name>"); printf (fmt, (int)set_nr); printf ("</name>\n");
						}
						else if (line->segment[i]->label)
							printf ("\t\t<name>%s</name>\n", line->segment[i]->label);
						else
							printf ("\t\t<name>%s %ld</name>\n", name[fmode], set_nr);
						if (GMT_parse_segment_item (line->segment[i]->header, "-D", description)) printf ("\t\t<description>%s</description>\n", description);
						printf ("\t\t<styleUrl>#GMT%ld</styleUrl>\n", index);
						printf ("\t\t<%s>\n", feature[fmode]);
						print_altmode ((int)extrude, (int)fmode, (int)altmode);
						if (fmode == POLYGON) {
							printf ("\t\t\t<outerBoundaryIs>\n\t\t\t\t<LinearRing>\n");
							if (line->segment[i]->min[GMT_X] < 180.0 && line->segment[i]->max[GMT_X] > 180.0) {
								/* GE cannot handle polygons crossing the dateline; warn for now */
								fprintf (stderr, "%s: Warning: A polygon is straddling the Dateline.  Google Earth will wrap this the wrong way\n", GMT_program);
								fprintf (stderr, "%s: Split the polygon into an East and West part and plot them as separate polygons.\n", GMT_program);
							}
						}
						printf ("%s<coordinates>\n", indent[fmode]);
					}

					/* Print the information for this point */
					if (fmode < LINE) {
						printf ("\t\t<Placemark>\n");
						if (no_label) { /* Nothing */ }
						else if (format_label) {
							printf ("\t\t\t<name>"); printf (fmt, (int)pnt_nr); printf ("</name>\n");
						}
						else if (line->segment[i]->label && line->segment[i]->n_rows > 1)
							printf ("\t\t\t<name>%s %ld</name>\n", line->segment[i]->label, j);
						else if (line->segment[i]->label)
							printf ("\t\t\t<name>%s</name>\n", line->segment[i]->label);
						else
							printf ("\t\t\t<name>%s %ld</name>\n", name[fmode], pnt_nr);
						if (line->segment[i]->n_rows == 1 && GMT_parse_segment_item (line->segment[i]->header, "-D", description)) printf ("\t\t<description>%s</description>\n", description);
						if (n_ext_cols) {
							printf ("\t\t\t<ExtendedData>\n");
							for (k = 0; k < n_ext_cols; k++) {
								printf ("\t\t\t\t<Data name = \"%s\">\n", D[k].name);
								printf ("\t\t\t\t\t<value>");
								ascii_output_one (line->segment[i]->coord[D[k].col][j], D[k].col);
								printf ("</value>\n\t\t\t\t</Data>\n");
							}
							printf ("\t\t\t</ExtendedData>\n");
						}
						if (fmode == SPAN) {
							printf ("\t\t\t<TimeSpan>\n");
							if (!GMT_is_dnan(line->segment[i]->coord[t1_col][j])) {
								printf ("\t\t\t\t<begin>");
								ascii_output_one (line->segment[i]->coord[t1_col][j], t1_col);
								printf ("</begin>\n");
							}
							if (!GMT_is_dnan(line->segment[i]->coord[t2_col][j])) {
								printf ("\t\t\t\t<end>");
								ascii_output_one (line->segment[i]->coord[t2_col][j], t2_col);
								printf ("</end>\n");
							}
							printf ("\t\t\t</TimeSpan>\n");
						}
						else if (fmode == EVENT) {
							printf ("\t\t\t<TimeStamp>\n\t\t\t\t<when>\n");
							ascii_output_one (line->segment[i]->coord[t1_col][j], t1_col);
							printf ("</when>\n\t\t\t</TimeStamp>\n");
						}
						printf ("\t\t\t<styleUrl>#GMT%ld</styleUrl>\n", index);
						printf ("\t\t\t<%s>\n\t", feature[fmode]);
						print_altmode ((int)extrude, FALSE, (int)altmode);
						printf ("\t\t\t\t<coordinates>");
						ascii_output_one (out[GMT_X], GMT_X);	printf (",");
						ascii_output_one (out[GMT_Y], GMT_Y);	printf (",");
						ascii_output_one (out[GMT_Z], GMT_Z);
						printf ("</coordinates>\n");
						printf ("\t\t\t</%s>\n", feature[fmode]);
						printf ("\t\t</Placemark>\n");
					}
					else {
						if (GMT_is_fnan(out[GMT_Z])) out[GMT_Z] = 0.0;	/* Google Earth can not handle lines at NaN altitude */
						printf ("%s\t", indent[fmode]);
						ascii_output_one (out[GMT_X], GMT_X);	printf (",");
						ascii_output_one (out[GMT_Y], GMT_Y);	printf (",");
						ascii_output_one (out[GMT_Z], GMT_Z);	printf ("\n");
					}
					pnt_nr++;
				}

				/* End of segment */
				if (pnt_nr == 0)
					set_nr--;
				else if (fmode < LINE)
					printf ("\t</Folder>\n");
				else {
					printf ("%s</coordinates>\n", indent[fmode]);
					if (fmode == POLYGON) printf ("\t\t\t\t</LinearRing>\n\t\t\t</outerBoundaryIs>\n");
					printf ("\t\t</%s>\n", feature[fmode]);
					printf ("\t</Placemark>\n");
				}
				set_nr++;
			}
			if (set_nr > 0 && !nofile) printf ("\t</Folder>\n");
			GMT_free_table (line);
		}
		if (fp != GMT_stdin) GMT_fclose (fp);
	}
	if (use_folder) printf ("</%s>\n", Document[KML_FOLDER]);
	if (GMT_ps.last_page) {
		printf ("</%s>\n", Document[KML_DOCUMENT]);
		printf ("</kml>\n");
	}

	GMT_free ((void *)D);
	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

GMT_LONG check_lon_lat (double *lon, double *lat, double west, double east, double south, double north)
{
	if (*lat < south || *lat > north) return (TRUE);
	if (*lon < west) *lon += 360.0;
	if (*lon > east) *lon -= 360.0;
	if (*lon < west) return (TRUE);
	return (FALSE);
}

void print_altmode (int extrude, int fmode, int altmode)
{
	char *RefLevel[5] = {"clampToGround", "relativeToGround", "absolute", "relativeToSeaFloor", "clampToSeaFloor"};
	if (extrude) printf ("\t\t\t<extrude>1</extrude>\n");
	if (fmode) printf ("\t\t\t<tesselate>1</tesselate>\n");
	if (altmode == KML_GROUND_REL || altmode == KML_ABSOLUTE) printf ("\t\t\t<altitudeMode>%s</altitudeMode>\n", RefLevel[altmode]);
	if (altmode == KML_SEAFLOOR_REL || altmode == KML_SEAFLOOR) printf ("\t\t\t<gx:altitudeMode>%s</gx:altitudeMode>\n", RefLevel[altmode]);
}

void place_region_tag (double west, double east, double south, double north, double alt_min, double alt_max, GMT_LONG lod_min, GMT_LONG lod_max, double fade_min, double fade_max)
{
	if (GMT_360_RANGE (west, east)) { west = -180.0; east = +180.0;}
	printf ("\t<Region>\n\t\t<LatLonAltBox>\n");
	printf ("\t\t\t<north>");		ascii_output_one (north, GMT_Y);
	printf ("</north>\n\t\t\t<south>");	ascii_output_one (south, GMT_Y);
	printf ("</south>\n");
	printf ("\t\t\t<east>");		ascii_output_one (east, GMT_X);
	printf ("</east>\n\t\t\t<west>");	ascii_output_one (west, GMT_X);
	printf ("</west>\n");
	if (alt_max > alt_min) {
		printf ("\t\t\t<minAltitude>%g</minAltitude>\n", alt_min);
		printf ("\t\t\t<maxAltitude>%g</maxAltitude>\n", alt_max);
	}
	printf ("\t\t</LatLonAltBox>\n");
	if (lod_max != lod_min) {
		printf ("\t\t<Lod>\n");
		printf ("\t\t\t<minLodPixels>%ld</minLodPixels>\n", lod_min);
		printf ("\t\t\t<maxLodPixels>%ld</maxLodPixels>\n", lod_max);
		if (fade_min > 0.0 || fade_max > 0.0) {
			printf ("\t\t\t<minFadeExtent>%g</minFadeExtent>\n", fade_min);
			printf ("\t\t\t<maxFadeExtent>%g</maxFadeExtent>\n", fade_max);
		}
		printf ("\t\t</Lod>\n");
	}
	printf ("\t</Region>\n");
}

GMT_LONG ascii_output_one (double x, GMT_LONG col)
{	/* Used instead of GMT_asfii_output_one since Windoze has trouble with GMT_stdout and stdout mixing */
	char text[GMT_LONG_TEXT];

	GMT_ascii_format_one (text, x, GMT_io.out_col_type[col]);
	return (printf ("%s", text));
}
