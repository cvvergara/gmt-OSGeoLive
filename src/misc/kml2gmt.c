/*--------------------------------------------------------------------
 *	$Id: kml2gmt.c 10173 2014-01-01 09:52:34Z pwessel $
 *
 *	Copyright (c) 2009-2014 by P. Wessel
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
 * kml2gmt is a reformatter that takes KML files and extracts GMT tables;
 * it is the opposite of gmt2kml.
 *
 * Author:	Paul Wessel
 * Date:	19-MAY-2009
 * Version:	1.0
 */
 
#include "gmt.h"

#define POINT			0
#define LINE			1
#define POLYGON			2

int main (int argc, char **argv)
{
	GMT_LONG i, n_files = 0, ncol = 2, start, fmode = POINT;
	GMT_LONG scan = TRUE, first = TRUE;

	GMT_LONG error = FALSE;
	
	char line[BUFSIZ], header[BUFSIZ], name[BUFSIZ], description[BUFSIZ], *file = NULL;

	double out[3];

	FILE *fp = NULL;

	GMT_LONG ascii_output_one (double x, GMT_LONG col);

	argc = (int)GMT_begin (argc, argv);

	GMT_io.in_col_type[GMT_X] = GMT_io.out_col_type[GMT_X] = GMT_IS_LON;
	GMT_io.in_col_type[GMT_Y] = GMT_io.out_col_type[GMT_Y] = GMT_IS_LAT;
	memset ((void *)header, 0, BUFSIZ);
	memset ((void *)name, 0, BUFSIZ);
	memset ((void *)description, 0, BUFSIZ);
	
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {

				/* Common parameters */

				case 'V':
				case 'b':
				case 'm':
				case ':':
				case '\0':
					error += GMT_parse_common_options (argv[i], NULL, NULL, NULL, NULL);
					break;

				/* Supplemental parameters */

				case 'Z':
					ncol = 3;
					break;
				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else {
			file = argv[i];
			n_files++;
		}
	}

	if (error || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "kml2gmt %s - Extract GMT data from a Google Earth KML file\n\n", GMT_VERSION);
		fprintf (stderr, "usage: kml2gmt [<infile>] [-V] [%s] [%s] [%s] > GMTdata.txt\n", GMT_m_OPT, GMT_t_OPT, GMT_bo_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\tinfile is the Google Earth KML file.\n");
		fprintf (stderr, "\t  If no file(s) is given, standard input is read.\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "-Z Output the z-column from the KML file [Only lon,lat is output]\n");
		GMT_explain_option ('V');
		GMT_explain_option ('m');
		GMT_explain_option (':');
		GMT_explain_option ('o');
		GMT_explain_option ('n');
		GMT_explain_option ('.');

		exit (EXIT_FAILURE);
	}

	if (n_files > 1) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Only one file can be processed at the time\n", GMT_program);
		error++;
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
        if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = 2;
	GMT_io.multi_segments[GMT_OUT] = TRUE;	/* So we can write segments */
	if (error) exit (EXIT_FAILURE);

	if (GMT_io.binary[GMT_IN] && gmtdefs.verbose) {
		char *type[2] = {"double", "single"};
		fprintf (stderr, "%s: Expects %ld-column %s-precision binary data\n", GMT_program, GMT_io.ncol[GMT_IN], type[GMT_io.single_precision[GMT_IN]]);
	}
#ifdef SET_IO_MODE
	GMT_setmode (GMT_OUT);
#endif
	/* Now we are ready to take on some input values */

	strcpy (gmtdefs.d_format, "%.12g");	/* Get enough decimals */
	
	if (n_files == 0) {	/* Just read standard input */
		fp = stdin;
		if (gmtdefs.verbose) fprintf (stderr, "%s: Reading from standard input\n", GMT_program);
#ifdef SET_IO_MODE
		GMT_setmode (GMT_IN);
#endif
		printf ("# %s: KML read from standard input\n", GMT_program);
	}
	else if ((fp = fopen (file, GMT_io.r_mode)) == NULL) {
		fprintf (stderr, "%s: Cannot open file %s\n", GMT_program, file);
		exit (EXIT_FAILURE);
	}
	else {
		if (gmtdefs.verbose) fprintf (stderr, "%s: Processing %s\n", GMT_program, file);
		printf ("# %s: KML read from %s\n", GMT_program, file);
	}

	while (fgets (line, BUFSIZ, fp)) {
		if (strstr (line, "<Placemark")) scan = TRUE;
		if (strstr (line, "</Placemark")) scan = FALSE;
		if (!scan) continue;
		if (strstr (line, "<Point")) fmode = POINT;
		if (strstr (line, "<LineString")) fmode = LINE;
		if (strstr (line, "<Polygon")) fmode = POLYGON;
		if (strstr (line, "<name>")) {
			for (i = 0; i < (GMT_LONG)strlen (line) && line[i] != '>'; i++);	/* Find end of <name> */
			start = i + 1;
			for (i = start; i < (GMT_LONG)strlen (line) && line[i] != '<'; i++);	/* Find start of </name> */
			line[i] = '\0';
			strcpy (name, &line[start]);
			GMT_chop (name);
			if (first) printf ("# %s\n", &line[start]);
			first = FALSE;
		}
		if (strstr (line, "<description>")) {
			for (i = 0; i < (GMT_LONG)strlen (line) && line[i] != '>'; i++);	/* Find end of <description> */
			start = i + 1;
			for (i = start; i < (GMT_LONG)strlen (line) && line[i] != '<'; i++);	/* Find start of </description> */
			line[i] = '\0';
			strcpy (description, &line[start]);
			GMT_chop (description);
			if (first) printf ("# %s\n", &line[start]);
			first = FALSE;
		}
		if (name[0] || description[0]) {
			sprintf (GMT_io.segment_header, "%c", GMT_io.EOF_flag[GMT_OUT]);
			if (name[0]) { strcat (GMT_io.segment_header, " -L\""); strcat (GMT_io.segment_header, name); strcat (GMT_io.segment_header, "\""); }
			if (description[0]) { strcat (GMT_io.segment_header, " -D\""); strcat (GMT_io.segment_header, description); strcat (GMT_io.segment_header, "\""); }
			strcat (GMT_io.segment_header, "\n");
		}
		
		if (!strstr (line, "<coordinates>")) continue;
		/* We get here when the line says coordinates */
		if (fmode == POINT) {	/* Process the single point */
			for (i = 0; i < (GMT_LONG)strlen (line) && line[i] != '>'; i++);		/* Find end of <coordinates> */
			sscanf (&line[i+1], "%lg,%lg,%lg", &out[GMT_X], &out[GMT_Y], &out[GMT_Z]);
			if (gmtdefs.xy_toggle[GMT_OUT]) d_swap (out[GMT_X], out[GMT_Y]);		/* Output lat/lon instead of lon/lat */
			ascii_output_one (out[GMT_X], GMT_X);	printf ("%s", gmtdefs.field_delimiter);
			ascii_output_one (out[GMT_Y], GMT_Y);	
			if (ncol == 3) { printf ("%s", gmtdefs.field_delimiter); ascii_output_one (out[GMT_Z], GMT_Z);}
			printf ("\n");
		}
		else {
			if (GMT_io.segment_header[0]) printf ("%s", GMT_io.segment_header);
			name[0] = description[0] = 0;
			while (fscanf (fp, "%lg,%lg,%lg", &out[GMT_X], &out[GMT_Y], &out[GMT_Z])) {
				if (gmtdefs.xy_toggle[GMT_OUT]) d_swap (out[GMT_X], out[GMT_Y]);		/* Output lat/lon instead of lon/lat */
				ascii_output_one (out[GMT_X], GMT_X);	printf ("%s", gmtdefs.field_delimiter);
				ascii_output_one (out[GMT_Y], GMT_Y);	
				if (ncol == 3) { printf ("%s", gmtdefs.field_delimiter); ascii_output_one (out[GMT_Z], GMT_Z);}
				printf ("\n");
			}
		}
	}
	
	if (fp != stdin) fclose (fp);
	
	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

GMT_LONG ascii_output_one (double x, GMT_LONG col)
{	/* Used instead of GMT_asfii_output_one since Windoze has trouble with GMT_stdout and stdout mixing */
	char text[GMT_LONG_TEXT];

	GMT_ascii_format_one (text, x, GMT_io.out_col_type[col]);
	return (printf ("%s", text));
}
