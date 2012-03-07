/*--------------------------------------------------------------------
 *	$Id: psmegaplot.c,v 1.21 2011/07/11 19:22:06 guru Exp $
 *
 *      Copyright (c) 1999-2011 by P. Wessel
 *      See LICENSE.TXT file for copying and redistribution conditions.
 *
 *      This program is free software; you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation; version 2 or any later version.
 *
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *
 *      Contact info: www.soest.hawaii.edu/wessel
 *--------------------------------------------------------------------*/
/*
 * psmegaplot allows a regular GMT produced postscript file to be magnified 
 * and split up into pieces that are plotted out to make up a jigsaw puzzle.
 * The source PS file is assumed to have any translate/offset/scale statements
 * reversed at the end.
 *
 * Author:	Paul Wessel
 * Date:	21-MAY-1991-1998
 * Version:	3.0 PW: File can have showpage; we skip everything after Trailer
 */
 
#include "gmt.h"

#define XL 575	/* Actual plotsize on Letter size paper, origin is <18,8> (varies from printer to printer) */
#define YL 775

int main (int argc, char **argv)
{
	GMT_LONG x0, y0, i, j, nx, ny, page;
	double scale = 0.0;
	FILE *fp = NULL;
	GMT_LONG first, crop_marks = FALSE, error = FALSE, done;
	char ifile[100], buffer[BUFSIZ];

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				case 'C' :
					crop_marks = TRUE;
					break;
				case 'S' :
					scale = atof (&argv[i][2]);
					break;
				case '\0' :
					GMT_give_synopsis_and_exit = TRUE;
					break;
				default :
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else
			strcpy (ifile, argv[i]);
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf(stderr, "psmegaplot %s - Make postersize plot using tiling\n\n", GMT_VERSION);
		fprintf(stderr, "usage : psmegaplot psfile -S<scale> [-C]\n");

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf(stderr, "	-S sets scale, must be > 1.0\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf(stderr, "	-C means plot crop marks at corners\n");
		exit (EXIT_FAILURE);
	}

	if (scale <= 1.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -S option:  scale must be larger than 1\n", GMT_program);
		error++;
	}
	if (error) exit (EXIT_FAILURE);

	if ((fp = fopen(ifile, "r")) == NULL) {
		fprintf(stderr, "psmegaplot: Could not open file %s\n", ifile);
		exit (EXIT_FAILURE);
	}

	ny = nx = (GMT_LONG)ceil (scale);
	first = TRUE;
	page = 0;
	printf("%%!PS\n");
	for (i = 0; i < nx; i++) {
		x0 = XL * i;
		for (j = 0; j < ny; j++) {
			y0 = YL * j;
			printf ("%%%%Page: poster %ld\n\n", ++page);
			if (crop_marks) {
				printf ("1 setlinewidth\n");
				printf ("18 9 moveto 0 -1 rlineto 1 0 rlineto\n");
				printf ("574 0 rmoveto 1 0 rlineto 0 1 rlineto\n");
				printf ("0 774 rmoveto 0 1 rlineto -1 0 rlineto\n");
				printf ("-574 0 rmoveto -1 0 rlineto 0 -1 rlineto stroke\n");
			}
			printf ("%ld %ld translate\n", -x0, -y0);
			printf ("%.3f %.3f scale\n", scale, scale);
			rewind (fp);
			done = FALSE;
			while (!done && fgets (buffer, BUFSIZ, fp)) {	/* Read 1 line at the time */
				if (!first && strstr (buffer, "dict begin")) continue;
				if (!first && strstr (buffer, "gsave")) continue;
				if (strstr (buffer, "setpagedevice")) continue;
				if (strstr (buffer, "%%BoundingBox")) continue;
				if (strstr (buffer, "%%Trailer")) {
					done = TRUE;
					continue;
				}
				if (buffer[0] == '/' && first)
					printf("%s", buffer);
				else if (buffer[0] != '/')
					printf("%s", buffer);
			}
			first = FALSE;
			printf ("%.3f %.3f scale\n", 1./scale, 1./scale);
			printf ("%ld %ld translate showpage\n", x0, y0);
		}
	}
	printf ("end\n");
	fclose(fp);
	fprintf (stderr, "psmegaplot: Produced %ld pages\n", page);

	exit (EXIT_SUCCESS);
}
