/*--------------------------------------------------------------------
 *	$Id: makepattern.c 9923 2012-12-18 20:45:53Z pwessel $
 *
 *      Copyright (c) 1999-2013 by P. Wessel
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
/* makepattern is a utility for creating color imagepatterns to be used
 * with psimage or as pattern fill -Gi<size>/<image>.  It accepts either
 * an icon file OR a 1-bit standard Sun rasterfile and
 * assigns a foreground and background color where the bitimage is 1 or
 * 0, respectively, and writes out an 8-bit standard Sun rasterfile.
 *
 *--------------------------------------------------------------------*/
/*
  
   Author: 	P. Wessel
   Date:	19 MAR, 1998.
   version:	1.1
		1.2, 12-APR-1999 (DOS-compatible)

*/

#include "gmt.h"
#include "pslib.h"

/* Structure for Sun rasterfile */

int main (int argc, char **argv)
{
	GMT_LONG pad = FALSE, icon = FALSE;

	int i, j, k, kk, nx, ny, extra, mx, step, error = 0;

	unsigned char map[6], *rgb = NULL;
#ifdef _WIN32
	char *ofile = NULL;
#endif

	unsigned int p;
	unsigned short int *pattern = VNULL;

	char *fname = CNULL;

	struct GMT_FILL fg, bg;

	struct imageinfo h;

	unsigned short int *get_icon (char *file, int *nx, int *ny, int invert);

	FILE *fp = NULL, *fpo = NULL;

	GMT_init_fill (&fg, 0, 0, 0);		/* Default foreground is black */
	GMT_init_fill (&bg, 255, 255, 255);	/* Default background is white */

	argc = (int)GMT_begin (argc, argv);

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {


				/* Common parameters */

				case '\0':
					error += (int)GMT_parse_common_options (argv[i], 0, 0, 0, 0);
					break;

				/* Supplemental parameters */

				case 'C':
					switch (argv[i][2]) {
						case 'f':
							if (GMT_getfill (&argv[i][3], &fg)) {
								GMT_fill_syntax ('C', " ");
								error++;
							}
							break;
						case 'b':
							if (GMT_getfill (&argv[i][3], &bg)) {
								GMT_fill_syntax ('C', " ");
								error++;
							}
							break;
						default:
							error++;
					}
					break;
#ifdef _WIN32
				case 'G':
					ofile = &argv[i][2];
					break;
#endif
				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else {
			fname = argv[i];
			if ((fp = fopen (fname, "rb")) == NULL) {
				fprintf(stderr,"%s: Cannot open file %s\n", GMT_program, fname);
				exit (EXIT_FAILURE);
			}
		}
	}

#ifdef _WIN32
	if (!ofile) {
		fprintf(stderr,"%s: No output file given\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if ((fpo = fopen (ofile, "wb")) == NULL) {
		fprintf(stderr,"%s: Cannot create file %s\n", GMT_program, ofile);
		exit (EXIT_FAILURE);
	}
#else
	fpo = stdout;
#endif
	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "makepattern %s - make color pattern from b/w pattern\n\n", GMT_VERSION);
#ifdef _WIN32
		fprintf(stderr,"usage:	makepattern 1bit.ras OR iconfile -Cf<rgb> -Cb<rgb> -Gimage.ras\n\n");
#else
		fprintf(stderr,"usage:	makepattern 1bit.ras OR iconfile -Cf<rgb> -Cb<rgb> > image.ras\n\n");
#endif

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

			fprintf(stderr,"\tGive the name of a 1-bit standard Sun rasterfile or iconfile\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf(stderr,"\t-Cf  Set foreground r/g/b color [0/0/0]\n");
		fprintf(stderr,"\t-Cb  Set background r/g/b color [255/255/255]\n");
#ifdef _WIN32
		fprintf(stderr,"\t-G  Name of output Sun rasterfile\n");
#endif
		exit (EXIT_FAILURE);
	}

	if (!fp) {
		fprintf(stderr,"%s: Specify either a 1-bit rasterfile or an iconfile\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	if (fp) {	/* Got a 1-bit Sun rasterfile */

		if (ps_read_rasheader (fp, &h, 0, 7)) {
			fprintf (stderr, "%s: Trouble reading Sun rasterfile header!\n", GMT_program);
			exit (EXIT_FAILURE);
		}
		if (h.magic != RAS_MAGIC) {	/* Try to read as icon instead */
			fclose (fp);
			pattern = get_icon (fname, &nx, &ny, FALSE);
			icon = TRUE;
		}
		else {

			if (h.depth != 1) {
				fprintf (stderr, "%s: File is not a 1-bit Sun rasterfile!\n", GMT_program);
				exit (EXIT_FAILURE);
			}

			nx = h.width;
			ny = h.height;
			mx = h.length/2;
			pattern = (unsigned short int *) GMT_memory (VNULL, (size_t)mx, sizeof (unsigned short int), GMT_program);
			if (fread ((void *)pattern, sizeof (unsigned short int), (size_t)mx, fp) != (size_t)mx) {
				fprintf (stderr, "%s: Trouble reading 1-bit Sun rasterfile image!\n", GMT_program);
				exit (EXIT_FAILURE);
			}
			fclose (fp);
		}
	}

	/* Fill out the 2-entry color table with chosen back and foreground colors */

	map[0] = bg.rgb[0];	map[1] = fg.rgb[0];
	map[2] = bg.rgb[1];	map[3] = fg.rgb[1];
	map[4] = bg.rgb[2];	map[5] = fg.rgb[2];

	/* Allocated byte array for 8-bit image; each row must be even # of bytes, pad if necessary */

	mx = nx;
	if (nx%2) {
		pad = TRUE;
		mx++;	/* Must be even number of bytes */
	}
	rgb = (unsigned char *) GMT_memory (VNULL, (size_t)(mx * ny), sizeof (unsigned char), GMT_program);

	mx = nx / 16;			/* Number of full 16-bit shorts */
	extra = nx - mx * 16;		/* Remainder of bits in the last short */
	step = (extra) ? mx + 1 : mx;	/* Number of shorts per row */

	for (j = kk = 0; j < ny; j++) {		/* For each row in image */

		for (i = 0; i < mx; i++) {	/* For each chunk of full 16 bits */

			for (k = 0; k < 16; k++, kk++) {	/* Deal with each bit */

				p = (32768 >> k);
				rgb[kk] = (((unsigned int)(pattern[j*step+i]) & p) > 0);

			}
		}
		if (extra) {	/* Deal with remainder of bits in last short */
			i++;
			for (k = 0; k < extra; k++, kk++) {
				p = (32768 >> k);
				rgb[kk] = (((unsigned int)(pattern[j*step+i]) & p) > 0);
			}
		}
		if (pad) kk++;
	}

	/* Fill out header structure for rasterfile */

	mx = nx;
	if (mx%2) mx++;	/* Must be even number of bytes */

	h.magic = RAS_MAGIC;
	h.width = nx;
	h.height = ny;
	h.depth = 8;
	h.length = mx * ny;
	h.type  = RT_STANDARD;
	h.maptype = RMT_EQUAL_RGB;
	h.maplength = 6;

	if (ps_write_rasheader (fp, &h, 0, 7)) {
		fprintf (stderr, "%s: Trouble writing Sun rasterfile header!\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	if (fwrite ((void *)map, sizeof (unsigned char), (size_t)h.maplength, fpo) != (size_t)h.maplength) {
		fprintf (stderr, "%s: Trouble writing Sun rasterfile colormap!\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	if (fwrite ((void *)rgb, sizeof (unsigned char), (size_t)h.length, fpo) != (size_t)h.length) {
		fprintf (stderr, "%s: Trouble writing Sun rasterfile image!\n", GMT_program);
		exit (EXIT_FAILURE);
	}

#ifdef _WIN32
	fclose (fpo);
#endif

	GMT_free ((void *)rgb);
	if (icon)
		free ((void *)pattern);
	else
		GMT_free ((void *)pattern);

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}


unsigned short int *get_icon (char *file, int *nx, int *ny, int invert)
{
	int i, n_items, n_int, j, k, n_read, n_alloc, n_lines, last;
	unsigned int integer;
	unsigned short int *pattern;
	char t[8][7], line[80], width[12], height[12], *not_used = NULL;
	FILE *fp;

	if ((fp = fopen (file, "r")) == NULL) {
		fprintf (stderr, "%s: Cannot open file %s\n", GMT_program, file);
		return ((unsigned short int *)NULL);
	}

	not_used = fgets (line, 80, fp);
	line[strlen(line)-1] = 0;
	sscanf (&line[3], "%*[^,], %[^,], %[^,]", width, height);
	i = 0;	while (width[i] != '=') i++;
	*nx = atoi (&width[i+1]);
	i = 0;	while (height[i] != '=') i++;
	*ny = atoi (&height[i+1]);
	not_used = fgets (line, 80, fp);
	line[strlen(line)-1] = 0;
	n_items = (*nx) * (*ny) / 4;
	n_int = n_items / 4;
	n_alloc = (int)ceil (n_int / 8.0) * n_int;
	pattern = (unsigned short int *) malloc ((size_t)n_alloc);
	n_lines = n_items / 32;
	for (i = k = 0; i < n_lines; i++) {
		not_used = fgets (line, 80, fp);
		line[strlen(line)-1] = 0;
		last = (i == (n_lines - 1));
		if (last) {
			n_read = sscanf (&line[1], "%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,]", t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7]);
			if (n_read == 4) {
				sscanf (&line[1], "%[^,], %[^,], %[^,], %[^,]", t[0], t[1], t[2], t[3]);	/* New style icon file */
				not_used = fgets (line, 80, fp);
				line[strlen(line)-1] = 0;
				sscanf (&line[1], " %[^,], %[^,], %[^,], %[^,]", t[4], t[5], t[6], t[7]);	/* New style icon file */
			}
		}
		else {
			n_read = sscanf (&line[1], "%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],", t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7]);
			if (n_read == 4) {
				sscanf (&line[1], "%[^,], %[^,], %[^,], %[^,],", t[0], t[1], t[2], t[3]);	/* New style icon file */
				not_used = fgets (line, 80, fp);
				line[strlen(line)-1] = 0;
				sscanf (&line[1], " %[^,], %[^,], %[^,], %[^,],", t[4], t[5], t[6], t[7]);	/* New style icon file */
			}
		}
		for (j = 0; j < 8; j++, k++) {
			sscanf (&t[j][2], "%x", &integer);
			pattern[k] = integer;
		}
	}
	fclose (fp);
	if (invert) for (k = 0; k < n_int; k++) pattern[k] = ~pattern[k];
	return (pattern);
}
