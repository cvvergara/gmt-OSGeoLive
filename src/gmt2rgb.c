/*--------------------------------------------------------------------
 *	$Id: gmt2rgb.c 9923 2012-12-18 20:45:53Z pwessel $
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
 * gmt2rgb reads either (1) an 8, 24, or 32 bit Sun rasterfile and writes out the
 * red, green, and blue components in three separate grid files, or (2) a z grd
 * file and a cpt file and compute r, g, b and write these out instead.
 *
 * Author:	Paul Wessel
 * Date:	17-SEP-2001
 * Version:	4.1.2
 *
 */

#include "gmt.h"
#include "pslib.h"

struct GMT2RGB_CTRL {
	struct C {	/* -C<cptfile> */
		GMT_LONG active;
		char *file;
	} C;
	struct F {	/* -F */
		GMT_LONG active;
	} F;
	struct G {	/* -G<nametemplate> */
		GMT_LONG active;
		char *name;
	} G;
	struct I {	/* -Idx[/dy] */
		GMT_LONG active;
		double xinc, yinc;
	} I;
	struct L {	/* -L<layer> */
		GMT_LONG active;
		char layer;
	} L;
	struct W {	/* -W<width/height>[/<n_bytes>] */
		GMT_LONG active;
		GMT_LONG nx, ny;	/* Dimension of image */
		GMT_LONG size;	/* Number of bytes per pixels */
	} W;
};

int main (int argc, char **argv)
{
	GMT_LONG i, error = 0, entry, one_or_zero, pos, nm, k, k3;
	int irgb[3];
	GMT_LONG guess = FALSE;
	float *z = NULL;
	double w, e, s, n;
	struct GRD_HEADER grd;
	struct imageinfo header;
	char rgb[3] = {'r', 'g', 'b'}, *comp[3] = {"red", "green", "blue"};
	char *file = CNULL, grdfile[BUFSIZ], ptr[GMT_TEXT_LEN];
	unsigned char *picture = NULL;
	struct GMT2RGB_CTRL *Ctrl = NULL;
	
	void guess_width (char *file, GMT_LONG byte_per_pixel, GMT_LONG *raw_nx, GMT_LONG *raw_ny);
	unsigned char *loadraw (char *file, struct imageinfo *header, GMT_LONG byte_per_pixel, GMT_LONG nx, GMT_LONG ny);
	void *New_gmt2rgb_Ctrl (), Free_gmt2rgb_Ctrl (struct GMT2RGB_CTRL *C);

	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct GMT2RGB_CTRL *)New_gmt2rgb_Ctrl ();	/* Allocate and initialize a new control structure */

	/* Check and interpret the command line arguments */

	w = e = s = n = 0.0;
	
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch(argv[i][1]) {

				/* Common parameters */

				case 'R':
				case 'V':
				case '\0':
					error += GMT_parse_common_options (argv[i], &w, &e, &s, &n);
					break;
				case 'C':
					Ctrl->C.file = strdup (&argv[i][2]);
					Ctrl->C.active = TRUE;
					break;
				case 'F':
					Ctrl->F.active = TRUE;
					break;
				case 'G':
					Ctrl->G.name = strdup (&argv[i][2]);
					break;
				case 'I':
					GMT_getinc (&argv[i][2], &Ctrl->I.xinc, &Ctrl->I.yinc);
					Ctrl->I.active = TRUE;
					break;
				case 'L':
					Ctrl->L.layer = argv[i][2];
					Ctrl->L.active = TRUE;
					break;
				case 'W':
					Ctrl->W.active = TRUE;
					guess = TRUE;
					entry = pos = 0;
					while ((GMT_strtok (&argv[i][2], "/", &pos, ptr))) {
						if (ptr[0] != '=') {
							switch (entry) {
								case 0:
									Ctrl->W.nx = atoi (ptr);
									guess = FALSE;
									break;
								case 1:
									Ctrl->W.ny = atoi (ptr);
									break;
								case 2:
									Ctrl->W.size = atoi (ptr);
									break;
								default:
									break;
							}
						}
						entry++;
					}
					break;

				/* Options not recognized */

				default:
					error = TRUE;
					fprintf (stderr, "GMT SYNTAX ERROR:  Unrecognized option -%c\n", argv[i][1]);
					break;
			}
		}
		else
			file = argv[i];
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr,"%s %s - Write r/g/b grid files from a grid file, a raw RGB file, or SUN rasterfile\n\n", GMT_program, GMT_VERSION);
		fprintf (stderr,"usage: %s <rasterfile|rgbfile||grdfile> [-C<cptfile>] [-F] [-G<nametemplate>] [%s]\n", GMT_program, GMT_Id_OPT);
		fprintf (stderr,"\t[-L<layer>] [%s] [-V] [-W<width/height>[/<n_bytes>]]\n\n", GMT_Rgeo_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr,"\t<infile> can be one of three different intput files:\n");
		fprintf (stderr,"\t  (1) An 8, 24, or 32-bit Sun rasterfile.  Use -I, -R, and -F to change the\n");
		fprintf (stderr,"\t      the default values of dx = dy = 1 and region 1/ncols/1/nrows.\n");
		fprintf (stderr,"\t  (2) A regular z grid file.  Use -C to provide a cpt file with which\n");
		fprintf (stderr,"\t      to convert z to r/g/b triplets. -R, -I, and -F are ignored.\n");
		fprintf (stderr,"\t  (3) A RGB or RGBA raw rasterfile. Since raw rasterfiles have no header, you have to\n");
		fprintf (stderr,"\t      give the image dimensions via -W\n");
		fprintf (stderr,"\t      However, you may take the chance of letting the program try to\n");
		fprintf (stderr,"\t      guess the image dimensions (slow, via FFT spectrum).\n");
		fprintf (stderr,"\n\tOPTIONS:\n");
		fprintf (stderr,"\t-C color palette file to convert z to rgb.  If given, we assume a z grid file is provided,\n");
		fprintf (stderr,"\t   else we will try to read a Sun rasterfile.\n");
		fprintf (stderr, "\t-F will force pixel registration [Default is grid registration].\n");
		fprintf (stderr,"\t-G Give outputfile name template for the three red, green, blue grid files.\n");
		fprintf (stderr,"\t   The template MUST contain the format code %%c which will be replaced with r, g, and b\n");
		fprintf (stderr,"\t   [Default is gmt2rgb_%%c.grd].\n");
		fprintf (stderr, "\t-I specifies grid size(s).  Append m (or c) to <dx> and/or <dy> for minutes (or seconds).\n");
		fprintf (stderr, "\t-L Only output the given layer (r, g, or b) [Default output all three].\n");
		GMT_explain_option ('R');
		GMT_explain_option ('V');
		fprintf (stderr, "\t-W sets the size of the raw raster file. By default an RGB file (which has 3 bytes/pixel)\n");
		fprintf (stderr, "\t   is assumed. For RGBA files use n_bytes = 4.\n");
		fprintf (stderr, "\t   Use -W for guessing the image size of a RGB raw file, and -W=/=/4\n");
		fprintf (stderr, "\t   if the raw image is of the RGBA type. Notice that this might be a\n");
		fprintf (stderr, "\t   bit slow because the guessing algorithm makes uses of FFTs.\n");
		exit (EXIT_FAILURE);
	}

	/* Check that the options selected are mutually consistent */

	GMT_check_lattice (&Ctrl->I.xinc, &Ctrl->I.yinc, &Ctrl->F.active, &Ctrl->I.active);

	if (!Ctrl->C.active) {
		if (!file) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify input raster file\n", GMT_program);
			error++;
		}
		if (!Ctrl->I.active && !Ctrl->W.active) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify -Idx/dy\n", GMT_program);
			error++;
		}
		if (Ctrl->I.active && (Ctrl->I.xinc == 0.0 || Ctrl->I.yinc == 0.0)) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR:  increments must be positive\n", GMT_program);
			error++;
		}
		if (Ctrl->W.size != 3 && Ctrl->W.size != 4) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR: byte_per_pixel must be either 3 or 4\n", GMT_program);
			error++;
		}
		if (guess)
			guess_width (file, Ctrl->W.size, &Ctrl->W.nx, &Ctrl->W.ny);

		if (Ctrl->W.active && Ctrl->W.nx <= 0) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR:  Witdth of raw raster file must be a positive integer. Not %ld\n", GMT_program, Ctrl->W.nx);
			error++;
		}
		if (Ctrl->W.active && Ctrl->W.ny <= 0) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR:  Height of raw raster file must be a positive integer. Not %ld\n", GMT_program, Ctrl->W.ny);
			error++;
		}
	}
	else {
		if (!file) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify input z grid file\n", GMT_program);
			error++;
		}
		if (!Ctrl->C.file) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify cpt file\n", GMT_program);
			error++;
		}
	}
	if (!strstr (Ctrl->G.name, "%c")) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  output template must contain %%c\n", GMT_program);
		error++;
	}
	if (Ctrl->I.active && !strchr ("rgb", Ctrl->L.layer)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  -L layer must be one of r, g, or b\n", GMT_program);
		error++;
	}
	if (error) exit (EXIT_FAILURE);

	GMT_grd_init (&grd, argc, argv, FALSE);

	if (Ctrl->C.active) {
		GMT_read_cpt (Ctrl->C.file);
		GMT_err_fail (GMT_read_grd_info (file, &grd), file);
		nm = GMT_get_nm (grd.nx, grd.ny);
		
		z = (float *) GMT_memory (VNULL, (size_t)nm, sizeof (float), GMT_program);
		for (i = 0; i < 3; i++) {	/* Do the r, g, and b channels */
			if (Ctrl->L.active && Ctrl->L.layer != rgb[i]) continue;
			if (gmtdefs.verbose) fprintf (stderr, "%s: Processing the %s components\n", GMT_program, comp[i]);
			GMT_err_fail (GMT_read_grd (file, &grd, z, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), file);
			sprintf (grdfile, Ctrl->G.name, rgb[i]);
			sprintf (grd.remark, "Grid of %s components in the 0-255 range", comp[i]);
			for (k = 0; k < nm; k++) {
				GMT_get_rgb_from_z (z[k], irgb);
				z[k] = (float)irgb[i];
			}
			GMT_err_fail (GMT_write_grd (grdfile, &grd, z, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), grdfile);
		}
	}
	else {
		if (GMT_access (file, R_OK)) {
			fprintf (stderr, "%s: Cannot find/open/read file %s\n", GMT_program, file);
			exit (EXIT_FAILURE);
		}

		if (!Ctrl->W.active)
			picture = ps_load_image (file, &header);
		else
			picture = loadraw (file, &header, Ctrl->W.size, Ctrl->W.nx, Ctrl->W.ny);

		if (!picture) {
			fprintf (stderr, "%s: Trouble loading/converting Sun rasterfile!\n", GMT_program);
			exit (EXIT_FAILURE);
		}
		if (header.depth < 8) {
			fprintf (stderr, "%s: Sun rasterfile must be at least 8 bits deep\n", GMT_program);
			exit (EXIT_FAILURE);
		}

		if (Ctrl->F.active) {
			grd.node_offset = 1;
			one_or_zero = 0;
		}
		else {
			grd.node_offset = 0;
			one_or_zero = 1;
		}
		if (!Ctrl->I.active) {
			if (gmtdefs.verbose) fprintf (stderr, "%s: Assign default dx = dy = 1\n", GMT_program);
			Ctrl->I.xinc = Ctrl->I.yinc = 1.0;
		}
		if (w == e && s == n) {	/* R not given, provide default */
			w = s = 0.0;
			e = (double)(header.width - one_or_zero);
			n = (double)(header.height - one_or_zero);
			if (gmtdefs.verbose) fprintf (stderr, "%s: Assign default -R0/%g/0/%g\n", GMT_program, e, n);
		}

		grd.nx = (int)(irint ((e-w)/Ctrl->I.xinc) + one_or_zero);
		grd.ny = (int)(irint ((n-s)/Ctrl->I.yinc) + one_or_zero);
		if (Ctrl->W.active && !Ctrl->I.active) {		/* This isn't correct because it doesn't deal with -F */
			grd.nx = (int)Ctrl->W.nx;
			grd.ny = (int)Ctrl->W.ny;
			Ctrl->I.xinc = (e-w)/(Ctrl->W.nx - one_or_zero);
			Ctrl->I.yinc = (n-s)/(Ctrl->W.ny - one_or_zero);
		}
		if (header.width != grd.nx) {
			fprintf (stderr, "%s: Sun rasterfile width and -R -I do not match (%d versus %d)  Need -F?\n", GMT_program, header.width, grd.nx);
			exit (EXIT_FAILURE);
		}
		if (header.height != grd.ny) {
			fprintf (stderr, "%s: Sun rasterfile height and -R -I do not match (%d versus %d)  Need -F?\n", GMT_program, header.height, grd.ny);
			exit (EXIT_FAILURE);
		}
		grd.x_min = w;	grd.x_max = e;
		grd.y_min = s;	grd.y_max = n;
		grd.x_inc = Ctrl->I.xinc;	grd.y_inc = Ctrl->I.yinc;
		nm = GMT_get_nm (grd.nx, grd.ny);

		GMT_err_fail (GMT_grd_RI_verify (&grd, 1), Ctrl->G.name);

		if (gmtdefs.verbose) fprintf (stderr, "%s: nx = %d  ny = %d\n", GMT_program, grd.nx, grd.ny);

		z = (float *) GMT_memory (VNULL, (size_t)nm, sizeof (float), GMT_program);

		for (i = 0; i < 3; i++) {	/* Do the r, g, and b channels */
			if (Ctrl->L.active && Ctrl->L.layer != rgb[i]) continue;
			if (gmtdefs.verbose) fprintf (stderr, "%s: Processing the %s components\n", GMT_program, comp[i]);
			sprintf (grdfile, Ctrl->G.name, rgb[i]);
			sprintf (grd.remark, "Grid of %s components in the 0-255 range", comp[i]);
			k3 = i;
			for (k = 0; k < nm; k++) {
				if (header.depth == 8)	/* Gray ramp */
					z[k] = (float)picture[k];
				else {				/* 24-bit image */
					z[k] = (float)picture[k3];
					k3 += 3;
				}
			}
			GMT_err_fail (GMT_write_grd (grdfile, &grd, z, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), grdfile);
		}
		GMT_free ((void *)picture);
	}
	GMT_free ((void *)z);

	Free_gmt2rgb_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

unsigned char *loadraw (char *file, struct imageinfo *header, GMT_LONG byte_per_pixel, GMT_LONG nx, GMT_LONG ny) {
	/* loadraw reads a raw binary grb or rgba rasterfile of depth 24, or 32 into memory */

	GMT_LONG j, i, nm;
	unsigned char *buffer;

	FILE *fp;

	if ((fp = GMT_fopen (file, "rb")) == NULL) {
		fprintf (stderr, "%s: Cannot open rasterfile %s!\n", GMT_program, file);
		exit (EXIT_FAILURE);
	}

	/* Lets pretend that the raw file is a sunraster file. This way the gmt2rgb code
	   can be used with very little changes */
	header->depth = 24;
	header->width = (int)nx;
	header->height = (int)ny;
	nm = ((size_t)nx) * ((size_t)ny) * ((size_t)byte_per_pixel);
	header->length = (int)nm;

	buffer = (unsigned char *) GMT_memory (VNULL, (size_t)nm, sizeof (unsigned char), GMT_program);
	if (GMT_fread ((void *)buffer, (size_t)1, (size_t)nm, fp) != (size_t)nm) {
		if (byte_per_pixel == 3)
			fprintf (stderr, "%s: Trouble reading raw 24-bit rasterfile!\n", GMT_program);
		if (byte_per_pixel == 4)
			fprintf (stderr, "%s: Trouble reading raw 32-bit rasterfile!\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	if (byte_per_pixel == 4) {		/* RGBA */
		for (i = 3, j = 4; j < nm; i += 3, j += 4) {
			buffer[i] = buffer[j];
			buffer[i+1] = buffer[j+1];
			buffer[i+2] = buffer[j+2];
		}
	}

	GMT_fclose (fp);
	return (buffer);
}

void guess_width (char *file, GMT_LONG byte_per_pixel, GMT_LONG *raw_nx, GMT_LONG *raw_ny) {
	unsigned char *buffer = NULL;
	float	*work = NULL, *datac = NULL, *img_pow = NULL, pow_max = -FLT_MAX, pm;
	GMT_LONG k = 0, j, inc, i, l, even, narray, img_size, n_pix;
	int rgb[3];
	FILE *fp = NULL;

	if ((fp = GMT_fopen (file, "rb")) == NULL) {
		fprintf (stderr, "%s: Cannot open rasterfile %s!\n", GMT_program, file);
		exit (EXIT_FAILURE);
	}

	GMT_fseek (fp, 0L, SEEK_END);
	img_size = GMT_ftell (fp);
	GMT_fseek (fp, 0L, SEEK_SET);

	n_pix = img_size / byte_per_pixel;

	buffer = (unsigned char *) GMT_memory (VNULL, (size_t)img_size, sizeof (unsigned char), GMT_program);
	datac = (float *) GMT_memory (VNULL, (size_t)2*n_pix, sizeof (float), GMT_program);
	work = (float *) GMT_memory (VNULL, (size_t)2*n_pix, sizeof(float), GMT_program);
	img_pow = (float *) GMT_memory (VNULL, (size_t)n_pix/2, sizeof (float), GMT_program);
	memset ((char *)work, 0, (size_t)(2*n_pix * sizeof(float)));

	if (GMT_fread ((void *)buffer, (size_t)1, (size_t)img_size, fp) != (size_t)img_size) {
		if (byte_per_pixel == 3)
			fprintf (stderr, "%s: Trouble_ reading raw 24-bit rasterfile!\n", GMT_program);
		if (byte_per_pixel == 4)
			fprintf (stderr, "%s: Trouble_ reading raw 32-bit rasterfile!\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	inc = (byte_per_pixel == 3) ? 3: 4;
	for (j = 0; j < img_size; j += inc) {
		rgb[0] = buffer[j];
		rgb[1] = buffer[j+1];
		rgb[2] = buffer[j+2];
		/* Convert rgb to gray using the GMT_YIQ transformation */
		datac[k] = (float) GMT_YIQ(rgb);
		k += 2;
	}

	narray = n_pix;
	GMT_fourt (datac, &narray, 1, -1, 1, work);

	/* Now compute the image's power spectrum */
	for (k = 0, j = 0; k < n_pix; k+= 2, j++) {
		img_pow[j] = (datac[k]*datac[k] + datac[k+1]*datac[k+1]) / n_pix; /* I*I-conj = power */
	}

	/* I'll assume that searching on one fifth of the spectrum is enough to find the
	   line frequency. */
	for (k = 5; k < n_pix/10; k++) {
		if (img_pow[k] > pow_max) {
			pow_max = img_pow[k];	 j = k+1;
		}
	}

	/* That's the way it should be but, I don't know why, the result is transposed. Instead
	   of the number of lines I get number of columns. This is very weard and smels	BUG */
	/* *raw_nx = j;		*raw_ny = irint((float)n_pix / raw_nx);*/

	/* So be it */
	*raw_ny = j;		*raw_nx = irint((float)n_pix / (*raw_ny));

	if ((*raw_nx) * (*raw_ny) != n_pix) {
		/* Let's make another attempt to find the right nx * ny combination. The idea is that we
	   	failed by a little, so we'll look arround the approximate solution adding 1 to nx and
	   	subtracting 1 to ny. Then we revert (subtract 1 to nx and add 1 to ny). Next apply the
	   	same test with an offset of 2, and so on until the offset is 10. */
		fprintf (stderr, "%s WARNING: first test based on FFT failed to guess image dimensions.\n\tI'll do now a second try\t", GMT_program);
		k = 1;		pm = 1;		l = 1;
		while (k < 41) {
			i = *raw_ny + (int)irint (copysign((double)l, (double)pm));
			pm *= -1.;
			j = (*raw_nx) + (int)irint (copysign((double)l, (double)pm));
			if (i*j == n_pix) {	/* Got a good candidate */
				*raw_ny = i;	*raw_nx = j;
				fprintf (stderr, "... SUCESS (W = %ld, H = %ld)\n", *raw_nx, *raw_ny);
				break;
			}
			even = (k%2 == 0) ? 1: 0;
			if (even) l++;
			k++;
		}
	}
	else
		if (gmtdefs.verbose) fprintf (stderr, "File %s has %ld Lines and %ld Cols\n", file, *raw_ny, *raw_nx);

	/* If both attempts failed */
	if ((*raw_nx) * (*raw_ny) != n_pix) {
		fprintf (stderr, "FAILURE while guessing image dimensions (W = %ld, H = %ld)\n", *raw_nx, *raw_ny);
		exit (EXIT_FAILURE);
	}

	GMT_fclose (fp);
	GMT_free ((void *)buffer);
	GMT_free ((void *)datac);
	GMT_free ((void *)work);
	GMT_free ((void *)img_pow);
}

void *New_gmt2rgb_Ctrl () {	/* Allocate and initialize a new control structure */
	struct GMT2RGB_CTRL *C;
	
	C = (struct GMT2RGB_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GMT2RGB_CTRL), "New_gmt2rgb_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	
	C->G.name = strdup ("gmt2rgb_%c.grd");
	C->W.size = 3;	/* 3 bytes per pixel */
		
	return ((void *)C);
}

void Free_gmt2rgb_Ctrl (struct GMT2RGB_CTRL *C) {	/* Deallocate control structure */
	if (C->C.file) free ((void *)C->C.file);	
	if (C->G.name) free ((void *)C->G.name);	
	GMT_free ((void *)C);	
}
