/*--------------------------------------------------------------------
 *	$Id: psimage.c,v 1.56 2011/07/08 21:27:06 guru Exp $
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
 * psimage reads a 1, 8, 24, or 32 bit Sun rasterfile and plots it on the page
 *
 * Author:	Paul Wessel
 * Date:	28-JUN-2000
 * Version:	4
 *
 */

#include "gmt.h"
#include "pslib.h"

struct PSIMAGE_CTRL {
	struct GMT_CONTOUR contour;
	struct C {	/* -C<xpos>/<ypos>[/<justify>] */
		GMT_LONG active;
		double x, y;
		char justify[3];
	} C;
	struct E {	/* -E<dpi> */
		GMT_LONG active;
		double dpi;
	} E;
	struct F {	/* -F<pen> */
		GMT_LONG active;
		struct GMT_PEN pen;
	} F;
	struct G {	/* -G[f|b|t]<rgb> */
		GMT_LONG active;
		int f_rgb[3];
		int b_rgb[3];
		int t_rgb[3];
	} G;
	struct I {	/* -I */
		GMT_LONG active;
	} I;
	struct M {	/* -M */
		GMT_LONG active;
	} M;
	struct N {	/* -N<nx>/<ny> */
		GMT_LONG active;
		GMT_LONG nx, ny;
	} N;
	struct W {	/* -W[-]<width>[/<height>] */
		GMT_LONG active;
		GMT_LONG interpolate;
		double width, height;
	} W;
};

int main (int argc, char **argv)
{
	GMT_LONG i, j, n, justify, PS_interpolate = 1, PS_transparent = 1;

	GMT_LONG error = FALSE;

	double x, y;

	char file[BUFSIZ], txt_a[GMT_LONG_TEXT], txt_b[GMT_LONG_TEXT], letter;

	unsigned char *picture = NULL, *buffer = NULL;

	struct imageinfo header;

	struct PSIMAGE_CTRL *Ctrl = NULL;

	void *New_psimage_Ctrl (), Free_psimage_Ctrl (struct PSIMAGE_CTRL *C);
	
	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct PSIMAGE_CTRL *)New_psimage_Ctrl ();	/* Allocate and initialize a new control structure */

	/* Check and interpret the command line arguments */

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch(argv[i][1]) {

				/* Common parameters */

				case 'K':
				case 'O':
				case 'P':
				case 'U':
				case 'V':
				case 'X':
				case 'x':
				case 'Y':
				case 'y':
				case 'c':
				case ':':
				case '\0':
					error += GMT_parse_common_options (argv[i], 0, 0, 0, 0);
					break;

				/* Supplemental parameters */

				case 'C':
					Ctrl->C.active = TRUE;
					n = sscanf (&argv[i][2], "%[^/]/%[^/]/%2s", txt_a, txt_b, Ctrl->C.justify);
					if (n < 2 || n > 3) {
						fprintf (stderr, "%s ERROR: Syntax is -C<xpos>/<ypos>[/<justify>]\n", GMT_program);
						error++;
					}
					Ctrl->C.x = GMT_convert_units (txt_a, GMT_INCH);
					Ctrl->C.y = GMT_convert_units (txt_b, GMT_INCH);
					if (n == 2) strcpy (Ctrl->C.justify, "LB");	/* Default positioning */
					break;
				case 'E':	/* Specify image dpi */
					Ctrl->E.active = TRUE;
					Ctrl->E.dpi = atof (&argv[i][2]);
					break;
				case 'F':	/* Specify frame pen */
					Ctrl->F.active = TRUE;
					if (argv[i][2] && GMT_getpen (&argv[i][2], &Ctrl->F.pen)) {
						GMT_pen_syntax ('F', " ");
						error++;
					}
					break;
				case 'G':
					Ctrl->G.active = TRUE;
					letter = (GMT_colorname2index (&argv[i][2]) >= 0) ? 'x' : argv[i][2];	/* If we have -G<colorname>, the x is used to bypass the case F|f|B|b switching below */
					switch (letter) {
						case 'F':
						case 'f':
							/* Set color for foreground pixels */
							if (argv[i][3] == '-' && argv[i][4] == '\0')
								Ctrl->G.f_rgb[0] = -1;
							else if (GMT_getrgb (&argv[i][3], Ctrl->G.f_rgb)) {
								GMT_rgb_syntax ('G', " ");
								error++;
							}
							break;
						case 'B':
						case 'b':
							/* Set color for background pixels */
							if (argv[i][3] == '-' && argv[i][4] == '\0')
								Ctrl->G.b_rgb[0] = -1;
							else if (GMT_getrgb (&argv[i][3], Ctrl->G.b_rgb)) {
								GMT_rgb_syntax ('G', " ");
								error++;
							}
							break;
						case 'T':
						case 't':
							/* Set transparent color */
							if (GMT_getrgb (&argv[i][3], Ctrl->G.t_rgb)) {
								GMT_rgb_syntax ('G', " ");
								error++;
							}
							break;
						default:	/* Gave either -G<r/g/b>, -G-, or -G<colorname>; all treated as -Gf */
							if (argv[i][2] == '-' && argv[i][3] == '\0')
								Ctrl->G.f_rgb[0] = -1;
							else if (GMT_getrgb (&argv[i][2], Ctrl->G.f_rgb)) {
								GMT_rgb_syntax ('G', " ");
								error++;
							}
							break;
					}
					break;
				case 'I':
					Ctrl->I.active = TRUE;
					break;
				case 'M':
					Ctrl->M.active = TRUE;
					break;
				case 'N':
					Ctrl->N.active = TRUE;
					n = sscanf (&argv[i][2], "%" GMT_LL "d/%" GMT_LL "d", &Ctrl->N.nx, &Ctrl->N.ny);
					if (n == 1) Ctrl->N.ny = Ctrl->N.nx;
					if (n < 1) {
						error++;
						fprintf (stderr, "%s: GMT SYNTAX ERROR -N option:  Must values for replication\n", GMT_program);
					}
					break;
				case 'W':
					Ctrl->W.active = TRUE;
					n = sscanf (&argv[i][2], "%[^/]/%s", txt_a, txt_b);
					Ctrl->W.width = GMT_convert_units (txt_a, GMT_INCH);
					if (n == 2) Ctrl->W.height = GMT_convert_units (txt_b, GMT_INCH);
					if (Ctrl->W.width < 0.0) {
						Ctrl->W.width = -Ctrl->W.width;
						Ctrl->W.interpolate = TRUE;
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
			strcpy (file, argv[i]);
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr,"psimage %s - To plot image file on maps\n\n", GMT_VERSION);
		fprintf (stderr,"usage: psimage <imagefile> [-E<dpi> or -W[-]<width>[/<height>]] [-C<xpos>/<ypos>[/<justify>]]\n");
		fprintf (stderr, "\t[-F<pen>] [-G[b|f|t]<color>] [-I] [-K] [-M] [-N<nx>[/<ny>]] [-O] [-P] [%s]\n", GMT_U_OPT);
		fprintf (stderr, "\t[-V] [%s] [%s] [%s]\n\n", GMT_X_OPT, GMT_Y_OPT, GMT_c_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr,"\t<imagefile> is an EPS file or a 1, 8, 24, or 32-bit Sun rasterfile.\n");
		fprintf (stderr, "\t-E sets image dpi (dots per inch), OR\n");
		fprintf (stderr, "\t-W sets the width (and height) of the image.  If <height> = 0\n");
		fprintf (stderr, "\t   then the original aspect ratio is maintained.  If <width> < 0\n");
		fprintf (stderr, "\t   then we use absolute value and interpolate image in PostScript.\n");
		fprintf (stderr,"\n\tOPTIONS:\n");
		fprintf (stderr,"\t-C sets the lower left position on the map for raster image [0/0].\n");
		fprintf (stderr,"\t   Optionally, append justification (see pstext for codes)\n");
		GMT_pen_syntax ('F', "draws a frame around the image with the given pen.");
		fprintf (stderr,"\t-Gb and -Gf (1-bit images only) sets the background and foreground color,\n");
		fprintf (stderr,"\t   respectively. Set <color> = - for transparency [Default is black and white]\n");
		fprintf (stderr,"\t-Gt (not for 1-bit images) indicate which color to be made transparent\n");
		fprintf (stderr,"\t   [Default no transparency].\n");
		fprintf (stderr,"\t-I invert 1-bit images (does not affect 8 or 24-bit images).\n");
		GMT_explain_option ('K');
		fprintf (stderr,"\t-M Force color -> monochrome image using GMT_YIQ-transformation.\n");
		fprintf (stderr,"\t-N Replicate image <nx> by <ny> times [Default is no replication]\n");
		GMT_explain_option ('O');
		GMT_explain_option ('P');
		GMT_explain_option ('U');
		GMT_explain_option ('V');
		GMT_explain_option ('X');
		GMT_explain_option ('c');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	/* If not done previously, set foreground to black, background to white */

	if (Ctrl->G.f_rgb[0] == -2) { Ctrl->G.f_rgb[0] = Ctrl->G.f_rgb[1] = Ctrl->G.f_rgb[2] = 0; }
	if (Ctrl->G.b_rgb[0] == -2) { Ctrl->G.b_rgb[0] = Ctrl->G.b_rgb[1] = Ctrl->G.b_rgb[2] = 255; }

	/* Check that the options selected are mutually consistent */

	if (!file[0]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify input raster file\n", GMT_program);
		error++;
	}
	if (Ctrl->W.width <= 0.0 && Ctrl->E.dpi <= 0.0) {
		fprintf (stderr, "%s: Must specify image width (-W) or dpi (-E)\n", GMT_program);
		error++;
	}
	if (Ctrl->N.active && (Ctrl->N.nx < 1 || Ctrl->N.ny < 1)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -N option:  Must specify positive values for replication\n", GMT_program);
		error++;
	}
	if (Ctrl->G.f_rgb[0] < 0 && Ctrl->G.b_rgb[0] < 0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -G option:  Only one of fore/back-ground can be transparent for 1-bit images\n", GMT_program);
		error++;
	}
	if (error) exit (EXIT_FAILURE);

	if (access (file, R_OK)) {
		fprintf (stderr, "%s: Cannot find/open/read file %s\n", GMT_program, file);
		exit (EXIT_FAILURE);
	}

	PS_interpolate = (Ctrl->W.interpolate) ? -1 : +1;
	
	picture = ps_load_image (file, &header);

	if (!picture) {
		fprintf (stderr, "%s: Trouble loading image file!\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	if (Ctrl->M.active) ps_rgb_to_mono (picture, &header);

	/* Add transparent color at beginning, if requested */
	if (Ctrl->G.t_rgb[0] < 0)
		PS_transparent = 1;
	else if (header.depth >= 8) {
		PS_transparent = -1;
		j = header.depth / 8;
		n = j * (header.width * header.height + 1);
		buffer = (unsigned char *) ps_memory (VNULL, (size_t)n, sizeof (unsigned char));
		for (i = 0; i < j; i++) buffer[i] = Ctrl->G.t_rgb[i];
		memcpy (&(buffer[j]), picture, n);
		ps_free ((void *)picture);
		picture = buffer;
	}
	else
		fprintf (stderr, "%s: Can only do transparent color for 8- or 24-bit images. -Gt ignored\n", GMT_program);

	if (!project_info.x_off_supplied && GMT_ps.overlay) GMT_ps.x_origin = 0.0;	/* Since map_setup is not called here */
	if (!project_info.y_off_supplied && GMT_ps.overlay) GMT_ps.y_origin = 0.0;

	GMT_plotinit (argc, argv);

	if (Ctrl->E.dpi > 0.0) Ctrl->W.width = (double) header.width / Ctrl->E.dpi;
	if (Ctrl->W.height == 0.0) Ctrl->W.height = header.height * Ctrl->W.width / header.width;
	justify = GMT_just_decode (Ctrl->C.justify, 12);
	Ctrl->C.x -= 0.5 * ((justify-1)%4) * Ctrl->W.width;
	Ctrl->C.y -= 0.5 * (justify/4) * Ctrl->W.height;

	for (j = 0; j < Ctrl->N.ny; j++) {
		y = Ctrl->C.y + j * Ctrl->W.height;
		if (Ctrl->N.ny > 1 && gmtdefs.verbose) fprintf (stderr, "%s: Replicating image %ld times for row %ld\n", GMT_program, Ctrl->N.nx, j);
		for (i = 0; i < Ctrl->N.nx; i++) {
			x = Ctrl->C.x + i * Ctrl->W.width;
			if (header.depth == 0)
				ps_epsimage (x, y, Ctrl->W.width, Ctrl->W.height, picture, (GMT_LONG)header.length, header.width, header.height, (GMT_LONG)header.xorigin, (GMT_LONG)header.yorigin);
			else if (header.depth == 1)
				/* Invert is opposite from what is expected. This is to match the behaviour of -Gp */
				ps_bitimage (x, y, Ctrl->W.width, Ctrl->W.height, picture, header.width, header.height, !Ctrl->I.active, Ctrl->G.f_rgb, Ctrl->G.b_rgb);
			else
				GMT_color_image (x, y, Ctrl->W.width, Ctrl->W.height, picture, PS_transparent * header.width, header.height, PS_interpolate * header.depth);
		}
	}

 	if (Ctrl->F.active) {
 		GMT_setpen (&Ctrl->F.pen);
 		ps_rect (Ctrl->C.x, Ctrl->C.y, Ctrl->C.x + (Ctrl->N.nx * Ctrl->W.width), Ctrl->C.y + (Ctrl->N.ny * Ctrl->W.height), GMT_no_rgb, TRUE);
 	}

	GMT_plotend ();

	ps_free ((void *)picture);

	Free_psimage_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}
void *New_psimage_Ctrl () {	/* Allocate and initialize a new control structure */
	struct PSIMAGE_CTRL *C;
	
	C = (struct PSIMAGE_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct PSIMAGE_CTRL), "New_psimage_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	GMT_init_pen (&C->F.pen, GMT_PENWIDTH);
	strcpy (C->C.justify, "LB");
	C->G.f_rgb[0] = C->G.b_rgb[0] = C->G.t_rgb[0] = -2;
	C->N.nx = C->N.ny = 1;	
	return ((void *)C);
}

void Free_psimage_Ctrl (struct PSIMAGE_CTRL *C) {	/* Deallocate control structure */
	GMT_free ((void *)C);	
}
