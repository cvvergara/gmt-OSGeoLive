/*--------------------------------------------------------------------
 *	$Id: psscale.c 10173 2014-01-01 09:52:34Z pwessel $
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
 * psscale draws a grayscale or colorscale either vertically or
 * horizontally.  psscale will interpolate colors if the lower
 * and upper rgb values for an interval are different.  The
 * resolution for this interpolation can be set on the command line.
 * If the scale is to be used with illuminated 2-D and 3-D plots
 * then options for intensities are available.
 *
 * Author:	Paul Wessel
 * Created:	10-MAY-1991
 * Modified:	01-AUG-1998
 *		3.3   13-APR-1999. PW: Fixed bad -X -Y behavior
 *		3.3.2 15-SEP-1999. PW: Fixed bad loop for -Z
 *		3.3.4 21-JAN-2000. PW: Now paint rectangles for discrete cpt files
 *		      08-FEB-2000. PW: Also understand patterns for discrete cpt files
 *		      10-NOV-2003. PW: Enhanced -E to set back or foreground only
 *		      08-MAR-2006. PW: Enhanced -I to set low/high intensity
 *		      23-MAY-2006. PW: Added -Q to handle logarithmic scales
 *		      04-OCT-2007. PW: Can now handle -B1p for log scales
 * Version:	4
 *
 */

#include "gmt.h"
#include "pslib.h"

struct PSSCALE_CTRL {
	struct A {	/* -A */
		GMT_LONG active;
		GMT_LONG mode;
	} A;
	struct C {	/* -C<cptfile> */
		GMT_LONG active;
		char *file;
	} C;
	struct D {	/* -D<xpos/ypos/length/width[h]> */
		GMT_LONG active;
		GMT_LONG horizontal;
		double x, y, width, length;
	} D;
	struct E {	/* -E[b|f][<length>] */
		GMT_LONG active;
		GMT_LONG mode;
		double length;
	} E;
	struct I {	/* -I[<intens>|<min_i>/<max_i>] */
		GMT_LONG active;
		double min, max;
	} I;
	struct M {	/* -M */
		GMT_LONG active;
	} M;
	struct N {	/* -N<dpi> */
		GMT_LONG active;
		GMT_LONG dpi;
	} N;
	struct L {	/* -L[i][<gap>] */
		GMT_LONG active;
		GMT_LONG interval;
		double spacing;
	} L;
	struct Q {	/* -Q */
		GMT_LONG active;
	} Q;
	struct S {	/* -S */
		GMT_LONG active;
	} S;
	struct Z {	/* -Z<zfile> */
		GMT_LONG active;
		char *file;
	} Z;
};

int main(int argc, char **argv)
{
	GMT_LONG error = FALSE, save_unix_time, must_shift_back = FALSE, B_set = FALSE;

	char flag, line[BUFSIZ], txt_a[GMT_LONG_TEXT], txt_b[GMT_LONG_TEXT];
	char txt_c[GMT_LONG_TEXT], txt_d[GMT_LONG_TEXT], text[GMT_LONG_TEXT];

	GMT_LONG i, j, n;

	double max_intens[2], dz, *z_width = NULL, x_origin = 0.0, y_origin = 0.0;
	double start_val, stop_val;

	FILE *fp = NULL;

	struct PSSCALE_CTRL *Ctrl = NULL;

	void GMT_draw_colorbar (double length, double width, double *z_width, GMT_LONG bit_dpi, GMT_LONG flip, GMT_LONG B_set, GMT_LONG equi, 
	   GMT_LONG horizontal, GMT_LONG logscl, GMT_LONG intens, double *max_intens, GMT_LONG skip_lines, GMT_LONG extend, double e_length, double gap, GMT_LONG interval_annot, GMT_LONG monochrome);

	void *New_psscale_Ctrl (), Free_psscale_Ctrl (struct PSSCALE_CTRL *C);
	
	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct PSSCALE_CTRL *)New_psscale_Ctrl ();	/* Allocate and initialize a new control structure */

	/* Check and interpret the command line arguments */

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch(argv[i][1]) {

				/* Common parameters */

				case 'B':
					if (! (argv[i][2] == ':' || argv[i][2] == '/')) B_set = TRUE;
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
				case '\0':
					error += GMT_parse_common_options (argv[i], 0, 0, 0, 0);
					break;

				/* Supplemental parameters */

				case 'A':
					Ctrl->A.active = TRUE;
					if (!argv[i][2]) Ctrl->A.mode |= 3;
					for (j = 2; argv[i][j]; j++) {
						switch (argv[i][j]) {
							case 'a': Ctrl->A.mode |= 1; break;
							case 'l': Ctrl->A.mode |= 2; break;
							case 'c': Ctrl->A.mode |= 4; break;
						}
					}
					break;
				case 'C':
					Ctrl->C.active = TRUE;
					Ctrl->C.file = strdup (&argv[i][2]);
					break;
				case 'D':
					Ctrl->D.active = TRUE;
					n = (GMT_LONG)strlen(argv[i]) - 1;
					flag = argv[i][n];
					if (flag == 'h' || flag == 'H') {
						Ctrl->D.horizontal = TRUE;
						argv[i][n] = 0;	/* Temporarily remove the flag */
					}
					sscanf (&argv[i][2], "%[^/]/%[^/]/%[^/]/%s", txt_a, txt_b, txt_c, txt_d);
					Ctrl->D.x   = GMT_convert_units (txt_a, GMT_INCH);
					Ctrl->D.y   = GMT_convert_units (txt_b, GMT_INCH);
					Ctrl->D.length = GMT_convert_units (txt_c, GMT_INCH);
					Ctrl->D.width  = GMT_convert_units (txt_d, GMT_INCH);
					if (Ctrl->D.horizontal) argv[i][n] = flag;	/* Restore the flag */
					break;
				case 'E':
					Ctrl->E.active = TRUE;
					j = 2;
					if (argv[i][2] == 'b') {
						Ctrl->E.mode = 1;
						j = 3;
					}
					else if (argv[i][2] == 'f') {
						Ctrl->E.mode = 2;
						j = 3;
					}
					else
						Ctrl->E.mode = 3;
					if (argv[i][j]) Ctrl->E.length = GMT_convert_units (&argv[i][j], GMT_INCH);
					break;
				case 'I':
					Ctrl->I.active = TRUE;
					if (argv[i][2]) {
						j = sscanf (&argv[i][2], "%[^/]/%s", txt_a, txt_b);
						if (j == 1) {
							Ctrl->I.max = atof (txt_a);
							Ctrl->I.min = -Ctrl->I.max;
						}
						else {
							Ctrl->I.min = atof (txt_a);
							Ctrl->I.max = atof (txt_b);
						}
					} 
					break;
				case 'M':
					Ctrl->M.active = TRUE;
					break;
				case 'N':
					Ctrl->N.active = TRUE;
					if (argv[i][2]) Ctrl->N.dpi = atoi (&argv[i][2]);
					break;
				case 'L':
					Ctrl->L.active = TRUE;
					j = 2;
					if (argv[i][2] == 'i') Ctrl->L.interval = TRUE, j = 3;
					if (argv[i][j]) Ctrl->L.spacing = GMT_convert_units (&argv[i][j], GMT_INCH);
					break;
				case 'Q':
					Ctrl->Q.active = TRUE;
					break;
				case 'S':
					Ctrl->S.active = TRUE;
					break;
				case 'Z':
					Ctrl->Z.active = TRUE;
					Ctrl->Z.file = strdup (&argv[i][2]);
					break;

				/* Options not recognized */

				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else
			fprintf (stderr, "%s: Warning: Ignoring filename %s\n", GMT_program, argv[i]);
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "%s %s - To create grayscale or colorscale for maps\n\n", GMT_program, GMT_VERSION);
		fprintf(stderr, "usage: psscale -D<xpos/ypos/length/width>[h] [-A[a|l|c]] [-C<cpt_file>] [-E[b|f][<length>]] [%s] [-I[<max_intens>|<low_i>/<high_i>]\n", GMT_B_OPT);
		fprintf(stderr, "\t[-K] [-L[i][<gap>]] [-M] [-N<dpi>] [-O] [-P] [-Q] [-S] [%s] [-V] [%s]\n\t[%s] [-Z<zfile>] [%s]\n\n", GMT_U_OPT, GMT_X_OPT, GMT_Y_OPT, GMT_c_OPT);
		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\t-D set mid-point position and length/width for scale.\n");
		fprintf (stderr, "\t   Give negative length to reverse the scalebar.\n");
		fprintf (stderr, "\t   Append h for horizontal scale.\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-A Place the desired annotations/labels on the other side of the colorscale instead.\n");
		fprintf (stderr, "\t   Append a or l to move only the annotations or labels to the other side.\n");
		fprintf (stderr, "\t   Append c to plot vertical labels as columns.\n");
		fprintf (stderr, "\t-B Set scale annotation interval and label. Use y-label to set unit label.\n");
		fprintf (stderr, "\t   If no annotation interval is set it is taken from the cpt file.\n");
		fprintf (stderr, "\t-C Color palette file. If not set, stdin is read.\n");
		fprintf (stderr, "\t   By default all color changes are annotated (but see -B).  To use a subset,\n");
		fprintf (stderr, "\t   add an extra column to the cpt-file with a L, U, or B\n");
		fprintf (stderr, "\t   to annotate Lower, Upper, or Both color segment boundaries.\n");
		fprintf (stderr, "\t   If a categorical CPT file is given the -Li is set automatically.\n");
		fprintf (stderr, "\t-E add sidebar triangles for back- and foreground colors.\n");
		fprintf (stderr, "\t   Specify b(ackground) or f(oreground) to get one only [Default is both].\n");
		fprintf (stderr, "\t   Optionally, append triangle height [Default is half the barwidth].\n");
		fprintf (stderr, "\t-I add illumination for +-<max_intens> or <low_i> to <high_i> [-1.0/1.0].\n");
		fprintf (stderr, "\t   Alternatively, specify <lower>/<upper> intensity values.\n");
		GMT_explain_option ('K');
		fprintf (stderr, "\t-L For equal-sized color rectangles. -B interval cannot be used.\n");
		fprintf (stderr, "\t   Append i to annotate the interval range instead of lower/upper.\n");
		fprintf (stderr, "\t   If <gap> is appended, we separate each rectangle by <gap> units and center each\n");
		fprintf (stderr, "\t   lower (z0) annotation on the rectangle.  Ignored if not a discrete cpt table.\n");
		fprintf (stderr, "\t   If -I is used then each rectangle will have the illuminated constant color.\n");
		fprintf (stderr, "\t-M force monochrome colorbar using GMT_YIQ transformation.\n");
		fprintf (stderr, "\t-N effective dots-per-inch for color scale [300].\n");
		GMT_explain_option ('O');
		GMT_explain_option ('P');
		fprintf (stderr, "\t-Q Plot colorbar using logarithmic scale and annotate powers of 10 [Default is linear].\n");
		fprintf (stderr, "\t-S Skip drawing color boundary lines on color scale [Default draws lines].\n");
		GMT_explain_option ('U');
		GMT_explain_option ('V');
		GMT_explain_option ('X');
		fprintf (stderr, "\t-Z give colorbar-width (in %s) per color entry.\n", GMT_unit_names[gmtdefs.measure_unit]);
		fprintf (stderr, "\t   By default, width of entry is scaled to color range,\n");
		fprintf (stderr, "\t   i.e., z = 0-100 gives twice the width as z = 100-150.\n");
		GMT_explain_option ('c');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	/* Check that the options selected are mutually consistent */

	if (!Ctrl->D.active) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR: -D is required and must be specified\n", GMT_program);
		error++;
	}
	else {
		if (fabs (Ctrl->D.length) < GMT_SMALL ) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR -D option: scale length must be nonzero\n", GMT_program);
			error++;
		}
		if (Ctrl->D.width <= 0.0) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR -D option: scale width must be positive\n", GMT_program);
			error++;
		}
	}
	if (Ctrl->N.dpi < 1) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -N option: dpi must be > 0\n", GMT_program);
		error++;
	}
	if (Ctrl->L.active && B_set) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -L option: Cannot be used with -B option.\n", GMT_program);
		error++;
	}
	      
	if (error) exit (EXIT_FAILURE);


	if (gmtdefs.verbose) {
		if (Ctrl->C.file)
			fprintf (stderr, "%s: Reading CPT file %s.", GMT_program, Ctrl->C.file);
		else
			fprintf (stderr, "%s: Reading standard input.", GMT_program);
	}

	GMT_read_cpt (Ctrl->C.file);

#ifdef GMT_CPT2	
	if (GMT_categorical) {
		Ctrl->L.active = Ctrl->L.interval = TRUE;
		if (gmtdefs.verbose) fprintf (stderr, "%s: CPT is for categorical data.", GMT_program);
	}
#endif
	
	if (gmtdefs.verbose) fprintf (stderr, "  CPT range from %g to %g\n",
		GMT_lut[0].z_low, GMT_lut[GMT_n_colors-1].z_high);

	if (Ctrl->Q.active) {	/* Take log of all z values */
		for (i = 0; i < GMT_n_colors; i++) {
			if (GMT_lut[i].z_low <= 0.0 || GMT_lut[i].z_high <= 0.0) {
				fprintf (stderr, "%s: GMT SYNTAX ERROR -Q option: All z-values must be positive for logarithmic scale\n", GMT_program);
				exit (EXIT_FAILURE);
			}
			GMT_lut[i].z_low = d_log10 (GMT_lut[i].z_low);
			GMT_lut[i].z_high = d_log10 (GMT_lut[i].z_high);
			GMT_lut[i].i_dz = 1.0 / (GMT_lut[i].z_high - GMT_lut[i].z_low);
		}
	}
	
	if (Ctrl->E.mode && Ctrl->E.length == 0.0) Ctrl->E.length = Ctrl->D.width * 0.5;
	max_intens[0] = Ctrl->I.min;
	max_intens[1] = Ctrl->I.max;

	z_width = (double *) GMT_memory (VNULL, (size_t)GMT_n_colors, sizeof (double), GMT_program);

	if (Ctrl->Z.file && (fp = GMT_fopen (Ctrl->Z.file, "r")) == NULL) {
		fprintf (stderr, "%s: Unable to open file %s\n", GMT_program, Ctrl->Z.file);
		exit (EXIT_FAILURE);
	}
	else if (Ctrl->Z.file) {	/* Opened successfully, now read */
		i = 0;
		while (GMT_fgets (line, BUFSIZ, fp)) {
			if (i == GMT_n_colors) {
				fprintf (stderr, "%s: -Z file %s has more slices than -C file %s!\n", GMT_program, Ctrl->Z.file, Ctrl->C.file);
				exit (EXIT_FAILURE);
			}
			GMT_chop (line);					/* Rid the world of CR/LF */
			if (line[0] == '#' || line[0] == '\0') continue;	/* Skip comments and blank lines */
			n = sscanf (line, "%lf", &z_width[i++]);
			if (n != 1) {
				fprintf (stderr, "%s: Error reading file %s near line %ld\n", GMT_program, Ctrl->Z.file, i);
				exit (EXIT_FAILURE);
			}
		}
		GMT_fclose (fp);
		if (i < GMT_n_colors) {
			fprintf (stderr, "%s: -Z file %s has fewer slices than -C file %s!\n", GMT_program, Ctrl->Z.file, Ctrl->C.file);
			exit (EXIT_FAILURE);
		}
	}
	else if (Ctrl->L.active) {
		dz = fabs (Ctrl->D.length) / GMT_n_colors;
		for (i = 0; i < GMT_n_colors; i++) z_width[i] = dz;
	}
	else {
		for (i = 0; i < GMT_n_colors; i++) z_width[i] = fabs (Ctrl->D.length) * (GMT_lut[i].z_high - GMT_lut[i].z_low) / (GMT_lut[GMT_n_colors-1].z_high - GMT_lut[0].z_low);
	}

/*-----------------start of kludge-------------------------------------------*/

	/* Because psscale uses -D to position things we need to make some
	 * changes so that BoundingBox and others are set ~correctly */

	if (Ctrl->Q.active && B_set) {
		sprintf (text, "X%gil/%gi", Ctrl->D.length, Ctrl->D.width);
		start_val = pow (10.0, GMT_lut[0].z_low);
		stop_val  = pow (10.0, GMT_lut[GMT_n_colors-1].z_high);
	}
	else {
		sprintf (text, "X%gi/%gi", Ctrl->D.length, Ctrl->D.width);
		start_val = GMT_lut[0].z_low;
		stop_val  = GMT_lut[GMT_n_colors-1].z_high;
	}
	GMT_parse_J_option (text);	/* Fake linear projection */
	GMT_err_fail (GMT_map_setup (start_val, stop_val, 0.0, Ctrl->D.width), "");

	/* We must do any origin translation manually in psscale */

	if (! (GMT_ps.x_origin == 0.0 && GMT_ps.y_origin == 0.0)) {	/* Must shift */
		must_shift_back = GMT_ps.absolute;	/* TRUE if absolute */
		x_origin = GMT_ps.x_origin;
		y_origin = GMT_ps.y_origin;
	}

	if (Ctrl->D.horizontal) {
		GMT_ps.x_origin = Ctrl->D.x - 0.5 * fabs (Ctrl->D.length);
		GMT_ps.y_origin = Ctrl->D.y - Ctrl->D.width;
		frame_info.side[1] = frame_info.side[2] = frame_info.side[3] = 0;
	}
	else {
		GMT_ps.x_origin = Ctrl->D.x;
		GMT_ps.y_origin = Ctrl->D.y - 0.5 * fabs (Ctrl->D.length);
		frame_info.side[0] = frame_info.side[2] = frame_info.side[3] = 0;
		d_swap (z_project.xmin, z_project.ymin);
		d_swap (z_project.xmax, z_project.ymax);
	}
	if (frame_info.axis[0].item[0].interval == 0.0) frame_info.plot = TRUE;
	GMT_ps.absolute = TRUE;	/* Ensure absolute coordinates */

	save_unix_time = GMT_ps.unix_time;
	GMT_ps.unix_time = FALSE;
	
/*-----------------end of kludge-------------------------------------------*/

	GMT_plotinit (argc, argv);
	GMT_ps.unix_time = save_unix_time;

	ps_transrotate (x_origin, y_origin, 0.0);

	ps_setpaint (gmtdefs.basemap_frame_rgb);

	if (GMT_ps.unix_time) {	/* Extra shift/unshift since psscale does things differently */
		ps_transrotate (-GMT_ps.x_origin, -GMT_ps.y_origin, 0.0);
		GMT_timestamp (GMT_ps.unix_time_pos[0], GMT_ps.unix_time_pos[1], GMT_ps.unix_time_just, GMT_ps.unix_time_label);
		ps_transrotate (GMT_ps.x_origin, GMT_ps.y_origin, 0.0);
	}

	GMT_draw_colorbar (Ctrl->D.length, Ctrl->D.width, z_width, Ctrl->N.dpi, Ctrl->A.mode, B_set, Ctrl->L.active, Ctrl->D.horizontal, Ctrl->Q.active, Ctrl->I.active, max_intens, Ctrl->S.active, Ctrl->E.mode, Ctrl->E.length, Ctrl->L.spacing, Ctrl->L.interval, Ctrl->M.active);

	if (must_shift_back) ps_transrotate (-x_origin, -y_origin, 0.0);

	GMT_plotend ();

	GMT_ps.absolute = must_shift_back;	/* Restore previous setting */

	GMT_free ((void *)z_width);

	Free_psscale_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void GMT_draw_colorbar (double length, double width, double *z_width, GMT_LONG bit_dpi, GMT_LONG flip, GMT_LONG B_set, GMT_LONG equi, GMT_LONG horizontal, GMT_LONG logscl, GMT_LONG intens, double *max_intens, GMT_LONG skip_lines, GMT_LONG extend, double e_length, double gap, GMT_LONG interval_annot, GMT_LONG monochrome)
{
	GMT_LONG i, ii, id, j, nb, ndec = -1, dec, p_val, depth;
	GMT_LONG nx = 0, ny = 0, nm, barmem, k, justify, l_justify, this_just, use_labels = 0, Label_justify;
	int rgb[3], rrggbb[3];
	GMT_LONG reverse, all = TRUE, use_image, center = FALSE, const_width = TRUE, do_annot;
	char format[GMT_LONG_TEXT], text[GMT_LONG_TEXT], test[GMT_LONG_TEXT], unit[GMT_LONG_TEXT], label[GMT_LONG_TEXT];
	unsigned char *bar = 0, *tmp = 0;
	double off, annot_off, label_off, len, len2, size, x0, x1, dx, xx, dir, y_base, y_annot, y_label;
	double z, xleft, xright, inc_i, inc_j, start_val, stop_val;
	double get_z (double x, double *width, GMT_LONG n);
	double xp[4], yp[4];
	struct GMT_FILL *f;
	void fix_format (char *unit, char *format);

	ps_setfont (gmtdefs.annot_font[0]);

	gmtdefs.annot_offset[0] = fabs (gmtdefs.annot_offset[0]);	/* No 'inside' annotations allowed in colorbar */

	/* Find max decimals needed */

	start_val = GMT_lut[0].z_low;
	stop_val  = GMT_lut[GMT_n_colors-1].z_high;
	if (B_set) {
		if (logscl) {	/* Must set original start/stop values for axis */
			start_val = pow (10.0, GMT_lut[0].z_low);
			stop_val  = pow (10.0, GMT_lut[GMT_n_colors-1].z_high);
		}
		ndec = GMT_get_format (frame_info.axis[0].item[0].interval, frame_info.axis[0].unit, frame_info.axis[0].prefix, format);
	}
	else {
		for (i = 0; i < GMT_n_colors; i++) {
			if (GMT_lut[i].label) use_labels++;
			if (GMT_lut[i].annot & 1) {
				if ((dec = GMT_get_format (GMT_lut[i].z_low, CNULL, CNULL, text)) > ndec) {
					strcpy (format, text);
					ndec = dec;
				}
			}
			if (GMT_lut[i].annot & 2) {
				if ((dec = GMT_get_format (GMT_lut[i].z_high, CNULL, CNULL, text)) > ndec) {
					strcpy (format, text);
					ndec = dec;
				}
			}
			if (GMT_lut[i].annot) all = FALSE;
		}
	}
	if (equi && use_labels == GMT_n_colors)
		all = use_labels = TRUE;	/* Only use optional text labels for equal length scales */
	else
		use_labels = FALSE;

	/* Check if steps in color map have constant width */
	for (i = 1; i < GMT_n_colors && const_width; i++)
		if (fabs(z_width[i] - z_width[0]) > GMT_SMALL) const_width = FALSE;
	
	if (ndec == 0) {	/* Not -B and no decimals are needed */
		strcpy (format, gmtdefs.d_format);
		fix_format (frame_info.axis[0].unit, format);	/* Add units if needed */
	}

	len = gmtdefs.tick_length;	/* +ve means draw on the outside of bar */
	len2 = 0.5 * len;
	xleft = x0 = 0.0;

	reverse = (length < 0.0);
	length = fabs (length);
	xright = length;
	use_image = (!GMT_cpt_pattern && gap <= 0.0 && (equi || const_width || GMT_continuous));

	if ((gap >= 0.0 || interval_annot) && !GMT_continuous) {	/* Want to center annotations for discrete colortable, using lower z0 value */
		center = TRUE;
		if (gap > 0.0) skip_lines = TRUE;
		gap *= 0.5;
		if (interval_annot) {
			sprintf (text, "%s - %s", format, format);
			strcpy (format, text);
		}
	}
	if (gap < 0.0) gap = 0.0;

	if (use_image || intens) {	/* Make bitimage for colorbar using bit_dpi */
		nx = (GMT_continuous) ? irint (length * bit_dpi) : GMT_n_colors;
		ny = (intens) ? irint (width * bit_dpi) : 1;
		nm = nx * ny;
		inc_i = length / nx;
		inc_j = (ny > 1) ? (max_intens[1] - max_intens[0]) / (ny - 1) : 0.0;
		barmem = (monochrome || GMT_gray) ? nm : 3 * nm;
		bar = (unsigned char *) GMT_memory (VNULL, (size_t)barmem, sizeof (char), GMT_program);

		/* Load bar image */

		for (i = 0; i < nx; i++) {
			z = (GMT_continuous) ? get_z ((i+0.5) * inc_i, z_width, GMT_n_colors) : GMT_lut[i].z_low;
			GMT_get_rgb_from_z (z, rrggbb);
			ii = (reverse) ? nx - i - 1 : i;
			for (j = 0; j < ny; j++) {
				for (k = 0; k < 3; k++) rgb[k] = rrggbb[k];
				k = j * nx + ii;
				if (intens) GMT_illuminate (max_intens[1] - j * inc_j, rgb);
				if (GMT_gray)	/* All gray, pick red */
					bar[k] =  (unsigned char) rgb[0];
				else if (monochrome)	/* Convert to gray using the GMT_YIQ transformation */
					bar[k] =  (unsigned char) GMT_YIQ (rgb);
				else {
					k *= 3;
					bar[k++] = (unsigned char) rgb[0];
					bar[k++] = (unsigned char) rgb[1];
					bar[k++] = (unsigned char) rgb[2];
				}
			}
		}
	}

	GMT_setpen (&gmtdefs.frame_pen);

	unit[0] = label[0] = 0;
	/* Defeat the auto-repeat of axis info */
	if (!strcmp (frame_info.axis[0].label, frame_info.axis[1].label)) frame_info.axis[1].label[0] = 0;
	/* Save label and unit, because we are going to switch them off in frame_info and do it ourselves */
	strcpy (label, frame_info.axis[0].label);
	strcpy (unit, frame_info.axis[1].label);
	frame_info.axis[0].label[0] = frame_info.axis[1].label[1] = 0;

	if (flip & 1) {
		justify = l_justify = (horizontal) ? 2 : 7;
		y_base = width;
		dir = 1.0;
	}
	else {
		justify = (horizontal) ? 10 : 7;
		l_justify = (horizontal) ? 10 : 5;
		y_base = 0.0;
		dir = -1.0;
	}
	if (flip & 2) {
		Label_justify = 2;
		y_label = width + gmtdefs.label_offset;
	}
	else {
		Label_justify = 10;
		y_label = -gmtdefs.label_offset;
	}

	/* Current point (0,0) is now at lower left location of bar */

	depth = (monochrome || GMT_gray) ? 8 : 24;
	if (horizontal) {
		if (use_image) {	/* Must plot as image */
			GMT_color_image (0.0, 0.0, length, width, bar, nx, ny, depth);
		}
		else {			/* Plot as rectangles */
			x0 = x1 = 0.0;
			yp[0] = yp[1] = 0.0;	yp[2] = yp[3] = width;
			for (i = 0; i < GMT_n_colors; i++) {
				ii = (reverse) ? GMT_n_colors - i - 1 : i;
				x1 += z_width[ii];
				if ((f = GMT_lut[ii].fill)) {	/* Using pattern fills */
					xp[0] = xp[3] = x0 + gap;	xp[1] = xp[2] = x1 - gap;
					GMT_fill (xp, yp, (GMT_LONG)4, f, center);	/* Outline drawn separately below if desired */
				}
				else if (intens) {
					nb = (GMT_gray || monochrome) ? 1 : 3;
					tmp = (unsigned char *) GMT_memory (VNULL, (size_t)ny*nb, sizeof (char), GMT_program);
					for (j = 0, k = i*nb; j < ny*nb; k+=(nx-1)*nb) {
						tmp[j++] = bar[k++];
						tmp[j++] = bar[k++];
						tmp[j++] = bar[k++];
					}
					GMT_color_image (x0 + gap, 0.0, z_width[ii] - 2*gap, width, tmp, 1, ny, depth);
					GMT_free ((void *)tmp);
					ps_rect (x0+gap, 0.0, x1-gap, width, GMT_no_rgb, center);
				}
				else {
					memcpy ((void *)rgb, (void *)GMT_lut[ii].rgb_low, (size_t)(3 * sizeof (int)));
					if (intens) GMT_illuminate (max_intens[1], rgb);
					if (monochrome) rgb[0] = rgb[1] = rgb[2] = GMT_YIQ (rgb);
					ps_rect (x0+gap, 0.0, x1-gap, width, rgb, center);
				}
				x0 = x1;
			}
		}

		annot_off = ((len > 0.0 && !center) ? len : 0.0) + gmtdefs.annot_offset[0];
		label_off = annot_off + gmtdefs.annot_font_size[0] * GMT_u2u[GMT_PT][GMT_INCH] + gmtdefs.label_offset;
		y_annot = y_base + dir * annot_off;
		if ((flip & 1) == (flip & 2) / 2) y_label = y_base + dir * label_off;

		if (extend & (reverse + 1)) {	/* Add color triangle on left side */
			xp[0] = xp[1] = xleft - gap;	xp[2] = xp[0] - e_length;
			yp[0] = width;	yp[1] = 0.0;	yp[2] = 0.5 * width;
			id = (reverse) ? GMT_FGD : GMT_BGD;
			if ((f = GMT_bfn[id].fill))
				GMT_fill (xp, yp, (GMT_LONG)3, f, TRUE);
			else {
				memcpy ((void *)rgb, (void *)GMT_bfn[id].rgb, (size_t)(3 * sizeof (int)));
				if (monochrome) rgb[0] = rgb[1] = rgb[2] = GMT_YIQ (rgb);
				ps_polygon (xp, yp, (GMT_LONG)3, rgb, TRUE);
			}
		}
		if (extend & (2 - reverse)) {	/* Add color triangle on right side */
			xp[0] = xp[1] = xright + gap;	xp[2] = xp[0] + e_length;
			yp[0] = width;	yp[1] = 0.0;	yp[2] = 0.5 * width;
			id = (reverse) ? GMT_BGD : GMT_FGD;
			if ((f = GMT_bfn[id].fill))
				GMT_fill (xp, yp, (GMT_LONG)3, f, TRUE);
			else {
				memcpy ((void *)rgb, (void *)GMT_bfn[id].rgb, (size_t)(3 * sizeof (int)));
				if (monochrome) rgb[0] = rgb[1] = rgb[2] = GMT_YIQ (rgb);
				ps_polygon (xp, yp, (GMT_LONG)3, rgb, TRUE);
			}
		}

		if (gap == 0.0) {
			ps_segment (xleft, 0.0, xleft + length, 0.0);
			ps_segment (xleft, width, xleft + length, width);
			ps_segment (xleft, 0.0, xleft, width);
			ps_segment (xright, 0.0, xright, width);
		}

		if (B_set) {	/* Used -B */
			GMT_xy_axis (xleft, y_base, length, start_val, stop_val, &frame_info.axis[0], !(flip & 1), frame_info.side[0]-1);
			if ((dx = GMT_get_map_interval (0, GMT_GRID_UPPER)) > 0.0) {
				GMT_setpen (&gmtdefs.grid_pen[0]);
				GMT_linearx_grid (GMT_lut[0].z_low, GMT_lut[GMT_n_colors-1].z_high, 0.0, width, dx);
			}
		}
		else {
			/* First draw gridlines, unless skip_lines is TRUE */

			if (!skip_lines) {
				GMT_setpen (&gmtdefs.grid_pen[0]);
				x1 = xleft;
				for (i = 0; i < GMT_n_colors; i++) {
					xx = (reverse) ? xright - x1 : x1;
					ps_segment (xx, 0.0, xx, width);
					x1 += z_width[i];
				}
				xx = (reverse) ? xleft : xright;
				ps_segment (xx, 0.0, xx, width);
			}

			/* Then annotate and draw tickmarks */

			GMT_setpen (&gmtdefs.tick_pen);
			x1 = xleft;
			if (center) x1 += 0.5 * z_width[0];
			for (i = 0; i < GMT_n_colors; i++) {
				xx = (reverse) ? xright - x1 : x1;
				if (!center) ps_plot (xx, y_base, PSL_PEN_MOVE);
				if (all || (GMT_lut[i].annot & 1)) {	/* Annotate this */
					if (!center) ps_plotr (0.0, dir * len, PSL_PEN_DRAW_AND_STROKE);
					this_just = justify;
					do_annot = TRUE;
					if (use_labels && GMT_lut[i].label) {
						strcpy (text, GMT_lut[i].label);
						this_just = l_justify;
					}
					else if (center && interval_annot)
						sprintf (text, format, GMT_lut[i].z_low, GMT_lut[i].z_high);
					else if (logscl) {
						p_val = irint (GMT_lut[i].z_low);
						if (GMT_IS_ZERO (GMT_lut[i].z_low - (double)p_val))
							sprintf (text, "10@+%ld@+", (GMT_LONG)irint (GMT_lut[i].z_low));
						else
							do_annot = FALSE;
					}
					else
						sprintf (text, format, GMT_lut[i].z_low);
					if (do_annot) ps_text (xx, y_annot, gmtdefs.annot_font_size[0], text, 0.0, -this_just, 0);
				}
				else
					if (!center) ps_plotr (0.0, dir * len2, PSL_PEN_DRAW_AND_STROKE);
				x1 += z_width[i];
			}
			if (!center && !use_labels) {
				xx = (reverse) ? xleft : xright;
				ps_plot (xx, y_base, PSL_PEN_MOVE);
				i = GMT_n_colors-1;
				if (all || (GMT_lut[i].annot & 2)) {
					ps_plotr (0.0, dir * len, PSL_PEN_DRAW_AND_STROKE);
					this_just = justify;
					do_annot = TRUE;
					if (logscl) {
						p_val = irint (GMT_lut[i].z_high);
						if (GMT_IS_ZERO (GMT_lut[i].z_high - (double)p_val))
							sprintf (text, "10@+%ld@+", p_val);
						else
							do_annot = FALSE;
					}
					else
						sprintf (text, format, GMT_lut[i].z_high);
					if (do_annot) ps_text (xx, y_annot, gmtdefs.annot_font_size[0], text, 0.0, -this_just, 0);
				}
				else
					ps_plotr (0.0, dir * len2, PSL_PEN_DRAW_AND_STROKE);
			}
		}
		if (label[0]) {	/* Add label */
			ps_setfont (gmtdefs.label_font);
			ps_text (xleft + 0.5 * length, y_label, gmtdefs.label_font_size, label, 0.0, Label_justify, 0);
		}
		if (unit[0]) {	/* Add unit label */
			ps_setfont (gmtdefs.annot_font[0]);
			ps_text (xright + e_length + gmtdefs.annot_offset[0], 0.5 * width, gmtdefs.annot_font_size[0], unit, 0.0, 5, 0);
		}
	}
	else {	/* Vertical scale */
		ps_transrotate (width, 0.0, 90.0);
		if (use_image) {	/* Must plot with image */
			GMT_color_image (0.0, 0.0, length, width, bar, nx, ny, depth);
		}
		else {			/* Plot as rectangles */
			x0 = x1 = 0.0;
			yp[0] = yp[1] = 0.0;	yp[2] = yp[3] = width;
			for (i = 0; i < GMT_n_colors; i++) {
				ii = (reverse) ? GMT_n_colors - i - 1 : i;
				x1 += z_width[ii];
				if ((f = GMT_lut[ii].fill)) {	/* Using pattern fills */
					ps_transrotate (x0 + gap, 0.0, -90.0);

					xp[0] = xp[3] = -width;	xp[1] = xp[2] = 0.0;
					yp[0] = yp[1] = 0.0;	yp[2] = yp[3] = x1 - x0 - 2.0 * gap;
					GMT_fill (xp, yp, (GMT_LONG)4, f, center);	/* Outline drawn separately below if desired */
					ps_rotatetrans (-(x0 + gap), 0.0, 90.0);
				}
				else if (intens) {
					nb = (GMT_gray || monochrome) ? 1 : 3;
					tmp = (unsigned char *) GMT_memory (VNULL, (size_t)ny*nb, sizeof (char), GMT_program);
					for (j = 0, k = i*nb; j < ny*nb; k+=(nx-1)*nb) {
						tmp[j++] = bar[k++];
						tmp[j++] = bar[k++];
						tmp[j++] = bar[k++];
					}
					GMT_color_image (x0 + gap, 0.0, z_width[ii] - 2*gap, width, tmp, 1, ny, depth);
					GMT_free ((void *)tmp);
					ps_rect (x0 + gap, 0.0, x1 - gap, width, GMT_no_rgb, center);
				}
				else {
					memcpy ((void *)rgb, (void *)GMT_lut[ii].rgb_low, (size_t)(3 * sizeof (int)));
					if (intens) GMT_illuminate (max_intens[1], rgb);
					if (monochrome) rgb[0] = rgb[1] = rgb[2] = GMT_YIQ (rgb);
					ps_rect (x0 + gap, 0.0, x1 - gap, width, rgb, center);
				}
				x0 = x1;
			}
		}

		if (center && interval_annot) {
			sprintf (text, "%ld - %ld", (GMT_LONG) (floor (GMT_lut[0].z_low)), (GMT_LONG) (ceil (GMT_lut[0].z_high)));
			sprintf (test, "%ld - %ld", (GMT_LONG) (floor (GMT_lut[GMT_n_colors-1].z_low)), (GMT_LONG) (ceil (GMT_lut[GMT_n_colors-1].z_high)));
			off = ((MAX ((GMT_LONG)strlen (text), (GMT_LONG)strlen (test)) + 2*ndec) * GMT_DEC_SIZE - 0.4 + 
				((ndec > 0) ? 2*GMT_PER_SIZE : 0.0))
				* gmtdefs.annot_font_size[0] * GMT_u2u[GMT_PT][GMT_INCH];
		}
		else {
			sprintf (text, "%ld", (GMT_LONG) (floor (GMT_lut[0].z_low)));
			sprintf (test, "%ld", (GMT_LONG) (ceil (center ? GMT_lut[GMT_n_colors-1].z_low : GMT_lut[GMT_n_colors-1].z_high)));
			off = ((MAX ((GMT_LONG)strlen (text), (GMT_LONG)strlen (test)) + ndec) * GMT_DEC_SIZE +
				((ndec > 0) ? GMT_PER_SIZE : 0.0))
				* gmtdefs.annot_font_size[0] * GMT_u2u[GMT_PT][GMT_INCH];
		}

		annot_off = ((len > 0.0 && !center) ? len : 0.0) + gmtdefs.annot_offset[0] + off;
		label_off = annot_off + gmtdefs.label_offset;
		if (use_labels || (flip & 1) || logscl) annot_off -= off;
		y_annot = y_base + dir * annot_off;
		if ((flip & 1) == (flip & 2) / 2) y_label = y_base + dir * label_off;

		if (extend & (reverse + 1)) {	/* Add color triangle at bottom */
			xp[0] = xp[1] = xleft - gap;	xp[2] = xp[0] - e_length;
			yp[0] = width;	yp[1] = 0.0;	yp[2] = 0.5 * width;
			id = (reverse) ? GMT_FGD : GMT_BGD;
			if ((f = GMT_bfn[id].fill))
				GMT_fill (xp, yp, (GMT_LONG)3, f, TRUE);
			else {
				memcpy ((void *)rgb, (void *)GMT_bfn[id].rgb, (size_t)(3 * sizeof (int)));
				if (monochrome) rgb[0] = rgb[1] = rgb[2] = GMT_YIQ (rgb);
				ps_polygon (xp, yp, (GMT_LONG)3, rgb, TRUE);
			}
		}
		if (extend & (2 - reverse)) {	/* Add color triangle at top */
			xp[0] = xp[1] = xright + gap;	xp[2] = xp[0] + e_length;
			yp[0] = width;	yp[1] = 0.0;	yp[2] = 0.5 * width;
			id = (reverse) ? GMT_BGD : GMT_FGD;
			if ((f = GMT_bfn[id].fill))
				GMT_fill (xp, yp, (GMT_LONG)3, f, TRUE);
			else {
				memcpy ((void *)rgb, (void *)GMT_bfn[id].rgb, (size_t)(3 * sizeof (int)));
				if (monochrome) rgb[0] = rgb[1] = rgb[2] = GMT_YIQ (rgb);
				ps_polygon (xp, yp, (GMT_LONG)3, rgb, TRUE);
			}
		}
		if (gap == 0.0) {
			ps_segment (xleft, 0.0, xleft + length, 0.0);
			ps_segment (xleft, width, xleft + length, width);
			ps_segment (xleft, 0.0, xleft, width);
			ps_segment (xright, 0.0, xright, width);
		}
		if (B_set) {	/* Used -B. Must kludge by copying x-axis and scaling to y since we must use GMT_xy_axis to draw a rotated x-axis... */
			PFL tmp;
			if ((dx = GMT_get_map_interval (0, GMT_GRID_UPPER)) > 0.0) {
				GMT_setpen (&gmtdefs.grid_pen[0]);
				GMT_linearx_grid (GMT_lut[0].z_low, GMT_lut[GMT_n_colors-1].z_high, 0.0, width, dx);
			}
			ps_transrotate (0.0, 0.0, -90.0);
			memcpy ((void *)&frame_info.axis[1], (void *)&frame_info.axis[0], sizeof (struct GMT_PLOT_AXIS));
			d_swap (project_info.x_scale, project_info.y_scale);
			d_swap (project_info.x0, project_info.y0);
			l_swap (project_info.xyz_projection[0], project_info.xyz_projection[1]);
			tmp = GMT_x_forward; GMT_y_forward = GMT_x_forward; GMT_x_forward = tmp;
			for (i = 0; i < 5; i++) frame_info.axis[1].item[i].parent = 1;
			GMT_xy_axis (-y_base, 0.0, length, start_val, stop_val, &frame_info.axis[1], (flip & 1), frame_info.side[0]-1);
			ps_rotatetrans (0.0, 0.0, 90.0);
		}
		else {
			if (!skip_lines) {	/* First draw gridlines */
				GMT_setpen (&gmtdefs.grid_pen[0]);
				x1 = xleft;
				for (i = 0; i < GMT_n_colors; i++) {
					xx = (reverse) ? xright - x1 : x1;
					ps_segment (xx, 0.0, xx, width);
					x1 += z_width[i];
				}
				xx = (reverse) ? xleft : xright;
				ps_segment (xx, 0.0, xx, width);
			}

			/* Then annotate and draw tickmarks */

			GMT_setpen (&gmtdefs.tick_pen);
			x1 = xleft;
			if (center) x1 += 0.5 * z_width[0];
			for (i = 0; i < GMT_n_colors; i++) {
				xx = (reverse) ? xright - x1 : x1;
				if (!center) ps_plot (xx, y_base, PSL_PEN_MOVE);
				if (all || (GMT_lut[i].annot & 1)) {
					if (!center) ps_plotr (0.0, dir * len, PSL_PEN_DRAW_AND_STROKE);
					this_just = justify;
					do_annot = TRUE;
					if (use_labels && GMT_lut[i].label) {
						strcpy (text, GMT_lut[i].label);
						this_just = l_justify;
					}
					else if (center && interval_annot)
						sprintf (text, format, GMT_lut[i].z_low, GMT_lut[i].z_high);
					else if (logscl) {
						p_val = irint (GMT_lut[i].z_low);
						if (GMT_IS_ZERO (GMT_lut[i].z_low - (double)p_val))
							sprintf (text, "10@+%ld@+", p_val);
						else
							do_annot = FALSE;
						this_just = l_justify;
					}
					else
						sprintf (text, format, GMT_lut[i].z_low);
					if (do_annot) ps_text (xx, y_annot, gmtdefs.annot_font_size[0], text, -90.0, -this_just, 0);
				}
				else
					if (!center) ps_plotr (0.0, dir * len2, PSL_PEN_DRAW_AND_STROKE);
				x1 += z_width[i];
			}
			if (!center && !use_labels) {
				xx = (reverse) ? xleft : xright;
				ps_plot (xx, y_base, PSL_PEN_MOVE);
				i = GMT_n_colors-1;
				if (all || (GMT_lut[i].annot & 2)) {
					ps_plotr (0.0, dir * len, PSL_PEN_DRAW_AND_STROKE);
					this_just = justify;
					do_annot = TRUE;
					if (logscl) {
						p_val = irint (GMT_lut[i].z_high);
						if (GMT_IS_ZERO (GMT_lut[i].z_high - (double)p_val))
							sprintf (text, "10@+%ld@+", p_val);
						else
							do_annot = FALSE;
						this_just = l_justify;
					}
					else
						sprintf (text, format, GMT_lut[i].z_high);
					if (do_annot) ps_text (xx, y_annot, gmtdefs.annot_font_size[0], text, -90.0, -this_just, 0);
				}
				else
					ps_plotr (0.0, dir * len2, PSL_PEN_DRAW_AND_STROKE);
			}
		}

		if (label[0]) {	/* Add label */
			ps_setfont (gmtdefs.label_font);
			if (strchr (label, '@') || strchr (label, '(') || !(flip & 4)) { /* Must set text along-side color bar */
				ps_text (xleft + 0.5 * length, y_label, gmtdefs.label_font_size, label, 0.0, Label_justify, 0);
			}
			else {	/* Set label text as column (AARRGGHH) */
				y_label += 0.5 * ((flip & 2) - 1) * gmtdefs.label_font_size * GMT_u2u[GMT_PT][GMT_INCH];
				size = 0.9 * gmtdefs.label_font_size * GMT_u2u[GMT_PT][GMT_INCH];
				x0 = 0.5 * (length + ((GMT_LONG)strlen (label) -1) * size);
				text[1] = 0;
				for (i = 0; i < (GMT_LONG)strlen (label); i++) {
					x1 = x0 - i * size;
					text[0] = label[i];
					ps_text (x1, y_label, gmtdefs.label_font_size, text, -90.0, 6, 0);
				}
			}
		}
		if (unit[0]) {	/* Add unit label */
			ps_setfont (gmtdefs.annot_font[0]);
			ps_text (xright + gmtdefs.annot_offset[0] + e_length, 0.5 * width, gmtdefs.annot_font_size[0], unit, -90.0, 2, 0);
		}
		ps_rotatetrans (-width, 0.0, -90.0);

	}
	if (use_image || intens) GMT_free ((void *)bar);
}

double get_z (double x, double *width, GMT_LONG n)
{
	GMT_LONG i = 0;
	double tmp;

	tmp = width[0];
	while (i < n && x > tmp) tmp += width[++i];
	if (i == n) return (GMT_lut[GMT_n_colors-1].z_high);
	return (GMT_lut[i].z_low + (x - tmp + width[i]) * (GMT_lut[i].z_high - GMT_lut[i].z_low) / width[i]);
}

void fix_format (char *unit, char *format)
{
	GMT_LONG i, j;
	char text[GMT_TEXT_LEN], new_format[BUFSIZ];

	/* Check if annotation units should be added */

	if (unit && unit[0]) {	/* Must append the unit string */
		if (!strchr (unit, '%'))	/* No percent signs */
			strncpy (text, unit, (size_t)GMT_TEXT_LEN);
		else {
			for (i = j = 0; i < (GMT_LONG)strlen (unit); i++) {
				text[j++] = unit[i];
				if (unit[i] == '%') text[j++] = unit[i];
			}
			text[j] = 0;
		}
		if (text[0] == '-')	/* No space between annotation and unit */
			sprintf (new_format, "%s%s", format, &text[1]);
		else		/* 1 space between annotation and unit */
			sprintf (new_format, "%s %s", format, text);
		strcpy (format, new_format);
	}
}

void *New_psscale_Ctrl () {	/* Allocate and initialize a new control structure */
	struct PSSCALE_CTRL *C;
	
	C = (struct PSSCALE_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct PSSCALE_CTRL), "New_psscale_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	C->N.dpi = 300;
	C->I.min = -1.0;
	C->I.max = +1.0;
	C->L.spacing = -1.0;
	return ((void *)C);
}

void Free_psscale_Ctrl (struct PSSCALE_CTRL *C) {	/* Deallocate control structure */
	if (C->C.file) free ((void *)C->C.file);
	if (C->Z.file) free ((void *)C->Z.file);
	GMT_free ((void *)C);	
}
