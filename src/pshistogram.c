/*--------------------------------------------------------------------
 *	$Id: pshistogram.c 9923 2012-12-18 20:45:53Z pwessel $
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
 * pshistogram.c -- a program for plotting histograms
 *
 *
 * Author:	Walter H. F. Smith
 * Date:	15 February, 1988
 *
 * Updated to v2.0 5-May-1991 Paul Wessel
 * Updated to v3.0 1-Jan-1995 Paul Wessel
 * Updated to v3.1 13-Jun-1998 Paul Wessel
 * Updated to v3.3 23-Apr-1999 Paul Wessel
 *		     Now uses two passes to get accurate mean and stdev
 * 27-JAN-2005: Added -T option to select column of interest
 * Version:	4.1.2 Free of global variables, and now uses double precision
 *		4.2.1 -C is now -F (-C will still work) and -C<cpt> paints bars
 */

#include "gmt.h"
#include "pslib.h"

struct PSHISTOGRAM_CTRL {
	struct A {	/* -A */
		GMT_LONG active;
	} A;
	struct C {	/* -C<cpt> */
		GMT_LONG active;
		char *file;
	} C;
	struct E {	/* -Eazimuth/elevation */
		GMT_LONG active;
		double azimuth, elevation;
	} E;
	struct F {	/* -F */
		GMT_LONG active;
	} F;
	struct G {	/* -Gfill */
		GMT_LONG active;
		struct GMT_FILL fill;
	} G;
	struct I {	/* -I[o] */
		GMT_LONG active;
		GMT_LONG mode;
	} I;
	struct L {	/* -L<pen> */
		GMT_LONG active;
		struct GMT_PEN pen;
	} L;
	struct Q {	/* -Q */
		GMT_LONG active;
	} Q;
	struct S {	/* -S */
		GMT_LONG active;
	} S;
	struct T {	/* -T<col> */
		GMT_LONG active;
		GMT_LONG col;
	} T;
	struct W {	/* -W<width> */
		GMT_LONG active;
		double inc;
	} W;
	struct Z {	/* -Z<type> */
		GMT_LONG active;
		GMT_LONG mode;
	} Z;
};

#define	PSHISTOGRAM_COUNTS 0
#define PSHISTOGRAM_FREQ_PCT 1
#define PSHISTOGRAM_LOG_COUNTS 2
#define PSHISTOGRAM_LOG_FREQ_PCT 3
#define PSHISTOGRAM_LOG10_COUNTS 4
#define PSHISTOGRAM_LOG10_FREQ_PCT 5

struct PSHISTOGRAM_INFO {	/* Control structure for pshistogram */
	double	yy0, yy1;
	GMT_LONG	*boxh;
	GMT_LONG	n_boxes;
	GMT_LONG	n_counted;
	double	box_width;
	double	west, east, south, north;
	GMT_LONG	center_box, cumulative;
	GMT_LONG	hist_type;
};

int main (int argc, char **argv)
{
	GMT_LONG error = FALSE, automatic = FALSE;

	char buffer[BUFSIZ], format[BUFSIZ];
	
	GMT_LONG i, n_alloc = GMT_CHUNK, n;
	GMT_LONG n_expected_fields, n_fields, n_read;
	
	double	*data = NULL, stats[6];

	double tmp, x_min, x_max, *in = NULL;
	
	FILE *fp_in = NULL;
	
	struct PSHISTOGRAM_INFO F;
	struct PSHISTOGRAM_CTRL *Ctrl = NULL;
	
	GMT_LONG fill_boxes (struct PSHISTOGRAM_INFO *F, double *x, GMT_LONG n);
	GMT_LONG plot_boxes (struct PSHISTOGRAM_INFO *F, GMT_LONG stairs, GMT_LONG flip_to_y, GMT_LONG draw_outline, struct GMT_PEN *pen, struct GMT_FILL *fill, GMT_LONG cpt);
	GMT_LONG get_loc_scl (double *x, GMT_LONG n, double *stats);
	void *New_pshistogram_Ctrl (), Free_pshistogram_Ctrl (struct PSHISTOGRAM_CTRL *C);

	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct PSHISTOGRAM_CTRL *)New_pshistogram_Ctrl ();	/* Allocate and initialize a new control structure */

	memset ((void *)&F, 0, sizeof (struct PSHISTOGRAM_INFO));
	F.hist_type = PSHISTOGRAM_COUNTS;

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {

				/* Common parameters */

				case 'B':
				case 'H':
				case 'J':
				case 'K':
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
				case '\0':
					error += GMT_parse_common_options (argv[i], &F.west, &F.east, &F.south, &F.north);
					break;

				/* Supplemental parameters */

				case 'A':
					Ctrl->A.active = TRUE;
					break;
				case 'C':
					if (argv[i][2]) {	/* Gave an argument */
						Ctrl->C.file = strdup (&argv[i][2]);
						Ctrl->C.active = TRUE;
					}
					else {/* Backwards compatible with old -C (now -F) */
						Ctrl->F.active = TRUE;
						fprintf (stderr, "%s: Warning: -C for centering is obsolete, please use -F\n", GMT_program);
					}
					break;
				case 'E':
					Ctrl->E.active = TRUE;
					sscanf (&argv[i][2], "%lf/%lf", &Ctrl->E.azimuth, &Ctrl->E.elevation);
					break;
				case 'F':
					Ctrl->F.active = TRUE;
					break;
				case 'G':
					Ctrl->G.active = TRUE;
					if (GMT_getfill (&argv[i][2], &Ctrl->G.fill)) {
						GMT_fill_syntax ('G', " ");
						error++;
					}
					break;
				case 'L':		/* Set line attributes */
					Ctrl->L.active = TRUE;
					if (GMT_getpen (&argv[i][2], &Ctrl->L.pen)) {
						GMT_pen_syntax ('L', " ");
						error++;
					}
					break;
				case 'I':
					Ctrl->I.active = TRUE;
					if (argv[i][2] == 'o') Ctrl->I.mode = 1;
					if (argv[i][2] == 'O') Ctrl->I.mode = 2;
					break;
				case 'Q':
					Ctrl->Q.active = TRUE;
					break;
				case 'S':
					Ctrl->S.active = TRUE;
					break;
				case 'T':
					Ctrl->T.active = TRUE;
					Ctrl->T.col = atoi (&argv[i][2]);
					break;
				case 'W':
					Ctrl->W.active = TRUE;
					Ctrl->W.inc = atof (&argv[i][2]);
					break;
				case 'Z':
					Ctrl->Z.active = TRUE;
					Ctrl->Z.mode = atoi (&argv[i][2]);
					break;
				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else {
			if ((fp_in = GMT_fopen(argv[i], GMT_io.r_mode)) == NULL) {
				fprintf (stderr, "%s: Cannot open file %s\n", GMT_program, argv[i]);
				exit (EXIT_FAILURE);
			}
		}
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr,"pshistogram %s - Calculate and plot histograms\n\n", GMT_VERSION);
		fprintf (stderr, "usage: pshistogram [file] %s -W<width> [%s] [-C<cpt>] [-Eaz/el] [-F] [-G<fill>]\n", GMT_Jx_OPT, GMT_B_OPT);
		fprintf (stderr, "\t[%s] [-I[o|O]] [-K] [-L<pen>] [-O] [-P] [-Q] [%s]\n", GMT_Ho_OPT, GMT_Rx_OPT);
		fprintf (stderr, "\t[-S] [-T<col>] [%s] [-V] [%s] [%s]\n", GMT_U_OPT, GMT_X_OPT, GMT_Y_OPT);
		fprintf (stderr, "\t[-Z[type]] [%s] [%s] [%s]\n\n", GMT_c_OPT, GMT_bi_OPT, GMT_f_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);
		fprintf (stderr, "\t-Jx|X for linear projection.  Scale in %s/units (or width in %s)\n", GMT_unit_names[gmtdefs.measure_unit], GMT_unit_names[gmtdefs.measure_unit]);
		fprintf (stderr, "\t    Use / to specify separate x/y scaling.\n");
		fprintf (stderr, "\t    If -JX is used then give axes lengths in %s rather than scales.\n", GMT_unit_names[gmtdefs.measure_unit]);

		fprintf (stderr, "\t-W sets the bin width.\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		GMT_explain_option ('b');
		fprintf (stderr, "\t-A will plot horizontal bars [Default is vertical].\n");
		fprintf (stderr, "\t-C Use cpt-file to assign fill to bars based on the mid x-value.\n");
		fprintf (stderr, "\t-E set azimuth and elevation of viewpoint for 3-D perspective [180/90].\n");
		fprintf (stderr, "\t-F will center the bins.\n");
		GMT_fill_syntax ('G', "Select color/pattern for columns.");
		GMT_explain_option ('H');
		fprintf (stderr, "\t-I will inquire about min/max x and y.  No plotting is done.\n");
		fprintf (stderr, "\t   Append o to output the resulting x, y data.\n");
		fprintf (stderr, "\t   Append O to output all resulting x, y data even with y=0.\n");
		GMT_explain_option ('K');
		GMT_pen_syntax ('L', "specify pen to draw histogram.");
		GMT_explain_option ('O');
		GMT_explain_option ('P');
		fprintf (stderr, "\t-Q plot a cumulative histogram.\n");
		GMT_explain_option ('r');
		fprintf (stderr, "\t   If neither -R nor -I are set, w/e/s/n will be based on input data.\n");
		fprintf (stderr, "\t-S Draws a stairs-step diagram [Default is bar histogram].\n");
		fprintf (stderr, "\t-T Append column to use (0 is first) [0].\n");
		GMT_explain_option ('U');
		GMT_explain_option ('V');
		GMT_explain_option ('X');
		fprintf (stderr, "\t-Z to choose type of vertical axis.  Select from\n");
		fprintf (stderr, "\t   0 - Counts [Default].\n");
		fprintf (stderr, "\t   1 - Frequency percent.\n");
		fprintf (stderr, "\t   2 - Log (1+counts).\n");
		fprintf (stderr, "\t   3 - Log (1+frequency percent).\n");
		fprintf (stderr, "\t   4 - Log10 (1+counts).\n");
		fprintf (stderr, "\t   5 - Log10 (1+frequency percent).\n");
		GMT_explain_option ('c');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf (stderr, "\t   Default is 2 input columns.\n");
		GMT_explain_option ('f');
		GMT_explain_option ('.');

		exit (EXIT_FAILURE);
	}

	if (!Ctrl->I.active && !GMT_IS_LINEAR) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -J option:  Only linear projection supported.\n", GMT_program);
		error++;
	}
	if (Ctrl->E.elevation <= 0.0 || Ctrl->E.elevation > 90.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -E option:  Elevation must be in 0-90 range\n", GMT_program);
		error++;
	}
	if (Ctrl->Z.mode < PSHISTOGRAM_COUNTS || Ctrl->Z.mode > PSHISTOGRAM_LOG10_FREQ_PCT) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -Z option:  histogram type must be in 0-5 range\n", GMT_program);
		error++;
	}
	if (Ctrl->W.inc <= 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -W option:  bin width must be nonzero\n", GMT_program);
		error++;
	}
	/* Now must specify either fill color with -G or outline pen with -L; JLL */
 
	if (!(Ctrl->C.active || Ctrl->I.active || Ctrl->G.active || Ctrl->L.active)) {
		fprintf (stderr, "%s: Must specify either fill (-G) or lookup colors (-C), outline pen attributes (-L), or both.\n", GMT_program);
		error++;
	}
	if (Ctrl->C.active && Ctrl->G.active) {
		fprintf (stderr, "%s: Cannot specify both fill (-G) and lookup colors (-C).\n", GMT_program);
		error++;
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
        if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = 2;
        if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] < 0) {
                fprintf (stderr, "%s: GMT SYNTAX ERROR.  Must specify number of columns in binary input data (-bi)\n", GMT_program);
		error++;
	}
        if (GMT_io.binary[GMT_IN] && Ctrl->T.col >= GMT_io.ncol[GMT_IN]) {
                fprintf (stderr, "%s: GMT SYNTAX ERROR.  Chosen column (%ld) exceeds the columns (%ld) in binary input data (-bi)\n", \
		GMT_program, Ctrl->T.col+1, GMT_io.ncol[GMT_IN]);
		error++;
	}
        if (Ctrl->T.col < 0) {
                fprintf (stderr, "%s: GMT SYNTAX ERROR.  Chosen column (%ld) must be >= 0\n", GMT_program, Ctrl->T.col);
		error++;
	}

	if (error) exit (EXIT_FAILURE);

	if (GMT_io.binary[GMT_IN] && gmtdefs.verbose) {
		char *type[2] = {"double", "single"};
		fprintf (stderr, "%s: Expects %ld-column %s-precision binary data\n", GMT_program, GMT_io.ncol[GMT_IN], type[GMT_io.single_precision[GMT_IN]]);
	}

	F.hist_type = Ctrl->Z.mode;
	F.box_width = Ctrl->W.inc;
	z_project.view_azimuth = Ctrl->E.azimuth;
	z_project.view_elevation = Ctrl->E.elevation;
	F.cumulative = Ctrl->Q.active;
	F.center_box = Ctrl->F.active;
	if (Ctrl->I.active && Ctrl->I.mode == 0) gmtdefs.verbose = TRUE;
	if (!Ctrl->I.active && !project_info.region_supplied) automatic = TRUE;

	if (Ctrl->C.active)  GMT_read_cpt (Ctrl->C.file);

	if (fp_in == NULL) {
		fp_in = GMT_stdin;
#ifdef SET_IO_MODE
		GMT_setmode (GMT_IN);
#endif
	}

	data = (double *) GMT_memory (VNULL, (size_t)n_alloc , sizeof (double), GMT_program);

	n = n_read = 0;
	x_min = DBL_MAX;	x_max = -DBL_MAX;
	n_expected_fields = (GMT_io.ncol[GMT_IN]) ? GMT_io.ncol[GMT_IN] : GMT_MAX_COLUMNS;

	if (GMT_io.io_header[GMT_IN]) for (i = 0; i < GMT_io.n_header_recs; i++) GMT_fgets (buffer, BUFSIZ, fp_in);

	while ((n_fields = GMT_input (fp_in, &n_expected_fields, &in)) >= 0 && !(GMT_io.status & GMT_IO_EOF)) {

		n_read++;

		if (GMT_io.status & GMT_IO_MISMATCH) {
			fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%d) fields near line %ld (skipped)\n", GMT_program, n_fields, 2, n_read);
			continue;
		}
		if (Ctrl->T.col >= n_fields) {
			fprintf (stderr, "%s: Expected at least %ld fields but only found %ld near line %ld (skipped)\n", GMT_program, Ctrl->T.col+1, n_fields, n_read);
			continue;
		}
		data[n] = in[Ctrl->T.col];
		if (!GMT_is_dnan (data[n])) n++;

		if (n == n_alloc) {
			n_alloc <<= 1;
			data = (double *) GMT_memory ((void *)data, (size_t) n_alloc, sizeof (double), GMT_program);
		}
		x_min = MIN (x_min, data[n-1]);
		x_max = MAX (x_max, data[n-1]);
	}

	if (n == 0) {
		fprintf (stderr, "%s: Fatal error, read only 0 points.\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	if (gmtdefs.verbose) fprintf (stderr,"%s: %ld points read\n", GMT_program, n);

	data = (double *) GMT_memory ((void *)data, (size_t)n , sizeof (double), GMT_program);
	if (fp_in != GMT_stdin) GMT_fclose (fp_in);

	get_loc_scl (data, n, stats);

	if (gmtdefs.verbose) {
		sprintf (format, "%%s: Extreme values of the data :\t%s\t%s\n", gmtdefs.d_format, gmtdefs.d_format);
		fprintf (stderr, format, GMT_program, data[0], data[n-1]);
		sprintf (format, "%%s: Locations: L2, L1, LMS; Scales: L2, L1, LMS\t%s\t%s\t%s\t%s\t%s\t%s\n", gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format);
		fprintf (stderr, format, GMT_program, stats[0], stats[1], stats[2], stats[3], stats[4], stats[5]);
	}

	if (F.west == F.east) {	/* Set automatic x range [ and tickmarks] */
		if (frame_info.axis[0].item[0].interval == 0.0) {
			tmp = pow (10.0, floor (d_log10 (x_max-x_min)));
			if (((x_max-x_min) / tmp) < 3.0) tmp *= 0.5;
		}
		else
			tmp = frame_info.axis[0].item[0].interval;
		F.west = floor (x_min / tmp) * tmp;
		F.east = ceil (x_max / tmp) * tmp;
		if (frame_info.axis[0].item[0].interval == 0.0) {
			frame_info.axis[0].item[0].interval = frame_info.axis[0].item[4].interval = tmp;
			frame_info.axis[0].item[0].parent = 0;
			frame_info.axis[0].item[0].active = TRUE;
			frame_info.plot = TRUE;
		}
	}

	if (fill_boxes (&F, data, n)) {
		fprintf (stderr, "%s: Fatal error during box fill.\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	if (gmtdefs.verbose) {
		sprintf (format, "\n%%s: min/max values are :\t%s\t%s\t%s\t%s\n", gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format);
		fprintf (stderr, format, GMT_program, x_min, x_max, F.yy0, F.yy1);
	}

	if (Ctrl->I.active) {	/* Only info requested, quit before plotting */
		if (Ctrl->I.mode) {
			GMT_LONG ibox;
			double xx, yy;
			sprintf (format, "%s\t%s\n", gmtdefs.d_format, gmtdefs.d_format);
			for (ibox = 0; ibox < F.n_boxes; ibox++) {
				if (Ctrl->I.mode == 1 && F.boxh[ibox] == 0) continue;
				xx = F.west + ibox * F.box_width;
				if (F.center_box) xx -= (0.5 * F.box_width);
				if (F.hist_type == PSHISTOGRAM_LOG_COUNTS)
					yy = d_log1p( (double)F.boxh[ibox]);
				else if (F.hist_type == PSHISTOGRAM_LOG10_COUNTS)
					yy = d_log101p( (double)F.boxh[ibox]);
				else if (F.hist_type == PSHISTOGRAM_FREQ_PCT)
					yy = (100.0 * F.boxh[ibox]) / F.n_counted;
				else if (F.hist_type == PSHISTOGRAM_LOG_FREQ_PCT)
					yy = d_log1p( 100.0 * F.boxh[ibox] / F.n_counted );
				else if (F.hist_type == PSHISTOGRAM_LOG10_FREQ_PCT)
					yy = d_log101p( 100.0 * F.boxh[ibox] / F.n_counted );
				else
					yy = (double)F.boxh[ibox];
				fprintf (stdout, format, xx, yy);
			}
		}
		GMT_free ((void *) data);
		GMT_free ((void *) F.boxh);
		exit (EXIT_SUCCESS);
	}

	if (automatic) {	/* Set up s/n based on 'clever' rounding up of the minmax values */
		project_info.region = TRUE;
		F.south = 0.0;
		if (frame_info.axis[1].item[0].interval == 0.0) {
			tmp = pow (10.0, floor (d_log10 (F.yy1)));
			if ((F.yy1 / tmp) < 3.0) tmp *= 0.5;
		}
		else
			tmp = frame_info.axis[1].item[0].interval;
		F.north = ceil (F.yy1 / tmp) * tmp;
		if (frame_info.axis[1].item[0].interval == 0.0) {	/* Tickmarks not set */
			frame_info.axis[1].item[0].interval = frame_info.axis[1].item[4].interval = tmp;
			frame_info.axis[1].item[0].parent = 1;
			frame_info.axis[1].item[0].active = TRUE;
			frame_info.plot = TRUE;
		}
		if (project_info.pars[0] == 0.0) {	/* Must give default xscale */
			project_info.pars[0] = gmtdefs.x_axis_length / (F.east - F.west);
			project_info.projection = GMT_LINEAR;
		}
		if (project_info.pars[1] == 0.0) {	/* Must give default yscale */
			project_info.pars[1] = gmtdefs.y_axis_length / F.north;
			project_info.projection = GMT_LINEAR;
		}
	}

	if (automatic && gmtdefs.verbose) {
		sprintf (format, "%%s: Use w/e/s/n = %s/%s/%s/%s and x-tick/y-tick = %s/%s\n", gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format);
		fprintf (stderr, format, GMT_program, F.west, F.east, F.south, F.north, frame_info.axis[0].item[0].interval, frame_info.axis[1].item[0].interval);
	}

	if (Ctrl->A.active) {
		char buffer[GMT_LONG_TEXT];
		d_swap (frame_info.axis[0].item[0].interval, frame_info.axis[1].item[0].interval);
		d_swap (frame_info.axis[0].item[4].interval, frame_info.axis[1].item[4].interval);
		d_swap (frame_info.axis[0].item[5].interval, frame_info.axis[1].item[5].interval);
		strcpy (buffer, frame_info.axis[0].label);
		strcpy(frame_info.axis[0].label, frame_info.axis[1].label);
		strcpy (frame_info.axis[1].label, buffer);
		GMT_err_fail (GMT_map_setup (F.south, F.north, F.west, F.east), "");
	}
	else
		GMT_err_fail (GMT_map_setup (F.west, F.east, F.south, F.north), "");

	GMT_plotinit (argc, argv);

	if (project_info.three_D) ps_transrotate (-z_project.xmin, -z_project.ymin, 0.0);

	GMT_map_clip_on (GMT_no_rgb, 3);
	if ( plot_boxes (&F, Ctrl->S.active, Ctrl->A.active, Ctrl->L.active, &Ctrl->L.pen, &Ctrl->G.fill, Ctrl->C.active) ) {
		fprintf (stderr, "%s: Fatal error during box plotting.\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	GMT_map_clip_off();

	GMT_map_basemap ();

	if (project_info.three_D) ps_rotatetrans (z_project.xmin, z_project.ymin, 0.0);
	GMT_plotend ();

	GMT_free ((void *) data);
	GMT_free ((void *) F.boxh);

	Free_pshistogram_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

GMT_LONG fill_boxes (struct PSHISTOGRAM_INFO *F, double *data, GMT_LONG n) {

	double	add_half = 0.0;
	GMT_LONG	b0, b1, i, ibox, count_sum;

	F->n_boxes = (GMT_LONG)ceil( ((F->east - F->west) / F->box_width) + 0.5);

	if (F->center_box) {
		F->n_boxes++;
		add_half = 0.5;
	}

	if (F->n_boxes <= 0) return (-1);

	F->boxh = (GMT_LONG *) GMT_memory (VNULL, (size_t)F->n_boxes , sizeof (GMT_LONG), GMT_program);

	F->n_counted = 0;

	/* First fill boxes with counts  */

	for (i = 0; i < n; i++) {
		ibox = (GMT_LONG)floor( ( (data[i] - F->west) / F->box_width) + add_half);
		if (ibox < 0 || ibox >= F->n_boxes) continue;
		F->boxh[ibox] ++;
		F->n_counted++;
	}

	if (F->cumulative) {
		count_sum = 0;
		for (ibox = 0; ibox < F->n_boxes; ibox++) {
			count_sum += F->boxh[ibox];
			F->boxh[ibox] = (GMT_LONG)count_sum;
		}
		b0 = 0;
		b1 = count_sum;
	}
	else {
		b0 = F->n_counted;
		b1 = 0;
		for (ibox = 0; ibox < F->n_boxes; ibox++) {
			if (b0 > F->boxh[ibox]) b0 = F->boxh[ibox];
			if (b1 < F->boxh[ibox]) b1 = F->boxh[ibox];
		}
	}

	/* Now find out what the min max y will be  */

	if (b0) {
		if (F->hist_type == PSHISTOGRAM_LOG_COUNTS)
			F->yy0 = d_log1p((double)b0);
		else if (F->hist_type == PSHISTOGRAM_LOG10_COUNTS)
			F->yy0 = d_log101p( (double)b0);
		else if (F->hist_type == PSHISTOGRAM_FREQ_PCT)
			F->yy0 = (100.0 * b0) / F->n_counted;
		else if (F->hist_type == PSHISTOGRAM_LOG_FREQ_PCT)
			F->yy0 = d_log1p( 100.0 * b0 / F->n_counted );
		else if (F->hist_type == PSHISTOGRAM_LOG10_FREQ_PCT)
			F->yy0 = d_log101p( 100.0 * b0 / F->n_counted );
		else
			F->yy0 = (double)b0;
	}
	else {
		F->yy0 = 0.0;
	}
	if (b1) {
		if (F->hist_type == PSHISTOGRAM_LOG_COUNTS)
			F->yy1 = d_log1p( (double)b1);
		else if (F->hist_type == PSHISTOGRAM_LOG10_COUNTS)
			F->yy1 = d_log101p( (double)b1);
		else if (F->hist_type == PSHISTOGRAM_FREQ_PCT)
			F->yy1 = (100.0 * b1) / F->n_counted;
		else if (F->hist_type == PSHISTOGRAM_LOG_FREQ_PCT)
			F->yy1 = d_log1p( 100.0 * b1 / F->n_counted );
		else if (F->hist_type == PSHISTOGRAM_LOG10_FREQ_PCT)
			F->yy1 = d_log101p( 100.0 * b1 / F->n_counted );
		else
			F->yy1 = (double)b1;
	}
	else {
		F->yy1 = 0.0;
	}
	return (0);
}

GMT_LONG plot_boxes (struct PSHISTOGRAM_INFO *F, GMT_LONG stairs, GMT_LONG flip_to_y, GMT_LONG draw_outline, struct GMT_PEN *pen, struct GMT_FILL *fill, GMT_LONG cpt)
{
	GMT_LONG	i, ibox, first = TRUE, index;
	int rgb[3];
	double	x[4], y[4], xx, yy, xval, *px, *py;
	struct GMT_FILL *f;

	if (draw_outline) GMT_setpen (pen);

	if (flip_to_y) {
		px = y;
		py = x;
	}
	else {
		px = x;
		py = y;
	}

	for (ibox = 0; ibox < F->n_boxes; ibox++) {
		if (stairs || F->boxh[ibox]) {
			x[0] = F->west + ibox * F->box_width;
			if (F->center_box) x[0] -= (0.5 * F->box_width);
			x[1] = x[0] + F->box_width;
			if (x[0] < F->west) x[0] = F->west;
			if (x[1] > F->east) x[1] = F->east;
			xval = 0.5 * (x[0] + x[1]);	/* Used for cpt lookup */
			x[2] = x[1];
			x[3] = x[0];
			y[0] = y[1] = F->south;
			if (F->hist_type == PSHISTOGRAM_LOG_COUNTS)
				y[2] = d_log1p( (double)F->boxh[ibox]);
			else if (F->hist_type == PSHISTOGRAM_LOG10_COUNTS)
				y[2] = d_log101p( (double)F->boxh[ibox]);
			else if (F->hist_type == PSHISTOGRAM_FREQ_PCT)
				y[2] = (100.0 * F->boxh[ibox]) / F->n_counted;
			else if (F->hist_type == PSHISTOGRAM_LOG_FREQ_PCT)
				y[2] = d_log1p( 100.0 * F->boxh[ibox] / F->n_counted );
			else if (F->hist_type == PSHISTOGRAM_LOG10_FREQ_PCT)
				y[2] = d_log101p( 100.0 * F->boxh[ibox] / F->n_counted );
			else
				y[2] = (double)F->boxh[ibox];
			y[3] = y[2];

			for (i = 0; i < 4; i++) {
				GMT_geo_to_xy (px[i], py[i], &xx, &yy);
				if (project_info.three_D) GMT_xyz_to_xy (xx, yy, project_info.z_level, &xx, &yy);
				px[i] = xx;	py[i] = yy;
			}

			if (stairs) {
				if (first) {
					first = FALSE;
					ps_plot (px[0], py[0], PSL_PEN_MOVE);
				}
				ps_plot (px[3], py[3], PSL_PEN_DRAW);
				ps_plot (px[2], py[2], PSL_PEN_DRAW);
			}
			else if (cpt) {
				index = GMT_get_rgb_from_z (xval, rgb);
				if ((index >= 0 && (f = GMT_lut[index].fill)) || (index < 0 && (f = GMT_bfn[index+3].fill)))	/* Pattern */
					GMT_fill (px, py, (GMT_LONG)4, f, draw_outline);
				else
					ps_patch (px, py, (GMT_LONG)4, rgb, draw_outline);
			}
			else
				GMT_fill (px, py, (GMT_LONG)4, fill, draw_outline);
		}
	}

	if (stairs) ps_plot (px[1], py[1], PSL_PEN_DRAW_AND_STROKE);

	return (0);
}

GMT_LONG get_loc_scl (double *data, GMT_LONG n, double *stats)
{
	/* Returns stats[] = L2, L1, LMS location, L2, L1, LMS scale  */

	GMT_LONG	i, j;
	GMT_LONG n_multiples;
	double	dx;

	if (n < 3) return (-1);

	qsort ((void *)data, (size_t)n, sizeof (double), GMT_comp_double_asc);

	/* Get median */
	j = n/2;
	stats[1] = (n%2) ? data[j] : (0.5 * (data[j] + data[j-1]));

	/* Get mode */

	GMT_mode (data, n, j, 0, 0, &n_multiples, &stats[2]);
	if (n_multiples > 0 && gmtdefs.verbose) fprintf (stderr, "%s: WARNING: %ld multiple modes found\n", GMT_program, n_multiples);

	/* Get MAD for L1 */

	GMT_getmad (data, n, stats[1], &stats[4]);

	/* Get LMSscale for mode */

	GMT_getmad (data, n, stats[2], &stats[5]);

	/* Calculate mean and stdev in two passes to minimize risk of overflow */

	stats[0] = stats[3] = 0.0;
	for (i = 0; i < n; i++) stats[0] += data[i];	/* Sum up the data */
	stats[0] /= n;	/* This is the mean value */
	for (i = 0; i < n; i++) {
		dx = data[i] - stats[0];
		stats[3] += (dx * dx);
	}
	stats[3] = sqrt (stats[3] / (n - 1));

	return (0);
}
	
void *New_pshistogram_Ctrl () {	/* Allocate and initialize a new control structure */
	struct PSHISTOGRAM_CTRL *C;
	
	C = (struct PSHISTOGRAM_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct PSHISTOGRAM_CTRL), "New_pshistogram_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	C->E.azimuth = 180.0;
	C->E.elevation = 90.0;
	GMT_init_fill (&C->G.fill, -1, 0, 0);	/* Do not fill is default */
	GMT_init_pen (&C->L.pen, GMT_PENWIDTH);
		
	return ((void *)C);
}

void Free_pshistogram_Ctrl (struct PSHISTOGRAM_CTRL *C) {	/* Deallocate control structure */
	if (C->C.file) free ((void *)C->C.file);	
	GMT_free ((void *)C);	
}
