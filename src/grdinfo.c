/*--------------------------------------------------------------------
 *	$Id: grdinfo.c 9923 2012-12-18 20:45:53Z pwessel $
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
 * grdinfo reads one or more grid file and [optionally] prints out various
 * statistics like mean/standard deviation and median/scale
 *
 * Author:	Paul Wessel
 * Date:	20-SEP-2001
 * Version:	4
 *		12-APR-2006: SHould be 64-bit ready
 */

#include "gmt.h"
#include "grdinfo.h"	/* Array with gridtype names */

struct GRDINFO_CTRL {
	struct C {	/* -C */
		GMT_LONG active;
	} C;
	struct F {	/* -F */
		GMT_LONG active;
	} F;
	struct I {	/* -Idx[/dy] */
		GMT_LONG active;
		GMT_LONG status;
		double xinc, yinc;
	} I;
	struct M {	/* -M */
		GMT_LONG active;
	} M;
	struct L {	/* -L[1|2] */
		GMT_LONG active;
		GMT_LONG norm;
	} L;
	struct T {	/* -T<dz> */
		GMT_LONG active;
		double inc;
	} T;
};

int main (int argc, char **argv)
{
	GMT_LONG error = FALSE, slow = FALSE;

	GMT_LONG nfiles = 0, k, i, n_grds = 0, col, row;
	GMT_LONG nm = 0, n_nan = 0, n = 0, ij, ij_min, ij_max;

	float *a = NULL;

	double x_min = DBL_MAX, y_min = DBL_MAX, z_min = DBL_MAX, x_max = -DBL_MAX, y_max = -DBL_MAX, z_max = -DBL_MAX;
	double global_xmin, global_xmax, global_ymin, global_ymax, global_zmin, global_zmax;
	double mean = 0.0, median = 0.0, sum2 = 0.0, stdev = 0.0, scale = 0.0, rms = 0.0, x;

	char format[BUFSIZ], text[GMT_TEXT_LEN], buffer[BUFSIZ];
	char *type[2] = { "Gridline", "Pixel"};

	struct GRD_HEADER grd;
	struct GRDINFO_CTRL *Ctrl = NULL;

	void *New_grdinfo_Ctrl (), Free_grdinfo_Ctrl (struct GRDINFO_CTRL *C);

	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct GRDINFO_CTRL *)New_grdinfo_Ctrl ();	/* Allocate and initialize a new control structure */

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				/* Common parameters */

				case 'V':
				case 'f':
				case '\0':
					error += GMT_parse_common_options (argv[i], 0, 0, 0, 0);
					break;

				/* Supplemental parameters */

				case 'C':
					Ctrl->C.active = TRUE;
					break;
				case 'D':	/* Left for backwards compatibility, use -f instead */
					GMT_io.in_col_type[0] = GMT_io.out_col_type[0] = GMT_IS_LON;
					GMT_io.in_col_type[1] = GMT_io.out_col_type[1] = GMT_IS_LAT;
					break;
				case 'F':
					Ctrl->F.active = TRUE;
					break;
				case 'I':
					Ctrl->I.active = TRUE;
					if (argv[i][2] == '\0')	/* No argus given, we want to output the -I string */
						Ctrl->I.status = 0;
					else if (argv[i][2] == '-' && argv[i][3] == '\0')	/* Dash given, we want to output the actual -R string */
						Ctrl->I.status = 1;
					else {	/* Report -R to nearest given multiple increment */
						Ctrl->I.status = 2;
						if (GMT_getinc (&argv[i][2], &Ctrl->I.xinc, &Ctrl->I.yinc)) {
							GMT_inc_syntax ('I', 1);
							error = TRUE;
						}
					}
					break;
				case 'M':
					Ctrl->M.active = TRUE;
					break;
				case 'L':
					Ctrl->L.active = TRUE;
					if (argv[i][2] == 0 || argv[i][2] == '2')
						Ctrl->L.norm |= 2;
					else if (argv[i][2] == '1')
						Ctrl->L.norm |= 1;
					break;
				case 'T':
					Ctrl->T.active = TRUE;
					Ctrl->T.inc = atof (&argv[i][2]);
					break;
				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else
			nfiles ++;
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "grdinfo %s - Extract information from netCDF grid files\n\n", GMT_VERSION);
		fprintf (stderr, "usage: grdinfo <grdfiles> [-C] [-F] [-I[<dx>[/<dy>]]] [-L[0|1|2]] [-M]\n");
		fprintf (stderr, "	[-T<dz>] [%s]\n", GMT_f_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\t<grdfiles> may be one or more netCDF grid files.\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-C formats report in fields on a single line using the format\n");
		fprintf (stderr, "\t   file w e s n z0 z1 dx dy nx ny [x0 y0 x1 y1] [med scale] [mean std rms] [n_nan]\n");
		fprintf (stderr, "\t   (-M gives [x0 y0 x1 y1] and [n_nan]; -L1 gives [med scale]; -L2 gives [mean std rms]).\n");
		fprintf (stderr, "\t-F reports domain in world mapping format [Default is generic].\n");
		fprintf (stderr, "\t-I returns textstring -Rw/e/s/n to nearest multiple of dx/dy.\n");
		fprintf (stderr, "\t   If -C is set then rounding off will occur but no -R string is issued.\n");
		fprintf (stderr, "\t   If no argument is given then the -I<xinc>/<yinc> string is issued.\n");
		fprintf (stderr, "\t   If -I- is given then the grid's -R string is issued.\n");
		fprintf (stderr, "\t-L0 reports range of data by actually reading them (not from header).\n");
		fprintf (stderr, "\t-L1 reports median and L1-scale of data set.\n");
		fprintf (stderr, "\t-L[2] reports mean, standard deviation, and rms of data set.\n");
		fprintf (stderr, "\t-M searches for the global min and max locations (x0,y0) and (x1,y1).\n");
		fprintf (stderr, "\t-T given increment dz, return global -Tzmin/zmax/dz in multiples of dz.\n");
		GMT_explain_option ('V');
		GMT_explain_option ('f');
		exit (EXIT_FAILURE);
	}

	if (nfiles == 0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR: Must specify one or more input files\n", GMT_program);
		error++;
	}
	if (Ctrl->T.active && Ctrl->T.inc <= 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -T: Must specify a positive increment\n", GMT_program);
		error++;
	}
	if (Ctrl->I.active && Ctrl->I.status == 2 && (Ctrl->I.xinc <= 0.0 || Ctrl->I.yinc <= 0.0)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -I: Must specify a positive increment(s)\n", GMT_program);
		error++;
	}
	if ((Ctrl->I.active || Ctrl->T.active) && slow) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -M: Not compatible with -I or -T\n", GMT_program);
		error++;
	}
	if (Ctrl->T.active && Ctrl->I.active) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR: Only one of -I -T can be specified\n", GMT_program);
		error++;
	}

	if (error) exit (EXIT_FAILURE);

	a = (float *) GMT_memory (VNULL, (size_t)1, sizeof (float), GMT_program);

	global_xmin = global_ymin = global_zmin = +DBL_MAX;
	global_xmax = global_ymax = global_zmax = -DBL_MAX;

	slow = (Ctrl->M.active || Ctrl->L.active);

	for (k = 1; k < argc; k++) {	/* Loop over arguments, skip options */

		if (argv[k][0] == '-') continue;

		GMT_grd_init (&grd, argc, argv, FALSE);

		if (GMT_err_pass (GMT_read_grd_info (argv[k], &grd), argv[k])) continue;

		if (gmtdefs.verbose) fprintf (stderr, "%s: Processing file %s\n", GMT_program, grd.name);

		n_grds++;

		if (slow) {	/* Must determine the location of global min and max values */

		 	nm = GMT_get_nm (grd.nx, grd.ny);
			a = (float *) GMT_memory ((void *)a, nm, sizeof (float), GMT_program);
			if (GMT_err_pass (GMT_read_grd (argv[k], &grd, a, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), grd.name)) continue;

			z_min = DBL_MAX;	z_max = -DBL_MAX;
			mean = median = sum2 = 0.0;
			ij_min = ij_max = n = 0;
			for (ij = 0; ij < nm; ij++) {
				if (GMT_is_fnan (a[ij])) continue;
				if (a[ij] < z_min) {
					z_min = a[ij];
					ij_min = ij;
				}
				if (a[ij] > z_max) {
					z_max = a[ij];
					ij_max = ij;
				}
				/* Use Welford (1962) algorithm to compute mean and corrected sum of squares */
				n++;
				x = a[ij] - mean;
				mean += x / n;
				sum2 += x * (a[ij] - mean);
			}

			n_nan = nm - n;
			col = (GMT_LONG)(ij_min % grd.nx);
			row = (GMT_LONG)(ij_min / grd.nx);
			x_min = GMT_i_to_x (col, grd.x_min, grd.x_max, grd.x_inc, grd.xy_off, grd.nx);
			y_min = GMT_j_to_y (row, grd.y_min, grd.y_max, grd.y_inc, grd.xy_off, grd.ny);
			col = (GMT_LONG)(ij_max % grd.nx);
			row = (GMT_LONG)(ij_max / grd.nx);
			x_max = GMT_i_to_x (col, grd.x_min, grd.x_max, grd.x_inc, grd.xy_off, grd.nx);
			y_max = GMT_j_to_y (row, grd.y_min, grd.y_max, grd.y_inc, grd.xy_off, grd.ny);
		}

		if (Ctrl->L.norm & 1) {	/* Calculate the median and L1 scale */
			qsort ((void *)a, nm, sizeof (float), GMT_comp_float_asc);
			median = (n%2) ? a[n/2] : 0.5*(a[n/2-1] + a[n/2]);
			for (ij = 0; ij < n; ij++) a[ij] = (float)fabs (a[ij] - median);
			qsort ((void *)a, n, sizeof (float), GMT_comp_float_asc);
			scale = (n%2) ? 1.4826 * a[n/2] : 0.7413 * (a[n/2-1] + a[n/2]);
		}
		if (Ctrl->L.norm & 2) {	/* Calculate the mean, standard deviation, and rms */
			x = (double)n;
			stdev = (n > 1) ? sqrt(sum2 / (x-1)) : GMT_d_NaN;
			rms = (n > 0) ? sqrt (sum2 / x + mean * mean) : GMT_d_NaN;
			mean = (n > 0) ? mean : GMT_d_NaN;
		}

		/* OK, time to report results */

		if (Ctrl->I.active && Ctrl->I.status == 1) {
			GMT_fputs ("-R", GMT_stdout);
			GMT_ascii_output_one (GMT_stdout, grd.x_min, 0);	GMT_fputs ("/", GMT_stdout);
			GMT_ascii_output_one (GMT_stdout, grd.x_max, 0);	GMT_fputs ("/", GMT_stdout);
			GMT_ascii_output_one (GMT_stdout, grd.y_min, 1);	GMT_fputs ("/", GMT_stdout);
			GMT_ascii_output_one (GMT_stdout, grd.y_max, 1);	GMT_fputs ("\n", GMT_stdout);
		} else if (Ctrl->I.active && Ctrl->I.status == 0) {
			GMT_fputs ("-I", GMT_stdout);
			GMT_ascii_output_one (GMT_stdout, grd.x_inc, 2);	GMT_fputs ("/", GMT_stdout);
			GMT_ascii_output_one (GMT_stdout, grd.y_inc, 2);
			GMT_fputs ("\n", GMT_stdout);
		} else if (Ctrl->C.active && !Ctrl->I.active) {
			sprintf (buffer, "%s\t", grd.name);			GMT_fputs (buffer, GMT_stdout);
			GMT_ascii_output_one (GMT_stdout, grd.x_min, 0);	GMT_fputs ("\t", GMT_stdout);
			GMT_ascii_output_one (GMT_stdout, grd.x_max, 0);	GMT_fputs ("\t", GMT_stdout);
			GMT_ascii_output_one (GMT_stdout, grd.y_min, 1);	GMT_fputs ("\t", GMT_stdout);
			GMT_ascii_output_one (GMT_stdout, grd.y_max, 1);	GMT_fputs ("\t", GMT_stdout);
			GMT_ascii_output_one (GMT_stdout, grd.z_min, 2);	GMT_fputs ("\t", GMT_stdout);
			GMT_ascii_output_one (GMT_stdout, grd.z_max, 2);	GMT_fputs ("\t", GMT_stdout);
			GMT_ascii_format_one (text, grd.x_inc, GMT_io.out_col_type[0]);
			if (isalpha ((int)text[strlen(text)-1])) text[strlen(text)-1] = '\0';	/* Chop of trailing WESN flag here */
			GMT_fputs (text, GMT_stdout);	GMT_fputs ("\t", GMT_stdout);
			GMT_ascii_format_one (text, grd.y_inc, GMT_io.out_col_type[1]);
			if (isalpha ((int)text[strlen(text)-1])) text[strlen(text)-1] = '\0';	/* Chop of trailing WESN flag here */
			GMT_fputs (text, GMT_stdout);	GMT_fputs ("\t", GMT_stdout);
			GMT_ascii_output_one (GMT_stdout, (double)grd.nx, 2);	GMT_fputs ("\t", GMT_stdout);	GMT_ascii_output_one (GMT_stdout, (double)grd.ny, 2);

			if (Ctrl->M.active) {
				GMT_fputs ("\t", GMT_stdout);	GMT_ascii_output_one (GMT_stdout, x_min, 0);
				GMT_fputs ("\t", GMT_stdout);	GMT_ascii_output_one (GMT_stdout, y_min, 1);
				GMT_fputs ("\t", GMT_stdout);	GMT_ascii_output_one (GMT_stdout, x_max, 0);
				GMT_fputs ("\t", GMT_stdout);	GMT_ascii_output_one (GMT_stdout, y_max, 1);
			}
			if (Ctrl->L.norm & 1) {
				GMT_fputs ("\t", GMT_stdout);	GMT_ascii_output_one (GMT_stdout, median, 2);
				GMT_fputs ("\t", GMT_stdout);	GMT_ascii_output_one (GMT_stdout, scale, 2);
			}
			if (Ctrl->L.norm & 2) {
				GMT_fputs ("\t", GMT_stdout);	GMT_ascii_output_one (GMT_stdout, mean, 2);
				GMT_fputs ("\t", GMT_stdout);	GMT_ascii_output_one (GMT_stdout, stdev, 2);
				GMT_fputs ("\t", GMT_stdout);	GMT_ascii_output_one (GMT_stdout, rms, 2);
			}
			if (Ctrl->M.active) {sprintf (buffer, "\t%d", (int)n_nan);	GMT_fputs(buffer, GMT_stdout);}
			GMT_fputs ("\n", GMT_stdout);
		}
		else if (!(Ctrl->T.active || (Ctrl->I.active && Ctrl->I.status == 2))) {
			sprintf (buffer, "%s: Title: %s\n", grd.name, grd.title);	GMT_fputs(buffer, GMT_stdout);
			sprintf (buffer, "%s: Command: %s\n", grd.name, grd.command);	GMT_fputs(buffer, GMT_stdout);
			sprintf (buffer, "%s: Remark: %s\n", grd.name, grd.remark);	GMT_fputs(buffer, GMT_stdout);
			if (grd.node_offset == 0 || grd.node_offset == 1)
				{sprintf (buffer, "%s: %s node registration used\n", grd.name, type[grd.node_offset]);	GMT_fputs(buffer, GMT_stdout);}
			else
				{sprintf (buffer, "%s: Unknown registration! Probably not a GMT grid\n", grd.name);	GMT_fputs(buffer, GMT_stdout);}
			if (grd.type >= 0 && grd.type < N_GRD_TYPES)
				{sprintf (buffer, "%s: Grid file format: %c%c (# %ld) %s\n", grd.name, (int)GMT_grdformats[grd.type][0], (int)GMT_grdformats[grd.type][1], grd.type, GMT_grd_type[grd.type-1]);	GMT_fputs(buffer, GMT_stdout);}
			else
				{sprintf (buffer, "%s: Unrecognized grid file format! Probably not a GMT grid\n", grd.name);	GMT_fputs(buffer, GMT_stdout);}
			if (Ctrl->F.active) {
				if ((fabs (grd.x_min) < 500.0) && (fabs (grd.x_max) < 500.0) && (fabs (grd.y_min) < 500.0) && (fabs (grd.y_max) < 500.0)) {
					sprintf (buffer, "%s: x_min: %.7f\n", grd.name, grd.x_min);	GMT_fputs(buffer, GMT_stdout);
					sprintf (buffer, "%s: x_max: %.7f\n", grd.name, grd.x_max);	GMT_fputs(buffer, GMT_stdout);
					sprintf (buffer, "%s: x_inc: %.7f\n", grd.name, grd.x_inc);	GMT_fputs(buffer, GMT_stdout);
					sprintf (buffer, "%s: name: %s\n", grd.name, grd.x_units);	GMT_fputs(buffer, GMT_stdout);
					sprintf (buffer, "%s: nx: %d\n", grd.name, grd.nx);		GMT_fputs(buffer, GMT_stdout);
					sprintf (buffer, "%s: y_min: %.7f\n", grd.name, grd.y_min);	GMT_fputs(buffer, GMT_stdout);
					sprintf (buffer, "%s: y_max: %.7f\n", grd.name, grd.y_max);	GMT_fputs(buffer, GMT_stdout);
					sprintf (buffer, "%s: y_inc: %.7f\n", grd.name, grd.y_inc);	GMT_fputs(buffer, GMT_stdout);
					sprintf (buffer, "%s: name: %s\n", grd.name, grd.y_units);	GMT_fputs(buffer, GMT_stdout);
					sprintf (buffer, "%s: ny: %d\n", grd.name, grd.ny);		GMT_fputs(buffer, GMT_stdout);
				}
				else {
					sprintf (buffer, "%s: x_min: %.2f\n", grd.name, grd.x_min);	GMT_fputs(buffer, GMT_stdout);
					sprintf (buffer, "%s: x_max: %.2f\n", grd.name, grd.x_max);	GMT_fputs(buffer, GMT_stdout);
					sprintf (buffer, "%s: x_inc: %.2f\n", grd.name, grd.x_inc);	GMT_fputs(buffer, GMT_stdout);
					sprintf (buffer, "%s: name: %s\n", grd.name, grd.x_units);	GMT_fputs(buffer, GMT_stdout);
					sprintf (buffer, "%s: nx: %d\n", grd.name, grd.nx);		GMT_fputs(buffer, GMT_stdout);
					sprintf (buffer, "%s: y_min: %.2f\n", grd.name, grd.y_min);	GMT_fputs(buffer, GMT_stdout);
					sprintf (buffer, "%s: y_max: %.2f\n", grd.name, grd.y_max);	GMT_fputs(buffer, GMT_stdout);
					sprintf (buffer, "%s: y_inc: %.2f\n", grd.name, grd.y_inc);	GMT_fputs(buffer, GMT_stdout);
					sprintf (buffer, "%s: name: %s\n", grd.name, grd.y_units);	GMT_fputs(buffer, GMT_stdout);
					sprintf (buffer, "%s: ny: %d\n", grd.name, grd.ny);		GMT_fputs(buffer, GMT_stdout);
				}
			}
			else {
				GMT_ascii_format_one (text, grd.x_inc, GMT_io.out_col_type[0]);
				if (isalpha ((int)text[strlen(text)-1])) text[strlen(text)-1] = '\0';	/* Chop of trailing WESN flag here */
				sprintf (buffer, "%s: x_min: ", grd.name);			GMT_fputs(buffer, GMT_stdout);
				GMT_ascii_output_one (GMT_stdout, grd.x_min, 0);
				GMT_fputs (" x_max: ", GMT_stdout);
				GMT_ascii_output_one (GMT_stdout, grd.x_max, 0);
				sprintf (buffer, " x_inc: %s", text);				GMT_fputs(buffer, GMT_stdout);
				sprintf (buffer, " name: %s nx: %d\n", grd.x_units, grd.nx);	GMT_fputs(buffer, GMT_stdout);
				sprintf (buffer, "%s: y_min: ", grd.name);			GMT_fputs(buffer, GMT_stdout);
				GMT_ascii_output_one (GMT_stdout, grd.y_min, 1);
				GMT_fputs (" y_max: ", GMT_stdout);
				GMT_ascii_output_one (GMT_stdout, grd.y_max, 1);
				GMT_ascii_format_one (text, grd.y_inc, GMT_io.out_col_type[1]);
				if (isalpha ((int)text[strlen(text)-1])) text[strlen(text)-1] = '\0';	/* Chop of trailing WESN flag here */
				sprintf (buffer, " y_inc: %s", text);				GMT_fputs(buffer, GMT_stdout);
				sprintf (buffer, " name: %s ny: %d\n", grd.y_units, grd.ny);	GMT_fputs(buffer, GMT_stdout);
			}

			if (Ctrl->M.active) {
				if (z_min == -DBL_MAX) z_min = GMT_d_NaN;
				if (z_max == +DBL_MAX) z_max = GMT_d_NaN;
				sprintf (buffer, "%s: z_min: ", grd.name);	GMT_fputs(buffer, GMT_stdout);
				GMT_ascii_output_one (GMT_stdout, z_min, 2);
				GMT_fputs (" at x = ", GMT_stdout);
				GMT_ascii_output_one (GMT_stdout, x_min, 0);
				GMT_fputs (" y = ", GMT_stdout);
				GMT_ascii_output_one (GMT_stdout, y_min, 1);
				GMT_fputs (" z_max: ", GMT_stdout);
				GMT_ascii_output_one (GMT_stdout, z_max, 2);
				GMT_fputs (" at x = ", GMT_stdout);
				GMT_ascii_output_one (GMT_stdout, x_max, 0);
				GMT_fputs (" y = ", GMT_stdout);
				GMT_ascii_output_one (GMT_stdout, y_max, 1);
				GMT_fputs ("\n", GMT_stdout);
			}
			else if (Ctrl->F.active) {
				sprintf (buffer, "%s: zmin: %g\n", grd.name, grd.z_min);	GMT_fputs(buffer, GMT_stdout);
				sprintf (buffer, "%s: zmax: %g\n", grd.name, grd.z_max);	GMT_fputs(buffer, GMT_stdout);
			 	sprintf (buffer, "%s: name: %s\n", grd.name, grd.z_units);	GMT_fputs(buffer, GMT_stdout);
			}
			else {
				sprintf (buffer, "%s: z_min: ", grd.name);	GMT_fputs(buffer, GMT_stdout);
				GMT_ascii_output_one (GMT_stdout, grd.z_min, 2);
				GMT_fputs (" z_max: ", GMT_stdout);
				GMT_ascii_output_one (GMT_stdout, grd.z_max, 2);
				sprintf (buffer, " name: %s\n", grd.z_units);	GMT_fputs(buffer, GMT_stdout);
			}

			GMT_ascii_format_one (text, grd.z_add_offset, GMT_io.out_col_type[2]);
			if (isalpha ((int)text[strlen(text)-1])) text[strlen(text)-1] = '\0';	/* Chop of trailing WESN flag here */
			sprintf (format, "%s: scale_factor: %s add_offset: %%s\n", grd.name, gmtdefs.d_format);
			sprintf (buffer, format, grd.z_scale_factor, text);	GMT_fputs(buffer, GMT_stdout);
			if (n_nan) fprintf (stdout, "%s: %ld nodes set to NaN\n", grd.name, n_nan);
			if (Ctrl->L.norm & 1) {
				sprintf (buffer, "%s: median: ", grd.name);	GMT_fputs(buffer, GMT_stdout);	GMT_ascii_output_one (GMT_stdout, median, 2);
				GMT_fputs (" scale: ", GMT_stdout);		GMT_ascii_output_one (GMT_stdout, scale, 2);
				GMT_fputs ("\n", GMT_stdout);
			}
			if (Ctrl->L.norm & 2) {
				sprintf (buffer, "%s: mean: ", grd.name);	GMT_fputs(buffer, GMT_stdout);	GMT_ascii_output_one (GMT_stdout, mean, 2);
				GMT_fputs (" stdev: ", GMT_stdout);		GMT_ascii_output_one (GMT_stdout, stdev, 2);
				GMT_fputs (" rms: ", GMT_stdout);		GMT_ascii_output_one (GMT_stdout, rms, 2);
				GMT_fputs ("\n", GMT_stdout);
			}
		}
		else {
			if (grd.z_min < global_zmin) global_zmin = grd.z_min;
			if (grd.z_max > global_zmax) global_zmax = grd.z_max;
			if (grd.x_min < global_xmin) global_xmin = grd.x_min;
			if (grd.x_max > global_xmax) global_xmax = grd.x_max;
			if (grd.y_min < global_ymin) global_ymin = grd.y_min;
			if (grd.y_max > global_ymax) global_ymax = grd.y_max;
		}
	}

	if (global_zmin == -DBL_MAX) global_zmin = GMT_d_NaN;	/* Never got set */
	if (global_zmax == +DBL_MAX) global_zmax = GMT_d_NaN;

	if (Ctrl->C.active && (Ctrl->I.active && Ctrl->I.status == 2)) {
		global_xmin = floor (global_xmin / Ctrl->I.xinc) * Ctrl->I.xinc;
		global_xmax = ceil  (global_xmax / Ctrl->I.xinc) * Ctrl->I.xinc;
		global_ymin = floor (global_ymin / Ctrl->I.yinc) * Ctrl->I.yinc;
		global_ymax = ceil  (global_ymax / Ctrl->I.yinc) * Ctrl->I.yinc;
		sprintf (buffer, "%d\t", (int)n_grds);	GMT_fputs(buffer, GMT_stdout);
		GMT_ascii_output_one (GMT_stdout, global_xmin, 0);	GMT_fputs ("\t", GMT_stdout);
		GMT_ascii_output_one (GMT_stdout, global_xmax, 0);	GMT_fputs ("\t", GMT_stdout);
		GMT_ascii_output_one (GMT_stdout, global_ymin, 1);	GMT_fputs ("\t", GMT_stdout);
		GMT_ascii_output_one (GMT_stdout, global_ymax, 1);	GMT_fputs ("\t", GMT_stdout);
		GMT_ascii_output_one (GMT_stdout, global_zmin, 1);	GMT_fputs ("\t", GMT_stdout);
		GMT_ascii_output_one (GMT_stdout, global_zmax, 1);	GMT_fputs ("\n", GMT_stdout);
	}
	else if (Ctrl->T.active) {
		global_zmin = floor (global_zmin / Ctrl->T.inc) * Ctrl->T.inc;
		global_zmax = ceil  (global_zmax / Ctrl->T.inc) * Ctrl->T.inc;
		GMT_fputs ("-T", GMT_stdout);
		GMT_ascii_output_one (GMT_stdout, global_zmin, 2);
		GMT_fputs ("/", GMT_stdout);
		GMT_ascii_output_one (GMT_stdout, global_zmax, 2);
		GMT_fputs ("/", GMT_stdout);
		GMT_ascii_output_one (GMT_stdout, Ctrl->T.inc, 2);
		GMT_fputs ("\n", GMT_stdout);
	}
	else if ((Ctrl->I.active && Ctrl->I.status == 2)) {
		global_xmin = floor (global_xmin / Ctrl->I.xinc) * Ctrl->I.xinc;
		global_xmax = ceil  (global_xmax / Ctrl->I.xinc) * Ctrl->I.xinc;
		global_ymin = floor (global_ymin / Ctrl->I.yinc) * Ctrl->I.yinc;
		global_ymax = ceil  (global_ymax / Ctrl->I.yinc) * Ctrl->I.yinc;
		GMT_fputs ("-R", GMT_stdout);
		GMT_ascii_output_one (GMT_stdout, global_xmin, 0);	GMT_fputs ("/", GMT_stdout);
		GMT_ascii_output_one (GMT_stdout, global_xmax, 0);	GMT_fputs ("/", GMT_stdout);
		GMT_ascii_output_one (GMT_stdout, global_ymin, 1);	GMT_fputs ("/", GMT_stdout);
		GMT_ascii_output_one (GMT_stdout, global_ymax, 1);	GMT_fputs ("\n", GMT_stdout);
	}

	GMT_free ((void *)a);

	Free_grdinfo_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void *New_grdinfo_Ctrl () {	/* Allocate and initialize a new control structure */
	struct GRDINFO_CTRL *C;

	C = (struct GRDINFO_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GRDINFO_CTRL), "New_grdinfo_Ctrl");

	/* Initialize values whose defaults are not 0/FALSE/NULL */

	return ((void *)C);
}

void Free_grdinfo_Ctrl (struct GRDINFO_CTRL *C) {	/* Deallocate control structure */
	GMT_free ((void *)C);
}
