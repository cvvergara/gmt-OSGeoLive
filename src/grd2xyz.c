/*--------------------------------------------------------------------
 *	$Id: grd2xyz.c 10173 2014-01-01 09:52:34Z pwessel $
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
 * grd2xyz.c reads a grid file and prints out the x,y,z values to
 * standard output.
 *
 * Author:	Paul Wessel
 * Date:	3-JAN-1991
 * Version:	4
 */

#include "gmt.h"

#include "gmt_parse_z_io.h"	/* To define the Z structure used for parsing */

struct GRD2XYZ_CTRL {
	struct E {	/* -E[f][<nodata>] */
		GMT_LONG active;
		GMT_LONG floating;
		double nodata;
	} E;
	struct N {	/* -N<nodata> */
		GMT_LONG active;
		double value;
	} N;
	struct S {	/* -S[r] */
		GMT_LONG active;
		GMT_LONG reverse;
	} S;
	struct W {	/* -W[<weight>] */
		GMT_LONG active;
		double weight;
	} W;
	struct Z Z;
};

int main (int argc, char **argv)
{
	GMT_LONG error = FALSE, b_only = FALSE, first = TRUE;

	GMT_LONG i, j, k, nx, ny, n_files = 0, n_out;
	GMT_LONG ij, nm, gmt_ij, n_suppressed = 0, n_total = 0;

	char  buffer[GMT_LONG_TEXT]; 
	float *z = NULL;
	
	double w, e, s, n, *x = NULL, *y = NULL, out[4], d_value;

	struct GRD_HEADER grd;
	struct GMT_Z_IO io;
	struct GRD2XYZ_CTRL *Ctrl = NULL;

	void *New_grd2xyz_Ctrl (), Free_grd2xyz_Ctrl (struct GRD2XYZ_CTRL *C);
	
	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct GRD2XYZ_CTRL *)New_grd2xyz_Ctrl ();	/* Allocate and initialize a new control structure */
	
	w = e = s = n = 0.0;

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				/* Common parameters */

				case 'b':
					b_only = TRUE;
				case 'H':
				case 'R':
				case 'V':
				case ':':
				case 'f':
				case '\0':
					error += GMT_parse_common_options (argv[i], &w, &e, &s, &n);
					break;

				/* Supplemental options */

				case 'E':
					Ctrl->E.active = TRUE;
					if (argv[i][2] == 'f') Ctrl->E.floating = TRUE;
					if (argv[i][2+Ctrl->E.floating]) Ctrl->E.nodata = atof (&argv[i][2+Ctrl->E.floating]);
					break;
				case 'L':	/* For backwards compatibility only; use -f instead */
					GMT_io.out_col_type[0] = GMT_IS_LON;
					GMT_io.out_col_type[1] = GMT_IS_LAT;
					break;
				case 'N':
					if (!argv[i][2]) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -N option:  Must specify value or NaN\n", GMT_program);
						error++;
					}
					else {
						Ctrl->N.value = (argv[i][2] == 'N' || argv[i][2] == 'n') ? GMT_d_NaN : atof (&argv[i][2]);
						Ctrl->N.active = TRUE;
					}
					break;
				case 'Z':
					Ctrl->Z.active = TRUE;
					error += GMT_parse_z_io (&argv[i][2], &Ctrl->Z);
					break;
				case 'S':
					Ctrl->S.active = TRUE;
					if (argv[i][2] == 'r') Ctrl->S.reverse = TRUE;
					break;
				case 'W':
					Ctrl->W.active = TRUE;
					Ctrl->W.weight = (argv[i][2]) ? atof (&argv[i][2]) : 1.0;
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
		fprintf (stderr, "grd2xyz %s - Converting netCDF grdfile(s) to ASCII xyz data\n\n", GMT_VERSION);
		fprintf( stderr, "usage: grd2xyz <grdfiles> [-E[f][<nodata>]] [%s] [-N<nodata>] [%s] [-S[r]] [-V]\n", GMT_Ho_OPT, GMT_Rgeo_OPT);
		fprintf( stderr, "\t[-W[<weight>]] [-Z[<flags>]] [%s] [%s] [%s] > xyzfile\n", GMT_t_OPT, GMT_bo_OPT, GMT_f_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\n\t<grdfiles> is one or more grid files to convert\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-E Write ESRI ArcInfo ASCII interchange format.  Only one grid file can be specified.\n");
		fprintf (stderr, "\t   Optionally append f for floating point output [Default is integer].\n");
		fprintf (stderr, "\t   Optionally append nodata value to represent NaNs [-9999].\n");
		fprintf (stderr, "\t-H Write 1 ASCII header record [Default is no header].\n");
		fprintf (stderr, "\t-N Replace z-values that equal NaN with this value [Default writes NaN].\n");
		GMT_explain_option ('R');
		fprintf (stderr, "\t-S Suppress output for nodes whose z equals NaN [Default prints all nodes].\n");
		fprintf (stderr, "\t   Append r to reverse the suppression (only output NaN nodes).\n");
		GMT_explain_option ('V');
		fprintf (stderr, "\t-W Write xyzw using supplied weight (or 1 if not given) [Default is xyz].\n");
		fprintf (stderr, "\t-Z sets exact specification of resulting 1-column output z-table.\n");
		fprintf (stderr, "\t   If data is in row format, state if first row is at T(op) or B(ottom);\n");
		fprintf (stderr, "\t     Then, append L or R to indicate starting point in row.\n");
		fprintf (stderr, "\t   If data is in column format, state if first columns is L(left) or R(ight);\n");
		fprintf (stderr, "\t     Then, append T or B to indicate starting point in column.\n");
		fprintf (stderr, "\t   Append x if gridline-registered, periodic data in x without repeating column at xmax.\n");
		fprintf (stderr, "\t   Append y if gridline-registered, periodic data in y without repeating row at ymax.\n");
		fprintf (stderr, "\t   Specify one of the following data types (all binary except a):\n");
		fprintf (stderr, "\t     a  Ascii\n");
		fprintf (stderr, "\t     c  signed 1-byte character\n");
		fprintf (stderr, "\t     u  unsigned 1-byte character\n");
		fprintf (stderr, "\t     h  signed short 2-byte integer\n");
		fprintf (stderr, "\t     H  unsigned short 2-byte integer\n");
		fprintf (stderr, "\t     i  signed 4-byte integer\n");
		fprintf (stderr, "\t     I  unsigned 4-byte integer\n");
		fprintf (stderr, "\t     l  long (4- or 8-byte) integer\n");
		fprintf (stderr, "\t     f  4-byte floating point single precision\n");
		fprintf (stderr, "\t     d  8-byte floating point double precision\n");
		fprintf (stderr, "\t   [Default format is scanline orientation in ascii representation: -ZTLa].\n");
		GMT_explain_option (':');
		GMT_explain_option ('o');
		GMT_explain_option ('n');
		GMT_explain_option ('f');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	if (n_files == 0) {
		fprintf (stderr, "%s: SYNTAX ERROR:  Must specify at least one input file\n", GMT_program);
		error++;
	}

	if (n_files > 1 && Ctrl->E.active) {
		fprintf (stderr, "%s: SYNTAX ERROR:  -E can only handle one input file\n", GMT_program);
		error++;
	}

	if (Ctrl->Z.active && Ctrl->E.active) {
		fprintf (stderr, "%s: SYNTAX ERROR:  -E is not compatible with -Z\n", GMT_program);
		error++;
	}
	if (b_only && Ctrl->Z.active) GMT_io.binary[GMT_OUT] = FALSE;

	GMT_init_z_io (Ctrl->Z.format, Ctrl->Z.repeat, Ctrl->Z.swab, Ctrl->Z.skip, Ctrl->Z.type, &io);

	if ((GMT_io.binary[GMT_OUT] || io.binary) && GMT_io.io_header[GMT_OUT]) {
		fprintf (stderr, "%s: SYNTAX ERROR.  Binary output data cannot have header -H\n", GMT_program);
		error++;
	}

	if (error) exit (EXIT_FAILURE);

	if (Ctrl->Z.active && io.binary) GMT_io.binary[GMT_OUT] = TRUE;
	
	if (b_only && Ctrl->Z.active) fprintf (stderr, "%s: Warning.  -Z overrides -bo\n", GMT_program);
	if (b_only && Ctrl->E.active) fprintf (stderr, "%s: Warning.  -E overrides -bo\n", GMT_program);

#ifdef SET_IO_MODE
		GMT_setmode (GMT_OUT);
#endif

	n_out = (Ctrl->W.active) ? 4 : 3;
	out[3] = Ctrl->W.weight;
	for (k = 1; k < argc; k++) {
		if (argv[k][0] == '-') continue;	/* Skip the options */

		GMT_err_fail (GMT_read_grd_info (argv[k], &grd), argv[k]);

		if (gmtdefs.verbose) fprintf (stderr, "%s: Working on file %s\n", GMT_program, argv[k]);

		if (e > w && n > s) {	/* Subset */
			GMT_err_fail (GMT_adjust_loose_wesn (&w, &e, &s, &n, &grd), "");	/* Make sure w,e,s,n matches header spacing */
			nx = GMT_get_n (w, e, grd.x_inc, grd.node_offset);
			ny = GMT_get_n (s, n, grd.y_inc, grd.node_offset);
			nm = GMT_get_nm (nx, ny);
		}
		else
			nm = GMT_get_nm (grd.nx, grd.ny);

		z = (float *) GMT_memory (VNULL, (size_t) nm, sizeof (float), GMT_program);
		GMT_err_fail (GMT_read_grd (argv[k], &grd, z, w, e, s, n, GMT_pad, FALSE), argv[k]);

		n_total += nm;

		GMT_err_fail (GMT_set_z_io (&io, &grd), argv[k]);

		if (Ctrl->Z.active) {
			if (GMT_io.io_header[GMT_OUT] && !io.binary) {sprintf (buffer, "%s\n", grd.z_units);	GMT_fputs(buffer, GMT_stdout);}

			for (ij = 0; ij < io.n_expected; ij++) {
				(io.get_gmt_ij) (&io, ij, &gmt_ij);
				d_value = z[gmt_ij];
				if (Ctrl->S.active && (GMT_is_dnan (d_value) + Ctrl->S.reverse) == 1) {
					n_suppressed++;
					continue;
				}
				if ((io.x_missing && io.gmt_i == io.x_period) || (io.y_missing && io.gmt_j == 0)) continue;
				if (Ctrl->N.active && GMT_is_dnan (d_value)) d_value = Ctrl->N.value;
				(io.write_item) (GMT_stdout, d_value);
			}
		}
		else if (Ctrl->E.active) {
			double slop;
			slop = 1.0 - (grd.x_inc / grd.y_inc);
			if (!GMT_IS_ZERO (slop)) {
				fprintf (stderr, "%s: ERROR: x_inc must equal y_inc when writing to ESRI format\n", GMT_program);
				exit (EXIT_FAILURE);
			}
			sprintf (buffer, "ncols %d\nnrows %d\n", grd.nx, grd.ny); GMT_fputs(buffer, GMT_stdout);
			if (grd.node_offset) {	/* Pixel format */
				GMT_fputs ("xllcorner ", GMT_stdout);
				sprintf (buffer, gmtdefs.d_format, grd.x_min);	GMT_fputs(buffer, GMT_stdout);
				GMT_fputs ("\nyllcorner ", GMT_stdout);
				sprintf (buffer, gmtdefs.d_format, grd.y_min);	GMT_fputs(buffer, GMT_stdout);
			}
			else {	/* Gridline format */
				GMT_fputs ("xllcenter ", GMT_stdout);
				sprintf (buffer, gmtdefs.d_format, grd.x_min);	GMT_fputs(buffer, GMT_stdout);
				GMT_fputs ("\nyllcenter ", GMT_stdout);
				sprintf (buffer, gmtdefs.d_format, grd.y_min);	GMT_fputs(buffer, GMT_stdout);
			}
			GMT_fputs ("\ncellsize ", GMT_stdout);
			sprintf (buffer, gmtdefs.d_format, grd.x_inc);		GMT_fputs(buffer, GMT_stdout);
			sprintf (buffer, "\nnodata_value %ld\n", (GMT_LONG)irint (Ctrl->E.nodata));	GMT_fputs(buffer, GMT_stdout);
			for (j = 0; j < grd.ny; j++) {	/* Scanlines, starting in the north (ymax) */
				ij = GMT_IJ (j, 0, grd.nx);
				for (i = 0; i < grd.nx; i++, ij++) {
					if (GMT_is_fnan (z[ij]))
						{sprintf (buffer, "%ld", (GMT_LONG)irint (Ctrl->E.nodata));	GMT_fputs(buffer, GMT_stdout);}
					else if (Ctrl->E.floating)
						{sprintf (buffer, gmtdefs.d_format, z[ij]);			GMT_fputs(buffer, GMT_stdout);}
					else
						{sprintf (buffer, "%ld", (GMT_LONG)irint ((double)z[ij]));	GMT_fputs(buffer, GMT_stdout);}
					if (i < (grd.nx-1)) GMT_fputs (" ", GMT_stdout);
				}
				GMT_fputs ("\n", GMT_stdout);
			}
		}
		else {

			x = (double *) GMT_memory (VNULL, (size_t) grd.nx, sizeof (double), GMT_program);
			y = (double *) GMT_memory (VNULL, (size_t) grd.ny, sizeof (double), GMT_program);

			/* Compute grid node positions once only */

			for (j = 0; j < grd.ny; j++) y[j] = GMT_j_to_y (j, grd.y_min, grd.y_max, grd.y_inc, grd.xy_off, grd.ny);
			for (i = 0; i < grd.nx; i++) x[i] = GMT_i_to_x (i, grd.x_min, grd.x_max, grd.x_inc, grd.xy_off, grd.nx);

			if (GMT_io.io_header[GMT_OUT] && first) {
				if (!grd.x_units[0]) strcpy (grd.x_units, "x");
				if (!grd.y_units[0]) strcpy (grd.y_units, "y");
				if (!grd.z_units[0]) strcpy (grd.z_units, "z");
				if (gmtdefs.xy_toggle[GMT_IN])
					{sprintf (buffer, "%s%s%s%s%s", grd.y_units, gmtdefs.field_delimiter, grd.x_units, gmtdefs.field_delimiter, grd.z_units);	GMT_fputs(buffer, GMT_stdout);}
				else
					{sprintf (buffer, "%s%s%s%s%s", grd.x_units, gmtdefs.field_delimiter, grd.y_units, gmtdefs.field_delimiter, grd.z_units);	GMT_fputs(buffer, GMT_stdout);}
				if (Ctrl->W.active) {
					GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
					GMT_fputs ("weight\n", GMT_stdout);
				}
				else
					GMT_fputs ("\n", GMT_stdout);
				first = FALSE;
			}

			for (j = ij = 0; j < grd.ny; j++) for (i = 0; i < grd.nx; i++, ij++) {
				out[GMT_Z] = z[ij];
				if (Ctrl->S.active && (GMT_is_dnan (out[GMT_Z]) + Ctrl->S.reverse) == 1) {
					n_suppressed++;
					continue;
				}
				out[GMT_X] = x[i];	out[GMT_Y] = y[j];
				if (Ctrl->N.active && GMT_is_dnan (out[GMT_Z])) out[GMT_Z] = Ctrl->N.value;
				GMT_output (GMT_stdout, n_out, out);
			}
			GMT_free ((void *)x);
			GMT_free ((void *)y);
		}

		GMT_free ((void *)z);
	}

	if (gmtdefs.verbose) fprintf (stderr, "%s: %ld values extracted\n", GMT_program, n_total - n_suppressed);
	if (n_suppressed && gmtdefs.verbose) {
		if (Ctrl->S.reverse)
			fprintf (stderr, "%s: %ld finite values suppressed\n", GMT_program, n_suppressed);
		else
			fprintf (stderr, "%s: %ld NaN values suppressed\n", GMT_program, n_suppressed);
	}

	Free_grd2xyz_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void *New_grd2xyz_Ctrl () {	/* Allocate and initialize a new control structure */
	struct GRD2XYZ_CTRL *C;
	
	C = (struct GRD2XYZ_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GRD2XYZ_CTRL), "New_grd2xyz_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	
	C->E.nodata = -9999.0;
	C->W.weight = 1.0;
	C->Z.type = 'a';
	C->Z.format[0] = 'T';	C->Z.format[1] = 'L';
		
	return ((void *)C);
}

void Free_grd2xyz_Ctrl (struct GRD2XYZ_CTRL *C) {	/* Deallocate control structure */
	GMT_free ((void *)C);	
}
