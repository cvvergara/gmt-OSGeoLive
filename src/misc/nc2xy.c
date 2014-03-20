/*--------------------------------------------------------------------
 *	$Id: nc2xy.c 10173 2014-01-01 09:52:34Z pwessel $
 *
 *	Copyright (c) 2006-2014 by R. Scharroo
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
 * nc2xy.c reads a netCDF file and prints out the x,y values to
 * standard output.
 *
 * Author:	Remko Scharroo
 * Date:	28-Nov-2005
 * Version:	1
 */

#define NVAR 10
#define GMT_WITH_NO_PS
#include "gmt.h"

int main (int argc, char **argv)
{
	GMT_LONG error = FALSE, suppress = FALSE, reverse = FALSE, skip, plus_const_col = FALSE;
	int ncid, varid[NVAR], dimid = 0, ndims;
	GMT_LONG i, k, n_files = 0, n_out;
	int n_dimid = -1;
	size_t n_total = 0, j, n_suppressed = 0;
	double w, e, s, n, out[NVAR+1], add_offset[NVAR], scale_factor[NVAR], missing_value[NVAR];
	double cte_col = 0.0;
	char varnm[NVAR][GMT_TEXT_LEN], long_name[GMT_LONG_TEXT], units[GMT_LONG_TEXT];
	struct GMT_TIME_SYSTEM time_system;

	argc = (int)GMT_begin (argc, argv);

	w = e = s = n = 0.0;

	for (i = 0; i < NVAR; i++) {
		varnm[i][0] = '\0';
	}

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				/* Common parameters */

				case 'f':
				case 'b':
				case 'V':
				case '\0':
					error += GMT_parse_common_options (argv[i], &w, &e, &s, &n);
					break;

				/* Supplemental options */

				case 'C':
					cte_col = atof(&argv[i][2]);
					plus_const_col = TRUE;
					break;
				case 'F':
					sscanf (&argv[i][2], "%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]/%[^/]", varnm[0], varnm[1], varnm[2], varnm[3], varnm[4], varnm[5], varnm[6], varnm[7], varnm[8], varnm[9]);
					break;
				case 'S':
					suppress = TRUE;
					if (argv[i][2] == 'r') reverse = TRUE;
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
		fprintf (stderr, "nc2xy %s - Converting netCDF column file(s) to ASCII xy data\n\n", GMT_VERSION);
		fprintf( stderr, "usage: nc2xy <files> [-F<var1>/<var2>/...] [-S[r]] [-V] [%s] [-bo] > xyfile\n", GMT_f_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\n\t<files> is one or more netCDF files to convert\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-F Specify variables to be extracted (up to %d)\n", NVAR);
		fprintf (stderr, "\t-S Suppress records with NaN values [Default prints all nodes]\n");
		fprintf (stderr, "\t   Append r to reverse the suppression (only output records with NaNs)\n");
		GMT_explain_option ('V');
		GMT_explain_option ('f');
		GMT_explain_option ('o');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	if (n_files == 0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify at least one input file\n", GMT_program);
		error++;
	}

	if (error) exit (EXIT_FAILURE);

#ifdef SET_IO_MODE
	GMT_setmode (GMT_OUT);
#endif

	for (k = 1; k < argc; k++) {
		if (argv[k][0] == '-') continue;	/* Skip the options */

		/* Open netCDF file */

		if (!strcmp (argv[k],"=")) {
			fprintf (stderr, "%s: %s [%s]\n", GMT_program, GMT_strerror(GMT_GRDIO_NC_NO_PIPE), argv[k]);
			exit (EXIT_FAILURE);
		}
		GMT_err_fail (nc_open (argv[k], NC_NOWRITE, &ncid), argv[k]);

		if (gmtdefs.verbose) fprintf (stderr, "%s: Working on file %s\n", GMT_program, argv[k]);

		/* Determine IDs of requested variables; take first two if none given */

		for (n_out = 0; n_out < NVAR && varnm[n_out][0]; n_out++)
			GMT_err_fail (nc_inq_varid (ncid, varnm[n_out], &varid[n_out]), argv[k]);

		if (n_out == 0) {
			n_out = 2;
			varid[0] = 0;
			varid[1] = 1;
		}

		/* Get further variable sizes and attributes */

		for (i = 0; i < n_out; i++) {
			/* Check that variable is 1-dimensional */
			GMT_err_fail (nc_inq_varndims (ncid, varid[i], &ndims), argv[k]);
			if (ndims != 1) {
				fprintf (stderr, "%s: Variable %s is not 1-dimensional\n", GMT_program, varnm[i]);
				exit (EXIT_FAILURE);
			}
			/* Check that dimension is same as the others */
			GMT_err_fail (nc_inq_vardimid (ncid, varid[i], &dimid), argv[k]);
			if (n_dimid < 0) {
				n_dimid = dimid;
				GMT_err_fail (nc_inq_dimlen (ncid, dimid, &n_total), argv[k]);
			}
			else if (n_dimid != dimid) {
				fprintf (stderr, "%s: Variable %s has different dimension than others\n", GMT_program, varnm[i]);
				exit (EXIT_FAILURE);
			}
			/* Check all attributes */
			if (nc_get_att_double (ncid, varid[i], "scale_factor", &scale_factor[i])) scale_factor[i] = 1.0;
			if (nc_get_att_double (ncid, varid[i], "add_offset", &add_offset[i])) add_offset[i] = 0.0;
			if (nc_get_att_double (ncid, varid[i], "_FillValue", &missing_value[i]) &&
			    nc_get_att_double (ncid, varid[i], "missing_value", &missing_value[i])) missing_value[i] = GMT_d_NaN;
			if (GMT_nc_get_att_text (ncid, varid[i], "long_name", long_name, (size_t)GMT_LONG_TEXT)) long_name[0] = 0;
			if (GMT_nc_get_att_text (ncid, varid[i], "units", units, (size_t)GMT_LONG_TEXT)) units[0] = 0;
			GMT_str_tolower (long_name); GMT_str_tolower (units);

			/* Scan for geographical or time units */
			if (!strcmp (long_name, "longitude") || strstr (units, "degrees_e"))
				GMT_io.in_col_type[i] = GMT_IS_LON;
			else if (!strcmp (long_name, "latitude") || strstr (units, "degrees_n"))
				GMT_io.in_col_type[i] = GMT_IS_LAT;
			else if (!strcmp (long_name, "time") || !strcmp (varnm[i], "time")) {
				GMT_io.in_col_type[i] = GMT_IS_RELTIME;
				memcpy ((void *)&time_system, (void *)&gmtdefs.time_system, sizeof (struct GMT_TIME_SYSTEM));
				if (GMT_get_time_system (units, &time_system) || GMT_init_time_system_structure (&time_system))
					fprintf (stderr, "%s: Warning: Time units [%s] in NetCDF file not recognised, defaulting to gmtdefaults.\n", GMT_program, units);
				/* Determine scale between data and internal time system, as well as the offset (in internal units) */
				scale_factor[i] = time_system.scale * gmtdefs.time_system.i_scale;
				add_offset[i] = (time_system.rata_die - gmtdefs.time_system.rata_die) + (time_system.epoch_t0 - gmtdefs.time_system.epoch_t0);
				add_offset[i] *= GMT_DAY2SEC_F * gmtdefs.time_system.i_scale;
			}
			else if (GMT_io.in_col_type[i] == GMT_IS_UNKNOWN)
				GMT_io.in_col_type[i] = GMT_IS_FLOAT;
		}

		/* Load data record by record and scale as required */

		for (j = 0; j < n_total; j++) {
			skip = FALSE;
			for (i = 0; i < n_out; i++) {
				GMT_err_fail (nc_get_var1_double (ncid, varid[i], &j, &out[i]), argv[k]);
				if (out[i] == missing_value[i])
					out[i] = GMT_d_NaN;
				else
					out[i] = out[i] * scale_factor[i] + add_offset[i];
				if (suppress && GMT_is_dnan (out[i])) {
					skip = TRUE;
					continue;
				}
			}
			if (plus_const_col) out[i] = cte_col;	/* Extra column with a const value on output */
			if (skip == reverse)
				GMT_output (GMT_stdout, n_out+plus_const_col, out);
			else
				n_suppressed++;
		}
	}

	if (gmtdefs.verbose) fprintf (stderr, "%s: %ld values extracted\n", GMT_program, (GMT_LONG)(n_total - n_suppressed));
	if (n_suppressed && gmtdefs.verbose) {
		if (reverse)
			fprintf (stderr, "%s: %ld finite values suppressed\n", GMT_program, (GMT_LONG)n_suppressed);
		else
			fprintf (stderr, "%s: %ld NaN values suppressed\n", GMT_program, (GMT_LONG)n_suppressed);
	}

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}
