/*--------------------------------------------------------------------
 *	$Id: grdreformat.c 10173 2014-01-01 09:52:34Z pwessel $
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
 * grdreformat.c reads a grid file in one format and outputs it in another
 *
 * Author:	Paul Wessel
 * Date:	3-JAN-1991-2011
 * Version:	4
 */

#include "gmt.h"

struct GRDREFORMAT_CTRL {
	struct N {	/* -N */
		GMT_LONG active;
	} N;
};

int main (int argc, char **argv)
{
#include "grdreformat.h"	/* Defines N_GRDTXT_LINES and char *grd_formats[N_GRDTXT_LINES] array used in usage message */
	GMT_LONG error = FALSE, global = FALSE, no_header;

	GMT_LONG i, nfiles = 0, nx, ny, type[2];

	size_t nm;
	
	double w, e, s, n;

	float *z = NULL;

	char *file[2], fname[2][BUFSIZ];

	struct GRD_HEADER grd;
	struct GRDREFORMAT_CTRL *Ctrl = NULL;

	void *New_grdreformat_Ctrl (), Free_grdreformat_Ctrl (struct GRDREFORMAT_CTRL *C);

	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct GRDREFORMAT_CTRL *)New_grdreformat_Ctrl ();	/* Allocate and initialize a new control structure */
	
	file[0] = file[1] = CNULL;
	w = e = s = n = 0.0;

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				/* Common parameters */

				case 'R':
				case 'V':
				case 'f':
				case '\0':
					error += GMT_parse_common_options (argv[i], &w, &e, &s, &n);
					break;

				case 'N':
					Ctrl->N.active = TRUE;
					break;
				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else if (nfiles < 2) {
			file[nfiles] = argv[i];
			nfiles++;
		}
		else
			nfiles++;
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {	/* Display usage */
		fprintf (stderr, "grdreformat %s - Converting between different grid file formats\n\n", GMT_VERSION);
		fprintf( stderr, "usage: grdreformat ingrdfile[=id[/scale/offset]] outgrdfile[=id[/scale/offset]] [-N]\n\t[%s] [%s] [-V]\n", GMT_Rgeo_OPT, GMT_f_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\tingrdfile is the grid file to convert.\n");
		fprintf (stderr, "\toutgrdfile is the new converted grid file.\n");
		fprintf( stderr, "\tscale and offset, if given, will multiply data by scale and add offset.\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-N  Do NOT write the header (native grids only - ignored otherwise).\n");
		fprintf (stderr, "\t\t  Useful when creating files to be used by grdraster.\n");
		GMT_explain_option ('r');
		GMT_explain_option ('f');
		GMT_explain_option ('V');

		fprintf (stderr, "\n	The following grid file formats are supported:\n\n");
#ifdef USE_GDAL
		for (i = 0; i < N_GRDTXT_LINES; i++) fprintf (stderr, "\t%s\n", grd_formats[i]);
#else
		for (i = 0; i < N_GRDTXT_LINES; i++) if (i != 33) fprintf (stderr, "\t%s\n", grd_formats[i]);
#endif
		exit (EXIT_FAILURE);
	}

	if (nfiles != 2) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR: Must specify both input and output file names.\n", GMT_program);
		error++;
	}

	if (error) exit (EXIT_FAILURE);

	GMT_grd_init (&grd, argc, argv, FALSE);
	no_header = (Ctrl->N.active) ? 64 : 0;
	GMT_err_fail (GMT_grd_get_format (file[0], &grd, TRUE), file[0]);
	type[0] = grd.type;
	strcpy (fname[0], grd.name);
	GMT_err_fail (GMT_grd_get_format (file[1], &grd, FALSE), file[1]);
	type[1] = grd.type;
	strcpy (fname[1], grd.name);

	if (type[1] == 20) {	/* Golden Surfer format 7 (double) is read-only */
		fprintf (stderr, "%s: Grid format sd (Golden Software Surfer format 7 (double)) is read-only!\n", GMT_program);
		exit (EXIT_FAILURE);
	}
#ifdef USE_GDAL
	if (type[1] == 22) {	/* GDAL format is read-only */
		fprintf (stderr, "%s: Grid format gd (GDAL) is read-only!\n", GMT_program);
		exit (EXIT_FAILURE);
	}
#endif	
	if (gmtdefs.verbose) {
		if (file[0][0] == '=') strcpy (fname[0], "<stdin>");
		if (file[1][0] == '=') strcpy (fname[1], "<stdout>");
		fprintf (stderr, "%s: Translating file %s (format = %ld) to file %s (format = %ld)\n", GMT_program, fname[0], type[0], fname[1], type[1]);
		if (no_header && GMT_grdformats[type[1]][0] != 'c' && GMT_grdformats[type[1]][0] != 'n') fprintf (stderr, "%s: No grd header will be written\n", GMT_program);
	}

	GMT_err_fail (GMT_read_grd_info (file[0], &grd), fname[0]);

	if (e > w && n > s) {
		global = GMT_grd_is_global (&grd);
		if (!global && (w < grd.x_min || e > grd.x_max)) error = TRUE;
		if (s < grd.y_min || n > grd.y_max) error = TRUE;
		if (error) {
			fprintf (stderr, "%s: Subset exceeds data domain!\n", GMT_program);
			exit (EXIT_FAILURE);
		}
		nx = GMT_get_n (w, e, grd.x_inc, grd.node_offset);
		ny = GMT_get_n (s, n, grd.y_inc, grd.node_offset);
		nm = GMT_get_nm (nx, ny);
	}
	else
		nm = GMT_get_nm (grd.nx, grd.ny);

	z = (float *) GMT_memory (VNULL, nm, sizeof (float), GMT_program);
	GMT_err_fail (GMT_read_grd (file[0], &grd, z, w, e, s, n, GMT_pad, FALSE), fname[0]);

	grd.type = type[1];

	GMT_grd_init (&grd, argc, argv, TRUE);

	GMT_err_fail (GMT_write_grd (file[1], &grd, z, 0.0, 0.0, 0.0, 0.0, GMT_pad, no_header), fname[1]);

	GMT_free ((void *)z);

	Free_grdreformat_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void *New_grdreformat_Ctrl () {	/* Allocate and initialize a new control structure */
	struct GRDREFORMAT_CTRL *C;
	
	C = (struct GRDREFORMAT_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GRDREFORMAT_CTRL), "New_grdreformat_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	
	return ((void *)C);
}

void Free_grdreformat_Ctrl (struct GRDREFORMAT_CTRL *C) {	/* Deallocate control structure */
	GMT_free ((void *)C);	
}
