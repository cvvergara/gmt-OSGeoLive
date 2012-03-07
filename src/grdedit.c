/*--------------------------------------------------------------------
 *	$Id: grdedit.c,v 1.66 2011/07/08 21:27:06 guru Exp $
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
 * grdedit reads an existing grid file and takes command line
 * arguments to redefine some of the grdheader parameters:
 *
 *	x_min/x_max OR x_inc (the other is recomputed)
 *	y_min/y_max OR y_inc (the other is recomputed)
 *	z_scale_factor/z_add_offset
 *	x_units/y_units/z_units
 *	title/command/remark
 *
 * Author:	Paul Wessel
 * Date:	9-SEP-1998
 * Version:	4
 * Modified:	19-FEB-2000 PW: Added -V option
 *		24-FEB-2006 PW: Added -T option
 */
 
#include "gmt.h"

struct GRDEDIT_CTRL {
	struct A {	/* -A */
		GMT_LONG active;
	} A;
	struct D {	/* -D<xname>/<yname>/<zname>/<scale>/<offset>/<title>/<remark> */
		GMT_LONG active;
		char *information;
	} D;
	struct E {	/* -E */
		GMT_LONG active;
	} E;
	struct N {	/* N<xyzfile> */
		GMT_LONG active;
		char *file;
	} N;
	struct S {	/* -S */
		GMT_LONG active;
	} S;
	struct T {	/* -T */
		GMT_LONG active;
	} T;
};

int main (int argc, char **argv)
{
	GMT_LONG error = FALSE;

	GMT_LONG i, j, n_expected_fields, n_fields;
	
	GMT_LONG n_data, k, nm;
	
	float *a = NULL;

	double w, e, s, n, shift_amount = 0.0, f, *in = NULL;

	char *grdfile = NULL, buffer[BUFSIZ];
	char *registration[2] = {"gridline", "pixel"};

	FILE *fp = NULL;

	struct GRD_HEADER grd;
	struct GRDEDIT_CTRL *Ctrl = NULL;

	void *New_grdedit_Ctrl (), Free_grdedit_Ctrl (struct GRDEDIT_CTRL *C);
	
	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct GRDEDIT_CTRL *)New_grdedit_Ctrl ();	/* Allocate and initialize a new control structure */
	
	grdfile = CNULL;
	w = e = s = n = 0.0;

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				/* Common parameters */

				case 'R':
				case 'H':
				case 'V':
				case ':':
  				case 'b':
  				case 'f':
				case '\0':
					error += GMT_parse_common_options (argv[i], &w, &e, &s, &n);
					break;

				/* Supplemental parameters */

				case 'A':
					Ctrl->A.active = TRUE;
					break;
				case 'D':
					Ctrl->D.active = TRUE;
					Ctrl->D.information = strdup (&argv[i][2]);
					break;
				case 'E':
					Ctrl->E.active = TRUE;
					break;
				case 'N':
					Ctrl->N.active = TRUE;
					Ctrl->N.file = strdup (&argv[i][2]);
					break;
				case 'S':
					Ctrl->S.active = TRUE;
					break;
				case 'T':
					Ctrl->T.active = TRUE;
					break;
				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else
			grdfile = argv[i];
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "grdedit %s - Modifying the header or content of a 2-D grid file\n\n", GMT_VERSION);
		fprintf (stderr, "usage: grdedit grdfile [-A] [%s]\n", GMT_GRDEDIT);
		fprintf (stderr, "\t[-E] [%s] [-N<xyzfile>] [%s] [-S] [-T]\n", GMT_H_OPT, GMT_Rgeo_OPT);
		fprintf (stderr, "\t[-V] [%s] [%s] [%s]\n", GMT_t_OPT, GMT_bi_OPT, GMT_f_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\tgrdfile is file to be modified\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-A will adjust dx/dy to be compatible with the files domain or -R\n");
		fprintf (stderr, "\t-D to enter information.  Specify '=' to get default value\n");
		fprintf (stderr, "\t-E to tranpose the entire grid (this will exchange x and y)\n");
		GMT_explain_option ('H');
		fprintf (stderr, "\t-N <file> has new xyz values to replace existing grid nodes\n");
		GMT_explain_option ('R');
		fprintf (stderr, "\t-S For global grids of 360 degree longitude range.\n");
		fprintf (stderr, "\t   Will rotate entire grid to coincide with new borders in -R\n");
		fprintf (stderr, "\t-T Toggle header from grid-line to pixel-registered grid or vice versa.\n");
		fprintf (stderr, "\t   This shrinks -R by 0.5*{dx,dy} going from pixel to grid-line registration\n");
		fprintf (stderr, "\t   and expands -R by 0.5*{dx,dy} going from grid-line to pixel registration.\n");
		GMT_explain_option ('V');
		GMT_explain_option (':');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		GMT_explain_option ('f');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	if (Ctrl->S.active && Ctrl->A.active) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -S option:  Incompatible with -A\n", GMT_program);
		error++;
	}
	if (Ctrl->E.active && (Ctrl->A.active || Ctrl->D.active || Ctrl->N.active || Ctrl->S.active || Ctrl->T.active)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -E option:  Incompatible with -A, -D, -N, -S, OR -T\n", GMT_program);
		error++;
	}
	if (Ctrl->S.active && Ctrl->T.active) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -S option:  Incompatible with -T\n", GMT_program);
		error++;
	}
	if (Ctrl->S.active && Ctrl->N.active) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -S option:  Incompatible with -N\n", GMT_program);
		error++;
	}
	if (Ctrl->S.active && !project_info.region_supplied) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -S option:  Must also specify -R\n", GMT_program);
		error++;
	}
	if (Ctrl->S.active && !GMT_360_RANGE (w,e)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -S option:  -R longitudes must span exactly 360 degrees\n", GMT_program);
		error++;
	}
	if (!grdfile) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify input file\n", GMT_program);
		error++;
	}
	if (Ctrl->N.active) {
		if (!Ctrl->N.file) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR -N option:  Must specify name of xyz file\n", GMT_program);
			error++;
		}
		if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
			error++;
		}
		if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = 3;
		if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] < 2) {
			fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data (-bi) must have at least 2 columns\n", GMT_program);
			error++;
		}
	}
	if (error) exit (EXIT_FAILURE);

	if (!strcmp (grdfile,  "=")) {
		fprintf (stderr, "%s: Piping of grid file not supported!\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	GMT_grd_init (&grd, argc, argv, TRUE);
	GMT_err_fail (GMT_read_grd_info (grdfile, &grd), grdfile);

	if ((grd.type == 6 || grd.type == 20) && Ctrl->T.active) {
		fprintf (stderr, "%s: Toggling registrations not possible for Surfer grid formats\n", GMT_program);
		fprintf (stderr, "%s: (Use grdreformat to convert to GMT default format and work on that file)\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	
	if (Ctrl->S.active && !GMT_grd_is_global (&grd)) {
		fprintf (stderr, "%s: Shift only allowed for global grids\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	if (gmtdefs.verbose) fprintf (stderr, "%s: Editing parameters for file %s\n", GMT_program, grdfile);

	/* Decode grd information given, if any */

	if (Ctrl->D.active) {
		double scale_factor, add_offset;
		if (gmtdefs.verbose) fprintf (stderr, "%s: Decode and change attributes in file %s\n", GMT_program, grdfile);
		scale_factor = grd.z_scale_factor;
		add_offset = grd.z_add_offset;
		GMT_decode_grd_h_info (Ctrl->D.information, &grd);
		if (scale_factor != grd.z_scale_factor || add_offset != grd.z_add_offset) {
			grd.z_min = (grd.z_min - add_offset) / scale_factor * grd.z_scale_factor + grd.z_add_offset;
			grd.z_max = (grd.z_max - add_offset) / scale_factor * grd.z_scale_factor + grd.z_add_offset;
		}
	}

	if (Ctrl->S.active) {
		shift_amount = w - grd.x_min;
		if (gmtdefs.verbose) fprintf (stderr, "%s: Shifting longitudes in file %s by %g degrees\n", GMT_program, grdfile, shift_amount);
		nm = GMT_get_nm (grd.nx, grd.ny);
		a = (float *) GMT_memory (VNULL, (size_t)nm, sizeof (float), GMT_program);
		GMT_err_fail (GMT_read_grd (grdfile, &grd, a, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), grdfile);
		GMT_grd_shift (&grd, a, shift_amount);
		GMT_err_fail (GMT_write_grd (grdfile, &grd, a, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), grdfile);
		GMT_free ((void *)a);
	}
	else if (Ctrl->N.active) {
		char *not_used = NULL;
		if (gmtdefs.verbose) fprintf (stderr, "%s: Replacing nodes using xyz values from file %s\n", GMT_program, Ctrl->N.file);
		if (GMT_io.binary[GMT_IN] && gmtdefs.verbose) {
			char *type[2] = {"double", "single"};
			fprintf (stderr, "%s: Expects %ld-column %s-precision binary data\n", GMT_program, GMT_io.ncol[GMT_IN], type[GMT_io.single_precision[GMT_IN]]);
		}
		nm = GMT_get_nm (grd.nx, grd.ny);
		a = (float *) GMT_memory (VNULL, (size_t)nm, sizeof (float), GMT_program);
		GMT_err_fail (GMT_read_grd (grdfile, &grd, a, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), grdfile);
		if ((fp = GMT_fopen (Ctrl->N.file, GMT_io.r_mode)) == NULL) {
			fprintf (stderr, "%s: Could not open file %s\n", GMT_program, Ctrl->N.file);
			exit (EXIT_FAILURE);
		}
		if (GMT_io.io_header[GMT_IN]) for (i = 0; i < GMT_io.n_header_recs; i++) not_used = GMT_fgets (buffer, BUFSIZ, fp);

		n_expected_fields = (GMT_io.ncol[GMT_IN]) ? GMT_io.ncol[GMT_IN] : 3;

		n_data = 0;
		while ((n_fields = GMT_input (fp, &n_expected_fields, &in)) >= 0 && !(GMT_io.status & GMT_IO_EOF)) {	/* Not yet EOF */
			n_data++;
			if (GMT_io.status & GMT_IO_MISMATCH) {
				fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld (skipped)\n", GMT_program, n_fields,  n_expected_fields, n_data);
				continue;
			}
			if (GMT_y_is_outside (in[GMT_Y],  grd.y_min, grd.y_max)) continue;	/* Outside y-range */
			if (GMT_x_is_outside (&in[GMT_X], grd.x_min, grd.x_max)) continue;	/* Outside x-range */
			j = GMT_y_to_j(in[GMT_Y], grd.y_min, grd.y_inc, grd.xy_off, grd.ny);
			if (j < 0 || j >= grd.ny) continue;
			i = GMT_x_to_i(in[GMT_X], grd.x_min, grd.x_inc, grd.xy_off, grd.nx);
			if (i < 0 || i >= grd.nx) continue;
			k = GMT_IJ (j, i, grd.nx);
			a[k] = (float)in[GMT_Z];
			if (GMT_io.in_col_type[0] == GMT_IS_LON && GMT_360_RANGE (grd.x_max, grd.x_min) && !grd.node_offset) {
				/* Possibly need to replicate e/w value */
				if (i == 0) {k = GMT_IJ (j, grd.nx-1, grd.nx); a[k] = (float)in[GMT_Z]; }
				if (i == (grd.nx-1)) {k = GMT_IJ (j, 0, grd.nx); a[k] = (float)in[GMT_Z]; }
			}
		}
		if (fp != GMT_stdin) GMT_fclose (fp);

		GMT_err_fail (GMT_write_grd (grdfile, &grd, a, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), grdfile);
		GMT_free ((void *)a);
	}
	else if (Ctrl->E.active) {	/* Transpose the matrix and exchange x and y info */
		struct GRD_HEADER grd_tr;
		GMT_LONG row, col;
		GMT_LONG ij, ij_tr;
		float *a_tr;
		
		a = (float *) GMT_memory (VNULL, (size_t)(grd.nx * grd.ny), sizeof (float), GMT_program);
		GMT_err_fail (GMT_read_grd (grdfile, &grd, a, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), grdfile);
		memcpy ((void *)&grd_tr, (void *)&grd, sizeof (struct GRD_HEADER));	/* First make a copy of header */
		grd_tr.nx = grd.ny;	/* Then exchange x and y values */
		grd_tr.x_min = grd.y_min;
		grd_tr.x_max = grd.y_max;
		grd_tr.x_inc = grd.y_inc;
		strcpy (grd_tr.x_units, grd.y_units);
		grd_tr.ny = grd.nx;
		grd_tr.y_min = grd.x_min;
		grd_tr.y_max = grd.x_max;
		grd_tr.y_inc = grd.x_inc;
		strcpy (grd_tr.y_units, grd.x_units);
		/* Now transpose the matrix */

		nm = GMT_get_nm (grd.nx, grd.ny);
		a_tr = (float *) GMT_memory (VNULL, (size_t)nm, sizeof (float), GMT_program);
		for (row = 0; row < grd.ny; row++) {
			for (col = 0; col < grd.nx; col++) {
				ij = GMT_IJ (row, col, grd.nx);
				ij_tr = GMT_IJ (col, row, grd.ny);
				a_tr[ij_tr] = a[ij];
			}
		}
		GMT_err_fail (GMT_write_grd (grdfile, &grd_tr, a_tr, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), grdfile);
		GMT_free ((void *)a_tr);
		GMT_free ((void *)a);
	}
	else {
		if (project_info.region_supplied) {
			if (gmtdefs.verbose) fprintf (stderr, "%s: Reset region in file %s to %g/%g/%g/%g\n",
				GMT_program, grdfile, w, e, s, n);
			grd.x_min = w;	grd.x_max = e;
			grd.y_min = s;	grd.y_max = n;
			Ctrl->A.active = TRUE;	/* Must ensure -R -I compatibility */
		}
		if (Ctrl->A.active) {
			grd.x_inc = GMT_get_inc (grd.x_min, grd.x_max, grd.nx, grd.node_offset);
			grd.y_inc = GMT_get_inc (grd.y_min, grd.y_max, grd.ny, grd.node_offset);
			if (gmtdefs.verbose) fprintf (stderr, "%s: Reset grid-spacing in file %s to %g/%g\n",
				GMT_program, grdfile, grd.x_inc, grd.y_inc);
		}
		if (Ctrl->T.active) {	/* Grid-line <---> Pixel toggling of the header */
			f = (grd.node_offset) ? -0.5 : +0.5;
			grd.x_min -= f * grd.x_inc;	grd.x_max += f * grd.x_inc;
			grd.y_min -= f * grd.y_inc;	grd.y_max += f * grd.y_inc;
			grd.node_offset = 1 - grd.node_offset;
			if (gmtdefs.verbose) fprintf (stderr, "%s: Toggled registration mode in file %s from %s to %s\n", GMT_program, grdfile, registration[1-grd.node_offset], registration[grd.node_offset]);
			if (gmtdefs.verbose) fprintf (stderr, "%s: Reset region in file %s to %g/%g/%g/%g\n", GMT_program, grdfile, grd.x_min, grd.x_max, grd.y_min, grd.y_max);
		}
		GMT_err_fail (GMT_update_grd_info (grdfile, &grd), grdfile);
	}

	if (gmtdefs.verbose) fprintf (stderr, "%s: File %s updated.\n", GMT_program, grdfile);

	Free_grdedit_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void *New_grdedit_Ctrl () {	/* Allocate and initialize a new control structure */
	struct GRDEDIT_CTRL *C;
	
	C = (struct GRDEDIT_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GRDEDIT_CTRL), "New_grdedit_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */

	return ((void *)C);
}

void Free_grdedit_Ctrl (struct GRDEDIT_CTRL *C) {	/* Deallocate control structure */
	if (C->D.information) free ((void *)C->D.information);	
	if (C->N.file) free ((void *)C->N.file);	
	GMT_free ((void *)C);	
}

