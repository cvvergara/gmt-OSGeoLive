/*--------------------------------------------------------------------
 *	$Id: grdcut.c 10173 2014-01-01 09:52:34Z pwessel $
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
 * grdcut.c reads a grid file and writes a portion within it
 * to a new file.
 *
 * Author:	Walter Smith
 * Date:	5 august, 1988
 * Version:	4
 * 		PW: 4/20/09 Added -Z option
 *
 */
 
#include "gmt.h"

struct GRDCUT_CTRL {
	struct G {	/* -G<output_grdfile> */
		GMT_LONG active;
		char *file;
	} G;
	struct Z {	/* -Z[min/max] */
		GMT_LONG active;
		GMT_LONG mode;	/* 1 means NaN */
		double min, max;
	} Z;
};

#define NAN_IS_INSIDE	0
#define NAN_IS_OUTSIDE	1

int main (int argc, char **argv)
{
	GMT_LONG	error = FALSE;

	GMT_LONG	i, k, nx_old, ny_old, nx_new, ny_new, one_or_zero;
	
	GMT_LONG	nm;

	char *grd_in, format[BUFSIZ], za[GMT_TEXT_LEN], zb[GMT_TEXT_LEN];

	float	*grd = NULL;

	double	w_new = 0.0, e_new = 0.0, s_new = 0.0, n_new = 0.0;
	double	w_old, e_old, s_old, n_old;

	struct GRD_HEADER header, test_header;
	struct GRDCUT_CTRL *Ctrl = NULL;

	void *New_grdcut_Ctrl (), Free_grdcut_Ctrl (struct GRDCUT_CTRL *C);
	
	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct GRDCUT_CTRL *) New_grdcut_Ctrl ();	/* Allocate and initialize a new control structure */

	grd_in = CNULL;
	za[0] = zb[0] = '\0';

	/* Check and interpret the command line arguments */

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch(argv[i][1]) {
				/* Common parameters */

				case 'R':
				case 'V':
				case 'f':
				case '\0':
					error += GMT_parse_common_options (argv[i], &w_new, &e_new, &s_new, &n_new);
					break;

	 			case 'G':
					Ctrl->G.active = TRUE;
					Ctrl->G.file = strdup (&argv[i][2]);
					break;
					
	 			case 'Z':
					Ctrl->Z.active = TRUE;
					k = 2;
					if (argv[i][k] == 'n') {
						Ctrl->Z.mode = NAN_IS_OUTSIDE;
						k = 3;
					}
					if (sscanf (&argv[i][k], "%[^/]/%s", za, zb) == 2) {
						if (!(za[0] == '-' && za[1] == '\0')) error += GMT_verify_expectations (GMT_io.in_col_type[GMT_Z], GMT_scanf_arg (za, GMT_io.in_col_type[GMT_Z], &Ctrl->Z.min), za);
						if (!(zb[0] == '-' && zb[1] == '\0')) error += GMT_verify_expectations (GMT_io.in_col_type[GMT_Z], GMT_scanf_arg (zb, GMT_io.in_col_type[GMT_Z], &Ctrl->Z.max), zb);
					}
					break;
	
				default:		/* Options not recognized */
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else
	 		grd_in = argv[i];
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr,"grdcut %s - Extract subsets from grid files\n\n", GMT_VERSION);
		fprintf (stderr, "usage: grdcut <input_grd> -G<output_grd> %s [-V] [-Z[n][min/max]] [%s]\n", GMT_Rgeo_OPT, GMT_f_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\t<input_grd> is file to extract a subset from.\n");
		fprintf (stderr, "\t-G specifies output grid file.\n");
		GMT_explain_option ('R');
		fprintf (stderr, "\t   Obviously, the WESN you specify must be within the WESN of the input file.\n");
		fprintf (stderr, "\t   If in doubt, run grdinfo first and check range of old file.\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		GMT_explain_option ('V');
		fprintf (stderr, "\t-Z Specify a range and determine the corresponding rectangular region so that\n");
		fprintf (stderr, "\t   all values outside this region are outside the range [-inf/+inf].\n");
		fprintf (stderr, "\t   Use -Zn to consider NaNs outside as well [Default just ignores NaNs].\n");
		GMT_explain_option ('f');
		exit (EXIT_FAILURE);
	}

	/* Check that the options selected make sense */

	if ((project_info.region_supplied + Ctrl->Z.active) != 1) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify either the -R or the -Z options\n", GMT_program);
		error++;
	}
	if (!Ctrl->G.file) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -G option:  Must specify output file\n", GMT_program);
		error++;
	}
	if (!grd_in) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify input file\n", GMT_program);
		error++;
	}
	if (error) exit (EXIT_FAILURE);

	GMT_err_fail (GMT_read_grd_info (grd_in, &header), grd_in);

	if (Ctrl->Z.active) {	/* Must determine new region via -Z */
		GMT_LONG i0, i1, j0, j1, j, ij;
		
		nm = GMT_get_nm (header.nx, header.ny);
		grd = (float *) GMT_memory (VNULL, (size_t)nm, sizeof (float), GMT_program);
		GMT_err_fail (GMT_read_grd (grd_in, &header, grd, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), grd_in);
		for (i = 0, i0 = -1; i0 == -1 && i < header.nx; i++) {	/* Scan from xmin towards xmax */
			for (j = 0, ij = i; i0 == -1 && j < header.ny; j++, ij += header.ny) {
				if (GMT_is_fnan (grd[ij])) {
					if (Ctrl->Z.mode == NAN_IS_OUTSIDE) i0 = i;	/* Must stop since this value defines the inner box */
				}
				else if (grd[ij] >= Ctrl->Z.min && grd[ij] <= Ctrl->Z.max)
					i0 = i;
			}
		}
		if (i0 == -1) {
			fprintf (stderr, "%s: The sub-region is empty - no file written\n", GMT_program);
			GMT_free ((void *)grd);
			exit (EXIT_FAILURE);
		}
		for (i = header.nx-1, i1 = -1; i1 == -1 && i > i0; i--) {	/* Scan from xmax towards xmin */
			for (j = 0, ij = i; i1 == -1 && j < header.ny; j++, ij += header.ny) {
				if (GMT_is_fnan (grd[ij])) {
					if (Ctrl->Z.mode == NAN_IS_INSIDE) i1 = i;	/* Must stop since this value defines the inner box */
				}
				else if (grd[ij] >= Ctrl->Z.min && grd[ij] <= Ctrl->Z.max)
					i1 = i;
			}
		}
		for (j = 0, j0 = -1; j0 == -1 && j < header.ny; j++) {	/* Scan from ymin towards ymax */
			for (i = i0, ij = GMT_IJ(j,i0,header.nx); j0 == -1 && i < i1; i++, ij++) {
				if (GMT_is_fnan (grd[ij])) {
					if (Ctrl->Z.mode == NAN_IS_INSIDE) j0 = j;	/* Must stop since this value defines the inner box */
				}
				else if (grd[ij] >= Ctrl->Z.min && grd[ij] <= Ctrl->Z.max)
					j0 = j;
			}
		}
		for (j = header.ny-1, j1 = -1; j1 == -1 && j >= j0; j--) {	/* Scan from ymax towards ymin */
			for (i = i0, ij = GMT_IJ(j,i0,header.nx); j1 == -1 && i < i1; i++, ij++) {
				if (GMT_is_fnan (grd[ij])) {
					if (Ctrl->Z.mode == NAN_IS_INSIDE) j1 = j;	/* Must stop since this value defines the inner box */
				}
				else if (grd[ij] >= Ctrl->Z.min && grd[ij] <= Ctrl->Z.max)
					j1 = j;
			}
		}
		if (i0 == 0 && j0 == 0 && i1 == (header.nx-1) && j1 == (header.ny-1)) {
			fprintf (stderr, "%s: Your -Z limits produced no subset - output grid is identical to input grid\n", GMT_program);
			w_new = header.x_min;	e_new = header.x_max;
			s_new = header.y_min;	n_new = header.y_max;
		}
		else {	/* Adjust boundaries inwards */
			w_new = header.x_min + i0 * header.x_inc;
			e_new = header.x_max - (header.nx - 1 - i1) * header.x_inc;
			s_new = header.y_min + (header.ny - 1 - j1) * header.y_inc;
			n_new = header.y_max - j0 * header.y_inc;
		}
		GMT_free ((void *)grd);
	}
	
	if (s_new < header.y_min || s_new > header.y_max) error = TRUE;
	if (n_new < header.y_min || n_new > header.y_max) error = TRUE;

	if (GMT_io.in_col_type[GMT_X] == GMT_IS_LON) {	/* Geographic data */
		if (w_new < header.x_min && e_new < header.x_min) {
			header.x_min -= 360.0;
			header.x_max -= 360.0;
		}
		if (w_new > header.x_max && e_new > header.x_max) {
			header.x_min += 360.0;
			header.x_max += 360.0;
		}
		if (!GMT_grd_is_global (&header) && (w_new < header.x_min || e_new > header.x_max)) error = TRUE;
	}
	else if (w_new < header.x_min || e_new > header.x_max)
		error = TRUE;

	if (error) {
		fprintf (stderr, "%s: Subset exceeds data domain!\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	/* Make sure output grid is kosher */

	GMT_adjust_loose_wesn (&w_new, &e_new, &s_new, &n_new, &header);

	test_header.x_min = w_new;	test_header.x_max = e_new;	test_header.x_inc = header.x_inc;
	test_header.y_min = s_new;	test_header.y_max = n_new;	test_header.y_inc = header.y_inc;
	GMT_err_fail (GMT_grd_RI_verify (&test_header, 1), Ctrl->G.file);

	/* OK, so far so good. Check if new wesn differs from old wesn by integer dx/dy */

	if (GMT_minmaxinc_verify (header.x_min, w_new, header.x_inc, GMT_SMALL) == 1) {
		fprintf (stderr, "%s: Old and new x_min do not differ by N * dx\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (GMT_minmaxinc_verify (e_new, header.x_max, header.x_inc, GMT_SMALL) == 1) {
		fprintf (stderr, "%s: Old and new x_max do not differ by N * dx\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (GMT_minmaxinc_verify (header.y_min, s_new, header.y_inc, GMT_SMALL) == 1) {
		fprintf (stderr, "%s: Old and new y_min do not differ by N * dy\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (GMT_minmaxinc_verify (n_new, header.y_max, header.y_inc, GMT_SMALL) == 1) {
		fprintf (stderr, "%s: Old and new y_max do not differ by N * dy\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	GMT_grd_init (&header, argc, argv, TRUE);

	w_old = header.x_min;	e_old = header.x_max;
	s_old = header.y_min;	n_old = header.y_max;
	nx_old = header.nx;	ny_old = header.ny;
	one_or_zero = (header.node_offset) ? 0 : 1;
	nx_new = irint ((e_new - w_new) / header.x_inc) + one_or_zero;
	ny_new = irint ((n_new - s_new) / header.y_inc) + one_or_zero;

	nm = GMT_get_nm (nx_new, ny_new);
	
	grd = (float *) GMT_memory (VNULL, (size_t)nm, sizeof (float), GMT_program);
	GMT_err_fail (GMT_read_grd (grd_in, &header, grd, w_new, e_new, s_new, n_new, GMT_pad, FALSE), grd_in);

	if (gmtdefs.verbose) {
		sprintf (format, "\t%s\t%s\t%s\t%s\t%s\t%s\t%%ld\t%%ld\n", gmtdefs.d_format, gmtdefs.d_format,
			gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format);
		fprintf (stderr, "%s: File spec:\tW E S N dx dy nx ny:\n", GMT_program);
		fprintf (stderr, "%s: Old:", GMT_program);
		fprintf (stderr, format, w_old, e_old, s_old, n_old, header.x_inc, header.y_inc, nx_old, ny_old);
		fprintf (stderr, "%s: New:", GMT_program);
		fprintf (stderr, format, w_new, e_new, s_new, n_new, header.x_inc, header.y_inc, nx_new, ny_new);
	}

	GMT_err_fail (GMT_write_grd (Ctrl->G.file, &header, grd, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), Ctrl->G.file);

	GMT_free ((void *) grd);

	Free_grdcut_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void *New_grdcut_Ctrl () {	/* Allocate and initialize a new control structure */
	struct GRDCUT_CTRL *C;
	
	C = (struct GRDCUT_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GRDCUT_CTRL), "New_grdcut_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
			
	C->Z.min = -DBL_MAX;	C->Z.max = DBL_MAX;			/* No limits on z-range */
	return ((void *)C);
}

void Free_grdcut_Ctrl (struct GRDCUT_CTRL *C) {	/* Deallocate control structure */
	if (C->G.file) free ((void *)C->G.file);	
	GMT_free ((void *)C);	
}
