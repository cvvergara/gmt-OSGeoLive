/*--------------------------------------------------------------------
 *	$Id: grdpaste.c 10173 2014-01-01 09:52:34Z pwessel $
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
 * grdpaste.c reads two grid files and writes a new file with
 * the first two pasted together along their common edge.
 *
 * Author:	Walter Smith
 * Date:	5 august, 1988
 * Updated to v2.0 20-May-1991 Paul Wessel
 * Updated to v3.1 23-May-1998 Paul Wessel
 * 3.3.5 23-JUN-2000 Paul Wessel
 * Version:	4
 *
 */
 
#include "gmt.h"

struct GRDPASTE_CTRL {
	struct G {	/* -G<output_grdfile> */
		GMT_LONG active;
		char *file;
	} G;
};

int main (int argc, char **argv)
{
	GMT_LONG	error = FALSE;

	char *grd_a = NULL, *grd_b = NULL, format[BUFSIZ];

	GMT_LONG	i, way, one_or_zero, n_in = 0;
	
	GMT_LONG nm;

	float	*c = NULL;

	double x_noise, y_noise, w, e;

	struct GRD_HEADER ha, hb, hc;
	struct GRDPASTE_CTRL *Ctrl = NULL;

	void *New_grdpaste_Ctrl (), Free_grdpaste_Ctrl (struct GRDPASTE_CTRL *C);
	
	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct GRDPASTE_CTRL *) New_grdpaste_Ctrl ();	/* Allocate and initialize a new control structure */

	grd_a = grd_b = CNULL;

	/* Check and interpret the command line arguments */

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch(argv[i][1]) {
				/* Common parameters */

				case 'V':
				case 'f':
				case '\0':
					error += GMT_parse_common_options (argv[i], 0, 0, 0, 0);
					break;

				/* Supplemental parameters */

	 			case 'G':
					Ctrl->G.active = TRUE;
					Ctrl->G.file = strdup (&argv[i][2]);
					break;
				default:		/* Options not recognized */
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else if (n_in == 0) {
	 		grd_a = argv[i];
			n_in++;
		}
		else if (n_in == 1) {
	 		grd_b = argv[i];
			n_in++;
		}
		else {
			error = TRUE;
			fprintf (stderr, "%s: GMT SYNTAX ERROR:  Only two files may be pasted\n", GMT_program);
		}
	}

	if (GMT_give_synopsis_and_exit || argc == 1) {
		fprintf (stderr,"grdpaste %s - Join two grid files along common edge\n\n", GMT_VERSION);
		fprintf (stderr, "usage: grdpaste <file_a> <file_b> -G<outfile> [-V] [%s]\n\n", GMT_f_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "	 where file_a and file_b are to be combined into outfile.\n");
		fprintf (stderr, "	 file_a and file_b must have same dx,dy and one edge in common.\n");
		fprintf (stderr, "	 If in doubt, run grdinfo first and check your files.\n");
		fprintf (stderr, "	 Use grdcut and/or grdsample to adjust files as necessary.\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		GMT_explain_option ('V');
		GMT_explain_option ('f');
		exit (EXIT_FAILURE);
	}

	if (!grd_a || !grd_b) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify two input files\n", GMT_program);
		error++;
	}
	if (!Ctrl->G.file) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -G:  Must specify output file\n", GMT_program);
		error++;
	}
	if (error) exit (EXIT_FAILURE);

	/* Check that the options selected make sense */


	/* Try to find a common side to join on  */

	GMT_grd_init (&hc, argc, argv, FALSE);
	GMT_err_fail (GMT_read_grd_info (grd_a, &ha), grd_a);
	GMT_err_fail (GMT_read_grd_info (grd_b, &hb), grd_b);

	if (ha.node_offset != hb.node_offset) error = TRUE;
	if ((ha.z_scale_factor != hb.z_scale_factor) || (ha.z_add_offset != hb.z_add_offset)) {
		fprintf (stderr, "%s: Scale/offset not compatible!\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	if ( fabs (ha.x_inc - hb.x_inc) < 1.0e-6 && fabs (ha.y_inc - hb.y_inc) < 1.0e-6) {
		hc.x_inc = ha.x_inc;
		hc.y_inc = ha.y_inc;
	}
	else {
		fprintf (stderr, "%s:  Grid intervals do not match!\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	one_or_zero = 1 - ha.node_offset;
	x_noise = GMT_SMALL * hc.x_inc;
	y_noise = GMT_SMALL * hc.y_inc;

	w = hb.x_min;	e = hb.x_max;	/* Get copies */
	if (GMT_io.in_col_type[GMT_IN] == GMT_IS_LON) {	/* Must be careful in determining a match */
		w -= 360.0;	e -= 360.0;	/* Wind to left first */
		while ((ha.x_min - w) > x_noise) {	/* Wind back to match grid A if possible */
			w += 360.0;
			e += 360.0;
		}
		hb.x_min = w;	hb.x_max = e;	/* Possibly update lon range */
	}
	if ( fabs (ha.x_min - hb.x_min) < x_noise && fabs (ha.x_max - hb.x_max) < x_noise ) {

		if (fabs (ha.y_max - hb.y_min) < y_noise) {
			way = 1;
			hc.nx = ha.nx;
			hc.ny = (int)(ha.ny + hb.ny - one_or_zero);
			hc.x_min = ha.x_min;
			hc.x_max = ha.x_max;
			hc.y_min = ha.y_min;
			hc.y_max = hb.y_max;
		}
		else if (fabs (ha.y_min - hb.y_max) < y_noise) {
			way = 2;
			hc.nx = ha.nx;
			hc.ny = (int)(ha.ny + hb.ny - one_or_zero);
			hc.x_min = ha.x_min;
			hc.x_max = ha.x_max;
			hc.y_min = hb.y_min;
			hc.y_max = ha.y_max;
		}
		else {
			fprintf (stderr, "%s:  Grids do not share a common edge!\n", GMT_program);
			exit (EXIT_FAILURE);
		}
	}
	else if ( fabs (ha.y_min - hb.y_min) < y_noise && fabs (ha.y_max - hb.y_max) < y_noise ) {

		if (fabs (ha.x_min - hb.x_max) < x_noise) {
			way = 3;
			hc.nx = (int)(ha.nx + hb.nx - one_or_zero);
			hc.ny = ha.ny;
			hc.x_min = hb.x_min;
			hc.x_max = ha.x_max;
			hc.y_min = ha.y_min;
			hc.y_max = ha.y_max;
		}
		else if (fabs (ha.x_max - hb.x_min) < x_noise) {
			way = 4;
			hc.nx = (int)(ha.nx + hb.nx - one_or_zero);
			hc.ny = ha.ny;
			hc.x_min = ha.x_min;
			hc.x_max = hb.x_max;
			hc.y_min = ha.y_min;
			hc.y_max = ha.y_max;
		}
		else {
			fprintf (stderr, "%s:  Grids do not share a common edge!\n", GMT_program);
			exit (EXIT_FAILURE);
		}
	}
	else {
		fprintf (stderr, "%s:  Grids do not share a common edge!\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	if (GMT_io.in_col_type[GMT_IN] == GMT_IS_LON && hc.x_max > 360.0) {	/* Take out 360.0 */
		hc.x_min -= 360.0;
		hc.x_max -= 360.0;
	}
	/* Now we can do it  */

	sprintf (format, "%%s: \t%s\t%s\t%s\t%s\t%s\t%s\t%%ld\t%%ld\n", gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format);
	if (gmtdefs.verbose) {
		fprintf (stderr, "File spec:\tW E S N dx dy nx ny:\n");
		fprintf (stderr, format, grd_a, ha.x_min, ha.x_max, ha.y_min, ha.y_max, ha.x_inc, ha.y_inc, ha.nx, ha.ny);
		fprintf (stderr, format, grd_b, hb.x_min, hb.x_max, hb.y_min, hb.y_max, hb.x_inc, hb.y_inc, hb.nx, hb.ny);
		fprintf (stderr, format, Ctrl->G.file, hc.x_min, hc.x_max, hc.y_min, hc.y_max, hc.x_inc, hc.y_inc, hc.nx, hc.ny);
	}

	nm = GMT_get_nm (hc.nx, hc.ny);
	c = (float *) GMT_memory (VNULL, (size_t)nm, sizeof (float), GMT_program);
	hc.node_offset = ha.node_offset;
	hc.z_scale_factor = ha.z_scale_factor;
	hc.z_add_offset = ha.z_add_offset;

	switch (way) {
		case 1:
			GMT_pad[3] = hb.ny - one_or_zero;
			GMT_err_fail (GMT_read_grd (grd_a, &ha, c, ha.x_min, ha.x_max, ha.y_min, ha.y_max, GMT_pad, FALSE), grd_a);
			GMT_pad[3] = 0;	GMT_pad[2] = ha.ny - one_or_zero;
			GMT_err_fail (GMT_read_grd (grd_b, &hb, c, hb.x_min, hb.x_max, hb.y_min, hb.y_max, GMT_pad, FALSE), grd_b);
			break;
		case 2:
			GMT_pad[2] = hb.ny - one_or_zero;
			GMT_err_fail (GMT_read_grd (grd_a, &ha, c, ha.x_min, ha.x_max, ha.y_min, ha.y_max, GMT_pad, FALSE), grd_a);
			GMT_pad[2] = 0;	GMT_pad[3] = ha.ny - one_or_zero;
			GMT_err_fail (GMT_read_grd (grd_b, &hb, c, hb.x_min, hb.x_max, hb.y_min, hb.y_max, GMT_pad, FALSE), grd_b);
			break;
		case 3:
			GMT_pad[0] = hb.nx - one_or_zero;
			GMT_err_fail (GMT_read_grd (grd_a, &ha, c, ha.x_min, ha.x_max, ha.y_min, ha.y_max, GMT_pad, FALSE), grd_a);
			GMT_pad[0] = 0;	GMT_pad[1] = ha.nx - one_or_zero;
			GMT_err_fail (GMT_read_grd (grd_b, &hb, c, hb.x_min, hb.x_max, hb.y_min, hb.y_max, GMT_pad, FALSE), grd_b);
			break;
		case 4:
			GMT_pad[1] = hb.nx - one_or_zero;
			GMT_err_fail (GMT_read_grd (grd_a, &ha, c, ha.x_min, ha.x_max, ha.y_min, ha.y_max, GMT_pad, FALSE), grd_a);
			GMT_pad[1] = 0;	GMT_pad[0] = ha.nx - one_or_zero;
			GMT_err_fail (GMT_read_grd (grd_b, &hb, c, hb.x_min, hb.x_max, hb.y_min, hb.y_max, GMT_pad, FALSE), grd_b);
			break;
	}

	GMT_pad[0] = GMT_pad[1] = GMT_pad[2] = GMT_pad[3] = 0;
	
	GMT_err_fail (GMT_write_grd (Ctrl->G.file, &hc, c, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), Ctrl->G.file);

	GMT_free ((void *) c);

	Free_grdpaste_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void *New_grdpaste_Ctrl () {	/* Allocate and initialize a new control structure */
	struct GRDPASTE_CTRL *C;
	
	C = (struct GRDPASTE_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GRDPASTE_CTRL), "New_grdpaste_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
			
	return ((void *)C);
}

void Free_grdpaste_Ctrl (struct GRDPASTE_CTRL *C) {	/* Deallocate control structure */
	if (C->G.file) free ((void *)C->G.file);	
	GMT_free ((void *)C);	
}
