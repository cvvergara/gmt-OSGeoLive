/*--------------------------------------------------------------------
 *	$Id: grdclip.c 10173 2014-01-01 09:52:34Z pwessel $
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
   grdclip read a grid file and sets all values < the user-supplied
   lower limit to the value <below>, and all values > the user-
   supplied upper limit to the value <above>.  above/below can
   be any number including NaN.

   WHF Smith, 22 Aug 88
   Modified for v2.0 4/26/91 P.Wessel
   Modified for v3.1 4/29/98 P.Wessel
   Modified for v3.3.5 06/29/2000 P.Wessel
   Version:	4
 
 */

#include "gmt.h"

struct GRDCLIP_CTRL {
	struct G {	/* -G<output_grdfile> */
		GMT_LONG active;
		char *file;
	} G;
	struct SA {	/* -Sa<high/above> */
		GMT_LONG active;
		float high, above;
	} Sa;
	struct SB {	/* -Sb<low/below> */
		GMT_LONG active;
		float low, below;
	} Sb;
};

int main (int argc, char **argv)
{
	GMT_LONG	error;
	char *infile = CNULL, format[BUFSIZ], txt[50];
	GMT_LONG i;
	GMT_LONG k, nxy, n_above, n_below, n;
	float	*data = NULL;
	struct GRD_HEADER header;
	struct GRDCLIP_CTRL *Ctrl = NULL;

	void *New_grdclip_Ctrl (), Free_grdclip_Ctrl (struct GRDCLIP_CTRL *C);
	
	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct GRDCLIP_CTRL *) New_grdclip_Ctrl ();	/* Allocate and initialize a new control structure */
	
	error = FALSE;
	n_above = n_below = 0;

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				/* Common parameters */

				case 'V':
				case '\0':
					error += GMT_parse_common_options (argv[i], 0, 0, 0, 0);
					break;

				/* Supplemental parameters */

				case 'S':
					if (argv[i][2] == 'a') {
						Ctrl->Sa.active = TRUE;
						n = sscanf (&argv[i][3], "%f/%s", &Ctrl->Sa.high, txt);
						if (n != 2) {
							fprintf (stderr, "%s: GMT SYNTAX ERROR -Sa option.  Correct syntax:\n", GMT_program);
							fprintf (stderr, "\t-Sa<high>/<above>, <above> may be set to NaN\n");
							error = TRUE;
						}
						else 
							Ctrl->Sa.above = (txt[0] == 'N' || txt[0] == 'n') ? GMT_f_NaN : (float)atof (txt);
					}
					else if (argv[i][2] == 'b') {
						Ctrl->Sb.active = TRUE;
						n = sscanf (&argv[i][3], "%f/%s", &Ctrl->Sb.low, txt);
						if (n != 2) {
							fprintf (stderr, "%s: GMT SYNTAX ERROR -Sb option.  Correct syntax:\n", GMT_program);
							fprintf (stderr, "\t-Sb<low>/<below>, <below> may be set to NaN\n");
							error = TRUE;
						}
						else
							Ctrl->Sb.below = (txt[0] == 'N' || txt[0] == 'n') ? GMT_f_NaN : (float)atof (txt);
					}
					else {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -S option.  Correct syntax:\n", GMT_program);
						fprintf (stderr, "\t-Sa<high>/<above> or -Sb<low>/<below>\n");
						error = TRUE;
					}
					break;
				case 'G':
					Ctrl->G.active = TRUE;
					Ctrl->G.file = strdup (&argv[i][2]);
					break;

				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else
			infile = argv[i];
	}

	if (GMT_give_synopsis_and_exit || argc == 1) {
		fprintf (stderr, "grdclip %s - Clipping of range in grid files\n\n", GMT_VERSION);
		fprintf (stderr, "usage: grdclip <input_grdfile> -G<output_grdfile>\n");
		fprintf (stderr, "\t[-Sa<high/above>] [-Sb<low/below>] [-V]\n");

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\tYou must choose at least one -S option.\n");
		fprintf (stderr, "\t-G name of output grid.\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-Sa will set all data > high to the <above> value.\n");
		fprintf (stderr, "\t-Sb will set all data < low to the <below> value.\n");
		fprintf (stderr, "\t    <above> and <below> can be any number, including NaN.\n");
		GMT_explain_option ('V');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	if (!Ctrl->G.file) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -G option:.  Must specify output file\n", GMT_program);
		error++;
	}
	if (!infile) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify input file\n", GMT_program);
		error++;
	}
	if (!(Ctrl->Sa.active || Ctrl->Sb.active)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -S option:  Must specify at least one of -Sa -Sb\n", GMT_program);
		error++;
	}

	if (error) exit (EXIT_FAILURE);

	GMT_err_fail (GMT_read_grd_info (infile, &header), infile);

	GMT_grd_init (&header, argc, argv, TRUE);

	nxy = GMT_get_nm (header.nx, header.ny);
	
	data = (float *) GMT_memory (VNULL, (size_t)nxy, sizeof (float), GMT_program);

	GMT_err_fail (GMT_read_grd (infile, &header, data, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), infile);

	for (k = 0; k < nxy; k++) {
		if (Ctrl->Sa.active && (!GMT_is_fnan (data[k])) && data[k] > Ctrl->Sa.high) {
			data[k] = Ctrl->Sa.above;
			n_above++;
		}
		else if (Ctrl->Sb.active && (!GMT_is_fnan (data[k])) && data[k] < Ctrl->Sb.low) {
			data[k] = Ctrl->Sb.below;
			n_below++;
		}
	}

	GMT_err_fail (GMT_write_grd (Ctrl->G.file, &header, data, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), Ctrl->G.file);

	GMT_free ((void *) data);

	if (gmtdefs.verbose) {
		sprintf (format, "%s set to %s\n", gmtdefs.d_format, gmtdefs.d_format);
		if (Ctrl->Sb.active) {
			fprintf (stderr, "%s: %ld values < ", GMT_program, n_below);
			fprintf (stderr, format, (double)Ctrl->Sb.low, (double)Ctrl->Sb.below);
		}
		if (Ctrl->Sa.active) {
			fprintf (stderr, "%s: %ld values > ", GMT_program, n_above);
			fprintf (stderr, format, (double)Ctrl->Sa.high, (double)Ctrl->Sa.above);
		}
	}

	Free_grdclip_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void *New_grdclip_Ctrl () {	/* Allocate and initialize a new control structure */
	struct GRDCLIP_CTRL *C;
	
	C = (struct GRDCLIP_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GRDCLIP_CTRL), "New_grdclip_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
			
	return ((void *)C);
}

void Free_grdclip_Ctrl (struct GRDCLIP_CTRL *C) {	/* Deallocate control structure */
	if (C->G.file) free ((void *)C->G.file);	
	GMT_free ((void *)C);	
}
