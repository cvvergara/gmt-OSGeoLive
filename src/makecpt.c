/*--------------------------------------------------------------------
 *	$Id: makecpt.c 9923 2012-12-18 20:45:53Z pwessel $
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
 *
 * Read an existing cpt table and desired z grid and produce
 * a GMT cpt file.  Can be inverted [-I] or made continuous [-Z].
 * Discrete color jumps in cpt tables are handled correctly.
 * Default color table is still rainbow.
 *
 * Author:	Walter H.f. Smith & P. Wessel
 * Date:	22-SEP-2000
 * Version:	4
 */

#include "gmt.h"

struct MAKECPT_CTRL {
	struct C {	/* -C<cpt> */
		GMT_LONG active;
		char *file;
	} C;
	struct D {	/* -D */
		GMT_LONG active;
	} D;
	struct I {	/* -I */
		GMT_LONG active;
	} I;
	struct M {	/* -M */
		GMT_LONG active;
	} M;
	struct N {	/* -N */
		GMT_LONG active;
	} N;
	struct T {	/* -T<z0/z1/dz> */
		GMT_LONG active;
		GMT_LONG cpt;
		double low, high, inc;
		char *file;
	} T;
	struct Q {	/* -Q[i|o */
		GMT_LONG active;
		GMT_LONG mode;
	} Q;
	struct Z {	/* -Z */
		GMT_LONG active;
	} Z;
};

int main(int argc, char **argv)
{
	GMT_LONG	i, nz;

	GMT_LONG	error = FALSE;

	double	*z = NULL;

	FILE	*fpc = NULL, *fpl = NULL;

	char	buffer[BUFSIZ], CPT_lis[BUFSIZ], CPT_file[BUFSIZ];

	struct MAKECPT_CTRL *Ctrl = NULL;

	void *New_makecpt_Ctrl (), Free_makecpt_Ctrl (struct MAKECPT_CTRL *C);

	argc = (int)GMT_begin(argc, argv);

	Ctrl = (struct MAKECPT_CTRL *)New_makecpt_Ctrl ();	/* Allocate and initialize a new control structure */
	
	/* Get list of available color tables in $GMT_SHAREDIR */

	GMT_getsharepath ("conf", "gmt_cpt", ".conf", CPT_lis);
	if ((fpc = fopen (CPT_lis, "r")) == NULL) {
		fprintf (stderr, "%s: ERROR: Cannot open file %s\n", GMT_program, CPT_lis);
		exit (EXIT_FAILURE);
	}

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch(argv[i][1]) {
				case 'V':
				case '\0':
					error += GMT_parse_common_options(argv[i], 0, 0, 0, 0);
					break;
				case 'C':
					Ctrl->C.active = TRUE;
					Ctrl->C.file = strdup (&argv[i][2]);
					break;
				case 'D':
					Ctrl->D.active = TRUE;
					break;
				case 'I':
					Ctrl->I.active = TRUE;
					break;
				case 'M':
					Ctrl->M.active = TRUE;
					break;
				case 'N':
					Ctrl->N.active = TRUE;
					break;
				case 'T':
					Ctrl->T.active = TRUE;
					if (!access (&argv[i][2], R_OK))
						Ctrl->T.file = strdup (&argv[i][2]);
					else
						sscanf (&argv[i][2], "%lf/%lf/%lf", &Ctrl->T.low, &Ctrl->T.high, &Ctrl->T.inc);
					break;
				case 'Q':
					Ctrl->Q.active = TRUE;
					if (argv[i][2] == 'o')	/* Input data is z, but take log10(z) before interpolation colors */
						Ctrl->Q.mode = 2;
					else			/* Input is log10(z) */
						Ctrl->Q.mode = 1;
					break;
				case 'Z':
					Ctrl->Z.active = TRUE;
					break;
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
		fprintf (stderr, "makecpt %s - Make GMT color palette tables\n\n", GMT_VERSION);
		fprintf (stderr, "usage:  makecpt [-C<table>] [-D] [-I] [-M] [-N] [-Q[i|o]] [-T<z0/z1/dz> | -T<file>] [-V] [-Z]\n");
		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-C Specify a colortable [Default is rainbow]:\n");
		fprintf (stderr, "\t   [Default min/max values for -T are given in brackets].\n");
		fprintf (stderr, "\t   ---------------------------------\n");
		while (fgets (buffer, BUFSIZ, fpc)) if (!(buffer[0] == '#' || buffer[0] == 0)) fprintf (stderr, "\t   %s", buffer);
		fclose (fpc);
		fprintf (stderr, "\t   ---------------------------------\n");
		fprintf (stderr, "\t-D Set back- and foreground color to match the bottom/top limits in the cpt file [Default uses color table].\n");
		fprintf (stderr, "\t-I Reverses the sense of the color table as well as back- and foreground color.\n");
		fprintf (stderr, "\t-M Use GMT defaults to set back-, foreground, and NaN colors [Default uses color table].\n");
		fprintf (stderr, "\t-N Do not write back-, foreground, and NaN colors [Default will].\n");
		fprintf (stderr, "\t-Q Assign a logarithmic colortable [Default is linear].\n");
		fprintf (stderr, "\t   -Qi: z-values are actually log10(z). Assign colors and write z [Default].\n");
		fprintf (stderr, "\t   -Qo: z-values are z, but take log10(z), assign colors and write z.\n");
		fprintf (stderr, "\t        If -T<z0/z1/dz> is given, dz is 1, 2, or 3 (as in logarithmic annotations).\n");
		fprintf (stderr, "\t-T Give start, stop, and increment for colorscale in z-units, or filename with custom z-values.\n");
		fprintf (stderr, "\t   If not given, the range in the master cptfile is used.\n");
		GMT_explain_option ('V');
		fprintf (stderr, "\t-Z Create a continuous color palette [Default is discontinuous, i.e., constant color intervals].\n");
		exit (EXIT_FAILURE);
	}

	fclose (fpc);

	if (Ctrl->T.active && !Ctrl->T.file && (Ctrl->T.low >= Ctrl->T.high || Ctrl->T.inc <= 0.0)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -T option:  Give start < stop and inc > 0.0\n", GMT_program);
		error++;
	}
	if (!Ctrl->C.active) {	/* Set default table */
		Ctrl->C.active = TRUE;
		Ctrl->C.file = strdup ("rainbow");
	}	
	error += GMT_set_cpt_path (CPT_file, Ctrl->C.file);

	if (Ctrl->T.file && (fpl = GMT_fopen (Ctrl->T.file, "r")) == NULL) {
		fprintf (stderr, "%s: Error: Could not open file %s\n", GMT_program, Ctrl->T.file);
		error++;
	}

	if (error) exit (EXIT_FAILURE);

	/* OK, we can now do the resampling */

	if (Ctrl->N.active) GMT_cpt_flags += 1;	/* bit 0 controls if BFN will be written out */
	if (Ctrl->D.active) GMT_cpt_flags += 2;	/* bit 1 controls if BF will be set to equal bottom/top rgb value */
	if (Ctrl->M.active) GMT_cpt_flags += 4;	/* bit 2 controls if BFN is determined by parameters */

	GMT_read_cpt (CPT_file);

	/* Set up arrays */

	if (Ctrl->T.file) {
		int n_alloc = GMT_SMALL_CHUNK;
		z = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof(double), GMT_program);
		nz = 0;
		while (GMT_fgets (buffer, BUFSIZ, fpl)) {
			if (GMT_is_a_blank_line (buffer)) continue;	/* Skip blank lines or # comments */
			GMT_chop (buffer);
			z[nz] = atof (buffer);
			nz++;
			if (nz == n_alloc) {
				n_alloc <<= 1;
				z = (double *) GMT_memory ((void *)z, (size_t)n_alloc, sizeof(double), GMT_program);
			}
		}
		GMT_fclose (fpl);
		if (nz == 0) {
			fprintf (stderr, "%s: Error: No intervals in file %s\n", GMT_program, Ctrl->T.file);
			exit (EXIT_FAILURE);
		}
		z = (double *) GMT_memory ((void *)z, (size_t)nz, sizeof(double), GMT_program);
	}
	else if (Ctrl->T.active && Ctrl->Q.mode == 2) {	/* Establish a log10 grid */
		if (!(Ctrl->T.inc == 1.0 || Ctrl->T.inc == 2.0 || Ctrl->T.inc == 3.0)) {
			fprintf (stderr, "%s: Error: For -Qo logarithmic spacing, dz must be 1, 2, or 3\n", GMT_program);
			exit (EXIT_FAILURE);
		}
		if (Ctrl->T.low <= 0.0) {
			fprintf (stderr, "%s: Error: For -Qo logarithmic spacing, z_start must be > 0\n", GMT_program);
			exit (EXIT_FAILURE);
		}
		nz = GMT_log_array (Ctrl->T.low, Ctrl->T.high, Ctrl->T.inc, &z);
	}
	else if (Ctrl->T.active) {	/* Establish linear grid */
		nz = irint ((Ctrl->T.high - Ctrl->T.low) / Ctrl->T.inc) + 1;
		z = (double *) GMT_memory (VNULL, (size_t)nz, sizeof(double), GMT_program);

		for (i = 0; i < nz; i++) z[i] = Ctrl->T.low + i * Ctrl->T.inc;	/* Desired z values */
	}
	else {	/* Just copy what was in the cpt file */
		nz = GMT_n_colors + 1;
		z = (double *) GMT_memory (VNULL, (size_t)nz, sizeof(double), GMT_program);
		if (Ctrl->I.active) {
			/* Reverse the intervals (only relavant for non-equidistant color maps) */
			for (i = 0; i < nz-1; i++) z[i] = GMT_lut[0].z_low + GMT_lut[GMT_n_colors-1].z_high - GMT_lut[GMT_n_colors-1-i].z_high;
		}
		else {
			for (i = 0; i < nz-1; i++) z[i] = GMT_lut[i].z_low;
		}
		z[i] = GMT_lut[i-1].z_high;
	}

	if (Ctrl->Q.mode == 2) for (i = 0; i < nz; i++) z[i] = d_log10 (z[i]);	/* Make log10(z) values for interpolation step */

	/* Write to GMT_stdout.  */

	fprintf (stdout, "#\tcpt file created by: %s", GMT_program);
	for (i = 1; i < argc; i++) fprintf (stdout, " %s", argv[i]);
	fprintf (stdout, "\n");
#ifdef WIN32
	fflush(stdout);	/* For some crazy reason if we don't do this, on Windows the comment line is written at the end of file */
#endif
	GMT_sample_cpt (z, nz, Ctrl->Z.active, Ctrl->I.active, Ctrl->Q.mode);

	GMT_free ((void *)z);

	Free_makecpt_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void *New_makecpt_Ctrl () {	/* Allocate and initialize a new control structure */
	struct MAKECPT_CTRL *C;
	
	C = (struct MAKECPT_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct MAKECPT_CTRL), "New_makecpt_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	return ((void *)C);
}

void Free_makecpt_Ctrl (struct MAKECPT_CTRL *C) {	/* Deallocate control structure */
	if (C->C.file) free ((void *)C->C.file);	
	if (C->T.file) free ((void *)C->T.file);	
	GMT_free ((void *)C);	
}
