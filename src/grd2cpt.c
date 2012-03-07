/*--------------------------------------------------------------------
 *	$Id: grd2cpt.c,v 1.65 2011/07/08 21:27:06 guru Exp $
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
 * grd2cpt reads a 2d binary gridded grdfile and creates a continuous-color-
 * palette cpt file, with a non-linear histogram-equalized mapping between
 * hue and data value.  (The linear mapping can be made with makecpt.)
 *
 * Creates a cumulative distribution function f(z) describing the data
 * in the grdfile.  f(z) is sampled at z values supplied by the user
 * [with -S option] or guessed from the sample mean and standard deviation.
 * f(z) is then found by looping over the grd array for each z and counting
 * data values <= z.  Once f(z) is found then a master cpt table is resampled
 * based on a normalized f(z).
 *
 * Author:	Walter H. F. Smith
 * Date:	12-JAN-1994
 * Revised:	PW: 12-MAY-1998, for GMT 3.1
 *		PW: 08-MAR-1998, for GMT 3.2 to allow use of master cptfiles
 *		PW: 08-JUL-2000, for GMT 3.3.5
 *		JL: 27-APR-2003, added a -R option
 *		SE: 17-SEP-2003, added a -T option
 * Version:	4
 *
 */

#include "gmt.h"

struct GRD2CPT_CTRL {
	struct C {	/* -C<cpt> */
		GMT_LONG active;
		char *file;
	} C;
	struct D {	/* -D */
		GMT_LONG active;
	} D;
	struct E {	/* -E<nlevels> */
		GMT_LONG active;
		int levels;
	} E;
	struct I {	/* -I */
		GMT_LONG active;
	} I;
	struct L {	/* -L<min_limit>/<max_limit> */
		GMT_LONG active;
		double min, max;
	} L;
	struct M {	/* -M */
		GMT_LONG active;
	} M;
	struct N {	/* -N */
		GMT_LONG active;
	} N;
	struct Q {	/* -Q[i|o */
		GMT_LONG active;
		GMT_LONG mode;
	} Q;
	struct S {	/* -S<z_start>/<z_stop>/<z_inc> */
		GMT_LONG active;
		GMT_LONG cpt;
		double low, high, inc;
		char *file;
	} S;
	struct T {	/* -T<kind> */
		GMT_LONG active;
		GMT_LONG kind; /* -1 symmetric +-zmin, +1 +-zmax, -2 = +-Minx(|zmin|,|zmax|), +2 = +-Max(|zmin|,|zmax|), 0 = min to max [Default] */
	} T;
	struct Z {	/* -Z */
		GMT_LONG active;
	} Z;
};

int main (int argc, char **argv)
{
	GMT_LONG error = FALSE, global = FALSE;

#define mgrd 1000

	char *grdfile[mgrd], CPT_lis[BUFSIZ], CPT_file[BUFSIZ], format[BUFSIZ], kind;

	GMT_LONG j, ngrd = 0, i, nxy, nxyg, nfound, ngood;

	float	*zdata = NULL;

	double *z = NULL, mean, sd, w, e, s, n;

	struct GRD_HEADER grd, new;
	struct CDF_CPT {
		double	z;	/* Data value  */
		double	f;	/* Cumulative distribution function f(z)  */
	} *cdf_cpt = NULL;
	struct GRD2CPT_CTRL *Ctrl = NULL;

	FILE	*fpc = NULL;

	void *New_grd2cpt_Ctrl (), Free_grd2cpt_Ctrl (struct GRD2CPT_CTRL *C);

	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct GRD2CPT_CTRL *)New_grd2cpt_Ctrl ();	/* Allocate and initialize a new control structure */

	w = e = s = n = 0.0;

	/* Get list of available color tables in $GMT_SHAREDIR */

	GMT_getsharepath ("conf", "gmt_cpt", ".conf", CPT_lis);
	if ((fpc = fopen (CPT_lis, "r")) == NULL) {
		fprintf (stderr, "%s: ERROR: Cannot open file %s\n", GMT_program, CPT_lis);
		exit (EXIT_FAILURE);
	}

	for (i = 1; !error && i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {

				/* Common parameters */

				case 'V':
				case 'R':
				case '\0':
					error += GMT_parse_common_options (argv[i], &w, &e, &s, &n);
					break;

				/* Supplemental parameters */

				case 'C':	/* Get cpt table */
					Ctrl->C.active = TRUE;
					Ctrl->C.file = strdup (&argv[i][2]);
					break;

				case 'D':
					Ctrl->D.active = TRUE;
					break;

				case 'E':
					Ctrl->E.active = TRUE;
					if (sscanf(&argv[i][2], "%d", &Ctrl->E.levels) != 1) {
						fprintf(stderr,"%s: GMT SYNTAX ERROR -E option:  Cannot decode value\n", GMT_program);
						error++;
					}
					break;

				case 'I':
					Ctrl->I.active = TRUE;
					break;

				case 'L':
					Ctrl->L.active = TRUE;
					if (sscanf(&argv[i][2], "%lf/%lf", &Ctrl->L.min, &Ctrl->L.max) != 2) {
						fprintf(stderr,"%s: GMT SYNTAX ERROR -L option:  Cannot decode limits\n", GMT_program);
						error++;
					}
					break;

				case 'M':
					Ctrl->M.active = TRUE;
					break;

				case 'N':
					Ctrl->N.active = TRUE;
					break;

				case 'Q':
					Ctrl->Q.active = TRUE;
					if (argv[i][2] == 'o')	/* Input data is z, but take log10(z) before interpolation colors */
						Ctrl->Q.mode = 2;
					else			/* Input is log10(z) */
						Ctrl->Q.mode = 1;
					break;

				case 'S':
					Ctrl->S.active = TRUE;
					if (sscanf(&argv[i][2], "%lf/%lf/%lf", &Ctrl->S.low, &Ctrl->S.high, &Ctrl->S.inc) != 3) {
						fprintf(stderr,"%s: GMT SYNTAX ERROR -S option:  Cannot decode values\n", GMT_program);
						error++;
					}
					break;
				case 'T':
					Ctrl->T.active = TRUE;
					kind = '\0';
					if (sscanf(&argv[i][2], "%c", &kind) != 1) {
						fprintf(stderr,"%s: GMT SYNTAX ERROR -T option:  Cannot decode option\n", GMT_program);
						error++;
					}
					switch (kind) {
						case '-':	/* Symmetric with |zmin| range */
							Ctrl->T.kind = -1;
							break;
						case '+':	/* Symmetric with |zmax| range */
							Ctrl->T.kind = +1;
							break;
						case '_':	/* Symmetric with min(|zmin|,|zmax|) range */
							Ctrl->T.kind = -2;
							break;
						case '=':	/* Symmetric with max(|zmin|,|zmax|) range */
							Ctrl->T.kind = +2;
							break;
						default:	/* plain min to max */
							fprintf(stderr,"%s: GMT SYNTAX ERROR -T option:  Must append modifier -, +, _, or =\n", GMT_program);
							error++;
							break;
					}
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
		else {
			if (ngrd > mgrd) {
				fprintf (stderr, "%s: GMT USAGE ERROR:  Cannot handle more than %i grids\n", GMT_program, mgrd);
				exit(EXIT_FAILURE);
			}

			grdfile[ngrd] = argv[i];

			if (ngrd == 0) {
				if (GMT_err_pass (GMT_read_grd_info (grdfile[0], &grd), grdfile[0])) error++;
			}
			else {
				if (GMT_err_pass (GMT_read_grd_info (grdfile[ngrd], &new), grdfile[ngrd])) error++;
				if (new.x_min != grd.x_min || new.x_max != grd.x_max || new.y_min != grd.y_min || new.y_max != grd.y_max || new.nx != grd.nx || new.ny != grd.ny) {
					fprintf (stderr, "%s: GMT USAGE ERROR:  Grid %s does not have the same dimensions and size as grid %s\n", GMT_program, grdfile[ngrd], grdfile[0]);
					exit (EXIT_FAILURE);
				}
			}
			ngrd++;
		}
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "grd2cpt %s - Make a linear or histogram-equalized color palette table from a grdfile\n\n", GMT_VERSION);
		fprintf (stderr, "usage: grd2cpt <grdfiles> [-C<table>] [-D] [-E<nlevels> [-I] [-L<min_limit>/<max_limit>]\n");
		fprintf (stderr, "\t[-M] [-N] [-Q[i|o]] [%s] [-S<z_start>/<z_stop>/<z_inc>] [-T<-|+|=|_>] [-V] [-Z]\n", GMT_Rgeo_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\t<grdfiles> names of one or more 2-D binary data sets\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-C  Specify a colortable [Default is rainbow]:\n");
		fprintf (stderr, "\t   [Original z-range is given in brackets]\n");
		fprintf (stderr, "\t   ---------------------------------\n");
		while (fgets (format, BUFSIZ, fpc)) if (!(format[0] == '#' || format[0] == 0)) fprintf (stderr, "\t   %s", format);
		fclose (fpc);
		fprintf (stderr, "\t   ---------------------------------\n");
		fprintf (stderr, "\t-D Set back- and foreground color to match the bottom/top limits in the cpt file [Default uses color table].\n");
		fprintf (stderr, "\t-E nlevels equidistant color levels from zmin to zmax.\n");
		fprintf (stderr, "\t-I Reverses the sense of the color table as well as back- and foreground color.\n");
		fprintf (stderr, "\t-M Use GMT defaults to set back-, foreground, and NaN colors [Default uses color table].\n");
		fprintf (stderr, "\t-N Do not write back-, foreground, and NaN colors [Default will].\n");
		fprintf (stderr, "\t-L Limit the range of the data [Default uses actual min,max of data].\n");
		fprintf (stderr, "\t-Q assign a logarithmic colortable [Default is linear]\n");
		fprintf (stderr, "\t   -Qi: z-values are actually log10(z). Assign colors and write z. [Default]\n");
		fprintf (stderr, "\t   -Qo: z-values are z, but take log10(z), assign colors and write z.\n");
		GMT_explain_option ('R');
		fprintf (stderr, "\t-S Sample points should Step from z_start to z_stop by z_inc [Default guesses some values].\n");
		fprintf (stderr, "\t-T Force color tables to be symmetric about 0. Append one modifier:\n");
		fprintf (stderr, "\t   - for values symmetric about zero from -|zmin| to +|zmin|\n");
		fprintf (stderr, "\t   + for values symmetric about zero from -|zmax| to +|zmax|\n");
		fprintf (stderr, "\t   _ for values symmetric about zero -+min(|zmin|,|zmax|)\n");
		fprintf (stderr, "\t   = for values symmetric about zero -+max(|zmin|,|zmax|)\n");
		GMT_explain_option ('V');
		fprintf (stderr, "\t-Z Create a continuous color palette [Default is discontinuous, i.e., constant color intervals]\n");
		exit (EXIT_FAILURE);
	}

	fclose (fpc);

	if (!Ctrl->C.active) {	/* Must assign default table [rainbow] */
		Ctrl->C.active = TRUE;
		Ctrl->C.file = strdup ("rainbow");
	}
	if (!Ctrl->E.active) Ctrl->E.levels =  11;	/* Default number of levels */

	/* Open the specified master color table */

	error += GMT_set_cpt_path (CPT_file, Ctrl->C.file);

	if (ngrd < 1) {
		fprintf(stderr,"%s: GMT USAGE ERROR:  No grid name(s) specified.\n", GMT_program);
		error++;
	}
	if (Ctrl->L.active && Ctrl->L.min >= Ctrl->L.max) {
		fprintf(stderr,"%s: GMT SYNTAX ERROR -L option:  min_limit must be less than max_limit.\n", GMT_program);
		error++;
	}
	if (Ctrl->S.active && (Ctrl->S.high <= Ctrl->S.low || Ctrl->S.inc <= 0.0)) {
		fprintf (stderr,"%s: GMT SYNTAX ERROR -S option:  Bad arguments\n", GMT_program);
		error++;
	}
	if (Ctrl->S.active && (Ctrl->T.active || Ctrl->E.active)) {
		fprintf (stderr,"%s: GMT USAGE ERROR -S option:  Cannot be combined with -E nor -T option.\n", GMT_program);
		error++;
	}

	if (error) exit (EXIT_FAILURE);

	if (Ctrl->N.active) GMT_cpt_flags += 1;	/* bit 0 controls if BFN will be written out */
	if (Ctrl->D.active) GMT_cpt_flags += 2;	/* bit 1 controls if BF will be set to equal bottom/top rgb value */
	if (Ctrl->M.active) GMT_cpt_flags += 4;	/* bit 2 controls if BFN is determined by parameters */

	GMT_read_cpt (CPT_file);

	if (e > w && n > s) {
		global = (fabs (grd.x_max - grd.x_min) == 360.0);
		if (!global && (w < grd.x_min || e > grd.x_max)) error = TRUE;
		if (s < grd.y_min || n > grd.y_max) error = TRUE;
		if (error) {
			fprintf (stderr, "%s: GMT ERROR: Subset exceeds data domain!\n", GMT_program);
			exit (EXIT_FAILURE);
		}
		nxy = GMT_get_nm (GMT_get_n (w, e, grd.x_inc, grd.node_offset), GMT_get_n (s, n, grd.y_inc, grd.node_offset));
	}
	else
		nxy = GMT_get_nm (grd.nx, grd.ny);

	nxyg = GMT_get_nm (nxy, ngrd);

	zdata = (float *) GMT_memory (VNULL, (size_t) nxyg, sizeof (float), GMT_program);

	for (i = 0; i < ngrd; i++) {
		GMT_err_fail (GMT_read_grd_info (grdfile[i], &grd), grdfile[i]);
		GMT_err_fail (GMT_read_grd (grdfile[i], &grd, &(zdata[i*nxy]), w, e, s, n, GMT_pad, FALSE), grdfile[i]);
	}

	/* Loop over the file and find NaNs.  If set limits, may create more NaNs  */
	nfound = 0;
	mean = sd = 0.0;
	if (Ctrl->L.active) {
		/* Loop over the grdfile, and set anything outside the limiting values to NaN.  */

		grd.z_min = Ctrl->L.min;
		grd.z_max = Ctrl->L.max;
		for (i = 0; i < nxyg; i++) {
			if (GMT_is_fnan (zdata[i]))
				nfound++;
			else {
				if (zdata[i] < Ctrl->L.min || zdata[i] > Ctrl->L.max) {
					nfound++;
					zdata[i] = GMT_f_NaN;
				}
				else {
					mean += zdata[i];
					sd += zdata[i] * zdata[i];
				}
			}
		}
	}
	else {
		Ctrl->L.min = grd.z_max;	/* This is just to double check grd.z_min, grd.z_max  */
		Ctrl->L.max = grd.z_min;
		for (i = 0; i < nxyg; i++) {
			if (GMT_is_fnan (zdata[i]))
				nfound++;
			else {
				if (zdata[i] < Ctrl->L.min) Ctrl->L.min = zdata[i];
				if (zdata[i] > Ctrl->L.max) Ctrl->L.max = zdata[i];
				mean += zdata[i];
				sd += zdata[i] * zdata[i];
			}
		}
		grd.z_min = Ctrl->L.min;
		grd.z_max = Ctrl->L.max;
	}
	ngood = nxyg - nfound;	/* This is the number of non-NaN points for the cdf function  */
	mean /= ngood;
	sd /= ngood;
	sd = sqrt(sd - mean*mean);
	if (gmtdefs.verbose) {
		sprintf(format,"%%s:  Mean and S.D. of data are %s %s\n", gmtdefs.d_format, gmtdefs.d_format);
		fprintf(stderr, format, GMT_program, mean, sd);
	}

	/* Now the zdata are ready.  Decide how to make steps in z.  */
	if (Ctrl->S.active) {
		/* Use predefined levels and interval */
		Ctrl->E.levels = (grd.z_min < Ctrl->S.low) ? 1 : 0;
		Ctrl->E.levels += (int)floor((Ctrl->S.high - Ctrl->S.low)/Ctrl->S.inc) + 1;
		if (grd.z_max > Ctrl->S.high) Ctrl->E.levels++;
		cdf_cpt = (struct CDF_CPT *)GMT_memory (VNULL, (size_t)Ctrl->E.levels, sizeof(struct CDF_CPT), GMT_program);
		if (grd.z_min < Ctrl->S.low) {
			cdf_cpt[0].z = grd.z_min;
			cdf_cpt[1].z = Ctrl->S.low;
			i = 2;
		}
		else {
			cdf_cpt[0].z = Ctrl->S.low;
			i = 1;
		}
		j = (grd.z_max > Ctrl->S.high) ? Ctrl->E.levels - 1 : Ctrl->E.levels;
		while (i < j) {
			cdf_cpt[i].z = cdf_cpt[i-1].z + Ctrl->S.inc;
			i++;
		}
		if (j == Ctrl->E.levels-1) cdf_cpt[j].z = grd.z_max;
	}

	else if (Ctrl->T.active || Ctrl->E.active) {
		/* Make a equaldistant color map from grd.z_min to grd.z_max */
		double start, range;

		switch (Ctrl->T.kind) {
			case -1:
				start = -fabs((double)grd.z_min);
				break;
			case 1:
				start = -fabs((double)grd.z_max);
				break;
			case -2:
				start = -MIN(fabs((double)grd.z_min), fabs((double)grd.z_max));
				break;
			case 2:
				start = -MAX(fabs((double)grd.z_min), fabs((double)grd.z_max));
				break;
			default:
				start = grd.z_min;
				break;
		}
		range = (Ctrl->T.kind) ? 2.0 * fabs (start) : grd.z_max - grd.z_min;
		Ctrl->S.inc = range / (double)(Ctrl->E.levels - 1);
		cdf_cpt = (struct CDF_CPT *)GMT_memory (VNULL, (size_t)Ctrl->E.levels, sizeof(struct CDF_CPT), GMT_program);
		for (i = 0; i < Ctrl->E.levels; i++) {
			cdf_cpt[i].z = start + i * Ctrl->S.inc;
		}
	}

	else {
		/* This is completely ad-hoc.  It chooses z based on steps of 0.1 for a Gaussian CDF:  */
		cdf_cpt = (struct CDF_CPT *)GMT_memory (VNULL, (size_t)Ctrl->E.levels, sizeof(struct CDF_CPT), GMT_program);
		if ((mean - 1.28155*sd) <= grd.z_min || (mean + 1.28155*sd) >= grd.z_max) {
			mean = 0.5 * (grd.z_min + grd.z_max);
			sd = (grd.z_max - mean) / 1.5;
			if (sd <= 0.0) {
				fprintf (stderr, "%s:  ERROR.  Min and Max data values are equal.\n", GMT_program);
				exit (EXIT_FAILURE);
			}
		}	/* End of stupid bug fix  */

		cdf_cpt[0].z = grd.z_min;
		cdf_cpt[1].z = mean - 1.28155 * sd;
		cdf_cpt[2].z = mean - 0.84162 * sd;
		cdf_cpt[3].z = mean - 0.52440 * sd;
		cdf_cpt[4].z = mean - 0.25335 * sd;
		cdf_cpt[5].z = mean;
		cdf_cpt[6].z = mean + 0.25335 * sd;
		cdf_cpt[7].z = mean + 0.52440 * sd;
		cdf_cpt[8].z = mean + 0.84162 * sd;
		cdf_cpt[9].z = mean + 1.28155 * sd;
		cdf_cpt[10].z = grd.z_max;
	}

	/* Get here when we are ready to go.  cdf_cpt[].z contains the sample points.  */

	if (gmtdefs.verbose) sprintf (format, "%%s: z = %s and CDF(z) = %s\n", gmtdefs.d_format, gmtdefs.d_format);
	for (j = 0; j < Ctrl->E.levels; j++) {
		if (cdf_cpt[j].z <= grd.z_min)
			cdf_cpt[j].f = 0.0;
		else if (cdf_cpt[j].z >= grd.z_max)
			cdf_cpt[j].f = 1.0;
		else {
			nfound = 0;
			for (i = 0; i < nxyg; i++) {
				if (!GMT_is_fnan (zdata[i]) && zdata[i] <= cdf_cpt[j].z) nfound++;
			}
			cdf_cpt[j].f = (double)(nfound-1)/(double)(ngood-1);
		}
		if (gmtdefs.verbose) fprintf (stderr, format, GMT_program, cdf_cpt[j].z, cdf_cpt[j].f);
	}

	/* Now the cdf function has been found.  We now resample the chosen cptfile  */

	/* Write to GMT_stdout.  */

	sprintf (format, "#\tcpt file created by: %s", GMT_program);	GMT_fputs (format, GMT_stdout);
	for (i = 1; i < argc; i++) {sprintf (format, " %s", argv[i]);	GMT_fputs (format, GMT_stdout);}
	GMT_fputs ("\n", GMT_stdout);


	z = (double *) GMT_memory (VNULL, (size_t)Ctrl->E.levels, sizeof (double), GMT_program);
	for (i = 0; i < Ctrl->E.levels; i++) z[i] = cdf_cpt[i].z;
	if (Ctrl->Q.mode == 2) for (i = 0; i < Ctrl->E.levels; i++) z[i] = d_log10 (z[i]);	/* Make log10(z) values for interpolation step */

	GMT_sample_cpt (z, -Ctrl->E.levels, Ctrl->Z.active, Ctrl->I.active, Ctrl->Q.mode);	/* -ve to keep original colors */

	GMT_free ((void *)cdf_cpt);
	GMT_free ((void *)zdata);
	GMT_free ((void *)z);

	Free_grd2cpt_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void *New_grd2cpt_Ctrl () {	/* Allocate and initialize a new control structure */
	struct GRD2CPT_CTRL *C;

	C = (struct GRD2CPT_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GRD2CPT_CTRL), "New_grd2cpt_Ctrl");

	/* Initialize values whose defaults are not 0/FALSE/NULL */
	return ((void *)C);
}

void Free_grd2cpt_Ctrl (struct GRD2CPT_CTRL *C) {	/* Deallocate control structure */
	if (C->C.file) free ((void *)C->C.file);
	GMT_free ((void *)C);
}
