/*--------------------------------------------------------------------
 *	$Id: sample1d.c,v 1.66 2011/07/08 21:27:06 guru Exp $
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
 * sample1d reads a 1-D dataset, and resamples the values on (1) user
 * supplied time values or (2) equidistant time-values based on <timestart> and
 * <dt>, both supplied at the command line. Choose among linear, cubic
 * spline, and Akima's spline.  sample1d will handle multiple column files,
 * user must choose which column contains the independent, monotonically
 * increasing variable.
 *
 * Author:	Paul Wessel
 * Date:	05-JUL-2000
 * Version:	4
 *
 */
 
#include "gmt.h"

#define INT_1D	0	/* Regular 1-D interpolation */
#define INT_2D	1	/* Spatial 2-D path interpolation */

struct SAMPLE1D_CTRL {
#ifdef DEBUG
	struct A {	/* -A[m|p] */
		int mode;
	} A;
#endif
	struct F {	/* -Fl|a|c */
		GMT_LONG active;
		GMT_LONG mode;
	} F;
	struct I {	/* -I<inc>[e|E|k|K|m|M|n|N|c|C|d|D] */
		GMT_LONG active;
		GMT_LONG mode;
		double inc;
		char unit;
	} I;
	struct T {	/* -T<time_col> */
		GMT_LONG active;
		GMT_LONG col;
	} T;
	struct N {	/* -N<knotfile> */
		GMT_LONG active;
		char *file;
	} N;
	struct S {	/* -S<xstart> */
		GMT_LONG active;
		double start;
	} S;
};

int main (int argc, char **argv)
{
	GMT_LONG result, n_col = 0, n_read, n_fields, n_expected_fields, n_files = 0, fno, n_args, n_req;
	GMT_LONG i, j, k, n, m, n_alloc = 0, m_alloc, rows = 1, m_supplied = 0;

	GMT_LONG error = FALSE, *nan_flag = VNULL, nofile = TRUE, done = FALSE;

	double *t_supplied_out = VNULL, *t_out, *ttime, *data, **col = VNULL, **out = VNULL;
	double tt, low_t, high_t, *in = NULL, *dout = VNULL;
#ifdef DEBUG
	GMT_LONG proj_type = 0;
	double d_scale;
	PFD distance_func;
#endif
	char line[BUFSIZ], type[3], *not_used = NULL;

	FILE *fp = NULL, *fpt = NULL;


	struct SAMPLE1D_CTRL *Ctrl = NULL;

	void *New_sample1d_Ctrl (), Free_sample1d_Ctrl (struct SAMPLE1D_CTRL *C);
	
	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct SAMPLE1D_CTRL *)New_sample1d_Ctrl ();	/* Allocate and initialize a new control structure */
	
	type[0] = 'l';	type[1] = 'a';	type[2] = 'c';

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {

				/* Common parameters */

				case 'H':
				case 'M':
				case 'V':
				case 'b':
				case 'f':
				case 'm':
				case '\0':
					error += (GMT_LONG)GMT_parse_common_options (argv[i], 0, 0, 0, 0);
					break;

				/* Supplemental parameters */

				case 'F':
					Ctrl->F.active = TRUE;
					switch (argv[i][2]) {
						case 'l':
							Ctrl->F.mode = 0;
							break;
						case 'a':
							Ctrl->F.mode = 1;
							break;
						case 'c':
							Ctrl->F.mode = 2;
							break;
						case 'n':
							Ctrl->F.mode = 3;
							break;
						default:	/* Use GMT defaults */
							break;
					}
					break;
				case 'I':
					Ctrl->I.active = TRUE;
					Ctrl->I.inc = atof (&argv[i][2]);
#ifdef DEBUG
					n = strlen (argv[i]) - 1;
					if (strchr ("ekmnd", argv[i][n])) {
						Ctrl->I.unit = argv[i][n];
						Ctrl->I.mode = INT_2D;
					}
#endif
					break;
				case 'T':
					Ctrl->T.active = TRUE;
					Ctrl->T.col = atoi (&argv[i][2]);
					break;
				case 'N':
					Ctrl->N.file = strdup (&argv[i][2]);
					Ctrl->N.active = TRUE;
					break;
				case 'S':
					Ctrl->S.active = TRUE;
					GMT_scanf_arg (&argv[i][2], GMT_IS_UNKNOWN, &Ctrl->S.start);
					break;

				/* For backward compatibility for now */

				case 'L':
					Ctrl->F.active = TRUE;
					Ctrl->F.mode = 0;
					break;
				case 'A':
#ifdef DEBUG
					if (argv[i][2] == 'm')
						Ctrl->A.mode = 1;
					else if (argv[i][2] == 'p')
						Ctrl->A.mode = 2;
					else {	/* Old backwards compatible -A option */
#endif
						Ctrl->F.active = TRUE;
						Ctrl->F.mode = 1;
#ifdef DEBUG
					}
#endif
					break;
				case 'C':
					Ctrl->F.active = TRUE;
					Ctrl->F.mode = 2;
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

	if (argc == 1 || GMT_give_synopsis_and_exit) {	/* Display usage */
		fprintf (stderr,"sample1d %s - Resampling of 1-D data sets\n\n", GMT_VERSION);
		fprintf (stderr, "usage: sample1d <infile(s)> [-Fl|a|c|n] [%s] [-I<t_inc>] [-N<knotfile>]\n", GMT_H_OPT);
		fprintf (stderr, "\t[-S<xstart>] [-T<time_col>] [-V] [%s] [%s] [%s]\n\n", GMT_b_OPT, GMT_f_OPT, GMT_m_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\t<infile> is one or more multicolumn ASCII (or binary, see -b) tables. [Default is standard input]\n");
		fprintf (stderr, "\tThe independent variable (see -T) must be monotonically in/de-creasing\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-F sets the interpolation mode.  Choose from:\n");
		fprintf (stderr, "\t   l Linear interpolation\n");
		fprintf (stderr, "\t   a Akima spline interpolation\n");
		fprintf (stderr, "\t   c Cubic spline interpolation\n");
		fprintf (stderr, "\t   n No interpolation (nearest point)\n");
		fprintf (stderr, "\t   [Default is -F%c]\n", type[Ctrl->F.mode]);
		GMT_explain_option ('H');
		fprintf (stderr, "\t-I <x_inc> sets equidistant grid interval [x1 - x0]\n");
#ifdef DEBUG
		fprintf (stderr, "\t   Append e|k|m|n|d to indicate that the first two columns contain\n");
		fprintf (stderr, "\t   longitude, latitude and you wish to resample this path using great\n");
		fprintf (stderr, "\t   circle segment with a nominal spacing of <t_inc> in those units.\n");
		fprintf (stderr, "\t   See -Am|p to only sample along meridians and parallels.\n");
#endif
		fprintf (stderr, "\t-N <knotfile> is an ASCII table with the desired time positions in column 0\n");
		fprintf (stderr, "\t   Overrides the -I and -S settings.  If none of -I, -S, and -N is set\n");
		fprintf (stderr, "\t   then <tstart> = first input point, <t_inc> = (t[1] - t[0])\n");
		fprintf (stderr, "\t-S <xstart> sets the first output point [first multiple of x_inc in range]\n");
		fprintf (stderr, "\t-T gives column number of the independent variable (time) [Default is 0 (first)]\n");
		GMT_explain_option ('V');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf (stderr, "\t   Default is 2 input columns\n");
		GMT_explain_option ('o');
		GMT_explain_option ('n');
		GMT_explain_option ('f');
		GMT_explain_option ('m');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	if (Ctrl->T.col < 0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -T option: Column number cannot be negative\n", GMT_program);
		error++;
	}
	if (Ctrl->N.active && Ctrl->I.active) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR: Specify only one of -N and -S\n", GMT_program);
		error++;
	}
	if (Ctrl->I.active && Ctrl->I.inc <= 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -I option: Must specify positive increment\n", GMT_program);
		error++;
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
        if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = 2;
	n_req = (Ctrl->T.col >= 2) ? Ctrl->T.col + 1 : 2;
        if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] < n_req) {
                fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data (-bi) must have at least %ld columns\n", GMT_program, n_req);
		error++;
	}
	if (Ctrl->N.active && (fpt = GMT_fopen (Ctrl->N.file, "r")) == NULL) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -N. Cannot open file %s\n", GMT_program, Ctrl->N.file);
		error++;
	}

#ifdef DEBUG
	if (Ctrl->I.mode) error += GMT_get_dist_scale (Ctrl->I.unit, &d_scale, &proj_type, &distance_func);
	if (Ctrl->I.mode) Ctrl->I.inc /= d_scale;			/* Convert increment to spherical degrees */
#endif
	if (error) exit (EXIT_FAILURE);

	if (GMT_io.binary[GMT_IN] && gmtdefs.verbose) {
		char *type[2] = {"double", "single"};
		fprintf (stderr, "%s: Expects %ld-column %s-precision binary data\n", GMT_program, GMT_io.ncol[GMT_IN], type[GMT_io.single_precision[GMT_IN]]);
	}

#ifdef SET_IO_MODE
	GMT_setmode (GMT_OUT);
#endif

	gmtdefs.interpolant = Ctrl->F.mode;
	GMT_io.skip_if_NaN[GMT_X] = GMT_io.skip_if_NaN[GMT_Y] = FALSE;	/* Turn off default GMT NaN-handling for (x,y) which is not the case here */
	GMT_io.skip_if_NaN[Ctrl->T.col] = TRUE;				/* ... But disallow NaN in "time" column */
	
	t_out = (double *)NULL;

	if (Ctrl->N.active) {	/* read file with abscissae */
		m_alloc = GMT_CHUNK;
		t_supplied_out = (double *) GMT_memory (VNULL, (size_t)m_alloc, sizeof (double), GMT_program);
		m = 0;
		if (GMT_io.io_header[GMT_IN]) for (i = 0; i < GMT_io.n_header_recs; i++) not_used = GMT_fgets (line, BUFSIZ, fpt);
		n_expected_fields = (GMT_io.ncol[GMT_IN]) ? GMT_io.ncol[GMT_IN] : GMT_MAX_COLUMNS;
		n_fields = GMT_input (fpt, &n_expected_fields, &in);
		while (! (GMT_io.status & GMT_IO_EOF)) {	/* Not yet EOF */
			while (GMT_io.status & GMT_IO_SEGMENT_HEADER) {
				n_fields = GMT_input (fpt, &n_expected_fields, &in);
			}
			if (GMT_io.status & GMT_IO_EOF) continue;	/* At EOF */
			t_supplied_out[m] = in[GMT_X];
			m++;
			if (m == m_alloc) {	/* Get more memory */
				m_alloc <<= 1;
				t_supplied_out = (double *) GMT_memory ((void *)t_supplied_out, (size_t)m_alloc, sizeof (double), GMT_program);
			}
			n_fields = GMT_input (fpt, &n_expected_fields, &in);
		}
		GMT_fclose (fpt);
		m_supplied = m;
		t_supplied_out = (double *) GMT_memory ((void *)t_supplied_out, (size_t)m_supplied, sizeof (double), GMT_program);
		t_out = (double *) GMT_memory (VNULL, (size_t)m_supplied, sizeof (double), GMT_program);
                if (gmtdefs.verbose) fprintf (stderr, "%s: Read %ld knots from file\n", GMT_program, m_supplied);
	}

	if (n_files > 0)
		nofile = FALSE;
	else
		n_files = 1;

	n_args = (argc > 1) ? argc : 2;

	for (fno = 1; !done && fno < n_args; fno++) {	/* Loop over input files, if any */
		if (!nofile && argv[fno][0] == '-') continue;

		if (nofile) {	/* Just read standard input */
			fp = GMT_stdin;
			done = TRUE;
#ifdef SET_IO_MODE
			GMT_setmode (GMT_IN);
#endif
		}
		else if ((fp = GMT_fopen (argv[fno], GMT_io.r_mode)) == NULL) {
			fprintf (stderr, "%s: Cannot open file %s\n", GMT_program, argv[fno]);
			continue;
		}

		if (!nofile && gmtdefs.verbose) fprintf (stderr, "%s: Working on file %s\n", GMT_program, argv[fno]);

		if (GMT_io.io_header[GMT_IN]) {
			for (i = 0; i < GMT_io.n_header_recs; i++) {
				not_used = GMT_fgets (line, BUFSIZ, fp);
				if (GMT_io.io_header[GMT_OUT]) GMT_fputs (line, GMT_stdout);
			}
		}

		n_expected_fields = (GMT_io.ncol[GMT_IN]) ? GMT_io.ncol[GMT_IN] : GMT_MAX_COLUMNS;
		n_read = n_col = 0;
		n_fields = GMT_input (fp, &n_expected_fields, &in);

		while (! (GMT_io.status & GMT_IO_EOF)) {	/* Not yet EOF */

			while (GMT_io.status & GMT_IO_SEGMENT_HEADER) {
				GMT_write_segmentheader (GMT_stdout, n_expected_fields);
				n_fields = GMT_input (fp, &n_expected_fields, &in);
			}
			if (GMT_io.status & GMT_IO_EOF) continue;	/* At EOF */

			n = 0;

			while (! (GMT_io.status & (GMT_IO_SEGMENT_HEADER | GMT_IO_EOF))) {	/* Keep going until FALSE or = 2 segment header */
				n_read++;
				if (GMT_io.status & GMT_IO_MISMATCH) {
					fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld\n", GMT_program, n_fields, n_expected_fields, n_read);
					continue;
				}

				if (n_col == 0) {	/* Allocate memory first time around */
					n_col = n_expected_fields;
					if (Ctrl->T.col >= n_col) {
						fprintf (stderr, "%s: time_col = %ld exceeds range of columns (%ld)!\n", GMT_program, Ctrl->T.col, n_col);
						exit (EXIT_FAILURE);
					}
					dout = (double *) GMT_memory (VNULL, (size_t)n_col, sizeof (double), GMT_program);
					col = (double **) GMT_memory (VNULL, (size_t)n_col, sizeof (double *), GMT_program);
					out = (double **) GMT_memory (VNULL, (size_t)n_col, sizeof (double *), GMT_program);
					nan_flag = (GMT_LONG *) GMT_memory (VNULL, (size_t)n_col, sizeof (GMT_LONG), GMT_program);
					n_alloc = GMT_CHUNK;
					for (j = 0; j < n_col; j++) col[j] = (double *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (double), GMT_program);
				}

				for (j = 0; j < n_col; j++) {
					if (GMT_is_dnan (in[j])) nan_flag[j] = TRUE;
					col[j][n] = in[j];
				}
				n++;

				if (n == n_alloc) {	/* Get more memory */
					n_alloc <<= 1;
					for (j = 0; j < n_col; j++) col[j] = (double *) GMT_memory ((void *)col[j], (size_t)n_alloc, sizeof (double), GMT_program);
				}

				n_fields = GMT_input (fp, &n_expected_fields, &in);
			}

			/* Here we have one segment to work on */

			/* If we didn't get input abscissa, now's the time to generate them */

#ifdef DEBUG
			if (Ctrl->I.active && Ctrl->I.mode == INT_2D) {	/* Special spatial interpolation */
				GMT_LONG np;
				double *lon, *lat, xyout[2];
				lon = (double *) GMT_memory (VNULL, n, sizeof (double), GMT_program);
				lat = (double *) GMT_memory (VNULL, n, sizeof (double), GMT_program);
				memcpy ((void *)lon, (void *)col[GMT_X], n*sizeof(double));
				memcpy ((void *)lat, (void *)col[GMT_Y], n*sizeof(double));
				np = GMT_fix_up_path (&lon, &lat, n, Ctrl->I.inc, Ctrl->A.mode);
				for (i = 0; i < np; i++) {
					xyout[GMT_X] = lon[i];	xyout[GMT_Y] = lat[i];
					GMT_output (GMT_stdout, 2, xyout);
				}
				rows += n + 1;
				GMT_free ((void *)lon);
				GMT_free ((void *)lat);
				continue;
			}
#endif
			if (Ctrl->N.active) {	/* Get relevant t_out segment */
				low_t  = MIN (col[Ctrl->T.col][0], col[Ctrl->T.col][n-1]);
				high_t = MAX (col[Ctrl->T.col][0], col[Ctrl->T.col][n-1]);
				for (i = m = 0; i < m_supplied; i++) {
					if (t_supplied_out[i] < low_t || t_supplied_out[i] > high_t) continue;
					t_out[m++] = t_supplied_out[i];
				}
				if (m == 0) fprintf (stderr, "%s: Warning: No output points for range %g to %g\n", GMT_program, col[Ctrl->T.col][0], col[Ctrl->T.col][n-1]);
			}
			else {	/* Generate evenly spaced grid */
				if (!Ctrl->I.active) Ctrl->I.inc = col[Ctrl->T.col][1] - col[Ctrl->T.col][0];
				if (Ctrl->I.active && (col[Ctrl->T.col][1] - col[Ctrl->T.col][0]) < 0.0 && Ctrl->I.inc > 0.0) Ctrl->I.inc = -Ctrl->I.inc;	/* For monotonically decreasing data */
				if (!Ctrl->S.active) {
					if (Ctrl->I.inc > 0.0) {
						Ctrl->S.start = floor (col[Ctrl->T.col][0] / Ctrl->I.inc) * Ctrl->I.inc;
						if (Ctrl->S.start < col[Ctrl->T.col][0]) Ctrl->S.start += Ctrl->I.inc;
					}
					else {
						Ctrl->S.start = ceil (col[Ctrl->T.col][0] / (-Ctrl->I.inc)) * (-Ctrl->I.inc);
						if (Ctrl->S.start > col[Ctrl->T.col][0]) Ctrl->S.start += Ctrl->I.inc;
					}
				}
				m = m_alloc = irint (fabs((col[Ctrl->T.col][n-1] - Ctrl->S.start) / Ctrl->I.inc)) + 1;
				t_out = (double *) GMT_memory ((void *)t_out, (size_t)m_alloc, sizeof (double), GMT_program);
				t_out[0] = Ctrl->S.start;
				i = 1;
				if (Ctrl->I.inc > 0.0) {
					while (i < m && (tt = Ctrl->S.start + i * Ctrl->I.inc) <= col[Ctrl->T.col][n-1]) {
						t_out[i] = tt;
						i++;
					}
				}
				else {
					while (i < m && (tt = Ctrl->S.start + i * Ctrl->I.inc) >= col[Ctrl->T.col][n-1]) {
						t_out[i] = tt;
						i++;
					}
				}
				m = i;
				if (fabs (t_out[m-1]-col[Ctrl->T.col][n-1]) < GMT_SMALL) {	/* Fix roundoff */
					t_out[m-1] = col[Ctrl->T.col][n-1];
				}
			}

			if (nan_flag[Ctrl->T.col]) {
				fprintf (stderr, "%s: Independent column has NaN's!\n", GMT_program);
				exit (EXIT_FAILURE);
			}

			for (j = 0; m && j < n_col; j++) {

				if (j == Ctrl->T.col) continue;	/* Skip the time column */

				out[j] = (double *) GMT_memory (VNULL, (size_t)m, sizeof (double), GMT_program);

				if (nan_flag[j] && !gmtdefs.nan_is_gap) {	/* NaN's present, need "clean" time and data columns */

					ttime = (double *) GMT_memory (VNULL, (size_t)n, sizeof (double), GMT_program);
					data = (double *) GMT_memory (VNULL, (size_t)n, sizeof (double), GMT_program);
					for (i = k = 0; i < n; i++) {
						if ( GMT_is_dnan (col[j][i]) ) continue;
						ttime[k] = col[Ctrl->T.col][i];
						data[k] = col[j][i];
						k++;
					}
					result = GMT_intpol (ttime, data, k, m, t_out, out[j], Ctrl->F.mode);
					GMT_free ((void *)ttime);
					GMT_free ((void *)data);
				}
				else
					result = GMT_intpol (col[Ctrl->T.col], col[j], n, m, t_out, out[j], Ctrl->F.mode);

				if (result != 0) {
					fprintf (stderr, "%s: Error from GMT_intpol near row %ld!\n", GMT_program, rows+result+1);
					exit (EXIT_FAILURE);
				}
			}

			out[Ctrl->T.col] = t_out;

			for (i = 0; i < m; i++) {
				for (j = 0; j < n_col; j++) dout[j] = out[j][i];
				GMT_output (GMT_stdout, n_col, dout);
			}
			for (j = 0; m && j < n_col; j++) if (j != Ctrl->T.col) GMT_free ((void *)out[j]);

			rows += n + 1;
		}

		if (fp != GMT_stdin) GMT_fclose(fp);
	}

	for (j = 0; j < n_col; j++) GMT_free ((void *)col[j]);

	GMT_free ((void *)t_out);
	GMT_free ((void *)out);
	GMT_free ((void *)col);
	if (dout) GMT_free ((void *)dout);
	if (nan_flag) GMT_free ((void *)nan_flag);
	if (Ctrl->N.active) GMT_free ((void *)t_supplied_out);
	
	Free_sample1d_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void *New_sample1d_Ctrl () {	/* Allocate and initialize a new control structure */
	struct SAMPLE1D_CTRL *C;
	
	C = (struct SAMPLE1D_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct SAMPLE1D_CTRL), "New_sample1d_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	C->F.mode = gmtdefs.interpolant;
		
	return ((void *)C);
}

void Free_sample1d_Ctrl (struct SAMPLE1D_CTRL *C) {	/* Deallocate control structure */
	if (C->N.file) free ((void *)C->N.file);	
	GMT_free ((void *)C);	
}
