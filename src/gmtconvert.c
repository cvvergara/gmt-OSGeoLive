/*--------------------------------------------------------------------
 *	$Id: gmtconvert.c 9923 2012-12-18 20:45:53Z pwessel $
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
 * gmtconvert.c
 * reads a GMT table and writes it out in a different format.
  
   Author: 	P. Wessel
   Date:	13-JUL-2000.
   version:	4
   Update:	21-APR-2004 PW: Now has built-in cut/paste capabilities

*/

#define GMT_WITH_NO_PS
#include "gmt.h"

struct GMTCONVERT_CTRL {
	struct A {	/* -A */
		GMT_LONG active;
	} A;
	struct D {	/* -D[<template>] */
		GMT_LONG active;
		char *name;
	} D;
	struct E {	/* -E */
		GMT_LONG active;
		GMT_LONG mode;
	} E;
	struct F {	/* -F<cols> */
		GMT_LONG active;
		GMT_LONG *col;
	} F;
	struct L {	/* -L */
		GMT_LONG active;
	} L;
	struct I {	/* -I */
		GMT_LONG active;
	} I;
	struct N {	/* -N */
		GMT_LONG active;
	} N;
	struct S {	/* -S[~]\"search string\" */
		GMT_LONG active;
		GMT_LONG inverse;
		char *pattern;
	} S;
};

int main (int argc, char **argv)
{
	GMT_LONG pos, start, stop, fno, n_files = 0, n_out = 0, n_tot_cols = 0;
	GMT_LONG i, j, k, n_total_read = 0, n_hdr_delay = 0, n_out_seg = 0;
	GMT_LONG n_rev_alloc = 0, n_rev = 0, n_rev_recs = 0, error = 0, n_seg = 0, n_in_seg = 0, n_nan, out_seg;
	
	GMT_LONG on, ascii_in_bin_out, match = FALSE;

	double	*in = NULL, *out = VNULL, *R_out = (double *)NULL;

	char record_str[BUFSIZ], ptr[BUFSIZ], header[BUFSIZ], filename[BUFSIZ];
	
	FILE *fp_out = NULL;

	struct GMTCONVERT_CTRL *Ctrl = NULL;

	void *New_gmtconvert_Ctrl (), Free_gmtconvert_Ctrl (struct GMTCONVERT_CTRL *C);
	
	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct GMTCONVERT_CTRL *)New_gmtconvert_Ctrl ();	/* Allocate and initialize a new control structure */
	
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {


				/* Common parameters */

				case 'H':
				case 'M':
				case 'V':
				case 'b':
				case 'f':
				case 'g':
				case 'm':
				case ':':
				case '\0':
					error += GMT_parse_common_options (argv[i], 0, 0, 0, 0);
					break;

				/* Supplemental parameters */

				case 'A':
					Ctrl->A.active = TRUE;
					break;
				case 'D':               /* Write each segment to a separate output file */
					if (argv[i][2]) Ctrl->D.name = strdup (&argv[i][2]);
					Ctrl->D.active = TRUE;
					break;
				case 'E':
					Ctrl->E.active = TRUE;
					if (argv[i][2] == 'f')		/* Get first point only */
						Ctrl->E.mode = 1;
					else if (argv[i][2] == 'l')	/* Get last point only */
						Ctrl->E.mode = 2;
					else				/* Get first and point only */
						Ctrl->E.mode = 3;
					break;
				case 'F':
					pos = 0;
					Ctrl->F.active = TRUE;
					while ((GMT_strtok (&argv[i][2], ",", &pos, ptr))) {	/* Process all tokens */
						if (strchr (ptr, '-'))	/* Range of columns given. e.g., 7-9 */
							sscanf (ptr, "%" GMT_LL "d-%" GMT_LL "d", &start, &stop);
						else if (isdigit ((int)ptr[0]))	/* Just a single column, e.g., 13 */
							start = stop = atoi (ptr);
						else {				/* Badness! */
							fprintf (stderr, "%s: GMT SYNTAX ERROR.  -F: Give cols or col-ranges separated by commas!\n", GMT_program);
							error++;
						}

						for (j = start; j <= stop; j++) Ctrl->F.col[n_out++] = j;
					}
					break;
				case 'L':
					Ctrl->L.active = TRUE;
					break;
				case 'I':
					Ctrl->I.active = TRUE;
					break;
				case 'N':
					Ctrl->N.active = TRUE;
					break;
				case 'S':
					Ctrl->S.active = TRUE;
					k = (argv[i][2] == '\\' && strlen (argv[i]) > 3 && argv[i][3] == '~') ? 3 : 2;	/* Special escape if pattern starts with ~ */
					if (argv[i][2] == '~') Ctrl->S.inverse = TRUE;
					Ctrl->S.pattern = strdup (&argv[i][k+Ctrl->S.inverse]);
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
		fprintf (stderr, "gmtconvert %s - Convert format, Paste, and/or Extract columns from table data\n\n", GMT_VERSION);
		fprintf(stderr,"usage:	gmtconvert [files] [-A] [-D[<template>]] [-E[f|l]] [-F<cols>] [%s] [-I]\n", GMT_H_OPT);
		fprintf(stderr,"\t[-L] [-N] [-S[~]\"search string\"] [-V] [%s]\n\t[%s] [%s] [%s] [%s]\n\n", GMT_t_OPT, GMT_b_OPT, GMT_f_OPT, GMT_g_OPT, GMT_m_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-A Input files should be pAsted, not concatenated [Default is concatenated].\n");
		fprintf (stderr, "\t-D Writes individual segments to separate files [Default writes one multisegment file to stdout].\n");
		fprintf (stderr, "\t   Append file name template which MUST contain a C-format specified for an integer (e.g., %%d)\n");
		fprintf (stderr, "\t   [Default uses gmtconvert_segment_%%d.d].\n");
		fprintf (stderr, "\t-E Extract first and last point per segment only [Default outputs all points].\n");
		fprintf (stderr, "\t   Append f for first only or l for last only\n");
		fprintf (stderr, "\t-F Give comma-separated list of desired columns or ranges (0 is first column) [Default is all].\n");
		GMT_explain_option ('H');
		fprintf (stderr, "\t-I Invert order of rows, i.e., outputting in reverse order.\n");
		fprintf (stderr, "\t-L Output all multisegment headers only, no data records (requires -m and ASCII data).\n");
		fprintf (stderr, "\t-N Skip output records where all fields == NaN [Default writes all records].\n");
		fprintf (stderr, "\t-S Only output segments whose segment headers contains the pattern \"string\".\n");
		fprintf (stderr, "\t   Use -S~\"string\" to output segment that DO NOT contain this pattern.\n");
		fprintf (stderr, "\t   If your pattern begins with ~, escape it with \\~.\n");
		fprintf (stderr, "\t   (-S requires -m and ASCII data) [All segments].\n");
		GMT_explain_option ('V');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		GMT_explain_option ('o');
		GMT_explain_option ('f');
		GMT_explain_option ('g');
		GMT_explain_option ('m');
		exit (EXIT_FAILURE);
	}

	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Must specify number of columns in binary input data (-bi)\n", GMT_program);
		error++;
	}
	if (GMT_io.binary[GMT_IN] && (Ctrl->L.active || Ctrl->S.active)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  -L or -S requires ASCII input data\n", GMT_program);
		error++;
	}
	if (!GMT_io.multi_segments[GMT_IN] && (Ctrl->L.active || Ctrl->S.active)) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  -L or -S requires multisegment files\n", GMT_program);
		error++;
	}
	if (Ctrl->D.active && Ctrl->D.name && !strstr (Ctrl->D.name, "%")) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  output template must contain %%d\n", GMT_program);
		error++;
	}
	if (error) exit (EXIT_FAILURE);

	if (GMT_io.binary[GMT_IN] && gmtdefs.verbose) {
		char *type[2] = {"double", "single"};
		fprintf (stderr, "%s: Expects %ld-column %s-precision binary data\n", GMT_program, GMT_io.ncol[GMT_IN], type[GMT_io.single_precision[GMT_IN]]);
	}

#ifdef SET_IO_MODE
	GMT_setmode (GMT_OUT);
#endif
	on = (Ctrl->L.active) ? FALSE : TRUE;
	ascii_in_bin_out = (!GMT_io.binary[GMT_IN] && GMT_io.binary[GMT_OUT]);
	fp_out = (Ctrl->D.active) ? NULL : GMT_stdout;
	if (Ctrl->D.active && !Ctrl->D.name) Ctrl->D.name = strdup ("gmtconvert_segment_%d.d");
	out_seg = -1;

	/* When NAN_RECORDS = pass and no multi-segment output, turn off default NaN handling */
	if (gmtdefs.nan_is_gap && !GMT_io.multi_segments[GMT_OUT]) {
		for (i = 0; i < 2; i++) GMT_io.skip_if_NaN[i] = FALSE;
	}

	if (Ctrl->A.active) {	/* Do paste and probably cut as well */
		GMT_LONG *n_expected_fields = NULL, *n_fields = NULL, all_set = 0;
		GMT_LONG n_alloc = GMT_SMALL_CHUNK;
		double *val = NULL;
		FILE **fp = NULL;

		fp = (FILE **) GMT_memory (VNULL, (size_t)n_files, sizeof (FILE *), GMT_program);
		n_expected_fields = (GMT_LONG *) GMT_memory (VNULL, n_files, sizeof (GMT_LONG), GMT_program);
		n_fields = (GMT_LONG *) GMT_memory (VNULL, n_files, sizeof (GMT_LONG), GMT_program);

		for (fno = 1, k = 0; fno < argc; fno++) {	/* Loop over input files, if any */
			if (argv[fno][0] == '-') continue;
			if ((fp[k] = GMT_fopen (argv[fno], GMT_io.r_mode)) == NULL) {
				fprintf (stderr, "%s: Cannot open file %s\n", GMT_program, argv[fno]);
				continue;
			}
			k++;
		}

		if (GMT_io.binary[GMT_IN]) {	/* If binary input then we know what to store */
			n_tot_cols = GMT_io.ncol[GMT_IN] * n_files;
			val = (double *) GMT_memory (VNULL, sizeof (double), (size_t)n_tot_cols, GMT_program);
			all_set = 2;
			if (n_out == 0) {	/* No cut action, just paste */
				n_out = n_tot_cols;
				for (k = 0; k < n_tot_cols; k++) Ctrl->F.col[k] = k;
			}
			out = (double *) GMT_memory (VNULL, sizeof (double), (size_t)n_out, GMT_program);
		}
		else {			/* Must find out by examining the input ascii data */
			val = (double *) GMT_memory (VNULL, sizeof (double), (size_t)n_alloc, GMT_program);
			all_set = 0;
		}

		if (!GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {	/* Only ascii data can have header records */
			for (i = 0; i < GMT_io.n_header_recs; i++) {
				for (k = 0; k < n_files; k++) {
					GMT_fgets (record_str, BUFSIZ, fp[k]);
					if (k < (n_files-1)) record_str[strlen(record_str)-1] = '\0';
					if (k == 0)
						strcpy (header, record_str);
					else
						strcat (header, &record_str[1]);
				}
				n_total_read ++;
				if (!GMT_io.binary[GMT_OUT] && GMT_io.io_header[GMT_OUT]) GMT_fputs (header, GMT_stdout);
			}
		}

		for (k = j = 0; k < n_files; k++) {	/* Get first record from each file */
			n_expected_fields[k] = (GMT_io.ncol[GMT_IN]) ? GMT_io.ncol[GMT_IN] : GMT_MAX_COLUMNS; /* GMT_MAX_COLUMNS -> Ascii table, must first find the number of fields */
			n_fields[k] = GMT_input (fp[k], &n_expected_fields[k], &in);
			if (! GMT_io.status & GMT_IO_SEGMENT_HEADER) {	/* Good, a data record, update number of cols etc */
				all_set = 1;
				n_tot_cols += n_expected_fields[k];
				if (n_tot_cols > n_alloc) {
					n_alloc <<= 1;
					val = (double *) GMT_memory ((char *)val, (size_t)n_alloc, sizeof (double), GMT_program);
				}
				for (i = 0; i < n_fields[k]; i++) val[j++] = in[i];	/* Copy the input values */
			}
			else {	/* No, just a multisegment header */
				if (k < (n_files-1)) GMT_io.segment_header[strlen(GMT_io.segment_header)-1] = '\0';
				if (k == 0)
					strcpy (header, GMT_io.segment_header);
				else
					strcat (header, &GMT_io.segment_header[1]);
			}
		}
		n_total_read ++;
		if (all_set == 1) {	/* Means we know how many total columns and can now allocate memory */
			all_set = 2;
			if (n_out == 0) {	/* No cut action, just paste and return all columns */
				n_out = n_tot_cols;
				for (k = 0; k < n_tot_cols; k++) Ctrl->F.col[k] = k;
			}
			out = (double *) GMT_memory (VNULL, sizeof (double), (size_t)n_out, GMT_program);
			if (n_tot_cols < n_out) {
				fprintf (stderr, "%s: Columns requested (%ld) exceed actual (%ld) - aborting\n", GMT_program, n_out, n_tot_cols);
				exit (EXIT_FAILURE);
			}
		}
		while (!GMT_REC_IS_EOF) {	/* Not yet EOF */

			while (GMT_REC_IS_NEW_SEGMENT && !GMT_REC_IS_EOF) {	/* As long as we have read a segment header */
				strcpy (GMT_io.segment_header, header);
				if (ascii_in_bin_out && n_out == 0)
					n_hdr_delay++;
				else if (Ctrl->D.active) {
					out_seg++;
					sprintf (filename, Ctrl->D.name, (int)out_seg);
					if (fp_out != NULL) GMT_fclose (fp_out);
					if ((fp_out = GMT_fopen (filename, GMT_io.w_mode)) == NULL ) {
						fprintf (stderr, "%s: Error creating file %s\n", GMT_program, filename);
						exit (EXIT_FAILURE);
					}
				}
				else
					GMT_write_segmentheader (fp_out, n_out);
				for (k = 0; k < n_files; k++) {	/* Get next record */
					n_fields[k] = GMT_input (fp[k],  &n_expected_fields[k], &in);
					if (!GMT_REC_IS_NEW_SEGMENT) {	/* Data record */
						if (! all_set) {	/* Very first data record */
							all_set = 1;
							n_tot_cols += n_expected_fields[k];
							if (n_tot_cols > n_alloc) {
								n_alloc <<= 1;
								val = (double *) GMT_memory ((char *)val, (size_t)n_alloc, sizeof (double), GMT_program);
							}
						}
						for (i = 0; i < n_fields[k]; i++) val[j++] = in[i];	/* Copy the input values */
					}
					else {	/* Segment header again */
						if (k < (n_files-1)) GMT_io.segment_header[strlen(GMT_io.segment_header)-1] = '\0';
						if (k == 0)
							strcpy (header, GMT_io.segment_header);
						else
							strcat (header, &GMT_io.segment_header[1]);
					}
				}
				n_total_read ++;
			}

			if (all_set == 1) {	/* Well, we still have not allocated space etc */
				all_set = 2;
				if (n_out == 0) {	/* No cut action, just paste */
					n_out = n_tot_cols;
					for (k = 0; k < n_tot_cols; k++) Ctrl->F.col[k] = k;
					if (n_hdr_delay) {	/* Delayed issue of segment header record(s) */
						for (k = 0; k < n_hdr_delay; k++) GMT_write_segmentheader (GMT_stdout, n_out);
						n_hdr_delay = 0;
					}
				}
				out = (double *) GMT_memory (VNULL, (size_t)n_out, sizeof (double), GMT_program);
				if (n_tot_cols < n_out) {
					fprintf (stderr, "%s: Columns requested (%ld) exceed actual (%ld) - aborting\n", GMT_program, n_out, n_tot_cols);
					exit (EXIT_FAILURE);
				}
			}
			while (!GMT_REC_IS_NEW_SEGMENT && !GMT_REC_IS_EOF) {	/* Keep going until FALSE or = 2 segment header */
				if (GMT_REC_IS_ERROR) {
					fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld (skipped)\n", GMT_program, n_fields[0], n_expected_fields[0], n_total_read);
					continue;
				}

				/* Now shuffle and output the chosen columns */

				for (i = n_nan = 0; i < n_out; i++) {
					out[i] = val[Ctrl->F.col[i]];
					if (GMT_is_dnan (out[i])) n_nan++;
				}

				if (Ctrl->I.active) {	/* Collect all rows an write in reverse order */
					if (n_rev_recs == n_rev_alloc) {
						n_rev_alloc = (n_rev_recs == 0) ? GMT_CHUNK : (n_rev_alloc << 1);
						R_out = (double *) GMT_memory ((void *)R_out, (size_t)(n_rev_alloc * n_out), sizeof (double), GMT_program);
					}
					if (!Ctrl->N.active || n_nan < n_out) {
						for (i = 0; i < n_out; i++) R_out[n_rev++] = out[i];
						n_rev_recs++;
					}
				}
				else if (!Ctrl->N.active || n_nan < n_out)
					GMT_output (fp_out, n_out, out);

				for (k = j = 0; k < n_files; k++) {	/* Get next record */
					n_fields[k] = GMT_input (fp[k], &n_expected_fields[k], &in);
					if (GMT_REC_IS_SEG_HEADER) {	/* Another segment header */
						if (k < (n_files-1)) GMT_io.segment_header[strlen(GMT_io.segment_header)-1] = '\0';
						if (k == 0)
							strcpy (header, GMT_io.segment_header);
						else
							strcat (header, &GMT_io.segment_header[1]);
					}
					else 
						for (i = 0; i < n_fields[k]; i++) val[j++] = in[i];	/* Copy the input values */
				}
				n_total_read ++;
			}
		}

		for (k = 0; k < n_files; k++) GMT_fclose(fp[k]);
		GMT_free ((void *)fp);
		GMT_free ((void *)n_expected_fields);
		GMT_free ((void *)n_fields);
		GMT_free ((void *)val);

		if (gmtdefs.verbose) fprintf(stderr, "gmtconvert: %ld records passed (input cols = %ld; output cols = %ld)\n", n_total_read, n_tot_cols, n_out);
	}
	else {	/* Just cutting out columns */
		GMT_LONG n_expected_fields, n_fields, n_args;
		GMT_LONG nofile = TRUE, done = FALSE;
		FILE *fp = NULL;

		if (GMT_io.binary[GMT_IN] && n_out == 0) {
			n_out = GMT_io.ncol[GMT_IN];
			for (k = 0; k < n_out; k++) Ctrl->F.col[k] = k;
			out = (double *) GMT_memory (VNULL, (size_t)n_out, sizeof (double), GMT_program);
		}
		else if (n_out)	/* Set by -F */
			out = (double *) GMT_memory (VNULL, (size_t)n_out, sizeof (double), GMT_program);

		if (n_files > 0)
			nofile = FALSE;
		else
			n_files = 1;

		n_args = (argc > 1) ? argc : 2;	/* So it will go through the loop at least once (when we pipe) */

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

			if (!GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
				for (i = 0; i < GMT_io.n_header_recs; i++) {
					GMT_fgets (record_str, BUFSIZ, fp);
					n_total_read ++;
					if (!GMT_io.binary[GMT_OUT] && GMT_io.io_header[GMT_OUT]) GMT_fputs (record_str, GMT_stdout);
				}
			}

			n_expected_fields = (GMT_io.ncol[GMT_IN]) ? GMT_io.ncol[GMT_IN] : GMT_MAX_COLUMNS; /* GMT_MAX_COLUMNS -> Ascii table, must first find the number of fields */
			n_fields = GMT_input (fp, &n_expected_fields, &in);
			if (GMT_IO_GAP_CHECKING) GMT_write_segmentheader (GMT_stdout, n_out);
			n_seg = n_in_seg = 0;
			while (!GMT_REC_IS_EOF) {	/* Not yet EOF */

				while (GMT_REC_IS_NEW_SEGMENT && !GMT_REC_IS_EOF) {
					if (Ctrl->E.active && (Ctrl->E.mode & 2) && n_in_seg > 0) GMT_fputs (record_str, GMT_stdout);		/* Previous segments last point */
					if (Ctrl->S.active) {	/* See if this header has text matching our search string */
						match = (strstr (GMT_io.segment_header, Ctrl->S.pattern) != NULL);	/* TRUE if we matched */
						on = (Ctrl->S.inverse == !match);	/* If match then we turn output ON */
					}
					if (ascii_in_bin_out && n_out == 0)
						n_hdr_delay++;	/* Unable to write header yet since we dont know the field count... */
					else if ((Ctrl->L.active || on) && Ctrl->D.active) {
						out_seg++;
						sprintf (filename, Ctrl->D.name, (int)out_seg);
						if (fp_out != NULL) GMT_fclose (fp_out);
						if ((fp_out = GMT_fopen (filename, GMT_io.w_mode)) == NULL ) {
							fprintf (stderr, "%s: Error creating file %s\n", GMT_program, filename);
							exit (EXIT_FAILURE);
						}
					}
					else if (Ctrl->L.active || on)
						GMT_write_segmentheader (GMT_stdout, n_out);
					n_seg++;
					if (on) n_out_seg++;
					n_in_seg = 0;
					n_fields = GMT_input (fp,  &n_expected_fields, &in);
					n_total_read ++;
				}
				if (GMT_REC_IS_GAP) {
					GMT_write_segmentheader (GMT_stdout, n_out);
					GMT_io.status = 0;	/* Done with gap */
				}

				while (!GMT_REC_IS_LINE_BREAK && !GMT_REC_IS_EOF) {	/* Keep going until FALSE or = 2 segment header */
					if (GMT_REC_IS_ERROR) {
						fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld\n", GMT_program, n_fields, n_expected_fields, n_total_read);
						exit (EXIT_FAILURE);
					}
					if (n_out == 0) {
						n_out = n_fields;
						if (n_hdr_delay) {	/* Delayed issue of segment header record(s) */
							for (k = 0; k < n_hdr_delay; k++) GMT_write_segmentheader (GMT_stdout, n_out);
							n_hdr_delay = 0;
						}
						for (k = 0; k < n_out; k++) Ctrl->F.col[k] = k;
						out = (double *) GMT_memory (VNULL, (size_t)n_out, sizeof (double), GMT_program);
					}
					if (n_expected_fields < n_out) {
						fprintf (stderr, "%s: Columns requested (%ld) exceed actual (%ld) - aborting\n", GMT_program, n_expected_fields, n_out);
						exit (EXIT_FAILURE);
					}

					/* Now output the chosen columns */

					for (i = n_nan = 0; i < n_out; i++) {
						out[i] = in[Ctrl->F.col[i]];
						if (GMT_is_dnan (out[i])) n_nan++;
					}

					if (Ctrl->I.active) {	/* Collect all rows an write in reverse order */
						if (n_rev_recs == n_rev_alloc) {
							n_rev_alloc = (n_rev_recs == 0) ? GMT_CHUNK : (n_rev_alloc << 1);
							R_out = (double *) GMT_memory ((void *)R_out, (size_t)(n_rev_alloc * n_out), sizeof (double), GMT_program);
						}
						if (!Ctrl->N.active || n_nan < n_out) {	/* OK to use this record */
							for (i = 0; i < n_out; i++) R_out[n_rev++] = out[i];
							n_rev_recs++;
						}
					}
					else if (Ctrl->E.active) {
						if (n_in_seg == 0 && (Ctrl->E.mode & 1)) GMT_fputs (GMT_io.current_record, GMT_stdout);
						strcpy (record_str, GMT_io.current_record);
						n_in_seg++;
					}
					else if (on)
						if (!Ctrl->N.active || n_nan < n_out) GMT_output (fp_out, n_out, out);

					n_total_read ++;

					n_fields = GMT_input (fp, &n_expected_fields, &in);
				}
			}
			if (Ctrl->E.active && (Ctrl->E.mode & 2) && n_in_seg > 0) GMT_fputs (record_str, GMT_stdout);

			if (fp != GMT_stdin) GMT_fclose(fp);
		}
		if (gmtdefs.verbose) fprintf(stderr, "%s: %ld records passed (input cols = %ld; output cols = %ld)\n", GMT_program, n_total_read, n_expected_fields, n_out);
		if (gmtdefs.verbose && Ctrl->S.active) fprintf (stderr, "%s: Extracted %ld from a total of %ld segments\n", argv[0], n_out_seg, n_seg);
	}

	if (Ctrl->I.active) {	/* Time to write out rows in reverse order */
		for (j = n_rev_recs - 1; j >= 0; j--) GMT_output (fp_out, n_out, &R_out[j*n_out]);
		GMT_free ((void *)R_out);
	}

	if (Ctrl->D.active && fp_out != NULL) GMT_fclose (fp_out);
	
	GMT_free ((void *)out);

	Free_gmtconvert_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void *New_gmtconvert_Ctrl () {	/* Allocate and initialize a new control structure */
	struct GMTCONVERT_CTRL *C;
	
	C = (struct GMTCONVERT_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GMTCONVERT_CTRL), "New_gmtconvert_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	
	C->F.col = (GMT_LONG *) GMT_memory (VNULL, (size_t)GMT_MAX_COLUMNS, sizeof (GMT_LONG), GMT_program);
	
	return ((void *)C);
}

void Free_gmtconvert_Ctrl (struct GMTCONVERT_CTRL *C) {	/* Deallocate control structure */
	if (C->D.name) free ((void *)C->D.name);	
	if (C->F.col) GMT_free ((void *)C->F.col);	
	if (C->S.pattern) free ((void *)C->S.pattern);	
	GMT_free ((void *)C);	
}
