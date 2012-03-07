/*--------------------------------------------------------------------
*    $Id: gmtdp.c,v 1.14 2011/07/11 19:22:06 guru Exp $
*
*	Copyright (c) 1991-2011 by P. Wessel and W. H. F. Smith
*	See LICENSE.TXT file for copying and redistribution conditions.
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
* gmtdp applies the Douglas-Peucker algorithm to simplify a line
* segment given a tolerance.  The algorithm is based on the paper
* Douglas, D. H., and T. K. Peucker, Algorithms for the reduction
*   of the number of points required to represent a digitized line
*   of its caricature, Can. Cartogr., 10, 112-122, 1973.
* The implementation of this algorithm has been kindly provided by
* Dr. Gary J. Robinson, Environmental Systems Science Centre,
* University of Reading, Reading, UK (gazza@mail.nerc-essc.ac.uk)
*
* Assembler:   Joaquim Luis
* Date:        24-Aug-2008
*/

#include "gmt.h"

GMT_LONG Douglas_Peucker_geog (double x_source[], double y_source[], GMT_LONG n_source, double band, GMT_LONG index[]);

int main (int argc, char **argv) {
	GMT_LONG i, j, k, m, np, n_files = 0, n_args, error = 0;
	GMT_LONG n, *index = NULL;
	GMT_LONG n_alloc = GMT_SMALL_CHUNK;
	struct GMT_DATASET *D;
	double tolerance = 1.0, *x = NULL, *y = NULL;
	static double out[GMT_MAX_COLUMNS];
	char format[BUFSIZ];
	FILE *fp = NULL;
	
	argc = (int)GMT_begin (argc, argv);
	memset ((void *)format, 0, BUFSIZ);
	
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {

				/* Common parameters */

				case 'H':
				case 'M':
				case 'V':
				case ':':
				case 'b':
				case 'f':
				case 'm':
				case '\0':
					error += GMT_parse_common_options (argv[i], NULL, NULL, NULL, NULL);
					break;

				/* Supplemental parameters */

				case 'T':
					tolerance = atof(&argv[i][2]);
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
		fprintf (stderr, "gmtdp %s - Line reduction using the Douglas-Peucker algorithm\n\n", GMT_VERSION);
		fprintf (stderr, "usage: gmtdp <infiles> [%s] -T<tolerance>\n", GMT_H_OPT);
		fprintf (stderr, "\t[-V[l]] [%s] [%s] [%s]\n\n", GMT_t_OPT, GMT_b_OPT, GMT_m_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\tinfiles (in ASCII or binary) have 2 or more columns with (x,y) or (y,x) in first columns.\n");
		fprintf (stderr, "\t  If no file(s) is given, standard input is read.\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-T tolerance is maximum mismatch in km.\n");
		GMT_explain_option ('H');
		GMT_explain_option ('V');
		GMT_explain_option (':');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf (stderr, "\t   Default is 2 input columns\n");
		GMT_explain_option ('o');
		GMT_explain_option ('n');
		GMT_explain_option ('m');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	if (n_files == 0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  No input files specified\n", GMT_program);
		error++;
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = 2;
	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] < 2) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data (-bi) must have at least %d columns\n", GMT_program, 2);
		error++;
	}
	if (error) exit (EXIT_FAILURE);

	/* Now we are ready to take on some input values */

	n_args = (argc > 1) ? argc : 2;

	/* Allocate memory and read in all the files; each file can have many lines (-M) */
	
	D = (struct GMT_DATASET *) GMT_memory (VNULL, (size_t)1, sizeof (struct GMT_DATASET), GMT_program);
	D->table = (struct GMT_TABLE **) GMT_memory (VNULL, (size_t)n_alloc, sizeof (struct GMT_TABLE *), GMT_program);
	for (k = 1, i = 0; k < argc; k++) {
		if (argv[k][0] == '-') continue;

		if (gmtdefs.verbose) fprintf (stderr, "%s: Reading file %s\n", GMT_program, argv[k]);

		GMT_import_table ((void *)argv[k], GMT_IS_FILE, &D->table[i], 0.0, FALSE, FALSE, TRUE);
		i++;
		if (i == n_alloc) {
			n_alloc <<= 1;
			D->table = (struct GMT_TABLE **) GMT_memory ((void *)D->table, (size_t)n_alloc, sizeof(struct GMT_TABLE *), GMT_program);
		}
	}
	if (i < n_alloc) D->table = (struct GMT_TABLE **) GMT_memory ((void *)D->table, (size_t)i, sizeof(struct GMT_TABLE *), GMT_program);
	D->n_tables = i;
	
	fp = GMT_stdout;
	for (k = 0; k < D->n_tables; k++) {
		for (j = 0; j < D->table[k]->n_segments; j++) {
			np = D->table[k]->segment[j]->n_rows;

			index = (GMT_LONG *) GMT_memory(VNULL, (size_t)np, sizeof(GMT_LONG), GMT_program);
			x = (double *) GMT_memory(VNULL, (size_t)np, sizeof(double), GMT_program);
			y = (double *) GMT_memory(VNULL, (size_t)np, sizeof(double), GMT_program);
			for (n = 0; n < np; n++) {
				x[n] = D->table[k]->segment[j]->coord[0][n];
				y[n] = D->table[k]->segment[j]->coord[1][n];
			}

			np = Douglas_Peucker_geog (x, y, np, tolerance, index);

			if (GMT_io.multi_segments[GMT_OUT]) {
				if (D->table[k]->segment[j]->header) strcpy (GMT_io.segment_header, D->table[k]->segment[j]->header);
				GMT_write_segmentheader (fp, D->table[k]->segment[j]->n_columns);
			}

			for (n = 0; n < np; n++) {
				out[0] = x[index[n]];		out[1] = y[index[n]];
				for (m = 2; m < D->table[k]->segment[j]->n_columns; m++)
					out[m] = D->table[k]->segment[j]->coord[m][index[n]];

				GMT_output (fp, D->table[k]->segment[j]->n_columns, out);
			}
		}
	}
	
	GMT_free ((void *)index);
	GMT_free ((void *)x);
	GMT_free ((void *)y);
	GMT_free_dataset (D);
	return (EXIT_SUCCESS);
}

/* -------------------------------------------------------------------------------------------- */
#define sqr(x) ((x)*(x))

/* Stack-based Douglas Peucker line simplification routine */
/* returned value is the number of output points */

GMT_LONG Douglas_Peucker_geog (double x_source[], double y_source[], GMT_LONG n_source, double band, GMT_LONG index[]) {
/* x/y_source	Input coordinates, n_source of them */
/* band;		tolerance in km */
/* index[]	output co-ordinates indices */

	GMT_LONG	n_stack, n_dest, start, end, i, sig;
	GMT_LONG	*sig_start, *sig_end;	/* indices of start&end of working section */

	double dev_sqr, max_dev_sqr, band_sqr;
	double  x12, y12, d12, x13, y13, d13, x23, y23, d23;

	/* check for simple cases */

	if ( n_source < 3 ) {     /* one or two points */
		for ( i = 0; i < n_source; i++) index[i] = i;
		return (n_source);
	}

	/* more complex case. initialise stack */

	sig_start = (GMT_LONG *) GMT_memory (NULL, n_source, sizeof (int), "Douglas_Peucker");
	sig_end   = (GMT_LONG *) GMT_memory (NULL, n_source, sizeof (int), "Douglas_Peucker");
	
	band *= 360.0 / (2.0 * M_PI * 6371.007181);	/* Now in degrees */
	band_sqr = sqr(band);

	n_dest = 0;

	sig_start[0] = 0;
	sig_end[0] = n_source-1;

	n_stack = 1;

	/* while the stack is not empty  ... */

	while ( n_stack > 0 ) {
		/* ... pop the top-most entries off the stacks */

		start = sig_start[n_stack-1];
		end = sig_end[n_stack-1];

		n_stack--;

		if ( end - start > 1 ) { /* any intermediate points ? */
			/* ... yes, so find most deviant intermediate point to
			   either side of line joining start & end points */

			x12 = x_source[end] - x_source[start];
			if (fabs (x12) > 180.0) x12 = 360.0 - fabs (x12);
			y12 = y_source[end] - y_source[start];
			x12 *= cosd (0.5 * (y_source[end] + y_source[start]));
			d12 = sqr(x12) + sqr(y12);

			for ( i = start + 1, sig = start, max_dev_sqr = -1.0; i < end; i++ ) {
				x13 = x_source[i] - x_source[start];
				if (fabs (x13) > 180.0) x13 = 360.0 - fabs (x13);
				y13 = y_source[i] - y_source[start];

				x23 = x_source[i] - x_source[end];
				if (fabs (x23) > 180.0) x23 = 360.0 - fabs (x23);
				y23 = y_source[i] - y_source[end];

				x13 *= cosd (0.5 * (y_source[i] + y_source[end]));
				x23 *= cosd (0.5 * (y_source[i] + y_source[end]));

				d13 = sqr(x13) + sqr(y13);
				d23 = sqr(x23) + sqr(y23);

				if ( d13 >= ( d12 + d23 ) )
					dev_sqr = d23;
				else if ( d23 >= ( d12 + d13 ) )
					dev_sqr = d13;
				else
					dev_sqr =  sqr( x13 * y12 - y13 * x12 ) / d12;

				if ( dev_sqr > max_dev_sqr  ) {
					sig = i;
					max_dev_sqr = dev_sqr;
				}
			}

			if ( max_dev_sqr < band_sqr ) {  /* is there a sig.  intermediate point ? */
				/* ... no, so transfer current start point */
				index[n_dest] = start;
				n_dest++;
			}
			else {
				/* ... yes, so push two sub-sections on stack for further processing */

				n_stack++;

				sig_start[n_stack-1] = sig;
				sig_end[n_stack-1] = end;

				n_stack++;

				sig_start[n_stack-1] = start;
				sig_end[n_stack-1] = sig;
			}
		}
		else {
			/* ... no intermediate points, so transfer current start point */
			index[n_dest] = start;
			n_dest++;
		}
	}


	/* transfer last point */

	index[n_dest] = n_source-1;
	n_dest++;

	GMT_free ((void *)sig_start);
	GMT_free ((void *)sig_end);

	return(n_dest);
}
