/*--------------------------------------------------------------------
 *	$Id: grdhisteq.c 10173 2014-01-01 09:52:34Z pwessel $
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
 * read a grid file and find the values which divide its range
 * into n_cell number of quantiles.
 *
 * Author:	W.H.F. Smith
 * Date: 	31 May 1990
 *
 * Modified:	12 June, 1990 by whf smith, adding [-Q] option for
 *	quadratic scaling.  Some rgb color systems consider that
 *	if black = (0,0,0) and white = (1,1,1) or (255,255,255),
 *	then a neutral gray "halfway" between black and while should
 *	be set to gray = (0.75,0.75,0.75) or (191,191,191).  If so,
 *	let 0 <= x <= 1 be the desired gradation between black and
 *	white (the intensity factor used by the coloring program.
 *	Then the gray tone level 0 <= y <= 1 is given by:
 *		y = 2*x - x**2.
 *	Using the -Q option will find the data values which divide
 *	the data range into <n_cells> values of y; default linear
 *	scaling will find the values for <n_cells> divisions of x.
 *
 * Updated to v2.0 15-May-1991 Paul Wessel
 * Updated to v3.1 14-Jun-1998 Paul Wessel
 * Updated to v3.3.5 14-Jun-2000 Paul Wessel
 * Version:	4
 */
 
#include "gmt.h"

struct GRDHISTEQ_CTRL {
	struct C {	/* -C<n_cells>*/
		GMT_LONG active;
		GMT_LONG value;
	} C;
	struct D {	/* -D */
		GMT_LONG active;
	} D;
	struct G {	/* -G<file> */
		GMT_LONG active;
		char *file;
	} G;
	struct N {	/* -N[<norm>] */
		GMT_LONG active;
		double norm;
	} N;
	struct Q {	/* -Q */
		GMT_LONG active;
	} Q;
};

struct	INDEXED_DATA {
	float x;
	GMT_LONG i;
};

struct	CELL {
	float low;
	float high;
};

int main (int argc, char **argv)
{
	GMT_LONG	error = FALSE;
	
	char *infile = CNULL;

	GMT_LONG i;
		
	struct GRDHISTEQ_CTRL *Ctrl = NULL;

	void do_usual (char *infile, char *outfile, GMT_LONG n_cells, GMT_LONG quadratic, GMT_LONG dump_intervals, int argc, char **argv);
	void do_gaussian (char *infile, char *outfile, double norm, int argc, char **argv);
	void *New_grdhisteq_Ctrl (), Free_grdhisteq_Ctrl (struct GRDHISTEQ_CTRL *C);

	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct GRDHISTEQ_CTRL *)New_grdhisteq_Ctrl ();	/* Allocate and initialize a new control structure */

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {

				/* Common parameters */

				case 'V':
				case '\0':
					error += GMT_parse_common_options (argv[i], 0, 0, 0, 0);
					break;

				/* Supplemental parameters */

				case 'C':
					Ctrl->C.active = TRUE;
					Ctrl->C.value = atoi(&argv[i][2]);
					break;
				case 'D':
					Ctrl->D.active = TRUE;
					break;
				case 'G':
					Ctrl->G.active = TRUE;
					Ctrl->G.file = strdup (&argv[i][2]);
					break;
				case 'N':
					Ctrl->N.active = TRUE;
					if (argv[i][2]) Ctrl->N.norm = atof (&argv[i][2]);
					break;
				case 'Q':
					Ctrl->Q.active = TRUE;
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

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr,"grdhisteq %s - Histogram equalization for grid files\n\n", GMT_VERSION);
		fprintf (stderr, "usage: grdhisteq <infile> -G<outfile> [-C<n_cells> -D -N[<norm>] -Q -V]\n");
		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);
		fprintf (stderr, "\t-C<n_cells> sets how many cells (divisions) of data range to make.\n");
		fprintf (stderr, "\t-D dump level information to stdout.\n");
		fprintf (stderr, "\t-G<outfile> will create an equalized output grid file.\n");
		fprintf (stderr, "\t-N use with -G to make an output grid file with standard normal scores.\n");
		fprintf (stderr, "\t   Append <norm> to normalize the scores to <-1,+1>.\n");
		fprintf (stderr, "\t-Q to use quadratic intensity scaling [Default is linear].\n");
		GMT_explain_option ('V');
		exit (EXIT_FAILURE);
	}

	if (!infile) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify input file\n", GMT_program);
		error++;
	}
	if (Ctrl->N.active && !Ctrl->G.file) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -N option:  Must also specify output file with -G\n", GMT_program);
		error++;
	}
	if (!Ctrl->N.active && Ctrl->C.value <= 0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -C option:  n_cells must be positive\n", GMT_program);
		error++;
	}
	if (error) exit (EXIT_FAILURE);

	if (!strcmp (infile, "=")) {
		fprintf (stderr, "%s: Piping of input grid file not supported!\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	if (Ctrl->N.active)
		do_gaussian (infile, Ctrl->G.file, Ctrl->N.norm, argc, argv);
	else
		do_usual (infile, Ctrl->G.file, Ctrl->C.value, Ctrl->Q.active, Ctrl->D.active, argc, argv);

	Free_grdhisteq_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void do_usual (char *infile, char *outfile, GMT_LONG n_cells, GMT_LONG quadratic, GMT_LONG dump_intervals, int argc, char **argv)
{
	char format[BUFSIZ];
	GMT_LONG	last_cell, n_cells_m1 = 0, current_cell, i, j, nxy, nxy_0;
	float	*data = NULL;
	double	delta_cell, target;
	struct GRD_HEADER header;
	struct	CELL *cell = NULL;
	
	float get_cell (float x, struct CELL *cell, GMT_LONG n_cells_m1, GMT_LONG last_cell);

	sprintf (format, "%s%s%s%s%%ld\n", gmtdefs.d_format, gmtdefs.field_delimiter, gmtdefs.d_format, gmtdefs.field_delimiter);

	GMT_err_fail (GMT_read_grd_info (infile, &header), infile);

	GMT_grd_init (&header, argc, argv, TRUE);

	nxy_0 = GMT_get_nm (header.nx, header.ny);
	
	data = (float *) GMT_memory (VNULL, (size_t)nxy_0, sizeof (float), GMT_program);
	GMT_err_fail (GMT_read_grd (infile, &header, data, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), infile);

	cell = (struct CELL *) GMT_memory (VNULL, (size_t)n_cells, sizeof(struct CELL), GMT_program);

	/* Sort the data and find the division points:  */

	qsort ((void *)data, (size_t)nxy_0, sizeof(float), GMT_comp_float_asc);
	nxy = nxy_0;
	while (nxy > 0 && GMT_is_fnan (data[nxy-1])) nxy--;	/* Only deal with real numbers */

	last_cell = n_cells/2;
	n_cells_m1 = n_cells - 1;

	current_cell = i = 0;
	delta_cell = ((double)nxy) / ((double)n_cells);

	while (current_cell < n_cells) {

		if (current_cell == (n_cells - 1) ) {
			j = nxy - 1;
		}
		else if (quadratic) {	/* Use y = 2x - x**2 scaling  */

			target = ( (double) (current_cell + 1) ) / ( (double) n_cells);
			j = (GMT_LONG)floor(nxy * (1.0 - sqrt(1.0 - target)));
		}
		else {	/* Use simple linear scale  */

			j = (GMT_LONG)(floor( (current_cell + 1) * delta_cell)) - 1;
		}

		cell[current_cell].low = data[i];
		cell[current_cell].high = data[j];

		if (dump_intervals) fprintf (stdout, format, data[i], data[j], current_cell);

		i = j;
		current_cell++;
	}

	if (outfile) {
		GMT_err_fail (GMT_read_grd (infile, &header, data, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), infile);

		for (i = 0; i < nxy_0; i++) data[i] = (GMT_is_fnan (data[i])) ? GMT_f_NaN : get_cell (data[i], cell, n_cells_m1, last_cell);

		GMT_err_fail (GMT_write_grd (outfile, &header, data, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), outfile);
	}

	GMT_free ((void *) data);
	GMT_free ((void *) cell);
}

float get_cell (float x, struct CELL *cell, GMT_LONG n_cells_m1, GMT_LONG last_cell)
{
	GMT_LONG	low, high, i;

	low = 0;
	high = n_cells_m1;
	i = last_cell;

	do {
		if (cell[i].low <= x && cell[i].high >= x) {
			last_cell = i;
			return ( (float)i);
		}
		else if (cell[low].low <= x && cell[low].high >= x) {
			return ( (float)low);
		}
		else if (cell[high].low <= x && cell[high].high >= x) {
			return ( (float)high);
		}
		else if (cell[i].low > x) {
			high = i;
			i = (low + high) / 2;
		}
		else if (cell[i].high < x) {
			low = i;
			i = (low + high) / 2;
		}
	} while (TRUE);
	return (0.0);	/* Cannot get here - just used to quiet compiler */
}

void do_gaussian (char *infile, char *outfile, double norm, int argc, char **argv)
{
	GMT_LONG	i, j, nxy, nxy_0;
	float	*data = NULL;
	double	dnxy;
	struct	GRD_HEADER header;
	struct	INDEXED_DATA *indexed_data = NULL;
	int	compare_indexed_floats(const void *point_1, const void *point_2);
	int	compare_indices(const void *point_1, const void *point_2);

	GMT_err_fail (GMT_read_grd_info (infile, &header), infile);

	GMT_grd_init (&header, argc, argv, TRUE);

	nxy_0 = header.nx * header.ny;
	data = (float *) GMT_memory (VNULL, (size_t)nxy_0, sizeof (float), GMT_program);
	GMT_err_fail (GMT_read_grd (infile, &header, data, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), infile);

	indexed_data = (struct INDEXED_DATA *) GMT_memory (VNULL, (size_t)nxy_0, sizeof (struct INDEXED_DATA), GMT_program);

	for (i = j = 0, nxy = nxy_0; i < nxy_0; i++) {
		if (GMT_is_fnan (data[i])) {	/* Put NaNs in the back */
			nxy--;
			indexed_data[nxy].i = i;
			indexed_data[nxy].x = data[i];
		}
		else {
			indexed_data[j].i = i;
			indexed_data[j].x = data[i];
			j++;
		}
	}

	/* Sort on data value  */

	qsort ((void *)indexed_data, (size_t)nxy, sizeof(struct INDEXED_DATA), compare_indexed_floats);

	dnxy = 1.0 / (nxy + 1);

	if (norm != 0.0) norm /= fabs (GMT_zcrit ((double)dnxy));	/* Normalize by abs(max score) */

	for (i = 0; i < nxy; i++) {
		indexed_data[i].x = (float)GMT_zcrit ((double)((i + 1) * dnxy));
		if (norm != 0.0) indexed_data[i].x *= (float)norm;
	}

	/* Sort on data index  */

	qsort ((void *)indexed_data, (size_t)nxy_0, sizeof(struct INDEXED_DATA), compare_indices);

	for (i = 0; i < nxy_0; i++) data[i] = indexed_data[i].x;

	GMT_err_fail (GMT_write_grd (outfile, &header, data, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), outfile);

	GMT_free ((void *) indexed_data);
	GMT_free ((void *) data);
}

int compare_indexed_floats (const void *point_1, const void *point_2)
{
	if (((struct INDEXED_DATA *)point_1)->x < ((struct INDEXED_DATA *)point_2)->x)
		return (-1);
	else if (((struct INDEXED_DATA *)point_1)->x > ((struct INDEXED_DATA *)point_2)->x)
		return (1);
	else 
		return (0);
}

int compare_indices (const void *point_1, const void *point_2)
{
	if (((struct INDEXED_DATA *)point_1)->i < ((struct INDEXED_DATA *)point_2)->i)
		return (-1);
	else if (((struct INDEXED_DATA *)point_1)->i > ((struct INDEXED_DATA *)point_2)->i)
		return (1);
	else 
		return (0);
}

void *New_grdhisteq_Ctrl () {	/* Allocate and initialize a new control structure */
	struct GRDHISTEQ_CTRL *C;
	
	C = (struct GRDHISTEQ_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GRDHISTEQ_CTRL), "New_grdhisteq_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	return ((void *)C);
}

void Free_grdhisteq_Ctrl (struct GRDHISTEQ_CTRL *C) {	/* Deallocate control structure */
	if (C->G.file) free ((void *)C->G.file);	
	GMT_free ((void *)C);	
}
