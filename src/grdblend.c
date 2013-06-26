/*--------------------------------------------------------------------
 *    $Id: grdblend.c 9923 2012-12-18 20:45:53Z pwessel $
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
 * grdblend reads any number of grid files that may partly overlap and
 * creates a blend of all the files given certain criteria.  Each input
 * grid is considered to have an "outer" and "inner" region.  The outer
 * region is the extent of the grid; the inner region is provided as
 * input in the form of a -Rw/e/s/n statement.  Finally, each grid can
 * be assigned its separate weight.  This information is given to the
 * program in ASCII format, one line per grid file; each line looks like
 *
 * grdfile	-Rw/e/s/n	weight
 *
 * The blending will use a 2-D cosine taper between the inner and outer
 * regions.  The output at any node is thus a weighted average of the
 * values from any grid that occupies that grid node.  Because the out-
 * put grid can be really huge (say global grids at fine resolution),
 * all grid input/output is done row by row so memory should not be a
 * limiting factor in making large grid.
 *
 * Author:	Paul Wessel
 * Date:	06-DEC-2001
 * Version:	4
 */

#include "gmt.h"

struct GRDBLEND_CTRL {
	struct G {	/* -G<grdfile> */
		GMT_LONG active;
		char *file;
	} G;
	struct I {	/* -Idx[/dy] */
		GMT_LONG active;
		double xinc, yinc;
	} I;
	struct N {	/* -N<nodata> */
		GMT_LONG active;
		double nodata;
	} N;
	struct Q {	/* -Q */
		GMT_LONG active;
	} Q;
	struct Z {	/* -Z<scale> */
		GMT_LONG active;
		double scale;
	} Z;
	struct W {	/* -W */
		GMT_LONG active;
	} W;
};

struct GRDBLEND_INFO {	/* Structure with info about each input grid file */
	struct GMT_GRDFILE G;				/* I/O structure for grid files, including grd header */
	GMT_LONG in_i0, in_i1, out_i0, out_i1;		/* Indices of outer and inner x-coordinates (in output grid coordinates) */
	GMT_LONG in_j0, in_j1, out_j0, out_j1;		/* Indices of outer and inner y-coordinates (in output grid coordinates) */
	GMT_LONG offset;				/* grid offset when the grid extends beyond north */
	long skip;					/* Byte offset to skip in native binary files */
	GMT_LONG outside;				/* TRUE if the current output row is outside the range of this grid */
	GMT_LONG invert;					/* TRUE if weight was given as negative and we want to taper to zero INSIDE the grid region */
	GMT_LONG open;					/* TRUE if file is currently open */
	char file[GMT_LONG_TEXT];			/* Name of grid file */
	double weight, wt_y, wxr, wxl, wyu, wyd;	/* Various weighting factors used for cosine-taper weights */
	double w_in, e_in, s_in, n_in;			/* Boundaries of inner region */
	float *z;					/* Row vector holding the current row from this file */
};

#define N_NOT_SUPPORTED	6

GMT_LONG found_unsupported_format (struct GRD_HEADER *h, char *text, char *file);

int main (int argc, char **argv)
{
	GMT_LONG i, col, pcol, row, nx_360 = 0, k, kk, m, n_blend, err, error = 0;
	GMT_LONG n_fill, n_tot, type;
	double wt_x, w, wt;
	float *z = NULL, no_data_f;
	char mode[2] = {'w', 'W'};
	GMT_LONG wrap_x;
	FILE *fp = NULL;
	struct GRDBLEND_INFO *blend = NULL;
	struct GMT_GRDFILE S;
	struct GRDBLEND_CTRL *Ctrl = NULL;

	void sync_input_rows (GMT_LONG row, struct GRDBLEND_INFO *blend, GMT_LONG n_blend, double half);
	GMT_LONG init_blend_job (FILE *fp, struct GRD_HEADER *h, struct GRDBLEND_INFO **blend);
	void *New_grdblend_Ctrl (), Free_grdblend_Ctrl (struct GRDBLEND_CTRL *C);

	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct GRDBLEND_CTRL *)New_grdblend_Ctrl ();	/* Allocate and initialize a new control structure */
	
	GMT_grd_init (&S.header, argc, argv, FALSE);

	n_fill = n_tot = 0;

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {

				/* Common parameters */

				case 'R':
				case 'V':
				case 'f':
				case ':':
				case '\0':
					error += GMT_parse_common_options (argv[i], &S.header.x_min, &S.header.x_max, &S.header.y_min, &S.header.y_max);
					break;

				/* Supplemental parameters */

				case 'I':
					Ctrl->I.active = TRUE;
					if (GMT_getinc (&argv[i][2], &Ctrl->I.xinc, &Ctrl->I.yinc)) {
						GMT_inc_syntax ('I', 1);
						error = TRUE;
					}
					break;
				case 'G':
					Ctrl->G.file = strdup (&argv[i][2]);
					Ctrl->G.active = TRUE;
					break;
				case 'N':
					if (!argv[i][2]) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -N option:  Must specify value or NaN\n", GMT_program);
						error++;
					}
					else {
						Ctrl->N.nodata = (argv[i][2] == 'N' || argv[i][2] == 'n') ? GMT_d_NaN : atof (&argv[i][2]);
					}
					Ctrl->N.active = TRUE;
					break;
				case 'Q':
					Ctrl->Q.active = TRUE;
					break;

				case 'W':
					Ctrl->W.active = TRUE;
					break;

				case 'Z':
					Ctrl->Z.active = TRUE;
					Ctrl->Z.scale = atof (&argv[i][2]);
					break;

				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else if ((fp = GMT_fopen (argv[i], "r")) == NULL) {
			fprintf (stderr, "%s: Cannot open file %s\n", GMT_program, argv[i]);
			exit (EXIT_FAILURE);
		}
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr,"grdblend %s - Blend several partially over-lapping grid files onto one grid\n\n", GMT_VERSION);
		fprintf (stderr, "usage: grdblend [<blendfile>] -G<grdfile> %s\n", GMT_I_OPT);
		fprintf (stderr, "\t%s [-N<nodata>] [-Q] [-Z<scale>] [-V] [-W]\n\t[%s]\n", GMT_Rgeo_OPT, GMT_f_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\t<blendfile> is an ASCII file (or stdin) with blending parameters for each input grid.\n");
		fprintf (stderr, "\t   Each record has three items:  filename -Rw/e/s/n weight.\n");
		fprintf (stderr, "\t   Relative weights are <weight> inside the given -R and cosine taper to 0 at actual grid -R.\n");
		fprintf (stderr, "\t   Give filename - weight if inner region should equal the actual region.\n");
		fprintf (stderr, "\t   Give a negative weight to invert the sense of the taper (i.e., |<weight>| outside given R.\n");
		fprintf (stderr, "\t   Each grid must be in netCDF or native binary format.\n");
		fprintf (stderr, "\t-G <grdfile> is the name of the output 2-D binary data set.\n");
		fprintf (stderr, "\t   Only netCDF or native binary grid formats are supported.\n");
		GMT_inc_syntax ('I', 0);
		GMT_explain_option ('R');
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-N Set value for nodes without constraints [Default is NaN].\n");
		fprintf (stderr, "\t-Q grdraster-compatible output without leading grd header [Default writes GMT grid file].\n");
		fprintf (stderr, "\t   Output grid must be in native binary format (i.e., not netCDF).\n");
		fprintf (stderr, "\t-Z Multiply z-values by this scale before writing to file [1].\n");
		GMT_explain_option ('V');
		fprintf (stderr, "\t-W Write out weights only (only applies to a single input file) [make blend grid].\n");
		GMT_explain_option ('f');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	GMT_check_lattice (&Ctrl->I.xinc, &Ctrl->I.yinc, NULL, &Ctrl->I.active);

	if (!project_info.region_supplied) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -R:  Must specify region\n", GMT_program);
		error++;
	}
	if (Ctrl->I.xinc <= 0.0 || Ctrl->I.yinc <= 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -I:  Must specify positive dx, dy\n", GMT_program);
		error++;
	}
	if (!Ctrl->G.file) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -G:  Must specify output file\n", GMT_program);
		error++;
	}
	if ((err = GMT_grd_get_format (Ctrl->G.file, &S.header, FALSE)) != GMT_NOERROR){
		fprintf (stderr, "%s: GMT SYNTAX ERROR: %s [%s]\n", GMT_program, GMT_strerror(err), Ctrl->G.file);
		exit (EXIT_FAILURE);
	}
	/* Only allow netcdf (both v3 and new) and native binary output */
	if (found_unsupported_format (&S.header, "GMT SYNTAX ERROR -G", Ctrl->G.file)) exit (EXIT_FAILURE);
	type = GMT_grdformats[S.header.type][0];
	if (Ctrl->Q.active && (type == 'c' || type == 'n')) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -Q:  Only native binary grids can skip writing the header.\n", GMT_program);
		error++;
	}

	if (error) exit (EXIT_FAILURE);

	S.header.x_inc = Ctrl->I.xinc;
	S.header.y_inc = Ctrl->I.yinc;
	
	GMT_RI_prepare (&S.header);	/* Ensure -R -I consistency and set nx, ny */

	if (fp == NULL) {	/* No file given, use standard input */
		fp = GMT_stdin;
		if (gmtdefs.verbose) fprintf (stderr, "%s: Reads blend parameters from standard input\n", GMT_program);
	}

	/* Process blend parameters and populate blend structure and open input files and seek to first row inside the output grid */

	no_data_f = (float)Ctrl->N.nodata;
	n_blend = init_blend_job (fp, &S.header, &blend);

	if (fp != GMT_stdin) GMT_fclose (fp);

	if (Ctrl->W.active && n_blend > 1) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -W:  Only applies when there is a single input grid file\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	/* Initialize header structure for output blend grid */

	S.header.z_min = S.header.z_max = 0.0;

	GMT_err_fail (GMT_grd_RI_verify (&S.header, 1), Ctrl->G.file);

	n_tot = GMT_get_nm (S.header.nx, S.header.ny);

	z = (float *) GMT_memory (VNULL, (size_t)S.header.nx, sizeof (float), GMT_program);	/* Memory for one output row */

	if (!Ctrl->Q.active) GMT_err_fail (GMT_write_grd_info (Ctrl->G.file, &S.header), Ctrl->G.file);

	GMT_err_fail (GMT_open_grd (Ctrl->G.file, &S, mode[Ctrl->Q.active]), Ctrl->G.file);
	
	if (gmtdefs.verbose && Ctrl->Z.active) {
		fprintf (stderr, "%s: Output data will be scaled by %g\n", GMT_program, Ctrl->Z.scale);
	}

	S.header.z_min = DBL_MAX;	/* These will be updated in the loop below */
	S.header.z_max = -DBL_MAX;
	wrap_x = (GMT_io.out_col_type[GMT_X] == GMT_IS_LON);	/* Periodic geographic grid */
	if (wrap_x) nx_360 = irint (360.0 / S.header.x_inc);

	for (row = 0; row < S.header.ny; row++) {	/* For every output row */

		memset ((void *)z, 0, (size_t)(S.header.nx * sizeof (float)));	/* Start from scratch */

		sync_input_rows (row, blend, n_blend, S.header.xy_off);	/* Wind each input file to current record and read each of the overlapping rows */

		for (col = 0; col < S.header.nx; col++) {	/* For each output node on the current row */

			w = 0.0;
			for (k = m = 0; k < n_blend; k++) {	/* Loop over every input grid */
				if (blend[k].outside) continue;					/* This grid is currently outside the s/n range */
				if (wrap_x) {	/* Special testing for periodic x coordinates */
					pcol = col + nx_360;
					while (pcol > blend[k].out_i1) pcol -= nx_360;
					if (pcol < blend[k].out_i0) continue;	/* This grid is currently outside the w/e range */
				}
				else {	/* Not periodic */
					if (col < blend[k].out_i0 || col > blend[k].out_i1) continue;	/* This grid is currently outside the xmin/xmax range */
					pcol = col;
				}
				kk = pcol - blend[k].out_i0;					/* kk is the local column variable for this grid */
				if (GMT_is_fnan (blend[k].z[kk])) continue;			/* NaNs do not contribute */
				if (pcol <= blend[k].in_i0)					/* Left cosine-taper weight */
					wt_x = 0.5 * (1.0 - cos ((pcol - blend[k].out_i0 + S.header.xy_off) * blend[k].wxl));
				else if (pcol >= blend[k].in_i1)					/* Right cosine-taper weight */
					wt_x = 0.5 * (1.0 - cos ((blend[k].out_i1 - pcol + S.header.xy_off) * blend[k].wxr));
				else								/* Inside inner region, weight = 1 */
					wt_x = 1.0;
				wt = wt_x * blend[k].wt_y;					/* Actual weight is 2-D cosine taper */
				if (blend[k].invert) wt = blend[k].weight - wt;			/* Invert the sense of the tapering */
				z[col] += (float)(wt * blend[k].z[kk]);				/* Add up weighted z*w sum */
				w += wt;							/* Add up the weight sum */
				m++;								/* Add up the number of contributing grids */
			}

			if (m) {		/* OK, at least one grid contributed to an output value */
				if (!Ctrl->W.active) {		/* Want output z blend */
					z[col] = (float)((w == 0.0) ? 0.0 : z[col] / w);	/* Get weighted average z */
					if (Ctrl->Z.active) z[col] *= (float)Ctrl->Z.scale;		/* Apply the global scale here */
				}
				else		/* Get the weight only */
					z[col] = (float)w;				/* Only interested in the weights */
				n_fill++;						/* One more cell filled */
				if (z[col] < S.header.z_min) S.header.z_min = z[col];	/* Update the extrema for output grid */
				if (z[col] > S.header.z_max) S.header.z_max = z[col];
			}
			else			/* No grids covered this node, defaults to the no_data value */
				z[col] = no_data_f;
		}
		GMT_err_fail (GMT_write_grd_row (&S, 0, z), Ctrl->G.file);

		if (gmtdefs.verbose && row%10 == 0)  fprintf (stderr, "%s: Processed row %7ld of %d\r", GMT_program, row, S.header.ny);

	}
	if (gmtdefs.verbose)  fprintf (stderr, "%s: Processed row %7ld\n", GMT_program, row);

	GMT_close_grd (&S);	/* Close the output gridfile */
	if (!Ctrl->Q.active) GMT_err_fail (GMT_update_grd_info (Ctrl->G.file, &S.header), Ctrl->G.file);
	GMT_free ((void *)z);

	for (k = 0; k < n_blend; k++) if (blend[k].open) {
		GMT_free ((void *)blend[k].z);
		GMT_close_grd (&blend[k].G);	/* Close all input grd files still open */
	}

	if (gmtdefs.verbose) {
		char empty[GMT_TEXT_LEN];
		fprintf (stderr, "%s: Blended grid size of %s is %d x %d\n", GMT_program, Ctrl->G.file, S.header.nx, S.header.ny);
		if (n_fill == n_tot)
			fprintf (stderr, "%s: All nodes assigned values\n", GMT_program);
		else {
			if (GMT_is_fnan (no_data_f))
				strcpy (empty, "NaN");
			else
				sprintf (empty, "%g", no_data_f);
			fprintf (stderr, "%s: %ld nodes assigned values, %ld set to %s\n", GMT_program, (GMT_LONG)n_fill, (GMT_LONG)(n_tot - n_fill), empty);
		}
	}

	GMT_free ((void *)blend);

	Free_grdblend_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

GMT_LONG found_unsupported_format (struct GRD_HEADER *h, char *text, char *file)
{
	GMT_LONG i;
	static char *not_supported[N_NOT_SUPPORTED] = {"rb", "rf", "sf", "sd", "af", "gd"};
	for (i = 0; i < N_NOT_SUPPORTED; i++) {	/* Only allow netcdf (both v3 and new) and native binary output */
		if (h->type == GMT_grd_format_decoder (not_supported[i])) {
			fprintf (stderr, "grdblend: %s:  Grid format type %s for file %s is not supported in grdblend\n", text, not_supported[i], file);
			return (1);
		}
	}
	return (0);
}
GMT_LONG init_blend_job (FILE *fp, struct GRD_HEADER *h, struct GRDBLEND_INFO **blend) {
	GMT_LONG n = 0, nr, one_or_zero = 0, n_alloc = 0, type;
	struct GRDBLEND_INFO *B = NULL;
	char line[BUFSIZ], r_in[GMT_LONG_TEXT], *sense[2] = {"normal", "inverse"};
	void decode_R (char *string, double *w, double *e, double *s, double *n);

	GMT_set_meminc (GMT_SMALL_CHUNK);
	while (GMT_fgets (line, BUFSIZ, fp)) {	/* Read each input file for grid information */
		GMT_chop (line);
		if (line[0] == '#' || line[0] == '\0') continue;	/* Skip comment lines or blank lines */

		if (n == n_alloc) n_alloc = GMT_alloc_memory ((void **)&B, n, n_alloc, sizeof (struct GRDBLEND_INFO), GMT_program);
		memset ((void *)&B[n], 0, sizeof (struct GRDBLEND_INFO));	/* Initialize memory */
		nr = sscanf (line, "%s %s %lf", B[n].file, r_in, &B[n].weight);
		if (nr != 3) {
			fprintf (stderr, "%s: Read error for blending parameters near row %ld\n", GMT_program, n);
			exit (EXIT_FAILURE);
		}
		GMT_err_fail (GMT_read_grd_info (B[n].file, &B[n].G.header), B[n].file);	/* Read header structure */
		if (found_unsupported_format (&B[n].G.header, "Error for input grid", B[n].file)) exit (EXIT_FAILURE);
		if (!strcmp (r_in, "-")) {	/* Set inner = outer region */
			B[n].w_in = B[n].G.header.x_min;	B[n].e_in = B[n].G.header.x_max;
			B[n].s_in = B[n].G.header.y_min;	B[n].n_in = B[n].G.header.y_max;
		}
		else	/* Must decode the -R string */
			decode_R (&r_in[2], &B[n].w_in, &B[n].e_in, &B[n].s_in, &B[n].n_in);	/* Decode inner region */

		/* Skip the file if its outer region does not lie within the final grid region */
		if (h->x_min > B[n].e_in || h->x_max < B[n].w_in || h->y_min > B[n].n_in || h->y_max < B[n].s_in) {
			fprintf (stderr, "%s: Warning: File %s entirely outside final grid region (skipped)\n", GMT_program, B[n].file);
			continue;
		}

		/* Various sanity checking - e.g., all grid of same registration type and grid spacing */

		if (n == 0) {
			h->node_offset = B[n].G.header.node_offset;
			one_or_zero = !h->node_offset;
			GMT_RI_prepare (h);	/* Ensure -R -I consistency and set nx, ny */
		}
		if (h->node_offset != B[n].G.header.node_offset){
			fprintf (stderr, "%s: File %s has different registration than the first file given\n", GMT_program, B[n].file);
			exit (EXIT_FAILURE);
		}
		if (!GMT_IS_ZERO (B[n].G.header.x_inc - h->x_inc)){
			fprintf (stderr, "%s: File %s has different x-increment (%g) than the chosen output increment (%g)\n", GMT_program, B[n].file, B[n].G.header.x_inc, h->x_inc);
			exit (EXIT_FAILURE);
		}
		if (!GMT_IS_ZERO (B[n].G.header.y_inc - h->y_inc)){
			fprintf (stderr, "%s: File %s has different y-increment (%g) than the chosen output increment (%g)\n", GMT_program, B[n].file, B[n].G.header.y_inc, h->y_inc);
			exit (EXIT_FAILURE);
		}
		if (B[n].weight < 0.0) {	/* Negative weight means invert sense of taper */
			B[n].weight = fabs (B[n].weight);
			B[n].invert = TRUE;
		}

		GMT_err_fail (GMT_open_grd (B[n].file, &B[n].G, 'r'), B[n].file);	/* Open the grid for incremental row reading */

		/* Here, i0, j0 is the very first col, row to read, while i1, j1 is the very last col, row to read .
		 * Weights at the outside i,j should be 0, and reach 1 at the edge of the inside block */

		/* The following works for both pixel and grid-registered grids since we are here using the i,j to measure the width of the
		 * taper zone in units of dx, dy. */
		 
		B[n].out_i0 = irint ((B[n].G.header.x_min - h->x_min) / h->x_inc);
		B[n].in_i0  = irint ((B[n].w_in - h->x_min) / h->x_inc) - 1;
		B[n].in_i1  = irint ((B[n].e_in - h->x_min) / h->x_inc) + one_or_zero;
		B[n].out_i1 = irint ((B[n].G.header.x_max - h->x_min) / h->x_inc) - B[n].G.header.node_offset;
		B[n].out_j0 = irint ((h->y_max - B[n].G.header.y_max) / h->y_inc);
		B[n].in_j0  = irint ((h->y_max - B[n].n_in) / h->y_inc) - 1;
		B[n].in_j1  = irint ((h->y_max - B[n].s_in) / h->y_inc) + one_or_zero;
		B[n].out_j1 = irint ((h->y_max - B[n].G.header.y_min) / h->y_inc) - B[n].G.header.node_offset;

		B[n].wxl = M_PI * h->x_inc / (B[n].w_in - B[n].G.header.x_min);
		B[n].wxr = M_PI * h->x_inc / (B[n].G.header.x_max - B[n].e_in);
		B[n].wyu = M_PI * h->y_inc / (B[n].G.header.y_max - B[n].n_in);
		B[n].wyd = M_PI * h->y_inc / (B[n].s_in - B[n].G.header.y_min);

		if (B[n].out_j0 < 0) {	/* Must skip to first row inside the present -R */
			type = GMT_grdformats[B[n].G.header.type][0];
			if (type == 'c')	/* Old-style, 1-D netcdf grid */
				B[n].offset = B[n].G.header.nx * GMT_abs (B[n].out_j0);
			else if (type == 'n')	/* New, 2-D netcdf grid */
				B[n].offset = B[n].out_j0;
			else
				B[n].skip = (long)(B[n].G.n_byte * GMT_abs (B[n].out_j0));	/* do the fseek when we are ready to read first row */
		}

		/* Allocate space for one entire row */

		B[n].z = (float *) GMT_memory (VNULL, (size_t)(B[n].G.header.nx), sizeof (float), GMT_program);
		if (gmtdefs.verbose) fprintf (stderr, "%s: Blend file %s in %g/%g/%g/%g with %s weight %g [%ld-%ld]\n",
			GMT_program, B[n].G.header.name, B[n].w_in, B[n].e_in, B[n].s_in, B[n].n_in, sense[B[n].invert], B[n].weight, B[n].out_j0, B[n].out_j1);

		GMT_close_grd (&B[n].G);

		n++;
	}

	n_alloc = GMT_alloc_memory ((void **)&B, 0, n, sizeof (struct GRDBLEND_INFO), GMT_program);
	*blend = B;
	GMT_reset_meminc ();

	return (n);
}

void sync_input_rows (GMT_LONG row, struct GRDBLEND_INFO *B, GMT_LONG n_blend, double half) {
	GMT_LONG k;

	for (k = 0; k < n_blend; k++) {	/* Get every input grid ready for the new row */
		if (row < B[k].out_j0 || row > B[k].out_j1) {	/* Either done with grid or haven't gotten to this range yet */
			B[k].outside = TRUE;
			if (B[k].open) {
				GMT_close_grd (&B[k].G);	/* Done with this file */
				B[k].open = FALSE;
				GMT_free ((void *)B[k].z);
			}
			continue;
		}
		B[k].outside = FALSE;
		if (row <= B[k].in_j0)		/* Top cosine taper weight */
			B[k].wt_y = 0.5 * (1.0 - cos ((row - B[k].out_j0 + half) * B[k].wyu));
		else if (row >= B[k].in_j1)	/* Bottom cosine taper weight */
			B[k].wt_y = 0.5 * (1.0 - cos ((B[k].out_j1 - row + half) * B[k].wyd));
		else				/* We are inside the inner region; y-weight = 1 */
			B[k].wt_y = 1.0;
		B[k].wt_y *= B[k].weight;

		if (!B[k].open) {
			GMT_err_fail (GMT_open_grd (B[k].file, &B[k].G, 'r'), B[k].file);	/* Open the grid for incremental row reading */
			if (B[k].skip) GMT_fseek (B[k].G.fp, B[k].skip, SEEK_CUR);	/* Position for native binary files */
			B[k].G.start[0] += B[k].offset;					/* Start position for netCDF files */
			B[k].open = TRUE;
		}

		GMT_err_fail (GMT_read_grd_row (&B[k].G, 0, B[k].z), B[k].file);	/* Get one row from this file */
	}
}

void decode_R (char *string, double *w, double *e, double *s, double *n) {
	double *p[4];
	GMT_LONG i, pos, error = 0;
	char text[BUFSIZ];

	/* Needed to decode the inner region -Rw/e/s/n string */

	p[0] = w;	p[1] = e;	p[2] = s;	p[3] = n;

	i = pos = 0;
	while (!error && (GMT_strtok (string, "/", &pos, text))) {
		error += GMT_verify_expectations (GMT_io.in_col_type[i/2], GMT_scanf_arg (text, GMT_io.in_col_type[i/2], p[i]), text);
		i++;
	}
	if (error || (i != 4) || GMT_check_region (*p[0], *p[1], *p[2], *p[3])) {
		GMT_syntax ('R');
	}
}

void *New_grdblend_Ctrl () {	/* Allocate and initialize a new control structure */
	struct GRDBLEND_CTRL *C;
	
	C = (struct GRDBLEND_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct GRDBLEND_CTRL), "New_grdblend_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	
	C->N.nodata = GMT_d_NaN;
	C->Z.scale = 1.0;
	
	return ((void *)C);
}

void Free_grdblend_Ctrl (struct GRDBLEND_CTRL *C) {	/* Deallocate control structure */
	if (C->G.file) free ((void *)C->G.file);	
	GMT_free ((void *)C);	
}
