/*--------------------------------------------------------------------
 *    $Id: minmax.c,v 1.75 2011/07/08 21:27:06 guru Exp $
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
 * minmax.c will read ascii or binary tables and report the
 * extreme values for all columns
 *
 * Author:	Paul Wessel
 * Date:	20-JUN-2000
 * Revised:	20-FEB-2001 BCH-J: Added -D option for dateline discontinuity
 *		13-MAR-2006 -D made obsolete by improved checking of range
 * Version:	4.1
 */

#include "gmt.h"

struct MINMAX_CTRL {	/* All control options for this program (except common args) */
	/* active is TRUE if the option has been activated */
	struct C {	/* -C */
		GMT_LONG active;
	} C;
	struct E {	/* -E<L|l|H|h><col> */
		GMT_LONG active;
		GMT_LONG abs;
		GMT_LONG mode;
		GMT_LONG col;
	} E;
	struct I {	/* -Idx[/dy[/<dz>..] */
		GMT_LONG active;
		GMT_LONG ncol;
		double inc[GMT_MAX_COLUMNS];
	} I;
	struct S {	/* -S[x|y] */
		GMT_LONG active;
		GMT_LONG xbar, ybar;
	} S;
	struct T {	/* -T<dz>[/<col>] */
		GMT_LONG active;
		double inc;
		GMT_LONG col;
	} T;
};

int main (int argc, char **argv)
{
	GMT_LONG  error = FALSE, nofile = TRUE, done = FALSE, got_stuff = FALSE, first, give_r_string = FALSE;
	GMT_LONG  dateline = FALSE, brackets = FALSE, work_on_abs_value, special = FALSE, quad[4] = {FALSE, FALSE, FALSE, FALSE};

	char line[BUFSIZ], file[BUFSIZ], chosen[BUFSIZ], buffer[BUFSIZ], *not_used = NULL;

	GMT_LONG i, j, ncol, n_files = 0, fno, n_args, n_read, n_fields, n_expected_fields, quad_no, n;
	
	double *xyzmin = NULL, *xyzmax = NULL, west, east, south, north, low, high, value, e_min = DBL_MAX, e_max = -DBL_MAX, *in = NULL;
	double xmin1 = 360.0, xmin2 = 360.0, xmax1 = -360.0, xmax2 = -360.0;

	FILE *fp = NULL;

	struct MINMAX_CTRL *Ctrl = NULL;

	GMT_LONG strip_blanks_and_output (double x, GMT_LONG col);
	void *New_minmax_Ctrl (), Free_minmax_Ctrl (struct MINMAX_CTRL *C);

	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct MINMAX_CTRL *) New_minmax_Ctrl ();		/* Allocate and initialize defaults in a new control structure */
	
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				/* Common parameters */
                      
				case 'H':
				case 'M':
				case ':':
  				case 'b':
  				case 'f':
				case 'm':
				case '\0':
					error += GMT_parse_common_options (argv[i], 0, 0, 0, 0);
					break;

				/* Supplemental parameters */

				case 'C':
					Ctrl->C.active = TRUE;
					break;
                                case 'D':	/* Backwards compatible */
                                        dateline = TRUE;
                                        break;
                                case 'E':
                                        Ctrl->E.active = TRUE;
					switch (argv[i][2]) {
						case 'L':
							Ctrl->E.abs = TRUE;
						case 'l':
							Ctrl->E.mode = -1;
							break;
						case 'H':
							Ctrl->E.abs = TRUE;
						case 'h':
							Ctrl->E.mode = +1;
							break;
						default:
							error ++;
							fprintf (stderr, "%s: GMT SYNTAX ERROR -E. Flags are L|l|H|h\n", GMT_program);
							break;
					}
					Ctrl->E.col = atoi (&argv[i][3]);
                                        break;
				case 'I':
					Ctrl->I.active = TRUE;
					if (argv[i][2] == 'p') special = TRUE;
					j = (special) ? 3 : 2;
					Ctrl->I.ncol = GMT_getincn (&argv[i][j], Ctrl->I.inc, GMT_MAX_COLUMNS);
					break;
				case 'L':	/* Obsolete, but backward compatibility prevails [use -f instead] */
					GMT_io.in_col_type[0] = GMT_io.out_col_type[0] = GMT_IS_LON;
					GMT_io.in_col_type[1] = GMT_io.out_col_type[1] = GMT_IS_LAT;
					fprintf (stderr, "%s: Option -L is obsolete (but is processed correctly).  Please use -fg instead\n", GMT_program);
					break;
				case 'T':
					Ctrl->T.active = TRUE;
					j = sscanf (&argv[i][2], "%lf/%" GMT_LL "d", &Ctrl->T.inc, &Ctrl->T.col);
					if (j == 1) Ctrl->T.col = 0;
					break;
				case 'S':
					Ctrl->S.active = TRUE;
					j = 2;
					while (argv[i][j]) {
						if (argv[i][j] == 'x') Ctrl->S.xbar = TRUE;
						if (argv[i][j] == 'y') Ctrl->S.ybar = TRUE;
						j++;
					}
					if (j == 2) Ctrl->S.xbar = Ctrl->S.ybar = TRUE;
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

	if (error || GMT_give_synopsis_and_exit) {	/* Because it's ok to give no arguments */
		fprintf (stderr, "minmax %s - Find extreme values in ASCII tables\n\n", GMT_VERSION);
		fprintf (stderr, "usage: minmax [files] [-C] [-E<L|l|H|h><col>] [%s] [-I[p]<dx>[/<dy>[/<dz>..]]\n", GMT_H_OPT);
		fprintf (stderr, "\t[-S[x][y]] [-T<dz>[/<col>]] [%s] [%s] [%s] [%s]\n", GMT_t_OPT, GMT_bi_OPT, GMT_f_OPT, GMT_m_OPT);
             
              	if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);
 
		fprintf (stderr, "\t-C formats the min and max into separate columns\n");
		fprintf (stderr, "\t-E Return the record with extreme value in specified column <col>\n");
		fprintf (stderr, "\t   Specify l or h for min or max value, respectively.  Upper case L or H\n");
		fprintf (stderr, "\t   means we operate instead on the absolute values of the data.\n");
		GMT_explain_option ('H');
		fprintf (stderr, "\t-I returns textstring -Rw/e/s/n to nearest multiple of dx/dy (assumes 2+ col data)\n");
		fprintf (stderr, "\t   If -C is set then no -R string is issued.  Instead, the number of increments\n");
		fprintf (stderr, "\t   given determines how many columns are rounded off to the nearest multiple.\n");
		fprintf (stderr, "\t   If only one increment is given we also use it for the second column (for backwards compatibility).\n");
		fprintf (stderr, "\t   To override this behaviour, use -Ip<dx>.\n");
		fprintf (stderr, "\t-S adds extra space for error bars. Useful together with -I.\n");
		fprintf (stderr, "\t   -Sx leaves space for horizontal error bar using value in third (2) column.\n");
		fprintf (stderr, "\t   -Sy leaves space for vertical error bar using value in third (2) column.\n");
		fprintf (stderr, "\t   -S or -Sxy leaves space for both error bars using values in third&fourth (2&3) columns.\n");
		fprintf (stderr, "\t-T returns textstring -Tzmin/zmax/dz to nearest multiple of the given dz\n");
		fprintf (stderr, "\t   Calculations are based on the first (0) column only.  Append /<col> to use another column\n");
		GMT_explain_option (':');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf (stderr, "\t   Default is 2 input columns\n");
		GMT_explain_option ('f');
		GMT_explain_option ('m');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}
	
	GMT_check_lattice (&Ctrl->I.inc[0], &Ctrl->I.inc[1], NULL, &Ctrl->I.active);
	if (Ctrl->I.active && !special && Ctrl->I.ncol == 1) {		/* Special case of dy = dx if not given */
		Ctrl->I.inc[1] = Ctrl->I.inc[0];
		Ctrl->I.ncol = 2;
	}
	if (Ctrl->I.active && !Ctrl->C.active && Ctrl->I.ncol < 2) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  -Ip requires -C\n", GMT_program);
		error++;
	}
	if (Ctrl->I.active && Ctrl->T.active) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Only one of -I and -T can be specified\n", GMT_program);
		error++;
	}
	if (Ctrl->T.active && Ctrl->T.inc <= 0.0 ) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -T option.  Must specify a positive increment\n", GMT_program);
		error++;
	}
	if (Ctrl->I.active) {
		for (i = 0; i < Ctrl->I.ncol; i++) {
			if (Ctrl->I.inc[i] <= 0.0) {
				fprintf (stderr, "%s: GMT SYNTAX ERROR -I option.  Must specify positive increment(s)\n", GMT_program);
				error++;
			}
		}
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = 2;
	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] < 1) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data (-bi) must have at least 1 column\n", GMT_program);
		error++;
	}
	if (GMT_io.in_col_type[0] != GMT_IS_LON && dateline) {
	  	fprintf (stderr, "%s: GMT SYNTAX ERROR -D option: requires -fg\n", GMT_program);
		error++;
	}
	if (Ctrl->E.active && Ctrl->E.col < 0) {
	  	fprintf (stderr, "%s: GMT SYNTAX ERROR -E option: requires a positive column\n", GMT_program);
		error++;
	}
	if (error) exit (EXIT_FAILURE);

	if (GMT_io.binary[GMT_IN] && gmtdefs.verbose) {
		char *type[2] = {"double", "single"};
		fprintf (stderr, "%s: Expects %ld-column %s-precision binary data\n", GMT_program, GMT_io.ncol[GMT_IN], type[GMT_io.single_precision[GMT_IN]]);
	}

	if (n_files > 0)
		nofile = FALSE;
	else
		n_files = 1;

	n_args = (argc > 1) ? argc : 2;
	west = south = DBL_MAX;	east = north = -DBL_MAX;

	xyzmin = (double *) GMT_memory (VNULL, (size_t)1, sizeof (double), GMT_program);
	xyzmax = (double *) GMT_memory (VNULL, (size_t)1, sizeof (double), GMT_program);

	give_r_string = (Ctrl->I.active && !Ctrl->C.active);
	brackets = !Ctrl->C.active;
	work_on_abs_value = (Ctrl->E.active && Ctrl->E.abs);
	if (GMT_io.in_col_type[0] == GMT_IS_LON) {	/* Must check that output format won't mess things up by printing west > east */
		if (!strcmp (gmtdefs.output_degree_format, "D")) {
			if (dateline)
				strcpy (gmtdefs.output_degree_format, "-D");
			else
				strcpy (gmtdefs.output_degree_format, "+D");
			GMT_err_fail (GMT_geo_C_format (gmtdefs.output_degree_format, &GMT_io.geo), "");
			if (gmtdefs.verbose) fprintf (stderr, "%s: Warning: OUTPUT_DEGREE_FORMAT reset from D to %s to ensure east > west\n", GMT_program, gmtdefs.output_degree_format);
		}
		else if (!strcmp (gmtdefs.output_degree_format, "ddd:mm:ss")) {
			strcpy (gmtdefs.output_degree_format, "ddd:mm:ssF");
			GMT_err_fail (GMT_geo_C_format (gmtdefs.output_degree_format, &GMT_io.geo), "");
			if (gmtdefs.verbose) fprintf (stderr, "%s: Warning: OUTPUT_DEGREE_FORMAT reset from ddd:mm:ss to %s to ensure east > west\n", GMT_program, gmtdefs.output_degree_format);
		}
	}

	n_expected_fields = (GMT_io.ncol[GMT_IN]) ? GMT_io.ncol[GMT_IN] : GMT_MAX_COLUMNS;

	for (fno = 1; !done && fno < n_args; fno++) {     /* Loop over input files, if any */
		if (!nofile && argv[fno][0] == '-') continue;


		if (nofile) {   /* Just read standard input */
			fp = GMT_stdin;
			done = TRUE;
			strcpy (file, "<stdin>");
			if (gmtdefs.verbose) fprintf (stderr, "%s: Reading from standard input\n", GMT_program);
#ifdef SET_IO_MODE
			GMT_setmode (GMT_IN);
#endif
		}
		else {
			strcpy (file, argv[fno]);
			if ((fp = GMT_fopen (file, GMT_io.r_mode)) == NULL) {
				fprintf (stderr, "%s: Cannot open file %s\n", GMT_program, file);
				continue;
			}
		}

		if (GMT_io.io_header[GMT_IN]) for (i = 0; i < GMT_io.n_header_recs; i++) not_used = GMT_fgets (line, BUFSIZ, fp);

		n = ncol = n_read = 0;
		n_expected_fields = (GMT_io.ncol[GMT_IN]) ? GMT_io.ncol[GMT_IN] : GMT_MAX_COLUMNS;
		first = TRUE;

		n_fields = GMT_input (fp, &n_expected_fields, &in);

		while (! (GMT_io.status & GMT_IO_EOF)) {	/* Not yet EOF */

			n_read++;

			while ((GMT_io.status & GMT_IO_SEGMENT_HEADER) && !(GMT_io.status & GMT_IO_EOF)) {
					n_fields = GMT_input (fp, &n_expected_fields, &in);
			}
			if ((GMT_io.status & GMT_IO_EOF)) continue;	/* At EOF */

			if (GMT_io.status & GMT_IO_MISMATCH) {
				fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld (skipped)\n", GMT_program, n_fields, n_expected_fields, n_read);
				continue;
			}

			if (first) {	/* First time, allocate # of columns */

				ncol = n_expected_fields;
				if (Ctrl->S.active && 2 + Ctrl->S.xbar + Ctrl->S.ybar > ncol) {
					fprintf (stderr, "%s: Not enough columns to support the -S option\n", GMT_program);
					exit (EXIT_FAILURE);
				}
				if (Ctrl->E.active && Ctrl->E.col >= ncol) {
	  				fprintf (stderr, "%s: GMT SYNTAX ERROR -E option: Chosen column exceeds column range (0-%ld)\n", GMT_program, ncol-1);
					exit (EXIT_FAILURE);
				}
				if (Ctrl->T.active && Ctrl->T.col >= ncol) {
					fprintf (stderr, "%s: GMT SYNTAX ERROR -T option: Chosen column exceeds column range (0-%ld)\n", GMT_program, ncol-1);
					exit (EXIT_FAILURE);
				}
				if (Ctrl->T.active) ncol = Ctrl->T.col+1;
				if (give_r_string) ncol = 2;

				/* Now we know # of columns, so allocate memory */

				xyzmin = (double *) GMT_memory ((void *)xyzmin, (size_t)ncol, sizeof (double), GMT_program);
				xyzmax = (double *) GMT_memory ((void *)xyzmax, (size_t)ncol, sizeof (double), GMT_program);

				for (i = 0; i < ncol; i++) {	/* Initialize */
					xyzmin[i] = +DBL_MAX;
					xyzmax[i] = -DBL_MAX;
				}
				xmin1 = xmin2 = 360.0;
				xmax1 = xmax2 = -360.0;
				if (Ctrl->I.active && ncol < 2 && !Ctrl->C.active) Ctrl->I.active = FALSE;
				first = FALSE;
			}

			/* Decode all fields and update minmax arrays */

			while (dateline && in[GMT_X] > 180.0) in[GMT_X] -= 360.0;

			if (Ctrl->E.active && !GMT_is_dnan (in[Ctrl->E.col])) {
				value = (work_on_abs_value) ? fabs (in[Ctrl->E.col]) : in[Ctrl->E.col];
				if (Ctrl->E.mode == -1 && value < e_min) {
					e_min = value;
					strcpy (chosen, GMT_io.current_record);
				}
				else if (Ctrl->E.mode == +1 && value > e_max) {
					e_max = value;
					strcpy (chosen, GMT_io.current_record);
				}
			}
			else if (!Ctrl->E.active) {
				for (i = 0; i < ncol; i++) {
					if (GMT_is_dnan (in[i])) continue;
					if (i == 0 && GMT_io.in_col_type[i] == GMT_IS_LON) {
						while (in[i] < 0.0) in[i] += 360.0;	/* Start off with everything in 0-360 range */
						xmin1 = MIN (in[i], xmin1);
						xmax1 = MAX (in[i], xmax1);
						quad_no = (int)floor (in[i]/90.0);	/* Yields quadrants 0-3 */
						if (quad_no == 4) quad_no = 0;		/* When in[i] == 360.0 */
						quad[quad_no] = TRUE;
						while (in[i] > 180.0) in[i] -= 360.0;	/* Switch to -180+/180 range */
						xmin2 = MIN (in[i], xmin2);
						xmax2 = MAX (in[i], xmax2);
					}
					else if ((i == 0 && Ctrl->S.xbar) || (i == 1 && Ctrl->S.ybar)) {
						/* Add/subtract value from error bar column */
						j = (i == 1 && Ctrl->S.xbar) ? 3 : 2;
						value = fabs(in[j]);
						if (in[i] - value < xyzmin[i]) xyzmin[i] = in[i] - value;
						if (in[i] + value > xyzmax[i]) xyzmax[i] = in[i] + value;
					}
					else {
						if (in[i] < xyzmin[i]) xyzmin[i] = in[i];
						if (in[i] > xyzmax[i]) xyzmax[i] = in[i];
					}
				}
			}

			n++;

			n_fields = GMT_input (fp, &n_expected_fields, &in);
		}
		if (fp != GMT_stdin) GMT_fclose (fp);

		if (!got_stuff) got_stuff = (n > 0);	/* We were able to open and read at least 1 record from a file */

		if (GMT_io.in_col_type[0] == GMT_IS_LON) {	/* Now determine longitude range */
			int n_quad;
			n_quad = (int)(quad[0] + quad[1] + quad[2] + quad[3]);	/* How many quadrants had data */
			if (quad[0] && quad[3]) {	/* Longitudes on either side of Greenwich only, must use -180/+180 notation */
				xyzmin[0] = xmin2;
				xyzmax[0] = xmax2;
				GMT_io.geo.range = 2;	/* Override this setting explicitly */
			}
			else if (quad[1] && quad[2]) {	/* Longitudes on either side of the date line, must user 0/360 notation */
				xyzmin[0] = xmin1;
				xyzmax[0] = xmax1;
				GMT_io.geo.range = 0;	/* Override this setting explicitly */
			}
			else if (n_quad == 2 && ((quad[0] && quad[2]) || (quad[1] && quad[3]))) {	/* Funny quadrant gap, pick shortest longitude extent */
				if ((xmax1 - xmin1) < (xmax2 - xmin2)) {	/* 0/360 more compact */
					xyzmin[0] = xmin1;
					xyzmax[0] = xmax1;
					GMT_io.geo.range = 0;	/* Override this setting explicitly */
				}
				else {						/* -180/+180 more compact */
					xyzmin[0] = xmin2;
					xyzmax[0] = xmax2;
					GMT_io.geo.range = 2;	/* Override this setting explicitly */
				}
			}
			else {						/* Either will do, use default settings */
				xyzmin[0] = xmin1;
				xyzmax[0] = xmax1;
			}
			if (xyzmin[0] > xyzmax[0]) xyzmin[0] -= 360.0;
			if (xyzmin[0] < 0.0 && xyzmax[0] < 0.0) xyzmin[0] += 360.0, xyzmax[0] += 360.0;
		}
		if (give_r_string && got_stuff) {
			west  = MIN (west, xyzmin[0]);
			east  = MAX (east, xyzmax[0]);
			south = MIN (south, xyzmin[1]);
			north = MAX (north, xyzmax[1]);
		}
		else if (Ctrl->T.active && got_stuff) {
			west  = MIN (west, xyzmin[Ctrl->T.col]);
			east  = MAX (east, xyzmax[Ctrl->T.col]);
		}
		else if (!Ctrl->E.active && n > 0) {
			if (!Ctrl->C.active) {sprintf (buffer, "%s: N = %ld\t", file, n);	GMT_fputs (buffer, GMT_stdout);}
			for (i = 0; i < ncol; i++) {
				if (xyzmin[i] == DBL_MAX) {
					low = high = GMT_d_NaN;
				}
				else if (i < Ctrl->I.ncol) {	/* Special treatment for x and y if selected */
					low  = (Ctrl->I.active) ? floor (xyzmin[i] / Ctrl->I.inc[i]) * Ctrl->I.inc[i] : xyzmin[i];
					high = (Ctrl->I.active) ? ceil (xyzmax[i] / Ctrl->I.inc[i]) * Ctrl->I.inc[i] : xyzmax[i];
				}
				else {
					low = xyzmin[i];
					high = xyzmax[i];
				}
				if (brackets) GMT_fputs ("<", GMT_stdout);
				GMT_ascii_output_one (GMT_stdout, low, i);
				(Ctrl->C.active) ? GMT_fputs ("\t", GMT_stdout) : GMT_fputs ("/", GMT_stdout);
				GMT_ascii_output_one (GMT_stdout, high, i);
				if (brackets) GMT_fputs (">", GMT_stdout);
				if (i < (ncol - 1)) GMT_fputs ("\t", GMT_stdout);
			}
			GMT_fputs ("\n", GMT_stdout);
		}
	}
	if (got_stuff) {
		if (give_r_string) {
			west  = floor (west / Ctrl->I.inc[0]) * Ctrl->I.inc[0];
			east  = ceil (east / Ctrl->I.inc[0]) * Ctrl->I.inc[0];
			south = floor (south / Ctrl->I.inc[1]) * Ctrl->I.inc[1];
			north = ceil (north / Ctrl->I.inc[1]) * Ctrl->I.inc[1];
			if (east < west) east += 360.0;
			if (GMT_IS_ZERO (east - 360.0) && west > 180.0) GMT_io.geo.range = 2;	/* Override this setting explicitly to get negative lons since 360 otherwise will become 0 */
			GMT_fputs ("-R", GMT_stdout);
			strip_blanks_and_output (west, 0);	GMT_fputs ("/", GMT_stdout);
			strip_blanks_and_output (east, 0);	GMT_fputs ("/", GMT_stdout);
			strip_blanks_and_output (south, 1);	GMT_fputs ("/", GMT_stdout);
			strip_blanks_and_output (north, 1);	GMT_fputs ("\n", GMT_stdout);
		}
		else if (Ctrl->T.active) {	/* -T option */
			west  = floor (west / Ctrl->T.inc) * Ctrl->T.inc;
			east  = ceil (east / Ctrl->T.inc) * Ctrl->T.inc;
			GMT_fputs ("-T", GMT_stdout);
			strip_blanks_and_output (west, 0);	GMT_fputs ("/", GMT_stdout);
			strip_blanks_and_output (east, 0);	GMT_fputs ("/", GMT_stdout);
			strip_blanks_and_output (Ctrl->T.inc, 0);	GMT_fputs ("\n", GMT_stdout);
		}
		else if (Ctrl->E.active)
			GMT_fputs (chosen, GMT_stdout);
	}
	else if (!got_stuff)
		fprintf (stderr, "%s: No input data found!\n", GMT_program);

	GMT_free ((void *)xyzmin);
	GMT_free ((void *)xyzmax);

	Free_minmax_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

GMT_LONG strip_blanks_and_output (double x, GMT_LONG col)
{
	/* Alternative to GMT_ascii_output_one that strips off leading blanks first */

	GMT_LONG k;
	char text[GMT_TEXT_LEN];

	GMT_ascii_format_one (text, x, GMT_io.out_col_type[col]);
	for (k = 0; text[k] && text[k] == ' '; k++);
	return (GMT_fputs (&text[k], GMT_stdout));
}

void *New_minmax_Ctrl () {	/* Allocate and initialize a new control structure */
	struct MINMAX_CTRL *C;
	
	C = (struct MINMAX_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct MINMAX_CTRL), "New_minmax_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	
	C->E.col = -1;
	
	return ((void *)C);
}

void Free_minmax_Ctrl (struct MINMAX_CTRL *C) {	/* Deallocate control structure */
	GMT_free ((void *)C);	
}
