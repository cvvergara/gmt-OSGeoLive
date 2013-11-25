/*--------------------------------------------------------------------
 *    $Id: blockmean.c 10014 2013-04-14 00:57:05Z pwessel $
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
   blockmean.c
   reads x, y, data, [weight] on GMT_stdin or file and writes out one value
   per cell, where cellular region is bounded by West East South North
   and cell dimensions are delta_x, delta_y.
   
   Latest method uses a hash table and linked lists if the region will
   not simply fit in core memory.
      
   Author:      Walter H. F. Smith
   Version:     3.0, testing hash tables
   Date:        4 October, 1988
   Modified:	26 April 1991 by WHFS for gmt v2.0
   Modified:	3 Jan 1995 by PW for gmt 3.0
   Modified:	3 May 1998 by PW for gmt 3.1
   Modified:	18 Oct 1999 by PW to add -S
   Modified:	3.3.5: 10 Jul 2000 by PW to add -L
   Version:	3.4: 01-MAR-2001 by PW, Use -F instead of -N, and add -C
   Version:	4: 01-AUG-2001 by PW, Added -f
   Version:	4.1: 14-SEP-2005 by PW, Added enhanced -I
   Version	4.1.2: 24-MAR-2006 by PW: No longer defines global variables
   Version	4.1.2: 4-APR-2006 by PW: Added -E for also returning standard
   		deviation, min, and max value per block.  Removed linked list
		approach in favor of full array, but only allocate space for
		items actually used.
			Also implemented size_t counters to be 64-bit compatible.
*/

#define BLOCKMEAN

#include "gmt.h"
#include "block_subs.h"

#define BLK_Z	0
#define BLK_W	1
#define BLK_S	0
#define BLK_L	1
#define BLK_H	2

struct BLK_PAIR {
	double a[2];
};

struct BLK_SLH {
	double a[3];
};

int main (int argc, char **argv)
{
	GMT_LONG	error = FALSE, nofile = TRUE, done = FALSE, first = TRUE, use_xy, duplicate_col;

	FILE *fp = NULL;

	double	weight, weighted_z, *in = NULL, half_dx, wesn[4], out[7], iw;

	GMT_LONG	i, j, n_expected_fields, n_fields, n_req, n_out, w_col;
	GMT_LONG	n_files = 0, fno, n_args, ij, n_blocks;
	GMT_LONG	n_cells_filled, n_read, n_lost, n_pitched, *np = NULL;
	
	char	modifier, line[BUFSIZ], format[BUFSIZ];
	
	struct GRD_HEADER h;

	struct	BLK_PAIR *xy = NULL, *zw = NULL;
	struct	BLK_SLH *slh = NULL;
	struct BLOCKMEAN_CTRL *Ctrl = NULL;
	
	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct BLOCKMEAN_CTRL *) New_blockmean_Ctrl ();	/* Allocate and initialize a new control structure */

	GMT_grd_init (&h, argc, argv, FALSE);
	
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
              
				/* Common parameters */
                      
				case 'H':
				case 'R':
				case 'V':
				case ':':
				case 'b':
				case 'f':
				case '\0':
					error += GMT_parse_common_options (argv[i], &h.x_min, &h.x_max, &h.y_min, &h.y_max);
					break;
                              
				/* Supplemental parameters */
                              
				case 'C':	/* Report center of block instead */
					Ctrl->C.active = TRUE;
					break;
				case 'E':	/* Extended report with standard deviation, min, and max in cols 4-6 */
					Ctrl->E.active = TRUE;
					break;
				case 'I':
					Ctrl->I.active = TRUE;
					if (GMT_getinc (&argv[i][2], &Ctrl->I.xinc, &Ctrl->I.yinc)) {
						GMT_inc_syntax ('I', 1);
						error = TRUE;
					}
					break;
				case 'L':	/* Obsolete, but backward compatibility prevails [use -f instead] */
					GMT_io.in_col_type[GMT_X] = GMT_io.out_col_type[GMT_X] = GMT_IS_LON;
					GMT_io.in_col_type[GMT_Y] = GMT_io.out_col_type[GMT_Y] = GMT_IS_LAT;
					fprintf (stderr, "%s: Option -L is obsolete (but is processed correctly).  Please use -f instead\n", GMT_program);
					break;
				case 'N':	/* Backward compatible with 3.3.6 */
				case 'F':
					Ctrl->F.active = TRUE;
					break;
				case 'S':
					Ctrl->S.active = TRUE;
					switch (argv[i][2]) {
						case 'w':	/* Report weight sum */
							Ctrl->S.mode = 2;
							break;
						default:	/* Report z sum */
							Ctrl->S.mode = 1;
							break;
					}
					break;
				case 'W':
					Ctrl->W.active = TRUE;
					if ( (modifier = argv[i][2]) == 'i' || modifier == 'I')
						Ctrl->W.weighted[GMT_IN] = TRUE;
					else if (modifier == 'O' || modifier == 'o')
						Ctrl->W.weighted[GMT_OUT] = TRUE;
					else
						Ctrl->W.weighted[GMT_IN] = Ctrl->W.weighted[GMT_OUT] = TRUE;
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
		fprintf (stderr, "blockmean %s - Block averaging by L2 norm\n\n", GMT_VERSION);
		fprintf (stderr, "usage: blockmean [infile(s)] %s %s\n", GMT_I_OPT, GMT_Rgeo_OPT);
		fprintf (stderr, "\t[-C] [-E] [-F] [%s] [-S[w|z]] [-V] [-W[i][o]] [%s] [%s]\n", GMT_H_OPT, GMT_t_OPT, GMT_b_OPT);
		fprintf (stderr, "\t[%s]\n\n", GMT_f_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		GMT_inc_syntax ('I', 0);
		GMT_explain_option ('R');
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-C Output center of block and mean z-value [Default outputs (mean x, mean y) location].\n");
		fprintf (stderr, "\t-E Extend output with st.dev (s), low (l), and high (h) value per block, i.e.,\n");
		fprintf (stderr, "\t   output (x,y,z,s,l,h[,w]) [Default outputs (x,y,z[,w]); see -W regarding w.\n");
		fprintf (stderr, "\t-F Offsets registration so block edges are on gridlines (pixel reg.) [Default: grid reg.].\n");
		GMT_explain_option ('H');
		fprintf (stderr, "\t-Sz report block sums rather than mean values [Default is mean values].\n");
		fprintf (stderr, "\t   -Sw reports weight sums instead of data sums.\n");
		GMT_explain_option ('V');
		fprintf (stderr, "\t-W sets Weight options.\n");
		fprintf (stderr, "\t   -Wi reads Weighted Input (4 cols: x,y,z,w) but writes only (x,y,z[,s,l,h]) Output.\n");
		fprintf (stderr, "\t   -Wo reads unWeighted Input (3 cols: x,y,z) but reports sum (x,y,z[,s,l,h],w) Output.\n");
		fprintf (stderr, "\t   -W with no modifier has both weighted Input and Output; Default is no weights used.\n");
		GMT_explain_option (':');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf (stderr, "\t   Default is 3 columns (or 4 if -W is set).\n");
		GMT_explain_option ('o');
		GMT_explain_option ('n');
		GMT_explain_option ('f');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	GMT_check_lattice (&Ctrl->I.xinc, &Ctrl->I.yinc, &Ctrl->F.active, &Ctrl->I.active);

	if (!project_info.region_supplied) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify -R option\n", GMT_program);
		error++;
	}
	if (Ctrl->I.xinc <= 0.0 || Ctrl->I.yinc <= 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -I option.  Must specify positive increment(s)\n", GMT_program);
		error++;
	}

	if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
	n_req = (Ctrl->W.weighted[GMT_IN]) ? 4 : 3;
	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = n_req;
	if (GMT_io.binary[GMT_IN] && n_req > GMT_io.ncol[GMT_IN]) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR.  binary input data must have at least %ld columns\n", GMT_program, n_req);
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
	
	h.x_inc = Ctrl->I.xinc;
	h.y_inc = Ctrl->I.yinc;
	h.node_offset = (int)Ctrl->F.active;
	
	GMT_RI_prepare (&h);	/* Ensure -R -I consistency and set nx, ny */

	duplicate_col = (GMT_360_RANGE (h.x_min, h.x_max) && h.node_offset == 0);	/* E.g., lon = 0 column should match lon = 360 column */
	half_dx = 0.5 * h.x_inc;
	n_blocks = GMT_get_nm (h.nx, h.ny);
	
	if ((zw = (struct BLK_PAIR *) GMT_memory (VNULL, n_blocks, sizeof (struct BLK_PAIR), GMT_program)) == NULL) {
		fprintf (stderr, "%s: ERROR: Unable to allocate memory for %ld blocks [%ld bytes for zw].\n",
			GMT_program, (GMT_LONG)n_blocks, (GMT_LONG)sizeof (struct BLK_PAIR));
		exit (EXIT_FAILURE);
	}
	if (!Ctrl->C.active && (xy = (struct BLK_PAIR *) GMT_memory (VNULL, n_blocks, sizeof (struct BLK_PAIR), GMT_program)) == NULL) {
		fprintf (stderr, "%s: ERROR: Unable to allocate memory for %ld blocks [%ld bytes for xy].\n",
			GMT_program, (GMT_LONG)n_blocks, (GMT_LONG)sizeof (struct BLK_PAIR));
		exit (EXIT_FAILURE);
	}
	if (Ctrl->E.active && (slh = (struct BLK_SLH *) GMT_memory (VNULL, n_blocks, sizeof (struct BLK_SLH), GMT_program)) == NULL) {
		fprintf (stderr, "%s: ERROR: Unable to allocate memory for %ld blocks [%ld bytes for slh].\n",
			GMT_program, (GMT_LONG)n_blocks, (GMT_LONG)sizeof (struct BLK_SLH));
		exit (EXIT_FAILURE);
	}
	if ((Ctrl->W.weighted[GMT_IN] && Ctrl->E.active) && (np = (GMT_LONG *) GMT_memory (VNULL, n_blocks, sizeof (GMT_LONG), GMT_program)) == NULL) {
		fprintf (stderr, "%s: ERROR: Unable to allocate memory for %ld blocks [%ld bytes for counter].\n",
			GMT_program, (GMT_LONG)n_blocks, (GMT_LONG)sizeof (GMT_LONG));
		exit (EXIT_FAILURE);
	}
	if (gmtdefs.verbose) {
		double mem;
		char unit = 'M';
		sprintf (format, "%%s: W: %s E: %s S: %s N: %s nx: %%ld ny: %%ld\n", gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format);
		fprintf (stderr, format, GMT_program, h.x_min, h.x_max, h.y_min, h.y_max, h.nx, h.ny);
		mem = (double)sizeof (struct BLK_PAIR);
		if (!Ctrl->C.active) mem += (double)sizeof (struct BLK_PAIR);
		if (Ctrl->E.active)  mem += (double)sizeof (struct BLK_SLH);
		if (Ctrl->W.weighted[GMT_IN] && Ctrl->E.active) mem += (double)sizeof (GMT_LONG);
		mem *= (double)n_blocks;
		mem /= (1024.0 * 1024.0);	/* Report Mbytes */
		if (mem > 1000.0) {		/* Report Gbytes */
			mem /= 1024.0;
			unit = 'G';
		}
		fprintf (stderr, "%s: Using a total of %.3g %cb for all arrays.\n", GMT_program, mem, unit);
	}

	use_xy = !Ctrl->C.active;	/* If not -C then we must keep track of x,y locations */
	GMT_set_xy_domain (wesn, &h);	/* May include some padding if gridline-registered */

	n_read = n_pitched = 0;

	n_expected_fields = (GMT_io.binary[GMT_IN]) ? GMT_io.ncol[GMT_IN] : 3 + Ctrl->W.weighted[GMT_IN];

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
				GMT_fgets (line, BUFSIZ, fp);
				GMT_chop (line);
				if (first && GMT_io.io_header[GMT_OUT]) {
					(Ctrl->W.weighted[GMT_OUT] && !(Ctrl->W.weighted[GMT_IN])) ? sprintf (format, "%s weights\n", line) : sprintf (format, "%s\n", line);
					GMT_fputs(format, GMT_stdout);
				}
			}
			first = FALSE;
		}

		while ((n_fields = GMT_input (fp, &n_expected_fields, &in)) >= 0 && !(GMT_io.status & GMT_IO_EOF)) {	/* Not yet EOF */

			n_read++;

			if (GMT_io.status & GMT_IO_MISMATCH) {
				fprintf (stderr, "%s: Mismatch between actual (%ld) and expected (%ld) fields near line %ld (skipped)\n", GMT_program, n_fields,  n_expected_fields, n_read);
				continue;
			}

			if (GMT_is_dnan (in[GMT_Z])) continue;	/* Skip when z = NaN */

			if (GMT_y_is_outside (in[GMT_Y],  wesn[2], wesn[3])) continue;	/* Outside y-range */
			if (GMT_x_is_outside (&in[GMT_X], wesn[0], wesn[1])) continue;	/* Outside x-range */
			if (duplicate_col && (wesn[1]-in[GMT_X] < half_dx)) {	/* Only compute mean values for the west column and not the repeating east column with lon += 360 */
				in[GMT_X] -= 360.0;	/* Make this point be considered for the western block mean value */
			}

			/* Get i and j indices of this block */
			
			i = GMT_x_to_i (in[GMT_X], h.x_min, h.x_inc, h.xy_off, h.nx);
			if ( i < 0 || i >= h.nx ) continue;
			j = GMT_y_to_j (in[GMT_Y], h.y_min, h.y_inc, h.xy_off, h.ny);
			if ( j < 0 || j >= h.ny ) continue;
			
			/* OK, this point is inside and will be used */
			
			weight = (Ctrl->W.weighted[GMT_IN]) ? in[3] : 1.0;
			weighted_z = in[GMT_Z] * weight;
			ij = GMT_IJ (j, i, h.nx);		/* 64-bit safe 1-D index */
			if (use_xy) {
				xy[ij].a[BLK_X] += (in[GMT_X]*weight);
				xy[ij].a[BLK_Y] += (in[GMT_Y]*weight);
			}
			if (Ctrl->E.active) {	/* Add up sum(w*z^2) and n for weighted stdev and keep track of min,max */
				slh[ij].a[BLK_S] += (weighted_z * in[GMT_Z]);
				if (Ctrl->W.weighted[GMT_IN]) np[ij]++;
				if (zw[ij].a[BLK_W] == 0) {	/* Initialize low,high */
					slh[ij].a[BLK_L] = +DBL_MAX;
					slh[ij].a[BLK_H] = -DBL_MAX;
				}
				if (in[GMT_Z] < slh[ij].a[BLK_L]) slh[ij].a[BLK_L] = in[GMT_Z];
				if (in[GMT_Z] > slh[ij].a[BLK_H]) slh[ij].a[BLK_H] = in[GMT_Z];
			}
			zw[ij].a[BLK_W] += weight;
			zw[ij].a[BLK_Z] += weighted_z;
			n_pitched++;
		}
		if (fp != GMT_stdin) GMT_fclose(fp);

	}

	n_out = ((Ctrl->W.weighted[GMT_OUT]) ? 4 : 3) + 3 * Ctrl->E.active;
	w_col = n_out - 1;
	n_cells_filled = 0;

	if (gmtdefs.verbose) fprintf (stderr, "%s: Calculating block means\n", GMT_program);

	for (ij = 0; ij < n_blocks; ij++) {

		if (zw[ij].a[BLK_W] == 0.0) continue;

		n_cells_filled++;
		if (Ctrl->W.weighted[GMT_OUT]) out[w_col] = zw[ij].a[BLK_W];
		iw = 1.0 / zw[ij].a[BLK_W];
		if (use_xy) {
			out[GMT_X] = xy[ij].a[BLK_X] * iw;
			out[GMT_Y] = xy[ij].a[BLK_Y] * iw;
		}
		else {	/* Use block center */
			i = ij % ((GMT_LONG)h.nx);
			j = ij / ((GMT_LONG)h.nx);
			out[GMT_X] = GMT_i_to_x (i, h.x_min, h.x_max, h.x_inc, h.xy_off, h.nx);
			out[GMT_Y] = GMT_j_to_y (j, h.y_min, h.y_max, h.y_inc, h.xy_off, h.ny);
		}
		if (Ctrl->S.mode) {
			out[GMT_Z] = (Ctrl->S.mode == 2) ? zw[ij].a[BLK_W] : zw[ij].a[BLK_Z];
		}
		else
			out[GMT_Z] = zw[ij].a[BLK_Z] * iw;
		if (Ctrl->E.active) {
			if (Ctrl->W.weighted[GMT_IN]) {
				out[3] = (np[ij] > 1) ? d_sqrt ((zw[ij].a[BLK_W] * slh[ij].a[BLK_S] - zw[ij].a[BLK_Z] * zw[ij].a[BLK_Z]) \
				/ (zw[ij].a[BLK_W] * zw[ij].a[BLK_W] * ((np[ij] - 1.0)/np[ij]))) : GMT_d_NaN;
			}
			else {
				out[3] = (zw[ij].a[BLK_W] > 1.0) ? d_sqrt ((zw[ij].a[BLK_W] * slh[ij].a[BLK_S] - zw[ij].a[BLK_Z] * zw[ij].a[BLK_Z]) \
				/ (zw[ij].a[BLK_W] * (zw[ij].a[BLK_W] - 1.0))) : GMT_d_NaN;
			}
			out[4] = slh[ij].a[BLK_L];
			out[5] = slh[ij].a[BLK_H];
		}
		GMT_output (GMT_stdout, n_out, out);
	}

	GMT_free ((void *)zw);
	if (use_xy) GMT_free ((void *)xy);
	if (Ctrl->E.active) {
		GMT_free ((void *)slh);
		if (Ctrl->W.weighted[GMT_IN]) GMT_free ((void *)np);
	}
	
	n_lost = n_read - n_pitched;
	if (gmtdefs.verbose) fprintf(stderr,"%s: N read: %ld N used: %ld N outside_area: %ld N cells filled: %ld\n",
		GMT_program, n_read, n_pitched, n_lost, n_cells_filled);

	Free_blockmean_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}
#include "block_subs.c"
