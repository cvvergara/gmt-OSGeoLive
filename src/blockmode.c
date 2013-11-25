/*--------------------------------------------------------------------
 *    $Id: blockmode.c 10014 2013-04-14 00:57:05Z pwessel $
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
   blockmode.c
   Takes lon, lat, data, [weight] on GMT_stdin or file and writes out one value
   per cell, where cellular region is bounded by West East South North and
   cell dimensions are delta_x, delta_y.
      
   Author: 	Walter H. F. Smith
   Date:	28 October, 1998
   Version:	1st version based on GMT3.1b blockmedian.c, but with changes
   		to struct BLK_DATA and to sort comparison routines to facilitate
   		mode coding.
   Modified:	3.3.5: 10-JUL-2000 PW: Added -L option
   Version:	3.4 01-MAR-2001 by PW, Replace -N with -F, and added -C
   Version:	4 01-MAR-2003 by PW
   Version:	4.1.x: 14-SEP-2005 by PW, Added enhanced -I
   		4-APR-2006 by PW: Added -E for LMS scale, low, and high value
			Also implemented size_t counters to be 64-bit compatible.
*/

#define BLOCKMODE

#include "gmt.h"
#include "block_subs.h"

int main (int argc, char **argv)
{
	GMT_LONG	error = FALSE, nofile = TRUE, done = FALSE, first = TRUE, mode_xy, duplicate_col;

	/* Default value for go_quickly = FALSE for backward compatibility with 3.0  */

	FILE *fp = NULL;

	double	*in = NULL, out[7], wesn[4], i_n_in_cell, weight, half_dx, *z_tmp = NULL;

	GMT_LONG i, j, ix, iy, fno, n_files = 0, n_args, n_req;
	GMT_LONG n_expected_fields, n_fields, n_out, w_col;
	GMT_LONG index, first_in_cell, first_in_new_cell, n_lost, n_read;
	GMT_LONG n_cells_filled, n_in_cell, nz;
	GMT_LONG n_alloc = 0, nz_alloc = 0, n_pitched;

	char	modifier, buffer[BUFSIZ], format[BUFSIZ];

	struct GRD_HEADER h;
	struct BLK_DATA *data = NULL;
	struct BLOCKMODE_CTRL *Ctrl = NULL;
	
	double weighted_mode (struct BLK_DATA *d, double wsum, GMT_LONG n, GMT_LONG k);
	
	argc = (int)GMT_begin (argc, argv);

	Ctrl = (struct BLOCKMODE_CTRL *) New_blockmode_Ctrl ();	/* Allocate and initialize a new control structure */

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
                              
				case 'C':
					Ctrl->C.active = TRUE;
					break;
				case 'E':
					Ctrl->E.active = TRUE;		/* Extended report with standard deviation, min, and max in cols 4-6 */
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
				case 'Q':
					Ctrl->Q.active = TRUE;
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
		fprintf (stderr, "blockmode %s - Block averaging by mode estimates\n\n", GMT_VERSION);
		fprintf (stderr, "usage: blockmode [infile(s)] %s %s\n", GMT_I_OPT, GMT_Rgeo_OPT);
		fprintf (stderr, "\t[-C] [-E] [-F] [%s] [-Q] [-V] [-W[i][o] ] [%s] [%s]\n", GMT_H_OPT, GMT_t_OPT, GMT_b_OPT);
		fprintf (stderr, "\t[%s]\n\n", GMT_f_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		GMT_inc_syntax ('I', 0);
		GMT_explain_option ('R');
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-C Output center of block and mode z-value [Default is mode location (but see -Q)].\n");
		fprintf (stderr, "\t-E Extend output with LMS scale (s), low (l), and high (h) value per block, i.e.,\n");
		fprintf (stderr, "\t   output (x,y,z,s,l,h[,w]) [Default outputs (x,y,z[,w]); see -W regarding w.\n");
		fprintf (stderr, "\t-F Offsets registration so block edges are on gridlines (pixel reg.) [Default: grid reg.].\n");
		GMT_explain_option ('H');
		fprintf (stderr, "\t-Q Quicker; get mode z and mean x,y [Default gets mode x, mode y, mode z].\n");
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

	if (Ctrl->C.active && Ctrl->Q.active) {
		fprintf (stderr, "%s: GMT WARNING:  -C overrides -Q\n", GMT_program);
		Ctrl->Q.active = FALSE;
	}
	if (!project_info.region_supplied) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify -R option\n", GMT_program);
		error++;
	}
	if (Ctrl->I.xinc <= 0.0 || Ctrl->I.yinc <= 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -I option.  Must specify positive increment(s)\n", GMT_program);
		error = TRUE;
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
	
	mode_xy = !Ctrl->C.active;
	h.x_inc = Ctrl->I.xinc;
	h.y_inc = Ctrl->I.yinc;
	h.node_offset = (int)Ctrl->F.active;
	GMT_RI_prepare (&h);	/* Ensure -R -I consistency and set nx, ny */
	h.xy_off = 0.5 * h.node_offset;	/* Use to calculate mean location of block */
	duplicate_col = (GMT_360_RANGE (h.x_min, h.x_max) && h.node_offset == 0);	/* E.g., lon = 0 column should match lon = 360 column */
	half_dx = 0.5 * h.x_inc;

	if (gmtdefs.verbose) {
		sprintf (format, "%%s: W: %s E: %s S: %s N: %s nx: %%ld ny: %%ld\n", gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format);
		fprintf (stderr, format, GMT_program, h.x_min, h.x_max, h.y_min, h.y_max, h.nx, h.ny);
	}

	n_read = n_pitched = 0;

	GMT_set_xy_domain (wesn, &h);	/* May include some padding if gridline-registered */
		
	/* Read the input data  */

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
				GMT_fgets (buffer, BUFSIZ, fp);
				GMT_chop (buffer);
				if (first && GMT_io.io_header[GMT_OUT]) {
					(Ctrl->W.weighted[GMT_OUT] && !(Ctrl->W.weighted[GMT_IN])) ? sprintf (format, "%s weights\n", buffer) : sprintf (format, "%s\n", buffer);
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

			ix = GMT_x_to_i (in[GMT_X], h.x_min, h.x_inc, h.xy_off, h.nx);
			if ( ix < 0 || ix >= h.nx ) continue;
			iy = GMT_y_to_j (in[GMT_Y], h.y_min, h.y_inc, h.xy_off, h.ny);
			if ( iy < 0 || iy >= h.ny ) continue;

			index = GMT_IJ (iy, ix, h.nx);		/* 64-bit safe 1-D index */

			if (n_pitched == n_alloc) n_alloc = GMT_alloc_memory ((void **)&data, n_pitched, n_alloc, sizeof (struct BLK_DATA), GMT_program);
			data[n_pitched].i = index;
			if (mode_xy) {
				data[n_pitched].a[BLK_X] = in[GMT_X];
				data[n_pitched].a[BLK_Y] = in[GMT_Y];
			}
			data[n_pitched].a[BLK_Z] = in[GMT_Z];
			data[n_pitched].a[BLK_W] = (Ctrl->W.weighted[GMT_IN]) ? in[3] : 1.0;

			n_pitched++;
		}
		if (fp != GMT_stdin) GMT_fclose(fp);

	}
	n_alloc = GMT_alloc_memory ((void **)&data, 0, n_pitched, sizeof (struct BLK_DATA), GMT_program);

	if (n_read == 0) {	/* Blank/empty input files */
		if (gmtdefs.verbose) fprintf (stderr, "%s: No data records found; no output produced", GMT_program);
		exit (EXIT_SUCCESS);
	}
	if (n_pitched == 0) {	/* No points inside region */
		if (gmtdefs.verbose) fprintf (stderr, "%s: No data points found inside the region; no output produced", GMT_program);
		exit (EXIT_SUCCESS);
	}
	n_lost = n_read - n_pitched;
	
	if (gmtdefs.verbose) fprintf(stderr,"%s: N read: %ld N used: %ld N outside_area: %ld\n", GMT_program, n_read, (GMT_LONG)n_pitched, n_lost);

	/* Ready to go. */

	n_out = (Ctrl->W.weighted[GMT_OUT]) ? 4 : 3;
	if (Ctrl->E.active) n_out += 3;
	w_col = n_out - 1;	/* Weights always reported in last output column */

	/* Sort on index and Z value */

	qsort((void *)data, n_pitched, sizeof (struct BLK_DATA), BLK_compare_index_z);

	/* Find n_in_cell and write appropriate output  */

	first_in_cell = n_cells_filled = nz = 0;
	while (first_in_cell < n_pitched) {
		weight = data[first_in_cell].a[BLK_W];
		if (Ctrl->E.active) {
			if (nz == nz_alloc) nz_alloc = GMT_alloc_memory ((void **)&z_tmp, nz, nz_alloc, sizeof (double), GMT_program);
			z_tmp[0] = data[first_in_cell].a[BLK_Z];
			nz = 1;
		}
		if (Ctrl->C.active) {
			j = data[first_in_cell].i / ((GMT_LONG)h.nx);
			i = data[first_in_cell].i % ((GMT_LONG)h.nx);
			out[GMT_X] = GMT_i_to_x (i, h.x_min, h.x_max, h.x_inc, h.xy_off, h.nx);
			out[GMT_Y] = GMT_j_to_y (j, h.y_min, h.y_max, h.y_inc, h.xy_off, h.ny);
		}
		else {
			out[GMT_X] = data[first_in_cell].a[BLK_X];
			out[GMT_Y] = data[first_in_cell].a[BLK_Y];
		}
		first_in_new_cell = first_in_cell + 1;
		while ( (first_in_new_cell < n_pitched) && (data[first_in_new_cell].i == data[first_in_cell].i) ) {
			weight += data[first_in_new_cell].a[BLK_W];	/* Summing up weights */
			if (mode_xy) {
				out[GMT_X] += data[first_in_new_cell].a[BLK_X];
				out[GMT_Y] += data[first_in_new_cell].a[BLK_Y];
			}
			if (Ctrl->E.active) {	/* Must get a temporary copy of the sorted z array */
				if (nz == nz_alloc) nz_alloc = GMT_alloc_memory ((void **)&z_tmp, nz, nz_alloc, sizeof (double), GMT_program);
				z_tmp[nz] = data[first_in_new_cell].a[BLK_Z];
				nz++;
			}
			first_in_new_cell++;
		}
		n_in_cell = first_in_new_cell - first_in_cell;
		if (n_in_cell > 2) {	/* data are already sorted on z; get z mode  */
			out[GMT_Z] = weighted_mode (&data[first_in_cell], weight, n_in_cell, 2);
			if (Ctrl->Q.active) {
				i_n_in_cell = 1.0 / n_in_cell;
				out[GMT_X] *= i_n_in_cell;
				out[GMT_Y] *= i_n_in_cell;
			}
			else if (mode_xy) {
				qsort((void *)&data[first_in_cell], (size_t)n_in_cell, sizeof (struct BLK_DATA), BLK_compare_x);
				out[GMT_X] = weighted_mode (&data[first_in_cell], weight, n_in_cell, 0);

				qsort((void *)&data[first_in_cell], (size_t)n_in_cell, sizeof (struct BLK_DATA), BLK_compare_y);
				out[GMT_Y] = weighted_mode (&data[first_in_cell], weight, n_in_cell, 1);
			}
		}
		else if (n_in_cell == 2) {
			if (data[first_in_cell].a[BLK_W] > data[first_in_cell+1].a[BLK_W]) {
				out[GMT_Z] = data[first_in_cell].a[BLK_Z];
				if (Ctrl->Q.active) {
					out[GMT_X] *= 0.5;
					out[GMT_Y] *= 0.5;
				}
				else if (mode_xy) {
					out[GMT_X] = data[first_in_cell].a[BLK_X];
					out[GMT_Y] = data[first_in_cell].a[BLK_Y];
				}
			}
			else if (data[first_in_cell].a[BLK_W] < data[first_in_cell+1].a[BLK_W]) {
				out[GMT_Z] = data[first_in_cell+1].a[BLK_Z];
				if (Ctrl->Q.active) {
					out[GMT_X] *= 0.5;
					out[GMT_Y] *= 0.5;
				}
				else if (mode_xy) {
					out[GMT_X] = data[first_in_cell+1].a[BLK_X];
					out[GMT_Y] = data[first_in_cell+1].a[BLK_Y];
				}
			}
			else {
				if (mode_xy) {	/* Need average location */
					out[GMT_X] *= 0.5;
					out[GMT_Y] *= 0.5;
				}
				out[GMT_Z] = 0.5 * (data[first_in_cell].a[BLK_Z] + data[first_in_cell+1].a[BLK_Z]);
			}
		}
		else
			out[GMT_Z] = data[first_in_cell].a[BLK_Z];

		if (Ctrl->E.active) {
			out[4] = z_tmp[0];	/* Low value */
			out[5] = z_tmp[nz-1];	/* High value */
			/* Turn z_tmp into absolute deviations from the mode (out[GMT_Z]) */
			if (nz > 1) {
				for (index = 0; index < nz; index++) z_tmp[index] = fabs (z_tmp[index] - out[GMT_Z]);
				qsort ((void *)z_tmp, (size_t)nz, sizeof (double), GMT_comp_double_asc);
				out[3] = (nz%2) ? z_tmp[nz/2] : 0.5 * (z_tmp[(nz-1)/2] + z_tmp[nz/2]);
				out[3] *= 1.4826;	/* This will be LMS MAD-based scale */
			}
			else
				out[3] = GMT_d_NaN;
		}
		if (Ctrl->W.weighted[GMT_OUT]) out[w_col] = weight;
		
		GMT_output (GMT_stdout, n_out, out);

		n_cells_filled++;
		first_in_cell = first_in_new_cell;
	}

	if (gmtdefs.verbose) fprintf(stderr,"%s: N_cells_filled: %ld\n", GMT_program, n_cells_filled);

	GMT_free ((void *)data);
	GMT_free ((void *)z_tmp);

	Free_blockmode_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

double	weighted_mode (struct BLK_DATA *d, double wsum, GMT_LONG n, GMT_LONG k)
{
	/* Estimate mode by finding a maximum in the estimated
	   pdf of weighted data.  Estimate the pdf as the finite 
	   difference of the cumulative frequency distribution 
	   over points from i to j.  This has the form top/bottom,
	   where top is the sum of the weights from i to j, and
	   bottom is (data[j] - data[i]).  Strategy is to start
	   with i=0, j=n-1, and then move i or j toward middle
	   while j-i > n/2 and bottom > 0.  At end while, midpoint
	   of range from i to j is the mode estimate.  Choose 
	   to move either i or j depending on which one will
	   cause greatest increase in pdf estimate.  If a tie,
	   move both.
	   
	   Strictly, the pdf estimated this way would need to be
	   scaled by (1/wsum), but this is constant so we don't
	   use it here, as we are seeking a relative minimum.  
	   
	   I assumed n > 2 when I wrote this.  */
	   
	double	top, topj, topi, bottomj, bottomi, pj, pi;
	GMT_LONG	i, j, nh;

	i = 0;
	j = n - 1;
	nh = n / 2;
	top = wsum;

	while ((j-i) > nh) {
		topi = top - d[i].a[BLK_W];
		topj = top - d[j].a[BLK_W];
		bottomi = d[j].a[k] - d[i+1].a[k];
		bottomj = d[j-1].a[k] - d[i].a[k];

		if (bottomj == 0.0) {
			return (d[j-1].a[k]);
		}
		else if (bottomi == 0.0) {
			return (d[i+1].a[k]);
		}
		else {
			pi = topi/bottomi;
			pj = topj/bottomj;
			if (pi > pj) {
				i++;
				top = topi;
			}
			else if (pi < pj) {
				j--;
				top = topj;
			}
			else {
				top -= (d[i].a[BLK_W] + d[j].a[BLK_W]);
                                i++;
                                j--;
			}
		}
	}
	return (0.5*(d[j].a[k] + d[i].a[k]));
}

#include "block_subs.c"
