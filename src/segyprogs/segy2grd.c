/*--------------------------------------------------------------------
 *	$Id: segy2grd.c 9923 2012-12-18 20:45:53Z pwessel $
 *
 *	Copyright (c) 1991-2013 by T. Henstock
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
 * segy2grd was modified from xyz2grd by T. J. Henstock, June 2002
 * to read a segy file and generate a
 * corresponding grd-file.
 *
 * Author:	Tim Henstock (then@noc.soton.ac.uk)
 * Date:	30-JUN-2002
 * Version:	3.4.1
 */
 
#include "gmt.h"
#include "segy_io.h"

int main (int argc, char **argv)
{
	GMT_LONG error = FALSE, got_input = FALSE, read_cont = FALSE, count = FALSE, average = TRUE;
	GMT_LONG nodata_set = FALSE, read_cdp = FALSE, read_offset = FALSE;
#ifdef WORDS_BIGENDIAN
	GMT_LONG swap_bytes = FALSE;
#else
	GMT_LONG swap_bytes = TRUE;
#endif

	GMT_LONG i, ij, ii, jj, nm, n_read = 0, n_filled = 0, n_used = 0, one_or_zero, *flag = NULL;
	GMT_LONG n_empty = 0, n_stuffed = 0, n_bad = 0, n_confused = 0;
	GMT_LONG check, n_traces=10000, n_samp=0, n_sampr=0, byte_x = 0, ix, isamp, ij0;

	GMT_LONG pixel = FALSE;
	double w, e, s, n, dx = 0.0, dy = 0.0, xy_off, idy, no_data_d;
	double s_int = 0.0, x_coord_scale=1.0, x0, yval;

	float no_data_f, *a = NULL;

	char *grdfile = NULL, line[BUFSIZ], input[BUFSIZ];

	FILE *fpi = NULL;

	struct GRD_HEADER grd;

/* SEGY parameters */
	char reelhead[3200];
	float *data = NULL;
	SEGYHEAD *header = NULL;
	char *head = NULL;
	int head2;
	SEGYREEL binhead;



	argc = (int)GMT_begin (argc, argv);

	grdfile = CNULL;
	input[0] = 0;
	w = e = s = n = 0.0;
	no_data_f = GMT_f_NaN;
	no_data_d = GMT_d_NaN;

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				/* Common parameters */

				case 'R':
				case 'V':
				case '\0':
					error += GMT_parse_common_options (argv[i], &w, &e, &s, &n);
					break;

				case 'A':
					if (argv[i][2] == 'n')
						count = TRUE;
					else if (argv[i][2] == '\0' || argv[i][2] == 'z')
						average = FALSE;
					else {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -A option:  Select -An or -A[z]\n", GMT_program);
						error++;
					}
					break;
				case 'D':
					strcpy (input, &argv[i][2]);
					got_input = TRUE;
					break;
				case 'F':
					pixel = TRUE;
					break;
				case 'G':
					grdfile = &argv[i][2];
					break;
				case 'I':
					GMT_getinc (&argv[i][2], &dx, &dy);
					break;
				case 'N':
					if (!argv[i][2]) {
						fprintf (stderr, "%s: GMT SYNTAX ERROR -N option:  Must specify value or NaN\n", GMT_program);
						error++;
					}
					else {
						no_data_d = (argv[i][2] == 'N' || argv[i][2] == 'n') ? GMT_d_NaN : atof (&argv[i][2]);
						no_data_f = (float)no_data_d;
						nodata_set = TRUE;
					}
					break;
				case 'X': /* -X and -Y can be changed in gmt routines to lower case...*/
				case 'x':
					x_coord_scale = atof (&argv[i][2]);
					break;
				case 'Y':
				case 'y':
					s_int = atof (&argv[i][2]);
					fprintf(stderr," %s: Overriding sample interval s_int = %f\n", GMT_program, s_int);
					break;
				case 'L':
					n_sampr = atoi (&argv[i][2]);
					break;
				case 'M':
					n_traces = atoi (&argv[i][2]);
					break;
				/* variable spacing */
				case 'S':
					switch(argv[i][2]){
						case 'o':
							read_offset = TRUE;
							break;
						case 'c':
							read_cdp = TRUE;
							break;
						case 'b':
							byte_x = atoi (&argv[i][3]);
							break;
					}
					break;
				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else if ((fpi = fopen (argv[i], "rb")) == NULL) {
			fprintf (stderr, "%s: Cannot find segy file %s\n", GMT_program,argv[i]);
			exit (EXIT_FAILURE);
		}
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "segy2grd %s - Converting segy data to a GMT grid file\n\n", GMT_VERSION);
		fprintf (stderr, "usage: segy2grd <segyfile> -G<grdfile> %s %s\n", GMT_Id_OPT, GMT_Rgeo_OPT);
		fprintf (stderr, "\t[-A[n|z]] [%s] [-F] \n", GMT_GRDEDIT);
		fprintf (stderr, "\t[-N<nodata>] [-X<x-scale>] [-Y<s_int>] [-V] \n");

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\tsegyfile(s) is an IEEE floating point SEGY file. Traces are all assumed to start at 0 time/depth\n");
		fprintf (stderr, "\t-G to name the output grid file.\n");
		fprintf (stderr, "\t-I specifies grid size(s). \n");
		GMT_explain_option ('R');
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-A (or -Az): Add multiple entries at the same node.\n");
		fprintf (stderr, "\t   Append n (-An): Count number of multiple entries per node instead.\n");
		fprintf (stderr, "\t   [Default (no -A option) will compute mean values]\n");
		fprintf (stderr, "\t-D to enter header information.  Specify '=' to get default value\n");
		fprintf (stderr, "\t-F will force pixel registration [Default is grid registration]\n");
		fprintf (stderr, "\t-N set value for nodes without corresponding input sample [Default is NaN]\n");
		fprintf (stderr,"\t-S<header> to set variable spacing\n");
		fprintf (stderr,"\t\t<header> is c for cdp, o for offset, b<number> for 4-byte float starting at byte number\n");
		fprintf (stderr,"\t\tIf -S not set, assumes even spacing of samples at dx, dy supplied with -I\n");
		fprintf (stderr,"\t-L<nsamp> to override number of samples\n");
		fprintf (stderr,"\t-X<x-scale> applies scalar x-scale to coordinates in trace header to match the coordinates specified in -R\n");
		fprintf (stderr,"\t-Y<s_int> specifies sample interval as s_int if incorrect in the SEGY file\n");
		fprintf (stderr,"\t-M<ntraces> to fix number of traces. Default reads all traces.\n\t\t-M0 will read number in binary header, -Mn will attempt to read only n traces.\n");
		GMT_explain_option ('V');
		GMT_explain_option ('.');
		exit (EXIT_FAILURE);
	}

	GMT_check_lattice (&dx, &dy, &pixel, NULL);

	read_cont = (!read_cdp && !read_offset && !byte_x);
	if (!project_info.region_supplied) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify -R option\n", GMT_program);
		error++;
	}
	if (dx <= 0.0 || dy <= 0.0) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR -I option.  Must specify positive increment(s)\n", GMT_program);
		error++;
	}
	if (!grdfile) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR option -G:  Must specify output file\n", GMT_program);
		error++;
	}
	if ((read_cdp && read_offset) || (read_cdp && byte_x) || (read_offset && byte_x)){
		fprintf(stderr, "%s: GMT SYNTAX ERROR option -S:  Must only specify one of cdp, offset, byte\n", GMT_program);
		error++;
	}

	if (error) exit (EXIT_FAILURE);

	GMT_grd_init (&grd, argc, argv, FALSE);

	/* if (!project_info.region) d_swap (s, e); */  /* Got w/s/e/n, make into w/e/s/n */

	/* Decode grd information given, if any */

	if (got_input) GMT_decode_grd_h_info (input, &grd);

	grd.node_offset = (int)pixel;
	one_or_zero = 1 - grd.node_offset;
	xy_off = 0.5 * grd.node_offset;

	grd.nx = (int)(irint ((e-w)/dx) + one_or_zero);
	grd.ny = (int)(irint ((n-s)/dy) + one_or_zero);
	grd.x_min = w;	grd.x_max = e;
	grd.y_min = s;	grd.y_max = n;
	grd.x_inc = dx;	grd.y_inc = dy;

	GMT_err_fail (GMT_grd_RI_verify (&grd, 1), grdfile);

	if (gmtdefs.verbose) fprintf (stderr, "%s: nx = %d  ny = %d\n", GMT_program, grd.nx, grd.ny);

	nm = GMT_get_nm (grd.nx, grd.ny);

	a = (float *) GMT_memory (VNULL, (size_t)nm, sizeof (float), GMT_program);
	flag = (GMT_LONG *) GMT_memory (VNULL, (size_t)nm, sizeof (GMT_LONG), GMT_program);

	idy = 1.0 / dy;
	ij = -1;	/* Will be incremented to 0 or set first time around */

/* read in reel headers from segy file */
	if (fpi == NULL) fpi = GMT_stdin;
	if ((check = get_segy_reelhd (fpi, reelhead)) != TRUE) exit(1);
	if ((check = get_segy_binhd (fpi, &binhead)) != TRUE) exit(1);

	if(swap_bytes){
/* this is a little-endian system, and we need to byte-swap ints in the reel header - we only
use a few of these*/
		if (gmtdefs.verbose) fprintf(stderr, "%s: swapping bytes for ints in the headers\n",GMT_program);
		binhead.num_traces = GMT_swab2(binhead.num_traces);
		binhead.nsamp = GMT_swab2(binhead.nsamp);
		binhead.dsfc = GMT_swab2(binhead.dsfc);
		binhead.sr = GMT_swab2(binhead.sr);
	}


/* set parameters from the reel headers */
	if (!n_traces)
		n_traces = binhead.num_traces;

	if (gmtdefs.verbose) fprintf(stderr, "%s: Number of traces in header is %ld\n", GMT_program, n_traces);


	if (!n_sampr){/* number of samples not overridden*/
		n_sampr = binhead.nsamp;
		fprintf(stderr,"%s: Number of samples per trace is %ld\n", GMT_program, n_sampr);
	}
	else if ((n_sampr != binhead.nsamp) && (binhead.nsamp))
		fprintf(stderr,"%s: warning nsampr input %ld, nsampr in header %d\n", GMT_program, n_sampr,  binhead.nsamp);

	if (!n_sampr){ /* no number of samples still - a problem! */
		fprintf(stderr, "%s: Error, number of samples per trace unknown\n", GMT_program);
		exit(EXIT_FAILURE);
	}

	if(gmtdefs.verbose) 
		fprintf(stderr, "%s: Number of samples for reel is %ld\n", GMT_program, n_sampr);

	if(binhead.dsfc != 5) fprintf(stderr, "pssegy: WARNING data not in IEEE format\n");

	if (!s_int){
		s_int = (double) binhead.sr; /* sample interval of data (microseconds) */
		s_int /= 1000000.0;
		fprintf(stderr,"%s: Sample interval is %f s\n", GMT_program, s_int);
	}
	else if ((s_int != binhead.sr) && (binhead.sr)) /* value in header overridden by input */
		fprintf(stderr, "%s: Warning s_int input %f, s_int in header %f\n", GMT_program, s_int, (float)binhead.sr);

	if (!s_int){ /* still no sample interval at this point is a problem! */
		fprintf(stderr, "%s: Error, no sample interval in reel header\n", GMT_program);
		exit(EXIT_FAILURE);
	}
	if (read_cont && (s_int != dy)){
		fprintf(stderr, "%s: Warning, grid spacing != sample interval, setting sample interval to grid spacing\n", GMT_program);
		s_int = dy;
	}

	if(dy < s_int)
		fprintf(stderr, "%s: Warning, grid spacing < sample interval, expect gaps in output....\n", GMT_program);



/* starts reading actual data here....... */

	if (read_cont) {	/* old-style segy2grd */
		ix=0;
		for (ij=0; ij<nm; ij++) a[ij] = no_data_f;
		if(grd.nx<n_traces){
			fprintf(stderr,"%s: Warning, number of traces in header > size of grid. Reading may be truncated\n", GMT_program);
			n_traces = grd.nx;
		}
		while ((ix<n_traces) && (header = get_segy_header(fpi))) {
			if (swap_bytes){
/* need to permanently byte-swap number of samples in the trace header */
				header->num_samps=GMT_swab4(header->num_samps);
			 	header->sampleLength=GMT_swab2(header->sampleLength);
			}

			data = (float *) get_segy_data(fpi, header); /* read a trace */
			/* get number of samples in _this_ trace or set to number in reel header */
			if ( !(n_samp = samp_rd(header)) ) n_samp = n_sampr;

			ij0 = (GMT_LONG)(s*idy);
			if (n_samp - ij0 > grd.ny) n_samp = grd.ny + ij0;

			if(swap_bytes){ /* need to swap the order of the bytes in the data even though assuming IEEE format */
				int *intdata = (int *) data;
				for (isamp=0; isamp<n_samp; isamp++){
					intdata[isamp]=GMT_swab4(intdata[isamp]);
				}
			}

			for (ij = ij0; ij < n_samp ; ij++) {  /* n*idy is index of first sample to be included in the grid */
				a[ix + grd.nx*(grd.ny+ij0-ij-1)] = data[ij];
			}

		free (data);
		free (header);
		ix++;
		}
	}
	else {	/* Get trace data and position by headers */
		ix=0;
		while ((ix<n_traces) && (header = get_segy_header(fpi))){   /* read traces one by one */
			if (read_offset){ /* plot traces by offset, cdp, or input order */
				int32_t offset = ((swap_bytes)? GMT_swab4(header->sourceToRecDist): header->sourceToRecDist);
				x0 = (double) offset;
			}
			else if (read_cdp){
				int32_t cdpval = ((swap_bytes)? GMT_swab4(header->cdpEns): header->cdpEns);
				x0 = (double) cdpval;
			}
			else if (byte_x){ /* ugly code - want to get value starting at byte_x of header into a double... */
				head = (char *) header;
			memcpy(&head2, &head[byte_x], 4); /* edited to fix bug where 8bytes were copied from head.
                                                Caused by casting to a long directly from char array*/ 
			x0 = (double) ((swap_bytes)? GMT_swab4(head2): head2);
			}
			else
				x0 = (1.0 + (double) ix);

			x0 *= x_coord_scale;

			if (swap_bytes){
/* need to permanently byte-swap some things in the trace header 
do this after getting the location of where traces are plotted in case the general byte_x case
overlaps a defined header in a strange way */
				 header->sourceToRecDist=GMT_swab4(header->sourceToRecDist);
			 	header->sampleLength=GMT_swab2(header->sampleLength);
				 header->num_samps=GMT_swab4(header->num_samps);
			}

			data = (float *) get_segy_data(fpi, header); /* read a trace */
			/* get number of samples in _this_ trace (e.g. OMEGA has strange ideas about SEGY standard)
			or set to number in reel header */
			if ( !(n_samp = samp_rd(header)) ) n_samp = n_sampr;

			if(swap_bytes){ /* need to swap the order of the bytes in the data even though assuming IEEE format */
				int *intdata = (int *) data;
				for (isamp=0; isamp<n_samp; isamp++){
					intdata[isamp]=GMT_swab4(intdata[isamp]);
				}
			}

			if (!(x0 < w || x0 > e)){	/* inside x-range */
				/* find horizontal grid pos of this trace */
				ii = GMT_x_to_i (x0, grd.x_min, grd.x_inc, xy_off, grd.nx);
				if (ii == grd.nx) ii--, n_confused++;
				for (isamp = 0; isamp< n_samp; isamp++){
					yval = isamp*s_int;
					if (!(yval < s || yval > n)){	/* inside y-range */
						jj = GMT_y_to_j (yval, grd.y_min, grd.y_inc, xy_off, grd.ny);
						if (jj == grd.ny) jj--, n_confused++;
						ij = jj * grd.nx + ii;
						a[ij] += data[isamp];	/* Add up incase we must average */
						flag[ij]++;
						n_used++;
					}
				}
			}
			free(data);
			free(header);
			ix++;
		}

		for (ij = 0; ij < nm; ij++) {	/* Check if all nodes got one value only */
			if (flag[ij] == 1) {
				if (count) a[ij] = 1.0;
				n_filled++;
			}
			else if (flag[ij] == 0) {
				n_empty++;
				a[ij] = no_data_f;
			}
			else {	/* More than 1 value went to this node */
				if (count)
					a[ij] = (float)flag[ij];
				else if (average)
					a[ij] /= (float)flag[ij];
				n_filled++;
				n_stuffed++;
			}
		}

		if (gmtdefs.verbose) {
			sprintf (line, "%s\n", gmtdefs.d_format);
			fprintf (stderr, "%s:  n_read: %ld  n_used: %ld  n_filled: %ld  n_empty: %ld set to ", GMT_program,
				n_read, n_used, n_filled, n_empty);
			(GMT_is_dnan (no_data_d)) ? fprintf (stderr, "NaN\n") : fprintf (stderr, line, no_data_d);
			if (n_bad) fprintf (stderr, "%s: %ld records unreadable\n", GMT_program, n_bad);
			if (n_stuffed) fprintf (stderr, "%s: Warning - %ld nodes had multiple entries that were averaged\n", GMT_program, n_stuffed);
			if (n_confused) fprintf (stderr, "%s: Warning - %ld values gave bad indices: Pixel vs gridline confusion?\n", GMT_program, n_confused);
		}
	}

	GMT_err_fail (GMT_write_grd (grdfile, &grd, a, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), grdfile);

	GMT_free ((void *)a);
	GMT_free ((void *)flag);

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}
