/*-----------------------------------------------------------------
 *	$Id: x2sys_report.c 9923 2012-12-18 20:45:53Z pwessel $
 *
 *      Copyright (c) 1999-2013 by P. Wessel
 *      See LICENSE.TXT file for copying and redistribution conditions.
 *
 *      This program is free software; you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation; version 2 or any later version.
 *
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *
 *      Contact info: gmt.soest.hawaii.edu
 *--------------------------------------------------------------------*/
/* x2sys_report will read the crossover data base and report on the statistics
 * of the crossovers for each track and overall.
 *
 * Author:	Paul Wessel
 * Date:	20-SEPT-2008
 * Version:	1.0, based on the spirit of the old x_system code x_report
 *		but completely rewritten from the ground up.
 *
 */

#include "x2sys.h"

struct COE_REPORT {	/* Holds summary info for each track */
	int nx;	/* Total number of COE for this track */
	double mean, stdev, rms;
	double sum, sum2, W;
	double d_max;	/* Length of track in distance units */
};

struct COE_ADJUST {	/* Holds adjustment spline knots */
	double d;	/* Distance along track */
	double c;	/* value to subtract from observation at crossover */
};

struct COE_ADJLIST {	/* Array with the growing arrays of COE_ADJUST per track */
	struct COE_ADJUST *K;
	int n, n_alloc;
};

int main (int argc, char **argv)
{
	char *TAG = CNULL, *column = CNULL, *track = CNULL, *ignorelist = CNULL, *dbase = NULL;
	char **trk_name = NULL, *correction_table = CNULL;
	struct X2SYS_INFO *s = NULL;
	struct X2SYS_BIX B;
	struct X2SYS_COE_PAIR *P = NULL;
	struct COE_REPORT *R = NULL;
	struct MGD77_CORRTABLE **CORR = NULL;
	GMT_LONG error = FALSE;
	GMT_LONG internal = FALSE;	/* FALSE if only external xovers are needed */
	GMT_LONG external = TRUE;	/* FALSE if only internal xovers are needed */
	GMT_LONG single = FALSE;		/* TRUE if a single track is specified -S */
	GMT_LONG apply_corrections = FALSE;	/* TRUE if -L is used */
	GMT_LONG adjust = FALSE;		/* TRUE and with -L will create an adjustment spline knots file per track */
	int i, k, n, coe_kind, nx_min = 0, n_use, n_tracks;
	GMT_LONG p, np, nx, Tnx = 0;
	double wesn[4], sum, sum2, sum_w, Tsum, Tsum2, COE, sign, scale, corr[2] = {0.0, 0.0};
	double Tmean, Tstdev, Trms;
	
	argc = (int)GMT_begin (argc, argv);
	
	for (i = (int)strlen(argv[0]); i >= 0 && argv[0][i] != '/'; i--);
	X2SYS_program = &argv[0][i+1];	/* Name without full path */
	memset ((void *)wesn, 0, 4*sizeof(double));	/* 0 means no box specified */
	
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				case 'V':
				case 'R':
				case '\0':
					error += GMT_parse_common_options (argv[i], &wesn[0], &wesn[1], &wesn[2], &wesn[3]);
					break;
				case 'A':
					adjust = TRUE;
					break;
				case 'C':
					column = &argv[i][2];
					break;
				case 'I':
					ignorelist = &argv[i][2];
					break;
				case 'L':	/* Crossover correction table */
					correction_table = &argv[i][2];
					apply_corrections = TRUE;
					break;
				case 'N':
					nx_min = atoi (&argv[i][2]);
					break;
				case 'Q':	/* Specify internal or external only */
					if (argv[i][2] == 'e') {external = TRUE; internal = FALSE;}
					if (argv[i][2] == 'i') {external = FALSE; internal = TRUE;}
					break;
				case 'S':
					track = &argv[i][2];
					single = TRUE;
					break;
				case 'T':
					TAG = &argv[i][2];
					break;
				default:
					error = TRUE;
					fprintf (stderr, "%s: Unrecognized option -%c\n", GMT_program, argv[i][1]);
					break;
			}
		}
		else
			dbase = argv[i];
	}

	if (argc == 1 || error || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "x2sys_report %s - Report statistics from crossover data base\n\n", X2SYS_VERSION);
		fprintf (stderr, "usage: x2sys_report -C<column> -T<TAG> [<COEdbase>] [-A] [-I<ignorelist>] [-L[<corrtable.txt>]]\n");
		fprintf (stderr, "\t [-N<nx_min>] [-Qe|i] [-S<track>] [%s] [-V]\n\n", GMT_Rgeo_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\t-C <column> is the name of the data column whose crossovers we want.\n");
		fprintf (stderr, "\t-T <TAG> is the system tag for the data set.\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t<COEdbase> File with crossover error data base [stdin]\n");
		fprintf (stderr, "\t-A Create adjustment splines per track to redistribute COEs between tracks\n");
		fprintf (stderr, "\t   according to their relative weight.\n");
		fprintf (stderr, "\t-I List of tracks to ignore [Use all tracks]s\n");
		fprintf (stderr, "\t-L Subtract systematic corrections from the data. If no correction file is given,\n");
		fprintf (stderr, "\t   the default file <TAG>_corrections.txt in $X2SYS_HOME/<TAG> is assumed.\n");
		fprintf (stderr, "\t-N Only output results for tracks with more than <nx_min> crossovers [0, i.e., report all tracks]\n");
		fprintf (stderr, "\t-Q Append e or i for external or internal crossovers [Default is external]\n");
		GMT_explain_option ('R');
		fprintf (stderr, "\t   [Default region is the entire data domain]\n");
		fprintf (stderr, "\t-S Return only crossovers involving this track [Use all tracks]\n");
		GMT_explain_option ('V');
		exit (EXIT_FAILURE);
	}

	if (!column) {
		fprintf (stderr, "%s: ERROR: Must use -C to specify observation of interest\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	
	/* Initialize system via the tag */
	
	x2sys_err_fail (x2sys_set_system (TAG, &s, &B, &GMT_io), TAG);

	/* Verify that the chosen column is known to the system */
	
	if (column) x2sys_err_fail (x2sys_pick_fields (column, s), "-C");
	if (s->n_out_columns != 1) {
		fprintf (stderr, "%s: ERROR: -C must specify a single column name\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	
	/* Select internal, external, or both */
	
	coe_kind = 0;
	if (internal) coe_kind |= 1;
	if (external) coe_kind |= 2;
	if (coe_kind == 0) coe_kind = 2;	/* External */
	
	/* Read the entire data base; note the -I, R and -S options are applied during reading */
	
	if (gmtdefs.verbose) fprintf (stderr, "%s: Read crossover database %s...\n", GMT_program, dbase);
	np = x2sys_read_coe_dbase (s, dbase, ignorelist, wesn, column, coe_kind, track, &P, &nx, &n_tracks);
	if (gmtdefs.verbose) fprintf (stderr, "Found %ld pairs and a total of %ld crossover records.\n", np, nx);

	if (np == 0 && nx == 0) {	/* End here since nothing was allocated */
		x2sys_end (s);
		GMT_end (argc, argv);
		exit (EXIT_SUCCESS);
	}
	
	R = (struct COE_REPORT *) GMT_memory (VNULL, n_tracks, sizeof (struct COE_REPORT), GMT_program);
	trk_name = (char **) GMT_memory (VNULL, n_tracks, sizeof (char *), GMT_program);
	
	for (p = 0; p < np; p++) {	/* Get track name list */
		for (k = 0; k < 2; k++) trk_name[P[p].id[k]] = P[p].trk[k];
	}
	
	if (apply_corrections) {	/* Load an ephemeral correction table */
		x2sys_get_corrtable (s, correction_table, n_tracks, trk_name, column, NULL, NULL, &CORR);
	}

	Tsum = Tsum2 = 0.0;
	for (p = 0; p < np; p++) {	/* For each pair of tracks that generated crossovers */
		for (i = n = 0, sum = sum2 = 0.0; i < P[p].nx; i++) {
			if (apply_corrections) {
				corr[0] = MGD77_Correction_Rec (CORR[P[p].id[0]][COE_Z].term, P[p].COE[i].data[0], NULL);
				corr[1] = MGD77_Correction_Rec (CORR[P[p].id[1]][COE_Z].term, P[p].COE[i].data[1], NULL);
			}
			COE = (P[p].COE[i].data[0][COE_Z] - corr[0]) - (P[p].COE[i].data[1][COE_Z] - corr[1]);
			if (GMT_is_dnan(COE)) continue;
			sum += COE;	sum2 += COE * COE;	n++;
			Tsum += COE;	Tsum2 += COE * COE;	Tnx++;
		}
		for (k = 0, sign = 1.0; n && k < 2; k++) {
			R[P[p].id[k]].nx += n;
			R[P[p].id[k]].sum += sign * sum;
			R[P[p].id[k]].sum2 += sum2;
			R[P[p].id[k]].d_max = P[p].dist[k];	/* Just copy over max distance for later */
			sign = -1.0;
		}
	}
	
	for (k = n_use = 0, sum_w = 0.0; k < n_tracks; k++) {	/* Calculate normalized weights */
		if (R[k].nx <= nx_min) continue;	/* Not enough COEs */
		if (R[k].nx && R[k].sum2 > 0.0) {
			R[k].W = R[k].nx / R[k].sum2;
			sum_w += R[k].W;
			n_use++;
		}
		else
			R[k].W = GMT_d_NaN;
	}
	scale = (n_use && sum_w > 0.0) ? n_use / sum_w : 1.0;
	
	/* Time to issue output */
	
	fprintf (stdout, "# Tag: %s %s\n", TAG, column);
	fprintf (stdout, "# Command: %s", GMT_program);
	if (!dbase) fprintf (stdout, " [stdin]");
	for (k = 1; k < argc; k++) fprintf (stdout, " %s", argv[k]);
	fprintf (stdout, "\n#track\tN\tmean\tstdev\trms\tweight[%d]\n", n_use);
	Tmean = (Tnx) ? Tsum / Tnx : GMT_d_NaN;
	Tstdev = (Tnx > 1) ? sqrt ((Tnx * Tsum2 - Tsum * Tsum) / (Tnx * (Tnx - 1.0))) : GMT_d_NaN;
	Trms = (Tnx) ? sqrt (Tsum2 / Tnx) : GMT_d_NaN;
	printf ("TOTAL\t%ld\t%g\t%g\t%g\t1\n", Tnx, Tmean, Tstdev, Trms);
	for (k = 0; k < n_tracks; k++) {	/* For each track that generated crossovers */
		if (R[k].nx <= nx_min) continue;			/* Not enough COEs */
		if (!GMT_is_dnan (R[k].W)) R[k].W *= scale;
		R[k].mean = (R[k].nx) ? R[k].sum / R[k].nx : GMT_d_NaN;
		R[k].stdev = (R[k].nx > 1) ? sqrt ((R[k].nx * R[k].sum2 - R[k].sum * R[k].sum) / (R[k].nx * (R[k].nx - 1.0))) : GMT_d_NaN;
		R[k].rms = (R[k].nx) ? sqrt (R[k].sum2 / R[k].nx) : GMT_d_NaN;
		printf ("%s\t%d\t%g\t%g\t%g\t%g\n", trk_name[k], R[k].nx, R[k].mean, R[k].stdev, R[k].rms, R[k].W);
	}
	
	if (adjust) {	/* Create track adjustment spline files for each track */
		int n_out, n1;
		char file[BUFSIZ];
		double out[2], z[2], z_ij;
		FILE *fp;
		struct COE_ADJLIST *adj;
		int comp_structs (const void *point_1, const void *point_2);

		GMT_io.out_col_type[GMT_X] = GMT_io.out_col_type[GMT_Y] = GMT_IS_FLOAT;	/* Since we will write (dist, COE) pairs */
		
		adj = (struct COE_ADJLIST *) GMT_memory (VNULL, n_tracks, sizeof (struct COE_ADJLIST), GMT_program);
		for (p = 0; p < np; p++) {	/* For each pair of tracks that generated crossovers */
			for (i = n = 0; i < P[p].nx; i++) {	/* For each COE between this pair */
				for (k = 0; k < 2; k++) {
					z[k] = P[p].COE[i].data[k][COE_Z];
					if (apply_corrections) z[k] -= MGD77_Correction_Rec (CORR[P[p].id[k]][COE_Z].term, P[p].COE[i].data[k], NULL);
				}
				z_ij = (z[0] * R[P[p].id[0]].W + z[1] * R[P[p].id[1]].W) / (R[P[p].id[0]].W + R[P[p].id[1]].W);	/* Desired z-value at crossover */
				if (GMT_is_dnan(z_ij)) continue;
				for (k = 0; k < 2; k++) {
					if (R[P[p].id[k]].nx <= nx_min) continue;	/* This track will not have enough total COE to be used in the end */
					if (adj[P[p].id[k]].n >= adj[P[p].id[k]].n_alloc) {	/* So first time both are zero and we allocate first */
						if (adj[P[p].id[k]].n_alloc) adj[P[p].id[k]].n_alloc <<= 1; else adj[P[p].id[k]].n_alloc = GMT_SMALL_CHUNK;
						adj[P[p].id[k]].K = (struct COE_ADJUST *) GMT_memory ((void *)adj[P[p].id[k]].K, adj[P[p].id[k]].n_alloc, sizeof (struct COE_ADJUST), GMT_program);
					}
					adj[P[p].id[k]].K[adj[P[p].id[k]].n].d = P[p].COE[i].data[k][COE_D];	/* Distance along current track at COE */
					adj[P[p].id[k]].K[adj[P[p].id[k]].n].c = z_ij - z[k];			/* Observed difference at COE to be adjusted */
					adj[P[p].id[k]].n++;
				}
			}
		}
		/* Time to create the track adjustment spline files */
		for (k = 0; k < n_tracks; k++) {	/* For each track that generated crossovers */
			if (adj[k].n <= nx_min) continue;			/* Not enough COEs */
			adj[k].K[adj[k].n].d = adj[k].K[adj[k].n].c = 0.0;	/* Add in anchor point (0,0) */
			adj[k].n++;
			adj[k].K[adj[k].n].d = R[k].d_max;			/* Add in anchor point (d_max,0) */
			adj[k].K[adj[k].n].c = 0.0;
			adj[k].n++;
			
			qsort((void *)adj[k].K, (size_t)adj[k].n, sizeof(struct COE_ADJUST), comp_structs);
			sprintf (file, "%s/%s/%s.%s.adj", X2SYS_HOME, TAG, trk_name[k], column);
			if ((fp = GMT_fopen (file, "w")) == NULL) {
				fprintf (stderr, "%s: Unable to create file %s!\n", GMT_program, file);
				exit (EXIT_FAILURE);
			}
			n1 = adj[k].n - 1;
			out[1] = 0.0;
			for (i = n_out = 0; i <= n1; i++) {
				out[0] = adj[k].K[i].d;
				out[1] += adj[k].K[i].c;
				n_out++;
				if (i == n1 || out[0] < adj[k].K[i+1].d) {	/* Time to output */
					out[1] /= n_out;
					GMT_output (fp, 2, out);
					out[1] = 0.0;
					n_out = 0;
				}
			}
			GMT_fclose (fp);
			GMT_free ((void *)adj[k].K);
		}
		GMT_free ((void *)adj);
	}
	
	/* Done, free up data base array, etc */
	
	x2sys_free_coe_dbase (P, np);
	GMT_free ((void *)trk_name);
	GMT_free ((void *)R);

	if (apply_corrections) {
		MGD77_Free_Correction (CORR, n_tracks);
	}
	x2sys_end (s);

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

int comp_structs (const void *point_1, const void *point_2) { /* Sort ADJ structure on distance */
        if ( ((struct COE_ADJUST *)point_1)->d < ((struct COE_ADJUST *)point_2)->d)
                return(-1);
        else if ( ((struct COE_ADJUST *)point_1)->d > ((struct COE_ADJUST *)point_2)->d)
                return(+1);
        else
                return(0);
}
