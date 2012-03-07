/*-----------------------------------------------------------------
 *	$Id: x2sys_list.c,v 1.38 2011/07/11 19:22:07 guru Exp $
 *
 *      Copyright (c) 1999-2011 by P. Wessel
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
/* x2sys_list will read the crossover data base and output a subset of
 * the crossovers in a format determined by the options.
 *
 * Author:	Paul Wessel
 * Date:	20-SEPT-2008
 * Version:	1.0, based on the spirit of the old x_system code x_list
 *		but completely rewritten from the ground up.
 *
 */

#include "x2sys.h"

#define GMT_T	3	/* Just used to indicate abs time formatting */
#define LETTERS "acdhiInNtTvwxyz"

void dump_ascii_cols (double *val, int col, int n, GMT_LONG first);

int main (int argc, char **argv)
{
	char *TAG = CNULL, *column = CNULL, *fflags = CNULL, *track = CNULL, *ignorelist = CNULL;
	char *dbase = NULL, *correction_table = CNULL, *weight_list = NULL, **trk_name = NULL, **weight_name = NULL;
	char buffer[BUFSIZ];
	struct X2SYS_INFO *s = NULL;
	struct X2SYS_BIX B;
	struct X2SYS_COE_PAIR *P = NULL;
	GMT_LONG error = FALSE, mixed = FALSE, check_for_NaN = FALSE, both = FALSE, first, z_is_requested = FALSE;
	GMT_LONG internal = TRUE;	/* FALSE if only external xovers are needed */
	GMT_LONG external = TRUE;	/* FALSE if only internal xovers are needed */
	GMT_LONG single = FALSE;		/* TRUE if a single track is specified -S */
	GMT_LONG symm_check = FALSE;	/* TRUE if -Y is used */
	GMT_LONG nx_check = FALSE;	/* TRUE if -N is used */
	GMT_LONG apply_corrections = FALSE;	/* TRUE if -L is used */
	int i, j, k, coe_kind, one, two, n_items, n_out, nx_min = 0, n_tracks, n_weights = 0, *trk_nx = NULL;
	GMT_LONG p, np_use = 0, nx_use = 0, np, m, nx, id;
	double wesn[4], val[2], out[128], corr[2] = {0.0, 0.0}, sec_2_unit = 1.0, w_k, w;
	double fixed_weight = 1.0, *weights = NULL, asymm_max = 1.0, *trk_symm = NULL;
	struct MGD77_CORRTABLE **CORR = NULL;
	
	argc = (int)GMT_begin (argc, argv);
	
	for (i = (int)strlen(argv[0]); i >= 0 && argv[0][i] != '/'; i--);
	X2SYS_program = &argv[0][i+1];	/* Name without full path */
	memset ((void *)wesn, 0, 4*sizeof(double));	/* 0 means no box specified */
	
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				case 'M':
				case 'R':
				case 'V':
				case 'm':
				case 'o':
				case '\0':
					error += GMT_parse_common_options (argv[i], &wesn[0], &wesn[1], &wesn[2], &wesn[3]);
					break;

				case 'A':
					asymm_max = atof (&argv[i][2]);
					symm_check = TRUE;
					break;
				case 'C':
					column = &argv[i][2];
					break;
				case 'F':
					fflags = &argv[i][2];
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
					nx_check = TRUE;
					break;
				case 'Q':	/* Specify internal or external only */
					if (argv[i][2] == 'e') internal = FALSE;
					if (argv[i][2] == 'i') external = FALSE;
					break;
				case 'S':
					if (argv[i][2] == '+') {	/* Print info relative to both cruises */
						both  = TRUE;
						track = &argv[i][3];
					}
					else
						track = &argv[i][2];
					single = TRUE;
					break;
				case 'T':
					TAG = &argv[i][2];
					break;
				case 'W':
					weight_list = &argv[i][2];
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
		fprintf (stderr, "x2sys_list %s - Extract subset from crossover data base\n\n", X2SYS_VERSION);
		fprintf (stderr, "usage: x2sys_list -C<column> -T<TAG> [<COEdbase>] [-A<asymm_max] [-F<flags>] [-I<ignorelist>] [-L[<corrtable.txt>]]\n");
		fprintf (stderr, "\t [-N<nx_min>] [-Qe|i] [-S[+]<track>] [%s] [-V] [-W<weight>] [%s]\n\n", GMT_Rgeo_OPT, GMT_mo_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\t-C <column> is the name of the data column whose crossovers we want.\n");
		fprintf (stderr, "\t-T <TAG> is the system tag for the data set.\n");
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t<COEdbase> File with crossover error data base [stdin].\n");
		fprintf (stderr, "\t-A Return only crossovers whose distribution in time [or dist if no time]\n");
		fprintf (stderr, "\t   are fairly symmetric about the mid-point. Specify max abs value for\n");
		fprintf (stderr, "\t   asymmetry = (n_right - n_left)/(nright + n_left) [1, i.e., use all tracks].\n");
		fprintf (stderr, "\t-F Specify any combination of %s in the order of desired output:\n", LETTERS);
		fprintf (stderr, "\t   a Angle (<= 90) between the two tracks at the crossover.\n");
		fprintf (stderr, "\t   c Crossover error in chosen observable (see -C).\n");
		fprintf (stderr, "\t   d Distance along tracks at the crossover.\n");
		fprintf (stderr, "\t   h Heading along tracks at the crossover.\n");
		fprintf (stderr, "\t   i Signed time interval between the two tracks at the crossover.\n");
		fprintf (stderr, "\t   I Unsigned time interval between the two tracks at the crossover.\n");
		fprintf (stderr, "\t   n Names of the two tracks.\n");
		fprintf (stderr, "\t   N Id numbers of the two tracks.\n");
		fprintf (stderr, "\t   t Absolute time along tracks at the crossover.\n");
		fprintf (stderr, "\t   T Time since start of track along tracks at the crossover.\n");
		fprintf (stderr, "\t   v Speed along tracks at the crossover.\n");
		fprintf (stderr, "\t   w weight assigned to the crossover.\n");
		fprintf (stderr, "\t   x x-coordinate of the crossover.\n");
		fprintf (stderr, "\t   y y-coordinate of the crossover.\n");
		fprintf (stderr, "\t   z Observed values (see -C) along tracks at the crossover.\n");
		fprintf (stderr, "\t   Unless -S is specified, d,h,n,t,T,v,z will yield two columns.\n");
		fprintf (stderr, "\t-I List of tracks to ignore [Use all tracks].\n");
		fprintf (stderr, "\t-L Subtract systematic corrections from the data. If no correction file is given,\n");
		fprintf (stderr, "\t   the default file <TAG>_corrections.txt in $X2SYS_HOME/<TAG> is assumed.\n");
		fprintf (stderr, "\t-N Only output results for tracks with more than <nx_min> crossovers [Use all tracks].\n");
		fprintf (stderr, "\t-Q Append e or i for external or internal crossovers [Default is both].\n");
		GMT_explain_option ('R');
		fprintf (stderr, "\t   [Default region is the entire data domain].\n");
		fprintf (stderr, "\t-S Return only crossovers involving this track [Use all tracks].\n");
		fprintf (stderr, "\t   Prepend a '+' to make it print info relative to both tracks [Default is selected track].\n");
		GMT_explain_option ('V');
		fprintf (stderr, "\t-W If argument can be opened as a file then we expect a List of tracks and their\n");
		fprintf (stderr, "\t   relative weights; otherwise the argument is the constant weight for all tracks [1].\n");
		GMT_explain_option ('o');
		GMT_explain_option ('m');
		exit (EXIT_FAILURE);
	}

	if (symm_check && (asymm_max <= 0.0 || asymm_max > 1.0)) {
		fprintf (stderr, "%s: ERROR: -S: Asymmetry must be in the range 0-1.\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (GMT_io.multi_segments[GMT_OUT] && GMT_io.binary[GMT_OUT]) {
		fprintf (stderr, "%s: ERROR: Cannot use -M with binary output.\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	if (!fflags) {
		fprintf (stderr, "%s: ERROR: Must use -F to specify output items.\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	n_items = (int)strlen (fflags);
	for (i = 0; i < n_items; i++) {
		if (!strchr (LETTERS, (int)fflags[i])) {
			fprintf (stderr, "%s: ERROR -F: Unknown item %c.\n", GMT_program, fflags[i]);
			exit (EXIT_FAILURE);			
		}
		if (fflags[i] == 'n') mixed = TRUE;		/* Both numbers and text - cannot use binary output */
		if (fflags[i] == 'c' || fflags[i] == 'z') check_for_NaN = TRUE;	/* Do not output records where the crossover or values are NaN */
	}
	if (mixed && GMT_io.binary[GMT_OUT]) {
		fprintf (stderr, "%s: ERROR: Cannot use -Fn with binary output.\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	
	sec_2_unit = gmtdefs.time_system.i_scale;	/* Save conversion from secs to TIME_UNIT before MGD77_Init switches to UNIX time system (seconds) */
	
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
	if (coe_kind == 0) coe_kind = 3;	/* Both */
	
	/* Read the entire data base; note the -I, R and -S options are applied during reading */
	
	if (gmtdefs.verbose) fprintf (stderr, "%s: Read crossover database %s...\n", GMT_program, dbase);
	np = x2sys_read_coe_dbase (s, dbase, ignorelist, wesn, column, coe_kind, track, &P, &nx, &n_tracks);
	if (gmtdefs.verbose) fprintf (stderr, "Found %ld pairs and a total of %ld crossover records.\n", np, nx);

	if (np == 0 && nx == 0) {	/* End here since nothing was allocated */
		x2sys_end (s);
		GMT_end (argc, argv);
		exit (EXIT_SUCCESS);
	}
	
	if (weight_list) {	/* Got file with weights for each track, OR a fixed weight [1] */
		if (x2sys_read_weights (weight_list, &weight_name, &weights, &n_weights) != X2SYS_NOERROR) fixed_weight = atof (weight_list);
	}
	
	/* Must count to see total number of COE per track */
	
	trk_nx = (int *) GMT_memory (VNULL, n_tracks, sizeof (int), GMT_program);
	trk_name = (char **) GMT_memory (VNULL, n_tracks, sizeof (char *), GMT_program);
	for (p = 0; p < np; p++) {	/* For each pair of tracks that generated crossovers */
		for (k = 0; k < 2; k++) {
			trk_nx[P[p].id[k]] += P[p].nx;
			trk_name[P[p].id[k]] = P[p].trk[k];
		}
	}
	
	/* Initialize column output types */
	
	one = 0;	two = 1;	/* Normal track order */
	if (!both)			/* It may already been set by -S+ */
		both = !single;		/* Two columns for many output choices */
	n_out = 1 + (int)both;		/* Number of column output for some cols if single is not specified */

	GMT_io.out_col_type[GMT_X] = (s->geographic) ? GMT_IS_LON : GMT_IS_FLOAT;
	GMT_io.out_col_type[GMT_Y] = (s->geographic) ? GMT_IS_LAT : GMT_IS_FLOAT;
	GMT_io.out_col_type[GMT_T] = GMT_IS_ABSTIME;
	
	for (i = j = 0; !mixed && i < n_items; i++, j++) {	/* Overwrite the above settings */
		switch (fflags[i]) {	/* acdhintTvxyz */
			case 'a':	/* Angle between tracks */
				GMT_io.out_col_type[j] = GMT_IS_FLOAT;
				break;
			case 'c':	/* Crossover value */
				GMT_io.out_col_type[j] = GMT_IS_FLOAT;
				z_is_requested = TRUE;
				break;
			case 'd':	/* Distance along track */
				GMT_io.out_col_type[j] = GMT_IS_FLOAT;
				if (both) GMT_io.out_col_type[++j] = GMT_IS_FLOAT;
				break;
			case 'h':	/* Heading along track */
				GMT_io.out_col_type[j] = GMT_IS_FLOAT;
				if (both) GMT_io.out_col_type[++j] = GMT_IS_FLOAT;
				break;
			case 'I':	/* Time interval (unsigned) */
			case 'i':	/* Time interval (signed) */
				GMT_io.out_col_type[j] = GMT_IS_FLOAT;
				break;
			case 'n':	/* Names of the track(s) [need this case to fall through] */
			case 'N':	/* ID numbers of tracks */
				GMT_io.out_col_type[j] = GMT_IS_FLOAT;
				if (both) GMT_io.out_col_type[++j] = GMT_IS_FLOAT;
				break;
			case 't':	/* Time along track */
				GMT_io.out_col_type[j] = GMT_IS_ABSTIME;
				if (both) GMT_io.out_col_type[++j] = GMT_IS_ABSTIME;
				break;
			case 'T':	/* Time along track since beginning of the first year of the track */
				GMT_io.out_col_type[j] = GMT_IS_FLOAT;
				if (both) GMT_io.out_col_type[++j] = GMT_IS_FLOAT;
				break;
			case 'v':	/* Speed along track */
				GMT_io.out_col_type[j] = GMT_IS_FLOAT;
				if (both) GMT_io.out_col_type[++j] = GMT_IS_FLOAT;
				break;
			case 'w':	/* Crossover composite weight */
				GMT_io.out_col_type[j] = GMT_IS_FLOAT;
				break;
			case 'x':	/* x coordinate of crossover */
				GMT_io.out_col_type[j] = (s->geographic) ? GMT_IS_LON : GMT_IS_FLOAT;
				break;
			case 'y':	/* y coordinate of crossover */
				GMT_io.out_col_type[j] = (s->geographic) ? GMT_IS_LAT : GMT_IS_FLOAT;
				break;
			case 'z':	/* Observed value along track */
				GMT_io.out_col_type[j] = GMT_IS_FLOAT;
				if (both) GMT_io.out_col_type[++j] = GMT_IS_FLOAT;
				z_is_requested = TRUE;
				break;
		}
	}

	if (apply_corrections && !check_for_NaN) {	/* Correction table would not be needed for output */
		fprintf (stderr, "%s: Warning: Correction table not needed for chosen output (corrections ignored).\n", GMT_program);
		apply_corrections = FALSE;
	}

	if (apply_corrections) {	/* Load an ephemeral correction table */
		x2sys_get_corrtable (s, correction_table, n_tracks, trk_name, column, NULL, NULL, &CORR);
	}
	
	if (symm_check) {
		int *x_side[2], half;
		double mid[2];
		for (j = 0; j < 2; j++) x_side[j] = (int *) GMT_memory (VNULL, n_tracks, sizeof (int), GMT_program);
		 trk_symm = (double *) GMT_memory (VNULL, n_tracks, sizeof (double), GMT_program);
		for (p = 0; p < np; p++) {	/* For each pair of tracks that generated crossovers */
			for (j = 0; j < 2; j++) {	/* Set mid-point for these two tracks */
				if (GMT_is_dnan (P[p].start[j]))
					mid[j] = P[p].dist[j];
				else
					mid[j] = 0.5 * (P[p].stop[j] + P[p].start[j]);
			}
			for (k = 0; k < P[p].nx; k++) {
				for (j = 0; j < 2; j++) {
					if (GMT_is_dnan (P[p].COE[k].data[j][COE_T]))
						half = (P[p].COE[k].data[j][COE_D] > mid[j]);
					else
						half = (P[p].COE[k].data[j][COE_T] > mid[j]);
					x_side[half][P[p].id[j]]++;
				}
			}
		}
		/* Compute symmetry */
		for (k = 0; k < n_tracks; k++)  trk_symm[k] = (x_side[1][k] - x_side[0][k]) / (x_side[1][k] + x_side[0][k]);
		for (j = 0; j < 2; j++) GMT_free ((void *)x_side[j]);
	}
	/* Time to issue output */
	
	if (!GMT_io.binary[GMT_OUT]) {	/* Write 3 header records */
		sprintf (buffer, "# Tag: %s %s\n", TAG, column);	GMT_fputs (buffer, GMT_stdout);
		sprintf (buffer, "# Command: %s ", GMT_program);	GMT_fputs (buffer, GMT_stdout);
		if (!dbase) GMT_fputs (" [stdin]", GMT_stdout);
		for (k = 1; k < argc; k++) {GMT_fputs (" ", GMT_stdout);	GMT_fputs (argv[k], GMT_stdout);}
		sprintf (buffer, "\n#");
		GMT_fputs (buffer, GMT_stdout);
		for (i = j = 0; i < n_items; i++, j++) {	/* Overwrite the above settings */
			if (i > 0) GMT_fputs ("\t", GMT_stdout);
			switch (fflags[i]) {	/* acdhintTvxyz */
				case 'a':	/* Angle between tracks */
					GMT_fputs ("angle", GMT_stdout);
					break;
				case 'c':	/* Crossover value */
					sprintf (buffer, "%s_x", column);
					GMT_fputs (buffer, GMT_stdout);
					break;
				case 'd':	/* Distance along track */
					(both) ? sprintf (buffer, "dist_1\tdist_2") : sprintf (buffer, "dist");
					GMT_fputs (buffer, GMT_stdout);
					break;
				case 'h':	/* Heading along track */
					(both) ? sprintf (buffer, "head_1\thead_2") : sprintf (buffer, "head");
					GMT_fputs (buffer, GMT_stdout);
					break;
				case 'I':	/* Time interval (unsigned) */
					GMT_fputs ("u_tint", GMT_stdout);
					break;
				case 'i':	/* Time interval (signed) */
					GMT_fputs ("s_tint", GMT_stdout);
					break;
				case 'n':	/* Names of the track(s) [need this case to fall through] */
					(both) ? sprintf (buffer, "track_1\ttrack_2") : sprintf (buffer, "track");
					GMT_fputs (buffer, GMT_stdout);
					break;
				case 'N':	/* ID numbers of tracks */
					(both) ? sprintf (buffer, "ID_1\tID_2") : sprintf (buffer, "ID");
					GMT_fputs (buffer, GMT_stdout);
					break;
				case 't':	/* Time along track */
					(both) ? sprintf (buffer, "t_1\tt_2") : sprintf (buffer, "t");
					GMT_fputs (buffer, GMT_stdout);
					break;
				case 'T':	/* Time along track since beginning of the first year of the track */
					(both) ? sprintf (buffer, "T_1\tT_2\n") : sprintf (buffer, "T");
					GMT_fputs (buffer, GMT_stdout);
					break;
				case 'v':	/* Speed along track */
					(both) ? sprintf (buffer, "vel_1\tvel_2") : sprintf (buffer, "vel");
					GMT_fputs (buffer, GMT_stdout);
					break;
				case 'w':	/* Composite weight of crossover */
					GMT_fputs ("weight", GMT_stdout);
					break;
				case 'x':	/* x coordinate of crossover */
					(s->geographic) ? GMT_fputs ("lon", GMT_stdout) : GMT_fputs ("x", GMT_stdout);
					break;
				case 'y':	/* y coordinate of crossover */
					(s->geographic) ? GMT_fputs ("lat", GMT_stdout) : GMT_fputs ("y", GMT_stdout);
					break;
				case 'z':	/* Observed value along track */
					if (both) {
						sprintf (buffer, "%s_1\t%s_2", column, column);
						GMT_fputs (buffer, GMT_stdout);
					}
					else
						GMT_fputs (column, GMT_stdout);
					break;
			}
		}
		GMT_fputs ("\n", GMT_stdout);
	}
	
	for (p = 0; p < np; p++) {	/* For each pair of tracks that generated crossovers */
		if (nx_check && (trk_nx[P[p].id[0]] < nx_min || trk_nx[P[p].id[1]] < nx_min)) continue;			/* Not enough COEs */
		if (symm_check && (fabs( trk_symm[P[p].id[0]]) > asymm_max || fabs( trk_symm[P[p].id[1]]) > asymm_max)) continue;	/* COEs not distributed symmatrically */
		np_use++;
		nx_use += P[p].nx;
		if (GMT_io.multi_segments[GMT_OUT]) {
			sprintf (buffer, "%c %s - %s nx = %d\n", GMT_io.EOF_flag[GMT_OUT], P[p].trk[0], P[p].trk[1], P[p].nx);
			GMT_fputs (buffer, GMT_stdout);
		}
		if (gmtdefs.verbose) fprintf (stderr, "%s: Crossovers from %s minus %s [%d].\n", GMT_program, P[p].trk[0], P[p].trk[1], P[p].nx);
		if (single) {	/* May have to flip which is track one and two */
			two = !strcmp (track, P[p].trk[0]);
			one = 1 - two;
		}
		
		for (k = 0; k < P[p].nx; k++) {	/* For each crossover */
			/* First check if this record has a non-NaN COE */
			if (check_for_NaN && (GMT_is_dnan (P[p].COE[k].data[one][COE_Z]) || GMT_is_dnan (P[p].COE[k].data[two][COE_Z]))) continue;
			
			for (i = j = 0, first = TRUE; i < n_items; i++, j++) {
				switch (fflags[i]) {	/* acdhintTvwxyz */
					case 'a':	/* Angle between tracks */
						val[0] = fabs (P[p].COE[k].data[0][COE_H] - P[p].COE[k].data[1][COE_H]);
						while (val[0] >= 180.0) val[0] -= 180.0;
						out[j] = (val[0] > 90.0) ? 180.0 - val[0] : val[0];
						if (mixed) dump_ascii_cols (val, GMT_Z, 1, first);
						break;
					case 'c':	/* Crossover value */
					 	if (apply_corrections) {
							corr[one] = MGD77_Correction_Rec (CORR[P[p].id[one]][COE_Z].term, P[p].COE[k].data[one], NULL);
							corr[two] = MGD77_Correction_Rec (CORR[P[p].id[two]][COE_Z].term, P[p].COE[k].data[two], NULL);
						}
						out[j] = val[0] = (P[p].COE[k].data[one][COE_Z] - corr[one]) - (P[p].COE[k].data[two][COE_Z] - corr[two]);
						if (mixed) dump_ascii_cols (val, GMT_Z, 1, first);
						break;
					case 'd':	/* Distance along track */
						out[j] = val[0] = P[p].COE[k].data[one][COE_D];
						if (both) out[++j] = val[1] = P[p].COE[k].data[two][COE_D];
						if (mixed) dump_ascii_cols (val, GMT_Z, n_out, first);
						break;
					case 'h':	/* Heading along track */
						out[j] = val[0] = P[p].COE[k].data[one][COE_H];
						if (both) out[++j] = val[1] = P[p].COE[k].data[two][COE_H];
						if (mixed) dump_ascii_cols (val, GMT_Z, n_out, first);
						break;
					case 'i':	/* Time interval in current TIME_UNIT */
						out[j] = val[0] = sec_2_unit * (P[p].COE[k].data[one][COE_T] - P[p].COE[k].data[two][COE_T]);
						if (mixed) dump_ascii_cols (val, GMT_Z, 1, first);
						break;
					case 'I':	/* Time interval in current TIME_UNIT */
						out[j] = val[0] = sec_2_unit * fabs (P[p].COE[k].data[one][COE_T] - P[p].COE[k].data[two][COE_T]);
						if (mixed) dump_ascii_cols (val, GMT_Z, 1, first);
						break;
					case 'n':	/* Names of the track(s) */
						if (both) {
							if (!first) GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
							GMT_fputs (P[p].trk[0], GMT_stdout);
							GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
							GMT_fputs (P[p].trk[1], GMT_stdout);
						}
						else {
							if (!first) GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
							if (!track)	/* Can we have a single and not track??? */
								GMT_fputs (P[p].trk[0], GMT_stdout);
							else if (strcmp(track, P[p].trk[0]))	/* Print name of other cruise */
								GMT_fputs (P[p].trk[0], GMT_stdout);
							else
								GMT_fputs (P[p].trk[1], GMT_stdout);
						}
						break;
					case 'N':	/* ID numbers of tracks */
						out[j] = val[0] = (double)P[p].id[one];
						if (both) out[++j] = val[1] = (double)P[p].id[two];
						if (mixed) dump_ascii_cols (val, GMT_Z, n_out, first);
						break;
					case 't':	/* Time along track */
						out[j] = val[0] = P[p].COE[k].data[one][COE_T];
						if (both) out[++j] = val[1] = P[p].COE[k].data[two][COE_T];
						if (mixed) dump_ascii_cols (val, GMT_T, n_out, first);
						break;
					case 'T':	/* Time along track since beginning of the track */
						out[j] = val[0] = sec_2_unit * (P[p].COE[k].data[one][COE_T] - P[p].start[one]);
						if (both) out[++j] = val[1] = sec_2_unit * (P[p].COE[k].data[two][COE_T] - P[p].start[two]);
						if (mixed) dump_ascii_cols (val, GMT_Z, n_out, first);
						break;
					case 'v':	/* Speed along track */
						out[j] = val[0] = P[p].COE[k].data[one][COE_V];
						if (both) out[++j] = val[1] = P[p].COE[k].data[two][COE_V];
						if (mixed) dump_ascii_cols (val, GMT_Z, n_out, first);
						break;
					case 'w':	/* Weight for this crossover */
						if (weights) {	/* Weightfile was given; compute composite weight for this COE */
							for (m = 0, w_k = 0.0; m < 2; m++) {
								if ((id = x2sys_find_track (P[p].trk[m], weight_name, n_weights)) == -1) {
									fprintf (stderr, "%s: No weights found for track %s - using weight = 1.\n", GMT_program, P[p].trk[m]);
									w = 1.0;
								}
								else
									w = weights[id];
								w_k += 1.0 / (w*w);
							}
							out[j] = sqrt (2.0/w_k);
						}
						else
							out[j] = fixed_weight;
						if (mixed) {val[0] = out[j];	dump_ascii_cols (val, GMT_Z, 1, first);}
						break;
					case 'x':	/* x coordinate of crossover */
						out[j] = val[0] = P[p].COE[k].data[0][COE_X];
						if (mixed) dump_ascii_cols (val, GMT_X, 1, first);
						break;
					case 'y':	/* y coordinate of crossover */
						out[j] = val[0] = P[p].COE[k].data[0][COE_Y];
						if (mixed) dump_ascii_cols (val, GMT_Y, 1, first);
						break;
					case 'z':	/* Observed value along track */
						if (apply_corrections) corr[one] = MGD77_Correction_Rec (CORR[P[p].id[one]][COE_Z].term, P[p].COE[k].data[one], NULL);
						out[j] = val[0] = P[p].COE[k].data[one][COE_Z] - corr[one];
						if (both) {
							if (apply_corrections) corr[two] = MGD77_Correction_Rec (CORR[P[p].id[two]][COE_Z].term, P[p].COE[k].data[two], NULL);
							out[++j] = val[1] = P[p].COE[k].data[two][COE_Z] - corr[two];
						}
						if (mixed) dump_ascii_cols (val, GMT_Z, n_out, first);
						break;
				}
				first = FALSE;
			}
			if (mixed)
				GMT_fputs ("\n", GMT_stdout);
			else
				GMT_output (GMT_stdout, j, out);
		}
	}
	if (gmtdefs.verbose) fprintf (stderr, "%s: Output %ld pairs and a total of %ld crossover records.\n", GMT_program, np_use, nx_use);
	
	/* Done, free up data base array */
	
	x2sys_free_coe_dbase (P, np);
	GMT_free ((void *)trk_nx);
	if (symm_check) GMT_free ((void *) trk_symm);

	if (apply_corrections) MGD77_Free_Correction (CORR, n_tracks);
	GMT_free ((void *)trk_name);
	x2sys_end (s);

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

void dump_ascii_cols (double *val, int col, int n, GMT_LONG first)
{	/* Short-hand to dump n = 1 or 2 numerical values in chosen format.
	 * col is used to set the format, and first is TRUE for first item per record.
	 */
	int i;
	for (i = 0; i < n; i++) {
		if (!first) GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
		GMT_ascii_output_one (GMT_stdout, val[i], col);
		first = FALSE;
	}
}
