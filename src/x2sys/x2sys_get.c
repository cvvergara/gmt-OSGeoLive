/*-----------------------------------------------------------------
 *	$Id: x2sys_get.c 10173 2014-01-01 09:52:34Z pwessel $
 *
 *      Copyright (c) 1999-2014 by P. Wessel
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
/* x2sys_get will read the track index database and report all the tracks
 * that matches the specified geographical or data-type criteria given
 * on the command line.
 *
 * Author:	Paul Wessel
 * Date:	14-JUN-2004
 * Version:	1.1, based on the spirit of the old mgg code
 *		31-MAR-2006: Changed -X to -L to avoid GMT -X clash
 *		06-DEC-2007: -L did not honor -F -N settings
 *
 */

#include "x2sys.h"

int find_leg (char *name, struct X2SYS_BIX *B, int n)
{	/* Return track id # for this leg */
	int i;

	for (i = 0; i < n; i++) if (!strcmp (name, B->head[i].trackname)) return (i);
	return (-1);
}

int main (int argc, char **argv)
{
	char *TAG = NULL, *fflags = NULL, *nflags = NULL, *fx = NULL;
	char *y_match = NULL, *n_match = NULL, line[BUFSIZ], *p = NULL;
	
	double x, y;

	struct X2SYS_INFO *s = NULL;
	struct X2SYS_BIX B;
	struct X2SYS_BIX_TRACK *track = NULL;

	GMT_LONG error = FALSE, no_suffix, center = FALSE, x_setup = FALSE;
	GMT_LONG y_ok, n_ok, first, *include = NULL, swap = FALSE, quick = FALSE;
	GMT_LONG internal = TRUE, external = TRUE, extension = FALSE, dump = FALSE;

	int combo = 0, n_tracks_found, n_tracks, bit, *in_bin_flag = NULL;
	int id1, id2, item, n_flags = 0, n_pairs, missing = 0, *ids_in_bin = NULL;
	GMT_LONG i, j, k, kk, c, start_j, start_i, stop_j, stop_i, ij;
	
	unsigned int *matrix = NULL;

	double west = 0.0, east = 0.0, south = 0.0, north = 0.0;
	
	FILE *fp = NULL, *fpd = NULL;

	argc = (int)GMT_begin (argc, argv);
	
	for (i = strlen(argv[0]); i >= 0 && argv[0][i] != '/'; i--);
	X2SYS_program = &argv[0][i+1];	/* Name without full path */

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				case 'V':
				case 'R':
				case '\0':
					error += GMT_parse_common_options (argv[i], &west, &east, &south, &north);
					break;
				case 'C':
					center = TRUE;
					break;
				case 'D':
					quick = TRUE;	/* Just list track names and not the full report  */
					break;
				case 'E':
					extension = TRUE;
					break;
				case 'F':
					fflags = &argv[i][2];
					break;
				case 'L':
					if (argv[i][2]) fx = &argv[i][2];
					x_setup = TRUE;
					break;
				case 'N':
					nflags = &argv[i][2];
					break;
				case 'T':
					TAG = &argv[i][2];
					break;
				case 'Q':	/* Specify internal or external only */
					if (argv[i][2] == 'e') internal = FALSE;
					if (argv[i][2] == 'i') external = FALSE;
					break;
				case 's':
					swap = TRUE;	/* Undocumented swap option for index.b reading */
					break;
				case 'd':
					dump = TRUE;	/* Undocumented dump option for index.b debugging */
					break;
				default:
					error = TRUE;
					fprintf (stderr, "%s: Unrecognized option -%c\n", GMT_program, argv[i][1]);
					break;
			}
		}
	}

	if (argc == 1 || error || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "x2sys_get %s - Get track listing from track index database\n\n", X2SYS_VERSION);
		fprintf (stderr, "usage: x2sys_get -T<TAG> [-C] [-D] [-E] [-F<fflags>] [-L[list]] [-N<nflags>] [-Qe|i] [%s] [-V]\n\n", GMT_Rgeo_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "	-C reports center of each tile with tracks instead [Default is track files].\n");
		fprintf (stderr, "	-D only reports the track names and not the report on each field.\n");
		fprintf (stderr, "	-E Append the filename extensions in the report [just track name].\n");
		fprintf (stderr, "	-F is comma-separated list of column names that ALL must be present [Default is any field].\n");
		fprintf (stderr, "	-L Setup mode: Return all pairs of cruises that might intersect given.\n");
		fprintf (stderr, "	   the bin distribution.  Optionally, give file with a list of cruises.\n");
		fprintf (stderr, "	   Then, only pairs with at least one cruise from the list is output.\n");
		fprintf (stderr, "	-N is comma-separated list of column names that ALL must be missing.\n");
		fprintf (stderr,"	-Q Append e for external or i for internal crossovers [Default is both].\n");
		fprintf (stderr,"	   Requires -L and affects the pairs that are returned.\n");
		GMT_explain_option ('R');
		fprintf (stderr, "	[Default region is the entire data domain]\n");
		GMT_explain_option ('V');
		exit (EXIT_FAILURE);
	}

	x2sys_err_fail (x2sys_set_system (TAG, &s, &B, &GMT_io), TAG);

	if (s->geographic) {
		GMT_io.out_col_type[0] = GMT_IS_LON;
		GMT_io.out_col_type[1] = GMT_IS_LAT;
		GMT_io.geo.range = s->geodetic;
	}
	else
		GMT_io.out_col_type[0] = GMT_io.out_col_type[1] = GMT_IS_FLOAT;
		
	if (west == east && south == north) {	/* Set default region */
		west = B.x_min;		east = B.x_max;
		south = B.y_min;	north = B.y_max;
	}

	if (fflags) x2sys_err_fail (x2sys_pick_fields (fflags, s), "-F");
	for (i = combo = 0; i < s->n_out_columns; i++) combo |= (1 << s->out_order[i]);

	if (nflags) {	/* Generate bit pattern that matches the columns that should be missing */
		x2sys_err_fail (x2sys_pick_fields (nflags, s), "-N");
		for (i = missing = 0; i < s->n_out_columns; i++) missing |= (1 << s->out_order[i]);
	}
	
	x2sys_bix_init (&B, FALSE);

	/* Read existing track-information from <ID>_tracks.d file */

	x2sys_err_fail (x2sys_bix_read_tracks (s, &B, 1, &n_tracks), "");

	/* Read geographical track-info from <ID>_index.b file */

	x2sys_err_fail (x2sys_bix_read_index (s, &B, swap), "");

	no_suffix = !extension;
	
	if (x_setup) {
		n_flags = (int)ceil (n_tracks / 32.0);
		include = (GMT_LONG *) GMT_memory (VNULL, (size_t)n_tracks, sizeof (GMT_LONG), GMT_program);
		if (fx) {
			if ((fp = fopen (fx, "r")) == NULL) {
				fprintf (stderr, "%s: ERROR: -L unable to open file %s\n", GMT_program, fx);
				exit (EXIT_FAILURE);
			}
			while (fgets (line, BUFSIZ, fp)) {
				GMT_chop (line);	/* Get rid of [CR]LF */
				if (line[0] == '#' || line[0] == '\0') continue;
				if ((p = strchr (line, '.'))) line[(int)(p-line)] = '\0';	/* Remove extension */
				k = find_leg (line, &B, n_tracks);	/* Return track id # for this leg */
				if (k == -1) {
					fprintf (stderr, "%s: Warning: Leg %s not in the data base\n", GMT_program, line);
					continue;
				}
				include[k] = TRUE;
			}
			fclose (fp);
		}
		else {	/* Use all */
			for (k = 0; k < n_tracks; k++) include[k] = TRUE;
		}
		matrix = (unsigned int *) GMT_memory (VNULL, (size_t)(n_tracks * n_flags + n_tracks / 32), sizeof (unsigned int), GMT_program);
		ids_in_bin = (int *) GMT_memory (VNULL, (size_t)n_tracks, sizeof (int), GMT_program);
	}
	else {	/* y_match will be TRUE for tracks that pass the -F test in at least one bin; same for n_match with the -N test */
		y_match = (char *) GMT_memory (VNULL, (size_t)n_tracks, sizeof (char), GMT_program);
		n_match = (char *) GMT_memory (VNULL, (size_t)n_tracks, sizeof (char), GMT_program);
	}
	in_bin_flag = (int *) GMT_memory (VNULL, (size_t)n_tracks, sizeof (int), GMT_program);
	
	/* Ok, now we can start finding the tracks requested */

	x2sys_err_fail (x2sys_bix_get_ij (west, south, &start_i, &start_j, &B, &j), "");
	x2sys_err_fail (x2sys_bix_get_ij (east, north, &stop_i, &stop_j, &B, &j), "");
	if (B.periodic && stop_i < start_i) stop_i += B.nx_bin;	/* Deal with longitude periodicity */

	if (dump) {
		fpd = fopen ("x2sys_index_dump.txt", "w");
		fprintf (fpd, "# Complete dump of x2sys_index.b for this region and tag.\n");
		fprintf (fpd, "# Segment header shows index number, the center location, and number of tracks crossing this bin.\n");
		fprintf (fpd, "# Each track entry shows track name, integer flag, and Y|N for all observables.  The columns are:\n# Track\t\tflag");
		for (c = 0; c < s->n_fields; c++) fprintf (fpd, "\t%s", s->info[c].name);
		fprintf (fpd, "\n#\n");
	}
	for (j = start_j; j <= stop_j; j++) {
		for (i = start_i; i <= stop_i; i++) {
			ij = j * B.nx_bin + (i % B.nx_bin);	/* Since i may exceed nx_bin due to longitude periodicity */
			if (dump) fprintf (fpd, "> Index %ld lon = %g lat = %g n_tracks = %d\n", ij, B.x_min + ((ij % B.nx_bin) + 0.5) * B.bin_x, B.y_min + ((ij / B.nx_bin) + 0.5) * B.bin_y, B.base[ij].n_tracks);
			if (B.base[ij].n_tracks == 0) continue;

			for (k = kk = 0, track = B.base[ij].first_track->next_track, first = TRUE; first && k < B.base[ij].n_tracks; k++, track = track->next_track) {
				in_bin_flag[track->track_id] |= track->track_flag;	/* Build the total bit flag for this cruise INSIDE the region only */
				if (dump) {	/* Write out all info for each track in this bin */
					fprintf (fpd, "%s\t%d", B.head[track->track_id].trackname, track->track_flag);
					for (c = 0, bit = 1; c < s->n_fields; c++, bit <<= 1) {
						(track->track_flag & bit) ? fprintf (fpd, "\tY") : fprintf (fpd, "\tN");
					}
					fprintf (fpd, "\n");
				}
					
				if (x_setup) {	/* Just build integer list of track ids for this bin */
					y_ok = (fflags) ? ((track->track_flag & combo) == combo) : TRUE;	/* Each cruise must have the requested fields */
					n_ok = (nflags) ? ((track->track_flag & missing) == missing) : TRUE;	/* Each cruise must have the missing fields */
					if (y_ok && n_ok) ids_in_bin[kk++] = track->track_id;
				}
				else {
					/* -F is straightforward: If at least one bin has all required cols then we flag the track to be reported */
					y_ok = (fflags) ? ((track->track_flag & combo) == combo && y_match[track->track_id] == 0) : TRUE;
					/* -N is less straightforward: We will skip it if any bin has any of the columns that all should be missing */
					n_ok = (nflags) ? ((track->track_flag & missing) != 0 && n_match[track->track_id] == 0) : FALSE;
					if (n_ok) n_match[track->track_id] = 1;	/* Track FAILED to satisfy the -N criteria (it had the data); we will skip tracks whose n_match is TRUE */
					if (y_ok) y_match[track->track_id] = 1;	/* Track satisfied the -F criteria; include this track */
					if (y_ok && !n_ok && center && first) {
						x = B.x_min + ((ij % B.nx_bin) + 0.5) * B.bin_x;
						y = B.y_min + ((ij / B.nx_bin) + 0.5) * B.bin_y;
						GMT_ascii_output_one (GMT_stdout, x, 0);
						GMT_fputs (gmtdefs.field_delimiter, GMT_stdout);
						GMT_ascii_output_one (GMT_stdout, y, 1);
						GMT_fputs ("\n", GMT_stdout);
						first = FALSE;
					}
				}
			}
			if (x_setup) {	/* Set bits for every possible pair, but exclude pairs not involving legs given */
				for (id1 = 0; id1 < kk; id1++) {
					for (id2 = id1 + 1; id2 < kk; id2++) {	/* Loop over all pairs */
						if (!(include[ids_in_bin[id1]] || include[ids_in_bin[id2]])) continue;	/* At last one leg must be from our list (if given) */
						item = ids_in_bin[id2] / 32;
						bit = ids_in_bin[id2] % 32;
						matrix[ids_in_bin[id1]*n_flags+item] |= (1 << bit);
						item = ids_in_bin[id1] / 32;
						bit = ids_in_bin[id1] % 32;
						matrix[ids_in_bin[id2]*n_flags+item] |= (1 << bit);
					}
				}
			}
				
		}
	}
	if (dump) {
		fclose (fpd);
		fprintf (stderr, "%s: Index debug dump for selected region written to x2sys_index_dump.txt\n", GMT_program);
	}
	if (x_setup) {
		for (id1 = n_pairs = 0; id1 < n_tracks; id1++) {
			if (internal && include[id1]) printf ("%s\t%s\n", B.head[id1].trackname, B.head[id1].trackname);
			if (!external) continue;	/* Skip all external comparisons */
			for (id2 = id1 + 1; id2 < n_tracks; id2++) {
				item = id2 / 32;
				bit = id2 % 32;
				if (!(matrix[id1*n_flags+item] & (1 << bit))) continue;	/* Pair not selected */
				n_pairs++;
				/* OK, print out pair, with lega alphabetically lower than legb */
				if (strcmp (B.head[id1].trackname, B.head[id2].trackname) < 0)
					printf ("%s\t%s\n", B.head[id1].trackname, B.head[id2].trackname);
				else
					printf ("%s\t%s\n", B.head[id2].trackname, B.head[id1].trackname);
			}
		}
		GMT_free ((void *) matrix);
		GMT_free ((void *) include);
		GMT_free ((void *) ids_in_bin);
		if (gmtdefs.verbose ) fprintf (stderr, "%s: Found %d pairs for crossover consideration\n", GMT_program, n_pairs);
	}
	else if (!center) {
		for (k = n_tracks_found = 0; k < n_tracks; k++) if (y_match[k] == 1 && n_match[k] == 0) n_tracks_found++;
		if (n_tracks_found) {
			if (gmtdefs.verbose ) fprintf (stderr, "%s: Found %d tracks\n", GMT_program, n_tracks_found);
	
			if (!quick) {
				printf ("# Search command: %s", GMT_program);
				for (i = 1; i < argc; i++) printf (" %s", argv[i]);
				printf ("\n#track_ID\t");
				for (c = 0; c < (s->n_fields-1); c++) printf ("%s\t", s->info[c].name);
				printf ("%s\n", s->info[s->n_fields-1].name);
			}
			for (k = 0; k < n_tracks; k++) {
				if (y_match[k] == 0 || n_match[k] == 1) continue;	/* Track that failed the -F test or the -N test are skipped */
				if (no_suffix) {
					for (i = strlen (B.head[k].trackname) - 1; i > 0 && B.head[k].trackname[i] != '.'; i--);
					if (i) B.head[k].trackname[i] = '\0';
				}
				printf ("%s", B.head[k].trackname);
				if (extension) printf (".%s", s->suffix);
				if (!quick) {
					for (c = 0, bit = 1; c < s->n_fields; c++, bit <<= 1) {
						(in_bin_flag[k] & bit) ? printf ("\tY") : printf ("\tN");
					}
				}
				printf ("\n");
			}
		}
		else
			if (gmtdefs.verbose ) fprintf (stderr, "%s: Search found no tracks\n", GMT_program);
	}
	
	GMT_free ((void *)y_match);
	GMT_free ((void *)n_match);
	GMT_free ((void *)in_bin_flag);
	x2sys_end (s);

	if (gmtdefs.verbose) fprintf (stderr, "%s completed successfully\n", GMT_program);

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}
