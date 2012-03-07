/*-----------------------------------------------------------------
 *	$Id: x2sys_put.c,v 1.47 2011/07/11 19:22:07 guru Exp $
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
/* x2sys_put will read one track bin file (.tbf) (from x2sys_binlist)
 * and update the track index database. The tbf file has the index
 * numbers of all the x by x degree boxes traversed by each data set.
 * The bin dimension and extent are determined by the tbf header.
 * This bin info is added to the track index database. This database
 * is then used by program x2sys_get that locates all the data track
 * within a given area that optionally contain specific data columns.
 *
 * Author:	Paul Wessel
 * Date:	14-JUN-2004
 * Version:	1.1, based on the spirit of the old mgg code
 *
 */

#include "x2sys.h"

int main (int argc, char **argv)
{
	char *TAG = CNULL;

	struct X2SYS_INFO *s = NULL;
	struct X2SYS_BIX B;

	struct X2SYS_BIX_TRACK_INFO *this_info = NULL, *new_info = NULL;
	struct X2SYS_BIX_TRACK *this_track = NULL;

	char track[GMT_TEXT_LEN], line[BUFSIZ], *c_unused = NULL;
	char track_file[BUFSIZ], index_file[BUFSIZ], old_track_file[BUFSIZ], old_index_file[BUFSIZ];
	char track_path[BUFSIZ], index_path[BUFSIZ], old_track_path[BUFSIZ], old_index_path[BUFSIZ];

	GMT_LONG error = FALSE, replace = FALSE, delete = FALSE, found_it, skip, swap = FALSE;

	FILE *fp = NULL, *fbin = NULL, *ftrack = NULL;

	int index, id, bin, free_id, max_flag, flag;
	int i, last_id, bit, total_flag;
	size_t s_unused = 0;

	int x2sys_bix_remove_track (int track_id, struct X2SYS_BIX *B);
	struct X2SYS_BIX_TRACK_INFO * x2sys_bix_find_track (char *track, GMT_LONG *found_it, struct X2SYS_BIX *B);

	argc = (int)GMT_begin (argc, argv);
	
	for (i = (int)strlen(argv[0]); i >= 0 && argv[0][i] != '/'; i--);
	X2SYS_program = &argv[0][i+1];	/* Name without full path */

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				case 'V':
				case '\0':
					error += GMT_parse_common_options (argv[i], NULL, NULL, NULL, NULL);
					break;
				case 'D':	/* Remove all traces of these tracks from the database */
					delete = TRUE;
					break;
				case 'F':	/* Force update of existing tracks if new data are found */
					replace = TRUE;
					break;
				case 'T':
					TAG = &argv[i][2];
					break;
				case 'S':
					swap = TRUE;	/* Swap option for index.b reading [Obsolete but left for backwardness] */
					break;
				default:
					error = TRUE;
					fprintf (stderr, "%s: Unrecognized option -%c\n", GMT_program, argv[i][1]);
					break;
			}
		}
		else if ((fp = GMT_fopen (argv[i], "r")) == NULL) {
			fprintf (stderr, "%s: ERROR: Could not open file %s\n", GMT_program, argv[i]);
			exit (EXIT_FAILURE);
		}
	}

	if (argc == 1 || error || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "x2sys_put %s - Update track index database from track bin file\n\n", X2SYS_VERSION);
		fprintf (stderr, "usage: x2sys_put [<info.tbf>] -T<TAG> [-D] [-F] [-V]\n\n");
		fprintf (stderr, "	<info.tbf> is one track bin file from x2sys_binlist [Default reads stdin]\n");
		fprintf (stderr, "	-T <TAG> is the system tag for this compilation\n");

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf(stderr,"	-D will remove the listed tracks  [Default will add to database]\n");
		fprintf(stderr,"	-F will force updates to earlier entries for a track with new information\n");
		fprintf(stderr,"	   [Default refuses to process tracks already in the database]\n");
#ifdef OBSOLETE
		fprintf (stderr, "	-S Byte swap binary records during read [no swapping]\n");
#endif
		GMT_explain_option ('V');
		exit (EXIT_FAILURE);
	}

	if (delete && replace) {
		fprintf (stderr, "%s: ERROR: Only specify one of -D and -F\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	
	if (replace) delete = TRUE;	/* Ironic, given previous if-test, but that is how the logic below works */
	
	if (fp == NULL) fp = GMT_stdin;	/* No file given; read stdin instead */

	c_unused = GMT_fgets (line, BUFSIZ, fp);	/* Got the first record from the track binindex file */
	if (strncmp (&line[2], TAG, strlen(TAG))) {	/* Hard check to see if the TAG matches what we says it should be */
		fprintf (stderr,"%s: The TAG specified (%s) does not match the one in the .tbf file (%s)\n", GMT_program, TAG, &line[2]);
		exit (EXIT_FAILURE);
	}

	/* Open TAG file and set the operational parameters */

	x2sys_err_fail (x2sys_set_system (TAG, &s, &B, &GMT_io), TAG);

	for (i = max_flag = 0, bit = 1; i < s->n_fields; i++, bit <<= 1) max_flag |= bit;

	x2sys_bix_init (&B, FALSE);

	/* Read existing track-information from <ID>_tracks.d file */

	 x2sys_err_fail (x2sys_bix_read_tracks (s, &B, 0, &last_id), "");
	 last_id--;

	/* Read geographical track-info from <ID>_index.b file */

	x2sys_err_fail (x2sys_bix_read_index (s, &B, swap), "");

	/* Ok, now we can start reading new info */

#ifdef DEBUG
	GMT_memtrack_off (GMT_mem_keeper);
#endif
	c_unused = GMT_fgets (line, BUFSIZ, fp);
	while (line[0] == '>') {	/* Next segment */
		sscanf (line, "> %s", track);
		
		/* Determine if this track is already in the data base */
		
		this_info = x2sys_bix_find_track (track, &found_it, &B);	/* Returns found_it = TRUE if found */
		
		/* In either case, this_info now points to the previous track so that this_info->next is the found track OR
		 * it is the point after which a new track should be inserted */
		
		free_id = 0;	/* Default is to add a new track entry */
		if (found_it) {	/* This track already exists in the database */
			if (delete) {	/* Here we wish to delete it (and possibly replace the contents) */
				if (gmtdefs.verbose) fprintf (stderr, "%s: Removing existing information for track: %s\n", GMT_program, track);
				free_id = x2sys_bix_remove_track (this_info->next_info->track_id, &B);
				fprintf (stderr, "%s: track %s removed\n", GMT_program, track);
				this_info->next_info = this_info->next_info->next_info;
				skip = !replace;	/* If we are not replacing the info then we skip the new info */
			}
			else {	/* Refuse to process tracks already in the database without the delete[and replace] options set */
				if (gmtdefs.verbose) fprintf (stderr, "%s: Track already in database (skipped): %s\n", GMT_program, track);
				skip = TRUE;
			}
		}
		else if (delete) {	/* Here we did not found the track: Give message and go back and read next track information */
			fprintf (stderr, "%s: track %s was not found in the database!\n", GMT_program, track);
			skip = TRUE;
		}
		else	/* Get here when we wish to add a new track not in the database */
			skip = FALSE;

		if (skip) {	/* Just wind past this segment */
			c_unused = GMT_fgets (line, BUFSIZ, fp);
			while (line[0] != '>' && (GMT_fgets (line, BUFSIZ, fp) != NULL));
		}
		else {	/* Read the tbf information for this track */

			if (gmtdefs.verbose) fprintf (stderr, "%s: Adding track: %s\n", GMT_program, track);

			/* If a track is replaced, then use the same id_no, else increment to get a new one */

			id = (free_id) ? free_id : ++last_id;
			if (!free_id) {	/* Must create a new entry */
				new_info = x2sys_bix_make_entry (track, id, 0);
				new_info->next_info = this_info->next_info;
				this_info->next_info = new_info;
				this_info = new_info;
			}
			else	/* Reuse the previously deleted entry */
				this_info = this_info->next_info;

			total_flag = 0;
			while (GMT_fgets (line, BUFSIZ, fp) && line[0] != '>') {
				i = sscanf (line, "%*s %*s %d %d", &index, &flag);
				if (i != 2) {	/* Could not decode the index and the flag entries */
					fprintf (stderr, "%s: Error processing record for track %s [%s]\n", GMT_program, track, line);
					exit (EXIT_FAILURE);
				}
				else if (flag > max_flag) {
					fprintf (stderr, "%s: data flag (%d) exceed maximum (%d) for track %s!\n", GMT_program, flag, max_flag, track);
					exit (EXIT_FAILURE);
				}
				if (B.base[index].n_tracks == 0) {	/* First track to cross this bin */
					B.base[index].first_track = x2sys_bix_make_track (0, 0);
					B.base[index].last_track = B.base[index].first_track;
				}
				B.base[index].last_track->next_track = x2sys_bix_make_track (id, flag);
				B.base[index].last_track = B.base[index].last_track->next_track;
				B.base[index].n_tracks++;
				total_flag |= flag;	/* Accumulate flags for the entire track */
			}
			this_info->flag = total_flag;	/* Store the track flags here */
		}
	}
	if (fp != GMT_stdin) GMT_fclose (fp);
#ifdef DEBUG
	GMT_memtrack_on (GMT_mem_keeper);
#endif

	/* Done, now we must rewrite the <ID>_index.b and <ID>_tracks.d files */

	sprintf (track_file, "%s/%s_tracks.d", TAG, TAG);
	sprintf (index_file, "%s/%s_index.b",  TAG, TAG);
	x2sys_path (track_file, track_path);
	x2sys_path (index_file, index_path);

	sprintf (old_track_file, "%s_old", track_file);
	sprintf (old_index_file, "%s_old", index_file);
	x2sys_path (old_track_file, old_track_path);
	x2sys_path (old_index_file, old_index_path);

	remove (old_track_path);	/* First delete old files */
	if (rename (track_path, old_track_path)) {
		fprintf (stderr, "%s: Rename failed for %s\t%s. Aborting %d!\n", GMT_program, track_path, old_track_path, i);
		exit (EXIT_FAILURE);
	}
	remove (old_index_path);	/* First delete old files */
	if (rename (index_path, old_index_path)) {
		fprintf (stderr, "%s: Rename failed for %s. Aborts!\n", GMT_program, index_path);
		exit (EXIT_FAILURE);
	}

	if ((ftrack = fopen (track_path, "w")) == NULL) {
		fprintf (stderr, "%s: Failed to create %s. Aborts!\n", GMT_program, track_path);
		exit (EXIT_FAILURE);
	}
	if ((fbin = fopen (index_path, "wb")) == NULL) {
		fprintf (stderr, "%s: Failed to create %s. Aborts!\n", GMT_program, index_path);
		exit (EXIT_FAILURE);
	}
	fprintf (ftrack,"# %s\n", TAG);
	for (this_info = B.head->next_info; this_info; this_info = this_info->next_info)
		fprintf (ftrack,"%s %d %d\n",this_info->trackname, this_info->track_id, this_info->flag);

	fclose (ftrack);
	chmod (track_file, (mode_t)S_RDONLY);

	for (bin = 0; bin < B.nm_bin; bin++) {
		if (B.base[bin].n_tracks == 0) continue;

		s_unused = fwrite ((void *)(&bin), (size_t)4, (size_t)1, fbin);
		s_unused = fwrite ((void *)(&B.base[bin].n_tracks), (size_t)4, (size_t)1, fbin);
		for (this_track = B.base[bin].first_track->next_track; this_track; this_track = this_track->next_track) {
			s_unused = fwrite ((void *)(&this_track->track_id), (size_t)4, (size_t)1, fbin);
			s_unused = fwrite ((void *)(&this_track->track_flag), (size_t)4, (size_t)1, fbin);
		}
	}
	fclose (fbin);
	chmod (index_file, (mode_t)S_RDONLY);

	if (gmtdefs.verbose) fprintf (stderr, "%s completed successfully\n", GMT_program);

	x2sys_end (s);

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

int x2sys_bix_remove_track (int track_id, struct X2SYS_BIX *B)
{
	/* Remove all traces of the track with id track_id from structure tree for all bins */

	struct X2SYS_BIX_TRACK *track, *skip_track;
	int bin;

	for (bin = 0; bin < B->nm_bin; bin++) {
		if (B->base[bin].n_tracks == 0) continue;	/* No tracks crossed this bin */

		for (track = B->base[bin].first_track; track->next_track && track->next_track->track_id != track_id; track = track->next_track);	/* Finds the track or end-of-list */

		if (track->next_track) {	/* Ok, found it. Lets remove it from this bin's list */
			skip_track = track->next_track;	/* These 3 lines sets the next points to skip the item to be removed */
			track->next_track = skip_track->next_track;
			skip_track->next_track = NULL;
			B->base[bin].n_tracks--;	/* One less entry for this bin */
			if (!track->next_track) B->base[bin].last_track = track;	/* Update the last track in case we just removed it */
			GMT_free ((void *) skip_track);	/* Remove memory associated with the track to be removed */
			if (B->base[bin].n_tracks == 0) GMT_free ((void *)B->base[bin].first_track);	/* OK, that was the only track in this bin, apparently */
		}
	}
	return (track_id);
}

struct X2SYS_BIX_TRACK_INFO * x2sys_bix_find_track (char *track, GMT_LONG *found_it, struct X2SYS_BIX *B)
{	/* Looks for given track in data base and if found returns pointer to the track before it and sets found_it to TRUE.
	 * I.e., the track is actually this_info->next_info.  If not found set found_it to FALSE and return pointer where
	 * this track should be inserted */
	
	struct X2SYS_BIX_TRACK_INFO *this_info;
	for (this_info = B->head; this_info->next_info && strcmp (this_info->next_info->trackname, track) < 0; this_info = this_info->next_info);
	*found_it = (this_info->next_info != (struct X2SYS_BIX_TRACK_INFO *)NULL && !strcmp (this_info->next_info->trackname, track));
	return (this_info);
}
