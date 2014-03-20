/*--------------------------------------------------------------------
 *	$Id: binlegs.c 10173 2014-01-01 09:52:34Z pwessel $
 *
 *    Copyright (c) 1991-2014 by P. Wessel and W. H. F. Smith
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/*
 * BINLEGS - A Program that will read a list of *.bix
 * files and update the legindex files. The bix file
 * has the index numbers of all the 1 by 1 degree boxes crossed
 * by this leg. This information is added to the legindex
 * database. This database is then used by program GMTLEGS that
 * will do the job of locating all legs within a given area (same
 * as Brown-book LOC, but now all the legs listed are really
 * inside the area).  The list of legs given to binlegs will have
 * the .bix suffix appended before opening the index file
 *
 * Author:	Paul Wessel
 * Date:	21-MAY-1987
 * Revised:	13-SEP-1989	Added option to remove legs and specify path
 * 		2-JAN-1991	Replaced fscanf with fgets
 * 		13-FEB-1991	Allow for legnames up to 16 chars
 * Version:	2.0 21-Jun-1991
 *		3.1 15-OCT-1998 POSIX compliant
 *		3.2 10-MAR-1999 New GMT_SHAREDIR/mgg directory place for data
 *
 */

#include "gmt.h"

#ifdef WIN32
#define _chmod(path,mode) chmod(path,mode)
extern int _chmod (const char *path, int mode);
#else
#include <sys/stat.h>
#endif

#define S_RDONLY 0000444

struct GRID {
	int bix;
	int n_legs;
	struct LEG *first_leg, *last_leg;
};

struct LEG {
	int leg_id;		/* Low 4 bytes carries gmt-flag, rest has id_no << 4 */
	struct LEG *next_leg;
};

struct INFO {
	char legname[16];
	int leg_id;
	int gmt;
	struct INFO *next_info;
};

int main (int argc, char **argv)
{

	struct INFO *info_head, *this_info, *new_info, *make_info(char *name, int id_no, int flag);
	struct LEG *this_leg, *make_leg(int id, int flag);

	char file[BUFSIZ], leg[16], lname[16], line[BUFSIZ], *c_unused = NULL;
	char leg_file[BUFSIZ], index_file[BUFSIZ], old_leg_file[BUFSIZ], old_index_file[BUFSIZ];
	char error = FALSE, replace = FALSE, delete = FALSE, found_it, verbose = FALSE;

	FILE *fp = NULL, *fbin = NULL, *fleg = NULL, *find = NULL;

	int no_of_legs, index, id, bin, free_id, max, flag, i, j, year, last_id;
	size_t not_used = 0;
	int remove_leg(int leg_id, struct GRID *grid);
	struct GRID grid[64800];

	argc = (int)GMT_begin (argc, argv);

	if (argc == 1) error = TRUE;
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			for (j = 1; argv[i][j]; j++) {
				switch (argv[i][1]) {
					case 'R':	/* Replace existing info for the legs, then append new info */
						replace = TRUE;
						break;
					case 'D':	/* Remove all traces of these legs from the info files */
						delete = TRUE;
						replace = TRUE;
						break;
					case 'V':
						verbose = TRUE;
						break;
					default:
						error = TRUE;
						fprintf (stderr, "SYNTAX ERROR:  Unrecognized option -%c\n", argv[i][1]);
						break;
				}
			}
		}
		else
			fp = fopen(argv[i],"r");
	}

	if (error || fp == NULL) {
		fprintf(stderr, "Usage: binlegs leglistfile [-D] [-R] [-V]\n");
		fprintf(stderr,"	-D will remove info for the listed legs\n");
		fprintf(stderr,"	-R will replace old info. [Default is append].\n");
		fprintf(stderr,"	-V (verbose) will print update-messages. [Default is silent]\n");
		exit (EXIT_FAILURE);
	}

	info_head = make_info ("-", 0, 0);
	this_info = info_head;
	for (bin = 0; bin < 64800; bin++) grid[bin].n_legs = 0;

	GMT_getsharepath ("mgg", "gmt_legs", ".d", leg_file);
	GMT_getsharepath ("mgg", "gmt_index", ".b", index_file);

	if ((fleg = fopen (leg_file,"r")) == NULL) {
		fprintf (stderr,"Could not find %s\n", leg_file);
		exit (EXIT_FAILURE);
	}
	if ((fbin = fopen (index_file,"rb")) == NULL) {
		fprintf (stderr,"Could not open %s\n", index_file);
		exit (EXIT_FAILURE);
	}

	last_id  = -1;

	/* Read existing leg-information from gmt_legs.d file */

	while (fgets (line, BUFSIZ, fleg)) {
		sscanf (line, "%s %d %d",lname, &id, &flag);
		this_info->next_info = make_info (lname, id, flag);
		this_info = this_info->next_info;
		if (id > last_id) last_id = id;
	}
	fclose (fleg);

	/* Read geographical leg-info from gmt_index.b file */

	while ((fread ((void *)(&index), (size_t)4, (size_t)1, fbin)) == 1) {
		not_used = fread ((void *)(&no_of_legs), (size_t)4, (size_t)1, fbin);
		if (grid[index].n_legs == 0) {
			grid[index].first_leg = make_leg (0, 0);
			grid[index].last_leg = grid[index].first_leg;
		}
		for (i = 0; i < no_of_legs; i++) {
			not_used = fread ((void *)(&id), (size_t)4, (size_t)1, fbin);
			flag = id & 15;
			id >>= 4;
			grid[index].last_leg->next_leg = make_leg (id, flag);
			grid[index].last_leg = grid[index].last_leg->next_leg;
			grid[index].n_legs++;
		}
	}
	fclose (fbin);

	/* Ok, now we can start reading new info */

	while (fgets (line, BUFSIZ, fp)) {
		sscanf (line, "%s", leg);
		for (this_info = info_head; this_info->next_info &&
			strcmp (this_info->next_info->legname, leg) < 0;
			this_info = this_info->next_info);
		free_id = 0;
		found_it = (this_info->next_info != (struct INFO *)NULL);
		if (found_it && !strcmp (this_info->next_info->legname, leg)) {
			if (replace) {
				if (verbose) fprintf (stderr, "binlegs: Removing old info for leg: %s\n", leg);
				free_id = remove_leg (this_info->next_info->leg_id, grid);
				if (delete) this_info->next_info = this_info->next_info->next_info;
			}
			else {
				if (verbose) fprintf (stderr, "binlegs: Leg %s already in system! (Skipped)\n",leg);
				continue;
			}
		}
		if (delete ) {	/* Give message and go back and read next leg */
			if (found_it)
				fprintf (stderr, "binlegs: Leg %s removed\n", leg);
			else
				fprintf (stderr, "binlegs: Leg %s not in system!\n", leg);
			continue;
		}

		/* Open the bix file */

		sprintf (file, "%s.bix", leg);
		if ((find = fopen (file,"r")) == NULL) {
			fprintf (stderr, "Could not find file: %s\n",file);
			continue;
		}
		if (verbose) fprintf (stderr, "binlegs: Adding leg: %s\n",leg);
		c_unused = fgets (line, BUFSIZ, find);
		sscanf (line, "%d %d", &year, &max);

		/* If a leg is replaced, then use the same id_no, else increment to get a new one */

		id = (free_id) ? free_id : ++last_id;
		if (!free_id) {
			new_info = make_info (leg, id, 0);
			new_info->next_info = this_info->next_info;
			this_info->next_info = new_info;
			this_info = new_info;
		}
		else
			this_info = this_info->next_info;
		while (fgets (line, BUFSIZ, find)) {
			sscanf (line, "%d %d", &index, &flag);
			if (flag > 7) {
				fprintf(stderr,"binlegs: gmt-flag no good! (Leg=%s)\n", leg);
				exit (EXIT_FAILURE);
			}
			if (grid[index].n_legs == 0) {
				grid[index].first_leg = make_leg (0, 0);
				grid[index].last_leg = grid[index].first_leg;
			}
			grid[index].last_leg->next_leg = make_leg (id, flag);
			grid[index].last_leg = grid[index].last_leg->next_leg;
			grid[index].n_legs++;
		}
		fclose (find);
		this_info->gmt = max;
	}
	fclose (fp);

	/* Done, now we must rewrite the gmt_index.b and gmt_legs.d files */

	sprintf (old_leg_file, "%s_old", leg_file);
	sprintf (old_index_file, "%s_old", index_file);

	if (rename (leg_file, old_leg_file)) {
		fprintf (stderr, "binlegs: Rename failed for %s. Aborts!\n", leg_file);
		exit (EXIT_FAILURE);
	}
	if (rename (index_file, old_index_file)) {
		fprintf (stderr, "binlegs: Rename failed for %s. Aborts!\n", index_file);
		exit (EXIT_FAILURE);
	}

	if ((fleg = fopen (leg_file, "w")) == NULL) {
		fprintf (stderr, "binlegs: Failed to create %s. Aborts!\n", leg_file);
		exit (EXIT_FAILURE);
	}
	if ((fbin = fopen (index_file, "wb")) == NULL) {
		fprintf (stderr, "binlegs: Failed to create %s. Aborts!\n", index_file);
		exit (EXIT_FAILURE);
	}
	for (this_info = info_head->next_info; this_info;
		this_info = this_info->next_info)
		fprintf (fleg,"%s %d %d\n",this_info->legname,
					  this_info->leg_id,
					  this_info->gmt);
	fclose (fleg);
	chmod (leg_file, (mode_t)S_RDONLY);

	for (bin = 0; bin < 64800; bin++) {
		if (grid[bin].n_legs > 0) {
			not_used = fwrite ((void *)(&bin), (size_t)4, (size_t)1, fbin);
			not_used = fwrite ((void *)(&grid[bin].n_legs), (size_t)4, (size_t)1, fbin);
			for (this_leg = grid[bin].first_leg->next_leg; this_leg; this_leg = this_leg->next_leg) not_used = fwrite ((void *)(&this_leg->leg_id), (size_t)4, (size_t)1, fbin);
		}
	}
	fclose (fbin);
	chmod (index_file, (mode_t)S_RDONLY);
	if (verbose) fprintf (stderr, "binlegs completed successfully\n");

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}

struct INFO *make_info (char *name, int id_no, int flag)
{
	struct INFO *new;
	new = (struct INFO *) GMT_memory (VNULL, (size_t)1, sizeof (struct INFO), "binlegs");
	strcpy (new->legname, name);
	new->leg_id = id_no;
	new->gmt = flag;
	new->next_info = 0;
	return (new);
}

struct LEG *make_leg (int id, int flag)
{
	struct LEG *new;
	new = (struct LEG *) GMT_memory (VNULL, (size_t)1, sizeof (struct LEG), "binlegs");
	new->leg_id = (id << 4) + flag;
	new->next_leg = 0;
	return (new);
}

int remove_leg (int leg_id, struct GRID *grid)
{
	/* Remove all traces of the leg with id leg_id from structures */
	struct LEG *leg, *skip_leg;
	int bin;

	for (bin = 0; bin < 64800; bin++) {
		if (grid[bin].n_legs == 0) continue;	/* No legs crossed this bin */

		for (leg = grid[bin].first_leg; leg->next_leg && leg->next_leg->leg_id == leg_id; leg = leg->next_leg);	/* Finds the leg or end-of-list */

		if (leg->next_leg) {	/* Ok, found it. Lets remove it from the list */
			skip_leg = leg->next_leg;
			leg->next_leg = skip_leg->next_leg;
			skip_leg->next_leg = 0;
			grid[bin].n_legs--;
			if (leg->next_leg == NULL) grid[bin].last_leg = leg;
			GMT_free ((void *) skip_leg);
			if (grid[bin].n_legs == 0) GMT_free ((void *)grid[bin].first_leg);
		}
	}
	return (leg_id);
}
