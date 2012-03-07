/*-----------------------------------------------------------------
 *	$Id: x2sys_merge.c,v 1.12 2011/07/11 19:22:07 guru Exp $
 *
 *      Copyright (c) 1999-2011 by J. Luis
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
/* x2sys_merge will read two crossovers data base and output the contents
 * of the main one updated with the COEs in the second one. The second
 * file should only contain updated COEs relatively to the first one.
 * That is, it MUST NOT contain any new two tracks intersections. 
 *
 * Author:	Joaquim Luis
 * Date:	09-Juin-2009
 *
 */

#include "gmt.h"

int main (int argc, char **argv) {

	int  i, j, k, n_alloc, n_base, n_merge, *map_base_start = NULL, *map_base_end = NULL, *map_merge_start = NULL, *map_merge_end = NULL;
	int  merge_start;
	char line[BUFSIZ], *dbase1 = CNULL, *dbase2 = CNULL, **pairs_base = NULL, **pairs_merge = NULL, *c_not_used = NULL;
	FILE *fp_base = NULL, *fp_merge = NULL;

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				case 'A':
					dbase1 = &argv[i][2];
					break;
				case 'M':
					dbase2 = &argv[i][2];
					break;
				default:
					GMT_give_synopsis_and_exit = 1;
					break;
			}
		}
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "x2sys_merge - Merge an updated COEs table (smaller) into a main one (bigger)\n\n");
		fprintf (stderr, "usage: x2sys_merge -A<main_COEdbase> -M<new_COEdbase>\n");

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\t-A<main_COEdbase> File with the main crossover error data base.\n");
		fprintf (stderr, "\t-M<new_COEdbase> File with the new crossover error data base.\n");
		fprintf (stderr, "\t  The new COEs will replace the old ones present in <main_COEdbase>.\n\n");
		fprintf (stderr, "\t  Result is printed to stdout.\n");
		exit (EXIT_FAILURE);
	}

	if (!dbase1) {
		fprintf (stderr, "%s: ERROR: Missing Base COEs database file. -A is mandatory\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	if (!dbase2) {
		fprintf (stderr, "%s: ERROR: Missing Updating COEs database file. -M is mandatory\n", GMT_program);
		exit (EXIT_FAILURE);
	}

	if ((fp_base = fopen (dbase1, "r")) == NULL) {
		fprintf (stderr, "%s: ERROR: Unable to open crossover file %s\n", GMT_program, dbase1);
		exit (EXIT_FAILURE);
	}

	if ((fp_merge = fopen (dbase2, "r")) == NULL) {
		fprintf (stderr, "%s: ERROR: Unable to open crossover file %s\n", GMT_program, dbase2);
		exit (EXIT_FAILURE);
	}

	n_alloc = GMT_CHUNK;
	map_base_start = (int *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (int), GMT_program);
	map_base_end = (int *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (int), GMT_program);
	pairs_base = (char **) GMT_memory (VNULL, (size_t)n_alloc, sizeof (char *), GMT_program);

	map_merge_start = (int *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (int), GMT_program);
	map_merge_end = (int *) GMT_memory (VNULL, (size_t)n_alloc, sizeof (int), GMT_program);
	pairs_merge = (char **) GMT_memory (VNULL, (size_t)n_alloc, sizeof (char *), GMT_program);

	/* Read in the main COEs dbase and store the pair track names */
	n_base = 0;		k = 1;
	while (fgets (line, BUFSIZ, fp_base)) {	/*  */
		if (line[0] == '>') {
			map_base_start[n_base] = k;
			if (n_base) map_base_end[n_base-1] = k - 1;
			pairs_base[n_base] = (char *) GMT_memory (VNULL, (size_t)24, sizeof (char), GMT_program);
			strncpy(pairs_base[n_base], &line[2], 19);
			n_base++;
			if (n_base == n_alloc) {
				n_alloc <<= 1;
				map_base_start = (int *) GMT_memory ((void *)map_base_start, (size_t)n_alloc, sizeof (int), GMT_program);
				map_base_end = (int *) GMT_memory ((void *)map_base_end, (size_t)n_alloc, sizeof (int), GMT_program);
				pairs_base = (char **) GMT_memory ((void *)pairs_base, (size_t)n_alloc, sizeof (char *), GMT_program);
			}
		}

		k++;
	}
	map_base_end[n_base - 1] = k - 1;	/* This one was not yet assigned */
	rewind(fp_base);

	/* Read in the updated COEs dbase and store the pair track names */
	n_alloc = GMT_CHUNK;
	n_merge = 0;		k = 1;
	while (fgets (line, BUFSIZ, fp_merge)) {	/*  */
		if (line[0] == '>') {
			map_merge_start[n_merge] = k;
			if (n_merge) map_merge_end[n_merge-1] = k - 1;
			pairs_merge[n_merge] = (char *) GMT_memory (VNULL, (size_t)24, sizeof (char), GMT_program);
			strncpy(pairs_merge[n_merge], &line[2], 19);
			n_merge++;
			if (n_merge == n_alloc) {
				n_alloc <<= 1;
				map_merge_start = (int *) GMT_memory ((void *)map_merge_start, (size_t)n_alloc, sizeof (int), GMT_program);
				map_merge_end = (int *) GMT_memory ((void *)map_merge_end, (size_t)n_alloc, sizeof (int), GMT_program);
				pairs_merge = (char **) GMT_memory ((void *)pairs_merge, (size_t)n_alloc, sizeof (char *), GMT_program);
			}
		}

		k++;
	}
	map_merge_end[n_merge - 1] = k - 1;	/* This one was not yet assigned */
	rewind(fp_merge);


	/* Jump comment lines in both files and osition the file poiter into the first data line */
	k = 0;
	while (fgets (line, BUFSIZ, fp_merge) && line[0] == '#') k++;	/* Jump the comment lines in the to-merge file */
	rewind(fp_merge);
	for (i = 0; i < k; i++) c_not_used = fgets (line, BUFSIZ, fp_merge);

	k = 0;
	while (fgets (line, BUFSIZ, fp_base) && line[0] == '#') {
		fprintf (stdout, "%s", line);
		k++;
	}
	rewind(fp_base);
	for (i = 0; i < k; i++) c_not_used = fgets (line, BUFSIZ, fp_base);

	/* Do the merging. COEs present in file dbase2 replace their pairs in dbase1 */
	merge_start = 0;
	for (i = 0; i < n_base; i++) {
		for (j = merge_start; j < n_merge; j++) {
			if (!strcmp(pairs_base[i], pairs_merge[j])) {		 /* Update these COEs */
				for (k = map_merge_start[j]; k <= map_merge_end[j]; k++) {
					c_not_used = fgets (line, BUFSIZ, fp_merge);
					fprintf (stdout, "%s", line);
				}
				for (k = map_base_start[i]; k <= map_base_end[i]; k++)	/* Advance also in the base file */
					c_not_used = fgets (line, BUFSIZ, fp_base);

				merge_start = j + 1;
				break;		/* Since we found this to update, no need to continue seeking for another repetition */
			}
			else if (j == (n_merge - 1)) {	/* Not equal. So do not to update, just recopy */
				for (k = map_base_start[i]; k <= map_base_end[i]; k++) {
					c_not_used = fgets (line, BUFSIZ, fp_base);
					fprintf (stdout, "%s", line);
				}
			}
		}
		if (merge_start == n_merge) {	/* Copy the rest of dbase1 file and stop */
			while (fgets (line, BUFSIZ, fp_base))
				fprintf (stdout, "%s", line);
			break;			/* Not very elegant way of stopping, but we are done */
		}
	}

	fclose(fp_base);
	fclose(fp_merge);

	for (i = 0; i < n_base; i++) GMT_free ((void *)pairs_base[i]);
	for (i = 0; i < n_merge; i++) GMT_free ((void *)pairs_merge[i]);
	GMT_free((void *) map_base_start);
	GMT_free((void *) map_base_end);
	GMT_free((void *) map_merge_start);
	GMT_free((void *) map_merge_end);

	exit (EXIT_SUCCESS);
}
