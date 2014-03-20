/*--------------------------------------------------------------------
 *	$Id: gmtset.c 10173 2014-01-01 09:52:34Z pwessel $
 *
 *	Copyright (c) 1991-2014 by P. Wessel and W. H. F. Smith
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
 * gmtset will set the specified options to the argument that follows them.
 *
 * Author:	Paul Wessel
 * Date:	10-JUN-2000
 * Version:	4.1.x
 */
 
#include "gmt.h"

int main (int argc, char **argv)
{
	GMT_LONG i, j;
	char *file = CNULL;

	/* SPECIAL INITIALIZATION SINCE BMT_begin IS NOT USED HERE !! */
#ifdef DEBUG
	GMT_memtrack_init (&GMT_mem_keeper);
#endif
	GMT_set_home ();
	GMT_init_fonts (&GMT_N_FONTS);
	GMT_begin_io ();
	if (argc == 2 && argv[1][0] == '-' && argv[1][1] == 0) GMT_give_synopsis_and_exit = TRUE;

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "gmtset %s - To set individual default parameters\n\n", GMT_VERSION);
		fprintf (stderr, "usage: gmtset [-G<defaultsfile>] PARAMETER1 [=] value1 PARAMETER2 [=] value2 PARAMETER3 [=] value3 ...\n");
		fprintf (stderr, "\tFor available PARAMETERS, see gmtdefaults man page\n");

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-G sets name of specific .gmtdefaults4 file to modify\n");
		fprintf (stderr, "\t   [Default looks for file in current directory.  If not found,\n");
		fprintf (stderr, "\t   it looks in the home directory, if not found it uses GMT defaults.\n");
		fprintf (stderr, "\t   The modified defaults are written to the current directory].\n");

		exit (EXIT_FAILURE);
	}

	for (i = strlen(argv[0]); i >= 0 && argv[0][i] != '/'; i--);
	GMT_program = &argv[0][i+1];	/* Name without full path */

	for (i = 1, j = 0; i < argc && j == 0; i++) {
		if (!strncmp (argv[i], "-G", (size_t)2)) {
			file = &argv[i][2];
			j = i;
		}
	}

	if (j) {
		for (i = j + 1; i < argc; i++, j++) argv[j] = argv[i];	/* Remove the -G string */
		argc--;
	}

	GMT_getdefaults (file);

	GMT_setdefaults (argc, argv);

	GMT_putdefaults (file);

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}
