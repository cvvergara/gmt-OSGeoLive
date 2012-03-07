/*--------------------------------------------------------------------
 *	$Id: gmtdefaults.c,v 1.45 2011/03/03 21:02:51 guru Exp $
 *
 *	Copyright (c) 1991-2011 by P. Wessel and W. H. F. Smith
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
 * gmtdefaults will list the users default settings for the GMT-SYSTEM or
 * (by using the -D option), get the GMT-SYSTEM's default settings.
 *
 * Author:	Paul Wessel
 * Date:	12-JUL-2000
 * Version:	4.0
 */
 
#define GMT_WITH_NO_PS
#include "gmt.h"

int main (int argc, char **argv)
{
	GMT_LONG i, get = 0;

	GMT_LONG error = FALSE, get_sys_defaults = FALSE, get_user_defaults = FALSE;

	char *path = NULL;

	/* SPECIAL INITIALIZATION SINCE BMT_begin IS NOT USED HERE !! */
#ifdef DEBUG
	GMT_memtrack_init (&GMT_mem_keeper);
#endif
	for (i = strlen(argv[0]); i >= 0 && argv[0][i] != '/'; i--);
	GMT_program = &argv[0][i+1];	/* Name without full path */
	GMT_set_home ();
	GMT_init_fonts (&GMT_N_FONTS);

	GMT_begin_io ();

	for (i = 1; !error && i < argc; i++) {
		if (argv[i][0] != '-') continue;
		switch (argv[i][1]) {
			case '\0':
				error += GMT_parse_common_options (argv[i], 0, 0, 0, 0);
				break;
			case 'D':	/* Get GMT defaults settings */
				get_sys_defaults = TRUE;
				switch (argv[i][2]) {
					case 'S':	/* SI version */
					case 's':
						get = 1;
						break;
					case 'U':	/* US version */
					case 'u':
						get = 2;
						break;
					default:	/* Version chosen in gmt.conf */
						get = 0;
						break;
				}
				break;
			case 'L':	/* List the user's current GMT defaults settings */
				get_user_defaults = TRUE;
				break;
			default:
				error = TRUE;
                                GMT_default_error (argv[i][1]);
                                exit (EXIT_FAILURE);
				break;
		}
	}

	if (argc == 1 || GMT_give_synopsis_and_exit) {
		fprintf (stderr, "gmtdefaults %s - List GMT default parameters to stdout\n\n", GMT_VERSION);
		fprintf (stderr, "usage: gmtdefaults [-D[s|u] | -L]\n\n");
		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);
		fprintf (stderr, "\t-D prints the default settings for the GMT system\n");
		fprintf (stderr, "\t   Append s to see the SI version of defaults\n");
		fprintf (stderr, "\t   Append u to see the US version of defaults\n");
		fprintf (stderr, "\t-L prints the users current GMT default settings\n");
		exit (EXIT_FAILURE);
	}

	if ((get_user_defaults + get_sys_defaults) != 1) {
		fprintf (stderr, "%s: GMT SYNTAX ERROR: Must specify one of -D or -L\n", GMT_program);
		error++;
	}

	if (error) exit (EXIT_FAILURE);

	if (get_user_defaults) {
		GMT_getdefaults (CNULL);
	}
	else {
		GMT_getdefpath (get, &path);
		GMT_getdefaults (path);
#ifdef WIN32
		/* We have a "race" condition here. Plain free() crashes Windows but *nix doesn't like GMT_free */
		GMT_free ((void *)path);
#else
		free ((void *)path);
#endif
	}

	GMT_putdefaults ("-");

	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}
