/*--------------------------------------------------------------------
 *	$Id: gmtpath.c 9923 2012-12-18 20:45:53Z pwessel $
 *
 *    Copyright (c) 1991-2013 by P. Wessel and W. H. F. Smith
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/*
 * gmtpath takes legid(s) as argument and returns the full path
 * to where this data file(s) can be found.
 *
 * Paul Wessel
 * 11/13/87
 */
 
#include "gmt.h"
#include "gmt_mgg.h"

int main (int argc, char **argv)
{
	char path[BUFSIZ];
	int i, error = FALSE;

	argc = (int)GMT_begin (argc, argv);
	
	gmtmggpath_init();

	for (i = 1; !error && i < argc; i++) {
		if (argv[i][0] == '-') {
			error = TRUE;
			continue;
		}
		if (!gmtmggpath_func (path, argv[i]))
			printf ("%s\n", path);
		else
			fprintf(stderr, "gmtpath: File %s not found\n", argv[i]);
	}
	if (error || argc == 1) {
		fprintf (stderr, "usage: gmtpath leg1 leg2 leg3 ...\n");
		exit (EXIT_FAILURE);
	}

	gmtmgg_end ();
	GMT_end (argc, argv);

	exit (EXIT_SUCCESS);
}
