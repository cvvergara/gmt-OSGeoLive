/*--------------------------------------------------------------------
 *	$Id: mgd77path.c,v 1.28 2011/07/11 19:22:04 guru Exp $
 *
 *    Copyright (c) 2004-2011 by P. Wessel
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/*
 * mgd77path accepts MGD77 cruise names and returns the full system
 * path to the file(s).
 *
 * Author:	Paul Wessel
 * Date:	26-AUG-2004
 * Version:	1.0 Based on the old gmtpath.c
 *
 *
 */
 
#include "mgd77.h"

int main (int argc, char **argv)
{
	GMT_LONG i, n_cruises = 0, n_paths;
	
	GMT_LONG error = FALSE, dirs = FALSE, strip = FALSE;
	
	char path[BUFSIZ], **list = NULL;
	
	struct MGD77_CONTROL M;

	argc = (int)GMT_begin (argc, argv);		/* Initialize GMT Machinery */
	
	MGD77_Init (&M);			/* Initialize MGD77 Machinery */

	for (i = 1; !error && i < argc; i++) {	/* Process each cruise */
		if (argv[i][0] != '-') continue;
		switch (argv[i][1]) {									
			case 'V':
			case '\0':
				error += GMT_parse_common_options (argv[i], NULL, NULL, NULL, NULL);
				break;
		
			case 'I':
				MGD77_Process_Ignore (argv[i][1], &argv[i][2]);
				break;
		
			case 'D':	/* Show list of directories with MGD77 files */
				dirs = TRUE;
				break;

			case 'P':	/* Show list of paths to MGD77 files */
				dirs = FALSE;
				if (argv[i][2] == '-') strip = TRUE;
				break;

			default:		/* Options not recognized */
				error = TRUE;
				break;
		}
	}

	if (GMT_give_synopsis_and_exit || argc == 1) {	/* Display usage */
		fprintf(stderr,"mgd77path %s - Return paths to MGD77 cruises and directories\n\n", MGD77_VERSION);
		fprintf(stderr,"usage: mgd77path <cruise(s)> -D|P[-] [-I<code>] [-V]\n\n");
         
		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);
              
		MGD77_Cruise_Explain ();
		fprintf(stderr,"	-D Instead of full cruise paths just list all directories with MGD77 files\n");
		fprintf(stderr,"	-P List full cruise paths [Default].  Append - to only get cruise names\n");
		fprintf(stderr,"	OPTIONS:\n\n");
		fprintf(stderr,"	-I Ignore certain data file formats from consideration. Append combination of act to ignore\n");
		fprintf(stderr,"	   (a) MGD77 ASCII, (c) MGD77+ netCDF, or (t) plain table files. [Default ignores none]\n");
		fprintf(stderr,"	-V verbose, report number of cruises returned\n");
		exit (EXIT_FAILURE);
	}

	if (dirs) {	/* Just list the current active MGD77 data directories and exit */
		printf ("# Currently, your $MGD77_HOME is set to: %s\n", M.MGD77_HOME);
		printf ("# $MGD77_HOME/mgd77_paths.txt contains these directories:\n");
		for (i = 0; i < M.n_MGD77_paths; i++) printf ("%s\n", M.MGD77_datadir[i]);
		exit (EXIT_SUCCESS);
	}

	n_paths = MGD77_Path_Expand (&M, argv, argc, &list);	/* Get list of requested IDs */

	if (n_paths == 0) {
		fprintf (stderr, "%s: No cruises found\n", GMT_program);
		exit (EXIT_FAILURE);
	}
	
	for (i = 0; i < n_paths; i++) {		/* Process each ID */
 		if (MGD77_Get_Path (path, list[i], &M))
   			fprintf (stderr, "%s : Cannot find cruise %s\n", GMT_program, list[i]);
		else if (strip) {
			printf ("%s\n", list[i]);
			n_cruises++;
		}
		else {
			printf ("%s\n", path);
			n_cruises++;
		}
	}
	
	if (gmtdefs.verbose) fprintf (stderr, "%s: Found %ld cruises\n", GMT_program, n_cruises);
	
	MGD77_Path_Free ((int)n_paths, list);
	MGD77_end (&M);
	GMT_end (argc, argv);
	
	exit (EXIT_SUCCESS);
}
